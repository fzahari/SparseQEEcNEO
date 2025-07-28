@testset "QEE Methods Tests" begin
    @testset "Qubit Mapping Optimization" begin
        configs = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0,0,0,0,0]), 1.0),
            Configuration("S(1→5)", Dict("e" => Int16[0,1,0,0,1,0,0,0]), 0.5),
            Configuration("S(2→7)", Dict("e" => Int16[1,0,0,0,0,0,1,0]), 0.3)
        ]
        
        orbital_to_qubit, n_qubits = SparseQEEcNEO.QEEMethods.optimize_qubit_mapping(configs)
        
        # Should map only used orbitals
        @test n_qubits == 4  # Orbitals 1,2,5,7 are used
        @test haskey(orbital_to_qubit, 1)
        @test haskey(orbital_to_qubit, 2)
        @test haskey(orbital_to_qubit, 5)
        @test haskey(orbital_to_qubit, 7)
        @test !haskey(orbital_to_qubit, 3)  # Not used
    end
    
    @testset "Resource Estimation - VQE" begin
        configs = [Configuration("C$i", Dict("e" => Int16[1,1,0,0]), rand()) 
                  for i in 1:50]
        
        resources = SparseQEEcNEO.QEEMethods.calculate_resource_requirements(
            configs, "vqe"
        )
        
        @test resources.n_configs == 50
        @test resources.n_qubits > 0
        @test resources.n_parameters > 0
        @test resources.circuit_depth > 0
        @test resources.n_measurements > 0
        @test resources.method == "vqe"
    end
    
    @testset "Resource Estimation - QPE" begin
        configs = [Configuration("C$i", Dict("e" => Int16[1,1,0,0]), rand()) 
                  for i in 1:20]
        
        resources = SparseQEEcNEO.QEEMethods.calculate_resource_requirements(
            configs, "qpe"
        )
        
        @test resources.method == "qpe"
        @test resources.n_parameters == 0  # QPE has no variational parameters
        @test resources.circuit_depth > resources.n_qubits  # QPE needs deep circuits
    end
    
    @testset "Initial State Generation" begin
        ref_config = Configuration("HF", Dict("e" => Int16[1,1,0,0,1,0,0,0]), 1.0)
        orbital_to_qubit = Dict(1=>1, 2=>2, 3=>3, 4=>4, 5=>5)
        
        circuit = SparseQEEcNEO.QEEMethods.generate_initial_state(
            ref_config, orbital_to_qubit
        )
        
        @test length(circuit) == 3  # Three occupied orbitals
        @test "X(1)" in circuit
        @test "X(2)" in circuit
        @test "X(5)" in circuit
    end
    
    @testset "Entanglement Structure Analysis" begin
        configs = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0),
            Configuration("S(1→3)", Dict("e" => Int16[0,1,1,0]), 0.5),
            Configuration("S(2→4)", Dict("e" => Int16[1,0,0,1]), 0.3),
            Configuration("D(12→34)", Dict("e" => Int16[0,0,1,1]), 0.2)
        ]
        
        entanglement = SparseQEEcNEO.QEEMethods.analyze_entanglement_structure(configs)
        
        @test entanglement.n_active_orbitals == 4
        @test entanglement.total_connections > 0
        @test 0 <= entanglement.connectivity_fraction <= 1
        @test size(entanglement.connectivity_matrix, 1) == 4
    end
    
    @testset "Quantum Advantage Estimation" begin
        # Small system - no advantage
        advantage_small = SparseQEEcNEO.QEEMethods.estimate_quantum_advantage(10, 4)
        
        @test advantage_small.classical_operations > 0
        @test advantage_small.quantum_operations > 0
        @test advantage_small.classical_memory_gb > 0
        
        # Large system - potential advantage
        advantage_large = SparseQEEcNEO.QEEMethods.estimate_quantum_advantage(1000, 20)
        
        @test advantage_large.advantage_ratio > advantage_small.advantage_ratio
        @test advantage_large.classical_memory_gb > advantage_small.classical_memory_gb
    end
    
    @testset "Circuit Depth Estimation" begin
        instructions = [
            "X(1)",
            "X(2)",
            "CNOT(1,2)",
            "RZ(0.5,1)",
            "CNOT(2,3)",
            "RY(0.3,2)"
        ]
        
        depth = SparseQEEcNEO.QEEMethods.estimate_circuit_depth(instructions)
        
        @test depth >= 3  # Minimum based on dependencies
        @test depth <= length(instructions)  # Maximum if all sequential
    end
end
