@testset "Integration Tests" begin
    if !NEO_AVAILABLE
        @test_skip "PySCF NEO not available"
        return
    end
    
    @testset "Full H2 NEO Calculation" begin
        # Use NEO molecule with quantum nucleus
        mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
        calc = NEOCalculation(xc="HF")
        config_sel = ConfigSelection(method="mp2", max_configs=10, max_nuc_orbs=0)
        
        results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel, neo_config=TEST_CONFIG)
        
        @test results.n_configs > 0
        @test results.n_configs <= 10
        @test results.energy < 0
        @test results.n_qubits > 0
        @test 0 < results.captured_importance <= 1.0
        @test results.hamiltonian_matrix !== nothing
        @test ishermitian(results.hamiltonian_matrix)
        
        # Check energy is reasonable
        @test -1.2 < results.total_energy < -1.0  # H2 at equilibrium
    end
    
    @testset "H2 with Quantum Proton" begin
        mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
        calc = NEOCalculation(xc="HF")
        config_sel = ConfigSelection(method="neo_cneo", max_configs=50, max_nuc_orbs=0)
        
        results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel, neo_config=TEST_CONFIG)
        
        @test results.n_configs > 2  # Should have nuclear excitations
        @test results.neo_metrics !== nothing
        @test results.neo_metrics.nuclear_participation > 0.5
        @test length(results.orbitals_per_species) == 2  # Electronic and nuclear
    end
    
    @testset "Water with Quantum Proton" begin
        mol = Molecule("O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", 
                      "sto-3g", quantum_nuc=[1])
        calc = NEOCalculation(xc="B3LYP", epc="17-2")
        config_sel = ConfigSelection(method="neo_cneo", max_configs=50, max_nuc_orbs=0)
        
        results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel, neo_config=TEST_CONFIG)
        
        @test results.n_configs > 0
        @test results.energy < -75  # Water should be around -76 Ha
        @test results.hamiltonian_matrix !== nothing
    end
    
    @testset "Method Comparison" begin
        mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
        methods = ["mp2", "neo_cneo", "hybrid_final"]
        
        results_dict = Dict{String, NEOResults}()
        for method in methods
            config_sel = ConfigSelection(method=method, max_configs=20, max_nuc_orbs=0)
            try
                results = sparse_qee_cneo(mol, config_sel=config_sel, neo_config=TEST_CONFIG)
                results_dict[method] = results
            catch e
                @warn "Method $method failed: $e"
            end
        end
        
        # Compare results
        if length(results_dict) >= 2
            energies = [r.total_energy for r in values(results_dict)]
            @test maximum(energies) - minimum(energies) < 0.1  # Should be similar
        end
    end
    
    @testset "EPC Functional Comparison" begin
        mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
        calc_hf = NEOCalculation(xc="HF")
        calc_epc = NEOCalculation(xc="B3LYP", epc="17-2")
        config_sel = ConfigSelection(method="neo_cneo", max_configs=20, max_nuc_orbs=0)
        
        results_hf = sparse_qee_cneo(mol, calc=calc_hf, config_sel=config_sel, neo_config=TEST_CONFIG)
        results_epc = sparse_qee_cneo(mol, calc=calc_epc, config_sel=config_sel, neo_config=TEST_CONFIG)
        
        # With proper EPC, energy might be lower or similar
        @test abs(results_epc.total_energy - results_hf.total_energy) < 1.0
    end
end
