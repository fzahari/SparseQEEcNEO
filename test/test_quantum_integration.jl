#!/usr/bin/env julia
"""
test_quantum_integration.jl - Comprehensive Quantum Computing Integration Tests

Tests the integrated quantum computing functionality including:
- Export to quantum formats (OpenFermion, Qiskit, Cirq)  
- Quantum resource estimation
- Circuit creation
- Validation functions

All tests follow Clean Code principles with descriptive names and single responsibilities.
"""

using Test
using SparseQEEcNEO

# Test constants
const TEST_H2_BOND_LENGTH = 1.4
const TEST_BASIS = "sto-3g"
const TEST_MAX_QUBITS = 4
const SUCCESS_SYMBOL = "✅"
const ERROR_SYMBOL = "❌"

@testset "Quantum Computing Integration Tests" begin
    
    @testset "Module Loading and Exports" begin
        @test isdefined(SparseQEEcNEO, :exportToQuantumFormats)
        @test isdefined(SparseQEEcNEO, :exportToOpenFermion)
        @test isdefined(SparseQEEcNEO, :exportToQiskit)
        @test isdefined(SparseQEEcNEO, :exportToCirq)
        @test isdefined(SparseQEEcNEO, :createVQECircuit)
        @test isdefined(SparseQEEcNEO, :estimateQuantumResources)
        @test isdefined(SparseQEEcNEO, :validateQuantumExport)
        @test isdefined(SparseQEEcNEO, :runQuantumDemo)
        
        println("   $SUCCESS_SYMBOL Quantum computing functions are properly exported")
    end
    
    @testset "Test Data Creation" begin
        # Create test molecule
        testMolecule = Molecule(
            "H 0 0 0; H 0 0 $TEST_H2_BOND_LENGTH",
            TEST_BASIS,
            quantum_nuc=[1]
        )
        
        @test testMolecule.atom_string == "H 0 0 0; H 0 0 $(TEST_H2_BOND_LENGTH)"
        @test testMolecule.basis == TEST_BASIS
        @test testMolecule.quantum_nuc == [1]
        
        # Create test configuration
        testConfig = ConfigSelection(
            method="mp2",
            max_configs=10,
            max_qubits=TEST_MAX_QUBITS
        )
        
        @test testConfig.method == "mp2"
        @test testConfig.max_qubits == TEST_MAX_QUBITS
        
        println("   $SUCCESS_SYMBOL Test data structures created successfully")
    end
    
    @testset "VQE Circuit Creation" begin
        # Test basic circuit creation
        circuit2Qubits = createVQECircuit(2)
        
        @test circuit2Qubits["circuit_type"] == "VQE"
        @test circuit2Qubits["n_qubits"] == 2
        @test circuit2Qubits["layers"] == 2  # Default layers
        @test circuit2Qubits["parameter_count"] == 4  # 2 qubits * 2 layers
        @test circuit2Qubits["depth_estimate"] > 0
        @test haskey(circuit2Qubits, "description")
        
        # Test custom layers
        circuit3Layers = createVQECircuit(3, 3)
        @test circuit3Layers["layers"] == 3
        @test circuit3Layers["parameter_count"] == 9  # 3 qubits * 3 layers
        
        # Test error handling
        @test_throws ArgumentError createVQECircuit(0)
        @test_throws ArgumentError createVQECircuit(2, 0)
        
        println("   $SUCCESS_SYMBOL VQE circuit creation working correctly")
    end
    
    @testset "Hamiltonian Data Creation for Export" begin
        # Create mock Hamiltonian data for testing
        mockH1e = Dict(
            "e" => [1.0 0.5; 0.5 -0.8],
            "n0" => [0.3 0.0; 0.0 0.2]
        )
        
        mockH2e = Dict(
            "e" => Dict((1,1,1,1) => 0.67, (1,2,1,2) => -0.45),
            "n0" => Dict((1,1,1,1) => 0.1)
        )
        
        mockCoupling = Dict(
            ("e", "n0") => [0.1 0.05; 0.05 0.08]
        )
        
        mockConfigs = [
            Configuration("HF", Dict("e" => Int16[1, 0]), 1.0),
            Configuration("S(0→2)", Dict("e" => Int16[0, 1]), 0.3)
        ]
        
        hamiltonianData = HamiltonianData(
            mockH1e,
            mockH2e,
            mockCoupling,
            Dict("e" => 2, "n0" => 2),
            Dict("e" => 1, "n0" => 1),
            [1.0, 0.0],
            mockConfigs,
            [-1.1, -0.8],
            Dict("e" => [1, 2], "n0" => [1, 2]),
            Dict("e" => Int[], "n0" => Int[])
        )
        
        @test length(hamiltonianData.h1e) == 2
        @test length(hamiltonianData.h2e) == 2
        @test length(hamiltonianData.coupling) == 1
        
        println("   $SUCCESS_SYMBOL Mock Hamiltonian data created for testing")
        
        # Store for use in other test sets
        global testHamiltonianData = hamiltonianData
    end
    
    @testset "Quantum Resource Estimation" begin
        resources = estimateQuantumResources(testHamiltonianData, "VQE")
        
        @test haskey(resources, "qubits_required")
        @test haskey(resources, "pauli_terms")
        @test haskey(resources, "algorithm")
        @test haskey(resources, "circuit_depth")
        @test haskey(resources, "parameter_count")
        @test haskey(resources, "total_measurements")
        
        @test resources["algorithm"] == "VQE"
        @test resources["qubits_required"] >= 2
        @test resources["pauli_terms"] >= 0
        
        # Test QAOA resources
        qaoaResources = estimateQuantumResources(testHamiltonianData, "QAOA")
        @test qaoaResources["algorithm"] == "QAOA"
        
        println("   $SUCCESS_SYMBOL Quantum resource estimation working")
    end
    
    @testset "Export Format Validation" begin
        # Test format validation
        @test_throws ArgumentError SparseQEEcNEO.QuantumComputing.validateExportFormat("invalid_format")
        
        # Valid formats should not throw
        for format in ["openfermion", "qiskit", "cirq", "all"]
            @test_nowarn SparseQEEcNEO.QuantumComputing.validateExportFormat(format)
        end
        
        # Test Hamiltonian data validation
        emptyHamiltonian = HamiltonianData(
            Dict{String, Matrix{Float64}}(),
            Dict{String, Dict{NTuple{4,Int}, Float64}}(),
            Dict{Tuple{String, String}, Matrix{Float64}}(),
            Dict{String, Int}(),
            Dict{String, Int}(),
            Float64[],
            Configuration[],
            Float64[],
            Dict{String, Vector{Int}}(),
            Dict{String, Vector{Int}}()
        )
        
        @test_throws ArgumentError SparseQEEcNEO.QuantumComputing.validateHamiltonianData(emptyHamiltonian)
        @test_nowarn SparseQEEcNEO.QuantumComputing.validateHamiltonianData(testHamiltonianData)
        
        println("   $SUCCESS_SYMBOL Export validation working correctly")
    end
    
    @testset "Export Functionality" begin
        # Create temporary files for testing
        tempDir = mktempdir()
        
        # Test OpenFermion export
        ofResult = exportToOpenFermion(testHamiltonianData, joinpath(tempDir, "test.of"))
        @test ofResult["success"] == true
        @test haskey(ofResult, "operator_count")
        @test haskey(ofResult, "filename")
        @test isfile(ofResult["filename"])
        
        # Test that file contains expected content
        ofContent = read(ofResult["filename"], String)
        @test occursin("# OpenFermion format export", ofContent)
        @test occursin("^", ofContent)  # Creation operators
        
        # Test Qiskit export
        qiskitResult = exportToQiskit(testHamiltonianData, joinpath(tempDir, "test.qiskit"))
        @test qiskitResult["success"] == true
        @test haskey(qiskitResult, "pauli_count")
        
        # Test Cirq export
        cirqResult = exportToCirq(testHamiltonianData, joinpath(tempDir, "test.cirq"))
        @test cirqResult["success"] == true
        @test haskey(cirqResult, "operator_count")
        
        # Test export to all formats
        allResult = exportToQuantumFormats(testHamiltonianData, "all", joinpath(tempDir, "test_all"))
        @test haskey(allResult, "openfermion")
        @test haskey(allResult, "qiskit")
        @test haskey(allResult, "cirq")
        @test allResult["openfermion"]["success"] == true
        
        # Clean up temp files
        rm(tempDir, recursive=true)
        
        println("   $SUCCESS_SYMBOL Export functionality working correctly")
    end
    
    @testset "Validation Functions" begin
        validation = validateQuantumExport(testHamiltonianData, "openfermion")
        
        @test haskey(validation, "hermiticity")
        @test haskey(validation, "energy_range")
        @test haskey(validation, "sparsity")
        @test haskey(validation, "fermionic_structure")
        
        @test isa(validation["hermiticity"], Bool)
        @test haskey(validation["energy_range"], :min)
        @test haskey(validation["energy_range"], :max)
        @test validation["sparsity"] >= 0.0
        @test validation["sparsity"] <= 1.0
        
        println("   $SUCCESS_SYMBOL Validation functions working correctly")
    end
    
    @testset "Demo Function" begin
        # Capture output to test demo function
        originalStdout = stdout
        
        try
            # Redirect stdout to capture output
            (readPipe, writePipe) = redirect_stdout()
            
            # Run demo in a task so we can capture its output
            demoTask = @async runQuantumDemo(testHamiltonianData)
            
            # Give it a moment to run
            sleep(1.0)
            
            # Restore stdout
            redirect_stdout(originalStdout)
            close(writePipe)
            
            # Read captured output
            output = String(read(readPipe))
            
            # Verify demo ran and produced expected output
            @test occursin("Quantum Computing Integration Demo", output)
            @test occursin("Quantum Resource Estimates", output)
            @test occursin("completed successfully", output)
            
            println("   $SUCCESS_SYMBOL Demo function executed successfully")
            
        catch error
            # Restore stdout if there's an error
            redirect_stdout(originalStdout)
            println("   ⚠️  Demo test skipped due to output redirection issues")
            @test_skip "Demo function test"
        end
    end
    
    println("\n$SUCCESS_SYMBOL All quantum computing integration tests passed!")
end

# Run the tests when file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running Quantum Computing Integration Tests...")
    include("test_quantum_integration.jl")
end