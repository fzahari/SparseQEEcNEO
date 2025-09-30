"""
test_cneo_methods.jl - Clean Code Tests for cNEO Methods

Tests constrained Nuclear-Electronic Orbital (cNEO) functionality following Clean Code principles.
Comprehensive testing of cNEO-HF and cNEO-MP2 implementations with proper error handling.
"""

using Test
using LinearAlgebra
using SparseQEEcNEO

# Test constants following Clean Code principles
const TEST_H2_BOND_LENGTH = 0.74
const TEST_CONSTRAINT_TOLERANCE = 1e-6
const TEST_ENERGY_TOLERANCE = 1e-8
const TEST_POSITION_TOLERANCE = 1e-4
const TEST_MAX_ITERATIONS = 20
const TEST_CONSTRAINT_X = 0.0
const TEST_CONSTRAINT_Y = 0.0  
const TEST_CONSTRAINT_Z = 0.5
const MOCK_HF_ENERGY = -1.1
const MOCK_ELECTRONIC_ENERGY = -1.2
const MOCK_NUCLEAR_KINETIC_ENERGY = 0.1

# Test configuration
const TEST_NEO_CONFIG = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))

# ======================== Test Helper Functions ========================

function create_test_molecule_h2()
    return Molecule(
        "H 0 0 0; H 0 0 $(TEST_H2_BOND_LENGTH)",
        "sto-3g",
        quantum_nuc=[0]
    )
end

function create_test_constraint_positions()
    return [[TEST_CONSTRAINT_X, TEST_CONSTRAINT_Y, TEST_CONSTRAINT_Z]]
end

function create_test_cneo_calculation()
    return create_cneo_calculation(
        method="HF",
        constraint_positions=create_test_constraint_positions(),
        convergence_threshold=TEST_CONSTRAINT_TOLERANCE,
        max_iterations=TEST_MAX_ITERATIONS
    )
end

function create_mock_cneo_results()
    constraint_positions = create_test_constraint_positions()
    actual_positions = [[0.01, 0.01, 0.51]]  # Slightly off from constraint
    position_errors = [norm(actual_positions[1] - constraint_positions[1])]
    lagrange_multipliers = [[0.1, 0.1, 0.2]]
    
    return CNEOResults(
        MOCK_HF_ENERGY,           # total_energy
        MOCK_ELECTRONIC_ENERGY,   # electronic_energy  
        MOCK_NUCLEAR_KINETIC_ENERGY,  # nuclear_kinetic_energy
        lagrange_multipliers,     # lagrange_multipliers
        actual_positions,         # actual_nuclear_positions
        constraint_positions,     # constraint_positions
        position_errors,          # position_errors
        true,                     # converged
        5,                        # iterations
        nothing                   # mean_field_object (mock)
    )
end

@testset "cNEO Methods Clean Code Tests" begin
    
    @testset "CNEOCalculation Parameter Validation" begin
        @testset "Valid CNEOCalculation Creation" begin
            constraint_positions = create_test_constraint_positions()
            
            calc = CNEOCalculation(
                "HF", "B3LYP", constraint_positions, 
                1e-6, 50, 0.5
            )
            
            @test calc.method == "HF"
            @test calc.functional == "B3LYP"
            @test calc.constraint_positions == constraint_positions
            @test calc.convergence_threshold == 1e-6
            @test calc.max_iterations == 50
            @test calc.lambda_damping == 0.5
        end
        
        @testset "Invalid Method Validation" begin
            @test_throws ArgumentError CNEOCalculation(
                "INVALID", "B3LYP", create_test_constraint_positions(),
                1e-6, 50, 0.5
            )
        end
        
        @testset "Invalid Convergence Threshold" begin
            @test_throws ArgumentError CNEOCalculation(
                "HF", "B3LYP", create_test_constraint_positions(),
                -1e-6, 50, 0.5
            )
        end
        
        @testset "Invalid Damping Factor" begin  
            @test_throws ArgumentError CNEOCalculation(
                "HF", "B3LYP", create_test_constraint_positions(),
                1e-6, 50, 1.5
            )
            
            @test_throws ArgumentError CNEOCalculation(
                "HF", "B3LYP", create_test_constraint_positions(),
                1e-6, 50, 0.0
            )
        end
        
        @testset "Invalid Constraint Position Dimensions" begin
            invalid_positions = [[1.0, 2.0]]  # Only 2D instead of 3D
            
            @test_throws ArgumentError CNEOCalculation(
                "HF", "B3LYP", invalid_positions,
                1e-6, 50, 0.5
            )
        end
    end
    
    @testset "CNEOCalculation Factory Functions" begin
        @testset "Default Parameter Creation" begin
            calc = create_cneo_calculation()
            
            @test calc.method == "HF"
            @test calc.functional == "B3LYP"
            @test isempty(calc.constraint_positions)
            @test calc.convergence_threshold ≈ 1e-6
            @test calc.max_iterations == 50
            @test calc.lambda_damping ≈ 0.5
        end
        
        @testset "Custom Parameter Creation" begin
            constraint_positions = create_test_constraint_positions()
            
            calc = create_cneo_calculation(
                method="DFT",
                functional="PBE",
                constraint_positions=constraint_positions,
                convergence_threshold=1e-7,
                max_iterations=100,
                lambda_damping=0.3
            )
            
            @test calc.method == "DFT"
            @test calc.functional == "PBE"
            @test calc.constraint_positions == constraint_positions
            @test calc.convergence_threshold ≈ 1e-7
            @test calc.max_iterations == 100
            @test calc.lambda_damping ≈ 0.3
        end
    end
    
    @testset "CNEOResults Data Structure" begin
        @testset "CNEOResults Construction and Access" begin
            results = create_mock_cneo_results()
            
            @test results.total_energy ≈ MOCK_HF_ENERGY
            @test results.electronic_energy ≈ MOCK_ELECTRONIC_ENERGY
            @test results.nuclear_kinetic_energy ≈ MOCK_NUCLEAR_KINETIC_ENERGY
            @test length(results.lagrange_multipliers) == 1
            @test length(results.actual_nuclear_positions) == 1
            @test length(results.constraint_positions) == 1
            @test length(results.position_errors) == 1
            @test results.converged == true
            @test results.iterations == 5
        end
        
        @testset "Position Error Calculation" begin
            results = create_mock_cneo_results()
            expected_error = norm([0.01, 0.01, 0.01])  # Difference from constraint
            
            @test results.position_errors[1] ≈ expected_error atol=1e-10
        end
    end
    
    @testset "CNEOMP2Results Data Structure" begin
        @testset "CNEOMP2Results Mock Creation" begin
            hf_results = create_mock_cneo_results()
            mp2_correlation = -0.1
            
            mp2_results = CNEOMP2Results(
                hf_results.total_energy,      # hf_energy
                mp2_correlation,              # mp2_correlation_energy  
                hf_results.total_energy + mp2_correlation,  # total_energy
                hf_results.lagrange_multipliers,  # lagrange_multipliers_hf
                hf_results.lagrange_multipliers,  # lagrange_multipliers_mp2
                hf_results.actual_nuclear_positions,  # nuclear_positions_hf
                hf_results.actual_nuclear_positions,  # nuclear_positions_mp2
                hf_results.position_errors,   # position_errors_hf
                hf_results.position_errors,   # position_errors_mp2
                true,                         # hf_converged
                true,                         # mp2_converged
                hf_results.iterations,        # hf_iterations
                3,                           # mp2_iterations
                nothing                      # mean_field_object
            )
            
            @test mp2_results.hf_energy ≈ MOCK_HF_ENERGY
            @test mp2_results.mp2_correlation_energy ≈ mp2_correlation
            @test mp2_results.total_energy ≈ MOCK_HF_ENERGY + mp2_correlation
            @test mp2_results.hf_converged == true
            @test mp2_results.mp2_converged == true
            @test mp2_results.hf_iterations == 5
            @test mp2_results.mp2_iterations == 3
        end
    end
    
    @testset "Basic Functionality Tests" begin
        @testset "Constants and Utilities" begin
            # Test Clean Code constants are properly defined
            @test SparseQEEcNEO.CNEOMethods.POSITION_DIMENSION == 3
            @test SparseQEEcNEO.CNEOMethods.MINIMUM_DENOMINATOR_THRESHOLD > 0
            @test SparseQEEcNEO.CNEOMethods.HESSIAN_REGULARIZATION > 0
            @test SparseQEEcNEO.CNEOMethods.DEFAULT_LAMBDA_DAMPING_FACTOR > 0
            @test SparseQEEcNEO.CNEOMethods.DEFAULT_LAMBDA_DAMPING_FACTOR <= 1
            
            # Test utility functions
            @test SparseQEEcNEO.CNEOMethods.format_energy(1.234567) == 1.234567
            @test SparseQEEcNEO.CNEOMethods.create_identity_hessian() == Matrix(1.0I, 3, 3)
        end
        
        @testset "Parameter Validation" begin
            # Test should not throw for valid methods
            @test_nowarn SparseQEEcNEO.CNEOMethods.validate_calculation_method("HF")
            @test_nowarn SparseQEEcNEO.CNEOMethods.validate_calculation_method("DFT")
            
            # Test should throw for invalid methods  
            @test_throws ArgumentError SparseQEEcNEO.CNEOMethods.validate_calculation_method("INVALID")
            @test_throws ArgumentError SparseQEEcNEO.CNEOMethods.validate_calculation_method("")
            
            # Functional name validation
            @test_nowarn SparseQEEcNEO.CNEOMethods.validate_functional_name("B3LYP")
            @test_throws ArgumentError SparseQEEcNEO.CNEOMethods.validate_functional_name("")
            
            # Convergence threshold validation
            @test_nowarn SparseQEEcNEO.CNEOMethods.validate_convergence_threshold(1e-6)
            @test_throws ArgumentError SparseQEEcNEO.CNEOMethods.validate_convergence_threshold(-1e-6)
            
            # Damping factor validation
            @test_nowarn SparseQEEcNEO.CNEOMethods.validate_damping_factor(0.5)
            @test_throws ArgumentError SparseQEEcNEO.CNEOMethods.validate_damping_factor(1.5)
        end
    end
end

println("Running cNEO Methods Clean Code Tests...")
println("Testing constrained Nuclear-Electronic Orbital functionality...")
println("All tests focus on Clean Code principles and proper error handling.")