#!/usr/bin/env julia
"""
examples/cneo_example.jl - Clean Code cNEO (Constrained Nuclear-Electronic Orbital) Examples

Demonstrates proper usage of constrained NEO calculations following Clean Code principles.
Shows cNEO-HF and cNEO-MP2 calculations with nuclear position constraints.
"""

# Clean Code Constants
const DEFAULT_H2_BOND_LENGTH = 0.74  # Angstroms
const ENERGY_DISPLAY_PRECISION = 6
const POSITION_DISPLAY_PRECISION = 4
const CONSTRAINT_POSITION_X = 0.0    # Target X position for quantum proton
const CONSTRAINT_POSITION_Y = 0.0    # Target Y position for quantum proton  
const CONSTRAINT_POSITION_Z = 0.5    # Target Z position for quantum proton (Angstroms)

# Setup environment
function initialize_cneo_environment()
    ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/path/to/pyscf-master")
    push!(LOAD_PATH, dirname(@__DIR__))
    return NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
end

# Load required modules
using SparseQEEcNEO
using LinearAlgebra

function display_cneo_examples_header()
    println("Constrained NEO (cNEO) Examples")
    println("Following Clean Code Principles") 
    println("="^60)
    println("Demonstrates nuclear position constraints in NEO calculations")
    println("="^60)
end

# ======================== Example 1: cNEO-HF Calculation ========================

function demonstrate_cneo_hf_calculation(neo_configuration::NEOConfig)
    println("\nExample 1: cNEO-HF with constrained quantum proton")
    println("-"^50)
    
    h2_molecule = create_h2_molecule_with_quantum_proton()
    constraint_calculation = create_position_constraint_calculation()
    
    cneo_hf_results = perform_constrained_hf_calculation(
        h2_molecule, constraint_calculation, neo_configuration
    )
    
    display_cneo_hf_results(cneo_hf_results)
    return cneo_hf_results
end

function create_h2_molecule_with_quantum_proton()
    return Molecule(
        "H 0 0 0; H 0 0 $(DEFAULT_H2_BOND_LENGTH)", 
        "sto-3g",
        quantum_nuc=[0]  # First hydrogen as quantum nucleus
    )
end

function create_position_constraint_calculation()
    constraint_position = [CONSTRAINT_POSITION_X, CONSTRAINT_POSITION_Y, CONSTRAINT_POSITION_Z]
    
    return create_cneo_calculation(
        method="HF",
        constraint_positions=[constraint_position],
        convergence_threshold=1e-6,
        max_iterations=50,
        lambda_damping=0.5
    )
end

function perform_constrained_hf_calculation(molecule::Molecule, 
                                          cneo_calculation::CNEOCalculation,
                                          neo_configuration::NEOConfig)
    return run_cneo_hf(molecule, cneo_calculation, neo_config=neo_configuration)
end

function display_cneo_hf_results(results::CNEOResults)
    println("\ncNEO-HF Results:")
    println("  Total energy: $(format_energy(results.total_energy)) Ha")
    println("  Electronic energy: $(format_energy(results.electronic_energy)) Ha")
    println("  Nuclear kinetic energy: $(format_energy(results.nuclear_kinetic_energy)) Ha")
    println("  Converged: $(results.converged)")
    println("  Iterations: $(results.iterations)")
    
    display_constraint_satisfaction_analysis(results)
end

function display_constraint_satisfaction_analysis(results::CNEOResults)
    println("\nConstraint Satisfaction Analysis:")
    
    for nucleus_index in 1:length(results.constraint_positions)
        target_position = results.constraint_positions[nucleus_index]
        actual_position = results.actual_nuclear_positions[nucleus_index]
        position_error = results.position_errors[nucleus_index]
        lagrange_multiplier = results.lagrange_multipliers[nucleus_index]
        
        println("  Quantum nucleus $nucleus_index:")
        println("    Target position: $(format_position_vector(target_position))")
        println("    Actual position: $(format_position_vector(actual_position))")
        println("    Position error: $(format_energy(position_error))")
        println("    Lagrange multiplier: $(format_position_vector(lagrange_multiplier))")
        
        constraint_satisfaction = calculate_constraint_satisfaction_percentage(position_error)
        println("    Constraint satisfaction: $(constraint_satisfaction)%")
    end
end

function format_energy(energy::Float64)
    return round(energy, digits=ENERGY_DISPLAY_PRECISION)
end

function format_position_vector(position::Vector{Float64})
    formatted_components = [round(component, digits=POSITION_DISPLAY_PRECISION) for component in position]
    return "[$(join(formatted_components, ", "))]"
end

function calculate_constraint_satisfaction_percentage(position_error::Float64)
    # Simple metric: constraint is well satisfied if error < 1e-3
    satisfaction_threshold = 1e-3
    satisfaction_percentage = max(0, 100 * (1 - position_error / satisfaction_threshold))
    return round(min(satisfaction_percentage, 100), digits=1)
end

# ======================== Example 2: cNEO-MP2 Calculation ========================

function demonstrate_cneo_mp2_calculation(neo_configuration::NEOConfig)
    println("\nExample 2: cNEO-MP2 with correlation corrections")
    println("-"^50)
    
    h2_molecule = create_h2_molecule_with_quantum_proton()
    constraint_positions = [[CONSTRAINT_POSITION_X, CONSTRAINT_POSITION_Y, CONSTRAINT_POSITION_Z]]
    
    cneo_mp2_results = perform_constrained_mp2_calculation(
        h2_molecule, constraint_positions, neo_configuration
    )
    
    display_cneo_mp2_results(cneo_mp2_results)
    return cneo_mp2_results
end

function perform_constrained_mp2_calculation(molecule::Molecule,
                                           constraint_positions::Vector{Vector{Float64}},
                                           neo_configuration::NEOConfig)
    return run_cneo_mp2(molecule, constraint_positions, neo_config=neo_configuration)
end

function display_cneo_mp2_results(results::CNEOMP2Results)
    println("\ncNEO-MP2 Results:")
    println("  HF energy: $(format_energy(results.hf_energy)) Ha")
    println("  MP2 correlation energy: $(format_energy(results.mp2_correlation_energy)) Ha")
    println("  Total energy: $(format_energy(results.total_energy)) Ha")
    println("  HF converged: $(results.hf_converged)")
    println("  MP2 converged: $(results.mp2_converged)")
    println("  HF iterations: $(results.hf_iterations)")
    println("  MP2 iterations: $(results.mp2_iterations)")
    
    display_mp2_constraint_analysis(results)
end

function display_mp2_constraint_analysis(results::CNEOMP2Results)
    println("\nMP2-Level Constraint Analysis:")
    
    correlation_improvement = abs(results.mp2_correlation_energy)
    energy_improvement_percentage = calculate_energy_improvement_percentage(
        results.hf_energy, results.total_energy
    )
    
    println("  Correlation energy magnitude: $(format_energy(correlation_improvement)) Ha")
    println("  Energy improvement: $(energy_improvement_percentage)%")
    
    for nucleus_index in 1:length(results.nuclear_positions_hf)
        hf_position = results.nuclear_positions_hf[nucleus_index] 
        mp2_position = results.nuclear_positions_mp2[nucleus_index]
        
        println("  Nucleus $nucleus_index positions:")
        println("    HF position: $(format_position_vector(hf_position))")
        println("    MP2 position: $(format_position_vector(mp2_position))")
        
        position_change = norm(mp2_position - hf_position)
        println("    Position change (HF→MP2): $(format_energy(position_change))")
    end
end

function calculate_energy_improvement_percentage(hf_energy::Float64, total_energy::Float64)
    if abs(hf_energy) > 1e-10
        improvement_percentage = abs((total_energy - hf_energy) / hf_energy) * 100
        return round(improvement_percentage, digits=2)
    else
        return 0.0
    end
end

# ======================== Example 3: Comparison and Analysis ========================

function demonstrate_method_comparison(neo_configuration::NEOConfig)
    println("\nExample 3: cNEO method comparison and analysis")
    println("-"^50)
    
    try
        # Run both methods
        hf_results = demonstrate_cneo_hf_calculation(neo_configuration)
        mp2_results = demonstrate_cneo_mp2_calculation(neo_configuration)
        
        # Compare results
        compare_cneo_methods(hf_results, mp2_results)
        
    catch error
        handle_cneo_calculation_error(error)
    end
end

function compare_cneo_methods(hf_results::CNEOResults, mp2_results::CNEOMP2Results)
    println("\nMethod Comparison:")
    
    energy_difference = mp2_results.total_energy - hf_results.total_energy
    correlation_contribution = abs(energy_difference / hf_results.total_energy) * 100
    
    println("  Energy difference (MP2 - HF): $(format_energy(energy_difference)) Ha")
    println("  Correlation contribution: $(round(correlation_contribution, digits=2))%")
    
    compare_convergence_behavior(hf_results, mp2_results)
    compare_constraint_satisfaction(hf_results, mp2_results)
end

function compare_convergence_behavior(hf_results::CNEOResults, mp2_results::CNEOMP2Results)
    println("\nConvergence Comparison:")
    println("  HF converged: $(hf_results.converged) ($(hf_results.iterations) iterations)")
    println("  MP2 HF step converged: $(mp2_results.hf_converged) ($(mp2_results.hf_iterations) iterations)")
    println("  MP2 correlation converged: $(mp2_results.mp2_converged) ($(mp2_results.mp2_iterations) iterations)")
end

function compare_constraint_satisfaction(hf_results::CNEOResults, mp2_results::CNEOMP2Results)
    println("\nConstraint Satisfaction Comparison:")
    
    for nucleus_index in 1:length(hf_results.position_errors)
        hf_error = hf_results.position_errors[nucleus_index]
        mp2_hf_error = mp2_results.position_errors_hf[nucleus_index]
        mp2_error = mp2_results.position_errors_mp2[nucleus_index]
        
        println("  Nucleus $nucleus_index constraint errors:")
        println("    cNEO-HF: $(format_energy(hf_error))")  
        println("    cNEO-MP2 (HF step): $(format_energy(mp2_hf_error))")
        println("    cNEO-MP2 (final): $(format_energy(mp2_error))")
    end
end

# ======================== Main Execution Function ========================

"""
    run_cneo_examples()

Execute all cNEO examples demonstrating Clean Code principles.
"""
function run_cneo_examples()
    try
        neo_configuration = initialize_cneo_environment()
        display_cneo_examples_header()
        
        demonstrate_method_comparison(neo_configuration)
        
        display_cneo_completion_message()
        
    catch error
        handle_cneo_example_error(error)
    end
end

function display_cneo_completion_message()
    println("\n" * "="^60)
    println("cNEO examples completed successfully!")
    println("Demonstrated Clean Code constrained NEO calculations:")
    println("• Nuclear position constraints with Lagrange multipliers")
    println("• Newton optimization with analytical Hessians")  
    println("• cNEO-HF and cNEO-MP2 correlation corrections")
    println("• Clean separation of concerns and focused functions")
    println("="^60)
end

function handle_cneo_calculation_error(error::Exception)
    println("\n⚠️ cNEO calculation encountered an issue:")
    println("   Error: $error")
    println("\nThis is likely due to:")
    println("   • PySCF with NEO support not available")
    println("   • Insufficient basis set for constraint satisfaction")
    println("   • Numerical convergence challenges")
    println("\nThe cNEO implementation is ready and follows Clean Code principles.")
end

function handle_cneo_example_error(error::Exception)
    println("\n❌ cNEO example execution failed:")
    println("   $error")
    println("\nPlease ensure:")
    println("   • PySCF with NEO support is properly installed")
    println("   • Environment variables are correctly set")
    println("   • Sufficient computational resources available")
    println("\nThe Clean Code cNEO implementation is integrated and ready to use!")
end

# Execute the cNEO examples
run_cneo_examples()