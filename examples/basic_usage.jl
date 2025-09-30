#!/usr/bin/env julia
"""
examples/basic_usage.jl - Clean Code examples of Sparse QEE-cNEO usage

Demonstrates proper usage patterns following Clean Code principles:
- Clear function names and structure
- Single responsibility principle
- Meaningful variable names
- Proper error handling
- Comprehensive documentation
"""

# Clean Code Constants
const DEFAULT_H2_BOND_LENGTH = 0.74  # Angstroms
const MAXIMUM_DISPLAYED_CONFIGURATIONS = 5
const ENERGY_DISPLAY_PRECISION = 6
const PERCENTAGE_DISPLAY_PRECISION = 1

# Setup environment
function initialize_environment()
    ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/path/to/pyscf-master")
    push!(LOAD_PATH, dirname(@__DIR__))
    return NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
end

# Load required modules
using SparseQEEcNEO
using LinearAlgebra

function display_examples_header()
    println("Sparse QEE-cNEO Basic Usage Examples")
    println("Following Clean Code Principles")
    println("="^60)
end

# Initialize
neo_configuration = initialize_environment()
display_examples_header()

# ======================== Example 1: Classical H2 Calculation ========================

function demonstrate_classical_h2_calculation(neo_configuration::NEOConfig)
    println("\nExample 1: H2 molecule (electronic only)")
    println("-"^40)
    
    h2_molecule = create_classical_h2_molecule()
    calculation_results = perform_classical_calculation(h2_molecule, neo_configuration)
    display_classical_results(calculation_results)
    
    return calculation_results
end

function create_classical_h2_molecule()
    return Molecule(
        "H 0 0 0; H 0 0 $(DEFAULT_H2_BOND_LENGTH)", 
        "sto-3g"
    )
end

function perform_classical_calculation(molecule::Molecule, neo_configuration::NEOConfig)
    return sparse_qee_cneo(molecule, neo_config=neo_configuration)
end

function display_classical_results(results::NEOResults)
    println("Classical H2 Results:")
    println("  Energy: $(format_energy(results.energy)) Ha")
    println("  MP2 correlation: $(format_energy(results.mp2_correlation)) Ha")
    println("  Total energy: $(format_energy(results.total_energy)) Ha")
    println("  Configurations: $(results.n_configs)")
    println("  Qubits needed: $(results.n_qubits)")
end

function format_energy(energy::Float64)
    return round(energy, digits=ENERGY_DISPLAY_PRECISION)
end

# ======================== Example 2: Quantum Nuclear H2 Calculation ========================

function demonstrate_quantum_nuclear_h2_calculation(neo_configuration::NEOConfig)
    println("\nExample 2: H2 with quantum proton")
    println("-"^40)
    
    quantum_h2_molecule = create_quantum_nuclear_h2_molecule()
    neo_calculation_parameters = create_neo_calculation_parameters()
    configuration_selection = create_neo_configuration_selection()
    
    quantum_results = perform_quantum_nuclear_calculation(
        quantum_h2_molecule, neo_calculation_parameters, 
        configuration_selection, neo_configuration
    )
    
    display_quantum_nuclear_results(quantum_results)
    return quantum_results
end

function create_quantum_nuclear_h2_molecule()
    return Molecule(
        "H 0 0 0; H 0 0 $(DEFAULT_H2_BOND_LENGTH)", 
        "sto-3g", 
        quantum_nuc=[0]  # First hydrogen as quantum nucleus
    )
end

function create_neo_calculation_parameters()
    return NEOCalculation(xc="B3LYP", epc="17-2")
end

function create_neo_configuration_selection()
    return ConfigSelection(
        method="neo_cneo", 
        max_configs=50, 
        max_nuc_orbs=0
    )
end

function perform_quantum_nuclear_calculation(molecule::Molecule, calculation::NEOCalculation,
                                           config_selection::ConfigSelection, neo_config::NEOConfig)
    return sparse_qee_cneo(
        molecule, 
        calc=calculation, 
        config_sel=config_selection, 
        neo_config=neo_config
    )
end

function display_quantum_nuclear_results(results::NEOResults)
    println("Quantum Nuclear H2 Results:")
    println("  Energy: $(format_energy(results.energy)) Ha")
    println("  Total energy (with EPC): $(format_energy(results.total_energy)) Ha")
    println("  Configurations: $(results.n_configs)")
end

function compare_calculation_energies(classical_results::NEOResults, quantum_results::NEOResults)
    energy_difference = quantum_results.energy - classical_results.energy
    println("\nEnergy Comparison:")
    println("  Classical H2: $(format_energy(classical_results.energy)) Ha")
    println("  Quantum nuclear H2: $(format_energy(quantum_results.energy)) Ha")
    println("  Difference (NEO - classical): $(format_energy(energy_difference)) Ha")
end

# ======================== Main Execution Function ========================

"""
    run_clean_code_examples()

Execute all Clean Code examples demonstrating Sparse QEE-cNEO functionality.
Follows Clean Code principles with clear separation of concerns.
"""
function run_clean_code_examples()
    try
        # Example 1: Classical H2 calculation
        classical_results = demonstrate_classical_h2_calculation(neo_configuration)
        
        # Example 2: Quantum nuclear H2 calculation
        quantum_results = demonstrate_quantum_nuclear_h2_calculation(neo_configuration)
        
        # Compare the two calculations
        compare_calculation_energies(classical_results, quantum_results)
        
        # Example 3: Method comparison
        demonstrate_configuration_method_comparison(quantum_results.configs[1].occupations, neo_configuration)
        
        # Example 4: Hamiltonian analysis
        demonstrate_hamiltonian_analysis(quantum_results)
        
        # Example 5: Data persistence
        demonstrate_results_saving(quantum_results)
        
        display_completion_message()
        
    catch error
        handle_example_error(error)
    end
end

function demonstrate_configuration_method_comparison(base_molecule_data::Dict, neo_configuration::NEOConfig)
    println("\nExample 3: Comparing configuration methods")
    println("-"^40)
    
    configuration_methods = ["mp2", "neo_cneo", "hybrid_final"]
    method_results = compare_configuration_methods(configuration_methods, neo_configuration)
    
    display_method_comparison_results(method_results)
end

function compare_configuration_methods(methods::Vector{String}, neo_configuration::NEOConfig)
    quantum_molecule = create_quantum_nuclear_h2_molecule()
    neo_calculation = create_neo_calculation_parameters()
    results_dict = Dict{String, NEOResults}()
    
    for method_name in methods
        config_selection = ConfigSelection(
            method=method_name, 
            max_configs=50, 
            max_nuc_orbs=0
        )
        
        method_results = perform_quantum_nuclear_calculation(
            quantum_molecule, neo_calculation, config_selection, neo_configuration
        )
        
        results_dict[method_name] = method_results
    end
    
    return results_dict
end

function display_method_comparison_results(method_results::Dict{String, NEOResults})
    for (method_name, results) in method_results
        importance_percentage = round(results.captured_importance * 100, digits=PERCENTAGE_DISPLAY_PRECISION)
        
        println("$method_name method:")
        println("  Configurations: $(results.n_configs)")
        println("  Importance captured: $importance_percentage%")
        println("  Qubits required: $(results.n_qubits)")
        println()
    end
end

function demonstrate_hamiltonian_analysis(calculation_results::NEOResults)
    println("Example 4: Working with Hamiltonians")
    println("-"^40)
    
    hamiltonian_matrix = calculation_results.hamiltonian_matrix
    hamiltonian_data = calculation_results.hamiltonian_data
    
    display_hamiltonian_properties(hamiltonian_matrix, hamiltonian_data)
    analyze_hamiltonian_eigenvalues(hamiltonian_matrix)
end

function display_hamiltonian_properties(hamiltonian_matrix, hamiltonian_data)
    println("Hamiltonian Properties:")
    println("  Matrix size: $(size(hamiltonian_matrix))")
    println("  Hermitian: $(ishermitian(hamiltonian_matrix))")
    println("  Components: $(keys(hamiltonian_data.h1e))")
end

function analyze_hamiltonian_eigenvalues(hamiltonian_matrix)
    eigenvalues = eigvals(Hermitian(hamiltonian_matrix))
    ground_state_energy = format_energy(eigenvalues[1])
    
    println("  Ground state energy: $ground_state_energy Ha")
    
    if length(eigenvalues) > 1
        first_excited_energy = format_energy(eigenvalues[2])
        energy_gap = format_energy(eigenvalues[2] - eigenvalues[1])
        
        println("  First excited state: $first_excited_energy Ha")
        println("  HOMO-LUMO gap: $energy_gap Ha")
    end
end

function demonstrate_results_saving(calculation_results::NEOResults)
    println("\nExample 5: Saving calculation results")
    println("-"^40)
    
    save_calculation_results(calculation_results)
    display_important_configurations(calculation_results.configs)
end

function save_calculation_results(results::NEOResults)
    output_filename = "h2_neo_hamiltonian.h5"
    
    save_hamiltonian(
        output_filename, 
        results.hamiltonian_data, 
        results.hamiltonian_matrix
    )
    
    println("Hamiltonian successfully saved to: $output_filename")
end

function display_important_configurations(configurations::Vector)
    println("\nMost important configurations:")
    
    max_configs_to_show = min(MAXIMUM_DISPLAYED_CONFIGURATIONS, length(configurations))
    
    for (index, configuration) in enumerate(configurations[1:max_configs_to_show])
        formatted_weight = round(configuration.weight, digits=4)
        println("  $index. $(configuration.name) - weight: $formatted_weight")
    end
end

function display_completion_message()
    println("\n" * "="^60)
    println("Clean Code examples completed successfully!")
    println("All functions demonstrate proper separation of concerns")
    println("and follow Robert C. Martin's Clean Code principles.")
    println("="^60)
end

function handle_example_error(error::Exception)
    println("\n❌ Example execution failed with error:")
    println("   $error")
    println("\nPlease ensure PySCF with NEO support is properly installed.")
end

# Execute the examples
run_clean_code_examples()
