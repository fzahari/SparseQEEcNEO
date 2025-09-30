"""
SparseQEEcNEO.jl - Main Module for Sparse QEE-cNEO Implementation
Version 6.0 - Complete implementation with Hamiltonian construction
"""
module SparseQEEcNEO

using PyCall
using LinearAlgebra
using SparseArrays
using Printf
using Statistics
using Base.Threads
using HDF5
using JSON

# Include all submodules in correct order
# IMPORTANT: Order matters due to dependencies
include("types.jl")
include("pyscf_interface.jl")
include("epc_functionals.jl")
include("nuclear_methods.jl")  # Must be before configuration_generation
include("configuration_generation.jl")  # Uses NuclearMethods
include("importance_analysis.jl")
include("qee_methods.jl")
include("hamiltonian_construction.jl")
# Note: cNEO modules are available as separate files but not included in main module due to dependencies

# Import from submodules
using .Types
using .PySCFInterface
using .EPCFunctionals
using .NuclearMethods
using .ConfigurationGeneration
using .ImportanceAnalysis
using .QEEMethods
using .HamiltonianConstruction

# Make sure this function is imported
import .NuclearMethods: apply_nuclear_orbital_truncation!

# Export main types and functions
export NEOConfig, Molecule, NEOCalculation, ConfigSelection, NEOResults
export Configuration, CompressedConfig, EPC_PARAMS
export sparse_qee_cneo, run_test_suite, run_demo_mode
export setup_pyscf, build_neo_molecule, run_neo_meanfield
export calculate_epc_energy, run_neo_mp2, run_neo_casci
export generate_configurations, select_important_configurations
export calculate_importance_metrics
export truncate_nuclear_orbitals, apply_nuclear_orbital_truncation!
export construct_second_quantized_hamiltonian, save_hamiltonian, load_hamiltonian
export HamiltonianData, analyze_hamiltonian_properties
export export_hamiltonian_openfermion

# Note: cNEO functionality available in advanced_examples/cneo/ directory

# ======================== Main Interface Function ========================

"""
    sparse_qee_cneo(mol::Molecule; kwargs...)

Main function for Sparse QEE-cNEO calculations with enhanced features.

# Arguments
- `mol::Molecule`: Molecular system with quantum nuclei specification
- `calc::NEOCalculation`: Calculation parameters (default: NEOCalculation())
- `config_sel::ConfigSelection`: Configuration selection parameters
- `neo_config::NEOConfig`: PySCF configuration

# Returns
- `NEOResults`: Complete results including configurations, energies, metrics, and Hamiltonian
"""
function sparse_qee_cneo(mol::Molecule; 
                        calc::NEOCalculation = NEOCalculation(),
                        config_sel::ConfigSelection = ConfigSelection(),
                        neo_config::NEOConfig = NEOConfig())
    
    performance_timer = create_performance_timer()
    
    try
        # Step 1: Initialize quantum chemistry environment
        meanfield_result = initialize_quantum_environment(mol, calc, neo_config)
        
        # Step 2: Apply nuclear orbital truncation if specified
        truncation_result = apply_orbital_truncation_if_needed(meanfield_result, config_sel)
        
        # Step 3: Calculate correlation corrections if requested
        correlation_result = calculate_correlation_corrections(
            meanfield_result, truncation_result, config_sel
        )
        
        # Step 4: Generate and select configurations
        configuration_result = generate_and_select_configurations(
            meanfield_result, correlation_result, config_sel
        )
        
        # Step 5: Construct Hamiltonian matrix
        hamiltonian_result = construct_hamiltonian_matrix(
            meanfield_result, configuration_result, config_sel
        )
        
        # Step 6: Apply additional energy corrections
        energy_corrections = calculate_energy_corrections(meanfield_result, calc)
        
        # Step 7: Assemble final results
        final_results = assemble_calculation_results(
            meanfield_result, correlation_result, configuration_result, 
            hamiltonian_result, energy_corrections, performance_timer
        )
        
        print_results_summary(final_results, mol, calc, config_sel)
        return final_results
        
    catch error
        handle_calculation_error(error, mol, calc)
        rethrow(error)
    end
end

# ======================== Clean Code Helper Functions ========================

# Constants for Clean Code compliance
const DEFAULT_ELECTRON_COUNT = 2
const DEFAULT_ORBITAL_COUNT = 4
const MEMORY_CONVERSION_FACTOR_MB = 1024^2
const MP2_WEIGHT_ENHANCEMENT_FACTOR = 10.0
const HAMILTONIAN_SPARSITY_THRESHOLD = 1e-12
const IMPORTANCE_DISPLAY_PRECISION = 3
const ENERGY_DISPLAY_PRECISION = 6

"""
    PerformanceTimer

Track computation time and memory usage during calculation.
"""
struct PerformanceTimer
    start_time::Float64
    initial_memory::Float64
end

function create_performance_timer()
    return PerformanceTimer(
        time(), 
        Base.gc_live_bytes() / MEMORY_CONVERSION_FACTOR_MB
    )
end

function calculate_elapsed_metrics(timer::PerformanceTimer)
    elapsed_time = time() - timer.start_time
    memory_used = Base.gc_live_bytes() / MEMORY_CONVERSION_FACTOR_MB - timer.initial_memory
    return elapsed_time, memory_used
end

"""
    MeanfieldResult

Encapsulate results from quantum chemistry mean-field calculation.
"""
struct MeanfieldResult
    meanfield_object::Any
    neo_molecule::Any
    pyscf_module::Any
end

"""
    TruncationResult

Results from nuclear orbital truncation process.
"""
struct TruncationResult
    was_truncated::Bool
    truncation_info::Dict{String, Any}
end

"""
    CorrelationResult

Results from correlation correction calculations (MP2, etc.).
"""
struct CorrelationResult
    mp2_energy::Float64
    t2_amplitudes::Any
end

"""
    ConfigurationResult

Results from configuration generation and selection.
"""
struct ConfigurationResult
    selected_configurations::Vector
    importance_data::Any
    final_method_used::String
    qubit_count::Int
end

"""
    HamiltonianResult

Results from Hamiltonian construction.
"""
struct HamiltonianResult
    hamiltonian_data::Any
    hamiltonian_matrix::Any
end

# ======================== Step 1: Initialize Quantum Environment ========================

function initialize_quantum_environment(mol::Molecule, calc::NEOCalculation, neo_config::NEOConfig)
    pyscf_module, has_neo = setup_pyscf(neo_config)
    
    validate_neo_availability(has_neo)
    
    neo_molecule = build_neo_molecule(mol, pyscf_module)
    meanfield_object = run_neo_meanfield(neo_molecule, calc, pyscf_module)
    
    return MeanfieldResult(meanfield_object, neo_molecule, pyscf_module)
end

function validate_neo_availability(has_neo::Bool)
    if !has_neo
        throw(ArgumentError("NEO module not available in PySCF. Please install PySCF with NEO support."))
    end
end

# ======================== Step 2: Orbital Truncation ========================

function apply_orbital_truncation_if_needed(meanfield_result::MeanfieldResult, config_sel::ConfigSelection)
    if should_apply_truncation(config_sel)
        was_truncated, truncation_info = apply_nuclear_orbital_truncation!(
            meanfield_result.meanfield_object, 
            meanfield_result.neo_molecule, 
            config_sel
        )
        return TruncationResult(was_truncated, truncation_info)
    else
        return TruncationResult(false, Dict{String, Any}())
    end
end

function should_apply_truncation(config_sel::ConfigSelection)
    return config_sel.max_nuc_orbs > 0
end

# ======================== Step 3: Correlation Corrections ========================

function calculate_correlation_corrections(meanfield_result::MeanfieldResult, 
                                         truncation_result::TruncationResult,
                                         config_sel::ConfigSelection)
    if should_calculate_mp2_correlation(config_sel)
        return attempt_mp2_calculation(meanfield_result, truncation_result, config_sel)
    else
        return CorrelationResult(0.0, nothing)
    end
end

function should_calculate_mp2_correlation(config_sel::ConfigSelection)
    return config_sel.method in ["mp2", "mp2_enhanced"]
end

function attempt_mp2_calculation(meanfield_result::MeanfieldResult,
                                truncation_result::TruncationResult,
                                config_sel::ConfigSelection)
    if should_skip_mp2_due_to_truncation(truncation_result, config_sel)
        @info "Skipping NEO-MP2 due to orbital truncation"
        return CorrelationResult(0.0, nothing)
    end
    
    try
        mp2_object, correlation_energy, t2_data = run_neo_mp2(
            meanfield_result.meanfield_object, 
            meanfield_result.neo_molecule
        )
        
        t2_amplitudes = extract_t2_amplitudes(t2_data, config_sel)
        
        return CorrelationResult(correlation_energy, t2_amplitudes)
        
    catch error
        handle_mp2_failure(error, truncation_result)
        return CorrelationResult(0.0, nothing)
    end
end

function should_skip_mp2_due_to_truncation(truncation_result::TruncationResult, config_sel::ConfigSelection)
    return truncation_result.was_truncated && !config_sel.force_mp2_with_truncation
end

function handle_mp2_failure(error::Exception, truncation_result::TruncationResult)
    if truncation_result.was_truncated
        @warn "NEO-MP2 failed with truncated orbitals (expected): $error"
        @info "Continuing without MP2 correlation"
    else
        @warn "NEO-MP2 failed: $error"
    end
end

# ======================== Step 4: Configuration Generation ========================

function generate_and_select_configurations(meanfield_result::MeanfieldResult,
                                          correlation_result::CorrelationResult,
                                          config_sel::ConfigSelection)
    initial_configurations = generate_configurations(
        meanfield_result.meanfield_object,
        meanfield_result.neo_molecule,
        config_sel,
        correlation_result.t2_amplitudes
    )
    
    importance_data = calculate_importance_metrics(
        initial_configurations, 
        meanfield_result.neo_molecule, 
        config_sel
    )
    
    optimized_result = optimize_method_selection(
        meanfield_result, correlation_result, initial_configurations, 
        importance_data, config_sel
    )
    
    selected_configurations, qubit_count = select_important_configurations(
        optimized_result.configurations, optimized_result.config_selection
    )
    
    return ConfigurationResult(
        selected_configurations,
        optimized_result.importance_data,
        optimized_result.config_selection.method,
        qubit_count
    )
end

struct OptimizedConfigurationResult
    configurations::Vector
    importance_data::Any
    config_selection::ConfigSelection
end

function optimize_method_selection(meanfield_result::MeanfieldResult,
                                 correlation_result::CorrelationResult,
                                 configurations::Vector,
                                 importance_data::Any,
                                 config_sel::ConfigSelection)
    if should_attempt_method_switching(config_sel, importance_data)
        return attempt_method_switching(
            meanfield_result, correlation_result, configurations, 
            importance_data, config_sel
        )
    else
        return OptimizedConfigurationResult(configurations, importance_data, config_sel)
    end
end

function should_attempt_method_switching(config_sel::ConfigSelection, importance_data::Any)
    return config_sel.auto_switch_method && 
           importance_data.total_importance < config_sel.importance_threshold_switch
end

function attempt_method_switching(meanfield_result::MeanfieldResult,
                                correlation_result::CorrelationResult,
                                original_configurations::Vector,
                                original_importance::Any,
                                original_config_sel::ConfigSelection)
    @warn "Low importance ($(round(original_importance.total_importance, digits=IMPORTANCE_DISPLAY_PRECISION))), switching methods"
    
    alternative_methods = ["neo_cneo", "neo_enhanced", "hybrid_final"]
    best_result = OptimizedConfigurationResult(original_configurations, original_importance, original_config_sel)
    
    for alternative_method in alternative_methods
        if alternative_method == original_config_sel.method
            continue
        end
        
        improved_result = try_alternative_method(
            meanfield_result, correlation_result, alternative_method, original_config_sel
        )
        
        if improved_result.importance_data.total_importance > best_result.importance_data.total_importance
            best_result = improved_result
            @info "Switched to method: $alternative_method"
            break
        end
    end
    
    return best_result
end

function try_alternative_method(meanfield_result::MeanfieldResult,
                              correlation_result::CorrelationResult,
                              method_name::String,
                              original_config_sel::ConfigSelection)
    alternative_config_sel = deepcopy(original_config_sel)
    alternative_config_sel.method = method_name
    
    alternative_configurations = generate_configurations(
        meanfield_result.meanfield_object,
        meanfield_result.neo_molecule,
        alternative_config_sel,
        correlation_result.t2_amplitudes
    )
    
    alternative_importance = calculate_importance_metrics(
        alternative_configurations,
        meanfield_result.neo_molecule,
        alternative_config_sel
    )
    
    return OptimizedConfigurationResult(
        alternative_configurations, alternative_importance, alternative_config_sel
    )
end

# ======================== Step 5: Hamiltonian Construction ========================

function construct_hamiltonian_matrix(meanfield_result::MeanfieldResult,
                                    configuration_result::ConfigurationResult,
                                    config_sel::ConfigSelection)
    hamiltonian_data, hamiltonian_matrix = construct_second_quantized_hamiltonian(
        meanfield_result.meanfield_object,
        meanfield_result.neo_molecule,
        configuration_result.selected_configurations,
        config_sel
    )
    
    return HamiltonianResult(hamiltonian_data, hamiltonian_matrix)
end

# ======================== Step 6: Energy Corrections ========================

function calculate_energy_corrections(meanfield_result::MeanfieldResult, calc::NEOCalculation)
    if should_apply_epc_correction(calc)
        return calculate_epc_energy(
            meanfield_result.meanfield_object, 
            meanfield_result.neo_molecule, 
            calc.epc
        )
    else
        return 0.0
    end
end

function should_apply_epc_correction(calc::NEOCalculation)
    return calc.epc != "none" && calc.epc != ""
end

# ======================== Step 7: Results Assembly ========================

function assemble_calculation_results(meanfield_result::MeanfieldResult,
                                    correlation_result::CorrelationResult,
                                    configuration_result::ConfigurationResult,
                                    hamiltonian_result::HamiltonianResult,
                                    epc_energy::Float64,
                                    performance_timer::PerformanceTimer)
    
    elapsed_time, memory_used = calculate_elapsed_metrics(performance_timer)
    
    orbital_analysis = analyze_orbital_structure(
        meanfield_result.meanfield_object, meanfield_result.neo_molecule
    )
    
    savings_metrics = calculate_computational_savings(
        configuration_result.selected_configurations,
        configuration_result.qubit_count,
        orbital_analysis.total_orbitals
    )
    
    base_energy = extract_base_energy(meanfield_result.meanfield_object)
    total_energy = base_energy + correlation_result.mp2_energy + epc_energy
    
    return NEOResults(
        configuration_result.selected_configurations,
        base_energy,
        correlation_result.mp2_energy,
        total_energy,
        length(configuration_result.selected_configurations),
        configuration_result.qubit_count,
        configuration_result.importance_data.total_importance,
        savings_metrics.determinant_savings,
        savings_metrics.qubit_savings,
        elapsed_time,
        memory_used,
        configuration_result.final_method_used,
        orbital_analysis.orbitals_per_species,
        configuration_result.importance_data.neo_metrics,
        false, # orbital_truncation - will be updated by truncation logic
        hamiltonian_result.hamiltonian_data,
        hamiltonian_result.hamiltonian_matrix
    )
end

struct OrbitalAnalysis
    total_orbitals::Int
    electron_count::Int
    orbitals_per_species::Vector{Int}
end

function analyze_orbital_structure(meanfield_object::Any, neo_molecule::Any)
    if has_electronic_component(meanfield_object)
        electronic_component = meanfield_object.components["e"]
        electron_count = extract_electron_count(electronic_component)
        orbital_count = extract_orbital_count(electronic_component)
        
        orbitals_per_species = [orbital_count]
        
        if has_nuclear_components(neo_molecule, meanfield_object)
            nuclear_orbital_count = extract_nuclear_orbital_count(meanfield_object)
            orbitals_per_species = [orbital_count, nuclear_orbital_count]
        end
        
        return OrbitalAnalysis(orbital_count, electron_count, orbitals_per_species)
    else
        return OrbitalAnalysis(DEFAULT_ORBITAL_COUNT, DEFAULT_ELECTRON_COUNT, [DEFAULT_ORBITAL_COUNT])
    end
end

function has_electronic_component(meanfield_object::Any)
    return pybuiltin("hasattr")(meanfield_object, "components") && 
           haskey(meanfield_object.components, "e")
end

function extract_electron_count(electronic_component::Any)
    if pybuiltin("hasattr")(electronic_component, "nelectron")
        return convert(Int, electronic_component.nelectron)
    elseif pybuiltin("hasattr")(electronic_component, "mo_occ")
        mo_occ = collect(electronic_component.mo_occ)
        return convert(Int, sum(mo_occ))
    else
        return DEFAULT_ELECTRON_COUNT
    end
end

function extract_orbital_count(electronic_component::Any)
    if pybuiltin("hasattr")(electronic_component, "mo_occ")
        return length(electronic_component.mo_occ)
    else
        return DEFAULT_ORBITAL_COUNT
    end
end

function has_nuclear_components(neo_molecule::Any, meanfield_object::Any)
    return pybuiltin("hasattr")(neo_molecule, "nuc_num") && neo_molecule.nuc_num > 0 &&
           pybuiltin("hasattr")(meanfield_object, "components") && haskey(meanfield_object.components, "n0")
end

function extract_nuclear_orbital_count(meanfield_object::Any)
    nuclear_component = meanfield_object.components["n0"]
    if pybuiltin("hasattr")(nuclear_component, "mo_occ")
        return length(nuclear_component.mo_occ)
    else
        return 0
    end
end

struct ComputationalSavings
    determinant_savings::Float64
    qubit_savings::Int
end

function calculate_computational_savings(selected_configs::Vector, qubit_count::Int, total_orbitals::Int)
    determinant_savings = 2^total_orbitals / max(length(selected_configs), 1)
    qubit_savings = total_orbitals - qubit_count
    
    return ComputationalSavings(determinant_savings, qubit_savings)
end

function extract_base_energy(meanfield_object::Any)
    if pybuiltin("hasattr")(meanfield_object, "e_tot")
        return convert(Float64, meanfield_object.e_tot)
    else
        return 0.0
    end
end

# ======================== Error Handling ========================

function handle_calculation_error(error::Exception, mol::Molecule, calc::NEOCalculation)
    @error "Calculation failed for molecule with $(length(mol.quantum_nuc)) quantum nuclei"
    @error "Method: $(calc.xc), EPC: $(calc.epc)"
    @error "Error details: $error"
end

# ======================== Results Summary ========================

function print_results_summary(results::NEOResults, mol::Molecule, calc::NEOCalculation, 
                              config_sel::ConfigSelection = ConfigSelection())
    
    println("\n" * "="^60)
    println("SPARSE QEE-CNEO CALCULATION RESULTS")
    println("="^60)
    
    print_method_information(calc, results)
    print_energy_information(results)
    print_configuration_information(results, mol)
    print_hamiltonian_information(results)
    print_performance_information(results)
    print_neo_specific_metrics(results)
    
    println("="^60)
end

function print_method_information(calc::NEOCalculation, results::NEOResults)
    println("\nMethod Information:")
    method_display = create_method_display_string(calc, results)
    println("  Method: $method_display")
    
    if should_display_epc_functional(calc)
        println("  EPC functional: $(calc.epc)")
    end
end

function create_method_display_string(calc::NEOCalculation, results::NEOResults)
    base_method = calc.xc
    configuration_method = results.method_used == "mp2" ? "MP2" : results.method_used
    return "$base_method+$configuration_method"
end

function should_display_epc_functional(calc::NEOCalculation)
    return calc.epc != "none" && !isempty(calc.epc)
end

function print_energy_information(results::NEOResults)
    println("\nEnergy Information:")
    println("  Base energy: $(round(results.energy, digits=ENERGY_DISPLAY_PRECISION)) Ha")
    
    if has_mp2_correlation(results)
        println("  MP2 correlation: $(round(results.mp2_correlation, digits=ENERGY_DISPLAY_PRECISION)) Ha")
    end
    
    println("  Total energy: $(round(results.total_energy, digits=ENERGY_DISPLAY_PRECISION)) Ha")
end

function has_mp2_correlation(results::NEOResults)
    return results.mp2_correlation != 0.0
end

function print_configuration_information(results::NEOResults, mol::Molecule)
    println("\nConfiguration Information:")
    
    electron_count = extract_electron_count_from_orbitals(results.orbitals_per_species)
    println("  Electrons: $electron_count")
    println("  Quantum nuclei: $(length(mol.quantum_nuc))")
    println("  Selected configurations: $(results.n_configs)")
    println("  Qubits required: $(results.n_qubits)")
    println("  Captured importance: $(format_percentage(results.captured_importance))%")
    println("  Computational savings:")
    println("    Determinant reduction: $(round(results.det_savings, digits=1))×")
    println("    Qubit reduction: $(results.qubit_savings)")
    println("  Orbitals per species: $(results.orbitals_per_species)")
    
    if results.orbital_truncation
        println("  ⚠️  Nuclear orbitals were truncated")
    end
end

function extract_electron_count_from_orbitals(orbitals_per_species::Vector{Int})
    # Estimate electron count from orbital count (rough approximation)
    return length(orbitals_per_species) > 0 ? orbitals_per_species[1] ÷ 2 : DEFAULT_ELECTRON_COUNT
end

function format_percentage(value::Float64)
    return round(value * 100, digits=1)
end

function print_hamiltonian_information(results::NEOResults)
    if has_hamiltonian_matrix(results)
        println("\nHamiltonian Information:")
        print_hamiltonian_basic_properties(results.hamiltonian_matrix)
        print_hamiltonian_analysis(results.hamiltonian_matrix)
    end
end

function has_hamiltonian_matrix(results::NEOResults)
    return results.hamiltonian_matrix !== nothing
end

function print_hamiltonian_basic_properties(hamiltonian_matrix::Any)
    matrix_size = size(hamiltonian_matrix)
    is_hermitian = ishermitian(hamiltonian_matrix)
    hermitian_status = is_hermitian ? "✓" : "✗"
    
    println("  Matrix dimensions: $matrix_size")
    println("  Hermitian: $hermitian_status")
end

function print_hamiltonian_analysis(hamiltonian_matrix::Any)
    if is_matrix_analyzable(hamiltonian_matrix)
        properties = analyze_hamiltonian_properties(hamiltonian_matrix)
        println("  Ground state energy: $(round(properties.ground_state_energy, digits=ENERGY_DISPLAY_PRECISION)) Ha")
        println("  Matrix sparsity: $(format_percentage(properties.sparsity))%")
    end
end

function is_matrix_analyzable(hamiltonian_matrix::Any)
    return size(hamiltonian_matrix, 1) > 0
end

function print_performance_information(results::NEOResults)
    println("\nPerformance Metrics:")
    println("  Computation time: $(round(results.computation_time, digits=IMPORTANCE_DISPLAY_PRECISION)) seconds")
    println("  Memory usage: $(round(results.memory_used, digits=1)) MB")
end

function print_neo_specific_metrics(results::NEOResults)
    if has_neo_metrics(results)
        println("\nNEO-Specific Metrics:")
        print_neo_importance_metrics(results.neo_metrics)
        print_neo_participation_metrics(results.neo_metrics)
    end
end

function has_neo_metrics(results::NEOResults)
    return results.neo_metrics !== nothing
end

function print_neo_importance_metrics(neo_metrics::Any)
    if has_standard_importance_field(neo_metrics)
        standard_importance = round(neo_metrics.standard_importance, digits=IMPORTANCE_DISPLAY_PRECISION)
        println("  Standard importance: $standard_importance")
    end
end

function has_standard_importance_field(neo_metrics::Any)
    return hasfield(typeof(neo_metrics), :standard_importance)
end

function print_neo_participation_metrics(neo_metrics::Any)
    nuclear_participation = format_percentage(neo_metrics.nuclear_participation)
    coupling_contribution = round(neo_metrics.coupling_contribution, digits=IMPORTANCE_DISPLAY_PRECISION)
    
    println("  Nuclear participation: $nuclear_participation%")
    println("  Coupling contribution: $coupling_contribution")
end

# ======================== Hamiltonian Analysis ========================

"""
    HamiltonianProperties

Structure to hold analyzed properties of a Hamiltonian matrix.
"""
struct HamiltonianProperties
    ground_state_energy::Float64
    energy_gap::Float64
    sparsity::Float64
    condition_number::Float64
    eigenvalues::Vector{Float64}
end

"""
    analyze_hamiltonian_properties(hamiltonian_matrix::Matrix)

Analyze key properties of the Hamiltonian matrix.

# Returns
- `HamiltonianProperties`: Struct containing analysis results
"""
function analyze_hamiltonian_properties(hamiltonian_matrix::Matrix)
    if is_empty_matrix(hamiltonian_matrix)
        return create_empty_hamiltonian_properties()
    end
    
    eigenvalues = calculate_hamiltonian_eigenvalues(hamiltonian_matrix)
    ground_state_energy = extract_ground_state_energy(eigenvalues)
    energy_gap = calculate_energy_gap(eigenvalues)
    matrix_sparsity = calculate_matrix_sparsity(hamiltonian_matrix)
    condition_number = calculate_condition_number(hamiltonian_matrix)
    
    return HamiltonianProperties(
        ground_state_energy,
        energy_gap,
        matrix_sparsity,
        condition_number,
        eigenvalues
    )
end

function is_empty_matrix(hamiltonian_matrix::Matrix)
    return size(hamiltonian_matrix, 1) == 0
end

function create_empty_hamiltonian_properties()
    return HamiltonianProperties(
        0.0,  # ground_state_energy
        0.0,  # energy_gap
        1.0,  # sparsity (empty matrix is 100% sparse)
        1.0,  # condition_number
        Float64[]  # eigenvalues
    )
end

function calculate_hamiltonian_eigenvalues(hamiltonian_matrix::Matrix)
    return try
        eigenvals_hermitian = eigvals(Hermitian(hamiltonian_matrix))
        sort!(eigenvals_hermitian)
        eigenvals_hermitian
    catch
        eigenvals_general = eigvals(hamiltonian_matrix)
        sort!(eigenvals_general)
        eigenvals_general
    end
end

function extract_ground_state_energy(eigenvalues::Vector{Float64})
    return length(eigenvalues) > 0 ? eigenvalues[1] : 0.0
end

function calculate_energy_gap(eigenvalues::Vector{Float64})
    return length(eigenvalues) > 1 ? eigenvalues[2] - eigenvalues[1] : 0.0
end

function calculate_matrix_sparsity(hamiltonian_matrix::Matrix)
    zero_elements_count = count(abs.(hamiltonian_matrix) .< HAMILTONIAN_SPARSITY_THRESHOLD)
    total_elements = length(hamiltonian_matrix)
    return zero_elements_count / total_elements
end

function calculate_condition_number(hamiltonian_matrix::Matrix)
    return try
        cond(hamiltonian_matrix)
    catch
        Inf
    end
end

# ======================== Test Suite ========================

"""
    run_test_suite(config::NEOConfig; run_modular_tests=true)

Run comprehensive test suite for the Sparse QEE-cNEO implementation.
"""
function run_test_suite(config::NEOConfig = NEOConfig(); run_modular_tests::Bool = true)
    println("\n" * "="^60)
    println("Enhanced cNEO Sparse QEE Tests (v6.0)")
    println("Complete implementation with Hamiltonian construction")
    println("="^60 * "\n")
    
    # Check PySCF availability
    pyscf, has_neo = setup_pyscf(config)
    
    if !has_neo
        println("Running in demo mode without NEO calculations...")
        run_demo_mode()
        return Dict("demo" => "completed")
    end
    
    # Run test cases
    test_results = Dict{String, Any}()
    
    # Test 1: H2 with quantum proton
    println("\nTest 1: H2 with quantum proton")
    println("-" * 40)
    mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
    try
        res = sparse_qee_cneo(mol_h2, neo_config=config)
        test_results["H2"] = res
        
        # Verify Hamiltonian
        if res.hamiltonian_matrix !== nothing
            println("  Hamiltonian constructed: $(size(res.hamiltonian_matrix))")
        end
        println("✓ Test passed")
    catch e
        println("✗ Test failed: $e")
        test_results["H2"] = e
    end
    
    # Test 2: Water with quantum proton
    println("\nTest 2: Water with quantum proton")
    println("-" * 40)
    mol_water = Molecule("O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", 
                        "6-31g", quantum_nuc=[1])
    calc = NEOCalculation(xc="B3LYP", epc="17-2")
    config_sel = ConfigSelection(method="neo_cneo", max_configs=100)
    
    try
        res = sparse_qee_cneo(mol_water, calc=calc, config_sel=config_sel, 
                             neo_config=config)
        test_results["Water"] = res
        println("✓ Test passed")
    catch e
        println("✗ Test failed: $e")
        test_results["Water"] = e
    end
    
    # Test 3: Method comparison
    println("\nTest 3: Method comparison for HCN")
    println("-" * 40)
    mol_hcn = Molecule("H 0 0 -1.064; C 0 0 0; N 0 0 1.156", 
                      "sto-3g", quantum_nuc=[0])
    
    for method in ["mp2", "neo_cneo", "hybrid_final"]
        println("\n  Testing method: $method")
        config_sel = ConfigSelection(method=method, max_configs=50)
        try
            res = sparse_qee_cneo(mol_hcn, config_sel=config_sel, neo_config=config)
            test_results["HCN_$method"] = res
            println("  ✓ Configurations: $(res.n_configs), Importance: $(round(res.captured_importance, digits=3))")
        catch e
            println("  ✗ Failed: $e")
            test_results["HCN_$method"] = e
        end
    end
    
    # Modular tests
    if run_modular_tests
        println("\nTest 4: Modular component tests")
        println("-" * 40)
        test_results["modular"] = run_modular_component_tests(config)
    end
    
    # Summary
    println("\n" * "="^60)
    println("Test Summary:")
    passed = count(v -> !(v isa Exception) for (k,v) in test_results)
    total = length(test_results)
    println("Passed: $passed/$total")
    println("="^60)
    
    return test_results
end

# ======================== Modular Component Tests ========================

"""
    run_modular_component_tests(config::NEOConfig)

Run tests for individual components in isolation.
"""
function run_modular_component_tests(config::NEOConfig)
    results = Dict{String, Any}()
    
    println("\n  a) Testing electronic component only...")
    mol_elec = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")  # No quantum nuclei
    try
        res = sparse_qee_cneo(mol_elec, neo_config=config)
        results["electronic_only"] = res
        println("    ✓ Electronic-only calculation successful")
    catch e
        results["electronic_only"] = e
        println("    ✗ Failed: $e")
    end
    
    println("\n  b) Testing nuclear methods...")
    mol_nuc = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0, 1])
    config_sel = ConfigSelection(
        method="neo_cneo",
        max_configs=20,
        debug_nuclear=true
    )
    try
        res = sparse_qee_cneo(mol_nuc, config_sel=config_sel, neo_config=config)
        results["nuclear_methods"] = res
        println("    ✓ Nuclear methods successful")
    catch e
        results["nuclear_methods"] = e
        println("    ✗ Failed: $e")
    end
    
    println("\n  c) Testing EPC functionals...")
    mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
    for epc in ["17-1", "17-2", "18-1", "18-2"]
        calc = NEOCalculation(xc="B3LYP", epc=epc)
        config_sel = ConfigSelection(method="mp2", max_configs=10)
        try
            res = sparse_qee_cneo(mol_h2, calc=calc, config_sel=config_sel, neo_config=config)
            results["epc_$epc"] = res.total_energy
            println("    ✓ EPC $epc: $(round(res.total_energy, digits=6)) Ha")
        catch e
            results["epc_$epc"] = e
            println("    ✗ EPC $epc failed")
        end
    end
    
    println("\n  d) Testing configuration compression...")
    config_sel = ConfigSelection(
        method="mp2",
        max_configs=200,
        use_compression=true
    )
    try
        res = sparse_qee_cneo(mol_elec, config_sel=config_sel, neo_config=config)
        compressed = count(c -> isa(c, CompressedConfig), res.configs)
        results["compression"] = compressed
        println("    ✓ Compressed $compressed configurations")
    catch e
        results["compression"] = e
        println("    ✗ Compression failed")
    end
    
    return results
end

# ======================== Demo Mode ========================

function run_demo_mode()
    println("\n" * "="^60)
    println("Running cNEO Demo Mode")
    println("="^60 * "\n")
    
    println("This demonstrates the implementation structure without calculations.\n")
    
    println("1. Configuration Generation Methods:")
    println("   - mp2: MP2-based selection")
    println("   - mp2_enhanced: MP2 with t2 amplitudes")
    println("   - casci: Complete active space")
    println("   - neo_cneo: cNEO-enhanced generation")
    println("   - neo_enhanced: Enhanced NEO approach")
    println("   - hybrid_final: Combined approach")
    
    println("\n2. EPC Functionals Available:")
    for (name, params) in EPC_PARAMS
        println("   - $name: a=$(params["a"]), b=$(params["b"]), c=$(params["c"])")
    end
    
    println("\n3. Expected Configuration Types:")
    println("   - Electronic: E(i→a), D(ij→ab)")
    println("   - Nuclear: N1(i→a)")
    println("   - Coupled: C(E(i→a)+N1(j→b))")
    
    println("\n4. Hamiltonian Construction:")
    println("   - One-body integrals (h1e)")
    println("   - Two-body integrals (h2e) - sparse storage")
    println("   - Electron-nuclear coupling terms")
    println("   - Active space determination")
    println("   - HDF5 storage format")
    
    println("\n5. Modular Testing Approach:")
    println("   - Test each component independently")
    println("   - Electronic-only calculations")
    println("   - Nuclear method validation")
    println("   - EPC functional comparison")
    println("   - Configuration compression")
    println("   - Hamiltonian construction")
    
    println("\n" * "="^60)
    println("Install PySCF with NEO for full functionality")
    println("="^60)
end

end # module
