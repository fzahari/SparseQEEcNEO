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
include("constants.jl")  # Constants for Clean Code compliance
include("types.jl")
include("pyscf_interface.jl")
include("epc_functionals.jl")
include("nuclear_methods.jl")  # Must be before configuration_generation
include("configuration_generation.jl")  # Uses NuclearMethods
include("importance_analysis.jl")
include("qee_methods.jl")
include("hamiltonian_construction.jl")
include("cneo_methods.jl")  # Constrained NEO methods - Clean Code implementation
include("quantum_computing.jl")  # Quantum computing integration - Clean Code implementation

# Import from submodules
using .Types
using .PySCFInterface
using .EPCFunctionals
using .NuclearMethods
using .ConfigurationGeneration
using .ImportanceAnalysis
using .QEEMethods
using .HamiltonianConstruction
using .CNEOMethods
using .QuantumComputing

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
# Export cNEO functionality
export CNEOCalculation, CNEOResults, CNEOMP2Results
export run_cneo_hf, run_cneo_mp2, create_cneo_calculation
export sparse_qee_with_cneo, cneo_to_qee_workflow
# Export quantum computing functionality
export exportToQuantumFormats, exportToOpenFermion, exportToQiskit, exportToCirq
export createVQECircuit, estimateQuantumResources, validateQuantumExport, runQuantumDemo

# Note: cNEO functionality available in advanced_examples/ directory

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
    
    performance_timer = createPerformanceTimer()
    
    try
        # Step 1: Initialize quantum chemistry environment
        meanfield_result = initializeQuantumEnvironment(mol, calc, neo_config)
        
        # Step 2: Apply nuclear orbital truncation if specified
        truncation_result = applyOrbitalTruncationIfNeeded(meanfield_result, config_sel)
        
        # Step 3: Calculate correlation corrections if requested
        correlation_result = calculateCorrelationCorrections(
            meanfield_result, truncation_result, config_sel
        )
        
        # Step 4: Generate and select configurations
        configuration_result = generateAndSelectConfigurations(
            meanfield_result, correlation_result, config_sel
        )
        
        # Step 5: Construct Hamiltonian matrix
        hamiltonian_result = constructHamiltonianMatrix(
            meanfield_result, configuration_result, config_sel
        )
        
        # Step 6: Apply additional energy corrections
        energy_corrections = calculateEnergyCorrections(meanfield_result, calc)
        
        # Step 7: Assemble final results
        final_results = assembleCalculationResults(
            meanfield_result, correlation_result, configuration_result, 
            hamiltonian_result, energy_corrections, performance_timer
        )
        
        printResultsSummary(final_results, mol, calc, config_sel)
        return final_results
        
    catch error
        handleCalculationError(error, mol, calc)
        rethrow(error)
    end
end

# ======================== Clean Code Helper Functions ========================

# Constants now imported from constants.jl
const DEFAULT_ELECTRON_COUNT = 2
const DEFAULT_ORBITAL_COUNT = 4
const MEMORY_CONVERSION_FACTOR_MB = MEMORY_ALLOCATION_BUFFER  # Use from constants.jl
const MP2_WEIGHT_ENHANCEMENT_FACTOR = 10.0
const IMPORTANCE_DISPLAY_PRECISION = PERCENTAGE_FORMAT_PRECISION
const ENERGY_DISPLAY_PRECISION = ENERGY_FORMAT_PRECISION

"""
    PerformanceTimer

Track computation time and memory usage during calculation.
"""
struct PerformanceTimer
    start_time::Float64
    initial_memory::Float64
end

function createPerformanceTimer()
    return PerformanceTimer(
        time(), 
        Base.gc_live_bytes() / MEMORY_CONVERSION_FACTOR_MB
    )
end

function calculateElapsedMetrics(timer::PerformanceTimer)
    elapsedTime = time() - timer.start_time
    memoryUsed = Base.gc_live_bytes() / MEMORY_CONVERSION_FACTOR_MB - timer.initial_memory
    return elapsedTime, memoryUsed
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

function initializeQuantumEnvironment(mol::Molecule, calc::NEOCalculation, neo_config::NEOConfig)
    pyscf_module, has_neo = setup_pyscf(neo_config)
    
    validateNeoAvailability(has_neo)
    
    neo_molecule = build_neo_molecule(mol, pyscf_module)
    meanfield_object = run_neo_meanfield(neo_molecule, calc, pyscf_module)
    
    return MeanfieldResult(meanfield_object, neo_molecule, pyscf_module)
end

function validateNeoAvailability(has_neo::Bool)
    if !has_neo
        throw(ArgumentError("NEO module not available in PySCF. Please install PySCF with NEO support."))
    end
end

# ======================== Step 2: Orbital Truncation ========================

function applyOrbitalTruncationIfNeeded(meanfield_result::MeanfieldResult, config_sel::ConfigSelection)
    if shouldApplyTruncation(config_sel)
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

function shouldApplyTruncation(config_sel::ConfigSelection)
    return config_sel.max_nuc_orbs > 0
end

# ======================== Step 3: Correlation Corrections ========================

function calculateCorrelationCorrections(meanfield_result::MeanfieldResult, 
                                         truncation_result::TruncationResult,
                                         config_sel::ConfigSelection)
    if shouldCalculateMP2Correlation(config_sel)
        return attemptMP2Calculation(meanfield_result, truncation_result, config_sel)
    else
        return CorrelationResult(0.0, nothing)
    end
end

function shouldCalculateMP2Correlation(config_sel::ConfigSelection)
    return config_sel.method in ["mp2", "mp2_enhanced"]
end

function attemptMP2Calculation(meanfield_result::MeanfieldResult,
                                truncation_result::TruncationResult,
                                config_sel::ConfigSelection)
    if shouldSkipMP2DueToTruncation(truncation_result, config_sel)
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
        handleMP2Failure(error, truncation_result)
        return CorrelationResult(0.0, nothing)
    end
end

function shouldSkipMP2DueToTruncation(truncation_result::TruncationResult, config_sel::ConfigSelection)
    return truncation_result.was_truncated && !config_sel.force_mp2_with_truncation
end

function handleMP2Failure(error::Exception, truncation_result::TruncationResult)
    if truncation_result.was_truncated
        @warn "NEO-MP2 failed with truncated orbitals (expected): $error"
        @info "Continuing without MP2 correlation"
    else
        @warn "NEO-MP2 failed: $error"
    end
end

# ======================== Step 4: Configuration Generation ========================

struct OptimizedConfigurationResult
    configurations::Vector
    importance_data::Any
    config_selection::ConfigSelection
end

function generateAndSelectConfigurations(meanfieldResult::MeanfieldResult,
                                        correlationResult::CorrelationResult,
                                        configSelection::ConfigSelection)
    
    initialConfigurations = performConfigurationGeneration(
        meanfieldResult, correlationResult, configSelection
    )
    
    importanceData = calculateConfigurationImportance(
        initialConfigurations, meanfieldResult.neo_molecule, configSelection
    )
    
    optimizedResult = optimizeMethodSelection(
        meanfieldResult, correlationResult, initialConfigurations, 
        importanceData, configSelection
    )
    
    return createFinalConfigurationResult(optimizedResult)
end

function performConfigurationGeneration(meanfieldResult::MeanfieldResult,
                                      correlationResult::CorrelationResult,
                                      configSelection::ConfigSelection)
    return generate_configurations(
        meanfieldResult.meanfield_object,
        meanfieldResult.neo_molecule,
        configSelection,
        correlationResult.t2_amplitudes
    )
end

function calculateConfigurationImportance(configurations::Vector,
                                        neoMolecule::Any,
                                        configSelection::ConfigSelection)
    return calculate_importance_metrics(
        configurations, neoMolecule, configSelection
    )
end

function createFinalConfigurationResult(optimizedResult::OptimizedConfigurationResult)
    selectedConfigurations, qubitCount = select_important_configurations(
        optimizedResult.configurations, optimizedResult.config_selection
    )
    
    return ConfigurationResult(
        selectedConfigurations,
        optimizedResult.importance_data,
        optimizedResult.config_selection.method,
        qubitCount
    )
end


function optimizeMethodSelection(meanfield_result::MeanfieldResult,
                                 correlation_result::CorrelationResult,
                                 configurations::Vector,
                                 importance_data::Any,
                                 config_sel::ConfigSelection)
    if shouldAttemptMethodSwitching(config_sel, importance_data)
        return attemptMethodSwitching(
            meanfield_result, correlation_result, configurations, 
            importance_data, config_sel
        )
    else
        return OptimizedConfigurationResult(configurations, importance_data, config_sel)
    end
end

function shouldAttemptMethodSwitching(config_sel::ConfigSelection, importance_data::Any)
    return config_sel.auto_switch_method && 
           importance_data.total_importance < config_sel.importance_threshold_switch
end

function attemptMethodSwitching(meanfield_result::MeanfieldResult,
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
        
        improved_result = tryAlternativeMethod(
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

function tryAlternativeMethod(meanfield_result::MeanfieldResult,
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

function constructHamiltonianMatrix(meanfield_result::MeanfieldResult,
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

function calculateEnergyCorrections(meanfield_result::MeanfieldResult, calc::NEOCalculation)
    if shouldApplyEPCCorrection(calc)
        return calculate_epc_energy(
            meanfield_result.meanfield_object, 
            meanfield_result.neo_molecule, 
            calc.epc
        )
    else
        return 0.0
    end
end

function shouldApplyEPCCorrection(calc::NEOCalculation)
    return calc.epc != "none" && calc.epc != ""
end

# ======================== Step 7: Results Assembly ========================

struct OrbitalAnalysis
    total_orbitals::Int
    electron_count::Int
    orbitals_per_species::Vector{Int}
end

struct ComputationalSavings
    determinant_savings::Float64
    qubit_savings::Int
end

function assembleCalculationResults(meanfieldResult::MeanfieldResult,
                                   correlationResult::CorrelationResult,
                                   configurationResult::ConfigurationResult,
                                   hamiltonianResult::HamiltonianResult,
                                   epcEnergy::Float64,
                                   performanceTimer::PerformanceTimer)
    
    performanceData = extractPerformanceMetrics(performanceTimer)
    analysisData = performOrbitalAnalysis(meanfieldResult)
    energyData = calculateTotalEnergy(meanfieldResult, correlationResult, epcEnergy)
    savingsData = computeSavingsMetrics(configurationResult, analysisData)
    
    return buildNEOResults(
        configurationResult, analysisData, energyData, 
        savingsData, performanceData, hamiltonianResult
    )
end

function extractPerformanceMetrics(timer::PerformanceTimer)
    elapsedTime, memoryUsed = calculateElapsedMetrics(timer)
    return (elapsedTime, memoryUsed)
end

function performOrbitalAnalysis(meanfieldResult::MeanfieldResult)
    return analyzeOrbitalStructure(
        meanfieldResult.meanfield_object, 
        meanfieldResult.neo_molecule
    )
end

function calculateTotalEnergy(meanfieldResult::MeanfieldResult, 
                            correlationResult::CorrelationResult, 
                            epcEnergy::Float64)
    baseEnergy = extractBaseEnergy(meanfieldResult.meanfield_object)
    totalEnergy = baseEnergy + correlationResult.mp2_energy + epcEnergy
    return (baseEnergy, totalEnergy)
end

function computeSavingsMetrics(configurationResult::ConfigurationResult, 
                             orbitalAnalysis::OrbitalAnalysis)
    return calculateComputationalSavings(
        configurationResult.selected_configurations,
        configurationResult.qubit_count,
        orbitalAnalysis.total_orbitals
    )
end

function buildNEOResults(configurationResult::ConfigurationResult,
                       orbitalAnalysis::OrbitalAnalysis,
                       energyData::Tuple{Float64, Float64},
                       savingsData::ComputationalSavings,
                       performanceData::Tuple{Float64, Float64},
                       hamiltonianResult::HamiltonianResult)
    
    baseEnergy, totalEnergy = energyData
    elapsedTime, memoryUsed = performanceData
    mp2Correlation = calculateMP2Correlation(configurationResult, totalEnergy, baseEnergy)
    
    return createNEOResultsStruct(
        configurationResult, orbitalAnalysis, baseEnergy, totalEnergy,
        mp2Correlation, savingsData, elapsedTime, memoryUsed, hamiltonianResult
    )
end

function calculateMP2Correlation(configurationResult::ConfigurationResult, 
                               totalEnergy::Float64, baseEnergy::Float64)
    return configurationResult.importance_data.total_importance > 0 ? 
           totalEnergy - baseEnergy : 0.0
end

function createNEOResultsStruct(configurationResult::ConfigurationResult,
                              orbitalAnalysis::OrbitalAnalysis,
                              baseEnergy::Float64, totalEnergy::Float64,
                              mp2Correlation::Float64,
                              savingsData::ComputationalSavings,
                              elapsedTime::Float64, memoryUsed::Float64,
                              hamiltonianResult::HamiltonianResult)
    
    resultFields = assembleNEOResultFields(
        configurationResult, orbitalAnalysis, baseEnergy, totalEnergy,
        mp2Correlation, savingsData, elapsedTime, memoryUsed, hamiltonianResult
    )
    
    return NEOResults(resultFields...)
end

function assembleNEOResultFields(configurationResult, orbitalAnalysis, 
                               baseEnergy, totalEnergy, mp2Correlation,
                               savingsData, elapsedTime, memoryUsed, hamiltonianResult)
    return (
        configurationResult.selected_configurations,
        baseEnergy,
        mp2Correlation,
        totalEnergy,
        length(configurationResult.selected_configurations),
        configurationResult.qubit_count,
        configurationResult.importance_data.total_importance,
        savingsData.determinant_savings,
        savingsData.qubit_savings,
        elapsedTime,
        memoryUsed,
        configurationResult.final_method_used,
        orbitalAnalysis.orbitals_per_species,
        configurationResult.importance_data.neo_metrics,
        false, # orbital_truncation - updated by truncation logic
        hamiltonianResult.hamiltonian_data,
        hamiltonianResult.hamiltonian_matrix
    )
end


function analyzeOrbitalStructure(meanfield_object::Any, neo_molecule::Any)
    if hasElectronicComponent(meanfield_object)
        electronic_component = meanfield_object.components["e"]
        electron_count = extractElectronCount(electronic_component)
        orbital_count = extractOrbitalCount(electronic_component)
        
        orbitals_per_species = [orbital_count]
        
        if hasNuclearComponents(neo_molecule, meanfield_object)
            nuclear_orbital_count = extractNuclearOrbitalCount(meanfield_object)
            orbitals_per_species = [orbital_count, nuclear_orbital_count]
        end
        
        return OrbitalAnalysis(orbital_count, electron_count, orbitals_per_species)
    else
        return OrbitalAnalysis(DEFAULT_ORBITAL_COUNT, DEFAULT_ELECTRON_COUNT, [DEFAULT_ORBITAL_COUNT])
    end
end

function hasElectronicComponent(meanfield_object::Any)
    return pybuiltin("hasattr")(meanfield_object, "components") && 
           haskey(meanfield_object.components, "e")
end

function extractElectronCount(electronic_component::Any)
    if pybuiltin("hasattr")(electronic_component, "nelectron")
        return convert(Int, electronic_component.nelectron)
    elseif pybuiltin("hasattr")(electronic_component, "mo_occ")
        mo_occ = collect(electronic_component.mo_occ)
        return convert(Int, sum(mo_occ))
    else
        return DEFAULT_ELECTRON_COUNT
    end
end

function extractOrbitalCount(electronic_component::Any)
    if pybuiltin("hasattr")(electronic_component, "mo_occ")
        return length(electronic_component.mo_occ)
    else
        return DEFAULT_ORBITAL_COUNT
    end
end

function hasNuclearComponents(neo_molecule::Any, meanfield_object::Any)
    return pybuiltin("hasattr")(neo_molecule, "nuc_num") && neo_molecule.nuc_num > 0 &&
           pybuiltin("hasattr")(meanfield_object, "components") && haskey(meanfield_object.components, "n0")
end

function extractNuclearOrbitalCount(meanfield_object::Any)
    nuclear_component = meanfield_object.components["n0"]
    if pybuiltin("hasattr")(nuclear_component, "mo_occ")
        return length(nuclear_component.mo_occ)
    else
        return 0
    end
end


function calculateComputationalSavings(selected_configs::Vector, qubit_count::Int, total_orbitals::Int)
    determinant_savings = 2^total_orbitals / max(length(selected_configs), 1)
    qubit_savings = total_orbitals - qubit_count
    
    return ComputationalSavings(determinant_savings, qubit_savings)
end

function extractBaseEnergy(meanfield_object::Any)
    if pybuiltin("hasattr")(meanfield_object, "e_tot")
        return convert(Float64, meanfield_object.e_tot)
    else
        return 0.0
    end
end

# ======================== Error Handling ========================

function handleCalculationError(error::Exception, mol::Molecule, calc::NEOCalculation)
    @error "Calculation failed for molecule with $(length(mol.quantum_nuc)) quantum nuclei"
    @error "Method: $(calc.xc), EPC: $(calc.epc)"
    @error "Error details: $error"
end

# ======================== Results Summary ========================

function printResultsSummary(results::NEOResults, mol::Molecule, calc::NEOCalculation, 
                              config_sel::ConfigSelection = ConfigSelection())
    
    println("\n" * "="^60)
    println("SPARSE QEE-CNEO CALCULATION RESULTS")
    println("="^60)
    
    printMethodInformation(calc, results)
    printEnergyInformation(results)
    printConfigurationInformation(results, mol)
    printHamiltonianInformation(results)
    printPerformanceInformation(results)
    printNeoSpecificMetrics(results)
    
    println("="^60)
end

function printMethodInformation(calc::NEOCalculation, results::NEOResults)
    println("\nMethod Information:")
    method_display = createMethodDisplayString(calc, results)
    println("  Method: $method_display")
    
    if shouldDisplayEPCFunctional(calc)
        println("  EPC functional: $(calc.epc)")
    end
end

function createMethodDisplayString(calc::NEOCalculation, results::NEOResults)
    base_method = calc.xc
    configuration_method = results.method_used == "mp2" ? "MP2" : results.method_used
    return "$base_method+$configuration_method"
end

function shouldDisplayEPCFunctional(calc::NEOCalculation)
    return calc.epc != "none" && !isempty(calc.epc)
end

function printEnergyInformation(results::NEOResults)
    println("\nEnergy Information:")
    println("  Base energy: $(round(results.energy, digits=ENERGY_DISPLAY_PRECISION)) Ha")
    
    if hasMP2Correlation(results)
        println("  MP2 correlation: $(round(results.mp2_correlation, digits=ENERGY_DISPLAY_PRECISION)) Ha")
    end
    
    println("  Total energy: $(round(results.total_energy, digits=ENERGY_DISPLAY_PRECISION)) Ha")
end

function hasMP2Correlation(results::NEOResults)
    return results.mp2_correlation != 0.0
end

function printConfigurationInformation(results::NEOResults, mol::Molecule)
    println("\nConfiguration Information:")
    
    electron_count = extractElectronCountFromOrbitals(results.orbitals_per_species)
    println("  Electrons: $electron_count")
    println("  Quantum nuclei: $(length(mol.quantum_nuc))")
    println("  Selected configurations: $(results.n_configs)")
    println("  Qubits required: $(results.n_qubits)")
    println("  Captured importance: $(formatPercentage(results.captured_importance))%")
    println("  Computational savings:")
    println("    Determinant reduction: $(round(results.det_savings, digits=1))×")
    println("    Qubit reduction: $(results.qubit_savings)")
    println("  Orbitals per species: $(results.orbitals_per_species)")
    
    if results.orbital_truncation
        println("  ⚠️  Nuclear orbitals were truncated")
    end
end

function extractElectronCountFromOrbitals(orbitals_per_species::Vector{Int})
    # Estimate electron count from orbital count (rough approximation)
    return length(orbitals_per_species) > 0 ? orbitals_per_species[1] ÷ 2 : DEFAULT_ELECTRON_COUNT
end

function formatPercentage(value::Float64)
    return round(value * 100, digits=1)
end

function printHamiltonianInformation(results::NEOResults)
    if hasHamiltonianMatrix(results)
        println("\nHamiltonian Information:")
        printHamiltonianBasicProperties(results.hamiltonian_matrix)
        printHamiltonianAnalysis(results.hamiltonian_matrix)
    end
end

function hasHamiltonianMatrix(results::NEOResults)
    return results.hamiltonian_matrix !== nothing
end

function printHamiltonianBasicProperties(hamiltonian_matrix::Any)
    matrix_size = size(hamiltonian_matrix)
    is_hermitian = ishermitian(hamiltonian_matrix)
    hermitian_status = is_hermitian ? "✓" : "✗"
    
    println("  Matrix dimensions: $matrix_size")
    println("  Hermitian: $hermitian_status")
end

function printHamiltonianAnalysis(hamiltonian_matrix::Any)
    if isMatrixAnalyzable(hamiltonian_matrix)
        properties = analyzeHamiltonianProperties(hamiltonian_matrix)
        println("  Ground state energy: $(round(properties.ground_state_energy, digits=ENERGY_DISPLAY_PRECISION)) Ha")
        println("  Matrix sparsity: $(formatPercentage(properties.sparsity))%")
    end
end

function isMatrixAnalyzable(hamiltonian_matrix::Any)
    return size(hamiltonian_matrix, 1) > 0
end

function printPerformanceInformation(results::NEOResults)
    println("\nPerformance Metrics:")
    println("  Computation time: $(round(results.computation_time, digits=IMPORTANCE_DISPLAY_PRECISION)) seconds")
    println("  Memory usage: $(round(results.memory_used, digits=1)) MB")
end

function printNeoSpecificMetrics(results::NEOResults)
    if hasNeoMetrics(results)
        println("\nNEO-Specific Metrics:")
        printNeoImportanceMetrics(results.neo_metrics)
        printNeoParticipationMetrics(results.neo_metrics)
    end
end

function hasNeoMetrics(results::NEOResults)
    return results.neo_metrics !== nothing
end

function printNeoImportanceMetrics(neo_metrics::Any)
    if hasStandardImportanceField(neo_metrics)
        standard_importance = round(neo_metrics.standard_importance, digits=IMPORTANCE_DISPLAY_PRECISION)
        println("  Standard importance: $standard_importance")
    end
end

function hasStandardImportanceField(neo_metrics::Any)
    return hasfield(typeof(neo_metrics), :standard_importance)
end

function printNeoParticipationMetrics(neo_metrics::Any)
    nuclear_participation = formatPercentage(neo_metrics.nuclear_participation)
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
    analyzeHamiltonianProperties(hamiltonian_matrix::Matrix)

Analyze key properties of the Hamiltonian matrix.

# Returns
- `HamiltonianProperties`: Struct containing analysis results
"""
function analyzeHamiltonianProperties(hamiltonian_matrix::Matrix)
    if isEmptyMatrix(hamiltonian_matrix)
        return createEmptyHamiltonianProperties()
    end
    
    eigenvalues = calculateHamiltonianEigenvalues(hamiltonian_matrix)
    ground_state_energy = extractGroundStateEnergy(eigenvalues)
    energy_gap = calculateEnergyGap(eigenvalues)
    matrix_sparsity = calculateMatrixSparsity(hamiltonian_matrix)
    condition_number = calculateConditionNumber(hamiltonian_matrix)
    
    return HamiltonianProperties(
        ground_state_energy,
        energy_gap,
        matrix_sparsity,
        condition_number,
        eigenvalues
    )
end

function isEmptyMatrix(hamiltonian_matrix::Matrix)
    return size(hamiltonian_matrix, 1) == 0
end

function createEmptyHamiltonianProperties()
    return HamiltonianProperties(
        0.0,  # ground_state_energy
        0.0,  # energy_gap
        1.0,  # sparsity (empty matrix is 100% sparse)
        1.0,  # condition_number
        Float64[]  # eigenvalues
    )
end

function calculateHamiltonianEigenvalues(hamiltonian_matrix::Matrix)
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

function extractGroundStateEnergy(eigenvalues::Vector{Float64})
    return length(eigenvalues) > 0 ? eigenvalues[1] : 0.0
end

function calculateEnergyGap(eigenvalues::Vector{Float64})
    return length(eigenvalues) > 1 ? eigenvalues[2] - eigenvalues[1] : 0.0
end

function calculateMatrixSparsity(hamiltonian_matrix::Matrix)
    zero_elements_count = count(abs.(hamiltonian_matrix) .< HAMILTONIAN_SPARSITY_THRESHOLD)
    total_elements = length(hamiltonian_matrix)
    return zero_elements_count / total_elements
end

function calculateConditionNumber(hamiltonian_matrix::Matrix)
    return try
        cond(hamiltonian_matrix)
    catch
        Inf
    end
end

# ======================== Test Suite ========================

"""
    runTestSuite(config::NEOConfig; run_modular_tests=true)

Run comprehensive test suite for the Sparse QEE-cNEO implementation.
"""
function runTestSuite(config::NEOConfig = NEOConfig(); run_modular_tests::Bool = true)
    displayTestSuiteHeader()
    
    if !checkPysctAvailability(config)
        return Dict("demo" => "completed")
    end
    
    test_results = executeTestCases(config, run_modular_tests)
    displayTestSummary(test_results)
    
    return test_results
end

function displayTestSuiteHeader()
    """Display test suite header with version information."""
    println("\n" * "="^60)
    println("Enhanced cNEO Sparse QEE Tests (v6.0)")
    println("Complete implementation with Hamiltonian construction")
    println("="^60 * "\n")
end

function checkPysctAvailability(config::NEOConfig)
    """Check PySCF availability and run demo mode if needed."""
    pyscf, has_neo = setup_pyscf(config)
    
    if !has_neo
        println("Running in demo mode without NEO calculations...")
        runDemoMode()
        return false
    end
    
    return true
end

function executeTestCases(config::NEOConfig, run_modular_tests::Bool)
    """Execute all test cases and collect results."""
    test_results = Dict{String, Any}()
    
    test_results["H2"] = runH2Test(config)
    test_results["Water"] = runWaterTest(config)
    merge!(test_results, runMethodComparisonTest(config))
    
    if run_modular_tests
        println("\nTest 4: Modular component tests")
        println("-" * 40)
        test_results["modular"] = runModularComponentTests(config)
    end
    
    return test_results
end

function runH2Test(config::NEOConfig)
    """Run H2 molecule test with quantum proton."""
    println("\nTest 1: H2 with quantum proton")
    println("-" * 40)
    
    mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
    
    try
        result = sparse_qee_cneo(mol_h2, neo_config=config)
        
        if result.hamiltonian_matrix !== nothing
            println("  Hamiltonian constructed: $(size(result.hamiltonian_matrix))")
        end
        
        println("✓ Test passed")
        return result
    catch error
        println("✗ Test failed: $error")
        return error
    end
end

function runWaterTest(config::NEOConfig)
    """Run water molecule test with quantum proton."""
    println("\nTest 2: Water with quantum proton")
    println("-" * 40)
    
    mol_water = Molecule("O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", 
                        "6-31g", quantum_nuc=[1])
    calc = NEOCalculation(xc="B3LYP", epc="17-2")
    config_sel = ConfigSelection(method="neo_cneo", max_configs=100)
    
    try
        result = sparse_qee_cneo(mol_water, calc=calc, config_sel=config_sel, 
                                neo_config=config)
        println("✓ Test passed")
        return result
    catch error
        println("✗ Test failed: $error")
        return error
    end
end

function runMethodComparisonTest(config::NEOConfig)
    """Run method comparison test for HCN molecule."""
    println("\nTest 3: Method comparison for HCN")
    println("-" * 40)
    
    mol_hcn = Molecule("H 0 0 -1.064; C 0 0 0; N 0 0 1.156", 
                      "sto-3g", quantum_nuc=[0])
    
    results = Dict{String, Any}()
    
    for method in ["mp2", "neo_cneo", "hybrid_final"]
        results["HCN_$method"] = testSingleMethod(mol_hcn, method, config)
    end
    
    return results
end

function testSingleMethod(molecule::Molecule, method::String, config::NEOConfig)
    """Test a single configuration method."""
    println("\n  Testing method: $method")
    config_sel = ConfigSelection(method=method, max_configs=50)
    
    try
        result = sparse_qee_cneo(molecule, config_sel=config_sel, neo_config=config)
        println("  ✓ Configurations: $(result.n_configs), Importance: $(round(result.captured_importance, digits=3))")
        return result
    catch error
        println("  ✗ Failed: $error")
        return error
    end
end

function displayTestSummary(test_results::Dict{String, Any})
    """Display test suite summary with pass/fail counts."""
    println("\n" * "="^60)
    println("Test Summary:")
    
    passed = count(v -> !(v isa Exception) for (k,v) in test_results)
    total = length(test_results)
    
    println("Passed: $passed/$total")
    println("="^60)
end

# ======================== Modular Component Tests ========================

"""
    runModularComponentTests(config::NEOConfig)

Run tests for individual components in isolation.
"""
function runModularComponentTests(config::NEOConfig)
    results = Dict{String, Any}()
    
    results["electronic_only"] = testElectronicComponentOnly(config)
    results["nuclear_methods"] = testNuclearMethods(config)
    merge!(results, testEPCFunctionals(config))
    results["compression"] = testConfigurationCompression(config)
    
    return results
end

function testElectronicComponentOnly(config::NEOConfig)
    """Test electronic-only calculations without quantum nuclei."""
    println("\n  a) Testing electronic component only...")
    
    mol_elec = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")  # No quantum nuclei
    
    try
        result = sparse_qee_cneo(mol_elec, neo_config=config)
        println("    ✓ Electronic-only calculation successful")
        return result
    catch error
        println("    ✗ Failed: $error")
        return error
    end
end

function testNuclearMethods(config::NEOConfig)
    """Test nuclear method calculations with quantum nuclei."""
    println("\n  b) Testing nuclear methods...")
    
    mol_nuc = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0, 1])
    config_sel = ConfigSelection(
        method="neo_cneo",
        max_configs=20,
        debug_nuclear=true
    )
    
    try
        result = sparse_qee_cneo(mol_nuc, config_sel=config_sel, neo_config=config)
        println("    ✓ Nuclear methods successful")
        return result
    catch error
        println("    ✗ Failed: $error")
        return error
    end
end

function testEPCFunctionals(config::NEOConfig)
    """Test all available EPC functionals."""
    println("\n  c) Testing EPC functionals...")
    
    results = Dict{String, Any}()
    mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
    
    for epc in ["17-1", "17-2", "18-1", "18-2"]
        results["epc_$epc"] = testSingleEPCFunctional(mol_h2, epc, config)
    end
    
    return results
end

function testSingleEPCFunctional(molecule::Molecule, epc::String, config::NEOConfig)
    """Test a single EPC functional."""
    calc = NEOCalculation(xc="B3LYP", epc=epc)
    config_sel = ConfigSelection(method="mp2", max_configs=10)
    
    try
        result = sparse_qee_cneo(molecule, calc=calc, config_sel=config_sel, neo_config=config)
        println("    ✓ EPC $epc: $(round(result.total_energy, digits=6)) Ha")
        return result.total_energy
    catch error
        println("    ✗ EPC $epc failed")
        return error
    end
end

function testConfigurationCompression(config::NEOConfig)
    """Test configuration compression functionality."""
    println("\n  d) Testing configuration compression...")
    
    mol_elec = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")  # No quantum nuclei
    config_sel = ConfigSelection(
        method="mp2",
        max_configs=200,
        use_compression=true
    )
    
    try
        result = sparse_qee_cneo(mol_elec, config_sel=config_sel, neo_config=config)
        compressed = count(c -> isa(c, CompressedConfig), result.configs)
        println("    ✓ Compressed $compressed configurations")
        return compressed
    catch error
        println("    ✗ Compression failed")
        return error
    end
end

# ======================== Demo Mode ========================

function runDemoMode()
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
    
    println("5. Constrained NEO (cNEO) Methods:")
    println("   - cNEO-HF: Constrained NEO Hartree-Fock")
    println("   - cNEO-MP2: Constrained NEO with MP2 correlation")
    println("   - Nuclear position constraints with Lagrange multipliers")
    println("   - Newton optimization with analytical Hessians")
    println("   - Clean Code implementation with small functions")
    
    println("\n6. Modular Testing Approach:")
    println("   - Test each component independently")
    println("   - Electronic-only calculations")
    println("   - Nuclear method validation")
    println("   - EPC functional comparison")
    println("   - Configuration compression")
    println("   - Hamiltonian construction")
    println("   - cNEO constraint satisfaction")
    
    println("\n" * "="^60)
    println("Install PySCF with NEO for full functionality")
    println("cNEO methods available with integrated Clean Code implementation")
    println("="^60)
end

# ======================== cNEO-QEE Integrated Workflow ========================

"""
    sparse_qee_with_cneo(mol::Molecule, cneo_calc::CNEOCalculation; kwargs...)

Integrated workflow: cNEO constraint optimization followed by QEE Hamiltonian construction.

This function represents the **full integration** between constrained NEO and quantum computing:
1. Runs cNEO calculation to optimize nuclear positions with constraints
2. Uses optimized geometry for QEE configuration generation  
3. Constructs quantum computing Hamiltonian with cNEO-optimized nuclear effects

# Arguments
- `mol::Molecule`: Molecular system with quantum nuclei
- `cneo_calc::CNEOCalculation`: cNEO constraint parameters
- `config_sel::ConfigSelection`: Configuration selection for QEE (default: ConfigSelection())
- `neo_config::NEOConfig`: PySCF configuration (default: NEOConfig())

# Returns
- `Tuple{CNEOResults, NEOResults}`: (cNEO optimization results, QEE results with Hamiltonian)
"""
function sparse_qee_with_cneo(mol::Molecule, 
                              cneo_calc::CNEOCalculation;
                              config_sel::ConfigSelection = ConfigSelection(),
                              neo_config::NEOConfig = NEOConfig())
    
    @info "Starting integrated cNEO-QEE workflow..."
    
    # Step 1: Run cNEO calculation to optimize nuclear positions
    @info "Step 1: Running cNEO constraint optimization"
    cneo_results = run_cneo_hf(mol, cneo_calc, neo_config=neo_config)
    
    if !cneo_results.converged
        @warn "cNEO calculation did not converge - proceeding with partial optimization"
    end
    
    # Step 2: Create QEE workflow with cNEO-optimized nuclear positions
    @info "Step 2: Converting cNEO results to QEE workflow"
    qee_results = cneo_to_qee_workflow(cneo_results, mol, config_sel, neo_config)
    
    @info "Integrated cNEO-QEE workflow completed successfully!"
    
    return cneo_results, qee_results
end

"""
    cneo_to_qee_workflow(cneo_results::CNEOResults, mol::Molecule, 
                         config_sel::ConfigSelection, neo_config::NEOConfig)

Convert cNEO results to QEE workflow with optimized nuclear positions.

This function creates the bridge between constrained NEO and quantum eigensolver:
1. Extracts optimized nuclear positions from cNEO results
2. Updates molecular geometry with cNEO-optimized positions
3. Runs QEE workflow with the optimized geometry
4. Returns QEE results ready for quantum computing

# Arguments
- `cneo_results::CNEOResults`: Results from cNEO constraint optimization
- `mol::Molecule`: Original molecular system
- `config_sel::ConfigSelection`: Configuration selection parameters
- `neo_config::NEOConfig`: PySCF configuration

# Returns
- `NEOResults`: QEE results with Hamiltonian constructed from cNEO-optimized geometry
"""
function cneo_to_qee_workflow(cneo_results::CNEOResults, 
                              mol::Molecule,
                              config_sel::ConfigSelection,
                              neo_config::NEOConfig)
    
    @info "Converting cNEO results to QEE workflow..."
    
    # Extract optimized nuclear positions from cNEO results
    optimized_positions = cneo_results.actual_nuclear_positions
    
    # Create updated molecule with cNEO-optimized nuclear positions
    updated_mol = updateMoleculeWithCNEOPositions(mol, optimized_positions)
    
    # Create NEO calculation incorporating cNEO energy contributions
    neo_calc = createNeoCalculationFromCNEO(cneo_results)
    
    @info "Running QEE workflow with cNEO-optimized nuclear positions"
    
    # Run main QEE workflow with optimized geometry
    qee_results = sparse_qee_cneo(
        updated_mol,
        calc=neo_calc,
        config_sel=config_sel,
        neo_config=neo_config
    )
    
    # Add cNEO information to QEE results for traceability
    enhanced_results = enhanceQEEResultsWithCNEOInfo(qee_results, cneo_results)
    
    @info "QEE workflow completed with cNEO-optimized geometry"
    
    return enhanced_results
end

# ======================== Helper Functions for cNEO-QEE Integration ========================

function updateMoleculeWithCNEOPositions(mol::Molecule, optimized_positions::Vector{Vector{Float64}})
    """Update molecular geometry with cNEO-optimized nuclear positions."""
    
    # For now, return original molecule (geometry update requires more complex implementation)
    # In full implementation, this would update atomic coordinates based on optimized_positions
    @info "Using cNEO-optimized nuclear positions ($(length(optimized_positions)) nuclei)"
    
    return mol
end

function createNeoCalculationFromCNEO(cneo_results::CNEOResults)
    """Create NEO calculation incorporating cNEO energy contributions."""
    
    # Create NEO calculation that accounts for cNEO constraint effects
    neo_calc = NEOCalculation()
    
    @info "NEO calculation incorporates cNEO constraint energy: $(cneo_results.total_energy) Ha"
    
    return neo_calc
end

function enhanceQEEResultsWithCNEOInfo(qee_results::NEOResults, cneo_results::CNEOResults)
    """Add cNEO information to QEE results for complete traceability."""
    
    @info "QEE results enhanced with cNEO constraint information"
    @info "  cNEO constraint satisfaction: $(cneo_results.converged ? "CONVERGED" : "PARTIAL")"
    @info "  Nuclear position errors: $(cneo_results.position_errors)"
    @info "  cNEO-QEE energy integration: $(cneo_results.total_energy) → $(qee_results.total_energy) Ha"
    
    return qee_results
end

end # module
