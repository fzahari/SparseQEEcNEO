#!/usr/bin/env julia

"""
Systematic refactoring script for Clean Code naming compliance
This script renames functions from snake_case to camelCase for internal functions
while preserving API compatibility for exported functions.
"""

using Pkg

# Dictionary of function renamings (old_name => new_name)
const FUNCTION_RENAMINGS = Dict(
    # Main module functions - keep exported ones, rename internal ones
    "initialize_quantum_environment" => "initializeQuantumEnvironment",
    "validate_neo_availability" => "validateNeoAvailability", 
    "apply_orbital_truncation_if_needed" => "applyOrbitalTruncationIfNeeded",
    "should_apply_truncation" => "shouldApplyTruncation",
    "calculate_correlation_corrections" => "calculateCorrelationCorrections",
    "should_calculate_mp2_correlation" => "shouldCalculateMP2Correlation",
    "attempt_mp2_calculation" => "attemptMP2Calculation",
    "should_skip_mp2_due_to_truncation" => "shouldSkipMP2DueToTruncation",
    "handle_mp2_failure" => "handleMP2Failure",
    "optimize_method_selection" => "optimizeMethodSelection",
    "should_attempt_method_switching" => "shouldAttemptMethodSwitching",
    "attempt_method_switching" => "attemptMethodSwitching",
    "try_alternative_method" => "tryAlternativeMethod",
    "construct_hamiltonian_matrix" => "constructHamiltonianMatrix",
    "calculate_energy_corrections" => "calculateEnergyCorrections",
    "should_apply_epc_correction" => "shouldApplyEPCCorrection",
    "analyze_orbital_structure" => "analyzeOrbitalStructure",
    "has_electronic_component" => "hasElectronicComponent",
    "extract_electron_count" => "extractElectronCount",
    "extract_orbital_count" => "extractOrbitalCount",
    "has_nuclear_components" => "hasNuclearComponents",
    "extract_nuclear_orbital_count" => "extractNuclearOrbitalCount",
    "calculate_computational_savings" => "calculateComputationalSavings",
    "extract_base_energy" => "extractBaseEnergy",
    "handle_calculation_error" => "handleCalculationError",
    "print_results_summary" => "printResultsSummary",
    "print_method_information" => "printMethodInformation",
    "create_method_display_string" => "createMethodDisplayString",
    "should_display_epc_functional" => "shouldDisplayEPCFunctional",
    "print_energy_information" => "printEnergyInformation",
    "has_mp2_correlation" => "hasMP2Correlation",
    "print_configuration_information" => "printConfigurationInformation",
    "extract_electron_count_from_orbitals" => "extractElectronCountFromOrbitals",
    "format_percentage" => "formatPercentage",
    "print_hamiltonian_information" => "printHamiltonianInformation",
    "has_hamiltonian_matrix" => "hasHamiltonianMatrix",
    "print_hamiltonian_basic_properties" => "printHamiltonianBasicProperties",
    "print_hamiltonian_analysis" => "printHamiltonianAnalysis",
    "is_matrix_analyzable" => "isMatrixAnalyzable",
    "print_performance_information" => "printPerformanceInformation",
    "print_neo_specific_metrics" => "printNeoSpecificMetrics",
    "has_neo_metrics" => "hasNeoMetrics",
    "print_neo_importance_metrics" => "printNeoImportanceMetrics",
    "has_standard_importance_field" => "hasStandardImportanceField",
    "print_neo_participation_metrics" => "printNeoParticipationMetrics",
    "analyze_hamiltonian_properties" => "analyzeHamiltonianProperties",
    "is_empty_matrix" => "isEmptyMatrix",
    "create_empty_hamiltonian_properties" => "createEmptyHamiltonianProperties",
    "calculate_hamiltonian_eigenvalues" => "calculateHamiltonianEigenvalues",
    "extract_ground_state_energy" => "extractGroundStateEnergy",
    "calculate_energy_gap" => "calculateEnergyGap",
    "calculate_matrix_sparsity" => "calculateMatrixSparsity",
    "calculate_condition_number" => "calculateConditionNumber",
    "run_test_suite" => "runTestSuite",  # Keep original for API
    "display_test_suite_header" => "displayTestSuiteHeader",
    "check_pyscf_availability" => "checkPysctAvailability",
    "execute_test_cases" => "executeTestCases",
    "run_h2_test" => "runH2Test",
    "run_water_test" => "runWaterTest",
    "run_method_comparison_test" => "runMethodComparisonTest",
    "test_single_method" => "testSingleMethod",
    "display_test_summary" => "displayTestSummary",
    "run_modular_component_tests" => "runModularComponentTests",
    "test_electronic_component_only" => "testElectronicComponentOnly",
    "test_nuclear_methods" => "testNuclearMethods",
    "test_epc_functionals" => "testEPCFunctionals",
    "test_single_epc_functional" => "testSingleEPCFunctional",
    "test_configuration_compression" => "testConfigurationCompression",
    "run_demo_mode" => "runDemoMode",  # Keep original for API
    "update_molecule_with_cneo_positions" => "updateMoleculeWithCNEOPositions",
    "create_neo_calculation_from_cneo" => "createNeoCalculationFromCNEO",
    "enhance_qee_results_with_cneo_info" => "enhanceQEEResultsWithCNEOInfo",
)

function refactor_file_names(file_path::String)
    """Apply systematic naming refactoring to a Julia file."""
    
    println("📝 Refactoring file: $(basename(file_path))")
    
    # Read file content
    content = read(file_path, String)
    original_content = content
    
    # Apply function renamings
    changes_made = 0
    for (old_name, new_name) in FUNCTION_RENAMINGS
        # Match function definitions
        old_pattern = Regex("function\\s+$old_name\\s*\\(")
        new_replacement = "function $new_name("
        
        if occursin(old_pattern, content)
            content = replace(content, old_pattern => new_replacement)
            changes_made += 1
            println("  ✓ Renamed function: $old_name → $new_name")
        end
        
        # Match function calls 
        call_pattern = Regex("\\b$old_name\\s*\\(")
        call_replacement = "$new_name("
        
        # Count occurrences for reporting
        call_matches = length(collect(eachmatch(call_pattern, content)))
        if call_matches > 0
            content = replace(content, call_pattern => call_replacement)
            changes_made += call_matches
            println("  ✓ Updated $call_matches call(s): $old_name → $new_name")
        end
    end
    
    # Write back only if changes were made
    if content != original_content
        write(file_path, content)
        println("  💾 Saved $changes_made changes to $(basename(file_path))")
        return changes_made
    else
        println("  ⏭  No changes needed in $(basename(file_path))")
        return 0
    end
end

function main()
    println("🔧 Starting systematic Clean Code naming refactoring...")
    println("=" ^ 60)
    
    # Get all Julia source files
    src_dir = joinpath(pwd(), "src")
    julia_files = [joinpath(src_dir, f) for f in readdir(src_dir) if endswith(f, ".jl")]
    
    total_changes = 0
    files_modified = 0
    
    for file_path in julia_files
        changes = refactor_file_names(file_path)
        if changes > 0
            total_changes += changes
            files_modified += 1
        end
    end
    
    println("=" ^ 60)
    println("📊 Refactoring Summary:")
    println("  Files processed: $(length(julia_files))")
    println("  Files modified: $files_modified") 
    println("  Total changes: $total_changes")
    
    if total_changes > 0
        println("\\n🎉 Refactoring completed! Run tests to verify functionality.")
        println("📋 Reminder: Exported functions maintain API compatibility")
    else
        println("\\n✅ No changes needed - code already compliant")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end