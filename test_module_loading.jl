#!/usr/bin/env julia
"""
SparseQEEcNEO Module Loading Test

Tests individual components of SparseQEEcNEO to verify module loading.
Follows Clean Code principles with small, focused functions.
"""

# Constants
const MODULE_FILES = [
    "types.jl",
    "pyscf_interface.jl", 
    "epc_functionals.jl",
    "nuclear_methods.jl",
    "configuration_generation.jl",
    "importance_analysis.jl",
    "qee_methods.jl",
    "hamiltonian_construction.jl"
]

const SUCCESS_SYMBOL = "✓"
const ERROR_SYMBOL = "✗"
const SOURCE_DIR = "src/"

# Functions - each does one thing and has a descriptive name
function display_test_header()
    println("=== SparseQEEcNEO Module Loading Test ===")
end

function test_basic_dependencies()
    println("\n1. Testing basic dependencies...")
    try
        @eval using PyCall, LinearAlgebra, SparseArrays, Printf, Statistics, HDF5, JSON
        println("   $SUCCESS_SYMBOL All basic dependencies loaded")
        return true
    catch error
        println("   $ERROR_SYMBOL Error loading dependencies: $error")
        return false
    end
end

function load_single_module(module_filename)
    try
        include(SOURCE_DIR * module_filename)
        println("   $SUCCESS_SYMBOL $module_filename loaded successfully")
        return true
    catch error
        println("   $ERROR_SYMBOL Error loading $module_filename: $error")
        return false
    end
end

function test_module_file_loading()
    println("\n2. Testing individual module files...")
    
    for module_file in MODULE_FILES
        if !load_single_module(module_file)
            return false
        end
    end
    return true
end

function test_module_accessibility()
    println("\n3. Testing module accessibility...")
    
    module_tests = [
        (:Types, "Types module accessible"),
        (:PySCFInterface, "PySCFInterface module accessible"),
        (:EPCFunctionals, "EPCFunctionals module accessible"),
        (:NuclearMethods, "NuclearMethods module accessible"),
        (:ConfigurationGeneration, "ConfigurationGeneration module accessible"),
        (:ImportanceAnalysis, "ImportanceAnalysis module accessible"),
        (:QEEMethods, "QEEMethods module accessible"),
        (:HamiltonianConstruction, "HamiltonianConstruction module accessible")
    ]
    
    try
        for (module_symbol, message) in module_tests
            @eval using .$module_symbol
            println("   $SUCCESS_SYMBOL $message")
        end
        return true
    catch error
        println("   $ERROR_SYMBOL Error using modules: $error")
        return false
    end
end

function create_test_molecule()
    return Molecule(
        "H 0 0 0; H 0 0 1.4",
        "sto-3g",
        quantum_nuc = [1]
    )
end

function test_type_creation()
    println("\n4. Testing type creation...")
    
    try
        molecule = create_test_molecule()
        println("   $SUCCESS_SYMBOL Molecule type created successfully")
        
        neo_calculation = NEOCalculation()
        println("   $SUCCESS_SYMBOL NEOCalculation type created successfully")
        
        configuration_selection = ConfigSelection()
        println("   $SUCCESS_SYMBOL ConfigSelection type created successfully")
        
        neo_config = NEOConfig()
        println("   $SUCCESS_SYMBOL NEOConfig type created successfully")
        
        return true
    catch error
        println("   $ERROR_SYMBOL Error creating types: $error")
        return false
    end
end

function display_test_results(all_tests_passed)
    if all_tests_passed
        println("\n5. All tests passed! $SUCCESS_SYMBOL")
        println("Individual components work correctly.")
        println("Module loading is functioning properly.")
    else
        println("\n5. Some tests failed! $ERROR_SYMBOL")
        exit(1)
    end
end

# Main execution
function run_module_loading_tests()
    display_test_header()
    
    dependencies_ok = test_basic_dependencies()
    if !dependencies_ok
        exit(1)
    end
    
    modules_ok = test_module_file_loading()
    if !modules_ok
        exit(1)
    end
    
    accessibility_ok = test_module_accessibility()
    if !accessibility_ok
        exit(1)
    end
    
    types_ok = test_type_creation()
    
    all_tests_passed = dependencies_ok && modules_ok && accessibility_ok && types_ok
    display_test_results(all_tests_passed)
end

# Execute tests
run_module_loading_tests()
