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

# Note: cNEO functionality available as separate modules in src/cneo_hf.jl and src/cneo_mp2.jl

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
    
    start_time = time()
    initial_memory = Base.gc_live_bytes() / 1024^2
    
    # Setup PySCF
    pyscf, has_neo = setup_pyscf(neo_config)
    
    if !has_neo
        error("NEO module not available in PySCF. Please install PySCF with NEO support.")
    end
    
    # Build NEO molecule
    mol_neo = build_neo_molecule(mol, pyscf)
    
    # Run mean-field calculation
    mf = run_neo_meanfield(mol_neo, calc, pyscf)
    
    # Apply nuclear orbital truncation if needed
    truncated = false
    truncation_info = Dict{String, Any}()
    
    if config_sel.max_nuc_orbs > 0
        truncated, truncation_info = apply_nuclear_orbital_truncation!(mf, mol_neo, config_sel)
    end
    
    # Calculate MP2 if needed
    mp2_correlation = 0.0
    t2_amplitudes = nothing
    
    if config_sel.method in ["mp2", "mp2_enhanced"]
        if !truncated || config_sel.force_mp2_with_truncation
            try
                mp, ecorr, t2 = run_neo_mp2(mf, mol_neo)
                mp2_correlation = ecorr
                t2_amplitudes = extract_t2_amplitudes(t2, config_sel)
            catch e
                if truncated
                    @warn "NEO-MP2 failed with truncated orbitals (expected): $e"
                    @info "Continuing without MP2 correlation"
                else
                    @warn "NEO-MP2 failed: $e"
                end
            end
        else
            @info "Skipping NEO-MP2 due to orbital truncation"
        end
    end
    
    # Generate configurations based on method
    configs = generate_configurations(mf, mol_neo, config_sel, t2_amplitudes)
    
    # Calculate importance metrics
    importance_data = calculate_importance_metrics(configs, mol_neo, config_sel)
    
    # Auto-switch method if importance too low
    if config_sel.auto_switch_method && importance_data.total_importance < config_sel.importance_threshold_switch
        @warn "Low importance ($(round(importance_data.total_importance, digits=3))), switching methods"
        
        # Try alternative methods
        for alt_method in ["neo_cneo", "neo_enhanced", "hybrid_final"]
            if alt_method != config_sel.method
                config_sel_new = deepcopy(config_sel)
                config_sel_new.method = alt_method
                configs_new = generate_configurations(mf, mol_neo, config_sel_new, t2_amplitudes)
                importance_new = calculate_importance_metrics(configs_new, mol_neo, config_sel_new)
                
                if importance_new.total_importance > importance_data.total_importance
                    configs = configs_new
                    importance_data = importance_new
                    config_sel = config_sel_new
                    @info "Switched to method: $alt_method"
                    break
                end
            end
        end
    end
    
    # Select configurations based on importance
    selected_configs, n_qubits = select_important_configurations(configs, config_sel)
    
    # Construct second-quantized Hamiltonian
    ham_data, H_matrix = construct_second_quantized_hamiltonian(
        mf, mol_neo, selected_configs, config_sel
    )
    
    # Apply EPC correction if requested
    epc_energy = 0.0
    if calc.epc != "none" && calc.epc != ""
        epc_energy = calculate_epc_energy(mf, mol_neo, calc.epc)
    end
    
    # Calculate final metrics
    elapsed_time = time() - start_time
    memory_used = Base.gc_live_bytes() / 1024^2 - initial_memory
    
    # Determine electron count and savings
    n_electrons = 2  # Default
    n_orbitals = 4   # Default
    
    if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "e")
        mf_e = mf.components["e"]
        if pybuiltin("hasattr")(mf_e, "nelectron")
            n_electrons = convert(Int, mf_e.nelectron)
        elseif pybuiltin("hasattr")(mf_e, "mo_occ")
            mo_occ = collect(mf_e.mo_occ)
            n_electrons = convert(Int, sum(mo_occ))
        end
        if pybuiltin("hasattr")(mf_e, "mo_occ")
            n_orbitals = length(mf_e.mo_occ)
        end
    end
    
    n_qubits_full = n_orbitals
    det_savings = 2^n_qubits_full / max(length(selected_configs), 1)
    qubit_savings = n_qubits_full - n_qubits
    
    # Get orbitals per species
    orbitals_per_species = Int[n_orbitals]  # Default just electronic
    if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
        if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "n0")
            mf_n = mf.components["n0"]
            if pybuiltin("hasattr")(mf_n, "mo_occ")
                n_nuc_orbs = length(mf_n.mo_occ)
                orbitals_per_species = Int[n_orbitals, n_nuc_orbs]
            end
        end
    end
    
    # Create results - use positional arguments
    results = NEOResults(
        selected_configs,
        pybuiltin("hasattr")(mf, "e_tot") ? convert(Float64, mf.e_tot) : 0.0,
        mp2_correlation,
        (pybuiltin("hasattr")(mf, "e_tot") ? convert(Float64, mf.e_tot) : 0.0) + mp2_correlation + epc_energy,
        length(selected_configs),
        n_qubits,
        importance_data.total_importance,
        det_savings,
        qubit_savings,
        elapsed_time,
        memory_used,
        config_sel.method,
        orbitals_per_species,
        importance_data.neo_metrics,
        truncated,
        ham_data,
        H_matrix
    )
    
    # Print summary
    print_results_summary(results, mol, calc, config_sel, n_electrons)
    
    return results
end

# ======================== Results Summary ========================

function print_results_summary(results::NEOResults, mol::Molecule, calc::NEOCalculation, 
                              config_sel::ConfigSelection, n_electrons::Int)
    println("\nResults:")
    println("  Method: $(calc.xc)+$(results.method_used == "mp2" ? "MP2" : results.method_used)")
    if calc.epc != "none"
        println("  EPC functional: $(calc.epc)")
    end
    println("  Energy: $(round(results.energy, digits=6)) Ha")
    if results.mp2_correlation != 0.0
        println("  MP2 correlation: $(round(results.mp2_correlation, digits=6)) Ha")
    end
    println("  Total energy: $(round(results.total_energy, digits=6)) Ha")
    println("  Electrons: $n_electrons")
    println("  Selected configs: $(results.n_configs)")
    println("  Qubits needed: $(results.n_qubits)")
    println("  Captured importance: $(round(results.captured_importance * 100, digits=1))%")
    println("  Determinant savings: $(round(results.det_savings, digits=1))x")
    println("  Qubit savings: $(results.qubit_savings)")
    println("  Orbitals per species: $(results.orbitals_per_species)")
    if results.orbital_truncation
        println("  * Nuclear orbitals were truncated")
    end
    
    # Hamiltonian info
    if results.hamiltonian_matrix !== nothing
        println("\nHamiltonian:")
        println("  Matrix size: $(size(results.hamiltonian_matrix))")
        println("  Hermitian: $(ishermitian(results.hamiltonian_matrix) ? "✓" : "✗")")
        
        # Quick analysis
        if size(results.hamiltonian_matrix, 1) > 0
            props = analyze_hamiltonian_properties(results.hamiltonian_matrix)
            println("  Ground state: $(round(props.ground_state_energy, digits=6)) Ha")
            println("  Sparsity: $(round(props.sparsity * 100, digits=1))%")
        end
    end
    
    println("\nPerformance:")
    println("  Computation time: $(round(results.computation_time, digits=3)) seconds")
    println("  Memory used: $(round(results.memory_used, digits=1)) MB")
    
    if results.neo_metrics !== nothing
        println("\nNEO-specific metrics:")
        if hasfield(typeof(results.neo_metrics), :standard_importance)
            println("  Standard importance: $(round(results.neo_metrics.standard_importance, digits=3))")
        end
        println("  Nuclear participation: $(round(results.neo_metrics.nuclear_participation * 100, digits=1))%")
        println("  Coupling contribution: $(round(results.neo_metrics.coupling_contribution, digits=3))")
    end
end

# ======================== Hamiltonian Analysis ========================

"""
    analyze_hamiltonian_properties(H::Matrix)

Analyze properties of the Hamiltonian matrix.
"""
function analyze_hamiltonian_properties(H::Matrix)
    # Handle empty matrix
    if size(H, 1) == 0
        return (
            ground_state_energy = 0.0,
            energy_gap = 0.0,
            sparsity = 1.0,
            condition_number = 1.0,
            eigenvalues = Float64[]
        )
    end
    
    # Calculate eigenvalues
    eigenvals = try
        eigvals(Hermitian(H))
    catch
        eigvals(H)
    end
    
    # Sort eigenvalues
    sort!(eigenvals)
    
    # Ground state energy
    ground_state = length(eigenvals) > 0 ? eigenvals[1] : 0.0
    
    # Energy gap
    gap = length(eigenvals) > 1 ? eigenvals[2] - eigenvals[1] : 0.0
    
    # Sparsity
    n_zero = count(abs.(H) .< 1e-12)
    sparsity = n_zero / length(H)
    
    # Condition number
    cond_num = try
        cond(H)
    catch
        Inf
    end
    
    return (
        ground_state_energy = ground_state,
        energy_gap = gap,
        sparsity = sparsity,
        condition_number = cond_num,
        eigenvalues = eigenvals
    )
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
