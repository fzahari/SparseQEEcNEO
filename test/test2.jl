# test2.jl - Comprehensive tests with fixes

ENV["PYTHONPATH"] = "/Users/federicozahariev/Work/Programs/QEE_Split_Grouping/pyscf-master"
push!(LOAD_PATH, dirname(@__DIR__))
using SparseQEEcNEO
using LinearAlgebra  # For eigvals

config = NEOConfig(pyscf_path="/Users/federicozahariev/Work/Programs/QEE_Split_Grouping/pyscf-master")

println("="^80)
println("Comprehensive SparseQEEcNEO Tests")
println("="^80)

# Test 1: H2 with different methods
println("\n1. Testing different configuration selection methods for H2")
println("-"^60)

mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
calc = NEOCalculation(xc="HF", epc="none")

methods = ["mp2", "neo_cneo", "neo_enhanced", "hybrid_final"]
results_by_method = Dict()

for method in methods
    print("Testing $method... ")
    config_sel = ConfigSelection(
        method=method,
        max_configs=100,
        max_nuc_orbs=0  # Skip truncation
    )
    
    try
        results = sparse_qee_cneo(mol_h2, calc=calc, config_sel=config_sel, neo_config=config)
        results_by_method[method] = results
        println("✓ E = $(round(results.total_energy, digits=6)) Ha, $(results.n_configs) configs")
    catch e
        println("✗ Failed: $e")
    end
end

# Test 2: EPC functionals (skip "none" for B3LYP)
println("\n2. Testing EPC functionals")
println("-"^60)

epc_results = Dict()
for epc in ["17-1", "17-2", "18-1", "18-2"]
    print("Testing EPC $epc with B3LYP... ")
    calc_epc = NEOCalculation(xc="B3LYP", epc=epc)
    config_sel = ConfigSelection(method="neo_cneo", max_configs=50, max_nuc_orbs=0)
    
    try
        results = sparse_qee_cneo(mol_h2, calc=calc_epc, config_sel=config_sel, neo_config=config)
        epc_results[epc] = results.total_energy
        println("✓ E = $(round(results.total_energy, digits=6)) Ha")
    catch e
        println("✗ Failed: $e")
    end
end

# Test 3: Water molecule with quantum proton
println("\n3. Testing water with quantum proton")
println("-"^60)

mol_water = Molecule("O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", "6-31g", quantum_nuc=[1])
calc_water = NEOCalculation(xc="B3LYP", epc="17-2")
config_sel_water = ConfigSelection(
    method="neo_cneo",
    max_configs=100,
    max_nuc_orbs=0,
    include_doubles=true
)

try
    results_water = sparse_qee_cneo(mol_water, calc=calc_water, config_sel=config_sel_water, neo_config=config)
    println("✓ Water results:")
    println("  Energy: $(round(results_water.total_energy, digits=6)) Ha")
    println("  Configurations: $(results_water.n_configs)")
    println("  Importance: $(round(results_water.captured_importance * 100, digits=1))%")
catch e
    println("✗ Failed: $e")
end

# Test 4: Hamiltonian properties
println("\n4. Analyzing Hamiltonian properties")
println("-"^60)

if haskey(results_by_method, "neo_cneo")
    results = results_by_method["neo_cneo"]
    if results.hamiltonian_matrix !== nothing
        H = results.hamiltonian_matrix
        eigenvals = eigvals(Hermitian(H))
        println("NEO-CNEO Hamiltonian analysis:")
        println("  Size: $(size(H))")
        println("  First 5 eigenvalues: $(eigenvals[1:min(5,end)])")
        println("  Ground state: $(minimum(eigenvals))")
        println("  Gap: $(length(eigenvals) > 1 ? eigenvals[2] - eigenvals[1] : 0.0)")
        
        # Save Hamiltonian
        filename = "h2_neo_cneo_hamiltonian.h5"
        save_hamiltonian(filename, results.hamiltonian_data, H)
        println("  Saved to: $filename")
    end
end

# Test 5: Configuration analysis
println("\n5. Configuration analysis")
println("-"^60)

for (method, results) in results_by_method
    if results isa NEOResults
        println("\n$method configurations:")
        for (i, config) in enumerate(results.configs[1:min(5, end)])
            println("  $i. $(config.name) - weight: $(round(config.weight, digits=4))")
        end
    end
end

# Summary of EPC energies
if length(epc_results) > 0
    println("\n6. EPC Energy Summary")
    println("-"^60)
    println("EPC Type | Total Energy (Ha)")
    println("---------|-----------------")
    for (epc, energy) in sort(collect(epc_results), by=x->x[2])
        println("$(rpad(epc,8)) | $(round(energy, digits=6))")
    end
end

println("\n" * "="^80)
println("Tests completed!")
println("="^80)
