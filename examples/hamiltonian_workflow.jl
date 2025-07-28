#!/usr/bin/env julia
"""
examples/hamiltonian_workflow.jl - Complete workflow demonstrating Hamiltonian construction
"""

# Setup
ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/path/to/pyscf-master")
push!(LOAD_PATH, dirname(@__DIR__))

using SparseQEEcNEO
using LinearAlgebra
using Printf
using Statistics

config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))

println("Sparse QEE-cNEO Hamiltonian Workflow Examples")
println("="^70)

# ======================== Example 1: Basic Hamiltonian Construction ========================
println("\nExample 1: Basic H2 Hamiltonian construction")
println("-"^50)

# Create molecule
mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")

# Run calculation
println("Running calculation...")
config_sel = ConfigSelection(method="mp2", max_configs=10, max_nuc_orbs=0)
results = sparse_qee_cneo(mol_h2, config_sel=config_sel, neo_config=config)

# Access Hamiltonian
H = results.hamiltonian_matrix
ham_data = results.hamiltonian_data

println("\nHamiltonian properties:")
println("  Matrix size: $(size(H))")
println("  Number of configurations: $(length(ham_data.configs))")
println("  Active orbitals: $(ham_data.active_orbitals)")
println("  Frozen orbitals: $(ham_data.frozen_orbitals)")

# ======================== Example 2: NEO Hamiltonian with Coupling ========================
println("\n\nExample 2: NEO Hamiltonian with electron-nuclear coupling")
println("-"^50)

# Water with quantum proton
mol_water = Molecule(
    "O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0",
    "6-31g",
    quantum_nuc=[1]  # Second H is quantum
)

calc = NEOCalculation(xc="B3LYP", epc="17-2")
config_sel = ConfigSelection(method="neo_cneo", max_configs=100, max_nuc_orbs=0)

println("Running NEO calculation...")
results_neo = sparse_qee_cneo(mol_water, calc=calc, config_sel=config_sel, neo_config=config)

# Analyze coupling
ham_neo = results_neo.hamiltonian_data
println("\nCoupling analysis:")
for ((c1, c2), V) in ham_neo.coupling
    println("  Coupling $c1-$c2: size $(size(V))")
    println("    Max element: $(round(maximum(abs.(V)), digits=6))")
    println("    Mean element: $(round(mean(abs.(V)), digits=6))")
end

# ======================== Example 3: Hamiltonian Matrix Analysis ========================
println("\n\nExample 3: Detailed Hamiltonian matrix analysis")
println("-"^50)

function analyze_hamiltonian_detailed(H::Matrix, configs)
    println("\nMatrix structure analysis:")
    
    # Sparsity
    n_zero = count(abs.(H) .< 1e-12)
    sparsity = n_zero / length(H) * 100
    println("  Sparsity: $(round(sparsity, digits=1))%")
    
    # Block structure
    println("\n  Block structure (by configuration type):")
    
    # Group configurations by type
    ref_indices = Int[]
    single_indices = Int[]
    double_indices = Int[]
    nuclear_indices = Int[]
    coupled_indices = Int[]
    
    for (i, config) in enumerate(configs)
        if occursin("HF", config.name) || occursin("REF", config.name)
            push!(ref_indices, i)
        elseif occursin("S(", config.name) || occursin("E(", config.name)
            push!(single_indices, i)
        elseif occursin("D(", config.name)
            push!(double_indices, i)
        elseif match(r"N\d+", config.name) !== nothing
            push!(nuclear_indices, i)
        elseif occursin("C(", config.name)
            push!(coupled_indices, i)
        end
    end
    
    # Analyze blocks
    if !isempty(ref_indices) && !isempty(single_indices)
        H_rs = H[ref_indices, single_indices]
        println("    Ref-Single block: max = $(round(maximum(abs.(H_rs)), digits=6))")
    end
    
    if !isempty(single_indices) && length(single_indices) > 1
        H_ss = H[single_indices, single_indices]
        println("    Single-Single block: max = $(round(maximum(abs.(H_ss)), digits=6))")
    end
    
    if !isempty(nuclear_indices) && length(nuclear_indices) > 1
        H_nn = H[nuclear_indices, nuclear_indices]
        println("    Nuclear-Nuclear block: max = $(round(maximum(abs.(H_nn)), digits=6))")
    end
    
    # Eigenvalue analysis
    println("\n  Eigenvalue analysis:")
    eigenvals = eigvals(Hermitian(H))
    println("    Ground state: $(round(eigenvals[1], digits=8)) Ha")
    println("    First 5 eigenvalues: $([round(e, digits=6) for e in eigenvals[1:min(5,end)]])")
    
    # Gap analysis
    if length(eigenvals) > 1
        gaps = diff(eigenvals)
        println("    HOMO-LUMO gap: $(round(gaps[1], digits=6)) Ha")
        println("    Max gap: $(round(maximum(gaps), digits=6)) Ha")
    end
    
    return eigenvals
end

eigenvals = analyze_hamiltonian_detailed(results_neo.hamiltonian_matrix, results_neo.configs)

# ======================== Example 4: Export to Quantum Computing Format ========================
println("\n\nExample 4: Exporting Hamiltonian for quantum computing")
println("-"^50)

# Export to OpenFermion format
filename_of = "water_neo_openfermion.txt"
export_hamiltonian_openfermion(ham_neo, filename_of)

println("Exported to OpenFermion format: $filename_of")

# Show sample of exported file
if isfile(filename_of)
    println("\nFirst 10 lines of OpenFermion file:")
    lines = readlines(filename_of)
    for i in 1:min(10, length(lines))
        println("  $i: $(lines[i])")
    end
    if length(lines) > 10
        println("  ... ($(length(lines) - 10) more lines)")
    end
end

# ======================== Example 5: Hamiltonian Storage ========================
println("\n\nExample 5: Hamiltonian storage analysis")
println("-"^50)

# Estimate storage requirements
storage = SparseQEEcNEO.HamiltonianConstruction.estimate_storage_requirements(ham_neo)

println("Storage requirements:")
println("  One-body integrals: $(round(storage.one_body_mb, digits=2)) MB")
println("  Two-body integrals: $(round(storage.two_body_mb, digits=2)) MB")
println("  Coupling terms: $(round(storage.coupling_mb, digits=2)) MB")
println("  Total: $(round(storage.total_mb, digits=2)) MB")

# ======================== Example 6: Effective Hamiltonian ========================
println("\n\nExample 6: Constructing effective Hamiltonians")
println("-"^50)

# Select subset of configurations based on importance
important_indices = Int[]
cumulative_weight = 0.0
total_weight = sum(c.weight for c in results_neo.configs)

for (i, config) in enumerate(results_neo.configs)
    push!(important_indices, i)
    cumulative_weight += config.weight
    
    if cumulative_weight / total_weight > 0.95
        break
    end
end

# Extract submatrix
H_full = results_neo.hamiltonian_matrix
H_eff = H_full[important_indices, important_indices]

println("Effective Hamiltonian:")
println("  Original size: $(size(H_full))")
println("  Reduced size: $(size(H_eff))")
println("  Configurations retained: $(length(important_indices))")
println("  Importance captured: $(round(cumulative_weight/total_weight * 100, digits=1))%")

# Compare eigenvalues
eigenvals_full = eigvals(Hermitian(H_full))
eigenvals_eff = eigvals(Hermitian(H_eff))

println("\nGround state comparison:")
println("  Full: $(round(eigenvals_full[1], digits=8)) Ha")
println("  Effective: $(round(eigenvals_eff[1], digits=8)) Ha")
println("  Error: $(round(abs(eigenvals_full[1] - eigenvals_eff[1]) * 1000, digits=3)) mHa")

# ======================== Summary ========================
println("\n" * "="^70)
println("Hamiltonian workflow examples completed!")
println("="^70)

println("\nKey takeaways:")
println("1. Hamiltonians are automatically constructed during sparse_qee_cneo")
println("2. NEO Hamiltonians include electron-nuclear coupling terms")
println("3. Hamiltonians can be analyzed for sparsity and block structure")
println("4. Export to quantum computing formats (OpenFermion)")
println("5. Effective Hamiltonians can be constructed from important configs")
println("6. HDF5 format enables efficient storage and retrieval")
