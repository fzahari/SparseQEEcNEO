#!/usr/bin/env julia
"""
examples/basic_usage.jl - Simple examples of Sparse QEE-cNEO usage
"""

# Setup environment
ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/path/to/pyscf-master")

# Load the module
push!(LOAD_PATH, dirname(@__DIR__))
using SparseQEEcNEO
using LinearAlgebra

println("Sparse QEE-cNEO Basic Usage Examples")
println("="^60)

# Create NEO configuration
config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))

# ======================== Example 1: Simple H2 ========================
println("\nExample 1: H2 molecule (electronic only)")
println("-"^40)

# Create H2 molecule without quantum nuclei
mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")

# Run with default settings (HF)
results = sparse_qee_cneo(mol_h2, neo_config=config)

println("Results:")
println("  Energy: $(round(results.energy, digits=6)) Ha")
println("  MP2 correlation: $(round(results.mp2_correlation, digits=6)) Ha")
println("  Total energy: $(round(results.total_energy, digits=6)) Ha")
println("  Configurations: $(results.n_configs)")
println("  Qubits needed: $(results.n_qubits)")

# ======================== Example 2: H2 with Quantum Proton ========================
println("\nExample 2: H2 with quantum proton")
println("-"^40)

# Create H2 with first H as quantum
mol_h2_neo = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])

# Use B3LYP with EPC functional
calc = NEOCalculation(xc="B3LYP", epc="17-2")
config_sel = ConfigSelection(method="neo_cneo", max_configs=50, max_nuc_orbs=0)
results_neo = sparse_qee_cneo(mol_h2_neo, calc=calc, config_sel=config_sel, neo_config=config)

println("Results with quantum proton:")
println("  Energy: $(round(results_neo.energy, digits=6)) Ha")
println("  Total energy (with EPC): $(round(results_neo.total_energy, digits=6)) Ha")
println("  Configurations: $(results_neo.n_configs)")

# Compare energies
println("\nEnergy difference (NEO - classical): $(round(results_neo.energy - results.energy, digits=6)) Ha")

# ======================== Example 3: Different Methods ========================
println("\nExample 3: Comparing configuration methods")
println("-"^40)

methods = ["mp2", "neo_cneo", "hybrid_final"]
method_results = Dict()

for method in methods
    config_sel = ConfigSelection(method=method, max_configs=50, max_nuc_orbs=0)
    res = sparse_qee_cneo(mol_h2_neo, calc=calc, config_sel=config_sel, neo_config=config)
    method_results[method] = res
    
    println("$method:")
    println("  Configurations: $(res.n_configs)")
    println("  Importance: $(round(res.captured_importance * 100, digits=1))%")
    println("  Qubits: $(res.n_qubits)")
end

# ======================== Example 4: Accessing Hamiltonian ========================
println("\nExample 4: Working with Hamiltonians")
println("-"^40)

# Get Hamiltonian from results
H = results_neo.hamiltonian_matrix
ham_data = results_neo.hamiltonian_data

println("Hamiltonian information:")
println("  Size: $(size(H))")
println("  Hermitian: $(ishermitian(H))")
println("  Components: $(keys(ham_data.h1e))")

# Calculate eigenvalues
eigenvals = eigvals(Hermitian(H))
println("  Ground state: $(round(eigenvals[1], digits=6)) Ha")
if length(eigenvals) > 1
    println("  First excited: $(round(eigenvals[2], digits=6)) Ha")
    println("  Gap: $(round(eigenvals[2] - eigenvals[1], digits=6)) Ha")
end

# ======================== Example 5: Save Results ========================
println("\nExample 5: Saving results")
println("-"^40)

# Save Hamiltonian
filename = "h2_neo_hamiltonian.h5"
save_hamiltonian(filename, ham_data, H)
println("Hamiltonian saved to: $filename")

# Show first few configurations
println("\nFirst 5 configurations:")
for (i, config) in enumerate(results_neo.configs[1:min(5, end)])
    println("  $i. $(config.name) - weight: $(round(config.weight, digits=4))")
end

println("\n" * "="^60)
println("Basic usage examples completed!")
