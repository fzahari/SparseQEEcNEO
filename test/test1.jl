# test1.jl - Test H2 with quantum proton

# Setup path to PySCF
ENV["PYTHONPATH"] = "/Users/federicozahariev/Work/Programs/QEE_Split_Grouping/pyscf-master"

# Load the package
push!(LOAD_PATH, dirname(@__DIR__))  # Add parent directory to load path
using SparseQEEcNEO

# Now you can use the types and functions
mol_h2 = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
config = NEOConfig(pyscf_path="/Users/federicozahariev/Work/Programs/QEE_Split_Grouping/pyscf-master")
calc = NEOCalculation(xc="HF", epc="none")
config_sel = ConfigSelection(
    method="mp2", 
    max_configs=50,
    max_nuc_orbs=50  # High enough to avoid truncation
)

println("Running H2 calculation with quantum proton...")
results_h2 = sparse_qee_cneo(mol_h2, calc=calc, config_sel=config_sel, neo_config=config)

println("\nH2 Results Summary:")
println("Energy: $(results_h2.total_energy) Ha")
println("Configurations: $(results_h2.n_configs)")
println("MP2 correlation: $(results_h2.mp2_correlation) Ha")
