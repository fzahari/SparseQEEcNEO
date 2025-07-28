================================================================================
                         SparseQEEcNEO.jl Documentation
                                  Version 6.0
================================================================================

OVERVIEW
--------
Sparse Quantum Eigensolver with constrained Nuclear-Electronic Orbital (cNEO) 
support for quantum computing applications.

FEATURES
--------
- Nuclear-Electronic Orbital (NEO) calculations with quantum nuclei
- Second-quantized Hamiltonian construction for quantum computing
- Multiple configuration selection methods (MP2, CASCI, NEO-enhanced)
- Electron-proton correlation (EPC) functionals
- Export to OpenFermion format

INSTALLATION
------------
1. Clone repository:
   git clone https://github.com/yourusername/SparseQEEcNEO
   cd SparseQEEcNEO

2. Install Julia dependencies:
   julia --project=.
   julia> using Pkg; Pkg.instantiate()

3. Install PySCF with NEO (optional but recommended):
   git clone https://github.com/corinwagen/pyscf
   cd pyscf && git checkout neo
   python setup.py install
   export PYSCF_PATH="/path/to/pyscf"

QUICK START
-----------
using SparseQEEcNEO

# H2 with quantum proton
mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
results = sparse_qee_cneo(mol)

println("Energy: $(results.total_energy) Ha")
println("Configurations: $(results.n_configs)")

# Save Hamiltonian
save_hamiltonian("h2.h5", results.hamiltonian_data, results.hamiltonian_matrix)

# Export for quantum computing
export_hamiltonian_openfermion(results.hamiltonian_data, "h2_quantum.txt")

DIRECTORY STRUCTURE
-------------------
SparseQEEcNEO/
├── Project.toml                    # Package dependencies
├── README.md                       # Main documentation
├── LICENSE
├── src/                           # Source code
│   ├── SparseQEEcNEO.jl          # Main module
│   ├── types.jl                  # Type definitions
│   ├── pyscf_interface.jl        # PySCF interface
│   ├── epc_functionals.jl        # EPC functionals
│   ├── nuclear_methods.jl        # Nuclear methods
│   ├── configuration_generation.jl
│   ├── importance_analysis.jl
│   ├── qee_methods.jl
│   └── hamiltonian_construction.jl
├── test/                          # Test suite
│   ├── runtests.jl
│   └── test_*.jl
├── examples/                      # Usage examples
│   ├── basic_usage.jl
│   ├── hamiltonian_workflow.jl
│   ├── batch_calculations.jl
│   └── visualization_demo.jl
└── tools/                         # Utilities
    ├── check_dependencies.jl
    ├── analyze_configurations.jl
    ├── visualize_hamiltonian.jl
    └── benchmark_suite.jl

USAGE EXAMPLES
--------------

Water with Quantum Proton:
--------------------------
# Water with second H as quantum
mol = Molecule("O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", 
               "6-31g", quantum_nuc=[1])

# B3LYP with EPC functional
calc = NEOCalculation(xc="B3LYP", epc="17-2")

# NEO-specific configuration selection
config_sel = ConfigSelection(method="neo_cneo", max_configs=200)

results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel)

Method Comparison:
------------------
methods = ["mp2", "neo_cneo", "hybrid_final"]

for method in methods
    config_sel = ConfigSelection(method=method, max_configs=100)
    res = sparse_qee_cneo(mol, config_sel=config_sel)
    println("$method: E = $(res.total_energy) Ha, $(res.n_configs) configs")
end

CONFIGURATION METHODS
---------------------
Method          Description                 Use Case
------          -----------                 --------
mp2             MP2-based selection         Weakly correlated systems
neo_cneo        NEO with coupling           Systems with quantum nuclei
hybrid_final    Combined approach           Maximum accuracy

KEY TYPES
---------

Molecule:
    Molecule(geometry, basis; charge=0, spin=0, quantum_nuc=[], nuc_basis="pb4d")

NEOCalculation:
    NEOCalculation(; xc="HF", epc="none", conv_tol=1e-8, max_cycle=100)

ConfigSelection:
    ConfigSelection(; method="mp2", max_configs=1000, max_nuc_orbs=0, 
                    use_neo_importance=true, importance_cutoff=1e-6)

MAIN API
--------
sparse_qee_cneo(mol; calc, config_sel, neo_config) → NEOResults
    Main entry point for calculations

save_hamiltonian(filename, ham_data, H_matrix)
    Save Hamiltonian to HDF5 file

load_hamiltonian(filename) → (ham_data, H_matrix)
    Load Hamiltonian from HDF5 file

export_hamiltonian_openfermion(ham_data, filename)
    Export to OpenFermion format

RUNNING TESTS
-------------
# All tests
include("test/runtests.jl")

# Specific tests
include("test/test_hamiltonian_construction.jl")

TOOLS
-----

Check Dependencies:
    julia tools/check_dependencies.jl

Analyze Configurations:
    julia tools/analyze_configurations.jl hamiltonian.h5

Benchmark Performance:
    julia tools/benchmark_suite.jl

Visualize Hamiltonian:
    julia tools/visualize_hamiltonian.jl hamiltonian.h5

PERFORMANCE TIPS
----------------
1. Memory: Use max_configs to limit configurations
2. Speed: Set JULIA_NUM_THREADS for parallelization
3. Nuclear Space: Use max_nuc_orbs to truncate nuclear orbitals
4. Importance: Set importance_cutoff to ignore small contributions

TROUBLESHOOTING
---------------
PySCF not found: 
    Set PYSCF_PATH environment variable

Memory errors: 
    Reduce max_configs or enable use_compression=true

Convergence issues: 
    Increase max_cycle in NEOCalculation

Debug mode: 
    Set ENV["JULIA_DEBUG"] = "SparseQEEcNEO"

THEORY
------
The NEO method treats selected nuclei (usually protons) quantum mechanically 
alongside electrons, capturing:
- Zero-point energy
- Quantum tunneling
- Isotope effects
- Proton-coupled electron transfer

EPC functionals (17-1/2, 18-1/2) provide electron-proton correlation corrections.

REFERENCES
----------
1. Hammes-Schiffer et al., Chem. Rev. 2020, 120, 4222 (NEO Theory)
2. Brorsen et al., J. Chem. Phys. 2017, 147, 114113 (EPC Functionals)

LICENSE
-------
See LICENSE file for details.

================================================================================
                              End of Documentation
================================================================================
