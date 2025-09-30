# SparseQEEcNEO.jl

**Sparse Quantum Eigensolver with constrained Nuclear-Electronic Orbital (cNEO) support for quantum computing applications.**

[![Julia](https://img.shields.io/badge/julia-v1.6+-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Features

- **Nuclear-Electronic Orbital (NEO) calculations** with quantum nuclei treatment
- **Second-quantized Hamiltonian construction** for quantum computing
- **Multiple configuration selection methods** (MP2, CASCI, NEO-enhanced)
- **Electron-proton correlation (EPC) functionals** (17-1/2, 18-1/2)
- **Export to OpenFermion format** for quantum algorithms
- **Constrained NEO (cNEO) methods** with Lagrange multiplier constraints

## Installation Requirements

⚠️ **IMPORTANT**: This package requires PySCF with NEO (Nuclear-Electronic Orbital) support from the specific fork at https://github.com/theorychemyang/pyscf. Standard PySCF installations will NOT work.

## Quick Setup

1. **Install PySCF with NEO support (required):**
   ```bash
   # Clone the specific PySCF fork with NEO support
   git clone https://github.com/theorychemyang/pyscf.git
   cd pyscf
   pip install -e .
   
   # Set environment variables
   export PYSCF_PATH="$(pwd)"
   export PYTHONPATH="$PYSCF_PATH:$PYTHONPATH"
   ```

2. **Clone and setup SparseQEEcNEO:**
   ```bash
   git clone https://github.com/fzahari/SparseQEEcNEO.jl.git
   cd SparseQEEcNEO
   source setup_env.sh  # Sets up PySCF paths and environment
   ```

3. **Install Julia dependencies:**
   ```julia
   julia --project=.
   julia> using Pkg; Pkg.instantiate()
   ```

4. **Verify setup:**
   ```bash
   julia tools/check_dependencies.jl
   ```

## Quick Example

```julia
using SparseQEEcNEO

# H2 molecule with quantum proton
mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
results = sparse_qee_cneo(mol)

println("Energy: $(results.total_energy) Ha")
println("Configurations: $(results.n_configs)")

# Save Hamiltonian for quantum computing
save_hamiltonian("h2.h5", results.hamiltonian_data, results.hamiltonian_matrix)
export_hamiltonian_openfermion(results.hamiltonian_data, "h2_quantum.txt")
```

## Constrained NEO (cNEO) Methods

The package includes integrated constrained NEO methods with position constraints on quantum nuclei:

```julia
using SparseQEEcNEO

# cNEO-HF with nuclear position constraint
mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
constraint_position = [0.0, 0.0, 0.5]  # Target position for quantum proton

# Create cNEO calculation
cneo_calc = create_cneo_calculation(
    method="HF",
    constraint_positions=[constraint_position],
    convergence_threshold=1e-6,
    max_iterations=50
)

# Run cNEO-HF calculation
cneo_results = run_cneo_hf(mol, cneo_calc)

println("cNEO-HF Energy: $(cneo_results.total_energy) Ha")
println("Constraint satisfied: $(cneo_results.converged)")
println("Position error: $(cneo_results.position_errors[1])")

# Run cNEO-MP2 with correlation
mp2_results = run_cneo_mp2(mol, [constraint_position])
println("cNEO-MP2 Energy: $(mp2_results.total_energy) Ha")
```

### cNEO Features

- **Nuclear Position Constraints**: Enforce specific positions for quantum nuclei
- **Lagrange Multiplier Optimization**: Newton method with analytical Hessians
- **cNEO-HF and cNEO-MP2**: Hartree-Fock and correlation-corrected methods
- **Clean Code Implementation**: Small, focused functions following Clean Code principles
- **Comprehensive Error Handling**: Robust validation and convergence monitoring

## Development

**📖 For detailed development information, see [WARP.md](WARP.md)**

The WARP.md file contains comprehensive guidance including:
- Environment setup commands
- Testing procedures
- Code architecture overview
- Development tools and utilities
- Configuration methods and examples

### Quick Commands

```bash
# Run examples
julia examples/basic_usage.jl

# Run tests
julia test/runtests.jl

# Run benchmarks
julia tools/benchmark_suite.jl
```

## Theory

The Nuclear-Electronic Orbital (NEO) method treats selected nuclei (typically protons) quantum mechanically alongside electrons, enabling accurate modeling of:

- **Zero-point energy** and nuclear quantum effects
- **Quantum tunneling** in chemical reactions  
- **Isotope effects** and nuclear dynamics
- **Proton-coupled electron transfer** processes

The package implements electron-proton correlation (EPC) functionals and constrained NEO methods for enhanced accuracy in quantum chemistry calculations.

## References

1. Hammes-Schiffer, S. et al. *Chem. Rev.* **2020**, *120*, 4222-4261. "Nuclear-Electronic Orbital Methods"
2. Brorsen, K. R. et al. *J. Chem. Phys.* **2017**, *147*, 114113. "Electron-Proton Correlation Functionals"

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**📖 For comprehensive documentation, development setup, and detailed examples, see [WARP.md](WARP.md)**
