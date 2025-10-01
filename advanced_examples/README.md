# Constrained Nuclear-Electronic Orbital (cNEO) Methods

This directory contains advanced implementations of constrained Nuclear-Electronic Orbital (cNEO) methods that enforce specific nuclear positions using Lagrange multiplier techniques.

## Overview

These modules implement sophisticated constrained NEO calculations:

- **`cneo_hf.jl`**: Constrained NEO-HF with Lagrange multiplier constraints
- **`cneo_mp2.jl`**: Constrained NEO-MP2 with correlation corrections
- **`test_cneo_hf.jl`**: Tests for cNEO-HF implementation
- **`test_cneo_mp2.jl`**: Tests for cNEO-MP2 implementation

## Features

### cNEO-HF (`cneo_hf.jl`)
- Newton's method for Lagrange multiplier optimization
- Analytical Hessian calculation using perturbation theory
- Nuclear position constraint enforcement
- Support for both HF and DFT functionals

### cNEO-MP2 (`cneo_mp2.jl`)
- Non-canonical MP2 implementation for constrained systems
- Iterative t-amplitude calculations
- MP2 density matrix construction
- Optional MP2-level constraint enforcement

## Theory

The cNEO method constrains quantum nuclear positions to specific values R₀ by adding Lagrange multiplier terms:

```
L = E[ψ] + Σᵢ λᵢ·(⟨ψ|r̂ᵢ|ψ⟩ - R₀ᵢ)
```

Where:
- `E[ψ]` is the electronic/nuclear energy functional  
- `λᵢ` are Lagrange multipliers for nucleus i
- `⟨ψ|r̂ᵢ|ψ⟩` is the expectation value of nuclear position
- `R₀ᵢ` is the desired constraint position

## Usage Example

```julia
# Load modules (adjust path as needed)
push!(LOAD_PATH, "path/to/SparseQEEcNEO")
include("cneo_hf.jl")
include("cneo_mp2.jl")

using .CNEOHF
using .CNEOMP2
using SparseQEEcNEO

# Define molecule with quantum nucleus
mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])

# Set constraint position for first hydrogen
constraint_pos = [[0.1, 0.0, 0.0]]  # Constrain to x=0.1

# Create calculation object
cneo_calc = CNEOCalculation(
    method="HF",
    constraint_positions=constraint_pos,
    convergence_threshold=1e-6
)

# Run cNEO-HF calculation
results = run_cneo_hf(mol, cneo_calc)

println("Constrained energy: $(results.energy) Ha")
println("Nuclear positions: $(results.nuclear_positions)")
println("Lagrange multipliers: $(results.lagrange_multipliers)")
```

## Status

⚠️ **Note**: These modules are not currently integrated into the main SparseQEEcNEO package due to circular dependency issues. They are provided as standalone advanced examples for researchers interested in constrained NEO methods.

## Future Integration

To integrate these modules into the main package, the following would need to be addressed:

1. **Resolve circular dependencies** between cNEO modules
2. **Unify API** with main sparse QEE workflow  
3. **Add to main exports** and documentation
4. **Integration tests** with main package functionality

## References

1. Hammes-Schiffer, S. "Nuclear-Electronic Orbital Methods" *Chem. Rev.* **2020**
2. Brorsen, K. R. "Multicomponent Quantum Chemistry" *J. Chem. Phys.* **2017**

## Testing

Run the tests with:

```bash
julia test_cneo_hf.jl
julia test_cneo_mp2.jl
```

Note: Tests require PySCF with NEO support to be properly configured.