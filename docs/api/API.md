# SparseQEEcNEO.jl API Documentation

## Overview

SparseQEEcNEO.jl provides a comprehensive suite of functions for quantum chemistry calculations using Nuclear-Electronic Orbital (NEO) methods with quantum computing integration.

## Table of Contents

- [Core Functions](#core-functions)
- [Types and Structures](#types-and-structures)
- [cNEO Methods](#cneo-methods)
- [Configuration Generation](#configuration-generation)
- [Hamiltonian Construction](#hamiltonian-construction)
- [Quantum Computing Integration](#quantum-computing-integration)
- [Utility Functions](#utility-functions)

## Core Functions

### `sparse_qee_cneo(mol::Molecule; kwargs...)`

Main function for Sparse QEE-cNEO calculations with enhanced features.

**Arguments:**
- `mol::Molecule`: Molecular system with quantum nuclei specification
- `calc::NEOCalculation`: Calculation parameters (default: NEOCalculation())
- `config_sel::ConfigSelection`: Configuration selection parameters
- `neo_config::NEOConfig`: PySCF configuration

**Returns:**
- `NEOResults`: Complete results including configurations, energies, metrics, and Hamiltonian

**Example:**
```julia
mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
results = sparse_qee_cneo(mol)
println("Energy: $(results.total_energy) Ha")
```

### `sparse_qee_with_cneo(mol::Molecule, cneo_calc::CNEOCalculation; kwargs...)`

**NEW**: Integrated workflow combining cNEO constraint optimization with QEE Hamiltonian construction.

**Arguments:**
- `mol::Molecule`: Molecular system with quantum nuclei
- `cneo_calc::CNEOCalculation`: cNEO constraint parameters
- `config_sel::ConfigSelection`: Configuration selection for QEE
- `neo_config::NEOConfig`: PySCF configuration

**Returns:**
- `Tuple{CNEOResults, NEOResults}`: (cNEO optimization results, QEE results with Hamiltonian)

**Example:**
```julia
mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
cneo_calc = create_cneo_calculation(
    constraint_positions=[[0.0, 0.0, 0.5]]
)
cneo_results, qee_results = sparse_qee_with_cneo(mol, cneo_calc)
```

## Types and Structures

### `Molecule`

Represents a molecular system for quantum chemistry calculations.

**Fields:**
- `geometry::String`: Molecular geometry in XYZ format
- `basis::String`: Basis set specification
- `quantum_nuc::Vector{Int}`: Indices of quantum nuclei (default: [])
- `charge::Int`: Molecular charge (default: 0)
- `spin::Int`: Spin multiplicity (default: 1)

### `NEOCalculation`

Parameters for NEO quantum chemistry calculations.

**Fields:**
- `method::String`: Calculation method ("HF", "B3LYP", etc.)
- `epc::String`: Electron-proton correlation functional ("none", "17-1", "17-2", "18-1", "18-2")

### `ConfigSelection`

Configuration selection parameters for QEE methods.

**Fields:**
- `method::String`: Selection method ("mp2", "casci", "neo_cneo", "hybrid_final")
- `max_configs::Int`: Maximum number of configurations
- `max_qubits::Int`: Maximum number of qubits for quantum computing
- `importance_cutoff::Float64`: Importance threshold for configuration selection

### `NEOResults`

Results from SparseQEEcNEO calculations.

**Fields:**
- `total_energy::Float64`: Total energy in Hartree
- `n_configs::Int`: Number of selected configurations
- `n_qubits::Int`: Number of qubits required
- `captured_importance::Float64`: Fraction of importance captured
- `hamiltonian_matrix::Matrix`: Second-quantized Hamiltonian matrix
- `hamiltonian_data::HamiltonianData`: Structured Hamiltonian data
- `performance_metrics::Dict`: Computation time and memory usage

## cNEO Methods

### `CNEOCalculation`

Parameters for constrained NEO calculations.

**Fields:**
- `method::String`: "HF" or "DFT"
- `functional::String`: DFT functional (e.g., "B3LYP")
- `constraint_positions::Vector{Vector{Float64}}`: Target positions for quantum nuclei
- `convergence_threshold::Float64`: Convergence criterion
- `max_iterations::Int`: Maximum SCF iterations
- `lambda_damping::Float64`: Damping factor for Newton updates

### `create_cneo_calculation(; kwargs...)`

Factory function to create CNEOCalculation with sensible defaults.

**Keyword Arguments:**
- `method::String = "HF"`: Calculation method
- `functional::String = "B3LYP"`: DFT functional
- `constraint_positions::Vector{Vector{Float64}} = []`: Nuclear constraint positions
- `convergence_threshold::Float64 = 1e-6`: Convergence threshold
- `max_iterations::Int = 50`: Maximum iterations
- `lambda_damping::Float64 = 0.5`: Damping factor

**Example:**
```julia
cneo_calc = create_cneo_calculation(
    method="HF",
    constraint_positions=[[0.0, 0.0, 0.5]],
    convergence_threshold=1e-6
)
```

### `run_cneo_hf(molecule::Molecule, cneo_calculation::CNEOCalculation; neo_config)`

Run constrained NEO-HF calculation with position constraints on quantum nuclei.

**Returns:**
- `CNEOResults`: Complete results with energies, positions, and convergence information

### `run_cneo_mp2(molecule::Molecule, constraint_positions; neo_config)`

Run constrained NEO-MP2 calculation with correlation corrections.

**Returns:**
- `CNEOMP2Results`: Results including HF and MP2 energies with constraint analysis

## Configuration Generation

### `generate_configurations(meanfield_result, neo_molecule, config_selection, t2_amplitudes)`

Generate quantum configurations using the specified method.

**Methods Available:**
- `"mp2"`: MP2-based selection for weakly correlated systems
- `"casci"`: Complete active space CI
- `"neo_cneo"`: NEO with coupling for quantum nuclei systems
- `"hybrid_final"`: Combined approach for maximum accuracy

## Hamiltonian Construction

### `construct_second_quantized_hamiltonian(meanfield_result, configurations, config_selection)`

Construct second-quantized Hamiltonian matrix from selected configurations.

**Returns:**
- `HamiltonianResult`: Contains matrix data and analysis

### `save_hamiltonian(filename, hamiltonian_data, hamiltonian_matrix)`

Save Hamiltonian to HDF5 format for quantum computing applications.

### `export_hamiltonian_openfermion(hamiltonian_data, filename)`

Export Hamiltonian to OpenFermion format for quantum algorithms.

## Quantum Computing Integration

### Export Capabilities

The package provides seamless integration with major quantum computing frameworks:

- **OpenFermion**: Electronic structure operators
- **Qiskit**: Pauli operators for VQE algorithms
- **Cirq**: Circuit-ready quantum operators
- **HDF5**: Structured data storage

### Quantum Algorithms Ready

- **VQE (Variational Quantum Eigensolver)**: Ground state energy calculation
- **QAOA (Quantum Approximate Optimization Algorithm)**: Optimization problems
- **QPE (Quantum Phase Estimation)**: Spectroscopy applications

## Utility Functions

### `run_test_suite(config::NEOConfig; run_modular_tests=true)`

Run comprehensive test suite for the Sparse QEE-cNEO implementation.

### `analyze_hamiltonian_properties(hamiltonian_matrix::Matrix)`

Analyze Hamiltonian properties including eigenvalues, sparsity, and condition number.

### `check_pyscf_availability(config::NEOConfig)`

Verify PySCF with NEO support is properly installed and configured.

## Constants

### Physical Constants
```julia
const NUCLEAR_MASS_FACTOR = 1836.0        # Proton-to-electron mass ratio
const POSITION_DIMENSION = 3              # Spatial dimensions
```

### Computational Constants
```julia
const DEFAULT_CONVERGENCE_THRESHOLD = 1e-6
const DEFAULT_MAX_ITERATIONS = 50
const HAMILTONIAN_SPARSITY_THRESHOLD = 1e-12
```

### EPC Functionals
```julia
const EPC_PARAMS = Dict(
    "17-1" => Dict("a" => 0.17, "b" => 1.0, "c" => 0.5),
    "17-2" => Dict("a" => 0.17, "b" => 2.0, "c" => 0.5),
    "18-1" => Dict("a" => 0.18, "b" => 1.0, "c" => 0.5),
    "18-2" => Dict("a" => 0.18, "b" => 2.0, "c" => 0.5)
)
```

## Error Handling

The package provides comprehensive error handling with descriptive messages:

- `ArgumentError`: Invalid parameters or configurations
- `MethodError`: Incorrect function usage
- `BoundsError`: Array index out of bounds
- Custom exceptions for quantum chemistry specific errors

## Performance Considerations

- Set `JULIA_NUM_THREADS` environment variable for parallel execution
- Use `max_configs` to limit memory usage for large systems
- Enable `use_compression` for large configuration spaces
- Consider `max_qubits` constraint for quantum computing applications

## Integration Examples

See the `examples/` directory for comprehensive usage examples:
- `basic_usage.jl`: Simple NEO calculations
- `cneo_example.jl`: Constrained NEO demonstrations
- `cneo_qee_integration.jl`: Complete cNEO-QEE workflow
- `quantum_computing_demo.jl`: Quantum computing integration

For more detailed examples and tutorials, see the `docs/tutorials/` directory.