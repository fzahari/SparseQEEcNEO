# Constants for SparseQEEcNEO.jl
# Organized by category for Clean Code compliance

# Physical Constants
const NUCLEAR_MASS_FACTOR = 1836.0        # Proton-to-electron mass ratio
const POSITION_DIMENSION = 3              # Spatial dimensions
const NUCLEAR_CHARGE_UNIT = 1.0          # Elementary charge
const MASS_ELECTRON = 1.0                # Electron mass in atomic units

# Chemical Constants
const HYDROGEN_MASS = 1.008              # Hydrogen atomic mass
const CARBON_MASS = 12.0                 # Carbon-12 atomic mass
const WATER_OH_DISTANCE = 0.957          # O-H bond distance in Angstrom
const WATER_HOH_ANGLE = 104.5           # H-O-H angle in degrees

# Computational Thresholds
const DEFAULT_CONVERGENCE_THRESHOLD = 1e-6
const DEFAULT_MAX_ITERATIONS = 50
const HAMILTONIAN_SPARSITY_THRESHOLD = 1e-12
const ENERGY_CONVERGENCE_THRESHOLD = 1e-8
const ORBITAL_OCCUPATION_THRESHOLD = 0.01
const NUMERICAL_PRECISION_CUTOFF = 1e-10

# Configuration Selection Parameters
const DEFAULT_MAX_CONFIGS = 1000
const DEFAULT_IMPORTANCE_CUTOFF = 0.01
const MAX_SINGLE_EXCITATIONS = 100
const MAX_DOUBLE_EXCITATIONS = 500
const CONFIGURATION_COMPRESSION_THRESHOLD = 0.99

# Quantum Computing Parameters
const DEFAULT_MAX_QUBITS = 30
const QUBIT_OVERHEAD_FACTOR = 2
const VQE_PARAMETER_SCALING = 4
const CIRCUIT_DEPTH_MULTIPLIER = 8
const MEASUREMENT_OVERHEAD = 1000

# EPC Functional Parameters
const EPC_17_1_PARAMS = Dict("a" => 0.17, "b" => 1.0, "c" => 0.5)
const EPC_17_2_PARAMS = Dict("a" => 0.17, "b" => 2.0, "c" => 0.5)
const EPC_18_1_PARAMS = Dict("a" => 0.18, "b" => 1.0, "c" => 0.5)
const EPC_18_2_PARAMS = Dict("a" => 0.18, "b" => 2.0, "c" => 0.5)

# cNEO Method Parameters
const DEFAULT_CNEO_CONVERGENCE = 1e-6
const DEFAULT_CNEO_MAX_ITERATIONS = 50
const DEFAULT_LAMBDA_DAMPING = 0.5
const NUCLEAR_POSITION_TOLERANCE = 0.1
const CONSTRAINT_SATISFACTION_THRESHOLD = 1e-5

# Memory and Performance Parameters
const MEMORY_ALLOCATION_BUFFER = 1024    # MB
const THREAD_POOL_SIZE = 4
const CHUNK_SIZE_CONFIGURATIONS = 100
const STORAGE_COMPRESSION_LEVEL = 6

# File I/O Parameters
const HDF5_COMPRESSION_LEVEL = 6
const JSON_INDENT_LEVEL = 2
const OUTPUT_PRECISION_DIGITS = 8

# Basis Set Parameters
const DEFAULT_BASIS_SET = "sto-3g"
const MINIMAL_BASIS_THRESHOLD = 5
const EXTENDED_BASIS_THRESHOLD = 20

# NEO-specific Parameters
const NUCLEAR_ORBITAL_CUTOFF = 0.5
const ELECTRON_NUCLEAR_COUPLING_THRESHOLD = 0.1
const NEO_ENHANCEMENT_FACTOR = 2.0
const NUCLEAR_KINETIC_SCALING = 0.8

# Test Suite Parameters
const TEST_ENERGY_TOLERANCE = 1e-4
const TEST_H2_BOND_LENGTH = 0.74        # Angstrom
const TEST_H2O_OH_DISTANCE = 0.957      # Angstrom
const BENCHMARK_TIMEOUT_SECONDS = 60

# Display and Formatting
const RESULTS_DISPLAY_WIDTH = 60
const ENERGY_FORMAT_PRECISION = 6
const PERCENTAGE_FORMAT_PRECISION = 2
const PROGRESS_BAR_WIDTH = 40

# Error Messages
const ERROR_PYSCF_NOT_FOUND = "PySCF with NEO support not found. Please check PYSCF_PATH environment variable."
const ERROR_INVALID_MOLECULE = "Invalid molecular specification. Please check geometry and basis set."
const ERROR_CONVERGENCE_FAILED = "SCF convergence failed. Consider adjusting convergence parameters."
const ERROR_MEMORY_EXCEEDED = "Memory limit exceeded. Consider reducing max_configs or enabling compression."

# Success Messages
const SUCCESS_CALCULATION_COMPLETE = "Calculation completed successfully"
const SUCCESS_HAMILTONIAN_CONSTRUCTED = "Hamiltonian matrix constructed successfully"
const SUCCESS_CONSTRAINTS_SATISFIED = "Nuclear constraints satisfied within tolerance"

# Method Names (for validation)
const VALID_NEO_METHODS = ["HF", "B3LYP", "PBE", "M06-2X"]
const VALID_CONFIG_METHODS = ["mp2", "casci", "neo_cneo", "hybrid_final"]
const VALID_EPC_FUNCTIONALS = ["none", "17-1", "17-2", "18-1", "18-2"]
const VALID_DFT_FUNCTIONALS = ["B3LYP", "PBE", "M06-2X", "wB97X-D"]

# Default File Extensions
const HDF5_EXTENSION = ".h5"
const JSON_EXTENSION = ".json"
const LOG_EXTENSION = ".log"
const CHECKPOINT_EXTENSION = ".chk"