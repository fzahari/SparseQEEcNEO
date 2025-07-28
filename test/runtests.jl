# test/runtests.jl - Complete test suite for SparseQEEcNEO

using Test
using LinearAlgebra
using Random
using PyCall

# Set up environment
ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/Users/federicozahariev/Work/Programs/QEE_Split_Grouping/pyscf-master")

# Include the module
include("../src/SparseQEEcNEO.jl")
using .SparseQEEcNEO

# Set random seed for reproducibility
Random.seed!(42)

# Global test configuration
const TEST_CONFIG = NEOConfig(
    pyscf_path=get(ENV, "PYSCF_PATH", "/Users/federicozahariev/Work/Programs/QEE_Split_Grouping/pyscf-master")
)

# Check if PySCF is available
const PYSCF_AVAILABLE, NEO_AVAILABLE = try
    setup_pyscf(TEST_CONFIG)
catch
    (false, false)
end

@testset "SparseQEEcNEO.jl Complete Test Suite" begin
    include("test_types.jl")
    include("test_epc_functionals.jl")
    include("test_configuration_generation.jl")
    include("test_importance_analysis.jl")
    include("test_hamiltonian_construction.jl")
    include("test_qee_methods.jl")
    include("test_nuclear_methods.jl")
    
    if PYSCF_AVAILABLE
        include("test_pyscf_interface.jl")
        include("test_integration.jl")
    else
        @warn "Skipping PySCF-dependent tests - PySCF not available"
    end
    
    include("test_performance.jl")
end
