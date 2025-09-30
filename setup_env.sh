#!/bin/bash
# SparseQEEcNEO environment setup

# Set PySCF path (update this to your PySCF installation)
export PYSCF_PATH="/path/to/pyscf-master"

# Set Python path
export PYTHONPATH="$PYSCF_PATH:$PYTHONPATH"

# Set number of threads for parallel execution
export JULIA_NUM_THREADS=4
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

# Fix OpenMP library conflict on macOS
export KMP_DUPLICATE_LIB_OK=TRUE

echo "SparseQEEcNEO environment configured"
echo "PYSCF_PATH: $PYSCF_PATH"
echo "Julia threads: $JULIA_NUM_THREADS"
