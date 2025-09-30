#!/bin/bash
# SparseQEEcNEO environment setup

# IMPORTANT: PySCF must be installed from https://github.com/theorychemyang/pyscf
# Standard PySCF installations will NOT work - NEO support is required

# Set PySCF path (update this to your PySCF installation with NEO support)
export PYSCF_PATH="/path/to/pyscf"  # Point to https://github.com/theorychemyang/pyscf clone

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
