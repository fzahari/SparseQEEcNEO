#!/usr/bin/env julia
"""
tools/check_dependencies.jl - Check and report on required dependencies
"""

using Pkg
using PyCall

function check_julia_packages()
    required_packages = [
        "PyCall", "LinearAlgebra", "SparseArrays", 
        "Printf", "Statistics", "HDF5", "JSON", "Test"
    ]
    
    optional_packages = [
        "Plots", "DataFrames", "CSV", "BenchmarkTools"
    ]
    
    println("Checking Julia packages...")
    println("\nRequired packages:")
    
    all_required_installed = true
    for pkg in required_packages
        installed = haskey(Pkg.project().dependencies, pkg) || 
                   haskey(Pkg.dependencies(), pkg)
        
        if installed
            println("  ✓ $pkg")
        else
            println("  ✗ $pkg - Install with: Pkg.add(\"$pkg\")")
            all_required_installed = false
        end
    end
    
    println("\nOptional packages (for examples/tools):")
    for pkg in optional_packages
        installed = haskey(Pkg.project().dependencies, pkg) || 
                   haskey(Pkg.dependencies(), pkg)
        
        if installed
            println("  ✓ $pkg")
        else
            println("  ○ $pkg - Install with: Pkg.add(\"$pkg\")")
        end
    end
    
    return all_required_installed
end

function check_python_packages()
    println("\nChecking Python packages...")
    
    py"""
    import subprocess
    import sys
    import importlib
    
    def check_package(package):
        try:
            importlib.import_module(package)
            return True, None
        except ImportError as e:
            return False, str(e)
    
    required = ['numpy', 'scipy', 'h5py']
    optional = ['matplotlib']
    
    print("\\nRequired Python packages:")
    all_installed = True
    for pkg in required:
        success, error = check_package(pkg)
        if success:
            print(f"  ✓ {pkg}")
        else:
            print(f"  ✗ {pkg} - Install with: {sys.executable} -m pip install {pkg}")
            all_installed = False
    
    print("\\nOptional Python packages:")
    for pkg in optional:
        success, error = check_package(pkg)
        if success:
            print(f"  ✓ {pkg}")
        else:
            print(f"  ○ {pkg} - Install with: {sys.executable} -m pip install {pkg}")
    """
    
    return py"all_installed"
end

function check_pyscf_neo()
    println("\nChecking PySCF with NEO...")
    
    pyscf_path = get(ENV, "PYSCF_PATH", "")
    if isempty(pyscf_path)
        println("  ⚠ PYSCF_PATH environment variable not set")
        println("  Set with: export PYSCF_PATH=/path/to/pyscf-master")
        return false
    end
    
    println("  PYSCF_PATH: $pyscf_path")
    
    if !isdir(pyscf_path)
        println("  ✗ Directory does not exist: $pyscf_path")
        return false
    end
    
    try
        py"""
        import sys
        sys.path.insert(0, $pyscf_path)
        
        try:
            import pyscf
            print(f"  ✓ PySCF found: {pyscf.__file__}")
            
            # Check for NEO
            has_neo = hasattr(pyscf, 'neo')
            if has_neo:
                print("  ✓ NEO module available")
                
                # Check NEO submodules
                neo_modules = ['hf', 'mp2', 'ks', 'mole']
                missing = []
                for mod in neo_modules:
                    if not hasattr(pyscf.neo, mod):
                        missing.append(mod)
                
                if missing:
                    print(f"  ⚠ Missing NEO submodules: {', '.join(missing)}")
                else:
                    print("  ✓ All NEO submodules present")
                    
                neo_available = True
            else:
                print("  ✗ NEO module not found")
                print("  Install NEO-enabled PySCF from: https://github.com/corinwagen/pyscf")
                neo_available = False
                
        except ImportError as e:
            print(f"  ✗ Failed to import PySCF: {e}")
            neo_available = False
        """
        
        return py"neo_available"
    catch e
        println("  ✗ Error checking PySCF: $e")
        return false
    end
end

function check_system_info()
    println("\nSystem Information:")
    println("  Julia version: $(VERSION)")
    println("  OS: $(Sys.KERNEL) ($(Sys.MACHINE))")
    println("  CPU threads: $(Sys.CPU_THREADS)")
    println("  WORD_SIZE: $(Sys.WORD_SIZE)")
    
    py"""
    import sys
    print(f"  Python version: {sys.version.split()[0]}")
    print(f"  Python executable: {sys.executable}")
    """
end

function generate_env_script()
    println("\nGenerating environment setup script...")
    
    script_content = """
    #!/bin/bash
    # SparseQEEcNEO environment setup
    
    # Set PySCF path (update this to your PySCF installation)
    export PYSCF_PATH="/path/to/pyscf-master"
    
    # Set Python path
    export PYTHONPATH="\$PYSCF_PATH:\$PYTHONPATH"
    
    # Set number of threads for parallel execution
    export JULIA_NUM_THREADS=4
    export OMP_NUM_THREADS=4
    export MKL_NUM_THREADS=4
    
    echo "SparseQEEcNEO environment configured"
    echo "PYSCF_PATH: \$PYSCF_PATH"
    echo "Julia threads: \$JULIA_NUM_THREADS"
    """
    
    filename = "setup_env.sh"
    open(filename, "w") do io
        write(io, script_content)
    end
    
    chmod(filename, 0o755)
    println("  Created $filename")
    println("  Edit this file to set your PySCF path, then run: source $filename")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    println("Sparse QEE-cNEO Dependency Check")
    println("="^60)
    
    check_system_info()
    
    julia_ok = check_julia_packages()
    python_ok = check_python_packages()
    neo_ok = check_pyscf_neo()
    
    println("\n" * "="^60)
    println("Summary:")
    
    if julia_ok && python_ok && neo_ok
        println("✓ All required dependencies are installed!")
        println("\nYou can run the examples with:")
        println("  julia examples/basic_usage.jl")
    else
        println("⚠ Some dependencies are missing or not configured")
        println("\nPlease install missing dependencies and ensure PYSCF_PATH is set correctly")
        
        if !neo_ok
            generate_env_script()
        end
    end
    
    println("="^60)
end
