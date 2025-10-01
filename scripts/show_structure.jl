#!/usr/bin/env julia

"""
Display clean directory structure of SparseQEEcNEO project
"""

function show_directory_structure()
    println("📁 SparseQEEcNEO.jl Project Structure")
    println("="^60)
    
    # Core directories and their contents
    directories = [
        ("📂 src/", "Core source code modules"),
        ("📂 test/", "Unit and integration tests"),
        ("📂 examples/", "Usage examples and demos"), 
        ("📂 advanced_examples/", "Advanced cNEO implementations"),
        ("📂 scripts/", "Development and analysis tools"),
        ("📂 docs/", "Documentation and API reference"),
        ("📂 benchmarks/", "Performance benchmarking"),
        ("📂 quantum_extension/", "Quantum computing integration"),
        ("📂 quantum_integration/", "Additional quantum features"),
        ("📂 notebooks/", "Jupyter notebooks"),
        ("📂 .github/", "GitHub workflows and templates")
    ]
    
    for (dir, desc) in directories
        if isdir(joinpath(pwd(), replace(dir, "📂 " => "", "/" => "")))
            println("$dir - $desc")
        end
    end
    
    println("\n📋 Key Files:")
    key_files = [
        ("README.md", "Project overview and installation"),
        ("Project.toml", "Julia package configuration"),
        ("LICENSE", "MIT License"),
        ("CITATION.bib", "Citation information"),
        ("CONTRIBUTING.md", "Development guidelines"),
        ("setup_env.sh", "Environment setup script")
    ]
    
    for (file, desc) in key_files
        if isfile(joinpath(pwd(), file))
            println("  📄 $file - $desc")
        end
    end
    
    println("\n🔍 Source Code Structure (src/):")
    src_files = [
        ("SparseQEEcNEO.jl", "Main module and API"),
        ("constants.jl", "Physical and computational constants"),
        ("types.jl", "Data structures and type definitions"),
        ("pyscf_interface.jl", "PySCF quantum chemistry interface"),
        ("epc_functionals.jl", "Electron-proton correlation functionals"),
        ("nuclear_methods.jl", "Nuclear orbital methods"),
        ("configuration_generation.jl", "Configuration selection algorithms"),
        ("importance_analysis.jl", "Configuration importance metrics"),
        ("qee_methods.jl", "Quantum eigensolver methods"),
        ("hamiltonian_construction.jl", "Second-quantized Hamiltonian"),
        ("cneo_methods.jl", "Constrained NEO implementation")
    ]
    
    for (file, desc) in src_files
        if isfile(joinpath(pwd(), "src", file))
            println("  📜 $file - $desc")
        end
    end
    
    println("\n🧪 Test Structure (test/):")
    test_files = readdir(joinpath(pwd(), "test"))
    test_files = sort([f for f in test_files if endswith(f, ".jl")])
    
    for file in test_files
        println("  🧪 $file")
    end
    
    println("\n💡 Examples Structure (examples/):")
    if isdir(joinpath(pwd(), "examples"))
        example_files = readdir(joinpath(pwd(), "examples"))
        example_files = sort([f for f in example_files if endswith(f, ".jl")])
        
        for file in example_files
            println("  💡 $file")
        end
    end
    
    println("\n🚀 Advanced Examples (advanced_examples/):")
    if isdir(joinpath(pwd(), "advanced_examples"))
        for subdir in readdir(joinpath(pwd(), "advanced_examples"))
            subdir_path = joinpath(pwd(), "advanced_examples", subdir)
            if isdir(subdir_path)
                println("  📁 $subdir/")
                files = [f for f in readdir(subdir_path) if endswith(f, ".jl")]
                for file in sort(files)
                    println("    📜 $file")
                end
            end
        end
    end
    
    println("\n🛠️  Development Scripts (scripts/):")
    if isdir(joinpath(pwd(), "scripts"))
        script_files = readdir(joinpath(pwd(), "scripts"))
        script_files = sort([f for f in script_files if endswith(f, ".jl")])
        
        for file in script_files
            println("  🛠️  $file")
        end
    end
    
    println("\n📊 Project Statistics:")
    # Count files by type
    total_jl = length([1 for (root, dirs, files) in walkdir(pwd()) 
                      for file in files if endswith(file, ".jl")])
    total_md = length([1 for (root, dirs, files) in walkdir(pwd()) 
                      for file in files if endswith(file, ".md")])
    
    println("  📜 Julia files (.jl): $total_jl")
    println("  📄 Markdown files (.md): $total_md")
    
    # Calculate total lines of code (rough estimate)
    total_lines = 0
    for (root, dirs, files) in walkdir(joinpath(pwd(), "src"))
        for file in files
            if endswith(file, ".jl")
                filepath = joinpath(root, file)
                total_lines += length(readlines(filepath))
            end
        end
    end
    println("  📊 Lines of code in src/: $total_lines")
    
    println("\n" * "="^60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    show_directory_structure()
end