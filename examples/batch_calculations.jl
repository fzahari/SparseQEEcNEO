#!/usr/bin/env julia
"""
examples/batch_calculations.jl - Run batch calculations and generate reports
"""

ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/path/to/pyscf-master")
push!(LOAD_PATH, dirname(@__DIR__))

using SparseQEEcNEO
using DataFrames
using CSV
using Printf

function run_batch_calculations(molecules, methods, output_dir="results")
    # Create output directory
    mkpath(output_dir)
    
    # NEO configuration
    config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
    
    # Results storage
    results = DataFrame(
        molecule = String[],
        method = String[],
        energy = Float64[],
        n_configs = Int[],
        n_qubits = Int[],
        importance = Float64[],
        time = Float64[],
        has_neo = Bool[]
    )
    
    # Run calculations
    for (mol_name, mol) in molecules
        has_quantum_nuc = !isempty(mol.quantum_nuc)
        
        for method in methods
            println("\nCalculating $mol_name with $method...")
            
            config_sel = ConfigSelection(method=method, max_configs=100, max_nuc_orbs=0)
            
            try
                t_start = time()
                res = sparse_qee_cneo(mol, config_sel=config_sel, neo_config=config)
                t_end = time()
                
                push!(results, (
                    molecule = mol_name,
                    method = method,
                    energy = res.total_energy,
                    n_configs = res.n_configs,
                    n_qubits = res.n_qubits,
                    importance = res.captured_importance,
                    time = t_end - t_start,
                    has_neo = has_quantum_nuc
                ))
                
                # Save Hamiltonian
                save_hamiltonian(
                    joinpath(output_dir, "$(mol_name)_$(method).h5"),
                    res.hamiltonian_data,
                    res.hamiltonian_matrix
                )
                
                println("  ✓ Success: E = $(round(res.total_energy, digits=6)) Ha")
                
            catch e
                println("  ✗ Failed: $e")
                push!(results, (
                    molecule = mol_name,
                    method = method,
                    energy = NaN,
                    n_configs = 0,
                    n_qubits = 0,
                    importance = 0.0,
                    time = 0.0,
                    has_neo = has_quantum_nuc
                ))
            end
        end
    end
    
    # Save results
    CSV.write(joinpath(output_dir, "results.csv"), results)
    
    # Generate summary report
    generate_report(results, output_dir)
    
    return results
end

function generate_report(results, output_dir)
    open(joinpath(output_dir, "summary_report.txt"), "w") do io
        println(io, "Sparse QEE-cNEO Batch Calculation Results")
        println(io, "="^60)
        println(io, "Generated: $(Dates.now())")
        println(io)
        
        # Summary by molecule
        println(io, "\nResults by Molecule:")
        println(io, "-"^60)
        
        for mol in unique(results.molecule)
            mol_results = results[results.molecule .== mol, :]
            println(io, "\n$mol:")
            
            for row in eachrow(mol_results)
                if !isnan(row.energy)
                    println(io, "  $(row.method): E = $(round(row.energy, digits=6)) Ha, " *
                              "$(row.n_configs) configs, $(round(row.time, digits=2))s")
                else
                    println(io, "  $(row.method): FAILED")
                end
            end
        end
        
        # Summary by method
        println(io, "\n\nResults by Method:")
        println(io, "-"^60)
        
        for method in unique(results.method)
            method_results = results[results.method .== method, :]
            successful = method_results[.!isnan.(method_results.energy), :]
            
            println(io, "\n$method:")
            println(io, "  Success rate: $(nrow(successful))/$(nrow(method_results))")
            
            if nrow(successful) > 0
                println(io, "  Average configs: $(round(mean(successful.n_configs), digits=1))")
                println(io, "  Average time: $(round(mean(successful.time), digits=2))s")
                println(io, "  Average importance: $(round(mean(successful.importance) * 100, digits=1))%")
            end
        end
        
        # NEO vs Classical comparison
        println(io, "\n\nNEO vs Classical Comparison:")
        println(io, "-"^60)
        
        neo_results = results[results.has_neo, :]
        classical_results = results[.!results.has_neo, :]
        
        println(io, "\nNEO calculations: $(nrow(neo_results))")
        println(io, "Classical calculations: $(nrow(classical_results))")
    end
    
    println("\nReport saved to $(joinpath(output_dir, "summary_report.txt"))")
end

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    # Define test molecules
    molecules = [
        ("H2", Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")),
        ("H2_NEO", Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])),
        ("HeH+", Molecule("He 0 0 0; H 0 0 1.4632", "sto-3g", charge=1)),
        ("HeH+_NEO", Molecule("He 0 0 0; H 0 0 1.4632", "sto-3g", charge=1, quantum_nuc=[1])),
        ("Water", Molecule("O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", "6-31g")),
        ("Water_NEO", Molecule("O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0", "6-31g", quantum_nuc=[1])),
        ("HCN", Molecule("H 0 0 -1.064; C 0 0 0; N 0 0 1.156", "sto-3g")),
        ("HCN_NEO", Molecule("H 0 0 -1.064; C 0 0 0; N 0 0 1.156", "sto-3g", quantum_nuc=[0]))
    ]
    
    methods = ["mp2", "neo_cneo", "hybrid_final"]
    
    results = run_batch_calculations(molecules, methods)
    println("\nBatch calculations complete! Results saved in 'results' directory.")
    
    # Print summary
    println("\nQuick Summary:")
    println("-"^40)
    for row in eachrow(results)
        if !isnan(row.energy)
            println("$(rpad(row.molecule, 12)) $(rpad(row.method, 12)): " *
                   "E = $(round(row.energy, digits=6)) Ha")
        end
    end
end
