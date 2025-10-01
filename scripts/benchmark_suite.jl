#!/usr/bin/env julia
"""
tools/benchmark_suite.jl - Performance benchmarking suite for Sparse QEE-cNEO
"""

ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/path/to/pyscf-master")
push!(LOAD_PATH, dirname(@__DIR__))

using SparseQEEcNEO
using Statistics
using Printf

# Optional packages
has_benchmarktools = try
    using BenchmarkTools
    true
catch
    false
end

has_dataframes = try
    using DataFrames
    using CSV
    true
catch
    false
end

# Simple timing macro if BenchmarkTools not available
if !has_benchmarktools
    macro benchmark(expr)
        quote
            times = Float64[]
            for _ in 1:5
                t0 = time()
                $(esc(expr))
                push!(times, time() - t0)
            end
            (times = times * 1e9, memory = 0)  # Convert to nanoseconds
        end
    end
end

function benchmark_configuration_generation()
    println("\nBenchmarking Configuration Generation")
    println("-"^50)
    
    # Test different system sizes
    n_orbitals_list = [4, 8, 16]
    methods = ["mp2", "neo_cneo"]
    
    results = []
    
    for n_orbs in n_orbitals_list
        # Create mock mean-field
        mf = create_mock_mf(n_orbs)
        mol = Molecule("Mock", "sto-3g")
        
        for method in methods
            config_sel = ConfigSelection(method=method, max_configs=100, max_nuc_orbs=0)
            
            # Time the function
            if has_benchmarktools
                b = @benchmark SparseQEEcNEO.ConfigurationGeneration.generate_configurations(
                    $mf, $mol, $config_sel, nothing
                )
                time_mean = mean(b.times) / 1e9  # Convert to seconds
                time_std = std(b.times) / 1e9
                memory = b.memory / 1024^2  # Convert to MB
            else
                times = Float64[]
                for _ in 1:5
                    t0 = time()
                    SparseQEEcNEO.ConfigurationGeneration.generate_configurations(
                        mf, mol, config_sel, nothing
                    )
                    push!(times, time() - t0)
                end
                time_mean = mean(times)
                time_std = std(times)
                memory = 0.0
            end
            
            # Get configurations for counting
            configs = SparseQEEcNEO.ConfigurationGeneration.generate_configurations(
                mf, mol, config_sel, nothing
            )
            
            result = (
                method = method,
                n_orbitals = n_orbs,
                time_mean = time_mean,
                time_std = time_std,
                memory = memory,
                n_configs = length(configs)
            )
            push!(results, result)
            
            println("  $method ($n_orbs orbitals): $(round(time_mean*1000, digits=2)) ms")
        end
    end
    
    return results
end

function benchmark_hamiltonian_construction()
    println("\nBenchmarking Hamiltonian Construction")
    println("-"^50)
    
    config_sizes = [10, 50, 100]
    n_orbs = 10
    
    results = []
    
    for n_configs in config_sizes
        # Create mock configurations
        configs = [
            Configuration("C$i", Dict("e" => create_random_occupation(n_orbs)), rand())
            for i in 1:n_configs
        ]
        
        # Create mock integrals
        h1e = Dict("e" => randn(n_orbs, n_orbs))
        h1e["e"] = (h1e["e"] + h1e["e"]') / 2
        h2e = Dict("e" => SparseQEEcNEO.HamiltonianConstruction.create_sparse_eri(n_orbs))
        coupling = Dict{Tuple{String,String}, Matrix{Float64}}()
        active = Dict("e" => collect(1:n_orbs))
        
        # Time the function
        times = Float64[]
        for _ in 1:5
            t0 = time()
            H = SparseQEEcNEO.HamiltonianConstruction.build_configuration_hamiltonian(
                configs, h1e, h2e, coupling, active
            )
            push!(times, time() - t0)
        end
        
        # Get actual Hamiltonian for analysis
        H = SparseQEEcNEO.HamiltonianConstruction.build_configuration_hamiltonian(
            configs, h1e, h2e, coupling, active
        )
        
        sparsity = count(abs.(H) .< 1e-12) / length(H)
        
        result = (
            n_configs = n_configs,
            n_orbitals = n_orbs,
            time_mean = mean(times),
            time_std = std(times),
            sparsity = sparsity
        )
        push!(results, result)
        
        println("  $n_configs configs: $(round(mean(times)*1000, digits=2)) ms, " *
                "sparsity: $(round(sparsity*100, digits=1))%")
    end
    
    return results
end

function benchmark_full_workflow()
    println("\nBenchmarking Full Workflow")
    println("-"^50)
    
    # Check if PySCF is available
    config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
    pyscf, has_neo = try
        SparseQEEcNEO.setup_pyscf(config)
    catch
        (nothing, false)
    end
    
    if !has_neo
        println("  PySCF NEO not available - skipping workflow benchmark")
        return []
    end
   
# Test molecules
    molecules = [
        ("H2", Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")),
        ("H2_NEO", Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])),
        ("HeH+", Molecule("He 0 0 0; H 0 0 1.4", "sto-3g", charge=1))
    ]
    
    methods = ["mp2", "neo_cneo"]
    
    results = []
    
    for (mol_name, mol) in molecules
        for method in methods
            println("  Benchmarking $mol_name with $method...")
            
            config_sel = ConfigSelection(method=method, max_configs=50, max_nuc_orbs=0)
            
            # Time full calculation
            t_start = time()
            mem_start = Base.gc_live_bytes()
            
            try
                res = sparse_qee_cneo(mol, config_sel=config_sel, neo_config=config)
                
                t_total = time() - t_start
                mem_peak = (Base.gc_live_bytes() - mem_start) / 1024^2
                
                result = (
                    molecule = mol_name,
                    method = method,
                    time_total = t_total,
                    memory_peak = mem_peak,
                    n_configs = res.n_configs,
                    n_qubits = res.n_qubits,
                    energy = res.total_energy
                )
                push!(results, result)
                
                println("    ✓ Time: $(round(t_total, digits=2))s, " *
                       "Configs: $(res.n_configs), " *
                       "Memory: $(round(mem_peak, digits=1)) MB")
            catch e
                println("    ✗ Failed: $e")
            end
        end
    end
    
    return results
end

# Helper functions
function create_mock_mf(n_orbs)
    mf = Dict()
    mf_e = Dict()
    mf_e["mo_occ"] = [i <= n_orbs÷2 ? 1.0 : 0.0 for i in 1:n_orbs]
    mf_e["mo_energy"] = collect(range(-2.0, 2.0, length=n_orbs))
    mf["components"] = Dict("e" => mf_e)
    mf["e_tot"] = -1.0
    return mf
end

function create_random_occupation(n_orbs, n_occupied=nothing)
    if n_occupied === nothing
        n_occupied = n_orbs ÷ 2
    end
    occ = zeros(Int16, n_orbs)
    occupied_indices = randperm(n_orbs)[1:n_occupied]
    occ[occupied_indices] .= 1
    return occ
end

function print_benchmark_summary(results_dict)
    println("\n" * "="^60)
    println("Benchmark Summary")
    println("="^60)
    
    if haskey(results_dict, "config_gen") && !isempty(results_dict["config_gen"])
        println("\nConfiguration Generation:")
        for r in results_dict["config_gen"]
            println("  $(r.method) ($(r.n_orbitals) orbs): " *
                   "$(round(r.time_mean*1000, digits=2)) ± " *
                   "$(round(r.time_std*1000, digits=2)) ms")
        end
    end
    
    if haskey(results_dict, "hamiltonian") && !isempty(results_dict["hamiltonian"])
        println("\nHamiltonian Construction:")
        for r in results_dict["hamiltonian"]
            println("  $(r.n_configs) configs: " *
                   "$(round(r.time_mean*1000, digits=2)) ± " *
                   "$(round(r.time_std*1000, digits=2)) ms " *
                   "($(round(r.sparsity*100, digits=1))% sparse)")
        end
    end
    
    if haskey(results_dict, "workflow") && !isempty(results_dict["workflow"])
        println("\nFull Workflow:")
        for r in results_dict["workflow"]
            println("  $(r.molecule) - $(r.method): " *
                   "$(round(r.time_total, digits=2))s, " *
                   "$(r.n_configs) configs")
        end
    end
end

function save_benchmark_results(results_dict, output_dir="benchmark_results")
    mkpath(output_dir)
    
    # Save text summary
    open(joinpath(output_dir, "benchmark_summary.txt"), "w") do io
        println(io, "SparseQEEcNEO Benchmark Results")
        println(io, "Generated: $(Dates.now())")
        println(io, "="^60)
        
        for (name, results) in results_dict
            if !isempty(results)
                println(io, "\n$name:")
                for r in results
                    println(io, "  $r")
                end
            end
        end
    end
    
    # Save as CSV if DataFrames available
    if has_dataframes
        for (name, results) in results_dict
            if !isempty(results)
                df = DataFrame(results)
                CSV.write(joinpath(output_dir, "$(name).csv"), df)
            end
        end
        println("\nResults saved to $output_dir/")
    end
end

# Main benchmark suite
function run_benchmark_suite(; save_results=true)
    println("\nSparse QEE-cNEO Performance Benchmark Suite")
    println("="^60)
    
    if !has_benchmarktools
        println("⚠ BenchmarkTools.jl not available - using simple timing")
        println("  Install with: using Pkg; Pkg.add(\"BenchmarkTools\")")
    end
    
    results_dict = Dict()
    
    # Run benchmarks
    println("\n1. Configuration Generation")
    results_dict["config_gen"] = benchmark_configuration_generation()
    
    println("\n2. Hamiltonian Construction")
    results_dict["hamiltonian"] = benchmark_hamiltonian_construction()
    
    println("\n3. Full Workflow")
    results_dict["workflow"] = benchmark_full_workflow()
    
    # Print summary
    print_benchmark_summary(results_dict)
    
    # Save results
    if save_results
        save_benchmark_results(results_dict)
    end
    
    return results_dict
end

# Run if called directly
if abspath(PROGRAM_FILE) == @__FILE__
    results = run_benchmark_suite()
    
    println("\n" * "="^60)
    println("Benchmarking complete!")
    
    # Performance recommendations
    println("\nPerformance Tips:")
    println("1. Set JULIA_NUM_THREADS for parallel execution")
    println("2. Use max_nuc_orbs to limit nuclear orbital space")
    println("3. Adjust importance_cutoff to reduce configurations")
    println("4. Enable use_compression for large systems")
end 
