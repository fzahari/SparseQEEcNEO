#!/usr/bin/env julia
"""
tools/analyze_configurations.jl - Analyze configuration space from saved Hamiltonians
"""

push!(LOAD_PATH, dirname(@__DIR__))
using SparseQEEcNEO
using Statistics
using LinearAlgebra

function analyze_configuration_space(configs)
    # Count types
    types = Dict{String, Int}()
    
    for config in configs
        type_key = get_config_type(config)
        types[type_key] = get(types, type_key, 0) + 1
    end
    
    # Weight distribution
    weights = [get_weight(config) for config in configs]
    
    println("\nConfiguration Space Analysis")
    println("="^50)
    println("Total configurations: $(length(configs))")
    
    println("\nConfiguration types:")
    for (type, count) in sort(collect(types), by=x->x[2], rev=true)
        percentage = round(count/length(configs)*100, digits=1)
        println("  $(rpad(type, 12)): $(lpad(count, 4)) ($(lpad(percentage, 5))%)")
    end
    
    println("\nWeight statistics:")
    println("  Maximum:     $(round(maximum(weights), digits=6))")
    println("  Minimum:     $(round(minimum(weights), digits=6))")
    println("  Mean:        $(round(mean(weights), digits=6))")
    println("  Median:      $(round(median(weights), digits=6))")
    println("  Std dev:     $(round(std(weights), digits=6))")
    
    # Cumulative importance
    sorted_weights = sort(weights, rev=true)
    cumsum_weights = cumsum(sorted_weights) / sum(sorted_weights)
    
    n_50 = findfirst(x -> x >= 0.5, cumsum_weights)
    n_90 = findfirst(x -> x >= 0.9, cumsum_weights)
    n_95 = findfirst(x -> x >= 0.95, cumsum_weights)
    n_99 = findfirst(x -> x >= 0.99, cumsum_weights)
    
    println("\nCumulative importance:")
    println("  50% from top $(lpad(n_50, 4)) configs ($(round(n_50/length(configs)*100, digits=1))%)")
    println("  90% from top $(lpad(n_90, 4)) configs ($(round(n_90/length(configs)*100, digits=1))%)")
    println("  95% from top $(lpad(n_95, 4)) configs ($(round(n_95/length(configs)*100, digits=1))%)")
    println("  99% from top $(lpad(n_99, 4)) configs ($(round(n_99/length(configs)*100, digits=1))%)")
    
    # Effective rank
    probs = weights / sum(weights)
    entropy = -sum(p * log(p) for p in probs if p > 0)
    effective_rank = exp(entropy)
    println("\nEffective rank: $(round(effective_rank, digits=2))")
    
    return types, weights, cumsum_weights
end

function get_config_type(config)
    name = config.name
    
    if occursin("HF", name) || occursin("REF", name)
        return "Reference"
    elseif occursin("C(", name)
        return "Coupled"
    elseif match(r"N\d+\(", name) !== nothing
        return "Nuclear"
    elseif occursin("S(", name) || occursin("E(", name)
        return "Single"
    elseif occursin("D(", name)
        return "Double"
    else
        return "Other"
    end
end

function get_weight(config)
    return config.weight
end

function analyze_hamiltonian_structure(H, configs)
    println("\nHamiltonian Structure Analysis")
    println("="^50)
    
    # Basic properties
    println("Matrix size: $(size(H))")
    println("Hermitian: $(ishermitian(H))")
    
    # Sparsity analysis
    n_zero = count(abs.(H) .< 1e-12)
    sparsity = n_zero / length(H) * 100
    println("Sparsity: $(round(sparsity, digits=1))%")
    
    # Norm analysis
    println("\nMatrix norms:")
    println("  Frobenius norm: $(round(norm(H), digits=6))")
    println("  Max element: $(round(maximum(abs.(H)), digits=6))")
    println("  Min non-zero: $(round(minimum(abs.(H[abs.(H) .> 1e-12])), sigdigits=3))")
    
    # Diagonal dominance
    diag_sum = sum(abs.(diag(H)))
    total_sum = sum(abs.(H))
    diag_fraction = diag_sum / total_sum * 100
    println("\nDiagonal contribution: $(round(diag_fraction, digits=1))%")
    
    # Block structure analysis
    println("\nBlock structure:")
    
    # Group configurations by type
    type_indices = Dict{String, Vector{Int}}()
    for (i, config) in enumerate(configs)
        type = get_config_type(config)
        if !haskey(type_indices, type)
            type_indices[type] = Int[]
        end
        push!(type_indices[type], i)
    end
    
    # Analyze blocks
    for (type1, indices1) in type_indices
        for (type2, indices2) in type_indices
            if !isempty(indices1) && !isempty(indices2)
                block = H[indices1, indices2]
                block_norm = norm(block)
                if block_norm > 1e-12
                    println("  $type1-$type2: norm = $(round(block_norm, digits=4))")
                end
            end
        end
    end
    
    # Eigenvalue analysis
    println("\nEigenvalue Analysis:")
    eigenvals = eigvals(Hermitian(H))
    
    println("  Ground state: $(round(eigenvals[1], digits=8)) Ha")
    if length(eigenvals) > 1
        println("  First excited: $(round(eigenvals[2], digits=8)) Ha")
        println("  HOMO-LUMO gap: $(round(eigenvals[2] - eigenvals[1], digits=6)) Ha")
    end
    
    if length(eigenvals) >= 10
        println("  10th eigenvalue: $(round(eigenvals[10], digits=6)) Ha")
    end
    
    # Spectral radius
    spectral_radius = maximum(abs.(eigenvals))
    println("  Spectral radius: $(round(spectral_radius, digits=6))")
    
    # Condition number
    if length(eigenvals) > 1
        condition = maximum(abs.(eigenvals)) / minimum(abs.(eigenvals[eigenvals .!= 0]))
        println("  Condition number: $(round(condition, sigdigits=3))")
    end
    
    return eigenvals
end

function export_configuration_data(configs, filename="configurations.txt")
    open(filename, "w") do io
        println(io, "# Configuration data from SparseQEEcNEO")
        println(io, "# Format: index name weight type")
        println(io, "#")
        
        for (i, config) in enumerate(configs)
            type = get_config_type(config)
            weight = get_weight(config)
            println(io, "$i $(config.name) $weight $type")
        end
    end
    
    println("\nConfiguration data exported to: $filename")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0
        println("Usage: julia analyze_configurations.jl <hamiltonian.h5> [output_prefix]")
        println("\nAnalyzes configuration space and Hamiltonian structure from saved file")
        exit(1)
    end
    
    filename = ARGS[1]
    output_prefix = length(ARGS) > 1 ? ARGS[2] : splitext(basename(filename))[1]
    
    println("Loading Hamiltonian from: $filename")
    
    # Load data
    ham_data, H = load_hamiltonian(filename)
    configs = ham_data.configs
    
    # Perform analyses
    types, weights, cumsum_weights = analyze_configuration_space(configs)
    eigenvals = analyze_hamiltonian_structure(H, configs)
    
    # Export data
    export_configuration_data(configs, "$(output_prefix)_configs.txt")
    
    # Save analysis results
    open("$(output_prefix)_analysis.txt", "w") do io
        println(io, "Configuration Space Analysis for $filename")
        println(io, "="^60)
        println(io, "Total configurations: $(length(configs))")
        println(io, "Hamiltonian size: $(size(H))")
        println(io, "Ground state energy: $(round(eigenvals[1], digits=8)) Ha")
        
        if length(eigenvals) > 1
            println(io, "HOMO-LUMO gap: $(round(eigenvals[2] - eigenvals[1], digits=6)) Ha")
        end
        
        println(io, "\nConfiguration types:")
        for (type, count) in sort(collect(types), by=x->x[2], rev=true)
            println(io, "  $type: $count")
        end
        
        println(io, "\nWeight distribution:")
        println(io, "  Max: $(round(maximum(weights), digits=6))")
        println(io, "  Min: $(round(minimum(weights), digits=6))")
        println(io, "  Mean: $(round(mean(weights), digits=6))")
    end
    
    println("\nAnalysis complete!")
    println("Results saved to:")
    println("  - $(output_prefix)_configs.txt")
    println("  - $(output_prefix)_analysis.txt")
end
