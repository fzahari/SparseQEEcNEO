#!/usr/bin/env julia
"""
tools/visualize_hamiltonian.jl - Create visualizations of Hamiltonian structure
"""

push!(LOAD_PATH, dirname(@__DIR__))
using SparseQEEcNEO
using LinearAlgebra

# Try to load Plots
has_plots = try
    using Plots
    gr()
    true
catch
    false
end

function text_visualization(H::Matrix, max_size=20)
    """Create text-based visualization for small matrices"""
    n = size(H, 1)
    
    if n > max_size
        println("Matrix too large for text display (size: $n×$n)")
        return
    end
    
    println("\nHamiltonian Matrix (rounded to 3 decimal places):")
    println("Row \\ Col" * join([lpad(j, 9) for j in 1:n]))
    
    for i in 1:n
        row_str = lpad(i, 4) * " |"
        for j in 1:n
            val = H[i,j]
            if abs(val) < 1e-12
                row_str *= "     .   "
            else
                row_str *= lpad(round(val, digits=3), 9)
            end
        end
        println(row_str)
    end
end

function create_visualizations(H::Matrix, configs, output_prefix)
    """Create various visualizations of the Hamiltonian"""
    
    println("\nCreating visualizations...")
    
    # Text-based analysis (always available)
    println("\n1. Matrix Statistics:")
    println("   Size: $(size(H))")
    println("   Sparsity: $(round(count(abs.(H) .< 1e-12) / length(H) * 100, digits=1))%")
    println("   Max element: $(round(maximum(abs.(H)), digits=6))")
    println("   Frobenius norm: $(round(norm(H), digits=6))")
    
    # Eigenvalue spectrum
    eigenvals = eigvals(Hermitian(H))
    n_show = min(20, length(eigenvals))
    
    println("\n2. Eigenvalue Spectrum (first $n_show):")
    for i in 1:n_show
        println("   State $i: $(round(eigenvals[i], digits=6)) Ha")
    end
    
    # Configuration analysis
    println("\n3. Configuration Analysis:")
    config_types = Dict{String, Int}()
    for config in configs
        type = occursin("HF", config.name) ? "Reference" :
               occursin("S(", config.name) || occursin("E(", config.name) ? "Single" :
               occursin("D(", config.name) ? "Double" :
               match(r"N\d+", config.name) !== nothing ? "Nuclear" :
               occursin("C(", config.name) ? "Coupled" : "Other"
        config_types[type] = get(config_types, type, 0) + 1
    end
    
    for (type, count) in sort(collect(config_types), by=x->x[2], rev=true)
        println("   $type: $count configurations")
    end
    
    # Small matrix text display
    if size(H, 1) <= 10
        text_visualization(H, 10)
    end
    
    # Graphical visualizations if Plots is available
    if has_plots
        println("\n4. Creating graphical visualizations...")
        
        # Hamiltonian heatmap
        p1 = heatmap(abs.(H),
                    color=:viridis,
                    title="Hamiltonian Matrix |H_ij|",
                    xlabel="Configuration j",
                    ylabel="Configuration i",
                    aspect_ratio=:equal,
                    colorbar_title="|H_ij|")
        savefig(p1, "$(output_prefix)_hamiltonian.png")
        
        # Sparsity pattern
        sparse_pattern = abs.(H) .> 1e-12
        p2 = heatmap(sparse_pattern,
                    color=[:white, :black],
                    title="Non-zero Pattern",
                    xlabel="Configuration j",
                    ylabel="Configuration i",
                    aspect_ratio=:equal,
                    colorbar=false)
        savefig(p2, "$(output_prefix)_sparsity.png")
        
        # Energy spectrum
        p3 = scatter(1:n_show, eigenvals[1:n_show],
                    xlabel="State Index",
                    ylabel="Energy (Ha)",
                    title="Energy Spectrum",
                    markersize=6,
                    label="Eigenvalues")
        savefig(p3, "$(output_prefix)_spectrum.png")
        
        # Configuration weights
        weights = [c.weight for c in configs]
        p4 = bar(1:min(30, length(weights)), 
                weights[1:min(30, length(weights))],
                xlabel="Configuration Index",
                ylabel="Weight",
                title="Configuration Weights (top 30)",
                label="Weight",
                color=:viridis)
        savefig(p4, "$(output_prefix)_weights.png")
        
        println("\nVisualizations saved:")
        println("   - $(output_prefix)_hamiltonian.png")
        println("   - $(output_prefix)_sparsity.png")
        println("   - $(output_prefix)_spectrum.png")
        println("   - $(output_prefix)_weights.png")
    else
        println("\n⚠ Plots.jl not available. Install with: using Pkg; Pkg.add(\"Plots\")")
        println("  Text-based analysis completed.")
    end
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 0
        println("Usage: julia visualize_hamiltonian.jl <hamiltonian.h5> [output_prefix]")
        println("\nCreates visualizations of Hamiltonian structure")
        exit(1)
    end
    
    filename = ARGS[1]
    output_prefix = length(ARGS) > 1 ? ARGS[2] : splitext(basename(filename))[1]
    
    println("Loading Hamiltonian from: $filename")
    
    # Load data
    ham_data, H = load_hamiltonian(filename)
    configs = ham_data.configs
    
    # Create visualizations
    create_visualizations(H, configs, output_prefix)
    
    println("\nVisualization complete!")
end
