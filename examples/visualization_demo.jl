#!/usr/bin/env julia
"""
examples/visualization_demo.jl - Visualization examples for Sparse QEE-cNEO
Note: This requires Plots.jl to be installed
"""

ENV["PYTHONPATH"] = get(ENV, "PYSCF_PATH", "/path/to/pyscf-master")
push!(LOAD_PATH, dirname(@__DIR__))

using SparseQEEcNEO
using LinearAlgebra

# Check if Plots is available
has_plots = try
    using Plots
    using ColorSchemes
    gr()
    default(size=(800, 600), dpi=150, titlefont=12, guidefont=11, tickfont=10, legendfont=10)
    true
catch
    false
end

if !has_plots
    println("Plots.jl not available. Install with: using Pkg; Pkg.add(\"Plots\")")
    println("Proceeding with text-based output only.")
end

config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))

println("Sparse QEE-cNEO Visualization Examples")
println("="^60)

# ======================== Example 1: Configuration Weight Distribution ========================
println("\nExample 1: Configuration weight distribution")
println("-"^40)

# Run calculation
mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
config_sel = ConfigSelection(method="neo_cneo", max_configs=100, max_nuc_orbs=0)
results = sparse_qee_cneo(mol, config_sel=config_sel, neo_config=config)

# Extract weights
weights = [c.weight for c in results.configs]
names = [c.name for c in results.configs]

println("Top 10 configurations:")
for i in 1:min(10, length(weights))
    println("  $(i). $(names[i]): $(round(weights[i], digits=4))")
end

if has_plots
    # Plot weight distribution
    p1 = bar(1:min(20, length(weights)), 
             weights[1:min(20, length(weights))],
             xlabel="Configuration Index",
             ylabel="Weight",
             title="Configuration Weights",
             label="Weight",
             color=:viridis)
    
    # Log scale plot
    p2 = scatter(1:length(weights),
                weights,
                xlabel="Configuration Index",
                ylabel="Weight (log scale)",
                title="Weight Distribution",
                yscale=:log10,
                label="Configurations",
                markersize=4,
                alpha=0.7)
    
    plot(p1, p2, layout=(1,2), size=(1200, 400))
    savefig("config_weights.png")
    println("\nVisualization saved to config_weights.png")
end

# ======================== Example 2: Hamiltonian Matrix Visualization ========================
println("\nExample 2: Hamiltonian matrix structure")
println("-"^40)

H = results.hamiltonian_matrix
println("Hamiltonian size: $(size(H))")
println("Sparsity: $(round(count(abs.(H) .< 1e-12) / length(H) * 100, digits=1))%")

if has_plots && size(H, 1) <= 100  # Only plot if reasonable size
    # Heatmap of absolute values
    p3 = heatmap(abs.(H),
                color=:viridis,
                title="Hamiltonian Matrix (|H_ij|)",
                xlabel="Config j",
                ylabel="Config i",
                aspect_ratio=:equal,
                colorbar_title="|H_ij|")
    
    # Sparsity pattern
    sparse_pattern = abs.(H) .> 1e-12
    p4 = heatmap(sparse_pattern,
                color=[:white, :black],
                title="Non-zero Pattern",
                xlabel="Config j",
                ylabel="Config i",
                aspect_ratio=:equal,
                colorbar=false)
    
    plot(p3, p4, layout=(1,2), size=(1000, 400))
    savefig("hamiltonian_structure.png")
    println("Visualization saved to hamiltonian_structure.png")
end

# ======================== Example 3: Energy Spectrum ========================
println("\nExample 3: Energy spectrum")
println("-"^40)

eigenvals = eigvals(Hermitian(H))
n_show = min(20, length(eigenvals))

println("First $n_show eigenvalues:")
for i in 1:n_show
    println("  State $i: $(round(eigenvals[i], digits=6)) Ha")
end

if has_plots
    # Energy level diagram
    p5 = scatter(ones(n_show), eigenvals[1:n_show],
                xlabel="",
                ylabel="Energy (Ha)",
                title="Energy Spectrum",
                markersize=8,
                label="Energy levels",
                xlims=(0.5, 1.5),
                xticks=false)
    
    # Add lines for each level
    for i in 1:n_show
        plot!([0.8, 1.2], [eigenvals[i], eigenvals[i]], 
              label=false, color=:black, alpha=0.5)
    end
    
    # Energy gaps
    gaps = diff(eigenvals[1:n_show])
    p6 = bar(1:length(gaps), gaps,
            xlabel="Gap Index",
            ylabel="Energy Gap (Ha)",
            title="Energy Gaps",
            label="ΔE",
            color=:orange)
    
    plot(p5, p6, layout=(1,2), size=(1000, 400))
    savefig("energy_spectrum.png")
    println("Visualization saved to energy_spectrum.png")
end

# ======================== Example 4: Method Comparison ========================
println("\nExample 4: Method comparison")
println("-"^40)

methods = ["mp2", "neo_cneo", "neo_enhanced"]
method_data = Dict()

for method in methods
    config_sel = ConfigSelection(method=method, max_configs=50, max_nuc_orbs=0)
    try
        res = sparse_qee_cneo(mol, config_sel=config_sel, neo_config=config)
        method_data[method] = (
            energy = res.total_energy,
            n_configs = res.n_configs,
            importance = res.captured_importance
        )
        println("$method: E = $(round(res.total_energy, digits=6)) Ha")
    catch e
        println("$method: Failed - $e")
    end
end

if has_plots && length(method_data) > 1
    # Comparison plots
    methods_plot = collect(keys(method_data))
    energies = [method_data[m].energy for m in methods_plot]
    n_configs = [method_data[m].n_configs for m in methods_plot]
    importance = [method_data[m].importance * 100 for m in methods_plot]
    
    p7 = bar(methods_plot, energies,
            ylabel="Energy (Ha)",
            title="Energy Comparison",
            label="Total Energy",
            color=:blues)
    
    p8 = bar(methods_plot, n_configs,
            ylabel="Configurations",
            title="Configuration Count",
            label="N configs",
            color=:greens)
    
    p9 = bar(methods_plot, importance,
            ylabel="Importance (%)",
            title="Captured Importance",
            label="Importance",
            color=:reds,
            ylims=(0, 105))
    
    plot(p7, p8, p9, layout=(1,3), size=(1400, 400))
    savefig("method_comparison.png")
    println("Visualization saved to method_comparison.png")
end

# ======================== Example 5: Configuration Space Analysis ========================
println("\nExample 5: Configuration space analysis")
println("-"^40)

# Analyze configuration types
config_types = Dict{String, Int}()
for config in results.configs
    type_key = if occursin("HF", config.name) || occursin("REF", config.name)
        "Reference"
    elseif occursin("C(", config.name)
        "Coupled"
    elseif match(r"N\d+", config.name) !== nothing
        "Nuclear"
    elseif occursin("S(", config.name) || occursin("E(", config.name)
        "Single"
    elseif occursin("D(", config.name)
        "Double"
    else
        "Other"
    end
    
    config_types[type_key] = get(config_types, type_key, 0) + 1
end

println("Configuration types:")
for (type, count) in sort(collect(config_types), by=x->x[2], rev=true)
    println("  $type: $count ($(round(count/length(results.configs)*100, digits=1))%)")
end

if has_plots && !isempty(config_types)
    # Pie chart of configuration types
    types = collect(keys(config_types))
    counts = [config_types[t] for t in types]
    
    pie(types, counts,
        title="Configuration Types",
        label=true,
        legend=:right)
    
    savefig("config_types.png")
    println("Visualization saved to config_types.png")
end

println("\n" * "="^60)
println("Visualization examples completed!")
if has_plots
    println("All visualizations saved as PNG files.")
else
    println("Install Plots.jl for graphical visualizations.")
end
