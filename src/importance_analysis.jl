"""
importance_analysis.jl - Methods for analyzing configuration importance
"""
module ImportanceAnalysis

using Statistics
using ..Types: Configuration, CompressedConfig, ConfigSelection

export calculate_importance_metrics, analyze_configuration_types
export calculate_neo_importance_metrics

# ======================== Main Importance Calculation ========================

function calculateImportanceMetrics(configs, molNeo, configSelection::ConfigSelection)
    """Calculate comprehensive importance metrics for configurations"""
    
    if isemptyConfigurationSet(configs)
        return createEmptyImportanceResult()
    end
    
    standardImportance = calculateStandardImportance(configs)
    configurationTypes = analyzeConfigurationTypes(configs)
    
    neoMetrics, totalImportance = processNEOMetrics(
        configs, molNeo, configSelection, configurationTypes, standardImportance
    )
    
    logImportanceResults(standardImportance, configs, configurationTypes, neoMetrics)
    
    return buildImportanceResult(totalImportance, neoMetrics, 
                               configurationTypes, standardImportance)
end

function isEmptyConfigurationSet(configs)
    return isempty(configs)
end

function createEmptyImportanceResult()
    return (
        total_importance = 0.0,
        neo_metrics = nothing,
        config_types = Dict{String, Int}()
    )
end

function calculateStandardImportance(configs)
    weights = [get_weight(c) for c in configs]
    totalWeight = sum(weights)
    normalizedWeights = weights ./ max(totalWeight, 1.0)
    return min(sum(normalizedWeights), 1.0)
end

function processNEOMetrics(configs, molNeo, configSelection, configTypes, standardImportance)
    if shouldCalculateNEOMetrics(configSelection, configs)
        neoMetrics = calculate_neo_importance_metrics(configs, molNeo, configTypes)
        enhancementFactor = calculate_neo_enhancement(neoMetrics)
        totalImportance = min(standardImportance * enhancementFactor, 1.0)
        return neoMetrics, totalImportance
    else
        return nothing, standardImportance
    end
end

function shouldCalculateNEOMetrics(configSelection, configs)
    return configSelection.use_neo_importance && hasNuclearConfigs(configs)
end

function getConfigName(config)
    return get_config_name(config)
end

function hasNuclearConfigs(configs)
    return has_nuclear_configs(configs)
end

function logImportanceResults(standardImportance, configs, configTypes, neoMetrics)
    @info "Importance analysis complete:" *
          "\n  Standard importance: $(round(standardImportance, digits=3))" *
          "\n  Total configurations: $(length(configs))" *
          "\n  Configuration types: $(configTypes)"
    
    if neoMetrics !== nothing
        logNEOMetrics(neoMetrics)
    end
end

function logNEOMetrics(neoMetrics)
    @info "NEO metrics:" *
          "\n  Nuclear participation: $(round(neoMetrics.nuclear_participation * 100, digits=1))%" *
          "\n  Coupling contribution: $(round(neoMetrics.coupling_contribution, digits=3))" *
          "\n  Enhancement factor: $(round(neoMetrics.enhancement_factor, digits=3))"
end

function buildImportanceResult(totalImportance, neoMetrics, configTypes, standardImportance)
    return (
        total_importance = totalImportance,
        neo_metrics = neoMetrics,
        config_types = configTypes,
        standard_importance = standardImportance
    )
end

# Backward compatibility wrapper
function calculate_importance_metrics(configs, molNeo, configSelection::ConfigSelection)
    return calculateImportanceMetrics(configs, molNeo, configSelection)
end

# ======================== Configuration Type Analysis ========================

function analyzeConfigurationTypes(configs)
    """Analyze and categorize configuration types"""
    typeCounts = initializeConfigurationTypeCounts()
    
    for config in configs
        configName = getConfigName(config)
        configType = determineConfigurationType(configName)
        typeCounts[configType] += 1
    end
    
    return typeCounts
end

function initializeConfigurationTypeCounts()
    return Dict{String, Int}(
        "reference" => 0,
        "single_elec" => 0,
        "double_elec" => 0,
        "single_nuc" => 0,
        "double_nuc" => 0,
        "coupled" => 0,
        "other" => 0
    )
end

function determineConfigurationType(configName::String)
    if isReferenceConfiguration(configName)
        return "reference"
    elseif isCoupledConfiguration(configName)
        return "coupled"
    elseif isSingleElectronicConfiguration(configName)
        return "single_elec"
    elseif isDoubleElectronicConfiguration(configName)
        return "double_elec"
    elseif isSingleNuclearConfiguration(configName)
        return "single_nuc"
    elseif isDoubleNuclearConfiguration(configName)
        return "double_nuc"
    else
        return "other"
    end
end

function isReferenceConfiguration(name::String)
    return occursin("HF", name) || occursin("REF", name)
end

function isCoupledConfiguration(name::String)
    return occursin("C(", name)
end

function isSingleElectronicConfiguration(name::String)
    return occursin("S(", name) || occursin("E(", name)
end

function isDoubleElectronicConfiguration(name::String)
    return occursin("D(", name)
end

function isSingleNuclearConfiguration(name::String)
    return match(r"N\d+\(\d+→\d+\)", name) !== nothing
end

function isDoubleNuclearConfiguration(name::String)
    return match(r"N\d+\(\d+\d+→\d+\d+\)", name) !== nothing
end

# Backward compatibility wrapper
function analyze_configuration_types(configs)
    return analyzeConfigurationTypes(configs)
end

# ======================== NEO-specific Importance ========================

function calculate_neo_importance_metrics(configs, mol_neo, config_types)
    """
    Calculate NEO-specific importance metrics
    """
    # Nuclear participation
    n_nuclear = config_types["single_nuc"] + config_types["double_nuc"] + config_types["coupled"]
    n_total = length(configs)
    nuclear_participation = n_nuclear / max(n_total, 1)
    
    # Coupling contribution
    n_coupled = config_types["coupled"]
    coupling_contribution = n_coupled / max(n_total, 1)
    
    # Weight analysis
    nuclear_weights = Float64[]
    coupled_weights = Float64[]
    electronic_weights = Float64[]
    
    for config in configs
        weight = get_weight(config)
        name = get_config_name(config)
        
        if occursin("C(", name)
            push!(coupled_weights, weight)
        elseif match(r"N\d+", name) !== nothing
            push!(nuclear_weights, weight)
        elseif occursin("S(", name) || occursin("D(", name) || occursin("E(", name)
            push!(electronic_weights, weight)
        end
    end
    
    # Calculate weight fractions
    total_weight = sum(get_weight(c) for c in configs)
    nuclear_weight_fraction = sum(nuclear_weights) / max(total_weight, 1e-10)
    coupled_weight_fraction = sum(coupled_weights) / max(total_weight, 1e-10)
    
    # Enhancement factor based on NEO characteristics
    enhancement_factor = calculate_enhancement_factor(
        nuclear_participation,
        coupling_contribution,
        nuclear_weight_fraction,
        coupled_weight_fraction
    )
    
    return (
        nuclear_participation = nuclear_participation,
        coupling_contribution = coupling_contribution,
        nuclear_weight_fraction = nuclear_weight_fraction,
        coupled_weight_fraction = coupled_weight_fraction,
        enhancement_factor = enhancement_factor,
        n_nuclear_configs = n_nuclear,
        n_coupled_configs = n_coupled,
        avg_nuclear_weight = isempty(nuclear_weights) ? 0.0 : mean(nuclear_weights),
        avg_coupled_weight = isempty(coupled_weights) ? 0.0 : mean(coupled_weights)
    )
end

function calculate_enhancement_factor(nuclear_part, coupling_cont, nuclear_wt, coupled_wt)
    """
    Calculate importance enhancement factor for NEO systems
    """
    # Base enhancement
    enhancement = 1.0
    
    # Nuclear participation enhancement
    if nuclear_part > 0.1
        enhancement *= (1.0 + nuclear_part)
    end
    
    # Coupling enhancement
    if coupling_cont > 0.05
        enhancement *= (1.0 + 2.0 * coupling_cont)
    end
    
    # Weight-based enhancement
    if nuclear_wt > 0.1
        enhancement *= (1.0 + 0.5 * nuclear_wt)
    end
    
    if coupled_wt > 0.05
        enhancement *= (1.0 + coupled_wt)
    end
    
    # Cap enhancement
    return min(enhancement, 10.0)
end

function calculate_neo_enhancement(neo_metrics)
    """
    Calculate total NEO enhancement factor
    """
    if neo_metrics === nothing
        return 1.0
    end
    
    return neo_metrics.enhancement_factor
end

# ======================== Helper Functions ========================

function get_weight(config)
    """
    Get weight from configuration (handles multiple types)
    """
    if isa(config, Configuration)
        return config.weight
    elseif isa(config, CompressedConfig)
        return Float64(config.weight)
    elseif isa(config, NamedTuple) && hasfield(typeof(config), :weight)
        return config.weight
    else
        return 0.0
    end
end

function get_config_name(config)
    """
    Get name from configuration (handles multiple types)
    """
    if isa(config, Configuration) || isa(config, CompressedConfig)
        return config.name
    elseif isa(config, NamedTuple) && hasfield(typeof(config), :name)
        return config.name
    else
        return "unknown"
    end
end

function has_nuclear_configs(configs)
    """
    Check if configurations include nuclear excitations
    """
    for config in configs
        name = get_config_name(config)
        if match(r"N\d+", name) !== nothing || occursin("C(", name)
            return true
        end
    end
    return false
end

# ======================== Advanced Analysis ========================

function analyze_importance_distribution(configs)
    """
    Analyze the distribution of importance weights
    """
    if isempty(configs)
        return nothing
    end
    
    weights = [get_weight(c) for c in configs]
    sort!(weights, rev=true)
    
    # Calculate cumulative importance
    cumulative = cumsum(weights) / sum(weights)
    
    # Find how many configs needed for different thresholds
    n_50 = findfirst(x -> x >= 0.5, cumulative)
    n_90 = findfirst(x -> x >= 0.9, cumulative)
    n_95 = findfirst(x -> x >= 0.95, cumulative)
    n_99 = findfirst(x -> x >= 0.99, cumulative)
    
    return (
        total_configs = length(configs),
        n_for_50_percent = something(n_50, length(configs)),
        n_for_90_percent = something(n_90, length(configs)),
        n_for_95_percent = something(n_95, length(configs)),
        n_for_99_percent = something(n_99, length(configs)),
        max_weight = weights[1],
        min_weight = weights[end],
        weight_ratio = weights[1] / max(weights[end], 1e-10)
    )
end

function calculate_effective_rank(configs)
    """
    Calculate effective rank (participation ratio) of configurations
    """
    if isempty(configs)
        return 0.0
    end
    
    weights = [get_weight(c) for c in configs]
    total_weight = sum(weights)
    
    if total_weight <= 0
        return 0.0
    end
    
    # Normalize weights
    probs = weights ./ total_weight
    
    # Calculate Shannon entropy
    entropy = -sum(p * log(p) for p in probs if p > 0)
    
    # Effective rank
    effective_rank = exp(entropy)
    
    return effective_rank
end

# ======================== Importance Filtering ========================

function filter_by_importance(configs, threshold::Float64)
    """
    Filter configurations by importance threshold
    """
    filtered = similar(configs, 0)
    
    for config in configs
        if get_weight(config) >= threshold
            push!(filtered, config)
        end
    end
    
    return filtered
end

function adaptive_importance_threshold(configs, target_configs::Int)
    """
    Find importance threshold to achieve target number of configurations
    """
    if length(configs) <= target_configs
        return 0.0
    end
    
    weights = sort([get_weight(c) for c in configs], rev=true)
    
    return weights[target_configs]
end

# ======================== Visualization Helpers ========================

function importance_summary_string(importance_data)
    """
    Create a summary string of importance analysis
    """
    summary = "Importance Analysis Summary\n"
    summary *= "="^40 * "\n"
    summary *= "Total importance: $(round(importance_data.total_importance, digits=3))\n"
    
    if hasfield(typeof(importance_data), :standard_importance)
        summary *= "Standard importance: $(round(importance_data.standard_importance, digits=3))\n"
    end
    
    if importance_data.neo_metrics !== nothing
        neo = importance_data.neo_metrics
        summary *= "\nNEO Metrics:\n"
        summary *= "  Nuclear participation: $(round(neo.nuclear_participation * 100, digits=1))%\n"
        summary *= "  Coupling contribution: $(round(neo.coupling_contribution * 100, digits=1))%\n"
        summary *= "  Enhancement factor: $(round(neo.enhancement_factor, digits=2))\n"
    end
    
    if hasfield(typeof(importance_data), :config_types)
        summary *= "\nConfiguration Types:\n"
        for (type, count) in importance_data.config_types
            if count > 0
                summary *= "  $type: $count\n"
            end
        end
    end
    
    return summary
end

end # module
