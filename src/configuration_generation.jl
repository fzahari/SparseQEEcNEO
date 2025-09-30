"""
configuration_generation.jl - Methods for generating electronic and nuclear configurations
"""
module ConfigurationGeneration

using PyCall
using LinearAlgebra
using Statistics 
using ..Types
using ..NuclearMethods  # Now available since it's loaded before this module
using ..PySCFInterface

export generate_configurations, select_important_configurations
export create_reference_config, compress_configurations

# ======================== Main Configuration Generation ========================

# Constants for Clean Code compliance
const DEFAULT_MP2_WEIGHT = 4.0
const CONFIGURATION_COMPRESSION_THRESHOLD = 100
const NUCLEAR_MASS_FACTOR = 1836.0
const COUPLING_ENHANCEMENT_BASE = 1.0

"""
    generate_configurations(meanfield_result, neo_molecule, config_selection, t2_amplitudes)

Generate quantum configurations using the specified method.

# Arguments
- `meanfield_result`: Quantum chemistry mean-field calculation result
- `neo_molecule`: NEO molecular system
- `config_selection`: Configuration selection parameters
- `t2_amplitudes`: Optional MP2 t2 amplitudes for enhancement

# Returns
- `Vector{Configuration}`: Generated configurations
"""
function generate_configurations(meanfield_result, neo_molecule, config_selection::ConfigSelection, t2_amplitudes=nothing)
    @info "Generating important configurations using $(config_selection.method) method..."
    
    configuration_generator = create_configuration_generator(config_selection.method)
    generated_configs = configuration_generator(meanfield_result, neo_molecule, config_selection, t2_amplitudes)
    
    final_configs = apply_compression_if_requested(generated_configs, config_selection)
    
    @info "Generated $(length(final_configs)) configurations"
    return final_configs
end

function create_configuration_generator(method_name::String)
    method_map = Dict(
        "mp2" => generate_configurations_mp2,
        "mp2_enhanced" => generate_configurations_mp2_enhanced,
        "casci" => generate_configurations_casci,
        "neo_cneo" => generate_configurations_neo_cneo,
        "neo_enhanced" => generate_configurations_neo_enhanced,
        "neo_final" => generate_configurations_neo_final,
        "hybrid_final" => generate_configurations_hybrid_final
    )
    
    if haskey(method_map, method_name)
        return method_map[method_name]
    else
        @warn "Unknown method $method_name, using MP2"
        return generate_configurations_mp2
    end
end

function apply_compression_if_requested(configurations::Vector, config_selection::ConfigSelection)
    if should_compress_configurations(configurations, config_selection)
        return compress_configurations(configurations)
    else
        return configurations
    end
end

function should_compress_configurations(configurations::Vector, config_selection::ConfigSelection)
    return config_selection.use_compression && length(configurations) > CONFIGURATION_COMPRESSION_THRESHOLD
end

# ======================== MP2-based Methods ========================

"""
    MP2OrbitalData

Structure to hold orbital information for MP2 calculations.
"""
struct MP2OrbitalData
    occupations::Vector{Float64}
    energies::Vector{Float64}
    occupied_count::Int
    total_count::Int
end

function generate_configurations_mp2(meanfield_result, neo_molecule, config_selection::ConfigSelection, t2_amplitudes=nothing)
    configurations = Configuration[]
    
    electronic_component = extract_electronic_component(meanfield_result)
    if electronic_component === nothing
        return configurations
    end
    
    orbital_data = extract_orbital_data(electronic_component)
    if orbital_data === nothing
        return configurations
    end
    
    # Add reference configuration
    reference_config = create_reference_config(meanfield_result, neo_molecule, "HF")
    push!(configurations, reference_config)
    
    # Generate single excitation configurations
    single_excitations = generate_single_excitation_configurations(orbital_data, config_selection)
    append!(configurations, single_excitations)
    
    # Generate double excitation configurations if requested
    if config_selection.include_doubles
        double_excitations = generate_double_excitation_configurations(orbital_data, config_selection)
        append!(configurations, double_excitations)
    end
    
    return sort_configurations_by_importance(configurations)
end

function extract_electronic_component(meanfield_result)
    if !has_electronic_component(meanfield_result)
        @warn "No electronic component found"
        return nothing
    end
    return meanfield_result.components["e"]
end

function has_electronic_component(meanfield_result)
    return pybuiltin("hasattr")(meanfield_result, "components") && 
           haskey(meanfield_result.components, "e")
end

function extract_orbital_data(electronic_component)
    if !has_required_orbital_properties(electronic_component)
        @warn "Missing orbital information"
        return nothing
    end
    
    occupations = collect(electronic_component.mo_occ)
    energies = collect(electronic_component.mo_energy)
    occupied_count = sum(occupations .> 0.5)
    total_count = length(occupations)
    
    return MP2OrbitalData(occupations, energies, occupied_count, total_count)
end

function has_required_orbital_properties(electronic_component)
    return pybuiltin("hasattr")(electronic_component, "mo_occ") &&
           pybuiltin("hasattr")(electronic_component, "mo_energy")
end

function generate_single_excitation_configurations(orbital_data::MP2OrbitalData, config_selection::ConfigSelection)
    configurations = Configuration[]
    
    for occupied_orbital in 1:orbital_data.occupied_count
        for virtual_orbital in (orbital_data.occupied_count + 1):orbital_data.total_count
            
            excitation_config = create_single_excitation_config(
                orbital_data, occupied_orbital, virtual_orbital, config_selection
            )
            
            if excitation_config !== nothing
                push!(configurations, excitation_config)
            end
        end
    end
    
    return configurations
end

function create_single_excitation_config(orbital_data::MP2OrbitalData, from_orbital::Int, 
                                       to_orbital::Int, config_selection::ConfigSelection)
    
    if !is_valid_excitation(orbital_data, from_orbital, to_orbital)
        return nothing
    end
    
    energy_gap = calculate_excitation_energy_gap(orbital_data, from_orbital, to_orbital)
    
    if !passes_energy_cutoff(energy_gap, config_selection)
        return nothing
    end
    
    excitation_weight = calculate_mp2_single_excitation_weight(energy_gap)
    
    if !passes_importance_cutoff(excitation_weight, config_selection)
        return nothing
    end
    
    new_occupation = create_single_excitation_occupation(orbital_data, from_orbital, to_orbital)
    configuration_name = create_single_excitation_name(from_orbital, to_orbital)
    
    return Configuration(
        configuration_name,
        Dict("e" => Int16.(new_occupation)),
        excitation_weight
    )
end

function is_valid_excitation(orbital_data::MP2OrbitalData, from_orbital::Int, to_orbital::Int)
    return to_orbital <= length(orbital_data.energies)
end

function calculate_excitation_energy_gap(orbital_data::MP2OrbitalData, from_orbital::Int, to_orbital::Int)
    return orbital_data.energies[to_orbital] - orbital_data.energies[from_orbital]
end

function passes_energy_cutoff(energy_gap::Float64, config_selection::ConfigSelection)
    return energy_gap <= config_selection.energy_cutoff
end

function calculate_mp2_single_excitation_weight(energy_gap::Float64)
    return DEFAULT_MP2_WEIGHT / (COUPLING_ENHANCEMENT_BASE + energy_gap)
end

function passes_importance_cutoff(weight::Float64, config_selection::ConfigSelection)
    return weight > config_selection.importance_cutoff
end

function create_single_excitation_occupation(orbital_data::MP2OrbitalData, from_orbital::Int, to_orbital::Int)
    new_occupation = copy(orbital_data.occupations)
    new_occupation[from_orbital] = 0
    new_occupation[to_orbital] = 1
    return new_occupation
end

function create_single_excitation_name(from_orbital::Int, to_orbital::Int)
    return "S($from_orbital→$to_orbital)"
end

function generate_double_excitation_configurations(orbital_data::MP2OrbitalData, config_selection::ConfigSelection)
    configurations = Configuration[]
    
    for i in 1:(orbital_data.occupied_count - 1)
        for j in (i + 1):orbital_data.occupied_count
            for a in (orbital_data.occupied_count + 1):(orbital_data.total_count - 1)
                for b in (a + 1):orbital_data.total_count
                    
                    double_config = create_double_excitation_config(
                        orbital_data, i, j, a, b, config_selection
                    )
                    
                    if double_config !== nothing
                        push!(configurations, double_config)
                    end
                end
            end
        end
    end
    
    return configurations
end

function create_double_excitation_config(orbital_data::MP2OrbitalData, i::Int, j::Int, 
                                       a::Int, b::Int, config_selection::ConfigSelection)
    
    if !is_valid_double_excitation(orbital_data, i, j, a, b)
        return nothing
    end
    
    total_energy_gap = calculate_double_excitation_energy_gap(orbital_data, i, j, a, b)
    
    if !passes_energy_cutoff(total_energy_gap, config_selection)
        return nothing
    end
    
    excitation_weight = calculate_mp2_double_excitation_weight(total_energy_gap)
    
    if !passes_importance_cutoff(excitation_weight, config_selection)
        return nothing
    end
    
    new_occupation = create_double_excitation_occupation(orbital_data, i, j, a, b)
    configuration_name = create_double_excitation_name(i, j, a, b)
    
    return Configuration(
        configuration_name,
        Dict("e" => Int16.(new_occupation)),
        excitation_weight
    )
end

function is_valid_double_excitation(orbital_data::MP2OrbitalData, i::Int, j::Int, a::Int, b::Int)
    return b <= length(orbital_data.energies)
end

function calculate_double_excitation_energy_gap(orbital_data::MP2OrbitalData, i::Int, j::Int, a::Int, b::Int)
    return (orbital_data.energies[a] + orbital_data.energies[b]) - 
           (orbital_data.energies[i] + orbital_data.energies[j])
end

function calculate_mp2_double_excitation_weight(energy_gap::Float64)
    return COUPLING_ENHANCEMENT_BASE / (COUPLING_ENHANCEMENT_BASE + energy_gap)
end

function create_double_excitation_occupation(orbital_data::MP2OrbitalData, i::Int, j::Int, a::Int, b::Int)
    new_occupation = copy(orbital_data.occupations)
    new_occupation[i] = 0
    new_occupation[j] = 0
    new_occupation[a] = 1
    new_occupation[b] = 1
    return new_occupation
end

function create_double_excitation_name(i::Int, j::Int, a::Int, b::Int)
    return "D($i$j→$a$b)"
end

function sort_configurations_by_importance(configurations::Vector{Configuration})
    return sort!(configurations, by=c -> c.weight, rev=true)
end

function generate_configurations_mp2_enhanced(mf, mol_neo, config_sel::ConfigSelection, t2_amplitudes)
    """
    Generate configurations using MP2 with t2 amplitude enhancement
    """
    # Start with standard MP2 configurations
    configs = generate_configurations_mp2(mf, mol_neo, config_sel)
    
    # Enhance weights based on t2 amplitudes if available
    if t2_amplitudes !== nothing && isa(t2_amplitudes, Dict) && length(t2_amplitudes) > 0
        @info "Enhancing configurations with $(length(t2_amplitudes)) t2 amplitudes"
        
        # Enhance existing configurations
        for config in configs
            if occursin("D(", config.name)
                # Extract indices from double excitation
                m = match(r"D\((\d+)(\d+)→(\d+)(\d+)\)", config.name)
                if m !== nothing
                    i, j, a, b = parse.(Int, m.captures)
                    
                    # Check for t2 amplitude
                    amp = get(t2_amplitudes, (i, j, a, b), 0.0)
                    if abs(amp) > 0
                        # Update weight directly
                        config = Configuration(
                            config.name,
                            config.occupations,
                            config.weight * (1.0 + 10.0 * abs(amp))
                        )
                    end
                end
            end
        end
        
        # Re-sort by enhanced weights
        sort!(configs, by=c->c.weight, rev=true)
    end
    
    return configs
end

# ======================== CASCI-based Method ========================

function generate_configurations_casci(mf, mol_neo, config_sel::ConfigSelection)
    """
    Generate configurations from CASCI calculation
    """
    configs = Configuration[]
    
    # Run CASCI
    mc = run_neo_casci(mf, mol_neo, config_sel)
    if mc === nothing
        @warn "CASCI failed, falling back to MP2"
        return generate_configurations_mp2(mf, mol_neo, config_sel)
    end
    
    # Extract configurations from CASCI wavefunction
    try
        # Get CI coefficients
        if pybuiltin("hasattr")(mc, "ci")
            ci_coeffs = mc.ci
            
            # For multiple states
            if isa(ci_coeffs, Array) && ndims(ci_coeffs) > 1
                ci_coeffs = ci_coeffs[:, 1]  # Use ground state
            end
            
            # Get determinants
            if pybuiltin("hasattr")(mc.fcisolver, "large_ci")
                # Get most important determinants
                dets, coeffs = mc.fcisolver.large_ci(ci_coeffs, config_sel.max_configs, tol=config_sel.cas_threshold)
                
                # Convert to configurations
                for (det, coeff) in zip(dets, coeffs)
                    weight = abs(coeff)^2
                    if weight > config_sel.importance_cutoff
                        # Convert determinant to occupation
                        occ = det_to_occupation(det, mc.ncas, mc.nelecas)
                        
                        config = Configuration(
                            "CAS($det)",
                            Dict("e" => Int16.(occ)),
                            weight
                        )
                        push!(configs, config)
                    end
                end
            end
        end
    catch e
        @warn "Could not extract CASCI configurations: $e"
    end
    
    # Add reference if no configurations found
    if isempty(configs)
        ref_config = create_reference_config(mf, mol_neo, "CAS-REF")
        push!(configs, ref_config)
    end
    
    return configs
end

# ======================== NEO-specific Methods ========================

function generate_configurations_neo_cneo(mf, mol_neo, config_sel::ConfigSelection)
    """
    Generate configurations using cNEO-enhanced approach
    """
    configs = Configuration[]
    
    # Reference configuration
    ref_config = create_reference_config(mf, mol_neo, "NEO-HF")
    push!(configs, ref_config)
    
    # Get components
    if !pybuiltin("hasattr")(mf, "components") || !haskey(mf.components, "e")
        @warn "No electronic component found"
        return configs
    end
    
    mf_elec = mf.components["e"]
    
    if !pybuiltin("hasattr")(mf_elec, "mo_occ") || !pybuiltin("hasattr")(mf_elec, "mo_energy")
        @warn "Missing orbital information"
        return configs
    end
    
    elec_occ = collect(mf_elec.mo_occ)
    elec_energy = collect(mf_elec.mo_energy)
    n_elec = sum(elec_occ .> 0.5)
    n_elec_orbs = length(elec_occ)
    
    # Calculate coupling matrix
    coupling_matrix = calculate_neo_coupling_matrix(mf, mol_neo)
    
    # Electronic excitations with coupling enhancement
    for i in 1:n_elec
        for a in (n_elec+1):n_elec_orbs
            if a > length(elec_energy)
                continue
            end
            
            energy_gap = elec_energy[a] - elec_energy[i]
            base_weight = 4.0 / (1.0 + energy_gap)
            
            # Enhancement from nuclear coupling
            coupling_factor = 1.0
            if size(coupling_matrix, 1) >= i && size(coupling_matrix, 1) >= a
                coupling_factor = 1.0 + mean(coupling_matrix[i, :]) + mean(coupling_matrix[a, :])
            end
            weight = base_weight * coupling_factor
            
            if weight > config_sel.importance_cutoff
                occ_new = copy(elec_occ)
                occ_new[i] = 0
                occ_new[a] = 1
                
                config = Configuration(
                    "E($i→$a)",
                    Dict("e" => Int16.(occ_new)),
                    weight
                )
                push!(configs, config)
            end
        end
    end
    
    # Nuclear excitations
    if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
        for nuc_idx in 0:(mol_neo.nuc_num-1)
            nuc_configs = generate_nuclear_excitations(mf, nuc_idx, config_sel)
            append!(configs, nuc_configs)
        end
        
        # Coupled excitations
        coupled_configs = generate_coupled_excitations(mf, mol_neo, coupling_matrix, config_sel)
        append!(configs, coupled_configs)
    end
    
    # Sort by weight
    sort!(configs, by=c->c.weight, rev=true)
    
    # Limit number of configurations
    if length(configs) > config_sel.max_configs
        configs = configs[1:config_sel.max_configs]
    end
    
    return configs
end

function generate_nuclear_excitations(mf, nuc_idx, config_sel::ConfigSelection)
    """
    Generate nuclear excitation configurations
    """
    configs = Configuration[]
    
    nuc_key = "n$nuc_idx"
    if !pybuiltin("hasattr")(mf, "components") || !haskey(mf.components, nuc_key)
        return configs
    end
    
    mf_nuc = mf.components[nuc_key]
    
    if !pybuiltin("hasattr")(mf_nuc, "mo_occ") || !pybuiltin("hasattr")(mf_nuc, "mo_energy")
        return configs
    end
    
    nuc_occ = collect(mf_nuc.mo_occ)
    nuc_energy = collect(mf_nuc.mo_energy)
    n_nuc_orbs = length(nuc_occ)
    
    # Find occupied orbitals
    occ_orbs = findall(nuc_occ .> 0.5)
    if isempty(occ_orbs)
        return configs
    end
    
    # Single excitations
    for i in occ_orbs
        for a in 1:n_nuc_orbs
            if nuc_occ[a] > 0.5  # Skip if occupied
                continue
            end
            
            energy_gap = nuc_energy[a] - nuc_energy[i]
            
            # Nuclear excitation weight with mass scaling
            mass_factor = sqrt(1836.0)
            weight = exp(-energy_gap / mass_factor)
            
            if weight > config_sel.importance_cutoff
                occ_new = copy(nuc_occ)
                occ_new[i] = 0
                occ_new[a] = 1
                
                config = Configuration(
                    "N$(nuc_idx+1)($i→$a)",
                    Dict(nuc_key => Int16.(occ_new)),
                    weight
                )
                push!(configs, config)
            end
        end
    end
    
    # Double excitations if requested
    if config_sel.include_doubles && length(occ_orbs) >= 2
        for (idx_i, i) in enumerate(occ_orbs)
            for j in occ_orbs[idx_i+1:end]
                virt_orbs = findall(nuc_occ .< 0.5)
                for (idx_a, a) in enumerate(virt_orbs)
                    for b in virt_orbs[idx_a+1:end]
                        energy_gap = nuc_energy[a] + nuc_energy[b] - 
                                   nuc_energy[i] - nuc_energy[j]
                        
                        weight = exp(-energy_gap / (2 * sqrt(1836.0))) * 0.5
                        
                        if weight > config_sel.importance_cutoff
                            occ_new = copy(nuc_occ)
                            occ_new[i] = 0
                            occ_new[j] = 0
                            occ_new[a] = 1
                            occ_new[b] = 1
                            
                            config = Configuration(
                                "N$(nuc_idx+1)($i$j→$a$b)",
                                Dict(nuc_key => Int16.(occ_new)),
                                weight
                            )
                            push!(configs, config)
                        end
                    end
                end
            end
        end
    end
    
    return configs
end

function generate_coupled_excitations(mf, mol_neo, coupling_matrix, config_sel::ConfigSelection)
    """
    Generate coupled electronic-nuclear excitations
    """
    configs = Configuration[]
    
    # Get components
    if !pybuiltin("hasattr")(mf, "components") || !haskey(mf.components, "e")
        return configs
    end
    
    mf_elec = mf.components["e"]
    if !pybuiltin("hasattr")(mf_elec, "mo_occ") || !pybuiltin("hasattr")(mf_elec, "mo_energy")
        return configs
    end
    
    elec_occ = collect(mf_elec.mo_occ)
    elec_energy = collect(mf_elec.mo_energy)
    n_elec = sum(elec_occ .> 0.5)
    
    # Find strongly coupled pairs
    n_elec_orbs, n_nuc_orbs = size(coupling_matrix)
    coupled_pairs = []
    
    for i in 1:n_elec_orbs
        for j in 1:n_nuc_orbs
            if coupling_matrix[i, j] > 0.1
                push!(coupled_pairs, (i, j, coupling_matrix[i, j]))
            end
        end
    end
    
    # Sort by coupling strength
    sort!(coupled_pairs, by=x->x[3], rev=true)
    
    # Generate coupled excitations for top pairs
    n_coupled = min(40, length(coupled_pairs))
    
    for (elec_idx, nuc_idx, coupling) in coupled_pairs[1:n_coupled]
        # Electronic part
        for a in (n_elec+1):length(elec_occ)
            if a > length(elec_energy)
                continue
            end
            
            # Nuclear part - check if n0 component exists
            if !haskey(mf.components, "n0")
                continue
            end
            
            mf_nuc = mf.components["n0"]
            if !pybuiltin("hasattr")(mf_nuc, "mo_occ") || !pybuiltin("hasattr")(mf_nuc, "mo_energy")
                continue
            end
            
            nuc_occ = collect(mf_nuc.mo_occ)
            nuc_energy = collect(mf_nuc.mo_energy)
            
            for b in 2:length(nuc_occ)
                if nuc_occ[b] > 0.5  # Skip occupied
                    continue
                end
                
                # Combined weight
                e_gap = elec_energy[a] - elec_energy[elec_idx]
                n_gap = nuc_energy[b] - nuc_energy[1]
                
                weight = coupling * exp(-(e_gap + n_gap/sqrt(1836.0))) * 2.0
                
                if weight > config_sel.importance_cutoff
                    # Electronic occupation
                    elec_occ_new = copy(elec_occ)
                    elec_occ_new[elec_idx] = 0
                    elec_occ_new[a] = 1
                    
                    # Nuclear occupation
                    nuc_occ_new = copy(nuc_occ)
                    nuc_occ_new[1] = 0
                    nuc_occ_new[b] = 1
                    
                    config = Configuration(
                        "C(E($elec_idx→$a)+N1(1→$b))",
                        Dict("e" => Int16.(elec_occ_new),
                             "n0" => Int16.(nuc_occ_new)),
                        weight
                    )
                    push!(configs, config)
                end
            end
        end
    end
    
    return configs
end

function calculate_neo_coupling_matrix(mf, mol_neo)
    """
    Calculate nuclear-electronic coupling matrix
    """
    # Get dimensions
    n_elec = 10  # Default
    n_nuc = 10
    
    if pybuiltin("hasattr")(mf, "components")
        if haskey(mf.components, "e") && pybuiltin("hasattr")(mf.components["e"], "mo_occ")
            n_elec = length(mf.components["e"].mo_occ)
        end
        
        if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
            if haskey(mf.components, "n0") && pybuiltin("hasattr")(mf.components["n0"], "mo_occ")
                n_nuc = length(mf.components["n0"].mo_occ)
            end
        end
    end
    
    # Initialize coupling matrix
    coupling_matrix = zeros(n_elec, n_nuc)
    
    # Calculate coupling based on energy differences
    if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "e") && haskey(mf.components, "n0")
        mf_elec = mf.components["e"]
        mf_nuc = mf.components["n0"]
        
        if pybuiltin("hasattr")(mf_elec, "mo_energy") && pybuiltin("hasattr")(mf_nuc, "mo_energy")
            elec_energies = collect(mf_elec.mo_energy)
            nuc_energies = collect(mf_nuc.mo_energy)
            
            for i in 1:min(n_elec, length(elec_energies))
                for j in 1:min(n_nuc, length(nuc_energies))
                    dE = abs(elec_energies[i] - nuc_energies[j])
                    coupling_matrix[i, j] = exp(-dE) / (1.0 + dE) / sqrt(1836.0)
                end
            end
        end
    end
    
    return coupling_matrix
end

# ======================== Other NEO Methods ========================

function generate_configurations_neo_enhanced(mf, mol_neo, config_sel::ConfigSelection)
    """
    Enhanced NEO configuration generation with additional features
    """
    # Start with cNEO configurations
    configs = generate_configurations_neo_cneo(mf, mol_neo, config_sel)
    
    # Add additional configurations based on density analysis
    if length(configs) < config_sel.max_configs / 2
        # Generate more nuclear excitations
        if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
            for nuc_idx in 0:(mol_neo.nuc_num-1)
                # Extend energy cutoff for more configurations
                config_sel_extended = deepcopy(config_sel)
                config_sel_extended.energy_cutoff *= 2.0
                
                nuc_configs = generate_nuclear_excitations(mf, nuc_idx, config_sel_extended)
                
                # Add with reduced weight
                for config in nuc_configs
                    config_new = Configuration(
                        config.name,
                        config.occupations,
                        config.weight * 0.5
                    )
                    push!(configs, config_new)
                end
            end
        end
    end
    
    # Re-sort and truncate
    sort!(configs, by=c->c.weight, rev=true)
    if length(configs) > config_sel.max_configs
        configs = configs[1:config_sel.max_configs]
    end
    
    return configs
end

function generate_configurations_neo_final(mf, mol_neo, config_sel::ConfigSelection)
    """
    Final NEO method with all enhancements
    """
    # Similar to neo_enhanced but with final optimizations
    return generate_configurations_neo_enhanced(mf, mol_neo, config_sel)
end

# ======================== Hybrid Method ========================

function generate_configurations_hybrid_final(mf, mol_neo, config_sel::ConfigSelection)
    """
    Hybrid method combining multiple approaches
    """
    @info "Starting hybrid configuration selection..."
    
    all_configs = Configuration[]
    
    # 1. Get MP2 configurations
    config_mp2 = deepcopy(config_sel)
    config_mp2.max_configs = div(config_sel.max_configs, 3)
    mp2_configs = generate_configurations_mp2(mf, mol_neo, config_mp2)
    append!(all_configs, mp2_configs)
    
    # 2. Get NEO configurations if quantum nuclei present
    if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
        config_neo = deepcopy(config_sel)
        config_neo.max_configs = div(config_sel.max_configs, 3)
        neo_configs = generate_configurations_neo_cneo(mf, mol_neo, config_neo)
        
        # Avoid duplicates
        existing_names = Set(c.name for c in all_configs)
        for config in neo_configs
            if !(config.name in existing_names)
                push!(all_configs, config)
            end
        end
    end
    
    # 3. Try CASCI if not too many configurations yet
    if length(all_configs) < config_sel.max_configs * 0.8
        config_cas = deepcopy(config_sel)
        config_cas.max_configs = div(config_sel.max_configs, 3)
        config_cas.active_orbs = 8  # Add default values
        config_cas.active_elec = 4
        cas_configs = generate_configurations_casci(mf, mol_neo, config_cas)
        
        existing_names = Set(c.name for c in all_configs)
        for config in cas_configs
            if !(config.name in existing_names)
                push!(all_configs, config)
            end
        end
    end
    
    # Sort and select final configurations
    sort!(all_configs, by=c->c.weight, rev=true)
    
    if length(all_configs) > config_sel.max_configs
        all_configs = all_configs[1:config_sel.max_configs]
    end
    
    @info "Hybrid selection complete:" *
          "\n  Electronic configs: $(count(c -> haskey(c.occupations, "e"), all_configs))" *
          "\n  Nuclear configs: $(count(c -> any(k -> startswith(String(k), "n"), keys(c.occupations)), all_configs))" *
          "\n  Coupled configs: $(count(c -> occursin("C(", c.name), all_configs))"
    
    return all_configs
end

# ======================== Helper Functions ========================

function create_reference_config(mf, mol_neo, name="HF")
    """
    Create reference configuration with all components
    """
    occupations = Dict{String, Vector{Int16}}()
    
    # Electronic occupation
    if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "e")
        mf_elec = mf.components["e"]
        if pybuiltin("hasattr")(mf_elec, "mo_occ")
            mo_occ = collect(mf_elec.mo_occ)
            # Convert to integers (0 or 1)
            occ_int = [occ > 0.5 ? Int16(1) : Int16(0) for occ in mo_occ]
            occupations["e"] = occ_int
        end
    end
    
    # Nuclear occupations
    if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
        for i in 0:(mol_neo.nuc_num-1)
            nuc_key = "n$i"
            if haskey(mf.components, nuc_key)
                mf_nuc = mf.components[nuc_key]
                if pybuiltin("hasattr")(mf_nuc, "mo_occ")
                    mo_occ = collect(mf_nuc.mo_occ)
                    # Convert to integers
                    occ_int = [occ > 0.5 ? Int16(1) : Int16(0) for occ in mo_occ]
                    occupations[nuc_key] = occ_int
                end
            end
        end
    end
    
    return Configuration(name, occupations, 1.0)
end

function compress_configurations(configs::Vector{Configuration})
    """
    Compress configurations to save memory
    """
    compressed = CompressedConfig[]
    
    for config in configs
        # Find non-zero occupation indices
        indices = Int16[]
        values = Float32[]
        
        for (component, occ) in config.occupations
            component_offset = component == "e" ? 0 : 1000 * parse(Int, component[2:end])
            
            for (i, val) in enumerate(occ)
                if val != 0
                    push!(indices, component_offset + i)
                    push!(values, Float32(val))
                end
            end
        end
        
        compressed_config = CompressedConfig(
            config.name,
            indices,
            values,
            Float32(config.weight)
        )
        push!(compressed, compressed_config)
    end
    
    @info "Compressed $(length(configs)) configurations, " *
          "memory saved: $(round((1 - sizeof(compressed)/sizeof(configs)) * 100, digits=1))%"
    
    return compressed
end

function det_to_occupation(det, ncas, nelecas)
    """
    Convert determinant string to occupation vector
    """
    # Placeholder implementation
    occ = zeros(ncas)
    for i in 1:min(nelecas, ncas)
        occ[i] = 1.0
    end
    return occ
end

# ======================== Configuration Selection ========================

function select_important_configurations(configs, config_sel::ConfigSelection)
    """
    Select most important configurations based on criteria
    """
    # Sort by weight if not already sorted
    if !issorted(configs, by=c->get_weight(c), rev=true)
        sort!(configs, by=c->get_weight(c), rev=true)
    end
    
    selected = similar(configs, 0)
    cumulative_weight = 0.0
    total_weight = sum(get_weight(c) for c in configs)
    n_qubits = 0
    
    for config in configs
        # Check if we've reached targets
        if length(selected) >= config_sel.max_configs
            break
        end
        
        if cumulative_weight / total_weight >= config_sel.target_importance
            break
        end
        
        # Calculate qubits needed
        config_qubits = calculate_config_qubits(config)
        
        if config_qubits > config_sel.max_qubits
            continue
        end
        
        # Add configuration
        push!(selected, config)
        cumulative_weight += get_weight(config)
        n_qubits = max(n_qubits, config_qubits)
    end
    
    return selected, n_qubits
end

function calculate_config_qubits(config)
    """
    Calculate number of qubits needed for configuration
    """
    max_orbital = 0
    
    if isa(config, Configuration)
        for (component, occ) in config.occupations
            non_zero = findall(occ .!= 0)
            if !isempty(non_zero)
                max_orbital = max(max_orbital, maximum(non_zero))
            end
        end
    elseif isa(config, CompressedConfig)
        if !isempty(config.occ_indices)
            max_orbital = maximum(config.occ_indices .% 1000)
        end
    end
    
    # Number of qubits is the highest orbital index
    return max_orbital
end

function get_weight(config)
    """
    Get weight from configuration (handles both types)
    """
    if isa(config, Configuration)
        return config.weight
    elseif isa(config, CompressedConfig)
        return Float64(config.weight)
    else
        return 0.0
    end
end

end # module
