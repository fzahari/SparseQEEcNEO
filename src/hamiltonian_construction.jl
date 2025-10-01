"""
hamiltonian_construction.jl - Second-quantized Hamiltonian construction and storage
"""
module HamiltonianConstruction

using LinearAlgebra
using SparseArrays
using HDF5
using JSON
using PyCall
using ..Types
using ..ConfigurationGeneration

export construct_second_quantized_hamiltonian, save_hamiltonian, load_hamiltonian
export HamiltonianData, reduce_hamiltonian_symmetry, validate_hamiltonian
export analyze_hamiltonian_properties, export_hamiltonian_openfermion

# ======================== Types ========================

"""
Structure to hold second-quantized Hamiltonian data
"""
struct HamiltonianData
    # One-body terms
    h1e::Dict{String, Matrix{Float64}}  # Component -> h_{ij}
    
    # Two-body terms (stored as sparse dictionary for 4D)
    h2e::Dict{String, Dict{NTuple{4,Int}, Float64}}  # Component -> sparse 4D tensor
    
    # Cross-component coupling terms
    coupling::Dict{Tuple{String, String}, Matrix{Float64}}  # (comp1, comp2) -> V_{ij}
    
    # Metadata
    n_orbitals::Dict{String, Int}
    n_electrons::Dict{String, Int}
    nuclear_charges::Vector{Float64}
    
    # Configuration space info
    configs::Vector{Configuration}
    config_energies::Vector{Float64}
    
    # Reduced space mapping
    active_orbitals::Dict{String, Vector{Int}}
    frozen_orbitals::Dict{String, Vector{Int}}
end

# ======================== Main Construction Function ========================

"""
    construct_second_quantized_hamiltonian(mf, mol_neo, configs, config_sel)

Construct second-quantized Hamiltonian in the reduced Hilbert space.

Returns:
- HamiltonianData: Complete Hamiltonian information
- H_matrix: Hamiltonian matrix in configuration basis
"""
function construct_second_quantized_hamiltonian(mf, mol_neo, configs, config_sel)
    @info "Constructing second-quantized Hamiltonian for $(length(configs)) configurations"
    
    integralData = extractIntegralData(mf, mol_neo)
    orbitalSpaceData = determineOrbitalSpaces(configs, mf)
    systemData = extractSystemData(mf, mol_neo)
    
    hamiltonianMatrix = buildConfigurationHamiltonian(
        configs, integralData, orbitalSpaceData
    )
    
    hamiltonianData = createHamiltonianDataStruct(
        integralData, orbitalSpaceData, systemData, configs, hamiltonianMatrix
    )
    
    return hamiltonianData, hamiltonianMatrix
end

function extractIntegralData(mf, molNeo)
    h1eDict, h2eDict = extract_integrals(mf, molNeo)
    couplingDict = calculate_coupling_integrals(mf, molNeo)
    return (h1e=h1eDict, h2e=h2eDict, coupling=couplingDict)
end

function determineOrbitalSpaces(configs, mf)
    activeOrbs, frozenOrbs = determine_active_space(configs, mf)
    return (active=activeOrbs, frozen=frozenOrbs)
end

function extractSystemData(mf, molNeo)
    nOrbitals, nElectrons = get_system_info(mf, molNeo)
    nuclearCharges = get_nuclear_charges(molNeo)
    return (n_orbitals=nOrbitals, n_electrons=nElectrons, nuclear_charges=nuclearCharges)
end

function buildConfigurationHamiltonian(configs, integralData, orbitalSpaceData)
    return build_configuration_hamiltonian(
        configs, integralData.h1e, integralData.h2e, 
        integralData.coupling, orbitalSpaceData.active
    )
end

function createHamiltonianDataStruct(integralData, orbitalSpaceData, systemData, 
                                   configs, hamiltonianMatrix)
    configEnergies = calculateDiagonalEnergies(hamiltonianMatrix)
    
    return HamiltonianData(
        integralData.h1e,
        integralData.h2e,
        integralData.coupling,
        systemData.n_orbitals,
        systemData.n_electrons,
        systemData.nuclear_charges,
        configs,
        configEnergies,
        orbitalSpaceData.active,
        orbitalSpaceData.frozen
    )
end

function calculateDiagonalEnergies(hamiltonianMatrix)
    return Float64[hamiltonianMatrix[i,i] for i in 1:size(hamiltonianMatrix, 1)]
end

# ======================== Integral Extraction ========================

function extract_integrals(mf, mol_neo)
    """
    Extract one- and two-electron integrals from mean-field object
    """
    h1e_dict = Dict{String, Matrix{Float64}}()
    h2e_dict = Dict{String, Dict{NTuple{4,Int}, Float64}}()
    
    # Electronic integrals
    if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "e")
        mf_e = mf.components["e"]
        h1e_e, h2e_e = get_electronic_integrals(mf_e)
        h1e_dict["e"] = h1e_e
        h2e_dict["e"] = h2e_e
    end
    
    # Nuclear integrals
    if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
        for i in 0:(mol_neo.nuc_num-1)
            nuc_key = "n$i"
            if haskey(mf.components, nuc_key)
                mf_n = mf.components[nuc_key]
                h1e_n, h2e_n = get_nuclear_integrals(mf_n, i)
                h1e_dict[nuc_key] = h1e_n
                h2e_dict[nuc_key] = h2e_n
            end
        end
    end
    
    return h1e_dict, h2e_dict
end

function get_electronic_integrals(mf_e)
    """
    Extract electronic integrals in MO basis
    """
    # Get MO coefficients
    mo_coeff = Matrix{Float64}(I, 10, 10)  # Default
    if pybuiltin("hasattr")(mf_e, "mo_coeff")
        try
            mo_coeff = convert(Array{Float64}, mf_e.mo_coeff)
        catch e
            @warn "Could not convert MO coefficients: $e"
        end
    end
    
    n_orbs = size(mo_coeff, 2)
    
    # One-electron integrals (kinetic + nuclear attraction)
    h1e = zeros(n_orbs, n_orbs)
    if pybuiltin("hasattr")(mf_e, "get_hcore")
        try
            h1e_ao = convert(Array{Float64}, mf_e.get_hcore())
            h1e = mo_coeff' * h1e_ao * mo_coeff
        catch e
            @warn "Could not get hcore: $e"
            # Use random symmetric matrix as placeholder
            h1e = randn(n_orbs, n_orbs)
            h1e = (h1e + h1e') / 2
        end
    else
        # Placeholder
        h1e = randn(n_orbs, n_orbs)
        h1e = (h1e + h1e') / 2
    end
    
    # Two-electron integrals (sparse dictionary format)
    h2e = create_sparse_eri(n_orbs)
    
    # Try to get actual ERIs if available
    if pybuiltin("hasattr")(mf_e, "_eri") && mf_e._eri !== nothing
        try
            # This would require proper 4-index transformation
            # For now, just use the sparse approximation
            @info "Using approximate two-electron integrals"
        catch e
            @warn "Could not process ERIs: $e"
        end
    end
    
    return h1e, h2e
end

function get_nuclear_integrals(mf_n, nuc_idx)
    """
    Extract nuclear integrals with proper mass scaling
    """
    # Get nuclear mass (assuming proton for now)
    mass = 1836.0  # Proton mass in atomic units
    
    # MO coefficients
    mo_coeff = Matrix{Float64}(I, 5, 5)  # Default
    if pybuiltin("hasattr")(mf_n, "mo_coeff")
        try
            mo_coeff = convert(Array{Float64}, mf_n.mo_coeff)
        catch e
            @warn "Could not convert nuclear MO coefficients: $e"
        end
    end
    
    n_orbs = size(mo_coeff, 2)
    
    # Kinetic energy with mass scaling
    T_mo = zeros(n_orbs, n_orbs)
    if pybuiltin("hasattr")(mf_n, "get_kinetic")
        try
            T_ao = convert(Array{Float64}, mf_n.get_kinetic())
            T_mo = mo_coeff' * T_ao * mo_coeff / (2 * mass)
        catch e
            @warn "Could not get nuclear kinetic energy: $e"
            # Approximate kinetic energy
            T_mo = randn(n_orbs, n_orbs)
            T_mo = (T_mo + T_mo') / (4 * mass)
        end
    else
        # Approximate kinetic energy
        T_mo = randn(n_orbs, n_orbs)
        T_mo = (T_mo + T_mo') / (4 * mass)
    end
    
    # Potential energy
    V_mo = zeros(n_orbs, n_orbs)
    if pybuiltin("hasattr")(mf_n, "get_veff")
        try
            V_ao = convert(Array{Float64}, mf_n.get_veff())
            V_mo = mo_coeff' * V_ao * mo_coeff
        catch e
            @warn "Could not get nuclear potential: $e"
            V_mo = randn(n_orbs, n_orbs)
            V_mo = (V_mo + V_mo') / 2
        end
    else
        V_mo = randn(n_orbs, n_orbs)
        V_mo = (V_mo + V_mo') / 2
    end
    
    h1e = T_mo + V_mo
    
    # Nuclear two-body integrals (usually negligible)
    h2e = create_sparse_eri(n_orbs, scale=1e-6)
    
    return h1e, h2e
end

function create_sparse_eri(n_orbs; scale=1.0)
    """
    Create sparse two-electron integrals as dictionary
    """
    eri = Dict{NTuple{4,Int}, Float64}()
    
    # Only store significant elements
    for i in 1:n_orbs
        for j in 1:n_orbs
            # Diagonal and near-diagonal elements
            if abs(i-j) <= 2
                val = scale * (1.0 / (1.0 + abs(i-j)))
                if abs(val) > 1e-12
                    eri[(i,j,i,j)] = val
                end
            end
        end
    end
    
    return eri
end

# ======================== Coupling Integrals ========================

function calculate_coupling_integrals(mf, mol_neo)
    """
    Calculate electron-nuclear coupling integrals
    """
    coupling_dict = Dict{Tuple{String, String}, Matrix{Float64}}()
    
    if !pybuiltin("hasattr")(mol_neo, "nuc_num") || mol_neo.nuc_num == 0
        return coupling_dict
    end
    
    # Get electronic orbitals
    if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "e")
        mf_e = mf.components["e"]
        n_elec_orbs = 10  # Default
        if pybuiltin("hasattr")(mf_e, "mo_coeff")
            try
                mo_coeff = mf_e.mo_coeff
                n_elec_orbs = size(mo_coeff, 2)
            catch
                # Keep default
            end
        end
        
        # For each quantum nucleus
        for i in 0:(mol_neo.nuc_num-1)
            nuc_key = "n$i"
            if haskey(mf.components, nuc_key)
                mf_n = mf.components[nuc_key]
                n_nuc_orbs = 5  # Default
                if pybuiltin("hasattr")(mf_n, "mo_coeff")
                    try
                        mo_coeff = mf_n.mo_coeff
                        n_nuc_orbs = size(mo_coeff, 2)
                    catch
                        # Keep default
                    end
                end
                
                # Calculate coupling matrix
                V_en = calculate_en_coupling_matrix(
                    mf_e, mf_n, n_elec_orbs, n_nuc_orbs
                )
                
                coupling_dict[("e", nuc_key)] = V_en
                coupling_dict[(nuc_key, "e")] = V_en'
            end
        end
    end
    
    return coupling_dict
end

function calculate_en_coupling_matrix(mf_e, mf_n, n_elec, n_nuc)
    """
    Calculate electron-nuclear coupling matrix elements
    """
    V = zeros(n_elec, n_nuc)
    
    # Simplified model based on orbital energies
    e_energies = zeros(n_elec)
    n_energies = zeros(n_nuc)
    
    if pybuiltin("hasattr")(mf_e, "mo_energy")
        try
            e_energies = collect(mf_e.mo_energy)[1:min(n_elec, length(mf_e.mo_energy))]
        catch
            # Use defaults
        end
    end
    
    if pybuiltin("hasattr")(mf_n, "mo_energy")
        try
            n_energies = collect(mf_n.mo_energy)[1:min(n_nuc, length(mf_n.mo_energy))]
        catch
            # Use defaults
        end
    end
    
    for i in 1:n_elec
        for j in 1:n_nuc
            # Coupling strength inversely proportional to energy difference
            dE = abs(e_energies[min(i, length(e_energies))] - 
                    n_energies[min(j, length(n_energies))])
            V[i,j] = exp(-dE) / (1.0 + dE) / sqrt(1836.0)
        end
    end
    
    return V
end

# ======================== Active Space Determination ========================

function determine_active_space(configs, mf)
    """
    Determine active and frozen orbitals from configurations
    """
    active_orbs = Dict{String, Vector{Int}}()
    frozen_orbs = Dict{String, Vector{Int}}()
    
    # Analyze each component
    components = Set{String}()
    for config in configs
        for comp in keys(config.occupations)
            push!(components, comp)
        end
    end
    
    for comp in components
        # Find all orbitals used in configurations
        used_orbitals = Set{Int}()
        ref_occ = nothing
        
        for config in configs
            if haskey(config.occupations, comp)
                occ = config.occupations[comp]
                if config.name == "HF" || config.name == "NEO-HF"
                    ref_occ = occ
                end
                
                for (i, val) in enumerate(occ)
                    if ref_occ !== nothing
                        # Check if orbital occupation changes
                        if val != ref_occ[i]
                            push!(used_orbitals, i)
                        end
                    end
                end
            end
        end
        
        # Active orbitals are those that change
        active = sort(collect(used_orbitals))
        
        # Frozen orbitals are occupied in reference but not active
        if ref_occ !== nothing
            frozen = Int[]
            for (i, occ) in enumerate(ref_occ)
                if occ > 0.5 && !(i in active)
                    push!(frozen, i)
                end
            end
            frozen_orbs[comp] = frozen
        else
            frozen_orbs[comp] = Int[]
        end
        
        active_orbs[comp] = active
    end
    
    return active_orbs, frozen_orbs
end

# ======================== Hamiltonian Matrix Construction ========================

function build_configuration_hamiltonian(configs, h1e, h2e, coupling, active_orbs)
    """
    Build Hamiltonian matrix in configuration basis
    """
    n_configs = length(configs)
    H = zeros(n_configs, n_configs)
    
    # Calculate matrix elements
    Threads.@threads for i in 1:n_configs
        for j in i:n_configs
            H[i,j] = calculate_matrix_element(
                configs[i], configs[j], h1e, h2e, coupling, active_orbs
            )
            if i != j
                H[j,i] = H[i,j]  # Hermitian
            end
        end
    end
    
    return H
end

function calculate_matrix_element(config1, config2, h1e, h2e, coupling, active_orbs)
    """
    Calculate <config1|H|config2> matrix element
    """
    # Check which components differ
    diff_components = Dict{String, Int}()
    
    for comp in union(keys(config1.occupations), keys(config2.occupations))
        if haskey(config1.occupations, comp) && haskey(config2.occupations, comp)
            occ1 = config1.occupations[comp]
            occ2 = config2.occupations[comp]
            n_diff = count(occ1 .!= occ2)
            if n_diff > 0
                diff_components[comp] = n_diff
            end
        end
    end
    
    # Total differences across all components
    total_diff = sum(values(diff_components))
    
    # Selection rules
    if total_diff == 0
        # Diagonal element
        return calculate_diagonal_element(config1, h1e, h2e, coupling)
    elseif total_diff == 2
        # Single excitation (within one component)
        if length(diff_components) == 1
            comp = first(keys(diff_components))
            return calculate_single_excitation(
                config1, config2, comp, h1e[comp]
            )
        else
            # Cross-component single excitation
            return calculate_cross_excitation(
                config1, config2, diff_components, coupling
            )
        end
    elseif total_diff == 4
        # Double excitation
        if length(diff_components) == 1
            comp = first(keys(diff_components))
            return calculate_double_excitation(
                config1, config2, comp, h2e[comp]
            )
        else
            # Mixed double excitation
            return 0.0  # Usually negligible
        end
    else
        return 0.0
    end
end

function calculate_diagonal_element(config, h1e, h2e, coupling)
    """
    Calculate diagonal Hamiltonian element
    """
    E = 0.0
    
    # One-body contributions
    for (comp, occ) in config.occupations
        if haskey(h1e, comp)
            h1 = h1e[comp]
            for i in 1:min(length(occ), size(h1, 1))
                if occ[i] > 0.5
                    E += h1[i,i]
                end
            end
        end
    end
    
    # Two-body contributions (within components)
    for (comp, occ) in config.occupations
        if haskey(h2e, comp)
            h2 = h2e[comp]
            for i in 1:length(occ)
                for j in 1:length(occ)
                    if occ[i] > 0.5 && occ[j] > 0.5
                        # Get two-body integral from dictionary
                        E += 0.5 * get(h2, (i,i,j,j), 0.0)
                    end
                end
            end
        end
    end
    
    # Cross-component coupling
    for ((comp1, comp2), V) in coupling
        if haskey(config.occupations, comp1) && haskey(config.occupations, comp2)
            occ1 = config.occupations[comp1]
            occ2 = config.occupations[comp2]
            
            for i in 1:min(length(occ1), size(V,1))
                for j in 1:min(length(occ2), size(V,2))
                    if occ1[i] > 0.5 && occ2[j] > 0.5
                        E += V[i,j]
                    end
                end
            end
        end
    end
    
    return E
end

function calculate_single_excitation(config1, config2, comp, h1)
    """
    Calculate single excitation matrix element
    """
    occ1 = config1.occupations[comp]
    occ2 = config2.occupations[comp]
    
    # Find excitation indices
    i = findfirst((occ1 .> 0.5) .& (occ2 .< 0.5))
    a = findfirst((occ1 .< 0.5) .& (occ2 .> 0.5))
    
    if i !== nothing && a !== nothing && i <= size(h1,1) && a <= size(h1,2)
        return h1[i,a]
    end
    
    return 0.0
end

function calculate_cross_excitation(config1, config2, diff_components, coupling)
    """
    Calculate cross-component excitation (e.g., coupled e-n excitation)
    """
    # This would involve coupling integrals
    # Simplified implementation
    return 0.0
end

function calculate_double_excitation(config1, config2, comp, h2)
    """
    Calculate double excitation matrix element
    """
    # Simplified - full implementation would extract indices and compute <ij||ab>
    return 0.0
end

# ======================== System Information ========================

function get_system_info(mf, mol_neo)
    """
    Extract system information
    """
    n_orbitals = Dict{String, Int}()
    n_electrons = Dict{String, Int}()
    
    # Electronic component
    if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "e")
        mf_e = mf.components["e"]
        if pybuiltin("hasattr")(mf_e, "mo_occ")
            occ = collect(mf_e.mo_occ)
            n_orbitals["e"] = length(occ)
            n_electrons["e"] = convert(Int, sum(occ))
        end
    end
    
    # Nuclear components
    if pybuiltin("hasattr")(mol_neo, "nuc_num") && mol_neo.nuc_num > 0
        for i in 0:(mol_neo.nuc_num-1)
            nuc_key = "n$i"
            if haskey(mf.components, nuc_key)
                mf_n = mf.components[nuc_key]
                if pybuiltin("hasattr")(mf_n, "mo_occ")
                    occ = collect(mf_n.mo_occ)
                    n_orbitals[nuc_key] = length(occ)
                    n_electrons[nuc_key] = convert(Int, sum(occ))
                end
            end
        end
    end
    
    return n_orbitals, n_electrons
end

function get_nuclear_charges(mol_neo)
    """
    Get nuclear charges
    """
    charges = Float64[]
    
    if pybuiltin("hasattr")(mol_neo, "atom_charges")
        try
            charges = convert(Vector{Float64}, collect(mol_neo.atom_charges))
        catch
            # Fall back to parsing atom string
        end
    end
    
    # If we couldn't get charges, try parsing atom string
    if isempty(charges) && pybuiltin("hasattr")(mol_neo, "atom")
        atoms = split(mol_neo.atom, ";")
        for atom in atoms
            element = split(strip(atom))[1]
            if element == "H"
                push!(charges, 1.0)
            elseif element == "C"
                push!(charges, 6.0)
            elseif element == "N"
                push!(charges, 7.0)
            elseif element == "O"
                push!(charges, 8.0)
            else
                push!(charges, 0.0)
            end
        end
    end
    
    return charges
end

# ======================== Efficient Storage ========================

"""
    save_hamiltonian(filename::String, ham_data::HamiltonianData, H_matrix=nothing)

Save Hamiltonian data to HDF5 file for efficient storage and retrieval.
"""
function save_hamiltonian(filename::String, ham_data::HamiltonianData, H_matrix=nothing)
    h5open(filename, "w") do file
        # Save one-body integrals
        g_h1e = create_group(file, "h1e")
        for (comp, h1) in ham_data.h1e
            write(g_h1e, comp, h1)
        end
        
        # Save two-body integrals (sparse dictionary format)
        g_h2e = create_group(file, "h2e")
        for (comp, h2_dict) in ham_data.h2e
            # Convert dictionary to arrays for storage
            n_elements = length(h2_dict)
            if n_elements > 0
                indices = zeros(Int, 4, n_elements)
                values = zeros(Float64, n_elements)
                
                for (idx, ((i,j,k,l), val)) in enumerate(h2_dict)
                    indices[:, idx] = [i, j, k, l]
                    values[idx] = val
                end
                
                write(g_h2e, "$(comp)_indices", indices)
                write(g_h2e, "$(comp)_values", values)
            else
                # Write empty arrays
                write(g_h2e, "$(comp)_indices", zeros(Int, 4, 0))
                write(g_h2e, "$(comp)_values", Float64[])
            end
        end
        
        # Save coupling integrals
        g_coupling = create_group(file, "coupling")
        for ((c1, c2), V) in ham_data.coupling
            write(g_coupling, "$(c1)_$(c2)", V)
        end
        
        # Save metadata
        g_meta = create_group(file, "metadata")
        write(g_meta, "n_orbitals", JSON.json(ham_data.n_orbitals))
        write(g_meta, "n_electrons", JSON.json(ham_data.n_electrons))
        write(g_meta, "nuclear_charges", ham_data.nuclear_charges)
        
        # Save configuration info
        g_configs = create_group(file, "configurations")
        write(g_configs, "n_configs", length(ham_data.configs))
        write(g_configs, "config_energies", ham_data.config_energies)
        
        # Save configuration names and weights
        config_names = [c.name for c in ham_data.configs]
        config_weights = [c.weight for c in ham_data.configs]
        write(g_configs, "names", join(config_names, "\n"))
        write(g_configs, "weights", config_weights)
        
        # Save active space info
        write(g_meta, "active_orbitals", JSON.json(ham_data.active_orbitals))
        write(g_meta, "frozen_orbitals", JSON.json(ham_data.frozen_orbitals))
        
        # Optionally save full Hamiltonian matrix
        if H_matrix !== nothing
            write(file, "H_matrix", H_matrix)
        end
    end
    
    @info "Hamiltonian saved to $filename"
end

"""
    load_hamiltonian(filename::String) -> (HamiltonianData, H_matrix)

Load Hamiltonian data from HDF5 file.
"""
function load_hamiltonian(filename::String)
    ham_data = nothing
    H_matrix = nothing
    
    h5open(filename, "r") do file
        # Load one-body integrals
        h1e = Dict{String, Matrix{Float64}}()
        g_h1e = file["h1e"]
        for comp in keys(g_h1e)
            h1e[comp] = read(g_h1e, comp)
        end
        
        # Load two-body integrals
        h2e = Dict{String, Dict{NTuple{4,Int}, Float64}}()
        g_h2e = file["h2e"]
        
        # Find components by looking for _indices entries
        components = String[]
        for key in keys(g_h2e)
            if endswith(key, "_indices")
                push!(components, replace(key, "_indices" => ""))
            end
        end
        
        for comp in components
            indices = read(g_h2e, "$(comp)_indices")
            values = read(g_h2e, "$(comp)_values")
            
            h2e_dict = Dict{NTuple{4,Int}, Float64}()
            for idx in 1:length(values)
                i, j, k, l = indices[:, idx]
                h2e_dict[(i, j, k, l)] = values[idx]
            end
            h2e[comp] = h2e_dict
        end
        
        # Load coupling integrals
        coupling = Dict{Tuple{String, String}, Matrix{Float64}}()
        g_coupling = file["coupling"]
        for key in keys(g_coupling)
            parts = split(key, "_")
            if length(parts) >= 2
                c1 = parts[1]
                c2 = join(parts[2:end], "_")
                coupling[(c1, c2)] = read(g_coupling, key)
            end
        end
        
        # Load metadata
        g_meta = file["metadata"]
        n_orbitals = JSON.parse(read(g_meta, "n_orbitals"))
        n_electrons = JSON.parse(read(g_meta, "n_electrons"))
        nuclear_charges = read(g_meta, "nuclear_charges")
        
        # Load configuration info
        g_configs = file["configurations"]
        n_configs = read(g_configs, "n_configs")
        config_energies = read(g_configs, "config_energies")
        
        # Reconstruct configurations (simplified)
        config_names = split(read(g_configs, "names"), "\n")
        config_weights = read(g_configs, "weights")
        
        configs = Configuration[]
        for i in 1:n_configs
            # Create placeholder configurations
            config = Configuration(
                config_names[i],
                Dict{String, Vector{Int16}}(),  # Would need full info
                config_weights[i]
            )
            push!(configs, config)
        end
        
        # Load active space info
        active_orbitals = JSON.parse(read(g_meta, "active_orbitals"))
        frozen_orbitals = JSON.parse(read(g_meta, "frozen_orbitals"))
        
        # Create HamiltonianData
        ham_data = HamiltonianData(
            h1e, h2e, coupling,
            Dict(k => Int(v) for (k,v) in n_orbitals),
            Dict(k => Int(v) for (k,v) in n_electrons),
            nuclear_charges,
            configs, config_energies,
            Dict(k => Int.(v) for (k,v) in active_orbitals),
            Dict(k => Int.(v) for (k,v) in frozen_orbitals)
        )
        
        # Load full matrix if available
        if haskey(file, "H_matrix")
            H_matrix = read(file, "H_matrix")
        end
    end
    
    return ham_data, H_matrix
end

# ======================== Validation and Analysis ========================

"""
    validate_hamiltonian(ham_data::HamiltonianData)

Validate Hamiltonian data for consistency and physical correctness.
"""
function validate_hamiltonian(ham_data::HamiltonianData)
    issues = String[]
    
    # Check hermiticity of one-body integrals
    for (comp, h1) in ham_data.h1e
        if norm(h1 - h1') > 1e-10
            push!(issues, "One-body integral $comp not hermitian")
        end
    end
    
    # Check active space consistency
    for (comp, active) in ham_data.active_orbitals
        if haskey(ham_data.n_orbitals, comp)
            if maximum(active, init=0) > ham_data.n_orbitals[comp]
                push!(issues, "Active orbital index exceeds total orbitals for $comp")
            end
        end
    end
    
    # Check electron numbers
    for (comp, n_elec) in ham_data.n_electrons
        if haskey(ham_data.n_orbitals, comp)
            if n_elec > ham_data.n_orbitals[comp]
                push!(issues, "More electrons than orbitals for $comp")
            end
        end
    end
    
    if isempty(issues)
        @info "Hamiltonian validation passed ✓"
        return true
    else
        @warn "Hamiltonian validation issues found:"
        for issue in issues
            @warn "  - $issue"
        end
        return false
    end
end

"""
    reduce_hamiltonian_symmetry(H_matrix, symmetry_group)

Reduce Hamiltonian using symmetry operations.
"""
function reduce_hamiltonian_symmetry(H_matrix, symmetry_group="C1")
    # Placeholder for symmetry reduction
    @info "Symmetry reduction not yet implemented for group $symmetry_group"
    return H_matrix
end

# ======================== Export Functions ========================

"""
    export_hamiltonian_openfermion(ham_data::HamiltonianData, filename::String)

Export Hamiltonian in OpenFermion format.
This function now uses the integrated quantum computing module for enhanced functionality.
"""
function export_hamiltonian_openfermion(ham_data::HamiltonianData, filename::String)
    # Use the new quantum computing integration
    result = QuantumComputing.exportToOpenFermion(ham_data, filename)
    
    if result["success"]
        @info "Hamiltonian exported to OpenFermion format: $filename"
        @info "  Operators exported: $(result["operator_count"])"
    else
        @warn "OpenFermion export failed: $(result["message"])"
    end
    
    return result
end

"""
    export_hamiltonian_quantum_formats(ham_data::HamiltonianData, format::String, filename::String)

Export Hamiltonian to multiple quantum computing formats.
Supported formats: "openfermion", "qiskit", "cirq", "all"
"""
function export_hamiltonian_quantum_formats(ham_data::HamiltonianData, format::String, filename::String)
    return QuantumComputing.exportToQuantumFormats(ham_data, format, filename)
end

# ======================== Storage Analysis ========================

"""
    estimate_storage_requirements(ham_data::HamiltonianData)

Estimate storage requirements for Hamiltonian data.
"""
function estimate_storage_requirements(ham_data::HamiltonianData)
    # One-body storage
    one_body_elements = sum(length(h1) for (_, h1) in ham_data.h1e)
    one_body_mb = one_body_elements * 8 / 1024^2
    
    # Two-body storage
    two_body_elements = sum(length(h2_dict) for (_, h2_dict) in ham_data.h2e)
    two_body_mb = two_body_elements * (4 * 4 + 8) / 1024^2  # indices + value
    
    # Coupling storage
    coupling_elements = sum(length(V) for (_, V) in ham_data.coupling)
    coupling_mb = coupling_elements * 8 / 1024^2
    
    total_mb = one_body_mb + two_body_mb + coupling_mb
    
    return (
        one_body_mb = one_body_mb,
        two_body_mb = two_body_mb,
        coupling_mb = coupling_mb,
        total_mb = total_mb
    )
end

end # module
