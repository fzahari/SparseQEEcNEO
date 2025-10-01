"""
epc_functionals.jl - Implementation of electron-proton correlation functionals
Based on Brorsen et al. J. Chem. Phys. 149, 044110 (2018)
"""

using Statistics

module EPCFunctionals

using PyCall
using LinearAlgebra
using Statistics
using ..Types

export calculate_epc_energy, get_electron_density, get_proton_density

# ======================== Main EPC Energy Function ========================

function calculate_epc_energy(mf, mol_neo, epc_type::String)
    """Calculate electron-proton correlation energy using specified functional"""
    
    if !isValidEPCType(epc_type)
        return 0.0
    end
    
    functionalParams = extractEPCParameters(epc_type)
    electronDensity, protonDensity = extractDensities(mf, mol_neo)
    
    if !areDensitiesValid(electronDensity, protonDensity)
        @warn "Invalid densities for EPC calculation"
        return 0.0
    end
    
    epcEnergy = computeEPCEnergy(electronDensity, protonDensity, epc_type, functionalParams)
    
    @info "EPC energy contribution: $epcEnergy Ha"
    return epcEnergy
end

function isValidEPCType(epcType::String)
    return epcType != "none" && haskey(EPC_PARAMS, epcType)
end

function extractEPCParameters(epcType::String)
    @info "Calculating EPC energy with functional: $epcType"
    params = EPC_PARAMS[epcType]
    return params["a"], params["b"], params["c"]
end

function extractDensities(mf, molNeo)
    electronDensity = get_electron_density(mf)
    protonDensity = get_proton_density(mf, molNeo)
    return electronDensity, protonDensity
end

function areDensitiesValid(electronDensity, protonDensity)
    return length(electronDensity) > 0 && length(protonDensity) > 0
end

function computeEPCEnergy(electronDensity, protonDensity, epcType::String, params)
    a, b, c = params
    
    if isEPC17Type(epcType)
        return calculate_epc17_energy(electronDensity, protonDensity, a, b, c)
    elseif isEPC18Type(epcType)
        return calculate_epc18_energy(electronDensity, protonDensity, a, b, c)
    else
        @warn "Unknown EPC functional: $epcType"
        return 0.0
    end
end

function isEPC17Type(epcType::String)
    return epcType in ["17-1", "17-2"]
end

function isEPC18Type(epcType::String)
    return epcType in ["18-1", "18-2"]
end

# ======================== EPC17 Functional ========================

function calculate_epc17_energy(rho_e, rho_p, a, b, c)
    """
    Calculate epc17 energy using geometric mean of inverse Wigner-Seitz radii
    
    E_epc17 = -∫ ρᵉ(r)ρᵖ(r) / D(r) dr
    where D(r) = a - b[ρᵉ(r)]^(1/2)[ρᵖ(r)]^(1/2) + c ρᵉ(r)ρᵖ(r)
    """
    
    # Ensure compatible dimensions
    n_points = min(length(rho_e), length(rho_p))
    if n_points == 0
        return 0.0
    end
    
    # Resize if needed
    if length(rho_e) != length(rho_p)
        rho_e = rho_e[1:n_points]
        rho_p = rho_p[1:n_points]
    end
    
    # Calculate energy density at each point
    energy_density = zeros(n_points)
    
    for i in 1:n_points
        if rho_e[i] > 1e-10 && rho_p[i] > 1e-10
            # Geometric mean term
            geom_mean = sqrt(rho_e[i] * rho_p[i])
            
            # Denominator D(r)
            D = a - b * sqrt(geom_mean) + c * rho_e[i] * rho_p[i]
            
            if D > 1e-10
                energy_density[i] = -rho_e[i] * rho_p[i] / D
            end
        end
    end
    
    # Integrate (simplified - should use proper grid integration)
    epc_energy = sum(energy_density) * (1.0 / n_points)
    
    return epc_energy
end

# ======================== EPC18 Functional ========================

function calculate_epc18_energy(rho_e, rho_p, a, b, c)
    """
    Calculate epc18 energy using arithmetic mean of inverse Wigner-Seitz radii
    
    Different correlation length approximation:
    β = ([ρᵉ]^(1/3) + [ρᵖ]^(1/3))/2
    """
    
    # Ensure compatible dimensions
    n_points = min(length(rho_e), length(rho_p))
    if n_points == 0
        return 0.0
    end
    
    # Resize if needed
    if length(rho_e) != length(rho_p)
        rho_e = rho_e[1:n_points]
        rho_p = rho_p[1:n_points]
    end
    
    # Calculate energy density
    energy_density = zeros(n_points)
    
    for i in 1:n_points
        if rho_e[i] > 1e-10 && rho_p[i] > 1e-10
            # Arithmetic mean of inverse Wigner-Seitz radii
            beta = (rho_e[i]^(1/3) + rho_p[i]^(1/3)) / 2
            
            # Modified denominator
            D = a - b * beta + c * beta^2
            
            if D > 1e-10
                energy_density[i] = -rho_e[i] * rho_p[i] / D
            end
        end
    end
    
    # Integrate
    epc_energy = sum(energy_density) * (1.0 / n_points)
    
    return epc_energy
end

# ======================== Density Extraction ========================

function get_electron_density(mf)
    """
    Extract electron density from mean-field object
    """
    try
        # Get electronic component
        if pybuiltin("hasattr")(mf, "components") && haskey(mf.components, "e")
            mf_e = mf.components["e"]
        else
            mf_e = mf
        end
        
        # Get density matrix
        if pybuiltin("hasattr")(mf_e, "make_rdm1")
            dm = mf_e.make_rdm1()
            
            # Convert to real array
            dm_real = real(convert(Array{Float64}, dm))
            
            # For now, return diagonal (orbital densities)
            # In full implementation, would transform to real space
            return diag(dm_real)
        else
            @warn "Cannot extract electron density - no make_rdm1 method"
            return Float64[]
        end
        
    catch e
        @warn "Error extracting electron density: $e"
        return Float64[]
    end
end

function get_proton_density(mf, mol_neo)
    """
    Extract proton density from NEO mean-field object
    """
    try
        # Check if there are quantum nuclei
        if !pybuiltin("hasattr")(mol_neo, "nuc_num") || mol_neo.nuc_num == 0
            return Float64[]
        end
        
        # Get first nuclear component (for single proton)
        nuc_key = "n0"
        if !haskey(mf.components, nuc_key)
            @warn "Nuclear component $nuc_key not found"
            return Float64[]
        end
        
        mf_n = mf.components[nuc_key]
        
        # Get density matrix
        if pybuiltin("hasattr")(mf_n, "make_rdm1")
            dm = mf_n.make_rdm1()
            
            # Convert to real array
            dm_real = real(convert(Array{Float64}, dm))
            
            # Return diagonal
            return diag(dm_real)
        else
            @warn "Cannot extract proton density - no make_rdm1 method"
            return Float64[]
        end
        
    catch e
        @warn "Error extracting proton density: $e"
        return Float64[]
    end
end

# ======================== Advanced EPC Features ========================

function calculate_epc_gradient(mf, mol_neo, epc_type::String)
    """
    Calculate gradient of EPC energy (for geometry optimization)
    Not implemented yet - placeholder for future development
    """
    @warn "EPC gradient calculation not yet implemented"
    return zeros(3)
end

function optimize_epc_parameters(reference_data, epc_type::String)
    """
    Optimize EPC parameters based on reference data
    Not implemented yet - placeholder for future development
    """
    @warn "EPC parameter optimization not yet implemented"
    return EPC_PARAMS[epc_type]
end

# ======================== Testing Functions ========================

function testEPCFunctionals()
    """
    Test EPC functional implementations with simple densities
    """
    println("\nTesting EPC Functionals")
    println("="^40)
    
    # Create test densities
    n_points = 10
    rho_e_test = ones(n_points) * 0.1  # Uniform electron density
    rho_p_test = exp.(-collect(0:n_points-1) / 2)  # Exponential proton density
    
    println("Test electron density: uniform 0.1")
    println("Test proton density: exponential decay")
    
    # Test each functional
    for (name, params) in EPC_PARAMS
        if name == "none"
            continue
        end
        
        a, b, c = params["a"], params["b"], params["c"]
        
        if name in ["17-1", "17-2"]
            energy = calculate_epc17_energy(rho_e_test, rho_p_test, a, b, c)
        elseif name in ["18-1", "18-2"]
            energy = calculate_epc18_energy(rho_e_test, rho_p_test, a, b, c)
        else
            continue
        end
        
        println("$name: E_epc = $energy Ha")
    end
    
    println("="^40)
end

end # module
