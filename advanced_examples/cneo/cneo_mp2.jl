module CNEOMP2

using LinearAlgebra
using PyCall
using ..Types
using ..PySCFInterface
# Note: CNEOHF functions will be imported at module level after both modules are loaded

export CNEOMP2Calculation, CNEOMP2Results, run_cneo_mp2

"""
Constrained NEO-MP2 calculation parameters
"""
struct CNEOMP2Calculation
    constraint_positions::Vector{Vector{Float64}}  # R₀ positions for each quantum nucleus
    constraint_at_mp2::Bool  # Whether to apply constraint at MP2 level
    convergence_threshold::Float64
    mp2_convergence_threshold::Float64
    max_iterations::Int
    max_mp2_iterations::Int
    lambda_damping::Float64
    use_canonical::Bool  # Use canonical or non-canonical MP2
end

# Constructor with defaults
function CNEOMP2Calculation(;
    constraint_positions::Vector{Vector{Float64}} = Vector{Vector{Float64}}(),
    constraint_at_mp2::Bool = false,  # By default, only constrain at HF level
    convergence_threshold::Float64 = 1e-6,
    mp2_convergence_threshold::Float64 = 1e-8,
    max_iterations::Int = 50,
    max_mp2_iterations::Int = 30,
    lambda_damping::Float64 = 0.5,
    use_canonical::Bool = false
)
    return CNEOMP2Calculation(
        constraint_positions, constraint_at_mp2, convergence_threshold,
        mp2_convergence_threshold, max_iterations, max_mp2_iterations,
        lambda_damping, use_canonical
    )
end

"""
Results from cNEO-MP2 calculation
"""
struct CNEOMP2Results
    hf_energy::Float64
    mp2_energy::Float64
    total_energy::Float64
    mp2_correlation::Float64
    lagrange_multipliers_hf::Vector{Vector{Float64}}
    lagrange_multipliers_mp2::Vector{Vector{Float64}}
    nuclear_positions_hf::Vector{Vector{Float64}}
    nuclear_positions_mp2::Vector{Vector{Float64}}
    converged::Bool
    iterations::Int
    mp2_iterations::Int
    mf::PyObject
end

"""
Build non-canonical Fock matrices for MP2
"""
function build_noncanonical_fock(mf::PyObject)
    pyscf = pyimport("pyscf")
    numpy = pyimport("numpy")
    
    # Get electronic component
    mf_elec = mf.components["e"]
    mo_coeff_e = mf_elec.mo_coeff
    mo_occ_e = mf_elec.mo_occ
    h1e = mf_elec.get_hcore()
    
    # Build electronic Fock matrix in AO basis
    dm_e = mf_elec.make_rdm1()
    vhf_e = mf_elec.get_veff(mf_elec.mol, dm_e)
    fock_e_ao = h1e + vhf_e
    
    # Transform to MO basis
    fock_e_mo = numpy.dot(mo_coeff_e.T, numpy.dot(fock_e_ao, mo_coeff_e))
    
    # Build nuclear Fock matrices
    fock_n_mo = Dict{String, PyObject}()
    
    for (key, mf_nuc) in mf.components
        if startswith(string(key), "n")
            mo_coeff_n = mf_nuc.mo_coeff
            mo_occ_n = mf_nuc.mo_occ
            h1n = mf_nuc.get_hcore()
            
            # Nuclear Fock only includes kinetic + external potential
            dm_n = mf_nuc.make_rdm1()
            
            # For NEO, nuclear Fock includes electron-nuclear attraction
            mol_nuc = mf_nuc.mol
            mol_elec = mf_elec.mol
            
            # Calculate electron density at nuclear grid points
            eri_en = pyscf.neo.eri.get_eri_ep(mol_elec, mol_nuc)
            ven = numpy.einsum("pqrs,rs->pq", eri_en, dm_e)
            
            fock_n_ao = h1n + ven
            fock_n_mo[string(key)] = numpy.dot(mo_coeff_n.T, numpy.dot(fock_n_ao, mo_coeff_n))
        end
    end
    
    return fock_e_mo, fock_n_mo
end

"""
Calculate MP2 t-amplitudes iteratively for non-canonical case
"""
function calculate_mp2_amplitudes_noncanonical(mf::PyObject, fock_e::PyObject, 
                                              fock_n::Dict{String, PyObject}, 
                                              max_iter::Int, conv_tol::Float64)
    pyscf = pyimport("pyscf")
    numpy = pyimport("numpy")
    
    # Get orbital information
    mf_elec = mf.components["e"]
    nocc_e = Int(sum(mf_elec.mo_occ) ÷ 2)  # Assuming closed shell
    nmo_e = size(mf_elec.mo_coeff, 2)
    nvir_e = nmo_e - nocc_e
    
    # Initialize t-amplitudes
    t2_ee = numpy.zeros((nocc_e, nocc_e, nvir_e, nvir_e))  # Electronic doubles
    t2_en = Dict{String, PyObject}()  # Electron-nuclear
    t2_nn = Dict{String, Dict{String, PyObject}}()  # Nuclear doubles
    
    # Get two-electron integrals
    mo_coeff_e = mf_elec.mo_coeff
    eri_mo = pyscf.ao2mo.kernel(mf_elec._eri, mo_coeff_e)
    eri_mo = eri_mo.reshape(nmo_e, nmo_e, nmo_e, nmo_e)
    
    # Initial guess from diagonal elements
    for i in 1:nocc_e, j in 1:nocc_e, a in 1:nvir_e, b in 1:nvir_e
        denom = (fock_e[i-1,i-1] + fock_e[j-1,j-1] - 
                 fock_e[nocc_e+a-1,nocc_e+a-1] - fock_e[nocc_e+b-1,nocc_e+b-1])
        if abs(denom) > 1e-10
            # <ij||ab> = <ij|ab> - <ij|ba>
            v_ijab = eri_mo[i-1,j-1,nocc_e+a-1,nocc_e+b-1] - 
                     eri_mo[i-1,j-1,nocc_e+b-1,nocc_e+a-1]
            t2_ee[i-1,j-1,a-1,b-1] = v_ijab / denom
        end
    end
    
    # Iterate to convergence
    converged = false
    iteration = 0
    
    while !converged && iteration < max_iter
        iteration += 1
        t2_old = numpy.copy(t2_ee)
        
        # Update amplitudes using residual equations
        for i in 1:nocc_e, j in 1:nocc_e, a in 1:nvir_e, b in 1:nvir_e
            # Calculate residual R_ijab
            r_ijab = eri_mo[i-1,j-1,nocc_e+a-1,nocc_e+b-1] - 
                     eri_mo[i-1,j-1,nocc_e+b-1,nocc_e+a-1]
            
            # Add Fock contributions
            for c in 1:nvir_e
                if c != a
                    r_ijab += fock_e[nocc_e+a-1,nocc_e+c-1] * t2_ee[i-1,j-1,c-1,b-1]
                end
                if c != b
                    r_ijab += fock_e[nocc_e+b-1,nocc_e+c-1] * t2_ee[i-1,j-1,a-1,c-1]
                end
            end
            
            for k in 1:nocc_e
                if k != i
                    r_ijab -= fock_e[k-1,i-1] * t2_ee[k-1,j-1,a-1,b-1]
                end
                if k != j
                    r_ijab -= fock_e[k-1,j-1] * t2_ee[i-1,k-1,a-1,b-1]
                end
            end
            
            # Update amplitude
            denom = (fock_e[i-1,i-1] + fock_e[j-1,j-1] - 
                     fock_e[nocc_e+a-1,nocc_e+a-1] - fock_e[nocc_e+b-1,nocc_e+b-1])
            if abs(denom) > 1e-10
                t2_ee[i-1,j-1,a-1,b-1] = r_ijab / denom
            end
        end
        
        # Check convergence
        error = numpy.linalg.norm(t2_ee - t2_old)
        converged = error < conv_tol
    end
    
    return t2_ee, converged, iteration
end

"""
Calculate MP2 energy from amplitudes
"""
function calculate_mp2_energy(mf::PyObject, t2_ee::PyObject)
    pyscf = pyimport("pyscf")
    numpy = pyimport("numpy")
    
    # Get orbital information
    mf_elec = mf.components["e"]
    mo_coeff_e = mf_elec.mo_coeff
    nocc_e = Int(sum(mf_elec.mo_occ) ÷ 2)
    nmo_e = size(mo_coeff_e, 2)
    
    # Get two-electron integrals in MO basis
    eri_mo = pyscf.ao2mo.kernel(mf_elec._eri, mo_coeff_e)
    eri_mo = eri_mo.reshape(nmo_e, nmo_e, nmo_e, nmo_e)
    
    # Calculate MP2 correlation energy
    emp2 = 0.0
    for i in 1:nocc_e, j in 1:nocc_e, a in 1:(nmo_e-nocc_e), b in 1:(nmo_e-nocc_e)
        v_ijab = eri_mo[i-1,j-1,nocc_e+a-1,nocc_e+b-1]
        v_ijba = eri_mo[i-1,j-1,nocc_e+b-1,nocc_e+a-1]
        
        emp2 += 0.25 * t2_ee[i-1,j-1,a-1,b-1] * (2*v_ijab - v_ijba)
    end
    
    return emp2
end

"""
Calculate MP2 density matrix
"""
function calculate_mp2_density(mf::PyObject, t2_ee::PyObject)
    numpy = pyimport("numpy")
    
    # Get orbital information
    mf_elec = mf.components["e"]
    nocc_e = Int(sum(mf_elec.mo_occ) ÷ 2)
    nmo_e = size(mf_elec.mo_coeff, 2)
    nvir_e = nmo_e - nocc_e
    
    # Initialize density matrix in MO basis
    dm_mo = numpy.zeros((nmo_e, nmo_e))
    
    # HF density contribution
    for i in 1:nocc_e
        dm_mo[i-1, i-1] = 2.0  # Double occupancy
    end
    
    # MP2 correction to density
    # Occupied-occupied block
    for i in 1:nocc_e, j in 1:nocc_e
        for a in 1:nvir_e, b in 1:nvir_e
            dm_mo[i-1, j-1] -= 0.5 * numpy.sum(
                t2_ee[i-1, :, a-1, b-1] * numpy.conj(t2_ee[j-1, :, a-1, b-1])
            )
        end
    end
    
    # Virtual-virtual block
    for a in 1:nvir_e, b in 1:nvir_e
        for i in 1:nocc_e, j in 1:nocc_e
            dm_mo[nocc_e+a-1, nocc_e+b-1] += 0.5 * numpy.sum(
                t2_ee[i-1, j-1, a-1, :] * numpy.conj(t2_ee[i-1, j-1, b-1, :])
            )
        end
    end
    
    return dm_mo
end

"""
Run constrained NEO-MP2 calculation
"""
function run_cneo_mp2(mol::Molecule, mp2_calc::CNEOMP2Calculation; 
                      neo_config::NEOConfig = NEOConfig())
    
    # Setup PySCF
    pyscf, has_neo = setup_pyscf(neo_config)
    if !has_neo
        error("NEO module not available")
    end
    
    # Build NEO molecule
    mol_neo = build_neo_molecule(mol, pyscf)
    
    # First run cNEO-HF
    @info "Running cNEO-HF..."
    mf = pyscf.neo.HF(mol_neo).density_fit()
    
    # Create cNEO-HF calculation object
    cneo_hf = CNEOHF.CNEOCalculation(
        method = "HF",
        constraint_positions = mp2_calc.constraint_positions,
        convergence_threshold = mp2_calc.convergence_threshold,
        max_iterations = mp2_calc.max_iterations,
        lambda_damping = mp2_calc.lambda_damping
    )
    
    # Solve for HF Lagrange multipliers
    lambdas_hf, converged_hf, iterations_hf = solve_lagrange_multipliers!(mf, cneo_hf)
    
    if !converged_hf
        @warn "cNEO-HF did not converge"
    end
    
    hf_energy = mf.e_tot
    @info "cNEO-HF energy: $hf_energy Ha"
    
    # Get HF nuclear positions
    nuclear_positions_hf = []
    nuc_keys = [k for k in keys(mf.components) if startswith(string(k), "n")]
    
    for i in 1:length(nuc_keys)
        nuc_key = "n$(i-1)"
        mf_nuc = mf.components[nuc_key]
        mol_nuc = mf_nuc.mol
        
        if hasproperty(mol_nuc, :atom_coords)
            origin = mol_nuc.atom_coords()[1, :]
        else
            origin = zeros(3)
        end
        
        pos = calculate_nuclear_position(mf_nuc, origin)
        push!(nuclear_positions_hf, pos)
    end
    
    # Now run MP2
    @info "Running MP2 correlation..."
    
    if mp2_calc.use_canonical
        # Use PySCF's canonical MP2
        mp2_obj = pyscf.neo.mp2.MP2(mf)
        emp2, t2 = mp2_obj.kernel()
        converged_mp2 = true
        mp2_iterations = 1
    else
        # Non-canonical MP2
        fock_e, fock_n = build_noncanonical_fock(mf)
        t2_ee, converged_mp2, mp2_iterations = calculate_mp2_amplitudes_noncanonical(
            mf, fock_e, fock_n, 
            mp2_calc.max_mp2_iterations, 
            mp2_calc.mp2_convergence_threshold
        )
        emp2 = calculate_mp2_energy(mf, t2_ee)
    end
    
    @info "MP2 correlation energy: $emp2 Ha"
    total_energy = hf_energy + emp2
    @info "Total cNEO-MP2 energy: $total_energy Ha"
    
    # If constraining at MP2 level, update densities and positions
    lambdas_mp2 = [zeros(3) for _ in 1:length(nuc_keys)]
    nuclear_positions_mp2 = copy(nuclear_positions_hf)
    
    if mp2_calc.constraint_at_mp2
        @info "Applying constraints at MP2 level..."
        # This would require iterating MP2 with density constraints
        # For now, we use the HF constraint positions
        @warn "MP2-level constraints not fully implemented yet"
    end
    
    # Create results
    results = CNEOMP2Results(
        hf_energy,
        emp2,
        total_energy,
        emp2,  # MP2 correlation
        lambdas_hf,
        lambdas_mp2,
        nuclear_positions_hf,
        nuclear_positions_mp2,
        converged_hf && converged_mp2,
        iterations_hf,
        mp2_iterations,
        mf
    )
    
    return results
end

end # module
