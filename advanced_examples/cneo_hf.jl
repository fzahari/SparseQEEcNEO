module CNEOHF

using LinearAlgebra
using PyCall
using ..Types
using ..PySCFInterface

export CNEOCalculation, CNEOResults, run_cneo_hf

"""
Constrained NEO-HF calculation parameters
"""
struct CNEOCalculation
    method::String  # "HF" or "DFT"
    functional::String  # DFT functional if method="DFT"
    constraint_positions::Vector{Vector{Float64}}  # R₀ positions for each quantum nucleus
    convergence_threshold::Float64
    max_iterations::Int
    lambda_damping::Float64  # Damping factor for Newton updates
end

# Constructor with defaults
function CNEOCalculation(;
    method::String = "HF",
    functional::String = "B3LYP",
    constraint_positions::Vector{Vector{Float64}} = Vector{Vector{Float64}}(),
    convergence_threshold::Float64 = 1e-6,
    max_iterations::Int = 50,
    lambda_damping::Float64 = 0.5
)
    return CNEOCalculation(method, functional, constraint_positions, 
                          convergence_threshold, max_iterations, lambda_damping)
end

"""
Results from cNEO-HF calculation
"""
struct CNEOResults
    energy::Float64
    electronic_energy::Float64
    nuclear_kinetic_energy::Float64
    lagrange_multipliers::Vector{Vector{Float64}}  # λ vectors for each nucleus
    nuclear_positions::Vector{Vector{Float64}}  # Actual expectation values
    converged::Bool
    iterations::Int
    mf::PyObject  # The converged mean-field object
end

"""
Calculate expectation value of nuclear position ⟨ϕⁿ|r̂|ϕⁿ⟩
"""
function calculate_nuclear_position(mf_nuc::PyObject, origin::Vector{Float64})
    # Get density matrix
    dm = mf_nuc.make_rdm1()
    
    # Get dipole integrals (position operators) in AO basis
    mol = mf_nuc.mol
    pyscf = pyimport("pyscf")
    
    # Dipole integrals relative to origin
    with_origin = pyscf.gto.mole.with_origin(mol, origin)
    dipole_ao = with_origin.intor_symmetric("int1e_r", comp=3)
    
    # Calculate expectation values
    position = zeros(3)
    for i in 1:3
        position[i] = real(tr(dm * dipole_ao[i]))
    end
    
    return position + origin
end

"""
Add constraint term λ·r to nuclear Hamiltonian
"""
function add_constraint_to_hamiltonian!(mf_nuc::PyObject, lambda_vec::Vector{Float64}, origin::Vector{Float64})
    # Get dipole integrals
    mol = mf_nuc.mol
    pyscf = pyimport("pyscf")
    
    with_origin = pyscf.gto.mole.with_origin(mol, origin)
    dipole_ao = with_origin.intor_symmetric("int1e_r", comp=3)
    
    # Build λ·r term
    lambda_dot_r = zeros(size(dipole_ao[1]))
    for i in 1:3
        lambda_dot_r += lambda_vec[i] * dipole_ao[i]
    end
    
    # Store original get_hcore if not already stored
    if !hasproperty(mf_nuc, :_original_get_hcore)
        mf_nuc._original_get_hcore = mf_nuc.get_hcore
    end
    
    # Create new get_hcore function that includes constraint
    py"""
    def make_constrained_hcore(original_hcore, constraint_term):
        def get_hcore(mol=None):
            h = original_hcore(mol)
            return h + constraint_term
        return get_hcore
    """
    
    mf_nuc.get_hcore = py"make_constrained_hcore"(mf_nuc._original_get_hcore, lambda_dot_r)
end

"""
Calculate gradient of Lagrangian: ∇L = ⟨ϕⁿ|r̂|ϕⁿ⟩ - R₀
"""
function lagrangian_gradient(mf_nuc::PyObject, R0::Vector{Float64}, origin::Vector{Float64})
    r_current = calculate_nuclear_position(mf_nuc, origin)
    return r_current - R0
end

"""
Calculate approximate Hessian using perturbation theory (Eq. 21 in paper)
"""
function lagrangian_hessian(mf_nuc::PyObject, origin::Vector{Float64})
    # Get orbital information
    mo_occ = get(mf_nuc, :mo_occ, nothing)
    mo_energy = get(mf_nuc, :mo_energy, nothing)
    mo_coeff = get(mf_nuc, :mo_coeff, nothing)
    
    if mo_occ === nothing || mo_energy === nothing || mo_coeff === nothing
        # Fallback to identity matrix if orbital info not available
        return Matrix(1.0I, 3, 3)
    end
    
    # Find occupied and virtual orbitals
    occ_idx = findall(mo_occ .> 0)
    vir_idx = findall(mo_occ .== 0)
    
    if isempty(occ_idx) || isempty(vir_idx)
        return Matrix(1.0I, 3, 3)
    end
    
    # Get dipole integrals in MO basis
    mol = mf_nuc.mol
    pyscf = pyimport("pyscf")
    numpy = pyimport("numpy")
    
    with_origin = pyscf.gto.mole.with_origin(mol, origin)
    dipole_ao = with_origin.intor_symmetric("int1e_r", comp=3)
    
    # Transform to MO basis
    dipole_mo = []
    for i in 1:3
        d_mo = numpy.dot(mo_coeff.T, numpy.dot(dipole_ao[i-1], mo_coeff))
        push!(dipole_mo, d_mo)
    end
    
    # Calculate Hessian using Eq. (21) from paper
    H = zeros(3, 3)
    
    for occ in occ_idx
        for vir in vir_idx
            de = mo_energy[occ-1] - mo_energy[vir-1]  # Python 0-based indexing
            if abs(de) > 1e-10  # Avoid division by zero
                for i in 1:3, j in 1:3
                    r_vo = dipole_mo[i][vir-1, occ-1]
                    r_ov = dipole_mo[j][occ-1, vir-1]
                    H[i,j] += 2 * real(r_vo * r_ov) / de
                end
            end
        end
    end
    
    # Ensure Hessian is negative definite (add small negative diagonal if needed)
    eigvals = eigvals(H)
    if any(eigvals .>= 0)
        H -= (maximum(eigvals) + 0.1) * I
    end
    
    return H
end

"""
Solve for Lagrange multipliers using Newton's method
"""
function solve_lagrange_multipliers!(mf::PyObject, cneo_calc::CNEOCalculation)
    # Get nuclear components
    nuc_keys = [k for k in keys(mf.components) if startswith(string(k), "n")]
    n_nuc = length(nuc_keys)
    
    if n_nuc != length(cneo_calc.constraint_positions)
        error("Number of constraint positions ($(length(cneo_calc.constraint_positions))) must match number of quantum nuclei ($n_nuc)")
    end
    
    # Initialize Lagrange multipliers
    lambdas = [zeros(3) for _ in 1:n_nuc]
    converged = false
    iteration = 0
    
    # Get nuclear basis origins
    origins = []
    for i in 1:n_nuc
        nuc_key = "n$(i-1)"
        mf_nuc = mf.components[nuc_key]
        mol_nuc = mf_nuc.mol
        # Get the position where nuclear basis is centered
        if hasproperty(mol_nuc, :atom_coords)
            coords = mol_nuc.atom_coords()
            push!(origins, coords[1, :])  # First atom is the quantum nucleus
        else
            push!(origins, zeros(3))
        end
    end
    
    @info "Starting cNEO-HF iterations..."
    
    while !converged && iteration < cneo_calc.max_iterations
        iteration += 1
        max_error = 0.0
        
        # Update each nuclear component with current Lagrange multipliers
        for i in 1:n_nuc
            nuc_key = "n$(i-1)"
            mf_nuc = mf.components[nuc_key]
            add_constraint_to_hamiltonian!(mf_nuc, lambdas[i], origins[i])
        end
        
        # Run SCF with current constraints
        mf.kernel()
        
        # Newton update for each nucleus
        for i in 1:n_nuc
            nuc_key = "n$(i-1)"
            mf_nuc = mf.components[nuc_key]
            
            # Calculate gradient
            grad = lagrangian_gradient(mf_nuc, cneo_calc.constraint_positions[i], origins[i])
            error = norm(grad)
            max_error = max(max_error, error)
            
            # Calculate Hessian
            H = lagrangian_hessian(mf_nuc, origins[i])
            
            # Newton update with damping
            try
                delta_lambda = -H \ grad
                lambdas[i] += cneo_calc.lambda_damping * delta_lambda
            catch e
                @warn "Hessian inversion failed, using gradient descent" exception=e
                lambdas[i] -= 0.01 * grad
            end
        end
        
        converged = max_error < cneo_calc.convergence_threshold
        
        if iteration % 5 == 0 || converged
            @info "cNEO iteration $iteration: max error = $(round(max_error, digits=8))"
        end
    end
    
    return lambdas, converged, iteration
end

"""
Run constrained NEO-HF calculation
"""
function run_cneo_hf(mol::Molecule, cneo_calc::CNEOCalculation; 
                     neo_config::NEOConfig = NEOConfig())
    
    # Setup PySCF
    pyscf, has_neo = setup_pyscf(neo_config)
    if !has_neo
        error("NEO module not available")
    end
    
    # Build NEO molecule
    mol_neo = build_neo_molecule(mol, pyscf)
    
    # Create mean-field object
    if cneo_calc.method == "HF"
        mf = pyscf.neo.HF(mol_neo).density_fit()
    else
        mf = pyscf.neo.KS(mol_neo, xc=cneo_calc.functional).density_fit()
    end
    
    # Run initial SCF without constraints
    mf.kernel()
    initial_energy = mf.e_tot
    
    @info "Initial NEO energy (unconstrained): $initial_energy Ha"
    
    # Solve for Lagrange multipliers
    lambdas, converged, iterations = solve_lagrange_multipliers!(mf, cneo_calc)
    
    if !converged
        @warn "cNEO-HF did not converge within $(cneo_calc.max_iterations) iterations"
    end
    
    # Get final nuclear positions
    nuc_positions = []
    origins = []
    nuc_keys = [k for k in keys(mf.components) if startswith(string(k), "n")]
    
    for i in 1:length(nuc_keys)
        nuc_key = "n$(i-1)"
        mf_nuc = mf.components[nuc_key]
        mol_nuc = mf_nuc.mol
        
        # Get origin
        if hasproperty(mol_nuc, :atom_coords)
            coords = mol_nuc.atom_coords()
            origin = coords[1, :]
        else
            origin = zeros(3)
        end
        
        pos = calculate_nuclear_position(mf_nuc, origin)
        push!(nuc_positions, pos)
    end
    
    # Calculate energy components
    electronic_energy = mf.energy_elec()[1]
    nuclear_kinetic = mf.e_tot - electronic_energy
    
    # Create results
    results = CNEOResults(
        mf.e_tot,
        electronic_energy,
        nuclear_kinetic,
        lambdas,
        nuc_positions,
        converged,
        iterations,
        mf
    )
    
    @info "cNEO-HF calculation complete:"
    @info "  Final energy: $(results.energy) Ha"
    @info "  Electronic energy: $(results.electronic_energy) Ha"
    @info "  Nuclear kinetic energy: $(results.nuclear_kinetic_energy) Ha"
    @info "  Converged: $(results.converged)"
    @info "  Iterations: $(results.iterations)"
    
    for i in 1:length(lambdas)
        @info "  Nucleus $i:"
        @info "    Constraint position: $(cneo_calc.constraint_positions[i])"
        @info "    Actual position: $(nuc_positions[i])"
        @info "    Lagrange multiplier: $(lambdas[i])"
        @info "    Position error: $(norm(nuc_positions[i] - cneo_calc.constraint_positions[i]))"
    end
    
    return results
end

end # module
