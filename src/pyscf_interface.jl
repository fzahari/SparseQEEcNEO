"""
pyscf_interface.jl - Interface to PySCF and NEO functionality
"""
module PySCFInterface

using PyCall
using Printf
using ..Types  # Use parent module, not include

export setup_pyscf, build_neo_molecule, run_neo_meanfield
export run_neo_mp2, run_neo_casci, extract_t2_amplitudes

# ======================== PySCF Setup ========================

function setup_pyscf(config::NEOConfig)
   py"""
   import sys
   
   # Add the PySCF path
   pyscf_path = $config.pyscf_path
   if pyscf_path not in sys.path:
       sys.path.insert(0, pyscf_path)
   
   # Import pyscf
   import pyscf
   print(f"Using PySCF from: {pyscf.__file__}")
   
   # Import NEO explicitly
   try:
       import pyscf.neo
       print("NEO module loaded successfully!")
       neo_available = True
   except ImportError as e:
       print(f"NEO module not available: {e}")
       neo_available = False
   """
   
   pyscf = pyimport("pyscf")
   
   # Check if NEO was imported successfully
   has_neo = py"neo_available"
   
   if has_neo
       # Import NEO module for use
       neo = pyimport("pyscf.neo")
       @info "✓ NEO module available and imported!"
   else
       @warn "NEO module not available"
   end
   
   return pyscf, has_neo
end

# ======================== NEO Molecule Building ========================

function build_neo_molecule(mol::Molecule, pyscf)
    # Import neo module
    neo = pyimport("pyscf.neo")
    
    # Convert quantum_nuc to Python list
    quantum_nuc_py = PyVector(mol.quantum_nuc)
    
    # Build NEO molecule
    mol_neo = neo.M(
        atom=mol.atom_string,
        basis=mol.basis,
        charge=mol.charge,
        spin=mol.spin,
        quantum_nuc=quantum_nuc_py,
        nuc_basis=mol.nuc_basis
    )
    
    @info "NEO molecule built with $(length(mol.quantum_nuc)) quantum nuclei"
    
    return mol_neo
end

# ======================== Mean-field Calculations ========================

function run_neo_meanfield(mol_neo, calc::NEOCalculation, pyscf)
    """
    Run NEO mean-field calculation (HF or DFT) using the same approach as the attached script
    """
    neo = pyscf.neo
    
    println("Running cNEO-$(calc.xc) calculation...")
    
    if calc.xc == "HF"
        mf = neo.HF(mol_neo)
    else
        # Use CDFT for DFT calculations
        mf = neo.CDFT(mol_neo, xc=calc.xc, epc=calc.epc)
    end
    
    # Set convergence criteria
    if pybuiltin("hasattr")(mf, "conv_tol")
        mf.conv_tol = 1e-8
    end
    if pybuiltin("hasattr")(mf, "max_cycle")
        mf.max_cycle = 100
    end
    if pybuiltin("hasattr")(mf, "diis")
        mf.diis = true
    end
    
    # Run calculation
    mf.kernel()
    
    if pybuiltin("hasattr")(mf, "converged") && mf.converged
        @info "cNEO-$(calc.xc) converged! Energy: $(mf.e_tot) Ha"
    else
        @warn "cNEO-$(calc.xc) did not converge fully"
    end
    
    # Debug nuclear components if available
    if pybuiltin("hasattr")(mf, "components")
        debug_neo_components(mf, mol_neo)
    end
    
    return mf
end

function debug_neo_components(mf, mol_neo)
    """
    Debug information for NEO components
    """
    if !pybuiltin("hasattr")(mol_neo, "nuc_num") || mol_neo.nuc_num == 0
        return
    end
    
    @info "NEO Component Analysis:"
    
    # Electronic component
    if haskey(mf.components, "e")
        mf_e = mf.components["e"]
        if pybuiltin("hasattr")(mf_e, "mo_occ")
            n_elec = sum(collect(mf_e.mo_occ) .> 0.5)
            @info "  Electronic: $n_elec occupied orbitals"
        end
    end
    
    # Nuclear components - using the actual quantum_nuc indices
    if pybuiltin("hasattr")(mol_neo, "quantum_nuc")
        for (idx, nuc_idx) in enumerate(mol_neo.quantum_nuc)
            nuc_key = "n$nuc_idx"
            if haskey(mf.components, nuc_key)
                mf_n = mf.components[nuc_key]
                if pybuiltin("hasattr")(mf_n, "mo_occ")
                    n_occ = sum(collect(mf_n.mo_occ) .> 0.5)
                    n_orbs = length(mf_n.mo_occ)
                    @info "  Nuclear $idx (index $nuc_idx): $n_occ occupied, $n_orbs total orbitals"
                    
                    if pybuiltin("hasattr")(mf_n, "mo_energy")
                        energies = collect(mf_n.mo_energy)[1:min(5, end)]
                        @info "    First orbital energies: $energies"
                    end
                end
            else
                @warn "  Nuclear component $nuc_key not found!"
            end
        end
    end
end

# ======================== Post-HF Methods ========================

function run_neo_mp2(mf, mol_neo, frozen=nothing)
    """
    Run NEO-MP2 calculation
    """
    pyscf = pyimport("pyscf")
    
    println("Running NEO-MP2 calculation...")
    
    # Create MP2 object
    if pybuiltin("hasattr")(pyscf.neo, "MP2")
        mp = pyscf.neo.MP2(mf)
        
        # Set frozen core if specified
        if frozen !== nothing
            mp.frozen = frozen
        end
        
        # Run MP2
        ecorr, t2 = mp.kernel()
        @info "NEO-MP2 correlation energy: $ecorr Ha"
        
    else
        # Fallback to electronic MP2 only
        @warn "NEO-MP2 not available, using electronic MP2 only"
        
        mp_e = pyscf.mp.MP2(mf.components["e"])
        if frozen !== nothing
            mp_e.frozen = frozen
        end
        
        ecorr, t2 = mp_e.kernel()
        @info "Electronic MP2 correlation energy: $ecorr Ha"
        
        mp = mp_e
    end
    
    return mp, ecorr, t2
end

function extract_t2_amplitudes(t2, config_sel::ConfigSelection)
    """
    Extract significant t2 amplitudes from MP2
    """
    if t2 === nothing
        return nothing
    end
    
    t2_dict = Dict{Tuple{Int,Int,Int,Int}, Float64}()
    threshold = config_sel.t2_threshold
    
    # Handle different t2 formats
    try
        if isa(t2, PyObject) && pybuiltin("hasattr")(t2, "shape")
            if length(t2.shape) == 4
                nocc, _, nvir, _ = t2.shape
                t2_array = convert(Array{Float64,4}, t2)
                
                for i in 1:nocc, j in 1:nocc, a in 1:nvir, b in 1:nvir
                    amp = t2_array[i,j,a,b]
                    if abs(amp) > threshold
                        t2_dict[(i,j,a,b)] = amp
                    end
                end
            end
        end
    catch e
        @warn "Could not extract t2 amplitudes: $e"
    end
    
    @info "Extracted $(length(t2_dict)) significant t2 amplitudes"
    return t2_dict
end

function run_neo_casci(mf, mol_neo, config_sel::ConfigSelection)
    """
    Run NEO-CASCI calculation
    """
    pyscf = pyimport("pyscf")
    
    println("Setting up CASCI calculation...")
    
    # Determine active space
    mf_e = mf.components["e"]
    n_elec_occ = sum(collect(mf_e.mo_occ) .> 0.5)
    n_orbs = length(mf_e.mo_occ)
    
    # Simple active space selection
    ncas = min(config_sel.active_orbs, n_orbs - n_elec_occ ÷ 2)
    nelecas = min(config_sel.active_elec, n_elec_occ)
    
    @info "Active space: ($nelecas e, $ncas o)"
    
    # Create CASCI object
    if pybuiltin("hasattr")(pyscf.neo, "CASCI")
        mc = pyscf.neo.CASCI(mf, ncas, nelecas)
    else
        # Fallback to electronic CASCI
        @warn "NEO-CASCI not available, using electronic CASCI"
        mc = pyscf.mcscf.CASCI(mf_e, ncas, nelecas)
    end
    
    # Set parameters
    mc.fcisolver.max_cycle = 100
    mc.fcisolver.conv_tol = 1e-8
    
    # Multiple states if requested
    if config_sel.cas_nstates > 1
        mc.nstate = config_sel.cas_nstates
    end
    
    # Run CASCI
    try
        mc.kernel()
        @info "CASCI calculation completed"
        return mc
    catch e
        @warn "CASCI failed: $e"
        return nothing
    end
end

end # module
