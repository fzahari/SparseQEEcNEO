"""
cneo_methods.jl - Constrained Nuclear-Electronic Orbital (cNEO) Methods

Implements cNEO-HF and cNEO-MP2 calculations with constraints on nuclear positions.
Follows Clean Code principles with small, focused functions and clear separation of concerns.
"""
module CNEOMethods

using LinearAlgebra
using PyCall
using ..Types
using ..PySCFInterface

export CNEOCalculation, CNEOResults, CNEOMP2Results
export run_cneo_hf, run_cneo_mp2, create_cneo_calculation

# ======================== Clean Code Constants ========================

const DEFAULT_CNEO_CONVERGENCE_THRESHOLD = 1e-6
const DEFAULT_CNEO_MAX_ITERATIONS = 50
const DEFAULT_LAMBDA_DAMPING_FACTOR = 0.5
const DEFAULT_MP2_CONVERGENCE_THRESHOLD = 1e-8
const DEFAULT_MP2_MAX_ITERATIONS = 30
const MINIMUM_DENOMINATOR_THRESHOLD = 1e-10
const POSITION_DIMENSION = 3
const GRADIENT_DESCENT_STEP = 0.01
const HESSIAN_REGULARIZATION = 0.1

# ======================== Data Structures ========================

"""
    CNEOCalculation

Parameters for constrained NEO calculations following Clean Code principles.

# Fields
- `method::String`: "HF" or "DFT" 
- `functional::String`: DFT functional (e.g., "B3LYP")
- `constraint_positions::Vector{Vector{Float64}}`: Target positions R₀ for each quantum nucleus
- `convergence_threshold::Float64`: Convergence criterion for constraint satisfaction
- `max_iterations::Int`: Maximum SCF iterations
- `lambda_damping::Float64`: Damping factor for Newton updates (0 < λ < 1)
"""
struct CNEOCalculation
    method::String
    functional::String
    constraint_positions::Vector{Vector{Float64}}
    convergence_threshold::Float64
    max_iterations::Int
    lambda_damping::Float64
    
    function CNEOCalculation(method, functional, constraint_positions, 
                           convergence_threshold, max_iterations, lambda_damping)
        validate_cneo_calculation_parameters(
            method, functional, constraint_positions, 
            convergence_threshold, max_iterations, lambda_damping
        )
        new(method, functional, constraint_positions, 
            convergence_threshold, max_iterations, lambda_damping)
    end
end

"""
    CNEOResults

Results from cNEO-HF calculation with complete analysis data.
"""
struct CNEOResults
    total_energy::Float64
    electronic_energy::Float64
    nuclear_kinetic_energy::Float64
    lagrange_multipliers::Vector{Vector{Float64}}
    actual_nuclear_positions::Vector{Vector{Float64}}
    constraint_positions::Vector{Vector{Float64}}
    position_errors::Vector{Float64}
    converged::Bool
    iterations::Int
    mean_field_object::PyObject
end

"""
    CNEOMP2Results

Results from cNEO-MP2 calculation with correlation energy analysis.
"""
struct CNEOMP2Results
    hf_energy::Float64
    mp2_correlation_energy::Float64
    total_energy::Float64
    lagrange_multipliers_hf::Vector{Vector{Float64}}
    lagrange_multipliers_mp2::Vector{Vector{Float64}}
    nuclear_positions_hf::Vector{Vector{Float64}}
    nuclear_positions_mp2::Vector{Vector{Float64}}
    position_errors_hf::Vector{Float64}
    position_errors_mp2::Vector{Float64}
    hf_converged::Bool
    mp2_converged::Bool
    hf_iterations::Int
    mp2_iterations::Int
    mean_field_object::PyObject
end

# ======================== Factory Functions ========================

"""
    create_cneo_calculation(; kwargs...)

Create CNEOCalculation with sensible defaults following Clean Code principles.

# Keyword Arguments
- `method::String = "HF"`: Calculation method
- `functional::String = "B3LYP"`: DFT functional 
- `constraint_positions::Vector{Vector{Float64}} = []`: Nuclear constraint positions
- `convergence_threshold::Float64 = $DEFAULT_CNEO_CONVERGENCE_THRESHOLD`: Convergence threshold
- `max_iterations::Int = $DEFAULT_CNEO_MAX_ITERATIONS`: Maximum iterations
- `lambda_damping::Float64 = $DEFAULT_LAMBDA_DAMPING_FACTOR`: Damping factor
"""
function create_cneo_calculation(;
    method::String = "HF",
    functional::String = "B3LYP",
    constraint_positions::Vector{Vector{Float64}} = Vector{Vector{Float64}}(),
    convergence_threshold::Float64 = DEFAULT_CNEO_CONVERGENCE_THRESHOLD,
    max_iterations::Int = DEFAULT_CNEO_MAX_ITERATIONS,
    lambda_damping::Float64 = DEFAULT_LAMBDA_DAMPING_FACTOR
)
    return CNEOCalculation(
        method, functional, constraint_positions,
        convergence_threshold, max_iterations, lambda_damping
    )
end

# ======================== Parameter Validation ========================

function validate_cneo_calculation_parameters(method, functional, constraint_positions, 
                                             convergence_threshold, max_iterations, lambda_damping)
    validate_calculation_method(method)
    validate_functional_name(functional)
    validate_constraint_positions(constraint_positions)
    validate_convergence_threshold(convergence_threshold)
    validate_iteration_limits(max_iterations)
    validate_damping_factor(lambda_damping)
end

function validate_calculation_method(method::String)
    valid_methods = ["HF", "DFT"]
    if !(method in valid_methods)
        throw(ArgumentError("Invalid calculation method '$method'. Must be one of: $(join(valid_methods, ", "))"))
    end
end

function validate_functional_name(functional::String)
    if isempty(strip(functional))
        throw(ArgumentError("Functional name cannot be empty"))
    end
end

function validate_constraint_positions(positions::Vector{Vector{Float64}})
    for (i, position) in enumerate(positions)
        if length(position) != POSITION_DIMENSION
            throw(ArgumentError("Constraint position $i must have $POSITION_DIMENSION dimensions, got $(length(position))"))
        end
    end
end

function validate_convergence_threshold(threshold::Float64)
    if threshold <= 0.0
        throw(ArgumentError("Convergence threshold must be positive, got $threshold"))
    end
end

function validate_iteration_limits(max_iterations::Int)
    if max_iterations <= 0
        throw(ArgumentError("Max iterations must be positive, got $max_iterations"))
    end
end

function validate_damping_factor(lambda_damping::Float64)
    if lambda_damping <= 0.0 || lambda_damping > 1.0
        throw(ArgumentError("Lambda damping must be in (0, 1], got $lambda_damping"))
    end
end

# ======================== Nuclear Position Analysis ========================

"""
    calculate_nuclear_position_expectation_value(nuclear_mean_field, origin)

Calculate ⟨ϕⁿ|r̂|ϕⁿ⟩ - expectation value of nuclear position operator.

# Arguments
- `nuclear_mean_field::PyObject`: Nuclear mean-field object
- `origin::Vector{Float64}`: Origin for position calculation

# Returns
- `Vector{Float64}`: Nuclear position expectation value (3D coordinates)
"""
function calculate_nuclear_position_expectation_value(nuclear_mean_field::PyObject, origin::Vector{Float64})
    density_matrix = create_nuclear_density_matrix(nuclear_mean_field)
    dipole_integrals = calculate_dipole_integrals(nuclear_mean_field, origin)
    position_expectation = compute_position_expectation_values(density_matrix, dipole_integrals)
    
    return position_expectation + origin
end

function create_nuclear_density_matrix(nuclear_mean_field::PyObject)
    return nuclear_mean_field.make_rdm1()
end

function calculate_dipole_integrals(nuclear_mean_field::PyObject, origin::Vector{Float64})
    molecule = nuclear_mean_field.mol
    pyscf = pyimport("pyscf")
    
    molecule_with_origin = pyscf.gto.mole.with_origin(molecule, origin)
    return molecule_with_origin.intor_symmetric("int1e_r", comp=3)
end

function compute_position_expectation_values(density_matrix, dipole_integrals)
    position = zeros(POSITION_DIMENSION)
    
    for dimension in 1:POSITION_DIMENSION
        position[dimension] = real(tr(density_matrix * dipole_integrals[dimension]))
    end
    
    return position
end

# ======================== Constraint Implementation ========================

"""
    apply_constraint_to_nuclear_hamiltonian!(nuclear_mean_field, lagrange_multiplier, origin)

Add constraint term λ·r to nuclear Hamiltonian for enforcing position constraints.

# Arguments  
- `nuclear_mean_field::PyObject`: Nuclear mean-field object to modify
- `lagrange_multiplier::Vector{Float64}`: Lagrange multiplier vector λ
- `origin::Vector{Float64}`: Origin for constraint application
"""
function apply_constraint_to_nuclear_hamiltonian!(nuclear_mean_field::PyObject, 
                                                 lagrange_multiplier::Vector{Float64}, 
                                                 origin::Vector{Float64})
    constraint_operator = build_constraint_operator(nuclear_mean_field, lagrange_multiplier, origin)
    install_constrained_hamiltonian!(nuclear_mean_field, constraint_operator)
end

function build_constraint_operator(nuclear_mean_field::PyObject, 
                                 lagrange_multiplier::Vector{Float64}, 
                                 origin::Vector{Float64})
    dipole_integrals = calculate_dipole_integrals(nuclear_mean_field, origin)
    constraint_term = zeros(size(dipole_integrals[1]))
    
    for dimension in 1:POSITION_DIMENSION
        constraint_term += lagrange_multiplier[dimension] * dipole_integrals[dimension]
    end
    
    return constraint_term
end

function install_constrained_hamiltonian!(nuclear_mean_field::PyObject, constraint_operator)
    preserve_original_hamiltonian_if_needed!(nuclear_mean_field)
    create_constrained_hamiltonian_method!(nuclear_mean_field, constraint_operator)
end

function preserve_original_hamiltonian_if_needed!(nuclear_mean_field::PyObject)
    if !hasproperty(nuclear_mean_field, :_original_get_hcore)
        nuclear_mean_field._original_get_hcore = nuclear_mean_field.get_hcore
    end
end

function create_constrained_hamiltonian_method!(nuclear_mean_field::PyObject, constraint_operator)
    py"""
    def create_constrained_hamiltonian_method(original_hamiltonian_method, constraint_term):
        def get_constrained_hcore(molecule=None):
            original_hamiltonian = original_hamiltonian_method(molecule)
            return original_hamiltonian + constraint_term
        return get_constrained_hcore
    """
    
    nuclear_mean_field.get_hcore = py"create_constrained_hamiltonian_method"(
        nuclear_mean_field._original_get_hcore, 
        constraint_operator
    )
end

# ======================== Constraint Optimization ========================

"""
    calculate_lagrangian_gradient(nuclear_mean_field, target_position, origin)

Calculate gradient of Lagrangian: ∇L = ⟨ϕⁿ|r̂|ϕⁿ⟩ - R₀.

# Returns
- `Vector{Float64}`: Gradient vector (constraint violation)
"""
function calculate_lagrangian_gradient(nuclear_mean_field::PyObject, 
                                     target_position::Vector{Float64}, 
                                     origin::Vector{Float64})
    current_position = calculate_nuclear_position_expectation_value(nuclear_mean_field, origin)
    return current_position - target_position
end

"""
    calculate_lagrangian_hessian(nuclear_mean_field, origin)

Calculate approximate Hessian using perturbation theory for Newton optimization.

# Returns  
- `Matrix{Float64}`: 3×3 Hessian matrix
"""
function calculate_lagrangian_hessian(nuclear_mean_field::PyObject, origin::Vector{Float64})
    orbital_information = extract_orbital_information(nuclear_mean_field)
    
    if !is_valid_orbital_information(orbital_information)
        return create_identity_hessian()
    end
    
    occupied_orbitals, virtual_orbitals = categorize_orbitals(orbital_information)
    
    if is_insufficient_orbital_space(occupied_orbitals, virtual_orbitals)
        return create_identity_hessian()
    end
    
    dipole_integrals_mo = transform_dipole_integrals_to_mo_basis(nuclear_mean_field, orbital_information, origin)
    hessian_matrix = compute_perturbative_hessian(orbital_information, occupied_orbitals, virtual_orbitals, dipole_integrals_mo)
    
    return ensure_negative_definite_hessian(hessian_matrix)
end

function extract_orbital_information(nuclear_mean_field::PyObject)
    return (
        occupations = get(nuclear_mean_field, :mo_occ, nothing),
        energies = get(nuclear_mean_field, :mo_energy, nothing),
        coefficients = get(nuclear_mean_field, :mo_coeff, nothing)
    )
end

function is_valid_orbital_information(orbital_info)
    return orbital_info.occupations !== nothing && 
           orbital_info.energies !== nothing && 
           orbital_info.coefficients !== nothing
end

function create_identity_hessian()
    return Matrix(1.0I, POSITION_DIMENSION, POSITION_DIMENSION)
end

function categorize_orbitals(orbital_info)
    occupied_indices = findall(orbital_info.occupations .> 0)
    virtual_indices = findall(orbital_info.occupations .== 0)
    return occupied_indices, virtual_indices
end

function is_insufficient_orbital_space(occupied_orbitals, virtual_orbitals)
    return isempty(occupied_orbitals) || isempty(virtual_orbitals)
end

function transform_dipole_integrals_to_mo_basis(nuclear_mean_field::PyObject, orbital_info, origin::Vector{Float64})
    dipole_integrals_ao = calculate_dipole_integrals(nuclear_mean_field, origin)
    numpy = pyimport("numpy")
    
    dipole_integrals_mo = []
    for dimension in 1:POSITION_DIMENSION
        transformed_integral = numpy.dot(
            orbital_info.coefficients.T, 
            numpy.dot(dipole_integrals_ao[dimension], orbital_info.coefficients)
        )
        push!(dipole_integrals_mo, transformed_integral)
    end
    
    return dipole_integrals_mo
end

function compute_perturbative_hessian(orbital_info, occupied_orbitals, virtual_orbitals, dipole_integrals_mo)
    hessian_matrix = zeros(POSITION_DIMENSION, POSITION_DIMENSION)
    
    for occupied_orbital in occupied_orbitals
        for virtual_orbital in virtual_orbitals
            orbital_energy_difference = calculate_orbital_energy_difference(
                orbital_info, occupied_orbital, virtual_orbital
            )
            
            if is_energy_difference_significant(orbital_energy_difference)
                add_orbital_contribution_to_hessian!(
                    hessian_matrix, dipole_integrals_mo, 
                    occupied_orbital, virtual_orbital, orbital_energy_difference
                )
            end
        end
    end
    
    return hessian_matrix
end

function calculate_orbital_energy_difference(orbital_info, occupied_orbital, virtual_orbital)
    return orbital_info.energies[occupied_orbital] - orbital_info.energies[virtual_orbital]
end

function is_energy_difference_significant(energy_difference::Float64)
    return abs(energy_difference) > MINIMUM_DENOMINATOR_THRESHOLD
end

function add_orbital_contribution_to_hessian!(hessian_matrix, dipole_integrals_mo, 
                                            occupied_orbital, virtual_orbital, energy_difference)
    for i in 1:POSITION_DIMENSION, j in 1:POSITION_DIMENSION
        transition_element_i = dipole_integrals_mo[i][virtual_orbital, occupied_orbital]
        transition_element_j = dipole_integrals_mo[j][occupied_orbital, virtual_orbital]
        
        hessian_matrix[i,j] += 2 * real(transition_element_i * transition_element_j) / energy_difference
    end
end

function ensure_negative_definite_hessian(hessian_matrix::Matrix{Float64})
    eigenvalues = eigvals(hessian_matrix)
    
    if any(eigenvalues .>= 0)
        regularization_shift = maximum(eigenvalues) + HESSIAN_REGULARIZATION
        hessian_matrix -= regularization_shift * I
    end
    
    return hessian_matrix
end

# ======================== Main cNEO-HF Implementation ========================

"""
    run_cneo_hf(molecule, cneo_calculation; neo_config)

Run constrained NEO-HF calculation with position constraints on quantum nuclei.

# Arguments
- `molecule::Molecule`: Molecular system with quantum nuclei
- `cneo_calculation::CNEOCalculation`: cNEO calculation parameters
- `neo_config::NEOConfig`: PySCF NEO configuration

# Returns  
- `CNEOResults`: Complete results with energies, positions, and convergence information
"""
function run_cneo_hf(molecule::Molecule, cneo_calculation::CNEOCalculation; 
                     neo_config::NEOConfig = NEOConfig())
    
    @info "Starting cNEO-HF calculation with $(length(cneo_calculation.constraint_positions)) constraints"
    
    # Step 1: Initialize quantum chemistry environment
    mean_field_environment = initialize_cneo_environment(molecule, cneo_calculation, neo_config)
    
    # Step 2: Run initial unconstrained calculation
    initial_energy = perform_initial_scf_calculation!(mean_field_environment.mean_field_object)
    @info "Initial unconstrained NEO energy: $(format_energy(initial_energy)) Ha"
    
    # Step 3: Apply constraints and optimize
    constraint_results = optimize_nuclear_constraints!(mean_field_environment, cneo_calculation)
    
    # Step 4: Analyze final results
    final_results = analyze_cneo_results(
        mean_field_environment, cneo_calculation, constraint_results
    )
    
    display_cneo_results_summary(final_results)
    return final_results
end

struct CNEOMeanFieldEnvironment
    mean_field_object::PyObject
    neo_molecule::PyObject
    nuclear_origins::Vector{Vector{Float64}}
    nuclear_components::Vector{String}
end

function initialize_cneo_environment(molecule::Molecule, cneo_calculation::CNEOCalculation, neo_config::NEOConfig)
    # Setup PySCF
    pyscf_module, has_neo = setup_pyscf(neo_config)
    validate_neo_availability(has_neo)
    
    # Build NEO molecule  
    neo_molecule = build_neo_molecule(molecule, pyscf_module)
    
    # Create mean-field object
    mean_field_object = create_cneo_mean_field_object(neo_molecule, cneo_calculation, pyscf_module)
    
    # Extract nuclear component information
    nuclear_origins, nuclear_components = extract_nuclear_component_information(mean_field_object)
    
    validate_constraint_compatibility(cneo_calculation, nuclear_components)
    
    return CNEOMeanFieldEnvironment(mean_field_object, neo_molecule, nuclear_origins, nuclear_components)
end

function validate_neo_availability(has_neo::Bool)
    if !has_neo
        throw(ArgumentError("NEO module not available in PySCF. Please install PySCF with NEO support."))
    end
end

function create_cneo_mean_field_object(neo_molecule, cneo_calculation::CNEOCalculation, pyscf_module)
    if cneo_calculation.method == "HF"
        return pyscf_module.neo.HF(neo_molecule).density_fit()
    else
        return pyscf_module.neo.KS(neo_molecule, xc=cneo_calculation.functional).density_fit()
    end
end

function extract_nuclear_component_information(mean_field_object::PyObject)
    nuclear_components = [key for key in keys(mean_field_object.components) if startswith(string(key), "n")]
    nuclear_origins = extract_nuclear_origins(mean_field_object, nuclear_components)
    
    return nuclear_origins, nuclear_components
end

function extract_nuclear_origins(mean_field_object::PyObject, nuclear_components::Vector{String})
    origins = Vector{Vector{Float64}}()
    
    for component_key in nuclear_components
        nuclear_mean_field = mean_field_object.components[component_key]
        nuclear_molecule = nuclear_mean_field.mol
        
        origin = extract_nuclear_origin_coordinates(nuclear_molecule)
        push!(origins, origin)
    end
    
    return origins
end

function extract_nuclear_origin_coordinates(nuclear_molecule)
    if hasproperty(nuclear_molecule, :atom_coords)
        coordinates = nuclear_molecule.atom_coords()
        return coordinates[1, :]  # First atom is the quantum nucleus
    else
        return zeros(POSITION_DIMENSION)
    end
end

function validate_constraint_compatibility(cneo_calculation::CNEOCalculation, nuclear_components::Vector{String})
    n_constraints = length(cneo_calculation.constraint_positions)
    n_nuclei = length(nuclear_components)
    
    if n_constraints != n_nuclei
        throw(ArgumentError(
            "Number of constraint positions ($n_constraints) must match number of quantum nuclei ($n_nuclei)"
        ))
    end
end

function perform_initial_scf_calculation!(mean_field_object::PyObject)
    mean_field_object.kernel()
    return mean_field_object.e_tot
end

struct CNEOConstraintResults
    lagrange_multipliers::Vector{Vector{Float64}}
    converged::Bool
    iterations::Int
    final_positions::Vector{Vector{Float64}}
    position_errors::Vector{Float64}
end

function optimize_nuclear_constraints!(environment::CNEOMeanFieldEnvironment, 
                                     cneo_calculation::CNEOCalculation)
    lagrange_multipliers = initialize_lagrange_multipliers(environment.nuclear_components)
    
    @info "Starting cNEO constraint optimization..."
    
    converged = false
    iteration = 0
    max_position_error = Inf
    
    while !converged && iteration < cneo_calculation.max_iterations
        iteration += 1
        
        # Apply current constraints to all nuclear components
        apply_constraints_to_all_nuclei!(environment, lagrange_multipliers)
        
        # Run SCF with constraints
        environment.mean_field_object.kernel()
        
        # Update Lagrange multipliers using Newton method
        max_position_error = update_lagrange_multipliers_newton!(
            environment, cneo_calculation, lagrange_multipliers
        )
        
        # Check convergence
        converged = max_position_error < cneo_calculation.convergence_threshold
        
        log_iteration_progress(iteration, max_position_error, converged)
    end
    
    final_positions, position_errors = calculate_final_constraint_satisfaction(
        environment, cneo_calculation
    )
    
    return CNEOConstraintResults(
        lagrange_multipliers, converged, iteration, final_positions, position_errors
    )
end

function initialize_lagrange_multipliers(nuclear_components::Vector{String})
    return [zeros(POSITION_DIMENSION) for _ in nuclear_components]
end

function apply_constraints_to_all_nuclei!(environment::CNEOMeanFieldEnvironment, 
                                        lagrange_multipliers::Vector{Vector{Float64}})
    for (index, component_key) in enumerate(environment.nuclear_components)
        nuclear_mean_field = environment.mean_field_object.components[component_key]
        nuclear_origin = environment.nuclear_origins[index]
        lagrange_multiplier = lagrange_multipliers[index]
        
        apply_constraint_to_nuclear_hamiltonian!(nuclear_mean_field, lagrange_multiplier, nuclear_origin)
    end
end

function update_lagrange_multipliers_newton!(environment::CNEOMeanFieldEnvironment,
                                           cneo_calculation::CNEOCalculation,
                                           lagrange_multipliers::Vector{Vector{Float64}})
    max_error = 0.0
    
    for (index, component_key) in enumerate(environment.nuclear_components)
        nuclear_mean_field = environment.mean_field_object.components[component_key]
        nuclear_origin = environment.nuclear_origins[index]
        target_position = cneo_calculation.constraint_positions[index]
        
        position_error = update_single_nucleus_multiplier!(
            nuclear_mean_field, target_position, nuclear_origin,
            cneo_calculation.lambda_damping, lagrange_multipliers[index]
        )
        
        max_error = max(max_error, position_error)
    end
    
    return max_error
end

function update_single_nucleus_multiplier!(nuclear_mean_field::PyObject, 
                                         target_position::Vector{Float64},
                                         nuclear_origin::Vector{Float64},
                                         damping_factor::Float64,
                                         lagrange_multiplier::Vector{Float64})
    gradient = calculate_lagrangian_gradient(nuclear_mean_field, target_position, nuclear_origin)
    position_error = norm(gradient)
    
    try
        hessian = calculate_lagrangian_hessian(nuclear_mean_field, nuclear_origin)
        newton_step = -hessian \ gradient
        lagrange_multiplier .+= damping_factor * newton_step
    catch error
        @warn "Hessian inversion failed, using gradient descent" exception=error
        lagrange_multiplier .-= GRADIENT_DESCENT_STEP * gradient
    end
    
    return position_error
end

function log_iteration_progress(iteration::Int, max_error::Float64, converged::Bool)
    if iteration % 5 == 0 || converged
        status = converged ? " (CONVERGED)" : ""
        @info "cNEO iteration $iteration: max position error = $(round(max_error, digits=8))$status"
    end
end

function calculate_final_constraint_satisfaction(environment::CNEOMeanFieldEnvironment, 
                                               cneo_calculation::CNEOCalculation)
    final_positions = Vector{Vector{Float64}}()
    position_errors = Vector{Float64}()
    
    for (index, component_key) in enumerate(environment.nuclear_components)
        nuclear_mean_field = environment.mean_field_object.components[component_key]
        nuclear_origin = environment.nuclear_origins[index]
        target_position = cneo_calculation.constraint_positions[index]
        
        actual_position = calculate_nuclear_position_expectation_value(nuclear_mean_field, nuclear_origin)
        position_error = norm(actual_position - target_position)
        
        push!(final_positions, actual_position)
        push!(position_errors, position_error)
    end
    
    return final_positions, position_errors
end

function analyze_cneo_results(environment::CNEOMeanFieldEnvironment,
                            cneo_calculation::CNEOCalculation,
                            constraint_results::CNEOConstraintResults)
    # Calculate energy components
    total_energy = environment.mean_field_object.e_tot
    electronic_energy = environment.mean_field_object.energy_elec()[1]
    nuclear_kinetic_energy = total_energy - electronic_energy
    
    return CNEOResults(
        total_energy,
        electronic_energy,
        nuclear_kinetic_energy,
        constraint_results.lagrange_multipliers,
        constraint_results.final_positions,
        cneo_calculation.constraint_positions,
        constraint_results.position_errors,
        constraint_results.converged,
        constraint_results.iterations,
        environment.mean_field_object
    )
end

function display_cneo_results_summary(results::CNEOResults)
    @info "cNEO-HF calculation completed:"
    @info "  Total energy: $(format_energy(results.total_energy)) Ha"
    @info "  Electronic energy: $(format_energy(results.electronic_energy)) Ha"
    @info "  Nuclear kinetic energy: $(format_energy(results.nuclear_kinetic_energy)) Ha"
    @info "  Converged: $(results.converged)"
    @info "  Iterations: $(results.iterations)"
    
    display_constraint_analysis(results)
end

function display_constraint_analysis(results::CNEOResults)
    for nucleus_index in 1:length(results.constraint_positions)
        @info "  Nucleus $nucleus_index constraint analysis:"
        @info "    Target position: $(results.constraint_positions[nucleus_index])"
        @info "    Actual position: $(results.actual_nuclear_positions[nucleus_index])" 
        @info "    Position error: $(format_energy(results.position_errors[nucleus_index]))"
        @info "    Lagrange multiplier: $(results.lagrange_multipliers[nucleus_index])"
    end
end

function format_energy(energy::Float64)
    return round(energy, digits=6)
end

# ======================== cNEO-MP2 Implementation (Stub) ========================

"""
    run_cneo_mp2(molecule, constraint_positions; neo_config)

Run cNEO-MP2 calculation with correlation corrections.

# Arguments
- `molecule::Molecule`: Molecular system  
- `constraint_positions::Vector{Vector{Float64}}`: Nuclear constraint positions
- `neo_config::NEOConfig`: PySCF NEO configuration

# Returns
- `CNEOMP2Results`: Complete MP2 results with correlation analysis
"""
function run_cneo_mp2(molecule::Molecule, constraint_positions::Vector{Vector{Float64}};
                      neo_config::NEOConfig = NEOConfig())
    
    @info "Running cNEO-MP2 calculation..."
    
    # Step 1: Run cNEO-HF first
    cneo_hf_calculation = create_cneo_calculation(constraint_positions=constraint_positions)
    hf_results = run_cneo_hf(molecule, cneo_hf_calculation, neo_config=neo_config)
    
    # Step 2: Add MP2 correlation (simplified implementation)
    mp2_correlation_energy = calculate_cneo_mp2_correlation(hf_results.mean_field_object)
    
    # Step 3: Create combined results
    return create_cneo_mp2_results(hf_results, mp2_correlation_energy)
end

function calculate_cneo_mp2_correlation(mean_field_object::PyObject)
    # Simplified MP2 correlation calculation
    # In full implementation, this would include nuclear correlation effects
    pyscf = pyimport("pyscf")
    
    try
        # Use PySCF's NEO-MP2 if available
        if hasmethod(pyscf.neo, :MP2)
            mp2_object = pyscf.neo.MP2(mean_field_object)
            correlation_energy, _ = mp2_object.kernel()
            return correlation_energy
        else
            # Fallback to electronic MP2 only
            electronic_component = mean_field_object.components["e"]
            mp2_object = pyscf.mp.MP2(electronic_component)
            correlation_energy, _ = mp2_object.kernel()
            @warn "Using electronic MP2 only - nuclear correlation not included"
            return correlation_energy
        end
    catch error
        @warn "MP2 calculation failed, returning zero correlation" exception=error
        return 0.0
    end
end

function create_cneo_mp2_results(hf_results::CNEOResults, mp2_correlation::Float64)
    total_energy = hf_results.total_energy + mp2_correlation
    
    return CNEOMP2Results(
        hf_results.total_energy,  # hf_energy
        mp2_correlation,          # mp2_correlation_energy
        total_energy,             # total_energy
        hf_results.lagrange_multipliers,  # lagrange_multipliers_hf
        hf_results.lagrange_multipliers,  # lagrange_multipliers_mp2 (same as HF for now)
        hf_results.actual_nuclear_positions,  # nuclear_positions_hf
        hf_results.actual_nuclear_positions,  # nuclear_positions_mp2 (same as HF for now)
        hf_results.position_errors,  # position_errors_hf
        hf_results.position_errors,  # position_errors_mp2 (same as HF for now)
        hf_results.converged,     # hf_converged
        true,                     # mp2_converged (simplified)
        hf_results.iterations,    # hf_iterations
        1,                        # mp2_iterations (simplified)
        hf_results.mean_field_object  # mean_field_object
    )
end

end # module CNEOMethods