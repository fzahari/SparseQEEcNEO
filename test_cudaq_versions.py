"""
Test suite for CUDA-Q versions
================================
Verifies that CUDA-Q implementations produce correct results.

These tests work without CUDA-Q installed (fall back to NumPy).
"""
import numpy as np
import pytest

# Test richerme_ion_analog_cudaq
def test_richerme_cudaq_basic():
    """Test basic gate synthesis with CUDA-Q version"""
    from richerme_ion_analog_cudaq import (
        n_body_string,
        n_body_string_arbitrary,
        UMQ,
        target_pauli_string_unitary,
        unitary_distance
    )

    # Test 1: ZXX gate synthesis
    t = 0.5
    U_synth = n_body_string(['Z', 'X', 'X'], t)
    U_target = target_pauli_string_unitary('ZXX', t)
    error = unitary_distance(U_target, U_synth)
    assert error < 1e-12, f"ZXX synthesis error: {error:.2e}"
    print(f"✓ ZXX synthesis: error = {error:.2e}")

    # Test 2: Arbitrary Pauli pattern
    U_zyz = n_body_string_arbitrary(['Z', 'Y', 'Z'], 0.3)
    U_target = target_pauli_string_unitary('ZYZ', 0.3)
    error = unitary_distance(U_target, U_zyz)
    assert error < 1e-12, f"ZYZ synthesis error: {error:.2e}"
    print(f"✓ ZYZ synthesis: error = {error:.2e}")

    # Test 3: UMQ unitarity
    n = 3
    chi = np.pi / 4
    U_umq = UMQ(n, chi)
    unitarity = np.linalg.norm(U_umq @ U_umq.conj().T - np.eye(2**n))
    assert unitarity < 1e-12, f"UMQ not unitary: {unitarity:.2e}"
    print(f"✓ UMQ unitarity: error = {unitarity:.2e}")


def test_richerme_cudaq_accessibility():
    """Test interaction graph accessibility functions"""
    from richerme_ion_analog_cudaq import (
        is_accessible,
        get_mode_weights_if_accessible,
        compute_sinusoidal_modes
    )

    N = 5
    B = compute_sinusoidal_modes(N)

    # Test orthonormality
    gram = B.T @ B
    ortho_error = np.linalg.norm(gram - np.eye(N))
    assert ortho_error < 1e-12, f"Modes not orthonormal: {ortho_error:.2e}"
    print(f"✓ Sinusoidal modes orthonormal: error = {ortho_error:.2e}")

    # Test all-to-all accessibility
    J_all_to_all = np.ones((N, N)) - np.eye(N)
    accessible = is_accessible(J_all_to_all, B)
    assert accessible, "All-to-all should be accessible"
    print(f"✓ All-to-all interaction is accessible")

    weights = get_mode_weights_if_accessible(J_all_to_all, B)
    assert weights is not None, "Should return mode weights"
    print(f"✓ Mode weights extracted: {weights}")


def test_h2_cudaq_basic():
    """Test H2 simulator CUDA-Q version"""
    from rich_sim_h2_cudaq import H2TimeEvolutionSimulator, TrappedIonSimulator

    ion_sim = TrappedIonSimulator(N=4)
    h2_sim = H2TimeEvolutionSimulator(ion_sim, use_hardware_gates=False)

    # Test Hamiltonian construction
    r = 0.74
    pauli_terms = h2_sim.build_h2_hamiltonian_pauli(r)
    assert len(pauli_terms) == 15, f"Expected 15 terms, got {len(pauli_terms)}"
    print(f"✓ H2 Hamiltonian has {len(pauli_terms)} terms")

    # Test Hamiltonian matrix
    H = h2_sim.pauli_to_matrix(pauli_terms)
    assert H.shape == (16, 16), f"Expected 16×16, got {H.shape}"
    hermitian_error = np.linalg.norm(H - H.conj().T)
    assert hermitian_error < 1e-12, f"Hamiltonian not Hermitian: {hermitian_error:.2e}"
    print(f"✓ H2 Hamiltonian is Hermitian: error = {hermitian_error:.2e}")

    # Test imaginary-time evolution (short run)
    imag_results = h2_sim.imaginary_time_evolution(r, beta_max=10.0, n_steps=100)
    assert 'ground_energy' in imag_results, "Missing ground_energy"
    assert 'energies' in imag_results, "Missing energies"
    assert len(imag_results['energies']) == 100, "Wrong number of steps"
    print(f"✓ Imaginary-time evolution completed")
    print(f"  Ground energy: {imag_results['ground_energy']:.6f} H")


def test_h2_cudaq_ansatz():
    """Test H2 hardware-efficient ansatz"""
    from rich_sim_h2_cudaq import H2TimeEvolutionSimulator, TrappedIonSimulator

    ion_sim = TrappedIonSimulator(N=4)
    h2_sim = H2TimeEvolutionSimulator(ion_sim, use_hardware_gates=False)

    # Test ansatz with XX-only gates, 2 layers, no nonlocal
    n_layers = 2
    use_full_gates = False
    use_nonlocal = False

    # Calculate expected parameters
    local_pairs = 3
    params_per_layer = 4 * 2 + local_pairs  # 4 Ry + 4 Rz + 3 XX
    n_params = params_per_layer * n_layers
    params = np.random.randn(n_params) * 0.1

    psi = h2_sim.hardware_efficient_ansatz(params, use_full_gates=use_full_gates,
                                            n_layers=n_layers, use_nonlocal=use_nonlocal)

    assert psi.shape == (16,), f"Expected shape (16,), got {psi.shape}"
    norm = np.linalg.norm(psi)
    assert abs(norm - 1.0) < 1e-12, f"State not normalized: |psi| = {norm}"
    print(f"✓ Ansatz generates normalized state: |psi| = {norm:.10f}")


def test_h2o_cudaq_basic():
    """Test H2O simulator CUDA-Q version"""
    from rich_sim_h2o_cudaq import H2O_QEE_Simulator, TrappedIonSimulator

    ion_sim = TrappedIonSimulator(N=10)
    h2o_sim = H2O_QEE_Simulator(ion_sim, use_hardware_gates=False)

    # Test QEE mapping
    assert h2o_sim.n_qubits_original == 14, "Should compress from 14 qubits"
    assert h2o_sim.n_qubits_qee == 10, "Should compress to 10 qubits"
    assert len(h2o_sim.qee_map['valid_states']) <= 1024, "Too many states"
    print(f"✓ QEE mapping: 14 → 10 qubits ({len(h2o_sim.qee_map['valid_states'])} states)")

    # Test Hamiltonian terms
    assert len(h2o_sim.hamiltonian_terms) > 0, "Should have Hamiltonian terms"
    print(f"✓ H2O Hamiltonian has {len(h2o_sim.hamiltonian_terms)} terms")

    # Test term grouping
    assert 'diagonal' in h2o_sim.grouped_terms, "Missing diagonal group"
    assert 'hopping' in h2o_sim.grouped_terms, "Missing hopping group"
    assert 'execution_order' in h2o_sim.grouped_terms, "Missing execution order"
    print(f"✓ Term grouping: {len(h2o_sim.grouped_terms['execution_order'])} operations")

    # Test evolution operator construction (just verify it builds)
    U = h2o_sim.build_evolution_operator(dt=0.01)
    assert U.shape == (1024, 1024), f"Expected 1024×1024, got {U.shape}"
    unitarity = np.linalg.norm(U @ U.conj().T - np.eye(1024))
    assert unitarity < 1e-10, f"Evolution operator not unitary: {unitarity:.2e}"
    print(f"✓ Evolution operator is unitary: error = {unitarity:.2e}")


def run_all_tests():
    """Run all tests"""
    print("=" * 70)
    print("Testing CUDA-Q Implementations (No CUDA-Q Required)")
    print("=" * 70)
    print()

    print("Test 1: Richerme CUDA-Q Basic Functionality")
    print("-" * 70)
    test_richerme_cudaq_basic()
    print()

    print("Test 2: Richerme CUDA-Q Accessibility")
    print("-" * 70)
    test_richerme_cudaq_accessibility()
    print()

    print("Test 3: H2 CUDA-Q Basic Functionality")
    print("-" * 70)
    test_h2_cudaq_basic()
    print()

    print("Test 4: H2 CUDA-Q Ansatz")
    print("-" * 70)
    test_h2_cudaq_ansatz()
    print()

    print("Test 5: H2O CUDA-Q Basic Functionality")
    print("-" * 70)
    test_h2o_cudaq_basic()
    print()

    print("=" * 70)
    print("All CUDA-Q tests passed!")
    print("=" * 70)


if __name__ == "__main__":
    run_all_tests()
