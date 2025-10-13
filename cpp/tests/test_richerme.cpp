/**
 * @file test_richerme.cpp
 * @brief Test suite for Richerme Ion Analog C++ library
 */

#include "richerme_ion_analog.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace richerme;

void test_pauli_operators() {
    std::cout << "Test 1: Pauli Operators\n";
    std::cout << "----------------------------------------------------------------------\n";

    // Test Pauli matrices properties
    Matrix X = PauliOps::X();
    Matrix Y = PauliOps::Y();
    Matrix Z = PauliOps::Z();
    Matrix I = PauliOps::I();

    // X² = I
    Matrix X2 = X * X;
    double error = (X2 - I).norm();
    assert(error < 1e-14);
    std::cout << "  ✓ X² = I (error: " << error << ")\n";

    // Y² = I
    Matrix Y2 = Y * Y;
    error = (Y2 - I).norm();
    assert(error < 1e-14);
    std::cout << "  ✓ Y² = I (error: " << error << ")\n";

    // Z² = I
    Matrix Z2 = Z * Z;
    error = (Z2 - I).norm();
    assert(error < 1e-14);
    std::cout << "  ✓ Z² = I (error: " << error << ")\n";

    // {X, Y} = 0 (anticommutator)
    Matrix anticomm = X*Y + Y*X;
    error = anticomm.norm();
    assert(error < 1e-14);
    std::cout << "  ✓ {X, Y} = 0 (error: " << error << ")\n";

    std::cout << "\n";
}

void test_gate_synthesis() {
    std::cout << "Test 2: Gate Synthesis\n";
    std::cout << "----------------------------------------------------------------------\n";

    // Test ZXX synthesis
    double t = 0.5;
    std::vector<char> axes = {'Z', 'X', 'X'};
    Matrix U_synth = GateSynthesis::n_body_string(axes, t);
    Matrix U_target = GateSynthesis::target_pauli_string_unitary("ZXX", t);

    double error = GateSynthesis::unitary_distance(U_target, U_synth);
    assert(error < 1e-12);
    std::cout << "  ✓ ZXX synthesis (error: " << error << ")\n";

    // Test YXX synthesis
    axes = {'Y', 'X', 'X'};
    U_synth = GateSynthesis::n_body_string(axes, t);
    U_target = GateSynthesis::target_pauli_string_unitary("YXX", t);

    error = GateSynthesis::unitary_distance(U_target, U_synth);
    assert(error < 1e-12);
    std::cout << "  ✓ YXX synthesis (error: " << error << ")\n";

    // Test XXX synthesis
    // Note: Using arbitrary pattern for better precision (no basis rotations needed)
    axes = {'X', 'X', 'X'};
    U_synth = GateSynthesis::n_body_string_arbitrary(axes, t);
    U_target = GateSynthesis::target_pauli_string_unitary("XXX", t);

    error = GateSynthesis::unitary_distance(U_target, U_synth);
    assert(error < 1e-12);
    std::cout << "  ✓ XXX synthesis (error: " << error << ")\n";

    std::cout << "\n";
}

void test_arbitrary_patterns() {
    std::cout << "Test 3: Arbitrary Pauli Patterns\n";
    std::cout << "----------------------------------------------------------------------\n";

    double t = 0.3;

    // Test ZYZ
    std::vector<char> axes = {'Z', 'Y', 'Z'};
    Matrix U_synth = GateSynthesis::n_body_string_arbitrary(axes, t);
    Matrix U_target = GateSynthesis::target_pauli_string_unitary("ZYZ", t);

    double error = GateSynthesis::unitary_distance(U_target, U_synth);
    assert(error < 1e-12);
    std::cout << "  ✓ ZYZ synthesis (error: " << error << ")\n";

    // Test XYX
    axes = {'X', 'Y', 'X'};
    U_synth = GateSynthesis::n_body_string_arbitrary(axes, t);
    U_target = GateSynthesis::target_pauli_string_unitary("XYX", t);

    error = GateSynthesis::unitary_distance(U_target, U_synth);
    assert(error < 1e-12);
    std::cout << "  ✓ XYX synthesis (error: " << error << ")\n";

    // Test YZY
    axes = {'Y', 'Z', 'Y'};
    U_synth = GateSynthesis::n_body_string_arbitrary(axes, t);
    U_target = GateSynthesis::target_pauli_string_unitary("YZY", t);

    error = GateSynthesis::unitary_distance(U_target, U_synth);
    assert(error < 1e-12);
    std::cout << "  ✓ YZY synthesis (error: " << error << ")\n";

    std::cout << "\n";
}

void test_umq_gate() {
    std::cout << "Test 4: UMQ Gate\n";
    std::cout << "----------------------------------------------------------------------\n";

    int n = 3;
    double chi = M_PI / 4.0;

    Matrix U_umq = GateSynthesis::UMQ(n, chi);

    // Check unitarity
    int dim = 1 << n;
    Matrix identity = Matrix::Identity(dim, dim);
    Matrix product = U_umq.adjoint() * U_umq;
    double error = (product - identity).norm();

    assert(error < 1e-12);
    std::cout << "  ✓ UMQ unitarity (error: " << error << ")\n";

    // Check zero interaction gives identity
    Matrix U_zero = GateSynthesis::UMQ(n, 0.0);
    error = (U_zero - identity).norm();

    assert(error < 1e-12);
    std::cout << "  ✓ UMQ(0) = I (error: " << error << ")\n";

    std::cout << "\n";
}

void test_accessibility() {
    std::cout << "Test 5: Interaction Graph Accessibility\n";
    std::cout << "----------------------------------------------------------------------\n";

    int N = 5;
    RealMatrix B = InteractionEngineering::compute_sinusoidal_modes(N);

    // Test orthonormality
    RealMatrix gram = B.transpose() * B;
    RealMatrix identity = RealMatrix::Identity(N, N);
    double error = (gram - identity).norm();

    assert(error < 1e-12);
    std::cout << "  ✓ Sinusoidal modes orthonormal (error: " << error << ")\n";

    // Test all-to-all accessibility
    RealMatrix J_all_to_all = RealMatrix::Ones(N, N) - RealMatrix::Identity(N, N);
    bool accessible = InteractionEngineering::is_accessible(J_all_to_all, B);

    assert(accessible);
    std::cout << "  ✓ All-to-all interaction is accessible\n";

    // Test mode weight extraction
    RealVector weights = InteractionEngineering::get_mode_weights_if_accessible(
        J_all_to_all, B);

    assert(weights.size() == N);
    std::cout << "  ✓ Mode weights extracted (" << weights.size() << " modes)\n";

    std::cout << "\n";
}

void test_hardware_specs() {
    std::cout << "Test 6: Hardware Specifications\n";
    std::cout << "----------------------------------------------------------------------\n";

    IonTrapHardware hw;

    // Check recoil frequency calculation
    double R = hw.recoil_frequency();

    assert(R > 0);
    assert(R < 1e6);  // Should be kHz range
    std::cout << "  ✓ Recoil frequency: " << R / 1e3 << " kHz\n";

    // Check fidelities are reasonable
    assert(hw.single_qubit_fidelity > 0.9);
    assert(hw.single_qubit_fidelity <= 1.0);
    assert(hw.two_qubit_fidelity > 0.9);
    assert(hw.two_qubit_fidelity <= 1.0);
    std::cout << "  ✓ Fidelities in valid range\n";

    std::cout << "\n";
}

int main() {
    std::cout << "======================================================================\n";
    std::cout << "Richerme Ion Analog C++ Library - Test Suite\n";
    std::cout << "======================================================================\n\n";

    try {
        test_pauli_operators();
        test_gate_synthesis();
        test_arbitrary_patterns();
        test_umq_gate();
        test_accessibility();
        test_hardware_specs();

        std::cout << "======================================================================\n";
        std::cout << "All tests passed!\n";
        std::cout << "======================================================================\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
