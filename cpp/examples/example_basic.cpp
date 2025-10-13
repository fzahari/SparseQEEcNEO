/**
 * @file example_basic.cpp
 * @brief Basic usage example for Richerme Ion Analog C++ library
 */

#include "richerme_ion_analog.h"
#include <iostream>
#include <iomanip>

using namespace richerme;

int main() {
    std::cout << "======================================================================\n";
    std::cout << "Richerme Ion Analog C++ Library - Basic Examples\n";
    std::cout << "======================================================================\n\n";

    // Example 1: Basic gate synthesis
    std::cout << "Example 1: Z⊗X⊗X Gate Synthesis\n";
    std::cout << "----------------------------------------------------------------------\n";

    double t = 0.5;
    std::vector<char> axes = {'Z', 'X', 'X'};

    Matrix U_synth = GateSynthesis::n_body_string(axes, t);
    Matrix U_target = GateSynthesis::target_pauli_string_unitary("ZXX", t);

    double error = GateSynthesis::unitary_distance(U_target, U_synth);

    std::cout << "  Evolution time: t = " << t << "\n";
    std::cout << "  Unitary shape: " << U_synth.rows() << "×" << U_synth.cols() << "\n";
    std::cout << "  Fidelity error: " << std::scientific << std::setprecision(2)
              << error << "\n";

    if (error < 1e-12) {
        std::cout << "  ✓ Synthesis achieved machine precision!\n";
    }
    std::cout << "\n";

    // Example 2: Arbitrary Pauli pattern
    std::cout << "Example 2: Arbitrary Z⊗Y⊗Z Pattern\n";
    std::cout << "----------------------------------------------------------------------\n";

    t = 0.3;
    axes = {'Z', 'Y', 'Z'};

    U_synth = GateSynthesis::n_body_string_arbitrary(axes, t);
    U_target = GateSynthesis::target_pauli_string_unitary("ZYZ", t);

    error = GateSynthesis::unitary_distance(U_target, U_synth);

    std::cout << "  Evolution time: t = " << t << "\n";
    std::cout << "  Fidelity error: " << std::scientific << error << "\n";

    if (error < 1e-12) {
        std::cout << "  ✓ Arbitrary pattern synthesis successful!\n";
    }
    std::cout << "\n";

    // Example 3: UMQ global entangling gate
    std::cout << "Example 3: UMQ Global Entangling Gate\n";
    std::cout << "----------------------------------------------------------------------\n";

    int n = 3;
    double chi = M_PI / 4.0;

    Matrix U_umq = GateSynthesis::UMQ(n, chi);

    // Check unitarity: U†U = I
    int dim = 1 << n;  // 2^n
    Matrix identity = Matrix::Identity(dim, dim);
    Matrix product = U_umq.adjoint() * U_umq;
    double unitarity_error = (product - identity).norm();

    std::cout << "  System size: " << n << " qubits\n";
    std::cout << "  Interaction strength: χ = π/4\n";
    std::cout << "  Unitarity check: " << std::scientific << unitarity_error << "\n";

    if (unitarity_error < 1e-12) {
        std::cout << "  ✓ UMQ is unitary!\n";
    }
    std::cout << "\n";

    // Example 4: Interaction graph accessibility
    std::cout << "Example 4: Interaction Graph Accessibility\n";
    std::cout << "----------------------------------------------------------------------\n";

    int N = 5;
    RealMatrix B = InteractionEngineering::compute_sinusoidal_modes(N);

    // All-to-all interaction
    RealMatrix J_all_to_all = RealMatrix::Ones(N, N) - RealMatrix::Identity(N, N);

    bool accessible = InteractionEngineering::is_accessible(J_all_to_all, B);

    std::cout << "  " << N << " ions, all-to-all interaction\n";
    std::cout << "  Accessible: " << (accessible ? "Yes" : "No") << "\n";

    if (accessible) {
        RealVector weights = InteractionEngineering::get_mode_weights_if_accessible(
            J_all_to_all, B);
        std::cout << "  Mode weights: [";
        for (int i = 0; i < weights.size(); ++i) {
            std::cout << std::fixed << std::setprecision(1) << weights(i);
            if (i < weights.size() - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "\n";

    // Example 5: Hardware specifications
    std::cout << "Example 5: 171Yb+ Hardware Specifications\n";
    std::cout << "----------------------------------------------------------------------\n";

    IonTrapHardware hw;

    std::cout << "  Hyperfine splitting: " << hw.hyperfine_splitting / 1e9
              << " GHz\n";
    std::cout << "  T2 coherence time: " << hw.T2_coherence << " s\n";
    std::cout << "  Single-qubit fidelity: "
              << std::fixed << std::setprecision(1) << hw.single_qubit_fidelity * 100
              << "%\n";
    std::cout << "  Two-qubit fidelity: " << hw.two_qubit_fidelity * 100 << "%\n";
    std::cout << "  Recoil frequency: "
              << std::scientific << std::setprecision(2) << hw.recoil_frequency() / 1e3
              << " kHz\n";
    std::cout << "\n";

    std::cout << "======================================================================\n";
    std::cout << "All examples completed successfully!\n";
    std::cout << "======================================================================\n";

    return 0;
}
