/**
 * @file richerme_ion_analog.h
 * @brief C++ header for Richerme Ion Analog quantum gate synthesis
 *
 * Quantum gate synthesis for trapped-ion analog quantum computers.
 * Uses CUDA-Q for GPU-accelerated operator algebra (no circuits).
 *
 * Based on:
 * - Richerme et al. (2025) - Multi-mode global driving
 * - Kyprianidis et al. (2024) - Interaction graph engineering
 *
 * @author Federico Zahariev
 * @date 2025-10-12
 */

#ifndef RICHERME_ION_ANALOG_H
#define RICHERME_ION_ANALOG_H

#include <complex>
#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>

namespace richerme {

// Type aliases for convenience
using Complex = std::complex<double>;
using Matrix = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;
using RealMatrix = Eigen::MatrixXd;
using RealVector = Eigen::VectorXd;

/**
 * @brief Hardware specifications for 171Yb+ trapped-ion system
 */
struct IonTrapHardware {
    double hyperfine_splitting = 12.6e9;        // Hz (12.6 GHz)
    double T2_coherence = 1.0;                  // seconds
    double T1_coherence = std::numeric_limits<double>::infinity();
    double single_qubit_fidelity = 0.998;       // 99.8%
    double two_qubit_fidelity = 0.97;           // 97%
    double omega_z = 2.0 * M_PI * 0.1e6;        // Hz (100 kHz axial)
    double omega_x = 2.0 * M_PI * 5.0e6;        // Hz (5 MHz radial)
    double omega_y = 2.0 * M_PI * 5.0e6;        // Hz (5 MHz radial)
    double wavelength = 369.5e-9;               // m (369.5 nm)
    double mass = 171.0 * 1.66054e-27;          // kg (171 amu)

    /**
     * @brief Calculate recoil frequency R = ℏ(Δk)²/(2m)
     * @return Recoil frequency in Hz
     */
    double recoil_frequency() const;
};

/**
 * @brief Pauli matrices and basic quantum operations
 */
class PauliOps {
public:
    /**
     * @brief Get Pauli matrix
     * @param pauli_char 'I', 'X', 'Y', or 'Z'
     * @return 2×2 Pauli matrix
     */
    static Matrix pauli(char pauli_char);

    /**
     * @brief Identity matrix
     */
    static Matrix I();

    /**
     * @brief Pauli X matrix
     */
    static Matrix X();

    /**
     * @brief Pauli Y matrix
     */
    static Matrix Y();

    /**
     * @brief Pauli Z matrix
     */
    static Matrix Z();
};

/**
 * @brief Quantum gate synthesis utilities
 */
class GateSynthesis {
public:
    /**
     * @brief Kronecker product of multiple matrices
     * @param ops Vector of matrices to tensor together
     * @return Kronecker product
     */
    static Matrix kronN(const std::vector<Matrix>& ops);

    /**
     * @brief Single Pauli operator in n-qubit space
     * @param n Number of qubits
     * @param pauli_char Pauli operator ('X', 'Y', 'Z')
     * @param qubit_idx Which qubit (0-indexed)
     * @return n-qubit operator
     */
    static Matrix pauli_operator(int n, char pauli_char, int qubit_idx);

    /**
     * @brief Sum of Pauli operators: P₀ + P₁ + ... + Pₙ₋₁
     * @param n Number of qubits
     * @param pauli_char Pauli operator ('X', 'Y', 'Z')
     * @return n-qubit operator
     */
    static Matrix sum_pauli(int n, char pauli_char);

    /**
     * @brief Matrix exponential: exp(-i * t * H)
     * @param H Hermitian operator
     * @param t Evolution time
     * @return Unitary operator
     */
    static Matrix expm(const Matrix& H, double t = 1.0);

    /**
     * @brief Single-qubit rotation around X-axis
     * @param n Number of qubits
     * @param qubit_idx Which qubit
     * @param theta Rotation angle
     * @return n-qubit unitary
     */
    static Matrix Rx(int n, int qubit_idx, double theta);

    /**
     * @brief Single-qubit rotation around Y-axis
     */
    static Matrix Ry(int n, int qubit_idx, double theta);

    /**
     * @brief Single-qubit rotation around Z-axis
     */
    static Matrix Rz(int n, int qubit_idx, double theta);

    /**
     * @brief Universal Multi-Qubit gate: UMQ(χ) = exp(-i(χ/4)(∑Xᵢ)²)
     * @param n Number of qubits
     * @param chi Interaction strength
     * @return Global MS-like gate
     */
    static Matrix UMQ(int n, double chi);

    /**
     * @brief Synthesize arbitrary Pauli string unitary
     * @param axes Vector of Pauli operators ['X', 'Y', 'Z']
     * @param t Evolution time
     * @return Unitary operator exp(-i t P₁⊗P₂⊗...⊗Pₙ)
     */
    static Matrix n_body_string_arbitrary(const std::vector<char>& axes, double t);

    /**
     * @brief Original UMQ-Rz-UMQ pattern (hardware-native)
     * @param axes Pauli string (first variable, rest X)
     * @param t Evolution time
     * @param flip_sign Whether to flip evolution sign
     * @return Unitary operator
     */
    static Matrix n_body_string(const std::vector<char>& axes, double t,
                                 bool flip_sign = false);

    /**
     * @brief Generate target Pauli string unitary for comparison
     * @param pauli_string String like "ZXX" or "XYZ"
     * @param t Evolution time
     * @return Ideal unitary
     */
    static Matrix target_pauli_string_unitary(const std::string& pauli_string, double t);

    /**
     * @brief Phase-invariant unitary distance
     * @param U First unitary
     * @param V Second unitary
     * @return Distance metric d(U,V) = min_φ ||U - e^(iφ)V|| / ||U||
     */
    static double unitary_distance(const Matrix& U, const Matrix& V);
};

/**
 * @brief Interaction graph engineering
 */
class InteractionEngineering {
public:
    /**
     * @brief Check if interaction matrix is accessible
     * @param J_desired Desired coupling matrix
     * @param B Mode matrix
     * @param tol Tolerance for off-diagonal elements
     * @return True if exactly realizable with global beams
     */
    static bool is_accessible(const RealMatrix& J_desired,
                             const RealMatrix& B,
                             double tol = 1e-10);

    /**
     * @brief Get mode weights if accessible
     * @param J_desired Desired coupling matrix
     * @param B Mode matrix
     * @return Mode weights (empty if not accessible)
     */
    static RealVector get_mode_weights_if_accessible(const RealMatrix& J_desired,
                                                     const RealMatrix& B);

    /**
     * @brief Result structure for mode weight optimization
     */
    struct OptimizationResult {
        RealVector weights;
        RealMatrix J_achieved;
        double infidelity;
        bool accessible;
        std::string method;
    };

    /**
     * @brief Optimize mode weights to approximate inaccessible graphs
     * @param J_desired Desired coupling matrix
     * @param B Mode matrix
     * @param method "least_squares" or "linear_program"
     * @return Optimization result
     */
    static OptimizationResult optimize_mode_weights(const RealMatrix& J_desired,
                                                    const RealMatrix& B,
                                                    const std::string& method = "least_squares");

    /**
     * @brief Compute sinusoidal modes for equispaced ions
     * @param N Number of ions
     * @return N×N mode matrix (Kyprianidis Eq. 18)
     */
    static RealMatrix compute_sinusoidal_modes(int N);

    /**
     * @brief Compute Ising couplings from multi-tone driving
     * @param B Mode matrix
     * @param omega Mode frequencies
     * @param mus Drive detunings
     * @param Omegas Rabi frequencies
     * @param R_recoil Recoil frequency
     * @return Coupling matrix J_ij (Richerme Eq. 4)
     */
    static RealMatrix Jij_from_multimode(const RealMatrix& B,
                                        const RealVector& omega,
                                        const RealVector& mus,
                                        const RealVector& Omegas,
                                        double R_recoil);
};

/**
 * @brief Main interface class for gate synthesis
 */
class RichermeIonAnalog {
public:
    /**
     * @brief Constructor
     * @param use_hardware_gates Whether to use hardware-native synthesis
     */
    explicit RichermeIonAnalog(bool use_hardware_gates = true);

    /**
     * @brief Get hardware specifications
     */
    const IonTrapHardware& get_hardware() const { return hardware_; }

    /**
     * @brief Synthesize gate using configured method
     */
    Matrix synthesize_gate(const std::vector<char>& axes, double t) const;

    /**
     * @brief Check interaction accessibility
     */
    bool check_accessibility(const RealMatrix& J_desired, const RealMatrix& B) const;

private:
    bool use_hardware_gates_;
    IonTrapHardware hardware_;
};

} // namespace richerme

#endif // RICHERME_ION_ANALOG_H
