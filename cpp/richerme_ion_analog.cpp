/**
 * @file richerme_ion_analog.cpp
 * @brief C++ implementation of Richerme Ion Analog library
 */

#include "richerme_ion_analog.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace richerme {

// ===== IonTrapHardware Implementation =====

double IonTrapHardware::recoil_frequency() const {
    const double hbar = 1.054571817e-34;  // J·s
    const double delta_k = 2.0 * M_PI / wavelength;
    return hbar * delta_k * delta_k / (2.0 * mass) / (2.0 * M_PI);
}

// ===== PauliOps Implementation =====

Matrix PauliOps::I() {
    Matrix id(2, 2);
    id << 1.0, 0.0,
          0.0, 1.0;
    return id;
}

Matrix PauliOps::X() {
    Matrix x(2, 2);
    x << 0.0, 1.0,
         1.0, 0.0;
    return x;
}

Matrix PauliOps::Y() {
    Matrix y(2, 2);
    y << 0.0, Complex(0, -1),
         Complex(0, 1), 0.0;
    return y;
}

Matrix PauliOps::Z() {
    Matrix z(2, 2);
    z << 1.0, 0.0,
         0.0, -1.0;
    return z;
}

Matrix PauliOps::pauli(char pauli_char) {
    switch (pauli_char) {
        case 'I': return I();
        case 'X': return X();
        case 'Y': return Y();
        case 'Z': return Z();
        default:
            throw std::invalid_argument("Invalid Pauli character");
    }
}

// ===== GateSynthesis Implementation =====

Matrix GateSynthesis::kronN(const std::vector<Matrix>& ops) {
    if (ops.empty()) {
        throw std::invalid_argument("Empty operator list");
    }

    Matrix result = ops[0];
    for (size_t i = 1; i < ops.size(); ++i) {
        int rows = result.rows() * ops[i].rows();
        int cols = result.cols() * ops[i].cols();
        Matrix temp(rows, cols);

        for (int r1 = 0; r1 < result.rows(); ++r1) {
            for (int c1 = 0; c1 < result.cols(); ++c1) {
                for (int r2 = 0; r2 < ops[i].rows(); ++r2) {
                    for (int c2 = 0; c2 < ops[i].cols(); ++c2) {
                        temp(r1 * ops[i].rows() + r2, c1 * ops[i].cols() + c2) =
                            result(r1, c1) * ops[i](r2, c2);
                    }
                }
            }
        }
        result = temp;
    }
    return result;
}

Matrix GateSynthesis::pauli_operator(int n, char pauli_char, int qubit_idx) {
    std::vector<Matrix> ops(n);
    for (int i = 0; i < n; ++i) {
        ops[i] = (i == qubit_idx) ? PauliOps::pauli(pauli_char) : PauliOps::I();
    }
    return kronN(ops);
}

Matrix GateSynthesis::sum_pauli(int n, char pauli_char) {
    int dim = 1 << n;  // 2^n
    Matrix sum = Matrix::Zero(dim, dim);
    for (int i = 0; i < n; ++i) {
        sum += pauli_operator(n, pauli_char, i);
    }
    return sum;
}

Matrix GateSynthesis::expm(const Matrix& H, double t) {
    // Use eigenvalue decomposition: exp(-i*t*H) = V * exp(-i*t*D) * V†
    Eigen::SelfAdjointEigenSolver<Matrix> solver(H);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed");
    }

    const auto& eigvals = solver.eigenvalues();
    const auto& eigvecs = solver.eigenvectors();

    // exp(-i * t * eigenvalue)
    Vector exp_eigs(eigvals.size());
    for (int i = 0; i < eigvals.size(); ++i) {
        exp_eigs(i) = std::exp(Complex(0, -t * eigvals(i)));
    }

    // Reconstruct: V * diag(exp_eigs) * V†
    Matrix result = eigvecs * exp_eigs.asDiagonal() * eigvecs.adjoint();
    return result;
}

Matrix GateSynthesis::Rx(int n, int qubit_idx, double theta) {
    Matrix H = 0.5 * pauli_operator(n, 'X', qubit_idx);
    return expm(H, theta);
}

Matrix GateSynthesis::Ry(int n, int qubit_idx, double theta) {
    Matrix H = 0.5 * pauli_operator(n, 'Y', qubit_idx);
    return expm(H, theta);
}

Matrix GateSynthesis::Rz(int n, int qubit_idx, double theta) {
    Matrix H = 0.5 * pauli_operator(n, 'Z', qubit_idx);
    return expm(H, theta);
}

Matrix GateSynthesis::UMQ(int n, double chi) {
    Matrix Sx = sum_pauli(n, 'X');
    Matrix H = 0.25 * (Sx * Sx);
    return expm(H, chi);
}

Matrix GateSynthesis::n_body_string_arbitrary(const std::vector<char>& axes, double t) {
    int n = axes.size();
    int dim = 1 << n;

    // Build basis rotation R
    Matrix R = Matrix::Identity(dim, dim);
    for (int i = 0; i < n; ++i) {
        if (axes[i] == 'Y') {
            R = Rz(n, i, -M_PI/2) * R;
        } else if (axes[i] == 'Z') {
            R = Ry(n, i, M_PI/2) * R;
        }
    }

    // Build XXX...X operator
    std::vector<Matrix> xx_ops(n, PauliOps::X());
    Matrix H_XXX = kronN(xx_ops);
    Matrix U_XXX = expm(H_XXX, t);

    // Conjugate: U = R† · XXX · R
    return R.adjoint() * U_XXX * R;
}

Matrix GateSynthesis::n_body_string(const std::vector<char>& axes, double t,
                                     bool flip_sign) {
    int n = axes.size();

    // Check pattern: first qubit variable, rest X
    for (size_t i = 1; i < axes.size(); ++i) {
        if (axes[i] != 'X') {
            throw std::invalid_argument("Original pattern expects X on qubits 1..n-1");
        }
    }

    char target_first = axes[0];
    double angle = flip_sign ? -2.0*t : 2.0*t;

    if (target_first == 'X') {
        Matrix Rpre = Ry(n, 0, -M_PI/2);
        Matrix Rpost = Ry(n, 0, M_PI/2);

        Matrix U = Rpre;
        U = UMQ(n, M_PI/2) * U;
        U = Rz(n, 0, angle) * U;
        U = UMQ(n, -M_PI/2) * U;
        U = Rpost * U;
        return U;

    } else if (target_first == 'Y') {
        Matrix Rpre = Rx(n, 0, -M_PI/2);
        Matrix Rpost = Rx(n, 0, M_PI/2);

        Matrix U = Rpre;
        U = UMQ(n, M_PI/2) * U;
        U = Rz(n, 0, angle) * U;
        U = UMQ(n, -M_PI/2) * U;
        U = Rpost * U;
        return U;

    } else {  // 'Z'
        Matrix Ucore = Rx(n, 0, -M_PI/2);
        Ucore = UMQ(n, M_PI/2) * Ucore;
        Ucore = Rz(n, 0, angle) * Ucore;
        Ucore = UMQ(n, -M_PI/2) * Ucore;
        Ucore = Rx(n, 0, M_PI/2) * Ucore;
        return Rx(n, 0, M_PI/2) * Ucore * Rx(n, 0, -M_PI/2);
    }
}

Matrix GateSynthesis::target_pauli_string_unitary(const std::string& pauli_string,
                                                   double t) {
    std::vector<Matrix> ops;
    for (char c : pauli_string) {
        ops.push_back(PauliOps::pauli(c));
    }
    Matrix H = kronN(ops);
    return expm(H, t);
}

double GateSynthesis::unitary_distance(const Matrix& U, const Matrix& V) {
    // Phase-invariant distance: d(U,V) = min_φ ||U - e^(iφ)V|| / ||U||
    Complex trace = (U.adjoint() * V).trace();
    double phase = std::arg(trace);

    Matrix diff = U - std::exp(Complex(0, phase)) * V;
    double norm_diff = diff.norm();
    double norm_U = U.norm();

    return norm_diff / norm_U;
}

// ===== InteractionEngineering Implementation =====

bool InteractionEngineering::is_accessible(const RealMatrix& J_desired,
                                          const RealMatrix& B,
                                          double tol) {
    RealMatrix C = B.transpose() * J_desired * B;

    // Check if off-diagonal elements are small
    double off_diag_norm = 0.0;
    for (int i = 0; i < C.rows(); ++i) {
        for (int j = 0; j < C.cols(); ++j) {
            if (i != j) {
                off_diag_norm += C(i, j) * C(i, j);
            }
        }
    }
    off_diag_norm = std::sqrt(off_diag_norm);

    return off_diag_norm < tol;
}

RealVector InteractionEngineering::get_mode_weights_if_accessible(
    const RealMatrix& J_desired, const RealMatrix& B) {

    if (!is_accessible(J_desired, B)) {
        return RealVector();  // Empty vector
    }

    RealMatrix C = B.transpose() * J_desired * B;
    return C.diagonal();
}

InteractionEngineering::OptimizationResult
InteractionEngineering::optimize_mode_weights(const RealMatrix& J_desired,
                                             const RealMatrix& B,
                                             const std::string& method) {
    int N = B.rows();

    // Build mode matrices J^(k) = b_k ⊗ b_k
    std::vector<RealMatrix> J_modes(N);
    for (int k = 0; k < N; ++k) {
        RealVector b_k = B.col(k);
        J_modes[k] = b_k * b_k.transpose();
    }

    // Least squares: solve for weights
    // Extract upper triangle elements
    std::vector<double> b_vec;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            b_vec.push_back(J_desired(i, j));
        }
    }

    int n_constraints = b_vec.size();
    Eigen::MatrixXd A(n_constraints, N);

    int row = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                A(row, k) = J_modes[k](i, j);
            }
            row++;
        }
    }

    Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size());

    // Solve least squares
    RealVector weights = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

    // Compute achieved J
    RealMatrix J_achieved = RealMatrix::Zero(N, N);
    for (int k = 0; k < N; ++k) {
        J_achieved += weights(k) * J_modes[k];
    }

    // Compute infidelity (Kyprianidis Eq. 12)
    RealMatrix J_des_tilde = J_desired;
    J_des_tilde.diagonal().setZero();

    RealMatrix J_exp_tilde = J_achieved;
    J_exp_tilde.diagonal().setZero();

    double inner = (J_exp_tilde.array() * J_des_tilde.array()).sum();
    double norm_exp = J_exp_tilde.norm();
    double norm_des = J_des_tilde.norm();

    double infidelity = 1.0;
    if (norm_exp * norm_des > 0) {
        infidelity = 0.5 * (1.0 - inner / (norm_exp * norm_des));
    }

    OptimizationResult result;
    result.weights = weights;
    result.J_achieved = J_achieved;
    result.infidelity = infidelity;
    result.accessible = is_accessible(J_desired, B);
    result.method = method;

    return result;
}

RealMatrix InteractionEngineering::compute_sinusoidal_modes(int N) {
    RealMatrix B(N, N);

    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < N; ++k) {
            if (k == 0) {
                // COM mode
                B(i, k) = 1.0 / std::sqrt(N);
            } else {
                // Higher modes
                double prefactor = std::sqrt(2.0 / N);
                double argument = (2*i + 1) * k * M_PI / (2.0 * N);
                B(i, k) = prefactor * std::cos(argument);
            }
        }
    }

    return B;
}

RealMatrix InteractionEngineering::Jij_from_multimode(const RealMatrix& B,
                                                      const RealVector& omega,
                                                      const RealVector& mus,
                                                      const RealVector& Omegas,
                                                      double R_recoil) {
    int N = B.rows();
    int M = mus.size();
    RealMatrix J = RealMatrix::Zero(N, N);

    // Richerme Eq. 4: J_ij = Σ_k Σ_m (Ω_m² R)/(μ_m² - ω_k²) · B_ik · B_jk
    for (int k = 0; k < N; ++k) {
        double c_k = 0.0;
        for (int m = 0; m < M; ++m) {
            c_k += (Omegas(m) * Omegas(m) * R_recoil) /
                   (mus(m) * mus(m) - omega(k) * omega(k));
        }

        // Add contribution J^(k) = c_k * (b_k ⊗ b_k)
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                J(i, j) += c_k * B(i, k) * B(j, k);
            }
        }
    }

    // Symmetrize and zero diagonal
    J = 0.5 * (J + J.transpose()).eval();
    J.diagonal().setZero();

    return J;
}

// ===== RichermeIonAnalog Implementation =====

RichermeIonAnalog::RichermeIonAnalog(bool use_hardware_gates)
    : use_hardware_gates_(use_hardware_gates) {
    // Initialize hardware parameters
}

Matrix RichermeIonAnalog::synthesize_gate(const std::vector<char>& axes, double t) const {
    if (use_hardware_gates_) {
        // Check if pattern is compatible with original UMQ-Rz-UMQ
        bool is_original_pattern = true;
        for (size_t i = 1; i < axes.size(); ++i) {
            if (axes[i] != 'X') {
                is_original_pattern = false;
                break;
            }
        }

        if (is_original_pattern) {
            return GateSynthesis::n_body_string(axes, t);
        } else {
            return GateSynthesis::n_body_string_arbitrary(axes, t);
        }
    } else {
        // Use arbitrary pattern (ideal gates)
        return GateSynthesis::n_body_string_arbitrary(axes, t);
    }
}

bool RichermeIonAnalog::check_accessibility(const RealMatrix& J_desired,
                                            const RealMatrix& B) const {
    return InteractionEngineering::is_accessible(J_desired, B);
}

} // namespace richerme
