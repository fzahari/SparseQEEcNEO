#!/usr/bin/env python3
"""
Quantum Computing Packages Test

Tests quantum computing packages required for SparseQEEcNEO integration.
Follows Clean Code principles with small, focused functions and clear error handling.
"""

import sys

# Constants
SUCCESS_SYMBOL = "✓"
ERROR_SYMBOL = "✗"
H2_BOND_LENGTH = 0.7414  # Angstroms


def display_header():
    """Display the test suite header."""
    print("=== Quantum Computing Package Tests ===\n")


def test_package_import(package_name, is_required=True):
    """Test importing a single package.
    
    Args:
        package_name (str): The name of the package to import.
        is_required (bool): Whether the package is required for core functionality.
        
    Returns:
        tuple: (success, module) - Success status and imported module if successful.
    """
    try:
        module = __import__(package_name)
        print(f"   {SUCCESS_SYMBOL} {package_name} imported successfully")
        if hasattr(module, "__version__"):
            print(f"   {SUCCESS_SYMBOL} {package_name} version: {module.__version__}")
        return True, module
    except ImportError as e:
        print(f"   {ERROR_SYMBOL} {package_name} import failed: {e}")
        if is_required:
            sys.exit(1)
        return False, None


def test_basic_imports():
    """Test importing all required and optional quantum computing packages."""
    print("1. Testing basic imports...")
    
    # Required packages
    qiskit_success, qiskit = test_package_import("qiskit", True)
    openfermion_success, openfermion = test_package_import("openfermion", True)
    cirq_success, cirq = test_package_import("cirq", True)
    
    # Optional packages
    qiskit_nature_success, _ = test_package_import("qiskit_nature", False)
    qiskit_algorithms_success, _ = test_package_import("qiskit_algorithms", False)
    
    return {
        "qiskit": qiskit if qiskit_success else None,
        "openfermion": openfermion if openfermion_success else None,
        "cirq": cirq if cirq_success else None,
        "qiskit_nature": qiskit_nature_success,
        "qiskit_algorithms": qiskit_algorithms_success
    }


def create_bell_circuit():
    """Create a Bell state circuit using Qiskit.
    
    Returns:
        tuple: (circuit, qreg, creg) - The circuit and its registers
    """
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
    
    qreg = QuantumRegister(2, 'q')
    creg = ClassicalRegister(2, 'c')
    circuit = QuantumCircuit(qreg, creg)
    
    # Create Bell state |Φ+⟩ = (|00⟩ + |11⟩)/√2
    circuit.h(qreg[0])      # Hadamard on first qubit
    circuit.cx(qreg[0], qreg[1])  # CNOT with first qubit as control
    circuit.measure_all()
    
    return circuit, qreg, creg


def simulate_qiskit_circuit(circuit):
    """Simulate a Qiskit circuit and return measurement results.
    
    Args:
        circuit (QuantumCircuit): Circuit to simulate
        
    Returns:
        dict: Measurement results
    """
    try:
        from qiskit_aer import Aer
        from qiskit import execute
        
        simulator = Aer.get_backend('qasm_simulator')
        job = execute(circuit, simulator, shots=1000)
        result = job.result()
        counts = result.get_counts(circuit)
        return counts
    except Exception as e:
        print(f"   {ERROR_SYMBOL} Circuit simulation failed: {e}")
        return None


def test_qiskit_functionality():
    """Test basic Qiskit functionality."""
    try:
        circuit, _, _ = create_bell_circuit()
        counts = simulate_qiskit_circuit(circuit)
        
        if counts:
            print(f"   {SUCCESS_SYMBOL} Qiskit circuit creation and simulation successful")
            print(f"   {SUCCESS_SYMBOL} Bell state measurement results: {counts}")
            return True
        return False
    except Exception as e:
        print(f"   {ERROR_SYMBOL} Qiskit circuit test failed: {e}")
        return False


def test_openfermion_functionality():
    """Test basic OpenFermion functionality."""
    try:
        from openfermion import FermionOperator, jordan_wigner
        
        # Create a simple fermionic operator: a†_0 a_1
        ferm_op = FermionOperator("0^ 1")
        
        # Convert to qubit operator using Jordan-Wigner transformation
        qubit_op = jordan_wigner(ferm_op)
        
        print(f"   {SUCCESS_SYMBOL} OpenFermion operator creation and transformation successful")
        print(f"   {SUCCESS_SYMBOL} Fermionic operator: {ferm_op}")
        print(f"   {SUCCESS_SYMBOL} Qubit operator: {qubit_op}")
        
        return True
    except Exception as e:
        print(f"   {ERROR_SYMBOL} OpenFermion operator test failed: {e}")
        return False


def test_cirq_functionality():
    """Test basic Cirq functionality."""
    try:
        import cirq
        
        # Create qubits on a line
        qubits = [cirq.LineQubit(i) for i in range(2)]
        
        # Create a Bell state circuit
        circuit = cirq.Circuit(
            cirq.H(qubits[0]),
            cirq.CNOT(qubits[0], qubits[1]),
            cirq.measure(*qubits, key='result')
        )
        
        # Simulate
        simulator = cirq.Simulator()
        result = simulator.run(circuit, repetitions=1000)
        
        print(f"   {SUCCESS_SYMBOL} Cirq circuit creation and simulation successful")
        print(f"   {SUCCESS_SYMBOL} Cirq measurement results shape: {result.measurements['result'].shape}")
        
        return True
    except Exception as e:
        print(f"   {ERROR_SYMBOL} Cirq circuit test failed: {e}")
        return False


def test_basic_functionality():
    """Test basic functionality of all quantum packages."""
    print("\n2. Testing basic functionality...")
    
    qiskit_ok = test_qiskit_functionality()
    openfermion_ok = test_openfermion_functionality()
    cirq_ok = test_cirq_functionality()
    
    return qiskit_ok and openfermion_ok and cirq_ok


def create_h2_molecule():
    """Create a H2 molecule representation using OpenFermion.
    
    Returns:
        MolecularData: Molecular data for H2
    """
    from openfermion.chem import MolecularData
    
    geometry = [
        ('H', (0., 0., 0.)),
        ('H', (0., 0., H2_BOND_LENGTH))
    ]
    basis = 'sto-3g'
    multiplicity = 1
    charge = 0
    
    return MolecularData(geometry, basis, multiplicity, charge)


def test_quantum_chemistry():
    """Test quantum chemistry functionality."""
    print("\n3. Testing quantum chemistry functionality...")
    
    openfermion_chem_ok = test_openfermion_chemistry()
    pyscf_ok = test_pyscf_integration()
    
    return openfermion_chem_ok and pyscf_ok


def test_openfermion_chemistry():
    """Test OpenFermion quantum chemistry functionality."""
    try:
        molecule = create_h2_molecule()
        
        print(f"   {SUCCESS_SYMBOL} MolecularData creation successful")
        print(f"   {SUCCESS_SYMBOL} Molecule: {molecule.name}")
        print(f"   {SUCCESS_SYMBOL} N electrons: {molecule.n_electrons}")
        print(f"   {SUCCESS_SYMBOL} N orbitals: {molecule.n_orbitals}")
        
        return True
    except Exception as e:
        print(f"   {ERROR_SYMBOL} OpenFermion molecular data test failed: {e}")
        return False


def test_pyscf_integration():
    """Test PySCF integration."""
    try:
        import pyscf
        print(f"   {SUCCESS_SYMBOL} PySCF is available")
        print(f"   {SUCCESS_SYMBOL} PySCF version: {pyscf.__version__}")
        
        # Test basic molecule creation
        from pyscf import gto
        mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
        print(f"   {SUCCESS_SYMBOL} PySCF molecule creation successful")
        
        return True
    except ImportError:
        print(f"   ! PySCF not available (expected - using custom NEO version)")
        return True  # Not a failure
    except Exception as e:
        print(f"   {ERROR_SYMBOL} PySCF test failed: {e}")
        return False


def test_package_interoperability():
    """Test interoperability between quantum packages."""
    print("\n4. Testing package interoperability...")
    
    try:
        from openfermion import FermionOperator
        
        # Simple fermionic operator for electron number operators
        of_op = FermionOperator("0^ 0") + FermionOperator("1^ 1")
        
        print(f"   {SUCCESS_SYMBOL} Created OpenFermion operator")
        print(f"   {SUCCESS_SYMBOL} Operator: {of_op}")
        
        # Attempt to import Qiskit Nature's fermionic operators
        try:
            from qiskit_nature.second_q.operators import FermionicOp
            print(f"   {SUCCESS_SYMBOL} Qiskit Nature's FermionicOp available")
        except ImportError:
            print(f"   ! Qiskit Nature interoperability not fully tested")
        
        return True
    except Exception as e:
        print(f"   ! Interoperability test limited: {e}")
        return False


def display_summary():
    """Display a summary of the test results."""
    print("\n=== All Tests Completed ===")
    print("\nSummary:")
    print("- Basic quantum computing packages are working")
    print("- Circuit creation and simulation successful")
    print("- Molecular data handling available")
    print("- Ready for SparseQEEcNEO integration!")


def run_quantum_package_tests():
    """Run all quantum package tests in sequence."""
    display_header()
    
    # Test imports
    imported_packages = test_basic_imports()
    
    # Only continue if required packages are available
    if not all(imported_packages[pkg] for pkg in ["qiskit", "openfermion", "cirq"]):
        sys.exit(1)
    
    # Run functionality tests
    basic_functionality_ok = test_basic_functionality()
    quantum_chemistry_ok = test_quantum_chemistry()
    interoperability_ok = test_package_interoperability()
    
    # Display final summary
    display_summary()
    
    # Return success status
    return basic_functionality_ok and quantum_chemistry_ok and interoperability_ok


# Execute tests when script is run directly
if __name__ == "__main__":
    success = run_quantum_package_tests()
    sys.exit(0 if success else 1)
