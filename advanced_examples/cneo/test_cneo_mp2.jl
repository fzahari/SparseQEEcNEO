using Test
using SparseQEEcNEO

@testset "cNEO-MP2 Tests" begin
    
    @testset "H2 cNEO-MP2" begin
        # Create H2 molecule with first H as quantum
        mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
        
        # Setup configuration
        config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
        
        # Run cNEO-HF with MP2
        calc = NEOCalculation(xc="HF")
        config_sel = ConfigSelection(method="mp2", max_configs=50, max_nuc_orbs=0)
        
        results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel, neo_config=config)
        
        # Tests
        @test results.energy < 0
        @test results.mp2_correlation < 0  # Should be negative
        @test results.total_energy < results.energy  # MP2 should lower energy
        @test results.n_configs > 0
        
        # Check that MP2 correlation is reasonable
        @test -0.1 < results.mp2_correlation < 0  # Typical range for small molecules
    end
    
    @testset "HF cNEO-MP2" begin
        # Create HF molecule with H as quantum
        mol = Molecule("H 0 0 0; F 0 0 0.92", "sto-3g", quantum_nuc=[0])
        
        config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
        calc = NEOCalculation(xc="HF")
        config_sel = ConfigSelection(method="mp2", max_configs=100, max_nuc_orbs=0)
        
        results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel, neo_config=config)
        
        @test results.energy < -99  # HF ground state
        @test results.mp2_correlation < 0
        @test results.total_energy < results.energy
        
        # Check Hamiltonian
        @test results.hamiltonian_matrix !== nothing
        @test ishermitian(results.hamiltonian_matrix)
    end
    
    @testset "MP2 vs NEO-MP2 Comparison" begin
        # Test with and without quantum nuclei
        mol_classical = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")  # No quantum nuclei
        mol_quantum = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
        
        config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
        calc = NEOCalculation(xc="HF")
        config_sel = ConfigSelection(method="mp2", max_configs=50, max_nuc_orbs=0)
        
        # Classical calculation
        results_classical = sparse_qee_cneo(mol_classical, calc=calc, config_sel=config_sel, neo_config=config)
        
        # Quantum calculation
        results_quantum = sparse_qee_cneo(mol_quantum, calc=calc, config_sel=config_sel, neo_config=config)
        
        # The energies should be different
        @test abs(results_classical.total_energy - results_quantum.total_energy) > 1e-6
        
        # Both should have MP2 correlation
        @test results_classical.mp2_correlation < 0
        @test results_quantum.mp2_correlation < 0
    end
    
    @testset "MP2 with Different Basis Sets" begin
        # Test MP2 with different basis sets
        basis_sets = ["sto-3g", "6-31g"]
        energies = Dict{String, Float64}()
        
        config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
        calc = NEOCalculation(xc="HF")
        config_sel = ConfigSelection(method="mp2", max_configs=50, max_nuc_orbs=0)
        
        for basis in basis_sets
            mol = Molecule("H 0 0 0; H 0 0 0.74", basis, quantum_nuc=[0])
            
            try
                results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel, neo_config=config)
                energies[basis] = results.total_energy
                
                @test results.mp2_correlation < 0
                @test results.n_configs > 0
            catch e
                @warn "Basis set $basis failed: $e"
            end
        end
        
        # 6-31g should give lower energy than sto-3g
        if haskey(energies, "sto-3g") && haskey(energies, "6-31g")
            @test energies["6-31g"] < energies["sto-3g"]
        end
    end
    
    @testset "MP2 Configuration Selection" begin
        # Test that MP2 method generates appropriate configurations
        mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
        
        config = NEOConfig(pyscf_path=get(ENV, "PYSCF_PATH", "/path/to/pyscf-master"))
        calc = NEOCalculation(xc="HF")
        
        # Test with different max_configs
        for max_configs in [10, 50, 100]
            config_sel = ConfigSelection(method="mp2", max_configs=max_configs, max_nuc_orbs=0)
            results = sparse_qee_cneo(mol, calc=calc, config_sel=config_sel, neo_config=config)
            
            @test results.n_configs <= max_configs
            @test results.n_configs > 0
            @test results.captured_importance > 0
        end
    end
end
