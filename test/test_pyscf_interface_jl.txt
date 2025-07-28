@testset "PySCF Interface Tests" begin
    config = NEOConfig()
    
    @testset "PySCF Setup" begin
        pyscf, has_neo = SparseQEEcNEO.PySCFInterface.setup_pyscf(config)
        
        if pyscf !== nothing
            @test pyscf !== nothing
            # Don't fail if NEO not available, just warn
            if !has_neo
                @test_skip "NEO module not available"
            else
                @test has_neo == true
            end
        else
            @test_skip "PySCF not available"
        end
    end
    
    # Only run remaining tests if PySCF is available
    pyscf, has_neo = SparseQEEcNEO.PySCFInterface.setup_pyscf(config)
    
    if has_neo
        @testset "NEO Molecule Building" begin
            mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
            mol_neo = SparseQEEcNEO.PySCFInterface.build_neo_molecule(mol, pyscf)
            
            @test pybuiltin("hasattr")(mol_neo, "nuc_num")
            @test mol_neo.nuc_num == 1
            @test pybuiltin("hasattr")(mol_neo, "quantum_nuc")
        end
        
        @testset "Mean-field Calculation" begin
            mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")
            mol_neo = SparseQEEcNEO.PySCFInterface.build_neo_molecule(mol, pyscf)
            
            calc = NEOCalculation(xc="HF")
            mf = SparseQEEcNEO.PySCFInterface.run_neo_meanfield(mol_neo, calc, pyscf)
            
            @test pybuiltin("hasattr")(mf, "e_tot")
            @test mf.e_tot < 0  # Energy should be negative
            @test pybuiltin("hasattr")(mf, "converged")
        end
        
        @testset "NEO Mean-field with Quantum Nuclei" begin
            mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
            mol_neo = SparseQEEcNEO.PySCFInterface.build_neo_molecule(mol, pyscf)
            
            calc = NEOCalculation(xc="B3LYP", epc="17-2")
            mf = SparseQEEcNEO.PySCFInterface.run_neo_meanfield(mol_neo, calc, pyscf)
            
            @test pybuiltin("hasattr")(mf, "components")
            @test haskey(mf.components, "e")
            @test haskey(mf.components, "n0")
        end
        
        @testset "MP2 Calculation" begin
            mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g")
            mol_neo = SparseQEEcNEO.PySCFInterface.build_neo_molecule(mol, pyscf)
            
            calc = NEOCalculation(xc="HF")
            mf = SparseQEEcNEO.PySCFInterface.run_neo_meanfield(mol_neo, calc, pyscf)
            
            mp, ecorr, t2 = SparseQEEcNEO.PySCFInterface.run_neo_mp2(mf, mol_neo)
            
            @test ecorr < 0  # Correlation energy should be negative
            @test t2 !== nothing
        end
    else
        @test_skip "All PySCF interface tests - NEO not available"
    end
end
