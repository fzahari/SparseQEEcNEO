@testset "Nuclear Methods Tests" begin
    @testset "Nuclear Orbital Coupling" begin
        # Create mock nuclear and electronic mean-field
        mf_nuc = Dict()
        mf_nuc["mo_energy"] = collect(range(-0.1, 0.1, length=5))
        
        mf_elec = Dict()
        mf_elec["mo_energy"] = collect(range(-2.0, 2.0, length=10))
        
        # Test coupling calculation
        coupling = SparseQEEcNEO.NuclearMethods.calculate_nuclear_orbital_coupling(
            mf_nuc, mf_elec
        )
        
        @test length(coupling) == 5
        @test all(0 <= c <= 1 for c in coupling)
        @test issorted(coupling, rev=true)  # Should decay with energy
    end
    
    @testset "Orbital Selection Scores" begin
        # Mock data
        mf_nuc = Dict()
        mf_nuc["mo_energy"] = [-0.5, -0.1, 0.0, 0.1, 0.5]
        
        mf_elec = Dict()
        mf_elec["mo_energy"] = collect(range(-2.0, 2.0, length=10))
        
        mo_energy = mf_nuc["mo_energy"]
        occ_orbital = 2  # Second orbital is occupied
        
        scores = SparseQEEcNEO.NuclearMethods.calculate_orbital_selection_scores(
            mf_nuc, mf_elec, occ_orbital, mo_energy
        )
        
        @test length(scores) == 5
        @test scores[occ_orbital] ≈ 1.0  # Occupied orbital should have max score
        @test all(0 <= s <= 1 for s in scores)
    end
    
    @testset "Nuclear Component Validation" begin
        # Create mock mean-field with nuclear components
        mf = Dict()
        mf["components"] = Dict()
        
        # Add electronic component
        mf_e = Dict()
        mf_e["mo_occ"] = [1.0, 1.0, 0.0, 0.0]
        mf_e["mo_energy"] = [-1.0, -0.5, 0.5, 1.0]
        mf["components"]["e"] = mf_e
        
        # Add nuclear component
        mf_n = Dict()
        mf_n["mo_occ"] = [1.0, 0.0, 0.0]
        mf_n["mo_energy"] = [-0.1, 0.0, 0.1]
        mf["components"]["n0"] = mf_n
        
        # Mock molecule
        mol_neo = Dict()
        mol_neo["nuc_num"] = 1
        
        # Test validation
        valid = SparseQEEcNEO.NuclearMethods.validate_nuclear_components(mf, mol_neo)
        @test valid == true
        
        # Test with missing component
        delete!(mf["components"], "n0")
        valid = SparseQEEcNEO.NuclearMethods.validate_nuclear_components(mf, mol_neo)
        @test valid == false
    end
    
    @testset "Nuclear Kinetic Energy" begin
        # Mock nuclear mean-field
        mf_nuc = Dict()
        
        # Mock kinetic energy integrals
        T = diagm([0.5, 1.0, 1.5])
        mf_nuc["kinetic"] = () -> T
        
        # Mock density matrix
        dm = diagm([1.0, 0.0, 0.0])
        mf_nuc["make_rdm1"] = () -> dm
        
        # Calculate kinetic energy
        kinetic = SparseQEEcNEO.NuclearMethods.calculate_nuclear_kinetic_energy(mf_nuc)
        
        @test kinetic ≈ 0.5  # T[1,1] * dm[1,1]
        @test kinetic > 0
    end
    
    @testset "Orbital Truncation" begin
        # Create mock mean-field
        mf = Dict()
        mf["components"] = Dict()
        
        # Electronic component
        mf_e = Dict()
        mf_e["mo_occ"] = ones(10)
        mf["components"]["e"] = mf_e
        
        # Nuclear component with many orbitals
        mf_n = Dict()
        mf_n["mo_occ"] = [1.0, zeros(19)...]  # 20 orbitals, first occupied
        mf_n["mo_energy"] = collect(range(-0.5, 1.5, length=20))
        mf_n["mo_coeff"] = Matrix{Float64}(I, 20, 20)
        mf["components"]["n0"] = mf_n
        
        # Mock molecule
        mol_neo = Dict()
        mol_neo["nuc_num"] = 1
        mol_neo["quantum_nuc"] = [0]
        
        # Test truncation
        config_sel = ConfigSelection(max_nuc_orbs=5)
        
        # The actual truncation modifies mf in-place
        # We test that it would select the right number
        scores = SparseQEEcNEO.NuclearMethods.calculate_orbital_selection_scores(
            mf_n, mf_e, 1, mf_n["mo_energy"]
        )
        
        @test length(scores) == 20
        @test scores[1] ≈ 1.0  # Occupied orbital should have highest score
        
        # Top 5 should include the occupied orbital
        top_5 = partialsortperm(scores, 1:5, rev=true)
        @test 1 in top_5
    end
end
