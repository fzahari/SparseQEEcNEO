@testset "EPC Functional Tests" begin
    @testset "EPC Parameters" begin
        # Test all EPC functionals have correct parameters
        for (name, params) in EPC_PARAMS
            @test haskey(params, "a")
            @test haskey(params, "b")
            @test haskey(params, "c")
            @test params["a"] >= 0
            @test params["b"] >= 0
            @test params["c"] >= 0
        end
        
        # Test specific values
        @test EPC_PARAMS["17-1"]["a"] ≈ 0.735
        @test EPC_PARAMS["17-2"]["a"] ≈ 0.760
        @test EPC_PARAMS["18-1"]["a"] ≈ 1.008
        @test EPC_PARAMS["18-2"]["a"] ≈ 0.895
    end
    
    @testset "EPC17 Energy Calculation" begin
        # Create test densities
        n_points = 100
        rho_e = 0.1 * ones(n_points)  # Uniform electron density
        rho_p = 0.01 * exp.(-collect(0:n_points-1) / 10)  # Exponential proton density
        
        # Test epc17-1
        params = EPC_PARAMS["17-1"]
        energy = SparseQEEcNEO.EPCFunctionals.calculate_epc17_energy(
            rho_e, rho_p, params["a"], params["b"], params["c"]
        )
        
        @test energy < 0  # Should be negative (attractive)
        @test isfinite(energy)
        
        # Test with zero densities
        energy_zero = SparseQEEcNEO.EPCFunctionals.calculate_epc17_energy(
            zeros(10), zeros(10), params["a"], params["b"], params["c"]
        )
        @test energy_zero ≈ 0.0
    end
    
    @testset "EPC18 Energy Calculation" begin
        # Create test densities
        n_points = 100
        rho_e = 0.1 * ones(n_points)
        rho_p = 0.01 * ones(n_points)
        
        # Test epc18-2
        params = EPC_PARAMS["18-2"]
        energy = SparseQEEcNEO.EPCFunctionals.calculate_epc18_energy(
            rho_e, rho_p, params["a"], params["b"], params["c"]
        )
        
        @test energy < 0
        @test isfinite(energy)
        
        # Compare epc17 vs epc18 for same densities
        params17 = EPC_PARAMS["17-2"]
        energy17 = SparseQEEcNEO.EPCFunctionals.calculate_epc17_energy(
            rho_e, rho_p, params17["a"], params17["b"], params17["c"]
        )
        
        # They should be different but same order of magnitude
        @test abs(energy - energy17) > 1e-6
        @test abs(energy / energy17) < 10
    end
    
    @testset "Density Extraction" begin
        # Create mock mean-field object
        mf = Dict()
        mf_e = Dict()
        mf_e["make_rdm1"] = () -> diagm([2.0, 2.0, 0.0, 0.0])
        mf["components"] = Dict("e" => mf_e)
        
        # Test electron density extraction
        rho_e = SparseQEEcNEO.EPCFunctionals.get_electron_density(mf)
        
        @test length(rho_e) == 4
        @test rho_e[1] ≈ 2.0
        @test rho_e[2] ≈ 2.0
        @test rho_e[3] ≈ 0.0
    end
    
    @testset "EPC Functional Integration" begin
        # Test that calculate_epc_energy handles all functional types
        mf = Dict()
        mol_neo = Molecule("H2", "sto-3g")
        
        # Mock densities
        SparseQEEcNEO.EPCFunctionals.get_electron_density = mf -> ones(10)
        SparseQEEcNEO.EPCFunctionals.get_proton_density = (mf, mol) -> ones(10)
        
        for epc_type in ["none", "17-1", "17-2", "18-1", "18-2"]
            energy = SparseQEEcNEO.EPCFunctionals.calculate_epc_energy(mf, mol_neo, epc_type)
            
            if epc_type == "none"
                @test energy ≈ 0.0
            else
                @test energy != 0.0
                @test isfinite(energy)
            end
        end
    end
end
