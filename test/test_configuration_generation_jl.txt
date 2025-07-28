@testset "Configuration Generation Tests" begin
    # Create mock mean-field object
    function create_mock_mf()
        mf = Dict()
        mf_e = Dict()
        mf_e["mo_occ"] = [1.0, 1.0, 0.0, 0.0]
        mf_e["mo_energy"] = [-1.0, -0.5, 0.5, 1.0]
        mf["components"] = Dict("e" => mf_e)
        mf["e_tot"] = -1.0
        return mf
    end
    
    @testset "MP2 Configuration Generation" begin
        mf = create_mock_mf()
        mol = Molecule("H2", "sto-3g")
        config_sel = ConfigSelection(method="mp2", max_configs=10)
        
        configs = SparseQEEcNEO.ConfigurationGeneration.generate_configurations_mp2(mf, mol, config_sel)
        
        @test length(configs) > 0
        @test configs[1].name == "HF"  # Reference should be first
        @test all(c.weight > 0 for c in configs)
        @test issorted([c.weight for c in configs], rev=true)  # Sorted by weight
    end
    
    @testset "Configuration Compression" begin
        configs = [
            Configuration("S(1→3)", Dict("e" => Int16[0,1,1,0,0,0]), 0.5),
            Configuration("S(2→4)", Dict("e" => Int16[1,0,0,1,0,0]), 0.3)
        ]
        
        compressed = SparseQEEcNEO.ConfigurationGeneration.compress_configurations(configs)
        
        @test length(compressed) == length(configs)
        @test all(isa(c, CompressedConfig) for c in compressed)
        @test compressed[1].weight ≈ 0.5
    end
    
    @testset "Reference Configuration" begin
        mf = create_mock_mf()
        mol = Molecule("H2", "sto-3g")
        
        ref = SparseQEEcNEO.ConfigurationGeneration.create_reference_config(mf, mol)
        
        @test ref.name == "HF"
        @test ref.weight == 1.0
        @test haskey(ref.occupations, "e")
        @test ref.occupations["e"] == Int16[1, 1, 0, 0]
    end
    
    @testset "Configuration Selection" begin
        configs = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0),
            Configuration("S1", Dict("e" => Int16[0,1,1,0]), 0.5),
            Configuration("S2", Dict("e" => Int16[1,0,0,1]), 0.3),
            Configuration("S3", Dict("e" => Int16[0,1,0,1]), 0.1)
        ]
        
        config_sel = ConfigSelection(max_configs=2, target_importance=0.8)
        selected, n_qubits = SparseQEEcNEO.ConfigurationGeneration.select_important_configurations(configs, config_sel)
        
        @test length(selected) == 2
        @test selected[1].name == "HF"
        @test selected[2].name == "S1"
        @test n_qubits >= 2
    end
end
