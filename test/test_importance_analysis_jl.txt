@testset "Importance Analysis Tests" begin
    @testset "Configuration Type Analysis" begin
        configs = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0),
            Configuration("S(1→3)", Dict("e" => Int16[0,1,1,0]), 0.5),
            Configuration("D(12→34)", Dict("e" => Int16[0,0,1,1]), 0.2),
            Configuration("N1(1→2)", Dict("n0" => Int16[0,1,0]), 0.3),
            Configuration("C(E(1→3)+N1(1→2))", 
                         Dict("e" => Int16[0,1,1,0], "n0" => Int16[0,1,0]), 0.4)
        ]
        
        type_counts = SparseQEEcNEO.ImportanceAnalysis.analyze_configuration_types(configs)
        
        @test type_counts["reference"] == 1
        @test type_counts["single_elec"] == 1
        @test type_counts["double_elec"] == 1
        @test type_counts["single_nuc"] == 1
        @test type_counts["coupled"] == 1
    end
    
    @testset "Importance Metrics" begin
        configs = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0),
            Configuration("S1", Dict("e" => Int16[0,1,1,0]), 0.5),
            Configuration("S2", Dict("e" => Int16[1,0,0,1]), 0.3)
        ]
        
        mol = Molecule("H2", "sto-3g")
        config_sel = ConfigSelection(use_neo_importance=false)
        
        importance = SparseQEEcNEO.ImportanceAnalysis.calculate_importance_metrics(
            configs, mol, config_sel
        )
        
        @test importance.total_importance ≈ 1.0  # Normalized
        @test importance.neo_metrics === nothing  # No NEO for this case
        @test haskey(importance.config_types, "reference")
    end
    
    @testset "NEO Importance Metrics" begin
        configs = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0),
            Configuration("N1(1→2)", Dict("n0" => Int16[0,1,0]), 0.3),
            Configuration("C(E(1→3)+N1(1→2))", 
                         Dict("e" => Int16[0,1,1,0], "n0" => Int16[0,1,0]), 0.4)
        ]
        
        mol = Molecule("H2", "sto-3g", quantum_nuc=[0])
        config_sel = ConfigSelection(use_neo_importance=true)
        
        importance = SparseQEEcNEO.ImportanceAnalysis.calculate_importance_metrics(
            configs, mol, config_sel
        )
        
        @test importance.neo_metrics !== nothing
        @test importance.neo_metrics.nuclear_participation > 0
        @test importance.neo_metrics.coupling_contribution > 0
        @test importance.neo_metrics.enhancement_factor >= 1.0
    end
    
    @testset "Importance Distribution" begin
        # Create configs with exponential weight decay
        n_configs = 100
        configs = [
            Configuration("C$i", Dict("e" => Int16[1,1,0,0]), exp(-i/10))
            for i in 1:n_configs
        ]
        
        dist = SparseQEEcNEO.ImportanceAnalysis.analyze_importance_distribution(configs)
        
        @test dist.total_configs == n_configs
        @test dist.n_for_50_percent < dist.n_for_90_percent
        @test dist.n_for_90_percent < dist.n_for_99_percent
        @test dist.max_weight ≈ exp(-1/10)
        @test dist.weight_ratio > 1
    end
    
    @testset "Effective Rank" begin
        # Uniform distribution - high effective rank
        configs_uniform = [
            Configuration("C$i", Dict("e" => Int16[1,1,0,0]), 1.0)
            for i in 1:10
        ]
        rank_uniform = SparseQEEcNEO.ImportanceAnalysis.calculate_effective_rank(configs_uniform)
        @test rank_uniform ≈ 10 atol=0.1
        
        # Single dominant config - low effective rank
        configs_dominant = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0]), 100.0),
            Configuration("S1", Dict("e" => Int16[0,1,1,0]), 0.1),
            Configuration("S2", Dict("e" => Int16[1,0,0,1]), 0.1)
        ]
        rank_dominant = SparseQEEcNEO.ImportanceAnalysis.calculate_effective_rank(configs_dominant)
        @test rank_dominant < 2
    end
    
    @testset "Configuration Filtering" begin
        configs = [
            Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0),
            Configuration("S1", Dict("e" => Int16[0,1,1,0]), 0.5),
            Configuration("S2", Dict("e" => Int16[1,0,0,1]), 0.1),
            Configuration("S3", Dict("e" => Int16[0,1,0,1]), 0.01)
        ]
        
        # Filter by importance
        filtered = SparseQEEcNEO.ImportanceAnalysis.filter_by_importance(configs, 0.1)
        @test length(filtered) == 3
        @test all(c.weight >= 0.1 for c in filtered)
        
        # Adaptive threshold
        threshold = SparseQEEcNEO.ImportanceAnalysis.adaptive_importance_threshold(configs, 2)
        @test threshold ≈ 0.5
    end
end
