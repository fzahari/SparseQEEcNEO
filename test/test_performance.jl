@testset "Performance Tests" begin
    @testset "Configuration Generation Scaling" begin
        times = Float64[]
        sizes = [10, 20, 50, 100]
        
        for n in sizes
            # Create mock configurations
            configs = Configuration[]
            
            t0 = time()
            for i in 1:n
                occ = zeros(Int16, 20)
                occ[1:10] .= 1
                if i > 1
                    # Create excitation
                    occ[rand(1:10)] = 0
                    occ[rand(11:20)] = 1
                end
                config = Configuration("test_$i", Dict("e" => occ), rand())
                push!(configs, config)
            end
            t1 = time()
            
            push!(times, t1 - t0)
        end
        
        # Check reasonable scaling
        @test all(t < 1.0 for t in times)  # Should be fast
        @test issorted(times)  # Should scale monotonically
    end
    
    @testset "Hamiltonian Construction Scaling" begin
        for n_configs in [10, 50, 100]
            configs = [
                Configuration("C$i", Dict("e" => Int16[1,1,0,0]), rand()) 
                for i in 1:n_configs
            ]
            
            h1e = Dict("e" => randn(4,4))
            h1e["e"] = (h1e["e"] + h1e["e"]') / 2
            h2e = Dict("e" => Dict{NTuple{4,Int}, Float64}())
            coupling = Dict{Tuple{String,String}, Matrix{Float64}}()
            active = Dict("e" => [1,2,3,4])
            
            t0 = time()
            H = SparseQEEcNEO.HamiltonianConstruction.build_configuration_hamiltonian(
                configs, h1e, h2e, coupling, active
            )
            t1 = time()
            
            @test size(H) == (n_configs, n_configs)
            @test t1 - t0 < 0.1 * n_configs  # Should scale linearly
            @test ishermitian(H)
        end
    end
    
    @testset "Memory Usage - Compression" begin
        # Test configuration compression efficiency
        n_orbs = 100
        n_configs = 1000
        
        # Uncompressed
        configs_full = Configuration[]
        for i in 1:n_configs
            occ = zeros(Int16, n_orbs)
            # Sparse occupation (only 5 occupied)
            for j in 1:5
                occ[rand(1:n_orbs)] = 1
            end
            push!(configs_full, Configuration("C$i", Dict("e" => occ), rand()))
        end
        mem_full = Base.summarysize(configs_full)
        
        # Compressed
        configs_comp = SparseQEEcNEO.ConfigurationGeneration.compress_configurations(configs_full)
        mem_comp = Base.summarysize(configs_comp)
        
        compression_ratio = mem_full / mem_comp
        @test compression_ratio > 2.0  # Should achieve good compression
        @test length(configs_comp) == length(configs_full)
        
        @info "Compression ratio: $(round(compression_ratio, digits=2))x"
    end
    
    if NEO_AVAILABLE
        @testset "Real Calculation Performance" begin
            mol = Molecule("H 0 0 0; H 0 0 0.74", "sto-3g", quantum_nuc=[0])
            config_sel = ConfigSelection(method="neo_cneo", max_configs=100, max_nuc_orbs=0)
            
            # Time the calculation
            t0 = time()
            results = sparse_qee_cneo(mol, config_sel=config_sel, neo_config=TEST_CONFIG)
            t1 = time()
            
            @test t1 - t0 < 30.0  # Should complete in reasonable time
            @test results.computation_time > 0
            @test results.memory_used > 0
            
            @info "H2 NEO calculation: $(round(t1-t0, digits=2))s, $(round(results.memory_used, digits=1)) MB"
        end
    end
end
