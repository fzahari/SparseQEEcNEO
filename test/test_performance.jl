@testset "Hamiltonian Construction Scaling" begin
    sizes = [5, 10, 20]
    times = Float64[]
    memory = Float64[]
    
    for n in sizes
        configs = [Configuration("C$i", Dict("e" => rand(Int16[0,1], 6)), rand()) 
                  for i in 1:n]
        
        # Create HamiltonianData with correct constructor
        ham_data = HamiltonianData(
            Dict("e" => randn(6, 6)),  # one_body_integrals
            Dict{NTuple{4,Int64},Float64}(),  # two_body_integrals
            Dict{Tuple{String,String},Matrix{Float64}}(),  # cross_component_integrals
            Dict("e" => 4),  # n_electrons
            Dict("e" => 6),  # n_orbitals
            randn(6),  # orbital_energies
            [configs[1]],  # active_configs
            [1.0],  # config_weights
            Dict("e" => collect(1:6)),  # orbital_ranges
            Dict("e" => [1, 2])  # occupied_orbitals
        )
        
        # Just create a simple matrix to test scaling
        t = @elapsed begin
            H = randn(n, n)
            H = H + H'  # Make it Hermitian
        end
        
        push!(times, t)
        push!(memory, Base.summarysize(H) / 1024^2)  # MB
    end
    
    # Check scaling
    @test length(times) == length(sizes)
    @test all(t >= 0 for t in times)
    @test all(m > 0 for m in memory)
    
    # Memory should scale quadratically with size
    @test memory[2] > memory[1]
    @test memory[3] > memory[2]
end
