@testset "Hamiltonian Construction Tests" begin
    @testset "Sparse ERI Creation" begin
        eri = SparseQEEcNEO.HamiltonianConstruction.create_sparse_eri(4)
        
        @test isa(eri, Dict{NTuple{4,Int}, Float64})
        @test length(eri) > 0
        @test all(v > 0 for v in values(eri))  # All positive
        @test haskey(eri, (1,1,1,1))  # Diagonal elements exist
    end
    
    @testset "Hamiltonian Data Structure" begin
        # Create minimal Hamiltonian data
        h1e = Dict("e" => randn(4,4))
        h1e["e"] = (h1e["e"] + h1e["e"]') / 2  # Symmetrize
        
        h2e = Dict("e" => SparseQEEcNEO.HamiltonianConstruction.create_sparse_eri(4))
        coupling = Dict(("e", "n0") => randn(4,2))
        
        configs = [Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0)]
        
        ham_data = HamiltonianData(
            h1e, h2e, coupling,
            Dict("e" => 4), Dict("e" => 2),
            [1.0], configs, [0.0],
            Dict("e" => [1,2]), Dict("e" => Int[])
        )
        
        @test haskey(ham_data.h1e, "e")
        @test haskey(ham_data.h2e, "e")
        @test length(ham_data.configs) == 1
    end
    
    @testset "Matrix Element Calculation" begin
        h1e = Dict("e" => diagm([1.0, 2.0, 3.0, 4.0]))
        h2e = Dict("e" => Dict{NTuple{4,Int}, Float64}())
        coupling = Dict{Tuple{String,String}, Matrix{Float64}}()
        
        config = Configuration("HF", Dict("e" => Int16[1,1,0,0]), 1.0)
        
        # Diagonal element
        E = SparseQEEcNEO.HamiltonianConstruction.calculate_diagonal_element(config, h1e, h2e, coupling)
        @test E ≈ 3.0  # 1.0 + 2.0 from occupied orbitals
    end
    
    @testset "Hamiltonian Validation" begin
        h1e = Dict("e" => randn(4,4))
        h1e["e"] = (h1e["e"] + h1e["e"]') / 2
        
        h2e = Dict("e" => Dict{NTuple{4,Int}, Float64}())
        coupling = Dict{Tuple{String,String}, Matrix{Float64}}()
        
        ham_data = HamiltonianData(
            h1e, h2e, coupling,
            Dict("e" => 4), Dict("e" => 2),
            Float64[], Configuration[], Float64[],
            Dict("e" => [1,2,3,4]), Dict("e" => Int[])
        )
        
        @test SparseQEEcNEO.HamiltonianConstruction.validate_hamiltonian(ham_data) == true
    end
end
