using EMD, Test
using Wavelets

@testset "EMD utilities" begin
    include("util.jl")
end

# deconstruct the signal into IMFs, then sum to prove reconstruction.
@testset "Partial Reconstruction" begin
    for sig in (:Doppler, :HeaviSine, :Blocks, :Bumps)
        x = testfunction(2^14, "$(sig)")
        imf = emd(x)
        @test isapprox(sum(imfs), x, rtol=1e-8)
    end
end