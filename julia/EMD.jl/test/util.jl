using Test, EMD.Util

in = [-4., 3., -2., 0., -4., 1., -3., 2., 2., 3.]
t = collect(1:length(in))

# test that the extrema extraction. using randi([-5, 5], 1, 10)
@testset "extrema" begin
    (minima, maxima) = extrminmax(in)
    @test minima == [3, 5, 7]
    @test maxima == [2, 4, 6]
    indzer = extrzeros(in)
    @test indzer == [1, 2, 4, 5, 6, 7]
end

# now check the boundary conditions with interpolation, (hardcoding nsym)
@testset "boundary conditions" begin
    (tmin, tmax, zmin, zmax) = boundarycheck(minima, maxima, t, in, in, 2)
    @test tmin == [-1, 1, 3, 5, 7, 13, 15]
    @test tmax == [-2, 0, 2, 4, 6, 10, 14]
    @test zmin == [-2, -4, -2, -4, -3, -3, -4]
    @test zmax == [0, 3, 3, 0, 1, 3, 1]
end

# testing mean/amp, which first requires boundary conditions, followed by
# extrema extraction.
@testset "mean/amplitude" begin
    (μ, nextr, nzer, amp) = meanamplitude(in)
    @test isapprox(μ, [-0.2293, 0.14191, -0.3121, -1.4473, -1.9597, -1.4225, -0.5689, 0.0649, 0.4105, 0.4945], rtol=1e-3)
    @test nextr == 6
    @test nzer == 6
    @test isapprox(amp, [3.7707, 2.8509, 1.6879, 1.4473, 2.0403, 2.4225, 2.4311, 2.4206, 2.4555, 2.5055], rtol=1e-3)
end



