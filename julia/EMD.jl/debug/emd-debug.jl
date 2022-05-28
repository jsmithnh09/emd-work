using Pkg
using Plots
Pkg.activate(".")
using EMD, DelimitedFiles

X = vec(readdlm(joinpath(@__DIR__, "signals", "doppler.txt")))
ref = reshape(readdlm(joinpath(@__DIR__, "signals", "imf4x4096.txt")), (4096, 4))


maxmodes = 0
maxiters = Int(2_000)

### Main EMD loop
k, curiter = 1, 0
r = copy(X)
maxmodes = 0
maxiters = Int(2_000)
imf = Vector{Vector{eltype(X)}}()
Xtol = √(eps(eltype(X))) * maximum(abs, X)


### first iteration.
m = copy(r)
(stopsift, μenv, _) = EMD.stopsifting(m, 0.05, 0.5, 0.05)

# check the noise floor flag.
noisefloor = Bool(maximum(abs, m) < Xtol)

## recurse until satisfied.
m -= μenv
(stopsift, μenv, _) = EMD.stopsifting(m, 0.05, 0.5, 0.05)

### if not stop sift, add.
push!(imf, m)


# for checking, specify index to compare.
ind = 1
dimf = abs.(imf[ind] - vec(ref[:,ind]))