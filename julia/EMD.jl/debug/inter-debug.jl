using Pkg
Pkg.activate(".")
using EMD
using DelimitedFiles: readdlm
using Dierckx: Spline1D
using Test

fpath = joinpath(@__DIR__, "signals")
t = collect(1:4096)

# manually confirmed the indices for minima/maxima matched MATLAB integers.
X = vec(readdlm(joinpath(fpath, "doppler.txt")))
(indmin, indmax) = EMD.extrminmax(X)
(jtmin, jtmax, jmmin, jmmax) = EMD.boundarycheck(indmin, indmax, t, X, X, 2)


tmin = vec(readdlm(joinpath(fpath, "tmin_iter1.txt")))
tmax = vec(readdlm(joinpath(fpath, "tmax_iter1.txt")))
mmin = vec(readdlm(joinpath(fpath, "mmin_iter1.txt")))
mmax = vec(readdlm(joinpath(fpath, "mmax_iter1.txt")))

all(jtmin .== tmin)
all(jtmax .== tmax)
all(jmmin .== mmin) # ERROR!!!
all(jmmax .== mmax)

envmin = vec(readdlm(joinpath(fpath, "envmin_iter1.txt")))
envmax = vec(readdlm(joinpath(fpath, "envmax_iter1.txt")))


## check if the interpolation method matches MATLAB.
jmin = Spline1D(tmin, mmin, k=3)(t)
jmax = Spline1D(tmax, mmax, k=3)(t)

mindiff = abs.(envmin - jmin)
maxdiff = abs.(envmax - jmax)

