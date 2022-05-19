"""
Empirical Mode Decomposition Package.
"""
module EMD

using Statistics: mean
using Dierckx: Spline1D

## include environment setup.
include("util.jl")
include("config.jl")
include("emdbase.jl")

EMD

end # module
