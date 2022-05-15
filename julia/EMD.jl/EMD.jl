"""
Empirical Mode Decomposition Package.
"""
module EMD

import DSP, Distributions, Interpolations
using Statistics: mean
using Dierckx: Spline1D

## include environment setup.
include("util.jl")
include("config.jl")

EMD

end # module