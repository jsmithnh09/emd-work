"""
Empirical Mode Decomposition Package.
"""
module EMD

import DSP, Distributions, Interpolations
using Statistics: mean
using Dierckx: Spline1D

## include environment setup.
include("config.jl")
include("util.jl")

EMD

end # module