module Config

using ..Util

### possible interpolation modes for EMD.
@enum InterpMode Linear Pchip Cubic Spline

"""
    EMDConfig(Opts...)
EMDConfig sets up the configuration for a decomposition. This
does not include the input signal, but purely the options surrounding
sifting, maximum iterations, the type of interpolation between extrema, etc.
"""

Base.@kwdef struct EMDConfig
    maxiters::Int = 2_000
    maxmodes::Int = typemax(Int)
    ndirs::Int = 4
    stop::Vector = [0.05, 0.5, 0.05]
    interp::InterpMode = Cubic
end



    

    



end # module