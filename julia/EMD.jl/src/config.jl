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

struct EMDConfig
    maxiters::Int
    maxmodes::Int
    ndirs::Int
    stop::Vector{AbstractFloat}
    interp::InterpMode
    function EMDConfig(
        maxiters=2_000, 
        maxmodes=typemax(Int), 
        ndirs=4, 
        stop=[0.05, 0.5, 0.05], 
        interp=Cubic
        )
        if length(stop) != 3
            throw(DimensionMismatch("stop vector must be 3-elements, instead got $(size(stop))"))
        end
        return new(maxiters, maxmodes, ndirs, stop, interp)
    end
end

### Show details about the config struct.
function Base.show(io::IO, cfg::EMDConfig)
    mmstr = (cfg.maxmodes == typemax(Int)) ? "∞" : "$(cfg.maxmodes)"
    print(io, "*** EMD Config ***\n") 
    print(io, "------------------\n")
    print(io, "maxmodes: ", mmstr, "\n")
    print(io, "stop: ", cfg.stop, "\n") 
    print(io, "interp: ", string(cfg.interp), "\n")
end

end # module