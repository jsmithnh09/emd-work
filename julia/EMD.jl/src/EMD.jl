"""
Empirical Mode Decomposition Package.
"""
module EMD

## include environment setup.
include("config.jl")
include("util.jl")

using ..Util, ..Config

export emd

"""
    imf = emd(x; kwargs...)

Given an input vector, Empirical Mode Decomposition is performed. Currently,
only cubic spline interpolation with default stoppage criterion is in place. The return
type `imf` will be a `Vector{Vector{Float64}}`.
"""
function emd(x::AbstractVector{T}; kwargs...) where {T <: AbstractFloat}
    
    # create the configuration object.
    cfg = EMDConfig(x, kwargs...)
    k, curiter = 1, 0
    r = copy(x)
    imf = Vector{Vector{T}}()

    while (!(stopemd(r)) && (k < cfg.maxmodes+1 || cfg.maxmodes == 0))
        m = copy(r)
        (stopsift, μenv, _) = stopsifting(m, cfg.stop[1], cfg.stop[2], cfg.stop[3])
        
        # sift loop
        while ((!stopsift) && (curiter < cfg.maxiters))
            m -= μenv
            (stopsift, μenv, _) = stopsifting(m, cfg.stop[1], cfg.stop[2], cfg.stop[3])
            curiter += 1
            
            # force sift stoppage in case there were too many iterations.
            if ((curiter == cfg.maxiters-1) && (curiter > 100))
                break
            end
        end
        push!(imf, m) # append IMF
        k += 1        # increment the mode iteration.
        r -= m        # extract mode from input.
        curiter = 0   # reset the sifting iteration.
    end
    if (any(r .> 0))
        push!(imf, r) # append residual if any fluctuations are still present.
    end
    imf
end # emd




EMD

end # module
