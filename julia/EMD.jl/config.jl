### possible interpolation modes for EMD.
@enum InterpMode Linear Pchip Cubic Spline

"""
    EMDConfig(Opts...)
EMDConfig sets up the configuration for a decomposition. This
does not include the input signal, but purely the options surrounding
sifting, maximum iterations, the type of interpolation between extrema, etc.
"""
struct EMDConfig{T<:AbstractFloat}
    maxiters::Int
    maxmodes::Int
    ndirs::Int
    sift::Int
    sift_h::Int
    taxis::StepRange
    stop::AbstractVector{T}
    x::AbstractVector{T}
    N::Int
    interp::InterpMode
    function EMDConfig(x::AbstractVector{T}) where T
        N = length(x)
        t = StepRange(1, 1, N)
        new{T}(Int(2_000), typemax(Int), Int(4), 0, 0, t, Array{T}([0.05, 0.5, 0.05]), x, N, Cubic)
    end
end

"""
    (tmin, tmax, zmin, zmax) = boundscheck(opts::EMDConfig, indmin::Array, indmax::Array, nbsym::Int)
Boundary check to define extrema beyond the input signal limits to prevent boundary issues.
Without mirror symmetry, we get ramping at either ends of the IMF signal.
"""
function boundscheck(indmin::AbstractVector, indmax::AbstractVector, opts::EMDConfig, z, nbsym::Int)
    length(indmin) + length(indmax) ≥ 3 || error("EMD: boundscheck: not enough extrema.")
    if (indmax[1] < indmin[1])
        if (opts.x[1] > opts.x[indmin[1]])
            # original EMD-m file used "end" keyword with minimum function??
        end
    end
end

"""
   stop = stopemd(opts::EMDConfig, imf::AbstractVector)
Returns a flag indicating if at least 3 extrema are present to continue 
the decomposition.
"""
function stopemd(otps::EMDConfig, imf::AbstractVector)
    (indmin, indmax) = extr(imf)
    Bool(length(indmin) + length(indmax) < 3)
end


function stopsift(opts::EMDConfig, imf::AbstractVector)
    (mean, npeaks, ndips, amp) = meanamplitude(opts, imf)
    Sx = abs.(mean) ./ amp
    S = mean(Sx)
end

"""
   oind = orthoindex(opts::EMDConfig, imf::AbstractVector)
`orthoindex` computes the index of orthogonality based on the input
signal `x` and the prospective mode function `imf`.
"""
function orthoindex(opts::EMDConfig, imf::AbstractVector)
    n = size(imf, 1)
    s = 0
    for i = 1:n
        for j = 1:n
            if i != j
                # we can improve this with `muladd`...
                s += abs(sum(imf[i,:]) .* conj.(imf[j,:]))/sum(x.^2)
            end
        end
    end
    0.5*s
end

"""
    (indmin, indmax) = extrminmax(x, t)
`extr` extracts the indices of extrema in value vector `x` over
domain `t`. min/max comparison attempts to match MATLAB behavior.
"""
function extrminmax(x::AbstractVector{T})
    dx = diff(x)
    m = length(x)
    d1 = dx[1:m-1]
    d2 = dx[2:m]
    bad = (dx .== 0)
    indmin = findall(d1.*d2<0 & d1<0) + 1
    indmax = findall(d1.*d2<0 & d1>0) + 1
    if any(bad)
        imax, imin = Int[], Int[]
        dd = diff(vcat(false, bad, false))
        head = findall(dd == 1)
        tail = findall(dd == -1)
        if (head[1] == 1)
            if length(head) > 1
                head = head[2:end]
                tail = tail[2:end]
            else
                head, tail = empty(head), empty(tail)
            end
        end
        if length(head) > 0
            if tail[end] == m
                if length(head) > 1
                    head = head[1:end-1]
                    tail = tail[1:end-1]
                else
                    head, tail = empty(head), empty(tail)
                end
            end
        end
        lh = length(head)
        if lh > 0
            for k = 1:lh
                # attempting to match rounding in MATLAB. 0.5 corner-case.
                half = (head[k] + tail[k]) / 2
                half = half == 0.5 ? Int(1) : round(Int, half)
                if d[head[k]-1] > 0
                    if d[tail[k]] < 0
                        append!(imax, half)
                    end
                else
                    if d[tail[k]] > 0
                        append!(imin, half)
                    end
                end
            end
        end
        if length(imax) > 0
            sort!(append!(indmax, imax))
        end
        if length(imin) > 0
            sort!(append!(indmin, imin))
        end
    end
    indmin, indmax
end
