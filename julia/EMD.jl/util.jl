module Util

"""
    y = fliplr(x)

alias to flip the input vector `x`.
"""
fliplr(x::AbstractVector) = x[end:-1:1]


"""
    (tmin, tmax, zmin, zmax) = boundarycheck(indmin::Array, indmax::Array, t::Array, z::Array, nbsym::Int)

Boundary check to define extrema beyond the input signal limits to prevent boundary issues.
Without mirror symmetry, we get ramping at either ends of the IMF signal.
"""
function boundarycheck(
    indmin::AbstractVector{T}, 
    indmax::AbstractVector{T}, 
    t::AbstractVector, 
    x::AbstractVector, 
    z::AbstractVector, 
    nbsym::Int
    ) where {T<:Int}
    length(indmin) + length(indmax) ≥ 3 || error("EMD: boundscheck: not enough extrema.")
    lx = length(x)
    if (indmax[1] > indmin[1])
        if (x[1] > x[indmin[1]])
            lmax = fliplr(indmax[2:min(end, nbsym+1)])
            lmin = fliplr(indmin[1:min(end, nbsym)])
            lsym = copy(indmax[1])
        else
            lmax = fliplr(indmax[1:min(end, nbsym)])
            lmin = vcat(fliplr(indmin[1:min(end, nbsym-1)]), 1)
            lsym = 1.
        end
    else
        if (x[1] < x[indmax[1]])
            lmax = fliplr(indmax[1:min(end, nbsym)])
            lmin = fliplr(indmin[2:min(end, nbsym+1)])
            lsym = copy(indmin[1])
        else
            lmax = vcat(fliplr(indmax[1:min(end, nbsym-1)]), 1)
            lmin = fliplr(indmin[1:min(end, nbsym)])
            lsym = 1.
        end
    end
    if (indmax[end] < indmin[end])
        if (x[end] < x[indmax[end]])
            rmax = fliplr(indmax[max(end-nbsym+1, 1):end])
            rmin = fliplr(indmin[max(end-nbsym, 1):end-1])
            rsym = copy(indmin[end])
        else
            rmax = vcat(lx, fliplr(indmax[max(end-nbsym+2, 1):end]))
            rmin = fliplr(indmin[max(end-nbsym+1, 1):end])
            rsym = copy(lx)
        end
    else
        if (x[end] > x[indmin[end]])
            rmax = fliplr(indmax[max(end-nbsym, 1):end-1])
            rmin = fliplr(indmin[max(end-nbsym+1, 1):end])
            rsym = copy(indmax[end])
        else
            rmax = fliplr(indmax[max(end-nbsym+1, 1):end])
            rmin = vcat(lx, fliplr(indmin[max(end-nbsym+2, 1):end]))
            rsym = copy(lx)
        end
    end
    tlmin = 2 .* t[lsym] .- t[lmin]
    tlmax = 2 .* t[lsym] .- t[lmax]
    trmin = 2 .* t[rsym] .- t[rmin]
    trmax = 2 .* t[rsym] .- t[rmax]
    if (tlmin[1] > t[1]) || (tlmax[1] > t[1])
        if (lsym == indmax[1])
            lmax = fliplr(indmax[1:min(end,nbsym)])
        else
            lmin = fliplr(indmin[1:min(end,nbsym)])
        end
        lsym = 1.
        tlmin = 2 .* t[lsym] .- t[lmin]
        tlmax = 2 .* t[lsym] .- t[lmax]
    end

    if (trmin[end] < t[lx]) || (trmax[end] < t[lx])
        if (rsym == indmax[end])
            rmax = fliplr(indmax[max(end-nbsym+1, 1):end])
        else
            rmin = fliplr(indmin[max(end-nbsym+1, 1):end])
        end
        rsym = copy(lx)
        trmin = 2 .* t[rsym] .- t[rmin]
        trmax = 2 .* t[rsym] .- t[rmax]
    end

    zlmax = z[lmax]
    zlmin = z[lmin]
    zrmax = z[rmax]
    zrmin = z[rmin]

    tmin = vcat(tlmin, t[indmin], trmin)
    tmax = vcat(tlmax, t[indmax], trmax)
    zmin = vcat(zlmin, z[indmin], zrmin)
    zmax = vcat(zlmax, t[indmax], zrmax)

    tmin, tmax, zmin, zmax         
end

"""
   stop = stopemd(opts::EMDConfig, imf::AbstractVector)

Returns a flag indicating if at least 3 extrema are present to continue 
the decomposition.
"""
function stopemd(imf::AbstractVector)
    (indmin, indmax) = extr(imf)
    Bool(length(indmin) + length(indmax) < 3)
end

"""
   oind = orthoindex(imf::AbstractVector)

`orthoindex` computes the index of orthogonality based on the input
signal `x` and the prospective mode function `imf`.
"""
function orthoindex(x::AbstractVector)
    n = length(x)
    s = 0.0
    for i = 1:n
        for j = 1:n
            if i != j
                # we can improve this with `muladd`...
                s += abs(sum(x[i,:]) .* conj.(x[j,:]))/sum(x.^2)
            end
        end
    end
    0.5*s
end

"""
    (indmin, indmax) = extrminmax(x)
    
`extr` extracts the indices of extrema in value vector `x` over
domain `t`. min/max comparison attempts to match MATLAB behavior.
"""
function extrminmax(x::AbstractVector)
    dx = diff(x)
    m = length(dx)
    d1 = dx[1:m-1]
    d2 = dx[2:m]
    bad = (dx .== 0)
    indmin = findall((d1.*d2 .< 0) .& (d1 .< 0)) .+ 1
    indmax = findall((d1.*d2 .< 0) .& (d1 .> 0)) .+ 1
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

"""
    indzers = extrzeros(x)

Extract indices of zero-crossing extrema.
"""
function extrzeros(x::AbstractVector)
    x1 = x[1:end-1]
    x2 = x[2:end]
    indzer = findall(x1.*x2 .< 0)
    iz = (x .== 0)
    if any(iz)
        zer = findall(x == 0)
        if any(diff(iz) == 1)
            dz = diff(vcat(false, iz, false))
            headz = findall(dz .== 1)
            tailz = findall(dz .== -1) .- 1
            zhalf = zhalf == 0.5 ? Int(1) : round(Int, (headz + tailz) / 2)
        else
            indz = iz
        end
        sort!(append!(indzer, indz))
    end
    indzer
end

end # module