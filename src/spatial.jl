#########################################################################
# Spatial indexing
#########################################################################
export spatialidxs
spatialidxs(A::AbDimArray) = spatialidxs(spacestructure(A), A)
function spatialidxs(::Grid1Deg, A)
    lats = (Lat(i) for i in 1:size(A, Lat))
    lons = (Lon(i) for i in 1:size(A, Lon))
    return Iterators.product(lons, lats)
end

function spatialidxs(::EqArea, A)
    return ((Coord(i),) for i in 1:size(A, Coord))
end

#########################################################################
# averaging functions over space or time
#########################################################################
export latmean, spacemean, zonalmean, spaceagg

"""
    zonalmean(a::DimensionalArray [, r])
Return the zonally-meand `a`.
Optionally do the mean for the data in range `r` of that dimension.
(`r` is fed into the dimension, `Dim(r)`, so it can be a range or an
arbitrary selector).
"""
zonalmean(a::AbDimArray) = dropagg(mean, a, Lon)
zonalmean(a::AbDimArray, r) = dropagg(mean, a[Lon(r)], Lon)

function zonalmean(a::AbDimArray{T, 2, <:Tuple{<:Coord, <: Time}}) where {T}
    idxs, lats = uniquelats(a)
    res = zeros(T, length(lats), size(a, 2))
    for (i, r) in enumerate(idxs)
        for j in 1:size(a, 2)
            res[i, j] = mean(view(a.data, r, j))
        end
    end
    return DimensionalArray(res, (Lat(lats), dims(a, 2)), label(a))
end
function zonalmean(a::AbDimArray{T, 1, <:Tuple{<:Coord}}) where {T}
    idxs, lats = uniquelats(a)
    res = zeros(T, length(lats))
    for (i, r) in enumerate(idxs)
        res[i] = mean(view(a.data, r))
    end
    return DimensionalArray(res, (Lat(lats),))
end

uniquelats(a::AbDimArray) = uniquelats(dims(a, Coord))
function uniquelats(c)
    @assert issorted(c, by = x -> x[2])
    idxs = Vector{UnitRange{Int}}()
    lats = eltype(eltype(c))[]
    iprev = 1
    for i in 2:length(c)
        if c[i][2] != c[i-1][2]
            push!(lats, c[i-1][2])
            push!(idxs, iprev:(i-1))
            iprev = i
        end
    end
    push!(lats, c[end][2])
    push!(idxs, iprev:length(c))
    return idxs, lats
end


"""
    latmean(a::DimensionalArray [, r])
Return the latitude-mean `a` (mean across dimension `Lat`).
Optionally do the mean for the data in range `r` of that dimension.

This function properly weights the mean by the cosine of the latitude.
"""
function latmean(a::AbDimArray, r = 1:size(a, Lat))
    a = a[Lat(r)]
    lw = _latweights(a)
    dropagg(sum, a ⊗ lw, Lat)
end
# Warning!!! `_latweights` divides by the weight sum, because it is intended to be
# used only with the `sum` function (for a)
_latweights(a::AbDimArray) = _latweights(dims(a, Lat))
function _latweights(a::Lat)
    we = cosd.(Array(a))
    we ./= sum(we)
    return DimensionalArray(we, (a,), "weights")
end

spacemean(a::AbDimArray) = spacemean(spacestructure(a), a)
spacemean(::Grid1Deg, a) = latmean(zonalmean(a))
spacemean(::EqArea, a) = dropagg(mean, a, Coord)

using StatsBase

"`spacemean(a, exw=nothing) = spaceagg(mean, a, exw)`"
spacemean(a, exw=nothing) = spaceagg(mean, a, exw)

"""
    spaceagg(f, a::DimensionalArray, w = nothing)
Aggregate array `a` using function `f` (e.g. `mean`) over all available space (i.e.
longitude and latitude) of `a`, weighting every part of `a` by its spatial area.
The function works for grid space layout as well as equal area.

`w` can be extra weights, to weight each spatial point with. They can either be
just an `AbDimArray` with same spatial dimensions as `a`, or the can be of exactly
same shape as `a`.
"""
spaceagg(f, a::AbDimArray, exw=nothing) = spaceagg(spacestructure(a), f, a, exw)
function spaceagg(::Grid1Deg, f, a, exw=nothing)
    # This assumes that lon is first dim and lat is second dim.
    w = repeat(cosd.(Array(dims(a, Lat)))', length(dims(a, Lon)))
    if hasdim(a, Time)
        r = zeros(length(dims(a, Time)))
        for i in 1:length(r)
            # The following codeblocks assume that `exw` can either have LonLat or LonLatTime
            if !isnothing(exw) && hasdim(exw, Time)
                w = w .* exw[Time(i)]
            elseif !isnothing(exw)
                w = w .* Array(exw)
            end
            r[i] = f(Array(view(a, Time(i))), pweights(w))
        end
        return DimensionalArray(r, (dims(a, Time),))
    else
        # TODO: We need code that also works if both `w` and `a` contain Time.
        if !isnothing(exw)
            w = w .* Array(exw)
        end
        return f(Array(a), pweights(w))
    end
end

spaceagg(::EqArea, f, a) = dropagg(f, a, Coord)
# I think the best scenario is to modify `dropagg` to take in weights.
function spaceagg(::EqArea, f, a, exw)
    error("TODO")
    w = pweights(Array(exw))
    dropagg(f, a, Coord)
end
