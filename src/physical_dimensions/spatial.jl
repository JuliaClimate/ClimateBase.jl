#=
Functionality regarding spatial coordinates of an array, i.e. Longitude and Lattitude
as well as equal area grid handling and anything else related to space.
=#
using Statistics, StatsBase
#########################################################################
# Spatial indexing
#########################################################################
export spatialidxs

"""
    spatialidxs(A::ClimArray) → idxs
Return an iterable that can be used to access all spatial points of `A` with the syntax
```julia
idxs = spatialidxs(A)
for i in idxs
    slice_at_give_space_point = A[i]
end
```
Works for standard grid as well as equal area.
"""
spatialidxs(A::AbDimArray) = spatialidxs(spacestructure(A), A)
function spatialidxs(::Grid, A)
    lons = (Lon(i) for i in 1:size(A, Lon))
    lats = (Lat(i) for i in 1:size(A, Lat))
    return Iterators.product(lons, lats)
end

function spatialidxs(::EqArea, A)
    return (Coord(i) for i in 1:size(A, Coord))
end


"""
    wrap_lon(x)
Wrap given longitude to -180 to 180 degrees.
"""
wrap_lon(x) = @. -180 + (360 + ((x+180) % 360)) % 360

#########################################################################
# averaging functions over space or time
#########################################################################
using StatsBase

export latmean, spacemean, zonalmean, spaceagg

# TODO: Document what it means to be Coord space.
# We expect that the coordinates are sorted by latitude

"""
    zonalmean(A::ClimArray [, r])
Return the zonal mean of `A`.
Optionally do the mean for the data in range `r` of the longitude
(`r` is fed into the dimension so it can be A range or an arbitrary selector).

Also works for equal area space.
"""
zonalmean(A::AbDimArray, args...) = zonalmean(spacestructure(A), A, args...)
zonalmean(::Grid, A::AbDimArray) = dropagg(mean, A, Lon)
zonalmean(::Grid, A::AbDimArray, r) = dropagg(mean, A[Lon(r)], Lon)

# Version of arbitrary dims (now assumes only 2, must fix)
# TODO: Extend this for `A` with more than 2 dimensions
function zonalmean(::EqArea, A::AbDimArray) where {T}
    idxs, lats = uniquelats(A)
    res = zeros(T, length(lats), size(A, 2))
    for (i, r) in enumerate(idxs)
        for j in 1:size(A, 2)
            res[i, j] = mean(view(A.data, r, j))
        end
    end
    return ClimArray(res, (Lat(lats), dims(A, 2)), label(A))
end
function zonalmean(::EqArea, A::AbDimArray{T, 1, <:Tuple{<:Coord}}) where {T}
    idxs, lats = uniquelats(A)
    res = zeros(T, length(lats))
    for (i, r) in enumerate(idxs)
        res[i] = mean(view(A.data, r))
    end
    return ClimArray(res, (Lat(lats),))
end

"""
    uniquelats(A::AbDimArray) → idxs, lats
Find the unique latitudes of `A`. Return the indices (vector of ranges) that each latitude
in `lats` covers, as well as the latitudes themselves.
"""
uniquelats(A::AbDimArray) = uniquelats(dims(A, Coord))
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
    latmean(A::ClimArray [, r])
Return the latitude-mean `A` (mean across dimension `Lat`).
Optionally do the mean for the data in range `r` of that dimension.

This function properly weights the mean by the cosine of the latitude.
"""
function latmean(A::AbDimArray, r = 1:size(A, Lat))
    A = A[Lat(r)]
    lw = _latweights(A)
    dropagg(sum, dimwise(*, A, lw), Lat)
end
# Warning!!! `_latweights` divides by the weight sum, because it is intended to be
# used only with the `sum` function (for A)
_latweights(A::AbDimArray) = _latweights(dims(A, Lat))
function _latweights(A::Lat)
    we = cosd.(Array(A))
    we ./= sum(we)
    return ClimArray(we, (A,))
end

using StatsBase

"""
    spacemean(A::ClimArray [, w]) = spaceagg(mean, A, w)
Average given `A` over its spatial coordinates.
Optionally provide statistical weights in `w`.
"""
spacemean(A, exw=nothing) = spaceagg(mean, A, exw)

"""
    spaceagg(f, A::ClimArray, w = nothing)
Aggregate `A` using function `f` (e.g. `mean`) over all available space (i.e.
longitude and latitude) of `A`, weighting every part of `A` by its spatial area.
The function works for grid as well as equal area space.

`w` can be extra weights, to weight each spatial point with. `w` can either be
just an `AbDimArray` with same space dimensions as `A`, or of exactly same shape as `A`.
"""
spaceagg(f, A::AbDimArray, exw=nothing) = spaceagg(spacestructure(A), f, A, exw)
function spaceagg(::Grid, f, A::AbDimArray, w=nothing)
    wtype = spaceweightassert(A, w)
    cosweights = repeat(cosd.(dims(A, Lat).val)', size(A, Lon))
    # in case Lat precedes Lon, transpose cosine weights
    dimindex(A, Lat) < dimindex(A, Lon) && (cosweights = cosweights')
    other = otherdims(A, (Lon, Lat))
    n = A.name == "" ? "" : A.name*", spatially aggregated with $(string(f))"
    R = ClimArray(zeros(eltype(A), size.(Ref(A), basenameof.(other)), other, n)
    # pre-calculate weights if possible
    if wtype == :no
        W = weights(cosweights)
    elseif wtype == :d2
        W = weights(cosweights .* Array(w))
    end
    # if A is only 2-dimensional, operation is trivial and single number
    if ndims(A) == 2
        return f(A, W)
    end
    # do the weighted average
    for i in otheridxs(A, (Lon, Lat))
        if wtype != :dany
            R[i] = f(view(A, i), W)
        else
            W = weights(view(w, i) .* cosweights)
            R[i] = f(view(A, i), W)
        end
    end
    return R
end

function spaceweightassert(A, w)
    if !isnothing(w)
        wdims = dims(w)
        if length(wdims) == 2
            @assert wdims == dims(A, (Lon, Lat))
            wtype = :d2
        else
            @assert wdims == dims(A)
            wtype = :dany
        end
    else
        wtype = :no
    end
    return wtype
end

spaceagg(::EqArea, f, A, ::Nothing) = dropagg(f, A, Coord)
# I think the best scenario is to modify `dropagg` to take in weights.
function spaceagg(::EqArea, f, A, exw)
    error("TODO")
    w = pweights(Array(exw))
    dropagg(f, A, Coord)
end

#########################################################################
# Hemispheric sum/difference
#########################################################################
export hemispheric_means, hemispheric_functions

"""
    hemispheric_functions(A::ClimArray) → north, south
Return two arrays `north, south`, by splitting `A` to its northern and southern hemispheres,
appropriately translating the latitudes of `south` so that both arrays have the same
latitudinal dimension (and thus can be compared and do opperations between them).
"""
hemispheric_functions(A) = hemispheric_functions(spacestructure(A), A)
function hemispheric_functions(::Grid, A)
    nh = A[Lat(Between(0,  90))]
    sh = A[Lat(Between(-90, 0))]
    # TODO: this can be a function "reverse dim"
    di = dimindex(sh, Lat)
    newdims = [dims(sh)...]
    newdims[di] = dims(nh, Lat)
    data = reverse(Array(sh); dims = di)
    sh = ClimArray(data, (newdims...,))
    return nh, sh
end

function hemispheric_functions(::EqArea, A)
    c = dims(A, Coord)
    @assert issorted(c, by = x -> x[2])
    shi, nhi = hemisphere_indices(c)
    nh = A[Coord(nhi)]
    sh = A[Coord(shi)]
    oldc = Array(dims(sh, Coord))
    si = sortperm(oldc, by = x -> x[2], rev = true)
    newc = [SVector(x[1], abs(x[2])) for x in oldc[si]]
    di = dimindex(sh, Coord)
    newdims = Any[dims(sh)...]
    newdims[di] = Coord(newc)
    data = reverse(Array(sh); dims = di)
    sh = ClimArray(data, (newdims...,))
    return nh, sh
end

"""
    hemispheric_means(A) → nh, sh
Return the (proper) averages of `A` over the northern and southern hemispheres.
Notice that this function explicitly does both zonal as well as meridional averaging.
Use [`hemispheric_functions`](@ref) to just split `A` into two hemispheres.
"""
hemispheric_means(A) = hemispheric_means(spacestructure(A), A)
function hemispheric_means(::Grid, A::AbDimArray)
    @assert hasdim(A, Lat)
    if hasdim(A, Lon)
        B = zonalmean(A)
    else
        B = A
    end
    nh = latmean(B[Lat(Between(0,  90))])
    sh = latmean(B[Lat(Between(-90, 0))])
    return nh, sh
end

function hemispheric_means(::EqArea, A::AbDimArray)
    shi, nhi = hemisphere_indices(c)
    nh = dropagg(mean, A[Coord(nhi)], Coord)
    sh = dropagg(mean, A[Coord(shi)], Coord)
    return nh, sh
end

latitudes(A) = latitudes(spacestructure(A), A)
latitudes(::Grid, A) = Array(dims(A, Lat))
latitudes(::EqArea, A) = unique!([x[2] for x in dims(A, Coord)])
