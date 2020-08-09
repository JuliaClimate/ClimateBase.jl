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
    spatialidxs(a::DimensionalArray) → idxs
Return an iterable that can be used to access all spatial points of `a` with the syntax
```julia
idxs = spatialidxs(a)
for i in idxs
    slice_at_give_space_point = a[i...]
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
    return ((Coord(i),) for i in 1:size(A, Coord))
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

# TODO: Document what it means to be Coord space. We expect that the coordinates
# are sorted by latitude

"""
    zonalmean(a::DimensionalArray [, r])
Return the zonal mean of `a`.
Optionally do the mean for the data in range `r` of the longitude
(`r` is fed into the dimension so it can be a range or an arbitrary selector).

Works for both grid as well as equal-area spaces.
"""
zonalmean(a::AbDimArray) = dropagg(mean, a, Lon)
zonalmean(a::AbDimArray, r) = dropagg(mean, a[Lon(r)], Lon)

# TODO: Extend this for `a` with more than 2 dimensions
# possible way is to permute dims so that Lat is first dim, and let the rest dims
# just be propagated?
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
    dropagg(sum, dimwise(*, a, lw), Lat)
end
# Warning!!! `_latweights` divides by the weight sum, because it is intended to be
# used only with the `sum` function (for a)
_latweights(a::AbDimArray) = _latweights(dims(a, Lat))
function _latweights(a::Lat)
    we = cosd.(Array(a))
    we ./= sum(we)
    return ClimArray(we, (a,))
end

using StatsBase

"""
    spacemean(a::DimensionalArray, w=nothing)
Average given `a` over its spatial coordinates.
Optionally provide statistical weights in `w`.
"""
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
function spaceagg(::Grid, f, a, exw=nothing)
    # This assumes that lon is first dim and lat is second dim.
    w = repeat(cosd.(Array(dims(a, Lat)))', length(dims(a, Lon)))
    # TODO: Extend this for abitrary matching weights.
    # We can use matching dims...?
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

#########################################################################
# Hemispheric sum/difference
#########################################################################
export hemispheric_means, hemispheric_functions

function even_odd_decomp(A)
    nh, sh = hemispheric_means(A)
    return (nh .+ sh)/2, (nh .- sh)/2
end

function even_odd_functions(A)
    nh, sh = hemispheric_functions(A)
    return (nh .+ sh)/2, (nh .- sh)/2
end

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
    sh = DimensionalArray(data, (newdims...,))
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
    sh = DimensionalArray(data, (newdims...,))
    return nh, sh
end

"""
    hemispheric_means(A) → nh, sh
Return the (proper) averages of `A` over the north and south hemispheres.
Notice that this function explicitly does both zonal as well as meridional averaging.
Use [`hemispheric_functions`](@ref) to just split `A` into two hemispheres.
"""
function hemispheric_means(A::AbDimArray)
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

latitudes(A) = latitudes(spacestructure(A), A)
latitudes(::Grid, A) = Array(dims(A, Lat))
latitudes(::EqArea, A) = unique!([x[2] for x in dims(A, Coord)])
