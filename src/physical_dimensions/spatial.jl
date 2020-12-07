#=
Functionality regarding spatial coordinates of an array, i.e. Longitude and Lattitude
as well as equal area grid handling and anything else related to space.
=#
using Statistics, StatsBase, StaticArrays
export SVector # for equal area grid
#########################################################################
# Spatial indexing
#########################################################################
export spatialidxs, lonlatfirst

"""
    spatialidxs(A::ClimArray) → idxs
Return an iterable that can be used to access all spatial points of `A` with the syntax
```julia
idxs = spatialidxs(A)
for i in idxs
    slice_at_give_space_point = A[i...]
end
```
Works for standard grid as well as equal area (`...` necessary because `i` is a `Tuple`).
"""
spatialidxs(A::AbDimArray) = spatialidxs(spacestructure(A), A)
function spatialidxs(::Grid, A)
    lons = (Lon(i) for i in 1:size(A, Lon))
    lats = (Lat(i) for i in 1:size(A, Lat))
    return Iterators.product(lons, lats)
end

"""
    lonlatfirst(A::ClimArray, args...) → B
Permute the dimensions of `A` to make a new array `B` that has first dimension longitude,
second dimension latitude, with the remaining dimensions of `A` following
(useful for most plotting functions). Optional extra dimensions
can be given as `args...`, specifying a specific order for the remaining dimensions.

Example:
```julia
B = lonlatfirst(A)
C = lonlatfirst(A, Time)
```
"""
function lonlatfirst(C, args...)
    permutedims(C, (Lon, Lat, args..., otherdims(C, (Lon, Lat, args...))...))
end

#########################################################################
# Periodicity of longitude
#########################################################################
export lon_distance, wrap_lon
"""
    wrap_lon(x)
Wrap given longitude to -180 to 180 degrees.
"""
wrap_lon(x) = @. -180 + (360 + ((x+180) % 360)) % 360

"""
    lon_distance(λ1, λ2, Δλ = 360) → δ
Calculate distance `δ` (also in degrees) between longitudes `λ1, λ2`, but taking into
account the periodic nature of longitude, which has period 360ᵒ.
"""
function lon_distance(x, y, p = eltype(x)(360))
    moddis = mod(abs(x - y), p)
    min(moddis, p - moddis)
end

#########################################################################
# averaging functions over space or time
#########################################################################
using StatsBase

export latmean, spacemean, zonalmean, spaceagg, uniquelats

"""
    zonalmean(A::ClimArray)
Return the zonal mean of `A`.
Optionally do the mean for the data in range `r` of the longitude
(`r` is fed into the dimension so it can be A range or an arbitrary selector).

Works for both grid and equal area space.
"""
zonalmean(A::AbDimArray, args...) = zonalmean(spacestructure(A), A, args...)
zonalmean(::Grid, A::AbDimArray) = dropagg(mean, A, Lon)

"""
    latmean(A::ClimArray [, r])
Return the latitude-mean `A` (mean across dimension `Lat`).
Optionally do the mean for the data in range `r` of that dimension.

This function properly weights the mean by the cosine of the latitude.
"""
function latmean(A::AbDimArray)
    lw = _latweights(A)
    if ndims(A) > 1
        return dropagg(sum, dimwise(*, A, lw), Lat)
    else
        return sum(A .* lw)
    end
end
latmean(A::AbstractDimArray, r) = latmean(A[Lat(r)])


# Warning!!! `_latweights` divides by the weight sum, because it is intended to be
# used only with the `sum` function (for A)
_latweights(A::AbDimArray) = _latweights(dims(A, Lat))
function _latweights(A::Lat)
    we = cosd.(A.val)
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
just an `AbDimArray` with same space as `A`, or of exactly same shape as `A`.
"""
spaceagg(f, A::AbDimArray, exw=nothing) = spaceagg(spacestructure(A), f, A, exw)
function spaceagg(::Grid, f, A::AbDimArray, w=nothing)
    wtype = spaceweightassert(A, w)
    cosweights = repeat(cosd.(dims(A, Lat).val)', size(A, Lon))
    # TODO: Extends so that this assertion is not necessary:
    @assert dimindex(A, Lon) < dimindex(A, Lat) "longitude must precede latitude"
    other = otherdims(A, (Lon, Lat))
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
    oidxs = otheridxs(A, (Lon(), Lat()))
    if wtype != :dany
        r = map(i -> f(view(A, i), W), oidxs)
    else
        # TODO: This multiplication .* here assumes that the lon-lat grid is the
        # first dimension.
        r = map(i -> f(view(A, i), weights(view(w, i) .* cosweights)), oidxs)
    end
    return ClimArray(r, other, A.name)
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

latitudes(A) = latitudes(spacestructure(A), A)
latitudes(::Grid, A) = dims(A, Lat).val
