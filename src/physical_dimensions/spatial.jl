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
Works for all types of space (`...` is necessary because `i` is a `Tuple`).
"""
spatialidxs(A::AbDimArray) = spatialidxs(spacestructure(A), A)
function spatialidxs(::LonLatGrid, A)
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
export lon_distance, wrap_lon, longitude_circshift
"""
    wrap_lon(x)
Wrap given longitude to -180 to 180 degrees.
"""
wrap_lon(x) = @. -180 + (360 + ((x+180) % 360)) % 360

"""
    lon_distance(λ1, λ2, Δλ = 360) → δ
Calculate distance `δ` (also in degrees) between longitudes `λ1, λ2`, but taking into
account the periodic nature of longitude, which has period `Δλ = 360ᵒ`.
"""
function lon_distance(x, y, p = eltype(x)(360))
    moddis = mod(abs(x - y), p)
    min(moddis, p - moddis)
end

"""
    longitude_circshift(X::ClimArray [, l]; wrap = true) → Y::ClimArray
Perform the same action as `Base.circshift`, but only for the longitudinal dimension
of `X` with shift `l`. If `wrap = true` the longitudes are wrapped to (-180, 180) degrees
using the modulo operation.

If `l` is not given, it is as much as necessary so that all longitudes > 180 are
wrapped.
"""
function longitude_circshift(X::ClimArray, l = nothing; wrap = true)
    if isnothing(l); l = count(≥(180), dims(X, Lon).val); end
    isnothing(l) && return X
    shifts = map(d -> d isa Lon ? l : 0, dims(X))
    shifted_data = circshift(X.data, shifts)
    shifted_lon = circshift(dims(X, Lon).val, l)
    if wrap; shifted_lon = wrap_lon.(shifted_lon); end
    shifted_lon = vector2range(shifted_lon)
    shifted_dim = Lon(shifted_lon; metadata = dims(X, Lon).metadata)
    new_dims = map(d -> d isa Lon ? shifted_dim : d, dims(X))
    return ClimArray(shifted_data, new_dims; name = X.name, attrib = X.attrib)
end

#########################################################################
# averaging functions over space or time
#########################################################################
using StatsBase

export latmean, spacemean, zonalmean, spaceagg, uniquelats

"""
    zonalmean(A::ClimArray [, W])
Return the zonal mean of `A`. Works for both [`LonLatGrid`](@ref) as well as
[`UnstructuredGrid`](@ref). Optionally provide statistical weights `W`.
These can be the same `size` as `A` or only having the same latitude structure as `A`.
"""
zonalmean(A::AbDimArray, W = nothing) = zonalmean(spacestructure(A), A, W)
zonalmean(::LonLatGrid, A::AbDimArray, W) = dropagg(mean, A, Lon, W)

"""
    latmean(A::ClimArray)
Return the latitude-mean `A` (mean across dimension `Lat`).
This function properly weights by the cosine of the latitude.
"""
function latmean(A::AbDimArray)
    lw = _latweights(A)
    if ndims(A) > 1
        return dropagg(sum, broadcast_dims(*, A, lw), Lat)
    else
        return sum(A .* lw)
    end
end

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
    spacemean(A::ClimArray [, W]) = spaceagg(mean, A, W)
Average given `A` over its spatial coordinates.
Optionally provide statistical weights in `W`.
"""
spacemean(A, exw=nothing) = spaceagg(mean, A, exw)

"""
    spaceagg(f, A::ClimArray, W = nothing)
Aggregate `A` using function `f` (e.g. `mean, std`) over all available space (i.e.
longitude and latitude) of `A`, weighting every part of `A` by its spatial area.

`W` can be extra weights, to weight each spatial point with. `W` can either be a
`ClimArray` with same spatial information as `A`, or having exactly same dimensions as `A`.
"""
spaceagg(f, A::AbDimArray, exw=nothing) = spaceagg(spacestructure(A), f, A, exw)
function spaceagg(::LonLatGrid, f, A::AbDimArray, w=nothing)
    wtype = spaceweightassert(A, w)
    cosweights = repeat(cosd.(dims(A, Lat).val)', size(A, Lon))
    if dimindex(A, Lon) > dimindex(A, Lat)
        error("At the moment this function assumes that Lon preceeds Lat, use `permutdims`.")
    end
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
        r = map(i -> f(view(A, i), weights(view(w, i) .* cosweights)), oidxs)
    end
    return ClimArray(r, other, A.name)
end

function spaceweightassert(A, w)
    if !isnothing(w)
        wdims = dims(w)
        # Check that longitude/latitude match even if the actual dimensions are not
        # identical (because they are e.g. not loaded from the same file)
        if length(wdims) == 2
            @assert Lat ∈ basetypeof.(wdims)
            @assert Lon ∈ basetypeof.(wdims)
            @assert val.(wdims) == val.(dims(A, (Lon, Lat)))
            wtype = :d2
        else
            @assert basetypeof.(wdims) == basetypeof.(dims(A))
            @assert size(w) == size(A)
            wtype = :dany
        end
    else
        wtype = :no
    end
    return wtype
end

#########################################################################
# Hemispheric
#########################################################################
export hemispheric_means, hemispheric_functions

"""
    hemispheric_functions(A::ClimArray) → north, south
Return two arrays `north, south`, by splitting `A` to its northern and southern hemispheres,
appropriately translating the latitudes of `south` so that both arrays have the same
latitudinal dimension (and thus can be compared and do opperations between them).
"""
hemispheric_functions(A) = hemispheric_functions(spacestructure(A), A)
function hemispheric_functions(::LonLatGrid, A)
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
    hemispheric_means(A [,W]) → nh, sh
Return the (proper) averages of `A` over the northern and southern hemispheres.
Notice that this function explicitly does both zonal as well as meridional averaging.
Use [`hemispheric_functions`](@ref) to just split `A` into two hemispheres.

Optionally provide weights `W` that need to have the same structure as [`spaceagg`](@ref).
"""
hemispheric_means(A, args...) = hemispheric_means(spacestructure(A), A, args...)
function hemispheric_means(::LonLatGrid, A::AbDimArray)
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
latitudes(::LonLatGrid, A) = dims(A, Lat).val
#########################################################################
# Tropics/extratropics
#########################################################################
export tropics_extratropics
"""
    tropics_extratropics(A::ClimArray; lower_lat=30, higher_lat=90) → tropics, extratropics
Separate the given array into two arrays: one having latitudes ℓ ∈ [-lower_lat, +lower_lat], and one
having [-higher_lat:-lower_lat, lower_lat:higher_lat].
"""
tropics_extratropics(A, args...; kwargs...) = 
tropics_extratropics(spacestructure(A), A, args...; kwargs...)

function tropics_extratropics(::LonLatGrid, A; lower_lat=30, higher_lat=90)
    tropics = A[Lat(Between(-lower_lat, lower_lat))]
    latdim = dims(A, Lat)
    extra_idxs_sh = DimensionalData.selectindices(latdim, Between(-higher_lat, -lower_lat))
    extra_idxs_nh = DimensionalData.selectindices(latdim, Between(lower_lat, higher_lat))
    extra_idxs = vcat(extra_idxs_sh, extra_idxs_nh)
    extratropics = A[Lat(extra_idxs)]
    return tropics, extratropics
end
