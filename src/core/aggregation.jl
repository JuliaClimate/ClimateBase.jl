#=
Data aggregation of any kind. Includes handling of missing values as well.
=#
#########################################################################
# missing values handling
#########################################################################
export missing_weights, missing_val

function nomissing(da::AbstractArray{Union{T,Missing},N}, check = any(ismissing, da)) where {T,N}
    check && error("array contains missing values")
    return Array{T,N}(da)
end
nomissing(da::AbstractArray{<:Union{Real, Dates.TimeType}}, args...) = da
nomissing(da::ClimArray) = ClimArray(nomissing(da.data), da.dims, da.refdims, da.name, da.attrib)


"""
    missing_weights(A::ClimArray, val = missing_val(A)) → B, W
Generate a new array `B` with values like `A`, but with `A`'s `missing` values replaced
with `val`. Also generate an array of weights, which has the value 0 when `A` had `missing`,
and the value `1` otherwise.

The output of this function should be used in conjunction with any of ClimateBase.jl
aggregating functions like `spacemean, timemean, ...`, when your data have `missing`
values which you want to _completely skip_ during the aggregation process.

This function returns `A, nothing` if `A` has no `missing` values.
"""
function missing_weights(A::ClimArray{Union{T, Missing}}, val = missing_val(A)) where {T}
    B = zeros(T, size(A))
    W = ones(T, size(A))
    missing_idxs = findall(ismissing, A)
    notmissing_idxs = findall(!ismissing, A)
    B[missing_idxs] .= val
    B[notmissing_idxs] .= A.data[notmissing_idxs]
    W[missing_idxs] .= 0
    return ClimArray(B, dims(A); name = A.name, attrib = A.attrib),
           ClimArray(W, dims(A); name = "weights_for_missing")
end
missing_weights(A::ClimArray{<:Number}, val = nothing) = A, nothing

"""
    missing_val(A)
Return the value that represents "missing" data in `A`, according to `A`'s metadata.
If `A` does not have the `_FillValue` metadata, return 0 instead.
"""
function missing_val(A)
    if A.attrib isa Dict
        return get(A.attrib, "_FillValue", 0)
    else
        return 0
    end
end


#########################################################################
# Aggregation of data, dropagg missings, dimensions, etc.
#########################################################################
export dropagg, nomissing, collapse, drop

"""
    dropagg(f, A::ClimArray, d [, W])
Apply statistics/aggregating function `f` (e.g. `sum` or `mean`) on array `A` across
dimension(s) `d` and drop the corresponding dimension(s) from the result
(Julia inherently keeps singleton dimensions).

If `A` is one dimensional, `dropagg` will return the single number of applying `f(A)`.

Optionally you can provide statistical weights in the form of a `W::ClimArray`.
`W` must either have same size as `A` or have only one dimension and satisfy
`length(W) == length(dims(A, d))` (i.e., a weight for each value of the dimension `d`).
The latter case can only work when `d` is a single dimension. See also 
[`missing_weights`](@ref) for (properly) dealing with data that have `missing` values.
"""
function dropagg(f, A, dims)
    length(size(A)) == 1 && return f(A)
    return dropdims(f(A; dims = dims); dims = dims)
end
# this method is necessary because of "reshaping" happening
# in DimensionalData.jl...
function dropagg(f, A::AbDimArray, dims)
    length(size(A)) == 1 && return f(A)
    r = dropdims(f(A; dims = dims); dims = dims)
    DimensionalData.rebuild(r, Array(r.data))
end

dropagg(f, A::ClimArray, d, ::Nothing) = dropagg(f, A, d)

function dropagg(f, A::ClimArray, d, W)
    odims = otherdims(A, d)
    oidxs = otheridxs(A, d)
    D = length(odims)
    if D == 0 # operation output is a single number
        return f(A, weights(W))
    elseif size(A) == size(W)
        r = map(i -> f(view(A, i), weights(view(W, i))), oidxs)
    elseif length(size(W)) == 1 && length(W) == length(dims(A, d))
        fw = weights(W)
        r = map(i -> f(view(A, i), fw), oidxs)
    else
        error("Given weights `W` have invalid form.")
    end
    return ClimArray(r, odims; name = A.name)
end


"""
    collapse(f, A, dim)
Reduce `A` towards dimension `dim` using the collapsing function `f` (e.g. `mean`).
This means that `f` is applied across all other dimensions of `A`, each of which are
subsequently dropped, leaving only the collapsed result of `A` vs. the remaining dimension.
"""
collapse(f, A, dim) = dropagg(f, A, otherdims(A, dim))

#########################################################################
# Other dimensions
#########################################################################
export otherdims
export otheridxs

"""
    otheridxs(A::ClimArray, Dim)
Return an iterator of indices, that when used can access all indices of `A` *except* those
belonging to dimension(s) `Dim`.

For example, if `A` has dims `(Lon, Lat, Time)` you can get all timeseries of `A`:
```julia
for i in otheridxs(A, Time())
    x = A[i...] # this is a timeseries (Vector) for each lon × lat combination
end
```
or all time+latitude slices with
```julia
for i in otheridxs(A, (Time(), Lat()))
    x = A[i...] # matrix of time × latitude slice for each longitude
end
```
(notice that splatting `i...` is necessary because `otheridxs` generates tuples)
"""
function otheridxs(A, D)
    z = otherdims(A, D)
    az = DimensionalData.basetypeof.(z)
    iters = [(Dim(i) for i in 1:size(A, Dim)) for Dim in az]
    return Iterators.product(iters...)
end

# This is the version @rafaqz suggested, but has worse performance...
# otheridxs(A, D) = map(identity, DimensionalData.dimwise_generators(otherdims(A, D)))


#########################################################################
# Dimensionwise
#########################################################################
#=
"""
    broadcast_dims(f, A, B)
Apply bivariate function `f` (e.g. multiplication `*`) across the (first) matching dimension
of `A` and `B`, where `B` has only one dimension.
For normal Julia arrays match is done on the first matching dimension size.

For example, if `A` has dimensions (Lon, Lat, Time) and `B` has dimension (Lat,)
(and of course both Lat are actually the same thing),
then the result would be an array `C` with same structure as `A` but with values
```julia
@. C[:, i, :] = f(A[:, i, :], B[i])
```
for every `i` across the matching dimension.

Useful when wanting to weight a field `A` based on its latitude or time.
"""
function broadcast_dims(f, A::AbDimArray, B::AbDimArray)
    @assert length(size(B)) == 1 "3rd argument of broadcast_dims must be a vector"
    m = find_matching_dim(A, B)
    # Array(data) necessary because reshapes somewhow fail to play well with the code...
    C = broadcast_dims(f, Array(data(A)), Array(data(B)), m)
    return DimensionalData.basetypeof(A)(C, dims(A))
end

broadcast_dims(f, A::AbstractArray, B, m::Int = find_matching_dim(A, B)) =
broadcast_dims!(copy(A), f, A, B, m)

function broadcast_dims!(C, f, A::AbstractArray{T, N}, B::AbstractVector, m::Int) where {T, N}
    shape = ntuple(i -> (i==m) ? length(B) : 1, N)
    b = reshape(B, shape)
    C .= f.(A, b)
    return C
end

function find_matching_dim(A::AbDimArray, B::AbDimArray)::Int
    dimB = dims(B)[1]
    m = findfirst(isequal(dimB), dims(A))
    isnothing(m) && error("A and B have no matching dimensions.")
    return m
end

function find_matching_dim(A::AbstractArray, B::AbstractVector)::Int
    m = findfirst(isequal(length(B)), size(A))
    isnothing(m) && error("A and B have no matching dimensions.")
    return m
end

=#
export broadcast_dims
