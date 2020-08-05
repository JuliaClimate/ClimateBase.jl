#=
Data aggregation of any kind
=#

#########################################################################
# Aggregation of data, dropagg missings, dimensions, etc.
#########################################################################
export dropagg, nomissing, collapse, drop

"""
    dropagg(f, A, dims)
Apply aggregating function `f` (e.g. `sum`) on array `A` across dimension(s) `dims` and drop
the corresponding dimension(s) from the result (Julia inherently keeps singleton dimensions).

If `A` is one dimensional, `dropagg` will return the single number of applying `f(A)`.
"""
function dropagg(f, A, dims)
    length(size(A)) == 1 && return f(A)
    return dropdims(f(A; dims = dims); dims = dims)
end

"""
    collapse(f, A, dim)
Reduce `A` towards dimension `dim` using the collapsing function `f` (e.g. `mean`).
This means that `f` is applied across all other dimensions of `A`, each of which are
subsequently dropped, leaving only the collapsed result of `A` vs. the remaining dimension.
"""
function collapse(f, A, dim)
    di = dimindex(A, dim)
    dimsA = dims(A)
    touse = setdiff(eachindex(dimsA), di)
    return dropagg(f, A; dims = ntuple(i -> dimsA[touse[i]], length(touse)))
end

dimindex(A, i::Int) = i

function nomissing(da::AbstractArray{Union{T,Missing},N}) where {T,N}
    any(ismissing, da) && error("array contains missing values")
    return Array{T,N}(da)
end
nomissing(a::AbstractArray, args...) = a
function nomissing(da::Array{Union{T,Missing},N}, value) where {T,N}
    return replace(da, missing => T(value))
end

#########################################################################
# Extend some DimensionlArray stuff
#########################################################################
export otherdims

nomissing(da::AbDimArray) = DimensionalData.rebuild(da, nomissing(da.data))

# this method is necessary because of "reshaping" happening
# in DimensionalData.jl...
function dropagg(f, A::AbDimArray, dims)
    length(size(A)) == 1 && return f(A)
    r = dropdims(f(A; dims = dims); dims = dims)
    DimensionalData.rebuild(r, Array(r.data))
end

function dimindex(A::AbDimArray, Dim)
    @assert hasdim(A, Dim)
    return findfirst(x -> typeof(x) <: Dim, dims(A))
end

Base.ones(A::AbDimArray) = DimensionalArray(ones(size(A)), dims(A))

function otherdims(A, D)
    x = dims(A)
    n = dimindex(A, D)
    ntuple(i -> i < n ? x[i] : x[i + 1], length(x) - 1)
end


# TODO: Better name
"""
    alongdimidxs(A::AbDimArray, D)
Return an iterator of indices, that when used to access `A` will provide the variation of
`A` along the given dimension `D`, for every other value of the remaining dimensions.

For example, if `A` has dims `(Lon, Lat, Time)` and you can get all timeseries of `A`
for every possible location doing
```julia
for i in alongdimidxs(A, Time)
    x = A[i...] # this is a timeseries (Vector)
end
```
"""
function alongdimidxs(A, D)
    z = otherdims(A, D)
    az = DimensionalData.basetypeof.(z)
    iters = [(Dim(i) for i in 1:size(A, Dim)) for Dim in az]
    return Iterators.product(iters...)
end

export alongdimidxs

#########################################################################
# Dimensionwise
#########################################################################
"""
    dimwise(f, A, B)
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
function dimwise(f, A::AbDimArray, B::AbDimArray)
    @assert length(size(B)) == 1 "3rd argument of dimwise must be a vector"
    m = find_matching_dim(A, B)
    # Array(data) necessary because reshapes somewhow fail to play well with the code...
    C = dimwise(f, Array(data(A)), Array(data(B)), m)
    return DimensionalArray(C, dims(A))
end

dimwise(f, A::AbstractArray, B, m::Int = find_matching_dim(A, B)) =
dimwise!(copy(A), f, A, B, m)

function dimwise!(C, f, A::AbstractArray{T, N}, B::AbstractVector, m::Int) where {T, N}
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

export dimwise
