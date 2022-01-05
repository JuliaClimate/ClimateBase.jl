"""
    transform_to_coord(A::ClimArray) → B
Transform given `A` to a new `B::ClimArray` so that the `Lon, Lat` dimensions in `A`
are transformed to a `Coord` dimension in `B`.
"""
function transform_to_coord(A; name = A.name, attrib = A.attrib)
    hasdim(A, Coord) && return A
    @assert hasdim(A, Lon) && hasdim(A, Lat)
    
    # First make coord dimension
    londim = dims(A, Lon)
    latdim = dims(A, Lat)
    lonlat = [SVector(l, lat) for lat in latdim for l in londim]
    coorddim = Coord(lonlat, (Lon, Lat))
    
    # Then, reshape A to have this dimension
    is = dimindex(A, (Lon, Lat))
    sizes = [size(A)...]
    sizes[is[1]] = sizes[is[1]]*sizes[is[2]]
    deleteat!(sizes, is[2])
    B = reshape(copy(A.data), sizes...)
    i = is[1]

    # Then, make new dimensions
    newdims = [dims(A)...]
    deleteat!(newdims, is)
    insert!(newdims, i, coorddim)
    X = ClimArray(B, Tuple(newdims))

    # Don't forget to sort lonlat coordinates
    si = sortperm(lonlat, by = reverse)
    X = X[Coord(si)]
    return ClimArray(X; name = Symbol(name), attrib)
end


#########################################################################
# Spatial aggregation functions
#########################################################################
spaceagg(::UnstructuredGrid, f, A, W = nothing) = dropagg(f, A, Coord, W)

function spatialidxs(::UnstructuredGrid, A)
    return ((Coord(i),) for i in 1:size(A, Coord))
end

function zonalmean(::UnstructuredGrid, A::AbDimArray, ::Nothing)
    idxs, lats = uniquelats(A)
    other = otherdims(A, Coord())
    r = zeros(eltype(A), (length(lats), size.(Ref(A), other)...))
    R = ClimArray(r, (Lat(lats), other...); name=A.name, attrib=A.attrib)
    for (i, r) in enumerate(idxs)
        for j in otheridxs(A, Coord())
            R[Lat(i), j...] = mean(view(A, Coord(r), j...))
        end
    end
    return R
end
function zonalmean(::UnstructuredGrid, A::AbDimArray{T, 1}, ::Nothing) where {T}
    idxs, lats = uniquelats(A)
    res = zeros(T, length(lats))
    for (i, r) in enumerate(idxs)
        res[i] = mean(view(A.data, r))
    end
    return ClimArray(res, (Lat(lats),); name=A.name, attrib=A.attrib)
end

# zonal mean with weights
function zonalmean(::UnstructuredGrid, A::ClimArray, W::AbstractArray)
	@assert size(A) == size(W)
    idxs, lats = uniquelats(A)
    other = otherdims(A, Coord())
    r0 = zeros(eltype(A), (length(lats), size.(Ref(A), other)...))
    R = ClimArray(r0, (Lat(lats), other...); name=A.name, attrib=A.attrib)
    for (i, r) in enumerate(idxs)
        for j in otheridxs(A, Coord())
            R[Lat(i), j...] = mean(view(A, Coord(r), j...), weights(view(W, Coord(r), j...)))
        end
    end
    return R
end
function zonalmean(::UnstructuredGrid, A::ClimArray{T, 1}, W::AbstractArray) where {T}
	@assert size(A) == size(W)
    idxs, lats = uniquelats(A)
    res = zeros(T, length(lats))
    for (i, r) in enumerate(idxs)
        res[i] = mean(view(A.data, r), weights(view(W.data, r)))
    end
    return ClimArray(res, (Lat(lats),); name=A.name, attrib=A.attrib)
end


"""
    uniquelats(A::ClimArray) → idxs, lats
    uniquelats(c::Vector{<:AbstractVector}) → idxs, lats
Find the unique latitudes of `A`. Return the indices (vector of ranges) that each latitude
in `lats` covers, as well as the latitudes themselves.
"""
uniquelats(A::AbDimArray) = uniquelats(gnv(dims(A, Coord)))
function uniquelats(c)
    @assert issorted(c; by = x -> x[2])
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

#########################################################################
# Special latitude splitting
#########################################################################
function hemispheric_functions(::UnstructuredGrid, A)
    c = gnv(dims(A, Coord))
    nhi, shi = hemisphere_indices(c)
    nh = A[Coord(nhi)]
    sh = A[Coord(shi)]
    oldc = gnv(dims(sh, Coord))
    si = sortperm(oldc, by = x -> x[2], rev = true)
    newc = [SVector(x[1], abs(x[2])) for x in oldc[si]]
    di = dimindex(sh, Coord)
    newdims = Any[dims(sh)...]
    newdims[di] = Coord(newc)
    data = reverse(Array(sh); dims = di)
    sh = ClimArray(data, (newdims...,))
    return nh, sh
end

function hemispheric_means(::UnstructuredGrid, A::AbDimArray)
    nhi, shi = hemisphere_indices(A)
    nh = dropagg(mean, A[Coord(nhi)], Coord)
    sh = dropagg(mean, A[Coord(shi)], Coord)
    return nh, sh
end

function hemispheric_means(::UnstructuredGrid, A::AbDimArray, W)
	@assert size(A) == size(W)
    nhi, shi = hemisphere_indices(A)
    nh = dropagg(mean, A[Coord(nhi)], Coord, W[Coord(nhi)])
    sh = dropagg(mean, A[Coord(shi)], Coord, W[Coord(shi)])
    return nh, sh
end

"""
	hemisphere_indices(coords) → nhi, shi
Return the indices of coordinates belonging to the north and south hemispheres.
"""
function hemisphere_indices(c)
    idxs, lats = uniquelats(c)
    i = findfirst(x -> x > 0, lats)
    shi = idxs[1][1]:idxs[i-1][end]
    nhi = idxs[i][1]:idxs[end][end]
    return nhi, shi
end

latitudes(::UnstructuredGrid, A) = unique!([x[2] for x in dims(A, Coord)])

function tropics_extratropics(::UnstructuredGrid, A; lower_lat=30)
    # TODO: Support `higher_lat` keyword
    c = gnv(dims(A, Coord))
    idxs, lats = uniquelats(c)
    i1 = findlast(x -> x < -lower_lat, lats)
    i2 = findfirst(x -> x > lower_lat, lats)
    # tropics indices (accounting for hemispheric only data as well)
    t1 = isnothing(i1) ? 0 : i1
    t2 = isnothing(i2) ? length(idxs)+1 : i2
    i_tropics = idxs[t1+1][1]:idxs[t2-1][end]
    # extratropics indices (accounting for hemispheric only data as well)
    i_sh_extra = isnothing(i1) ? (1:0) : idxs[1][1]:idxs[i1][end]
    i_nh_extra = isnothing(i2) ? (1:0) : idxs[i2][1]:idxs[end][end]
    i_extra = vcat(i_sh_extra, i_nh_extra)

    tropics = A[Coord(i_tropics)]
    extratropics = A[Coord(i_extra)]
    return tropics, extratropics
end

#########################################################################
# Extention of convenience indexing of `Coord`
#########################################################################
# The code here is exclusively a performance optimization that relies
# on the fact that we sort coordinates by latitude
function coord_latitudes_between(c, l1, l2)
    idxs, lats = uniquelats(c)
    i1 = searchsortedfirst(lats, l1)
    i2 = searchsortedlast(lats, l2)
    i1 > i2 && ((i1, i2) = (i2, i1)) # in case bounds are given in reverse order
    return idxs[i1][1]:idxs[i2][end]
end

# This modifies what happens on A[Coord(Lat(Between(x,y)))]
function DimensionalData.selectindices(c::Coord, sel::Tuple{<:Lat{ <: Between}})
    l1, l2 = sel[1].val.val
    return coord_latitudes_between(gnv(c), l1, l2) # this is Vector{Int}
end

# This modifies what happens on A[Coord(Lat(x..y))]
function DimensionalData.selectindices(c::Coord,
    sel::Tuple{<:Lat{ <: DimensionalData.LookupArrays.IntervalSets.AbstractInterval}})
    l1 = sel[1].val.left; l2 = sel.val.right
    return coord_latitudes_between(gnv(c), l1, l2) # this is Vector{Int}
end