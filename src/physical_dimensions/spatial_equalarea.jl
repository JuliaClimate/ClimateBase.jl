#########################################################################
# Spatial functions
#########################################################################
function spatialidxs(::UnstructuredGrid, A)
    return ((Coord(i),) for i in 1:size(A, Coord))
end

function zonalmean(::UnstructuredGrid, A::AbDimArray)
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
function zonalmean(::UnstructuredGrid, A::AbDimArray{T, 1}) where {T}
    idxs, lats = uniquelats(A)
    res = zeros(T, length(lats))
    for (i, r) in enumerate(idxs)
        res[i] = mean(view(A.data, r))
    end
    return ClimArray(res, (Lat(lats),); name=A.name, attrib=A.attrib)
end


"""
    uniquelats(A::AbDimArray) → idxs, lats
    uniquelats(c::Vector{<:AbstractVector}) → idxs, lats
Find the unique latitudes of `A`. Return the indices (vector of ranges) that each latitude
in `lats` covers, as well as the latitudes themselves.
"""
uniquelats(A::AbDimArray) = uniquelats(dims(A, Coord).val)
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

spaceagg(::UnstructuredGrid, f, A, ::Nothing) = dropagg(f, A, Coord)
# I think the best scenario is to modify `dropagg` to take in weights.
function spaceagg(::UnstructuredGrid, f, A, exw)
    error("TODO")
    w = pweights(Array(exw))
    dropagg(f, A, Coord)
end

function hemispheric_functions(::UnstructuredGrid, A)
    c = dims(A, Coord).val
    shi, nhi = hemisphere_indices(c)
    nh = A[Coord(nhi)]
    sh = A[Coord(shi)]
    oldc = dims(sh, Coord).val
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

function hemisphere_indices(c)
    idxs, lats = uniquelats(c)
    i = findfirst(x -> x > 0, lats)
    shi = idxs[1][1]:idxs[i][end]
    nhi = idxs[i+1][1]:idxs[end][end]
    return nhi, shi
end

latitudes(::UnstructuredGrid, A) = unique!([x[2] for x in dims(A, Coord)])

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
function DimensionalData.sel2indices(c::Coord, sel::Tuple{<:Lat{ <: Between}})
    l1, l2 = sel[1].val.val
    return coord_latitudes_between(c.val, l1, l2) # this is Vector{Int}
end
