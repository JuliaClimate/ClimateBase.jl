
function spatialidxs(::EqArea, A)
    return ((Coord(i),) for i in 1:size(A, Coord))
end

function zonalmean(::EqArea, A::AbDimArray)
    idxs, lats = uniquelats(A)
    other = otherdims(A, Coord())
    r = zeros(eltype(A), (length(lats), size.(Ref(A), other)...), other)
    R = ClimArray(r, (Lat(lats), other...), A.name)
    for (i, r) in enumerate(idxs)
        for j in otheridxs(A, Coord())
            R[Lat(i), j...] = mean(view(A, Coord(r), j...))
        end
    end
    return R
end
function zonalmean(::EqArea, A::AbDimArray{T, 1}) where {T}
    idxs, lats = uniquelats(A)
    res = zeros(T, length(lats))
    for (i, r) in enumerate(idxs)
        res[i] = mean(view(A.data, r))
    end
    return ClimArray(res, (Lat(lats),))
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
    sizehint!(lats, round(Int, sqrt(length(c))))
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
    sizehint!(lats, length(lats))
    return idxs, lats
end

spaceagg(::EqArea, f, A, ::Nothing) = dropagg(f, A, Coord)
# I think the best scenario is to modify `dropagg` to take in weights.
function spaceagg(::EqArea, f, A, exw)
    error("TODO")
    w = pweights(Array(exw))
    dropagg(f, A, Coord)
end

function hemispheric_functions(::EqArea, A)
    c = dims(A, Coord).val
    @assert issorted(c; by = x -> x[2])
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

function hemispheric_means(::EqArea, A::AbDimArray)
    shi, nhi = hemisphere_indices(c)
    nh = dropagg(mean, A[Coord(nhi)], Coord)
    sh = dropagg(mean, A[Coord(shi)], Coord)
    return nh, sh
end

latitudes(::EqArea, A) = unique!([x[2] for x in dims(A, Coord)])
