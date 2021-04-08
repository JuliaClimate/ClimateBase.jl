#########################################################################
# Loading equal area
#########################################################################
export SVector, ClimArray_eqarea

ClimArray_eqarea(s::String, args...) = ClimArray_eqarea(NCDataset(s), args...)
function ClimArray_eqarea(ds::NCDatasets.AbstractDataset, var::String, name = var)
    svar = string(var)
    cfvar = ds[svar]
    attrib = Dict(cfvar.attrib)
    A = cfvar |> Array
    # if haskey(ds, "ncells") # this is the equal area grid, so we make a Coord dimension
    #     # TODO: This code has not yet been generalized to arbitrary dimensions
    #     lon = ds["lon"] |> Array .|> wrap_lon
    #     lat = ds["lat"] |> Array
    #     time = ds["time"] |> Array
    #     lonlat = [SVector(lon[i], lat[i]) for i in 1:length(lon)]
    #     # here we sort lonlat and A in ascending latitude order,
    #     # because the CDO output has reverse or even totally unsorted order
    #     si = sortperm(lonlat, by = reverse)
    #     data = ClimArray(A[si, :], (Coord(lonlat[si]), Time(time));
    #     attrib = attrib, name = svar)

    if !haskey(ds, "reduced_points")
        error("""
        We expected a key `"reduced_points"` in the .nc file, which contains the information
        to reconstruct the Gaussian equal area points.
        """)
    end

    # TODO: I've noticed that this converts integer dimension (like pressure)
    # into Float64, but I'm not sure why...
    alldims = [NCDatasets.dimnames(cfvar)...]
    @assert "rgrid" ∈ alldims
    i = findfirst(x -> x == "rgrid", alldims)
    remainingdims = deleteat!(copy(alldims), i)
    actualdims = Any[create_dims(ds, remainingdims)...]

    lonlat = reduced_grid_to_points(ds["lat"], ds["reduced_points"])
    si = sortperm(lonlat, by = reverse)
    coord_text = "Gaussian equal area grid with $(length(ds["lat"])) points."
    coords = Coord(lonlat; metadata = Dict("grid" => coord_text))

    insert!(actualdims, i, coords)

    X = ClimArray(A, Tuple(actualdims))
    X = X[Coord(si)]
    if !any(ismissing, X)
        X = nomissing(X)
    end
    return ClimArray(X; name = Symbol(name), attrib)
end

function reduced_grid_to_points(lat, reduced_points)
    lonlat = SVector{2, Float32}[]
    for (i, θ) in enumerate(lat)
        n = reduced_points[i]
        dλ = Float32(360/n)
        for j in 0:n-1
            push!(lonlat, SVector(0 + dλ*j, θ))
        end
    end
    return lonlat
end

#########################################################################
# Spatial functions
#########################################################################
function spatialidxs(::GaussianEqualArea, A)
    return ((Coord(i),) for i in 1:size(A, Coord))
end

function zonalmean(::GaussianEqualArea, A::AbDimArray)
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
function zonalmean(::GaussianEqualArea, A::AbDimArray{T, 1}) where {T}
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

spaceagg(::GaussianEqualArea, f, A, ::Nothing) = dropagg(f, A, Coord)
# I think the best scenario is to modify `dropagg` to take in weights.
function spaceagg(::GaussianEqualArea, f, A, exw)
    error("TODO")
    w = pweights(Array(exw))
    dropagg(f, A, Coord)
end

function hemispheric_functions(::GaussianEqualArea, A)
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

function hemispheric_means(::GaussianEqualArea, A::AbDimArray)
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

latitudes(::GaussianEqualArea, A) = unique!([x[2] for x in dims(A, Coord)])

#########################################################################
# Extention of convenience indexing of `Coord`
#########################################################################
# TODO: this section can be made much more general, but with much more
# effort: https://github.com/rafaqz/DimensionalData.jl/issues/207
export coord_latitudes_between
function coord_latitudes_between(A::ClimArray, l1, l2)
    coord_latitudes_between(dims(A, Coord).val, l1, l2)
end

function coord_latitudes_between(c, l1, l2)
    idxs, lats = uniquelats(c)
    # Notice that lats is guaranteed sorted for gaussian equal area
    i1 = searchsortedfirst(lats, l1)
    i2 = searchsortedlast(lats, l2)
    i1 > i2 && ((i1, i2) = (i2, i1)) # in case bounds are given in reverse order
    return idxs[i1][1]:idxs[i2][end]
end

# This method is necessary so that I can do A[Coord(Lat(...))]
function DimensionalData._dims2indices(c::Coord, sel::Coord{ <: Lat{ <: Between}}, emptyval = Colon())
    l1, l2 = sel.val.val.val
    return coord_latitudes_between(c, l1, l2) # this is Vector{Int}
end

# This allows Between to access Coord directly and be translated to latitude
function DimensionalData.sel2indices(c::Coord, sel::Between{Tuple{X,Y}}) where {X<:Real, Y<:Real}
    l1, l2 = sel.val
    return coord_latitudes_between(c, l1, l2) # this is Vector{Int}
end
