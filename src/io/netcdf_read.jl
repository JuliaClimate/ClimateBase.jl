
#########################################################################
# Reading
#########################################################################
"""
    ncread(file, var [, selection]; kwargs...) → A
Load the variable `var` from the `file` and convert it into a [`ClimArray`](@ref)
with proper dimension mapping and also containing the variable attributes as a dictionary.
Dimension attributes are also given to the dimensions of `A`, if any exist.
See keywords below for specifications for unstructured grids.

`file` can be a string to a `.nc` file. Or, it can be an
`NCDataset`, which allows you to lazily combine different
`.nc` data (typically split by time), e.g.
```julia
using Glob # for getting all files
alldata = glob("toa_fluxes_*.nc")
file = NCDataset(alldata; aggdim = "time")
A = ClimArray(file, "tow_sw_all")
```

`var` is a `String` denoting which variable to load.
For `.nc` data containing groups `var` can also be a tuple `("group_name", "var_name")`
that loads a specific variable from a specific group.
In this case, the attributes of both the group and the CF-variable are attributed to
the created [`ClimArray`](@ref).

Optionally you can provide a `selection` for selecting a smaller part of the full array.
The `selection` must be a tuple of indices that compose the selection and you must specify
exactly as many ranges as the dimensions of the array and in the correct order.
For example, if `var` corresponds
to an array with three dimensions, such syntaxes are possible:
```julia
(:, :, 1:3)
(1:5:100, 1:1, [1,5,6])
```
The function [`ncsize`](@ref) can be useful for `selection`.

See also [`ncdetails`](@ref), [`nckeys`](@ref) and [`ncwrite`](@ref).

## Smart loading
The following things make loading data with `ncread` smarter than directly trying to use
NCDatasets.jl and then convert to some kind of dimensional container.
1. Data are directly transformed into `ClimArray`, conserving metadata and dimension names.
1. If there are no missing values in the data (according to CF standards), the
   returned array is automatically converted to a concrete type (i.e. `Union{Float32, Missing}`
   becomes `Float32`).
1. Dimensions that are ranges (i.e. sampled with constant step size) are automatically
   transformed to a standard Julia `Range` type (which makes sub-selecting faster).
1. Automatically deducing whether the spatial information is in an orthogonal
   grid or not, and creating a single `Coord` dimension if not.

## Keywords
* `name` optionally rename loaded array.
* `grid = nothing` optionally specify whether the underlying grid is `grid = OrthogonalSpace()`
  or `grid = CoordinateSpace()`, see [Types of spatial information](@ref).
  If `nothing`, we try to deduce automatically based on
  the names of dimensions and other keys of the `NCDataset`.
* `lon, lat`. These two keywords are useful in unstructured grid data where the grid
  information is provided in a *separate .nc file*. What we need is the user to
  provide vectors of the central longitude and central latitude of each grid point.
  This is done e.g. by
  ```julia
  ds = NCDataset("path/to/grid.nc");
  lon = Array(ds["clon"]);
  lat = Array(ds["clat"]);
  ```
  If `lon, lat` are given, `grid` is automatically assumed `CoordinateSpace()`.
"""
function ncread(path::Union{String, Vector{String}}, args...; kwargs...)
    NCDataset(path) do ds
        data = ncread(ds, args...; kwargs...)
        return data
    end
end

# TODO: Allow reading multiple variables at once. This has the performance benefit
# of not re-creating dimensions all the time.

function ncread(ds::NCDatasets.AbstractDataset, var, selection = nothing;
        name = var2name(var), grid = nothing, lon = nothing, lat = nothing,
    )
    if lon isa Vector && lat isa Vector
        gridtype = CoordinateSpace()
    elseif isnothing(grid)
        gridtype = autodetect_grid(ds)
    else
        gridtype = grid
    end

    if gridtype == CoordinateSpace()
        return ncread_unstructured(ds, var, name, lon, lat, selection)
    else
        return ncread_lonlat(ds, var, name, selection)
    end
end

function autodetect_grid(ds)
    if haskey(ds, "reduced_points") || haskey(ds, "clon") ||
        any(x -> x ∈ ds.dim, POSSIBLE_CELL_NAMES)
        # Common cases of coordinate spaces in NetCDF
        return CoordinateSpace()
    elseif haskey(ds, "lat") && length(size(ds["lat"])) > 1
        # Uncommon case of storing lon/lat as matrices, used in some weird ocean grids
        return CoordinateSpace()
    else
        return OrthogonalSpace()
    end
end

var2name(var) = string(var)
var2name(var::Tuple) = join(var, "_")

#########################################################################
# Reading: OrthogonalSpace and main reading functionality
#########################################################################
# Notice that this function properly loads even without any spatial coordinate
function ncread_lonlat(ds::NCDatasets.AbstractDataset, var, name, selection)
    cfvar = var isa String ? ds[var] : ds.group[var[1]][var[2]]
    sel = isnothing(selection) ? selecteverything(cfvar) : selection
    @assert length(sel) == length(size(cfvar))
    attrib = get_attributes_from_var(ds, cfvar, var)
    A = cfvar[sel...]
    dnames = Tuple(NCDatasets.dimnames(cfvar))
    if !any(ismissing, A)
        A = nomissing(A)
    end
    dimensions = create_dims(ds, dnames, sel)
    data = ClimArray(A, dimensions; name = Symbol(name), attrib = attrib)
    return data
end
selecteverything(cfvar) = map(i -> 1:i, size(cfvar))

get_attributes_from_var(ds, cfvar, var::String) = Dict(cfvar.attrib)
function get_attributes_from_var(ds, cfvar, var::Tuple)
    d1 = Dict(cfvar.attrib)
    d2 = Dict(ds.group[var[1]].attrib)
    return merge!(d2, d1)
end


"""
    create_dims(ds::NCDatasets.AbstractDataset, dnames, cfvar, sel = selecteverything(A))
Create a tuple of `Dimension`s from the `dnames` (tuple of strings).
"""
function create_dims(ds::NCDatasets.AbstractDataset, dnames, sel = selecteverything(A))
    true_dims = to_proper_dimensions(dnames)
    dim_values = extract_dim_values(ds, dnames, sel)
    # Some stupid datasets return a union{Missing} type for dimensions.
    # this is of course nonsense, a dimension cannot have "missing" values.
    if any(d -> Missing <: eltype(d), dim_values)
        dim_values = nomissing.(dim_values)
    end
    optimal_values = vector2range.(dim_values)
    attribs = [
        (haskey(ds, d) && ds[d].attrib isa NCDatasets.BaseAttributes) ?
            Dict(ds[d].attrib) :
            Dict()
        for d in dnames
    ]
    out = Dimension[]
    for i in 1:length(true_dims)
        push!(out, true_dims[i](optimal_values[i]; metadata = attribs[i]))
    end
    return (out...,)
end

function to_proper_dimensions(dnames)
    r = []
    for n in dnames
        if haskey(COMMONNAMES, n)
            push!(r, COMMONNAMES[n])
        else
            @warn """
            Dimension name "$n" not in common names.
            Please consider opening an issue/PR on GitHub proposing "$n" in common names!
            Making a generic dimension named `Dim{:$n}` for now...
            """
            push!(r, Dim{Symbol(n)})
        end
    end
    return (r...,)
end
export Dim # for generic dimensions this must be exported

function extract_dim_values(ds::NCDatasets.AbstractDataset, dnames, sel)
    dim_values = AbstractVector[]
    for (i, n) in enumerate(dnames)
        if haskey(ds, n)
            push!(dim_values, vec(ds[n][sel[i]]))
        else
            @warn """
            Dimension named "$n" does not have values in the dataset.
            Using the range `1:size(cfvar, i)` as the dimension values instead,
            where `i` is the dimension index.
            """
            push!(dim_values, sel[i])
        end
    end
    return dim_values
end


#########################################################################
# Reading: CoordinateSpace
#########################################################################
export SVector

function ncread_unstructured(
        ds::NCDatasets.AbstractDataset, var::String, name, lon, lat, selection
    )
    cfvar = var isa String ? ds[var] : ds.group[var[1]][var[2]]
    sel = isnothing(selection) ? selecteverything(cfvar) : selection
    attrib = get_attributes_from_var(ds, cfvar, var)

    # Here we generate the longitudes and latitudes based on whether we have
    # reduced points or not, and obtain the name of the coord dimension in the `ds`
    if isnothing(lon)
        lonlat, original_grid_dim = load_coordinate_points(ds)
    else
        lonlat = [SVector(lo, la) for (lo, la) in zip(lon, lat)]
        original_grid_dim = intersect(NCDatasets.dimnames(cfvar), POSSIBLE_CELL_NAMES)[1]
    end

    alldims = Any[NCDatasets.dimnames(cfvar)...]

    # Set up the remaining dimensions of the dataset
    if original_grid_dim isa Tuple
        # stupid case where `lon` data are `Matrix`

        # TODO: I'm not sure I've set up this to work correctly in the case
        # of `selection` given by user... But the `Matrix` case is so rare
        # I honestly don't care at the moment.

        A = cfvar[sel...]
        sizes = [size(A)...]
        # First, find where lon, lat dimensions are positioned
        is = findall(x -> x ∈ original_grid_dim, alldims)
        @assert length(is) == 2
        remainingdims = deleteat!(copy(alldims), is)
        actualdims = Any[create_dims(ds, remainingdims, A)...]
        # Then, reshape A accordingly
        sizes[is[1]] = sizes[is[1]]*sizes[is[2]]
        deleteat!(sizes, is[2])
        A = reshape(A, sizes...)
        i = is[1]
    else
        @assert original_grid_dim ∈ alldims
        i = findfirst(x -> x == original_grid_dim, alldims)
        remainingdims = deleteat!(copy(alldims), i)
        A = cfvar[sel...]
        remainingsel = Tuple(deleteat!([sel...], i))
        actualdims = Any[create_dims(ds, remainingdims, remainingsel)...]
    end

    # Make coordinate dimension
    lonlat = lonlat[sel[i]]
    si = sortperm(lonlat; by = reverse)
    coords = Coord(lonlat, (Lon, Lat))
    insert!(actualdims, i, coords)

    # Create array, sort by latitude and remove missings
    X = ClimArray(A, Tuple(actualdims))
    X = X[Coord(si)]
    if !any(ismissing, X)
        X = nomissing(X)
    end
    return ClimArray(X; name = Symbol(name), attrib)
end

function has_unstructured_key(ds)
    any(x -> haskey(ds.dim, x), POSSIBLE_CELL_NAMES) ||
    haskey(ds, "lon") || haskey(ds, "clon")
end

function load_coordinate_points(ds)
    if haskey(ds, "reduced_points")
        lonlat = reduced_grid_to_points(ds["lat"], ds["reduced_points"])
        original_grid_dim = "rgrid" # Specific to CDO Gaussian grid
    elseif haskey(ds, "lon") && length(size(ds["lon"])) == 2
        # This is the case of having a non-orthogonal grid, but
        # still saving the lon/lat information as matrices whose dimensions are lon/lat
        # (this is a rather stupid way to format things, unfortunately)
        lons = ds["lon"] |> Matrix |> vec
        lats = ds["lat"] |> Matrix |> vec
        original_grid_dim = ("lon", "lat")
        lonlat = [SVector(lo, la) for (lo, la) in zip(lons, lats)]
    elseif has_unstructured_key(ds)
        if haskey(ds, "lon")
            lons = ds["lon"] |> Vector
            lats = ds["lat"] |> Vector
            original_grid_dim = NCDatasets.dimnames(ds["lon"])[1]
        elseif haskey(ds, "clon")
            lons = ds["clon"] |> Vector
            lats = ds["clat"] |> Vector
            original_grid_dim = NCDatasets.dimnames(ds["clon"])[1]
        else
            error("""
            We didn't find key `"lon"` or `"clon"` that represents the longitude of each
            polygon in a non-orthogonal grid.
            """)
        end
        lonlat = [SVector(lo, la) for (lo, la) in zip(lons, lats)]
    else
        error("""
        We couldn't automatically identify the lon/lat values of cell centers.
        Please provide explicitly keywords `lon, lat` in `ncread`.
        """)
    end
    lonlat = convert_to_degrees(lonlat, ds)
    return lonlat, original_grid_dim
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

function convert_to_degrees(lonlat, ds)
    !(haskey(ds, "clat") || haskey(ds, "lat")) && return convert_to_degrees(lonlat)
    x = haskey(ds, "clat") ? ds["clat"] : ds["lat"]
    if get(x.attrib, "units", nothing) == "radian" || !any(ll -> abs(ll[2]) > π/2, lonlat)
        lonlat = [SVector(lo*180/π, la*180/π) for (lo, la) in lonlat]
    end
    return lonlat
end
function convert_to_degrees(lonlat)
    if !any(ll -> abs(ll[2]) > π/2, lonlat)
        lonlat = [SVector(lo*180/π, la*180/π) for (lo, la) in lonlat]
    end
    return lonlat
end
