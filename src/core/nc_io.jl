#=
Code related with input output (IO) of .nc files directly to/from ClimArrays
An initial version of parts of this code was taken from:
https://github.com/rafaqz/GeoData.jl
=#
using NCDatasets: NCDatasets, NCDataset
export NCDataset
export nckeys, ncdetails, globalattr
export ncread, ncwrite

dim_to_commonname(::Lat) = "lat"
dim_to_commonname(::Lon) = "lon"
dim_to_commonname(::Time) = "time"
dim_to_commonname(::Pre) = "level"
dim_to_commonname(D::Dim) = string(DimensionalData.name(D))

#########################################################################
# Utilities
#########################################################################
"""
    nckeys(file::String)
Return all keys of the `.nc` file in `file`.
"""
function nckeys(path::String)
    NCDataset(path) do ds
        return keys(ds)
    end
end
nckeys(a::NCDataset) = keys(a)

"""
    ncdetails(file::String, io = stdout)
Print details about the `.nc` file in `file` on `io`.
"""
function ncdetails(file::String, io = stdout)
    NCDataset(file) do ds
        show(io, MIME"text/plain"(), ds)
    end
end
ncdetails(ds::NCDataset, io = stdout) = show(io, MIME"text/plain"(), ds)

"""
    globalattr(file::String) → Dict
Return the global attributes of the .nc file.
"""
function globalattr(file::String)
    NCDataset(file) do ds
        return Dict(ds.attrib)
    end
end

#########################################################################
# Reading
#########################################################################
"""
    ncread(file, var; name, kwargs...) → A
Load the variable `var` from the `file` and convert it into a [`ClimArray`](@ref)
with proper dimension mapping and also containing the variable attributes as a dictionary.
Dimension attributes are also given to the dimensions of `A`, if any exist.
See keywords below for specifications for unstructured grids.

`file` can be a string to a `.nc` file. Or, it can be an
`NCDataset`, which allows you to lazily combine different
`.nc` data (typically split by time), e.g.
```julia
using Glob # for getting all files
alldata = glob("toa_fluxes_2020_*.nc")
file = NCDataset(alldata; aggdim = "time")
A = ClimArray(file, "tow_sw_all")
```

`var` is a `String` denoting which variable to load.
For `.nc` data containing groups `var` can also be a tuple `("group_name", "var_name")`
that loads a specific variable from a specific group.
In this case, the attributes of both the group and the CF-variable are attributed to
the created [`ClimArray`](@ref).

See also [`ncdetails`](@ref), [`nckeys`](@ref) and [`ncwrite`](@ref).

We do two performance improvements while loading the data:
1. If there are no missing values in the data (according to CF standards), the
   returned array is automatically converted to a concrete type (i.e. `Union{Float32, Missing}`
   becomes `Float32`).
2. Dimensions that are ranges (i.e. sampled with constant step size) are automatically
   transformed to a standard Julia `Range` type (which makes sub-selecting faster).

## Keywords
* `name` optionally rename loaded array.
* `grid = nothing` optionally specify whether the underlying grid is `grid = LonLatGrid()`
  or `grid = UnstructuredGrid()`, see [Types of spatial coordinates](@ref).
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
  If `lon, lat` are given, `grid` is automatically assumed `UnstructuredGrid()`.
"""
function ncread(path::Union{String, Vector{String}}, args...; kwargs...)
    NCDataset(path) do ds
        data = ncread(ds, args...; kwargs...)
        return data
    end
end

# TODO: Allow this function to take as input a tuple of indices, e.g. (:, :, 1:5)
# and only load this part, and correctly and instantly make it a ClimArray, which
# can solve "large memory" or "large data" problems. This funcionality
# must be sure to load the correct ranges of dimensions as well though!
#
# TODO: Allow reading multiple variables at once. This has the performance benefit
# of not re-creating dimensions all the time.

function ncread(ds::NCDatasets.AbstractDataset, var;
        name = var2name(var), grid = nothing, lon = nothing, lat = nothing,
    )
    if lon isa Vector && lat isa Vector
        gridtype = UnstructuredGrid()
    elseif isnothing(grid)
        gridtype = autodetect_grid(ds)
    else
        gridtype = grid
    end

    if gridtype == UnstructuredGrid()
        return ncread_unstructured(ds, var, name, lon, lat)
    else
        return ncread_lonlat(ds, var, name)
    end
end

function autodetect_grid(ds)
    if haskey(ds, "reduced_points") || haskey(ds.dim, "ncells") || haskey(ds, "clon")
        return UnstructuredGrid()
    else
        return LonLatGrid()
    end
end

var2name(var) = string(var)
var2name(var::Tuple) = join(var, "_")

#########################################################################
# Reading: LonLatGrid and main reading functionality
#########################################################################
# Notice that this function properly loads even without any spatial coordinate
function ncread_lonlat(ds::NCDatasets.AbstractDataset, var, name)
    cfvar = var isa String ? ds[var] : ds.group[var[1]][var[2]]
    attrib = get_attributes_from_var(ds, cfvar, var)
    A = cfvar |> Array
    dnames = Tuple(NCDatasets.dimnames(cfvar))
    if !any(ismissing, A)
        A = nomissing(A)
    end
    dimensions = create_dims(ds, dnames, A)
    data = ClimArray(A, dimensions; name = Symbol(name), attrib = attrib)
    return data
end

get_attributes_from_var(ds, cfvar, var::String) = Dict(cfvar.attrib)
function get_attributes_from_var(ds, cfvar, var::Tuple)
    d1 = Dict(cfvar.attrib)
    d2 = Dict(ds.group[var[1]].attrib)
    return merge!(d2, d1)
end


"""
    create_dims(ds::NCDatasets.AbstractDataset, dnames)
Create a tuple of `Dimension`s from the `dnames` (tuple of strings).
"""
function create_dims(ds::NCDatasets.AbstractDataset, dnames, A)
    true_dims = to_proper_dimensions(dnames)
    dim_values = extract_dim_values(ds, dnames, A)
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
            Making a generic dimension named `Dim{:$n}` for now...
            """
            push!(r, Dim{Symbol(n)})
        end
    end
    return (r...,)
end
export Dim # for generic dimensions this must be exported

function extract_dim_values(ds::NCDatasets.AbstractDataset, dnames, A)
    dim_values = AbstractVector[]
    for (i, n) in enumerate(dnames)
        if haskey(ds, n)
            push!(dim_values, Vector(ds[n]))
        else
            @warn """
            Dimension named "$n" does not have values in the dataset.
            Using the range `1:size(A, i)` as the dimension values instead,
            where `i` is the dimension index.
            """
            push!(dim_values, 1:size(A, i))
        end
    end
    return dim_values
end


#########################################################################
# Making vectors → ranges
#########################################################################
function vector2range(x::Vector{<:Real})
    length(x) < 3 && return x
    dx = x[2]-x[1]
    for i in 3:length(x)
        x[i]-x[i-1] ≠ dx && return x # if no constant step, return array as is
    end
    # do not check value equality, only difference equality
    return x[1]:dx:x[end]
end

function vector2range(t::Vector{<:Dates.AbstractTime})
    tsamp = temporal_sampling(t)
    period = tsamp2period(tsamp)
    isnothing(period) && return t
    t1 = period == :hourly ? t[1] : Date(t[1])
    tf = period == :hourly ? t[end] : Date(t[end])
    r = t1:period:tf
    return r == t ? r : t # final safety check to ensure equal values
end

vector2range(r::AbstractRange) = r


#########################################################################
# Reading: UnstructuredGrid
#########################################################################
export SVector

function ncread_unstructured(ds::NCDatasets.AbstractDataset, var::String, name, lon, lat)
    cfvar = var isa String ? ds[var] : ds.group[var[1]][var[2]]
    attrib = get_attributes_from_var(ds, cfvar, var)

    # Here we generate the longitudes and latitudes based on whether we have
    # reduced points or not, and obtain the name of the coord dimension in the `ds`
    if isnothing(lon)
        lonlat, original_grid_dim = load_coordinate_points(ds)
        lonlat = convert_to_degrees(lonlat, ds)
    else
        lonlat = [SVector(lo, la) for (lo, la) in zip(lon, lat)]
        original_grid_dim = intersect(NCDatasets.dimnames(cfvar), POSSIBLE_CELL_NAMES)[1]
        lonlat = convert_to_degrees(lonlat)
    end

    # TODO: I've noticed that this converts integer dimension (like pressure)
    # into Float64, but I'm not sure why...
    alldims = [NCDatasets.dimnames(cfvar)...]

    # Set up the remaining dimensions of the dataset
    @assert original_grid_dim ∈ alldims
    i = findfirst(x -> x == original_grid_dim, alldims)
    remainingdims = deleteat!(copy(alldims), i)
    A = Array(cfvar)
    actualdims = Any[create_dims(ds, remainingdims, A)...]

    # Make coordinate dimension
    si = sortperm(lonlat, by = reverse)
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

const POSSIBLE_CELL_NAMES = ["ncells", "cell", "rgrid", "grid"]

function has_unstructured_key(ds)
    any(x -> haskey(ds.dim, x), POSSIBLE_CELL_NAMES) ||
    haskey(ds, "lon") || haskey(ds, "clon")
end

function load_coordinate_points(ds)
    if haskey(ds, "reduced_points")
        lonlat = reduced_grid_to_points(ds["lat"], ds["reduced_points"])
        original_grid_dim = "rgrid" # TODO: Specific to CDO Gaussian grid
    elseif has_unstructured_key(ds)
        if haskey(ds, "lon")
            lons = ds["lon"] |> Array .|> wrap_lon
            lats = ds["lat"] |> Array
            original_grid_dim = NCDatasets.dimnames(ds["lon"])[1]
        elseif haskey(ds, "clon")
            lons = ds["clon"] |> Array .|> wrap_lon
            lats = ds["clat"] |> Array
            original_grid_dim = NCDatasets.dimnames(ds["clon"])[1]
        else
            error("""
            We didn't find key `"lon"` or `"clon"` that represents the longitude of each
            polygon in an unstructured grid.
            """)
        end
        lonlat = [SVector(lo, la) for (lo, la) in zip(lons, lats)]
    else
        error("""
        We didn't find any of the following keys: `"ncells", "cells", "reduced_points",
        "clon", "lon"`, at least one of which is mandatory for unstructured grid.
        If your data stores the "cell" information with a different name, please open
        an issue and let us know!
        """)
    end
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

#########################################################################
# Saving to .nc files
#########################################################################
const DEFAULT_ATTRIBS = Dict(
    "time" => Dict(
        "units" => "days since 0000-00-01 00:00:00",
        "standard_name" => "time"
    ),
    "lon" => Dict(
        "units" => "degrees_east",
        "standard_name" => "longitude",
        "valid_range" => Float32[-180.0, 360.0]
    ),
    "lat" => Dict(
        "units" => "degrees_north",
        "standard_name" => "latitude",
        "valid_range" => Float32[-90.0, 90.0]
    ),
    "level" => Dict(
        "units" => "millibars",
        "long_name" => "pressure_level",
    ),
)

"""
    ncwrite(file::String, Xs; globalattr = Dict())
Write the given `ClimArray` instances (any iterable of `ClimArray`s or a single `ClimArray`)
to a `.nc` file following CF standard conventions using NCDatasets.jl.
Optionally specify global attributes for the `.nc` file.

The metadata of the arrays in `Xs`, as well as their dimensions, are properly written
in the `.nc` file and any necessary type convertions happen automatically.

**WARNING**: We assume that any dimensions shared between the `Xs` are identical.

See also [`ncread`](@ref).
"""
function ncwrite(file::String, X::ClimArray; globalattr = Dict())
    ncwrite(file, (X,); globalattr)
end
function ncwrite(file::String, Xs; globalattr = Dict())

    # TODO: Fixing this is very easy. Simply make a `"ncells"` dimension, and then write
    # the `"lon"` and `"lat"` cfvariables to the nc file by decomposing the coordinates
    # into longitude and latitude.
    if any(X -> hasdim(X, Coord), Xs)
        error("""
        Outputing `UnstructuredGrid` coordinates to .nc files is not yet supported,
        but it is an easy fix, see source of `ncwrite`.
        """)
    end

    # ds = NCDataset(file, "c"; attrib = globalattr)
    NCDataset(file, "c"; attrib = globalattr) do ds
        for (i, X) in enumerate(Xs)
            n = string(X.name)
            if n == ""
                n = "x$i"
                @warn "$i-th ClimArray has no name, naming it $(n) instead."
            end
            println("processing variable $(n)...")
            add_dims_to_ncfile!(ds, dims(X))
            attrib = X.attrib
            if (isnothing(attrib) || attrib == DimensionalData.NoMetadata())
                attrib = Dict()
            end
            dnames = dim_to_commonname.(dims(X))
            data = Array(X)
            NCDatasets.defVar(ds, n, data, (dnames...,); attrib)
        end
        # close(ds)
    end
end

function add_dims_to_ncfile!(ds::NCDatasets.AbstractDataset, dimensions::Tuple)
    dnames = dim_to_commonname.(dimensions)
    for (i, d) ∈ enumerate(dnames)
        haskey(ds, d) && continue
        println("writing dimension $d...")
        v = dimensions[i].val
        # this conversion to DateTime is necessary because CFTime.jl doesn't support Date
        eltype(v) == Date && (v = DateTime.(v))
        l = length(v)
        NCDatasets.defDim(ds, d, l) # add dimension entry
        attrib = dimensions[i].metadata
        if (isnothing(attrib) || attrib == DimensionalData.NoMetadata()) && haskey(DEFAULT_ATTRIBS, d)
            @warn "Dimension $d has no attributes, adding default attributes (mandatory)."
            attrib = DEFAULT_ATTRIBS[d]
        end
        # write dimension values as a variable as well (mandatory)
        NCDatasets.defVar(ds, d, v, (d, ); attrib = attrib)
    end
end
