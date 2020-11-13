#=
Code related with input output (IO) of .nc files directly to/from ClimArrays
An initial version of parts of this code was taken from:
https://github.com/rafaqz/GeoData.jl
=#
using NCDatasets
export NCDataset
export nckeys, ncdetails
export climarrays_to_nc

dim_to_commonname(::Lat) = "lat"
dim_to_commonname(::Lon) = "lon"
dim_to_commonname(::Time) = "time"
dim_to_commonname(::Pre) = "level"
dim_to_commonname(D::Dim) = string(DimensionalData.name(D))

#########################################################################
# NCDatasets → DimensionalArray convertions and loading
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
    ClimArray(file::NCDataset, var::String, name = var) -> A
Load the variable `var` from the `file` and convert it
into a `ClimArray` which also contains the variable attributes as a dictionary.

Notice that `file` should be an `NCDataset`, which allows you to lazily combine different
`.nc` data (typically split by time), e.g.
```julia
alldata = ["toa_fluxes_2020_\$(i).nc" for i in 1:12]
file = NCDataset(alldata; aggdim = "time")
A = ClimArray(file, "tow_sw_all")
```
(of course you can just do `NCDataset("file.nc")` for single files).

We do two performance improvements while loading the data:
1. If there are no missing values in the data (according to CF standards), the
   returned array is automatically converted to a concrete type (i.e. `Union{Float32, Missing}`
   becomes `Float32`).
2. Dimensions that are ranges (i.e. sampled with constant step size) are automatically
   transformed to a standard Julia `Range` type (which makes sub-selecting faster).


At the moment, support for auto-loading equal area space types does not exist,
see [Types of spatial coordinates](@ref). But
you can easily transform them yourself into a `ClimArray` by doing e.g.:
```julia
file = NCDataset("some_file_with_eqarea.nc")
lons = file["lon"]
lats = file["lat"]
coords = [SVector(lo, la) for (lo, la) in zip(lons, lats)]
t = file["time"]
dimensions = (Coord(coords), Time(t))
data = file["actual_data_like_radiation"]
A = ClimArray(data, dimensions)
```
"""
function ClimArray(path::Union{String, Vector{String}}, args...; kwargs...)
    NCDataset(path) do ds
        data = ClimArray(ds, args...; kwargs...)
        return data
    end
end

# TODO: Allow this function to take as input a tuple of indices, e.g. (:, :, 1:5)
# and only load this part, and correctly and instantly make it a ClimArray, which
# can solve "large memory" or "large data" problems. This funcionality
# must be sure to load the correct ranges of dimensions as well though!

function ClimArray(ds::NCDatasets.AbstractDataset, var::String, name = var; eqarea = false)
    svar = string(var)
    cfvar = ds[svar]
    attrib = Dict(cfvar.attrib)
    A = cfvar |> Array
    if eqarea
        # TODO: This piece of code is specific to CDO output...
        if haskey(ds, "ncells") # this is the equal area grid, so we make a Coord dimension
            lon = ds["lon"] |> Array .|> wrap_lon
            lat = ds["lat"] |> Array
            time = ds["time"] |> Array
            lonlat = [SVector(lon[i], lat[i]) for i in 1:length(lon)]
            # here we sort lonlat and A in ascending latitude order,
            # because the CDO output has reverse or even totally unsorted order
            si = sortperm(lonlat, by = reverse)
            data = ClimArray(A[si, :], (Coord(lonlat[si]), Time(time));
            attrib = attrib, name = svar)
        elseif haskey(ds, "reduced_points")
            # TODO: This can be easily upgraded to arbitary dimensions via a simple
            # dimension replacement / permutation at the end
            lonlat = reduced_grid_to_points(ds["lat"], ds["reduced_points"])
            si = sortperm(lonlat, by = reverse)
            time = ds["time"] |> Array
            data = ClimArray(A[si, :], (Coord(lonlat[si]), Time(time));
            name = svar, attrib = attrib)
        else
            error("Don't know how to handle this equal area grid!")
        end
    else # standard variables
        dnames = Tuple(NCDatasets.dimnames(cfvar))
        data = ClimArray(A, create_dims(ds, dnames); name = Symbol(name), attrib = attrib)
    end
    if !any(ismissing, data)
        data = nomissing(data)
    end
    return data
end

"""
    create_dims(ds::NCDatasets.AbstractDataset, dnames)
Create a tuple of `Dimension`s from the `dnames` (tuple of strings).
"""
function create_dims(ds::NCDatasets.AbstractDataset, dnames)
    # true_dims = getindex.(Ref(COMMONNAMES), dnames)
    true_dims = to_proper_dimensions(dnames)
    dim_values = Array.(getindex.(Ref(ds), dnames))
    optimal_values = vector2range.(dim_values)
    attribs = [
        ds[d].attrib isa NCDatasets.BaseAttributes ? Dict(ds[d].attrib) : nothing
        for d in dnames
    ]
    out = []
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
            Dimension name "$n" not in common names. Strongly recommended to ask for
            adding this name to COMMONNAMES on github. Making generic dimension for now...
            """
            push!(r, Dim{Symbol(n)})
        end
    end
    return (r...,)
end

export Dim # for generic dimensions this must be exported

#########################################################################
# Making vectors → ranges
#########################################################################
function vector2range(x::Vector{<:Real})
    dx = x[2]-x[1]
    for i in 3:length(x)
        x[i]-x[i-1] ≠ dx && return x # if no constant step, return array as is
    end
    r = x[1]:dx:x[end]
    @assert r == x
    return r
end

function vector2range(t::Vector{<:DateTime})
    !sampled_less_than_date(t) && return vector2range(Date.(t))
    # TODO: implement hourly sampling here
    @warn "Hourly sampling not yet implemented."
    return t
end

function vector2range(t::Vector{<:Date})
    tsamp = temporal_sampling(t)
    period = tsamp2period(tsamp)
    r = t[1]:period:t[end]
    @assert r == t
    return r
end


#########################################################################
# Equal area related
#########################################################################
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
    climarrays_to_nc(file::String, Xs; globalattr = Dict())
Write the given `ClimArray` instances (any iterable of `ClimArray`s or a single `ClimArray`)
to a `.nc` file following CF standard conventions using NCDatasets.jl.
Optionally specify global attributes for the `.nc` file.

The metadata of the arrays in `Xs`, as well as their dimensions, are properly written
in the `.nc` file and any necessary type convertions happen automatically.

**WARNING**: We assume that any dimensions shared between the `Xs` are identical.
"""
function climarrays_to_nc(file::String, X::ClimArray; globalattr = Dict())
    climarrays_to_nc(file, (X,); globalattr)
end
function climarrays_to_nc(file::String, Xs; globalattr = Dict())
    ds = NCDataset(file, "c"; attrib = globalattr)
    for (i, fieldname) in enumerate(Xs)
        println("processing variable $fieldname...")
        W = climarray_from_xarray(xa, fieldname)
        println("converting to CERES format...")
        X = to_CERES_latitude(W)
        println("writing dimensions...")
        add_dims_to_ncfile!(ds, dims(X))
        println("writing the CF-variable...")
        attrib = X.attrib
        dnames = dim_to_commonname.(dims(X))
        data = Array(X)
        @show (fieldname, summary(data), dnames)
        defVar(ds, fieldname, data, (dnames...,); attrib)
    end
    close(ds)
end

function add_dims_to_ncfile!(ds::NCDatasets.AbstractDataset, dimensions::Tuple)
    dnames = dim_to_commonname.(dimensions)
    for (i, d) ∈ enumerate(dnames)
        haskey(ds, d) && continue
        v = dimensions[i].val
        # this conversion to DateTime is necessary because CFTime.jl doesn't support Date
        eltype(v) == Date && (v = DateTime.(v))
        l = length(v)
        defDim(ds, d, l) # add dimension entry
        attrib = dimensions[i].metadata
        if isnothing(attrib) && haskey(DEFAULT_ATTRIBS, d)
            @warn "Dimension $d has no attributes, adding default attributes (mandatory)."
            attrib = DEFAULT_ATTRIBS[d]
        end
        # write dimension values as a variable as well (mandatory)
        defVar(ds, d, v, (d, ); attrib = attrib)
    end
end
