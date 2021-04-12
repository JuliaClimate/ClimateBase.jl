#=
Code related with input output (IO) of .nc files directly to/from ClimArrays
An initial version of parts of this code was taken from:
https://github.com/rafaqz/GeoData.jl
=#
using NCDatasets
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
    ncread(file::Union{String,NCDataset}, var::String; name, grid) → A
Load the variable `var` from the `file` and convert it into a `ClimArray`
with proper dimension mapping and also containing the variable attributes as a dictionary.
Dimension attributes are also given to the dimensions of `A`, if any exist.

Notice that `file` can be an `NCDataset`, which allows you to lazily combine different
`.nc` data (typically split by time), e.g.
```julia
alldata = ["toa_fluxes_2020_\$(i).nc" for i in 1:12]
file = NCDataset(alldata; aggdim = "time")
A = ClimArray(file, "tow_sw_all")
```
(but you can also directly give the string to a single file `"file.nc"`
if data are contained in a single file).

We do two performance improvements while loading the data:
1. If there are no missing values in the data (according to CF standards), the
   returned array is automatically converted to a concrete type (i.e. `Union{Float32, Missing}`
   becomes `Float32`).
2. Dimensions that are ranges (i.e. sampled with constant step size) are automatically
   transformed to a standard Julia `Range` type (which makes sub-selecting faster).

## Keywords
* `name = var` optionally rename loaded array.
* `grid = nothing` optionally specify whether the underlying grid is [`LonLatGrid`](@ref)
  or [`UnstructuredGrid`](@ref).
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

function ncread(ds::NCDatasets.AbstractDataset, var::String, name = var)
    svar = string(var)
    cfvar = ds[svar]
    attrib = Dict(cfvar.attrib)
    A = cfvar |> Array
    dnames = Tuple(NCDatasets.dimnames(cfvar))
    if !any(ismissing, A)
        A = nomissing(A)
    end
    data = ClimArray(A, create_dims(ds, dnames); name = Symbol(name), attrib = attrib)
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
    length(x) < 3 && return x
    dx = x[2]-x[1]
    for i in 3:length(x)
        x[i]-x[i-1] ≠ dx && return x # if no constant step, return array as is
    end
    r = x[1]:dx:x[end]
    return r == x ? r : x
end

function vector2range(t::Vector{<:Dates.AbstractTime})
    tsamp = temporal_sampling(t)
    period = tsamp2period(tsamp)
    isnothing(period) && return t
    r = t[1]:period:t[end]
    return r == t ? r : t
end

vector2range(r::AbstractRange) = r

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
"""
function ncwrite(file::String, X::ClimArray; globalattr = Dict())
    ncwrite(file, (X,); globalattr)
end
function ncwrite(file::String, Xs; globalattr = Dict())
    ds = NCDataset(file, "c"; attrib = globalattr)
    # NCDataset("file.nc", "c"; attrib = globalattr) do ds
        for (i, X) in enumerate(Xs)
            n = string(X.name)
            if n == ""
                n = "x$i"
                @warn "$i-th ClimArray has no name, naming it $(n) instead."
            end
            println("processing variable $(n)...")
            add_dims_to_ncfile!(ds, dims(X))
            println("writing the CF-variable...")
            attrib = X.attrib
            isnothing(attrib) && (attrib = Dict())
            dnames = dim_to_commonname.(dims(X))
            data = Array(X)
            defVar(ds, n, data, (dnames...,); attrib)
        end
        close(ds)
    # end
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
        defDim(ds, d, l) # add dimension entry
        attrib = dimensions[i].metadata
        if (isnothing(attrib) || attrib == DimensionalData.NoMetadata()) && haskey(DEFAULT_ATTRIBS, d)
            @warn "Dimension $d has no attributes, adding default attributes (mandatory)."
            attrib = DEFAULT_ATTRIBS[d]
        end
        # write dimension values as a variable as well (mandatory)
        defVar(ds, d, v, (d, ); attrib = attrib)
    end
end
