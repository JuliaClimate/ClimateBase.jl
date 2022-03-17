#########################################################################
# Defaults
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

#########################################################################
# Write: standard data
#########################################################################
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
        Outputing `CoordinateSpace` coordinates to .nc files is not yet supported,
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
        attrib = DimensionalData.metadata(dimensions[i])
        if (isnothing(attrib) || attrib == DimensionalData.NoMetadata()) && haskey(DEFAULT_ATTRIBS, d)
            @warn "Dimension $d has no attributes, adding default attributes (mandatory)."
            attrib = DEFAULT_ATTRIBS[d]
        end
        # write dimension values as a variable as well (mandatory)
        NCDatasets.defVar(ds, d, v, (d, ); attrib = attrib)
    end
end
