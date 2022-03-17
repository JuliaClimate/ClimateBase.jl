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
    end
end

function add_dims_to_ncfile!(ds::NCDatasets.AbstractDataset, dimensions::Tuple)
    dnames = dim_to_commonname.(dimensions)
    dims_in_ds = [x[1] for x in ds.dim]
    for (d, dname) ∈ zip(dimensions, dnames)
        dname ∈ dims_in_ds && continue
        println("writing dimension $dname...")
        v = gnv(d); l = length(v)
        NCDatasets.defDim(ds, dname, l) # add dimension entry
        if d isa Coord
            # Define clon/clat variables with this dimension
            lons = getindex.(v, 1); lats = getindex.(v, 2)
            NCDatasets.defVar(ds, "clon", lons, (dname, ); attrib = DEFAULT_ATTRIBS["lon"])
            NCDatasets.defVar(ds, "clat", lats, (dname, ); attrib = DEFAULT_ATTRIBS["lat"])
        else
            # this conversion to DateTime is necessary because CFTime.jl doesn't support Date
            eltype(v) == Date && (v = DateTime.(v))
            attrib = DimensionalData.metadata(d)
            if (isnothing(attrib) || attrib == DimensionalData.NoMetadata()) &&
                    haskey(DEFAULT_ATTRIBS, dname)
                @warn "Dimension $dname has no attributes, adding default attributes."
                attrib = DEFAULT_ATTRIBS[dname]
            end
            # write dimension values as a variable as well (mandatory)
            NCDatasets.defVar(ds, dname, v, (dname, ); attrib = attrib)
        end
    end
end

#########################################################################
# Write: CoordinateSpace
#########################################################################
