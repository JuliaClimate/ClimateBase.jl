#=
Code related with loading .nc file data directly into a dimensional array
An initial version of parts of this code was taken from:
https://github.com/rafaqz/GeoData.jl
=#
using NCDatasets

#########################################################################
# NCDatasets → DimensionalArray convertions and loading
#########################################################################
"""
    create_dims(ds::NCDataset, dnames)
Create a tuple of `Dimension`s from the `dnames` (tuple of strings).
"""
function create_dims(ds::NCDataset, dnames)
    true_dims = getindex.(Ref(COMMONNAMES), dnames)
    dim_values = Array.(getindex.(Ref(ds), dnames))
    return dim_values .|> true_dims
end

function globalattributes(file::String)
    # TODO
end

function allkeys(path::String)
    NCDataset(path) do ds
        return keys(ds)
    end
end
allkeys(a::NCDataset) = keys(a)

# TODO: Make a wrapper type, or just a new type that subtypes AbDimArray
# and has only one extra field called attributes.

"""
    get_var_as_dimarray(file::String, var::String; eqarea = false)) -> A, attrib
Return the variable `var` from the `.nc` file located at path, by converting it
into a `DimensionalArray`. In addition return its attributes as a dictionary.

Notice that if there are no missing values in the data (according to CF standards), the
returned array is automatically converted to a concrete type (i.e. `Union{Float32, Missing}`
becomes `Float32`).

Keyword `eqarea` denotes if the underlying grid is structured or unstructured (equal area).
"""
function get_var_as_dimarray(path, var;
        eqarea = false, hasmissing = false)
    NCDataset(path) do ds
        svar = string(var)
        cfvar = ds[svar]
        A = cfvar |> Array
        if eqarea
            if haskey(ds, "ncells") # this is the equal area grid, so we make a Coord dimension
                lon = ds["lon"] |> Array .|> wrap_lon
                lat = ds["lat"] |> Array
                time = ds["time"] |> Array
                lonlat = [SVector(lon[i], lat[i]) for i in 1:length(lon)]
                # here we sort lonlat and A in ascending latitude order,
                # because the CDO output has reverse or even totally unsorted order
                si = sortperm(lonlat, by = x -> x[2])
                data = DimensionalArray(A[si, :], (Coord(lonlat[si]), Time(time)), svar)
            elseif haskey(ds, "reduced_points")
                lonlat = reduced_grid_to_points(ds["lat"], ds["reduced_points"])
                si = sortperm(lonlat, by = x -> x[2])
                time = ds["time"] |> Array
                data = DimensionalArray(A[si, :], (Coord(lonlat[si]), Time(time)), svar)
            end
        else # standard variables
            dnames = Tuple(NCDatasets.dimnames(cfvar))
            data = DimensionalArray(A, create_dims(ds, dnames), svar)
        end
        if !any(ismissing, data)
            data = nomissing(data)
        end
        return data, Dict(cfvar.attrib)
    end
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
