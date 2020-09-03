#=
Code related with loading .nc file data directly into a dimensional array
An initial version of parts of this code was taken from:
https://github.com/rafaqz/GeoData.jl
=#
using NCDatasets
export NCDataset
export nckeys, ncdetails
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

"""
    ClimArray(file::NCDataset, var::String) -> A
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

If there are no missing values in the data (according to CF standards), the
returned array is automatically converted to a concrete type (i.e. `Union{Float32, Missing}`
becomes `Float32`).

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
# can solve "large memory" or "large data" problems.
function ClimArray(ds::NCDataset, var::String; eqarea = false)
    svar = string(var)
    cfvar = ds[svar]
    attrib = Dict(cfvar.attrib)
    A = cfvar |> Array
    if eqarea
        # TODO: I have to re-work this code to be more general and allow other dimensions
        # as well!!!!
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
        data = ClimArray(A, create_dims(ds, dnames); name = svar, attrib = attrib)
    end
    if !any(ismissing, data)
        data = nomissing(data)
    end
    return data
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
