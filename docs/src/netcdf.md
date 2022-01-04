# NetCDF IO

ClimateBase.jl has support for `"file.nc" ⇆ ClimArray`.
Usually this is done using NCDatasets.jl, but see below for a function that translates a loaded `xarray` (from Python) into `ClimArray`.

## Read

To load a `ClimArray` directly from an `.nc` file do:
```@docs
ncread
```

Notice that (at the moment) we use a pre-defined mapping of common names to proper dimensions - please feel free to extend the following via a Pull Request:
```@example main
using ClimateBase # hide
ClimateBase.COMMONNAMES
```

Also, the following convenience functions are provided for examining the content of on-disk `.nc` files without loading all data on memory.
```@docs
nckeys
ncdetails
ncsize
globalattr
```

## Write
You can also write a bunch of `ClimArray`s directly into an `.nc` file with
```@docs
ncwrite
```

## xarray
You can use the following functions (which are not defined and exported in `ClimateBase` to avoid dependency on PyCall.jl) to load data using Python's `xarray`.
```julia
using ClimateBase, Dates
# This needs to numpy, xarray and dask installed from Conda
using PyCall
xr = pyimport("xarray")
np = pyimport("numpy")

function climarray_from_xarray(xa, fieldname, name = fieldname)
    w = getproperty(xa, Symbol(fieldname))
    raw_data = Array(w.values)
    dnames = collect(w.dims) # dimensions in string name
    dim_values, dim_attrs = extract_dimension_values_xarray(xa, dnames)
    @assert collect(size(raw_data)) == length.(dim_values)
    actual_dims = create_dims_xarray(dnames, dim_values, dim_attrs)
    ca = ClimArray(raw_data, actual_dims, name; attrib = w.attrs)
end

function extract_dimension_values_xarray(xa, dnames = collect(xa.dims))
    dim_values = []
    dim_attrs = Vector{Any}(fill(nothing, length(dnames)))
    for (i, d) in enumerate(dnames)
        dim_attrs[i] = getproperty(xa, d).attrs
        x = getproperty(xa, d).values
        if d ≠ "time"
            push!(dim_values, x)
        else
            # Dates need special handling to be transformed into `DateTime`.
            dates = [np.datetime_as_string(y)[1:19] for y in x]
            dates = DateTime.(dates)
            push!(dim_values, dates)
        end
    end
    return dim_values, dim_attrs
end

function create_dims_xarray(dnames, dim_values, dim_attrs)
    true_dims = ClimateBase.to_proper_dimensions(dnames)
    optimal_values = ClimateBase.vector2range.(dim_values)
    out = []
    for i in 1:length(true_dims)
        push!(out, true_dims[i](optimal_values[i]; metadata = dim_attrs[i]))
    end
    return (out...,)
end

# Load some data
xa = xr.open_mfdataset(ERA5_files_path)
X = climarray_from_xarray(xa, "w", "optional name")
```

## Ensemble types
A dedicated type representing ensembles has no reason to exist in ClimateBase.jl.
As the package takes advantage of standard Julia datastructures and syntax, those can be used to represent "ensembles". For example to do an "ensemble global mean" you can just do:
```julia
# load all data
E = [ClimArray("ensemble_$i.nc", "x") for i in 1:10]
# mean from all data
global_mean = mean(spacemean(X) for X in E)
```
where you see that the "ensemble" was represented just as a `Vector{ClimArray}`.
Of course, this requires that all data can fit into memory, but this is so far the only way ClimateBase.jl operates anyways.
