# Introduction
`ClimateBase` is a Julia package offering basic functionality for analyzing data that are typically in the form used by climate sciences.
These data are dimensional & spatiotemporal but the corresponding dimensions all need special handling.
For example the most common dimensions are longitude, latitude and time.

* longitude is by definition a periodic dimension
* latitude is a linear dimension. However because the coordinate system most often used in climate sciences is a grid of longitude × latitude (in equal degrees) the area element of space depends on latitude and this needs to be taken into account.
* time is a linear dimension *in principle*, but its values are `<: AbstractDateTime` instead of `<: Real`. The human calendar (where these values come from) is periodic but each period may not correspond to the same physical time, and this also needs to be taken into account.

`ClimateBase` is structured to deal with these intricacies, and in addition offer several functionalities commonly used, and sought after, by climate scientists.
It also serves as the base building block for `ClimateTools`, which offers more advanced functionalities.

At the moment the focus of `ClimateBase` is **not** on operating on data *on disk*. It is designed for in-memory climate data exploration and manipulation.

### Installation
This package is registered and you can install it with
```julia
using Pkg; Pkg.add("ClimateBase")
```
Make sure your installed version coincides with the one in this docs (see bottom left corner of this page).

## `ClimArray`: the core data structure
This project treats "climate data" as an instance of [`ClimArray`](@ref).
At the moment `ClimArray` is a subtype of `DimensionalArray` from DimensionalData.jl.
A brief introduction to DimensionalData.jl is copied here from its docs, because basic handling of a `ClimArray` comes from DimensionalData.jl.
DimensionalData.jl allows to dimensionally-index data by their values.

E.g. you can create an array with
```@example main
using ClimateBase, Dates
Time = ClimateBase.Ti # `Time` is more intuitive than `Ti`
lats = -90:5:90
lons = 0:10:359
t = Date(2000, 3, 15):Month(1):Date(2020, 3, 15)
# Here we wrap all dimension data into proper dimensions:
dimensions = (Lon(lons), Lat(lats), Time(t))
# where `Lon, Lat, Time` are `Dimension`s exported by ClimateBase
# combining the array data with dimensions makes a `ClimArray`:
data = rand(36, 37, 241) # same size as `dimensions`
A = ClimArray(data, dimensions)
```

You can then select a specific time slice at `Date(2011,5,15)` and a longitude interval between 0 and 30 degrees like so:
```@example main
B = A[Lon(Between(0, 30)), Time(At(Date(2011,5,15)))]
```

With `ClimArray` you can use convenience, physically-inspired functions that do automatic (and correct) weighting.
For example the latitudinal mean of `B` is simply
```@example main
C = latmean(B)
```
where in this averaging process each data point is weighted by the cosine of its latitude.

### Making a `ClimArray`
You can create a `ClimArray` yourself, or you can load data from an `.nc` file with CF-conventions, see [NetCDF IO](@ref).
```@docs
ClimArray(::AbstractArray, ::Tuple)
```
It is strongly recommended to use the dimensions we export (because we dispatch on them and use their information):
```@example
using ClimateBase, DimensionalData # hide
for D in ClimateBase.STANDARD_DIMS
    println(D, " (full name = $(DimensionalData.name(D)))")
end
```
We explicitly assume that `Lon, Lat` are measured in degrees and not radians or meters (extremely important for spatial averaging processes).

## NetCDF IO

ClimateBase.jl has support for `file.nc ⇆ ClimArray`.
Usually this is done using NCDatasets.jl, but see below for a function that translates a loaded `xarray` (from Python) into `ClimArray`.

### Read

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
globalattr
```

### Write
You can also write a bunch of `ClimArray`s directly into an `.nc` file with
```@docs
ncwrite
```

### xarray
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


## Crash-course to DimensionalData.jl
```@docs
DimensionalData
```
