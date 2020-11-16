# Introduction
`ClimateBase` is a Julia package offering basic functionality for analyzing data that are typically in the form used by climate sciences.
These data are dimensional & spatiotemporal but the corresponding dimensions all need special handling.
For example the most common dimensions are longitude, latitude and time.

* longitude is by definition a periodic dimension
* latitude is a linear dimension. However because the coordinate system most often used in climate sciences is a grid of longitude × latitude (in equal degrees) the area element of space depends on latitude and this needs to be taken into account.
* time is a linear dimension *in principle*, but its values are `<: AbstractDateTime` instead of `<: Real`. The human calendar (where these values come from) is periodic but each period may not correspond to the same physical time, and this also needs to be taken into account.

`ClimateBase` is structured to deal with these intricacies, and in addition offer several functionalities commonly used, and sought after, by climate scientists.
It also serves as the base building block for `ClimateTools`, which offers more advanced functionalities.

The focus of `ClimateBase` is **not** loading data, nor operating on data *on disk*. It is designed for in-memory climate data exploration and manipulation.
That being said, basic data loading functionality is offered in terms of `NCDatasets`, see below.

## `ClimArray`: the core data structure
This project treats "climate data" as a [`ClimArray`](@ref), which uses the DimensionalData.jl interface.
`ClimArray` is *almost* equivalent to `DimensionalArray`.
A (brief) introduction to DimensionalData.jl is copied here from its docs, because basic knowledge of how to handle a `ClimArray` is assumed in our docs.

DimensionalData.jl allows truly convenient handling of climate data, where it is important to be able to dimensionally-index data by their values.

E.g. you can create an array with
```@example main
using ClimateBase, Dates
Time = ClimateBase.Ti # more intuitive
lats = -90:5:90
lons = 0:10:359
t = Date(2000, 3, 15):Month(1):Date(2020, 3, 15)
dimensions = (Lon(lons), Lat(lats), Time(t))
A = ClimArray(rand(36, 37, 241), dimensions)
```
and then select a specific timeslice at `Date(2011,5,15)` and a longitude interval between 0 and 30 degrees like so:
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
ClimArray(::Union{String, Vector{String}})
```

Notice that (at the moment) we use a pre-defined mapping of common names to proper dimensions - please feel free to extend the following via a Pull Request:
```@example main
using ClimateBase # hide
ClimateBase.COMMONNAMES
```

Also, two convenience functions are provided for examining the content of on-disk `.nc` files without loading all data on memory.
```@docs
nckeys
ncdetails
```

### Write
You can also write a bunch of `ClimArray`s directly into an `.nc` file with
```@docs
climarrays_to_nc
```

### xarray
You can use the following functions (which are not defined and exported in `ClimateBase` to avoid dependency on PyCall.jl)
```julia
using ClimateBase, Dates
# This needs to numpy, xarray and dask installed from Conda
using PyCall
xr = pyimport("xarray")
np = pyimport("numpy")

"""
    climarray_from_xarray(xa, fieldname, name = Symbol(fieldname))
Load underlying field with given `fieldname` from the given xarray instance `xa`,
optionally providing a name for it. This `xa` can be loaded with commands like
```julia
using PyCall, ClimateBase
xr = pyimport("xarray")
ERA5_files = "some_file_name.nc"
xa = xr.open_mfdataset(ERA5_files)
X = climarray_from_xarray(xa, "w", "optional name")
```
"""
function climarray_from_xarray(xa, fieldname, name = Symbol(fieldname))
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
            # This date specification assumes up to day sampling (hence the 1:10)
            dates = [np.datetime_as_string(y)[1:10] for y in x]
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
```

## Temporal
Functions related with the `Time` dimension.
```@docs
timemean
timeagg
monthlyagg
yearlyagg
temporalrange
maxyearspan
temporal_sampling
time_in_days
```

## Spatial

### Types of spatial coordinates
Most of the time the spatial information of your data is in the form of a Longitude × Latitude grid. This is simply achieved via the existence of two dimensions (`Lon, Lat`) in your dimensional data array. Height, although representing physical space as well, is not considered part of the "spatial dimensions", and is treated as any other additional dimension.
This type of space is called `Grid`. It is assumed throughout that longitude and latitude are measured in **degrees**.

Another type of spatial coordinates is supported, and that is of **equal-area**, called `EqArea`.
There, the spatial dimension is instead given by a single `Vector` of coordinate locations, i.e. 2-element `SVector(longitude, latitude)`. The dimension of this vector is `Coord`.
Each point in this vector corresponds to a polygon (typically triangle or trapezoid) that covers equal amount of spatial area as any other point.
The actual limits of each polygon are not included in the dimension.
Typical examples of such equal area grids are reduced (or thinned) Gaussian grids or icosahedral-based grids.

Within ClimateBase.jl aims to work with either type of spatial coordinate system. Therefore, physically inspired averaging functions, like [`spacemean`](@ref) or [`zonalmean`](@ref), work for both types of spatial coordinates.
In addition, the function `spatialidxs` returns an iterator over the spatial coordinates of the data, and works for both types (grid or equal-area):
```@docs
spatialidxs
```

### Spatial aggregation
```@docs
zonalmean
latmean
spacemean
spaceagg
hemispheric_means
hemispheric_functions
lonlatfirst
```

## General aggregation
The physical averages of the previous section are done by taking advantage of a general aggregation syntax, which works with any aggregating function like `mean, sum, std`, etc.
```@docs
dropagg
collapse
```

## Timeseries Analysis
```@docs
sinusoidal_continuation
seasonal_decomposition
```

## Climate quantities
Functions that calculate climate-related quantities.
```@docs
insolation
surface_atmosphere_contributions
total_toa_albedo
```

## Crash-course to DimensionalData.jl
```@docs
DimensionalData
```
