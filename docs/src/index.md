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
ClimArray(::Union{String, Vector{String}})
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
climarrays_to_nc
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

## Temporal
Functions related with the `Time` dimension.
```@docs
timemean
timeagg
monthlyagg
yearlyagg
seasonalyagg
temporalrange
maxyearspan
temporal_sampling
time_in_days
```

## Spatial

### Spatial aggregation
All functions in this section work for both types of space, see [Types of spatial coordinates](@ref).
```@docs
zonalmean
latmean
spacemean
spaceagg
hemispheric_means
hemispheric_functions
lonlatfirst
```

### Types of spatial coordinates
Most of the time the spatial information of your data is in the form of a Longitude × Latitude grid. This is simply achieved via the existence of two dimensions (`Lon, Lat`) in your dimensional data array.
This type of space is called `LonLatGrid`. It is assumed throughout that longitude and latitude are measured in **degrees**.
Height, although representing physical space as well, is not considered part of the "spatial dimensions", and is treated as any other additional dimension.

Another type of spatial coordinates is supported, and that is of **equal-area**.
Currently only a single type, `GaussianEqualArea`, exists for this purpose, which represents coordinates in a Gaussian grid as shown here: https://en.wikipedia.org/wiki/Gaussian_grid.
In `GaussianEqualArea` the spatial information is instead given by single dimension whose elements are coordinate locations, i.e. 2-element `SVector(longitude, latitude)`.
The dimension name is `Coord`.
Each point in this dimension corresponds to a polygon (for `GaussianEqualArea` a trapezoid) that covers equal amount of spatial area as any other point.
The actual limits of each polygon are not included in the dimension for performance reasons.

ClimateBase.jl works with either type of spatial coordinate system.
Therefore, physically inspired averaging functions, like [`spacemean`](@ref) or [`zonalmean`](@ref), work for both types of spatial coordinates.
In addition, the function `spatialidxs` returns an iterator over the spatial coordinates of the data, and works for both types (grid or equal-area):
```@docs
spatialidxs
```

### Equal area creation

!!! warn
    Equal area functionality is currently in an **experimental phase**!
    You can try using the function `ClimArray_eqarea` to load a `ClimArray` with Gaussian grid directly from a `.nc` file. This function assumes that this grid was created with CDO, using e.g. `cdo remapbil,gea250 IN.nc OUT.nc`.


To manually make a Gaussian grid `ClimArray`, try the following approach:
```julia
file = NCDataset("some_file_with_eqarea.nc")
# the following lines make the coordinates of the equal area, which depends on
# how your .nc file is structured
lons = Array(file["lon"])
lats = Array(file["lat"])
coords = [SVector(lo, la) for (lo, la) in zip(lons, lats)]
# Sort points by their latitude (important!)
si = sortperm(coords, by = reverse)
# Load some remaining dimensions and create the proper `Dimension` tuple:
t = Array(file["time"])
dimensions = (Coord(coords), Time(t))
# Finally load the array data and make a ClimArray
data = Array(file["actual_data_like_radiation"])
data = data[si, :] # permute like the coordinate
A = ClimArray(data, dimensions)
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

## Plotting
Currently ClimateBase.jl does not have integrated plotting support. In the near future it will have this based on the upcoming GeoMakie.jl.

For now, you can use PyCall.jl, matplotlib, and the Python library cartopy.
In the file [`ClimateBase/plotting/python.jl`](https://github.com/JuliaClimate/ClimateBase.jl/tree/master/plotting/python.jl) we provide two functions that plot maps of `ClimArray` in arbitrary projections: `earthsurface` for `LonLatGrid` and `earthscatter` for `GaussianEqualArea`. You can incorporate these in your source code as a temporary solution.

## Crash-course to DimensionalData.jl
```@docs
DimensionalData
```
