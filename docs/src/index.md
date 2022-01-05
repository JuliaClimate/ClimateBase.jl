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
At the moment `ClimArray` is a subtype of `DimArray` from DimensionalData.jl.
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
gnv
```
It is strongly recommended to use the dimensions we export (because we dispatch on them and use their information):
```@example
using ClimateBase, DimensionalData # hide
for D in ClimateBase.STANDARD_DIMS
    println(D, " (full name = $(DimensionalData.name(D)))")
end
```
We explicitly assume that `Lon, Lat` are measured in degrees and not radians or meters (extremely important for spatial averaging processes).

## Crash-course to DimensionalData.jl
```@docs
DimensionalData
```

## Available selectors
| Selector                | Description                                                         |
| :---------------------- | :------------------------------------------------------------------ |
| [`At(x)`]               | get the index exactly matching the passed in value(s)               |
| [`Near(x)`]             | get the closest index to the passed in value(s)                     |
| [`Contains(x)`]         | get indices where the value x falls within an interval              |
| [`Where(f)`]            | filter the array axis by a function of the dimension index values.  |
| [`a..b`]                | get all indices between two values, inclusively.                    |
| [`OpenInterval(a, b)`]  | get all indices between `a` and `b`, exclusively.                   |
| [`Interval{A,B}(a, b)`] | get all indices between `a` and `b`, as `:closed` or `:open`.       |
