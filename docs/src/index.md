# Introduction
`ClimateBase` is a Julia package offering basic tools for analyzing data that are typically in the form used by climate sciences.
These data are dimensional spatiotemporal but the corresponding dimensions all need special handling. For example the most common dimensions are longitude, latitude and time.

* longitude is by definition a periodic dimension
* latitude is a linear dimension. However because the coordinate system most often used in climate sciences is a grid of longitude × latitude (in equal degrees) the area element of space depends on latitude and this needs to be taken into account.
* time is a linear dimension *in principle*, but its values are `<: AbstractDateTime` instead of `<: Real`. The human calendar (where these values come from) is periodic but each period may not correspond to the same physical time, and this also needs to be taken into account.

`ClimateBase` is structured to deal with these difficulties, and in addition offer several functionalities commonly used, and sought after, by climate scientists.
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

## Making a `ClimArray`
You can create a `ClimArray` yourself, or you can load data from an `.nc` file with CF-conventions, using `ClimArray`:
```@docs
ClimArray
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


## Types of spatial coordinates
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

## Physical averages
```@docs
zonalmean
latmean
spacemean
timemean
hemispheric_means
hemispheric_functions
```

## Aggregation
The physical averages of the previous section are done by taking advantage of a general aggregation syntax, which works with any aggregating function like `mean, sum, std`, etc.
```@docs
dropagg
collapse
```

### Statistical weighting
Functons like `timemean` and `spacemean` perform statistically-proper averaging by weighting each data point by its length in time or by the cosine of its latitude.
Functional versions that can explicitly take *extra* weights (e.g. if you want to weight your data with e.g. the ice area fraction) are provided:
```@docs
timeagg
spaceagg
```
**TODO: I need to improve the statistical weighting to work with all possible dimensional combinations.**

These functions internally use the following function that can be used also outside of statistical weighting:
```@docs
dimwise
```


## Time-handling
```@docs
maxyearspan
yearly
monthspan
temporal_sampling
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

---

There are some standard dimensions that we use within our functions, and it as assumed that if your data contains such dimensions (in concept) then they must match the appropriate type:
```@example
using ClimateBase, DimensionalData # hide
for D in ClimateBase.STANDARD_DIMS
    println(D, " (full name = $(DimensionalData.name(D)))")
end
```
We explicitly assume that `Lon, Lat` are measured in degrees and not radians or meters (extremely important for spatial averaging processes).
