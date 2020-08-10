# Introduction
Basic tools for dealing with climate (spatiotemporal) data.

This project treats "climate data" as a `ClimArray`, which uses the DimensionalData.jl interface and can be thought of as a syntactic equivalent to `DimensionalArray`.
A (very) brief introduction to DimensionalData.jl is copied here from its docs, because basic knowledge of how to handle a `ClimArray` is assumed in our docs.
DimensionalData.jl allows truly convenient handling of climate data, where it is important to be able to dimensionally-index data by their values. E.g. you can do
```@example main
using ClimateBase, Dates
Time = ClimateBase.Ti # more intuitive
lats = -90:5:90
lons = 0:10:359
t = Date(2000, 3, 15):Month(1):Date(2020, 3, 15)
dimensions = (Lon(lons), Lat(lats), Time(t))
A = ClimArray(rand(36, 37, 241), dimensions)
B = A[Lon(Between(0, 30)), Time(At(Date(2011,5,15)))]
```
and use convenience, physically-inspired functions that do automatic (and correct) statistical weighting, like
```@example main
C = latmean(B)
```
where in this averaging process each data point is weighted by the cosine of its latitude.

**Notice: at the moment the entirety of this package relies on doing operations in-memory. In the future, doing operations from-disk is planned.**

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
This type of space is called "grid". It is assumed throughout that the Longitude dimension precedes the Latitude dimension, and both are measured in **degrees**.

Another type of spatial coordinates is supported, and that is of **equal-area**. There, the spatial dimension is instead given by a single `Vector` of coordinate locations, i.e. 2-element `SVector(longitude, latitude)`.
Each point in this vector corresponds to a polygon (typically triangle or trapezoid) that covers equal amount of spatial area as any other point.
The actual limits of each polygon are not included in the dimension.
Physically inspired averaging functions, like [`spacemean`](@ref) below, work for both types of spatial coordinates.

**TODO: Finish this section.**

The function `spatialidxs` returns an iterator over the spatial coordinates of the data, and works for both types (grid or equal-area):
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
