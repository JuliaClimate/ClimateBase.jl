# Introduction
Basic tools for dealing with climate (spatiotemporal) data.

This project treats "climate data" as a `ClimArray`, which uses the DimensionalData.jl interface and can be thought of as a syntactic equivalent to `DimensionalArray`.
A (very) brief introduction to DimensionalData.jl is copied here from its docs, for clarity.

## Making a `ClimArray`
You can create a `ClimArray` yourself, or you can load data from an `.nc` file with CF-conventions, using `ClimArray`:
```@docs
ClimArray
```

Notice that (at the moment) we use a pre-defined mapping of common names to proper dimensions - please feel free to extend the following via a Pull Request:
```@example
using ClimateBase # hide
ClimateBase.COMMONNAMES
```

## Types of spatial coordinates
Most of the time the spatial information of your data is in the form of a Longitude × Latitude grid. This is simply achieved via the existence of two dimensions (`Lon, Lat`) in your dimensional data array. Height, although representing physical space as well, is not considered part of the "spatial dimensions", and is treated as any other additional dimension.

Another type of spatial coordinates is supported, and that is of **equal-area**. There, the spatial dimension is instead given by a single `Vector` of coordinate locations, i.e. 2-element `SVector(longitude, latitude)`. Each point in this vector corresponds to a polygon (typically triangle or trapezoid) that covers equal amount of spatial area as any other point.

**TODO: Finish this section.**

The function `spatialidxs` returns an iterator over the spatial coordinates of the data, and works for both types (grid or equal-area). It is used like so:
```julia
for i in spatialidxs(A)
    x = A[i...] # each possible location slice
end
```


## Aggregation
Physically-inspired aggregation:
```@docs
zonalmean
latmean
spacemean
timemean
hemispheric_means
hemispheric_functions
```

General aggregation:
```@docs
dropagg
collapse
```

## Time-handling
```@docs
maxyearspan
yearly
monthspan
temporal_sampling
```

## Statistical weighting

General:
```@docs
dimwise
```
**TODO:  complete this section**

## Timeseries Analysis
```@docs
sinusoidal_continuation
seasonal_decomposition
```

## Crash-course to DimensionalData.jl
```@docs
DimensionalData
```
---

There are some standard dimensions that we use within our functions, and it as assumed that if your data contains such dimensions (in concept) then they must match the appropriate type:
```@example
using ClimateBase, DimensionalData # hide
for D in ClimateBase.ALLDIMS
    println(D, " (full name = $(DimensionalData.name(D)))")
end
```
We explicitly assume that `Lon, Lat` are measured in degrees and not radians or meters (extremely important for spatial averaging processes).
