# Statistics

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
realtime_days
realtime_milliseconds
seasonality
sametimespan
```

## Spatial

All functions in this section work for both types of space, see [Types of spatial information](@ref).
```@docs
zonalmean
latmean
spacemean
spaceagg
hemispheric_means
hemispheric_functions
tropics_extratropics
lonlatfirst
longitude_circshift
```

### Types of spatial information
Spatial information (excluding height/pressure dimensions) in ClimateBase.jl exists in one of two forms:
```@docs
OrthogonalSpace
CoordinateSpace
```

ClimateBase.jl works with either type of spatial coordinate system.
Therefore, physically inspired averaging functions, like [`spacemean`](@ref) or [`zonalmean`](@ref), work for both types of spatial coordinates.
In addition, the function `spatialidxs` returns an iterator over the spatial coordinates of the data, and works for both types (grid or equal-area):
```@docs
spatialidxs
```

[`ncread`](@ref) tries to automatically deduce the correct space type and create the appropriate dimension.


## General aggregation
The physical averages of the previous section are done by taking advantage of a general aggregation syntax, which works with any aggregating function like `mean, sum, std`, etc.
```@docs
dropagg
collapse
```

## Missing data
When loading an array with [`ncread`](@ref), the values of the return array may contain missing values if the actual data contain missing values according to the CF-standards.
In other packages or other programming languages these missing values are handled "internally" and e.g. in statistical operations like `mean`, the statistics explicitly skip over missing values.
For example this is a typical workflow of creating an array, assigning `missing` to all values of an array over land, and then taking the `mean` of the array, which would be the "mean over ocean".

ClimateBase.jl _does not_ follow this approach for two reasons: 1) it does not comply with [Julia's `missing` propagation logic](https://docs.julialang.org/en/v1/manual/missing/), 2) using proper statistical weights gives more power to the user. As you have already seen in the documentation strings of e.g. [`timeagg`](@ref), [`spaceagg`](@ref) or [`dropagg`](@ref), you can provide explicit statistical weights of various forms.
This gives you more power, because in the case of `missing` your statistical weights can only be 0 (missing value) or 1 (non-missing value). As an example, "pixel" of your spatial grid will have ambiguous values if it is not 100% covered by ocean, and to do a _proper_ average over ocean you should instead provide weights `W` whose value is quite simply the ocean fraction of each pixel.

But what if you already have an array with `missing` values and you want to do what was described in the beginning, e.g. average by skipping the missings? Do not worry, we have you covered! Use the function [`missing_weights`](@ref)! See also [`sinusoidal_continuation`](@ref) if the missing values are only in a subset of your temporal coverage.

```@docs
missing_weights
missing_val
```
