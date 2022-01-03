# Advanced functionality
This page lists functionality that is currently in ClimateBase.jl but will be ported over to ClimateTools.jl once ClimateTools.jl has been properly updated to latest ClimateBase.jl.

## Vertical Interpolation
We offer 3 functions for vertical interpolation. Extrapolation results in missing values by default, but can also be linear (`extrapolation_bc = Line()`). Check the Interpolations.jl package for more information: https://juliamath.github.io/Interpolations.jl/latest/ For additional extrapolation boundary conditions: https://juliamath.github.io/Interpolations.jl/latest/extrapolation/

```@example main
D = ClimArray([1.0:11.0 2.0:12.0 3.0:2.0:23.0], (Hei(0.0:2000.0:20000.0), Ti(1:3)))
pressure_levels = [950, 850, 650, 350, 250, 150] .* 100.0
D_pre = interpolate_height2pressure(D, pressure_levels)
D_back = interpolate_pressure2height(D_pre, Vector(dims(D,Hei).val),extrapolation_bc=Line())

```

```@docs
interpolation2pressure
interpolate_height2pressure
interpolate_pressure2height
```



## Climate quantities
Functions that calculate climate-related quantities.
```@docs
insolation
surface_atmosphere_contributions
total_toa_albedo
```

## Timeseries Analysis
```@docs
sinusoidal_continuation
seasonal_decomposition
```
