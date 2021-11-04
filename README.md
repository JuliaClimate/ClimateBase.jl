# ClimateBase.jl
| **Documentation**   |  **Tests**     | **CodeCov** |
|:--------:|:---------------:|:--------:|
|[![](https://img.shields.io/badge/docs-online-blue.svg)](https://JuliaClimate.github.io/ClimateBase.jl/dev)| [![CI](https://github.com/JuliaClimate/ClimateBase.jl/workflows/CI/badge.svg)](https://github.com/JuliaClimate/ClimateBase.jl/actions) | [![codecov](https://codecov.io/gh/JuliaClimate/ClimateBase.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaClimate/ClimateBase.jl) 

`ClimateBase` is a Julia package offering tools to analyze and manipulate climate (spatiotemporal) data. These data are dimensional & spatiotemporal but the corresponding dimensions all need special handling. For example the most common dimensions are longitude, latitude and time.
`ClimateBase` is structured to deal with these intricacies, and in addition offer several functionalities commonly used, and sought after, by climate scientists.

In the future it will serve as the base building block for `ClimateTools`, which offers more advanced functionalities, such as bias correction, interpolation and shapefiles  spatial subsetting.

ClimateBase.jl supports Julia 1.5+.
