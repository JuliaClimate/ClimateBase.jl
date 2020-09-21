# ClimateBase.jl
| **Documentation**   |  **Tests**     |
|:--------:|:---------------:|
|[![](https://img.shields.io/badge/docs-online-blue.svg)](https://JuliaClimate.github.io/ClimateBase.jl/dev)| [![Build Status](https://travis-ci.org/JuliaClimate/ClimateBase.jl.svg?branch=master)](https://travis-ci.org/JuliaClimate/ClimateBase.jl)

`ClimateBase` is a Julia package offering basic functionality for analyzing data that are typically in the form used by climate sciences.
These data are dimensional & spatiotemporal but the corresponding dimensions all need special handling.
For example the most common dimensions are longitude, latitude and time.
`ClimateBase` is structured to deal with these intricacies, and in addition offer several functionalities commonly used, and sought after, by climate scientists.
It also serves as the base building block for `ClimateTools`, which offers more advanced functionalities.
