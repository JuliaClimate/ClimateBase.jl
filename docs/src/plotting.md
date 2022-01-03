# Plotting

## Native Julia plotting with GeoMakie.jl
*TODO* Plotting now works with GeoMakie.jl

## Python based plotting with `cartopy`
For now, you can use PyCall.jl, matplotlib, and the Python library cartopy.
In the file [`ClimateBase/plotting/python.jl`](https://github.com/JuliaClimate/ClimateBase.jl/tree/master/plotting/python.jl) we provide two functions that plot maps of `ClimArray` in arbitrary projections: `earthsurface` for `LonLatGrid` and `earthscatter` for `UnstructuredGrid`. You can incorporate these in your source code as a temporary solution.
