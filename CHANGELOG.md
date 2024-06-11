# 0.17.0

- New function `dailyagg`
- Many bug fixes in terms of compatibility with other packages
- Usage of Julia Package extensions for GeoMakie integration

# 0.16.3
- New functions `value_space, quantile_space`.
- `globalattr` has been renamed to `ncglobalattr`.

# 0.16
- Data with `Coord` dimension can now be saved as .nc files with `ncwrite`.
- When loading data with `Coord`, the longitude was automatically wrapped to -180 to 180 degrees. There was no reason for this, now no wrapping is done, the user may want to do it on their own after loading.

# 0.15
- Removed functions for calculating some decompositions of albedo to atmosphere/surface contributions. They were not basic enough for this package.


_Changelog for ClimateBase.jl is kept from version 0.14 onwards_