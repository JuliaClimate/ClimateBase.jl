# ClimateBase.jl

ClimateBase implement the basic type that are used by [ClimateTools](https://github.com/Balinus/ClimateTools.jl) and [ClimatePlots](https://github.com/JuliaClimate/ClimatePlots.jl).

## ClimGrid

The `ClimGrid` type is a in-memory representation of a CF-compliant netCDF file for a single variable.

```julia
struct ClimGrid
  data::AxisArray # labeled axis
  longrid::AbstractArray{N,2} where N # the longitude grid
  latgrid::AbstractArray{N,2} where N # the latitude grid
  msk::Array{N, 2} where N
  grid_mapping::Dict # bindings of native grid
  dimension_dict::Dict
  model::String
  frequency::String
  experiment::String
  run::String
  project::String # CORDEX, CMIP5, etc.
  institute::String
  filename::String
  dataunits::String
  latunits::String # of the coordinate variable
  lonunits::String # of the coordinate variable
  variable::String # Type of variable (i.e. can be the same as "var", but it is changed when calculating indices)
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)
  typeofcal::String # Calendar type
  timeattrib::Dict # Time attributes
  varattribs::Dict # Variable attributes
  globalattribs::Dict # Global attributes

end
```
