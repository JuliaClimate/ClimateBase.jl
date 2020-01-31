module ClimateBase

using AxisArrays
using ArgCheck
using Images
using Dates

# TYPES

"""
    ClimGrid{A <: AxisArray}
In-memory representation of Climate Forecast netCDF files.
struct ClimGrid\n
  data::AxisArray # Data\n
  longrid::AbstractArray{N,2} where N # the longitude grid\n
  latgrid::AbstractArray{N,2} where N # the latitude grid\n
  msk::Array{N, 2} where N # Data mask (NaNs and 1.0)\n
  grid_mapping::Dict#{String, Any} # bindings for native grid\n
  dimension_dict::Dict\n
  model::String\n
  frequency::String # Day, month, years\n
  experiment::String # Historical, RCP4.5, RCP8.5, etc.\n
  run::String\n
  project::String # CORDEX, CMIP5, etc.\n
  institute::String # UQAM, DMI, etc.\n
  filename::String # Path of the original file\n
  dataunits::String # Celsius, kelvin, etc.\n
  latunits::String # latitude coordinate unit\n
  lonunits::String # longitude coordinate unit\n
  variable::String # Type of variable (i.e. can be the same as "typeofvar", but it is changed when calculating indices)\n
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)\n
  typeofcal::String # Calendar type\n
  timeattrib::Dict # Time attributes (e.g. days since ... )\n
  varattribs::Dict # Variable attributes dictionary\n
  globalattribs::Dict # Global attributes dictionary\n
end\n
"""
struct ClimGrid{A <: AxisArray}
  data::A
  longrid::Array{N,T} where T where N
  latgrid::Array{N,T} where T where N
  msk::Array{N,T} where T where N
  grid_mapping::Dict # information of native grid
  dimension_dict::Dict
  timeattrib::Dict
  model::String
  frequency::String
  experiment::String
  run::String
  project::String
  institute::String
  filename::String
  dataunits::String
  latunits::String # of the coordinate variable
  lonunits::String # of the coordinate variable
  variable::String # Type of variable
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)
  typeofcal::String # Calendar type
  varattribs::Dict # Variable attributes
  globalattribs::Dict # Global attributes

end

"""
    ClimGrid(data; longrid=[], latgrid=[], msk=[], grid_mapping=Dict(), dimension_dict=Dict(), model="NA", frequency="NA", experiment="NA", run="NA", project="NA", institute="NA", filename="NA", dataunits="NA", latunits="NA", lonunits="NA", variable="NA", typeofvar="NA", typeofcal="NA", varattribs=Dict(), globalattribs=Dict())
Constructor of the ClimGrid function. Data is an AxisArray. Everything else is optional, but usually needed for further processing (mapping, interpolation, etc...).
struct ClimGrid\n
  data::AxisArray # Data \n
  longrid::AbstractArray{N,2} where N # the longitude grid \n
  latgrid::AbstractArray{N,2} where N # the latitude grid \n
  msk::Array{N, 2} where N # Data mask (NaNs and 1.0) \n
  grid_mapping::Dict#{String, Any} # bindings for native grid \n
  dimension_dict::Dict\n
  model::String\n
  frequency::String # Day, month, years\n
  experiment::String # Historical, RCP4.5, RCP8.5, etc.\n
  run::String\n
  project::String # CORDEX, CMIP5, etc.\n
  institute::String # UQAM, DMI, etc.\n
  filename::String # Path of the original file\n
  dataunits::String # Celsius, kelvin, etc.\n
  latunits::String # latitude coordinate unit\n
  lonunits::String # longitude coordinate unit\n
  variable::String # Type of variable (i.e. can be the same as "typeofvar", but it is changed when calculating indices)\n
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)\n
  typeofcal::String # Calendar type\n
  timeattrib::Dict # Time attributes (e.g. days since ... )\n
  varattribs::Dict # Variable attributes dictionary\n
  globalattribs::Dict # Global attributes dictionary\n
end\n
"""
function ClimGrid(data; longrid=[], latgrid=[], msk=[], grid_mapping=Dict(), dimension_dict=Dict(), timeattrib=Dict(), model="NA", frequency="NA", experiment="NA", run="NA", project="NA", institute="NA", filename="NA", dataunits="NA", latunits="NA", lonunits="NA", variable="NA", typeofvar="NA", typeofcal="NA", varattribs=Dict(), globalattribs=Dict())

    if isempty(dimension_dict)
        dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
    end

    if isempty(msk)
        msk = Array{Float64}(ones((size(data, 1), size(data, 2))))
    end

    ClimGrid(data, longrid, latgrid, msk, grid_mapping, dimension_dict, timeattrib, model, frequency, experiment, run, project, institute, filename, dataunits, latunits, lonunits, variable, typeofvar, typeofcal, varattribs, globalattribs)
end

include("functions.jl")

export ClimGrid
export periodmean
export temporalsubset
export get_timevec
export buildtimetype
export timeindex
export buildarray_climato
export buildarrayinterface
export buildarray_annual
export buildarray_resample


end # module
