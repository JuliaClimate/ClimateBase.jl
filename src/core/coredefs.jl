using DimensionalData

using DimensionalData: @dim, AbDimArray, hasdim, Dimension, IndependentDim
using Dates

Time = DimensionalData.Ti

export At, Between, Near # Selectors from DimensionalArrays.jl
export hasdim, AbDimArray, DimensionalArray
export get_var_as_dimarray, allkeys
export Time, Lon, Lat, dims, Coord, Height
export EqArea, Grid, spacestructure, wrap_lon

@dim Lon IndependentDim "Longitude" "lon"
@dim Lat IndependentDim "Latitude" "lat"
@dim Coord IndependentDim "Coordinates"
@dim Categ IndependentDim "Category"
@dim Height IndependentDim "Height" "height"

const COMMONNAMES = Dict(
    "lat" => Lat,
    "latitude" => Lat,
    "lon" => Lon,
    "long" => Lon,
    "longitude" => Lon,
    "time" => Time,
    "height" => Height,
)

# the trait EqArea is for equal area grids. Functions can use the `spacestructure` and
# dispatch on `EqArea` or other types while still being type-stable
struct EqArea end
struct Grid end
spacestructure(a::AbDimArray) = spacestructure(dims(a))
function spacestructure(dims)
    if hasdim(dims, Coord)
        EqArea()
    elseif hasdim(dims, Lon) || hasdim(dims, Lat)
        Grid()
    else
        nothing
    end
end
