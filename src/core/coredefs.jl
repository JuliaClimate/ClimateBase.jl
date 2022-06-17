##########################################################################################
# Basic imports and dimension definitions
##########################################################################################
using DimensionalData
using DimensionalData.Dimensions, DimensionalData.LookupArrays
using DimensionalData: basetypeof, broadcast_dims
using DimensionalData.Dimensions: setdims
using Dates

Time = DimensionalData.Ti
Tim = DimensionalData.Ti

AbDimArray = DimensionalData.AbstractDimArray

"""
    dimindex(A::ClimArray, d) → i
Get the index of dimension `d`.
"""
dimindex(A, i::Int) = i
dimindex(A, dim) = DimensionalData.dimnum(A, dim)

"""
    gnv(object) → x
Short for "get numeric value", this function will return the pure numeric value(s)
of the given object. Convenience function for quickly getting the numeric data of
dimensional arrays or dimensions.
"""
gnv(x) = x
gnv(x::Union{AbDimArray, LookupArray}) = parent(x)
gnv(x::Dimension) = parent(parent(x))

export At, (..), Between, Near # Selectors from DimensionalArrays.jl
export hasdim, dims, dimindex
export Time, Lon, Lat, dims, Coord, Hei, Pre, Ti, Tim
export CoordinateSpace, OrthogonalSpace, spacestructure
export DimensionalData # for accessing its functions
export setdims
export gnv

@dim Lon IndependentDim "Longitude"
@dim Lat IndependentDim "Latitude"
@dim Hei IndependentDim "Height"
@dim Pre IndependentDim "Pressure"

STANDARD_DIMS = (Lon, Lat, Time, Hei, Pre, Coord)

"""
    COMMONNAMES
A dictionary of common names of dimensions (as strings) to actual dimension types.
"""
const COMMONNAMES = Dict(
    "lat" => Lat,
    "latitude" => Lat,
    "y" => Lat,
    "yc" => Lat,
    "lon" => Lon,
    "long" => Lon,
    "x" => Lon,
    "xc" => Lon,
    "longitude" => Lon,
    "time" => Time,
    "t" => Time,
    "age" => Time,
    "height" => Hei,
    "altitude" => Hei,
    "pressure" => Pre,
    "level" => Pre,
)

##########################################################################################
# Space types
##########################################################################################

# the following traits for the the way space is configured. currently the options are
# CoordinateSpace, which is for points/coordinates
# while the OrthogonalSpace is for standard Longitude x Latitude dimensions.
# Dispatch on `CoordinateSpace` or other types is type-stable
# because it comes from the underlying dimension

abstract type SpaceType end

"""
Space information is represented by two orthogonal dimensions `Lon, Lat`,
one being longitude and the other being latitude.
"""
struct OrthogonalSpace <: SpaceType end

"""
Space information is represented by a single dimension `Coord`, whose
elements are coordinates, i.e. 2-element `SVector(longitude, latitude)`.
Each coordinate represents the center of an arbitrary polygon in space.
The actual limits of each polygon are not included in the dimension for performance reasons.

In statistical functions such as [`spaceagg`](@ref), it is assumed that entry of the
coordinates covers **an equal amount of area**. If this is not the case, you can simply
provide an additional weights vector which would correspond to the area covered.

This dimension also allows indexing by latitude ranges, e.g. you can do
```julia
A # some `ClimArray` with a `Coord` dimension
A[Coord(Lat(-30..30)))]
```

Most functions of ClimateBase.jl implicitly assume that the coordinates are
sorted by latidude. You can achieve this with the following code:
```julia
A # some `ClimArray` with a `Coord` dimension
coords = gnv(dims(A, Coord))
si = sortperm(coords; by = reverse)
A = A[Coord(si)]
```
**This is done automatically by [`ncread`](@ref).**
"""
struct CoordinateSpace <: SpaceType end

spacestructure(a::AbDimArray) = spacestructure(dims(a))
function spacestructure(dims)
    if hasdim(dims, Coord)
        CoordinateSpace()
    elseif hasdim(dims, Lon) || hasdim(dims, Lat)
        OrthogonalSpace()
    else
        error("Array does not have any spatial dimensions: `Lon`, `Lat`, or `Coord`.")
    end
end

##########################################################################################
# ClimArray definition and DimensionalData.jl extensions
##########################################################################################
# TODO: remove `refdims` field
export ClimArray
struct ClimArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Me} <: AbDimArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Symbol
    attrib::Me
end
ClimArray(A::AbDimArray) = ClimArray(A.data, A.dims, A.refdims, A.name, A.metadata)
ClimArray(A::ClimArray, dims::Tuple = A.dims; name = A.name, attrib = A.attrib) =
ClimArray(A.data, format(dims, A.data), A.refdims, Symbol(name), attrib)

"""
    ClimArray(A::Array, dims::Tuple; name = "", attrib = nothing)

`ClimArray` is a structure that contains numerical array data bundled with dimensional
information, a name and an `attrib` field (typically a dictionary) that holds general
attributes.
You can think of `ClimArray` as a in-memory representation of a CFVariable.

At the moment, a `ClimArray` is using `DimArray` from DimensionalData.jl, and
all basic handling of `ClimArray` is offered by `DimensionalData` (see below).

`ClimArray` is created by passing in standard array data `A` and a
tuple of dimensions `dims`. See [`ncread`](@ref) to automatically create a `ClimArray`
from a .nc file. For obtaining the raw numeric values of a `ClimArray` or any of its
dimensions, use the function [`gnv`](@ref).

## Example
```julia
using ClimateBase, Dates
Time = ClimateBase.Ti # more intuitive name for time dimension
lats = -90:5:90
lons = 0:10:359
t = Date(2000, 3, 15):Month(1):Date(2020, 3, 15)
# dimensional information:
dimensions = (Lon(lons), Lat(lats), Time(t))
data = rand(36, 37, 241) # numeric data
A = ClimArray(data, dimensions)
```
"""
ClimArray(A::AbstractArray, dims::Tuple; refdims=(), name="", attrib=nothing) =
    ClimArray(A, format(dims, A), refdims, Symbol(name), attrib)
ClimArray(A::AbstractArray, dims::Tuple, name; refdims=(), attrib=nothing) =
    ClimArray(A, format(dims, A), refdims, Symbol(name), attrib)

Base.parent(A::ClimArray) = A.data

DimensionalData.metadata(A::ClimArray) = A.attrib
DimensionalData.basetypeof(::ClimArray) = ClimArray

DimensionalData.rebuild(
    ::ClimArray, data, dims, refdims, name, metadata
) = ClimArray(data, format(dims, data), refdims, name, metadata)
DimensionalData.rebuild(
    A::ClimArray;
    data=gnv(A), dims=dims(A), refdims=refdims(A), name=name(A), metadata=metadata(A)
) = ClimArray(data, format(dims, data), refdims, name, metadata)


# The following basic methods allow indexing with tuples, (Time(5), Lon(3))
Base.getindex(A::ClimArray, i::Tuple) = A[i...]
Base.setindex!(A::ClimArray, x, i::Tuple) = setindex!(A, x, i...)
Base.view(A::ClimArray, i::Tuple) = view(A, i...)

# Convenience
Base.ones(A::AbDimArray) = basetypeof(A)(ones(eltype(A), size(A)), dims(A))
Base.zeros(A::AbDimArray) = basetypeof(A)(zeros(eltype(A), size(A)), dims(A))
