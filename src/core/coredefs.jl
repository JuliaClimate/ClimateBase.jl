##########################################################################################
# Basic imports and dimension definitions
##########################################################################################
using DimensionalData
using DimensionalData: @dim, hasdim, Dimension, IndependentDim
using DimensionalData: basetypeof
using Dates

Time = DimensionalData.Ti

AbDimArray = DimensionalData.AbstractDimArray

"""
    dimindex(A::ClimArray, d) → i
Get the index of dimension `d`.
"""
dimindex(A, i::Int) = i
dimindex(A, dim) = DimensionalData.dimnum(A, dim)

export At, Between, Near # Selectors from DimensionalArrays.jl
export hasdim, dims, dimindex
export Time, Lon, Lat, dims, Coord, Hei, Pre, Ti
export UnstructuredGrid, LonLatGrid, spacestructure
export DimensionalData # for accessing its functions

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
# UnstructuredGrid, which is for equal area "grids" (or better points/coordinates)
# while the LonLatGrid is for standard Longitude x Latitude dimensions.
# Dispatch on `UnstructuredGrid` or other types is type-stable
# because it comes from the underlying dimension

abstract type SpaceType end

"""
Space coordinates are represented by two orthogonal dimensions `Lon, Lat`,
one being longitude and the other being latitude.
"""
struct LonLatGrid <: SpaceType end

"""
Space coordinates are represented by a single dimension `Coord`, whose
elements are coordinate locations, i.e. 2-element `SVector(longitude, latitude)`.
Each coordinate represents an **equal area polygon** corresponding to the point in space.
The actual limits of each polygon are not included in the dimension for performance reasons.

This dimension allows indexing according to the underlying `Lon, Lat` representation,
e.g. you can do
```julia
A # some `ClimArray` with unstructured grid type.
A[Coord(Lon(Between(0, 30)), Lat(Between(-30, 30)))]
```

To use functions such as [`zonalmean`](@ref) or [`hemispheric_means`](@ref) with this grid,
you must first sort the `ClimArray` so that the latitudes
of its coordinates are sorted in ascending order. I.e.
```julia
A # some `ClimArray` with unstructured grid type.
coords = dims(A, Coord).val
si = sortperm(coords, by = reverse)
A = A[Coord(si)]
```
**This is done automatically by [`ncread`](@ref).**

!!! warn
    `UnstructuredGrid` functionality is currently in an **experimental phase**!
    Notice that non-equal area unstructured grids are not supported yet.
"""
struct UnstructuredGrid <: SpaceType end

spacestructure(a::AbDimArray) = spacestructure(dims(a))
function spacestructure(dims)
    if hasdim(dims, Coord)
        UnstructuredGrid()
    elseif hasdim(dims, Lon) || hasdim(dims, Lat)
        LonLatGrid()
    else
        error("Array does not have any spatial dimensions: `Lon`, `Lat`, or `Coord`.")
    end
end

##########################################################################################
# ClimArray definition and DimensionalData.jl extensions
##########################################################################################
export ClimArray
struct ClimArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Me} <: AbstractDimensionalArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Symbol
    attrib::Me
end
ClimArray(A::DimensionalArray) = ClimArray(A.data, A.dims, A.refdims, A.name, A.metadata)
ClimArray(A::ClimArray, dims::Tuple = A.dims; name = A.name, attrib = A.attrib) =
ClimArray(A.data, DimensionalData.formatdims(A.data, dims), A.refdims, Symbol(name), attrib)

"""
    ClimArray(A::Array, dims::Tuple; name = "", attrib = nothing)

`ClimArray` is a structure that contains numerical array data bundled with dimensional
information, a name and an `attrib` field (typically a dictionary) that holds general
attributes.
You can think of `ClimArray` as a in-memory representation of a CFVariable.

At the moment, a `ClimArray` is using `DimensionalArray` from DimensionalData.jl, and
all basic handling of `ClimArray` is offered by `DimensionalData` (see below).

`ClimArray` is created by passing in standard array data `A` and a
tuple of dimensions `dims`. See [`ncread`](@ref) to automatically create a `ClimArray`
from a .nc file.

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
    ClimArray(A, DimensionalData.formatdims(A, dims), refdims, Symbol(name), attrib)
ClimArray(A::AbstractArray, dims::Tuple, name; refdims=(), attrib=nothing) =
    ClimArray(A, DimensionalData.formatdims(A, dims), refdims, Symbol(name), attrib)

Base.parent(A::ClimArray) = A.data
Base.@propagate_inbounds Base.setindex!(A::ClimArray, x, I::Vararg{DimensionalData.StandardIndices}) =
    setindex!(A.data, x, I...)

DimensionalData.metadata(A::ClimArray) = A.attrib
DimensionalData.rebuild(A::ClimArray, data::Any, dims::Tuple=dims(A), refdims=DimensionalData.refdims(A),
name="", attrib=nothing) = ClimArray(data, dims, refdims, Symbol(name), attrib)
DimensionalData.basetypeof(::ClimArray) = ClimArray

DimensionalData.rebuild(
    A::ClimArray;
    data=data(A), dims=dims(A), refdims=refdims(A), name=name(A), metadata=metadata(A)
) = ClimArray(data, dims, refdims, name, metadata)


# The following basic methods allow indexing with tuples, (Time(5), Lon(3))
Base.getindex(A::ClimArray, i::Tuple) = A[i...]
Base.setindex!(A::ClimArray, x, i::Tuple) = setindex!(A, x, i...)
Base.view(A::ClimArray, i::Tuple) = view(A, i...)

# Convenience
Base.ones(A::AbDimArray) = basetypeof(A)(ones(eltype(A), size(A)), dims(A))
Base.zeros(A::AbDimArray) = basetypeof(A)(zeros(eltype(A), size(A)), dims(A))
