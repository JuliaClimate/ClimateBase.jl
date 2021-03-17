##########################################################################################
# Basic imports and dimension definitions
##########################################################################################
using DimensionalData
using DimensionalData: @dim, hasdim, Dimension, IndependentDim
using DimensionalData: basetypeof
using Dates

Time = DimensionalData.Ti

AbDimArray = DimensionalData.AbstractDimArray

export At, Between, Near # Selectors from DimensionalArrays.jl
export hasdim, DimensionalArray, dimnum
export Time, Lon, Lat, dims, Coord, Hei, Pre
export GaussianEqualArea, LonLatGrid, spacestructure
export DimensionalData # for accessing its functions

@dim Lon IndependentDim "Longitude"
@dim Lat IndependentDim "Latitude"
@dim Coord IndependentDim "Spatial Coordinates"
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
    "lon" => Lon,
    "long" => Lon,
    "longitude" => Lon,
    "time" => Time,
    "height" => Hei,
    "altitude" => Hei,
    "pressure" => Pre,
    "level" => Pre,
)

# the following traits for the the way space is configured. currently the options are
# GaussianEqualArea, which is for equal area "grids" (or better points/coordinates)
# while the LonLatGrid is for standard Longitude x Latitude dimensions.
# Dispatch on `GaussianEqualArea` or other types is type-stable
# because it comes from the underlying dimension
abstract type SpaceType end
struct GaussianEqualArea <: SpaceType end
struct LonLatGrid <: SpaceType end
spacestructure(a::AbDimArray) = spacestructure(dims(a))
function spacestructure(dims)
    if hasdim(dims, Coord)
        GaussianEqualArea()
    elseif hasdim(dims, Lon) || hasdim(dims, Lat)
        LonLatGrid()
    else
        nothing
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

`ClimArray` is created by passing in standard array data `A` and a tuple of dimensions `dims`.

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

##########################################################################################
# Pretty printing
##########################################################################################
# Remove reference dims from printing, and show attributes if any
function Base.show(io::IO, A::ClimArray)
    summary(io, A)
    print(io, "and")
    printstyled(io, " data: "; color=:green)
    dataA = data(A)
    print(io, summary(dataA), "\n")
    x = 2length(dims(A)) + attriblength(A.attrib) + 5
    custom_show(io, data(A), x)
end

attriblength(d::AbstractDict) = length(d)
attriblength(d) = 0

# Define summary
function Base.summary(io::IO, A::ClimArray)
    l = nameof(typeof(A))
    printstyled(io, nameof(typeof(A)); color=:blue)
    if A.name ≠ Symbol("")
        print(io, " (named ")
        printstyled(io, A.name; color=:blue)
        print(io, ")")
    end

    print(io, " with dimensions:\n")
    for d in dims(A)
        print(io, " ", d, "\n")
    end
    if !isnothing(A.attrib)
        printstyled(io, "attributes: "; color=:magenta)
        show(io, MIME"text/plain"(), A.attrib)
        print(io, '\n')
    end
end

# Thanks to Michael Abbott for the following function
function custom_show(io::IO, A::AbstractArray{T,0}, x) where T
    Base.show(IOContext(io, :compact => true, :limit => true), A)
end
function custom_show(io::IO, A::AbstractArray{T,1}, x) where T
    Base.show(IOContext(io, :compact => true, :limit => true, :displaysize => displaysize(io) .- (x, 0)), A)
end
function custom_show(io::IO, A::AbstractArray{T,2}, x) where T
    Base.print_matrix(IOContext(io, :compact => true, :limit => true, :displaysize => displaysize(io) .- (x, 0)), A)
end
function custom_show(io::IO, A::AbstractArray{T,N}, x) where {T,N}
    o = ones(Int, N-2)
    frame = A[:, :, o...]
    onestring = join(o, ", ")
    println(io, "[:, :, $(onestring)]")
    Base.print_matrix(
        IOContext(io, :compact => true, :limit=>true, :displaysize => displaysize(io) .- (x, 0)),
        frame)
    print(io, "\n[and ", prod(size(A,d) for d=3:N) - 1," more slices...]")
end
