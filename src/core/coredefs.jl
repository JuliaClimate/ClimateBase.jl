using DimensionalData

using DimensionalData: @dim, AbDimArray, hasdim, Dimension, IndependentDim
using DimensionalData: basenameof
using Dates

Time = DimensionalData.Ti

export At, Between, Near # Selectors from DimensionalArrays.jl
export hasdim, AbDimArray, DimensionalArray
export get_var_as_dimarray
export Time, Lon, Lat, dims, Coord, Hei, Ti
export EqArea, Grid, spacestructure, wrap_lon

@dim Lon IndependentDim "Longitude" "lon"
@dim Lat IndependentDim "Latitude" "lat"
@dim Coord IndependentDim "Coordinates (spatial)"
@dim Hei IndependentDim "Height" "height"
@dim Pre IndependentDim "Pressure" "pressure"

STANDARD_DIMS = (Lon, Lat, Time, Hei, Pre, Coord)

const COMMONNAMES = Dict(
    "lat" => Lat,
    "latitude" => Lat,
    "lon" => Lon,
    "long" => Lon,
    "longitude" => Lon,
    "time" => Time,
    "height" => Hei,
    "pressure" => Pre,
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

# ClimArray type
export ClimArray
struct ClimArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na<:AbstractString,Me} <: AbstractDimensionalArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    attrib::Me
end
ClimArray(A::DimensionalArray) = ClimArray(A.data, A.dims, A.refdims, A.name, A.metadata)

"""
    ClimArray(A, dims::Tuple; name = "", attrib = nothing)

`ClimArray` is a simple wrapper of a standard dimensional array from DimensionalData.jl
bundled with an extra `attrib` field (typically a dictionary) that holds general attributes.

`ClimArray` is created by passing in standard array data `A` and a tuple of dimensions `dims`.
"""
ClimArray(A::AbstractArray, dims::Tuple; refdims=(), name="", attrib=nothing) =
ClimArray(A, DimensionalData.formatdims(A, dims), refdims, name, attrib)
ClimArray(A::AbstractArray, dims::Tuple, name::String; refdims=(), attrib=nothing) =
ClimArray(A, DimensionalData.formatdims(A, dims), refdims, name, attrib)

Base.parent(A::ClimArray) = A.data
Base.@propagate_inbounds Base.setindex!(A::ClimArray, x, I::Vararg{DimensionalData.StandardIndices}) =
    setindex!(A.data, x, I...)

DimensionalData.metadata(A::ClimArray) = A.attrib
DimensionalData.rebuild(A::ClimArray, data::Any, dims::Tuple=dims(A), refdims=DimensionalData.refdims(A),
name="", attrib=nothing) = ClimArray(data, dims, refdims, name, attrib)
DimensionalData.basetypeof(::ClimArray) = ClimArray

# The following basic methods allow indexing with tuples, (Time(5), Lon(3))
Base.getindex(A::ClimArray, i::Tuple) = A[i...]
Base.setindex!(A::ClimArray, x, i::Tuple) = setindex!(A, x, i...)
Base.view(A::ClimArray, i::Tuple) = view(A, i...)

# Remove reference dims from printing, and show attributes if any
function Base.show(io::IO, A::ClimArray)
    l = nameof(typeof(A))
    printstyled(io, nameof(typeof(A)); color=:blue)
    if A.name != ""
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
    print(io, "and")
    printstyled(io, " data: "; color=:green)
    dataA = data(A)
    print(io, summary(dataA), "\n")
    DimensionalData.custom_show(io, data(A))
end
