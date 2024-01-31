module ClimateBase
export NCDataset

include("core/coredefs.jl")
include("core/prettyprint.jl")
include("core/aggregation.jl")

include("physical_dimensions/spatial.jl")
include("physical_dimensions/spatial_equalarea.jl")
include("physical_dimensions/temporal.jl")

include("io/vector2range.jl")
include("io/netcdf.jl")

include("interpolations/height_interpolation.jl")
include("interpolations/quantile_space.jl")

# All following will be moved to ClimateTools.jl once its updated
include("climate/solar.jl")
include("tsa/continuation.jl")
include("tsa/decomposition.jl")

include("exports.jl")
include("plotting/geomakie.jl")

end # module
