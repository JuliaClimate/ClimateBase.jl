module ClimateBase
export NCDataset

include("core/coredefs.jl")
include("core/prettyprint.jl")
include("core/nc_io.jl")
include("core/aggregation.jl")

include("physical_dimensions/spatial.jl")
include("physical_dimensions/spatial_equalarea.jl")
include("physical_dimensions/temporal.jl")

include("climate/solar.jl")
include("climate/albedo.jl")

include("tsa/continuation.jl")
include("tsa/decomposition.jl")

end # module
