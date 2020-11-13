module ClimateBase

# TODO: Be sure all exported names have docstrings
include("core/coredefs.jl")
include("core/nc_io.jl")
include("core/aggregation.jl")

include("physical_dimensions/spatial.jl")
include("physical_dimensions/temporal.jl")

include("climate/solar.jl")
include("climate/albedo.jl")

include("tsa/continuation.jl")
include("tsa/decomposition.jl")

end # module
