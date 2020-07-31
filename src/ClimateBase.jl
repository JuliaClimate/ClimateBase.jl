module ClimateBase

# TODO: Be sure all exported names have docstrings
include("aggregation.jl")
include("loading_nc.jl")
include("continuation.jl")

include("physical_dimensions/spatial.jl")
include("physical_dimensions/temporal.jl")

include("climate/solar.jl")
include("climate/albedo.jl")

end # module
