module ClimateBase
export NCDataset

include("core/coredefs.jl")
include("core/prettyprint.jl")
include("core/nc_io.jl")
include("core/aggregation.jl")
include("core/interpolation.jl")

include("physical_dimensions/spatial.jl")
include("physical_dimensions/spatial_equalarea.jl")
include("physical_dimensions/temporal.jl")

include("climate/solar.jl")

include("tsa/continuation.jl")
include("tsa/decomposition.jl")

include("exports.jl")

using Requires
function __init__()
    @require GeoMakie="db073c08-6b98-4ee5-b6a4-5efafb3259c6" begin
        include("plotting/geomakie.jl")
    end
end

end # module
