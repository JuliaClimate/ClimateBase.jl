using ClimateBase, Test, Dates
using Statistics
using ClimateBase.Interpolations
import StatsBase
Time = ClimateBase.Time
cd(@__DIR__)

# Create the artificial dimensional arrays
function monthly_insolation(t::TimeType, args...)
    d = ClimateBase.monthspan(t)
    mean(insolation(τ, args...) for τ in d)
end

lats = -86:4:86
lons = collect(0.5:10:360)
t = Date(2000, 3, 15):Month(1):Date(2020, 2, 15)
tdaily = Date(2000, 3, 15):Day(1):Date(2020, 3, 14)
thourly = DateTime(2000, 3, 15):Hour(1):DateTime(2000, 6, 14)

d = (Lon(lons), Lat(lats), Time(t))
A = zeros([length(x) for x in (lons, lats, t)]...)
B = copy(A)

# generate solar rad
for i in 1:length(lats)
    θ = lats[i]
    s = monthly_insolation.(t, θ)
    for j in 1:length(lons)
        c = cosd(lons[j])
        A[j, i, :] .= c.*s
        B[j, i, :] .= s
    end
end

A = ClimArray(A, d; name = "insolation")
B = ClimArray(B, d; attrib = Dict("a" => 2)) # on purpose without name

# Create similar arrays with Coord dimension
reduced_points = [round(Int, length(lons)*cosd(θ)) for θ in lats]
coords = ClimateBase.reduced_grid_to_points(lats, reduced_points)
sort!(coords; by = reverse) # critical!

coord_dim = Coord(coords, (Lon, Lat))
C = zeros(length(coords), length(t))
for (i, (lon, θ)) in enumerate(coords)
    C[i, :] .= monthly_insolation.(t, θ) .* cosd(lon)
end

C = ClimArray(C, (coord_dim, Time(t)); name = "has_coords", attrib = Dict("a" => 2))

# %%
@testset "ClimateBase.jl tests" begin
include("general_tests.jl")
include("temporal_tests.jl")
include("space_tests.jl")
include("space_coord_tests.jl")
include("io_tests.jl")
include("advanced_tests.jl")
end