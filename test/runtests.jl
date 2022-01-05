using ClimateBase, Test, Dates
using Statistics
using ClimateBase.Interpolations
import StatsBase
Time = ClimateBase.Time
cd(@__DIR__)

# Create the artificial dimensional array A that will be used in tests
function monthly_insolation(t::TimeType, args...)
    d = ClimateBase.monthspan(t)
    mean(insolation(τ, args...) for τ in d)
end

lats = -86:4:86
lons = collect(0.5:10:360)
t = Date(2000, 3, 15):Month(1):Date(2020, 2, 15)
tdaily = Date(2000, 3, 15):Day(1):Date(2020, 3, 14)

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

# %% 
include("general_tests.jl")
include("temporal_tests.jl")
include("space_tests.jl")
include("space_coord_tests.jl")
include("io_tests.jl")
include("advanced_tests.jl")