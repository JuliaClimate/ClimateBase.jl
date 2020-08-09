using ClimateBase, Test, Dates
using Statistics

# TODO: Report to DimemnsionalData.jl that `dropdim` does not retain type.

# Create the artificial dimensional array A that will be used in tests
function monthly_insolation(t::TimeType, args...)
    d = monthspan(t)
    mean(insolation(τ, args...) for τ in d)
end

lats = -90:5:90
lons = 0:10:359
t = Date(2000, 3, 15):Month(1):Date(2020, 3, 15)

d = (Lon(lons), Lat(lats), Time(t))
A = zeros([length(x) for x in (lons, lats, t)]...)
B = zeros([length(x) for x in (lons, lats, t)]...)

# generate solar rad
for i in 1:length(lats)
    θ = lats[i]
    s = [monthly_insolation(ti, θ) for ti in t]
    for j in 1:length(lons)
        c = cosd(lons[j])
        A[j, i, :] .= c.*s
        B[j, i, :] .= s
    end
end

A = ClimArray(A, d; name = "lon variation")
B = ClimArray(B, d; name = "lon constant")

@testset "Dropping dimensions" begin
    dt = timemean(A)
    # @test dt isa ClimArray
    @test !hasdim(dt, Time)
    @test hasdim(dt, Lon)
    dl = zonalmean(A)
    # @test dl isa ClimArray
    @test !hasdim(dl, Lon)
    @test hasdim(dl, Time)
    da = latmean(A)
    # @test da isa ClimArray
    @test !hasdim(da, Lat)
    @test hasdim(da, Lon)
    @test hasdim(da, Time)
    ds = spacemean(A)
    # @test ds isa ClimArray
    @test !hasdim(ds, Lon)
    @test !hasdim(ds, Lat)
    @test hasdim(ds, Time)
end

@testset "Proper physical weights" begin
    x = spacemean(B)
    @test timemean(x) != dropagg(mean, x, Time)
    x = B[:, :, 1]
    @test zonalmean(latmean(x)) ≈ latmean(zonalmean(x)) ≈ spacemean(x)
    @test spacemean(x) - mean(x) > 50
end

@testset "Insolation/Hemispheres" begin
    @test !any(x -> x < 0, B)
    x = timemean(B)
    x, y = hemispheric_functions(zonalmean(x))
    # @test x isa ClimArray
    x[Lat(1)] == y[Lat(1)]
    for j in 2:size(x, Lat)
        @test x[j] < x[j-1]
        @test y[j] < y[j-1]
    end
end
