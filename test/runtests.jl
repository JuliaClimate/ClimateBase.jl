using ClimateBase, Test, Dates
using Statistics

Time = ClimateBase.Time

# TODO: Further test spatial averaging by making one hemisphere 1 and other 0

# Create the artificial dimensional array A that will be used in tests
function monthly_insolation(t::TimeType, args...)
    d = ClimateBase.monthspan(t)
    mean(insolation(τ, args...) for τ in d)
end

lats = -86:4:86
lons = collect(0.5:10:360)
t = Date(2000, 3, 15):Month(1):Date(2020, 3, 15)
tdense = Date(2000, 3, 15):Day(1):Date(2020, 3, 15)

d = (Lon(lons), Lat(lats), Time(t))
A = zeros([length(x) for x in (lons, lats, t)]...)
B = copy(A)

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

A = ClimArray(A, d; name = "lon-variation")
B = ClimArray(B, d; name = "lon-constant")

# %%

@testset "Dropping dimensions" begin
    dt = timemean(A)
    @test dt isa ClimArray
    @test !hasdim(dt, Time)
    @test hasdim(dt, Lon)
    @test hasdim(dt, Lat)
    dl = zonalmean(A)
    @test dl isa ClimArray
    @test !hasdim(dl, Lon)
    @test hasdim(dl, Time)
    @test hasdim(dl, Lat)
    da = latmean(A)
    @test da isa ClimArray
    @test !hasdim(da, Lat)
    @test hasdim(da, Lon)
    @test hasdim(da, Time)
    ds = spacemean(A)
    @test ds isa ClimArray
    @test !hasdim(ds, Lon)
    @test !hasdim(ds, Lat)
    @test hasdim(ds, Time)
end

@testset "Default physical weights" begin
    x, y = hemispheric_means(B)
    @test timemean(x) != dropagg(mean, x, Time)
    d1 = timemean(x) - timemean(y)
    @test abs(d1) < 0.01
    d2 = dropagg(mean, x, Time) - dropagg(mean, y, Time)
    @test d2 < -0.1
    x = B[:, :, 1]
    @test zonalmean(latmean(x)) ≈ latmean(zonalmean(x)) ≈ spacemean(x)
    @test spacemean(x) - mean(x) > 50
    @test all(y -> abs(y) < 1e-8, spacemean(A))
end

@testset "Temporal weighting" begin
    # spatiotemporal weights:
    W = zero(A)
    W[Time(5)] .= 1.0
    res = timeagg(mean, A, W)
    @test all(res .≈ A[Time(5)])
    @test dims(A, Lon) == dims(res, Lon)
    res = timeagg(std, A, W)
    @test all(x -> isapprox(x, 0; atol = 1e-8), res)

    # just time weight
    w = zeros(length(t))
    w[5] = 1.0
    res = timeagg(mean, A, w)
    @test all(res .≈ A[Time(5)])
    @test dims(A, Lon) == dims(res, Lon)
end

@testset "Advance temporal manipulation" begin
    tdense = Date(2000, 3, 1):Day(1):Date(2020, 3, 31)
    tyearly = Date(2000, 3, 1):Year(1):Date(2020, 3, 31)
    mdates = unique!([(year(d), month(d)) for d in tdense])
    ydates = unique!([year(d) for d in tdense])
    tranges = temporalrange(tdense, Dates.month)
    yranges = temporalrange(tdense, Dates.year)
    @test length(tranges) == length(mdates)
    @test length(yranges) == length(ydates)

    @test temporal_sampling(t) == :monthly
    @test temporal_sampling(tdense) == :daily
    @test temporal_sampling(tyearly) == :yearly

    Bz = zonalmean(B)
    # Generate an array that has daily insolation
    C = ClimArray(zeros(length(lats), length(tdense)), (Lat(lats), Time(tdense)))

    # First version: just count number of days
    for j in 1:length(tdense)
        for i in 1:length(lats)
            C[i, j] = daysinmonth(tdense[j])
        end
    end
    Cm = monthlyagg(C)

    @test length(unique(Array(Cm))) == 4 # there are four unique number of days
    # test that each value when rounded to an integer is an integer (for first slice only
    # all remaining slices are the same)
    for e in Cm[:, 1]
        @test round(Int, e) == e
    end
    @test step(dims(Cm, Time).val) == Month(1)
    @test temporal_sampling(dims(Cm, Time).val) == :monthly

    for j in 1:length(tdense)
        for i in 1:length(lats)
            C[i, j] = daysinyear(tdense[j])
        end
    end
    Cy = yearlyagg(C)
    @test length(unique(Array(Cy))) == 2
    for e in Cy
        @test round(Int, e) == e
    end
    @test step(dims(Cy, Time).val) == Year(1)
    @test temporal_sampling(dims(Cy, Time).val) == :yearly


    # Second version: test actual physics
    for j in 1:length(tdense)
        for i in 1:length(lats)
            C[i, j] = insolation(tdense[j], lats[i])
        end
    end
    Cm = monthlyagg(C)
    @test all(Cm .≈ Bz)
end

@testset "Spatial weighting" begin
    W = zero(A)
    W[Lon(5)] .= 1.0
    res = spaceagg(mean, A, W)
    @test all(res .≈ latmean(A)[Lon(5)])
    @test dims(A, Time) == dims(res, Time)
    res = spaceagg(std, A, W)
    @test all(x -> x>0, res)

    W = zeros(length(lons), length(lats))
    W[5, :] .= 1.0
    W = ClimArray(W, dims(A, (Lon, Lat)))
    res = spaceagg(mean, A, W)
    @test all(res .≈ latmean(A)[Lon(5)])
    @test dims(A, Time) == dims(res, Time)
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
