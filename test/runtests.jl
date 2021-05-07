using ClimateBase, Test, Dates
using Statistics
Time = ClimateBase.Time
cd(@__DIR__)

# Create the artificial dimensional array A that will be used in tests
function monthly_insolation(t::TimeType, args...)
    d = ClimateBase.monthspan(t)
    mean(insolation(τ, args...) for τ in d)
end

lats = -86:4:86
lons = collect(0.5:10:360)
t = Date(2000, 3, 15):Month(1):Date(2020, 3, 15)
tdaily = Date(2000, 3, 15):Day(1):Date(2020, 3, 15)

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

A = ClimArray(A, d; name = "insolation")
B = ClimArray(B, d; attrib = Dict("a" => 2)) # on purpose without name

# %% General Tests
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
    ds = dropagg(mean, A, Time)
    @test dt isa ClimArray
    @test !hasdim(dt, Time)
    @test hasdim(dt, Lon)
    @test hasdim(dt, Lat)
end

@testset "dropagg with weighting" begin
    # spatiotemporal weights:
    W = zero(A)
    W[Time(5)] .= 1.0
    res = dropagg(mean, A, Time, W)
    @test all(res .≈ A[Time(5)])
    @test dims(A, Lon).val == dims(res, Lon).val
    res = dropagg(std, A, Time, W)
    @test all(x -> isapprox(x, 0; atol = 1e-8), res)

    # just time weight
    w = zeros(length(t))
    w[5] = 1.0
    res = dropagg(mean, A, Time, w)
    @test all(res .≈ A[Time(5)])
    @test dims(A, Lon).val == dims(res, Lon).val
end


# %% Time Tests
@testset "Temporal weighting" begin
    # spatiotemporal weights:
    W = zero(A)
    W[Time(5)] .= 1.0
    res = timeagg(mean, A, W)
    @test all(res .≈ A[Time(5)])
    @test dims(A, Lon).val == dims(res, Lon).val
    res = timeagg(std, A, W)
    @test all(x -> isapprox(x, 0; atol = 1e-8), res)

    # just time weight
    w = zeros(length(t))
    w[5] = 1.0
    res = timeagg(mean, A, w)
    @test all(res .≈ A[Time(5)])
    @test dims(A, Lon).val == dims(res, Lon).val

    # TODO: more tests needed here, e.g. for timeagg(mean, t, a, w)
end

@testset "Advanced temporal manipulation" begin
    tdaily = Date(2000, 3, 1):Day(1):Date(2020, 3, 31)
    tyearly = Date(2000, 3, 1):Year(1):Date(2020, 3, 31)
    thourly = DateTime(2000, 3, 1):Hour(1):DateTime(2001, 4, 15)
    mdates = unique!([(year(d), month(d)) for d in tdaily])
    ydates = unique!([year(d) for d in tdaily])
    tranges = temporalrange(tdaily, Dates.month)
    yranges = temporalrange(tdaily, Dates.year)
    @test length(tranges) == length(mdates)
    @test length(yranges) == length(ydates)

    @test temporal_sampling(t) == :monthly
    @test temporal_sampling(tdaily) == :daily
    @test temporal_sampling(tyearly) == :yearly

    # Test hourly stuff
    @test temporal_sampling(thourly) == :hourly
    @test temporal_sampling(collect(thourly)) == :hourly
    hmys = maxyearspan(thourly, :hourly)
    @test hmys < length(thourly)
    @test year(thourly[1]) == year(thourly[hmys])-1

    # Arbitrary random range
    tother = collect(thourly)
    tother[4] = DateTime(2000, 7, 2, 2, 2)
    @test temporal_sampling(tother) == :other

    # test yearly temporal weights (i.e., no weighting)
    X = ClimArray(rand(3,3), (Lon(1:3), Time(tyearly[1:3])))
    W = [0, 1, 0]
    @test vec(timemean(X).data) == vec(mean(X.data; dims = 2))
    @test timemean(X, W).data == X.data[:, 2]

    C = ClimArray(zeros(length(lats), length(tdaily)), (Lat(lats), Time(tdaily)))

    # First version: just count number of days
    for j in 1:length(tdaily)
        for i in 1:length(lats)
            C[i, j] = daysinmonth(tdaily[j])
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

    for j in 1:length(tdaily)
        for i in 1:length(lats)
            C[i, j] = daysinyear(tdaily[j])
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
    for j in 1:length(tdaily)
        for i in 1:length(lats)
            C[i, j] = insolation(tdaily[j], lats[i])
        end
    end
    Bz = zonalmean(B)
    Cm = monthlyagg(C)
    @test all(Cm .≈ Bz)

    Asea = seasonalyagg(A)
    tsea = dims(Asea, Time).val
    @test Base.step(tsea) == Month(3)
end

# %% Space tests
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

# %% IO tests
@testset "NetCDF file IO" begin
    globat = Dict("history" => "test")
    ncwrite("test.nc", (A, B); globalattr = globat)
    Aloaded = ncread("test.nc", "insolation")
    Bloaded = ncread("test.nc", "x2")

    @test A.data == Aloaded.data
    @test dims(Aloaded, Lon).metadata["units"] == "degrees_east"
    @test B.data == Bloaded.data
    @test string(Bloaded.name) == "x2"
    @test dims(Bloaded, Time).metadata["standard_name"] == "time"

    ds = NCDataset("test.nc")
    @test ds.attrib["history"] == "test"
    close(ds)
    rm("test.nc")
end
