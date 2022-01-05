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

# %% Pretty printing
@testset "pretty printing" begin
    P = sprint(show, MIME"text/plain"(), A)
    @test contains(P, "ClimArray")
    @test contains(P, "Lon")
    @test contains(P, "data")
    @test contains(P, "44")
end

# %% General Statistics Tests
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
    tyearly = Date(2000, 3, 1):Year(1):Date(2020, 3, 31)
    thourly = DateTime(2000, 3, 1):Hour(1):DateTime(2001, 4, 15)
    mdates = unique!([(year(d), month(d)) for d in tdaily])
    ydates = unique!([year(d) for d in tdaily])
    tranges = temporalrange(tdaily, Dates.month)
    yranges = temporalrange(tdaily, Dates.year)
    @testset "time sampling" begin
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
    end

    @testset "sametimespan" begin
        tm1 = Date(2000, 3, 15):Month(1):Date(2010, 3, 15)
        tm2 = Date(2001, 1, 1):Month(1):Date(2011, 1, 1)
        A1 = ClimArray(rand(length(tm1)), (Time(tm1),))
        A2 = ClimArray(rand(length(tm2)), (Time(tm2),))
        B1, B2 = sametimespan(A1, A2)
        @test size(B1) == size(B2)
        @test dims(B1, Time)[1] == Date(2001, 1, 15)
        @test dims(B1, Time)[end] == Date(2010, 3, 15)
        @test dims(B2, Time)[end] == Date(2010, 3, 1)
    end
    
    @testset "monthlyagg and co." begin
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
        @test size(Cm, Time) == size(B, Time) + 1
        # For times besides the first and last month, they sould be the same:
        for k in 5:length(t)-5
            @test all(Cm[Time(k)] .≈ Bz[Time(k)])
        end

        Asea = seasonalyagg(A)
        tsea = dims(Asea, Time).val
        @test Base.step(tsea) == Month(3)
    end

    @testset "interannual variability" begin
        C = spacemean(B)
        for y in unique(year.(t))[1:end-1]
            C[Time(At(Date(y, 3, 15)))] = 0
        end
        dates, vals = seasonality(C)
        @test length(dates) == 12
        @test all(v -> length(v) == 20, vals)
        @test month(dates[3]) == 3
        @test all(iszero, vals[3])
        @test all(!iszero, vals[2])
    end
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

@testset "tropics/extratropics" begin
    tropics, extratropics = tropics_extratropics(B)
    tlats = dims(tropics, Lat).val
    elats = dims(extratropics, Lat).val
    @test all(l -> -30 ≤ l ≤ 30, tlats)
    @test all(l -> abs(l) ≥ 30, elats)
    @test hasdim(tropics, Ti)

    tmean = spacemean(timemean(tropics))
    emean = spacemean(timemean(extratropics))
    amean = spacemean(timemean(B))
    @test amean ≈ (emean + tmean)/2 rtol = 5

    # Test for keyword options for latitude bounds
    tropics, extratropics = tropics_extratropics(B, lower_lat=20, higher_lat=78)
    tlats = dims(tropics, Lat).val
    elats = dims(extratropics, Lat).val
    @test all(l -> -20 ≤ l ≤ 20, tlats)
    @test all(l -> abs(l) ≥ 20, elats)
    @test maximum(abs.(elats)) <= 78
    @test hasdim(tropics, Ti) # probably redundant, already tested

end

# %% IO tests
@testset "NetCDF file IO" begin
    globat = Dict("history" => "test")
    ncwrite("test.nc", (A, B); globalattr = globat)
    Aloaded = ncread("test.nc", "insolation")
    Bloaded = ncread("test.nc", "x2")

    @test A.data == Aloaded.data
    @test DimensionalData.metadata(dims(Aloaded, Lon))["units"] == "degrees_east"
    @test B.data == Bloaded.data
    @test string(Bloaded.name) == "x2"
    @test DimensionalData.metadata(dims(Bloaded, Time))["standard_name"] == "time"

    Bpartly = ncread("test.nc", "x2", ([1,4,5], :, 1:3))
    @test Bpartly.data[:, :, :] == Bloaded.data[[1,4,5], :, 1:3]

    ds = NCDataset("test.nc")
    @test ds.attrib["history"] == "test"
    close(ds)
    rm("test.nc")
end

@testset "Missings handling" begin
    m = rand(Float64, length.((lons, lats)))
    midx = CartesianIndices(m)[1:2, :]
    L = length(midx)
    m[midx] .= -99.9
    M = ClimArray(m, (Lon(lons), Lat(lats)); attrib = Dict("_FillValue" => -99.9), name = "M")
    ncwrite("missing_test.nc", M)
    M2 = ncread("missing_test.nc", "M")
    for i in midx; @test ismissing(M2[i]); end
    # now test `missing_weights`
    B, W = missing_weights(M2)
    @test all(iszero, W.data[midx])
    @test all(isequal(-99.9), B.data[midx])
    mmean = mean(skipmissing(M2.data))
    bmean = mean(B, StatsBase.weights(W))
    @test bmean == mmean
    # test `missing_weights` application to zonal mean
    C = B[3:end, :]
    z1 = zonalmean(C)
    z2 = zonalmean(B, W)
    @test all(z1.data .≈ z2.data)
    rm("missing_test.nc")
end

@testset "Vertical interpolation" begin
    D = ClimArray([1.0:1.0:11.0 2.0:1.0:12.0 3.0:1.0:13.0], (Hei(0.0:2000.0:20000.0), Ti(1:3)))
    pressure_levels = [950.0,850.0,650.0,350.0,250.0,150.0] .* 100.0
    D_pre = interpolate_height2pressure(D, pressure_levels,extrapolation_bc=NaN)
    D_back = interpolate_pressure2height(D_pre, Vector(dims(D,Hei).val),extrapolation_bc=Line())
    my_dim = Dim{:My_Dim}
    E = ClimArray([1.0:1.0:10.0 2.0:1.0:11.0 3.0:2.0:21.0], (my_dim(1:10), Ti(1:3)))
    pressure = ClimArray([1000.0:-100.0:100.0 1000.0:-100.0:100.0 1000.0:-100.0:100.0] *100.0, (my_dim(1:10), Ti(1:3)))
    E_pre = interpolation2pressure(E, pressure, pressure_levels; vertical_coord=my_dim, extrapolation_bc=NaN )
    E_pre2 = interpolation2pressure(reverse(E,dims=my_dim), reverse(pressure,dims=my_dim), pressure_levels; vertical_coord=my_dim, extrapolation_bc=NaN, descending = false )

    @test hasdim(D_pre,Pre())
    @test gnv(dims(D_back,Hei)) == gnv(dims(D,Hei)) && gnv(dims(D_back,Ti)) == gnv(dims(D,Ti))
    @test E_pre.data == E_pre2.data
    @test E_pre[Pre(1)] < E_pre[Pre(2)] < E_pre[Pre(3)] < E_pre[Pre(4)] < E_pre[Pre(5)]
    @test D_pre[Pre(1)] < D_pre[Pre(2)] < D_pre[Pre(3)] < D_pre[Pre(4)] < D_pre[Pre(5)]
    @test D_back[Hei(1)] < D_back[Hei(2)] < D_back[Hei(3)] < D_back[Hei(4)] < D_back[Hei(5)]
    @test D[Hei(1)] < D_pre[Pre(3)] < D[Hei(11)]
end
