# %% Time Tests
@testset "Temporal Tests" begin 

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



end # tempooral tests