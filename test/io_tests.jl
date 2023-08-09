@testset "NetCDF IO" begin

@testset "vector2range" begin
    th = collect(thourly)
    th2 = ClimateBase.vector2range(th)
    @test th2 isa AbstractRange
    @test th2 == th
end

@testset "NetCDF standard tests" begin
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

@testset "NetCDF Coord tests" begin
    ncwrite("test.nc", C)
    n = string(DimensionalData.name(C))
    Cloaded = ncread("test.nc", "has_coords")
    @test size(Cloaded) == size(C)
    @test hasdim(Cloaded, Coord)
    @test gnv(dims(Cloaded, Coord)) == gnv(dims(C, Coord))
    @test gnv(Cloaded) == gnv(C)

    # TODO: Here we need a lot of tests for all the super weird different
    # ways that there are to save a Coord datafile... But we can't creat them.
    # So we need to upload files somewhere and load them here.
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
    @test bmean ≈ mmean # for some reason here we lose accuracy: 0.4992455417474702 vs 0.49924554174747
    # test `missing_weights` application to zonal mean
    C = B[3:end, :]
    z1 = zonalmean(C)
    z2 = zonalmean(B, W)
    @test all(z1.data .≈ z2.data)
    rm("missing_test.nc")
end

@testset "CFTime dates" begin
    using NCDatasets.CFTime: DateTime360Day
    cfdates = collect(DateTime360Day(1900,01,01):Day(1):DateTime360Day(1919,12,30))
    x = float.(month.(cfdates))
    X = ClimArray(x, (Tim(cfdates),); name = "x")

    @testset "temporal stats" begin
        @test temporal_sampling(cfdates) == :daily
        @test monthday_indices(cfdates) == 1:360:length(cfdates)
        trange = temporalranges(cfdates)
        for i in 1:length(trange)
            @test trange[i] == (1 + (i-1)*30):(i*30)
        end


        Y = monthlyagg(X)
        @test length(Y) == 20*12
        ty = gnv(dims(Y, Tim))
        @test temporal_sampling(ty) == :monthly
        @test step(ty) == Month(1)
        for (i, y) in enumerate(Y)
            @test y == mod1(i, 12)
        end
        Z = yearlyagg(X)
        @test length(Z) == 20
        # The mean of 1 to 12 is by definition 6.5
        @test all(isequal(6.5), Z)
        tz = gnv(dims(Z, Tim))
        @test temporal_sampling(tz) == :yearly
        @test step(tz) == Year(1)
    end
    @testset "Writing/Reading CFTime" begin
        ncwrite("cftime_test.nc", X)
        @test isfile("cftime_test.nc")
        X2 = ncread("cftime_test.nc", "x")
        @test eltype(dims(X2, Tim)) == DateTime360Day
        @test gnv(dims(X2, Tim)) == cfdates
    end

end

end # NetCDF tests