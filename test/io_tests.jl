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

end # NetCDF tests