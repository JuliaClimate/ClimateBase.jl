@testset "pretty printing" begin
    P = sprint(show, MIME"text/plain"(), A)
    @test contains(P, "ClimArray")
    @test contains(P, "Lon")
    @test contains(P, "data")
    @test contains(P, "44")
end

@testset "General Statistics" begin
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

end # general statistics