@testset "Space Tests" begin

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



end # Space Tests