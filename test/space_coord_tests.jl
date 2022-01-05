@testset "Spatial Coordinates Tests" begin
    @testset "Indexing and Basics" begin
        # Needs https://github.com/rafaqz/DimensionalData.jl/issues/357
        # subsel = C[Lat(0..8)]
        # @test hasdim(subsel, Coord)
        # @test size(subsel, Coord) < size(C, Coord)
        # @test all(0 ≤ c[2] ≤ 8 for c in gnv(dims(subsel, Coord)))

        @test length(collect(spatialidxs(C))) == size(C, Coord)

        Atran = transform_to_coord(A)
        @test hasdim(Atran, Coord)
        @test lons == sort!(unique!([c[1] for c in gnv(dims(Atran, Coord))]))

    end

    @testset "Zonal mean" begin
        zC = zonalmean(C)
        @test hasdim(zC, Lat)
        @test gnv(dims(zC, Lat)) == lats
        @test maximum(C) < 1 # because we aggregate over cosine Lon, the result is small
        # TODO: test with weights
    end

    @testset "Hemispheric stuff" begin
        # TODO: Would be useful to test this hemispheric stuff with coordinates
        # that include 0 latitude and that do not cover symmetric latitdues
        # even though it is super uncommon
        nhi, shi = ClimateBase.hemisphere_indices(gnv(dims(C, Coord)))
        @test length(nhi) == length(shi)
        @test sort!(vcat(nhi, shi)) == 1:size(C, Coord)

        Cnh, Csh = hemispheric_functions(C) # symmetric, non-zero latitudes
        @test hasdim(Cnh, Coord) && hasdim(Csh, Coord)
        @test size(Cnh, Coord) == size(Csh, Coord)

        # hemispheric_Functions makes latitudes the same by convention:
        @test gnv(dims(Cnh, Coord)) == gnv(dims(Csh, Coord))
        @test minimum(uniquelats(Cnh)[2]) ≥ 0
        @test maximum(uniquelats(Csh)[2]) ≤ 0

        Cnh, Csh = hemispheric_means(C)
        @test !hasdim(Cnh, Coord) && !hasdim(Csh, Coord)
        @test hasdim(Cnh, Time)
        # TODO: test actual numeric values
    end


end