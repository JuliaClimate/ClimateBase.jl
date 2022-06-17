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

@testset "quantile space" begin

end

@testset "value space" begin
    t = 0:0.01:5π
    x = cos.(t)
    y = sin.(t)
    rx = range(-1.0, nextfloat(1.0); length = 21)
    @testset "1D" begin
        using StatsBase
        ymeans, bin_indices = value_space(x, y; Arange = rx)
        weights = length.(bin_indices)
        @test sum(weights) == length(y)
        # trigonometric functions have conentrated value weight at the edges
        @test weights[1] > weights[length(weights)÷2]
        # sine must be very small when cosine is large
        @test ymeans[1] < weights[length(weights)÷2]
        @test mean(ymeans, Weights(weights)) ≈ StatsBase.mean(y)
    end
    @testset "2D" begin
        z = y
        zmeans, bin_indices = value_space(x, y, z; Arange = rx, Brange = rx)
        # Because cozine and size cannever be 1 at the same time, all corners
        # of z must be nan:
        L = length(rx)-1
        @test all(isnan, [zmeans[1,1], zmeans[1,L], zmeans[L,1], zmeans[L,L]])
        # and again because of where cosines and sines may have values at the same time
        # zmeans only is non NaN at specific locations one can derive (depends on step size)
        non_nan_i = 6:15
        @test all(!isnan, zmeans[1, non_nan_i])
        @test all(!isnan, zmeans[non_nan_i, end])
        # and because z is y, it must increase along the second dimension
        @test issorted(zmeans[1, non_nan_i])
    end
end
