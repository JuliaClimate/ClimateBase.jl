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
