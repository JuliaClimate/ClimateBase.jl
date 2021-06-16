const T0 = 288.15
const gamma = -6.5e-3
const P0 = 101325
const gravity = 9.8076
const Rd = 287.053

using Interpolations
export interpolation2pressure, interpolate_height2pressure, interpolate_pressure2height, pressure2height, height2pressure, hello, Line

"""
interpolation2pressure(A::ClimArray, pressure::ClimArray, pressure_levels::Vector; heightname=Hei(), extrapolation_bc=NaN )
Return an ClimArray where the vertical coordinate is pressure. Pressure values need to ascending or descending

`A`: `ClimArray` with some arbitrary height coordinate (height above ground, model levels, etc...) with the variable that should be interpolated
`pressure`: `ClimArray` with pressure values in Pascal and the same coordinates as `A`.
`pressure_levels`: Vector which contains the pressure levels in Pascal that A should be interpolated on
`vertical_coord` : name of the vertical coordinate. Should usually be `Hei` but can be model levels or something else that does not represent physical height
`extrapolation_bc`: extrapolation is set to `NaN` by default, but can be any value or linear extrapolation with `extrapolation_bc = Line()` For other extrapolation methods use the `Interpolations` package.
`order` : specifies whether pressure is ordered descendingly or ascendingly along the vertical coordinate
"""
function interpolation2pressure(A::ClimArray, pressure::ClimArray, pressure_levels::Vector; vertical_coord=Hei(), extrapolation_bc=NaN,descending=true )

    dims(A) == dims(pressure) || error("Pressure and Variable Array do not have the same dimensions")
    hasdim(A, vertical_coord) || error("Vertical dimension " + string(vertical_coord) +  " not found in array")

    # The interplation requires ascending pressure values:
    if descending == true
        A = reverse(A,dims=vertical_coord)
        pressure = reverse(pressure,dims=vertical_coord)
    end

    # construct output Array:
    pre = Pre(pressure_levels; metadata = Dict())
    dims_no_height = otherdims(A, vertical_coord)
    out_dims = (dims_no_height...,pre)
    dimension_lengths = length.(out_dims)

    int_array = ClimArray(zeros(eltype(A), Tuple(dimension_lengths)), out_dims ; name = A.name, attrib = A.attrib)
    for i in otheridxs(A, vertical_coord)
        #print(pressure[i].data)

        itp = LinearInterpolation(pressure[i],A[i],extrapolation_bc=extrapolation_bc)
        int_array[i] = itp(pressure_levels)

    end

    return int_array
end

"""
pressure2height(height)
Converts pressure [Pa] to height above ground [m].
The calculation is done under the assumption of hydrostatic equilibrium with a constant lapse rate of -6.5 K/km, starting from a surface temperature T0=288.15 K.
See also equations 39 and 40 in: Berberan-Santos, M. N., Bodunov, E. N., & Pogliani, L. (1997). On the barometric formula. American Journal of Physics, 65(5), 404-412.
"""
function pressure2height(pressure)
    return (T0/gamma) * ( (pressure/P0)^( - gamma*Rd/gravity ) - 1.0 )
end

"""
height2pressure(height)
Converts height above ground [m] to pressure [Pa].
The calculation is done under the assumption of hydrostatic equilibrium with a constant lapse rate of -6.5 K/km, starting from a surface temperature T0=288.15 K.
See also equations 39 and 40 in: Berberan-Santos, M. N., Bodunov, E. N., & Pogliani, L. (1997). On the barometric formula. American Journal of Physics, 65(5), 404-412.
"""
function height2pressure(height)
    return P0 * exp(- gravity/(Rd * gamma ) * log(1.0 + (height*gamma)/T0 ))
end


"""
interpolate_height2pressure(A::ClimArray, pressure_levels::Vector; extrapolation_bc=NaN )
Return a `ClimArray` where the vertical coordinate is pressure. Pressure levels are calculated with the `height2pressure` function, based on hydrostatic equilibrium, and then interpolated to the given levels.

`A`: `ClimArray` with height above ground [m] as height coordinate
`pressure_levels`: Vector which contains the pressure levels in Pascal that `A` should be interpolated on
`extrapolation_bc`: extrapolation is set to `NaN` by default, but can be any value or linear extrapolation with `extrapolation_bc = Line()` For other extrapolation methods use the `Interpolations` package.
"""
function interpolate_height2pressure(A::ClimArray,pressure_levels::Vector; extrapolation_bc=NaN)

    hasdim(A, Hei) || error("Hei() dimension not found in array")

    # The interplation requires ascending coordinates:
    # Pressure has to be ascending, therefore height has to be descending
    if issorted(dims(A,Hei)) ==  true
        A = reverse(A,dims=Hei())
    else
        error("Height is not an ascending coordinate. Try reversing with reverse().")
    end

    pressure = height2pressure.(dims(A,Hei()).val)

    # construct output Array:
    pre = Pre(pressure_levels; metadata = Dict())
    dims_no_height = otherdims(A, Hei())
    out_dims = (dims_no_height...,pre)
    dimension_lengths = length.(out_dims)

    int_array = ClimArray(zeros(eltype(A), Tuple(dimension_lengths)), out_dims ; name = A.name, attrib = A.attrib)

    for i in otheridxs(A, Hei())

        itp = LinearInterpolation(pressure,A[i],extrapolation_bc=extrapolation_bc)
        int_array[i] = itp(pressure_levels)

    end
    return int_array

end


"""
interpolate_pressure2height(A::ClimArray,heights::Vector; extrapolation_bc=NaN)
Return a `ClimArray` where the vertical coordinate is height above ground.  Height above ground is calculated with the `pressure2height` function, based on hydrostatic equilibrium, and then interpolated to the given heights.

`A`: ClimArray with pressure [Pa] as height coordinate
`heights`: Vector which contains the heights that `A` should be interpolated on in meters above ground
`extrapolation_bc`: extrapolation is set to `NaN` by default, but can be any value or linear extrapolation with `extrapolation_bc = Line()` For other extrapolation methods use the `Interpolations` package.
"""
function interpolate_pressure2height(A::ClimArray,heights::Vector; extrapolation_bc=NaN)

    hasdim(A, Pre) || error("Pre() dimension not found in array")

    height = pressure2height.(dims(A,Pre()).val)

    # construct output Array:
    hei = Hei(heights; metadata = Dict())
    dims_no_pre = otherdims(A, Pre())
    out_dims = (dims_no_pre...,hei)
    dimension_lengths = length.(out_dims)

    int_array = ClimArray(zeros(eltype(A), Tuple(dimension_lengths)), out_dims ; name = A.name, attrib = A.attrib)

    for i in otheridxs(A, Pre())

        itp = LinearInterpolation(height,A[i],extrapolation_bc=extrapolation_bc)
        int_array[i] = itp(heights)

    end
    return int_array

end
