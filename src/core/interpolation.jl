T0 = 288.15 
gamma = -6.5e-3
P0 = 101325
const gravity = 9.8076
const Rd = 287.053

using Interpolations

"""

interpolation2pressure(A::ClimArray, pressure::ClimArray, pressure_levels::Vector; heightname=Hei(), extrapolation_bc=NaN )
Return an ClimArray where the vertical coordinate is pressure.

A: ClimArray with some arbitrary height corrdinate (height above ground, model levels, etc...) with the variable that should be interpolated
pressure: ClimArray with pressure values in Pascal and the same coordinates as A.
pressure_levels: Vector which contains the pressure levels in Pascal that A should be interpolated on
extrapolation_bc: extrapolation is set to NaN by default, but can be e.g. linear with extrapolation_bc = Line()

"""

function interpolation2pressure(A::ClimArray, pressure::ClimArray, pressure_levels::Vector; heightname=Hei(), extrapolation_bc=NaN )

    dims(A) == dims(pressure) || error("Pressure and Variable Array do not have the same dimensions")
    
    # The interplation requires ascending coordinates:
    if issorted(dims(A,Hei)) ==  false
        A = reverse(A,dims=Hei())
        pressure = reverse(pressure,dims=Hei())
    end
        
    # construct output Array:
    pre = Pre(pressure_levels; metadata = Dict())
    dims_no_height = otherdims(A, heightname)
    out_dims = (dims_no_height...,pre)
    dimension_lengths = length.(out_dims)
    
    int_array = ClimArray(zeros(eltype(A), Tuple(dimension_lengths)), out_dims ; name = A.name, attrib = A.attrib)    
    for i in otheridxs(A, heightname)
        #print(pressure[i].data)
        
        itp = LinearInterpolation(pressure[i],A[i],extrapolation_bc=extrapolation_bc)
        int_array[i] = itp(pressure_levels)
            
    end
    
    return int_array
end



function pressure2height(pressure)
    return (T0/gamma) * ( (pressure/P0)^( - gamma*Rd/gravity ) - 1 )
end


function height2pressure(height)
    return P0 * exp(- gravity/(Rd * gamma ) * log(1 + (height*gamma)/T0 ))
end

"""

interpolate_height2pressure(A::ClimArray, pressure_levels::Vector; extrapolation_bc=NaN )
Return an ClimArray where the vertical coordinate is pressure.

A: ClimArray with height above ground [m] as height coordinate
pressure_levels: Vector which contains the pressure levels in Pascal that A should be interpolated on

extrapolation_bc: extrapolation is set to NaN by default, but can be e.g. linear with extrapolation_bc = Line()

"""



function interpolate_height2pressure(A::ClimArray,pressure_levels::Vector; extrapolation_bc=NaN)
    
    # check whether height is in A?
    
    # The interplation requires ascending coordinates:
    # Pressure has to be ascending, therefore height has to be descending
    if issorted(dims(A,Hei)) ==  true
        A = reverse(A,dims=Hei())
    else
        error("Height is not an ascending coordinate.")
    end
    
    print(A)
    
    pressure = height2pressure.(dims(A,Hei()).val)
    print(pressure)
    
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
Return an ClimArray where the vertical coordinate is pressure.

A: ClimArray with height above ground [m] as height coordinate
heights: Vector which contains the heights that A should be interpolated on in meters above ground

extrapolation_bc: extrapolation is set to NaN by default, but can be e.g. linear with extrapolation_bc = Line()

"""


function interpolate_pressure2height(A::ClimArray,heights::Vector; extrapolation_bc=NaN)
    
    # check whether pressure is in A?

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

#example use case:

#example_data = ClimArray([1.:1.:11. 2.:1.:12. 3.:2.:23.], (Hei(0.:2000.:20000.), Ti(1:3)))
#pressure_levels = [1000.,850.,700.,500.,400.,300.,200.,150.,100.] .* 100.
#example_data_pre = interpolate_height2pressure(example_data, pressure_levels,extrapolation_bc=NaN)
#example_data_back = interpolate_pressure2height(my_data_pre, Vector(dims(my_data,Hei).val),extrapolation_bc=Line())


