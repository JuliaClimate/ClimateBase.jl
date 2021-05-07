using Interpolations

"""

interpolation2pressure(A::ClimArray, Pressure::ClimArray, PressureLevels::Vector; HeightName=Hei() )
Return an ClimArray where the vertical coordinate is pressure.

A: ClimArray with some arbitrary height corrdinate (height above ground, model levels, etc...) with the variable that should be interpolated
Pressure: ClimArray with pressure values and the same coordinates as A.
PressureLevels: Array which contains the pressure levels that A should be interpolated on

"""

function interpolation2pressure(A::ClimArray, Pressure::ClimArray, PressureLevels::Vector; HeightName=Hei() )

    dims(A) == dims(Pressure) || error("Pressure and Variable Array do not have the same dimensions")
    # also check whether the coordinates are the same?!
    
    # construct output Array: There might be an easier way to do this
    pre = Pre(PressureLevels; metadata = Dict())
    DimsNoHeight = otherdims(A, HeightName)
    OutDimensions = (DimsNoHeight...,pre)
    DimensionLengths = []
    for dim in OutDimensions
        push!(DimensionLengths, length(dim))
    end
    
    IntArray = ClimArray(zeros(Tuple(DimensionLengths)), OutDimensions ; name = A.name, attrib = A.attrib)
    
    for i in otheridxs(A, HeightName)
        
        itp = LinearInterpolation(Pressure[i],A[i],extrapolation_bc=NaN) # Extrapolation should be NaNs
        IntArray[i] = itp(PressureLevels)
        
    end
    return IntArray
end

