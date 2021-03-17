#=
timeseries decompositions, i.e. separating some data into components
=#
export seasonal_decomposition
"""
    seasonal_decomposition(A::ClimArray, fs = [1, 2]) → seasonal, residual
Decompose `A` into a seasonal and residual components, where the seasonal contains the
periodic parts of `A`, with frequencies given in `fs`, and residual contains what's left.

`fs` is measured in 1/year. This function works even for non-equispaced time axis (e.g.
monthly averages) and uses LPVSpectral.jl and SignalDecomposition.jl.
"""
function seasonal_decomposition(A::AbDimArray, fs::Vector)
    @assert hasdim(A, Time)
    E = _numbertype(eltype(A))
    method = Sinusoidal(E.(fs ./ DAYS_IN_ORBIT))
    seasonal = DimensionalData.basetypeof(A)(copy(Array(A)), dims(A))
    residual = DimensionalData.basetypeof(A)(copy(Array(A)), dims(A))

    t = dims(A, Time).val
    truetime = time_in_days(t, E)
    for i in otheridxs(A, Time)
        y = Array(A[i...])
        sea, res = SignalDecomposition.decompose(truetime, y, method)
        seasonal[i...] .= sea
        residual[i...] .= res
    end
    return seasonal, residual
end

function seasonal_decomposition(A::AbDimArray{T, 1}, fs::Vector) where {T}
    @assert hasdim(A, Time)
    E = _numbertype(eltype(A))
    method = Sinusoidal(E.(fs ./ DAYS_IN_ORBIT))
    seasonal = DimensionalData.basetypeof(A)(copy(Array(A)), dims(A))
    residual = DimensionalData.basetypeof(A)(copy(Array(A)), dims(A))

    t = dims(A, Time).val
    truetime = time_in_days(t, E)
    y = Array(A)
    sea, res = SignalDecomposition.decompose(truetime, y, method)
    seasonal .= sea
    residual .= res
    return seasonal, residual
end
