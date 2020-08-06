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
seasonal_decomposition(A::AbDimArray, b) = seasonal_decomposition(dims(A, Time), A, b)

function seasonal_decomposition(t, A::AbDimArray, fs::Vector)
    @assert hasdim(A, Time)
    E = _numbertype(T)
    method = Sinusoidal(E.(fs ./ DAYS_IN_YEAR))
    seasonal = DimensionalData.basetypeof(A)(copy(Array(A)), dims(A); name = A.name*"seasonal")
    residual = DimensionalData.basetypeof(A)(copy(Array(A)), dims(A); name = A.name*"residual")

    # TODO: Fix this to use "real time"
    truetime = Float32.(cumsum(daysinmonth.(t)))
    for i in alongdimidxs(T, Time)
        y = Array(A[i...])
        sea, res = SignalDecomposition.decompose(truetime, y, method)
        seasonal[i...] .= sea
        residual[i...] .= res
    end
    return seasonal, residual
end
