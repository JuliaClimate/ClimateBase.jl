#=
Continuation of data, to e.g. fill missing values
=#
using SignalDecomposition
export sinusoidal_continuation

_numbertype(T) = _numbertype(eltype(T))
_numbertype(::Type{Union{T, Missing}}) where T = T
_numbertype(::Type{T}) where {T <: Real} = T

# TODO: perhaps it is worth considering a version that makes true time in seconds
# instead of months, which in the future would allow support of diurnal cycles
# TODO: YES! Make true time be seconds!

"""
    sinusoidal_continuation(T::ClimArray, fs = [1, 2]; Tmin = -Inf, Tmax = Inf)
Fill-in the missing values of spatiotemporal field `T`, by fitting sinusoidals
to the non-missing values, and then using the fitted sinusoidals for the missing values.

Frequencies are given per year (frequency 2 means 1/2 of a year).

`Tmin, Tmax` limits are used to clamp the result into this range (no clamping in the
default case).
"""
function sinusoidal_continuation(T, frequencies = [1.0, 2.0]; Tmin = -Inf, Tmax = Inf)
    E = _numbertype(T)
    lpv = Sinusoidal(E.(frequencies ./ DAYS_IN_ORBIT))
    fullT = copy(T)
    # TODO: this must be extended to a general "true time" function
    truetime = time_in_days(dims(T, Time).val, E)
    for i in otheridxs(T, Time)
        x = T[i...]
        any(ismissing, x) || continue # this timeseries needs no correction
        mi = findall(!ismissing, x)
        x0 = Array{Float32}(x[mi])
        t0 = truetime[mi]
        sea, res = decompose(t0, x0, lpv)
        c = cosines(truetime, lpv)
        imi = findall(ismissing, x)
        x[imi] .= c[imi]
        fullT[i...] .= x
    end
    if !isinf(Tmin) && !isinf(Tmin)
        fullT = clamp.(fullT, Tmin, Tmax)
    end
    return nomissing(fullT)
end

function cosines(t, lpv::Sinusoidal)
    s = fill(lpv.A[1], length(t))
    for i in 2:length(lpv.A)
        @. s += lpv.A[i] * cos(2π * lpv.fs[i-1] * t + lpv.φ[i])
    end
    return s
end
