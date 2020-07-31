#=
Continuation of data, to e.g. fill missing values
=#
using SignalDecomposition
export sinusoidal_continuation

"""
    sinusoidal_continuation(T::AbDimArray, fs = [1, 2]; Tmin = 0, Tmax = Inf)
Fill-in the missing values of spatiotemporal field `T`, by fitting sinusoidals
to the non-missing values, and then using the fitted sinusoidals for the missing values.

Frequencies are given per year (frequency 2 means 1/2 of a year).

`Tmin, Tmax` limits are used to clamp the result into this range.
"""
function sinusoidal_continuation(T, frequencies = [1.0, 2.0]; Tmin = 0, Tmax = Inf)
    lpv = Sinusoidal(Float32.(frequencies ./ DAYS_IN_YEAR))
    fullT = copy(T)
    truetime = Float32.(cumsum(daysinmonth.(dims(T, Time))))
    for i in spatialidxs(T)
        x = T[i...]
        any(ismissing, x) || continue # this space needs no correction
        mi = findall(!ismissing, x)
        x0 = Array{Float32}(x[mi])
        t0 = truetime[mi]
        sea, res = decompose(t0, x0, lpv)
        c = cosines(truetime, lpv)
        imi = findall(ismissing, x)
        x[imi] .= c[imi]
        fullT[i...] .= x
    end
    return clamp.(fullT, Tmin, Tmax)
end

function cosines(t, lpv::Sinusoidal)
    s = fill(lpv.A[1], length(t))
    for i in 2:length(lpv.A)
        @. s += lpv.A[i] * cos(2π * lpv.fs[i-1] * t + lpv.φ[i])
    end
    return s
end
