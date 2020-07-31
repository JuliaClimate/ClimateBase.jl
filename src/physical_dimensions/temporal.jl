#=
Handling of time in data as a physical quantity, and time-related data processing
=#
#########################################################################
# Datetime related
#########################################################################
# TODO: Identify monthly, daily, yearly or arbitrary spacing
# to simplify this identification is done exclusively on first 3 time points
# TODO: monthlymean funcion

using Dates
export yearly, maxyearspan, monthspan, daymonth, DAYS_IN_YEAR, monthamount
const DAYS_IN_YEAR = 365.26
millisecond2month(t) = Month(round(Int, t.value / 1000 / 60 / 60 / 24 / 30))
daymonth(t) = day(t), month(t)

maxyearspan(A::AbDimArray) = maxyearspan(Array(dims(A, Time)))
"""
    maxyearspan(t::AbstractVector) → i
Find the maximum index `i` of `t` so that the total time covered is a multiple
of years (12 months), assuming monthly spaced data.

    maxyearspan(A::AbDimArray) = maxyearspan(dims(A, Time))
"""
function maxyearspan(times)
    length(times) % 12 == 0 && return length(times)
    m = month(times[1])
    findlast(i -> month(times[i]) == m, 1:length(times)) - 1
end

"""
    yearly(times, date = times[1]) # TODO: bad name
Find the indices in `times` (which is a `Vector{Date}`) at which
the date in `times` gives the same day and month as `date`.
"""
function yearly(times, date = times[1])
    d1, m1 = daymonth(date)
    a = findall(i -> daymonth(times[i]) == (d1, m1), 1:length(times))
    return a
end

"""
    monthamount(x::AbstractArray{<:AbstractDateTime})
Convert `x` to the amount of months starting from `x[1]`.
"""
function monthamount(t)
    L = length(t)
    yd = [year(t[i]) - year(t[1]) for i in 1:L]
    md = [month(t[i]) - month(t[1]) for i in 1:L]
    [yd[i]*12 + md[i] for i in 1:L]
end

"""
    monthspan(t::TimeType)
Return a vector of `Date`s that span the entire month that `t` belongs in.
"""
function monthspan(t)
    m = mod1(month(t), 12)
    n = mod1(m+1, 12)
    y = year(t)
    u = m == 12 ? y+1 : y
    d = collect(Date(y, m, 1):Day(1):Date(u, n, 1))[1:end-1]
end

"""
    daycount(t::AbstractArray{<:TimeType}, T = Float32)
Convert a given date time array into measurement units of days:
a real-valued array which counts time in days, always increasing.
TODO: better description.
"""
function daycount(t::AbstractArray{<:TimeType}, T = Float32)
    ts = temporal_sampling(t)
    if ts == :monthly
        truetime = daysinmonth.(t)
        r = T.(cumsum(truetime))
    elseif
        ts == :yearly
        # TODO
    elseif ts == :daily
        # TODO
    end
    return r
end
daycount(t::AbstractArray{<:Real}) = t

#########################################################################
# temporal statistics
#########################################################################
export timemean, timeagg
# TODO: Make this work for general arrays that do not have time as only
# final dimension. Probably have a look at the "drop" code form DimensionalData.jl.

"`timemean(a) = timeagg(mean, a)`"
timemean(a) = timeagg(mean, a)
timemean(t, a) = timeagg(mean, t, a)


"""
    timeagg(f, a::DimArray, w = nothing)
Perform a proper temporal aggregation of the function `f` (e.g. `mean`)
on `a` (assuming monthly spaced data) where:
* Only full year spans of `a` are included.
* Each month in `a` is weighted with its length in days.

If you don't want these features, just do `dropagg(f, a, Time)`.

`w` are possible statistical weights that are used in conjuction to the monthly weighting.

    timeagg(f, t, a::Vector, w = nothing)
Same as above, but for arbitrary vector accompanied by time vector `t`.
"""
function timeagg(f, a::AbDimArray, exw = nothing)
    mys = maxyearspan(dims(a, Time))
    t = Array(dims(a, Time))[1:mys]
    w = if isnothing(exw)
        weights(daysinmonth.(t))
    else
        # TODO: Here I must use `dimwise` if `exw` has more than 1 dimensions
        weights(daysinmonth.(t) .* view(Array(exw), 1:mys))
    end
    y = zeros(eltype(a), Base.front(size(a)))
    @assert dims(a)[end] isa Time
    r = [1:e for e in Base.front(size(a))] # this line assumes that time is last dimension
    for i in Iterators.product(r...)
        if f == std
            y[i...] = f(view(Array(a), i..., 1:mys), w; corrected=true)
        else
            y[i...] = f(view(Array(a), i..., 1:mys), w)
        end
    end
    DimensionalArray(y, Base.front(dims(a)), a.name)
end

function timeagg(f, a::AbDimArray{T, 1}, exw = nothing) where {T} # version with only time dimension
    t = Array(dims(a, Time))
    return timeagg(f, t, a, exw)
end

function timeagg(f, t::Vector{<:TimeType}, a, exw = nothing) # version with just vectors
    mys = maxyearspan(t)
    t = t[1:mys]
    w = if isnothing(exw)
        weights(daysinmonth.(t))
    else
        weights(daysinmonth.(t) .* view(Array(exw), 1:mys))
    end
    y = f(Array(a)[1:mys], w)
end
