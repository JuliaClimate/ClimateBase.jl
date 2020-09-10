#=
Handling of time in data as a physical quantity, and time-related data processing
=#
using Statistics, StatsBase
#########################################################################
# Datetime related
#########################################################################
# TODO: Identify monthly, daily, yearly or arbitrary spacing
# to simplify this identification is done exclusively on first 3 time points
# TODO: monthlymean funcion

using Dates
export monthday_indices, maxyearspan, monthspan, daymonth, DAYS_IN_YEAR, monthamount
const DAYS_IN_YEAR = 365.26
millisecond2month(t) = Month(round(Int, t.value / 1000 / 60 / 60 / 24 / 30))
daymonth(t) = day(t), month(t)

maxyearspan(A::AbDimArray, tsamp = temporal_sampling(A)) =
maxyearspan(dims(A, Time).val, tsamp)

"""
    maxyearspan(t::AbstractVector) → i
Find the maximum index `i` of `t` so that the total time covered is a multiple
of years (12 months), assuming monthly spaced data.

    maxyearspan(A::AbDimArray) = maxyearspan(dims(A, Time))
"""
function maxyearspan(times, tsamp = temporal_sampling(times))
    # TODO: Make maxyearspan work with sampling
    length(times) % 12 == 0 && return length(times)
    m = month(times[1])
    findlast(i -> month(times[i]) == m, 1:length(times)) - 1
end

# TODO: bad name
"""
    monthday_indices(times, date = times[1])
Find the indices in `times` (which is a `Vector{Date}`) at which
the date in `times` gives the same day and month as `date`.
"""
function monthday_indices(times, date = times[1])
    d1, m1 = daymonth(date)
    a = findall(i -> daymonth(times[i]) == (d1, m1), 1:length(times))
end

"""
    monthamount(x::AbstractArray{<:AbstractDateTime})
Convert `x` to the amount of months (starting from 0).
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
function monthspan(t::TimeType)
    m = mod1(month(t), 12)
    n = mod1(m+1, 12)
    y = year(t)
    u = m == 12 ? y+1 : y
    d = collect(Date(y, m, 1):Day(1):Date(u, n, 1))[1:end-1]
end


# TODO: better description.
"""
    daycount(t::AbstractArray{<:TimeType}, T = Float32)
Convert a given date time array into measurement units of days:
a real-valued array which counts time in days, always increasing.
"""
function daycount(t::AbstractArray{<:TimeType}, T = Float32)
    ts = temporal_sampling(t)
    if ts == :monthly
        truetime = daysinmonth.(t)
        r = T.(cumsum(truetime))
    elseif ts == :yearly
        # TODO
    elseif ts == :daily
        # TODO
    end
    return r
end
daycount(t::AbstractArray{<:Real}) = t


export temporal_sampling
"""
    temporal_sampling(x) → symbol
Return the temporal sampling type of `x`, which is either an array of `Date`s or
a dimensional array (with `Time` dimension). Possible return values are:
- `:yearly`
- `:monthly`
- `:else`
where `:else` covers the cases of either irregular sampling or daily (or even shorter)
sampling. This function is used to perform proper temporal averages.
"""
temporal_sampling(a::AbDimArray) = temporal_sampling(Array(dims(a, Time)))
function temporal_sampling(t)
    # TODO: Actually code this.
    :monthly
end

#########################################################################
# temporal statistics
#########################################################################
export timemean, timeagg

"""
    timemean(A::ClimArray [, w]) = timeagg(mean, A, w)
Temporal average of `A`.
"""
timemean(A::ClimArray, w = nothing) = timeagg(mean, A, w)


"""
    timeagg(f, A::ClimArray, W = nothing)
Perform a proper temporal aggregation of the function `f` (e.g. `mean, std`)
on `A` (assuming monthly spaced data) where:
* Only full year spans of `A` are included (because most processes are affected by yearly cycle)
* Each month in `a` is weighted with its length in days.

If you don't want these features, just do [`dropagg`](@ref)`(f, A, Time)`.

`W` are possible statistical weights that are used in conjuction to the temporal weighting,
to weight each time point differently.
If they are not a vector (a weight for each time point), then they have to be a dimensional
array of same dimensional layout as `A` (a weight for each data point).

    timeagg(f, t::Vector, x::Vector, w = nothing)
Same as above, but for arbitrary vector `x` accompanied by time vector `t`.
"""
function timeagg(f, A::AbDimArray, w = nothing)
    t = dims(A, Time).val
    w isa AbDimArray && @assert dims(w) == dims(A)
    w isa Vector && @assert length(w) == length(t)
    tsamp = temporal_sampling(t)
    mys = maxyearspan(t, tsamp)
    tw = temporal_weights(t, tsamp)
    W = if isnothing(w)
        tw
    elseif w isa Vector
        tw .* w
    elseif w isa AbDimArray
        _w = dimwise(*, w, ClimArray(tw, (Time(t),)))
        @view _w[Time(1:mys)]
    end
    other = otherdims(A, Time)
    _A = @view A[Time(1:mys)]
    !(w isa AbDimArray) && (fw = weights(view(W, 1:mys)))
    if w isa AbDimArray
        r = map(i -> f(view(_A, i), weights(view(W, i))), otheridxs(A, Time()))
    else
        r = map(i -> f(view(_A, i), fw), otheridxs(A, Time()))
    end
    n = A.name == "" ? "" : A.name*", temporally aggregated with $(string(f))"
    R = ClimArray(r, other, n)
    return R
end

function timeagg(f, a::AbDimArray{T, 1}, exw = nothing) where {T} # version with only time dimension
    t = Array(dims(a, Time))
    return timeagg(f, t, a, exw)
end

function timeagg(f, t::Vector{<:TimeType}, a, exw = nothing) # version with just vectors
    mys = maxyearspan(t)
    # TODO: type stability
    t = @view t[1:mys]
    tw = temporal_weights(t)
    w = if isnothing(exw)
        weights(daysinmonth.(t))
    else
        weights(daysinmonth.(t) .* view(Array(exw), 1:mys))
    end
    y = f(Array(a)[1:mys], w)
end

function temporal_weights(t, tsamp = temporal_sampling(t))
    if tsamp == :monthly
        w = daysinmonth.(t)
    end
    return w
end

#########################################################################
# Monthly/yearly/daily means
#########################################################################
export monthlymean, temporalrange

"""
    monthlymean(A::ClimArray, f = mean) -> B
Create a new array where the temporal daily information has been aggregated to months
using the function `f`.
By convention, the dates of the new array always have day number of `15`.
"""
function monthlymean(A::ClimArray, f = mean)
    t0 = dims(A, Time) |> Array
    finaldate = Date(year(t0[end]), month(t0[end]), 16)
    startdate = Date(year(t0[1]), month(t0[1]), 15)
    t = startdate:Month(1):finaldate
    tranges = temporalrange(t0, Dates.month)
    other = otherdims(A, Time)
    n = A.name == "" ? "" : A.name*", monthly averaged"
    B = ClimArray(zeros(eltype(A), length.(other)..., length(t)), (other..., Time(t)), n)
    for i in 1:length(tranges)
        B[Time(i)] .= dropagg(f, view(A, Time(tranges[i])), Time)
    end
    return B
end

"""
    temporalrange(t::AbstractVector{<:TimeType}}, f = Dates.month) → r
Return a vector of ranges so that each range of indices are values of `t` that
belong in either the same month, year, or day, depending on `f`.
`f` can take the values: `Dates.year, Dates.month, Dates.day` (functions).
"""
function temporalrange(t::AbstractArray{<:TimeType}, f = Dates.month)
    @assert issorted(t) "Sorted time required."
    L = length(t)
    r = Vector{UnitRange{Int}}()
    i, x = 1, f(t[1]) # previous entries
    for j in 2:L
        y = f(t[j])
        x == y && continue
        push!(r, i:(j-1))
        i, x = j, y
    end
    push!(r, i:L) # final range not included in for loop
    return r
end
