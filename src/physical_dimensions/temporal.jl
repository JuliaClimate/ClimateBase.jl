#=
Handling of time in data as a physical quantity, and time-related data processing
=#
using Statistics, StatsBase
export monthday_indices, maxyearspan, daymonth, time_in_days
export temporal_sampling
export timemean, timeagg
export monthlyagg, yearlyagg, temporalrange, seasonalyagg, season
export DAYS_IN_ORBIT, HOURS_IN_ORBIT
#########################################################################
# Datetime related
#########################################################################
using Dates

const DAYS_IN_ORBIT = 365.26
const HOURS_IN_ORBIT = 365.26*24

daymonth(t) = day(t), month(t)

maxyearspan(A::AbstractDimArray, tsamp = temporal_sampling(A)) =
maxyearspan(dims(A, Time).val, tsamp)

"""
    temporal_sampling(x) → symbol
Return the temporal sampling type of `x`, which is either an array of `Date`s or
a dimensional array (with `Time` dimension).

Possible return values are:
- `:hourly`, where the temporal difference between entries is exactly 1 hour.
- `:daily`, where the temporal difference between dates is exactly 1 day.
- `:monthly`, where all dates have the same day, but different month.
- `:yearly`, where all dates have the same month+day, but different year.
- `:other`, which means that `x` doesn't fall to any of the above categories.
"""
temporal_sampling(A::AbstractDimArray) = temporal_sampling(dims(A, Time).val)
temporal_sampling(t::Dimension) = temporal_sampling(t.val)

function temporal_sampling(t::AbstractVector{<:TimeType})

    function issame(f)
        x0 = f(t[1])
        return all(i -> f(t[i]) == x0, 2:length(t))
    end
    samehour, sameday, samemonth, sameyear = map(issame, (hour, day, month, year))

    if !samehour
        # check if they have exactly 1 hour difference
        dh = getproperty.(diff(t), :value)
        return all(isequal(3600000), dh) ? :hourly : :other
    elseif !sameday && samemonth
        return :daily
    elseif !sameday && !samemonth && (day(t[2])-day(t[1])<0 || day(t[3])-day(t[2])<0)
        # this clause checks daily data where the days wrap over the end of the month!
        return :daily
    elseif sameday && !samemonth
        return :monthly
    elseif sameday && samemonth && !sameyear
        return :yearly
    else
        return :other
    end
end
temporal_sampling(t::AbstractVector) = :other
temporal_sampling(t::StepRange{<:Any,Month}) = :monthly
temporal_sampling(t::StepRange{<:Any,Year}) = :yearly
temporal_sampling(t::StepRange{<:Any,Day}) = :daily
temporal_sampling(t::StepRange{<:Any,Hour}) = :hourly
temporal_sampling(t::StepRange{<:Any,<:Any}) = :other

"return the appropriate subtype of `Dates.Period` or `nothing`."
function tsamp2period(tsamp)
    tsamp == :monthly && return Month(1)
    tsamp == :yearly && return Year(1)
    tsamp == :daily && return Day(1)
    tsamp == :hourly && return Hour(1)
    tsamp == :other && return nothing
end


"""
    maxyearspan(A::ClimArray) = maxyearspan(dims(A, Time))
    maxyearspan(t::Vector{<:DateTime}) → i
Find the maximum index `i` of `t` so that `t[1:i]` covers exact(*) multiples of years.

(*) For monthly spaced data `i` is a multiple of `12` while for daily data we find
the largest possible multiple of `DAYS_IN_ORBIT`.
"""
function maxyearspan(times, tsamp = temporal_sampling(times))
    l = length(times)
    if tsamp == :monthly
        l % 12 == 0 && return l
        m = month(times[1])
        x = findlast(i -> month(times[i]) == m, 1:l) - 1
        if x != 0
            return x
        elseif x == 0
            @warn "Caution: data does not cover a full year."
            return l
        end
    elseif tsamp == :yearly
        return length(times)
    elseif tsamp == :daily
        n_max = Int(l÷DAYS_IN_ORBIT)
        nb_years = findlast(n -> round(Int, n * DAYS_IN_ORBIT) ≤ l, 1:n_max)
        if nb_years != nothing
            return round(Int,nb_years * DAYS_IN_ORBIT)-1
        elseif nb_years == nothing
            @warn "Caution: data does not cover a full year."
            return l
        end
    elseif tsamp == :hourly
        n_max = Int(l÷HOURS_IN_ORBIT)
        nb_years = findlast(n -> round(Int, n * HOURS_IN_ORBIT) ≤ l, 1:n_max)
        if nb_years != nothing
            return round(Int,nb_years * HOURS_IN_ORBIT)-1
        elseif nb_years == nothing
            @warn "Caution: data does not cover a full year."
            return l
        end
    else
        error("maxyearspan: not implemented yet for $(tsamp)-sampled data")
    end
end


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
Return a vector of daily spaced `Date`s that span the entire month that `t` belongs in.
"""
function monthspan(t::TimeType)
    m = mod1(month(t), 12)
    n = mod1(m+1, 12)
    y = year(t)
    u = m == 12 ? y+1 : y
    d = collect(Date(y, m, 1):Day(1):Date(u, n, 1))[1:end-1]
end

"""
    time_in_days(t::AbstractArray{<:TimeType}, T = Float32)
Convert a given date time array into measurement units of days:
a real-valued array which counts time in days, always increasing (cumulative).
"""
function time_in_days(t::AbstractArray{<:TimeType}, T = Float32)
    ts = temporal_sampling(t)
    if ts == :monthly
        truetime = daysinmonth.(t)
        return r = T.(cumsum(truetime))
    elseif ts == :yearly
        truetime = daysinmonth.(t)
        return r = T.(cumsum(truetime))
    elseif ts == :daily
        return T.(1:length(t))
    else
        error("Don't know how to find days for sampling $(ts)")
    end
end

#########################################################################
# temporal statistics
#########################################################################
"""
    timemean(A::ClimArray [, w]) = timeagg(mean, A, w)
Temporal average of `A`, see [`timeagg`](@ref).
"""
timemean(A::ClimArray, w = nothing) = timeagg(mean, A, w)


"""
    timeagg(f, A::ClimArray, W = nothing)
Perform a proper temporal aggregation of the function `f` (e.g. `mean, std`) on `A` where:
* Only full year spans of `A` are included, see [`maxyearspan`](@ref)
  (because most processes are affected by yearly cycle,
  and averaging over an uneven number of cycles typically results in artifacts)
* Each month in `A` is weighted with its length in days (for monthly sampled data)

If you don't want these features, just do [`dropagg`](@ref)`(f, A, Time, W)`.

`W` are possible statistical weights that are used in conjuction to the temporal weighting,
to weight each time point differently.
If they are not a vector (a weight for each time point), then they have to be a dimensional
array of same dimensional layout as `A` (a weight for each data point).

See also [`monthlyagg`](@ref), [`yearlyagg`](@ref), [`seasonalyagg`](@ref).

    timeagg(f, t::Vector, x::Vector, w = nothing)
Same as above, but for arbitrary vector `x` accompanied by time vector `t`.
"""
function timeagg(f, A::AbDimArray, w = nothing)
    w isa AbDimArray && @assert dims(w) == dims(A)
    w isa Vector && @assert length(w) == size(A, Time)
    tsamp = temporal_sampling(A)
    r = if tsamp == :daily || tsamp == :other
        timeagg_daily(f, A, w)
    elseif tsamp == :monthly
        timeagg_monthly(f, A, w)
    elseif tsamp == :yearly
        timeagg_yearly(f, A, w)
    end
    R = ClimArray(r, otherdims(A, Time()), A.name)
end

function timeagg_yearly(f, A, w)
     if isnothing(w)
        dropagg(f, A, Time)
    elseif w isa AbDimArray
        map(i -> f(view(A, i), weights(view(W, i))), otheridxs(A, Time()))
    elseif w isa Vector
        fw = weights(w)
        map(i -> f(view(A, i), fw), otheridxs(A, Time()))
    end
end

function timeagg_monthly(f, A::AbDimArray, w)
    t = dims(A, Time).val
    mys = maxyearspan(t, :monthly)
    tw = daysinmonth.(t)
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
    if w isa AbDimArray
        r = map(i -> f(view(_A, i), weights(view(W, i))), otheridxs(A, Time()))
    else
        fw = weights(view(W, 1:mys))
        r = map(i -> f(view(_A, i), fw), otheridxs(A, Time()))
    end
end

function timeagg_daily(f, A::AbDimArray, w)
    t = dims(A, Time).val
    mys = maxyearspan(t)
    _A = view(A, Time(1:mys))
    if w isa AbDimArray
        W = view(w, Time(1:mys))
        r = map(i -> f(view(_A, i), weights(view(W, i))), otheridxs(A, Time()))
    elseif w isa Vector
        fw = weights(view(w, 1:mys))
        r = map(i -> f(view(_A, i), fw), otheridxs(A, Time()))
    elseif isnothing(w)
        r = dropagg(f, _A, Time)
    end
end


function timeagg(f, a::AbDimArray{T, 1}, w = nothing) where {T} # version with only time dimension
    t = dims(a, Time).val
    return timeagg(f, t, Array(a), w)
end

function timeagg(f, T::AbstractVector{<:TimeType}, a::Vector, w = nothing) # version with just vectors
    tsamp = temporal_sampling(T)
    mys = maxyearspan(T, tsamp)
    t = view(T, 1:mys)
    if tsamp == :monthly
        dimw = float.(daysinmonth.(t))
        !isnothing(w) && (dimw .*= view(w, 1:mys))
        return f(view(a, 1:mys), weights(dimw))
    else
        if isnothing(w)
            return f(view(a, 1:mys))
        else
            return f(view(a, 1:mys), weights(view(w, 1:mys)))
        end
    end
end


#########################################################################
# Monthly/yearly/daily/seasonal means
#########################################################################
"""
    monthlyagg(A::ClimArray, f = mean; mday = 15) -> B
Create a new array where the temporal information has been aggregated into months
using the function `f`.
The dates of the new array always have day number of `mday`.
"""
function monthlyagg(A::ClimArray, f = mean; mday = 15)
    t0 = dims(A, Time).val
    startdate = Date(year(t0[1]), month(t0[1]), mday)
    finaldate = Date(year(t0[end]), month(t0[end]), mday+1)
    t = startdate:Month(1):finaldate
    tranges = temporalrange(t0, Dates.month)
    return timegroup(A, f, t, tranges, "monthly")
end

"""
    yearlyagg(A::ClimArray, f = mean) -> B
Create a new array where the temporal information has been aggregated into years
using the function `f`.
By convention, the dates of the new array always have month and day number of `1`.
"""
function yearlyagg(A::ClimArray, f = mean)
    t0 = dims(A, Time).val
    startdate = Date(year(t0[1]), 1, 1)
    finaldate = Date(year(t0[end]), 2, 1)
    t = startdate:Year(1):finaldate
    tranges = temporalrange(t0, Dates.year)
    return timegroup(A, f, t, tranges, "yearly")
end

function timegroup(A, f, t, tranges, name)
    other = otherdims(A, Time)
    n = A.name == Symbol("") ? A.name : Symbol(A.name, ", $(name)")
    B = ClimArray(zeros(eltype(A), length.(other)..., length(t)), (other..., Time(t)), n)
    for i in 1:length(tranges)
        B[Time(i)] .= dropagg(f, view(A, Time(tranges[i])), Time)
    end
    return B
end

"""
    temporalrange(A::ClimArray, f = Dates.month) → r
    temporalrange(t::AbstractVector{<:TimeType}}, f = Dates.month) → r
Return a vector of ranges so that each range of indices are values of `t` that
belong in either the same month, year, day, or season, depending on `f`.
`f` can take the values: `Dates.year, Dates.month, Dates.day` or `season`
(all are functions).

Used in e.g. [`monthlyagg`](@ref), [`yearlyagg`](@ref) or [`seasonalyagg`](@ref).
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
temporalrange(A::AbstractDimArray, f = Dates.month) = temporalrange(dims(A, Time).val, f)


"""
    seasonalyagg(A::ClimArray, f = mean) -> B
Create a new array where the temporal information has been aggregated into seasons
using the function `f`.
By convention, seasons are represented as Dates spaced 3-months apart, where only the
months December, March, June and September are used to specify the date, with day 1.
"""
function seasonalyagg(A::ClimArray, f = mean)
    t0 = dims(A, Time).val
    startdate = to_seasonal_date(t0[1])
    finaldate = to_seasonal_date(t0[end])
    t = startdate:Month(3):finaldate
    tranges = temporalrange(t0, season)
    return timegroup(A, f, t, tranges, "seasonaly")
end

function to_seasonal_date(t)
    y, m = year(t), month(t)
    if m ∈ 3:5
        return Date(y, 3, 1)
    elseif m ∈ 6:8
        return Date(y, 6, 1)
    elseif m ∈ 9:11
        return Date(y, 9, 1)
    elseif m == 12
        return Date(y, 12, 1)
    elseif m ∈ 1:2
        return Date(y-1, 12, 1)
    end
end

"""
    season(date) → s::Int
Return the season of the given date, 1 for DJF, 2 for MAM, 3 for JJA, 4 for SON.
Complements functions like `Dates.year, Dates.month, Dates.day`.
"""
function season(t::Dates.AbstractTime)
    m = month(t)
    if m ∈ 3:5
        2
    elseif m ∈ 6:8
        3
    elseif m ∈ 9:11
        4
    else
        1
    end
end
