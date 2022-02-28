#=
Handling of time in data as a physical quantity, and time-related data processing
=#
using Statistics, StatsBase
#########################################################################
# Datetime related
#########################################################################
using Dates

const DAYS_IN_ORBIT = 365.26
const HOURS_IN_ORBIT = 365.26*24

"daymonth(t) = day(t), month(t)"
daymonth(t) = day(t), month(t)

maxyearspan(A::AbstractDimArray, tsamp = temporal_sampling(A)) =
maxyearspan(dims(A, Time).val, tsamp)

"""
    temporal_sampling(x) → symbol
Return the temporal sampling type of `x`, which is either an array of `Date`s or
a dimensional array (with `Time` dimension).

Possible return values are:
- `:hourly`, where the temporal difference between successive entries is exactly 1 hour.
- `:daily`, where the temporal difference between successive entries is exactly 1 day.
- `:monthly`, where all dates have the same day, but different month.
- `:yearly`, where all dates have the same month and day, but different year.
- `:other`, which means that `x` doesn't fall to any of the above categories.
"""
temporal_sampling(A::AbstractDimArray) = temporal_sampling(dims(A, Time).val)
temporal_sampling(t::Dimension) = temporal_sampling(t.val)

function temporal_sampling(t::AbstractVector{<:TimeType})
    function issame(f)
        f === hour && eltype(t) <: Date && return true
        x0 = f(t[1])
        return all(i -> f(t[i]) == x0, 2:length(t))
    end
    samehour, sameday, samemonth, sameyear = map(issame, (hour, day, month, year))
    constdailydiff = all(i -> day(t[i])-day(t[i-1]) ∈ (1, -30, -27, -28, -29), 2:length(t))

    if !samehour
        # check if they have exactly 1 hour difference
        dh = getproperty.(diff(t), :value)
        return all(isequal(3600000), dh) ? :hourly : :other
    elseif !sameday && samemonth
        return :daily
    elseif samehour && !sameday && !samemonth && constdailydiff
        # this clause checks daily data that span more than 1 month
        return :daily
    elseif sameday && !samemonth
        return :monthly
    elseif sameday && samemonth && !sameyear
        return :yearly
    else
        return :other
    end
end
temporal_sampling(::AbstractVector) = :other
temporal_sampling(::StepRange{<:Any,Month}) = :monthly
temporal_sampling(::StepRange{<:Any,Year}) = :yearly
temporal_sampling(::StepRange{<:Any,Day}) = :daily
temporal_sampling(::StepRange{<:Any,Hour}) = :hourly
temporal_sampling(::StepRange{<:Any,<:Any}) = :other

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
    realtime_days(t::AbstractVector{<:TimeType}, T = Float32)
Convert the given _sequential_ date time vector `t` in a vector in a format of "real time",
where time is represented by real numbers, increasing cumulatively, as is the case when
representing a timeseries `x(t)`.
As only differences matter in this form, the returned vector always starts from 0.
The measurement unit of time here is days.

For temporal sampling less than daily return `realtime_milliseconds(t) ./ (24*60*60*1000)`.

Example:
```juliarepl
julia> t = Date(2004):Month(1):Date(2004, 6)
Date("2004-01-01"):Month(1):Date("2004-06-01")

julia> realtime_days(t)
6-element Vector{Float32}:
   0.0
  29.0
  60.0
  90.0
 121.0
 151.0
```
"""
function realtime_days(t::AbstractArray{<:TimeType}, T = Float32)
    @assert issorted(t)
    ts = temporal_sampling(t)
    if ts == :monthly
        truetime = cumsum(daysinmonth.(t))
        return T.(truetime .- truetime[1])
    elseif ts == :yearly
        truetime = cumsum(daysinyear.(t))
        return T.(truetime .- truetime[1])
    elseif ts == :daily
        return T.(0:length(t)-1)
    else
        return T.(realtime_milliseconds(t) ./ 86400000)
    end
end
realtime_days(A) = realtime_days(dims(A, Ti).val)

"""
    realtime_milliseconds(t::AbstractArray{<:TimeType}, T = Float64)
Similar with [`realtime_days`](@ref), but now the measurement unit is millisecond.
For extra accuracy, direct differences in `t` are used.
"""
function realtime_milliseconds(t::AbstractArray{<:TimeType}, T = Float64)
    @assert issorted(t)
    r = cumsum([T(x.value) for x in diff(t)])
    pushfirst!(r, 0)
    return r
end
realtime_milliseconds(A) = realtime_milliseconds(dims(A, Ti).val)


"""
    sametimespan(Xs; mintime, maxtime) → Ys
Given a container of `ClimArray`s, return the same `ClimArray`s but now accessed in the
`Time` dimension so that they all span the same time interval.
Also works for dictionaries with values `ClimArray`s.

Optionally you can provide `Date` or `DatTime` values for the keywords `mintime, maxtime`
that can further limit the minimum/maximum time span accessed.

`sametimespan` takes into consideration the temporal sampling of the arrays for
better accuracy.
"""
function sametimespan(Xs; kwargs...)
   mint, maxt = findsametimespan(Xs; kwargs...)
   map(X -> X[Time(mint..maxt)], Xs)
end
function sametimespan(Xs::AbstractDict; kwargs...)
    mint, maxt = findsametimespan(values(Xs); kwargs...)
    return Dict(k => X[Time(mint..maxt)] for (k, X) in Xs)
end

function findsametimespan(Xs; maxtime = nothing, mintime = nothing)
    mint = maximum(minimum(dims(X, Time).val) for X in Xs)
    maxt = minimum(maximum(dims(X, Time).val) for X in Xs)
    mint = isnothing(mintime) ? mint : max(mint, mintime)
    maxt = isnothing(maxtime) ? maxt : min(maxt, maxtime)

    # Make an intelligent decision for monthly/yearly sampled data
    tsamps = temporal_sampling.(Xs)
    if all(isequal(:monthly), tsamps)
        mint = Date(year(mint), month(mint), 1)
        d = daysinmonth(maxt)
        maxt = Date(year(maxt), month(maxt), d)
    elseif all(isequal(:yearly), tsamps)
        mint = Date(year(mint), 1, 1)
        maxt = Date(year(maxt), 12, 31)
    end
    return mint, maxt
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
This is also done in the case where the time sampling is unknown.

`W` are possible statistical weights that are used in conjuction to the temporal weighting,
to weight each time point differently.
If they are not a vector (a weight for each time point), then they have to be a dimensional
array of same dimensional layout as `A` (a weight for each data point).

See also [`monthlyagg`](@ref), [`yearlyagg`](@ref), [`seasonalyagg`](@ref).

    timeagg(f, t::Vector, x::Vector, w = nothing)
Same as above, but for arbitrary vector `x` accompanied by time vector `t`.
"""
function timeagg(f, A::AbDimArray, w = nothing)
    !hasdim(A, Time) && error("Array does not have `Time` dimension!")
    w isa AbDimArray && @assert size(w) == size(A)
    w isa Vector && @assert length(w) == size(A, Time)
    tsamp = temporal_sampling(A)
    if tsamp == :other
        return dropagg(f, A, Time, w)
    end
    r = if tsamp == :daily
        timeagg_daily(f, A, w)
    elseif tsamp == :monthly
        timeagg_monthly(f, A, w)
    elseif tsamp == :yearly
        timeagg_yearly(f, A, w)
    end
    return ClimArray(r, otherdims(A, Time()); name = A.name)
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
        _w = broadcast_dims(*, w, ClimArray(tw, (Time(t),)))
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

#########################################################################
# Other advanced functions
#########################################################################
"""
    seasonality(t, x; y0 = year(t[1])) → dates, vals
Calculate the "seasonality" of a vector `x` defined with respect to a datetime
vector `t` and return `dates, vals`.
`dates` are all unique dates present in `t` *disregarding the year* (so only the
month and day are compared). The `dates` have as year entry `y0`.
`vals` is a vector of vectors, where `vals[i]` are all the values of `x` that have
day and month same as `dates[i]`. The elements of `vals` are sorted as encountered in `x`.

Typically one is interested in `mean.(vals)`, which actually is the seasonality,
and `std.(vals)` which is the interannual variability at each date.

    seasonality(A::ClimArray) → dates, vals
If given a `ClimArray`, then the array must have only one dimension (time).
"""
function seasonality(t, xs; y0 = year(t[1]))
    dates = unique(daymonth.(t))
    vals = [eltype(xs)[] for u in dates]
    for (τ, x) in zip(t, xs)
        i = findfirst(isequal(daymonth(τ)), dates)
        push!(vals[i], x)
    end
    outdates = [Date(y0, mon, day) for (day, mon) in dates]
    si = sortperm(outdates)
    return outdates[si], vals[si]
end

function seasonality(A::ClimArray; kwargs...)
    @assert length(dims(A)) == 1
    return seasonality(dims(A, Time).val, A.data; kwargs...)
end
