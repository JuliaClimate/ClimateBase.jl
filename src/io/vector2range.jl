#########################################################################
# Making vectors → ranges
#########################################################################
function vector2range(x::AbstractVector{<:Real})
    length(x) < 3 && return x
    dx = x[2]-x[1]
    for i in 3:length(x)
        x[i]-x[i-1] ≠ dx && return x # if no constant step, return array as is
    end
    # do not check value equality, only difference equality
    return x[1]:dx:x[end]
end

function vector2range(t::AbstractVector{<:Dates.AbstractTime})
    tsamp = temporal_sampling(t)
    period = tsamp2period(tsamp)
    isnothing(period) && return t
    t1 = tsamp == :hourly ? t[1] : Date(t[1])
    tf = tsamp == :hourly ? t[end] : Date(t[end])
    r = t1:period:tf
    return r == t ? r : t # final safety check to ensure equal values
end

vector2range(r) = r
