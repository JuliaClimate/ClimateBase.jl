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

function vector2range(t::AbstractVector{<:Z}) where {Z<:Dates.AbstractTime}
    tsamp = temporal_sampling(t)
    period = tsamp2period(tsamp)
    isnothing(period) && return t
    special_format = Z <: NCDatasets.CFTime.AbstractCFDateTime
    use_base_date = (tsamp == :hourly || special_format)
    ti = use_base_date ? t[1] : Date(t[1])
    tf = use_base_date ? t[end] : Date(t[end])
    r = ti:period:tf
    return r == t ? r : t # final safety check to ensure equal values
end

vector2range(r) = r
