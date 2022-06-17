import StatsBase
export quantile_space, value_space

###########################################################################################
# Quantile space
###########################################################################################
"""
    quantile_space(A, B; n = 50) → Aq, Bq, q
Express the array `B` into the quantile space of `A`. `B, A` must have the same indices.

The array `B` is binned according to the quantile values of the elements of `A`.
`n` bins are created in total and the bin edges `p` are returned.
The `i`-th bin contains data whose `A`-quantile values are ∈ [`q[i]`, `q[i+1]`).
The indices of these `A` values are used to group `B`.

In each bin, the binned values of `A, B` are averaged, resulting in `Aq, Bq`.

The amount of datapoints per quantile is by definition `length(A)/n`.
"""
function quantile_space(A, B; n = 50)
    @assert size(A) == size(B)
    bin_idxs, qs = quantile_decomposition(A, n)
    Aquant, Bquant = averages_from_indices(bin_idxs, A, B)
    return Aquant, Bquant, qs
end

function quantile_decomposition(A, n)
    Avec = vec(A)
    qs = range(0, 1; length = n+1)
    bin_width = Base.step(qs)
    A_cdf = StatsBase.ecdf(Avec)
    A_quantiles = A_cdf.(Avec)
    bin_idxs = [Int[] for _ in 1:length(qs)-1]
    for (i, q) in enumerate(A_quantiles)
        j = Int(q÷bin_width) + 1 # index of quantile bin
        push!(bin_idxs[j], i)
    end
    return bin_idxs, qs
end

function averages_from_indices(idxs, As...)
    map(A -> map(i -> StatsBase.mean(view(vec(A), i)), idxs), As)
end

# TODO:
# function quantile_space(A, B, C; n = 50)
# I am honestly not sure how to go about that.

###########################################################################################
# Quantile space
###########################################################################################
"""
    value_space(A, B; Arange) → Bmeans, bin_indices
Express the array `B` into the value space of `A`. `B, A` must have the same indices.
This means that `A` is binned according to the given `Arange`,
and the same indices that binned `A` are also used to bin `B`.
The `i`-th bin contains data whose `A` values ∈ [`Arange[i]`, `Arange[i+1]`).
In each bin, the binned values of `B` are averaged, resulting in `Bmeans`.

Elements of `A` that are not ∈ `Arange` are skipped.
The returned `bin_indices` are the indices in each bin (hence, the weights
for the means are just `length.(bin_indices)`)

By default `Arange = range(minimum(A), nextfloat(maximum(A)); length = 100)`.
"""
function value_space(A, B; Arange = _default_val_range(A))
    @assert size(A) == size(B)
    @assert issorted(Arange)
    bin_idxs = indices_in_values(A, Arange)
    Bmeans, = averages_from_indices(bin_idxs, B)
    return Bmeans, bin_idxs
end

_default_val_range(A) = range(minimum(A), nextfloat(maximum(A)); length = 100)

function indices_in_values(A, Arange)
    bin_idxs = [Int[] for _ in 1:length(Arange)-1]
    # `li` are linear indices of `A`. `bi` are bin indices of the values of `A`.
    for (li, a) in enumerate(vec(A))
        a > maximum(Arange) && continue
        bi = searchsortedlast(Arange, a)
        bi == 0 && continue
        push!(bin_idxs[bi], li)
    end
    return bin_idxs
end

# 2D version
"""
    value_space(A, B, C; Arange, Brange) → Cmeans, bin_indices
Express the array `C` into the joint value space of `A, B`.
`A, B, C` must have the same indices.

A 2D histogram is created based on the given ranges, the and elements of `C` are binned
according to the values of `A, B`. Then, the elements are averaged, which returns a matrix
`Cmeans`, defined over the joint space S = `Arange × Brange`.
`bin_indices` is also a matrix (with vector elements).
`Cmeans` is `NaN` for bins without any elements.
"""
function value_space(A, B, C;
        Arange = _default_val_range(A), Brange = _default_val_range(B)
    )
    @assert size(A) == size(B) == size(C)
    @assert issorted(Arange)
    bin_idxs_1D = indices_in_values(A, Arange)
    bin_idxs_2D = refine_indices_in_values(bin_idxs_1D, Arange, B, Brange)
    Cmeans, = averages_from_indices(bin_idxs_2D, C)
    return Cmeans, bin_idxs_2D
end


function refine_indices_in_values(bin_idxs_1D, Arange, B, Brange)
    vecB = vec(B)
    bin_idxs = Matrix{Vector{Int}}(undef, length(Arange)-1, length(Brange)-1)
    for k in eachindex(bin_idxs); bin_idxs[k] = Int[]; end
    for (bi, idxs) in enumerate(bin_idxs_1D) # iterate over bins of A
        # notice that `idxs` is a container of linear indices for B!
        for li in idxs
            b = vecB[li]
            b > maximum(Brange) && continue
            bj = searchsortedlast(Brange, b)
            bj == 0 && continue
            push!(bin_idxs[bi, bj], li)
        end
    end
    return bin_idxs
end
