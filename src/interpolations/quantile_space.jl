import StatsBase
export quantile_space

"""
    quantile_space(A, B; n = 50) → q, Aq, Bq
Express the array `B` into the quantile space of `A`. `B, A` must have the same indices.

The array `B` is binned according to the quantile values of the elements of `A`.
`n` bins are created in total and the bin edges `p` are returned.
The `i`-th bin contains data whose `A`-quantile values are ∈ [`q[i]`, `q[i+1]`).
The indices of these `A` values are used to group `B`.

In each bin, the binned values of `A, B` are averaged, resulting in `Aq, Bq`.
"""
function quantile_space(A, B; n = 50)
    @assert size(A) == size(B)
    bin_idxs, qs = quantile_decomposition(A, n)
    Aquant = map_to_quantile(A, bin_idxs)
    Bquant = map_to_quantile(B, bin_idxs)
    return qs, Aquant, Bquant
end

function quantile_decomposition(A, n)
    Avec = vec(A)
    qs = range(0, 1; length = n+1)
    bin_width = step(qs)
    A_cdf = StatsBase.ecdf(Avec)
    A_quantiles = A_cdf.(Avec)
    bin_idxs = [Int[] for _ in 1:length(qs)-1] # data indices in a given quantile bin
    for (i, q) in enumerate(A_quantiles)
        j = Int(q÷bin_width) + 1 # index of quantile bin
        push!(bin_idxs[j], i)
    end
    return bin_idxs, qs
end

function map_to_quantile(A, bin_idxs)
    map(idxs -> StatsBase.mean(view(vec(A), idxs)), bin_idxs)
end