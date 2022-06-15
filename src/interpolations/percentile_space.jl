import StatsBase
export percentile_space

"""
    percentile_space(A, B; n = 40) → p, Ap, Bp, weights
Express the given arrays into the percentile (or quantile) space of `A`.

The given data are binned according to the percentile values of the elements of `A`.
`n` bins are created in total and the bin edges `p` are returned.
The `i`-th bin contains data whose A-percentile values are ∈ [`p[i]`, `p[i+1]`).

In each bin, the binned values of `A, B` are averaged, resulting in `Ap, Bp`.
The returned `weights` simply contain how many elements were in each bin, so that
`mean(A) ≈ mean(Ap, Weights(weights))`.
"""
function percentile_space(A, B; n = 40)
    @assert size(A) == size(B)
    bin_idxs, weights, qs = quantile_decomposition(A, n)
    Aquant = map_to_quant(A, bin_idxs)
    Bquant = map_to_quant(B, bin_idxs)
    return qs, Aquant, Bquant, weights
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
    return bin_idxs, length.(bin_idxs), qs
end

function map_to_quant(A, bin_idxs)
    map(idxs -> StatsBase.mean(view(vec(A), idxs)), bin_idxs)
end