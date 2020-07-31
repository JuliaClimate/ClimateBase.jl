#=
Physics-inspired process simulation and data analysis related
to the albedo
=#
export surface_atmosphere_contributions, total_toa_albedo

"""
    surface_atmosphere_contributions(I, F_toa_⬆, F_s_⬆, F_s_⬇) → α_ATM, α_SFC
Calculate the atmospheric and surface **contributions** of the planetary albedo, so that
the TOA albedo is `α = α_ATM + α_SFC`, using the
simple 1-layer radiative transfer model by Donohoe & Battisti (2011) or G. Stephens (2015).
Stephens' formulas are incorrect and I have corrected them!
"""
function surface_atmosphere_contributions(I, F_toa_⬆, F_s_⬆, F_s_⬇)
    R = F_toa_⬆ ./ I   # planetary albedo (system reflectance)
    T = F_s_⬇ ./ I     # system transmisttance
    α = F_s_⬆ ./ F_s_⬇ # surface albedo

    # Formulas by Graeme, which are wrong!
    # τ = @. (1 - α*R) / (1 - (α^2) * (T^2)) # atmosphere transmittance
    # r = @. R - τ*α*T # atmospheric contribution to albedo

    # My calculations:
    r = @. (R - α*T^2) / (1 - (α*T)^2) # atmospheric contribution to albedo
    t = @. T*(1 - r*α)                 # atmospheric transmittance
    s = @. (α*t^2) / (1 - r*α)         # surface contribution to albedo
    return r, s
end

"""
    total_toa_albedo(a, s, t) → α
Combine given atmosphere albedo `a`, surface albedo `s` and atmosphere transmittance `t`
into a total top-of-the-atmosphere albedo `α` according to the model of Donohoe & Battisti (2011).
"""
total_toa_albedo(a, s, t) = a + s*t^2/(1-a*s)
