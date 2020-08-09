#########################################################################
# Insolation
#########################################################################
#=
Code originally from:
https://github.com/climate-machine/RRTMGP.jl/blob/master/src/rte/SolarZenithAngle.jl

with examples at:
https://climate-machine.github.io/RRTMGP.jl/latest/RTE/SolarZenithAngle.html

The physics and details about the algorithm
can be found in Tapio Schneider’s book Physics of the Earth's Climate, section 3.4
http://climate-dynamics.org/wp-content/uploads/2017/04/Climate_Book.pdf
=#
export insolation, monthly_insolation
"""
    insolation(t, ϕ; kwargs...)
Calculate daily averaged insolation in W/m² at given time and latitude `t, φ`.
`φ` is given in **degrees**, and `t` in **days** (real number or date).

Keywords:
```
Ya = DAYS_IN_YEAR # = 365.26 # days
t_VE = 76.0 # days of vernal equinox
S_0 = 1362.0 # W/m^2
γ=23.44
ϖ=282.95
e=0.017 # eccentricity
```
"""
insolation(t::TimeType, φ; kwargs...) = insolation(float(dayofyear(t)), φ; kwargs...)
function insolation(t::Real, ϕ;
        Ya = DAYS_IN_YEAR # = 365.26 # days
        t_VE = 76.0 # days since Jan 1
        S_0 = 1362.0 # W/m^2
        γ=23.44
        ϖ=282.95
        e=0.017
    )

    # convert inputs from degrees to radians
    ϕ = ϕ * π / 180
    γ = γ * π / 180
    ϖ = ϖ * π / 180

    # step 1, calculate the mean anomaly at vernal equinox
    β = sqrt(1 - e^2)
    M_VE = -ϖ + (e + e^3 / 4) * (1 + β) * sin(ϖ)

    # step 2, calculate the mean anomaly
    M = (2 * π * (t - t_VE)) / (Ya) + M_VE

    # step 3, calculate the true anomaly
    A = M + (2 * e - e^3 / 4) * sin(M)

    # step 4, calculate the distance to the sun
    d = (1 - e^2) / (1 + e * cos(A))

    # step 5, calculate the solar longitude
    L_s = A + ϖ

    # step 6, calculate the declination angle
    δ = asin(sin(γ) * sin(L_s))

    # step 7, calculate the sunrise/sunset angle
    T = tan(ϕ) * tan(δ)
    if T >= 1
        η_d = π
    elseif T <= -1
        η_d = 0.0
    else
        η_d = acos(-1 * T)
    end

    # step 8, calculate the daily averaged cos(zenith angle)
    c1 = η_d * sin(ϕ) * sin(δ)
    c2 = cos(ϕ) * cos(δ) * sin(η_d)
    cosbar = (1 / π) * (c1 + c2)

    # step 9, calculate the flux
    F = S_0 * (1 / d)^2 * cosbar
    return F
end
