"""
    tinker08MF(lnν, z, Δm)

Halo multiplicity function.

*Reference*: Equation (3, 5-8) and Table 2 of Tinker et al., ApJ, 688, 709 (2008)
- `tinker08MF(lnν, z, Δm)` corresponds to f(σ)/2 in this reference, where σ=√(δc/ν).

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R at a redshift `z`.
- `z::Real`: redshift.
- `Δm::Real`: overdensity within a spherical region of radius R, whose mean density is equal to Δm times the mean **mass** density of the Universe.

Tinker et al. (2008)'s multiplicity function is NOT normalized, i.e,

``∫_-∞^∞ dlnν tinker08MF(lnν, z, Δm) ≠ 1``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function tinker08MF(lnν::Real, z::Real, Δm::Real)
    Δ = [200, 300, 400, 600, 800, 1200, 1600, 2400, 3200]
    A0 = [0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260]
    a0 = [1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66]
    b0 = [2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41]
    c0 = [1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44]
    As = Spline1D(Δ, A0)
    as = Spline1D(Δ, a0)
    bs = Spline1D(Δ, b0)
    cs = Spline1D(Δ, c0)
    α = 10^(-(0.75 / log10(Δm / 75))^1.2)
    if z < 3
        A = As(Δm) * (1 + z)^-0.14
        a = as(Δm) * (1 + z)^-0.06
        b = bs(Δm) * (1 + z)^-α
    else
        A = As(Δm) * 4^-0.14
        a = as(Δm) * 4^-0.06
        b = bs(Δm) * 4^-α
    end
    c = cs(Δm)
    ν = exp(lnν)
    σ = 1.6865 / √ν
    mf = 0.5 * A * ((σ / b)^-a + 1) * exp(-c / σ^2)
end

"""
    tinker10MF(lnν, z, Δm)

Halo multiplicity function.

*Reference*: Equation (8-12) and Table 4 of Tinker et al., ApJ, 724, 878 (2010)
- `tinker10MF(lnν, z, Δm)` corresponds to f(σ)/2 in this reference, where σ=√(δc/ν).

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within
    a top-hat smoothing of scale R at a redshift `z`.
- `z::Real`: redshift.
- `Δm::Real`: overdensity within a spherical region of radius R, whose mean density is equal to Δm times the mean **mass** density of the Universe.

Tinker et al. (2010)'s multiplicity function is normalized at **z=0**, i.e,

``∫_-∞^∞ dlnν tinker10MF(lnν, z=0, Δm) = 1``

It is not normalized for other `z`.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function tinker10MF(lnν::Real, z::Real, Δm::Real)
    Δ = [200, 300, 400, 600, 800, 1200, 1600, 2400, 3200]
    α0 = [0.368, 0.363, 0.385, 0.389, 0.393, 0.365, 0.379, 0.355, 0.327]
    β0 = [0.589, 0.585, 0.544, 0.543, 0.564, 0.623, 0.637, 0.673, 0.702]
    γ0 = [0.864, 0.922, 0.987, 1.09, 1.20, 1.34, 1.50, 1.68, 1.81]
    ϕ0 = [-0.729, -0.789, -0.910, -1.05, -1.20, -1.26, -1.45, -1.50, -1.49]
    η0 =
        [-0.243, -0.261, -0.261, -0.273, -0.278, -0.301, -0.301, -0.319, -0.336]
    αs = Spline1D(Δ, α0)
    βs = Spline1D(Δ, β0)
    γs = Spline1D(Δ, γ0)
    ϕs = Spline1D(Δ, ϕ0)
    ηs = Spline1D(Δ, η0)
    if z < 3
        β = βs(Δm) * (1 + z)^0.20
        ϕ = ϕs(Δm) * (1 + z)^-0.08
        η = ηs(Δm) * (1 + z)^0.27
        γ = γs(Δm) * (1 + z)^-0.01
    else
        β = βs(Δm) * 4^0.20
        ϕ = ϕs(Δm) * 4^-0.08
        η = ηs(Δm) * 4^0.27
        γ = γs(Δm) * 4^-0.01
    end
    ν = exp(lnν)
    σ = 1.6865 / √ν
    mf = 0.5 * αs(Δm) * (1 + (β^2 * ν)^-ϕ) * ν^η * exp(-γ * ν / 2) * √ν
end
