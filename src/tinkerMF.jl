"""
    tinker08MF(lnν, z, Δm)

Halo multiplicity function.

*Reference*: Equation (3, 5-8) and Table 2 of Tinker et al., ApJ, 688, 709 (2008)

Note that `tinker08MF(lnν, z)` defined here corresponds to f(σ)/2 in this reference, where σ=√(δc/ν).

# Arguments
- `lnν`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), which is defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R at a redshift `z`.
- `z`: redshift.
- `Δm`: overdensity within a spherical region of radius R, whose mean density is equal to Δm times the mean **mass** density of the Universe.

Tinker et al. (2008)'s multiplicity function is NOT normalized, i.e,

``∫_-∞^∞ dlnν tinker08MF(lnν, z, Δm) ≠ 1``

The halo mass function, dn/dM, can be computed in the following way.

1. The comoving number density of halos per lnν is related to the mass function, dn/dM, as

``dM M dn/dM = ρm dlnν tinker08MF(lnν, z, Δm)``

where ρm is the mean mass density of the Universe **today**, as dn/dM is the comoving number density of halos.
Dividing both sides by dM, one finds

``dn/dlnM = ρm dlnν/dlnM tinker08MF(lnν, z, Δm) / M``

2. ν is more directly related to R, as ``ν = [δc/σ(R,z)]^2``. So, it is more convenient to write dn/dM as

``dn/dlnM = dlnR/dlnM dn/dlnR``

and compute dn/dlnR first, via

``dn/dlnR = ρm dlnν/dlnR tinker08MF(lnν, z, Δm) / M(R)``

Now, we can use ``M(R) = (4π/3)ρm R^3`` to obtain

``dn/dlnR = (3/4π) dlnν/dlnR tinker08MF(lnν, z, Δm) / R^3``

with lnν and dlnν/dlnR related to lnR via

- ``lnν = 2ln(δc) - ln[σ^2(R,z)]``

- ``dlnν/dlnR = -dln[σ^2(R)]/dlnR`` (which is independent of `z`)

3. Once dn/dlnR is obtained as a function of R, one can compute dn/dlnM using ``dlnR/dlnM = 1/3``, and obtain

``dn/dlnM = (1/3) dn/dlnR = dlnν/dlnR `tinker08MF(lnν, z, Δm)` / (4πR^3)``

with ``M(R) = (4π/3)ρm R^3``, and ``ρm = 2.775e11 (Ωm h^2) M⊙ Mpc^{-3}``.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function tinker08MF(lnν, z, Δm)
    Δ = [200, 300, 400, 600, 800, 1200, 1600, 2400, 3200]
    A0 = [0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260]
    a0 = [1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66]
    b0 = [2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41]
    c0 = [1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44]
    As = Spline1D(Δ, A0)
    as = Spline1D(Δ, a0)
    bs = Spline1D(Δ, b0)
    cs = Spline1D(Δ, c0)
    A = As(Δm) * (1 + z)^-0.14
    a = as(Δm) * (1 + z)^-0.06
    α = 10^(-(0.75 / log10(Δm / 75))^1.2)
    b = bs(Δm) * (1 + z)^-α
    c = cs(Δm)
    ν = exp(lnν)
    σ = 1.6865 / √ν
    0.5 * A * ((σ / b)^-a + 1) * exp(-c / σ^2)
end

"""
    tinker10MF(lnν, z, Δm)

Halo multiplicity function.

*Reference*: Equation (8-12) and Table 4 of Tinker et al., ApJ, 724, 878 (2010)

Note that `tinker10MF(lnν, z)` defined here corresponds to f(σ)/2 in this reference, where σ=√(δc/ν).

# Arguments
- `lnν`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), which is defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within
    a top-hat smoothing of scale R at a redshift `z`.
- `z`: redshift.
- `Δm`: overdensity within a spherical region of radius R, whose mean density is equal to Δm times the mean **mass** density of the Universe.

Tinker et al. (2010)'s multiplicity function is normalized at **z=0**, i.e,

``∫_-∞^∞ dlnν tinker10MF(lnν, z, Δm) = 1``

The halo mass function, dn/dM, can be computed in the following way.

1. The comoving number density of halos per lnν is related to the mass function, dn/dM, as

``dM M dn/dM = ρm dlnν tinker10MF(lnν, z, Δm)``

where ρm is the mean mass density of the Universe **today**, as dn/dM is the comoving number density of halos.
Dividing both sides by dM, one finds

``dn/dlnM = ρm dlnν/dlnM tinker10MF(lnν, z, Δm) / M``

2. ν is more directly related to R, as ``ν = [δc/σ(R,z)]^2``. So, it is more convenient to write dn/dM as

``dn/dlnM = dlnR/dlnM dn/dlnR``

and compute dn/dlnR first, via

``dn/dlnR = ρm dlnν/dlnR tinker10MF(lnν, z, Δm) / M(R)``

Now, we can use ``M(R) = (4π/3)ρm R^3`` to obtain

``dn/dlnR = (3/4π) dlnν/dlnR tinker10MF(lnν, z, Δm) / R^3``

with lnν and dlnν/dlnR related to lnR via

- ``lnν = 2ln(δc) - ln[σ^2(R,z)]``

- ``dlnν/dlnR = -dln[σ^2(R)]/dlnR`` (which is independent of `z`)

3. Once dn/dlnR is obtained as a function of R, one can compute dn/dlnM using ``dlnR/dlnM = 1/3``, and obtain

``dn/dlnM = (1/3) dn/dlnR = dlnν/dlnR `tinker10MF(lnν, z, Δm)` / (4πR^3)``

with ``M(R) = (4π/3)ρm R^3``, and ``ρm = 2.775e11 (Ωm h^2) M⊙ Mpc^{-3}``.

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function tinker10MF(lnν, z, Δm)
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
    0.5 * αs(Δm) * (1 + (β^2 * ν)^-ϕ) * ν^η * exp(-γ * ν / 2) * √ν
end
