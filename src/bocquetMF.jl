"""
    bocquetMFhy(lnν, z)

Halo multiplicity function for Δm = 200.

*Reference*: Equation (3, 4) and Table 2 of Bocquet et al., MNRAS, 456, 2361 (2016)
- For the parameters of "M200m" and "Hydro" in Table 2
- `bocquetMFhy(lnν, z)` corresponds to f(σ)/2 in this reference, where σ=√(δc/ν).

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R at a redshift `z`.
- `z::Real`: redshift.

Bocquet et al.'s multiplicity function is NOT normalized, i.e,

``∫_-∞^∞ dlnν bocquetMFhy(lnν, z) ≠ 1``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function bocquetMFhy(lnν::Real, z::Real)
    A0, a0, b0, c0 = 0.228, 2.15, 1.69, 1.30
    Az, az, bz, cz = 0.285, -0.058, -0.366, -0.045
    A = A0 * (1 + z)^Az
    a = a0 * (1 + z)^az
    b = b0 * (1 + z)^bz
    c = c0 * (1 + z)^cz
    ν = exp(lnν)
    σ = 1.6865 / √ν
    mf = 0.5 * A * ((σ / b)^-a + 1) * exp(-c / σ^2)
end

"""
    bocquetMFdm(lnν, z)

Halo multiplicity function for Δm = 200.

*Reference*: Equation (3, 4) and Table 2 of Bocquet et al., MNRAS, 456, 2361 (2016)
- For the parameters of "M200m" and "DMonly" in Table 2
- `bocquetMFhy(lnν, z)` corresponds to f(σ)/2 in this reference, where σ=√(δc/ν).

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R at a redshift `z`.
- `z::Real`: redshift.

Bocquet et al.'s multiplicity function is NOT normalized, i.e,

``∫_-∞^∞ dlnν bocquetMFdm(lnν, z) ≠ 1``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function bocquetMFdm(lnν::Real, z::Real)
    A0, a0, b0, c0 = 0.175, 1.53, 2.55, 1.19
    Az, az, bz, cz = -0.012, -0.040, -0.194, -0.021
    A = A0 * (1 + z)^Az
    a = a0 * (1 + z)^az
    b = b0 * (1 + z)^bz
    c = c0 * (1 + z)^cz
    ν = exp(lnν)
    σ = 1.6865 / √ν
    mf = 0.5 * A * ((σ / b)^-a + 1) * exp(-c / σ^2)
end
