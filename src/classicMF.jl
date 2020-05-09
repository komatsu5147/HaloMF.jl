"""
    psMF(lnν)

Press-Schechter halo multiplicity function.

*Reference*: Press & Schechter, 187, 425 (1974)
- See also Bond, Cole, Efstathiou & Kaiser, ApJ, 379, 440 (1991)
- See Sheth & Tormen, MNRAS, 329, 61 (2002) for the notation used here. `psMF(lnν)` corresponds to νf(ν) in this reference.

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R)]^2. Here, δc = 1.6865 and σ(R) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R.

Press & Schechter's multiplicity function is normalized, i.e,

``∫_-∞^∞ dlnν psMF(lnν) = 1``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function psMF(lnν::Real)
    ν = exp(lnν)
    mf = exp(-ν / 2) * √(ν / 2π)
end

"""
    stMF(lnν)

Sheth-Tormen halo multiplicity function.

*Reference*: Equation (2) of Sheth & Tormen, MNRAS, 329, 61 (2002)
- See also Sheth, Mo & Tormen, MNRAS, 323, 1 (2001)
- `stMF(lnν)` corresponds to νf(ν) in these references

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R)]^2. Here, δc = 1.6865 and σ(R) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R.

Sheth & Tormen's multiplicity function is normalized, i.e,

``∫_-∞^∞ dlnν stMF(lnν) = 1``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function stMF(lnν::Real)
    A, p, q = 0.3222, 0.3, 0.707
    ν = exp(lnν)
    mf = A * (1 + (q * ν)^-p) * exp(-q * ν / 2) * √(q * ν / 2π)
end

"""
    jenkinsMF(lnν)

Jenkins et al.'s halo multiplicity function.

*Reference*: Equation (B3) of Jenkins et al., MNRAS, 321, 372 (2001)
- `jenkinsMF(lnν)` corresponds to f(σ)/2 in this reference, where σ=√(δc/ν).
- We use the fitting function for Δm = 180. Here, Δm is an overdensity within a spherical region of radius R, whose mean density is equal to Δm times the mean **mass** density of the Universe.

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R)]^2. Here, δc = 1.6865 and σ(R) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R.

Jenkins et al.'s multiplicity function is NOT normalized, i.e,

``∫_-∞^∞ dlnν jenkinsMF(lnν) ≠ 1``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function jenkinsMF(lnν::Real)
    ν = exp(lnν)
    σ = 1.6865 / √ν
    mf = 0.5 * 0.301 * exp(-abs(log(1 / σ) + 0.64)^3.82)
end
