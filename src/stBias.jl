# Sheth-Tormen constants (same as in stMF)
const _st_q = 0.707
const _st_Ap = 0.3222
const _st_p = 0.3
const _st_δc = 1.6865

"""
    stBias(lnν)
    stBias1(lnν)

Sheth & Tormen's linear halo bias, b1(lnν).

*Reference*: Equations (68-70) of Cooray & Sheth, Phys. Rept., 372, 1 (2002) with typos corrected.

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R at a redshift `z`.

Sheth & Tormen's bias is normalized as

``∫_-∞^∞ dlnν stMF(lnν) stBias(lnν) = 1``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function stBias(lnν::Real)
    ν = exp(lnν)
    qν = _st_q * ν
    ϵ1 = (qν - 1) / _st_δc
    E1 = 2 * _st_p / _st_δc / (1 + qν^_st_p)
    return 1 + ϵ1 + E1
end
const stBias1 = stBias

"""
    stBias2(lnν)

Sheth & Tormen's second-order halo bias, b2(lnν).

*Reference*: Equations (68-70) of Cooray & Sheth, Phys. Rept., 372, 1 (2002) with typos corrected.

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν (see `stBias`).

Sheth & Tormen's second-order bias satisfies

``∫_-∞^∞ dlnν stMF(lnν) stBias2(lnν) = 0``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function stBias2(lnν::Real)
    ν = exp(lnν)
    qν = _st_q * ν
    a2 = -17 / 21
    ϵ1 = (qν - 1) / _st_δc
    ϵ2 = qν * (qν - 3) / _st_δc^2
    E1 = 2 * _st_p / _st_δc / (1 + qν^_st_p)
    E2 = ((1 + 2 * _st_p) / _st_δc + 2 * ϵ1) * E1
    return 2 * (1 + a2) * (ϵ1 + E1) + ϵ2 + E2
end

"""
    stBias3(lnν)

Sheth & Tormen's third-order halo bias, b3(lnν).

*Reference*: Equations (68-70) of Cooray & Sheth, Phys. Rept., 372, 1 (2002) with typos corrected.

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν (see `stBias`).

Sheth & Tormen's third-order bias satisfies

``∫_-∞^∞ dlnν stMF(lnν) stBias3(lnν) = 0``

This function is based on [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/).
"""
function stBias3(lnν::Real)
    ν = exp(lnν)
    qν = _st_q * ν
    a2 = -17 / 21
    a3 = 341 / 567
    ϵ1 = (qν - 1) / _st_δc
    ϵ2 = qν * (qν - 3) / _st_δc^2
    ϵ3 = qν * (qν^2 - 6 * qν + 3) / _st_δc^3
    E1 = 2 * _st_p / _st_δc / (1 + qν^_st_p)
    E2 = ((1 + 2 * _st_p) / _st_δc + 2 * ϵ1) * E1
    E3 = ((4 * (_st_p^2 - 1) + 6 * _st_p * qν) / _st_δc^2 + 3 * ϵ1^2) * E1
    return 6 * (a2 + a3) * (ϵ1 + E1) + 3 * (1 + 2 * a2) * (ϵ2 + E2) + ϵ3 + E3
end