"""
    pbsBias(lnν, MF)
    pbsBias1(lnν, MF)

Linear halo bias, b1(lnν), from the peak-background split (PBS) approximation.

*Reference*: Section 3.3 of Desjacques, Jeong & Schmidt, Phys. Rept., 733, 1 (2018)

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within
    a top-hat smoothing of scale R at a redshift `z`.
- `MF`(lnν): a function which returns a halo multiplicity function with the argument lnν.

The linear PBS bias is normalized as

``∫_-∞^∞ dlnν MF(lnν) pbsBias(MF) = 1``
"""
function pbsBias(lnν::Real, MF)
    δc = 1.6865
    b1 = 1 - ForwardDiff.derivative(MF, lnν) * 2 / δc / MF(lnν)
    return b1
end
const pbsBias1 = pbsBias

"""
    pbsBias2(lnν, MF)

Second-order PBS halo bias, b2(lnν), from the peak-background split (PBS) approximation.

*Reference*: Section 3.3 of Desjacques, Jeong & Schmidt, Phys. Rept., 733, 1 (2018)

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν (see `pbsBias`).
- `MF`(lnν): a function which returns a halo multiplicity function with the argument lnν.

The second-order PBS bias satisfies

``∫_-∞^∞ dlnν MF(lnν) pbsBias2(MF) = 0``
"""
function pbsBias2(lnν::Real, MF)
    δc = 1.6865
    a2 = -17 / 21
    dMFdlnν(x) = ForwardDiff.derivative(MF, x)
    b2 = (
        -2 * (1 + a2) * dMFdlnν(lnν) * 2 / δc / MF(lnν) +
        ForwardDiff.derivative(dMFdlnν, lnν) * 4 / δc^2 / MF(lnν) -
        dMFdlnν(lnν) * 2 / δc^2 / MF(lnν)
    )
    return b2
end