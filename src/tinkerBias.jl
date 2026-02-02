"""
    tinker10Bias(lnν, Δm)

Tinker et al.'s halo bias.

*Reference*: Equation (6) and Table 2 of Tinker et al., ApJ, 724, 878 (2010)

# Arguments
- `lnν::Real`: natural logarithm of a threshold, ν, i.e., `lnν` = log(ν), defined by ν ≡ [δc/σ(R,z)]^2. Here, δc = 1.6865 and σ(R,z) is the r.m.s. mass fluctuation within
    a top-hat smoothing of scale R at a redshift `z`.
- `Δm::Real`: overdensity within a spherical region of radius R, whose mean density is equal to Δm times the mean **mass** density of the Universe.

Tinker et al. (2010)'s bias is normalized as

``∫_-∞^∞ dlnν tinker10MF(lnν, z=0, Δm) tinker10Bias(lnν, Δm) = 1``

However, it is not normalized for otinker10MF(lnν, z, Δm) with other `z`.

"""
function tinker10Bias(lnν::Real, Δm::Real)
    y = log10(Δm)
    A = 1.0 + 0.24y * exp(-(4/y)^4)
    a = 0.44y - 0.88
    B, b, c = 0.183, 1.5, 2.4
    C = 0.019 + 0.107y + 0.19*exp(-(4/y)^4)
    ν = exp(lnν)
    ap, bp, cp = a/2, b/2, c/2
    bias = 1 - A * ν^ap / (ν^ap + 1.6865^a) + B * ν^bp + C * ν^cp
end
