# HaloMF

This package contains functions to return a **halo multiplicity function** (hence the name, "HaloMF"), which is the fundamental building block for computing the (comoving) number density of gravitationally collapsed structures, called *halos*, in the Universe.

The package contains
- `tinker08MF(lnν, z, Δm)`: Equation (3, 5-8) and Table 2 of [Tinker et al., ApJ, 688, 709 (2008)](https://iopscience.iop.org/article/10.1086/591439)
- `tinker10MF(lnν, z, Δm)`: Equation (8-12) and Table 4 of [Tinker et al., ApJ, 724, 878 (2010)](https://iopscience.iop.org/article/10.1088/0004-637X/724/2/878)
- `bocquetMFhy(lnν, z)`: Equation (3, 4) with parameters of "M200m, Hydro" in Table 2 of [Bocquet et al., MNRAS, 456, 2361 (2016)](https://academic.oup.com/mnras/article/456/3/2361/1085699)
- `bocquetMFdm(lnν, z)`: Equation (3, 4) with parameters of "M200m, DMonly" in Table 2 of [Bocquet et al., MNRAS, 456, 2361 (2016)](https://academic.oup.com/mnras/article/456/3/2361/1085699)

## Arguments

- `lnν::Real`: natural logarithm of a *threshold*, ν, i.e., `lnν` = log(ν), defined by ``ν ≡ [δc/σ(R,z)]^2``. Here, **δc = 1.6865** and σ(R,z) is the r.m.s. mass fluctuation within a top-hat smoothing of scale R at a redshift `z`.
    - Note that ν in Tinker et al.'s papers is defined as ν = δc/σ(R,z). Be careful about the power of 2. We follow the notation of Sheth & Tormen (2002; see below) of ν = [δc/σ(R,z)]^2.
- `z::Real`: redshift.
- `Δm::Real`: overdensity within a spherical region of radius R, whose mean density is equal to Δm times the mean **mass** density of the Universe, ``ρm(z) = ρm(z=0)(1+z)^3``.
   - The mass enclosed within R is given by ``M = (4π/3)ρm(z)Δm R^3``.

Note:
- Functions `tinker08MF` and `tinker10MF` interporate parameters of the halo multiplicity function in ``200 ≤ Δm ≤ 3200``.
    - Outside this region, the functions use the parameters of `Δm = 200` for `Δm < 200` and `Δm = 3200` for `Δm > 3200`.
- When you wish to compute a halo multiplicity function for an overdensity `Δc` with respect to the **critical** density of the Universe, ``ρc(z) = ρc(z=0)E^2(z)``, use the relation ``Δm = Δc E^z(z)/[Ωm(1+z)^3]`` with a desired value of Δc.
   - ``Ωm = ρm(z=0)/ρc(z=0)``
   - ``E^2(z) = Ωm(1+z)^3 + ΩΛ + (1 - Ωm - ΩΛ)(1+z)^2``.

## Classic multiplicity functions

The package als contains some of the classic halo multiplicity functions
- `psMF(lnν)`: [Press & Schechter, 187, 425 (1974)](http://articles.adsabs.harvard.edu/pdf/1974ApJ...187..425P)
- `stMF(lnν)`: Equation (2) of [Sheth & Tormen, MNRAS, 329, 61 (2002)](https://academic.oup.com/mnras/article/329/1/61/1112679)
- `jenkinsMF(lnν)`: Equation (B3) of [Jenkins et al., MNRAS, 321, 372 (2001)](https://academic.oup.com/mnras/article/321/2/372/980658)

These multiplicity functions are assumed to be "universal", in the sense that they depend only on ν and do not depend explicitly on `z`. This assumption was challenged by Tinker et al. (2008), hence the explicit dependence on `z` (see `tinker08MF` and `tinker10MF`).

## On normalization

Some of the multiplicity functions (`tinker10MF` for `z=0`, `psMF`, `stMF`) are normalized such that

``∫_-∞^∞ dlnν MF(lnν) = 1``

This normalization is equivalent to saying that all the mass in the Universe is contained in collapsed structures, i.e., halos. In terms of the (comoving) number density of halos per mass, dn/dM, one can write this normalization as

``∫_0^∞ dM M dn/dM = ρm``

where ρm is the mean mass density of the (present-day) Universe. This normalization is convenient mathematically but is not necessarily physical; thus, you do not have to pay too much attention to this. It is certainly useful for checking the code when the MF is normalized.

## Relation to the halo mass function

The halo mass function, dn/dM, can be computed from the halo multiplicity function, `MF`, in the following way.

1. The number density of halos per log(ν) is related to the mass function, dn/dM, as

``dM M dn/dM = ρm dlnν MF(lnν)``

where ρm is the mean mass density of the Universe. If dn/dM is the *comoving number density*, i.e., the number of halos per comoving volume (which is almost always assumed in the literature), ρm must be the mean mass density of the Universe **at the present epoch**, i.e., ρm(z=0).

Dividing both sides by dM, one finds

``dn/dlnM = ρm dlnν/dlnM MF(lnν) / M``

2. ν is more directly related to R, as ``ν = [δc/σ(R,z)]^2``. So, it is more convenient to write dn/dM as

``dn/dlnM = dlnR/dlnM dn/dlnR``

and compute dn/dlnR first, via

``dn/dlnR = ρm dlnν/dlnR MF(lnν) / M(R)``

Now, we can use ``M(R) = (4π/3)ρm R^3`` to obtain

``dn/dlnR = (3/4π) dlnν/dlnR MF(lnν) / R^3``

with lnν and dlnν/dlnR related to lnR by

- ``lnν = 2ln(δc) - ln[σ^2(R,z)]``

- ``dlnν/dlnR = -dln[σ^2(R)]/dlnR`` (which is independent of `z`)

3. Once dn/dlnR is obtained as a function of R, one can compute dn/dlnM as

``dn/dlnM = (1/3) dn/dlnR = dlnν/dlnR MF(lnν) / (4πR^3)``

with ``M(R) = (4π/3)ρm R^3``, and ``ρm = 2.775e11 (Ωm h^2) M⊙ Mpc^{-3}`` is the present-day mean mass density of the Universe.

## Example Juia code to compute dn/dM

This example code is avaiable in [examples/MassFunction.jl](https://github.com/komatsu5147/HaloMF.jl/blob/master/examples/MassFunction.jl).

Below you need to supply a linear matter power spectrum `pk(k)` and a function to compute variance of the mass density fluctuation `sigma2(R)`. If you do not have them available already, they can be found in [MatterPower.jl](https://github.com/komatsu5147/MatterPower.jl).

If you would like to generate a nice figure showing dn/dlnM as a function of M, take a look at [examples/PlotMassFunction.jl](https://github.com/komatsu5147/HaloMF.jl/blob/master/examples/PlotMassFunction.jl).
```
using HaloMF

# We use some functions in https://github.com/komatsu5147/MatterPower.jl
import MatterPower

# %% Specify a redshift
redshift = 0

# %% Define a function to return a linear matter power spectrum (in units of Mpc^3/h^3)
# as a function of the comoving wavenumber, k_ov_h, in units of h/Mpc.
# Here is an example using Einstein & Hu's analytical transfer function in MatterPower.jl

# Cosmological parameters
As, ns, kpivot = 2.097e-9, 0.9652, 0.05
Ωm, ΩΛ, Ωb = 0.315, 0.685, 0.049
Ωk = 1 - Ωm - ΩΛ
h0 = 0.674
ωm, fb = Ωm * h0^2, Ωb / Ωm

# Tabulate linear growth factor as a function of scale factor, a
sol = MatterPower.setup_growth(Ωm, ΩΛ)
a = 1 / (1 + redshift)
D1 = sol(a)[1]

pk(k_ov_h) =
   D1^2 *
   As * (k_ov_h * h0 / kpivot)^(ns - 1) *
   (2 * k_ov_h^2 * 2998^2 / 5 / Ωm)^2 *
   MatterPower.t_nowiggle(k_ov_h * h0, ωm, fb)^2 *
   2 * π^2 / k_ov_h^3

# %% Alternatively you may read in pre-computed data and define a spline function
# using Dierckx
# pk = Spline1D(tabulated_k_ov_h, tabulated_power_spectrum)

# %% Spefify a halo mass, Mh, in units of M⊙/h (M⊙ is the mass of Sun)
Mh = 1e14

# Compute the corresponding radius, Rh, in units of Mpc/h
# ρc is the present-day critical density of the Universe in units of h^2 M⊙/Mpc^3
ρc = 2.775e11
Rh = cbrt(Mh * 3 / 4π / ρc / Ωm)

# %% Compute variance of the mass density fluctuation and its derivative
# Here is an example using functions in MatterPower.jl
σ2 = MatterPower.sigma2(pk, Rh)
dlnσ2dlnRh = Rh * MatterPower.dsigma2dR(pk, Rh) / σ2

# %% Finally, the mass function per logarithmic mass interval, dn/dlnMh, in units of h^3 Mpc^-3
# Specify the desired value of the overdensity
Δm = 200
z = redshift
lnν = 2 * log(1.6865) - log(σ2)
MF = tinker08MF(lnν, z, Δm)
dndlnMh = -dlnσ2dlnRh * MF / 4π / Rh^3
```
## Acknowledgment

The functions provided in this package are adapted from [Cosmology Routine Library (CRL)](https://wwwmpa.mpa-garching.mpg.de/~komatsu/crl/), which is based upon work supported in part by NSF under Grant AST-0807649 and PHY-0758153, NASA under Grant NNX08AL43G, and Alfred P. Sloan Research Foundation via a Sloan Fellowship. This work is also supported in part by JSPS KAKENHI Grant Number JP15H05896.
