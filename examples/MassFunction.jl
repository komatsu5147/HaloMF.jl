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
