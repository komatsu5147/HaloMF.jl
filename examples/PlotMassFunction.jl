using HaloMF
using Plots, LaTeXStrings

# We use some functions in https://github.com/komatsu5147/MatterPower.jl
import MatterPower

# %% Specify a redshift
redshift = 0

# %% Define a function to return a linear matter power spectrum (in units of Mpc^3/h^3)
# as a function of the comoving wavenumber, kovh, in units of h/Mpc.
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

pk(kovh) =
   D1^2 *
   As *
   (kovh * h0 / kpivot)^(ns - 1) *
   (2 * kovh^2 * 2998^2 / 5 / Ωm)^2 *
   MatterPower.t_nowiggle(kovh * h0, ωm, fb)^2 *
   2 *
   π^2 / kovh^3

# %% Alternatively you may read in pre-computed data and define a spline function
# using Dierckx
# pk = Spline1D(tabulated_kovh, tabulated_power_spectrum)

# %% Spefify a range of halo masses, Mh, in units of M⊙/h (M⊙ is the mass of Sun)
lnMh = log(3e11):0.3:log(3e15)
Mh = exp.(lnMh)

# Compute the corresponding radius, Rh, in units of Mpc/h
# ρc is the present-day critical density of the Universe in units of h^2 M⊙/Mpc^3
ρc = 2.775e11
Rh = cbrt.(Mh * 3 / 4π / ρc / Ωm)

# %% Loop over mass
nmass = length(Mh)
dndlnMh = zeros(nmass)
for imass = 1:nmass
   # Compute variance of the mass density fluctuation and its derivative
   # Here is an example using functions in MatterPower.jl
   σ2 = MatterPower.sigma2(pk, Rh[imass])
   dlnσ2dlnRh = Rh[imass] * MatterPower.dsigma2dR(pk, Rh[imass]) / σ2

   # Finally, the mass function per logarithmic mass interval, dn/dlnMh, in units of h^3 Mpc^-3
   # Specify the desired value of the overdensity
   Δm = 200
   z = redshift
   lnν = 2 * log(1.6865) - log(σ2)
   MF = tinker08MF(lnν, z, Δm)
   dndlnMh[imass] = -dlnσ2dlnRh * MF / 4π / Rh[imass]^3
end

# %% Plot results and save to "dndlnm.pdf"
p = plot(
   Mh,
   dndlnMh,
   xaxis = :log,
   yaxis = :log,
   legend = :none,
   m = 2,
   ylab = L"dn/dlnM~[h^3~Mpc^{-3}]",
   xlab = L"M~[h^{-1}~M_\odot]",
)
savefig("dndlnm.pdf")
display(p)
