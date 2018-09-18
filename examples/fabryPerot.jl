#!/usr/bin/env julia

using PyPlot
close("all")

# Workig directory
path = "/home/leniac/GithubProjects/optics/thin_film_matrix/"
cd(path)

# Load modules
include("ThinFilmOptics.jl") # main calculation program
using Main.ThinFilmOptics: Spectra
include("RIdb.jl") # collection of refractive indexes data
using Main.RIdb: aluminum, air, bk7, chrome, dummy, glass, gold, silicon, silicontemperature, silver, sno2f, h2o, etoh
include("MixingRules.jl") # collection of mixing rules for dielectric functions
using Main.MixingRules: bruggemanspheres, looyengacylinders, looyengaspheres, lorentzlorenz, maxwellgarnettspheres, monecke, gedf, gem
include("indexprofileplot.jl")
include("photonicdispersionplot.jl")

# Wavelength range [nm]
λ = 200:1:1000. # [730.]
# Central/reference wavelength [nm]: must be inside the λ range
λ0 = 700.
# Angle of incidence [degrees], recommended to enter it as an array
θ = [0.] #0:1:10
# Refractive index of materials: elements represent the index position below in the structure (surface down to substrate)
nprofile = [1 2 3 2 3 2 3 2 3 3 3 3 2 3 2 3 2 3 4]
# Thickness vector of each layer [nm], incident and substrate are not included
d = 1/4 * ones.(lastindex(nprofile)-2,1)
# Flag that indicates whether the thicknesses input are optical (1) or geometrical (0), it must be specified for each layer as they can be combined
dflag = ones.(lastindex(nprofile)-2,1)
# Materials profile for each different layer: the nprofile keeps all the information about how this materials will be placed in the structure. Here there incident material is air, then number 2 in nprofile is the looyenga with 0.86 porosity, the number 3 in nprofile is looyenga with 0.54 porosity, the number 4 is silicon. The number of elements inside 'materials' variable must match the maximum value of nprofile
l1 = air(λ) # outermost medium
l2 = looyengaspheres(air(λ),silicon(λ),0.86)
l3 = looyengaspheres(air(λ),silicon(λ),0.54)
l4 = silicon(λ) # substrate medium
materials = [l1 l2 l3 l4]
# polarization, w=1 p-type (TM) wave, and w=0 s-type (TE) wave or a mix of both with any intermediate value
w = 1.
# calculation of the electromagnetic field profile: yes (1) or no (0)
emfflag = 1
# subdivision of each layer for the calculation of the EMF
layerdivision = 10
# calculation of the photonic dispersion: yes (1) or no (0). Only available for crystals without defects
pbgflag = 1

# call main script
results1 = Spectra(λ, λ0, nprofile, dflag, d, w, materials, θ, emfflag, layerdivision, pbgflag)

### Optional examples to plot results

# plot the R, T and A spectra
figure()
plot(λ, results1.R, label="Reflectance")
plot(λ, results1.T, label="Transmittance")
plot(λ, 1 .- (results1.T + results1.R), label="Absorbance")
legend(loc="best")
ax1 = gca()
ax1[:tick_params](which="both", direction="in", pad=10, labelsize=22) # ticks offset
# ax1[:yaxis][:set_ticks_position]("both")
# ax1[:xaxis][:set_ticks_position]("both")
axis("tight")
xlabel("Wavelength [nm]")
ylabel("Reflectance")

# plot the refractive index profile
indexprofileplot(λ, results1.nλ0, results1.d, nprofile, results1.emf, results1.multilayerdepth, θ, λ0)

# plot the EMF pattern
emfield = log10.(copy(dropdims(results1.emf, dims=2)')) # surface plots cannot handle Adjoint yet
figure()
surface = contourf(λ, vec(results1.multilayerdepth), emfield, 20)
cb1_tags = floor.(LinRange(minimum(emfield), maximum(emfield), 5))
cb1 = colorbar(surface, ticks=cb1_tags)
cb1[:set_label]("EMF intensity")
ax2 = gca()
ax2[:tick_params](which="both", direction="in", pad=10, labelsize=22) # ticks offset
ax2[:yaxis][:set_ticks_position]("both")
ax2[:xaxis][:set_ticks_position]("both")
axis("tight")
xlabel("Wavelength [nm]")
ylabel("Depth profile [nm]")

# plot the photonic dispersion
photonicdispersionplot(λ, θ, results1.d, results1.Rp, results1.Rs, results1.R, results1. kblochp, results1.kblochs, results1.kbloch)
