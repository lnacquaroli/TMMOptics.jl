#!/usr/bin/env julia
clearconsole()

# Workig directory
cd("/home/leniac/GithubProjects/optics/uploaded/")

# Load modules
# using PyPlot, PlotUtils, LaTeXStrings
include("ThinFilmOptics.jl") # main calculation program
using Main.ThinFilmOptics: Spectra
include("drawindexprofile.jl")
include("photonicdispersionplot.jl")

# Close figures
close("all")

# Wavelength range [nm]
λ = 200:1:1000. # [730.]
# Central/reference wavelength [nm]: must be inside the λ range
λ0 = 700.
# Angle of incidence [degrees], recommended to enter it as array with []
θ = [0.] #0:1:10
# Index of materials: elements represent the index position below in the structure (top to bottom, minimum of 3)
nprofile = [1 2 3 2 3 2 3 2 3 3 3 3 2 3 2 3 2 3 4]
# nprofile = [1 2 3]
# Thickness vector of each layer [nm], incident and substrate are not included
d = 1/4 * ones.(lastindex(nprofile)-2,1)
# Flag that indicates whether the thicknesses introduced are optical (1) or geometrical (0), it must be specified for each layer
dflag = ones.(lastindex(nprofile)-2,1)
# Materials profile for each different layer: the nprofile keeps all the information about how this materials will be placed in the structure. Here there incident material is air, then number 2 in nprofile is the looyenga with 0.86 porosity, the number 3 in nprofile is looyenga with 0.54 porosity, and the substrate (4 in nprofile) is silicon
l1 = :(air($(λ))) # outermost medium
l2 = :(looyengaspheres(air($(λ)),silicon($(λ)),0.86))
l3 = :(looyengaspheres(air($(λ)),silicon($(λ)),0.54))
l4 = :(silicon($(λ))) # substrate medium
materials = [l1 l2 l3 l4] # the number of elements inside this variable must match the maximum value of nprofile
# polarization, w=1 p-type (TM) wave, and w=0 s-type (TE) wave or a mix of both with any intermediate value (0.4)
w = 1.
# calculation of the electromagnetic field profile: 1, yes or 0, no
emfflag = 1
# subdivision of each layer for the calculation of the EMF
layerdivision = 10
# calculation of the photonic dispersion: 1, yes or 0, no (only available for crystals without defects)
pbgflag = 1.

# call program
results1 = Spectra(λ, λ0, nprofile, dflag, d, w, materials, θ, emfflag, layerdivision, pbgflag)

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
drawindexprofile(λ, results1.nλ0, results1.d, nprofile, results1.emf, results1.multilayerdepth, θ, λ0)

# plot the EMF pattern
emfield = log10.(copy(dropdims(results1.emf, dims=2)')) # surface plots cannot handle Adjoint yet
figure()
# surface = contourf(output1.multilayerdepth', λ_range, squeeze(output1.emf,2), 40)
surface = contourf(λ, vec(results1.multilayerdepth), emfield, 20)
# cb1_tags = floor.(LinRange(minimum(emfield), maximum(emfield), 5))
cb1 = colorbar(surface)#, ticks=cb1_tags)
cb1[:set_label]("EMF intensity")#(L"$K_0(qr)\exp(i2\omega t)$")
ax2 = gca()
ax2[:tick_params](which="both", direction="in", pad=10, labelsize=22) # ticks offset
ax2[:yaxis][:set_ticks_position]("both")
ax2[:xaxis][:set_ticks_position]("both")
axis("tight")
xlabel("Wavelength [nm]")
ylabel("Depth profile [nm]")

# # plot the admittance pattern
# ηsfield = copy(dropdims(results1.ηs, dims=2)') # surface plots cannot handle Adjoint yet
# figure()
# # surface = contourf(output1.multilayerdepth', λ_range, squeeze(output1.emf,2), 40)
# surface2 = contourf(λ, vec(results1.d), real.(ηsfield[2:end-1,:]), 100)
# # cb1_tags = floor.(LinRange(minimum(emfield), maximum(emfield), 5))
# cb2 = colorbar(surface2)#, ticks=cb1_tags)
# cb2[:set_label]("Admittance")#(L"$K_0(qr)\exp(i2\omega t)$")
# ax3 = gca()
# ax3[:tick_params](which="both", direction="in", pad=10, labelsize=22) # ticks offset
# ax3[:yaxis][:set_ticks_position]("both")
# ax3[:xaxis][:set_ticks_position]("both")
# axis("tight")
# xlabel("Wavelength [nm]")
# ylabel("Depth profile [nm]")
#
# plot the photonic dispersion
photonicdispersionplot(λ, θ, results1.d, results1.Rp, results1.Rs, results1.R, results1. kblochp, results1.kblochs, results1.kbloch)
