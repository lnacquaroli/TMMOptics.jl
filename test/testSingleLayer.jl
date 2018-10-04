# test dbr

function dummy(λ::AbstractArray{T,N}, x::S, y::S) where {T<:Number, N, S<:Number}
    Nc = ones.(length(λ)) .* (x + im*y)
end # EOF dummy(...)

# Define beam
λi = 400 # intial wavelength [nm]
λf = 1000 # final wavelength [nm]
λ = LinRange(λi, λf, λf-λi+1) # wavelength range [nm]
# λ = λi:λf
λ0 = 700. # reference wavelength
θ = [0.] # angle of incidence [degrees]
p = 1. # polatization (s-wave = 0. and p-wave = 1., or any intermediate)
beam = PlaneWave(λ, λ0, θ, p)

# Define layers
l0 = Geometrical(dummy(beam.λ, 1., 0.), 0.)
l1 = Geometrical(dummy(beam.λ, 1.5, 0.1), 70.)
l2 = Optical(dummy(beam.λ, 1.1, 0.1), 1/4)
l3 = Geometrical(dummy(beam.λ, 2.5, 0.6), 0.)

# Refractive index of materials: elements represent the index position below in the structure (surface down to substrate)
nseq = [l0 l1 l2 l3]

# calculation of the electromagnetic field profile: yes (1) or no (0)
emfflag = 1
# subdivision of each layer for the calculation of the EMF
h = 10
# calculation of the photonic dispersion: yes (1) or no (0)
pbgflag = 1

# call main script
results1 = thinfilmoptics(beam, nseq, emfflag, h, pbgflag)
