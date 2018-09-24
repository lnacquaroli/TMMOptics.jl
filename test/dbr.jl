# test dbr

function air(λ::AbstractArray{T,M}) where {T<:Number, M}
    N = 1.00029 * ones.(lastindex(λ)) + im .* zeros.(lastindex(λ))
end

#http://www.cvilaser.com/PublicPages/Pages/TechnicalTips.aspx
function glass(λ::AbstractArray{T,M}) where {T<:Number, M}
    N = ones.(length(λ)) .* (1.52 + im * 0)
end # EOF glass(...)

function dummy(λ::AbstractArray{T,N}, x::S, y::S) where {T<:Number, N, S<:Number}
    Nc = ones.(length(λ)) .* (x + im*y)
end # EOF dummy(...)

# Wavelength range [nm]
λ = 200:1:1000. # [730.]
# Central/reference wavelength [nm]: must be inside the λ range
λ0 = 700.
# Angle of incidence [degrees], recommended to enter it as an array
θ = [0.] #0:1:10
# Refractive index of materials: elements represent the index position below in the structure (surface down to substrate)
n = [1 2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3 4]
# Thickness vector of each layer [nm], incident and substrate are not included
d = 1/4 * ones.(lastindex(n)-2,1)
# Flag that indicates whether the thicknesses input are optical (1) or geometrical (0), it must be specified for each layer as they can be combined
dflag = ones.(lastindex(n)-2,1)
# Materials profile for each different layer: the n keeps all the information about how this materials will be placed in the structure. Here there incident material is air, then number 2 in n is the looyenga with 0.86 porosity, the number 3 in n is looyenga with 0.54 porosity, the number 4 is silicon. The number of elements inside 'materials' variable must match the maximum value of n
l1 = air(λ) # outermost medium
l2 = dummy(λ, 2.5, 0.0)
l3 = dummy(λ, 1.5, 0.1)
l4 = dummy(λ, 3.5, 0.5) # substrate medium
materials = [l1 l2 l3 l4]
# polarization, w=1 p-type (TM) wave, and w=0 s-type (TE) wave or a mix of both with any intermediate value
w = 1.
# calculation of the electromagnetic field profile: yes (1) or no (0)
emfflag = 1
# subdivision of each layer for the calculation of the EMF
h = 10
# calculation of the photonic dispersion: yes (1) or no (0). Only available for crystals without defects
pbgflag = 1

# call main script
results1 = Spectra(λ, λ0, n, dflag, d, w, materials, θ, emfflag, h, pbgflag)