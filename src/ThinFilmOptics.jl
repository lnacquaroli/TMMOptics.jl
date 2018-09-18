"""
Performs the calculation of reflectance, transmittance, absorptance, electric field and photonic band gap (for photonic crystals only) using the transfer matrix formalism.

Usage:
    results = Spectra(λ, λ0, nlayers, dflags, dinput, polarization, materials, θ, emfflag, emflayerdivision, pbgflag)

Input:
    λ: wavelength range [nm]
    λ0: central wavelength [nm]
    nlayers: index profile of the multilayer sequence
    dflags: flag to indicate the type of dinput input (optical o physic)
    dinput: array with d, physical dinput in nm
    polarization: polarization (1-->p, 0-->s, between 0 and 1 calculations are averaged)
    materials: materials inside the multilayer
    θ: angle of incidence [degrees]
    emfflag: flag that indicates if the calculation of electric field is required (1) or not (0). Default is 0.
    emflayerdivision: number of sublayers in which each layer will be divided to calculate the electric field. Default is 1.
    pbgflag: flag that indicates if photonic band gap must be calculated (1) or not (0). Only valid for periodic structures. Default is 0.

Output:
    results: structure with the following fields:
        Rp: p-wave complex reflectance
        Rs: s-wave complex reflectance
        R: polarization averaged reflectance
        Tp: p-wave complex transmittance
        Ts: s-wave complex transmittance
        T: polarization averaged transmittance
        rrp: p-wave complex reflection coefficient
        rrs: s-wave complex reflection coefficient
        ttp: p-wave complex transmission coefficient
        tts: s-wave complex transmission coefficient
        emfp: p-wave electric field distribution
        emfs: s-wave electric field distribution
        emf: polarization averaged electric field distribution
        d: thickness of each layer [nm]
        kblochp: p-wave Bloch dispersion wavevectors
        kblochs: s-wave Bloch dispersion wavevectors
        kbloch: polarization averaged Bloch dispersion wavevectors
        multilayerdepth: depth profile in the multilayer (taking into account the emflayerdivision for the emf plot) [nm]
        nλ0: refractive indexes profile at the central wavelength reflectance
        ηs: admittance of the whole structure for s-polarization
        ηp: admittance of the whole structure for p-polarization
author: lnacquaroli
"""

module ThinFilmOptics

using LinearAlgebra

# export the results structure
export Spectra

# define the output type
struct Spectra{N1<:Float64, N2<:Number, N3<:ComplexF64, N4<:Number, N5<:Number}
    Rp::Array{N1}; Rs::Array{N1}; R::Array{N1}; Tp::Array{N1}; Ts::Array{N1}; T::Array{N1};
    rrp::Array{N3}; rrs::Array{N3}; ttp::Array{N3}; tts::Array{N3};
    emfp::Array{N1}; emfs::Array{N1}; emf::Array{N1};
    d::Array{N2};
    kblochp::Array{N5}; kblochs::Array{N5}; kbloch::Array{N5};
    multilayerdepth::Array{N1};
    nλ0::Array{N4};
    ηs::Array{N3}; ηp::Array{N3}
end

function Spectra(λ::AbstractArray{A1,M1}, λ0::M7, nlayers::AbstractArray{A2,M2}, dflags::AbstractArray{A3,M3}, dinput::AbstractArray{A4,M4}, polarization::A5, materials::AbstractArray{A6,M5}, θ::AbstractArray{A7,M6}, emfflag::A9=0, emflayerdivision::A10=1, pbgflag::A11=0) where {A1<:Number, M1, A2<:Number, M2, A3<:Number, M3, A4<:Number, M4, A5<:Number, M5, A6<:Number, M6, A7<:Number, A9<:Number, A10<:Number, A11<:Number, M7<:Number}

    # Load usable variables, work in columns and radians
    λ = makecolumns(λ)
    θ = deg2rad.(θ)
    θ = makecolumns(θ)
    dflags = makecolumns(dflags)
    nlayers = makecolumns(nlayers)
    dinput = makecolumns(dinput)

    # To make a faster clean computation
    numberlayers = lastindex(nlayers)
    numberlambda = lastindex(λ)
    numberangles = lastindex(θ)

    # Check variables sizes
    @assert maximum(nlayers) == size(materials,2) "The number of refractive indexes in the materials variable should be at least equal to the maximum value inside the media profile"
    @assert numberlayers-2 == lastindex(dflags) "Number of elements of the thickness vector must be equal to that of sequence vector minus two"

    # Build the complex refractive index matrix for the whole structures
    indexesprofile = Array{ComplexF64,2}(undef, (numberlambda, numberlayers))
    for s = 1 : numberlayers
        indexesprofile[:,s] = materials[:,nlayers[s]]
    end # for s = 1 : numberlayers

    # Build the array of thickness depending on the input
    d = Array{Float64,1}(undef, numberlayers-2)
    # Find λ closest to λ0
    idxλ0 = findmin(abs.(λ .- λ0))[2][1]
    # If the input are optical thicknesses
    for s = 2 : numberlayers-1
        # convert fraction of central wavelength into physical thickness or pass the input
        d[s-1] = dflags[s-1] == 1 ? dinput[s-1] * λ0 / real.(indexesprofile[idxλ0,s]) : dinput[s-1]
    end # for s = 2 : numberlayers-1

    # Remove λ0 from λ and refractive index profile. Since dflags is a vector, we add λ0 to all indexesprofile.
    nλ0 = indexesprofile[idxλ0,:] # nice to return to plot the profile steps

    # Calculation of complex coefficients of reflection, transmission and emf
    Rp, Rs, R, Tp, Ts, T, rrp, rrs, ttp, tts, emfp, emfs, emf, ηs, ηp = reflectiontransmissionemf(indexesprofile, d, λ, θ, polarization, emflayerdivision, emfflag)

    # Provide the multilayer depth considering the emflayerdivision
    temp1 = d / emflayerdivision
    temp2 = temp1 * ones.(1,emflayerdivision) # outer product
    # Cumulative dinput vector: this gives the final depth profile
    multilayerdepth = cumsum([0; temp2[:]], dims=1)
    multilayerdepth = multilayerdepth[1:end-1] # remove last line just for multilayer structures

    # Calculation of photonic band gap for crystals without defects
    if (pbgflag == 1) & (numberlayers > 3)
        kblochp, kblochs, kbloch = photonicdispersion(λ, θ, indexesprofile[:,2:3], [d[1] d[2]], λ0, polarization)
    else
        kbloch = [0.]
        kblochs = [0.]
        kblochp = [0.]
    end

    # Return variables
    Spectra(Rp, Rs, R, Tp, Ts, T, rrp, rrs, ttp, tts, emfp, emfs, emf, d, kblochp, kblochs, kbloch, multilayerdepth, nλ0, ηs, ηp)

end # Spectra(...)

"""
Function that makes column-wise arrays.
"""
function makecolumns end
makecolumns(x::AbstractArray{T,N}) where {T<:Number,N} = size(x,2) > size(x,1) ? x[:] : x
makecolumns(x::UnitRange{T} where T<:Int64) = x
makecolumns(x::LinRange{T} where T<:Float64) = x
makecolumns(x::Float64) = x
makecolumns(x::Int64) = x

"""
Function that computes the transfer matrix method to calculate the reflection and transmission coefficients, their spectra and the EMF. In this function the electromagnetic field is calculated (emfflag = 1).
"""
function reflectiontransmissionemf(indexprofile::AbstractArray{U,P}, d::AbstractArray{V,M}, λ::AbstractArray{V,O}, θ::AbstractArray{V,Q}, w::W, nsl::S, emfflag::S2) where {U<:ComplexF64, P, V<:Float64, M, O, Q, W<:Number, S<:Int64, S2<:Number}

    numberlambda = lastindex(λ)
    numberangles = lastindex(θ)
    totalnumberlayers = lastindex(d) * nsl # total number of layers

    # Calculation of complex coefficients of reflection, transmission and emf
    tts = Array{ComplexF64,2}(undef, (numberlambda, numberangles))
    ttp = similar(tts)
    rrs = similar(tts)
    rrp = similar(tts)
    T = Array{Float64,2}(undef, (numberlambda, numberangles))
    R = similar(T)
    Ts = similar(T)
    Rs = similar(T)
    Tp = similar(T)
    Rp = similar(T)
    emfs = Array{Float64,3}(undef, (numberlambda, numberangles, totalnumberlayers))
    emfp = similar(emfs)
    emf = similar(emfs)
    ηs = Array{ComplexF64,3}(undef, (numberlambda, numberangles, size(indexprofile,2)))
    ηp = similar(ηs)

    for l = 1:numberlambda, a = 1:numberangles

        # Number of layers in the structure
        N = indexprofile[l,:]
        numlay = lastindex(N)

        # Calculation of the optical transfer matrix for all layers
        ηs[l,a,:], ηp[l,a,:], Ms, Mp, δ = tmatrix(N, d, λ[l], θ[a], numlay)

        # Compute the Fresnell coefficients
        rrs[l,a], tts[l,a] = fresnell(ηs[l,a,:], Ms)
        rrp[l,a], ttp[l,a] = fresnell(ηp[l,a,:], Mp)

        Rs[l,a] = real.(rrs[l,a]*conj.(rrs[l,a]))
        Rp[l,a] = real.(rrp[l,a]*conj.(rrp[l,a]))
        R[l,a] = w * Rp[l,a] + (1-w) * Rs[l,a]
        Ts[l,a] = real.(ηs[l,a,1] * ηs[l,a,end] * tts[l,a] * conj.(tts[l,a]))
        Tp[l,a] = real.(ηp[l,a,1] * ηp[l,a,end] * ttp[l,a] * conj.(ttp[l,a]))
        T[l,a] = w * Tp[l,a] + (1-w) * Ts[l,a]

        if emfflag == 1
            # Calculation of the inverse optical transfer for all layers
            Gs11, Gs12 = tmatrixinv(N, d, λ[l], θ[a], numlay, nsl, δ, ηs[l,a,:], Ms)
            Gp11, Gp12 = tmatrixinv(N, d, λ[l], θ[a], numlay, nsl, δ, ηp[l,a,:], Mp)
            # Field intensity respect to the incident beam
            emfs[l,a,:] = emfield(Gs11, Gs12, ηs[l,a,:], Ms)
            emfp[l,a,:] = emfield(Gp11, Gp12, ηp[l,a,:], Mp)
            emf[l,a,:] = w * emfp[l,a,:] + (1-w) * emfs[l,a,:]
        else
            emf = [0.]
            emfp = [0.]
            emfs = [0.]
        end

    end # for l = 1:numberlambda-1, a = 1:numberangles

    return Rp, Rs, R, Tp, Ts, T, rrp, rrs, ttp, tts, emfp, emfs, emf, ηs, ηp

end # function reflectiontransmissionemf(...)

"""
Computes the total transfer matrix, and admittance for the whole structure at each wavelenth and angle of incidence
"""
function tmatrix(N::AbstractArray{B1,C1}, d::AbstractArray{B2,C2}, λ::C3, θ::C4, numlay::C5) where {B1<:ComplexF64, C1, B2<:Number, C2, C3<:Float64, C4<:Number, C5<:Number}

    # Warm-up
    ϕ = Array{ComplexF64,1}(undef, numlay)
    ηs = similar(ϕ)
    ηp = similar(ϕ)
    δ = Array{ComplexF64,1}(undef, numlay)
    ϕ[1] = cos.(θ) # cosine Snell's law
    Ms = Matrix{ComplexF64}(1.0*I, 2,2)
    Mp = Matrix{ComplexF64}(1.0*I, 2,2)

    # According to polarization type, calculate the admittance of the medium
    ηp[1] = N[1] / ϕ[1]
    ηs[1] = N[1] * ϕ[1]

    # Calculation of the optical transfer matrix for each layer between the medium and the substrate
    for c = 2 : numlay-1

        # Compute angle of incidence inside each layer according to the cosine Snell law, to avoid cutoff of total internal reflection with complex angles
        # ϕ[c] = asin.( (N[c-1] / N[c] ) * sin.(ϕ[c-1]) )
        ϕ[c] = (1 - (N[c-1] / N[c])^2 * (1-ϕ[c-1]^2) )^0.5 # this is the cosine already

        # Phase shifts for each layer
        δ[c] = 2 * π * N[c] * d[c-1] * ϕ[c] / λ

        # Admittance of the layer c
        ηp[c] = N[c] / ϕ[c]
        ηs[c] = N[c] * ϕ[c]

        # Total transfer matrix
        Ms = Ms * Ω(δ[c], ηs[c])
        Mp = Mp * Ω(δ[c], ηp[c])

    end # for c = 2 : numlay-1

    # Compute the admittance of the substrate
    ϕ[numlay] = (1 - (N[numlay-1] / N[numlay])^2 * (1-ϕ[numlay-1]^2) )^0.5 # this is the cosine already
    ηp[numlay] = N[numlay] / ϕ[numlay]
    ηs[numlay] = N[numlay] * ϕ[numlay]

    return ηs, ηp, Ms, Mp, δ

end # tmatrix(...)

"""
Function that calculates the optical transfer matrix of a layer, φ: phase shift of the layer y: admittance of the layer, Τ: 2x2 optical tranfer matrix.
"""
function Ω(φ, y)
    tempii = cos.(φ)
    Τ = [tempii (-im ./ y .* sin.(φ)); (-im .* y .* sin.(φ)) tempii]
end # Ω(...)

"""
Compute the reflection and transmission coefficients given the admittance and transfer matrix of the whole structure per wavelenth and angle of incidence.
"""
function fresnell(y::AbstractArray{A1,B1}, A::AbstractArray{A2,B2}) where {A1<:ComplexF64, B1, A2<:ComplexF64, B2}
    # reflection coefficient
    rr = (y[1]*A[1,1]-A[2,1]+y[1]*y[end]*A[1,2]-y[end]*A[2,2]) ./ (y[1]*A[1,1]+A[2,1]+y[1]*y[end]*A[1,2]+y[end]*A[2,2])
    # Transmission coefficient
    tt = 2 ./ (y[1]*A[1,1]+A[2,1]+y[1]*y[end]*A[1,2]+y[end]*A[2,2])
    return rr, tt
end # fresnell(...)

"""
Computes the inverse total transfer matrix for the whole structure at each wavelenth and angle of incidence.
"""
function tmatrixinv(N::AbstractArray{B1,C1}, d::AbstractArray{B2,C2}, λ::C3, θ::C4, numlay::C5, nsl::C6, δ::AbstractArray{B3,C7}, η::AbstractArray{B4,C8}, Ma::AbstractArray{B5,C9}) where {B1<:ComplexF64, C1, B2<:Number, C2, C3<:Float64, C4<:Number, C5<:Number, C6<:Number, B3<:ComplexF64, C7, B4<:ComplexF64, C8, B5<:ComplexF64, C9}
    # Calculation of the matrix elements for the EM field
    M1 = Matrix{ComplexF64}(1.0*I, 2,2)
    Uaux1 = Matrix{ComplexF64}(1.0*I, 2,2)
    G11 = Array{ComplexF64,1}(undef, (numlay-2)*nsl)
    G12 = similar(G11)
    temp1 = δ/nsl
    for c = 2 : numlay-1, j = 1 : nsl
        idx1 = nsl * (c-2) + j
        Uaux1 = Ωinv(temp1[c], η[c]) * Uaux1
        M1 = Uaux1 * Ma
        G11[idx1] = M1[1,1]
        G12[idx1] = M1[1,2]
    end # for c = 2 : numlay-1, j = 1 : nsl
    return G11, G12
end # tmatrixinv(...)

"""
Function that calculates the inverse of optical transfer matrix of a layer, φ:  phase shift of the layer, y: admittance of the layer, Τ: 2x2 optical tranfer matrix.
"""
function Ωinv(φ, y)
    tempii = cos.(φ)
    Τ = [tempii (im ./ y .* sin.(φ)); (im .* y .* sin.(φ)) tempii]
end # Ωinv(...)

"""
Compute the electric field distribution.
"""
function emfield(G11::Array{E1,F1}, G12::Array{E1,F1}, η::Array{E1,F2}, Ma::AbstractArray{A2,B2}) where {E1<:ComplexF64, F1, F2, A2<:ComplexF64, B2}
    # Field intensity calculation E0+
    emf = abs2.( (G11+η[end].*G12)./ (0.25 .* (η[1]*Ma[1,1]+Ma[2,1]+η[1]*η[end]*Ma[1,2]+η[end]*Ma[2,2])) )
end

"""
Function that computes the photonic dispersion of ordered structures (crystals only) alternating two different dielectric layers (pbgflag = 1).
"""
function photonicdispersion(λ::AbstractArray{T,M}, θ::AbstractArray{U,N}, Nvec::AbstractArray{V,O}, esp::AbstractArray{W,P}, λ0::S, wave::S) where {T<:Number, M, U<:Number, N, V<:Number, O, W<:Number, P, S<:Number}

    numbertheta = lastindex(θ)
    numberlam = lastindex(λ)

    # find central wavelength
    # aux1 = findall(abs.(λ0-λ) == min(abs.(λ0-λ)))
    aux1 = findmin(abs.(λ .- λ0))[2][1]
    n1 = Nvec[aux1,1]
    n2 = Nvec[aux1,2]

    # frequency range
    f = 2 * pi ./ λ

    Ys1 = Array{ComplexF64,1}(undef,numbertheta)
    Ys2 = similar(Ys1) # use similar
    Yp1 = similar(Ys1)
    Yp2 = similar(Ys1)
    θ1 = similar(Ys1)
    kbloch = Array{ComplexF64}(undef, (numberlam, numbertheta))
    kblochs = similar(kbloch)
    kblochp = similar(kbloch)
    for a = 1 : numbertheta

        # admittances for p and s type
        Ys1 = n1*cos.(θ[a])
        Ys2 = n2*cos.(θ[a])
        Yp1 = n1/cos.(θ[a])
        Yp2 = n2/cos.(θ[a])

        # incidence angle of the second layer with snell's law of cosine
        # θ1[a] = asin.(n1*sin.(θ[a])/n2)
        θ1[a] = (1 - (n1 / n2)^2 * (1- θ[a]^2) )^0.5

        # Bloch wavevector
        kblochs[:,a] = acos.(0.5 .* (2 .* cos.(esp[1] .* f .* n1 .* cos.(θ[a]) .* cos.(f .* esp[2] .* n2 .* θ1[a]) - (Ys2.^2 + Ys1.^2) ./ Ys2 ./ Ys1 .* sin.(f .* esp[1] .* n1 .* cos.(θ[a])) .* sin.(f .* esp[2] .* n2 .* θ1[a]))))

        kblochp[:,a] = acos.(0.5 * (2 .* cos.(esp[1] .* f .* n1 .* cos.(θ[a])) .* cos.(f .* esp[2] .* n2 .* θ1[a]) - (Yp2.^2 + Yp1.^2) ./ Yp2 ./ Yp1 .* sin.(f .* esp[1] .* n1 .* cos.(θ[a])) .* sin.(f .* esp[2] .* n2 .* θ1[a])))

    end
    kbloch = kblochp * wave + kblochs * (1-wave)
    return kblochp, kblochs, kbloch

end # function photonicdispersion(...)

end # ThinFilmOptics
