module TMMOptics

# Wrap results structure
export Spectra

# Define the output type
struct Spectra{N1<:Float64, N2<:Number, N3<:ComplexF64, N4<:Number, N5<:Number}
    Rp::Array{N1}; Rs::Array{N1}; R::Array{N1}; Tp::Array{N1}; Ts::Array{N1}; T::Array{N1};
    ρp::Array{N3}; ρs::Array{N3}; τp::Array{N3}; τs::Array{N3};
    emfp::Array{N1}; emfs::Array{N1}; emf::Array{N1};
    d::Array{N2};
    κp::Array{N5}; κs::Array{N5}; κ::Array{N5};
    ℓ::Array{N1};
    nλ0::Array{N4};
    ηp::Array{N3}; ηs::Array{N3};
    δ::Array{N3}
end

function Spectra(λ::AbstractArray{A1,M1}, λ0::M7, n::AbstractArray{A2,M2}, dflags::AbstractArray{A3,M3}, dinput::AbstractArray{A4,M4}, w::A5, materials::AbstractArray{A6,M5}, θ::AbstractArray{A7,M6}, emfflag::A9=0, h::A10=1, pbgflag::A11=0) where {A1<:Number, M1, A2<:Number, M2, A3<:Number, M3, A4<:Number, M4, A5<:Number, M5, A6<:Number, M6, A7<:Number, A9<:Number, A10<:Number, A11<:Number, M7<:Number}

    # Work in columns and radians
    λ = makecolumns(λ)
    θ = deg2rad.(θ)
    θ = makecolumns(θ)
    dflags = makecolumns(dflags)
    n = makecolumns(n)
    dinput = makecolumns(dinput)

    # Length of few variables
    nLen = lastindex(n)
    λLen = lastindex(λ)
    θLen = lastindex(θ)

    # Check variables sizes
    maximum(n) == size(materials,2) || throw(DimensionMismatch("the number of media (index of refraction) should be at least equal to the maximum value inside the sequence of media"))
    lastindex(dflags) == lastindex(dinput) || throw(DimensionMismatch("lastindex(dflags) not equal to lastindex(d)"))
    nLen-2 == lastindex(dinput) || throw(DimensionMismatch("lastindex(n)-2 not equal to lastindex(dinput)"))
    # nLen-2 == lastindex(dflags) || throw(DimensionMismatch("lastindex(n)-2 not equal to lastindex(dflags)"))

    # Build the complex refractive index matrix for the whole structures
    nseq = Array{ComplexF64,2}(undef, (λLen, nLen))
    for s = 1 : nLen
        nseq[:,s] = materials[:,n[s]]
    end # for s = 1 : nLen

    # Find λ closest to λ0
    idxλ0 = findmin(abs.(λ .- λ0))[2][1]

    # Build the array of thickness depending on the input
    d = Array{Float64,1}(undef, nLen-2)
    # If the input are optical thicknesses
    for s = 2 : nLen-1
        # convert fraction of central wavelength into physical thickness or pass the input
        d[s-1] = dflags[s-1] == 1 ? dinput[s-1] * λ0 / real(nseq[idxλ0,s]) : dinput[s-1]
    end # for s = 2 : nLen-1

    # Nice to return to plot the profile steps
    nλ0 = nseq[idxλ0,:]

    # Calculation of complex coefficients of reflection, transmission and emf
    Rp, Rs, R, Tp, Ts, T, ρp, ρs, τp, τs, emfp, emfs, emf, ηp, ηs, δ = rtemf(nseq, d, λ, θ, w, h, emfflag)

    # Provide the multilayer depth considering the h division
    temp1 = d / h
    temp2 = temp1 * ones.(1,h) # outer product
    # Cumulative dinput vector: this gives the final depth profile
    ℓ = cumsum([0; temp2[:]], dims=1)
    ℓ = ℓ[1:end-1] # remove last line just for multilayer structures

    # Calculation of photonic band gap for crystals without defects
    if (pbgflag == 1) & (nLen > 3)
        κp, κs, κ = pbg(λ, θ, nseq[:,2:3], [d[1] d[2]], idxλ0, w)
    else
        κ = [0.]
        κs = [0.]
        κp = [0.]
    end

    # Return variables
    Spectra(Rp, Rs, R, Tp, Ts, T, ρp, ρs, τp, τs, emfp, emfs, emf, d, κp, κs, κ, ℓ, nλ0, ηp, ηs, δ[:,:,2:end-1])

end # Spectra(...)

"""
Make column-wise input.
"""
function makecolumns end
makecolumns(x::AbstractArray{T,N}) where {T<:Number,N} = size(x,2) > size(x,1) ? x[:] : x
makecolumns(x::UnitRange{T}) where {T<:Int64} = x
makecolumns(x::LinRange{T}) where {T<:Float64}= x
makecolumns(x::T) where {T<:Float64} = x
makecolumns(x::T) where {T<:Int64} = x

"""
Computes the reflection and transmission coefficients, and their spectra. The electromagnetic field is calculated if emfflag = 1.
"""
function rtemf(nseq::AbstractArray{U,P}, d::AbstractArray{V,M}, λ::AbstractArray{V,O}, θ::AbstractArray{V,Q}, w::W, h::S, emfflag::S2) where {U<:ComplexF64, P, V<:Float64, M, O, Q, W<:Number, S<:Int64, S2<:Number}

    λLen = lastindex(λ)
    θLen = lastindex(θ)
    hLen = lastindex(d) * h # total number of layers

    # Calculation of complex coefficients of reflection, transmission and emf
    τs = Array{ComplexF64,2}(undef, (λLen, θLen))
    τp = similar(τs)
    ρs = similar(τs)
    ρp = similar(τs)
    T = Array{Float64,2}(undef, (λLen, θLen))
    R = similar(T)
    Ts = similar(T)
    Rs = similar(T)
    Tp = similar(T)
    Rp = similar(T)
    emfs = Array{Float64,3}(undef, (λLen, θLen, hLen))
    emfp = similar(emfs)
    emf = similar(emfs)
    δ = Array{ComplexF64,3}(undef, (λLen, θLen, size(nseq,2)))
    ηs = similar(δ)
    ηp = similar(δ)

    for l = 1:λLen, a = 1:θLen

        # Number of layers in the structure
        N = nseq[l,:]
        NLen = lastindex(N)

        # Calculation of the optical transfer matrix for all layers
        ηs[l,a,:], ηp[l,a,:], Ωs, Ωp, δ[l,a,:] = tmatrix(N, d, λ[l], θ[a], NLen)

        # Compute the Fresnell coefficients
        ρs[l,a] = ρ(ηs[l,a,:], Ωs)
        ρp[l,a] = ρ(ηp[l,a,:], Ωp)
        τs[l,a] = τ(ηs[l,a,:], Ωs)
        τp[l,a] = τ(ηp[l,a,:], Ωp)

        Rs[l,a] = abs2.(ρs[l,a]) #real(ρs[l,a] * conj(ρs[l,a]))
        Rp[l,a] = abs2.(ρs[l,a]) #real(ρp[l,a] * conj(ρp[l,a]))
        R[l,a] = w * Rp[l,a] + (1-w) * Rs[l,a]
        Ts[l,a] = real(ηs[l,a,1] * ηs[l,a,end]) * abs2.(τs[l,a]) #real(ηs[l,a,1] * ηs[l,a,end] * τs[l,a] * conj(τs[l,a]))
        Tp[l,a] = real(ηp[l,a,1] * ηp[l,a,end]) * abs2.(τp[l,a]) #real(ηp[l,a,1] * ηp[l,a,end] * τp[l,a] * conj(τp[l,a]))
        T[l,a] = w * Tp[l,a] + (1-w) * Ts[l,a]

        if emfflag == 1
            # Calculation of the inverse optical transfer for all layers
            gs11, gs12 = G(N, d, λ[l], θ[a], NLen, h, δ[l,a,:], ηs[l,a,:], Ωs)
            gp11, gp12 = G(N, d, λ[l], θ[a], NLen, h, δ[l,a,:], ηp[l,a,:], Ωp)
            # Field intensity respect to the incident beam
            emfs[l,a,:] = F⃗(gs11, gs12, ηs[l,a,:], Ωs)
            emfp[l,a,:] = F⃗(gp11, gp12, ηp[l,a,:], Ωp)
            emf[l,a,:] = w * emfp[l,a,:] + (1-w) * emfs[l,a,:]
        else
            emf = [0.]
            emfp = [0.]
            emfs = [0.]
        end

    end # for l = 1:λLen-1, a = 1:θLen

    return Rp, Rs, R, Tp, Ts, T, ρp, ρs, τp, τs, emfp, emfs, emf, ηp, ηs, δ

end # function rtemf(...)

"""
Computes the total transfer matrix and admittance for the whole structure at each wavelenth and angle of incidence.
"""
function tmatrix(N::AbstractArray{B1,C1}, d::AbstractArray{B2,C2}, λ::C3, θ::C4, NLen::C5) where {B1<:ComplexF64, C1, B2<:Number, C2, C3<:Float64, C4<:Number, C5<:Number}
    # Warm-up
    ϕ = Array{ComplexF64,1}(undef, NLen)
    ηs = similar(ϕ)
    ηp = similar(ϕ)
    δ = Array{ComplexF64,1}(undef, NLen)
    ϕ[1] = cos(θ) # cosine Snell's law
    Ωs = complex.([1. 0.; 0. 1.]) #Matrix{ComplexF64}(1.0*I, 2,2)
    Ωp = complex.([1. 0.; 0. 1.]) #Matrix{ComplexF64}(1.0*I, 2,2)
    # Calculate the admittance of the first medium for both polarizations
    ηp[1] = N[1] / ϕ[1]
    ηs[1] = N[1] * ϕ[1]
    # Calculation of the optical transfer matrix for each layer between the medium and the substrate
    for c = 2 : NLen-1
        # Compute angle of incidence inside each layer according to the cosine Snell law, to avoid cutoff of total internal reflection with complex angles
        ϕ[c] = (1 - (N[c-1] / N[c])^2 * (1-ϕ[c-1]^2) )^0.5 # this is the cosine already
        # Phase shifts for each layer
        δ[c] = 2 * π * N[c] * d[c-1] * ϕ[c] / λ
        # Admittance of the layer c
        ηp[c] = N[c] / ϕ[c]
        ηs[c] = N[c] * ϕ[c]
        # Total transfer matrix
        Ωs = Ωs * Φ(δ[c], ηs[c])
        Ωp = Ωp * Φ(δ[c], ηp[c])
    end # for c = 2 : NLen-1
    # Compute the angle of incidence
    ϕ[NLen] = (1 - (N[NLen-1] / N[NLen])^2 * (1-ϕ[NLen-1]^2) )^0.5 # this is the cosine already
    # Admittance of the last layer
    ηp[NLen] = N[NLen] / ϕ[NLen]
    ηs[NLen] = N[NLen] * ϕ[NLen]
    return ηs, ηp, Ωs, Ωp, δ
end # tmatrix(...)

"""
Calculates the optical transfer matrix of a layer. φ: phase shift of the layer y: admittance of the layer, Τ: 2x2 optical tranfer matrix.
"""
Φ(φ::T1, s::T1) where {T1<:ComplexF64} = [cos(φ) (-im*sin(φ)/s); (-im*sin(φ)*s) cos(φ)]

"""
Computes the reflection and transmission coefficients given the admittance and transfer matrix of the whole structure per wavelenth and angle of incidence.
"""
# reflection coefficient
ρ(s::AbstractArray{A1,B1}, Ψ::AbstractArray{A2,B2}) where {A1<:ComplexF64, B1, A2<:ComplexF64, B2} = (s[1]*Ψ[1,1] - Ψ[2,1] + s[1]*s[end] * Ψ[1,2] - s[end]*Ψ[2,2]) / (s[1]*Ψ[1,1] + Ψ[2,1] + s[1]*s[end]*Ψ[1,2] + s[end]*Ψ[2,2])
# transmission coefficient
τ(s::AbstractArray{A1,B1}, Ψ::AbstractArray{A2,B2}) where {A1<:ComplexF64, B1, A2<:ComplexF64, B2} = 2 / (s[1]*Ψ[1,1] + Ψ[2,1] + s[1]*s[end]*Ψ[1,2] + s[end]*Ψ[2,2])

"""
Computes the inverse total transfer matrix for the whole structure at each wavelenth and angle of incidence.
"""
function G(N::AbstractArray{B1,C1}, d::AbstractArray{B2,C2}, λ::C3, θ::C4, NLen::C5, h::C6, δ::AbstractArray{B3,C7}, η::AbstractArray{B4,C8}, Ψ::AbstractArray{B5,C9}) where {B1<:ComplexF64, C1, B2<:Number, C2, C3<:Float64, C4<:Number, C5<:Number, C6<:Number, B3<:ComplexF64, C7, B4<:ComplexF64, C8, B5<:ComplexF64, C9}
    # Calculation of the matrix elements for the EM field
    temp0 = Array{ComplexF64,2}(undef,2,2)
    temp1 = complex.([1. 0.; 0. 1.])
    g11 = Array{ComplexF64,1}(undef, (NLen-2)*h)
    g12 = similar(g11)
    temp2 = δ/h
    for c = 2 : NLen-1, j = 1 : h
        idx1 = h * (c-2) + j
        temp1 = Ξ(temp2[c], η[c]) * temp1
        temp0 = temp1 * Ψ
        g11[idx1] = temp0[1,1]
        g12[idx1] = temp0[1,2]
    end # for c = 2 : NLen-1, j = 1 : h
    return g11, g12
end # G(...)

"""
Calculates the inverse of optical transfer matrix of a layer. φ:  phase shift of the layer, y: admittance of the layer, Τ: 2x2 optical tranfer matrix.
"""
Ξ(φ::T1, s::T1) where {T1<:ComplexF64} = [cos(φ) (im*sin(φ)/s); (im*sin(φ)*s) cos(φ)]

"""
Compute the electric field distribution.
"""
F⃗(g11::Array{E1,F1}, g12::Array{E1,F1}, s::Array{E1,F2}, Ψ::AbstractArray{A2,B2}) where {E1<:ComplexF64, F1, F2, A2<:ComplexF64, B2} = abs2.((g11+s[end]*g12)/(0.25*(s[1]*Ψ[1,1]+Ψ[2,1]+s[1]*s[end]*Ψ[1,2]+s[end]*Ψ[2,2])))

"""
Computes the photonic dispersion of ordered structures (crystals only) alternating two different dielectric layers (pbgflag = 1).
"""
function pbg(λ::AbstractArray{T,M}, θ::AbstractArray{U,N}, n::AbstractArray{V,O}, d::AbstractArray{W,P}, aux1::Y, w::S) where {T<:Number, M, U<:Number, N, V<:Number, O, W<:Number, P, Y<:Int64, S<:Number}
    θLen = lastindex(θ)
    λLen = lastindex(λ)
    # Indexes at chosen wavelentgth (reference)
    n1 = n[aux1,1]
    n2 = n[aux1,2]
    # Frequency range
    f = 2 * pi ./ λ
    # Admittances for p and s type
    cosθ = cos.(θ)
    ys1 = n1 .* cosθ
    ys2 = n2 .* cosθ
    yp1 = n1 ./ cosθ
    yp2 = n2 ./ cosθ
    # Angle of incidence of the second layer with snell's law of cosine
    θ1 = (1 .- (n1 ./ n2).^2 .* (1 .- cosθ.^2) ).^0.5
    # Warm-up
    κ = Array{ComplexF64}(undef, (λLen, θLen))
    κs = similar(κ)
    κp = similar(κ)
    for a = 1 : θLen
        # Bloch wavevector
        κs[:,a] = acos.(0.5 .* (2 .* cos.(d[1] .* f .* n1 .* cosθ[a] .* cos.(f .* d[2] .* n2 .* θ1[a]) - (ys2[a].^2 + ys1[a].^2) ./ ys2[a] ./ ys1[a] .* sin.(f .* d[1] .* n1 .* cosθ[a]) .* sin.(f .* d[2] .* n2 .* θ1[a]))))
        κp[:,a] = acos.(0.5 * (2 .* cos.(d[1] .* f .* n1 .* cosθ[a]) .* cos.(f .* d[2] .* n2 .* θ1[a]) - (yp2[a].^2 + yp1[a].^2) ./ yp2[a] ./ yp1[a] .* sin.(f .* d[1] .* n1 .* cosθ[a]) .* sin.(f .* d[2] .* n2 .* θ1[a])))
    end
    κ .= κp * w .+ κs * (1-w)
    return κp, κs, κ
end # function pbg(...)

end # TMMOptics
