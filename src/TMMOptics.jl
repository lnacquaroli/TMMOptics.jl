module TMMOptics

export thinfilmoptics, PlaneWave, Geometrical, Optical

"""
Definition of the beam type.
"""
abstract type LightSource end
struct PlaneWave{T1} <: LightSource where {T1<:Float64}
    λ; λ0::T1; θ; p::T1;
end

"""
Material type definition with geometrical (physical) and optical thickness.
"""
abstract type Material end
struct Geometrical{T1, T2} <: Material where {T1<:ComplexF64, T2<:Float64}
    n::Array{T1}; d::T2;
end
struct Optical{T1, T2} <: Material where {T1<:ComplexF64, T2<:Float64}
    n::Array{T1}; f::T2;
end

"""
Wrap the output into different types.
"""
abstract type Output end
struct Spectra{T1, T2} <: Output where {T1<:Float64, T2<:ComplexF64}
    Rp::Array{T1}; Rs::Array{T1}; Tp::Array{T1}; Ts::Array{T1};
    ρp::Array{T2}; ρs::Array{T2}; τp::Array{T2}; τs::Array{T2};
end
struct Field{T1} <: Output where {T1<:Float64}
    emfp::Array{T1}; emfs::Array{T1};
end
struct Bloch{T1} <: Output where {T1<:Number}
    κp::Array{T1}; κs::Array{T1};
end
struct Misc{T1, T2, T3} <: Output where {T1<: Float64, T2<:Number, T3<:Number}
    d::Array{T1}; ℓ::Array{T2}; nλ0::Array{T3};
end
struct AdmPhase{T1} <: Output where {T1<:ComplexF64}
    ηp::Array{T1}; ηs::Array{T1}; δ::Array{T1};
end
struct Results{T1, T2, T3, T4, T5} <:Output where {T1<:Spectra, T2<:Field, T3<:Bloch, T4<:Misc, T5<:AdmPhase}
    Spectra::T1; Field::T2; Bloch::T3; Misc::T4; AdmPhase::T5;
end

"""
Data conditioning and main calculations.
"""
function thinfilmoptics(Beam::T1, Layers::Array{T2,N2}, emfflag::T3=0, h::T3=1, pbgflag::T3=0) where {T1<:PlaneWave, T2<:Material, N2, T3<:Int64}
    # Work in columns
    λ::Vector{Float64} = vec(Beam.λ)
    θ::Vector{Float64} = vec(deg2rad.(Beam.θ)) # radians
    # Length of few variables
    nLen::Int64 = size(Layers,2)
    λLen::Int64 = lastindex(λ)
    θLen::Int64 = lastindex(θ)
    # Find λ closest to λ0
    idxλ0::Int64 = findmin(abs.(λ .- Beam.λ0))[2][1]
    # Build the sequence of index of refractions and the array of thickness depending on the input
    d, nseq, nλ0 = buildArrays(Layers, idxλ0, Beam.λ0, nLen, λLen)
    # Provide the multilayer depth considering the h division
    temp1::Vector{Float64} = d[2:end-1] / h
    temp2::Array{Float64,2} = temp1 * ones.(1,h) # outer product
    # Cumulative dinput vector: this gives the final depth profile
    ℓ = cumsum([0; temp2[:]], dims=1)
    # Remove last line from cumsum
    ℓ = ℓ[1:end-1]
    # call transfer matrix method
    tmmout = tmm(nseq, d, λ, θ, Beam.p, emfflag, h, nLen, λLen, θLen)
    # Calculation of photonic band gap for crystals without defects
    if (pbgflag == 1) & (nLen > 3)
        κp, κs = pbg(λ, θ, nseq[idxλ0,2], nseq[idxλ0,3], d[2], d[3], idxλ0, Beam.p)
    else
        κp = [0.]; κs = [0.]
    end
    # Return results
    Results(tmmout[1], tmmout[2], Bloch(κp, κs), Misc(d[2:end-1], ℓ, nλ0), tmmout[3])
end # thinfilmoptics(...)

"""Build the sequence of index of refractions and the array of thickness depending on the input"""
function buildArrays(layers::Array{T0,N0}, idxλ0::T1, λ0::T2, nLen::T1, λLen::T1) where {T0<:Material, N0, T1<:Int64, T2<:Float64}
    d = Array{Float64,1}(undef, nLen)
    nλ0 = Array{ComplexF64,1}(undef, nLen)
    nseq = Array{ComplexF64,2}(undef, (λLen, nLen))
    for s in LinearIndices(layers)
        # Refractive index
        nseq[:,s] = layers[s].n
        # Nice to return to plot the profile steps
        nλ0[s] = layers[s].n[idxλ0]
        # convert fraction of central wavelength into physical thickness or pass the input
        if typeof(layers[s]) == Optical{Complex{Float64},Float64}
            d[s] = layers[s].f * λ0 / real(nλ0[s])
        elseif typeof(layers[s]) == Geometrical{Complex{Float64},Float64}
            d[s] = layers[s].d
        end
    end # for s in LinearIndices(Layers)
    return d, nseq, nλ0
end # buildArrays()

"""
Computes the reflection and transmission coefficients, and their spectra. The electromagnetic field is calculated if emfflag = 1.
"""
function tmm(nseq::AbstractArray{T1,N1}, d::AbstractArray{T2,N2}, λ::AbstractArray{T3,N3}, θ::AbstractArray{T4,N4}, w::T2, emfflag::T5, h::T5, nLen::T5, λLen::T5, θLen::T5) where {T1<:ComplexF64, N1, T2<:Float64, N2, T3<:Number, N3, T4<:Number, N4, T5<:Int64}
    # Useful lengths
    hLen::Int64 = (lastindex(d)-2) * h
    # Warm-up
    τs = Array{ComplexF64,2}(undef, (λLen, θLen))
    τp = similar(τs); ρs = similar(τs); ρp = similar(τs)
    Ts = Array{Float64,2}(undef, (λLen, θLen))
    Tp = similar(Ts); Rs = similar(Ts); Rp = similar(Ts)
    emfs = Array{Float64,3}(undef, (λLen, θLen, hLen))
    emfp = similar(emfs);
    δ = Array{ComplexF64,3}(undef, (λLen, θLen, nLen))
    ηs = similar(δ); ηp = similar(ηs)
    # Calculation of complex coefficients of reflection, transmission and emf
    for l in LinearIndices(λ), a in LinearIndices(θ)
        # Calculation of the optical transfer matrix for all layers
        ηs[l,a,:], ηp[l,a,:], Ψs, Ψp, δ[l,a,:] = tmatrix(nseq[l,:], d, λ[l], θ[a], nLen)
        # Compute the Fresnell coefficients
        ρs[l,a] = ρ(ηs[l,a,1], ηs[l,a,end], Ψs)
        ρp[l,a] = ρ(ηp[l,a,1], ηp[l,a,end], Ψp)
        τs[l,a] = τ(ηs[l,a,1], ηs[l,a,end], Ψs)
        τp[l,a] = τ(ηp[l,a,1], ηp[l,a,end], Ψp)
        # Reflectance and transmittance
        Rs[l,a] = abs2(ρs[l,a])
        Rp[l,a] = abs2(ρp[l,a])
        Ts[l,a] = real(ηs[l,a,1] * ηs[l,a,end]) * abs2(τs[l,a])
        Tp[l,a] = real(ηp[l,a,1] * ηp[l,a,end]) * abs2(τp[l,a])
        # EMF
        if emfflag == 1
            # Calculation of the inverse optical transfer for all layers
            gs11, gs12 = G(nseq[l,:], d, λ[l], θ[a], nLen, h, δ[l,a,:], ηs[l,a,:], Ψs)
            gp11, gp12 = G(nseq[l,:], d, λ[l], θ[a], nLen, h, δ[l,a,:], ηp[l,a,:], Ψp)
            # Field intensity respect to the incident beam
            emfs[l,a,:] = F⃗(gs11, gs12, ηs[l,a,1], ηs[l,a,end], Ψs)
            emfp[l,a,:] = F⃗(gp11, gp12, ηp[l,a,1], ηp[l,a,end], Ψp)
        else
            emf = [0.]; emfp = [0.]
        end
    end # for l in LinearIndices(λ), a in LinearIndices(θ)
    # return results
    return (Spectra(Rp, Rs, Tp, Ts, ρp, ρs, τp, τs), Field(emfp, emfs), AdmPhase(ηp, ηs, δ[:,:,2:end-1]))
end # function tmm(...)

"""
Computes the total transfer matrix and admittance for the whole structure at each wavelenth and angle of incidence.
"""
function tmatrix(N::AbstractArray{T1,N1}, d::AbstractArray{T2,N2}, λ::T3, θ::T4, NLen::T5) where {T1<:ComplexF64, N1, T2<:Float64, N2, T3<:Number, T4<:Number, T5<:Int64}
    # Warm-up
    cosϕ = Array{ComplexF64,1}(undef, NLen)
    δ = Array{ComplexF64,1}(undef, NLen)
    Ω = complex.([1. 0.; 0. 1.])
    # Compute angle of incidence inside each layer according to the cosine Snell law, to avoid cutoff of total internal reflection with complex angles
    cosϕ[1] = cos(θ)
    for c = 2 : NLen
        cosϕ[c] = cosϑ(N[c-1], N[c], cosϕ[c-1])
    end # for c = 2 : NLen
    # Calculate the admittance of the first medium for both polarizations
    ηp = ζₚ.(N, cosϕ)
    ηs = ζₛ.(N, cosϕ)
    # Phase shift angles
    @. δ = 2 * π * N * d * cosϕ / λ
    # Total transfer matrix
    tempp = [[Ω]; Φ.(δ[2:end-1], ηp[2:end-1])]
    Ψp::Array{ComplexF64} = reduce(*, tempp)
    temps = [[Ω]; Φ.(δ[2:end-1], ηs[2:end-1])]
    Ψs::Array{ComplexF64} = reduce(*, temps)
    # return results
    return ηs, ηp, Ψs, Ψp, δ
end # tmatrix(...)

"""
Calculates the optical transfer matrix of a layer. φ: phase shift of the layer y: admittance of the layer, Τ: 2x2 optical tranfer matrix.
"""
Φ(φ::T1, s::T1) where {T1<:ComplexF64} = [cos(φ) (-im*sin(φ)/s); (-im*sin(φ)*s) cos(φ)]

"""
Computes the reflection and transmission coefficients given the admittance and transfer matrix of the whole structure per wavelenth and angle of incidence.
"""
# reflection coefficient
ρ(s1::T1, sf::T1, Ψ::AbstractArray{T1,N2}) where {T1<:ComplexF64, N2} = (s1*Ψ[1,1] - Ψ[2,1] + s1*sf*Ψ[1,2] - sf*Ψ[2,2]) / (s1*Ψ[1,1] + Ψ[2,1] + s1*sf*Ψ[1,2] + sf*Ψ[2,2])
# transmission coefficient
τ(s1::T1, sf::T1, Ψ::AbstractArray{T1,N2}) where {T1<:ComplexF64, N2} = 2 / (s1*Ψ[1,1] + Ψ[2,1] + s1*sf*Ψ[1,2] + sf*Ψ[2,2])

"""
Computes the inverse total transfer matrix for the whole structure at each wavelenth and angle of incidence.
"""
function G(N::AbstractArray{T1,N1}, d::AbstractArray{T2,N2}, λ::T3, θ::T4, NLen::T5, h::T5, δ::AbstractArray{T6,N7}, η::AbstractArray{T6,N8}, Ψ::AbstractArray{T6,N9}) where {T1<:ComplexF64, N1, T2<:Float64, N2, T3<:Number, T4<:Number, T5<:Int64, T6<:ComplexF64, N7, N8, N9}
    # Warm-up
    m0 = Array{ComplexF64,2}(undef,2,2)
    m1 = complex.([1. 0.; 0. 1.])
    g11 = Array{ComplexF64,1}(undef, (NLen-2)*h)
    g12 = similar(g11)
    # Divide the phase shift by h but keep η as is for each layer
    mδ::Array{ComplexF64} = δ/h
    for c = 2 : NLen-1, j = 1 : h
        k::Int64 = h * (c-2) + j
        m1 = Ξ(mδ[c], η[c]) * m1
        m0 = m1 * Ψ
        g11[k] = m0[1,1]
        g12[k] = m0[1,2]
    end # for c = 2 : nN-1, j = 1 : h
    # return results
    return g11, g12
end # G(...)

"""
Calculates the inverse of optical transfer matrix of a layer. φ:  phase shift of the layer, y: admittance of the layer, Τ: 2x2 optical tranfer matrix.
"""
Ξ(φ::T1, s::T1) where {T1<:ComplexF64} = [cos(φ) (im*sin(φ)/s); (im*sin(φ)*s) cos(φ)]

"""
Compute the electric field distribution.
"""
F⃗(g11::Array{T1,N1}, g12::Array{T1,N1}, s1::T1, sf::T1, Ψ::AbstractArray{T1,N2}) where {T1<:ComplexF64, N1, N2} = abs2.((g11 + sf*g12) / (0.25*(s1*Ψ[1,1] + Ψ[2,1] + s1*sf*Ψ[1,2] + sf*Ψ[2,2])))

"""
Computes the photonic dispersion of ordered structures (crystals only) alternating two different dielectric layers (pbgflag = 1).
"""
function pbg(λ::AbstractArray{T1,N1}, θ::AbstractArray{T2,N2}, n1::T3, n2::T3, d1::T4, d2::T4, idxλ0::T5, w::T6) where {T1<:Number, N1, T2<:Number, N2, T3<:ComplexF64, T4<:Float64, T5<:Int64, T6<:Float64}
    # angular frequency
    ω::Vector{Float64} = 2 * π ./ λ
    # Unit cell
    L::Float64 = d1 + d2
    # Angle of incidence of the second layer with Snell's law of cosine
    cosθ1::Vector{ComplexF64} = cos.(θ)
    cosθ2::Vector{ComplexF64} = cosϑ.(n1, n2, cosθ1)
    # Prefactor for Bloch wavevector
    factor_s = admFactor.(ζₛ.(n1, cosθ1), ζₛ.(n2, cosθ2))
    factor_p = admFactor.(ζₚ.(n1, cosθ1), ζₚ.(n2, cosθ2))
    # Warm-up
    cosκp = Array{ComplexF64,2}(undef, (lastindex(ω), lastindex(cosθ1)))
    cosκs = similar(cosκp)
    # Bloch wavevectors
    for a in LinearIndices(cosθ1), b in LinearIndices(ω)
        cosκp[b,a] = cosκ(d1, d2, n1, n2, cosθ1[a], cosθ2[a], ω[b], factor_p[a])
        cosκs[b,a] = cosκ(d1, d2, n1, n2, cosθ1[a], cosθ2[a], ω[b], factor_s[a])
    end # for a in LinearIndices(cosθ1), b in LinearIndices(ω)
    # return results of Bloch wavevectors
    return acos.(cosκp)./L, acos.(cosκs)./L
end # function pbg(...)

"""Prefactor for bloch wavevector"""
admFactor(s1::T1, s2::T2) where {T1<:Number, T2<:Number} = 0.5 * (s1^2 + s2^2) / s1 / s2

"""Bloch wavevector"""
cosκ(d1::T1, d2::T1, n1::T2, n2::T2, cosθ1::T3, cosθ2::T4, ω::T1, f::T5) where {T1<:Float64, T2<:ComplexF64, T3<:Number, T4<:Number, T5<:Number} = cos(d1*ω*n1*cosθ1) * cos(ω*d2*n2*cosθ2) - f * sin(ω*d1*n1*cosθ1) * sin(ω*d2*n2*cosθ2)

"""
Snell's law in cosine form. Returns de cosine already
"""
cosϑ(n1::T1, n2::T1, cosθ::T2) where {T1<:ComplexF64, T2<:Number} = sqrt(1 - (n1/n2)^2 * (1-cosθ^2) )

"""
Admittance of the medium for polarization.
"""
ζₚ(n::T1, cosθ::T2) where {T1<:ComplexF64,T2<:Number} = n / cosθ
ζₛ(n::T1, cosθ::T2) where {T1<:ComplexF64,T2<:Number} = n * cosθ

end # TMMOptics
