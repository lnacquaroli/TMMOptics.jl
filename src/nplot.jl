"""
Function that plots the multilayer structure, using PyPlot, PlotUtils and LaTeXStrings. Also using Printf: @sprintf to display some data.
author: lnacquaroli
"""

using PyPlot, PlotUtils, LaTeXStrings
using Printf: @sprintf

function nplot(λ::AbstractArray{T,M}, nλ0::AbstractArray{U,N}, d::AbstractArray{V,O}, n::AbstractArray{W,P}, emf::AbstractArray{X,Q}, τ::AbstractArray{Y,R}, θ::AbstractArray{Z,S}, λ0::A2) where {T<:Number, M, U<:Number, N, V<:Number, O, W<:Number, P, X<:Number, Q, Y<:Number, R, Z<:Number, S, A1<:Number, A2<:Number}

    # check input
    nn = lastindex(n)
    nd = lastindex(d)
    # if size(d,1) != size(nλ0,1)-2
    @assert nd == lastindex(nλ0)-2 "Variables dimensions disagree: length(nλ0)-2 = length(d)"

    # define the offset for the incident and substrate medium
    doffset = 200. # nm

    # Maps each certain to a certain color
    aux2 = unique(n)
    naux2 = lastindex(aux2)
    cm = cgrad(:Spectral)
    β = [cm[i] for i in LinRange(0,1,naux2)] # range(0,stop=1,length=naux2) should work also
    aux3 = [getfield(β[j],i) for i=1:4, j=1:lastindex(β)]

    βassigned = Array{Float64,2}(undef, (4,nn))
    for idx2 = 1 : naux2, idx1 = 1 : nn
        if aux2[idx2] == n[idx1]
            βassigned[:,idx1] = aux3[:,idx2]
        end
    end

    nλ0 = real.(nλ0)
    figure()
    ax_profile = gca()
    # draw incident medium
    ax_profile[:fill_between]([-doffset, 0], 0, nλ0[1], facecolor=βassigned[:,1], interpolate="true", alpha=0.8, edgecolor = "black")

    new_d = cumsum([0; d], dims=1)
    for idx = 1 : nd
        # draw successive layers
        ax_profile[:fill_between]([new_d[idx], new_d[idx+1]], 0, nλ0[idx+1], facecolor=βassigned[:,idx+1], interpolate="true", alpha=.6, edgecolor = "black")
    end
    # draw substrate medium
    ax_profile[:fill_between]([new_d[end], new_d[end]+doffset], 0, nλ0[end], facecolor=βassigned[:,end], interpolate="true", alpha=.6, edgecolor = "black")

    ax_profile[:tick_params](which="both", direction="in", pad=10, labelsize=22) # ticks offset
    ax_profile[:yaxis][:set_ticks_position]("both")
    ax_profile[:xaxis][:set_ticks_position]("both")
    # axis("tight")
    xlabel("Thickness profile [nm]")
    ylabel(L"\mathrm{Refractive index at }\,\, \lambda_0")

    if size(emf,1) > 1
        if (length(λ) == 1) | (length(θ) == 1)

            figure()
            ax_profile = gca()
            # draw incident medium
            ax_profile[:fill_between]([-doffset, 0], 0, nλ0[1], facecolor=βassigned[:,1], interpolate="true", alpha=0.8, edgecolor = "black")

            new_d = cumsum([0; d], dims=1)
            for idx = 1 : nd
                # draw successive layers
                ax_profile[:fill_between]([new_d[idx], new_d[idx+1]], 0, nλ0[idx+1], facecolor=βassigned[:,idx+1], interpolate="true", alpha=.6, edgecolor = "black")
            end
            # draw substrate medium
            ax_profile[:fill_between]([new_d[end], new_d[end]+doffset], 0, nλ0[end], facecolor=βassigned[:,end], interpolate="true", alpha=.6, edgecolor = "black")

            ax_profile[:tick_params](which="both", direction="in", pad=10, labelsize=22) # ticks offset
            ax_profile[:yaxis][:set_ticks_position]("both")
            ax_profile[:xaxis][:set_ticks_position]("both")
            # axis("tight")
            xlabel("Thickness profile [nm]")
            ylabel(L"\mathrm{Refractive index at }\,\, \lambda_0")

            # plot the emf for the reference wavelength
            aux1 = findmin(abs.(λ .- λ0))[2][1]
            plot(τ, emf[aux1, 1, :])
            title((@sprintf "EMF-highest resonance at %0.0f [nm]" λ[aux1]), y=1.02)

        end # if (length(λ) == 1) | (length(θ) == 1)
    end # size(emf,1) > 1

end # function nplot(...)
