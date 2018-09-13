"""
Function that plots the multilayer structure, using PyPlot, PlotUtils and LaTeXStrings. Also using Printf: @sprintf to display some data.
author: lnacquaroli
"""

using PyPlot, PlotUtils, LaTeXStrings
using Printf: @sprintf

function drawindexprofile(λ::AbstractArray{T,M}, nprofile::AbstractArray{U,N}, dprofile::AbstractArray{V,O}, nlayers::AbstractArray{W,P}, emf::AbstractArray{X,Q}, multilayerdepth::AbstractArray{Y,R}, θ::AbstractArray{Z,S}, λ0::A2) where {T<:Number, M, U<:Number, N, V<:Number, O, W<:Number, P, X<:Number, Q, Y<:Number, R, Z<:Number, S, A1<:Number, A2<:Number}

    # check input
    numelnlayers = lastindex(nlayers)
    numeldprofile = lastindex(dprofile)
    # if size(dprofile,1) != size(nprofile,1)-2
    @assert numeldprofile == lastindex(nprofile)-2 "Variables dimensions disagree: length(nprofile)-2 = length(dprofile)"

    # define the offset for the incident and substrate medium
    doffset = 200. # nm

    # Maps each certain to a certain color
    aux2 = unique(nlayers)
    numelaux2 = lastindex(aux2)
    cm = cgrad(:Spectral)
    colores = [cm[i] for i in LinRange(0,1,numelaux2)] # range(0,stop=1,length=numelaux2) should work also
    aux3 = [getfield(colores[j],i) for i=1:4, j=1:lastindex(colores)]

    assigned_color = Array{Float64,2}(undef, (4,numelnlayers))
    for idx2 = 1 : numelaux2, idx1 = 1 : numelnlayers
        if aux2[idx2] == nlayers[idx1]
            assigned_color[:,idx1] = aux3[:,idx2]
        end
    end

    nλ0 = real.(nprofile)
    figure()
    ax_profile = gca()
    # draw incident medium
    ax_profile[:fill_between]([-doffset, 0], 0, nλ0[1], facecolor=assigned_color[:,1], interpolate="true", alpha=0.8, edgecolor = "black")

    new_dprofile = cumsum([0; dprofile], dims=1)
    for idx = 1 : numeldprofile
        # draw successive layers
        ax_profile[:fill_between]([new_dprofile[idx], new_dprofile[idx+1]], 0, nλ0[idx+1], facecolor=assigned_color[:,idx+1], interpolate="true", alpha=.6, edgecolor = "black")
    end
    # draw substrate medium
    ax_profile[:fill_between]([new_dprofile[end], new_dprofile[end]+doffset], 0, nλ0[end], facecolor=assigned_color[:,end], interpolate="true", alpha=.6, edgecolor = "black")

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
            ax_profile[:fill_between]([-doffset, 0], 0, nλ0[1], facecolor=assigned_color[:,1], interpolate="true", alpha=0.8, edgecolor = "black")

            new_dprofile = cumsum([0; dprofile], dims=1)
            for idx = 1 : numeldprofile
                # draw successive layers
                ax_profile[:fill_between]([new_dprofile[idx], new_dprofile[idx+1]], 0, nλ0[idx+1], facecolor=assigned_color[:,idx+1], interpolate="true", alpha=.6, edgecolor = "black")
            end
            # draw substrate medium
            ax_profile[:fill_between]([new_dprofile[end], new_dprofile[end]+doffset], 0, nλ0[end], facecolor=assigned_color[:,end], interpolate="true", alpha=.6, edgecolor = "black")

            ax_profile[:tick_params](which="both", direction="in", pad=10, labelsize=22) # ticks offset
            ax_profile[:yaxis][:set_ticks_position]("both")
            ax_profile[:xaxis][:set_ticks_position]("both")
            # axis("tight")
            xlabel("Thickness profile [nm]")
            ylabel(L"\mathrm{Refractive index at }\,\, \lambda_0")

            # plot the emf for the reference wavelength
            aux1 = findmin(abs.(λ .- λ0))[2][1]
            plot(multilayerdepth, emf[aux1, 1, :])
            title((@sprintf "EMF-highest resonance at %0.0f [nm]" λ[aux1]), y=1.02)

        end # if (length(λ) == 1) | (length(θ) == 1)
    end # size(emf,1) > 1

end # function drawindexprofile(...)
