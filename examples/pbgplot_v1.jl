"""
Plots the photonic dispersion of the multilayer structure processed using PyPlot, PlotUtils and LaTeXStrings.
See examples folder for more details.
"""

using PyPlot, PlotUtils, LaTeXStrings
using Printf: @sprintf

function pbgplot(λ::AbstractArray{T1,N1}, θ::AbstractArray{T2,N2}, d::AbstractArray{T3,N3}, Rp::AbstractArray{T1,N4}, Rs::AbstractArray{T1,N4}, kblochp::AbstractArray{T5,N5}, kblochs::AbstractArray{T5,N5}) where {T1<:Number, N1, T2<:Number, N2, T3<:Float64, N3, N4, T5<:ComplexF64, N5}

    nθ = lastindex(θ)

    # periodicity
    Λ = d[1] + d[2]

    # frequency range normalized
    ω = Λ ./ λ

    # kbloch for one angle only
    temp0 = Λ ./ π
    kblochs = kblochs .* temp0
    kblochp = kblochp .* temp0

    # k normalized for angle-frequency dependence
    # logical matrices, used to select points which belong to the forbidden bands
    # maskp = abs.(kblochp) .>= 1.
    # masks = abs.(kblochs) .>= 1.
    maskp = imag.(kblochp) .== 0.
    masks = imag.(kblochs) .== 0.
    κp = deepcopy(kblochp)
    κp[.!maskp] .= 1.
    κp[maskp] .= 0.
    κs = deepcopy(kblochs)
    κs[.!masks] .= 1.
    κs[masks] .= 0.
    # parallel wavevector qz
    qz = sin.(deg2rad.(θ))
    # choose better colormap
    cmp = ColorMap(get_cmap("Blues"))
    cms = ColorMap(get_cmap("Reds"))

    if nθ == 1
        figure()
        title("p/TM-wave")
        subplot(2,2,1),plot(Rp, ω, linestyle="-",color="#386cb0")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel("Reflectance")
        subplot(2,2,2),plot(real.(kblochp), ω, linestyle="-",color="#66a61e")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel(L"K^{r}_{Bloch}\Lambda/\pi")
        subplot(2,2,3),plot(Rp, ω,linestyle="-",color="#386cb0")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel("Reflectance")
        subplot(2,2,4),plot(abs.(imag.(kblochp)), ω,linestyle="-",color="#e41a1c")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel(L"K^{i}_{Bloch}\Lambda/\pi")
        figure()
        title("s/TE-wave")
        subplot(2,2,1),plot(Rs, ω,linestyle="-",color="#386cb0")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel("Reflectance")
        subplot(2,2,2),plot(real.(kblochs), ω, linestyle="-",color="#66a61e")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel(L"K^{r}_{Bloch}\Lambda/\pi")
        subplot(2,2,3),plot(Rs, ω, linestyle="-",color="#386cb0")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel("Reflectance")
        subplot(2,2,4),plot(abs.(imag.(kblochs)), ω, linestyle="-",color="#e41a1c")
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel(L"K^{i}_{Bloch}\Lambda/\pi")
    else
        figure()
        title("p/TM-wave")
        contourf(qz, ω, κp, 60, cmap = cmp, alpha = 0.7)
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel(L"K_{Bloch}\Lambda/\pi")
        # shading flat
        figure()
        title("s/TE-wave")
        contourf(qz, ω, κs, 60, cmap = cms, alpha = 0.7)
        ylabel(L"\omega\Lambda/(2\pi)")
        xlabel(L"K_{Bloch}\Lambda/\pi")
    end

end # function pbgplot(...)
