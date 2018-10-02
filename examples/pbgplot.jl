"""
Plots the photonic dispersion of the multilayer structure processed
using PyPlot, PlotUtils and LaTeXStrings.
See examples folder for more details.
"""

using PyPlot, PlotUtils, LaTeXStrings
using Printf: @sprintf

function pbgplot(λ::AbstractArray{T1,N1}, θ::AbstractArray{T2,N2}, d::AbstractArray{T3,N3}, Rp::AbstractArray{T1,N4}, Rs::AbstractArray{T1,N4}, R::AbstractArray{T1,N4}, kblochp::AbstractArray{T5,N5}, kblochs::AbstractArray{T5,N5}, kbloch::AbstractArray{T5,N5}) where {T1<:Number, N1, T2<:Number, N2, T3<:Float64, N3, N4, T5<:ComplexF64, N5}

    nθ = lastindex(θ)

    # frequency range
    f = 2 * pi ./ λ

    # unit cell
    L = d[1] + d[2]

    if nθ == 1
        figure()
        subplot(2,2,1),plot(Rp,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (p-wave)")
        subplot(2,2,2),plot(real.(kblochp),f*L,linestyle="-",color="#66a61e")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{r}_{Bloch}\lambda\,\, \mathrm{(p-wave)}")
        subplot(2,2,3),plot(Rp,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (p-wave)")
        subplot(2,2,4),plot(abs.(imag.(kblochp)),f*L,linestyle="-",color="#e41a1c")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{i}_{Bloch}\lambda\,\, \mathrm{ (p-wave)}")
        figure()
        subplot(2,2,1),plot(Rs,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (s-wave)")
        subplot(2,2,2),plot(real.(kblochs),f*L,linestyle="-",color="#66a61e")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{r}_{Bloch}\lambda\,\, \mathrm{ (s-wave)}")
        subplot(2,2,3),plot(Rs,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (s-wave)")
        subplot(2,2,4),plot(abs.(imag.(kblochs)),f*L,linestyle="-",color="#e41a1c")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{i}_{Bloch}\lambda\,\, \mathrm{ (s-wave)}")
    else
        figure()
        contourf(f.*L,rad2deg.(θ),copy(real.(kblochp)'),40)
        xlabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        ylabel("Angle of incidence")
        title(L"K^{r}_{Bloch}\,\, \mathrm{(p-wave)}")
        figure()
        contourf(f.*L,rad2deg.(θ),copy(abs.(imag.(kblochp))'),40)
        xlabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        ylabel("Angle of incidence")
        title(L"K^{i}_{Bloch}\,\, \mathrm{(p-wave)}")
        # shading flat
        figure()
        contourf(f.*L,rad2deg.(θ),copy(real.(kblochs)'),40)
        xlabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        ylabel("Angle of incidence")
        title(L"K^{r}_{Bloch}\,\, \mathrm{(s-wave)}")
        # shading flat
        figure()
        contourf(f.*L,rad2deg.(θ),copy(abs.(imag.(kblochs))'),40)
        xlabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        ylabel("Angle of incidence")
        title(L"K^{i}_{Bloch}\,\, \mathrm{(s-wave)}")
    end

end # function pbgplot(...)
