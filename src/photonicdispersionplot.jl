"""
Example of function that plots the photonic dispersion, using PyPlot, PlotUtils and LaTeXStrings.
author: lnacquaroli
"""

using PyPlot, PlotUtils, LaTeXStrings
using Printf: @sprintf

function photonicdispersionplot(λ::AbstractArray{T,M}, θ::AbstractArray{U,N}, d::AbstractArray{W,P}, Rp::AbstractArray{X,Q}, Rs::AbstractArray{X,Q}, R::AbstractArray{X,Q}, kblochp::AbstractArray{A1,B1}, kblochs::AbstractArray{A1,B1}, kbloch::AbstractArray{A1,B1}) where {T<:Number, M, U<:Number, N, W<:Number, P, X<:Number, Q, A1<:ComplexF64, B1}

    numbertheta = lastindex(θ)

    # frequency range
    f = 2 * pi ./ λ

    # unit cell
    L = d[1] + d[2]

    if numbertheta == 1
        figure()
        subplot(2,2,2),plot(real.(kblochp),f*L,linestyle="-",color="#66a61e")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{r}_{Bloch}\lambda\,\, \mathrm{(p-wave)}")
        subplot(2,2,1),plot(Rp,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (p-wave)")
        subplot(2,2,4),plot(abs.(imag.(kblochp)),f*L,linestyle="-",color="#e41a1c")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{i}_{Bloch}\lambda\,\, \mathrm{ (p-wave)}")
        subplot(2,2,3),plot(Rp,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (p-wave)")
        figure()
        subplot(2,2,2),plot(real.(kblochs),f*L,linestyle="-",color="#66a61e")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{r}_{Bloch}\lambda\,\, \mathrm{ (s-wave)}")
        subplot(2,2,1),plot(Rs,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (s-wave)")
        subplot(2,2,4),plot(abs.(imag.(kblochs)),f*L,linestyle="-",color="#e41a1c")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel(L"K^{i}_{Bloch}\lambda\,\, \mathrm{ (s-wave)}")
        subplot(2,2,3),plot(Rs,f*L,linestyle="-",color="#386cb0")
        ylabel(L"\omega\lambda = (2\pi /\lambda)(d_{1}+d_{2})")
        xlabel("Reflectance (s-wave)")
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

end # function photonicdispersionplot(...)
