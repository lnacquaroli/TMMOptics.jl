# ThinFilmOptics.jl
 Simulation and optimization of thin film stacks using optical matrices. Performs the calculation of reflectance, transmittance, absorptance, electric field and photonic band gap (for photonic crystals only) using the transfer matrix formalism.

# Usage:
    results = Spectra(λ, λ0, n, dflags, dinput, w, materials, θ, emfflag, h, pbgflag)

Input:
    λ: wavelength range [nm]
    λ0: central wavelength [nm]
    n: index profile of the multilayer sequence
    dflags: flag to indicate the type of dinput input (optical o physic)
    dinput: array with d, optical thickness is fractional λ0 or geometrical input in nm
    w: wave polarization (1-->p, 0-->s, 0<w<1 calculations are averaged)
    materials: materials inside the multilayer
    θ: angle of incidence [degrees]
    emfflag: flag that indicates if the calculation of electric field is required (1) or not (0). Default is 0.
    h: number of sublayers in which each layer will be divided to calculate the electric field. Default is 1.
    pbgflag: flag that indicates if photonic band gap must be calculated (1) or not (0). Only valid for periodic structures. Default is 0.
Output:
    results: structure with the following fields:
        Rp: p-wave complex reflectance
        Rs: s-wave complex reflectance
        R: w averaged reflectance
        Tp: p-wave complex transmittance
        Ts: s-wave complex transmittance
        T: w averaged transmittance
        ρp: p-wave complex reflection coefficient
        ρs: s-wave complex reflection coefficient
        τp: p-wave complex transmission coefficient
        τs: s-wave complex transmission coefficient
        emfp: p-wave electric field distribution
        emfs: s-wave electric field distribution
        emf: w averaged electric field distribution
        d: thickness of each layer [nm]
        κp: p-wave Bloch dispersion wavevectors
        κs: s-wave Bloch dispersion wavevectors
        κ: w averaged Bloch dispersion wavevectors
        ℓ: depth profile in the multilayer (taking into account the h for the emf plot) [nm]
        nλ0: refractive indexes profile at the central wavelength reflectance
        δ: phase shift
