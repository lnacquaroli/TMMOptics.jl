# ThinFilmOptics.jl
 Simulation and optimization of thin film stacks using optical matrices. Performs the calculation of reflectance, transmittance, absorptance, electric field and photonic band gap (for photonic crystals only) using the transfer matrix formalism.

## Usage
    results = Spectra(λ, λ0, n, dflags, dinput, w, materials, θ, emfflag, h, pbgflag)

`Spectra` is the function called fom the main module `ThinFilmOptics.jl`. (See `results` also)

`λ` is the wavelenth range in nanometers. The accepted input values can be 1-d array or linear ranges. For instance:

    λ = 200:1000
    λ = [245. 300. 450. 868.]
    λ = [632.]
  
`λ0` is the reference wavelength in nanometers, especially useful when computing mutlilayer stacks. It is used to plot the profile of index of refraction allowing to visualize the optical structure processed. Accepts floating numbers. For instance:

    λ0 = 632.    

In progress...
