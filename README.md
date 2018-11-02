# TMMOptics.jl

[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](http://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/lnacquaroli/TMMOptics.jl.svg?branch=master)](https://travis-ci.com/lnacquaroli/TMMOptics.jl)

Simulation of the propagation of plane waves through arbitrary layered system of thin film stacks using optical matrices. Performs the calculation of reflectance, transmittance, absorptance, electromagnetic field and photonic band gap (for photonic crystals only) using the transfer matrix formalism. For further details see [https://arxiv.org/abs/1809.07708](https://arxiv.org/abs/1809.07708).

## Installation

This package is not yet registered. It can be installed in Julia with the following ([see further](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html#Adding-unregistered-packages-1)):
```julia
julia> ]
(v1.0) pkg> add https://github.com/lnacquaroli/TMMOptics.jl
```

`TMMOptics.jl` is compatible with Julia version 1.0 or later.

If you want to avoid reading any further you can jump to the examples posted inside the examples folder.

## Usage

The typical calling structure is as follow:

```julia
results = ThinFilmOptics(beam, nseq, emfflag=0, h=1, pbgflag=0)
```

Where `ThinFilmOptics` is the main function called from the module `TMMOptics.jl`.

## Input arguments

### LightSource type

`beam::PlaneWave` is a subtype that defines the type of light source used in the simulation. At the moment the only subtype accepted is `PlaneWave <: LightSource`, which is exported by `TMMOptics.jl`. To call this subtype we can write `beam = PlaneWave(λ, λ0, θ, p)`, where the parameters required are described below.

#### Wavelength range

`λ` defines the wavelength range in nanometers, over which the calculations are performed. The recommended input are 1-d arrays, linear or step ranges (although not explicitly defined). For instance:
```julia
λ = LinRange(200, 1000, 801) # this ensures floating input
λ = 200:1000
λ = [245. 300. 450. 868.]
λ = [632.]
```

#### Reference wavelength

`λ0::Float64` sets the reference wavelength in nanometers, especially useful when computing mutlilayer stacks. It is used to compute the geometrical (physical) thickness of each layer and also to plot the profile of index of refraction allowing to visualize the optical structure processed. Accepts floating numbers. For instance:
```julia
λ0 = 632.
```

#### Angle of incidence

`θ` defines the angle of incidence (in degrees) of the wave into the first (incident) medium. It accepts the same type of input as `λ`. For instance:
```julia
θ = LinRange(0, 90, 901) # this ensures floating input
θ = [0.] # zero degress --> normal incidence
θ = [45.] # 45. degrees
θ = 0:1:60
```

By setting `θ` and `λ` as arrays with more than one element each, the quantities that depends on them will be output as 2d-array or 3d-array.

#### Polarization of the wave

`p::Float64` indicates the polarization (linear polarization) of the wave. Accepted values are `p = 1` for the p-wave (TM), and `p = 0` for the s-wave (TE). In fact, the program computes both types at every run, however this parameter set the admittances calculations. The relevance of this parameter is subject to `θ != 0`, otherwise, there is no polarization effect. Floating values are recommended for this input. E.g.:
```julia
p = 1. # p-wave
p = 0. # s-wave
```

### Materials type

The `Material` supertype allows two subtypes exported by `TMMOptics.jl`. The first one is `Geometrical(n, d)`, which is used to indicate an input layer with its complex index of refraction `n`and geometrical (physical) thickness `d` in nanometers. The other subtype allowed is `Optical(n, f)`, which is used to construct a layer with refractive `n` and the fraction of optical thickness at `beam.λ0`. The last one is useful for multilayer stacks normally written in terms of quarter/half-wave layers. Both subtypes can be mixed when defining the structure.

`nseq::Array{Materials}` is an array setting the sequence of the materials involved in the simulation, containing the information about the index of refraction and the thickness and type of thickness that compose the multilayer stack. The order set here will be used to alternate the layers in the multilayer.

The first element in `nseq` indicates the incident medium while the last one the substrate. 

#### Geometrical subtype

For instance, in the case of simulation a single layer configuration:
```julia
l0 = Geometrical(air(beam.λ), 0.) # incident medium
l1 = Geometrical(sno2f(beam.λ), 150.)
l2 = Geometrical(silicon(beam.λ), 0.) # subtrate
nseq = [l0 l1 l2]
```

#### Optical subtype

For a multilayer stack ([Distributed Bragg Reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector)) alternating two different index of refraction materials 1 and 2, wrapped inside an incident medium 1 and a substrate 4, we can write:
```julia
l0 = Geometrical(air(beam.λ), 0.)
l1 = Optical(looyengaspheres(air(beam.λ),silicon(beam.λ),0.54), 1/4.)
l2 = Optical(looyengaspheres(air(beam.λ),silicon(beam.λ),0.86), 1/4.)
l3 = Geometrical(silicon(beam.λ), 0.)
nseq = [l0 l1 l2 l1 l2 l1 l2 l1 l2 l1 l2 l1 l2 l1 l2 l1 l2 l3]
```

where `l1` and `l2` are quarter-wave optical thickness at the reference wavelength `λ0`.

For all cases, the thicknesses of the first and last layers are not taken into account as the calculations are performed for semi-infinite incident/substrate media, so we just put a zero regardless of the subtype.

### Electromagnetic field (EMF) flag

`emfflag::Int64` is a flag that indicates whether to compute the field distribution inside the multilayer or not. If `emfflag = 1`, the elecromagnetic field is computed, while if `emfflag = 0` is not. This is useful especially when you just want the reflectance or transmission spectra information but not the field distribution, because the computation might be time consuming. The default value is `0`.

### Number of layers division for EMF calculation

`h::Int64` set the number of sub-layers the script will use to divide each layer defined in `nseq` to calculate the EMF. If you opted for `emfflag = 1` then you might want to set this to, say `10` to better visualize the field distribution. You can think as to have more numerical resolution when increasing `h`. The default value is `1`. If you choose `emfflag = 0`, then `h` is neglected.

### Photonic band gap flag

`pbgflag::Int64` indicates whether to calculate or not the photonic dispersion. It is only useful when computing DBR, i.e., when computing binary periodic systems that alternate two different refractive indexes. By default `pbgflag = 0`.

## Output

`results::Results <: Output` is a structure containing different subtypes of outputs depending on their relation (more or less). There are five subtypes included in the output:

### Spectra subtype of output

`results.Spectra::Results <: Output` contains information on the reflectance, transmittance and fresnell coefficients spectra. Inside this subtype there exist the following fields:

`Rp::Array{Float64}(lastindex(beam.λ), lastindex(beam.θ))`: p-wave reflectance.

`Rs::Array{Float64}(lastindex(beam.λ), lastindex(beam.θ))`: s-wave reflectance.

`Tp::Array{Float64}(lastindex(beam.λ), lastindex(beam.θ))`: p-wave transmittance.
  
`Ts::Array{Float64}(lastindex(beam.λ), lastindex(beam.θ))`: s-wave transmittance.
 
`ρp::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ))`: p-wave complex reflection coefficient.
  
`ρs::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ))`: s-wave complex reflection coefficient.
  
`τp::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ))`: p-wave complex transmission coefficient.
  
`τs::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ))`: s-wave complex transmission coefficient.

### Field subtype of output

`results.Field::Results <: Output` contains information on the electromagnetic field distribution calculated, with the following fields:

`emfp::Array{Float64}(lastindex(beam.λ), lastindex(beam.θ), lastindex(results.Misc.ℓ))`: p-wave electric field distribution.
  
`emfs::Array{Float64}(lastindex(beam.λ), lastindex(beam.θ), lastindex(results.Misc.ℓ))`: s-wave electric field distribution.

### Bloch subtype of output

`results.Bloch::Results <: Output` contains information on the Bloch wavevectors calculated for DBRs, with the following fields:

`κp::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ))`: p-wave Bloch dispersion wavevectors.

`κs::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ))`: s-wave Bloch dispersion wavevectors.

### Misc subtype of output

`results.Misc::Results <: Output` gives information about the geometrical (physical) thicknesses in nanometers (even if you input the optical thickness information, the output returns the physical one), the geometrical thickness taking into account the `h` number defined above, and the index of refraction of each layer defined in `nseq` at the wavelength of reference `λ0`. The fields are as follow:

`d::Array{Float64}(size(nseq,2)-2)`: geometrical (physical) thickness of each layer as ordered by `n` [nm].

`ℓ::Array{Float64}((size(nseq,2)-2)*h)`: Geometrical (physical) thickness of each layer, where each layer is divided into `h` sub-layers [nm].

`nλ0::Array{Float64}(size(nseq,2))`: profile of index of refraction profile computed at the wavelength reference `λ0`.

### AdmPhase subtype of output

`results.AdmPhase::Results <: Output` wraps the admittances for both polarizations and the phase shift angle (thickness) for the whole structure. The fields are as follow:

`ηp::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ), lastindex(size(nseq,2)))`: p-wave admittance of the structure.

`ηs::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ), lastindex(size(nseq,2)))`: s-wave admittance of the structure.

`δ::Array{ComplexF64}(lastindex(beam.λ), lastindex(beam.θ), lastindex(size(nseq,2))-2)`: phase shift (thickness) of the structure without the incident and substrate media.

## Examples of index of refraction functions used

The next two modules are optional since include functions with index of refraction information compiled for the purpose of the docs. They can also be used as provided for problems involving the same materials though.

#### RIdb.jl

Module containing a collection of functions with index of refration for the following materials: aluminum, air, bk7, chrome, dummy, glass, gold, silicon, silicontemperature, silver, sno2f, h2o, etoh. These functions accept as input arguments the wavelength range `λ`, and return the index of refraction as a complex floating number. Even for non-aborbent materials (where a list of zeros in the imaginary part is placed), the index works better with complex character. Users can use their own functions as well that output complex types though.

This module depends on [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) to return the index of refraction as a function of the input wavelength, which might differ from that of the experimental data. Also depends on [HDF5.jl](https://github.com/JuliaIO/HDF5.jl) since the data is compiled into a h5 file for simplicity. 

Use of this module as follow, e.g.:
```julia
include("pathToFile/RIdb.jl")
using Main.RIdb: aluminum, air, bk7, chrome, dummy, glass, gold, silicon, silicontemperature, silver, sno2f, h2o, etoh
λ = 200:1000.
n = silicon(λ)
t = 240 # temperature in C
n = silicontemperature(λ, t)
n = sno2f(λ)
```

Data taken from database: https://refractiveindex.info, http://www.ioffe.ru/SVA/NSM/nk/, https://github.com/ulfgri/nk. You have to check the range of validity for the wavelength of each index.

#### MixingRules.jl

Module containing a collection of functions with effective index of refraction mixing rules for binary systems, accepting two indexes of refraction with the same lengths and the volume fraction of one of them.

```julia
include("pathToFile/RIdb.jl")
using Main.RIdb: aluminum, air, bk7, chrome, dummy, glass, gold, silicon, silicontemperature, silver, sno2f, h2o, etoh
include("pathToFile/MixingRules.jl")
using Main.MixingRules: bruggemanspheres, looyengacylinders, looyengaspheres, lorentzlorenz, maxwellgarnettspheres, monecke, gedf, gem
λ = 200:1000.
f = 0.86 # volume fraction of air inside the host material
l1 = looyengaspheres(air(λ), silicon(λ), f)
l2 = looyengacylinder(air(λ), silicon(λ), f)
l3 = lorentzlorens(h2o(λ), etoh(λ), f)
l4 = maxwellgarnettspheres(air(λ), gold(λ), f)
```

Sources: [Theiss (1997)](https://www.sciencedirect.com/science/article/pii/S016757299600012X), [Torres-Costa et al. (2014)](https://www.sciencedirect.com/science/article/pii/B9780857097118500088), [Liu (2016](https://doi.org/10.1063/1.4943639), [Celzard et al. (2002)](https://doi.org/10.1016/S0008-6223(02)00196-3), [Urteaga et al. (2013)](https://dx.doi.org/10.1021/la304869y). See inside the module for more especific details.

## Examples of useful plotting functions

Included two types of plots so far are quite useful to visualize the type of structure being modeled. These are out of the main module as they can be customized by the user.

#### nplot.jl (using PyPlot)

Example of plotting the index of refraction profile selected in `n`, usually for `λ0` (but could be different). This is particularly advantageous to visualize the chosen stack and track possible mistakes. For instance:
```julia
nplot(beam.λ, beam.θ, beam.λ0, results.Misc.d, results.Misc.ℓ, results.Field.emfp, results.Misc.nλ0, nseq)
```

If `size(results.emfp,1) > 1` prompts another plot with the EMF at `λ0` overlapping the structure.

#### pbgplot.jl (using PyPlot)

Example of plotting the photonic dispersion of DBRs. For instance:
```julia
pbgplot(beam.λ, beam.θ, results.Misc.d, results.Spectra.Rp, results.Spectra.Rs, results.Bloch.κp, results.Bloch.κs)
```

## We welcome suggestions

If you have ideas and suggestions to improve TMMOptics in Julia, PRs and issues are welcomed.

## Other projects in optics (not necessarily in Julia)

* [https://github.com/ngedwin98/ABCDBeamTrace.jl](https://github.com/ngedwin98/ABCDBeamTrace.jl)

* [https://github.com/lbolla/EMpy](https://github.com/lbolla/EMpy)

* [https://github.com/ederag/geoptics](https://github.com/ederag/geoptics)

* [https://ricktu288.github.io/ray-optics/simulator/](https://ricktu288.github.io/ray-optics/simulator/)

* [https://github.com/brownjm/geometric-optics](https://github.com/brownjm/geometric-optics)

* [https://github.com/NGC604/thinerator](https://github.com/NGC604/thinerator)

* [https://github.com/marcus-cmc/Optical-Modeling](https://github.com/marcus-cmc/Optical-Modeling)

* [http://opticspy.org/](http://opticspy.org/)

* [https://github.com/Bitburg-chef/WavefrontOptics](https://github.com/Bitburg-chef/WavefrontOptics)

* [https://github.com/jlcodona/AOSim2](https://github.com/jlcodona/AOSim2)

* [https://github.com/stevengj/meep](https://github.com/stevengj/meep)

* [https://github.com/stevengj/mpb](https://github.com/stevengj/mpb)

## To do

* Rugate filters

