# ThinFilmOptics.jl

Simulation of the propagation of plane waves through arbitrary layered system of thin film stacks using optical matrices. Performs the calculation of reflectance, transmittance, absorptance, electric field and photonic band gap (for photonic crystals only) using the transfer matrix formalism.

ThinFilmOptics.jl is compatible with Julia version 1.0 or later.

If you want to avoid reading any further you can jump to the examples posted inside the examples folder.

## Usage

The typical calling structure is as follow:

```julia
results = Spectra(λ, λ0, n, dflags, dinput, w, materials, θ, emfflag=0, h=1, pbgflag=0)
```

## Input

#### Spectra

`Spectra` is the function called fom the main module `ThinFilmOptics.jl`. (See also `results`)

#### Wavelength range

`λ` is the wavelenth range in nanometers. The accepted input values can be 1-d array or linear ranges. For instance:
```julia
λ = 200:1000
λ = [245. 300. 450. 868.]
λ = [632.]
```
  
#### Reference wavelength

`λ0` is the reference wavelength in nanometers, especially useful when computing mutlilayer stacks. It is used to plot the profile of index of refraction allowing to visualize the optical structure processed. Accepts floating numbers. For instance:
```julia
λ0 = 632.
```

#### Layers profile

`n` is an 1d-array setting the sequence of index of refraction that composes the multilayer stack. The order set here will be used to alternate the layers of media selected below in materials. The first element indicates the incident medium while the last one the substrate. For instance, in the case of simulation a single layer configuration:
```julia
n = [1 2 3]
```
which indicates that the incident medium is layer number 1 and the substrate is layer 3. The second layer, number 2, is the layer of interest. For a multilayer stack ([Distributed Bragg Reflector](https://en.wikipedia.org/wiki/Distributed_Bragg_reflector)) alternating two different index of refraction materials 1 and 2, wrapped inside an incident medium 1 and a substrate 4, we can write:
```julia
n = [1 2 3 2 3 2 3 2 3 2 3 2 3 4]
```

#### Thickness profile

`d` is an 1d-array setting the thickness of each layer inside the multilayer stack, which accepts two types of input: the geometrical (physical) thickness in nanometers, or the fraction of optical thickness (product of the index of refraction at the wavelength of reference times the geometrical thickness) at the wavelength of reference `λ0`. The latter is useful when computing multilayer stacks and normally using quarter-wavelength or half-wavelength thicknesses. The elements of this array contains information about all the layers inside the multilayer, except the first and last ones (semi-infinite calculation). For instance, for the single layer example, we have
```julia
d = [648.]
```

while for the DBR, we can set the whole stack to have quarter-wavelength optical thickness as follow:
```julia
d = 1/4 * ones.(lastindex(n)-2,1)
```

It is also accepted to mixed both type of input here. We can set for instance:
```julia
d = [1/4 250. 1/2 980.]
```

for a four layer system. Remember that `d` does not include information about the incident and substrate.

#### Type of thickness profile

`dflags` tells the script to treat `d` as optical or geometrical thickness. This gives more flexibility to the computation and also allows to mix the `d` inputs. If the thickness input element is optical `dflags = 1`, or if it is geometrical `dflags = 0`. The number of elements of `dflags` must equal that of `d`. For instance, for the single layer example above:
```julia
dflags = [0.]
```

while for the DBR:
```julia
dflags = ones.(lastindex(n)-2,1)
```
or
```julia
dflags = ones.(lastindex(d),1)
```

In the case of mixing flags input, we can write:
```julia
dflags = [1 0 1 0]
```

#### Polarization of the wave

`w` indicates the polarization (linear polarization) of the wave. Accepted values are `w = 1` for the p-wave (TM), and `w = 0` for the s-wave (TE). An option to calculate quantities with intermediate values, `0 < w < 1` is accepted as well. The relevance of this parameter is subject to `θ != 0`, otherwise, there is no polarization effect. Floating values are recommended for this input. E.g.:
```julia
w = 1. # p-wave
w = 0. # s-wave
w = 0.5
w = 0.4 # mixes 60% of s-wave polarized light with 40% p-wave.
```

#### Index of refraction data

`materials` contains the index of refraction of the media used in the stack, while `n` keeps all the information about how these media will be placed in the whole structure. The idea is that we feed the script with the index if refraction as a function of `λ`. The `RIdb.jl` module holds a custom database with experimental data of index of refraction for different materials. The order of the index of refraction input in this variable are directly related to the element number inside the `n` profile. It is possible to create and input your own function with an index of refraction of your interest. For example in the single layer case we have:
```julia
l1 = air(λ) # outermost medium
l2 = chrome(λ) # active layer
l3 = silicon(λ) # substrate medium
materials = [l1 l2 l3]
```

where the first material (`air`, `l1`) is the incident medium and the third material (`silicon`, `l3`) is the substrate. The second material is the active layer `l2`. To understand the relation with `n`, `l1` will be placed in the structure where the number `1`is located inside `n`, while `l2` will be placed at `2`, and `l3` at `3`. When it comes to a DBR, alternating two indexes of refraction, we can write as follow:
```julia
l1 = air(λ) # outermost medium
l2 = looyengaspheres(air(λ),silicon(λ),0.86)
l3 = looyengaspheres(air(λ),silicon(λ),0.54)
l4 = silicon(λ) # substrate medium
materials = [l1 l2 l3 l4]
```

where the incident material is air and the substrate is silicon. The number `2` is a medium with index of refractive index in terms of Looyenga mixing rule with 0.86 porosity, and `3` with 0.54 porosity mixing air and silicon. Note here that `l2` will be placed inside the structure where there is a number `2` in `n`, `l3` will be placed where there is a `3` in `n`. This way of writing `n` and `materials` makes suitable to make cleaner input, especially when you want to compute multilayer stacks with lots of layers.

The number of materials included here should be at least equal to the maximum element in `n`. You can put more if you know where to place them, but not less `maximum(n) == size(materials,2)`.

#### Angle of incidence

`θ` is an 1d-array that set the angle of incidence of the wave into the first (incident) medium. It accepts the same type of input as `λ`, in degrees. For instance:
```julia
θ = [0.] # zero degress --> normal incidence
θ = [45.] # 45. degrees
θ = 0:1:60
```

By setting `θ` and `λ` as arrays with more than one element each, the quantities that depends on them will be output as 2d-array or 3d-array.

#### Electromagnetic field (EMF) flag

`emfflag` is a flag that indicates whether to compute the field distribution inside the multilayer or not. If `emfflag = 1`, the elecromagnetic field is computed, while if `emfflag = 0` is not. This is useful especially when you just want the reflectance or transmission spectra information but not the field distribution, because the computation might be time consuming. The default value is `0`.

#### Number of layers division

`h` is a floating number that tells in how many sub-layers the script will use to calculate the EMF. If you opted for `emfflag = 1` then you might want to set this to, say `10` to better visualize the field distribution. You can think as to have more numerical resolution. The default value is `1`. If you choose `emfflag = 0`, then `h` is neglected.

#### Photonic band gap flag

`pbgflag` indicates whether to calculate or not the photonic dispersion. It is only useful when computing DBR, i.e., when computing binary periodic systems that alternate two different refractive indexes. By default `pbgflag = 0`.

## Output

`results` is a structure with type `Spectra` that contains the following fields computed by the main program:

  `results.Rp`: p-wave complex reflectance. 2d-array with lengths of `λ` and `θ`
  
  `results.Rs`: s-wave complex reflectance. 2d-array with lengths of `λ` and `θ`
  
  `results.R`: w averaged reflectance. 2d-array with lengths of `λ` and `θ`
  
  `results.Tp`: p-wave complex transmittance. 2d-array with lengths of `λ` and `θ`
  
  `results.Ts`: s-wave complex transmittance. 2d-array with lengths of `λ` and `θ`
  
  `results.T`: w averaged transmittance. 2d-array with lengths of `λ` and `θ`
  
  `results.ρp`: p-wave complex reflection coefficient. 2d-array with lengths of `λ` and `θ`
  
  `results.ρs`: s-wave complex reflection coefficient. 2d-array with lengths of `λ` and `θ`
  
  `results.τp`: p-wave complex transmission coefficient. 2d-array with lengths of `λ` and `θ`
  
  `results.τs`: s-wave complex transmission coefficient. 2d-array with lengths of `λ` and `θ`
  
  `results.emfp`: p-wave electric field distribution. 3d-array with lengths of `λ`, `θ`, and `ℓ`
  
  `results.emfs`: s-wave electric field distribution. 3d-array with lengths of `λ`, `θ`, and `ℓ`
  
  `results.emf`: w averaged electric field distribution. 3d-array with lengths of `λ`, `θ`, and `ℓ`
  
  `results.d`: thickness of each layer [nm]. Geometrical (physical) thickness of each layer as ordered by `n`. 1d-array with length of `d`
  
  `results.κp`: p-wave Bloch dispersion wavevectors. 2d-array with lengths of `λ` and `θ`
  
  `results.κs`: s-wave Bloch dispersion wavevectors. 2d-array with lengths of `λ` and `θ`
  
  `results.κ`: w averaged Bloch dispersion wavevectors. 2d-array with lengths of `λ` and `θ`
  
  `results.ℓ`: depth profile in the multilayer stack taking into account `h` [nm]. Geometrical (physical) thickness of each layer, where each layer is divided into `h` sub-layers. 1d-array with length of `lastindex(d) * h`
  
  `results.nλ0`: profile of index of refraction profile computed at the wavelength reference `λ0`. 1d-array with length of `n`
  
  `results.δ`: phase shift of the whole structure. 3d-array with lengths `λ`, `θ`, and `d`


## Welcome ideas

If you have any ideas and suggestions to improve this in Julia, PRs and issues are welcomed.


## Other projects in optics (not necessarily in Julia)

[https://github.com/ngedwin98/ABCDBeamTrace.jl](https://github.com/ngedwin98/ABCDBeamTrace.jl)

[https://github.com/lbolla/EMpy](https://github.com/lbolla/EMpy)

[https://github.com/ederag/geoptics](https://github.com/ederag/geoptics)

[https://ricktu288.github.io/ray-optics/simulator/](https://ricktu288.github.io/ray-optics/simulator/)

[https://github.com/brownjm/geometric-optics](https://github.com/brownjm/geometric-optics)

[https://github.com/NGC604/thinerator](https://github.com/NGC604/thinerator)

[https://github.com/marcus-cmc/Optical-Modeling](https://github.com/marcus-cmc/Optical-Modeling)

[http://opticspy.org/](http://opticspy.org/)

[https://github.com/Bitburg-chef/WavefrontOptics](https://github.com/Bitburg-chef/WavefrontOptics)

[https://github.com/jlcodona/AOSim2](https://github.com/jlcodona/AOSim2)

[https://github.com/stevengj/meep](https://github.com/stevengj/meep)
