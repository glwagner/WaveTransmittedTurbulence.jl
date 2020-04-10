# WaveTransmittedTurbulence

Tools, setups, run scripts, and visualization scripts for reproducing the simulation data and plots in 

> Wagner et al., "Near inertial waves and turbulence driven by the growth of surface waves", _submitted to JPO in April 2020_

## Simulations in section 3

### Spinup

To run a large eddy simulation of free convection, type

```
julia --project simulations/spinup/run_free_convection.jl
```

from the top-level of this repository.

The above command will run a simulation of free convection at a default low resolution of `32 x 32 x 32` and default surface buoyancy flux `1e-9 m^2 / s^3`, which will complete on some laptops in a matter of minutes.
The simulation will execute on a GPU if one is available.

To recreate the 'spin up' simulations reported in section 3, type

```
julia --project simulations/spinup/run_free_convection.jl --buoyancy_flux 1e-9 --Nh 256 --Nz 256
```

and

```
julia --project simulations/spinup/run_free_convection.jl --buoyancy_flux 1e-8 --Nh 256 --Nz 256
```

Additional command line arguments are revealed by typing

```
julia --project simulations/spinup/run_free_convection.jl --help
```

### Turbulence and near-inertial waves griven by gradual surface wave growth

## Simulations in section 4

# Repository structure
