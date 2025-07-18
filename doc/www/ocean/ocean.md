(ref_ocean_overview)=
# Ocean Modelling Overview

The ocean component of ICON consists of three parts:

1. [The ocean model](ref_ocean_model) itself (referred to as ICON-O)
2. [The sea-ice module](ref_ocean_seaice)
3. and the [bio-geochemistry module HAMOCC](ref_ocean_biogeochem).

Additionally, an [ocean skin parameterisation](ref_ocean_skin) is available.

## Ocean model components

(ref_ocean_model)=
### Ocean model
The ocean model ICON-O ({term}`Korn 2018`; {term}`Korn et al. 2022`) solves the hydrostatic Boussinesq equations, the classical set of dynamical equations for global ocean dynamics.
ICON-O uses, as the atmosphere, an icosahedral–triangular C grid.
The global grid of ICON-O can be locally refined to create a “computational telescope” that zooms into a region of interest.
The numerics of ICON-O share similarities to the atmosphere component but also have important differences.
Both components use a mimetic discretization of discrete differential operators, but ICON-O uses the novel concept of Hilbert-space compatible reconstructions to calculate volume and tracer fluxes on the staggered ICON grid (see {term}`Korn 2018`).
The numerical scheme of ICON-O allows for generalized vertical coordinates: of particular importance here are the depth-based z-level coordinate and the so-called z-star coordinate where all water levels are stretched locally according to a factor depending on the sea-surface elevation and the water depth.

Several parameterizations for sub-grid processes are available in ICON-O.
Upper-ocean vertical mixing is realized by a turbulent kinetic energy scheme following {term}`Gaspar et al. 1990`.
Deep-ocean mixing can be either approximated by constant background values of a vertical diffusivity or a constant value for turbulent kinetic energy.
Furthermore, ICON-O has the option of applying the novel energetically consistent parameterization IDEMIX ({term}`Olbers and Eden 2013`) for internal-wave driven mixing (see {term}`Brüggemann et al. 2024` for details).
Horizontal dissipation of momentum is achieved by harmonic and biharmonic operators where Smagorinsky and Leith closures can be applied to specify the horizontal viscosities.
In ICON-O simulations that are too coarse to resolve mesoscale eddies, a structure-preserving formulation of the Gent and McWilliams parameterization can be applied (see {term}`Korn 2018` for details).

(ref_ocean_seaice)=
### Sea-ice model
The sea ice model is part of ICON-O. It consists of a dynamic and a thermodynamic component.
Sea ice thermodynamics describe freezing and melting by a single-category, zero-layer formulation ({term}`Semtner 1976`). The current sea ice dynamics are based on the sea ice dynamics component of the Finite-Element Sea Ice Model (FESIM) ({term}`Danilov et al. 2015`). The sea ice model solves the momentum equation for sea ice with an elastic–viscous–plastic (EVP) rheology.
Since ICON-O and FESIM use different variable staggering, a wrapper is needed to transfer variables between the ICON-O grid and the sea ice dynamics component (see {term}`Korn et al. 2022`).
A new sea ice dynamics model has been developed to bypass these limitations (see {term}`Mehlmann and Korn 2021`) and will be employed in the future.

(ref_ocean_biogeochem)=
### Ocean biogeochemistry model
The ocean biogeochemistry component is provided by HAMOCC6 ({term}`Ilyina et al. 2013`). It simulates at least 20 biogeochemical tracers in the water column, following an ex- tended nutrient, phytoplankton, zooplankton, and detritus approach, also including dissolved organic matter, as described in {term}`Six and Maier-Reimer 1996`.
It also simulates the upper sediment by 12 biologically active layers and a burial layer to represent the dissolution and decomposition of inorganic and organic matter as well as the diffusion of pore water constituents.
The co-limiting nutrients consist of phosphate, nitrate, silicate, and iron.
A fixed stoichiometry for all organic compounds is assumed.

(ref_ocean_skin)=
### Ocean warm layer and cold skin
The ocean surface layer often features a __warm layer__ caused by solar radiation penetrating up to 3 meters, resulting in a strong diurnal cycle. Also a __cold skin__ develops due to cooling effects of surface latent and sensible heat fluxes, influencing the top millimeter of the ocean.

The ICON implementation is based on the ideas of {term}`Zeng and Beljaars 2005` with modifications from {term}`Takaya, Bidlot, Beljaars and Janssen 2010`.
The following applications are currently supported:
- atmospheric forecasts
- forecasts coupled to an ocean
- full data assimilation (DA) with cycling of the variables `sst_warm_layer` and `sst_cold_cycle` including weak coupling to ocean DA.

The ocean surface layer parameterisation components - warm layer and cold_skin - can be turned on separately with the parameters `itype_oskin_warm` and `itype_oskin_cold` (0/1 meaning off/activated respectively).

## Model configurations

### Configuring ICON-O on Levante (DKRZ)

To set up an ICON-O model simulation, it is recommended to use the tool `make_target_runscript`.
Before, this script can be applied, a template for an ocean simulation needs to be copied into `./run`.
We suggest to begin with the `exp.ocean_omip` template:

```bash
cd run/
cp checksuite.ocean_internal/omip/exp.ocean_omip .
./make_target_runscript in_script=exp.ocean_omip \
                        in_script=exec.iconrun \
                        EXPNAME=exp.ocean_omip \
                        cpu_time=00:10:00 \
                        no_of_nodes=1 \
                        openmp_threads=4 \
                        account_no=<your-project-number>
```

This generate the `run/exp.ocean_omip.run`.
In this run script several replacements need to be done:

```none
RXBY -> Replace with desired grid type, there are templates for certain configurations like R2B4, R2B6, and R2B8
ZZZ  -> Replace with desired vertical grid, e.g. L72
NNN  -> Replace with desired vertical coordinate type: 0 for zlev and 1 for zstar
```

After this basic configuration, a first simulation can be started by:
```bash
sbatch exp.ocean_omip.run
```

Once the model is running, a directory is generated `./experiments/exp.ocean_omip`, where the model output is saved.

### Configuring HAMOCC

_to be added_

## Ocean diagnostics

With ICON-O not only the prognostic variables but also a multitude of other diagnostic variables can be saved.
In general, the output is configured within the `output_nml` namelist sections by blocks similar to this example:

```fortran
&output_nml
  filetype                        = 5
  output_filename                 = "${EXPNAME}_P1M_3d"
  output_start                    = "${start_date}"
  output_end                      = "${end_date}"
  output_interval                 = "P1M"
  filename_format                 = "<output_filename>_<datetime2>"
  operation                       = "mean"
  file_interval                   = "${file_interval}"
  mode                            = 1
  include_last                    = .FALSE.
  output_grid                     = .TRUE.
  ml_varlist                      = 'to','so','u','v','w','A_veloc_v','A_tracer_v_to','A_tracer_v_so',
                                    'normal_velocity','rho','rhopot','mass_flux',
                                    'heat_content_liquid_water','swrab','rsdoabsorb',
                                    'group:oce_eddy'
/
```

Here, the following definitions apply:

```
filetype: 5 for netcdf
output_filename: first part of the filename, e.g."${EXPNAME}_P1M_3d"
filename_format: structure of the filename, e.g. "<output_filename>_<datetime2>"
output_start: ???
output_end: ???
output_interval: Specify how often an output should be written. "P1M" means every month an output will be written.
file_interval: Specify the interval once a new file will be opened. "P1Y" means every year a new file will be created.
operation: Specify if a snapshot (delete the line) or a time average ("mean") will be written.
mode: ????
include_lat: ????
output_grid: If .TRUE. the grid will be written to the file
ml_varlist: a list of variable names specified by `add_var` in the model, please use single quotes (') to encapsulate each variable.
```

Most output variables are structured in groups.
The most important groups are:

```
oce_default         -> default model output
oce_essentials      -> ???
ice_default         -> default ice output
ocean_monitor       -> globally integrated quantities
ocean_moc           -> diagnostics regarding meridional overturning and heat transport
oce_eddy            -> diagnostics for products of eddy quantities like eddy kinetic energy
oce_layers          -> diagnostics for analysis in isopycnal space
oce_vmix_tke        -> diagnostics of the TKE mixing scheme
oce_vmix_iwe        -> diagnostics of the internal wave mixing scheme
oce_ts_budget       -> diagnostics for temperature and salinity budgets
```

Within the post-processing toolbox [pyicon](https://gitlab.dkrz.de/m300602/pyicon/-/tree/master/pyicon) there are several example python notebooks that illustrate how different diagnostics can be applied to study ocean dynamics.
