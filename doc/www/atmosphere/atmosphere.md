(ref_atmosphere)=
# Atmosphere Overview

The ICON atmosphere model predicts the spatio-temporal evolution
of the atmospheric state in terms of the prognostic variables
virtual potential temperature, 3D wind, total air density
and mass fractions of atmospheric water constituents and trace gases.
In addition, the model provides a comprehensive set of diagnostic quantities,
such as surface pressure, wind gusts or potential vorticity
just to name a few. An extensive, but still incomplete list of available
output variables is provided in Appendix A of the {term}`ICON Tutorial 2024`.

In mathematical terms, the ICON atmosphere model solves the fully compressible non-hydrostatic
Navier-Stokes equations on the sphere. The explicitly resolved scales of
motion are treated by the so called **[](ref_atmosphere_dycore)**.
The latter is accompanyied by a set of physical parameterizations which account for the effect motions
that fall below a chosen mesh size.
ICON offers two different physics packages which are known as the AES Physics Package and the NWP Physics Package.

The **[](ref_atmosphere_aes_physics)** was originally derived from the physics package of the ECHAM model and subsequently further developed.
Applications range from general circulation model on long time scales to km-scale, storm-resolving climate simulations.
The [](ref_atmosphere_aes_physics) can be chosen by setting the namelist parameter {term}`iforcing``=2`.

The **[](ref_atmosphere_nwp_physics)** parameterizations were chosen from different sources, most notably the [COSMO-Model](https://www.cosmo-model.org) and the [IFS model](https://www.ecmwf.int/en/forecasts/documentation-and-support/changes-ecmwf-model).
Originally developed for numerical weather prediction, the [](ref_atmosphere_nwp_physics) was extended for seamless application across scales.
The [](ref_atmosphere_nwp_physics) can be activated by setting {term}`iforcing``=3`.

(ref_atmosphere_dycore)=
## Dynamical Core

The dynamical core can be considered as the foundations of any numerical model.
It predicts the evolution in space and time of all
atmospheric motions which are resolvable on a given mesh.
To this end, the dynamical core solves the Euler equations which is a set
of partial differential equations describing adiabatic and inviscid flow.
There exist various approximative forms of the Euler equations in which
certain types of motion are filtered out that are difficult to handle numerically.
Well known examples are the hydrostatic form or the Boussinesq form,
both of them not supporting the propagation of acoustic waves.

The ICON dynamical core solves the fully compressible form of the Euler
equations, which does support acoustic waves.
Notable approximations to the exact Euler equations relate to the treatment
of the Earth as a spherical geoid and the shallow atmosphere approximation
(see e.g. {term}`Thuburn & White 2013`).
The latter may be deactivated with the Namelist switch
{term}`ldeepatmo``=.TRUE.`, leading to the so called deep atmosphere form of the
governing equations ({term}`Borchert et al. 2019`). Given a suitable set of physical
parameterizations, this allows for simulations on model domains
reaching all the way up to the lower thermosphere.

The governing set of equations is discretized on an Icosahedral-triangular C-grid
in horizontal directions, while in vertical direction a height-based terrain following
coordinate with Lorenz-type staggering of the prognostic variables is used.
The discrete numerical operators, such as the divergence,
gradient or laplace operator, are constructed using a mixture of finite
difference and finite volume discretizations of mostly second order-accuracy.
They are combined with a predictor-corrector type two-level time integration
scheme, leading to a discretization which is mass conserving, but not strictly
energy conserving.

For additional details on the dynamical core, the reader is referred to {term}`Zaengl et al. 2015`
and Chapter 3 of the {term}`ICON Tutorial 2024`.

(ref_atmosphere_tracer_transport)=
## Tracer Transport

The tracer transport module is an important building block of any weather
prediction or climate model.
It solves the tracer mass continuity equation,
in order to describe the redistribution of gaseous, liquid or solid
atmospheric constituents such as water vapour or rain water due to air
motion or gravitational settling (sedimentation).

In ICON, finite volume methods of second order accuracy in time and up to fourth order accuracy
in space are applied to construct mass conserving and mass consistent transport schemes.
If needed, these schemes can be combined with monotonicity or positivity preserving limiters.

More details on the tracer transport module can be found in the ICON reports
({term}`Reinert 2020`, {term}`Reinert & Zaengl 2021`) and Section 3.6 of the {term}`ICON Tutorial 2024`.

(ref_atmosphere_physics)=
## Physical Parameterizations

(ref_atmosphere_aes_physics)=
### AES Physics Package

_to be added_

(ref_atmosphere_nwp_physics)=
### NWP Physics Package

:::{table} Overview NWP physics package
:width: 65
:widths: auto
:align: center

| Parameterization                     | References    | Namelist Parameter |
| :-----------                         | :------------ | :------------      |
| **Radiation**                        | RRTM ({term}`Mlawer et al. 1997`, {term}`Barker et al. 2003`), ecRad ({term}`Hogan & Bozzo 2018`) | {term}`inwp_radiation` |
| **Non-Orographic Gravity Wave Drag** | {term}`Orr et al. 2010` | {term}`inwp_gwd` |
| **Sub-grid scale Orographic Drag**   | {term}`Lott & Miller 1997` | {term}`inwp_sso` |
| **Cloud Cover**                      | - | {term}`inwp_cldcover` |
| **Microphysics**                     | Single Moment ({term}`Doms et al. 2011`), Double Moment ({term}`Seifert & Beheng 2006`), SBM ({term}`Khain & Sednev 1996`, {term}`Khain et al. 2004`) | {term}`inwp_gscp` |
| **Convection**                       | {term}`Tiedtke 1989`, {term}`Bechtold et al. 2008` | {term}`inwp_convection` |
| **Turbulent Transfer**               | Prognostic TKE ({term}`Raschendorfer 2001`), 3D Smagorinsky ({term}`Smagorinsky 1963`, {term}`Lilly 1962`) | {term}`inwp_turb` |
| **Land**                             | See **[](ref_land_schemes)** | {term}`inwp_surface` |
:::

More detailed descriptions of some of above options are available here:

::::{grid} 1 2 2 3
:gutter: 1 1 1 2

:::{grid-item-card}
**Radiation (ecRad)**
^^^
[](ref_atmosphere_ecrad)  
[](ref_atmosphere_ecrad_redgrid)  
[](ref_atmosphere_ecrad_aerosol)
:::
::::

_to be extended_

You can find a brief overview on the NWP physics package in chapter 3 of the **{term}`ICON Tutorial 2024`**.

(ref_atmosphere_nwp_waves)=
## Waves (NWP)

_to be added_

### Glossary of Namelist Parameters

_Operational NWP setting marked by {material-regular}`settings;1em;pst-color-secondary`_

:::{glossary}
iforcing
  (`&run_nml`) Forcing of dynamics and transport by parameterized processes. 2: AES forcing, 3:{material-regular}`settings;1em;pst-color-secondary` NWP forcing

ldeepatmo
  (`dynamics_nml`) Switch for deep-atmosphere modification of non-hydrostatic atmosphere. Specific settings can be found in `&upatmo_nml`.

inwp_radiation
  (`&nwp_phy_nml`) Radiation parameterization. 1: RRTM radiation, 4:{material-regular}`settings;1em;pst-color-secondary` ecRad radiation

inwp_gwd
  (`&nwp_phy_nml`) 1:{material-regular}`settings;1em;pst-color-secondary` Orr et al. scheme

inwp_sso
  (`&nwp_phy_nml`) 1:{material-regular}`settings;1em;pst-color-secondary` Lott-Miller scheme

inwp_cldcover
  (`&nwp_phy_nml`) 1:{material-regular}`settings;1em;pst-color-secondary` Diagnostic PDF 5: All or nothing scheme (grid-scale clouds)

inwp_gscp
  (`&nwp_phy_nml`) 1:{material-regular}`settings;1em;pst-color-secondary` Single moment 2:{material-regular}`settings;1em;pst-color-secondary` Single moment incl. graupel 4: Double moment 8: Warm spectral bin microphysics

inwp_convection
  (`&nwp_phy_nml`) 1:{material-regular}`settings;1em;pst-color-secondary` Tiedtke-Bechtold

inwp_turb
  (`&nwp_phy_nml`) 1:{material-regular}`settings;1em;pst-color-secondary` Prognostic TKE (COSMO) 5: 3D Smagorinsky diffusion

inwp_surface
  (`&nwp_phy_nml`) 1:{material-regular}`settings;1em;pst-color-secondary` TERRA
:::
