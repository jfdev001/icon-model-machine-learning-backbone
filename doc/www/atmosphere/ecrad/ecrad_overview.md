```{eval-rst}
:orphan:
```

(ref_atmosphere_ecrad)=
# ecRad Overview

The NWP physics package of ICON uses the ecRad radiation scheme developed by [ECMWF](https://ecmwf.int).
This page aims to described the configuration and input options for ecRad that are available in ICON.
For a more detailed and complete description of the radiation scheme itself is provided in the [ECMWF Confluence Wiki](https://confluence.ecmwf.int/display/ECRAD).
The original source code is publicly available on [GitHub](https://github.com/ecmwf-ifs/ecrad), please note that ICON uses a modified version of this code.

(ref_atmosphere_ecrad_redgrid)=
# Reduced Radiation Grid

Radiation is one of the computationally expensive physical parameterizations.
There are several possibilities to decrease the computational cost of radiation, e.g. by reducing the temporal, spatial or spectral resolution.
By activating **{term}`lredgrid_phys``=.true.`** and specifying the corresponding grid file with **{term}`radiation_grid_filename`**, the radiation is calculated on a one grid level coarser domain to reduce the computational cost of the radiation by about a factor of 4.

**We highly recommend activating the reduced radiation grid option for the following reasons:**

- It is computationally cheaper by a factor of 4, usually without a degradation of the results. Since the radiation is treated as a slow physics process, it is not called every time step anyways. A coarser horizontal grid thus fits better to the advective time scale.
- Only for the reduced radiation grid, there is an additional option named **{term}`latm_above_top`**. This option adds an extra layer at the top to account for the incoming long-wave radiation. This reduces the biases at the model top significantly.
- For global domains, there is a load balancing for sunlit and shadowed parts of the earth for the reduced grid.

For a more detailed description of the reduced radiation grid implementation, see {term}`ICON Tutorial 2024`.

(ref_atmosphere_ecrad_aerosol)=
# Aerosol Input Options

## Tegen climatology

Climatological aerosol based on the {term}`Tegen et al. 1997` climatology can be selected by choosing **{term}`irad_aero``=6`**.
This options has the following characteristics:

- Optical thicknesses at the wavelength 550 nm of the 5 species **Sea Salt**, **Soil Dust**, **Sulfate**, **Organic Carbon** and **Black Carbon** are provided in the [external parameter file](ref_buildrun_external_param).
- The annual cycle is considered by providing monthly data which is linearly interpolated inside ICON to the target date.
- The original data is vertically integrated optical thickness. For the use in ecRad, an exponentially decaying, normalized vertical profile is added by ICON.
- The target variables optical thickness (SW/LW), single scattering albedo (SW) and asymmetry parameter (SW) at the radiation wavelength bands are derived based on lookup tables in the ICON code.

## CAMS climatology or CAMS forecast aerosol

ICON currently supports the use of either the CAMS 49R2 aerosol climatology, or CAMS forecast aerosol for direction aerosol-radiation interactions. Aerosol-cloud interactions using CAMS are not yet supported.

###	Using the 49R2 CAMS climatology

ICON supports use of the recent (December 2024) CAMS climatology, version 49R2. This option is activated with the namelist parameter **{term}`irad_aero``=7`**.

Older versions (43R3) are no longer supported. Aerosol mixing ratios are supplied monthly on 21 pressure surfaces.
Aerosol species affected by human activity have an additional dimension `epoch` covering thirteen 5-year long periods from 1955 to 2015.
The original climatology file can be downloaded from [this ECMWF webpage](https://aux.ecmwf.int/ecpds/home/radiation/aerosol_climatology/aerosol_cams_climatology_49r2_1951-2019_4D.nc).

At this point in time, ICON does **not** support the new `epoch` dimension of the climatology (i.e. anthropogenic change of aerosol over 5+ year periods).
The ICON repository contains the script {{ '[`make_camsclim_onICONgrid.sh`]({}/scripts/preprocessing/make_camsclim_onICONgrid.sh)'.format(base_url) }} to extract the latest (2015) epoch from the original data file and interpolate the climatology onto an ICON grid of the user’s choice.
Installation of CDO, NCO and python3 (numpy, xarray) tools is required to run this script.

###  Using CAMS forecasts

ICON can also use CAMS forecast aerosol fields on 137 model levels. This option is activated with the namelist parameter **{term}`irad_aero``=8`**.

CAMS forecast aerosol fields can be retrieved from ECMWF via MARS request. The ICON repository contains the script {{ '[`make_camsforc_onICONgrid.sh`]({}/scripts/preprocessing/make_camsforc_onICONgrid.sh)'.format(base_url) }} which then interpolates the CAMS forecast aerosol onto an ICON grid of the user’s choice.
The script header contains more information on how to retrieve CAMS forecast aerosol from MARS.
CDO, python3 and ecmwf-toolbox are required to run this script.

### Information relevant to both CAMS climatology and forecast

#### A note on remapping to the ICON grid

CDO remapping is used for the interpolation. None of the available remapping options are perfect:
- Conservative remapping (**remapcon**) leads to visible 'squares' corresponding to the original, coarser CAMS climatology in the interpolated fields, which translates into visible 'square' shapes in the clear-sky radiation also. This is obviously undesireable.
- Bicubic remapping (**remapbic**) leads to 'overshooting' features around high orography (Himalayas, Andes). Also undesireable.
- The option considered to be best at the moment (and **implemented by default**) is bilinear remapping (**remapbil**), which produces a reasonably smooth field without obvious overshooting.

However, the 'best' option may depend on the application: for regional simulations away from steep orography, **remapbic** may be more advantageous (smoother).

#### Vertical interpolation onto ICON model levels

ICON interpolates from the original CAMS vertical levels (21 pressure surfaces for climatology, 137 IFS model levels for forecast) to the ICON vertical model levels. Because the horizontal resolution of the original CAMS files can be much lower than the ICON resolution, the difference in surface pressure between ICON and CAMS can be in excess of 200hPa around steep orography. A straight interpolation between pressure levels would therefore neglect a significant amount of aerosol mass if the lowest 200hPa of the CAMS profile were ignored.
To avoid this, the CAMS profile is re-distributed between the ICON surface pressure and top of the atmosphere before interpolation, conserving total aerosol mass in the column.

#### Differences to previous CAMS versions

The previous CAMS climatology 43R3 provided aerosol as **layer mass**. This climatology now also exists (v2) in a format providing **mixing ratios** instead. Since layer mass will no longer be used from the CAMS side, reading in layer mass is no longer supported by ICON.

The 43R3 CAMS climatology also had some unrealistic features, such as high dust accumulations in the stratosphere. These artefacts have now been removed in the 49R2 climatology. Since it is unlikely that the older, flawed climatology will continue to be used, the option to read in the 43R3 climatology is also no longer supported. However, it is easy to make ICON compatible with this older version again, if so desired, by changing the number of CAMS levels (`nlev_cams`) in `mo_reader_cams.f90` from 21 back to 60, and using the more recent mixing ratio version of the [43R4 CAMS climatology which can be obtained from ECMWF](https://aux.ecmwf.int/ecpds/home/radiation/aerosol_climatology/aerosol_cams_climatology_43r3_v2_3D.nc).

#### Appropriate aerosol optical properties

The new 49R2 CAMS climatology was created using an updated version (CY48R1) of the CAMS aerosol model, which contains some significant changes relative to the older version used to create the 43R3 climatology. Therefore, the **new climatology should be used with an appropriate set of optical properties!** A tabulated list of the appropriate optical properties for each version of the climatology/forecasts can be found **[here](https://confluence.ecmwf.int/display/ECRAD/Aerosol-radiation+interactions+in+the+IFS)**.

For this implementation, the appropriate aerosol optical properties have been pre-selected for use with the CAMS 49R2 climatology. The pre-selected default for the CAMS forecasts is set to the most recent CAMS IFS cycle 49R1. If working with **older CAMS forecasts**, the **user has to adapt the selected properties** in `mo_nwp_ecrad_init.f90` (line 342 and following) according to the table linked above.

:::{admonition} Known limitations
:class: admonition-icontheme
Known limitations of the new 49R2 climatology include that the "far field" aerosol such as in the Arctic is too low.
The IFS is run with an additional artificial small constant background term to get the best results.
:::

# Glossary of Namelist Parameters

_Operational NWP setting marked by {material-regular}`settings;1em;pst-color-secondary`_

:::{glossary}
lredgrid_phys
  (`&grid_nml`) If set to `.TRUE.`{material-regular}`settings;1em;pst-color-secondary` radiation is calculated on a coarser grid (i.e. one grid level coarser).

radiation_grid_filename
  (`&grid_nml`) Filename of the grid to be used for the radiation model. Must only be specified for the base domain, since for child domains the grid of the respective parent domain serves as radiation grid. An empty string is required, if radiation is computed on the full (non-reduced) grid.

latm_above_top
  (`&nwp_phy_nml`) Adds an extra layer at the model top to account for the incoming long-wave radiation if set to `.TRUE.`{material-regular}`settings;1em;pst-color-secondary`.

irad_aero
  (`&radiation_nml`) Specify aerosol input for radiation. **0:** None, **3:** externally specified (e.g. [](ref_tools_comin)) **6:** {material-regular}`settings;1em;pst-color-secondary` Tegen climatology, **7:** CAMS 3D climatology, **8:** CAMS 3D forecasted, **9:** [](ref_atmosphere_art), **12:** tropospheric Kinne climatology (constant in time), **13:** tropospheric Kinne climatology (time-dependent), **14:** volcanic stratospheric aerosols for CMIP6 (time dependent), **15:** combination of 13 and 14, **18:** tropospheric natural Kinne climatology + volcanic stratospheric aerosols + anthropogenic 'simple plumes' (time-dependent), **19:** as 18 without volcanic stratospheric aerosols
:::
