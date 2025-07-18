```{eval-rst}
:orphan:
```

(ref_buildrun_recommconf)=
# Recommended Configurations

For testing purposes we provide a number of use cases. All use cases provide everything necessary for running ICON:

* [grids and external parameters](ref_buildrun_gridextpar),
* [initial and (where necessary) boundary conditions](ref_buildrun_icbc),
* namelist input to start the ICON run

_Please note that the preparation of further use cases is on-going._


:::{admonition} Support
:class: admonition-icontheme
Besides testing, these use cases serve as recommended configurations which are supported by the developers.
There will be no support from the ICON partners for configurations that are not listed below.
:::

## NWP Use Cases

The use cases for NWP are available on the [DKRZ Swiftbrowser](https://swiftbrowser.dkrz.de/public/dkrz_4d992e1b-f237-4258-a2bc-138ca6a1cf59/icon-model-use-cases/).
From there you can download tar-balls for every use case.
Every tar-ball contains a README with additional information on how to run the use case.

nwp-global-R02B06:
: This use case runs a global ICON application with a resolution of about 40 km (R02B06). It includes a nest (R02B07) over Europe.

nwp-local-R19B07:
: This use case runs a limited-area (local) ICON application over Germany with a resolution of about 2 km (R19B07).
  It is comparable to DWD's operational application ICON-D2.
