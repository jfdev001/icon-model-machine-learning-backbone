```{eval-rst}
:orphan:
```

```{image} logo_art_black.svg
:class: only-light
:width: 200px
:align: center
```

```{image} logo_art_white.svg
:class: only-dark
:width: 200px
:align: center
```

(ref_atmosphere_art)=
# Aerosols & Reactive Trace gases (ART)

## Introduction to the ART Module
The [ART](https://www.icon-art.kit.edu) (**A**erosols and **R**eactive **T**race gases) module is an extension of the ICON model that significantly enhances its capabilities by enabling the simulation of gases, aerosol particles, and related feedback processes in the atmosphere. It was developed and is provided by the Karlsruhe Institute of Technology (KIT).

## Directory Structure and Naming Conventions
When setting up the ICON model, the root directory is typically named "icon-model." For consistency, this directory is referred to as the "ICON folder". The ART module is incorporated into the ICON model by including the ART code within the icon-model directory, specifically under the path _icon-model/externals/art_, which will be referred to as the "ART folder." The integrated system of the ICON model and the ART module is referred to as the "ICON-ART model" throughout this guide.

## User Guide Details
:::topic
For detailed instructions and guidelines on utilizing the ART module, the comprehensive [**ICON-ART User guide**](https://www.icon-art.kit.edu/userguide/index.php?title=Main_Page) provides essential information and step-by-step procedures for effectively using the ART module within the ICON model.
:::

To begin, it is necessary to configure the required [input](https://www.icon-art.kit.edu/userguide/index.php?title=Input) files and parameters, which involves setting up the [namelist](https://www.icon-art.kit.edu/userguide/index.php?title=Namelist) that defines the configurations for the simulation. Once the model is operational, it will generate [output](https://www.icon-art.kit.edu/userguide/index.php?title=Output) parameters that can be printed to provide insights into various atmospheric processes and dynamics.

[Configuring aerosol dynamics](https://www.icon-art.kit.edu/userguide/index.php?title=AERODYN) and [atmospheric chemistry](https://www.icon-art.kit.edu/userguide/index.php?title=Atmospheric_Chemistry) settings are critical aspects of using the ART module, enabling the simulation of complex interactions within the atmosphere. After running simulations, the [postprocessing](https://www.icon-art.kit.edu/userguide/index.php?title=Postprocessing) phase assists in extracting and visualizing the data.

For those interested in customizing the ART module, the guide also covers [programming ART](https://www.icon-art.kit.edu/userguide/index.php?title=Programming_ART), offering insights into how to adapt the module to specific needs. [Tutorial examples](https://www.icon-art.kit.edu/userguide/index.php?title=Tutorial_Examples) are included to facilitate practical applications of the ART module, and a [training course PDF](https://www.icon-art.kit.edu/userguide/index.php?title=Training_Course) is available for comprehensive instruction on using ART effectively.
