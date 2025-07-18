# mkexp for ICON

MakeExperiments! (_mkexp_) is a tool for preparing experiments with MPI-M's earth system models. It provides a unified command line interface to perform experiments uisng ICON configurations maintained by the MPI-M.

Information on its use is provided in the [Quick Start](./doc/Quick_Start.md) guide.


## Basic structure

To set up an experiment, you create a config file eg.

```
EXP_TYPE = torus
EXP_OPTIONS =
ENVIRONMENT = ios
[jobs]
  [[run_start]]
    time_limit = 00:05:00
```

and run _mkexp_ on that file. The configuration is then gathered from default settings, experiment type, possibly options, the computing environment settings, and finally the config file itself. Latter settings amend the former.

The final config is applied to the requested job templates, creating one script for each. Templates are sought in the same locations as configs except for the environment. The latest template found is chosen to create the script. In addition, the environment template is included at the beginning of each script.

For technical details, see `utils/mkexp/doc/mkexp.pdf`.

### Default settings and templates (`mkexp/run/defaults`)

* Global default configuration for all experiments (`DEFAULT.config`)
* Default job templates (`DEFAULT.*.tmpl`)
* Ancillary templates as required by formatting or job templates

### Experiment types (`mkexp/run/types`)

* Experiment specific configuration (`<EXP_TYPE>.config`)
* Specialized job templates if required (`<EXP_TYPE>.*.tmpl`)

### Packaged Options (`mkexp/run/options`)

* Collections of optional config settings per topic (`<EXP_OPTION>.config`)

### Computing environment settings (`mkexp/run/environments`)

* Fallback environment settings if unknown or not specified (`DEFAULT.config`)
* Specific environment settings and job resource specifications (`<ENVIRONMENT>.config`)
* Environment template formatting information according to system needs (`<ENVIRONMENT>.tmpl`)
