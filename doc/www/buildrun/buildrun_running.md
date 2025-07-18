```{eval-rst}
:orphan:
```

(ref_buildrun_running)=
# Running ICON

The ICON model is typically run through a *runscript*, which sets up the working directory, populates it with all required input files (grid files, namelists, etc.), sets environment variables, runs the model, and postprocesses its output.

You can generate a runscript with one of the two tools described in this document. The first tool [`mkexp`](ref_buildrun_mkexp) is well documented and continuously maintained. The only downside is that it does not have support for many environments and experiments that you can run with ICON yet. The second one [`make_runscript`](ref_buildrun_make_runscript) is a legacy, poorly documented set of shell scripts, which, however, supports a lot of (mainly HPC) environments and experiments.

(ref_buildrun_mkexp)=
## Using `mkexp` to prepare ICON experiments

[MakeExperiments! (`mkexp`)](https://gitlab.dkrz.de/esmenv/mkexp) is a Python tool for preparing experiments with ICON. It helps users set up an experimental workflow, and generate the runscript needed to execute supported configurations. The tool presently supports setting up ICON configurations in the DKRZ environment (CPU and GPU) and is and can be adapted to other environments.

You can find the `mkexp` tool in the `utils/mkexp` directory of the ICON source code repository, which is managed as a git submodule.

### Requirements

Before you start using `mkexp`, make sure that your software environment meets the following requirements:

1. The tool needs the Python interface of the [MTIME](https://gitlab.dkrz.de/icon-libraries/libmtime) library. The most simple way to get the required Python module is to configure ICON with the `--enable-bundled-python=mtime` argument. For example:
    ```sh
    ./config/generic/gcc --enable-bundled-python=mtime
    ```
2. Additionally, you need to make sure that [Jinja2](https://pypi.org/project/Jinja2/) and [six](https://pypi.org/project/six/) are available in your Python environment.
3. If ICON is configured and built [out-of-source](ref_buildrun_configuration_oos), you need to set the `ICON_BUILD_DIR` environment variable to the absolute path to the root build directory of ICON:
    ```sh
    export ICON_BUILD_DIR=/path/to/icon/build/directory
    ```

### Steps to run an experiment

Running an experiment using `mkexp` generally consists of three steps briefly described below in this section. For more details, see the [documentation](https://gitlab.dkrz.de/esmenv/mkexp/-/blob/master/doc/mkexp.pdf).

**Step 1: Create the configuration file**

The configuration file overrides the generic parameters specified in `run/mkexp/defaults/DEFAULT.config`. The experiment- and environment-specific parameters can be overridden or set via the `EXP_TYPE` and `ENVIRONMENT` variables of the configuration file. The values of those variables correspond to the `.config` files in the `run/mkexp/types` and `run/mkexp/environments` directories. The configuration file of an experiment can override or set additional parameters via the `EXP_OPTIONS` variable (the value corresponds to the files in the `run/mkexp/options` directory or directly. The resulting set of parameters is used by `mkexp` to create the experiment runscript, set the run environment and configure the required directory structure for the run.

The most simple way to create the configuration file for an experiment is to copy one of the `/run/examples` (e.g. `/run/examples/bubble.config`) and adjust it to your needs. First, you need to come up with the experiment identifier (e.g. `exp_id`) and copy the example configuration to the {{ '[`run`]({}/run)'.format(base_url) }} directory under the name that corresponds to that identifier:

```sh
cd ./run
cp ./examples/bubble.config ./exp_id.config
```

Now you need to review and edit the contents of `exp_id.config`. If you intend to run the experiment on the DKRZ machine, make sure the `ACCOUNT` is set correctly (it should be set to the SLURM account that you normally submit jobs with). Running the experiment in an unknown environment requires more adjustments. For example, if you want to run the bubble experiment on your personal machine, we recommend considering the following changes:

```diff
@@ -5,4 +5,4 @@
 EXP_TYPE = torus
-ENVIRONMENT = levante
 ACCOUNT = mh0287
+INPUT_ROOT = $HOME/data

@@ -18,2 +18,5 @@ OUTPUT_INTERVAL = $ATMO_TIME_STEP

+# workaround for set-up info
+use_build_env =
+
 [namelists]
@@ -30,2 +33,4 @@ OUTPUT_INTERVAL = $ATMO_TIME_STEP
       output_grid = true
+    [[[parallel_nml]]]
+      num_io_procs = 0
     [[[output_nml atm_3d]]]
@@ -38,3 +43,5 @@ OUTPUT_INTERVAL = $ATMO_TIME_STEP
     nodes = 1
-    threads_per_task = 4
+    threads_per_task = 2
+    cpus_per_node = 4
+    hardware_threads = true
     time_limit = 00:05:00
```

In the macOS environment, you also need to override one of the diagnostic utilities:

```diff
@@ -37,2 +37,3 @@ OUTPUT_INTERVAL = $ATMO_TIME_STEP
   [[run]]
+    ldd = otool -L
     nodes = 1
```

**Step 2: Generate the scripts and workflow environment**

Execute [`mkexp`](https://gitlab.dkrz.de/esmenv/mkexp/-/blob/master/mkexp):

```sh
../utils/mkexp/mkexp exp_id.config
```

The command will create the required directory structure. For example:

```sh
Script directory: '/path/to/icon-srcdir/experiments/exp_id/scripts'
Data directory: '/work/your-account-number/your-user-number/master/experiments/exp_id/outdata'
Work directory: '/scratch/your-account-type/your-account-number/master/experiments/exp_id/work'
Log directory: '/work/your-account-number/your-user-number/master/experiments/exp_id/log'
```

The command will also create the `exp_id.run_start` runscript and place it to the `Script directory`. Review the script to make sure that the path to the grid file (`icon_grid_G.nc`) is set correctly for your environment.

**Step 3: Execute the runscript.**

For the last step, switch to the aforementioned `Script directory` and either execute or submit (if you are in the HPC environment) the runscript:

```sh
cd ../experiments/exp_id/scripts
sbatch exp_id.run_start
```

(ref_buildrun_make_runscript)=
## Using `make_runscript` to prepare ICON experiments

The `make_runscript` shell script is an ICON-specific tool for runscript generation. It takes the experiment template files from the {{ '[`run`]({}/run)'.format(base_url) }} directory, prepends the environment-specific shell snippets from the `run/create_target_header` file, adjusts the result based on how ICON is configured and built (see `run/collect.set-up.info.in`) and produces a shell script that is ready for the execution or submission.

### Requirements

The tool does not support the [out-of-source](ref_buildrun_configuration_oos) builds. To circumvent this, most of the existing [configure wrappers](ref_buildrun_configuration_wrappers) make the required files available in the build directory:

```sh
# Copy runscript-related files when building out-of-source:
if test $(pwd) != $(cd "${icon_dir}"; pwd); then
  echo "Copying runscript input files from the source directory..."
  rsync -uavz ${icon_dir}/run . --exclude='*.in' --exclude='.*' --exclude='standard_*' --exclude=mkexp
  ln -sf -t run/ ${icon_dir}/run/{standard_*,mkexp}
  rsync -uavz ${icon_dir}/externals . --exclude='.git' --exclude='*.f90' --exclude='*.F90' --exclude='*.c' --exclude='*.h' --exclude='*.Po' --exclude='tests' --exclude='*.mod' --exclude='*.o'
  rsync -uavz ${icon_dir}/make_runscripts .
  rsync -uavz ${icon_dir}/scripts .
  ln -sf ${icon_dir}/data
  ln -sf ${icon_dir}/vertical_coord_tables
fi
```

### Steps to run an experiment

To generate a runscript based on a particular template (e.g. `run/exp.atm_tracer_Hadley`), switch to the root build directory of ICON and run `make_runscript` while providing the name of the experiment (without the prefix) as an argument:

```sh
./make_runscripts atm_tracer_Hadley
```

Alternatively, you can generate runscripts for all existing experiments:

```sh
./make_runscripts --all
```

The generated runscripts are saved to the {{ '[`run`]({}/run)'.format(base_url) }} subdirectory of the build directory. The headers of the runscripts containing arguments for the HPC workload manager, e.g. [SLURM](https://slurm.schedmd.com/), might require additional manual adjustments regarding CPU time accounting, node allocation, etc.

Once the runscript is created and adjusted, switch to the {{ '[`run`]({}/run)'.format(base_url) }} subdirectory of the root build directory of ICON and either execute or submit (if you are in the HPC environment) the runscript:

```sh
cd ./run
sbatch ./exp.atm_tracer_Hadley.run
```

:::{admonition} checksuite.nwp specifics
:class: admonition-icontheme
Experiments in the `run/checksuite.nwp` directory are prepared to be run in the hybrid mode, which uses two icon binaries. To generate those experiments the user can specify the `-r <run directory to process>` argument:
```sh
./make_runscripts --all -r run/checksuite.nwp
```
In the case of the NEC-Aurora, the second binary is the host (or x86 scalar) binary. It can be found automatically if the two build directories are either `something/vector` and `something/host` or `something_else/VH` and `something_else/VE`. In either case, `./make_runscripts` must be called in the `vector` (or `VE`) build directory. Alternatively, the x86-host binary can be specified by the `-s <secondary build dir>` option:
```sh
./make_runscripts run_ICON_01_R3B9_lam -r run/checksuite.nwp -s ../host_gcc-9.1.0/
```
:::

:::{admonition} make_target_runscript
:class: admonition-icontheme
Alternatively, the users can employ the low-level tool `/run/make_target_runscript` for runscript generation offering more fine-grained control over certain parameters. For example, the wall clock limit and the number of allocated nodes can be injected into the runscript as follows:
```sh
cd ./run && ln -sf ./checksuite.ocean_internal/omip/exp.ocean_omip_long exp.ocean_omip_long
./make_target_runscript in_script=exp.ocean_omip_long in_script=exec.iconrun \
  EXPNAME=ocean_omip_long cpu_time=08:00:00 no_of_nodes=20
```
:::

(ref_buildrun_gridextpar)=
## Grids & External Parameters

### Grid Files

The ICON model receives information about the horizontal grid from so-called **grid files** in the [NetCDF format](https://www.unidata.ucar.edu/software/netcdf/).
These files store coordinates and topological index relations between cells, edges and vertices of the chosen domain.
A detailed description of the content of these grid files is provided in the _Necessary Input Data_ section of the **{term}`ICON Tutorial 2024`**.

The grid files for ICON usually follow the nomenclature `R<n>B<k>`, where `<n>` denotes the number of root divisions and `<k>` the number of subsequent bisections.
From `<n>` and `<k>` the resolution of the grid can be estimated by the formula:

```{math}
\Delta x \sim \frac{5050}{n \cdot 2^k} \quad km.
```

:::{admonition} Download Grid & External Parameter Data
:class: admonition-icontheme
A set of predefined grid and external parameter datasets is available at **[http://icon-downloads.mpimet.mpg.de/](http://icon-downloads.mpimet.mpg.de/)**.
:::

(ref_buildrun_external_param)=
### External Parameters (NWP)

_Please note that this description applies to the [](ref_atmosphere_nwp_physics)_.

External parameter datasets contain topological and climatological data that is assumed to be constant during a typical NWP integration.
These datasets are aggregated to a given ICON grid using the **[EXTPAR Software](http://www.cosmo-model.org/content/support/software/default.htm)**.
Like for the grid files, a more detailed description is given in the _Necessary Input Data_ section of the **{term}`ICON Tutorial 2024`**.

(ref_buildrun_icbc)=
## Initial & Boundary Data
