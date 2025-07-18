```{eval-rst}
:orphan:
```

(ref_buildrun_introduction)=
# Building ICON

The process of building ICON consists of two parts: *configuring* the options and compiler flags, and *building* the source code with those options and flags.

(ref_buildrun_configuration)=
# Configuration

The configuration step is done by calling the {{ '[`configure`]({}/configure)'.format(base_url) }} script with arguments specifying the location of libraries and tools required for building, as well as options enabling or disabling particular model features. For example:

```sh
./configure CC=mpicc FC=mpif90 LIBS='-lnetcdff -lnetcdf -llapack -lblas' --disable-ocean --disable-coupling
```

:::{admonition} Full list of configuration options
:class: admonition-icontheme
Users are recommended to get familiar with the full list of configuration options and their default values by running:
```sh
./configure --help
```
:::

The {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON is implemented using [Autoconf](https://www.gnu.org/software/autoconf/) and its interface should be familiar to those who have experience with [Autotools](https://en.wikipedia.org/wiki/GNU_Autotools)-based building systems.

The following sections provide information on some features and implementation details of the configuration process.

(ref_buildrun_configuration_icondep)=
## ICON dependencies

[Fig. 1](fig_icon_depgraph) shows a partial dependency graph of the model. A dependency can be either *mandatory* (i.e. the package is required regardless of the specified configure options) or *optional* (i.e. the package is required only if particular features of the model are enabled), and some of the dependencies are provided together with the source code of ICON as git submodules and referred to as *bundled packages* later in the text.

**_NOTE:_** The term *bundled package* does not apply to all packages listed in {{ '[`.gitmodules`]({}/.gitmodules)'.format(base_url) }}: some of them, e.g. [JSBACH](https://gitlab.dkrz.de/jsbach/jsbach) and [ICON-ART](https://www.icon-art.kit.edu), have circular dependencies with the ICON source code and therefore are treated as part of it.

```{image} ./icon-depgraph/icon-depgraph.svg
:alt: ICON dependency graph
:align: center
:name: fig_icon_depgraph
```
*Fig. 1. ICON dependency graph*

The list of packages required for successful configuration and building depends on the selected options. To make the configuration process more transparent, the {{ '[`configure`]({}/configure)'.format(base_url) }} script does not accept paths to the installation directories of the packages, which would be used to extend the corresponding compiler flags. Instead, paths to the header and library files of the packages must be provided by the user as compiler and linker flags, i.e. in the `FCFLAGS`, `CPPFLAGS`, and `LDFLAGS` arguments. Moreover, the script does not try to guess the list of libraries to be used, therefore all the `-l` linker flags need to be specified in the `LIBS` argument in the correct order. The recommended (topologically sorted) order for the `LIBS` argument is presented in [Table 1](tab_icon_depgraph), and the recommended order for the `FCFLAGS`, `CPPFLAGS`, `LDFLAGS` is the reversed one.

:::{table} *Table 1. Topologically sorted list of the ICON dependency graph*
:widths: auto
:align: center
:name: tab_icon_depgraph

| Package | Dependency condition<a name="f1back"> <sup><a href="#f1">1</a></sup> | Required flags<sup><a href="#f1">1</a></sup> |
| :---: | :---: | :---: |
| [FORTRAN-SUPPORT](https://gitlab.dkrz.de/icon-libraries/libfortran-support) | `--with-external-fortran-support` | `FCFLAGS='-I/path/to/fortran-support/include' LDFLAGS='-L/path/to/fortran-support/lib' LIBS='-lfortran-support'` |
| [ICON-TIXI](https://gitlab.dkrz.de/icon-libraries/libtixi) (a modified version of [TIXI](https://github.com/DLR-SC/tixi)) | `--enable-art --with-external-tixi` | `FCFLAGS='-I/path/to/tixi/include' LDFLAGS='-L/path/to/tixi/lib' LIBS='-licon_tixi'` |
| [XML2](https://gitlab.gnome.org/GNOME/libxml2/-/wikis/home) | `--enable-art` | `CPPFLAGS='-I/path/to/libxml2/include/libxml2' LDFLAGS='-L/path/to/libxml2/lib' LIBS='-lxml2'` |
| [YAC](https://gitlab.dkrz.de/dkrz-sw/yac) | `--enable-coupling --with-external-yac` | `FCFLAGS='-I/path/to/yac/include' LDFLAGS='-L/path/to/yac/lib' LIBS='-lyac'` |
| [FYAML](https://github.com/pantoniou/libfyaml) | `--enable-coupling` | `CPPFLAGS='-I/path/to/libxml2/include/libfyaml' LDFLAGS='-L/path/to/libfyaml/lib' LIBS='-lfyaml'` |
| [MTIME](https://gitlab.dkrz.de/icon-libraries/libmtime) (Fortran interface) | `--with-external-mtime` | `FCFLAGS='-I/path/to/mtime/include' LDFLAGS='-L/path/to/mtime/lib' LIBS='-lmtime'` |
| [MTIME](https://gitlab.dkrz.de/icon-libraries/libmtime) (C interface) | `--enable-coupling --without-external-yac --with-external-mtime` | `CPPFLAGS='-I/path/to/mtime/include' LDFLAGS='-L/path/to/mtime/lib' LIBS='-lmtime'` |
| [SERIALBOX](https://github.com/GridTools/serialbox) | `--enable-serialization` | `FCFLAGS='-I/path/to/serialbox2/include' LDFLAGS='-L/path/to/serialbox2/lib' LIBS='-lSerialboxFortran'` |
| [CDI](https://gitlab.dkrz.de/mpim-sw/libcdi) | `--with-external-cdi` | `FCFLAGS='-I/path/to/libcdi/include' LDFLAGS='-L/path/to/libcdi/lib' LIBS='-lcdi_f2003 -lcdi'` (or `LIBS='-lcdi_f2003 -lcdipio -lcdi'`) |
| [PPM](https://gitlab.dkrz.de/dkrz-sw/ppm) (C interface) | `--enable-cdi-pio --with-external-ppm --without-external-cdi` |  `CPPFLAGS='-I/path/to/ppm/include' LDFLAGS='-L/path/to/ppm/lib' LIBS='-lscalesppmcore'` |
| [ECCODES](https://confluence.ecmwf.int/display/ECC) (Fortran interface) | `--enable-emvorado` | `FCFLAGS='-I/path/to/eccodes/include' LDFLAGS='-L/path/to/eccodes/lib' LIBS='-leccodes_f90'` |
| [ECCODES](https://confluence.ecmwf.int/display/ECC) (C interface) | `--enable-grib2 --without-external-cdi` | `CPPFLAGS='-I/path/to/eccodes/include' LDFLAGS='-L/path/to/eccodes/lib' LIBS='-leccodes'` (or `LIBS='-lgrib_api'`) |
| [YAXT](https://gitlab.dkrz.de/dkrz-sw/yaxt) (Fortran interface) | `--enable-yaxt --with-external-yaxt` or `--enable-cdi-pio --with-external-yaxt` | `FCFLAGS='-I/path/to/yaxt/include' LDFLAGS='-L/path/to/yaxt/lib' LIBS='-lyaxt'` |
| [YAXT](https://gitlab.dkrz.de/dkrz-sw/yaxt) (C interface) | `--enable-cdi-pio --without-external-cdi --with-external-yaxt` or `--enable-coupling --with-external-yaxt` | `CPPFLAGS='-I/path/to/yaxt/include' LDFLAGS='-L/path/to/yaxt/lib' LIBS='-lyaxt_c'` |
| [SCT](https://gitlab.dkrz.de/dkrz-sw/sct) (Fortran interface) | `--enable-sct --with-external-sct` | `FCFLAGS='-I/path/to/sct/include' LDFLAGS='-L/path/to/sct/lib' LIBS='-lsct'` |
| [RTTOV](https://www.nwpsaf.eu/site/software/rttov/) | `--enable-rttov` | `FCFLAGS='-I/path/to/rttov/include -I/path/to/rttov/mod' LDFLAGS='-L/path/to/rttov/lib' LIBS='-lrttov_other -lrttov_emis_atlas -lrttov_brdf_atlas -lrttov_parallel -lrttov_coef_io -lrttov_hdf -lrttov_main'` |
| [LAPACK](http://www.netlib.org/lapack/) (or analogue) | mandatory | `LDFLAGS='-L/path/to/lapack/lib' LIBS='-llapack'` (depends on the implementation) |
| [BLAS](http://www.netlib.org/blas/) (or analogue) | mandatory | `LDFLAGS='-L/path/to/blas/lib' LIBS='-lblas'` (depends on the implementation) |
| [ECRAD](https://confluence.ecmwf.int/display/ECRAD/ECMWF+Radiation+Scheme+Home) | `--enable-ecrad --with-external-ecrad` | `FCFLAGS='-I/path/to/ecrad/include' LDFLAGS='-L/path/to/ecrad/lib' LIBS='-lradiation -lifsrrtm -lutilities -lifsaux'` |
| [RTE+RRTMGP](https://github.com/earth-system-radiation/rte-rrtmgp) | `--enable-rte-rrtmgp --with-external-rte-rrtmgp` | `FCFLAGS='-I/path/to/rte-rrtmgp/include' LDFLAGS='-L/path/to/rte-rrtmgp/lib' LIBS='-lrrtmgp -lrte'` |
| [NetCDF-Fortran](https://docs.unidata.ucar.edu/netcdf-fortran/current/) | mandatory | `FCFLAGS='-I/path/to/netcdf-fortran/include' LDFLAGS='-L/path/to/netcdf-fortran/lib' LIBS='-lnetcdff'` |
| [NetCDF-C](https://docs.unidata.ucar.edu/netcdf-c/current/) | `--without-external-cdi` or `--enable-coupling` | `CPPFLAGS='-I/path/to/netcdf/include' LDFLAGS='-L/path/to/netcdf/lib' LIBS='-lnetcdf'` |
| [HDF5](https://www.hdfgroup.org/solutions/hdf5/) (low- and high-level Fortran interfaces) | `--enable-emvorado` or `--enable-rttov` | `FCFLAGS='-I/path/to/hdf5/include' LDFLAGS='-L/path/to/hdf5/lib' LIBS='-lhdf5_hl_fortran -lhdf5_fortran'` |
| [HDF5](https://www.hdfgroup.org/solutions/hdf5/) (low-level C interface) | `--enable-sct --without-external-sct` | `CPPFLAGS='-I/path/to/hdf5/include' LDFLAGS='-L/path/to/hdf5/lib' LIBS='-lhdf5'` |
| [ZLIB](https://zlib.net) | `--enable-emvorado` | `LDFLAGS='-L/path/to/zlib/lib' LIBS='-lz'`<sup><a href="#f2">2</a></sup><a name="f2back"> |
| [AEC](https://gitlab.dkrz.de/k202009/libaec) | static linking | `LDFLAGS='-L/path/to/aec/lib' LIBS='-lsz -laec'` |
| [MPI](https://www.mpi-forum.org) (Fortran interface) | `--enable-mpi` | `FC='/path/to/mpi/bin/mpif90'` or `FCFLAGS='-I/path/to/mpi/include' LDFLAGS='-L/path/to/mpi/lib' LIBS='-lmpifort -lmpi'` (depends on the implementation) |
| [MPI](https://www.mpi-forum.org) (C interface) | `--enable-coupling` or `--enable-yaxt --without-external-yaxt` or `--enable-mpi --enable-sct --without-external-sct` or `--enable-cdi-pio --without-external-cdi` | `CC=/path/to/mpi/bin/mpicc` or `CPPFLAGS='-I/path/to/mpi/include' LDFLAGS='-L/path/to/mpi/lib' LIBS='-lmpi'` (depends on the implementation) |
| [ROCm](https://www.amd.com/en/products/software/rocm.html) | `--enable-gpu=openacc+hip` | `LDFLAGS='-L/path/to/rocm/lib' LIBS='-lamdhip64'` (depends on the platform) |
| [CUDA](https://developer.nvidia.com/cuda-zone) | `--enable-gpu=openacc+cuda` | `LDFLAGS='-L/path/to/cuda/lib' LIBS='-lcudart'` |
| [STDC++](https://isocpp.org/)<sup><a href="#f3">3</a></sup><a name="f3back"> | `--enable-gpu=openacc+cuda` or `--enable-gpu=openacc+hip` or `--enable-comin` | `LDFLAGS='-L/path/to/gcc/used/by/CUDACXX-or-HIPCXX-or-COMIN/lib' LIBS='-lstdc++'` (depends on the implementation) |
:::

1. <a name="f1"/>The dependency conditions and required flags are specified assuming that the shared versions of the libraries containing `RPATH` entries pointing to their dependencies are used (see section [](ref_buildrun_configuration_dynlibs)).<a href="#f1back">{octicon}`undo;1em;pst-color-secondary`</a>
2. <a name="f2"/>ZLIB is used via the `ISO_C_BINDING` interface and does not require additional preprocessor flags.<a href="#f2back">{octicon}`undo;1em;pst-color-secondary`</a>
3. <a name="f3"/>The provided standard C++ library must be compatible with the code generated by `CXX` and/or the host compiler of `CUDACXX`/`HIPCXX`.<a href="#f3back">{octicon}`undo;1em;pst-color-secondary`</a>

(ref_buildrun_configuration_bundled)=
## Bundled packages

As it was mentioned in the previous section, some of the packages are bundled together with the ICON source code. However, users can download and install those packages before configuring ICON and use them instead. This is controlled by the `--with-external-<package>` arguments of the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON. The arguments accept either `yes` or `no`. If the usage of an external version of a package is requested (i.e. `--with-external-<package>=yes`), the compiler and linker flags, i.e. `FCFLAGS`, `CPPFLAGS`, `LDFLAGS` and `LIBS`, are supposed to be extended accordingly (see [Table 1](tab_icon_depgraph)). The {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON fails if the flags are not set correctly or enable a version of the package that is known to be incompatible with ICON.

By default, the bundled versions of the packages are used. In this case, the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON runs the configure and CMake scripts of the packages and extends the compiler and linker flags automatically. The arguments that are passed to the configure and CMake scripts of the bundled packages are composed based on the arguments provided to the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON as follows:
- [Autotools](https://en.wikipedia.org/wiki/GNU_Autotools)-based packages:
  - by default, the arguments are passed unchanged, which means that if you need to give an additional argument to the configure script of a bundled package, you can specify it when calling the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON, even though, the argument is not listed in its help message;
  - arguments that potentially break the configuration and building of ICON are filtered out (see the expansion of the `ACX_SUBDIR_REMOVE_ARGS` macro in {{ '[`configure.ac`]({}/configure.ac)'.format(base_url) }}), for example, the libraries of the bundled packages must be linked statically, therefore the argument `--enable-shared` is never passed to the configure scripts that support it;
  - the list of arguments is extended with ones that enforce consistent building (see expansion to the `ACX_SUBDIR_APPEND_ARGS` macro in {{ '[`configure.ac`]({}/configure.ac)'.format(base_url) }}), for example, the configure scripts of the packages receive additional arguments `--disable-shared` and `--enable-static`;
  - compiler flags are modified as described in section [](ref_buildrun_configuration_flags).
- [CMake](https://cmake.org)-based packages:
  - the packages are configured out-of-source using the `CMAKE` command (see section [](ref_buildrun_configuration_compilers)) with the [`Unix Makefiles`](https://cmake.org/cmake/help/v3.18/generator/Unix%20Makefiles.html) specified as the generator;
  - the list of arguments is package-specific (see expansion to the `ACX_SUBDIR_APPEND_ARGS` macro in {{ '[`configure.ac`]({}/configure.ac)'.format(base_url) }}) but the shared libraries and tests are generally disabled (`-DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=OFF`);
  - compiler flags that are passed to the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON are modified as described in section [](ref_buildrun_configuration_flags) and passed as the respective `-DCMAKE_<LANG>_FLAGS` arguments;
  - compiler commands (see section [](ref_buildrun_configuration_compilers)) are passed via the respective [environment variables](https://cmake.org/cmake/help/latest/manual/cmake-env-variables.7.html#environment-variables-for-languages).

:::{admonition} Configuration Errors
:class: admonition-icontheme
Oftentimes, error and warning messages that are printed at configure time are emitted not by the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON but by the configure scripts of the bundled packages. Each of the scripts generates its own `config.log` file, which can help in better understanding of the reported issue. The log files are put to the corresponding subdirectories of the `./externals` directory residing in the root build directory of ICON.
:::

(ref_buildrun_configuration_compilers)=
## Compilers and tools

Compilers and tools to be used for building are read by the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON from the following environment variables (can be passed as the command-line arguments):
- `FC` &mdash; Fortran compiler command;
- `CC` &mdash; C compiler command;
- `CXX` &mdash; C++ compiler command (used by some of the bundled packages);
- `CUDACXX` &mdash; [CUDA C++ compiler](https://developer.nvidia.com/cuda-llvm-compiler) command (used only when the GPU support with CUDA is enabled);
- `HIPCXX` &mdash; [HIP C++ compiler](https://github.com/ROCm/llvm-project/tree/amd-staging/amd/hipcc) command (used only when the GPU support with HIP is enabled);
- `AR` &mdash; archiver command (used to create static libraries);
- `RANLIB` &mdash; archive indexer command (used to create static libraries);
- `PYTHON` &mdash; [Python](https://www.python.org/) interpreter command;
- `CMAKE` &mdash; [CMake](https://cmake.org/) command (used to configure some of the bundled packages);
- `FPP` &mdash; Fortran preprocessor command (used when explicit Fortran preprocessing is enabled, see section [](ref_buildrun_building_preprocessing) for more details), must treat the first positional command-line argument as the path to the input source file and print the result to the standard output stream;
- `SB2PP` &mdash; [Serialbox2](https://gridtools.github.io/serialbox/) preprocessor command (used when the Serialbox2 serialization is enabled, see section [](ref_buildrun_building_preprocessing) for more details);
- `MPI_LAUNCH` &mdash; interactive (synchronous) MPI launcher command (used by the bundled packages for configure-time checks).

If the variables are set, the {{ '[`configure`]({}/configure)'.format(base_url) }} script checks whether their values meet the requirements, otherwise, the script tries to guess suitable values for them. Thus, if you want to make sure that a particular command for a particular operation is used, you need to specify the corresponding variable explicitly. For example, the usage of NAG compiler is enforced with the following additional command-line argument of the {{ '[`configure`]({}/configure)'.format(base_url) }} script:

```sh
./configure FC=nagfor <other arguments>
```

(ref_buildrun_configuration_flags)=
## Compiler flags

The configure script supports several groups of compiler flags. Each group is associated with one of the following environment variables (can be passed as the command-line arguments):
- `FCFLAGS` &mdash; Fortran compiler flags to be used when configuring, compiling and linking ICON, as well as passed to the configure and CMake (via the `-DCMAKE_Fortran_FLAGS` argument) scripts of the bundled packages (in contrast to standard Autoconf-based scripts, the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON does not set `FCFLAGS` to `-g -O2` by default);
- `ICON_FCFLAGS` &mdash; Fortran compiler flags to be appended to `FCFLAGS` when configuring, compiling and linking ICON;
- `ICON_<NAME>_FCFLAGS` &mdash; Fortran compiler flags to be appended to `FCFLAGS` instead of `ICON_FCFLAGS` when compiling files of the [Fortran compile group](ref_buildrun_configuration_compilegroups) `<NAME>` (defaults to `ICON_FCFLAGS`, which can be overridden by setting the variable to an empty value: `ICON_OCEAN_FCFLAGS=`);
- `ICON_BUNDLED_FCFLAGS` &mdash; Fortran compiler flags to be appended to `FCFLAGS` when configuring the bundled packages (defaults to `ICON_FCFLAGS`, which can be overridden by setting the variable to an empty value: `ICON_BUNDLED_FCFLAGS=`);
- `ICON_RTE_RRTMGP_FCFLAGS`, `ICON_ECRAD_FCFLAGS`, etc. &mdash; Fortran compiler flags to be appended `FCFLAGS` when configuring the respective bundled packages (defaults to `ICON_BUNDLED_FCFLAGS`, which can be overridden by setting the variables to empty values: `ICON_RTE_RRTMGP_FCFLAGS=`, `ICON_ECRAD_FCFLAGS=`, etc.);
- `CFLAGS` &mdash; C compiler flags to be used when configuring and compiling ICON, as well as passed to the configure and CMake (via the `-DCMAKE_C_FLAGS` argument) scripts of the bundled packages (in contrast to standard Autoconf-based scripts, the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON does not set `FCFLAGS` to `-g -O2` by default);
- `CPPFLAGS` &mdash; C preprocessor flags to be used when configuring and compiling ICON, as well as passed to the configure and CMake (via the `-DCMAKE_C_FLAGS` and `-DCMAKE_CXX_FLAGS` arguments) scripts of the bundled packages;
- `ICON_CFLAGS` &mdash; C compiler flags to be appended to `CFLAGS` when configuring and compiling ICON;
- `ICON_BUNDLED_CFLAGS` &mdash; C compiler flags to be appended to `CFLAGS` when configuring the bundled packages (defaults to `ICON_CFLAGS`, which can be overridden by setting the variable to an empty value: `ICON_BUNDLED_CFLAGS=`);
- `ICON_CDI_CFLAGS`, `ICON_MTIME_CFLAGS`, etc. &mdash; C compiler flags to be appended to `CFLAGS` when configuring the respective bundled packages (defaults to `ICON_BUNDLED_CFLAGS`, which can be overridden by setting the variables to empty values: `ICON_CDI_CFLAGS=`, `ICON_MTIME_CFLAGS=`, etc.);
- `CXXFLAGS` &mdash; C++ compiler flags to be passed to the configure and CMake (via the `-DCMAKE_CXX_FLAGS` argument) scripts of the bundled packages (in contrast to standard Autoconf-based scripts, the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON does not set `CXXFLAGS` to `-g -O2` by default);
- `ICON_BUNDLED_CXXFLAGS` &mdash; C++ compiler flags to be appended to `CXXFLAGS` when configuring the bundled packages;
- `CUDAFLAGS` &mdash; CUDA C++ compiler flags to be used when configuring and compiling ICON;
- `HIPFLAGS` &mdash; HIP C++ compiler flags to be used when configuring and compiling ICON;
- `LDFLAGS` &mdash; common Fortran, C and C++ compiler flags to be used when configuring and linking ICON, as well as passed to the configure and CMake (via the `-DCMAKE_EXE_LINKER_FLAGS`, `-DCMAKE_MODULE_LINKER_FLAGS` and `-DCMAKE_SHARED_LINKER_FLAGS` arguments) scripts and of the bundled packages;
- `ICON_LDFLAGS` &mdash; Fortran compiler flags to be appended to `LDFLAGS` when configuring and linking ICON;
- `LIBS` &mdash; a list of libraries (see [Table 1](tab_icon_depgraph) for the recommended order) to be passed to the linker by the Fortran compiler when linking ICON and to the configure and CMake (via the `-DCMAKE_Fortran_STANDARD_LIBRARIES`, `-DCMAKE_C_STANDARD_LIBRARIES` and `-DCMAKE_CXX_STANDARD_LIBRARIES` arguments) scripts of the bundled packages.

The general recommendation to follow when composing the flags is:
1. Flags specifying search paths for header and module files, i.e. the `-I<path>` flags, should be specified as `FCFLAGS` and `CPPFLAGS`, depending on whether they need to be passed to Fortran or C compiler, respectively.
2. By default, other flags, e.g. the optimization ones, that are meant for Fortran and C compilers should be appended to `FCFLAGS` and `CFLAGS`, respectively.
3. Fortran and C compiler flags that need to be used when configuring, compiling and linking ICON but at the same time can break the configuration (a flag is too restrictive for the configure checks to pass, e.g. `-fimplicit-none` for Gfortran) or the functionality of the bundled packages (e.g. the optimization level required for ICON is too high and leads to errors in the functionality of the bundled packages) can be put to `ICON_FCFLAGS`, `ICON_CFLAGS` or `ICON_LDFLAGS`.
4. Special optimization flags for a selected set of Fortran source files of ICON can be put to `ICON_<NAME>_FCFLAGS`, where `<NAME>` is the name of a [Fortran compile group](ref_buildrun_configuration_compilegroups).
5. Fortran and C compiler flags that need to be used when compiling and linking the bundled packages but at the same time conflict with the flags required for ICON (e.g. you want to compile ICON with `-O3` flag but the bundled packages
need to be compiled with `-O2`) can be specified as `ICON_BUNDLED_FCFLAGS` and `ICON_BUNDLED_CFLAGS`, respectively.
6. If a set of Fortran (or C) compiler flags needs to be passed only to some particular bundled package, it can be specified in the respective variable, e.g. in `ICON_CDI_FCFLAGS` (or `ICON_CDI_CFLAGS`).

(ref_buildrun_configuration_compilegroups)=
### Fortran compile groups

Certain Fortran source files of ICON, including its components that do not have their own build system (e.g. JSBACH) might require special compiler flags. Either due to a decision made by the developers (e.g. the ocean component) or due to compiler bugs. Such cases can be covered with Fortran compile groups. For example, a separate set of flags for the ocean component can be specified at the configure time as follows:

```sh
./configure \
  --enable-fcgroup-OCEAN=src/hamocc:src/ocean:src/sea_ice \
  ICON_OCEAN_FCFLAGS=<special ocean component flags> \
  <other arguments>
```

Or alternatively:

```sh
./configure \
  --enable-fcgroup-OCEAN \
  ICON_OCEAN_PATH=src/hamocc:src/ocean:src/sea_ice \
  ICON_OCEAN_FCFLAGS=<special ocean component flags> \
  <other arguments>
```

(ref_buildrun_configuration_dynlibs)=
## Dynamic libraries

For each `-L<path>` flag found in the `LDFLAGS` and `LIBS` variables, the {{ '[`configure`]({}/configure)'.format(base_url) }} script of ICON generates an additional linker flag that puts the `<path>` on the list of runtime library search paths of the ICON executable. This allows for the automatic location of the required libraries by the *dynamic linker* at the runtime. The flags are appended to `LDFLAGS` at the configure time and their actual form depends on the Fortran compiler in use. By default, the flags are composed using the template `-Wl,-rpath -Wl,<path>` with currently the only exception for NAG compiler, which accepts the flags in the form `-Wl,-Wl,,-rpath -Wl,-Wl,,<path>`. If the `-rpath` flags generated by the {{ '[`configure`]({}/configure)'.format(base_url) }} script break the building or you perform a completely static linking, you can disable the feature with the `--disable-rpaths` argument.


:::{admonition} Notes on RPATH
:class: admonition-icontheme
- The GNU Linker (GNU ld) implements a feature called `new-dtags`. If this feature is enabled (usually by default), the linker treats the `-rpath <path>` flag differently: instead of setting the `DT_RPATH` attribute of the output shared library file to the `<path>`, it sets the `DT_RUNPATH` attribute of the file to the same value. This alters the way the dynamic loader locates the required dynamic libraries at the runtime: if the dynamic loader finds a `DT_RUNPATH` attribute, it **ignores** the value of the `DT_RPATH` attribute (if any), with the effect that the `LD_LIBRARY_PATH` environment variable is checked first and the paths in the `DT_RUNPATH` attribute are only searched afterwards. Moreover, the dynamic loader does not search `DT_RUNPATH` locations for transitive dependencies, unlike `DT_RPATH`. Therefore, it is important to keep in mind that:
  1. The ICON executable is not agnostic to the environment if it has been linked with the `new-dtags` feature enabled: the `LD_LIBRARY_PATH` environment variable can override the rpath entries set by the linker. The feature can be disabled by appending the `-Wl,--disable-new-dtags` flags to the `LDFLAGS` variable.
  2. If an immediate ICON dependency, e.g. `libnetcdf.so`, has at least one `DT_RUNPATH` entry but none of them points to a directory containing one of the libraries required by that dependency, e.g. `libsz.so`, the dynamic loader will not be able to locate the latter at the runtime even if the ICON executable has been linked without the `new-dtags` feature and an `-Wl,-rpath` flag pointing to the right location, e.g. `-Wl,--disable-new-dtags -Wl,-rpath -Wl,/path/to/libsz`. A possible workaround for this is to make the secondary dependency of ICON a primary one by *overlinking* (i.e. linking to a library, which is not used by the executable directly) to it, e.g. by adding the `-lsz` to the `LIBS` variable. This way, the dependency will become a non-transitive one and the dynamic loader will be able to locate it using either `DT_RUNPATH` or `DT_RPATH` entries of the ICON executable.
For more details, refer to the man pages of the linker (`man ld`) and the dynamic loader (`man ld.so`).

- Due to a significant maintenance overhead (see commit [36ab00dd]( https://gitlab.dkrz.de/icon/icon/-/commit/36ab00dd5641419e5afbb918b47cdf697c54b737#87db583be5c13c1f7b3c958b10e03d67b6a2ca06)), the generated `-rpath` flags are **not** passed to the configure scripts of the bundled packages. However, some of their checks imply running executables linked with the flags listed in the `LIBS` variable. To prevent those checks from false negative results, which oftentimes are reported with misleading messages, the dynamic loader needs to be able to locate all the libraries referenced in the `LIBS` variable. A way to achieve that is to list paths to the libraries in the `LD_LIBRARY_PATH` variable and export it, i.e. run `export LD_LIBRARY_PATH="<path1>:<path2>:..."` before running the {{ '[`configure`]({}/configure)'.format(base_url) }} script.
:::

:::{admonition} Libtool
:class: admonition-icontheme
Some of the bundled packages employ [Libtool](https://www.gnu.org/software/libtool/), which is known to be **not** fully compatible with some compilers. For example, the flags `-Wl,-Wl,,-rpath -Wl,-Wl,,<path>`, which are valid for NAG compiler, are incorrectly transformed by Libtool into `-Wl,-Wl -Wl,"" -Wl,-rpath -Wl,-Wl -Wl,"" -Wl,<path>`. A possible solution for this problem is to add the flags in the form understood by NAG compiler not to `LDFLAGS` but to `ICON_LDFLAGS`.
:::

(ref_buildrun_configuration_envs)=
## Configuration and building environments

It is important that both the configuration and the building stages are performed in the same environment, i.e. the environment variables that might influence the way the compilers and the linker work are set to the same values when running the {{ '[`configure`]({}/configure)'.format(base_url) }} script and when running the `make` command for building. For example, NAG compiler will not work if the environment variable `NAG_KUSARI_FILE` is not set properly. Keeping track of all the steps required to re-initialize the environment for the building stage, e.g. when the configuration stage has already been done but in another terminal session, might be challenging, especially in HPC environments offering multiple compilers and libraries.

To ensure consistency between configuration and building environments, you can set the `BUILD_ENV` argument of the {{ '[`configure`]({}/configure)'.format(base_url) }} script with a series of shell commands that initialize the necessary environment variables. If the `BUILD_ENV` variable is not empty, the {{ '[`configure`]({}/configure)'.format(base_url) }} script will run the commands it contains before running any checks. Additionally, the commands will be saved to the `Makefile`, so they will be executed each time the `make` command is launched for building. The shell script that is provided as value of the `BUILD_ENV` argument must be a one-liner ending with a semicolon (;) symbol, for example:

```sh
./configure BUILD_ENV='. /etc/profile.d/modules.sh; module purge; module load intel;' <other arguments>
```

Also, a proper implementation of the `BUILD_ENV` script allows for switching between multiple build directories (see section [](ref_buildrun_configuration_oos)) without having to re-initialize the environment accordingly.

(ref_buildrun_configuration_wrappers)=
## Configure wrappers

Real-case configuration commands might be rather long and complex. For example, below, you can find an example of the configuration command for [Levante@DKRZ](https://docs.dkrz.de/doc/levante/index.html):

```sh
./configure \
  AR=xiar \
  BUILD_ENV=". ./config/dkrz/module_switcher; \
             switch_for_module \
               intel-oneapi-compilers/2022.0.1-gcc-11.2.0 \
               openmpi/4.1.2-intel-2021.5.0;" \
  CC=mpicc \
  CFLAGS="-g -gdwarf-4 -qno-opt-dynamic-align -m64 -march=core-avx2 \
          -mtune=core-avx2 -fma -ip -pc64 -std=gnu99" \
  CPPFLAGS="-I/sw/spack-levante/hdf5-1.12.1-tvymb5/include \
            -I/sw/spack-levante/netcdf-c-4.8.1-2k3cmu/include \
            -I/sw/spack-levante/eccodes-2.21.0-3ehkbb/include \
            -I/usr/include/libxml2" \
  FC=mpif90 \
  FCFLAGS="-I/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/include \
           -m64 -march=core-avx2 -mtune=core-avx2 -g -gdwarf-4 -pc64 \
           -fp-model source" \
  ICON_BUNDLED_CFLAGS='-O2 -ftz' \
  ICON_CDI_CFLAGS='-O2 -ftz' \
  ICON_CFLAGS='-O3 -ftz' \
  ICON_ECRAD_FCFLAGS='-qno-opt-dynamic-align -no-fma -fpe0' \
  ICON_FCFLAGS="-DDO_NOT_COMBINE_PUT_AND_NOCHECK -O3 -ftz \
                -qoverride-limits -assume realloc_lhs \
                -align array64byte -fma -ip \
                -D__SWAPDIM -DOCE_SOLVE_OMP" \
  ICON_OCEAN_FCFLAGS="-O3 -assume norealloc_lhs -reentrancy threaded \
                      -qopt-report-file=stdout -qopt-report=0 \
                      -qopt-report-phase=vec" \
  ICON_YAC_CFLAGS='-O2 -ftz' \
  LDFLAGS="-L/sw/spack-levante/hdf5-1.12.1-tvymb5/lib \
           -L/sw/spack-levante/netcdf-c-4.8.1-2k3cmu/lib \
           -L/sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g/lib \
           -L/sw/spack-levante/netlib-lapack-3.9.1-rwhcz7/lib64 \
           -L/sw/spack-levante/eccodes-2.21.0-3ehkbb/lib64" \
  LIBS="-Wl,--disable-new-dtags -Wl,--as-needed -lxml2 -leccodes \
        -llapack -lblas -lnetcdff -lnetcdf -lhdf5" \
  MPI_LAUNCH=mpiexec \
  --enable-intel-consistency \
  --enable-vectorized-lrtm \
  --enable-parallel-netcdf \
  --enable-grib2 \
  --enable-fcgroup-OCEAN=src/hamocc:src/ocean:src/sea_ice \
  --enable-yaxt \
  --enable-art \
  --enable-ecrad \
  --disable-mpi-checks
```

Repeatedly composing such a command is exhausting and error-prone. Therefore, each team involved in the development of ICON is encouraged to implement and maintain *configure wrappers*, which would simplify the configuration stage for their users. The wrappers should be put in a subdirectory with a relevant name of the {{ '[`config`]({}/config)'.format(base_url) }} directory. Although there are no hard requirements on how the scripts should be implemented, we recommend considering the following features:

1. Ensure that the wrapper script passes command-line arguments to the {{ '[`configure`]({}/configure)'.format(base_url) }} script, allowing users to override default configuration options. For example:

    ```sh
    ./config/dkrz/levante.intel --enable-openmp
    ```
    Also, account for the case of calling the wrapper with the `--help` argument.
2. Account for out-of-source building:
    - the wrapper script should be able to locate the {{ '[`configure`]({}/configure)'.format(base_url) }} script when called from a directory other than the source root directory, for example:
        ```sh
        script_dir=$(cd "$(dirname "$0")"; pwd)
        icon_dir=$(cd "${script_dir}/../.."; pwd)
        "${icon_dir}/configure" <a list of predefined arguments> "$@"
        ```
    - the wrapper script should prepare the current working directory for the following [runscript generation](ref_buildrun_mkexp).
3. Prepend the `LIBS` variable with `-Wl,--as-needed` flag so that the actual list of the libraries the ICON executable depends on would include only those required for the particular configuration of the model (see the man pages of the linker for more details: `man ld`).
4. Allow for running the `make check` command (see section [](ref_buildrun_building_bundled) for more details) by extending the `LD_LIBRARY_PATH` environment variable in the `BUILD_ENV` script, instead of doing so in the wrapper script itself.

(ref_buildrun_configuration_wrappersgeneric)=
## Generic configure wrapper

The subfolder {{ '[`config/generic`]({}/config/generic)'.format(base_url) }} contains generic [configure wrappers](ref_buildrun_configuration_wrappers). These generic configure wrappers assume that all required [](ref_buildrun_configuration_wrappersgeneric_sl) are installed under the same prefix. The prefix defaults to `/opt/local` with a fallback to `/opt/homebrew` on macOS and to `/usr` on other platforms. The default values can be overridden by setting the environment variable `ICON_SW_PREFIX`:

```sh
export ICON_SW_PREFIX='/path/to/icon/prerequisites'
```

The following sections provide a list of software required for building and running ICON. Users can build and install (to the same prefix) the listed packages manually or use the [package managers](https://en.wikipedia.org/wiki/Package_manager) available for their platform. Basic instructions on how to do it on several popular platforms are provided in section [](ref_buildrun_hardware).

(ref_buildrun_configuration_wrappersgeneric_bt)=
### Building tools

- [GNU Make](https://www.gnu.org/software/make) v3.81+
- [CMake](https://cmake.org) v3.18+
- [Python](https://www.python.org) v3.9+
- [Perl](https://www.perl.org) v5.10+
- Interoperable C, CXX and Fortran compilers

(ref_buildrun_configuration_wrappersgeneric_sl)=
### Software libraries

- [MPICH](https://www.mpich.org), [OpenMPI](https://www.open-mpi.org) or any other [MPI](https://www.mpi-forum.org) implementation that provides compiler wrappers `mpicc`, `mpicxx` and `mpif90` for C, C++ and Fortran, respectively, as well as the job launcher `mpiexec`

:::{admonition} Notes on MPI libraries
:class: admonition-icontheme
- The job launcher of OpenMPI fails to run more MPI processes than the number of real processor cores available on the machine by default. That might lead to failures when configuring or running ICON. The solution to the problem is to run the configure wrapper with an additional argument `MPI_LAUNCH='mpiexec --oversubscribe'` (alternatively, you can set the `OMPI_MCA_rmaps_base_oversubscribe` environment variable to `1`).
- It is not rare that the latest versions (or the default versions available via the package managers) of OpenMPI and MPICH are affected with bugs that make the libraries unusable for ICON. A way to make sure that the MPI library does not have significant defects is to switch to the root source directory of ICON and run the following commands (do not forget the aforementioned extra arguments for the configure wrapper if you are using OpenMPI):
```sh
./config/generic/gcc --enable-yaxt --enable-cdi-pio --enable-coupling
make -j4 check-bundled TESTS=  # this step speeds up the next one but can be skipped
make check-bundled  # avoid running this step in parallel on a weak machine, i.e. omit the -j argument
test $? -eq 0 && echo "Everything is fine" || echo "Something went wrong"
```
- After that, you can clean up the source directory and reconfigure ICON the way you need.
:::

- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) with high-level interface (for <a href="#netcdf-c">NetCDF-С</a>), thread-safety (for <a href="#cdo">CDO</a>), and szlib filtering support (only C interface required, not a direct dependency of ICON)
- <a name="netcdf-c"/> [NetCDF-C](https://docs.unidata.ucar.edu/netcdf-c/current/) with NetCDF-4 support
- [NetCDF-Fortran](https://docs.unidata.ucar.edu/netcdf-fortran/current/)
- [BLAS](http://www.netlib.org/blas)
- [LAPACK](http://www.netlib.org/lapack)
- [ecCodes](https://confluence.ecmwf.int/display/ECC) with JPEG2000 and AEC support (only C interface required)
- [libfyaml](https://github.com/pantoniou/libfyaml)
- [Libxml2](http://www.xmlsoft.org)

See section [ICON dependencies](ref_buildrun_configuration_icondep) for more details.

(ref_buildrun_configuration_wrappersgeneric_ot)=
### Optional tools

- <a name="cdo"/> [CDO](https://code.mpimet.mpg.de/projects/cdo) for pre- and post-processing, also used by some of the [generated runscripts](ref_buildrun_mkexp)
- [rsync](https://rsync.samba.org/) for the generated runscipts in the case of [out-of-source building](ref_buildrun_configuration_oos)

(ref_buildrun_configuration_oos)=
## Out-of-source configuration (building)

The building system of ICON supports so-called *out-of-source* builds. This means that you can build ICON in a directory other than the *source* root directory. The main advantage of this is that you can easily switch between several different configurations (each in its own *build* directory) of the model while working on the same source code. Hence, it is possible to introduce changes into the source code and test them with different compilers, flags and options without having to re-configure the model or copy the updated source files to other directories.

The {{ '[`configure`]({}/configure)'.format(base_url) }} script (also when called via a configure wrapper) prepares the *current working directory* for the following building. The particular case of the current working directory being the source root directory is called *in-source* configuration.


:::{admonition} make distclean
:class: admonition-icontheme
It is not allowed to mix in-source and out-of-source builds: the source directory must be cleaned (`make distclean`) from the files generated as a result of a prior in-source configuration and building before an out-of-source configuration can take place.
:::

The following example shows how ICON can be configured on [Levante@DKRZ](https://docs.dkrz.de/doc/levante/index.html) in two different directories using Intel and GCC compiler toolchains (assuming that the source root directory of ICON is `/path/to/icon-srcdir`):

```sh
mkdir intel && cd intel
/path/to/icon-srcdir/config/dkrz/levante.intel
cd ..
mkdir gcc && cd gcc
/path/to/icon-srcdir/config/dkrz/levante.gcc
cd..
```

After executing the commands above, you will have two directories named `intel` and `gcc`, each ready for building ICON with Intel and GCC compilers, respectively.

(ref_buildrun_building)=
# Building

The building stage is done with [GNU make](https://www.gnu.org/software/make/) upon successful completion of the configuration stage. The oldest supported version of `make` is **3.81**, however, it has significant limitations and it is recommended to use version **4.1** or later.

The building step is done by running `make` command with an optional argument specifying the number of jobs to run simultaneously. For example:

```sh
make -j8
```

:::{admonition} Basic targets
:class: admonition-icontheme
Users are recommended to get familiar with the list of basic *targets*, i.e. supported subcommands, by running:
```sh
make help
```
:::

By default, `make` reads the instructions from a file called `Makefile`, which in the case of ICON is generated based on the {{ '[`Makefile.in`]({}/Makefile.in)'.format(base_url) }} template and is mainly responsible for the initialization of the building environment. The initialization is performed by the shell script provided to the {{ '[`configure`]({}/configure)'.format(base_url) }} script as the `BUILD_ENV` argument (see section [](ref_buildrun_configuration_envs)). The script is saved in the `Makefile`. Each time `make` is executed in the root of the ICON build directory, it reads the `Makefile`, runs the initialization script and [recursively](https://www.gnu.org/software/make/manual/html_node/Recursion.html) executes itself in the same directory but with another input makefile called `icon.mk`. The latter is generated based on the {{ '[`icon.mk.in`]({}/icon.mk.in)'.format(base_url) }} template and contains the instructions on how to perform the following building steps required to generate the ICON executable.


:::{admonition} Revert modifications of automatically generated files
:class: admonition-icontheme
All modifications of the makefiles and other files that are automatically generated at the configuration stage can be reverted by calling the `./config.status` script residing in the build directory.
:::

The following sections provide information on some features and implementation details of the building process.

:::{admonition} Verbose make
:class: admonition-icontheme
By default, `make` produces very minimalistic output to the standard output stream. This allows for better recognition of the warning messages emitted by the compilers and other tools. This can be altered either at the configuration stage by calling the {{ '[`configure`]({}/configure)'.format(base_url) }} script with an additional option `--disable-silent-rules` or at the building stage by calling `make` with an additional command-line argument `V`, which can be set either to `1` (enable verbose output) or to `0` (disable verbose output). For example, if you want to see the exact commands executed by `make`, you can run:
```sh
make V=1 <other arguments>
```
:::

(ref_buildrun_building_collection)=
## Source file collection

The list of source files that need to be compiled to produce the ICON executable is generated dynamically each time the `make` command is executed. This is done using `find` command (both [GNU](https://www.gnu.org/software/findutils/) and [BSD](https://www.freebsd.org/cgi/man.cgi?find(1)) versions are supported) with the assumption that the source files have the following filename extensions:

- `.f90` &mdash; Fortran source files, regardless of whether they contain Fortran preprocessor directives;
- `.inc` &mdash; Fortran header files included with the quoted form of the Fortran preprocessor `#include` directives, e.g. `#include "filename.inc"`;
- `.incf` &mdash; Fortran header files included with the Fortran `INCLUDE` statements, e.g. `include 'filename.incf'` (these files are not allowed to have Fortran preprocessor directives);
- `.c` &mdash; C source files;
- `.cu` &mdash; CUDA source files;
- `.hip.cc` &mdash; HIP source files.

The list of source files is a result of the recursive search for files that have the aforementioned extensions and reside in the `src` and `support` subdirectories of the source root directory of ICON. Additionally, depending on whether the corresponding components of the model were enabled at the configuration stage, the list is extended with Fortran source files from the `./externals/jsbach/src`, `./externals/dace_icon/src_for_icon`, `./externals/emvorado` and `./externals/art` subdirectories.

:::{admonition} Notes on source file collection
:class: admonition-icontheme
- In general, you can extend the source base of ICON just by adding the source files to the `src` subdirectory of the source root directory of ICON.
- Files and directories with names starting with a dot (`.`), as well as files containing the common preprocessor infix `.pp-` (see section [](ref_buildrun_building_preprocessing) for more details) in their names, are ignored.
:::

(ref_buildrun_building_preprocessing)=
## Preprocessing

Depending on the configuration, Fortran source files undergo one or more of the following preprocessing procedures.
1. Fortran source files residing in the `./externals/jsbach/src` are preprocessed with the [`dsl4jsb.py`](https://gitlab.dkrz.de/jsbach/jsbach/blob/master/scripts/dsl4jsb/dsl4jsb.py) script. This is done only if the JSBACH component has been enabled at the configuration stage (`--enable-jsbach`). Otherwise, the source files of the component are completely ignored. The output files of this procedure have the same name as the input files plus an additional infix `.pp-jsb`.
2. Depending on whether the *explicit Fortran preprocessing* is enabled (`--enable-explicit-fpp`), the results of the **actual** previous preprocessing step, together with Fortran source files that have not been preprocessed yet, are preprocessed with the standard Fortran preprocessor. The output files of this procedure have the same name as the input files plus an additional infix `.pp-fpp`.
3. If the *Serialbox2 serialization* is enabled (`--enable-serialization`), the results of the **actual** previous preprocessing step, together with Fortran source files that have not been preprocessed yet, are preprocessed with the [corresponding script](https://github.com/GridTools/serialbox/blob/master/src/serialbox-python/pp_ser/pp_ser.py) of the [Serialbox2](https://gridtools.github.io/serialbox/) toolkit. The output files of this procedure have the same name as the input files plus an additional infix `.pp-sb2`.

:::{admonition} Explicit Fortran preprocessing
:class: admonition-icontheme
The explicit Fortran preprocessing is enabled automatically when the Serialbox2 serialization is enabled. You can override this by disabling the preprocessing with the `--disable-explicit-fpp` option.
:::

The output directories of the preprocessing steps have the same layout as the directories containing their input files and the output files have the same prefixes and suffixes (i.e. extensions) as the corresponding input files. For example, if the source file of JSBACH `./externals/jsbach/src/base/mo_jsb_base.f90` is preprocessed by each of the preprocessing steps, the corresponding output files are saved as follows:
- `./externals/jsbach/src/base/mo_jsb_base.pp-jsb.f90` &mdash; JSBACH preprocessing output;
- `./externals/jsbach/src/base/mo_jsb_base.pp-jsb.pp-fpp.f90` &mdash; explicit Fortran preprocessing output;
- `./externals/jsbach/src/base/mo_jsb_base.pp-jsb.pp-fpp.pp-sb2.f90` &mdash; Serialbox2 preprocessing output.

:::{admonition} Language-specific preprocessing
:class: admonition-icontheme
Source files are additionally preprocessed with the corresponding standard language-specific (i.e. Fortran, C, CUDA, HIP) preprocessors as part of the compilation process. In contrast to the procedures described in this section, the compilation-time preprocessing is non-optional and run by the compilers *implicitly*.
:::

(ref_buildrun_building_deptrack)=
## Source dependency tracking

Before the compilation can take place, it is required to identify the exact list of source files that need to be compiled to produce the ICON executable. The actual content of the list depends not only on how the code is [configured](ref_buildrun_configuration), but possible modifications of the source code made since the last call of `make` need to be taken into account as well. Moreover, in the case of Fortran source files, the compilation order becomes important since a source file declaring a Fortran module must be compiled before any other source file using that module. Both tasks are accomplished as follows.

Once the [](ref_buildrun_building_preprocessing) is finished, all source files of ICON (including the enabled components) or their **final** preprocessed versions are processed with the dependency generator ({{ '[`depgen.py`]({}/utils/mkhelper/depgen.py)'.format(base_url) }}). The tool parses each source file, detects which header and module files are required for its successful compilation, and stores this information in the form of a makefile. The makefiles are then read by `make` and the dependency listing script {{ '[`deplist.py`]({}/utils/mkhelper/deplist.py)'.format(base_url) }}. The former makes sure that the source files are compiled in the right order, and the latter identifies the list of source files that need to be compiled.

The dependency generator ({{ '[`depgen.py`]({}/utils/mkhelper/depgen.py)'.format(base_url) }}) recognizes preprocessor `#include`, `#if` and the associated directives as well as Fortran `INCLUDE`, `USE` and `MODULE` statements. If the usage of a module or a header file is surrounded by the `#ifdef SOME_MACRO` and `#endif` directives, it will be put on the list of files required for the compilation only if the macro `SOME_MACRO` is defined. The list of macro definitions enabling various features of the model is generated at the [configuration](ref_buildrun_configuration) stage in the form of compiler flags, e.g. `-DSOME_MACRO -DSOME_OTHER_MACRO`, which are appended to `FCFLAGS`. This way, some Fortran modules do not need to be generated in particular configurations of the model and, therefore, the source files declaring them are not compiled.

:::{admonition} Undetectable dependencies
:class: admonition-icontheme
Two types of source dependencies cannot be detected by the dependency generator ({{ '[`depgen.py`]({}/utils/mkhelper/depgen.py)'.format(base_url) }}). Undetectable dependencies of the first type are related to [C/Fortran interoperability](http://fortranwiki.org/fortran/show/C+interoperability): if a Fortran source file contains a declaration of a binding to a function defined in a C source file, the dependency of the respective object files must be reflected in the recipe for target `c_binding.d` of the {{ '[`icon.mk.in`]({}/icon.mk.in)'.format(base_url) }} template. The second type of undetectable dependencies is associated with Fortran external procedures: when a Fortran source file contains a call to an external procedure, i.e. a function or a subroutine that is not part of any Fortran module. Dependencies of the second type are not supported.
:::

The dependency listing script {{ '[`deplist.py`]({}/utils/mkhelper/deplist.py)'.format(base_url) }} reads the dependency makefiles, builds a source dependency graph and traverses it starting with the vertex associated with the `src/drivers/icon.o` object file. Each vertex of the graph accessible from the starting one is printed to the output. The output is then filtered by `make` to generate the list of object files required for the ICON executable. This is the main purpose of the listing script. Additionally, the tool can run [](ref_buildrun_building_consistencychecks) described in the following subsection.

(ref_buildrun_building_consistencychecks)=
### Code consistency checks

Each failed code consistency check run by the dependency listing script {{ '[`deplist.py`]({}/utils/mkhelper/deplist.py)'.format(base_url) }} is reported to the standard error stream. The identified problems are expressed in terms of files and makefile dependencies and, therefore, require additional explanation provided in this subsection. Normally, the codebase is kept consistent and users do not see the messages described below until they introduce a modification to the source code that breaks the consistency. Currently, the dependency listing script {{ '[`deplist.py`]({}/utils/mkhelper/deplist.py)'.format(base_url) }} checks the source dependency graph for the following problems:

1. **Two or more Fortran source files declare modules with the same name.**

    This type of inconsistency is reported as follows:
    ```
    deplist.py: WARNING: target 'mod/some_module.mod.proxy' has more than one immediate prerequisite matching pattern '*.o':
    	some/dir/some_file.o
    	some/other/dir/some_other_file.o
    ```
    This means that the Fortran module `some_module` is declared twice. The first declaration is found in the file `some/dir/some_file.f90` and the second declaration is found in `some/other/dir/some_other_file.f90`.

2. **Two or more Fortran modules circularly depend on each other.**

    This type of inconsistency is reported as follows:
    ```
    deplist.py: WARNING: the dependency graph has a cycle:
    	src/drivers/icon.o
    	...
    	mod/some_module.mod.proxy
    	some/dir/some_file.o
    	mod/some_module_1.mod.proxy <- start of cycle
    	some/other/dir/some_file_1.o
    	mod/some_module_2.mod.proxy
    	some/other/dir2/some_file_2.o
    	mod/some_module_1.mod.proxy <- end of cycle
    ```
    This reads as that the module `some_module_1` (declared in `some/dir/some_file_1.f90`) uses module `some_module_2` (declared in `some/other/dir/some_file_2.f90`), which in turn uses `some_module_1`. Usually, this means that the compilation of `some/dir/some_file_1.f90` will fail.

3. **A Fortran module is used but not declared.**

    This problem is reported by the dependency listing script with the following message:
    ```
    deplist.py: WARNING: target 'mod/missing_module.mod.proxy' does not have an immediate prerequisite matching any of the patterns: '*.o'
    ```
    This means that the module `missing_module` is used in one of the source files but there is no Fortran source file in the ICON codebase that declares it. It might be the case, however, that the module is not missing but just not part of the ICON codebase, e.g. `mpi`, `sct`, `yaxt`, etc. Such modules are external to ICON and need to be explicitly specified as such in the file `depgen.f90.config` residing in the current build directory (the file is generated at configuration time based on a template file residing in the source directory. Therefore, to make the modifications persistent, you need to introduce them in the file `depgen.f90.config.in`.

4. **Two or more source files have the same basename.**

    The problem is reported as follows:
    ```
    deplist.py: WARNING: the dependency graph contains more than one target with basename 'some_file.o':
    	some/dir/some_file.o
    	some/other/dir/some_file.o
    ```
    This message reports about two (or more) source (not necessarily Fortran) files `some/dir/some_file.f90` and `some/other/dir/some_file.f90` that compile into objects with the same [basename](https://docs.python.org/3/library/os.path.html#os.path.basename). Although handled by the building system in most cases, having several source files with the same basename in a project is considered bad practice, potentially has negative side effects, and, therefore, is deprecated.

(ref_buildrun_building_cascadeprevention)=
### Compilation cascade prevention

It is important, especially for the development process, that the modifications of the source code done after the initial compilation trigger as few recompilations as possible. One of the basic features of `make` is to keep track of the file modification timestamps. Based on the information from the makefiles generated by the dependency generator {{ '[`depgen.py`]({}//utils/mkhelper/depgen.py)'.format(base_url) }}, the tool triggers recompilation of a source file only if the file itself or a header file it includes, or a Fortran module file it uses has been modified since the last execution. Unfortunately, most of the Fortran compilers update the module files even if their relevant contents do not change, i.e. the modification timestamp of a module file gets updated even if the declaration of the associated Fortran module in the source file remains the same. This leads to so-called *compilation cascades*.

Partially, this issue is circumvented in the building system of ICON as follows.
1. If a Fortran source file `filename.f90` uses a module `modulename`, the corresponding dependency makefile `filename.f90.d` (created by the dependency generator {{ '[`depgen.py`]({}//utils/mkhelper/depgen.py)'.format(base_url) }}) gets an entry declaring the dependency of the respective object file on a module *proxy* file:
    ```make
    filename.o: mod/modulename.mod.proxy
    ```
2. When the compilation of the file declaring the module `modulename` takes place for the first time, the original module file `mod/modulename.mod` generated by the compiler is backed up under the name `mod/modulename.mod.proxy`.
3. When `make` checks whether the object file `filename.o` needs to be updated (i.e. the source file `filename.f90` needs to be recompiled), it compares the potentially updated original module file `mod/modulename.mod` with the proxy file `mod/modulename.mod.proxy` and triggers the recompilation only if they are *significantly different*. The latter fact is determined in two steps: first, the files are compared for binary identity with the `cmp` command, second, if the first check shows the difference, the contents of the files are compared with the `/utils/mkhelper/fortmodcmp.py` script, which employs compiler-specific heuristics.
4. Each time the proxy file `mod/modulename.mod.proxy` is detected to be significantly different from the original module file `mod/modulename.mod`, it is replaced with a new copy of the latter.

The described mechanism helps to avoid compilation cascades in many cases. However, the structure of the module files generated by most of the compilers is usually not documented, which makes the comparison of the module files difficult. Thus, the redundant recompilations are not guaranteed to be eliminated entirely.

(ref_buildrun_building_bundled)=
## Building of the bundled packages

The building of the [bundled packages](ref_buildrun_configuration_bundled) is based on the makefiles generated by their configure and CMake scripts. The makefiles are put to the corresponding subdirectories of the `./externals` directory residing in the root build directory of ICON. The packages are built before the compilation of any Fortran source file of ICON takes place. This is done to make sure that the interface Fortran modules of the libraries are available in advance.

Once called in the build directory of ICON, `make` [recursively](https://www.gnu.org/software/make/manual/html_node/Recursion.html) runs itself in the build directories of the bundled packages. The list of targets passed to the instances of `make` running in the directories of the bundled packages depends on the list of targets specified by the user when calling `make` in the build directory of ICON. The targets `all`, `mostlyclean` (for Autotools-based packages), `clean`, `distclean`, or `check` (for Autotools-based packages) and `test` (for CMake-based packages) are preserved and passed over. All other targets are filtered out.

:::{admonition} Check bundled libraries
:class: admonition-icontheme
Many of the bundled packages have a collection of tests, which can be triggered by running `make check` command from the root build directory of ICON. The tests can help to identify potential runtime problems at an early stage and make sure that the core functionality of a package works as expected in the given software environment.
:::

(ref_buildrun_building_provenance)=
## Source provenance collection

Source provenance information is collected at the building stage and injected into the ICON executable. This information is saved at runtime in the output files of the model so that the latter can be matched with the exact version of ICON that was used to produce them. The information is collected automatically with the help of the {{ '[`pvcs.py`]({}/utils/pvcs.py)'.format(base_url) }} script. The script generated a source file `version.c` containing the URL of the git repository, the name of the git branch, and the hash of the git commit. The source file is then treated by `make` as part of the ICON codebase.
