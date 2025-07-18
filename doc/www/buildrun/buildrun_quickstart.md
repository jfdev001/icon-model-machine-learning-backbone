(ref_buildrun_quickstart)=
# Quick Start

## Obtaining the Code

ICON is simultaneously developed in several repositories. There is the *primary* repository [icon](https://gitlab.dkrz.de/icon/icon) and several *secondary* ones: [icon-nwp](https://gitlab.dkrz.de/icon/icon-nwp), [icon-mpim](https://gitlab.dkrz.de/icon/icon-mpim), etc. If unsure, use the primary repository. If you don't have access to any of those, choose the *public* [icon-model](https://gitlab.dkrz.de/icon/icon-model) repository.

Clone the ICON repository of choice and its submodules using the following command:

```sh
git clone --recursive https://gitlab.dkrz.de/icon/icon-model.git
```

Users _with a DKRZ GitLab account_ are recommended to clone the repository via the SSH protocol:

```sh
git clone --recursive git@gitlab.dkrz.de:icon/icon-model.git
```

## Configuring and Building

The configuration step is typically executed by running the {{ '[`configure`]({}/configure)'.format(base_url) }} script with command-line arguments. These arguments specify the locations of libraries and tools needed for building. Note that the list of arguments required for successful configuration can be lengthy and complex, so it's recommended to use a platform- or machine-specific *configure wrapper*, which automatically sets the necessary compiler and linker flags along with the recommended configuration options. You can find the configure wrappers in the respective subdirectories of the {{ '[`config`]({}/config)'.format(base_url) }} directory.

For example, to build ICON on [Levante@DKRZ](https://docs.dkrz.de/doc/levante/index.html) with OpenMP enabled using the Intel compiler, run the following command:

```sh
./config/dkrz/levante.intel --enable-openmp
```

Alternatively, you can create a directory and perform an *out-of-source* build:

```sh
mkdir build && cd build
/path/to/icon/config/dkrz/levante.intel --enable-openmp
```

Using an out-of-source build, you can build ICON in multiple configurations using the same copy of the source code.

The building step is done by running `make` command with an optional argument specifying the number of jobs to run simultaneously. Usually, 8 is a good choice. For example:

```sh
make -j8
```

The result of the building &mdash; the executable file of ICON &mdash; is saved to the `bin` subdirectory of the build directory.

:::{admonition} ICON on your system
:class: admonition-icontheme
If you want to build and run ICON on your personal computer, consider using the [**generic configure wrapper**](ref_buildrun_configuration_wrappersgeneric).
For **detailed information** on **[](ref_buildrun_configuration)** and **[](ref_buildrun_building)** please refer to the respective section.
:::

## Running ICON

To run ICON, you need to create a runscript that sets the required environment variables and calls the executable. One way to get started with running ICON is to use the [`mkexp`](https://gitlab.dkrz.de/esmenv/mkexp) utility to generate a runscript and experimental configuration for the "bubble" test.

```sh
../utils/mkexp/mkexp bubble.config
```

The command above generates a runscript called `../experiments/bubble/scripts/bubble.run_start`. To execute the experiment, navigate to the directory containing the script and submit it for the execution:

```sh
cd ../experiments/bubble/scripts
sbatch bubble.run_start
`````

The output will be directed to the `Work directory` identified after the execution of `mkexp`.

:::{admonition} Running ICON
:class: admonition-icontheme
For **detailed information** please refer to section [](ref_buildrun_running).
:::

# FAQ

1. <a name="faq1" href="#faq1">**I run the configure script without any arguments and it fails. What should I do?**</a>

    First, we recommend checking whether there is a suitable [](ref_buildrun_configuration_wrappers) in the {{ '[`config`]({}/config)'.format(base_url) }} directory that you could use instead of running the configure script directly. If that is not the case, you need at least to specify the `LIBS` variable telling the configure script which libraries to link the executables to. The content of the list depends on the configure options you specify (see [Table 1](tab_icon_depgraph)), for example:

    ```sh
    ./configure --disable-mpi --disable-coupling LIBS='-lnetcdff -lnetcdf -llapack -lblas'
    ```

    If the libraries reside in nonstandard directories, you might also need to specify the `FCFLAGS`, `CPPFLAGS`, and `LDFLAGS` variables to tell the script which directories need to be searched for header and library files (see section [](ref_buildrun_configuration)) for more details).

2. <a name="faq2" href="#faq2">**How can I reproduce the configuration of the model used in a Buildbot test?**</a>

    Scripts run by Buildbot for configuration and building of the model reside in the `config/buildbot` directory. You can run them manually on the corresponding machine.

3. <a name="faq3" href="#faq3">**I get an error message from the configure script starting with _"configure: error: unable to find sources of..."_. What does this mean?**</a>

    Most probably, you forgot to initialize and/or update git submodules. You can do that by switching to the *source* root directory of ICON and running the following command:

    ```sh
    git submodule update --init
    ```

4. <a name="faq4" href="#faq4">**I have problems configuring/building ICON. What is the most efficient way to ask for help?**</a>

    Whoever you ask for help will appreciate receiving the log files. You can generate a tarball with the log files by running the following commands from the root build directory of ICON:

    ```sh
    make V=1 2>&1 | tee make.log
    tar --transform 's:^:build-report/:' -czf build-report.tar.gz $(find . -name 'config.log' -o -name 'CMakeCache.txt') make.log
    ```

    The result of the commands above will be file `build-report.tar.gz`, which should be attached to the very first email describing your problem. Please, do not forget to specify the **repository** and the **branch** that you experience the issue with, preferably in the form of a URL (e.g. {{ '[https://gitlab.dkrz.de/icon/icon-model/-/tree/icon-{}-public](https://gitlab.dkrz.de/icon/icon-model/-/tree/icon-{}-public)'.format(release, release) }}).
