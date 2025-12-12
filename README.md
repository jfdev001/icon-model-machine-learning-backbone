# Building ICON with FTorch

The purpose of this repository is to illustrate how one can build the
[ICON](https://gitlab.dkrz.de/icon/icon-model) numerical weather prediction and climate model with support for
[FTorch](https://github.com/Cambridge-ICCS/FTorch) as an external library. The
goal is to provide boilerplate for including machine learning models built using PyTorch into
a Fortran codebase---in this case we focus on ICON. 

**IMPORTANT**: The first part of this readme focuses on building the *public*
version of ICON with FTorch. Be aware that for private versions of ICON (e.g.,
`icon-nwp/uaicon-iap-dev`) you will need to check [Setup FTorch with a
different version of ICON](#setup-ftorch-with-a-different-version-of-icon),
though the steps are mostly the same, so it is advisable to read the public
version setup first.

## Getting the code 

To get started, do the following:

```shell 
# Get fork of public ICON
git clone https://github.com/jfdev001/icon-model-machine-learning-backbone.git

# Cache the path to the icon model source directory
icon_model_src=$(readlink -f icon-model-machine-learning-backbone)

# enter icon src dir
cd icon-model-machine-learning-backbone 

# get external libraries
git submodule update --init 
```

## Setting up Torch

To build FTorch, you will need a working PyTorch or LibTorch installation.
To that end, we install a PyTorch[CPU] version since this is the 
most straightforward. Note, we have not yet tested ICON with PyTorch[GPU].

```shell 
# setup virtual environment
python -m venv .venv 

# Activate the venv 
. .venv/bin/activate 

# install torch dependency
pip install torch --index-url https://download.pytorch.org/whl/cpu

# Save the libtorch directory 
python_version=$(python --version | grep -o -E "[0-9]\.[0-9]*")
export Torch_DIR=$(pwd)/.venv/lib/python${python_version}/site-packages/
```

While the above steps use a `venv`, you could easily use a `conda` environment
instead. If you, for example, already have a working PyTorch installation, you
should be able to simply set `Torch_DIR` to the `site-packages` corresponding
to your `conda` environment with the PyTorch installation, e.g.,:

```shell
# here, my conda environment is named `ftorch`, so make sure and change 
# that to your env name... pay attention also to the python version since you 
# may have a different version (i.e., not 3.13)
export Torch_DIR=$HOME/miniconda3/envs/ftorch/lib/python3.13/site-packages
```

If in an HPC environment, you may refer to the Torch installation associated 
with `module load`, e.g., 

```shell
module load pytorch/2.5.1
export Torch_DIR=$(dirname $(python -c "import torch; print(torch.__file__)"))
```

You could also just use the C++ library LibTorch if you want a more lightweight
approach. The CPU only version is available
[here](https://docs.pytorch.org/cppdocs/installing.html). Remember, as in the
previous steps, you will need to `export Torch_DIR` with something like 

```shell 
export Torch_DIR=/path/to/libtorch
```

## On a local system using Docker

Once you have installed the PyTorch dependencies, you can configure and compile
ICON in a generic docker environment locally:

```shell 
# Pull docker image 
docker image pull iconmodel/icon-dev 

# Run the container in interactive mode, mounting appropriate volumes and
# setting environment variables 
sudo docker run \
    -v ${icon_model_src}:/home/icon/icon-src \
    -e Torch_DIR=/home/icon/libtorch \
    -v ${Torch_DIR}:/home/icon/libtorch \
    -it iconmodel/icon-dev

# in docker image... make build dir 
mkdir icon-build

# Configure ICON with FTorch 
cd icon-build
$HOME/icon-src/config/generic/gcc --enable-ftorch

# Compile ICON 
make -j8
```

The Docker environment does not easily facilitate running tests because one
needs access to grid files and other data. This section was written only to
work out any kinks involved in configuring and compiling ICON + FTorch, but it
has not undergone integration testing, though configuration and compilation
succeed. For further ICON usage documentation, see
[icon-model.org](https://www.icon-model.org/).

## On Levante 

This section is relevant for users of the [Levante HPC system at
DKRZ](https://docs.dkrz.de/doc/levante/index.html).

The same setup described in [Setting up Torch](#setting-up-torch) is necessary
before performing the below steps, which are nearly identical to the steps in
[On a local system using Docker](#on-a-local-system-using-docker).

```shell 
# Switch to a directory with sufficient space (e.g., datawave)
cd /work/bm1233/${USER}

# Make a directory for ICON source codes 
mkdir icon-srcs

# Make a build directory for out of source ICON builds 
mkdir -p icon-builds/icon-model-machine-learning-backbone

# Cache the path to the ICON build directory 
icon_build_dir=$(readlink -f icon-builds/icon-model-machine-learning-backbone)

# Get fork of public ICON
cd icon-srcs 
git clone https://github.com/jfdev001/icon-model-machine-learning-backbone.git

# Cache the path to the icon model source directory
icon_model_src=$(readlink -f icon-model-machine-learning-backbone)

# Configure ICON with FTorch
cd ${icon_build_dir}
${icon_model_src}/config/dkrz/levante.gcc --enable-ftorch

# Compile ICON 
make -j8

# Make the run scripts 
./make_runscripts atm_tracer_Hadley
```

On Levante, sufficient data exists to run a test case to verify that FTorch 
is linked properly to ICON:

```shell
# Run the test job (takes < 1 minute and uses only 1 compute node)
email=YOUR_EMAIL_HERE # modify this!
sbatch --mail-user=${email} --mail-type=ALL run/exp.atm_tracer_Hadley.run
```

## Configuring ICON with your own FTorch installation

The previous section relies on the ICON build system to compile FTorch. If you
prefer, you can build FTorch separately, but then you must manually provide
`FCFLAGS`,`LIBS` as well as `LIBDIR` during configuration. Assuming on DKRZ
Levante, you may instead do the following to configure ICON with your own
FTorch installation:

```shell
# change to build directory 
cd ${icon_build_dir}

# configure using an external FTorch installation...
# assumes you have already installed FTorch
# see https://github.com/Cambridge-ICCS/FTorch?tab=readme-ov-file#installation
ftorch_dir=</path/to/FTorch/installed>
${icon_model_src}/config/dkrz/levante.gcc \
    --enable-ftorch --with-external-ftorch \
    ftorch_FCFLAGS="-I${ftorch_dir}/include -I${ftorch_dir}/include/ftorch" \
    ftorch_LIBDIR="${ftorch_dir}/lib64" \
    ftorch_LIBS="-lftorch"

# compile as usual
make -j8

# Make the run scripts 
./make_runscripts atm_tracer_Hadley

# Run the test job (takes < 1 minute and uses only 1 compute node)
email=YOUR_EMAIL_HERE # modify this!
sbatch --mail-user=${email} --mail-type=ALL run/exp.atm_tracer_Hadley.run
```

## Internally coupling ICON with FTorch 

For coupling machine learning models with ICON, the most straightforward
approach is the internal coupling approach: that is, you will create and modify
source files in the ICON codebase itself. You could use a dedicated coupler
such as [ComIn](https://gmd.copernicus.org/articles/18/1001/2025/), note
however its implementation of MPI blocking during coupling, which could cause
significant performance penalties.

As an example of the internal coupling approach, you might create a file
containing your key machine learning model logic (i.e., this might be model
loading, any preprocessing etc.). You might define a new routine called `Unet`
in a file called
`src/upper_atmosphere/mo_upatmo_iondrag_machine_learning_model.f90`, which
could be subsequently imported in another file as shown below:

```fortran 
! @file src/upper_atmosphere/mo_upatmo_phy_iondrag.f90
MODULE mo_upatmo_phy_iondrag
    ...
    USE mo_upatmo_iondrag_machine_learning_model, ONLY: Unet
    ...
END MODULE mo_upatmo_phy_iondrag
```

**IMPORTANT**: Anytime you modify source files, you will have to recompile ICON.
Switch back to the build directory and call make like so:

```shell
cd ${icon_build_dir}
make -j8
```

If you add new files to the ICON source code or change configuration files like
`configure.ac`, you will have to reconfigure **and** recompile like so:

```shell
cd ${icon_model_src}
autoreconf -fvi
cd ${icon_build_dir}
${icon_model_src}/config/dkrz/levante.gcc --enable-ftorch && make -j8
```

Currently, to ensure the the test case in section [On Levante](#on-levante)
works, we have added a dummy call to an FTorch function in the
`src/io/shared/mo_output_event_handler.f90`. That call serves no other purpose.

For more FTorch examples, see
[FTorch/examples](https://github.com/Cambridge-ICCS/FTorch/tree/ef44cf6d70edef38003dec41ee7a1b496922a1a5/examples).

## Setup FTorch with a different version of ICON

In the current repository, we have already made the necessary additions to
facilitate building ICON with FTorch. *However*, if you are using a different
version of ICON (e.g., `icon-nwp`), you will need to either (a) update the
`externals/` and `configure.ac` in order to support FTorch for which we provide
instructions in [The automated approach](#the-automated-approach) or (b) make
the minimum set of modifications to your existing ICON configuration following
the instructions in [The manual approach](#the-manual-approach).

The example in this section assumes access to DKRZ's Levante; however, if
you simply substitute this example with your local or HPC system's
ICON, the instructions still apply.

To illustrate the necessary changes to integrate FTorch into a different
version of ICON, we select the `icon-nwp/uaicon-iap-dev` branch 
and the commit `881b54a9d12a17f8e5a5a6930f86d3c3ce5aa19c` corresponding 
to 2025-07-02. At the time of this writing (2025-12-05), that branch
still exists, but may not in the future due to refactoring efforts.

To get access to this version of ICON, do the following:

```shell
# Switch to a directory with sufficient space (e.g., datawave)
cd /work/bm1233/${USER}

# Make a build directory for out of source ICON builds 
# NOTE: If this directory already exists, it would be wise to delete it 
# since if you attempt to re-compile sources of potentially mismatched UAICON
# versions/config files, you will get errors
mkdir -p icon-builds/uaicon-with-ftorch

# Cache the path to the ICON build directory 
icon_build_dir=$(readlink -f icon-builds/uaicon-with-ftorch)

# Make a directory for ICON source codes 
mkdir icon-srcs

# Get the uaicon codebase -- this will fail if you do not
# have developer/report access to the DKRZ gitlab 
cd icon-srcs 
git clone git@gitlab.dkrz.de:icon/icon-nwp.git uaicon-with-ftorch

# Cache the path to the uaicon source code 
icon_model_src=$(readlink -f uaicon-with-ftorch)

# Change to the source dir and use a specific uaicon version 
cd uaicon-with-ftorch
git checkout -b uaicon-iap-dev origin/uaicon-iap-dev
git checkout -b uaicon-iap-dev-2025-07-02 881b54a9d12a17f8e5a5a6930f86d3c3ce5aa19c
```

At this point, if you do not have a version of PyTorch/LibTorch installed, you 
can follow the instructions in section [Setting up Torch](#setting-up-torch).

If you end up taking the `.venv` approach to setting up Torch and setup the 
`.venv` in your ICON project directory, make sure you add this to your
`.gitignore` with 

```shell
echo ".venv" >> .gitignore
```

### The automated approach

This approach automatically modifies `configure.ac` so that the ICON build 
system will compile FTorch for you.

To integrate the necessary changes into an arbitrary ICON project, add the
`setup_ftorch` script in the current repository to the other repository of
interest.

```shell 
# You must define your icon source and build directories... the below paths
# are dummy paths, please modify them
icon_model_src=/path/to/icon-srcs/my-icon
icon_build_dir=/path/to/icon-builds/my-icon

# Get the setup script and execute it 
wget https://raw.githubusercontent.com/jfdev001/icon-model-machine-learning-backbone/refs/heads/release-2025.04-public/setup_ftorch
./setup_ftorch
```

Then, as before, you build your version of ICON in a separate build directory, 
making sure to enable the FTorch components:

```shell
cd ${icon_build_dir}
${icon_model_src}/config/dkrz/levante.gcc --enable-ftorch && make -j8
```

Note that `setup_ftorch` will also add a dummy hello ftorch program to 
`src/io/shared/mo_output_event_handler.f90` in order to verify that FTorch 
works correctly. You can delete that dummy logic if you wish. It only calls 
a dummy ftorch program at the very beginning of the model start so it does not
add anything computationally expensive.

If you happen to be on Levante, you can call the below example run script 
and check the LOG file for the substring `hello ftorch` to verify your 
FTorch installation works:

```shell
cd ${icon_build_dir}
./make_runscripts atm_tracer_Hadley

# Run the test job (takes < 1 minute and uses only 1 compute node)
email=YOUR_EMAIL_HERE # modify this!
sbatch --mail-user=${email} --mail-type=ALL run/exp.atm_tracer_Hadley.run
```
### The manual approach

If do not wish to use the automated approach, we describe here the minimum
modifications you must make in order to compile ICON with FTorch support. This
approach is less generic than the automated approach. You can just modify one
of the `config/<target>/<config file>` (henceforth called ICON config files)
files directly and add the appropriate `FCFLAGS`, `LDFLAGS`, and `LIBS`
corresponding to FTorch.

For the manual approach, you will need to compile FTorch using the *same*
compilers you are using to compile ICON. So, you cannot compile FTorch with
Intel compilers and then compile ICON with GNU compilers. This means that you
need to pay attention to the compilers specified in the ICON config files. For
example, if you wish to build FTorch using the same compilers as
`config/dkrz/levante.gcc`, inspecting that file shows that the `MPI_ROOT`
variable can be used to determine what the corresponding C, C++, and Fortran
compilers are that you must use to build FTorch (search `MPI_ROOT`, `CC`,
`CXX`, and `FC`):

```shell
# @file levante.gcc
MPI_ROOT='/sw/spack-levante/openmpi-4.1.2-mnmady'
CC="${MPI_ROOT}/bin/mpicc"
CXX=$("${MPI_ROOT}/bin/mpicxx" -show | sed 's: .*$::')
FC="${MPI_ROOT}/bin/mpif90"
```

You must then use these compilers to build FTorch accordingly (make sure you
have a Torch installation available, see [Setting up Torch](#setting-up-torch)):

```shell
# You must define your icon source and build directories... the below paths
# are dummy paths, please modify them
icon_model_src=/path/to/icon-srcs/my-icon
icon_build_dir=/path/to/icon-builds/my-icon

# Change to your icon model project directory
cd ${icon_model_src}

# Get FTorch source code (placed in externals directory for convenience)
mkdir externals
cd externals 
git clone --depth 1 https://github.com/Cambridge-ICCS/FTorch.git

# Define compilers for FTorch based on levante.gcc 
MPI_ROOT="/sw/spack-levante/openmpi-4.1.2-mnmady"
CC="${MPI_ROOT}/bin/mpicc"
CXX="${MPI_ROOT}/bin/mpicxx"
FC="${MPI_ROOT}/bin/mpif90"

# Build FTorch -- assumes PyTorch installed in a venv in the ICON project dir
cd FTorch
mkdir build
cd build

python_version=$(python -c "import sys; print('.'.join(sys.version.split('.')[:2]))")
export Torch_DIR=${icon_model_src}/.venv/lib/python${python_version}/site-packages/torch

cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(pwd)/installed -DCMAKE_Fortran_COMPILER=${FC} -DCMAKE_C_COMPILER=${CC} -DCMAKE_CXX_COMPILER=${CXX}
cmake --build . --target install
```

We suggest that you make a modification to one of the files in the `src/`
directory such that the `USE ftorch` statement is included. This way, when 
you compile ICON in the next steps, you know that you can use FTorch in your 
codebase.

You must also modify the base configuration file that you wish to use to add
the appropriate flags to support FTorch. For example, you might update
`config/dkrz/levante.gcc` as follows (search `BUILD_ENV`, `FCFLAGS`, `LDFLAGS`,
and `LIBS`):

```shell
# @file levante.gcc
# FTorch configuration
FTORCH_ROOT="${icon_dir}/externals/FTorch/build/installed"
FTORCH_LIBDIR="${FTORCH_ROOT}/lib64"
FTORCH_LIBS="-lftorch"
FTORCH_LDFLAGS="-L${FTORCH_LIBDIR}"
FTORCH_INCLUDEDIR="${FTORCH_ROOT}/include"
FTORCH_MODULEDIR="${FTORCH_INCLUDEDIR}/ftorch"
FTORCH_FCFLAGS="-I${FTORCH_INCLUDEDIR} -I${FTORCH_MODULEDIR}"

BUILD_ENV="export LD_LIBRARY_PATH=\"${FTORCH_LIBDIR}:${FYAML_ROOT}/lib:${HDF5_ROOT}/lib:${NETCDF_ROOT}/lib:${NETCDFF_ROOT}/lib:${ECCODES_ROOT}/lib64:\${LD_LIBRARY_PATH}\"; export PATH=\"${HDF5_ROOT}/bin:\${PATH}\"; ${BLAS_LAPACK_BUILD_ENV}"

FCFLAGS="${FTORCH_FCFLAGS} -I${HDF5_ROOT}/include -I${NETCDFF_ROOT}/include -I${ECCODES_ROOT}/include -fmodule-private -fimplicit-none -fmax-identifier-length=63 -Wall -Wcharacter-truncation -Wconversion -Wunderflow -Wunused-parameter -Wno-surprising -fall-intrinsics -g -march=native -mpc64"

LDFLAGS="${FTORCH_LDFLAGS} -L${HDF5_ROOT}/lib -L${NETCDF_ROOT}/lib -L${NETCDFF_ROOT}/lib ${BLAS_LAPACK_LDFLAGS} -L${ECCODES_ROOT}/lib64 -L${FYAML_ROOT}/lib"

LIBS="-Wl,--disable-new-dtags -Wl,--as-needed ${FTORCH_LIBS} ${XML2_LIBS} ${FYAML_LIBS} ${ECCODES_LIBS} ${BLAS_LAPACK_LIBS} ${NETCDFF_LIBS} ${NETCDF_LIBS} ${HDF5_LIBS} ${ZLIB_LIBS} ${STDCPP_LIBS}"
```

Then, configure and compile your version of ICON in a separate build directory:

```shell
cd ${icon_build_dir}
${icon_model_src}/config/dkrz/levante.gcc && make -j8
```

Note that we do not pass the `--enable-ftorch` flag since we assume the
`configure.ac` is the default `configure.ac` which does not have the
appropriate modifications to support this flag.

One last addition: if you are using an FTorch version on or after commit
`ddfd35a5b0b0874f05aa7345f06ae25a8c35a4a8`, FTorch has a pkg config file. This
means you do not have to define the `LDFLAGS`, `LIBS`, and `FCFLAGS` manually,
but rather can rely on querying the flags that the authors of the FTorch
library have specified are necessary for compilation. Therefore,
`config/dkrz/levante.gcc` might instead look like the following:

```shell
# FTorch configuration
FTORCH_PKG_CONFIG="${icon_dir}/externals/FTorch/build/installed/lib64/pkgconfig"
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:${FTORCH_PKG_CONFIG}"
FTORCH_LIBS="$(pkg-config --libs-only-l ftorch)"
FTORCH_LDFLAGS="$(pkg-config --libs-only-L ftorch)"
# regex that returns the first -L${FTORCH_LIBDIR} by matching characters
# until 1 or more spaces is encountered... e.g., for
# -L/path/lib -Wl,-rpath,/path/lib, the regex returns /path/lib
FTORCH_LIBDIR=$(echo ${FTORCH_LDFLAGS} | grep -oP "(?<=-L)[^\s]+" | head -n 1)
FTORCH_FCFLAGS="$(pkg-config --cflags ftorch)"
```

As before, you would still need to modify the `BUILD_ENV`, `FCFLAGS`,
`LDFLAGS`, and `LIBS` variables in the ICON config file.

# Brief Description of ICON Build System

This section is only relevant if you'd like to know some of the internal 
details of the ICON build system. It is not strictly necessary to be aware 
of these in order to use ICON + FTorch.

ICON is built using GNU autotools (i.e., autoconf, automake, and libtool). The
main files one should inspect to see the dependencies and build system are
therefore `configure.ac` and `icon.mk.in`. ICON also provides example
configurations in the `config` directory for different platforms (e.g., DKRZ
Levante). The ICON build system, unlike like e.g., CMake, has custom helper
scripts to resolve dependencies that source/header files have with each other.
From a developer perspective, this is not that important since you don't have
to modify such scripts/helper functions, but it is something to be aware of.
Note, a subset of these helper functions allows ICON to build an external whose
build system is CMake (e.g., FTorch, ComIn, etc.).

With respect to ICON's "externals" (aka dependencies) you can differentiate
between two categories: (1) externals that are completely independent of ICON
and can therefore be used and integrated into arbitrary applications, and (2)
externals that are internally coupled to ICON due to their dependency on ICON
core data structures/functions (e.g., ICON-ART depends on ICON's grid data
structure). Category (2) externals **cannot** be used independently of the ICON
source code. 

As of this writing (2025-12-05), the following externals belong to category
(2):

* ART 
* JSBACH 
* EMVORADO (not available in ICON public release, but available in ICON-NWP)
* DACE  (not available in ICON public release, but available in ICON-NWP)

The remaining externals, belong to category (1). Note, category (2) externals 
were identified by inspecting `icon.mk.in`.

From a developer perspective you only need to be concerned with `configure.ac`
if you are adding configuration for a category (1) external. It is not relevant
for FTorch, but should you wish to add an external that is tightly coupled to
to ICON (e.g., something similar to ART), then you would have to modify both
`icon.mk.in` *and* `configure.ac`. 

Should you make any changes to any configuration scripts, make sure to call

```shell
# assuming in ICON repo
autoreconf --force --install --verbose 
```

# Contributing

This repository is intended to be used as a template for facilitating including
ML models into ICON via FTorch. For up to date ICON model developments, please
open issues and make pull requests in the [DKRZ ICON model
repository](https://gitlab.dkrz.de/icon/icon-model).

Please feel free, however, to open issues in the present repository related to
documentation clarity or problems building ICON with FTorch.
