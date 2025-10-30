# Building ICON with FTorch

The purpose of this repository is to illustrate how one can build the
[ICON](https://gitlab.dkrz.de/icon/icon-model) model with support for
[FTorch](https://github.com/Cambridge-ICCS/FTorch) as an external library. The
goal is to faciliate including machine learning models built using PyTorch into
a Fortran codebase---in this case we focus on ICON. 

**IMPORTANT**: The first part of this readme focuses on building the *public*
version of ICON with FTorch. Be aware that for private versions of ICON (e.g.,
`icon-nwp/uaicon-iap-dev`) you will need to check [Setup with a different
version of ICON](#Setup-with-a-different-version-of-ICON), though the steps 
are mostly the same, so it is advisable to read the public version setup first.

## General setup with a public version of ICON

To get started with a public version of ICON, do the following:

```shell 
# get fork of public ICON
git clone --branch jfrazier/integrate-ftorch --depth 1 \
    https://gitlab.dkrz.de/b383137/icon-model-with-ftorch.git 

# Cache the path to the icon model source directory
icon_model_src=$(readlink -f icon-model-with-ftorch)

# enter icon src dir
cd icon-model-with-ftorch 

# get external libraries
git submodule update --init 
```

In this repository, I have already made the necessary additions to the
`configure.ac` to facilitate building ICON with FTorch. Since this repository
is based on the public release of ICON, you may need to make adjustments
in `configure.ac` pertaining to FTorch into the `configure.ac` of **your** ICON
repository (i.e., you may be using `icon-nwp/uaicon-iap-dev` which might rely
on a different version of `configure.ac`). 

Now, to build FTorch, you will need a working PyTorch or LibTorch installation.
To that end, we install a PyTorch[CPU] version since this is the 
most straightforward. Note, I have not yet tested ICON with PyTorch[GPU].

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

Note that while the above steps use a `venv`, you could easily use 
a `conda` environment instead. If you, for example, already have 
a working PyTorch installation, you should be able to simply 
set `Torch_DIR` to the `site-packages` corresponding to your 
`conda` environment with the PyTorch installation, e.g.,:

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

## On a local system using Docker

Once you've got the PyTorch dependencies, you can configure and compile ICON in
a generic docker environment locally:

```shell 
# Pull docker image 
docker image pull iconmodel/icon-dev 

# Run the container in interactive mode, mounting appropriate volumes and
# setting environment variables 
sudo docker run -it \
    -v ${icon_model_src}:/home/icon/icon-src \
    -e Torch_DIR=/home/icon/libtorch \
    -v ${Torch_DIR}:/home/icon/libtorch

# in docker image... make build dir 
mkdir icon-build

# Configure ICON with FTorch 
cd icon-build
$HOME/icon-src/config/generic/gcc --without-external-ftorch

# Compile ICON 
make -j8
```

The Docker environment does not easily facilitate running tests because one
needs access to grid files and other data that are not easily accessible. This
section was written only to work out any kinks involved in configuring and
compiling ICON + FTorch, but it is technically untested in terms of
functionality (even though config and compilation succeeds).

## On Levante 

The same setup described in [General setup with a public version of ICON](#General-setup-with-a-public-version-of-ICON) is necessary before performing the below steps,
which are nearly identical to the steps in [On a local system using Docker](#On-a-local-system-using-Docker).

```shell 
# Switch to a directory with sufficient space (e.g., datawave)
cd /work/bm1233/${USER}

# Make a build directory for out of source ICON builds 
mkdir -p icon-builds/icon-model-with-ftorch

# Cache the path to the ICON build directory 
icon_build_dir=$(readlink -f icon-builds/icon-model-with-ftorch)

# Configure ICON with FTorch  (TODO: try with intel as well)
cd icon-build/icon-model-with-ftorch
${icon_model_src}/config/dkrz/levante.gcc --without-external-ftorch

# Compile ICON 
make -j8

# Make the run scripts 
./make_runscripts --all
```

On Levante, sufficient data exists to run a test case to verify that FTorch 
is linked properly to ICON:

```shell
# NOTE: if you want to get notifications on your jobs completion status,
# make sure to add the following fields to the below run script:
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH --mail-type=ALL

# Run the test job (takes < 1 minute and uses only 1 compute node)
${icon_build_dir}/run/exp.atm_tracer_Hadley.run
```

**IMPORTANT**: Any new runscripts you write need to modify the
`LD_LIBRARY_PATH` in order to ensure that ICON is able to find FTorch during
runtime. The runscripts generated by `make_runscripts --all` already have 
the appropriate modification. Make sure that 

```shell 
export LD_LIBRARY_PATH=${FTorch_DIR}:${LD_LIBRARY_PATH}
```

is somewhere in your runscript. `FTorch_DIR` is the path to the directory 
containing the file `libftorch.so`. Unless you have done something strange 
during installation, this should be 

```shell
FTorch_DIR=${icon_build_dir}/externals/FTorch/build
```

More details on internal coupling of ICON with FTorch in the next section.


## Internally coupling ICON with FTorch 

For coupling machine learning models with ICON, you will have to take the
internal coupling approach: that is, you will create and modify source files in
the ICON codebase itself. Using a dedicated coupler such as
[ComIn](https://gmd.copernicus.org/articles/18/1001/2025/) is not ideal due to
its implementation of MPI blocking during coupling, which could cause
significant performance penalties. More practically, ComIn may not work with a
non-standard version of ICON (e.g., `icon-nwp/uaicon-iap-dev`).

As an example of the internal coupling approach, you might create a file
containing your key machine learning model logic (i.e., this might be
model loading, any preprocessing etc.). A simple hello world file is provided 
in `src/upper_atmosphere/mo_upatmo_iondrag_machine_learning_model.f90`.
You might define a new routine in that file called `Unet`, which could 
be subsequently imported in another file as shown below:

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

If you add new files to the ICON source code or other change configuration
files like `configure.ac` or `run/create_target_header`, etc., you will have to
reconfigure **and** recompile like so:

```shell
cd ${icon_build_dir}
${icon_model_src}/config/dkrz/levante.gcc --without-external-ftorch && make -j8
```

Currently, to ensure the the test case in section [On Levante](#On-Levante)
works, I have added a dummy call to an FTorch function in the 
`src/io/shared/mo_output_event_handler.f90`. That call serves no other purpose.

For more FTorch examples, see
[FTorch/examples](https://github.com/Cambridge-ICCS/FTorch/tree/ef44cf6d70edef38003dec41ee7a1b496922a1a5/examples).

## Setup with a different version of ICON

TODO

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
build system is CMake (e.g., like FTorch, ComIn, etc.).

With respect to ICON's "externals" (aka dependencies... you could think of
these as "packages" if you are coming from a Python background), you can
differentiate between two categories: (1) externals that are completely
independent of ICON and can therefore be used and integrated into arbitrary
applications, and (2) externals that are internally coupled to ICON due to
their dependency on ICON core data structures/functions (e.g., ICON-ART depends
on ICON's grid data structure). Category (2) externals **cannot** be used
independently of the ICON source code. 

As of this writing (2025-10-28), the following externals belong to category
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
to ICON (e.g., something similar to ART), then you would have to modify 
both `icon.mk.in` and `configure.ac`.
Should you make any changes to any configuration scripts, make sure to call

```shell
# assuming in ICON repo
autoreconf --force --install --verbose 
```

# References 
* [iap/misc-issues#17: On which configuration files one needs to
change](https://igit.iap-kborn.de/icon-iap/misc-tasks/-/issues/17)

