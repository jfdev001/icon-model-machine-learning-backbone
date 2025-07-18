#!/bin/bash

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#############################################
#
#  Download 4D monthly aerosol climatology file
#  'aerosol_cams_climatology_49r2_1951-2019_4D.nc' from:
#  https://aux.ecmwf.int/ecpds/home/radiation/aerosol_climatology/
#  and rename it as 'aerosol_cams_climatology_49r2_1951-2019_4D_orig.nc'
#
#  The climatology contains aerosol mixing ratios with dimensions
#  (epoch, month, lev, lat, lon). Since CDO cannot cope with five
#  dimensions and requires the time variable to be named 'time',
#  this script renames 'month' into 'time', selects only the
#  last epoch of the dataset (2015) and drops the epoch dimension
#  before interpolating the fields onto the ICON grid.
#
#  References:
#  Bozzo et al. 2019 https://doi.org/10.5194/gmd-2019-149
#  Flemming et al. 2013 https://doi.org/10.5194/acp-17-1945-2017
#  Remy et al. 2019 https://doi.org/10.5194/gmd-12-4627-2019
#  More information on the climatology can also be found in the
#  file metadata.
#
#  How to run:
#  Copy this script, together with ICON grid file and original
#  climatology file into the same folder. Ensure that
#  CDO (version 2.1.0 or later), NCO (4.8.1 or later)  and
#  python3 (3.10.11 or later, incl. numpy, xarray) are available.
#  Execute the script.
#
#  Note: CDO versions more recent than 2.0.6 may produce HDF-
#  related warnings. These can be ignored.
#
#############################################

BaseName=aerosol_cams_climatology_49r2_1951-2019_4D
GridName=icon_grid_0026_R03B07_G # an example

# name of the original CAMS climatology file
origFile=${BaseName}_orig.nc
timeFile=${BaseName}_time.nc

# rename month dimension as 'time'
ncrename -O -d month,time ${origFile} ${timeFile}

# create netcdf file only containing correct time values
cat > create_time.py << EOF
import xarray as xr
import numpy as np

fileout='./timevar.nc'

time_values = np.array([1,2,3,4,5,6,7,8,9,10,11,12])
dummy_data = np.random.random(len(time_values))

ds = xr.Dataset(
    {
        "dummy_var": (["time"], dummy_data)  # Variable definition
    },
    coords={
        "time": time_values  # Define the time coordinate
    }
)
ds.encoding["unlimited_dims"] = {"time"}
ds.to_netcdf(fileout)
EOF

runpython=`python3 create_time.py`
$runpython

# append newly created time variable
ncks -m -A -C -v time ./timevar.nc ${timeFile}
# edit attributes
ncatted -O -a standard_name,time,o,c,'time' ${timeFile}
ncatted -O -a units,time,o,c,"months since 2001-1-15 24:00:00" ${timeFile}
ncatted -O -a calendar,time,o,c,"proleptic_gregorian" ${timeFile}
ncatted -O -a axis,time,o,c,"T" ${timeFile}


# extract last epoch (index 12)
ncks -d epoch,12,12 ${timeFile} ${BaseName}_epoch12.nc

# remove "epoch" dimension from file, leaving dimensions
# (time,lev,lat,lon)
cat > drop_dimension.py << EOF
import xarray as xr

filein='./aerosol_cams_climatology_49r2_1951-2019_4D_epoch12.nc'
fileout='./aerosol_cams_climatology_49r2_1951-2019_4D.nc'

ds=xr.open_dataset(filein,decode_times=False)

ds_without_epoch = ds.isel(epoch=0)
ds_without_epoch = ds_without_epoch.drop_vars("epoch")

ds_without_epoch.to_netcdf(fileout)
EOF

runpython=`python3 drop_dimension.py`
$runpython

# name of the CAMS climatology file after the above modifications
sourceFile=${BaseName}.nc

# name of the target grid file. File provided by user
TARGETGRID=${GridName}.nc

# name of the output file with the interpulated CAMS climatology on ICON grid
OFILE=icon_cams_clim_${GridName}.nc

# interpolate climatology file onto ICON grid
# NOTE: bicubic (remapbic) interpolation leads to overshooting features around high orography!
# NOTE: conservative (remapcon) interpolation leads to visible "edges" along lat/lon squares of
#       original grid
# remapbil appears to be the best interpolation option using cdo at this time.
# Pick your poison!
GRIDNUM=`cdo sinfov ${TARGETGRID} | grep nvertex=3 | awk '{print $1}'`
cdo -s -P 4 remapbil,"${TARGETGRID}:${GRIDNUM}"  ${sourceFile} t1.nc
cdo mul -gec,0.0 t1.nc t1.nc t2.nc
cdo add -mulc,0.0 -ltc,0.0 t1.nc t2.nc ${OFILE}
rm -rf t1.nc t2.nc

#cleanup
rm ${timeFile} ${BaseName}_epoch12.nc timevar.nc create_time.py drop_dimension.py

echo 'Done'
