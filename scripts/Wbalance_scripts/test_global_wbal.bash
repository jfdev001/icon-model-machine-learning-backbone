#!/bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#-----------------------------------------------------------------------------
# Script to check global water conservation of ICON-XPP
#   and to narrow down the problem in case of a water imbalance.
#
#                                                   Veronika Gayler, Feb. 2023
#-----------------------------------------------------------------------------
set -e
CDO="cdo -s"

#------------------------------------------------------------------------------
# Settings section
# ----------------
# Time periode to be examined (needs to match complete files in case of multi year output chunks)
yr1=1800      # First year
yr2=1829      # Last year
interval=10   # output file interval in years
exp=vga0424
dt_atm=900    # atmosphere time step in seconds
dt_oce=1800   # ocean time step in seconds
outdata=/work/mj0060/m220053/icon-esm/land-esm-dev/experiments/$exp/outdata
restdir=/work/mj0060/m220053/icon-esm/land-esm-dev/experiments/$exp/restart
fractions=/pool/data/ICON/grids/public/mpim/0043-0035/land/r0003/bc_land_frac_11pfts_1850.nc
grid_A=/pool/data/ICON/grids/public/mpim/0043/icon_grid_0043_R02B04_G.nc
grid_O=/pool/data/ICON/grids/public/mpim/0035/icon_grid_0035_R02B06_O.nc
#------------------------------------------------------------------------------

workdir=${outdata}/wbal-test_${yr1}-${yr2}
[[ -d ${workdir} ]] || mkdir ${workdir}
cd ${workdir}

# Get fractions
$CDO selvar,glac ${fractions} glac.nc
$CDO selvar,notsea ${fractions} notsea.nc
$CDO -addc,1 -mulc,-1 notsea.nc ocean_fract.nc
$CDO gtc,0 notsea.nc notsea_gt0.nc  # 1 in cells with land fraction
$CDO selvar,cell_area ${grid_A} area_A.nc
$CDO fldsum area_A.nc global_area.nc

$CDO selvar,cell_area ${grid_O} area_O.nc
$CDO selvar,cell_sea_land_mask ${grid_O} cell_sea_land_mask_O.nc
$CDO -gtc,0 -mulc,-1 cell_sea_land_mask_O.nc ocean_fract_O.nc
$CDO fldsum -mul area_O.nc ocean_fract_O.nc global_ocean_area.nc
$CDO div global_ocean_area.nc global_area.nc global_ocean_fraction.nc

# Further preparations
yearlist="" ; yr=$yr1
while [[ $yr -le $yr2 ]]; do
  yearlist="$yearlist $yr"
  (( yr = yr + $interval ))
done
(( nyrs = yr2 - yr1 + 1 ))
(( nsec = $nyrs * 31557600 ))  # 365.25 x 86400
ocean_area=$($CDO output global_ocean_area.nc) # [m2]


# Function to calculate sea level rise [m3 s-1] -> [mm over the years]
function get_sea_level_rise {

wbal_change_in_m3_per_sec=$1
wbal_change=$($CDO output $wbal_change_in_m3_per_sec)
sea_level_rise=$($CDO output -mulc,$nsec -mulc,1000 -div $wbal_change_in_m3_per_sec global_ocean_area.nc)

}
echo ""
echo "------------------------------------------------------------------------"
echo "  Water balance check for experiment $exp"
echo "    Taking into account simulation years $yr1 to $yr2"
echo "------------------------------------------------------------------------"

#
# 1. Diagnostics from the land monitoring
#
filelist_lmon=""
for yr in $yearlist; do
  filelist_lmon="$filelist_lmon ${outdata}/${exp}_lnd_mon_${yr}0101.nc"
done

# a) hd_water_error_gsum_box [m3/time step] -> [m3/s]
$CDO -s -timavg -divc,${dt_atm} -mergetime -apply,-selvar,hd_water_error_gsum_box [ ${filelist_lmon} ] wbal_change.nc
get_sea_level_rise wbal_change.nc
wbal_change=$($CDO output wbal_change.nc)
echo ""
echo "Diagnostics from land monitoring"
echo "  hd_water_error        [m3 s-1]: ${wbal_change} -> Sea level rise in $nyrs years [mm]: $sea_level_rise"

# b) hydro_weq_balance_err_gsum_box [m3/time step] -> [m3/s]
$CDO -s -timavg -divc,${dt_atm} -mergetime -apply,-selvar,hydro_weq_balance_err_gsum_box  [ ${filelist_lmon} ] wbal_change.nc
get_sea_level_rise wbal_change.nc
echo "  hydro_weq_balance_err [m3 s-1]: ${wbal_change} -> Sea level rise in $nyrs years [mm]: $sea_level_rise"

#
# 2. Diagnostics from the ocean monitoring
#
filelist_omon=""
for yr in $yearlist; do
  filelist_omon="$filelist_omon ${outdata}/${exp}_oce_mon_${yr}0101.nc"
done

echo ""
echo "Diagnostics from ocean monitoring"

# a) ssh_global: Sea surface hight [m] -> [mm]
$CDO -setunit,"mm" -mulc,1000 -yearavg -mergetime -apply,-selvar,ssh_global [ ${filelist_omon} ] ssh_global.nc
$CDO -sub -selyear,${yr2} ssh_global.nc -selyear,${yr1} ssh_global.nc ssh_global_rise.nc
ssh_global_rise=$($CDO output ssh_global_rise.nc)
(( nyrm1 = nyrs - 1 ))
echo "  ssh_global                  -> Sea level rise from $yr1 to $yr2, i.e. $nyrm1 years [mm]: $ssh_global_rise"

# b) Fresh water flux at ocean the surface - above ice [m/s] -> [m3/s]
$CDO -O -mergetime \
     -apply,-expr,'FrshFlux_above_ice=(FrshFlux_Precipitation_global+FrshFlux_Evaporation_global+FrshFlux_Runoff_global)' \
     [ ${filelist_omon} ] FrshFlux_above_ice.nc
$CDO -mul global_ocean_area.nc -timavg FrshFlux_above_ice.nc wbal_change.nc
get_sea_level_rise wbal_change.nc
echo "  Fresh water flux above ice [m3 s-1]: ${wbal_change} -> Sea level rise in $nyrs years [mm]: $sea_level_rise"

# c) Fresh water flux at ocean surface - below ice [m/s] -> [m3/s]
$CDO -O -mergetime \
     -apply,-expr,"FrshFlux_below_ice=(FrshFlux_Runoff_global+FrshFlux_TotalOcean_global+FrshFlux_VolumeIce_global+totalsnowfall_global/$dt_oce)" \
     [ ${filelist_omon} ] FrshFlux_below_ice.nc
$CDO -mul global_ocean_area.nc -timavg FrshFlux_below_ice.nc wbal_change.nc
get_sea_level_rise wbal_change.nc
echo "  Fresh water flux below ice [m3 s-1]: ${wbal_change} -> Sea level rise in $nyrs years [mm]: $sea_level_rise"

#
# Extract variables fo further diagnostics
#
filelist=""
for yr in $yearlist; do
  filelist="${filelist} ${outdata}/${exp}_lnd_basic_ml_${yr}0101.nc"
done
for var in hydro_discharge_box hydro_discharge_ocean_box hydro_evapotrans_box hydro_runoff_box hydro_drainage_box; do
  $CDO -O mergetime -apply,-selvar,$var [ ${filelist} ] ${exp}_lnd_basic_ml.$var.tmp
  # Get rid of 1st timestep (Dec 31 of previous year)
  $CDO selyear,$yr1/$yr2 ${exp}_lnd_basic_ml.$var.tmp ${exp}_lnd_basic_ml.$var.nc
  rm ${exp}_lnd_basic_ml.$var.tmp
done

filelist=""
for yr in $yearlist; do
  filelist="${filelist} ${outdata}/${exp}_atm_2d_ml_${yr}0101.nc"
done
for var in pr evspsbl; do
  $CDO -O mergetime -apply,-selvar,$var [ ${filelist} ] ${exp}_atm_2d_ml.$var.tmp
  # Select time period - and get rid of 1st timestep (Dec 31 of previous year)
  $CDO selyear,$yr1/$yr2 ${exp}_atm_2d_ml.$var.tmp ${exp}_atm_2d_ml.$var.nc
  rm ${exp}_atm_2d_ml.$var.tmp
done

#
# 2. Is global PME close to zero in multy year average?
# -----------------------------------------------------
echo ""
echo "Atmosphere:"
# Global PME on atmospheric grid [kg m-2 s-1]
$CDO -setvar,pme -add ${exp}_atm_2d_ml.pr.nc ${exp}_atm_2d_ml.evspsbl.nc pme_A.nc
$CDO -fldmean pme_A.nc pme_A.fldmean.nc
$CDO timavg -setunit,"m3 s-1" -mulc,0.001 -mul pme_A.fldmean.nc global_area.nc pme_A-mismatch_${yr1}-${yr2}.nc
get_sea_level_rise pme_A-mismatch_${yr1}-${yr2}.nc
wbal_change=$($CDO output pme_A-mismatch_${yr1}-${yr2}.nc)
echo "Global mean PME average from $yr1 to $yr2 (close to zero in equilibrium):"
echo "  wbal_change in [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

#
# 3. Same from ocean perspective: Is net fresh water flux to ocean close to zero?
#
echo ""
echo "Ocean:"
for var in FrshFlux_Runoff FrshFlux_Precipitation FrshFlux_Evaporation; do
  filelist=""
  for yr in $yearlist; do
    filelist="${filelist} ${outdata}/${exp}_oce_P1M_2d_${yr}0101.nc"
  done
  $CDO -O mergetime -apply,-selvar,${var} [ ${filelist} ] ${exp}_oce_P1M_2d_${var}.nc
done
# PME over ocean (on ocean grid) [m/s] -> [m3/s]
$CDO add ${exp}_oce_P1M_2d_FrshFlux_Precipitation.nc ${exp}_oce_P1M_2d_FrshFlux_Evaporation.nc pme_oce_O.nc
$CDO -timavg -setunit,"m3/s" -fldsum -mul pme_oce_O.nc area_O.nc pme_oce_O.fldsum.nc

# River inflow on ocean grid [m/s] -> [m3/s]
$CDO -timavg -setunit,"m3/s" -fldsum -mul ${exp}_oce_P1M_2d_FrshFlux_Runoff.nc area_O.nc runoff_O.fldsum.nc

# Global mean freshwater fluxes to the ocean shoud be close to zero in multi-year simulations.
$CDO add runoff_O.fldsum.nc pme_oce_O.fldsum.nc frshflux-mismatch_${yr1}-${yr2}.nc
get_sea_level_rise frshflux-mismatch_${yr1}-${yr2}.nc
wbal_change=$($CDO output frshflux-mismatch_${yr1}-${yr2}.nc)
echo "Global mean freshwater flux to the ocean (should be close to zero in equilibrium):"
echo "  wbal_change in [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

#
# 4. Does global land PME correspond to global ocean discharge?
# -------------------------------------------------------------
echo ""
echo "Land:"
# Global sum of PME over land [kg m-2 s-1] -> [m3/s]
$CDO add -mul ${exp}_atm_2d_ml.pr.nc notsea_gt0.nc \
              ${exp}_lnd_basic_ml.hydro_evapotrans_box.nc pme_land.nc
$CDO -setunit,"m3 s-1" -mulc,0.001 -fldsum -mul pme_land.nc \
       -mul notsea.nc area_A.nc   pme_land.fldsum.nc

# Global sum of ocean discharge [m3/s]
$CDO fldsum ${exp}_lnd_basic_ml.hydro_discharge_ocean_box.nc disch-oce.fldsum.nc

# In multi year mean, PME over land and global discharge should be identical
$CDO sub disch-oce.fldsum.nc pme_land.fldsum.nc pme-disch-oce-mismatch.nc
$CDO timavg pme-disch-oce-mismatch.nc pme-disch-oce-mismatch_${yr1}-${yr2}.nc
get_sea_level_rise pme-disch-oce-mismatch_${yr1}-${yr2}.nc
wbal_change=$($CDO output pme-disch-oce-mismatch_${yr1}-${yr2}.nc)
echo "Mismatch between global land PME and global discharge to ocean:"
echo "  wbal_change in [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

# What causes the Land water balance mismatch?
#
# 4a. Does global land PME corresponds to global runoff and drainage?
# ------------------------------------------------------------------
# Sum up runoff and drainage [kg m-2 s-1]
$CDO add ${exp}_lnd_basic_ml.hydro_runoff_box.nc ${exp}_lnd_basic_ml.hydro_drainage_box.nc \
    runoff-plus-drainage.nc

# In long simulations PME on land should correspond to runoff and drainage.
$CDO timavg -sub runoff-plus-drainage.nc pme_land.nc pme-run-drain-missmatch.nc
# Global sum [m3/s]
$CDO -setunit,"m3 s-1" -mulc,0.001 -fldsum \
     -mul pme-run-drain-missmatch.nc -mul notsea.nc area_A.nc \
          pme-run-drain-mismatch_${yr1}-${yr2}.nc
get_sea_level_rise pme-run-drain-mismatch_${yr1}-${yr2}.nc
wbal_change=$($CDO output pme-run-drain-mismatch_${yr1}-${yr2}.nc)
echo "  a) Mismatch between global land PME and runoff+drainage:"
echo "     wbal_change in [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

#
# 4b. Does global runoff and drainage correspond to global ocean discharge?
# ------------------------------------------------------------------------
# Global runoff and drainage (on land fraction) [kg m-2 s-1] -> [m3/s]
$CDO -mulc,0.001 -fldsum -mul runoff-plus-drainage.nc \
     -mul notsea.nc area_A.nc   runoff-plus-drainage.fldsum.nc

# Global ocean discharge [m3/s]
$CDO fldsum ${exp}_lnd_basic_ml.hydro_discharge_ocean_box.nc discharge_ocean.fldsum.nc

# In multi year mean, runoff and drainage over land should match global discharge
$CDO timavg -sub discharge_ocean.fldsum.nc runoff-plus-drainage.fldsum.nc \
       run-drain-disch-mismatch_${yr1}-${yr2}.nc
get_sea_level_rise run-drain-disch-mismatch_${yr1}-${yr2}.nc
wbal_change=$($CDO output run-drain-disch-mismatch_${yr1}-${yr2}.nc)
echo "  b) Mismatch between global runoff+drainage and ocean discharge:"
echo "     wbal_change in [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

#
# 4c. Does the mismatch correspond to jsbach total land water content change?
#
# Global land water budget change  [m3/nyrs] -> [m3/s]
echo "Global land water content change:"
(( yrn = yr2 + 1 ))   # next year: restart file at the end of yr2
$CDO -sub -fldsum -selvar,hydro_weq_land_box ${restdir}/${exp}_restart_atm_${yrn}0101.nc \
          -fldsum -selvar,hydro_weq_land_box ${restdir}/${exp}_restart_atm_${yr1}0101.nc \
     weq_land_change.nc
$CDO -divc,$nsec weq_land_change.nc weq_land_change_${yr1}-${yr2}.nc
get_sea_level_rise weq_land_change_${yr1}-${yr2}.nc
wbal_change=$($CDO output weq_land_change_${yr1}-${yr2}.nc)
echo "  Land (without HD) [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

# Global HD water budget change  [m3/nyrs] -> [m3/s]
$CDO sub -fldsum -selvar,hd_water_budget_box ${restdir}/${exp}_restart_atm_${yrn}0101.nc \
         -fldsum -selvar,hd_water_budget_box ${restdir}/${exp}_restart_atm_${yr1}0101.nc \
     hd_budget_change.nc
$CDO -divc,$nsec hd_budget_change.nc hd_budget_change_${yr1}-${yr2}.nc
get_sea_level_rise hd_budget_change_${yr1}-${yr2}.nc
wbal_change=$($CDO output hd_budget_change_${yr1}-${yr2}.nc)
echo "  HD reservoirs     [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

#
# 5. Freshwater fluxes to the ocean on atmosphere and ocean grid
# --------------------------------------------------------------
echo ""
echo "Mismatch due to remapping between atmosphere and ocean grids:"
#
# 5a PME
# -------
# PME on ocean grid [m/s] -> [m3/s]
$CDO add ${exp}_oce_P1M_2d_FrshFlux_Precipitation.nc ${exp}_oce_P1M_2d_FrshFlux_Evaporation.nc pme_oce_O.nc
$CDO setunit,"m3/s" -fldsum -mul pme_oce_O.nc area_O.nc pme_oce_O.fldsum.nc

# PME over ocean on atm grid
# As PME ocean is not in the basic atm output we need to calculate it:
#  pme_A = sea * pme_ocean + notsea * pme_land
#    => pme_ocean = (pme_A - notsea * pme_land) / sea
$CDO -div -sub pme_A.nc -mul pme_land.nc notsea.nc ocean_fract.nc pme_oce_A.nc
# Global sum [kg m-2 s-1] -> [m3/s]
$CDO setunit,"m3/s" -fldsum -mulc,0.001 \
       -mul pme_oce_A.nc -mul ocean_fract.nc area_A.nc pme_oce_A.fldsum.nc

# Mismatch O-A
$CDO -timavg -sub pme_oce_O.fldsum.nc pme_oce_A.fldsum.nc pme-O-A-mismatch_${yr1}-${yr2}.nc
get_sea_level_rise pme-O-A-mismatch_${yr1}-${yr2}.nc
wbal_change=$($CDO output pme-O-A-mismatch_${yr1}-${yr2}.nc)
echo "  PME          [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"

#
# 5b River discharge
# ------------------
# River inflow on ocean grid [m/s] -> [m3 s-1]
$CDO -setunit,"m3 s-1" -mul ${exp}_oce_P1M_2d_FrshFlux_Runoff.nc area_O.nc runoff_O.nc
$CDO fldsum runoff_O.nc runoff_O.fldsum.nc
# River discharge on atm grid [m3 s-1]
$CDO fldsum ${exp}_lnd_basic_ml.hydro_discharge_ocean_box.nc runoff_A.fldsum.nc
# Mismatch [m3 s-1]
$CDO timavg -sub runoff_O.fldsum.nc runoff_A.fldsum.nc runoff-mismatch_${yr1}-${yr2}.nc
get_sea_level_rise runoff-mismatch_${yr1}-${yr2}.nc
wbal_change=$($CDO output runoff-mismatch_${yr1}-${yr2}.nc)
echo "  Disch. Ocean [m3 s-1]: ${wbal_change} -> sea_level_rise in $nyrs years [mm]: ${sea_level_rise}"
