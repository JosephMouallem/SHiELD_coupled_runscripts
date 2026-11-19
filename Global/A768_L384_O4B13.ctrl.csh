#!/bin/tcsh
#SBATCH --output=./stdout/%x.%j
#SBATCH --job-name=C768_coupled
#SBATCH --clusters=c6
#SBATCH --time=1:59:00
#SBATCH --nodes=50
#SBATCH --exclusive
#SBATCH --account=bil-coastal-gfdl
##SBATCH --mail-user=joseph.mouallem@noaa.gov
##SBATCH --mail-type=ALL
##SBATCH --export=CLU=c6,DAYS=9,NAME=20240926.00Z,MEMO=LM4,LX=32,LY=32,NT=1,ALL

# Script to run the fully coupled SHiELD atmosphere-land-ocean-sea ice model.
# Configuration: C768 atmosphere, C384 LM4 land, and the OM5 b13 MOM6 ocean
# and SIS2 sea ice configuration (A768_L384_O4b13).
# Stages input data, writes runtime tables and namelists, and runs the model.
# Default experiment: 2024-09-26 00 UTC, integrated for 3 days. (Helene)
# Before use, review the Slurm allocation, account, paths, and run settings.
# @developper <mouallem@princeton.edu>
# version:
# 11/18/2026 - Initial Version

set echo

# Experiment settings
set DAYS = "3"
set NAME = "20240926.00Z"
set LX = "32"
set LY = "32"
set NT = "1"  # Threads per MPI rank
set ocean_ice_IC = ".F."
echo "Cluster: $SLURM_CLUSTER_NAME"

# Site paths: customize for your account and checkout.
set YourGroup = "bil-coastal-gfdl"

if (${SLURM_CLUSTER_NAME} == "c6") then
  set BASE_DIR   = "/gpfs/f6/${YourGroup}/scratch/${USER}/"  # Work directory
  set BUILD_DIR  = "${HOME}/shiemom_update/SHiELD_build/"  # Build directory
  set INPUT_DATA = "/gpfs/f6/bil-coastal-gfdl/proj-shared/gfdl_w/SHiELD_INPUT_DATA/"
endif

if (${SLURM_CLUSTER_NAME} == "c7") then
  set BASE_DIR   = "/gpfs/f7/${YourGroup}/scratch/${USER}/"  # Work directory
  set BUILD_DIR  = "${HOME}/c7/SHiELD_build/"  # Build directory
  set INPUT_DATA = "/gpfs/f7/bil-coastal-gfdl/proj-shared/gfdl_w/SHiELD_INPUT_DATA/"
endif

set RUN_DIR = "/ncrc/home2/${USER}/shiemom_dev/run/"
set COUPLED_FILES = ${INPUT_DATA}/Coupled_SHiELD/
# Used only when ocean_ice_IC = ".T."; edit for your initialization dataset.
set OCEAN_ICE_IC_DIR = "/gpfs/f6/bil-coastal-gfdl/scratch/${USER}/MOM_ICs_OBCs/temp/ICs/OM5/NK75"

if ( ! $?COMPILER ) then
  set COMPILER = "intel"
endif

set RELEASE = "test_fully_coupled"

unset echo
source ${BUILD_DIR}/site/environment.${COMPILER}.csh
set echo

set res = 768

set days = ${DAYS}
# Model build and experiment naming
set TYPE = "nh"         # choices:  nh, hydro
set MODE = "64bit"      # choices:  32bit, 64bit
set MONO = "non-mono"   # choices:  mono, non-mono
set CASE = "C$res"
set MEMO = "$SLURM_JOB_NAME.${days}d_1hr_diag.coldstart_ref"
set EXE = "x"
set HYPT = "off"         # choices:  on, off  (controls hyperthreading)
set COMP = "repro"       # choices:  debug, repro, prod
set NO_SEND = "no_send"    # choices:  send, no_send
set RESTART_RUN = "F"

# Run directory and executable
set WORKDIR    = ${BASE_DIR}/SHiELD_${RELEASE}/${NAME}.${TYPE}.${COMP}.${MODE}.${COMPILER}.${MONO}.${MEMO}/
set executable = ${BUILD_DIR}/Build/bin/SHiEMOM_LM4_${TYPE}.${COMP}.${MODE}.${COMPILER}.${EXE}

# Atmospheric input files
set ICDIR   = ${INPUT_DATA}/global.v201810/${CASE}_test/IC/${CASE}/${NAME}_IC/
set ICS  = ${ICDIR}/GFS_INPUT.tar
set FIX  = ${INPUT_DATA}/fix.v201810/
set GFS  = ${INPUT_DATA}/GFS_STD_INPUT.20160311.tar
set GRID = ${INPUT_DATA}/global.v201810/${CASE}/GRID/

set DIAG_TABLE  = ${RUN_DIR}/tables/diag_table_6species_tc_dp_ocean
set FIELD_TABLE = ${RUN_DIR}/tables/field_table_6species_atmland_tke # Must match the selected physics

# Timestamp helper for optional archiving
set TIME_STAMP = ${BUILD_DIR}/site/time_stamp.csh

# Grid and parallel decomposition
set npx = "769"
set npy = "769"
set npz = "79"

set layout_x = $LX
set layout_y = $LY
set layout_x_land = $layout_x
set layout_y_land = $layout_y

set atmos_npes = `expr $layout_x \* $layout_y \* 6`
set land_npes = `expr $layout_x_land \* $layout_y_land \* 6`

set io_layout = "1,1"
set nthreads = $NT
# Blocking factor for threading and physics performance
set blocksize = "36"

# Integration length; dt_atmos and dt_therm are in seconds
set months = "0"
set hours = "1"  # Added to DAYS; the default total is 73 hours
set seconds = "0"
set dt_atmos = "150"
set dt_therm = "300"  # Ocean thermodynamic time step

# Optional YAML field/data tables
set use_yaml = ".F."  # Requires field_table.yaml and data_table.yaml when enabled

# set the pre-conditioning of the solution
# =0 implies no pre-conditioning
# >0 means new adiabatic pre-conditioning
# <0 means older adiabatic pre-conditioning
set na_init = 1

# NCEP/NGGPS atmospheric initialization
set filtered_terrain = ".true."
set ncep_levs = "64"
set gfs_dwinds = ".true."

# GFS diagnostic output, accumulation reset, and surface cycling intervals (hours)
set fdiag = "1."
set fhzer = "1."
set fhcyc = "24."

# Select whether GFS physics calculates geopotential
set gfs_phil = ".false."

# Diagnostic and debug options
set no_dycore = ".false."
set dycore_only = ".false."
set chksum_debug = ".false."
set print_freq = "6"

if (${TYPE} == "nh") then
  # Nonhydrostatic configuration
  set make_nh = ".T."
  set hydrostatic = ".F."
  set phys_hydrostatic = ".F."
  set use_hydro_pressure = ".F."
  set consv_te = "1."
  # FV3 time splitting
  set k_split = "1"
  set n_split = "8"
else
  # Hydrostatic configuration
  set make_nh = ".F."
  set hydrostatic = ".T."
  set phys_hydrostatic = ".F."     # will be ignored in hydro mode
  set use_hydro_pressure = ".T."   # have to be .T. in hydro mode
  set consv_te = "0."
  # FV3 time splitting
  set k_split = "2"
  set n_split = "6"
endif

if (${MONO} == "mono" || ${MONO} == "monotonic") then
  # Monotonic configuration
  set d_con = "1."
  set do_vort_damp = ".false."
else
  # Nonmonotonic configuration
  set d_con = "1."
  set do_vort_damp = ".true."
endif

# Hyperthreading and launch command
if (${HYPT} == "on") then
  set hyperthread = ".true."
  set div = 2
else
  set hyperthread = ".false."
  set div = 1
endif

@ npes = ${atmos_npes}
@ skip = ${nthreads} / ${div}
set run_cmd = "srun --label --ntasks=$npes --cpus-per-task=$skip ./$executable:t"

setenv SLURM_CPU_BIND verbose

setenv MPICH_ENV_DISPLAY
setenv MPICH_MPIIO_CB_ALIGN 2
setenv MALLOC_MMAP_MAX_ 0
setenv MALLOC_TRIM_THRESHOLD_ 536870912
setenv NC_BLKSZ 1M
# Intel OpenMP stack size
setenv KMP_STACKSIZE 256m
setenv FI_CXI_RX_MATCH_MODE software # C6 workaround for receive-resource exhaustion

# Prepare the run directory and initial state
if (${RESTART_RUN} == "F") then

  \rm -rf $WORKDIR/rundir

  mkdir -p $WORKDIR/rundir
  cd $WORKDIR/rundir

  mkdir -p RESTART
  mkdir -p INPUT

  # Atmospheric initial conditions for NAME
  if ( -e ${ICDIR}/gfs_data.tile1.nc ) then
    cp -rf ${ICDIR}/* INPUT/
  else
    tar xf ${ICS}
  endif

  # Cold-start namelist settings
  set nggps_ic = ".T."
  set mountain = ".F."
  set external_ic = ".T."
  set warm_start = ".F."

  set input_filename = 'n' # MOM6 and SIS2 input mode
else

  cd $WORKDIR/rundir
  \rm -rf INPUT/*

  # Move restart data into INPUT
  mv RESTART/* INPUT/.

  # Restart namelist settings
  set make_nh = ".F."
  set nggps_ic = ".F."
  set mountain = ".T."
  set external_ic = ".F."
  set warm_start = ".T."
  set na_init = 0
  set input_filename = 'r' # MOM6 and SIS2 input mode

endif

cp $executable .

# GFS standard input data
tar xf ${GFS}

# Grid and orography data
cp -rf ${GRID}/oro* INPUT/.
cp -rf ${COUPLED_FILES}/GLOBAL/A768_L384_O4b13/C768_mosaic/* INPUT/.

# River input files
cp -rf ${COUPLED_FILES}/GLOBAL/A768_L384_O4b13/river/* INPUT/.

# Parse NAME (YYYYMMDD.HHZ) for the model and diagnostics

set y = `echo ${NAME} | cut -c1-4`
set m = `echo ${NAME} | cut -c5-6`
set d = `echo ${NAME} | cut -c7-8`
set h = `echo ${NAME} | cut -c10-11`
set echo
set curr_date = "${y},${m},${d},${h},0,0"

# Diagnostic table header and date placeholders
cat >! diag_table << EOF
${NAME}.${CASE}.${MODE}.${MONO}
$y $m $d $h 0 0
EOF

cat ${DIAG_TABLE} >> diag_table
sed -i "s/YYYY MM DD CYC/${y} ${m} ${d} ${h}/g" diag_table # Resolve the diagnostic reference date

cp ${FIELD_TABLE} field_table

cp $FIX/global_sfc_emissivity_idx.txt INPUT/sfc_emissivity_idx.txt
cp INPUT/aerosol.dat .
cp INPUT/co2historicaldata_*.txt .
cp INPUT/sfc_emissivity_idx.txt .
cp INPUT/solarconstant_noaa_an.txt .

# MOM6 and SIS2 input files
set MOM_INPUT_DIR = "${COUPLED_FILES}/GLOBAL/om5-config-b13/OM5/"
cp $MOM_INPUT_DIR/MOM_* .
cp $MOM_INPUT_DIR/SIS_* .
ln -sf $MOM_INPUT_DIR/INPUT/MOM_channels_global_025 INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/layer_coord.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/hycom1_75_800m.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/woa13_decav_ptemp_monthly_fulldepth_01.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/woa13_decav_s_monthly_fulldepth_01.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/tidal_amplitude_v20250916.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/geothermal_davies2013_v1_v20250916.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/diag_z3.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/diag_rho2_74layer.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/diag_rho0_60layer.nc INPUT/
ln -sf $MOM_INPUT_DIR/INPUT/seawifs-clim-1997-2010_v20250916.nc INPUT/

# LM4 land input files
set MOMSISLM_INPUT_DIR = "${COUPLED_FILES}/"
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/geohydrology.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/geohydrology_table_2a2n.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/ground_type.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/cover_type.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/biodata.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/landuse.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/soil_brdf.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/soil_type_hwsd_5minute.nc INPUT/soil_type.nc
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/c384_topo_rough_land_expanded/* INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/landfrac.nc INPUT/
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/RCP4.5_SSP2_2001-2100_states_v20180926.nc INPUT/states.nc
ln -sf $MOMSISLM_INPUT_DIR/common_LM4/RCP4.5_SSP2_2001-2100_transitions_v20180926.nc INPUT/transitions.nc

# remaped land data if restarting from remaped land/ocean files
###############################################################3
#set SPEAR_INPUT_DIR = "/gpfs/f6/bil-coastal-gfdl/scratch/Joseph.Mouallem/SPEAR_INPUT_C384/SPEAR_INPUT_IC_2024_old/"
#
#ln -sf $SPEAR_INPUT_DIR/ic* INPUT/
#ln -sf $SPEAR_INPUT_DIR/M* INPUT/
#ln -sf $SPEAR_INPUT_DIR/riv* INPUT/
#ln -sf $SPEAR_INPUT_DIR/topo* INPUT/
#ln -sf $SPEAR_INPUT_DIR/cg* INPUT/
#
#set SPEAR_INPUT_DIR = "/gpfs/f6/bil-coastal-gfdl/scratch/Joseph.Mouallem/SPEAR_INPUT_C384_soilwater_sensitivity/SPEAR_INPUT_IC_2024_old_remap/" # for remaped land files
##set SPEAR_INPUT_DIR = "/gpfs/f6/bil-coastal-gfdl/scratch/Joseph.Mouallem/SPEAR_INPUT_C384/SPEAR_INPUT_IC_2024_old_remap_halfheight/" # for remaped land files
#
#ln -sf $SPEAR_INPUT_DIR/la* INPUT/
#ln -sf $SPEAR_INPUT_DIR/gla* INPUT/
##ln -sf $SPEAR_INPUT_DIR/riv* INPUT/
##ln -sf $SPEAR_INPUT_DIR/topo* INPUT/
#ln -sf $SPEAR_INPUT_DIR/soi* INPUT/
#ln -sf $SPEAR_INPUT_DIR/veg* INPUT/
#ln -sf $SPEAR_INPUT_DIR/ca* INPUT/
##ln -sf $SPEAR_INPUT_DIR/cg* INPUT/
#ln -sf $SPEAR_INPUT_DIR/sn* INPUT/

#ln -sf $SPEAR_INPUT_DIR/ic* INPUT/
#ln -sf $SPEAR_INPUT_DIR/M* INPUT/

unset echo

# MOM6 and SIS2 initial conditions and parameter overrides

if ( ${ocean_ice_IC} == ".T." ) then

  ln -sf ${OCEAN_ICE_IC_DIR}/MOM6_IC_${y}${m}${d}${h}_OM5.nc INPUT/
  ln -sf ${OCEAN_ICE_IC_DIR}/SIS2_IC_${y}${m}${d}_OM5.nc INPUT/

  cat >! MOM_override <<EOF
#override DT=${dt_atmos}
#override DT_THERM=${dt_therm}
#override ALE_COORDINATE_CONFIG = "HYBRID:hycom1_75_800m.nc,sigma2,dz"
#override INIT_LAYERS_FROM_Z_FILE = True
#override TEMP_Z_INIT_FILE = "MOM6_IC_${y}${m}${d}${h}_OM5.nc"
#override SALT_Z_INIT_FILE = "MOM6_IC_${y}${m}${d}${h}_OM5.nc"
#override Z_INIT_FILE_PTEMP_VAR = "temp" ! default = "temp"
#override Z_INIT_FILE_SALT_VAR = "salt"   ! default = "salt"
#override Z_INIT_ALE_REMAPPING = True     !   [Boolean] default = False
#override Z_INIT_REMAP_GENERAL = False     !   [Boolean] default = False
#override Z_INIT_REMAP_OLD_ALG = False    !   [Boolean] default = True
#override Z_INIT_REMAP_FULL_COLUMN = True
#override TEMP_SALT_INIT_VERTICAL_REMAP_ONLY = True
#override REMAP_UV_USING_OLD_ALG = False
#override REMAP_AFTER_INITIALIZATION = True
#override DEPRESS_INITIAL_SURFACE = True
#override SURFACE_HEIGHT_IC_FILE = "MOM6_IC_${y}${m}${d}${h}_OM5.nc"
#override SURFACE_HEIGHT_IC_VAR = "ssh"
#override VELOCITY_CONFIG = "file"
#override VELOCITY_FILE = "MOM6_IC_${y}${m}${d}${h}_OM5.nc" ! barotropic adjusted currents
EOF

  cat >! SIS_override <<EOF
CONCENTRATION_INIT_CONFIG = "file"
ICE_THICKNESS_INIT_CONFIG = "file"
ICE_CONCENTRATION_FILE = "SIS2_IC_${y}${m}${d}_OM5.nc"
ICE_CONCENTRATION_IC_VAR = "aice"
ICE_THICKNESS_FILE = "SIS2_IC_${y}${m}${d}_OM5.nc"
ICE_THICKNESS_IC_VAR = "hm"
EOF

else # Use the default ocean initialization in MOM_input

  cat >! MOM_override <<EOF
#override DT=${dt_atmos}
#override DT_THERM=${dt_therm}
EOF

  cat >! SIS_override <<EOF
EOF

endif

# Coupled model namelist
cat > input.nml <<EOF

&MOM_input_nml
  output_directory = './',
  input_filename = ${input_filename}
  restart_input_dir = 'INPUT/',
  restart_output_dir = 'RESTART/',
  parameter_filename = 'MOM_input','MOM_override'
/

&SIS_input_nml
  output_directory = './',
  input_filename = ${input_filename}
  restart_input_dir = 'INPUT/',
  restart_output_dir = 'RESTART/',
  parameter_filename = 'SIS_input','SIS_override'
/

&fms_affinity_nml
  affinity=.F.
/

&glac_data_nml
  use_lm2_awc = .true.,
  rsa_exp_global = 10.,
  f_iso_cold = 0.92, 0.73
  f_vol_cold = 0.06, 0.08
  f_geo_cold = 0.00, 0.00
  f_iso_warm = 0.92, 0.73
  f_vol_warm = 0.06, 0.08
  f_geo_warm = 0.00, 0.00
  dat_emis_dry = 1.,
  dat_emis_sat = 1.,
  geothermal_heat_flux_constant = 0.065
/

&glac_nml
  conserve_glacier_mass = .true.
  lm2 = .true.
  albedo_to_use = 'brdf-params'
/

&harvesting_nml
  do_harvesting          =.true.
  grazing_intensity      = 18.25 ! roughly 5% of biomass per day
  grazing_residue        = 0.1
  grazing_frequency      = 'daily'
  min_lai_for_grazing    = 1.0
  frac_wood_wasted_harv  = 0.25
  frac_wood_wasted_clear = 0.25
  crop_seed_density      = 0.1
  waste_below_ground_wood = .TRUE.
/

&hlsp_nml
  max_num_topo_hlsps = 1
  num_vertclusters = 10
  use_geohydrodata = .true.
  strm_depth_penetration = 1.5
  use_hlsp_aspect_in_gwflow = .true.
  diagnostics_by_cluster = .true.
  init_wt_strmelev = .false.
  dammed_strm_bc = .true.
  limit_intertile_flow = .true.
  flow_ratio_limit = 0.25
/

&soil_data_nml
  geohydrology_to_use    = 'hill'
  use_comp_for_ic        = .true.
  use_comp_for_push      = .true.
  use_tau_fix            = .false.
  use_sat_fix            = .true.
  use_alt3_soil_hydraulics = .true.
  gw_scale_soil_depth  = 0.1666667
  gw_scale_perm        = 1.0
  gw_scale_relief      = 0.5
  aspect           = 0.01
  k0_macro_z       = 0.0
  k0_macro_x       = 1.0
  freeze_factor    = 1.5
  comp             = 1e-4
  num_l            = 20
  dz = 0.020, 0.040, 0.040, 0.050, 0.05, 0.100, 0.10, 0.200, 0.200, 0.2,
  0.40,  0.40,  0.40,  0.4,   0.4,   1.0,   1.0,  1.0,  1.5,   2.5
  soil_type_file          = "INPUT/soil_type.nc"
  input_cover_types       =  1,       2,       3,       4,       5,       6,       7,       8,       9,       10,      11,      12,      13,       14
  tile_names              =  'hec',   'sic',   'lic',   'sicl',  'cl',    'si',    'sil',   'sac',   'l',     'sacl',  'sal',  'ls',     's',     'u'
  dat_w_sat               =  0.468,   0.468,   0.468,   0.464,   0.465,   0.476,   0.476,   0.406,   0.439,   0.404,   0.434,   0.421,   0.339,   0.439
  dat_k_sat_ref           =  0.00097, 0.0013,  0.00097, 0.002,   0.0024,  0.0028,  0.0028,  0.0072,  0.0033,  0.0045,  0.0052,  0.014,   0.047,   0.0033
  dat_psi_sat_ref         = -0.47,   -0.32,   -0.47,   -0.62,   -0.26,   -0.76,   -0.76,   -0.098,  -0.35,   -0.13,   -0.14,   -0.036   -0.069,  -0.35
  dat_chb                 =  12.0,    10.0,    12.0,    8.7,     8.2,     5.3,     5.3,     11.0,    5.2,     6.8,     4.7,     4.3,     2.8,     5.2
  dat_heat_capacity_dry   =  1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6,   1.1e6
  dat_thermal_cond_dry    =  0.21,    0.21,    0.21,    0.21,    0.21,    0.21,    0.21,    0.21,    0.21,    0.21,    0.21,    0.21,    0.14,    0.21
  dat_thermal_cond_sat    =  1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     2.3,     1.5
  dat_thermal_cond_exp    =  6,       6,       6,       5,       5,       5,       5,       6,       5,       5,       5,       5,       3,       5
  dat_thermal_cond_scale  =  10,      10,      10,      0.5,     0.5,     0.5,     0.5,     10,      0.5,     0.5,     0.5,     0.5,     15,      0.5
  dat_thermal_cond_weight =  0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.7,     0.2,     0.7
  dat_emis_dry            =  1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0
  dat_emis_sat            =  1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0
  dat_tf_depr             =  2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0,     2.0
  dat_z0_momentum         =  0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01,    0.01
  Dpsi_min_const      = -1.e12
  K_min               = 1.e-12
  geothermal_heat_flux_constant = 0.065
  retro_a0N1 = .true.
/

&soil_nml
  init_wtdep    = 10.
  albedo_to_use = 'brdf-maps'
  uptake_to_use = 'darcy2d-linearized'
  uptake_oneway = .true.
  active_layer_drainage_acceleration = 1.
  gw_flux_max   = 0.001
  lrunf_from_div = .false.
  cold_infilt   = .false.
  use_depth_to_wt_4 = .true.
  layer_for_gw_switch = 16
  bottom_up_cold_infilt = .true.
  cold_depth = 3.0
  bwood_macinf = 1.
  allow_neg_wl = .true.
  prohibit_negative_water_div = .true.
  fix_neg_subsurface_wl = .False.
/

&river_nml
  dt_slow     = 150.
  land_area_called_cellarea = .true.
  all_big_outlet_ctn0 = .true.
  ave_DHG_exp = 0.55, 0.45, 0.0
  ave_DHG_coef = 7.5, 0.3, 0.444444444
  ave_AAS_exp = 0.26, 0.40, 0.34
  stop_on_mask_mismatch=.F.
/

&river_physics_nml
  lake_sfc_w_min = 10.
  prohibit_cold_ice_outflow = .FALSE.
/
&snow_data_nml
  emis_snow_max = 1.
  emis_snow_min = 1.
  z0_momentum = 0.01
  num_l = 5
  dz = 0.05, 0.2, 0.5, 0.2, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0
  f_iso_cold = 0.92, 0.58
  f_vol_cold = 0.06, 0.08
  f_geo_cold = 0.00, 0.00
  f_iso_warm = 0.77, 0.43
  f_vol_warm = 0.06, 0.08
  f_geo_warm = 0.00, 0.00
  depth_crit = 0.0125
  thermal_cond_ref = 0.09
  distinct_snow_on_glacier = .TRUE.
  f_iso_cold_on_glacier = 0.92, 0.68
  f_vol_cold_on_glacier = 0.06, 0.08
  f_geo_cold_on_glacier = 0.00, 0.00
  f_iso_warm_on_glacier = 0.77, 0.43
  f_vol_warm_on_glacier = 0.06, 0.08
  f_geo_warm_on_glacier = 0.00, 0.00
/

&snow_nml
  max_snow = 3000
  min_snow_mass = 1.e-7
  albedo_to_use = 'brdf-params'
  snow_density = 250.
  prevent_tiny_snow = .true.
/

&ice_spec_nml
  mcm_ice = .false.
  do_leads = .false.
  minimum_ice_concentration = 0.5
  minimum_ice_thickness = 0.1
  sst_degk = .true.
/

&icebergs_nml
  make_calving_reproduce=.TRUE.
  really_debug=.FALSE.
  debug=.FALSE.
  verbose=.FALSE.
  verbose_hrs=7200
  use_operator_splitting=.TRUE.
  bergy_bit_erosion_fraction=0.0
  sicn_shift=0.1
  parallel_reprod=.TRUE.
  traj_sample_hrs=0
  add_weight_to_ocean=.false.
  tidal_drift = 0.005
  coastal_drift = 0.001
/

&interpolator_nml
  retain_cm3_bug = .false.
/

&lake_data_nml
  lake_to_use = 'from-rivers'
  lake_specific_width = .false.
  large_lake_sill_width = 800.
  num_l = 20,
  dat_emis_dry = 1.,
  dat_emis_sat = 1.,
  dat_heat_capacity_ref = 0.,
  f_iso_liq = 0.025, 0.01
  f_iso_ice = 0.025, 0.01
  f_vol_liq = 0., 0.
  f_vol_ice = 0., 0.
  f_geo_liq = 0., 0.
  f_geo_ice = 0., 0.
  dry_lake_depth_frac = 0.5
  dat_z0_momentum(1)     = 1.4e-4
  dat_z0_momentum_ice(1) = 1.4e-4
  k_over_B_ice           = 0.25
  geothermal_heat_flux_constant = 0.065
/

&lake_nml
  albedo_to_use = 'brdf-params'
  float_ice_to_top = .true.
  lake_depth_max = 1.e10
  max_plain_slope = 0.0
  K_z_background  = 15.e-5
/

&land_conservation_nml
  do_check_conservation = .FALSE.
/

&cana_nml
  turbulence_to_use = "lm3v"
  init_co2 = 286.0e-6
  canopy_air_mass_for_tracers = 10.0
  allow_small_z0 = .true.
  sai_turb = .true.
/

&land_debug_nml
/

&land_model_nml
  layout   = $layout_x_land,$layout_y_land
  io_layout = $io_layout
  tau_snow_T_adj = 604800.
  prohibit_negative_canopy_water=.false.
  min_frac = 1e-8
/

&landuse_nml
  do_landuse_change = .F.
  data_type   = 'luh2'
  distribute_transitions = 'min-n-tiles'
  input_file  = 'INPUT/transitions.nc'
  state_file  = 'INPUT/states.nc'
  static_file = 'INPUT/landfrac.nc'
/

&vegn_data_nml
  vegn_to_use = 'uniform'
  K1=10,
  K2=0.1,
  fsc_liv=0.9,
  fsc_wood=0.45,
  c1(4) = 0.3
  c2(4) = 0.3
  Vmax = 2.0E-5, 2.0E-5, 2.0E-5, 2.0E-5, 1.50E-5,
  m_cond = 4., 9., 9., 7., 7.,
  alpha_phot = 0.05, 0.06, 0.06, 0.06, 0.06,
  gamma_resp = 0.03, 0.02, 0.02, 0.02, 0.02,
  tc_crit(0:2) = 3*273.16
  fact_crit_phen(0:4) = 0., 0., 0., 0., 0.
  fact_crit_fire(0:4) = 0., 0., 0., 0., 0.
  cnst_crit_phen(0:4) = 0.30, 0.15, 0.15, 0.30, 0.30
  cnst_crit_fire(0:4) = 0.15,  0.40, 0.15, 0.15, 0.15
  wet_leaf_dreg(0:4)  = .3, .3, .3, .3, .3
  ksi =  0, 0, 0, 0, 0,
  leaf_refl(0:4,1) = 0.11, 0.11, 0.10, 0.10, 0.10
  leaf_refl(0:4,2) = 0.58, 0.58, 0.45, 0.45, 0.39,
  dat_root_zeta(0:4) =   0.35212, 0.17039, 0.28909, 0.25813, 0.17039
  critical_root_density = 0.0,
  tau_drip_s     = 259200.0
  cmc_lai(0:4) = 0.02, 0.02, 0.02, 0.02, 0.02
  csc_lai(0:4) = 0.30, 0.30, 0.30, 0.30, 0.2
  dat_snow_crit = 2*1.e7, .1, 1.e7, .1
  t_transp_min = 268.
  srl(0:1) = 112.0e3, 150.0e3
  root_perm = 14*5e-7
  alpha(1,3) = 4
  leaf_age_tau(2) = 150
  smoke_fraction = 0.9, 0.9, 0.6, 0.6, 0.6
  tg_c3_thresh = 1
  phen_ev1 = 0.2
  phen_ev2 = 0.7
  cmc_eps = 0.01
  alpha(0:4,6) = 0.0, 0.0, 0.012, 0.012, 0.012
  treefall_disturbance_rate = 0.175, 0.185, 0.025, 0.0275, 0.027
  use_light_saber = .true.
  sai_cover = .true.
  sai_rad = .true.
  min_lai = 1.e-5   ! h1g, prevent numerical instability, 2017-08-24
/

&vegn_nml
  init_Tv = 288.0
  rad_to_use = 'two-stream'
  snow_rad_to_use = 'paint-leaves'
  do_cohort_dynamics   = .TRUE.,
  do_patch_disturbance = .TRUE.,
  do_phenology         = .TRUE.,
  do_biogeography      = .TRUE.,
  do_seed_transport    = .TRUE.,
  tau_smooth_ncm = 22.0
  rav_lit_bwood = 300.
/

&lscale_cond_nml
  do_evap = .true.
/

&sat_vapor_pres_nml
  construct_table_wrt_liq = .true.,
  construct_table_wrt_liq_and_ice = .true.,
  use_exact_qs = .true.
/

&monin_obukhov_nml
  neutral = .false. ! KGao: should not use true for coupled mode
/

&surface_flux_nml
  ncar_ocean_flux = .false. ! KGao: must be false for the coupled mode
  raoult_sat_vap = .true.
  do_iter_monin_obukhov = .true. ! KGao: turn this on for hwrf17 rough scheme
  niter_monin_obukhov = 2 ! KGao: 2 should work great
/
&xgrid_nml
  make_exchange_reproduce = .true.
  interp_method = 'first_order'
  monotonic_exchange = .false.
/

&atmos_model_nml
  blocksize = $blocksize
  chksum_debug = $chksum_debug
  dycore_only = $dycore_only
  fdiag = $fdiag
  fullcoupler_fluxes = 2
/

&fms_io_nml
  checksum_required   = .false.
  max_files_r = 100,
  max_files_w = 100,
/

&fms_nml
  clock_grain = 'ROUTINE',
  domains_stack_size = 300000000,
  print_memory_usage = .false.
/

&fv_grid_nml
  grid_file = 'INPUT/grid_spec.nc'
/

&fv_core_nml
  layout   = $layout_x,$layout_y
  io_layout = $io_layout
  npx      = $npx
  npy      = $npy
  ntiles   = 6
  npz    = $npz
  grid_type = -1
  make_nh = $make_nh
  fv_debug = .F.
  range_warn = .F.
  reset_eta = .F.
  n_sponge = 30
  nudge_qv = .T.
  rf_fast = .F.
  tau = 5.
  rf_cutoff = 7.5e2
  d2_bg_k1 = 0.15
  d2_bg_k2 = 0.02
  kord_tm = -9
  kord_mt =  9
  kord_wz =  9
  kord_tr =  9
  hydrostatic = $hydrostatic
  phys_hydrostatic = $phys_hydrostatic
  use_hydro_pressure = $use_hydro_pressure
  beta = 0.
  a_imp = 1.
  p_fac = 0.1
  k_split  = $k_split
  n_split  = $n_split
  nwat = 6
  na_init = $na_init
  d_ext = 0.0
  dnats = 1
  fv_sg_adj = 600
  d2_bg = 0.
  nord =  3
  dddmp = 0.2
  d4_bg = 0.15
  vtdm4 = 0.03
  delt_max = 0.002
  ke_bg = 0.
  do_vort_damp = $do_vort_damp
  external_ic = $external_ic
  gfs_phil = $gfs_phil
  nggps_ic = $nggps_ic
  mountain = $mountain
  ncep_ic = .F.
  d_con = $d_con
  hord_mt = 5
  hord_vt = 5
  hord_tm = 5
  hord_dp = -5
  hord_tr = -5
  adjust_dry_mass = .F.
  consv_te = $consv_te
  consv_am = .F.
  fill = .T.
  dwind_2d = .F.
  print_freq = $print_freq
  warm_start = $warm_start
  no_dycore = $no_dycore
  z_tracer = .T.
/

&integ_phys_nml
  do_inline_mp = .F.
  do_sat_adj = .F.
/

&coupler_nml
  months = $months
  days  = $days
  hours = $hours
  seconds = $seconds
  dt_atmos = $dt_atmos
  current_date =  $curr_date
  calendar = 'julian'
  atmos_nthreads = $nthreads
  use_hyper_thread = $hyperthread
  do_ocean=.True.
  dt_cpld = $dt_atmos
  do_flux=.True.
  do_land=.True.
  do_ice=.True.

  ! concurrent=.True.
  ! atmos_npes = $atmos_npes
  ! land_npes  = $land_npes
  !  atmos_npes = 1152
  !  ocean_npes=72

/

&external_ic_nml
  filtered_terrain = $filtered_terrain
  levp = $ncep_levs
  gfs_dwinds = $gfs_dwinds
  checker_tr = .F.
  nt_checker = 0
/

&gfs_physics_nml
  fhzero         = $fhzer
  ldiag3d        = .F.
  fhcyc          = $fhcyc
  nst_anl        = .true.
  use_ufo        = .true.
  pre_rad        = .false.
  ncld           = 5
  zhao_mic       = .false.
  pdfcld         = .false.
  fhswr          = 3600.
  fhlwr          = 3600.
  ialb           = 1
  iems           = 1
  IAER           = 111
  ico2           = 2
  isubc_sw       = 2
  isubc_lw       = 2
  isol           = 2
  lwhtr          = .true.
  swhtr          = .true.
  cnvgwd         = .true.
  do_deep        = .true.
  shal_cnv       = .true.
  cal_pre        = .false.
  redrag         = .true.
  dspheat        = .true.
  hybedmf        = .false.
  random_clds    = .false.
  trans_trac     = .true.
  cnvcld         = .false.
  imfshalcnv     = 2
  imfdeepcnv     = 2
  cdmbgwd        = 3.5, 0.25
  prslrd0        = 0.
  ivegsrc        = 1
  isot           = 1
  ysupbl         = .F.
  satmedmf       = .T.
  isatmedmf      = 1
  rlmx           = 150.
  cs0            = 0.25
  clam_deep      = 0.1
  clam_shal      = 0.1
  xkzminv        = 1.0
  xkzm_m         = 1.0
  xkzm_h         = 1.0
  cloud_gfdl     = .true.
  do_ocean       = .false.

/

! &ocean_nml
!     mld_option       = "obs"
!     ocean_option     = "MLM"
!     restore_method   = 2
!     mld_obs_ratio    = 1.
!     use_rain_flux    = .true.
!     sst_restore_tscale = 2.
!     start_lat        = -30.
!     end_lat          = 30.
!     Gam              = 0.2
!     use_old_mlm      = .true.
!     do_mld_restore   = .true.
!     mld_restore_tscale = 2.
!     stress_ratio     = 1.
!     eps_day          = 10.
!/

&gfdl_mp_nml
  do_sedi_heat = .true.
  rad_snow = .true.
  rad_graupel = .true.
  rad_rain = .true.
  const_vi = .false.
  const_vs = .false.
  const_vg = .false.
  const_vr = .false.
  vi_fac = 1.
  vs_fac = 1.
  vg_fac = 1.
  vr_fac = 1.
  vi_max = 1.
  vs_max = 2.
  vg_max = 12.
  vr_max = 12.
  qi_lim = 1.
  prog_ccn = .false.
  do_qa = .true.
  tau_l2v = 225.
  tau_v2l = 150.
  rthresh = 10.e-6  ! This is a key parameter for cloud water
  dw_land  = 0.16
  dw_ocean = 0.10
  ql_mlt = 1.0e-3
  qi0_crt = 8.0E-5
  qs0_crt = 1.0e-3
  tau_i2s = 1000.
  c_psaci = 0.05
  c_pgacs = 0.01
  rh_inc = 0.30
  rh_inr = 0.30
  ccn_l = 300.
  ccn_o = 100.
  c_paut = 0.5
  z_slope_liq  = .true.
  z_slope_ice  = .true.
  fix_negative = .false.
  icloud_f = 0
/

&cloud_diagnosis_nml
  ql0_max = 2.0e-3
  qi0_max = 2.0e-4
  ccn_o = 100.
  ccn_l = 300.
  qmin = 1.0e-12
  beta = 1.22
  rewflag = 1
  reiflag = 1
  rewmin = 5.0
  rewmax = 10.0
  reimin = 10.0
  reimax = 150.0
  liq_ice_combine = .true.
/

&diag_manager_nml
  prepend_date = .F.
  mix_snapshot_average_fields = .false.,
  max_input_fields = 800,
  max_output_fields = 1500
  max_axes = 100
  max_num_axis_sets = 50
  max_files = 40
  issue_oor_warnings = .false.
  max_axis_attributes=4
/

&interpolator_nml
  interp_method = 'conserve_great_circle'
/

&namsfc
  FNGLAC   = "$FIX/global_glacier.2x2.grb",
  FNMXIC   = "$FIX/global_maxice.2x2.grb",
  FNTSFC   = "$FIX/RTGSST.1982.2012.monthly.clim.grb",
  FNMLDC   = ""
  FNSNOC   = "$FIX/global_snoclim.1.875.grb",
  FNZORC   = "igbp",
  FNALBC   = "$FIX/global_snowfree_albedo.bosu.t1534.3072.1536.rg.grb",
  FNALBC2  = "$FIX/global_albedo4.1x1.grb",
  FNAISC   = "$FIX/CFSR.SEAICE.1982.2012.monthly.clim.grb",
  FNTG3C   = "$FIX/global_tg3clim.2.6x1.5.grb",
  FNVEGC   = "$FIX/global_vegfrac.0.144.decpercent.grb",
  FNVETC   = "$FIX/global_vegtype.igbp.t1534.3072.1536.rg.grb",
  FNSOTC   = "$FIX/global_soiltype.statsgo.t1534.3072.1536.rg.grb",
  FNSMCC   = "$FIX/global_soilmgldas.t1534.3072.1536.grb",
  FNMSKH   = "$FIX/seaice_newland.grb",
  FNTSFA   = "",
  FNACNA   = "",
  FNSNOA   = "",
  FNVMNC   = "$FIX/global_shdmin.0.144x0.144.grb",
  FNVMXC   = "$FIX/global_shdmax.0.144x0.144.grb",
  FNSLPC   = "$FIX/global_slope.1x1.grb",
  FNABSC   = "$FIX/global_mxsnoalb.uariz.t1534.3072.1536.rg.grb",
  LDEBUG   =.false.,
  FSMCL(2) = 99999
  FSMCL(3) = 99999
  FSMCL(4) = 99999
  FTSFS    = 90
  FAISS    = 99999
  FSNOL    = 99999
  FSICL    = 99999
  FTSFL    = 99999,
  FAISL    = 99999,
  FVETL    = 99999,
  FSOTL    = 99999,
  FvmnL    = 99999,
  FvmxL    = 99999,
  FSLPL    = 99999,
  FABSL    = 99999,
  FSNOS    = 99999,
  FSICS    = 99999,
/
EOF

if ( ${use_yaml} == ".T." ) then
  cat >> input.nml << EOF

&field_manager_nml
  use_field_table_yaml = $use_yaml
/

&data_override_nml
  use_data_table_yaml = $use_yaml
/
EOF
endif

# Run the model
${run_cmd} | tee fms.out
set run_status = $status
if ( $run_status != 0 ) then
  echo "ERROR: model launch or output capture failed (status $run_status)."
  exit $run_status
endif

# Site-specific postprocessing (requires scripts under HOME).
cd $WORKDIR/rundir
~/land_tools/decompland.csh && ~/bin/pp.csh .
~/land_tools/do_384land.csh

# Optional GFDL archiving: define gfdl_archive and SEND_FILE before enabling.
if ($NO_SEND == "send") then

#########################################################################
# generate date for file names
########################################################################

    set begindate = `$TIME_STAMP -bhf digital`
    if ( $begindate == "" ) set begindate = tmp`date '+%j%H%M%S'`

    set enddate = `$TIME_STAMP -ehf digital`
    if ( $enddate == "" ) set enddate = tmp`date '+%j%H%M%S'`
    set fyear = `echo $enddate | cut -c -4`

    cd $WORKDIR/rundir
    cat time_stamp.out

########################################################################
# save ascii output files
########################################################################

    if ( ! -d $WORKDIR/ascii ) mkdir $WORKDIR/ascii
    if ( ! -d $WORKDIR/ascii ) then
     echo "ERROR: $WORKDIR/ascii is not a directory."
     exit 1
    endif

    mkdir -p $WORKDIR/ascii/$begindate
    foreach out (`ls *.out *.results input*.nml *_table`)
      mv $out $WORKDIR/ascii/$begindate/
    end

    tar cvf - *\.out *\.results | gzip -c > $WORKDIR/ascii/$begindate.ascii_out.tgz

    #sbatch --export=source=$WORKDIR/ascii/$begindate.ascii_out.tgz,destination=gfdl:$gfdl_archive/ascii/$begindate.ascii_out.tgz,extension=null,type=ascii --output=$HOME/STDOUT/%x.o%j $SEND_FILE

########################################################################
# move restart files
########################################################################

    cd $WORKDIR

    if ( ! -d $WORKDIR/restart ) mkdir -p $WORKDIR/restart

    if ( ! -d $WORKDIR/restart ) then
      echo "ERROR: $WORKDIR/restart is not a directory."
      exit
    endif

    find $WORKDIR/rundir/RESTART -iname '*.res*' > $WORKDIR/rundir/file.restart.list.txt
    find $WORKDIR/rundir/RESTART -iname '*_data*' >> $WORKDIR/rundir/file.restart.list.txt
    set resfiles     = `wc -l $WORKDIR/rundir/file.restart.list.txt | awk '{print $1}'`

   if ( $resfiles > 0 ) then

      set dateDir = $WORKDIR/restart/$enddate
      set restart_file = $dateDir

      set list = `ls -C1 $WORKDIR/rundir/RESTART`
#      if ( $irun < $segmentsPerJob ) then
#        rm -r $workDir/INPUT/*.res*
#        foreach index ($list)
#          cp $workDir/RESTART/$index $workDir/INPUT/$index
#        end
#      endif

      if ( ! -d $dateDir ) mkdir -p $dateDir

      if ( ! -d $dateDir ) then
        echo "ERROR: $dateDir is not a directory."
        exit
      endif

      foreach index ($list)
        mv $WORKDIR/rundir/RESTART/$index $restart_file/$index
      end

      ln -sf $restart_file/* $WORKDIR/rundir/RESTART/
      #sbatch --export=source=$WORKDIR/restart/$enddate,destination=gfdl:$gfdl_archive/restart/$enddate,extension=tar,type=restart --output=$HOME/STDOUT/%x.o%j $SEND_FILE

   endif

########################################################################
# move history files
########################################################################

    cd $WORKDIR

    if ( ! -d $WORKDIR/history ) mkdir -p $WORKDIR/history
    if ( ! -d $WORKDIR/history ) then
      echo "ERROR: $WORKDIR/history is not a directory."
      exit 1
    endif

    set dateDir = $WORKDIR/history/$begindate
    if ( ! -d  $dateDir ) mkdir $dateDir
    if ( ! -d  $dateDir ) then
      echo "ERROR: $dateDir is not a directory."
      exit 1
    endif

    find $WORKDIR/rundir -maxdepth 1 -type f -regex '.*.nc'      -exec mv {} $dateDir \;
    find $WORKDIR/rundir -maxdepth 1 -type f -regex '.*.nc.....' -exec mv {} $dateDir \;

    cd $dateDir
      if ( ! -d ${begindate}_atmos ) mkdir -p ${begindate}_atmos
      mv *atmos.*.nc* ${begindate}_atmos
      mv ${begindate}_atmos ../.

    cd $WORKDIR/rundir

    sbatch --export=source=$WORKDIR/history/$begindate,destination=gfdl:$gfdl_archive/history/$begindate,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE

    sbatch --export=source=$WORKDIR/history/${begindate}_atmos,destination=gfdl:$gfdl_archive/history/${begindate}_atmos,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
endif
