#!/bin/tcsh
#SBATCH --output=./stdout/%x.%j
#SBATCH --job-name=Global_shiemom
#SBATCH --clusters=c6
#SBATCH --time=07:30:00
#SBATCH --nodes=106

# Script to run a global C48 SHiELD + 1deg MOM6
# Joseph.mouallem@noaa.gov

set echo

set res = 48
set ocean_mod = 'real'  #options: ideal or real

echo "Cluster: $SLURM_CLUSTER_NAME"

if (${SLURM_CLUSTER_NAME} == "c6") then
  set BASE_DIR    = "/gpfs/f6/bil-coastal-gfdl/proj-shared/Joseph.Mouallem/"           # Sim/work directory - to be set
  set BASE_DIR    = "/gpfs/f6/bil-coastal-gfdl/scratch/Joseph.Mouallem/"           # Sim/work directory - to be set
  set BUILD_DIR = "~${USER}/shiemom/SHiELD_build/"                  # Build directory - To be set
  set INPUT_DATA = "/gpfs/f6/bil-coastal-gfdl/proj-shared/gfdl_w/SHiELD_INPUT_DATA/"
endif

#if (${SLURM_CLUSTER_NAME} == "c5") then
#  set BASE_DIR    = "/gpfs/f5/gfdl_w/scratch/Joseph.Mouallem/"                     # Sim/work directory - to be set
#  set BUILD_DIR = "~${USER}/shiemom/SHiELD_build/"                  # Build directory - To be set
#  set INPUT_DATA = "/gpfs/f5/gfdl_w/proj-shared/SHiELD_INPUT_DATA/"
#endif

unset echo
source ${BUILD_DIR}/site/environment.intel.csh
set echo

set RUN_DIR = "`pwd`"                           # For the diag_table, field_table, momsis input

if ( ! $?COMPILER ) then
  set COMPILER = "intel"
endif

set RELEASE = "SHiEMOM/GLOBAL_ocean/"

# case specific details
set TYPE = "nh"          # choices:  nh, hydro
set MODE = "64bit"      # choices:  32bit, 64bit
set CASE = "C${res}"
set MONO = "non-mono"
set NAME = "20160801.00Z"
set MEMO = "$SLURM_JOB_NAME.res$res"
set MEMO = "res$res.ocean_$ocean_mod"
set PBL  = "TKE"        # choices:  TKE or YSU
set HYPT = "off"         # choices:  on, off  (controls hyperthreading)
set COMP = "debug"       # choices:  debug, repro, prod
set NO_SEND = "no_send"  # choices:  send, no_send
set EXE = "x"
set HYPT = "on"         # choices:  on, off  (controls hyperthreading)
if ( ! $?COMP ) then
  set COMP = "repro"       # choices:  debug, repro, prod
endif
set NO_SEND = "no_send"    # choices:  send, no_send
set RESTART_RUN = "F"

# directory structure
set WORKDIR    = ${BASE_DIR}/${RELEASE}/${NAME}.${CASE}.${TYPE}.${MODE}.${COMPILER}.${MONO}.${MEMO}/
set executable = ${BUILD_DIR}/Build/bin/SHiEMOM_${TYPE}.${COMP}.${MODE}.${COMPILER}.${EXE}

set ICDIR   = ${INPUT_DATA}/global.v201810/${CASE}/${NAME}_IC/
set ICS  = ${ICDIR}/GFS_INPUT.tar
set FIX  = ${INPUT_DATA}/fix.v201810/
set GFS  = ${INPUT_DATA}/GFS_STD_INPUT.20160311.tar
set GRID = ${INPUT_DATA}/global.v201810/${CASE}/GRID/
set FIX_bqx  = ${INPUT_DATA}/climo_data.v201807

# sending file to gfdl
set TIME_STAMP = ${USER}/Util/time_stamp.csh

set DIAG_TABLE  = ${RUN_DIR}/tables/diag_table_6species_tc_dp_ocean
set FIELD_TABLE = ${RUN_DIR}/tables/field_table_6species_atmland # will be changed later if tke-edmf is used

# changeable parameters
# dycore definitions
set npx = "49"
set npy = "49"
set npz = "79"
set layout_x = "12"
set layout_y = "12"
set io_layout = "1,1"
set nthreads = "4"

# blocking factor used for threading and general physics performance
set blocksize = "36"

# run length
set months = "0"
set days = "0"
set hours = "0"
set seconds = "900"
set dt_atmos = "450"
set dt_therm = "900"  # for ocean

#fms yaml
set use_yaml=".F." #if True, requires data_table.yaml and field_table.yaml

set na_init = 1
set rough = "hwrf17" # hwrf17; coare3.5; beljaars; charnock

# variables for controlling initialization of NCEP/NGGPS ICs
set filtered_terrain = ".true."
set ncep_levs = "64"
set gfs_dwinds = ".true."

# PBL related settings
switch ($PBL)
case "TKE":
    set FIELD_TABLE = ${FIELD_TABLE}_tke
    set satmedmf = ".true."
    set ysupbl = ".false."
    breaksw
case "YSU":
    set satmedmf = ".false."
    set ysupbl = ".true."
endsw

# variables for gfs diagnostic output intervals and time to zero out time-accumulated data
# set fdiag = "6.,12.,18.,24.,30.,36.,42.,48.,54.,60.,66.,72.,78.,84.,90.,96.,102.,108.,114.,120.,126.,132.,138.,144.,150.,156.,162.,168.,174.,180.,186.,192.,198.,204.,210.,216.,222.,228.,234.,240."
set fdiag = "1."
set fhzer = "1."
set fhcyc = "24."

# determines whether FV3 or GFS physics calculate geopotential
set gfs_phil = ".false."

# determine whether ozone production occurs in GFS physics
set ozcalc = ".true."

# set various debug options
set no_dycore = ".false."
set dycore_only = ".false."
set chksum_debug = ".false."
set print_freq = "6"

if (${TYPE} == "nh") then
  # non-hydrostatic options
  set make_nh = ".T."
  set hydrostatic = ".F."
  set phys_hydrostatic = ".F."     # can be tested
  set use_hydro_pressure = ".F."   # can be tested
  set consv_te = "1."
    # time step parameters in FV3
  set k_split = "2"
  set k_split = "4"
  set n_split = "6"
else
  # hydrostatic options
  set make_nh = ".F."
  set hydrostatic = ".T."
  set phys_hydrostatic = ".F."     # will be ignored in hydro mode
  set use_hydro_pressure = ".T."   # have to be .T. in hydro mode
  set consv_te = "0."
    # time step parameters in FV3
  set k_split = "2"
  set n_split = "6"
endif

if (${MONO} == "mono" || ${MONO} == "monotonic") then
  # monotonic options
  set d_con = "1."
  set do_vort_damp = ".false."
else
  # non-monotonic options
  set d_con = "1."
  set do_vort_damp = ".true."
endif

# variables for hyperthreading
if (${HYPT} == "on") then
  set hyperthread = ".true."
  set j_opt = "-j2"
  set div = 2
else
  set hyperthread = ".false."
  set j_opt = "-j1"
  set div = 1
endif

# when running with threads, need to use the following command
@ npes = ${layout_x} * ${layout_y} * 6
@ skip = ${nthreads} / ${div}
set run_cmd = "srun --ntasks=$npes --cpus-per-task=$skip ./$executable:t"

setenv SLURM_CPU_BIND verbose

setenv MPICH_ENV_DISPLAY
setenv MPICH_MPIIO_CB_ALIGN 2
setenv MALLOC_MMAP_MAX_ 0
setenv MALLOC_TRIM_THRESHOLD_ 536870912
setenv NC_BLKSZ 1M
# necessary for OpenMP when using Intel
setenv KMP_STACKSIZE 256m

if (${RESTART_RUN} == "F") then
  \rm -rf $WORKDIR/rundir
  mkdir -p $WORKDIR/rundir
  cd $WORKDIR/rundir

  mkdir -p RESTART
  mkdir -p INPUT

  # Date specific ICs
  if ( -e ${ICDIR}/gfs_data.tile1.nc ) then
    cp -rf ${ICDIR}/* INPUT/
  else
    tar xf ${ICS}
  endif

  # set variables in input.nml for initial run
  set nggps_ic = ".T."
  set mountain = ".F."
  set external_ic = ".T."
  set warm_start = ".F."

else

  cd $WORKDIR/rundir
  \rm -rf INPUT/*

  # move the restart data into INPUT/
  mv ${RESTART}/* INPUT/.

  # reset values in input.nml for restart run
  set make_nh = ".F."
  set nggps_ic = ".F."
  set mountain = ".T."
  set external_ic = ".F."
  set warm_start = ".T."
  set na_init = 0

endif

# copy over the other tables and executable
if ( ${use_yaml} == ".T." ) then
  #cp ${BUILD_AREA}/tables/data_table.yaml data_table.yaml
  cp ${BUILD_DIR}/tables/field_table_6species.yaml field_table.yaml
else
  #cp ${BUILD_AREA}/tables/data_table data_table
  cp ${BUILD_DIR}/tables/field_table_6species field_table
endif

cp ${BUILD_DIR}/tables/diag_table_no3d diag_table
cp $executable .

# GFS standard input data
tar xf ${GFS}

# Grid and orography data
cp -rf ${GRID}/* INPUT/.

# build the date for curr_date from NAME

#unset echo
set y = `echo ${NAME} | cut -c1-4`
set m = `echo ${NAME} | cut -c5-6`
set d = `echo ${NAME} | cut -c7-8`
set h = `echo ${NAME} | cut -c10-11`
set echo
set curr_date = "${y},${m},${d},${h},0,0"

# build the diag_table with the experiment name and date stamp
cat >! diag_table << EOF
${NAME}.${CASE}.${MODE}.${MONO}
$y $m $d $h 0 0 
EOF

cat ${DIAG_TABLE} >> diag_table
cp ${FIELD_TABLE} field_table
cp $FIX/global_sfc_emissivity_idx.txt INPUT/sfc_emissivity_idx.txt
cp INPUT/aerosol.dat .
cp INPUT/co2historicaldata_*.txt .
cp INPUT/sfc_emissivity_idx.txt .
cp INPUT/solarconstant_noaa_an.txt .

###########################################################################
###########################################################################
############################ ADD MOM6 input files
###########################################################################
###########################################################################

if ( ${ocean_mod} == "ideal" ) then # similar to the paper config, constant sst, u=v=0
  set MOM_INPUT_DIR = "${RUN_DIR}/MOMSIS_inputfiles/"
  cp $MOM_INPUT_DIR/MOM_* .
  cp $MOM_INPUT_DIR/SIS_* .
  ln -sf $MOM_INPUT_DIR/INPUT/layer_coord.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/hycom1_75_800m.nc INPUT/
  ln -sf /gpfs/f6/bil-coastal-gfdl/proj-shared/Joseph.Mouallem/shiemom_pdata/GLOBAL/C48_1deg/* INPUT/ # grids, mosaic, etc
endif

if ( ${ocean_mod} == "real" ) then # same as OM1_deg
  set MOM_INPUT_DIR = "/gpfs/f6/bil-coastal-gfdl/proj-shared/Joseph.Mouallem/shiemom_pdata/GLOBAL/MOM6-examples/ice_ocean_SIS2/OM_1deg/" # public github files 
  cp $MOM_INPUT_DIR/MOM_* .
  cp $MOM_INPUT_DIR/SIS_* .
  ln -sf $MOM_INPUT_DIR/INPUT/topo_edits_011818.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/MOM_channels_SPEAR INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/layer_coord.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/hycom1_75_800m.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/woa13_decav_ptemp_monthly_fulldepth_01.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/woa13_decav_s_monthly_fulldepth_01.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/vgrid_75_2m.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/KH_background_2d.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/tidal_amplitude.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/topog.nc INPUT/
  ln -sf $MOM_INPUT_DIR/INPUT/seawifs_1998-2006_smoothed_2X.nc INPUT/
  ln -sf /gpfs/f6/bil-coastal-gfdl/proj-shared/Joseph.Mouallem/shiemom_pdata/GLOBAL/C48_1deg/*mosaic* INPUT/
  ln -sf /gpfs/f6/bil-coastal-gfdl/proj-shared/Joseph.Mouallem/shiemom_pdata/GLOBAL/C48_1deg/grid_spec.nc INPUT/
  ln -sf /gpfs/f6/bil-coastal-gfdl/proj-shared/Joseph.Mouallem/shiemom_pdata/GLOBAL/C48_1deg/ocean_hgrid.nc INPUT/
endif

set input_filename = 'n' # for mom6/sis2

cat >! MOM_override <<EOF
#override DT=${dt_atmos}
#override DT_THERM=${dt_therm}
EOF

cat >! SIS_override <<EOF
EOF


cat >! input.nml <<EOF

&xgrid_nml
  !interp_method=second_order
  !make_exchange_reproduce=.true.
/

 &SIS_input_nml
       output_directory = './',
       input_filename = ${input_filename}
       restart_input_dir = 'INPUT/',
       restart_output_dir = 'RESTART/',
       parameter_filename = 'SIS_input',
                            'SIS_override' 
/

 &MOM_input_nml
       output_directory = '.',
       input_filename = ${input_filename}
       restart_input_dir = 'INPUT',
       restart_output_dir = 'RESTART',
       parameter_filename = 'MOM_input',
                             'MOM_override'
/

 &ice_albedo_nml
       t_range = 10. 
/

 &ice_model_nml
/

 &icebergs_nml
/

 &monin_obukhov_nml
       neutral = .false. ! KGao: should not use true for coupled mode
/

 &ocean_albedo_nml
       ocean_albedo_option = 2
/

 &ocean_rough_nml
       rough_scheme = $rough
/

 &sat_vapor_pres_nml
       construct_table_wrt_liq = .true.
       construct_table_wrt_liq_and_ice = .true.
/

 &surface_flux_nml
       ncar_ocean_flux = .false. ! KGao: must be false for the coupled mode
       raoult_sat_vap = .true.
       do_iter_monin_obukhov = .true. ! KGao: turn this on for hwrf17 rough scheme
       niter_monin_obukhov = 2 ! KGao: 2 should work great
/

 &fms_affinity_nml
       affinity=.F.
/


 &atmos_model_nml
     blocksize = $blocksize
     chksum_debug = $chksum_debug
     dycore_only = $dycore_only
     fdiag = $fdiag
     fullcoupler_fluxes=1
/

&diag_manager_nml
    prepend_date = .F.
    max_num_axis_sets = 100 
    max_files = 100
    max_axes = 240
/

 &fms_io_nml
       checksum_required = .false.
       max_files_r = 100,
       max_files_w = 100,
/

 &fms_nml
       clock_grain = 'ROUTINE',
       domains_stack_size = 3000000,
       print_memory_usage = .F.
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
       make_nh =  .F.
       fv_debug = .F.
       range_warn = .F.
       reset_eta = .F.
       nudge_qv = .T.
       d2_bg_k1 = 0.20
       d2_bg_k2 = 0.15
       kord_tm = -11 ! 2019: slightly better inversion structures
       kord_mt =  11
       kord_wz =  11
       kord_tr =  11
       fill = .T.     ! 2019: enabled filling negatives from remapping
       fill_gfs = .F. !2019: disabled filling negatives from GFS physics
       hydrostatic = .F.
       phys_hydrostatic = .F.
       use_hydro_pressure = .F.
       beta = 0.
       a_imp = 1.
       p_fac = 0.1
       k_split  = $k_split
       n_split  = $n_split
       nwat = 6 
       na_init =$na_init 
       d_ext = 0.0
       dnats = 2 ! 2019: improved efficiency by not advecting o3
       fv_sg_adj = 1800 ! 2019: full-domain weak 2dz damping
       n_sponge = $npz
       d2_bg = 0.
       nord =  3  
       dddmp = 0.1
       d4_bg = 0.14
       vtdm4 = 0.02
       do_vort_damp = .T.
       external_ic = $external_ic
       nggps_ic = $nggps_ic 
       hrrrv3_ic= .F.
       mountain = $mountain 
       ncep_ic = .F.
       d_con = 1.0 ! 2019: Full-strength dissipative heating
       hord_mt = 6
       hord_vt = 6
       hord_tm = 6
       hord_dp = 6
       hord_tr = -5 
       adjust_dry_mass = .F.
       consv_te = 0.
       consv_am = .F.
       dwind_2d = .F.
       print_freq = $print_freq
       warm_start = $warm_start
       no_dycore = $no_dycore

       rf_fast = .F.
       tau = 5.
       rf_cutoff = 50.e2

       delt_max = 0.002


/

 &integ_phys_nml
       do_inline_mp = .T.
       do_sat_adj = .F.
/

 &coupler_nml
       months = $months
       days  = $days
       hours = $hours
       seconds = $seconds
       dt_atmos = $dt_atmos
       dt_cpld = $dt_therm         !off with shield, on with shieldfull
       current_date =  $curr_date
       calendar = 'julian'
       atmos_nthreads = $nthreads
       use_hyper_thread = $hyperthread
       do_ocean=.T.               !off with shield, on with shieldfull
       do_flux=.T.
       do_land=.False.
       do_ice=.T.

     !  ice_npes=-1
     !  land_npes=-1
       do_debug=.true.
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
       ldiag3d        = .false.
       fhcyc          = $fhcyc
       nst_anl        = .true.
       use_ufo        = .true.
       pre_rad        = .false.
       ncld           = 5
       zhao_mic       = .false.
       pdfcld         = .true.
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
       do_deep        = .false.
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
       ysupbl         = $ysupbl
       satmedmf       = $satmedmf 
       isatmedmf      = 1
       !use_lup_only   = .T.
       rlmx           = 150.
       cs0            = 0.25 
       clam_deep      = 0.1 
       clam_shal      = 0.1 
       !c0s_deep       = 0.01
       !c0s_shal       = 0.01
       !c1_deep        = 1.
       !c1_shal        = 1.
       xkzminv        = 1.0
       xkzm_m         = 1.0
       xkzm_h         = 1.0
       cloud_gfdl     = .true.
       do_ocean       = .false.
/

 &gfdl_mp_nml
       use_rhc_revap = .true.
       rhc_revap = .8
       do_sedi_uv = .true.
       do_sedi_w = .true.
       do_sedi_heat = .false.
       rad_snow = .true.
       rad_graupel = .true.
       rad_rain = .true.
       const_vi = .F.
       const_vs = .F.
       const_vg = .F.
       const_vr = .F.
       vi_max = 0.6
       vs_max = 8.
       vg_max = 10.
       vr_max = 12.
       qi_lim = 2.
       prog_ccn = .false.
       do_qa = .true.
       tau_l2v = 180. !225.
       tau_v2l = 90.  !150.
       rthresh = 10.e-6  ! This is a key parameter for cloud water
       dw_land  = 0.16
       dw_ocean = 0.10
  !     ql_gen = 1.0e-3
       ql_mlt = 1.0e-3
       qi0_crt = 1.2e-4
       qi0_max = 2.0e-4 
       qs0_crt = 1.0e-2
       tau_i2s = 1000.
       c_psaci = 0.1
       c_pgacs = 0.001
       rh_inc = 0.30
       rh_inr = 0.30
   !    rh_ins = 0.30
       ccn_l = 270. !300.
       ccn_o = 90. !200.
       c_paut = 0.5
       c_pracw = 0.75
       z_slope_liq  = .true.
       z_slope_ice  = .true.
       fix_negative = .true.
       icloud_f = 1
       beta = 1.22
       rewflag = 1
       reiflag = 1 
       rewmin = 5.0
       rewmax = 10.0
       reimin = 10.0
       reimax = 150.0
       rermin = 10.0
       rermax = 10000.0
       resmin = 150.0
       resmax = 10000.0
       liq_ice_combine = .false.
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

# run the executable

   sleep 1
   ${run_cmd} | tee fms.out
   if ( $? != 0 ) then
	exit
   endif
    @ irun++

if ($NO_SEND == "no_send") then
  continue
endif

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

    foreach out (`ls *.out *.results *.nml *_table`)
      mv $out $begindate.$out
    end

    tar cvf - *\.out *\.results *\.nml *_table | gzip -c > $WORKDIR/ascii/$begindate.ascii_out.tgz

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

    cd $WORKDIR/rundir


 #   sbatch --export=source=$WORKDIR/history/$begindate,destination=gfdl:$gfdl_archive/history/$begindate,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE


end # while ( $irun <= $nruns )
