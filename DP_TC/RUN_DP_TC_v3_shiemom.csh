#!/bin/tcsh 
#SBATCH --output=./stdout/%x.o%j
#SBATCH --job-name=TC_DP_shiemom
#SBATCH --time=09:00:00
#SBATCH --cluster=c6
#SBATCH --nodes=9

# This script is optimized for double-periodic TC configuration 
# Joseph.Mouallem.noaa.gov or Kun.Gao@noaa.gov

# Joseph: upgraded to run shiemom (shield+mom6)
# will need to external mom6 files and grid files by BuildGrid.csh
# check slides for more info:
# https://docs.google.com/presentation/d/1DZ3vk96gJ_ETbtRxSWxbdSnP7xBgagBkslnJr-g_lMM/edit#slide=id.p

set echo

echo "Cluster: $SLURM_CLUSTER_NAME"

if (${SLURM_CLUSTER_NAME} == "c6") then
set BASE_DIR    = "/gpfs/f6/bil-coastal-gfdl/scratch/Joseph.Mouallem/test/"           # Sim/work directory - to be set
set BUILD_DIR = "~${USER}/shiemom_clean/pubrelease_c6/SHiELD_build/"                  # Build directory - To be set
set INPUT_DATA = "/gpfs/f6/bil-coastal-gfdl/proj-shared/gfdl_w/SHiELD_INPUT_DATA/"
set MOSAIC_DATA = "/gpfs/f6/bil-coastal-gfdl/proj-shared/Joseph.Mouallem/shiemom_pdata/INPUT/"
endif

if (${SLURM_CLUSTER_NAME} == "c5") then
set BASE_DIR    = "/gpfs/f5/gfdl_w/scratch/Joseph.Mouallem/test/"                     # Sim/work directory - to be set
set BUILD_DIR = "~${USER}/shiemom_clean/pubrelease_c5/SHiELD_build/"                  # Build directory - To be set
set INPUT_DATA = "/gpfs/f5/gfdl_w/proj-shared/SHiELD_INPUT_DATA/"
set MOSAIC_DATA = "/gpfs/f5/gfdl_w/world-shared/Joseph.Mouallem/shiemom_pdata/INPUT/"
endif

unset echo
source ${BUILD_DIR}/site/environment.intel.csh
set echo

set RUN_DIR = "`pwd`"                           # For the diag_table, field_table, momsis input

# --- KG: optional paramters for experiments
set sst = 303
set ini_storm = "big2" # "small", "med2", "big2"
set deglat = 20 
set days = 1 
set hord_option = "hord6"
set warm_start = ".F." # if true, set restart_dir
set res = "4km" # available options: 1km or 2km or 4km
set MEMO = ${hord_option}"_"${deglat}"N_"${sst}"K_"${ini_storm}"_"${days}"d_"${res}

set RELEASE = "SHiEMOM"

# case specific details
set TYPE = "nh"        
set MODE = "64bit"
set MONO = "non-mono" 
set GRID = "DP"
set CASE = "TC"
set DATE = "20240920.00Z" 
set PBL  = "YSU"        # choices:  TKE or YSU
set HYPT = "on"         # choices:  on, off  (controls hyperthreading)
set COMP = "repro"       # choices:  debug, repro, prod
set NO_SEND = "no_send" # choices:  send, no_send
set EXE  = "intel.x"

# directory structure
set WORKDIR    = ${BASE_DIR}/${RELEASE}/${DATE}.${CASE}.${TYPE}.${MODE}.${MONO}.${MEMO}/
set executable = ${BUILD_DIR}/Build/bin/SHiEMOM_${TYPE}.${COMP}.${MODE}.${EXE}

set RUNNAME = ${DATE}.${GRID}.${MEMO}

# input filesets
set FIXDIR  = ${INPUT_DATA}/fix_am/
set CLIMO_DATA = ${INPUT_DATA}/climo_data.v201807/

# sending file to gfdl
set gfdl_archive = /archive/${USER}/${RELEASE}/${DATE}.${GRID}.${MEMO}/
set SEND_FILE = ${USER}/Util/send_file_slurm.csh
set TIME_STAMP = ${USER}/Util/time_stamp.csh

set DIAG_TABLE  = ${RUN_DIR}/tables/diag_table_6species_tc_dp_ocean
set FIELD_TABLE = ${RUN_DIR}/tables/field_table_6species_atmland # will be changed later if tke-edmf is used

# Resolution
switch ($res) #assuming domain size=10deg, need to adjust timestep, move it here XXXX
case "1km":
   set npx = "1025"
   set npy = "1025"
   set k_split = "10"
   breaksw
case "2km":
   set npx = "513"
   set npy = "513"
   set k_split = "5"
   breaksw
case "4km":
   set npx = "257"
   set npy = "257"
   set k_split = "3"
endsw

set npz = "50"
set layout_x = "30" 
set layout_y = "30" 
set io_layout = "1,1"
set nthreads = "2"

# time step parameters
set n_split = "5"
set dt_atmos = "90"

# hord related settings
switch ($hord_option)
case "hord5":
   set hord_mt = 5
   set hord_vt = 5
   set hord_tm = 5
   set hord_dp = -5
   set hord_tr = -5
   breaksw
case "hord5_hordtr8":
   set hord_mt = 5
   set hord_vt = 5
   set hord_tm = 5
   set hord_dp = -5
   set hord_tr = 8
   breaksw
case "hord6":
   set hord_mt = 6
   set hord_vt = 6
   set hord_tm = 6
   set hord_dp = 6
   set hord_tr = -5
endsw

# initial storm related settings
switch ($ini_storm)
case "small":
   set rp_TC = 100000.
   set dp_TC = 1115.
   breaksw
case "med":
   set rp_TC = 150000.
   set dp_TC = 1115.
   breaksw
case "big":
   set rp_TC = 200000.
   set dp_TC = 1115.
case "med2":
   set rp_TC = 150000.
   set dp_TC = 1250. # same vmax as small
   breaksw
case "big2":
   set rp_TC = 200000.
   set dp_TC = 1350. # same vmax as small
endsw

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

# blocking factor used for threading and general physics performance
set blocksize = "32"

# run length
set months = "0"
set days = $days
set hours = "1"
set minutes = "0"
set seconds = "0"

set print_freq = "1"

# variables for gfs diagnostic output intervals and time to zero out time-accumulated data
set fdiag = "1."
set fhzer = "1."
set fhcyc = "0."

# set various debug options
set no_dycore = ".F."
set dycore_only = ".F." 

# variables for hyperthreading
if (${HYPT} == "on") then
   set hyperthread = ".true."
   set div = 2
else
   set hyperthread = ".false."
   set div = 1
   endif
@ skip = ${nthreads} / ${div}

# when running with threads, ned to use the following command
@ npes = ${layout_x} * ${layout_y} 
set run_cmd = "srun --ntasks=$npes --cpus-per-task=$skip ./$executable:t"

setenv SLURM_CPU_BIND verbose
setenv MPICH_ENV_DISPLAY
setenv MPICH_MPIIO_CB_ALIGN 2
setenv MALLOC_MMAP_MAX_ 0
setenv MALLOC_TRIM_THRESHOLD_ 536870912
setenv NC_BLKSZ 1M
# necessary for OpenMP when using Intel
setenv KMP_STACKSIZE 512m

rm -rf $WORKDIR/rundir
mkdir -p $WORKDIR/rundir
cd $WORKDIR/rundir
mkdir -p RESTART INPUT

cp ${RUN_DIR}/MOMSIS_INPUTFILES/MOM_input .
cp ${RUN_DIR}/MOMSIS_INPUTFILES/SIS_input .
cp ${RUN_DIR}/MOMSIS_INPUTFILES/MOM_override_${res} MOM_override
cp ${RUN_DIR}/MOMSIS_INPUTFILES/SIS_override_${res} SIS_override
ln -sf ${MOSAIC_DATA}/DP_${res}/* INPUT/

set irun = 1
set nruns = 1

while ( $irun <= $nruns )

if ( $irun == 1 ) then
 set warm_start = ".F."
 set input_filename = 'n'
else
   # move the restart data into INPUT/
   if ($NO_SEND == "no_send") then
    mv RESTART/* INPUT/.
   else
    ln -s $restart_file/* INPUT/.
   endif

   # reset values in input.nml for restart run
   set warm_start = ".T."
   set input_filename = 'r'
endif

# build the date for curr_date and diag_table from DATE
unset echo
set y = `echo ${DATE} | cut -c1-4`
set m = `echo ${DATE} | cut -c5-6`
set d = `echo ${DATE} | cut -c7-8`
set h = `echo ${DATE} | cut -c10-11`
set echo
set curr_date = "${y},${m},${d},${h},0,0"

# build the diag_table with the experiment name and date stamp
cat >! diag_table << EOF
${RUNNAME}
$y $m $d $h 0 0 
EOF

cat ${DIAG_TABLE} >> diag_table

# copy over the other tables and executable
cp ${RUN_DIR}/data_table data_table
cp ${FIELD_TABLE} field_table
cp $executable .

# GFS FIX data
ln -sf $FIXDIR/ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77 INPUT/global_o3prdlos.f77
ln -sf $FIXDIR/global_h2o_pltc.f77 INPUT/global_h2oprdlos.f77
ln -sf $FIXDIR/global_solarconstant_noaa_an.txt INPUT/solarconstant_noaa_an.txt
ln -sf $FIXDIR/global_sfc_emissivity_idx.txt INPUT/sfc_emissivity_idx.txt
ln -sf $FIXDIR/global_co2historicaldata_glob.txt INPUT/co2historicaldata_glob.txt
ln -sf $FIXDIR/co2monthlycyc.txt INPUT/co2monthlycyc.txt
foreach file ( $FIXDIR/fix_co2_proj/global_co2historicaldata_????.txt )
        ln -sf $file INPUT/`echo $file:t | sed s/global_co2historicaldata/co2historicaldata/g`
end
ln -sf $FIXDIR/global_climaeropac_global.txt INPUT/aerosol.dat
foreach file ( $FIXDIR/global_volcanic_aerosols_????-????.txt )
        ln -sf $file INPUT/`echo $file:t | sed s/global_volcanic_aerosols/volcanic_aerosols/g`
end

echo "ls working dir"
ls .

cat >! input.nml <<EOF

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
       rough_scheme = 'hwrf17' ! KGao: 'hwrf17' is highly recommended; options: 'hwrf17','coare3.5','beljaars','charnock'
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

 &topography_nml
       topog_file = 'INPUT/topog.nc'
/

 &fms_affinity_nml
       affinity=.F.
/

 &xgrid_nml
 ! xgrid_clocks_on=.T.
/

 &atmos_model_nml
       dycore_only = $dycore_only
       fdiag = $fdiag
       fullcoupler_fluxes = .true.
       debug = .false.
/

 &fms_io_nml
       checksum_required = .false.
       max_files_r = 100,
       max_files_w = 100,
/

 &fms_nml
       clock_grain = 'ROUTINE',
       domains_stack_size = 10000000,
       print_memory_usage = .F.
/

 &fv_core_nml
       layout   = $layout_x,$layout_y
       io_layout= $io_layout
       npx      = $npx
       npy      = $npy
       ntiles = 1 
       npz = $npz
       !npz_type = "emc"
       grid_type = 4
       deglat = $deglat
       !dx_const = 3250. ! not used if domain_deg is set 
       !dy_const = 3250. ! not used if domain_deg is set
       domain_deg = 10. 
       make_nh = .F. 
       fv_debug = .F.
       range_warn = .T.
       reset_eta = .F.
       nudge_qv = .F.
       rf_fast = .F.
       kord_tm = -9
       kord_mt =  9
       kord_wz =  9
       kord_tr =  9
       hydrostatic = .F. 
       phys_hydrostatic = .F.
       use_hydro_pressure = .F.
       beta = 0.
       a_imp = 1.
       p_fac = 0.1
       k_split  = $k_split
       n_split  = $n_split
       nwat = 6 !with field_table_6spe 
       na_init = 0 
       d_ext = 0.0
       dnats = 1 

       ! --> Damping options

       ! Divergence damping 
       nord =  3    ! order of divergence damping
       d4_bg = 0.12 ! 0.15 was used in T-SHiELD 
       d2_bg = 0.   ! background second-order divergence damping; active even nord != 0

       ! Vorticity damping (same order as divergence damping)
       do_vort_damp = .T.
       vtdm4 = 0.03

       ! Smagorinsky
       dddmp = 0.2 ! 2nd order Smagorinsky-type divergence damping

       ! Dissipative heating
       d_con = 1.0      ! Fraction of lost KE to be converted to heat
       delt_max = 0.002 ! a limiter ? 

       ! 2dz filter
       n_sponge = 11 
       fv_sg_adj = 300 

       ! Rayleigh damping 
       tau = 1. 
       rf_cutoff = 0. ! 50e2 in T-SHiELD

       ! 2nd order diffusion for the top 2 layers
       d2_bg_k1 = 0.20 
       d2_bg_k2 = 0.12

       ! <--- end of Damping options

       ke_bg = 0.
       external_ic = .F.
       gfs_phil = .F. 
       nggps_ic = .F.
       mountain = .F.
       ncep_ic = .F.
       hord_mt = $hord_mt 
       hord_vt = $hord_vt
       hord_tm = $hord_tm
       hord_dp = $hord_dp
       hord_tr = $hord_tr 
       adjust_dry_mass = .F.
       consv_te = 0.
       consv_am = .F.
       fill = .T.
       dwind_2d = .F.
       print_freq = $print_freq
       warm_start = $warm_start 
       no_dycore = $no_dycore 
       write_3d_diags = .T.
       is_ideal_case = .true.
/

&integ_phys_nml
       do_sat_adj = .F.
       do_inline_mp = .F.
/

&test_case_nml
       test_case = 20
       rp_TC = $rp_TC
       dp_TC = $dp_TC
       Ts_TC = $sst
/

 &coupler_nml
       months = $months
       days  = $days
       hours = $hours
       minutes = $minutes
       seconds = $seconds
       dt_atmos = $dt_atmos
       do_ocean=.T.               !off with shield, on with shieldfull
       do_ice=.T.               
       do_land=.F.             
       do_atmos=.T.           
       do_flux=.T.           
       dt_cpld = $dt_atmos        !off with shield, on with shieldfull
       current_date =  $curr_date
       calendar = 'NOLEAP'
      ! do_chksum= .T.
      ! do_debug= .T.
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
       IAER           = 000 
       ico2           = 0  
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
       ysupbl         = $ysupbl
       satmedmf       = $satmedmf 
       isatmedmf      = 1
       rlmx           = 150.
       cs0            = 0.25 
       xkzminv        = 1.0
       xkzm_m         = 1.0
       xkzm_h         = 1.0
       cloud_gfdl     = .true.
       do_ocean       = .false.
       fixed_date     = .true.
       fixed_solhr    = .true.
       fixed_sollat   = .true.
       sollat         = 25.
       daily_mean     = .true.
       Ts0            = $sst
/

 &gfdl_mp_nml
       delay_cond_evap = .true.
       do_cond_timescale = .true.
       tau_v2l = 150. !600.  ! not activated in T-SHiELD 
       use_rhc_revap = .true.
       !rhc_revap = .8
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
       tau_l2v = 180. 
       !tau_v2l = 150.  ! T-SHiELD uses 90.; TC20 used same
       rthresh = 10.e-6  ! This is a key parameter for cloud water
       dw_land  = 0.16
       dw_ocean = 0.10
       ql_mlt = 1.0e-3
       qi0_crt = 1.2e-4
       qi0_max = 2.0e-4
       qs0_crt = 1.0e-2
       tau_i2s = 1000.
       c_psaci = 0.1
       c_pgacs = 0.001
       rh_inc = 0.30
       rh_inr = 0.30
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

 &diag_manager_nml 
       prepend_date = .F.
       max_axes=100 !needed when adding ocean diag entries in diag_table
/

  &interpolator_nml
       interp_method = 'conserve_great_circle'
/

  &namsfc
/

EOF

# run the executable
sleep 1
${run_cmd} | tee fms.out

# if error occurs
if ( $? != 0 ) then 
   exit
endif

@ irun++
end

## handle data transfer below
#
#if ($NO_SEND == "send") then
#
##########################################################################
## generate date for file names
#########################################################################
#
#    set begindate = `$TIME_STAMP -bhf digital`
#    if ( $begindate == "" ) set begindate = tmp`date '+%j%H%M%S'`
#
#    set enddate = `$TIME_STAMP -ehf digital`
#    if ( $enddate == "" ) set enddate = tmp`date '+%j%H%M%S'`
#    set fyear = `echo $enddate | cut -c -4`
#
#    cd $WORKDIR/rundir
#    cat time_stamp.out
#
#########################################################################
## save ascii output files
#########################################################################
#
#    if ( ! -d $WORKDIR/ascii ) mkdir $WORKDIR/ascii
#    if ( ! -d $WORKDIR/ascii ) then
#     echo "ERROR: $WORKDIR/ascii is not a directory."
#     exit 1
#    endif
#
#    mkdir -p $WORKDIR/ascii/$begindate
#    foreach out (`ls *.out *.results input*.nml *_table`)
#      cp $out $WORKDIR/ascii/$begindate/
#    end
#
#    cd $WORKDIR/ascii/$begindate
#    tar cvf - *\.out *\.results | gzip -c > $WORKDIR/ascii/$begindate.ascii_out.tgz
#
#    sbatch --export=source=$WORKDIR/ascii/$begindate.ascii_out.tgz,destination=gfdl:$gfdl_archive/ascii/$begindate.ascii_out.tgz,extension=null,type=ascii --output=$HOME/STDOUT/%x.o%j $SEND_FILE
#
#########################################################################
## move restart files
#########################################################################
#
#    cd $WORKDIR
#
#    if ( ! -d $WORKDIR/restart ) mkdir -p $WORKDIR/restart
#
#    if ( ! -d $WORKDIR/restart ) then
#      echo "ERROR: $WORKDIR/restart is not a directory."
#      exit
#    endif
#
#    find $WORKDIR/rundir/RESTART -iname '*.res*' > $WORKDIR/rundir/file.restart.list.txt
#    find $WORKDIR/rundir/RESTART -iname '*_data*' >> $WORKDIR/rundir/file.restart.list.txt
#    set resfiles     = `wc -l $WORKDIR/rundir/file.restart.list.txt | awk '{print $1}'`
#
#   if ( $resfiles > 0 ) then
#
#      set dateDir = $WORKDIR/restart/$enddate
#      set restart_file = $dateDir
#
#      set list = `ls -C1 $WORKDIR/rundir/RESTART`
#
#      if ( ! -d $dateDir ) mkdir -p $dateDir
#
#      foreach index ($list)
#        mv $WORKDIR/rundir/RESTART/$index $restart_file/$index
#      end
#
#      ln -sf $restart_file/* $WORKDIR/rundir/RESTART/
#
#      sbatch --export=source=$WORKDIR/restart/$enddate,destination=gfdl:$gfdl_archive/restart/$enddate,extension=tar,type=restart --output=$HOME/STDOUT/%x.o%j $SEND_FILE
#
#   endif
#
#########################################################################
## move history files
#########################################################################
#
#    cd $WORKDIR
#
#    if ( ! -d $WORKDIR/history ) mkdir -p $WORKDIR/history
#    if ( ! -d $WORKDIR/history ) then
#      echo "ERROR: $WORKDIR/history is not a directory."
#      exit 1
#    endif
#
#    set dateDir = $WORKDIR/history/$begindate
#    if ( ! -d  $dateDir ) mkdir $dateDir
#    if ( ! -d  $dateDir ) then
#      echo "ERROR: $dateDir is not a directory."
#      exit 1
#    endif
#
#    find $WORKDIR/rundir -maxdepth 1 -type f -regex '.*.nc'      -exec mv {} $dateDir \;
#    find $WORKDIR/rundir -maxdepth 1 -type f -regex '.*.nc.....' -exec mv {} $dateDir \;
#
#    cd $dateDir
#    #  if ( ! -d ${begindate}_nggps3d ) mkdir -p ${begindate}_nggps3d
#    #  mv nggps3d*.nc* ${begindate}_nggps3d
#    #  mv ${begindate}_nggps3d ../.
#    #  if ( ! -d ${begindate}_tracer3d ) mkdir -p ${begindate}_tracer3d
#    #  mv tracer3d*.nc* ${begindate}_tracer3d
#    #  mv ${begindate}_tracer3d ../.
#    #  if ( ! -d ${begindate}_gfs_physics ) mkdir -p ${begindate}_gfs_physics
#    #  mv gfs_physics*.nc* ${begindate}_gfs_physics
#    #  mv ${begindate}_gfs_physics ../.
#    #  if ( ! -d ${begindate}_cloud3d ) mkdir -p ${begindate}_cloud3d
#    #   mv cloud3d*.nc* ${begindate}_cloud3d
#    #  mv ${begindate}_cloud3d ../.
#
#    cd $WORKDIR/rundir
#
#    sbatch --export=source=$WORKDIR/history/$begindate,destination=gfdl:$gfdl_archive/history/$begindate,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
#    #sbatch --export=source=$WORKDIR/history/${begindate}_nggps3d,destination=gfdl:$gfdl_archive/history/${begindate}_nggps3d,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
#    #sbatch --export=source=$WORKDIR/history/${begindate}_tracer3d,destination=gfdl:$gfdl_archive/history/${begindate}_tracer3d,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
#    #sbatch --export=source=$WORKDIR/history/${begindate}_gfs_physics,destination=gfdl:$gfdl_archive/history/${begindate}_gfs_physics,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
#    #sbatch --export=source=$WORKDIR/history/${begindate}_cloud3d,destination=gfdl:$gfdl_archive/history/${begindate}_cloud3d,extension=tar,type=history --output=$HOME/STDOUT/%x.o%j $SEND_FILE
#
#else # "NO_SEND"
#
#    cd $WORKDIR/rundir
#
#    set begindate = `$TIME_STAMP -bhf digital`
#    if ( $begindate == "" ) set begindate = tmp`date '+%j%H%M%S'`
#    set enddate = `$TIME_STAMP -ehf digital`
#    if ( $enddate == "" ) set enddate = tmp`date '+%j%H%M%S'`
#
#    mkdir -p $WORKDIR/ascii/$begindate
#    mv *.out *.results *.nml *_table $WORKDIR/ascii/$begindate
#
#    mkdir -p $WORKDIR/history/$begindate
#    mv *.nc $WORKDIR/history/$begindate
#
#    mkdir -p $WORKDIR/restart/$enddate
#    mv RESTART/* $WORKDIR/restart/$enddate
#    ln -sf $WORKDIR/restart/$enddate/* RESTART/
#
#endif
