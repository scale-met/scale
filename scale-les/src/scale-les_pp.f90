!-------------------------------------------------------------------------------
!> Program SCALE-LES pp
!!
!! @par Description
!!          This program is driver of preprocess tools
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scaleles_pp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use dc_log, only: &
     LogInit
  use gtool_file, only: &
     FileCloseAll
  use scale_precision
  use scale_stdio
  use scale_prof

  use scale_process, only: &
     PRC_setup,    &
     PRC_MPIstart, &
     PRC_MPIstop,  &
     PRC_MPIfinish
  use scale_const, only: &
     CONST_setup
  use scale_calendar, only: &
     CALENDAR_setup
  use scale_random, only: &
     RANDOM_setup
  use scale_time, only: &
     TIME_setup
  use scale_grid_index, only: &
     GRID_INDEX_setup
  use scale_land_grid_index, only: &
     LAND_GRID_INDEX_setup
  use scale_grid, only: &
     GRID_setup
  use scale_land_grid, only: &
     LAND_GRID_setup
  use scale_tracer, only: &
     TRACER_setup
  use scale_fileio, only: &
     FILEIO_setup
  use scale_comm, only: &
     COMM_setup
  use scale_topography, only: &
     TOPO_setup
  use scale_landuse, only: &
     LANDUSE_setup
  use scale_grid_real, only: &
     REAL_setup,   &
     REAL_update_Z
  use scale_gridtrans, only: &
     GTRANS_setup
  use scale_interpolation, only: &
     INTERP_setup
  use scale_stats, only: &
     STAT_setup
  use scale_history, only: &
     HIST_setup
  use mod_atmos_vars, only: &
     ATMOS_vars_setup
  use mod_atmos_vars_sf, only: &
     ATMOS_vars_sf_setup
  use scale_atmos_hydrostatic, only: &
     ATMOS_HYDROSTATIC_setup
  use scale_atmos_thermodyn, only: &
     ATMOS_THERMODYN_setup
  use scale_atmos_saturation, only: &
     ATMOS_SATURATION_setup
  use mod_land_vars, only: &
     LAND_vars_setup
  use mod_ocean_vars, only: &
     OCEAN_vars_setup
  use mod_cpl_vars, only: &
     CPL_vars_setup
  use mod_cnvboundary, only: &
     CONVERT_setup, &
     CONVERT
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  !########## Initial setup ##########

  ! setup standard I/O
  call IO_setup

  ! start MPI
  call PRC_MPIstart

  ! setup process
  call PRC_setup

  ! setup Log
  call LogInit(IO_FID_CONF, IO_FID_LOG, IO_L)

  ! setup constants
  call CONST_setup

  ! setup calendar
  call CALENDAR_setup

  ! setup random number
  call RANDOM_setup

  ! setup time
  call TIME_setup( setup_TimeIntegration = .false. )

  call PROF_rapstart('Initialize')

  ! setup horisontal/vertical grid index
  call GRID_INDEX_setup
  call LAND_GRID_INDEX_setup

  ! setup grid coordinates (cartesian,idealized)
  call GRID_setup
  call LAND_GRID_setup

  ! setup tracer index
  call TRACER_setup

  ! setup file I/O
  call FILEIO_setup

  ! setup mpi communication
  call COMM_setup

  ! setup topography
  call TOPO_setup
  ! setup land use category index/fraction
  call LANDUSE_setup
  ! setup grid coordinates (real world)
  call REAL_setup

  ! setup grid transfer metrics (uses in ATMOS_dynamics)
  call GTRANS_setup
  ! setup Z-ZS interpolation factor (uses in History)
  call INTERP_setup

  ! setup statistics
  call STAT_setup
  ! setup history I/O
  call HIST_setup



  ! setup atmos
  call ATMOS_HYDROSTATIC_setup
  call ATMOS_THERMODYN_setup
  call ATMOS_SATURATION_setup

  ! setup variable container: atmos
  call ATMOS_vars_setup
  ! setup variable container: atmos_sf
  call ATMOS_vars_sf_setup
  ! setup variable container: land
  call LAND_vars_setup
  ! setup variable container: ocean
  call OCEAN_vars_setup
  ! setup variable container: coupler
  call CPL_vars_setup

  ! setup preprocess converter
!  call CONVERT_setup

  call PROF_rapend('Initialize')

  !########## main ##########

  call PROF_rapstart('Main')

  ! execute mktopo
  call PROF_rapstart('Convert')
!  call CONVERT
  call PROF_rapend  ('Convert')

  call PROF_rapend('Main')

  !########## Finalize ##########

  call PROF_rapreport

  call FileCloseAll

  ! stop MPI
  call PRC_MPIfinish

  stop
  !=============================================================================
end program scaleles_pp
