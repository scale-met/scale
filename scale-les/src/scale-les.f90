!-------------------------------------------------------------------------------
!> Program SCALE-LES
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scaleles
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
     PRC_MPIfinish
  use scale_const, only: &
     CONST_setup
  use scale_calendar, only: &
     CALENDAR_setup
  use scale_random, only: &
     RANDOM_setup
  use scale_time, only: &
     TIME_setup,           &
     TIME_checkstate,      &
     TIME_advance,         &
     TIME_DOATMOS_step,    &
     TIME_DOLAND_step,     &
     TIME_DOOCEAN_step,    &
     TIME_DOCPL_calc,      &
     TIME_DOATMOS_restart, &
     TIME_DOLAND_restart,  &
     TIME_DOOCEAN_restart, &
     TIME_DOCPL_restart,   &
     TIME_DOend
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
     REAL_setup
  use scale_gridtrans, only: &
     GTRANS_setup
  use scale_interpolation, only: &
     INTERP_setup
  use scale_stats, only: &
     STAT_setup
  use scale_history, only: &
     HIST_setup, &
     HIST_write
  use scale_monitor, only: &
     MONIT_setup, &
     MONIT_write, &
     MONIT_finalize
  use mod_atmos, only: &
     ATMOS_setup, &
     ATMOS_step
  use mod_atmos_vars, only: &
     ATMOS_vars_setup,         &
     ATMOS_vars_restart_read,  &
     ATMOS_vars_restart_write, &
     ATMOS_vars_restart_check, &
     ATMOS_sw_restart,         &
     ATMOS_sw_check
  use mod_atmos_vars_sf, only: &
     ATMOS_vars_sf_setup,       &
     ATMOS_vars_sf_restart_read
  use scale_atmos_hydrostatic, only: &
     ATMOS_HYDROSTATIC_setup
  use scale_atmos_thermodyn, only: &
     ATMOS_THERMODYN_setup
  use scale_atmos_saturation, only: &
     ATMOS_SATURATION_setup
  use mod_land, only: &
     LAND_setup, &
     LAND_step
  use mod_land_vars, only: &
     LAND_vars_setup,         &
     LAND_vars_restart_read,  &
     LAND_vars_restart_write, &
     LAND_sw_restart
  use mod_ocean, only: &
     OCEAN_setup, &
     OCEAN_step
  use mod_ocean_vars, only: &
     OCEAN_vars_setup,         &
     OCEAN_vars_restart_read,  &
     OCEAN_vars_restart_write, &
     OCEAN_sw_restart
  use mod_cpl, only: &
     CPL_setup, &
     CPL_calc
  use mod_cpl_vars, only: &
     CPL_vars_setup,         &
     CPL_vars_restart_read,  &
     CPL_vars_restart_write, &
     CPL_sw_restart
  use mod_user, only: &
     USER_setup, &
     USER_step
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
  call TIME_setup( setup_TimeIntegration = .true. )

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
  ! setup monitor I/O
  call MONIT_setup

  ! setup common tools
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

  ! setup restart data: atmos
  call ATMOS_vars_restart_read
  ! setup restart data: atmos_sf
  call ATMOS_vars_sf_restart_read
  ! setup restart data: land
  call LAND_vars_restart_read
  ! setup restart data: ocean
  call OCEAN_vars_restart_read
  ! setup restart data: coupler
  call CPL_vars_restart_read

  ! setup atmosphere
  call ATMOS_setup

  ! setup land
  call LAND_setup

  ! setup ocean
  call OCEAN_setup

  ! setup coupler
  call CPL_setup

  ! setup user-defined procedure
  call USER_setup

  call PROF_rapend('Initialize')

  !########## main ##########

#ifdef _FIPP_
  call fipp_start
#endif
#ifdef _PAPI_
  call PROF_PAPI_rapstart
#endif

  if( IO_L ) write(IO_FID_LOG,*)
  if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
  call PROF_rapstart('Main Loop(Total)')

  do

    ! report current time
    call TIME_checkstate

    ! user-defined procedure
    call USER_step

    ! change to next state
    if ( TIME_DOATMOS_step ) call ATMOS_step
    if ( TIME_DOLAND_step  ) call LAND_step
    if ( TIME_DOOCEAN_step ) call OCEAN_step
    if ( TIME_DOCPL_calc   ) call CPL_calc

    ! time advance
    call TIME_advance

    ! history&monitor file output
    call HIST_write
    call MONIT_write('MAIN')

    ! restart output
    if ( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_write
    if ( LAND_sw_restart  .AND. TIME_DOLAND_restart  ) call LAND_vars_restart_write
    if ( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_write
    if ( CPL_sw_restart   .AND. TIME_DOCPL_restart   ) call CPL_vars_restart_write

    if ( TIME_DOend ) exit

  enddo

  call PROF_rapend('Main Loop(Total)')
  if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
  if( IO_L ) write(IO_FID_LOG,*)

#ifdef _FIPP_
  call fipp_stop
#endif
#ifdef _PAPI_
  call PROF_PAPI_rapstop
#endif

  !########## Finalize ##########

  ! check data
  if ( ATMOS_sw_check ) call ATMOS_vars_restart_check

  call PROF_rapreport
#ifdef _PAPI_
  call PROF_PAPI_rapreport
#endif

  call FileCloseAll

  call MONIT_finalize
  ! stop MPI
  call PRC_MPIfinish

  stop
  !=============================================================================
end program scaleles
