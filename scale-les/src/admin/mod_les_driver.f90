!-------------------------------------------------------------------------------
!> module SCALE-LES (a main routine of les model)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for LES-scale weather
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-12-12 (R.Yoshida)  [mod] from program of scaleles
!!
!<
!-------------------------------------------------------------------------------
module mod_les_driver
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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scale-les.h"
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: scaleles

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_MID), private, parameter :: MODELNAME = "SCALE-LES ver. "//VERSION

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine scaleles( &
       comm_world,       &
       intercomm_parent, &
       intercomm_child,  &
       cnf_fname         )
    use scale_process, only: &
       PRC_LOCAL_setup
    use scale_les_process, only: &
       PRC_setup
    use scale_const, only: &
       CONST_setup
    use scale_calendar, only: &
       CALENDAR_setup
    use scale_random, only: &
       RANDOM_setup
    use scale_grid_index, only: &
       GRID_INDEX_setup
    use scale_grid, only: &
       GRID_setup
    use scale_grid_nest, only: &
       NEST_setup
    use scale_land_grid_index, only: &
       LAND_GRID_INDEX_setup
    use scale_land_grid, only: &
       LAND_GRID_setup
    use scale_urban_grid_index, only: &
       URBAN_GRID_INDEX_setup
    use scale_urban_grid, only: &
       URBAN_GRID_setup
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
    use scale_les_statistics, only: &
       STAT_setup
    use scale_history, only: &
       HIST_setup, &
       HIST_write
    use scale_monitor, only: &
       MONIT_setup, &
       MONIT_write, &
       MONIT_finalize
    use scale_atmos_hydrostatic, only: &
       ATMOS_HYDROSTATIC_setup
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_setup
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_setup
    use scale_bulkflux, only: &
       BULKFLUX_setup
    use scale_roughness, only: &
       ROUGHNESS_setup

    use mod_admin_restart, only: &
       ADMIN_restart_setup
    use mod_admin_time, only: &
       ADMIN_TIME_setup,      &
       ADMIN_TIME_checkstate, &
       ADMIN_TIME_advance,    &
       TIME_DOATMOS_step,     &
       TIME_DOLAND_step,      &
       TIME_DOURBAN_step,     &
       TIME_DOOCEAN_step,     &
       TIME_DOATMOS_restart,  &
       TIME_DOLAND_restart,   &
       TIME_DOURBAN_restart,  &
       TIME_DOOCEAN_restart,  &
       TIME_DOresume,         &
       TIME_DOend
    use mod_atmos_admin, only: &
       ATMOS_admin_setup, &
       ATMOS_do
    use mod_atmos_vars, only: &
       ATMOS_vars_setup,                         &
       ATMOS_sw_restart => ATMOS_RESTART_OUTPUT, &
       ATMOS_vars_restart_write,                 &
       ATMOS_sw_check => ATMOS_RESTART_CHECK,    &
       ATMOS_vars_restart_check
    use mod_atmos_driver, only: &
       ATMOS_driver_setup,    &
       ATMOS_driver,           &
       ATMOS_driver_finalize
    use mod_ocean_admin, only: &
       OCEAN_admin_setup, &
       OCEAN_do
    use mod_ocean_vars, only: &
       OCEAN_vars_setup,                         &
       OCEAN_sw_restart => OCEAN_RESTART_OUTPUT, &
       OCEAN_vars_restart_write
    use mod_ocean_driver, only: &
       OCEAN_driver_setup, &
       OCEAN_driver
    use mod_land_admin, only: &
       LAND_admin_setup, &
       LAND_do
    use mod_land_vars, only: &
       LAND_vars_setup,                        &
       LAND_sw_restart => LAND_RESTART_OUTPUT, &
       LAND_vars_restart_write
    use mod_land_driver, only: &
       LAND_driver_setup, &
       LAND_driver
    use mod_urban_admin, only: &
       URBAN_admin_setup, &
       URBAN_do
    use mod_urban_vars, only: &
       URBAN_vars_setup,                         &
       URBAN_sw_restart => URBAN_RESTART_OUTPUT, &
       URBAN_vars_restart_write
    use mod_urban_driver, only: &
       URBAN_driver_setup, &
       URBAN_driver
    use mod_cpl_admin, only: &
       CPL_admin_setup
    use mod_cpl_vars, only: &
       CPL_vars_setup
    use mod_user, only: &
       USER_setup, &
       USER_step
    implicit none

    integer,               intent(in) :: comm_world
    integer,               intent(in) :: intercomm_parent
    integer,               intent(in) :: intercomm_child
    character(len=H_LONG), intent(in) :: cnf_fname

    integer :: myrank
    logical :: ismaster
    !---------------------------------------------------------------------------

    !########## Initial setup ##########

    ! setup standard I/O
    call IO_setup( MODELNAME, .true., cnf_fname )

    ! setup MPI
    call PRC_LOCAL_setup( comm_world, & ! [IN]
                          myrank,     & ! [OUT]
                          ismaster    ) ! [OUT]

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    call LogInit( IO_FID_CONF, IO_FID_LOG, IO_L )

    ! setup process
    call PRC_setup

    ! setup PROF
    call PROF_setup

    ! setup constants
    call CONST_setup

    ! setup calendar
    call CALENDAR_setup

    ! setup random number
    call RANDOM_setup

    ! setup time
    call ADMIN_TIME_setup( setup_TimeIntegration = .true. )

    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)

    ! setup horizontal/vertical grid coordinates (cartesian,idealized)
    call GRID_INDEX_setup
    call GRID_setup

    call LAND_GRID_INDEX_setup
    call LAND_GRID_setup

    call URBAN_GRID_INDEX_setup
    call URBAN_GRID_setup

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

    ! setup restart
    call ADMIN_restart_setup
    ! setup statistics
    call STAT_setup
    ! setup history I/O
    call HIST_setup
    ! setup monitor I/O
    call MONIT_setup

    ! setup nesting grid
    call NEST_setup ( intercomm_parent, intercomm_child )

    ! setup common tools
    call ATMOS_HYDROSTATIC_setup
    call ATMOS_THERMODYN_setup
    call ATMOS_SATURATION_setup

    call BULKFLUX_setup
    call ROUGHNESS_setup

    ! setup submodel administrator
    call ATMOS_admin_setup
    call OCEAN_admin_setup
    call LAND_admin_setup
    call URBAN_admin_setup
    call CPL_admin_setup

    ! setup variable container
    call ATMOS_vars_setup
    call OCEAN_vars_setup
    call LAND_vars_setup
    call URBAN_vars_setup
    call CPL_vars_setup

    ! setup submodel driver
    call ATMOS_driver_setup
    call OCEAN_driver_setup
    call LAND_driver_setup
    call URBAN_driver_setup

    call USER_setup

    call PROF_rapend('Initialize', 0)

    !########## main ##########

#ifdef _FIPP_
    call fipp_start
#endif
#ifdef _PAPI_
    call PROF_PAPI_rapstart
#endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START TIMESTEP ++++++'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)

    do

      ! report current time
      call ADMIN_TIME_checkstate

      if ( TIME_DORESUME ) then
         ! resume state from restart files
         call resume_state

         ! history&monitor file output
         call HIST_write ! if needed
         call MONIT_write('MAIN')
      end if


      ! time advance
      call ADMIN_TIME_advance

      ! user-defined procedure
      call USER_step

      ! change to next state
      if( OCEAN_do .AND. TIME_DOOCEAN_step ) call OCEAN_driver
      if( LAND_do  .AND. TIME_DOLAND_step  ) call LAND_driver
      if( URBAN_do .AND. TIME_DOURBAN_step ) call URBAN_driver
      if( ATMOS_do .AND. TIME_DOATMOS_step ) call ATMOS_driver

      ! history&monitor file output
      call HIST_write
      call MONIT_write('MAIN')

      ! restart output
      if( OCEAN_sw_restart .AND. TIME_DOOCEAN_restart ) call OCEAN_vars_restart_write
      if( LAND_sw_restart  .AND. TIME_DOLAND_restart  ) call LAND_vars_restart_write
      if( URBAN_sw_restart .AND. TIME_DOURBAN_restart ) call URBAN_vars_restart_write
      if( ATMOS_sw_restart .AND. TIME_DOATMOS_restart ) call ATMOS_vars_restart_write

      if( TIME_DOend ) exit

      if( IO_L ) call flush(IO_FID_LOG)

    enddo

    if( ATMOS_do ) call ATMOS_driver_finalize

    call PROF_rapend('Main_Loop', 0)
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
    if( ATMOS_sw_check ) call ATMOS_vars_restart_check

    call PROF_rapreport
#ifdef _PAPI_
    call PROF_PAPI_rapreport
#endif

    call MONIT_finalize

    call FileCloseAll

    return
  end subroutine scaleles


  subroutine resume_state
    use mod_atmos_driver, only: &
       ATMOS_driver_resume1,    &
       ATMOS_driver_resume2,    &
       ATMOS_SURFACE_SET
    use mod_ocean_driver, only: &
       OCEAN_driver_resume, &
       OCEAN_SURFACE_SET
    use mod_land_driver, only: &
       LAND_driver_resume, &
       LAND_SURFACE_SET
    use mod_urban_driver, only: &
       URBAN_driver_resume, &
       URBAN_SURFACE_SET
    use mod_atmos_vars, only: &
       ATMOS_vars_diagnostics, &
       ATMOS_vars_monitor,     &
       ATMOS_vars_restart_read
    use mod_ocean_vars, only: &
       OCEAN_vars_restart_read
    use mod_land_vars, only: &
       LAND_vars_restart_read
    use mod_urban_vars, only: &
       URBAN_vars_restart_read
    use mod_user, only: &
       USER_resume0, &
       USER_resume
    use mod_atmos_admin, only: &
       ATMOS_do
    use mod_ocean_admin, only: &
       OCEAN_do
    use mod_land_admin, only: &
       LAND_do
    use mod_urban_admin, only: &
       URBAN_do
    implicit none

    ! read restart data
    if( ATMOS_do ) call ATMOS_vars_restart_read
    if( OCEAN_do ) call OCEAN_vars_restart_read
    if( LAND_do  ) call LAND_vars_restart_read
    if( URBAN_do ) call URBAN_vars_restart_read

    ! setup user-defined procedure before setup of other components
    call USER_resume0

    ! calc diagnostics
    call ATMOS_vars_diagnostics

    ! first monitor
    call ATMOS_vars_monitor

    ! setup surface condition
    call ATMOS_SURFACE_SET( countup=.false. )
    call OCEAN_SURFACE_SET( countup=.false. )
    call LAND_SURFACE_SET ( countup=.false. )
    call URBAN_SURFACE_SET( countup=.false. )

    ! setup submodel driver
    call ATMOS_driver_resume1
    call OCEAN_driver_resume
    call LAND_driver_resume
    call URBAN_driver_resume
    call ATMOS_driver_resume2

    ! setup user-defined procedure
    call USER_resume

    return
  end subroutine resume_state

end module mod_les_driver
