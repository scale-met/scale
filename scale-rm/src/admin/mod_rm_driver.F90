!-------------------------------------------------------------------------------
!> module SCALE-RM (a main routine of regional model)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Limited area model for regional weather, regional climate, and large-Eddy Simulation (LES)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_rm_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scale-rm.h"
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: rm_driver

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
  character(len=H_MID), private, parameter :: MODELNAME = "SCALE-RM ver. "//VERSION

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine rm_driver( &
       comm_world,       &
       intercomm_parent, &
       intercomm_child,  &
       cnf_fname         )
    use scale_file, only: &
       FILE_Close_All
    use scale_prc, only: &
       PRC_abort,      &
       PRC_LOCAL_setup
    use scale_fpm, only: &
       FPM_alive,       &
       FPM_Polling,     &
       FPM_POLLING_FREQ
    use scale_prc_cartesC, only: &
       PRC_CARTESC_setup
    use scale_const, only: &
       CONST_setup
    use scale_calendar, only: &
       CALENDAR_setup
    use scale_random, only: &
       RANDOM_setup
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_setup
    use scale_atmos_grid_cartesC_index, only: &
       ATMOS_GRID_CARTESC_INDEX_setup, &
       IA, JA
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_setup, &
       DOMAIN_CENTER_Y => ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
       CY => ATMOS_GRID_CARTESC_CY, &
       DX, DY
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_setup,        &
       ATMOS_GRID_CARTESC_REAL_calc_areavol, &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT
    use scale_atmos_grid_cartesC_metric, only: &
       ATMOS_GRID_CARTESC_METRIC_setup, &
       ATMOS_GRID_CARTESC_METRIC_MAPF
    use scale_ocean_grid_cartesC_index, only: &
       OCEAN_GRID_CARTESC_INDEX_setup
    use scale_ocean_grid_cartesC, only: &
       OCEAN_GRID_CARTESC_setup
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_setup, &
       OCEAN_GRID_CARTESC_REAL_set_areavol
    use scale_land_grid_cartesC_index, only: &
       LAND_GRID_CARTESC_INDEX_setup
    use scale_land_grid_cartesC, only: &
       LAND_GRID_CARTESC_setup
    use scale_land_grid_cartesC_real, only: &
       LAND_GRID_CARTESC_REAL_setup, &
       LAND_GRID_CARTESC_REAL_set_areavol
    use scale_urban_grid_cartesC_index, only: &
       URBAN_GRID_CARTESC_INDEX_setup
    use scale_urban_grid_cartesC, only: &
       URBAN_GRID_CARTESC_setup
    use scale_urban_grid_cartesC_real, only: &
       URBAN_GRID_CARTESC_REAL_setup, &
       URBAN_GRID_CARTESC_REAL_set_areavol
    use scale_file_cartesC, only: &
       FILE_CARTESC_setup, &
       FILE_CARTESC_cleanup
    use scale_comm_cartesC, only: &
       COMM_setup , &
       COMM_cleanup
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_setup
    use scale_topography, only: &
       TOPOGRAPHY_setup
    use scale_landuse, only: &
       LANDUSE_setup
    use scale_statistics, only: &
       STATISTICS_setup
    use scale_time, only: &
       TIME_NOWDATE, &
       TIME_NOWMS,   &
       TIME_NOWSTEP, &
       TIME_DTSEC
    use scale_coriolis, only: &
       CORIOLIS_setup
    use scale_atmos_hydrostatic, only: &
       ATMOS_HYDROSTATIC_setup
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_setup
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_setup
    use scale_bulkflux, only: &
       BULKFLUX_setup
    use scale_file_external_input_cartesC, only: &
       FILE_EXTERNAL_INPUT_CARTESC_setup
    use scale_file_history, only: &
       FILE_HISTORY_write, &
       FILE_HISTORY_set_nowdate, &
       FILE_HISTORY_finalize
    use scale_file_history_cartesC, only: &
       FILE_HISTORY_CARTESC_setup
    use scale_monitor_cartesC, only: &
       MONITOR_CARTESC_setup
    use scale_monitor, only: &
       MONITOR_write, &
       MONITOR_finalize

    use mod_atmos_driver, only: &
       ATMOS_driver_tracer_setup
    use mod_admin_versioncheck, only: &
       ADMIN_versioncheck
    use mod_admin_time, only: &
       ADMIN_TIME_setup,      &
       ADMIN_TIME_checkstate, &
       ADMIN_TIME_advance,    &
       TIME_DOATMOS_step,     &
       TIME_DOLAND_step,      &
       TIME_DOURBAN_step,     &
       TIME_DOOCEAN_step,     &
       TIME_DOresume,         &
       TIME_DOend
    use mod_admin_restart, only: &
       ADMIN_restart_setup, &
       ADMIN_restart_write
    use mod_atmos_admin, only: &
       ATMOS_admin_setup, &
       ATMOS_do,          &
       ATMOS_PHY_MP_TYPE
    use mod_atmos_vars, only: &
       ATMOS_vars_setup,           &
       ATMOS_RESTART_CHECK,        &
       ATMOS_vars_restart_check,   &
       ATMOS_vars_history_setpres, &
       ATMOS_vars_history,         &
       ATMOS_vars_monitor
    use mod_atmos_driver, only: &
       ATMOS_driver_setup,                    &
       ATMOS_driver_calc_tendency,            &
       ATMOS_driver_calc_tendency_from_sflux, &
       ATMOS_driver_update,                   &
       ATMOS_driver_finalize
    use mod_atmos_phy_mp_vars, only: &
       QA_MP
    use mod_ocean_admin, only: &
       OCEAN_admin_setup, &
       OCEAN_do
    use mod_ocean_vars, only: &
       OCEAN_vars_setup,  &
       OCEAN_vars_history
    use mod_ocean_driver, only: &
       OCEAN_driver_setup,         &
       OCEAN_driver_calc_tendency, &
       OCEAN_driver_update
    use mod_land_admin, only: &
       LAND_admin_setup, &
       LAND_do
    use mod_land_vars, only: &
       LAND_vars_setup,  &
       LAND_vars_history
    use mod_land_driver, only: &
       LAND_driver_setup,         &
       LAND_driver_calc_tendency, &
       LAND_driver_update
    use mod_urban_admin, only: &
       URBAN_admin_setup, &
       URBAN_do,          &
       URBAN_land
    use mod_urban_vars, only: &
       URBAN_vars_setup,  &
       URBAN_vars_history
    use mod_urban_driver, only: &
       URBAN_driver_setup,         &
       URBAN_driver_calc_tendency, &
       URBAN_driver_update
    use mod_lake_admin, only: &
       LAKE_admin_setup, &
       LAKE_do
    use mod_cpl_admin, only: &
       CPL_admin_setup, &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_vars_setup
    use mod_user, only: &
       USER_tracer_setup,  &
       USER_setup,         &
       USER_calc_tendency, &
       USER_update
    implicit none

    integer,          intent(in) :: comm_world
    integer,          intent(in) :: intercomm_parent
    integer,          intent(in) :: intercomm_child
    character(len=*), intent(in) :: cnf_fname

    integer :: myrank
    integer :: fpm_counter
    logical :: ismaster
    logical :: sign_exit
    !---------------------------------------------------------------------------

    !########## Initial setup ##########

    ! setup standard I/O
    call IO_setup( MODELNAME, cnf_fname )

    ! setup MPI
    call PRC_LOCAL_setup( comm_world, & ! [IN]
                          myrank,     & ! [OUT]
                          ismaster    ) ! [OUT]

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )

    ! namelist compatibility check
    call ADMIN_versioncheck

    ! setup process
    call PRC_CARTESC_setup

    ! setup PROF
    call PROF_setup

    !###########################################################################
    ! profiler start
    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)

    ! setup constants
    call CONST_setup

    ! setup calendar
    call CALENDAR_setup

    ! setup random number
    call RANDOM_setup

    ! setup submodel administrator
    call ATMOS_admin_setup
    call OCEAN_admin_setup
    call LAND_admin_setup
    call URBAN_admin_setup
    call LAKE_admin_setup
    call CPL_admin_setup

    ! setup horizontal/vertical grid coordinates (cartesian,idealized)
    if ( ATMOS_do ) then
       call ATMOS_GRID_CARTESC_INDEX_setup
       call ATMOS_GRID_CARTESC_setup
    endif

    if ( OCEAN_do ) then
       call OCEAN_GRID_CARTESC_INDEX_setup
       call OCEAN_GRID_CARTESC_setup
    endif

    if ( LAND_do ) then
       call LAND_GRID_CARTESC_INDEX_setup
       call LAND_GRID_CARTESC_setup
    endif

    if ( URBAN_do ) then
       call URBAN_GRID_CARTESC_INDEX_setup
       call URBAN_GRID_CARTESC_setup
    endif

    ! setup tracer index
    call ATMOS_HYDROMETEOR_setup
    call ATMOS_driver_tracer_setup
    call USER_tracer_setup

    ! setup file I/O
    call FILE_CARTESC_setup

    ! setup mpi communication
    call COMM_setup

    ! setup topography
    call TOPOGRAPHY_setup
    ! setup land use category index/fraction
    call LANDUSE_setup( OCEAN_do, (.not. URBAN_land), LAKE_do )

    ! setup grid coordinates (real world)
    if ( ATMOS_do ) then
       call ATMOS_GRID_CARTESC_REAL_setup
       ! setup grid transfer metrics (uses in ATMOS_dynamics)
       call ATMOS_GRID_CARTESC_METRIC_setup
       call ATMOS_GRID_CARTESC_REAL_calc_areavol( ATMOS_GRID_CARTESC_METRIC_MAPF(:,:,:,:) )
    endif
    if ( OCEAN_do ) then
       call OCEAN_GRID_CARTESC_REAL_setup
       call OCEAN_GRID_CARTESC_REAL_set_areavol
    end if
    if ( LAND_do  ) then
       call LAND_GRID_CARTESC_REAL_setup
       call LAND_GRID_CARTESC_REAL_set_areavol
    end if
    if ( URBAN_do ) then
       call URBAN_GRID_CARTESC_REAL_setup
       call URBAN_GRID_CARTESC_REAL_set_areavol
    end if

    ! setup restart
    call ADMIN_restart_setup
    ! setup time
    call ADMIN_TIME_setup( setup_TimeIntegration = .true. )
    ! setup statistics
    call STATISTICS_setup
    ! setup history I/O
    call FILE_HISTORY_CARTESC_setup
    ! setup monitor I/O
    call MONITOR_CARTESC_setup( TIME_DTSEC, ATMOS_do, OCEAN_do, LAND_do, URBAN_do )
    ! setup external in
    call FILE_EXTERNAL_INPUT_CARTESC_setup

    ! setup nesting grid
    call COMM_CARTESC_NEST_setup ( QA_MP, ATMOS_PHY_MP_TYPE, intercomm_parent, intercomm_child )

    ! setup coriolis parameter
    call CORIOLIS_setup( IA, JA, REAL_LAT(:,:), CY(:), DOMAIN_CENTER_Y )

    ! setup common tools
    call ATMOS_HYDROSTATIC_setup
    call ATMOS_THERMODYN_setup
    call ATMOS_SATURATION_setup

    call BULKFLUX_setup( sqrt(DX**2+DY**2) )

    ! setup variable container
    if ( ATMOS_do ) call ATMOS_vars_setup
    if ( OCEAN_do ) call OCEAN_vars_setup
    if ( LAND_do  ) call LAND_vars_setup
    if ( URBAN_do ) call URBAN_vars_setup
    if ( CPL_sw   ) call CPL_vars_setup

    ! setup driver
    if ( ATMOS_do ) call ATMOS_driver_setup
    if ( OCEAN_do ) call OCEAN_driver_setup
    if ( LAND_do  ) call LAND_driver_setup
    if ( URBAN_do ) call URBAN_driver_setup

    call USER_setup

    call PROF_rapend('Initialize', 0)

    !########## main ##########

#ifdef FIPP
    call fipp_start
#endif
#ifdef PAPI
    call PROF_PAPI_rapstart
#endif

    LOG_NEWLINE
    LOG_PROGRESS(*) 'START TIMESTEP'
    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_Loop', 0)

    fpm_counter = 0
    do
      ! report current time
      call ADMIN_TIME_checkstate

      if ( TIME_DOresume ) then
         ! set state from restart files
         call restart_read

         ! history&monitor file output
         call MONITOR_write('MAIN', TIME_NOWSTEP)
         call FILE_HISTORY_write ! if needed
      endif

      ! time advance
      call ADMIN_TIME_advance
      call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

      ! change to next state
      if( OCEAN_do .AND. TIME_DOOCEAN_step ) call OCEAN_driver_update
      if( LAND_do  .AND. TIME_DOLAND_step  ) call LAND_driver_update
      if( URBAN_do .AND. TIME_DOURBAN_step ) call URBAN_driver_update
      if( ATMOS_do .AND. TIME_DOATMOS_step ) call ATMOS_driver_update( TIME_DOend )
                                             call USER_update
      ! restart & monitor output
      if ( ATMOS_do ) call ATMOS_vars_monitor
      call ADMIN_restart_write
      call MONITOR_write('MAIN', TIME_NOWSTEP)

      ! calc tendencies and diagnostices
      if( ATMOS_do .AND. TIME_DOATMOS_step ) call ATMOS_driver_calc_tendency( force = .false. )
      if( OCEAN_do .AND. TIME_DOOCEAN_step ) call OCEAN_driver_calc_tendency( force = .false. )
      if( LAND_do  .AND. TIME_DOLAND_step  ) call LAND_driver_calc_tendency( force = .false. )
      if( URBAN_do .AND. TIME_DOURBAN_step ) call URBAN_driver_calc_tendency( force = .false. )
      if( CPL_sw   .AND. TIME_DOATMOS_step ) call ATMOS_driver_calc_tendency_from_sflux( force = .false. )
                                             call USER_calc_tendency

      ! history file output
      !   Set hydrostatic pressure coordinate
      if ( ATMOS_do ) call ATMOS_vars_history_setpres
      if ( ATMOS_do ) call ATMOS_vars_history
      if ( OCEAN_do ) call OCEAN_vars_history
      if ( LAND_do  ) call LAND_vars_history
      if ( URBAN_do ) call URBAN_vars_history

      call FILE_HISTORY_write

      if( TIME_DOend ) exit

      if( IO_L ) call flush(IO_FID_LOG)

      ! FPM polling
      if ( FPM_alive .AND. FPM_POLLING_FREQ > 0 ) then
         if ( fpm_counter > FPM_POLLING_FREQ ) then
            sign_exit = .false.
            call FPM_Polling( .true., sign_exit )
            if ( sign_exit ) then
               LOG_ERROR("scalerm",*) 'receive stop signal'
               call PRC_abort
            endif
            fpm_counter = 0
         endif

         fpm_counter = fpm_counter + 1
      endif

    enddo

    call PROF_rapend('Main_Loop', 0)

    LOG_PROGRESS(*) 'END TIMESTEP'
    LOG_NEWLINE

#ifdef FIPP
    call fipp_stop
#endif
#ifdef PAPI
    call PROF_PAPI_rapstop
#endif

    !########## Finalize ##########

    call PROF_setprefx('FIN')

    call PROF_rapstart('All', 1)

    if( ATMOS_do ) call ATMOS_driver_finalize

    ! check data
    if( ATMOS_RESTART_CHECK ) call ATMOS_vars_restart_check

    call PROF_rapstart('Monit', 2)
    call MONITOR_finalize
    call PROF_rapend  ('Monit', 2)

    call PROF_rapstart('File', 2)
    call FILE_HISTORY_finalize
    ! clean up resource allocated for I/O
    call FILE_CARTESC_cleanup

    call COMM_cleanup

    call FILE_Close_All
    call PROF_rapend  ('File', 2)

    call PROF_rapend  ('All', 1)

    call PROF_rapreport
#ifdef PAPI
    call PROF_PAPI_rapreport
#endif

    return
  end subroutine rm_driver

  !-----------------------------------------------------------------------------
  subroutine restart_read
    use scale_atmos_grid_cartesC_index
    use scale_atmos_grid_cartesC, only: &
       CZ   => ATMOS_GRID_CARTESC_CZ,  &
       FZ   => ATMOS_GRID_CARTESC_FZ,  &
       FDZ  => ATMOS_GRID_CARTESC_FDZ, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ,  &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ,  &
       REAL_PHI => ATMOS_GRID_CARTESC_REAL_PHI, &
       AREA     => ATMOS_GRID_CARTESC_REAL_AREA
    use scale_time, only: &
       TIME_NOWDAYSEC
    use mod_admin_restart, only: &
       ADMIN_restart_read
    use mod_atmos_admin, only: &
       ATMOS_do
    use mod_atmos_driver, only: &
       ATMOS_driver_calc_tendency, &
       ATMOS_driver_calc_tendency_from_sflux, &
       ATMOS_SURFACE_SET
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics, &
       ATMOS_vars_history_setpres,  &
       ATMOS_vars_history,          &
       ATMOS_vars_monitor,          &
       DENS,                        &
       POTT,                        &
       TEMP,                        &
       PRES,                        &
       QV
    use mod_atmos_bnd_driver, only: &
       ATMOS_BOUNDARY_driver_set
    use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_UPDATE
    use mod_ocean_admin, only: &
       OCEAN_do
    use mod_ocean_driver, only: &
       OCEAN_driver_calc_tendency, &
       OCEAN_SURFACE_SET
    use mod_ocean_vars, only: &
       OCEAN_vars_history
    use mod_land_admin, only: &
       LAND_do
    use mod_land_driver, only: &
       LAND_driver_calc_tendency, &
       LAND_SURFACE_SET
    use mod_land_vars, only: &
       LAND_vars_history
    use mod_urban_admin, only: &
       URBAN_do
    use mod_urban_driver, only: &
       URBAN_driver_calc_tendency, &
       URBAN_SURFACE_SET
    use mod_urban_vars, only: &
       URBAN_vars_history
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_user, only: &
       USER_calc_tendency
    implicit none
    !---------------------------------------------------------------------------

    ! read restart data
    call ADMIN_restart_read

    if ( ATMOS_do ) then
       ! calc diagnostics
       call ATMOS_vars_calc_diagnostics
       call ATMOS_REFSTATE_update( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                   DENS(:,:,:), POTT(:,:,:), TEMP(:,:,:), PRES(:,:,:), QV(:,:,:), & ! [IN]
                                   CZ(:), FZ(:), FDZ(:), RCDZ(:),                                 & ! [IN]
                                   REAL_CZ(:,:,:), REAL_FZ(:,:,:), REAL_PHI(:,:,:), AREA(:,:),    & ! [IN]
                                   TIME_NOWDAYSEC,                                                & ! [IN]
                                   force = .true.                                                 )
       call ATMOS_BOUNDARY_driver_set
       call ATMOS_vars_history_setpres
    endif

    ! setup surface condition
    if( ATMOS_do ) call ATMOS_SURFACE_SET( countup=.false. )
    if( OCEAN_do ) call OCEAN_SURFACE_SET( countup=.false. )
    if( LAND_do  ) call LAND_SURFACE_SET ( countup=.false. )
    if( URBAN_do ) call URBAN_SURFACE_SET( countup=.false. )

    ! calc tendencies
    if( ATMOS_do ) call ATMOS_driver_calc_tendency           ( force=.true. )
    if( OCEAN_do ) call OCEAN_driver_calc_tendency           ( force=.true. )
    if( LAND_do  ) call LAND_driver_calc_tendency            ( force=.true. )
    if( URBAN_do ) call URBAN_driver_calc_tendency           ( force=.true. )
    if( CPL_sw   ) call ATMOS_driver_calc_tendency_from_sflux( force=.true. )
                   call USER_calc_tendency

    !########## History & Monitor ##########
    if( ATMOS_do ) call ATMOS_vars_history
    if( OCEAN_do ) call OCEAN_vars_history
    if( LAND_do  ) call LAND_vars_history
    if( URBAN_do ) call URBAN_vars_history

    call ATMOS_vars_monitor

    return
  end subroutine restart_read

end module mod_rm_driver
