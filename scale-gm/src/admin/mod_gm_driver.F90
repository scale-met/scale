!-------------------------------------------------------------------------------
!> module SCALE-GM (a main routine of global model)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Numerical model for global scale weather/climate
!!          based on the Non-hydrostatic ICosahedral Atmosphere Model (NICAM)
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_gm_driver
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
#include "scale-gm.h"
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: gm_driver

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
  character(len=H_MID), private, parameter :: MODELNAME = "SCALE-GM ver. "//VERSION

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine gm_driver( &
       comm_world,       &
       intercomm_parent, &
       intercomm_child,  &
       cnf_fname         )
    use scale_file, only: &
       FILE_Close_All
    use scale_prc, only: &
       PRC_abort,       &
       PRC_LOCAL_setup
    use scale_fpm, only: &
       FPM_alive,       &
       FPM_Polling,     &
       FPM_POLLING_FREQ
    use scale_prc_icoA, only: &
       PRC_ICOA_setup
    use scale_const, only: &
       CONST_setup,            &
       CONST_THERMODYN_TYPE,   &
       RADIUS => CONST_RADIUS, &
       PI     => CONST_PI
    use scale_calendar, only: &
       CALENDAR_setup
    use scale_random, only: &
       RANDOM_setup
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_setup
    use scale_atmos_grid_icoA_index, only: &
       ATMOS_GRID_icoA_INDEX_setup, &
       ADM_glevel
    use scale_atmos_grid_icoA, only: &
       ATMOS_GRID_icoA_setup
    use scale_ocean_grid_icoA_index, only: &
       OCEAN_GRID_icoA_INDEX_setup
    use scale_ocean_grid_icoA, only: &
       OCEAN_GRID_icoA_setup
    use scale_land_grid_icoA_index, only: &
       LAND_GRID_icoA_INDEX_setup
    use scale_land_grid_icoA, only: &
       LAND_GRID_icoA_setup
    use scale_urban_grid_icoA_index, only: &
       URBAN_GRID_icoA_INDEX_setup
    use scale_urban_grid_icoA, only: &
       URBAN_GRID_icoA_setup
    use scale_comm_icoA, only: &
       COMM_setup
    use scale_landuse, only: &
       LANDUSE_setup
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_setup
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_setup
    use scale_bulkflux, only: &
       BULKFLUX_setup
    use mod_fio, only: &
       FIO_setup
    use mod_grd, only: &
       GRD_setup
    use mod_gmtr, only: &
       GMTR_setup
    use mod_oprt, only: &
       OPRT_setup
    use mod_vmtr, only: &
       VMTR_setup
    use mod_time, only: &
       TIME_setup,     &
       TIME_report,    &
       TIME_advance,   &
       TIME_LSTEP_MAX
    use mod_extdata, only: &
       extdata_setup
    use mod_runconf, only: &
       runconf_setup
    use mod_prgvar, only: &
       prgvar_setup,            &
       restart_input_basename,  &
       restart_output_basename, &
       restart_input,           &
       restart_output
    use mod_dynamics, only: &
       dynamics_setup, &
       dynamics_step
    use mod_forcing_driver, only: &
       forcing_setup, &
       forcing_step
    use mod_history, only: &
       history_setup, &
       history_out,   &
       HIST_output_step0
    use mod_history_vars, only: &
       history_vars_setup, &
       history_vars
    use mod_embudget, only: &
       embudget_setup, &
       embudget_monitor
    use mod_gm_topography, only: &
       TOPO_setup
    use mod_atmos_vars, only: &
       ATMOS_vars_setup
    use mod_atmos_phy_driver, only: &
       ATMOS_phy_driver_tracer_setup, &
       ATMOS_phy_driver_setup, &
       ATMOS_phy_driver
    use mod_atmos_phy_sf_vars, only: &
       ATMOS_PHY_SF_vars_setup
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_setup
    use mod_time, only: &
       TIME_RES_ATMOS_PHY_RD, &
       TIME_DOATMOS_PHY_RD
    use mod_atmos_admin, only: &
       ATMOS_admin_setup, &
       ATMOS_do
    use mod_ocean_admin, only: &
       OCEAN_admin_setup, &
       OCEAN_do
    use mod_land_admin, only: &
       LAND_admin_setup, &
       LAND_do
    use mod_urban_admin, only: &
       URBAN_admin_setup, &
       URBAN_do
    use mod_lake_admin, only: &
       LAKE_admin_setup, &
       LAKE_do
    use mod_cpl_admin, only: &
       CPL_admin_setup, &
       CPL_sw
    implicit none

    integer,          intent(in) :: comm_world
    integer,          intent(in) :: intercomm_parent
    integer,          intent(in) :: intercomm_child
    character(len=*), intent(in) :: cnf_fname

    integer :: myrank
    integer :: fpm_counter
    logical :: ismaster
    logical :: sign_exit

    real(RP) :: dx, dy
    integer  :: n
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

    ! setup process
    call PRC_ICOA_setup

    ! setup PROF
    call PROF_setup

    !###########################################################################
    ! profiler start
    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)

    ! setup constants
    call CONST_setup
    CONST_THERMODYN_TYPE='EXACT'

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

    ! setup horizontal/vertical grid coordinates (icosahedral,idealized)
    if ( ATMOS_do ) then
       call ATMOS_GRID_icoA_INDEX_setup
       call ATMOS_GRID_icoA_setup
    endif

    if ( OCEAN_do ) then
       call OCEAN_GRID_icoA_INDEX_setup
       !call OCEAN_GRID_icoA_setup
    endif

    if ( LAND_do ) then
       call LAND_GRID_icoA_INDEX_setup
       !call LAND_GRID_icoA_setup
    endif

    if ( URBAN_do ) then
       call URBAN_GRID_icoA_INDEX_setup
       !call URBAN_GRID_icoA_setup
    endif

    !---< tracer setup >---
    call ATMOS_HYDROMETEOR_setup
    call ATMOS_phy_driver_tracer_setup

    !---< I/O module setup >---
    call FIO_setup

    !---< comm module setup >---
    call COMM_setup

    !---< grid module setup >---
    call GRD_setup

    !---< geometrics module setup >---
    call GMTR_setup

    !---< operator module setup >---
    call OPRT_setup

    !---< vertical metrics module setup >---
    call VMTR_setup

    !---< time module setup >---
    call TIME_setup

    !---< external data module setup >---
    call extdata_setup


    !---< forcing module setup >---
    call forcing_setup

    !---< nhm_runconf module setup >---
    call runconf_setup

    !---< topography module setup >---
    call TOPO_setup

    !---< landuse module setup >---
    call LANDUSE_setup( .false., .false., .false. )

    ! setup common tools
    call ATMOS_THERMODYN_setup
    call ATMOS_SATURATION_setup

    dx = sqrt( 4.0_RP * PI * RADIUS**2 / real(10*4**ADM_glevel+2,kind=RP) )
    dy = dx
    call BULKFLUX_setup( sqrt(dx**2+dy**2) )

    !---< prognostic variable module setup >---
    call prgvar_setup
    call restart_input( restart_input_basename )

    !---< scale variable module setup >---
    call ATMOS_vars_setup

    !---< dynamics module setup >---
    call dynamics_setup

    !---< surface module setup >---
    call atmos_phy_sf_driver_setup

    !---< physics module setup >---
    call ATMOS_phy_driver_setup
    TIME_RES_ATMOS_PHY_RD = 0
    TIME_DOATMOS_PHY_RD   = .true.

    !---< energy&mass budget module setup >---
    call embudget_setup

    !---< history module setup >---
    call history_setup

    !---< history variable module setup >---
    call history_vars_setup

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

    !--- history output at initial time
    if ( HIST_output_step0 ) then
       call TIME_advance( reverse=.true. )
       call TIME_report
       call history_vars
       call TIME_advance
       call history_out
    endif

    fpm_counter = 0

    do n = 1, TIME_LSTEP_MAX

       call TIME_report

       call PROF_rapstart('_Atmos',1)
       call dynamics_step
       call ATMOS_phy_driver
       call forcing_step
       call PROF_rapend  ('_Atmos',1)

       call PROF_rapstart('_History',1)
       call history_vars
       call embudget_monitor

       call TIME_advance

       call history_out

       call restart_output( restart_output_basename, n )

       call PROF_rapend  ('_History',1)

      if( IO_L ) call flush(IO_FID_LOG)

      ! FPM polling
      if ( FPM_alive ) then
      if ( FPM_POLLING_FREQ > 0 ) then
         if ( fpm_counter > FPM_POLLING_FREQ ) then
            sign_exit = .false.
            call FPM_Polling( .true., sign_exit )
            if ( sign_exit ) then
               LOG_ERROR("gm_driver",*) 'receive stop signal'
               call PRC_abort
            endif
            fpm_counter = 0
         endif
         fpm_counter = fpm_counter + 1
      endif
      endif

    enddo

    call PROF_rapend('Main_Loop', 0)

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END TIMESTEP ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

#ifdef FIPP
    call fipp_stop
#endif
#ifdef PAPI
    call PROF_PAPI_rapstop
#endif

    !########## Finalize ##########

    call PROF_setprefx('FIN')

    call PROF_rapstart('All', 1)

    call PROF_rapstart('File', 2)
    call FILE_Close_All
    call PROF_rapend  ('File', 2)

    call PROF_rapend  ('All', 1)

    call PROF_rapreport
#ifdef PAPI
    call PROF_PAPI_rapreport
#endif

    return
  end subroutine gm_driver

end module mod_gm_driver
