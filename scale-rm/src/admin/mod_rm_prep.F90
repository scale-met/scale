!-------------------------------------------------------------------------------
!> module SCALE-RM prep
!!
!! @par Description
!!          This program is driver of preprocess tools
!!          1) boundary data (e.g. topography, land use index)
!!          2) initial data for ideal/real test cases
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_rm_prep
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
  public :: rm_prep

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
  subroutine rm_prep( &
       comm_world,       &
       intercomm_parent, &
       intercomm_child,  &
       cnf_fname         )
    use scale_file, only: &
       FILE_Close_All
    use scale_prc, only: &
       PRC_LOCAL_setup
    use scale_prc_cartesC, only: &
       PRC_CARTESC_setup
    use scale_const, only: &
       CONST_setup
    use scale_calendar, only: &
       CALENDAR_setup
    use scale_random, only: &
       RANDOM_setup
    use scale_atmos_grid_cartesC_index, only: &
       ATMOS_GRID_CARTESC_INDEX_setup
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_setup
    use scale_comm_cartesC_nest, only: &
       COMM_CARTESC_NEST_setup
    use scale_ocean_grid_cartesC_index, only: &
       OCEAN_GRID_CARTESC_INDEX_setup
    use scale_ocean_grid_cartesC, only: &
       OCEAN_GRID_CARTESC_setup
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_setup
    use scale_land_grid_cartesC_index, only: &
       LAND_GRID_CARTESC_INDEX_setup
    use scale_land_grid_cartesC, only: &
       LAND_GRID_CARTESC_setup
    use scale_land_grid_cartesC_real, only: &
       LAND_GRID_CARTESC_REAL_setup
    use scale_urban_grid_cartesC_index, only: &
       URBAN_GRID_CARTESC_INDEX_setup
    use scale_urban_grid_cartesC, only: &
       URBAN_GRID_CARTESC_setup
    use scale_urban_grid_cartesC_real, only: &
       URBAN_GRID_CARTESC_REAL_setup
    use scale_file_cartesC, only: &
       FILE_CARTESC_setup, &
       FILE_CARTESC_cleanup
    use scale_comm_cartesC, only: &
       COMM_setup
    use scale_topography, only: &
       TOPO_setup
    use scale_landuse, only: &
       LANDUSE_setup
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_setup,   &
       ATMOS_GRID_CARTESC_REAL_update_Z
    use scale_atmos_grid_cartesC_metric, only: &
       ATMOS_GRID_CARTESC_METRIC_setup
    use scale_statistics, only: &
       STATISTICS_setup
    use scale_atmos_hydrostatic, only: &
       ATMOS_HYDROSTATIC_setup
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_setup
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_setup
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_setup
    use mod_atmos_driver, only: &
       ATMOS_driver_tracer_setup
    use mod_admin_restart, only: &
       ADMIN_restart_setup
    use mod_admin_versioncheck, only: &
       ADMIN_versioncheck
    use mod_admin_time, only: &
       ADMIN_TIME_setup
    use mod_atmos_admin, only: &
       ATMOS_admin_setup, &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_do
    use mod_atmos_vars, only: &
       ATMOS_vars_setup
    use mod_atmos_phy_mp_vars, only: &
       QA_MP
    use mod_ocean_admin, only: &
       OCEAN_admin_setup, &
       OCEAN_do
    use mod_ocean_vars, only: &
       OCEAN_vars_setup
    use mod_land_admin, only: &
       LAND_admin_setup, &
       LAND_do
    use mod_land_vars, only: &
       LAND_vars_setup
    use mod_urban_admin, only: &
       URBAN_admin_setup, &
       URBAN_do, &
       URBAN_land
    use mod_urban_vars, only: &
       URBAN_vars_setup
    use mod_lake_admin, only: &
       LAKE_admin_setup, &
       LAKE_do
    use mod_cpl_admin, only: &
       CPL_admin_setup, &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_vars_setup
    use mod_convert, only: &
       CONVERT_setup, &
       CONVERT
    use mod_mktopo, only: &
       MKTOPO_setup, &
       MKTOPO
    use mod_mkinit, only: &
       MKINIT_setup, &
       MKINIT
    use mod_user, only: &
       USER_tracer_setup, &
       USER_mkinit
    use mod_atmos_driver, only: &
       ATMOS_SURFACE_GET
    use mod_ocean_driver, only: &
       OCEAN_SURFACE_SET
    use mod_land_driver, only: &
       LAND_SURFACE_SET
    use mod_urban_driver, only: &
       URBAN_SURFACE_SET
    use mod_admin_restart, only: &
       ADMIN_restart_write
    use mod_admin_time, only: &
       TIME_DOATMOS_restart,  &
       TIME_DOLAND_restart,   &
       TIME_DOURBAN_restart,  &
       TIME_DOOCEAN_restart
    implicit none

    integer,          intent(in) :: comm_world
    integer,          intent(in) :: intercomm_parent
    integer,          intent(in) :: intercomm_child
    character(len=*), intent(in) :: cnf_fname

    integer :: myrank
    logical :: ismaster

    logical :: output
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


    ! profiler start
    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)


    ! setup constants
    call CONST_setup

    ! setup calendar
    call CALENDAR_setup

    ! setup random number
    call RANDOM_setup

    ! setup horizontal/vertical grid coordinates (cartesian,idealized)
    call ATMOS_GRID_CARTESC_INDEX_setup
    call ATMOS_GRID_CARTESC_setup

    call OCEAN_GRID_CARTESC_INDEX_setup
    call OCEAN_GRID_CARTESC_setup

    call LAND_GRID_CARTESC_INDEX_setup
    call LAND_GRID_CARTESC_setup

    call URBAN_GRID_CARTESC_INDEX_setup
    call URBAN_GRID_CARTESC_setup

    ! setup submodel administrator
    call ATMOS_admin_setup
    call OCEAN_admin_setup
    call LAND_admin_setup
    call URBAN_admin_setup
    call LAKE_admin_setup
    call CPL_admin_setup

    ! setup tracer index
    call ATMOS_HYDROMETEOR_setup
    call ATMOS_driver_tracer_setup
    call USER_tracer_setup

    ! setup file I/O
    call FILE_CARTESC_setup

    ! setup mpi communication
    call COMM_setup

    ! setup topography
    call TOPO_setup
    ! setup land use category index/fraction
    call LANDUSE_setup( OCEAN_do, (.not. URBAN_land), LAKE_do )
    ! setup grid coordinates (real world)
    call ATMOS_GRID_CARTESC_REAL_setup

    ! setup grid transfer metrics (uses in ATMOS_dynamics)
    call ATMOS_GRID_CARTESC_METRIC_setup

    call OCEAN_GRID_CARTESC_REAL_setup
    call LAND_GRID_CARTESC_REAL_setup
    call URBAN_GRID_CARTESC_REAL_setup

    ! setup restart
    call ADMIN_restart_setup
    ! setup time
    call ADMIN_TIME_setup( setup_TimeIntegration = .false. )
    ! setup statistics
    call STATISTICS_setup

    ! setup nesting grid
    call COMM_CARTESC_NEST_setup ( QA_MP, ATMOS_PHY_MP_TYPE, intercomm_parent, intercomm_child )


    ! setup common tools
    call ATMOS_HYDROSTATIC_setup
    call ATMOS_THERMODYN_setup
    call ATMOS_SATURATION_setup

    ! setup variable container
    if ( ATMOS_do ) call ATMOS_vars_setup
    if ( OCEAN_do ) call OCEAN_vars_setup
    if ( LAND_do  ) call LAND_vars_setup
    if ( URBAN_do ) call URBAN_vars_setup
    if ( CPL_sw   ) call CPL_vars_setup

    ! setup preprocess converter
    call CONVERT_setup

    ! setup mktopo
    call MKTOPO_setup

    ! setup mkinit
    call MKINIT_setup

    call PROF_rapend('Initialize',0)

    !########## main ##########

    call PROF_setprefx('MAIN')
    call PROF_rapstart('Main_prep',0)

    ! execute preprocess
    call PROF_rapstart('Convert',1)
    call CONVERT
    call PROF_rapend  ('Convert',1)

    ! execute mktopo
    call PROF_rapstart('MkTopo',1)
    call MKTOPO
    call PROF_rapend  ('MkTopo',1)

    ! re-setup
    call ATMOS_GRID_CARTESC_REAL_update_Z

    ! execute mkinit
    call PROF_rapstart('MkInit',1)
    call MKINIT( output )
    call USER_mkinit
    call PROF_rapend  ('MkInit',1)

    if ( output ) then
      call PROF_rapstart('MkInit_restart',1)

      ! setup surface condition
      call OCEAN_SURFACE_SET( countup = .false. )
      call LAND_SURFACE_SET ( countup = .false. )
      call URBAN_SURFACE_SET( countup = .false. )
      call ATMOS_SURFACE_GET

      ! output restart file
      TIME_DOOCEAN_restart = .TRUE.
      TIME_DOLAND_restart  = .TRUE.
      TIME_DOURBAN_restart = .TRUE.
      TIME_DOATMOS_restart = .TRUE.
      call ADMIN_restart_write

      call PROF_rapend  ('MkInit_restart',1)

   end if

    call PROF_rapend('Main_prep',0)

    !########## Finalize ##########

    ! setup file I/O
    call FILE_CARTESC_cleanup

    call FILE_Close_All

    call PROF_rapreport

    return
  end subroutine rm_prep

end module mod_rm_prep
