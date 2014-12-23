!-------------------------------------------------------------------------------
!> module OCEAN driver
!!
!! @par Description
!!          Ocean model driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_ocean_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index

  use scale_const, only: &
     I_SW  => CONST_I_SW, &
     I_LW  => CONST_I_LW
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_driver_setup
  public :: OCEAN_driver
  public :: OCEAN_SURFACE_GET
  public :: OCEAN_SURFACE_SET

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_driver_setup
    use mod_ocean_phy_driver, only: &
       OCEAN_PHY_driver_setup
!    use mod_ocean_frc_nudge, only: &
!       OCEAN_FRC_driver_setup
    use mod_ocean_vars, only: &
       OCEAN_vars_history
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN] / Origin[SCALE-LES]'

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    call OCEAN_PHY_driver_setup

!    if( OCEAN_FRC_sw ) call OCEAN_FRC_driver_setup

    !########## Set Surface Boundary to coupler ##########
    call OCEAN_SURFACE_SET( countup=.true. )

    !########## History & Monitor ##########
    if ( OCEAN_sw ) then
       call PROF_rapstart('OCN History', 1)
       call OCEAN_vars_history
       call PROF_rapend  ('OCN History', 1)
    endif

    return
  end subroutine OCEAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_driver
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_ocean_vars, only: &
       OCEAN_TEMP,         &
       OCEAN_SFC_TEMP,     &
       OCEAN_SFC_albedo,   &
       OCEAN_SFC_Z0M,      &
       OCEAN_SFC_Z0H,      &
       OCEAN_SFC_Z0E,      &
       OCEAN_TEMP_t,       &
       OCEAN_SFC_TEMP_t,   &
       OCEAN_SFC_albedo_t, &
       OCEAN_SFC_Z0M_t,    &
       OCEAN_SFC_Z0H_t,    &
       OCEAN_SFC_Z0E_t,    &
       OCEAN_vars_total,   &
       OCEAN_vars_history
    use mod_ocean_phy_driver, only: &
       OCEAN_PHY_driver
!    use mod_ocean_forcing, only: &
!       OCEAN_forcing
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    !########## Physics ##########
    if ( OCEAN_sw ) then
       call PROF_rapstart('OCN Physics', 1)
       call OCEAN_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('OCN Physics', 1)
    endif

    !########## Forcing ##########
!    if ( OCEAN_FORCE_sw ) then
!       call PROF_rapstart('OCN Forcing', 1)
!       call OCEAN_forcing
!       call PROF_rapend  ('OCN Forcing', 1)
!    endif

    !########## Update ##########
    do j = JS, JE
    do i = IS, IE
       OCEAN_TEMP      (i,j)      = OCEAN_TEMP      (i,j)      + OCEAN_TEMP_t      (i,j)      * dt
       OCEAN_SFC_TEMP  (i,j)      = OCEAN_SFC_TEMP  (i,j)      + OCEAN_SFC_TEMP_t  (i,j)      * dt
       OCEAN_SFC_albedo(i,j,I_LW) = OCEAN_SFC_albedo(i,j,I_LW) + OCEAN_SFC_albedo_t(i,j,I_LW) * dt
       OCEAN_SFC_albedo(i,j,I_SW) = OCEAN_SFC_albedo(i,j,I_SW) + OCEAN_SFC_albedo_t(i,j,I_SW) * dt
       OCEAN_SFC_Z0M   (i,j)      = OCEAN_SFC_Z0M   (i,j)      + OCEAN_SFC_Z0M_t   (i,j)      * dt
       OCEAN_SFC_Z0H   (i,j)      = OCEAN_SFC_Z0H   (i,j)      + OCEAN_SFC_Z0H_t   (i,j)      * dt
       OCEAN_SFC_Z0E   (i,j)      = OCEAN_SFC_Z0E   (i,j)      + OCEAN_SFC_Z0E_t   (i,j)      * dt
    enddo
    enddo

    call OCEAN_vars_total

    !########## Set Surface Boundary to coupler ##########
    call OCEAN_SURFACE_SET( countup=.true. )

    !########## reset tendencies ##########
    do j = JS, JE
    do i = IS, IE
       OCEAN_TEMP_t      (i,j)      = 0.0_RP
       OCEAN_SFC_TEMP_t  (i,j)      = 0.0_RP
       OCEAN_SFC_albedo_t(i,j,I_LW) = 0.0_RP
       OCEAN_SFC_albedo_t(i,j,I_SW) = 0.0_RP
       OCEAN_SFC_Z0M_t   (i,j)      = 0.0_RP
       OCEAN_SFC_Z0H_t   (i,j)      = 0.0_RP
       OCEAN_SFC_Z0E_t   (i,j)      = 0.0_RP
    enddo
    enddo

    !########## History & Monitor ##########
    call PROF_rapstart('OCN History', 1)
    call OCEAN_vars_history
    call PROF_rapend  ('OCN History', 1)

    return
  end subroutine OCEAN_driver

  !-----------------------------------------------------------------------------
  !> Get surface boundary from other model
  subroutine OCEAN_SURFACE_GET
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_ocean_vars, only: &
       ATMOS_TEMP,     &
       ATMOS_PRES,     &
       ATMOS_W,        &
       ATMOS_U,        &
       ATMOS_V,        &
       ATMOS_DENS,     &
       ATMOS_QV,       &
       ATMOS_PBL,      &
       ATMOS_SFC_PRES, &
       ATMOS_SFLX_LW,  &
       ATMOS_SFLX_SW,  &
       ATMOS_SFLX_prec
    use mod_cpl_vars, only: &
       CPL_getATM_OCN
    implicit none
    !---------------------------------------------------------------------------

    if ( OCEAN_sw ) then
       call CPL_getATM_OCN( ATMOS_TEMP     (:,:), & ! [OUT]
                            ATMOS_PRES     (:,:), & ! [OUT]
                            ATMOS_W        (:,:), & ! [OUT]
                            ATMOS_U        (:,:), & ! [OUT]
                            ATMOS_V        (:,:), & ! [OUT]
                            ATMOS_DENS     (:,:), & ! [OUT]
                            ATMOS_QV       (:,:), & ! [OUT]
                            ATMOS_PBL      (:,:), & ! [OUT]
                            ATMOS_SFC_PRES (:,:), & ! [OUT]
                            ATMOS_SFLX_LW  (:,:), & ! [OUT]
                            ATMOS_SFLX_SW  (:,:), & ! [OUT]
                            ATMOS_SFLX_prec(:,:)  ) ! [OUT]
    endif

    return
  end subroutine OCEAN_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine OCEAN_SURFACE_SET( countup )
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_ocean_vars, only: &
       OCEAN_SFC_TEMP,      &
       OCEAN_SFC_albedo,    &
       OCEAN_SFC_Z0M,       &
       OCEAN_SFC_Z0H,       &
       OCEAN_SFC_Z0E,       &
       OCEAN_SFLX_MW,       &
       OCEAN_SFLX_MU,       &
       OCEAN_SFLX_MV,       &
       OCEAN_SFLX_SH,       &
       OCEAN_SFLX_LH,       &
       OCEAN_SFLX_WH,       &
       OCEAN_SFLX_evap,     &
       OCEAN_U10,           &
       OCEAN_V10,           &
       OCEAN_T2,            &
       OCEAN_Q2
    use mod_cpl_vars, only: &
       CPL_putOCN
    implicit none

    ! arguments
    logical, intent(in) :: countup
    !---------------------------------------------------------------------------

    if ( OCEAN_sw ) then
       call CPL_putOCN( OCEAN_SFC_TEMP  (:,:),   & ! [IN]
                        OCEAN_SFC_albedo(:,:,:), & ! [IN]
                        OCEAN_SFC_Z0M   (:,:),   & ! [IN]
                        OCEAN_SFC_Z0H   (:,:),   & ! [IN]
                        OCEAN_SFC_Z0E   (:,:),   & ! [IN]
                        OCEAN_SFLX_MW   (:,:),   & ! [IN]
                        OCEAN_SFLX_MU   (:,:),   & ! [IN]
                        OCEAN_SFLX_MV   (:,:),   & ! [IN]
                        OCEAN_SFLX_SH   (:,:),   & ! [IN]
                        OCEAN_SFLX_LH   (:,:),   & ! [IN]
                        OCEAN_SFLX_WH   (:,:),   & ! [IN]
                        OCEAN_SFLX_evap (:,:),   & ! [IN]
                        OCEAN_U10       (:,:),   & ! [IN]
                        OCEAN_V10       (:,:),   & ! [IN]
                        OCEAN_T2        (:,:),   & ! [IN]
                        OCEAN_Q2        (:,:),   & ! [IN]
                        countup                  ) ! [IN]
    endif

    return
  end subroutine OCEAN_SURFACE_SET

end module mod_ocean_driver
