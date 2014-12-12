!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land model driver
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_land_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_driver_setup
  public :: LAND_driver
  public :: LAND_SURFACE_SET

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
  subroutine LAND_driver_setup
    use mod_land_phy_driver, only: &
       LAND_PHY_driver_setup
    use mod_land_vars, only: &
       LAND_vars_history
    use mod_land_admin, only: &
       LAND_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[LAND] / Origin[SCALE-LES]'

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET( setup=.true. )

    call LAND_PHY_driver_setup

    !########## History & Monitor ##########
    if ( LAND_sw ) then
       call PROF_rapstart('LND History', 1)
       call LAND_vars_history
       call PROF_rapend  ('LND History', 1)
    endif

    return
  end subroutine LAND_driver_setup

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_driver
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use mod_land_admin, only: &
       LAND_sw
    use mod_land_vars, only: &
       LAND_TEMP,        &
       LAND_WATER,       &
       LAND_TEMP_t,      &
       LAND_WATER_t,     &
       LAND_vars_total,  &
       LAND_vars_history
    use mod_land_phy_driver, only: &
       LAND_PHY_driver
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET( setup=.false. )

    !########## Physics ##########
    if ( LAND_sw ) then
       call PROF_rapstart('LND Physics', 1)
       call LAND_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('LND Physics', 1)
    endif

    !########## Update ##########
    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
       LAND_TEMP (k,i,j) = LAND_TEMP (k,i,j) + LAND_TEMP_t (k,i,j) * dt
       LAND_WATER(k,i,j) = LAND_WATER(k,i,j) + LAND_WATER_t(k,i,j) * dt
    enddo
    enddo
    enddo

    !########## Negative Fixer ##########
    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
       LAND_WATER(k,i,j) = max( LAND_WATER(k,i,j), 0.0_RP )
    enddo
    enddo
    enddo

    call LAND_vars_total

    !########## Put Surface Boundary to coupler ##########
    call LAND_SURFACE_SET( setup=.false. )

    !########## History & Monitor ##########
    call PROF_rapstart('LND History', 1)
    call LAND_vars_history
    call PROF_rapend  ('LND History', 1)

    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
       LAND_TEMP_t (k,i,j) = 0.0_RP
       LAND_WATER_t(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    return
  end subroutine LAND_driver

  !-----------------------------------------------------------------------------
  !> Get surface boundary
  subroutine LAND_SURFACE_GET( setup )
    use mod_land_vars, only: &
       LAND_SFC_TEMP,      &
       LAND_SFC_albedo,    &
       LAND_SFLX_heat,     &
       LAND_SFLX_prec,     &
       LAND_SFLX_evap
    use mod_cpl_admin, only: &
       CPL_sw_AtmLnd
    use mod_cpl_vars, only: &
       CPL_getLnd
    implicit none

    logical, intent(in) :: setup
    !---------------------------------------------------------------------------

    if ( CPL_sw_AtmLnd ) then
       call CPL_getLnd( LAND_SFC_TEMP  (:,:),   & ! [OUT]
                        LAND_SFC_albedo(:,:,:), & ! [OUT]
                        LAND_SFLX_heat (:,:),   & ! [OUT]
                        LAND_SFLX_prec (:,:),   & ! [OUT]
                        LAND_SFLX_evap (:,:)    ) ! [OUT]
    endif

    return
  end subroutine LAND_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine LAND_SURFACE_SET( setup )
    use scale_land_grid, only: &
       dz => GRID_LCDZ
    use mod_admin_restart, only: &
       RESTART_RUN
    use mod_land_vars, only: &
       I_WaterCritical,    &
       I_ThermalCond,      &
       I_Z0M,              &
       I_Z0H,              &
       I_Z0E,              &
       LAND_TEMP,          &
       LAND_WATER,         &
       LAND_SFC_TEMP,      &
       LAND_SFC_albedo,    &
       LAND_SFLX_heat,     &
       LAND_SFLX_prec,     &
       LAND_SFLX_evap,     &
       LAND_PROPERTY
    use mod_cpl_admin, only: &
       CPL_sw_AtmLnd
    use mod_cpl_vars, only: &
       CPL_putLnd_setup,   &
       CPL_putLnd_restart, &
       CPL_putLnd
    implicit none

    logical, intent(in) :: setup

    real(RP) :: LAND_DZ     (IA,JA)
    real(RP) :: LAND_TEMP_Z1(IA,JA)
    real(RP) :: LAND_BETA   (IA,JA)

    real(RP), parameter :: BETA_MAX = 1.0_RP
    !---------------------------------------------------------------------------

    if ( CPL_sw_AtmLnd ) then
       if ( setup ) then
          LAND_DZ(:,:) = dz(LKS)

          call CPL_putLnd_setup( LAND_TEMP_Z1(:,:),                  & ! [IN]
                                 LAND_BETA   (:,:),                  & ! [IN]
                                 LAND_SFC_TEMP  (:,:),               & ! [IN]
                                 LAND_SFC_albedo(:,:,:),             & ! [IN]
                                 LAND_PROPERTY  (:,:,I_ThermalCond), & ! [IN]
                                 LAND_DZ        (:,:),               & ! [IN]
                                 LAND_PROPERTY  (:,:,I_Z0M),         & ! [IN]
                                 LAND_PROPERTY  (:,:,I_Z0H),         & ! [IN]
                                 LAND_PROPERTY  (:,:,I_Z0E)          ) ! [IN]

          if ( RESTART_RUN ) then
             call CPL_putLnd_restart( LAND_SFLX_heat (:,:),   & ! [IN]
                                      LAND_SFLX_prec (:,:),   & ! [IN]
                                      LAND_SFLX_evap (:,:)    ) ! [IN]
          endif
       else
          LAND_TEMP_Z1(:,:) = LAND_TEMP(LKS,:,:)

          LAND_BETA   (:,:) = min( LAND_WATER(LKS,:,:) / LAND_PROPERTY(:,:,I_WaterCritical), BETA_MAX )

          call CPL_putLnd( LAND_TEMP_Z1(:,:), & ! [IN]
                           LAND_BETA   (:,:)  ) ! [IN]
       endif
    endif

    return
  end subroutine LAND_SURFACE_SET

end module mod_land_driver
