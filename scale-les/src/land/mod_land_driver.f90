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
  public :: LAND_driver_setup
  public :: LAND_driver
  public :: LAND_SURFACE_GET
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
    call LAND_SURFACE_GET

    call LAND_PHY_driver_setup

    !########## Set Surface Boundary to coupler ##########
    call LAND_SURFACE_SET( countup=.true. )

    !########## History & Monitor ##########
    if ( LAND_sw ) then
       call PROF_rapstart('LND_History', 1)
       call LAND_vars_history
       call PROF_rapend  ('LND_History', 1)
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
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_SFC_TEMP,     &
       LAND_SFC_albedo,   &
       LAND_TEMP_t,       &
       LAND_WATER_t,      &
       LAND_SFC_TEMP_t,   &
       LAND_SFC_albedo_t, &
       LAND_vars_total,   &
       LAND_vars_history
    use mod_land_phy_driver, only: &
       LAND_PHY_driver
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET

    !########## Physics ##########
    if ( LAND_sw ) then
       call PROF_rapstart('LND_Physics', 1)
       call LAND_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('LND_Physics', 1)
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
    do j = JS, JE
    do i = IS, IE
       LAND_SFC_TEMP  (i,j)      = LAND_SFC_TEMP  (i,j)      + LAND_SFC_TEMP_t  (i,j)      * dt
       LAND_SFC_albedo(i,j,I_LW) = LAND_SFC_albedo(i,j,I_LW) + LAND_SFC_albedo_t(i,j,I_LW) * dt
       LAND_SFC_albedo(i,j,I_SW) = LAND_SFC_albedo(i,j,I_SW) + LAND_SFC_albedo_t(i,j,I_SW) * dt
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

    !########## Set Surface Boundary to coupler ##########
    call LAND_SURFACE_SET( countup=.true. )

    !########## reset tendencies ##########
    do j = JS, JE
    do i = IS, IE
    do k = LKS, LKE
       LAND_TEMP_t (k,i,j) = 0.0_RP
       LAND_WATER_t(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       LAND_SFC_TEMP_t  (i,j)      = 0.0_RP
       LAND_SFC_albedo_t(i,j,I_LW) = 0.0_RP
       LAND_SFC_albedo_t(i,j,I_SW) = 0.0_RP
    enddo
    enddo

    !########## History & Monitor ##########
    call PROF_rapstart('LND_History', 1)
    call LAND_vars_history
    call PROF_rapend  ('LND_History', 1)

    return
  end subroutine LAND_driver

  !-----------------------------------------------------------------------------
  !> Get surface boundary from other model
  subroutine LAND_SURFACE_GET
    use mod_land_admin, only: &
       LAND_sw
    use mod_land_vars, only: &
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
       ATMOS_cosSZA,   &
       ATMOS_SFLX_prec
    use mod_cpl_vars, only: &
       CPL_getATM_LND
    implicit none
    !---------------------------------------------------------------------------

    if ( LAND_sw ) then
       call CPL_getATM_LND( ATMOS_TEMP     (:,:), & ! [OUT]
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
                            ATMOS_cosSZA   (:,:), & ! [OUT]
                            ATMOS_SFLX_prec(:,:)  ) ! [OUT]
    endif

    return
  end subroutine LAND_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine LAND_SURFACE_SET( countup )
    use mod_land_admin, only: &
       LAND_sw
    use mod_land_vars, only: &
       LAND_PROPERTY,   &
       I_Z0M,           &
       I_Z0H,           &
       I_Z0E,           &
       LAND_SFC_TEMP,   &
       LAND_SFC_albedo, &
       LAND_SFLX_MW,    &
       LAND_SFLX_MU,    &
       LAND_SFLX_MV,    &
       LAND_SFLX_SH,    &
       LAND_SFLX_LH,    &
       LAND_SFLX_GH,    &
       LAND_SFLX_evap,  &
       LAND_U10,        &
       LAND_V10,        &
       LAND_T2,         &
       LAND_Q2
    use mod_cpl_vars, only: &
       CPL_putLND
    implicit none

    ! arguments
    logical, intent(in) :: countup
    !---------------------------------------------------------------------------

    if ( LAND_sw ) then
       call CPL_putLND( LAND_SFC_TEMP  (:,:),       & ! [IN]
                        LAND_SFC_albedo(:,:,:),     & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0M), & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0H), & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0E), & ! [IN]
                        LAND_SFLX_MW   (:,:),       & ! [IN]
                        LAND_SFLX_MU   (:,:),       & ! [IN]
                        LAND_SFLX_MV   (:,:),       & ! [IN]
                        LAND_SFLX_SH   (:,:),       & ! [IN]
                        LAND_SFLX_LH   (:,:),       & ! [IN]
                        LAND_SFLX_GH   (:,:),       & ! [IN]
                        LAND_SFLX_evap (:,:),       & ! [IN]
                        LAND_U10       (:,:),       & ! [IN]
                        LAND_V10       (:,:),       & ! [IN]
                        LAND_T2        (:,:),       & ! [IN]
                        LAND_Q2        (:,:),       & ! [IN]
                        countup                     ) ! [IN]
    endif

    return
  end subroutine LAND_SURFACE_SET

end module mod_land_driver
