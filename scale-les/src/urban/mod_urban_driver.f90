!-------------------------------------------------------------------------------
!> module URBAN driver
!!
!! @par Description
!!          Urban module driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_urban_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_driver_setup
  public :: URBAN_driver
  public :: URBAN_SURFACE_GET
  public :: URBAN_SURFACE_SET

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
  subroutine URBAN_driver_setup
    use mod_urban_phy_driver, only: &
       URBAN_PHY_driver_setup
    use mod_urban_vars, only: &
       URBAN_vars_history
    use mod_urban_admin, only: &
       URBAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[URBAN] / Origin[SCALE-LES]'

    !########## Get Surface Boundary from coupler ##########
    call URBAN_SURFACE_GET

    call URBAN_PHY_driver_setup

    !########## Set Surface Boundary to coupler ##########
    call URBAN_SURFACE_SET( countup=.false. )

    !########## History & Monitor ##########
    if ( URBAN_sw ) then
       call PROF_rapstart('URB History', 1)
       call URBAN_vars_history
       call PROF_rapend  ('URB History', 1)
    endif

    return
  end subroutine URBAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Urban step
  subroutine URBAN_driver
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN
    use mod_urban_admin, only: &
       URBAN_sw
    use mod_urban_vars, only: &
       URBAN_TR_t,        &
       URBAN_TB_t,        &
       URBAN_TG_t,        &
       URBAN_TC_t,        &
       URBAN_QC_t,        &
       URBAN_UC_t,        &
       URBAN_TRL_t,       &
       URBAN_TBL_t,       &
       URBAN_TGL_t,       &
       URBAN_RAINR_t,     &
       URBAN_RAINB_t,     &
       URBAN_RAING_t,     &
       URBAN_ROFF_t,      &
       URBAN_TR,          &
       URBAN_TB,          &
       URBAN_TG,          &
       URBAN_TC,          &
       URBAN_QC,          &
       URBAN_UC,          &
       URBAN_TRL,         &
       URBAN_TBL,         &
       URBAN_TGL,         &
       URBAN_RAINR,       &
       URBAN_RAINB,       &
       URBAN_RAING,       &
       URBAN_ROFF,        &
       URBAN_vars_total,  &
       URBAN_vars_history
    use mod_urban_phy_driver, only: &
       URBAN_PHY_driver
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call URBAN_SURFACE_GET

    !########## Physics ##########
    if ( URBAN_sw ) then
       call PROF_rapstart('URB Physics', 1)
       call URBAN_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('URB Physics', 1)
    endif

    !########## Update ##########
    do j = JS, JE
    do i = IS, IE
    do k = UKS, UKE
       URBAN_TRL(k,i,j) = URBAN_TRL(k,i,j) + URBAN_TRL_t(k,i,j) * dt
       URBAN_TBL(k,i,j) = URBAN_TBL(k,i,j) + URBAN_TBL_t(k,i,j) * dt
       URBAN_TGL(k,i,j) = URBAN_TGL(k,i,j) + URBAN_TGL_t(k,i,j) * dt
    end do
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       URBAN_TR(i,j) = URBAN_TR(i,j) + URBAN_TR_t(i,j) * dt
       URBAN_TB(i,j) = URBAN_TB(i,j) + URBAN_TB_t(i,j) * dt
       URBAN_TG(i,j) = URBAN_TG(i,j) + URBAN_TG_t(i,j) * dt
       URBAN_TC(i,j) = URBAN_TC(i,j) + URBAN_TC_t(i,j) * dt
       URBAN_QC(i,j) = URBAN_QC(i,j) + URBAN_QC_t(i,j) * dt
       URBAN_UC(i,j) = URBAN_UC(i,j) + URBAN_UC_t(i,j) * dt

       URBAN_RAINR(i,j) = URBAN_RAINR(i,j) + URBAN_RAINR_t(i,j) * dt
       URBAN_RAINB(i,j) = URBAN_RAINB(i,j) + URBAN_RAINB_t(i,j) * dt
       URBAN_RAING(i,j) = URBAN_RAING(i,j) + URBAN_RAING_t(i,j) * dt
       URBAN_ROFF (i,j) = URBAN_ROFF (i,j) + URBAN_ROFF_t (i,j) * dt
    end do
    end do

    call URBAN_vars_total

    !########## Set Surface Boundary to coupler ##########
    call URBAN_SURFACE_SET( countup=.true. )

    !########## reset tendencies ##########
    do j = JS, JE
    do i = IS, IE
    do k = UKS, UKE
       URBAN_TRL_t(k,i,j) = 0.0_RP
       URBAN_TBL_t(k,i,j) = 0.0_RP
       URBAN_TGL_t(k,i,j) = 0.0_RP
    end do
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       URBAN_TR_t(i,j) = 0.0_RP
       URBAN_TB_t(i,j) = 0.0_RP
       URBAN_TG_t(i,j) = 0.0_RP
       URBAN_TC_t(i,j) = 0.0_RP
       URBAN_QC_t(i,j) = 0.0_RP
       URBAN_UC_t(i,j) = 0.0_RP

       URBAN_RAINR_t(i,j) = 0.0_RP
       URBAN_RAINB_t(i,j) = 0.0_RP
       URBAN_RAING_t(i,j) = 0.0_RP
       URBAN_ROFF_t (i,j) = 0.0_RP
    enddo
    enddo

    !########## History & Monitor ##########
    call PROF_rapstart('URB History', 1)
    call URBAN_vars_history
    call PROF_rapend  ('URB History', 1)

    return
  end subroutine URBAN_driver

  !-----------------------------------------------------------------------------
  !> Get surface boundary
  subroutine URBAN_SURFACE_GET
    use mod_urban_admin, only: &
       URBAN_sw
    use mod_urban_vars, only: &
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
       CPL_getATM_URB
    implicit none
    !---------------------------------------------------------------------------

    if ( URBAN_sw ) then
       call CPL_getATM_URB( ATMOS_TEMP     (:,:), & ! [OUT]
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
  end subroutine URBAN_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Set surface boundary to other model
  subroutine URBAN_SURFACE_SET( countup )
    use mod_urban_admin, only: &
       URBAN_sw
    use mod_urban_vars, only: &
       URBAN_SFC_TEMP,   &
       URBAN_SFC_albedo, &
       URBAN_SFLX_MW,    &
       URBAN_SFLX_MU,    &
       URBAN_SFLX_MV,    &
       URBAN_SFLX_SH,    &
       URBAN_SFLX_LH,    &
       URBAN_SFLX_GH,    &
       URBAN_SFLX_evap,  &
       URBAN_Z0M,        &
       URBAN_Z0H,        &
       URBAN_Z0E,        &
       URBAN_U10,        &
       URBAN_V10,        &
       URBAN_T2,         &
       URBAN_Q2
    use mod_cpl_vars, only: &
       CPL_putURB
    implicit none

    ! arguments
    logical, intent(in) :: countup
    !---------------------------------------------------------------------------

    if ( URBAN_sw ) then
       call CPL_putURB( URBAN_SFC_TEMP  (:,:),   & ! [IN]
                        URBAN_SFC_albedo(:,:,:), & ! [IN]
                        URBAN_Z0M       (:,:),   & ! [IN]
                        URBAN_Z0H       (:,:),   & ! [IN]
                        URBAN_Z0E       (:,:),   & ! [IN]
                        URBAN_SFLX_MW   (:,:),   & ! [IN]
                        URBAN_SFLX_MU   (:,:),   & ! [IN]
                        URBAN_SFLX_MV   (:,:),   & ! [IN]
                        URBAN_SFLX_SH   (:,:),   & ! [IN]
                        URBAN_SFLX_LH   (:,:),   & ! [IN]
                        URBAN_SFLX_GH   (:,:),   & ! [IN]
                        URBAN_SFLX_evap (:,:),   & ! [IN]
                        URBAN_U10       (:,:),   & ! [IN]
                        URBAN_V10       (:,:),   & ! [IN]
                        URBAN_T2        (:,:),   & ! [IN]
                        URBAN_Q2        (:,:),   & ! [IN]
                        countup                  ) ! [IN]
    endif

    return
  end subroutine URBAN_SURFACE_SET

end module mod_urban_driver
