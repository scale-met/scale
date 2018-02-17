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
  use scale_ocean_grid_cartesC_index

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
  public :: OCEAN_driver_resume
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
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN] / Origin[SCALE-RM]'

    call OCEAN_PHY_driver_setup

!    if( OCEAN_FRC_sw ) call OCEAN_FRC_driver_setup

    return
  end subroutine OCEAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine OCEAN_driver_resume
    use mod_ocean_phy_driver, only: &
       OCEAN_PHY_driver_resume
!    use mod_ocean_frc_nudge, only: &
!       OCEAN_FRC_driver_resume
    use mod_ocean_vars, only: &
       OCEAN_vars_history, &
       OCEAN_TEMP_t,       &
       OCEAN_SALT_t,       &
       OCEAN_UVEL_t,       &
       OCEAN_VVEL_t,       &
       OCEAN_SFC_TEMP_t,   &
       OCEAN_SFC_albedo_t, &
       OCEAN_SFC_Z0M_t,    &
       OCEAN_SFC_Z0H_t,    &
       OCEAN_SFC_Z0E_t
    use mod_ocean_admin, only: &
       OCEAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN] / Origin[SCALE-RM]'

    !########## Get Surface Boundary from coupler ##########
    call OCEAN_SURFACE_GET

    !########## initialize tendencies ##########
!OCL XFILL
    OCEAN_SFC_TEMP_t  (:,:)   = 0.0_RP
!OCL XFILL
    OCEAN_SFC_albedo_t(:,:,:) = 0.0_RP
!OCL XFILL
    OCEAN_SFC_Z0M_t   (:,:)   = 0.0_RP
!OCL XFILL
    OCEAN_SFC_Z0H_t   (:,:)   = 0.0_RP
!OCL XFILL
    OCEAN_SFC_Z0E_t   (:,:)   = 0.0_RP
!OCL XFILL
    OCEAN_TEMP_t      (:,:,:) = 0.0_RP
!OCL XFILL
    OCEAN_SALT_t      (:,:,:) = 0.0_RP
!OCL XFILL
    OCEAN_UVEL_t      (:,:,:) = 0.0_RP
!OCL XFILL
    OCEAN_VVEL_t      (:,:,:) = 0.0_RP

    ! setup each components
    call OCEAN_PHY_driver_resume
!    if( OCEAN_FRC_sw ) call OCEAN_FRC_driver_resume

    !########## Set Surface Boundary to coupler ##########
    call OCEAN_SURFACE_SET( countup=.true. )

    !########## History & Monitor ##########
    if ( OCEAN_sw ) then
       call PROF_rapstart('OCN_History', 1)
       call OCEAN_vars_history
       call PROF_rapend  ('OCN_History', 1)
    endif

    return
  end subroutine OCEAN_driver_resume

  !-----------------------------------------------------------------------------
  !> Ocean step
  subroutine OCEAN_driver
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_ocean_vars, only: &
       OCEAN_TEMP,         &
       OCEAN_SALT,         &
       OCEAN_UVEL,         &
       OCEAN_VVEL,         &
       OCEAN_SFC_TEMP,     &
       OCEAN_SFC_albedo,   &
       OCEAN_SFC_Z0M,      &
       OCEAN_SFC_Z0H,      &
       OCEAN_SFC_Z0E,      &
       OCEAN_TEMP_t,       &
       OCEAN_SALT_t,       &
       OCEAN_UVEL_t,       &
       OCEAN_VVEL_t,       &
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

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call PROF_rapstart('OCN_SfcExch', 2)
    call OCEAN_SURFACE_GET
    call PROF_rapend  ('OCN_SfcExch', 2)

    !########## Physics ##########
    if ( OCEAN_sw ) then
       call PROF_rapstart('OCN_Physics', 1)
       call OCEAN_PHY_driver( update_flag = .true. )
       call PROF_rapend  ('OCN_Physics', 1)
    endif

    !########## Forcing ##########
!    if ( OCEAN_FORCE_sw ) then
!       call PROF_rapstart('OCN_Forcing', 1)
!       call OCEAN_forcing
!       call PROF_rapend  ('OCN_Forcing', 1)
!    endif

    !########## Update ##########
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_SFC_TEMP  (i,j)      = OCEAN_SFC_TEMP  (i,j)      + OCEAN_SFC_TEMP_t  (i,j)      * dt
       OCEAN_SFC_albedo(i,j,I_LW) = OCEAN_SFC_albedo(i,j,I_LW) + OCEAN_SFC_albedo_t(i,j,I_LW) * dt
       OCEAN_SFC_albedo(i,j,I_SW) = OCEAN_SFC_albedo(i,j,I_SW) + OCEAN_SFC_albedo_t(i,j,I_SW) * dt
       OCEAN_SFC_Z0M   (i,j)      = OCEAN_SFC_Z0M   (i,j)      + OCEAN_SFC_Z0M_t   (i,j)      * dt
       OCEAN_SFC_Z0H   (i,j)      = OCEAN_SFC_Z0H   (i,j)      + OCEAN_SFC_Z0H_t   (i,j)      * dt
       OCEAN_SFC_Z0E   (i,j)      = OCEAN_SFC_Z0E   (i,j)      + OCEAN_SFC_Z0E_t   (i,j)      * dt
    enddo
    enddo

    do j = OJS, OJE
    do i = OIS, OIE
    do k = OKS, OKE
       OCEAN_TEMP(k,i,j) = OCEAN_TEMP(k,i,j) + OCEAN_TEMP_t(k,i,j) * dt
       OCEAN_SALT(k,i,j) = OCEAN_SALT(k,i,j) + OCEAN_SALT_t(k,i,j) * dt
       OCEAN_UVEL(k,i,j) = OCEAN_UVEL(k,i,j) + OCEAN_UVEL_t(k,i,j) * dt
       OCEAN_VVEL(k,i,j) = OCEAN_VVEL(k,i,j) + OCEAN_VVEL_t(k,i,j) * dt
    enddo
    enddo
    enddo

    call OCEAN_vars_total

    !########## Set Surface Boundary to coupler ##########
    call PROF_rapstart('OCN_SfcExch', 2)
    call OCEAN_SURFACE_SET( countup=.true. )
    call PROF_rapend  ('OCN_SfcExch', 2)

    !########## reset tendencies ##########
!OCL XFILL
    do j = OJS, OJE
    do i = OIS, OIE
       OCEAN_SFC_TEMP_t  (i,j)      = 0.0_RP
       OCEAN_SFC_albedo_t(i,j,I_LW) = 0.0_RP
       OCEAN_SFC_albedo_t(i,j,I_SW) = 0.0_RP
       OCEAN_SFC_Z0M_t   (i,j)      = 0.0_RP
       OCEAN_SFC_Z0H_t   (i,j)      = 0.0_RP
       OCEAN_SFC_Z0E_t   (i,j)      = 0.0_RP
    enddo
    enddo
!OCL XFILL
    do j = OJS, OJE
    do i = OIS, OIE
    do k = OKS, OKE
       OCEAN_TEMP_t(k,i,j) = 0.0_RP
       OCEAN_SALT_t(k,i,j) = 0.0_RP
       OCEAN_UVEL_t(k,i,j) = 0.0_RP
       OCEAN_VVEL_t(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    !########## History & Monitor ##########
    call PROF_rapstart('OCN_History', 1)
    call OCEAN_vars_history
    call PROF_rapend  ('OCN_History', 1)

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
       ATMOS_SFC_DENS, &
       ATMOS_SFC_PRES, &
       ATMOS_SFLX_LW,  &
       ATMOS_SFLX_SW,  &
       ATMOS_cosSZA,   &
       ATMOS_SFLX_prec
    use mod_cpl_vars, only: &
       CPL_getATM_OCN
    implicit none

    real(RP) :: ATMOS_SFLX_rad_dn(OIA,OJA,2,2)
    real(RP) :: ATMOS_SFLX_rain  (OIA,OJA)
    real(RP) :: ATMOS_SFLX_snow  (OIA,OJA)

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( OCEAN_sw ) then
       call CPL_getATM_OCN( ATMOS_TEMP       (:,:),     & ! [OUT]
                            ATMOS_PRES       (:,:),     & ! [OUT]
                            ATMOS_W          (:,:),     & ! [OUT]
                            ATMOS_U          (:,:),     & ! [OUT]
                            ATMOS_V          (:,:),     & ! [OUT]
                            ATMOS_DENS       (:,:),     & ! [OUT]
                            ATMOS_QV         (:,:),     & ! [OUT]
                            ATMOS_PBL        (:,:),     & ! [OUT]
                            ATMOS_SFC_DENS   (:,:),     & ! [OUT]
                            ATMOS_SFC_PRES   (:,:),     & ! [OUT]
                            ATMOS_SFLX_rad_dn(:,:,:,:), & ! [OUT]
                            ATMOS_cosSZA     (:,:),     & ! [OUT]
                            ATMOS_SFLX_rain  (:,:),     & ! [OUT]
                            ATMOS_SFLX_snow  (:,:)      ) ! [OUT]
    endif

!OCL XFILL
    do j = OJS, OJE
    do i = OIS, OIE
       ATMOS_SFLX_SW  (i,j) = ATMOS_SFLX_rad_dn(i,j,I_SW,1) + ATMOS_SFLX_rad_dn(i,j,I_SW,2) ! direct+diffuse
       ATMOS_SFLX_LW  (i,j) = ATMOS_SFLX_rad_dn(i,j,I_LW,1) + ATMOS_SFLX_rad_dn(i,j,I_LW,2) ! direct+diffuse

       ATMOS_SFLX_prec(i,j) = ATMOS_SFLX_rain(i,j) + ATMOS_SFLX_snow(i,j) ! liquid+ice
    enddo
    enddo

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
