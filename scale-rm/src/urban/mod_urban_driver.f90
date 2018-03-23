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
  use scale_urban_grid_cartesC_index

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
  public :: URBAN_driver_setup
  public :: URBAN_driver_resume
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
    use scale_prc, only: &
       PRC_abort
    use mod_urban_phy_driver, only: &
       URBAN_PHY_driver_setup
    use mod_urban_admin, only: &
       URBAN_do, &
       URBAN_DYN_TYPE
    use scale_urban_dyn_kusaka01, only: &
       URBAN_DYN_KUSAKA01_setup
    use mod_urban_vars, only: &
       URBAN_Z0M, &
       URBAN_Z0H, &
       URBAN_Z0E
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER Setup] / Categ[URBAN] / Origin[SCALE-RM]'

    if ( URBAN_do ) then

       select case ( URBAN_DYN_TYPE )
       case ( 'KUSAKA01' )
          call URBAN_DYN_KUSAKA01_setup( UIA, UIS, UIE, UJA, UJS, UJE, &
                                         URBAN_Z0M(:,:), URBAN_Z0H(:,:), URBAN_Z0E(:,:) ) ! [OUT]

       case ( 'CONST' )
          ! do nothing
       case default
          write(*,*) 'xxx LAND_DYN_TYPE is invalid: ', trim(URBAN_DYN_TYPE)
          call PRC_abort
       end select

       call URBAN_PHY_driver_setup

    end if

    return
  end subroutine URBAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine URBAN_driver_resume
    use mod_urban_phy_driver, only: &
       URBAN_PHY_driver_resume
    use mod_urban_vars, only: &
       URBAN_vars_history, &
       URBAN_TRL_t,        &
       URBAN_TBL_t,        &
       URBAN_TGL_t,        &
       URBAN_TR_t,         &
       URBAN_TB_t,         &
       URBAN_TG_t,         &
       URBAN_TC_t,         &
       URBAN_QC_t,         &
       URBAN_UC_t,         &
       URBAN_RAINR_t,      &
       URBAN_RAINB_t,      &
       URBAN_RAING_t,      &
       URBAN_ROFF_t
    use mod_admin_restart, only: &
       RESTART_RUN
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[URBAN] / Origin[SCALE-RM]'

    !########## Get Surface Boundary from coupler ##########
    call URBAN_SURFACE_GET

    !########## initialize tendency ##########
!OCL XFILL
    do j = UJS, UJE
    do i = UIS, UIE
    do k = UKS, UKE
       URBAN_TRL_t(k,i,j) = 0.0_RP
       URBAN_TBL_t(k,i,j) = 0.0_RP
       URBAN_TGL_t(k,i,j) = 0.0_RP
    end do
    end do
    end do

!OCL XFILL
    do j = UJS, UJE
    do i = UIS, UIE
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

    ! resume each component
    if ( .NOT. RESTART_RUN ) then ! tentative
       call URBAN_PHY_driver_resume
       call URBAN_driver
    end if

    !########## Set Surface Boundary to coupler ##########
    call URBAN_SURFACE_SET( countup=.true. )

!!$    !########## History & Monitor ##########
!!$    call PROF_rapstart('URB_History', 1)
!!$    call URBAN_vars_history
!!$    call PROF_rapend  ('URB_History', 1)

    return
  end subroutine URBAN_driver_resume

  !-----------------------------------------------------------------------------
  !> Urban step
  subroutine URBAN_driver
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use mod_urban_vars, only: &
       ATMOS_TEMP,      &
       ATMOS_PRES,      &
       ATMOS_W,         &
       ATMOS_U,         &
       ATMOS_V,         &
       ATMOS_DENS,      &
       ATMOS_QV,        &
       ATMOS_PBL,       &
       ATMOS_SFC_DENS,  &
       ATMOS_SFC_PRES,  &
       ATMOS_SFLX_LW,   &
       ATMOS_SFLX_SW,   &
       ATMOS_cosSZA,    &
       ATMOS_SFLX_rain, &
       ATMOS_SFLX_snow, &
       URBAN_TRL_t,       &
       URBAN_TBL_t,       &
       URBAN_TGL_t,       &
       URBAN_TR_t,        &
       URBAN_TB_t,        &
       URBAN_TG_t,        &
       URBAN_TC_t,        &
       URBAN_QC_t,        &
       URBAN_UC_t,        &
       URBAN_RAINR_t,     &
       URBAN_RAINB_t,     &
       URBAN_RAING_t,     &
       URBAN_ROFF_t,      &
       URBAN_SFC_TEMP,    &
       URBAN_SFC_albedo,  &
       URBAN_SFLX_MW,     &
       URBAN_SFLX_MU,     &
       URBAN_SFLX_MV,     &
       URBAN_SFLX_SH,     &
       URBAN_SFLX_LH,     &
       URBAN_SFLX_evap,   &
       URBAN_SFLX_GH,     &
       URBAN_Z0M,         &
       URBAN_Z0H,         &
       URBAN_Z0E,         &
       URBAN_U10,         &
       URBAN_V10,         &
       URBAN_T2,          &
       URBAN_Q2,          &
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
    use scale_urban_dyn_kusaka01, only: &
       URBAN_DYN_kusaka01
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_urban_grid_cartesC, only: &
       CDZ => URBAN_GRID_CARTESC_CDZ
    use scale_landuse, only: &
       LANDUSE_fact_urban
    use mod_urban_admin, only: &
       URBAN_DYN_TYPE
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE
    use scale_mapprojection, only: &
       BASE_LON => MAPPROJECTION_basepoint_lon, &
       BASE_LAT => MAPPROJECTION_basepoint_lat
    implicit none

    real(RP) :: LHV(UIA,UJA) ! latent heat of vaporization [J/kg]

    real(RP) :: LAT, LON
    integer  :: tloc     ! local time (1-24h)
    real(RP) :: dsec     ! second [s]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !########## Get Surface Boundary from coupler ##########
    call PROF_rapstart('URB_SfcExch', 2)
    call URBAN_SURFACE_GET
    call PROF_rapend  ('URB_SfcExch', 2)

    !########## Dynamics / Update variables ##########
    select case ( URBAN_DYN_TYPE )
    case ( 'CONST' )
       ! do nothing
    case ( 'KUSAKA01' )

       LAT = BASE_LAT / D2R
       LON = BASE_LON / D2R
       tloc = mod( (NOWDATE(4) + int(LON/15.0_RP)),24 )
       dsec = real( NOWDATE(5)*60.0_RP + NOWDATE(6), kind=RP ) / 3600.0_RP
       if( tloc == 0 ) tloc = 24

       call HYDROMETEOR_LHV( LHV(:,:), ATMOS_TEMP(:,:) )

       call URBAN_DYN_kusaka01( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                                ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                            & ! [IN]
                                ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                    & ! [IN]
                                ATMOS_DENS(:,:), ATMOS_QV(:,:), LHV(:,:),                    & ! [IN]
                                REAL_Z1(:,:), ATMOS_PBL(:,:),                                & ! [IN]
                                ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                    & ! [IN]
                                ATMOS_SFLX_LW(:,:,:), ATMOS_SFLX_SW(:,:,:),                  & ! [IN]
                                ATMOS_SFLX_rain(:,:), ATMOS_SFLX_snow(:,:),                  & ! [IN]
                                URBAN_TR(:,:), URBAN_TB(:,:), URBAN_TG(:,:),                 & ! [INOUT]
                                URBAN_TC(:,:), URBAN_QC(:,:), URBAN_UC(:,:),                 & ! [INOUT]
                                URBAN_TRL(:,:,:), URBAN_TBL(:,:,:), URBAN_TGL(:,:,:),        & ! [INOUT]
                                URBAN_RAINR(:,:), URBAN_RAINB(:,:), URBAN_RAING(:,:),        & ! [INOUT]
                                URBAN_ROFF(:,:),                                             & ! [INOUT]
                                CDZ(:),                                                      & ! [IN]
                                LANDUSE_fact_urban(:,:),                                     & ! [IN]
                                tloc, dsec, dt,                                              & ! [IN]
                                URBAN_SFC_TEMP(:,:),                                         & ! [OUT]
                                URBAN_SFC_albedo(:,:,I_LW), URBAN_SFC_albedo(:,:,I_SW),      & ! [OUT]
                                URBAN_SFLX_MW(:,:), URBAN_SFLX_MU(:,:), URBAN_SFLX_MV(:,:),  & ! [OUT]
                                URBAN_SFLX_SH(:,:), URBAN_SFLX_LH(:,:), URBAN_SFLX_GH(:,:),  & ! [OUT]
                                URBAN_Z0M(:,:), URBAN_Z0H(:,:), URBAN_Z0E(:,:),              & ! [OUT]
                                URBAN_U10(:,:), URBAN_V10(:,:), URBAN_T2(:,:), URBAN_Q2(:,:) ) ! [OUT]

!OCL XFILL
       do j = UJS, UJE
       do i = UIS, UIE
          URBAN_SFLX_evap(i,j) = URBAN_SFLX_LH(i,j) / LHV(i,j)
       end do
       end do

    end select

!!$    do j = UJS, UJE
!!$    do i = UIS, UIE
!!$       URBAN_TR(i,j) = URBAN_TR(i,j) + URBAN_TR_t(i,j) * dt
!!$       URBAN_TB(i,j) = URBAN_TB(i,j) + URBAN_TB_t(i,j) * dt
!!$       URBAN_TG(i,j) = URBAN_TG(i,j) + URBAN_TG_t(i,j) * dt
!!$       URBAN_TC(i,j) = URBAN_TC(i,j) + URBAN_TC_t(i,j) * dt
!!$       URBAN_QC(i,j) = URBAN_QC(i,j) + URBAN_QC_t(i,j) * dt
!!$       URBAN_UC(i,j) = URBAN_UC(i,j) + URBAN_UC_t(i,j) * dt
!!$
!!$       URBAN_RAINR(i,j) = URBAN_RAINR(i,j) + URBAN_RAINR_t(i,j) * dt
!!$       URBAN_RAINB(i,j) = URBAN_RAINB(i,j) + URBAN_RAINB_t(i,j) * dt
!!$       URBAN_RAING(i,j) = URBAN_RAING(i,j) + URBAN_RAING_t(i,j) * dt
!!$       URBAN_ROFF (i,j) = URBAN_ROFF (i,j) + URBAN_ROFF_t (i,j) * dt
!!$    end do
!!$    end do

    call URBAN_vars_total

    !########## Set Surface Boundary to coupler ##########
    call PROF_rapstart('URB_SfcExch', 2)
    call URBAN_SURFACE_SET( countup=.true. )
    call PROF_rapend  ('URB_SfcExch', 2)

    !########## reset tendencies ##########
!OCL XFILL
    do j = UJS, UJE
    do i = UIS, UIE
    do k = UKS, UKE
       URBAN_TRL_t(k,i,j) = 0.0_RP
       URBAN_TBL_t(k,i,j) = 0.0_RP
       URBAN_TGL_t(k,i,j) = 0.0_RP
    end do
    end do
    end do

!OCL XFILL
    do j = UJS, UJE
    do i = UIS, UIE
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

    !########## Physics ##########
    call PROF_rapstart('URB_Physics', 1)
    call URBAN_PHY_driver( update_flag = .true. )
    call PROF_rapend  ('URB_Physics', 1)

    !########## History & Monitor ##########
    call PROF_rapstart('URB_History', 1)
    call URBAN_vars_history
    call PROF_rapend  ('URB_History', 1)

    return
  end subroutine URBAN_driver

  !-----------------------------------------------------------------------------
  !> Get surface boundary
  subroutine URBAN_SURFACE_GET
    use mod_urban_admin, only: &
       URBAN_do
    use mod_urban_vars, only: &
       ATMOS_TEMP,      &
       ATMOS_PRES,      &
       ATMOS_W,         &
       ATMOS_U,         &
       ATMOS_V,         &
       ATMOS_DENS,      &
       ATMOS_QV,        &
       ATMOS_PBL,       &
       ATMOS_SFC_DENS,  &
       ATMOS_SFC_PRES,  &
       ATMOS_SFLX_LW,   &
       ATMOS_SFLX_SW,   &
       ATMOS_cosSZA,    &
       ATMOS_SFLX_rain, &
       ATMOS_SFLX_snow
    use mod_cpl_vars, only: &
       CPL_getATM_URB
    implicit none

    real(RP) :: ATMOS_SFLX_rad_dn(UIA,UJA,2,2)

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( URBAN_do ) then
       call CPL_getATM_URB( ATMOS_TEMP       (:,:),     & ! [OUT]
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
    do j = UJS, UJE
    do i = UIS, UIE
       ATMOS_SFLX_SW  (i,j,1) = ATMOS_SFLX_rad_dn(i,j,I_SW,1) ! direct
       ATMOS_SFLX_LW  (i,j,1) = ATMOS_SFLX_rad_dn(i,j,I_LW,1) ! direct
       ATMOS_SFLX_SW  (i,j,2) = ATMOS_SFLX_rad_dn(i,j,I_SW,2) ! diffuse
       ATMOS_SFLX_LW  (i,j,2) = ATMOS_SFLX_rad_dn(i,j,I_LW,2) ! diffuse
    enddo
    enddo

    return
  end subroutine URBAN_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Set surface boundary to other model
  subroutine URBAN_SURFACE_SET( countup )
    use mod_urban_admin, only: &
       URBAN_do
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

    if ( URBAN_do ) then
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
