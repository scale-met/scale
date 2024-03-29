!-------------------------------------------------------------------------------
!> module URBAN driver
!!
!! @par Description
!!          Urban module driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_urban_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_urban_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_driver_setup
  public :: URBAN_driver_finalize
  public :: URBAN_driver_calc_tendency
  public :: URBAN_driver_update
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
  real(RP), private, allocatable :: AH_URB (:,:,:) ! urban grid average of anthropogenic sensible heat [W/m2]
  real(RP), private, allocatable :: AHL_URB(:,:,:) ! urban grid average of anthropogenic latent heat [W/m2]
  real(RP), private              :: AH_TOFFSET     ! time offset for AH [Hour]
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use mod_urban_admin, only: &
       URBAN_do, &
       URBAN_DYN_TYPE, &
       URBAN_SFC_TYPE
    use scale_urban_dyn_kusaka01, only: &
       URBAN_DYN_KUSAKA01_setup
    use scale_landuse, only: &
       LANDUSE_fact_urban
    use mod_urban_vars, only: &
       URBAN_Z0M, &
       URBAN_Z0H, &
       URBAN_Z0E, &
       URBAN_ZD
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_driver_setup",*) 'Setup'

    if ( URBAN_do ) then

       allocate( AH_URB (UIA,UJA,1:24) )
       allocate( AHL_URB(UIA,UJA,1:24) )
       AH_URB  (:,:,:) = UNDEF
       AHL_URB (:,:,:) = UNDEF
       !$acc enter data create(AH_URB,AHL_URB)

       select case ( URBAN_DYN_TYPE )
       case ( 'KUSAKA01' )
          call URBAN_DYN_KUSAKA01_setup( UIA, UIS, UIE, UJA, UJS, UJE,                   & ! [IN]
                                         LANDUSE_fact_urban(:,:),                        & ! [IN]
                                         URBAN_Z0M(:,:), URBAN_Z0H(:,:), URBAN_Z0E(:,:), & ! [OUT]
                                         URBAN_ZD(:,:),                                  & ! [OUT]
                                         AH_URB(:,:,:), AHL_URB(:,:,:), AH_TOFFSET       ) ! [OUT]

          URBAN_SFC_TYPE = 'KUSAKA01'
       case default
          LOG_ERROR("URBAN_driver_setup",*) 'LAND_DYN_TYPE is invalid: ', trim(URBAN_DYN_TYPE)
          call PRC_abort
       end select

       select case ( URBAN_SFC_TYPE )
       case ( 'KUSAKA01' )
          ! do nothing
       case default
          LOG_ERROR("URBAN_driver_setup",*) 'LAND_SFC_TYPE is invalid: ', trim(URBAN_SFC_TYPE)
          call PRC_abort
       end select

    end if

    return
  end subroutine URBAN_driver_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine URBAN_driver_finalize
    use mod_urban_admin, only: &
       URBAN_do, &
       URBAN_DYN_TYPE, &
       URBAN_SFC_TYPE
    use scale_urban_dyn_kusaka01, only: &
       URBAN_DYN_kusaka01_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_driver_finalize",*) 'Finalize'

    if ( URBAN_do ) then

       select case ( URBAN_DYN_TYPE )
       case ( 'KUSAKA01' )
          call URBAN_DYN_kusaka01_finalize
       end select

       select case ( URBAN_SFC_TYPE )
       case ( 'KUSAKA01' )
       end select

       !$acc exit data delete(AH_URB,AHL_URB)
       deallocate( AH_URB  )
       deallocate( AHL_URB )

    end if

    return
  end subroutine URBAN_driver_finalize

  !-----------------------------------------------------------------------------
  !> Calclate tendency
  subroutine URBAN_driver_calc_tendency( force )
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_urban_grid_cartesC_real, only: &
       URBAN_GRID_CARTESC_REAL_VOL,    &
       URBAN_GRID_CARTESC_REAL_TOTVOL, &
       URBAN_GRID_CARTESC_REAL_AREA,   &
       URBAN_GRID_CARTESC_REAL_TOTAREA
    use scale_topography, only: &
       TanSL_X => TOPOGRAPHY_TanSL_X, &
       TanSL_Y => TOPOGRAPHY_TanSL_Y
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_ch
    use mod_atmos_phy_ch_driver, only: &
       ATMOS_PHY_CH_driver_URBAN_flux
    use mod_urban_vars, only: &
       ATMOS_TEMP,      &
       ATMOS_PRES,      &
       ATMOS_U,         &
       ATMOS_V,         &
       ATMOS_DENS,      &
       ATMOS_QV,        &
       ATMOS_SFC_DENS,  &
       ATMOS_SFC_PRES,  &
       ATMOS_SFLX_LW,   &
       ATMOS_SFLX_SW,   &
       ATMOS_cosSZA,    &
       ATMOS_SFLX_water, &
       ATMOS_SFLX_ENGI, &
       URBAN_TRL_t,     &
       URBAN_TBL_t,     &
       URBAN_TGL_t,     &
       URBAN_TR_t,      &
       URBAN_TB_t,      &
       URBAN_TG_t,      &
       URBAN_TC_t,      &
       URBAN_QC_t,      &
       URBAN_UC_t,      &
       URBAN_RAINR_t,   &
       URBAN_RAINB_t,   &
       URBAN_RAING_t,   &
       URBAN_SFC_TEMP,    &
       URBAN_SFC_albedo,  &
       URBAN_SFLX_MW,     &
       URBAN_SFLX_MU,     &
       URBAN_SFLX_MV,     &
       URBAN_SFLX_SH,     &
       URBAN_SFLX_LH,     &
       URBAN_SFLX_SHEX,   &
       URBAN_SFLX_LHEX,   &
       URBAN_SFLX_QVEX,   &
       URBAN_SFLX_QTRC,   &
       URBAN_SFLX_GH,     &
       URBAN_Z0M,         &
       URBAN_Z0H,         &
       URBAN_Z0E,         &
       URBAN_ZD,          &
       URBAN_AH,          &
       URBAN_AHL,         &
       URBAN_Ustar,       &
       URBAN_Tstar,       &
       URBAN_Qstar,       &
       URBAN_Wstar,       &
       URBAN_RLmo,        &
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
       URBAN_ROFF
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV, &
       ATMOS_HYDROMETEOR_dry,                    &
       I_QV
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN, &
       NOWDATE => TIME_NOWDATE
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_urban_grid_cartesC, only: &
       CDZ => URBAN_GRID_CARTESC_CDZ
    use scale_landuse, only: &
       LANDUSE_fact_urban, &
       exists_urban => LANDUSE_exists_urban
    use mod_urban_admin, only: &
       URBAN_SFC_TYPE
    use scale_urban_dyn_kusaka01, only: &
       URBAN_DYN_kusaka01
    implicit none
    logical, intent(in) :: force

    real(RP) :: TRL(UKA,UIA,UJA), TBL(UKA,UIA,UJA), TGL(UKA,UIA,UJA)
    real(RP) :: TR(UIA,UJA), TB(UIA,UJA), TG(UIA,UJA)
    real(RP) :: TC(UIA,UJA), QC(UIA,UJA), UC(UIA,UJA)
    real(RP) :: RAINR(UIA,UJA), RAINB(UIA,UJA), RAING(UIA,UJA)

    real(RP) :: LHV(UIA,UJA)        ! latent heat of vaporization [J/kg]

   ! real(RP) :: URBAN_SFLX_LHEX(UIA,UJA)

    integer  :: tloc, tloc_next     ! universal time (1-24h)
    real(RP) :: dsec                ! second [s]

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('URB_CalcTend', 1)

    !$acc data create(TRL,TBL,TGL,TR,TB,TG,TC,QC,UC,RAINR,RAINB,RAING,LHV)

    !########## Get Surface Boundary from coupler ##########
    call URBAN_SURFACE_GET

    !########## initialize tendency ##########
!OCL XFILL
    !$acc kernels
    do j = UJS, UJE
    do i = UIS, UIE
    do k = UKS, UKE
       URBAN_TRL_t(k,i,j) = 0.0_RP
       URBAN_TBL_t(k,i,j) = 0.0_RP
       URBAN_TGL_t(k,i,j) = 0.0_RP
    end do
    end do
    end do
    !$acc end kernels

!OCL XFILL
    !$acc kernels
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
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$omp parallel do
    !$acc kernels
    do iq = 1, QA
    do j  = UJS, UJE
    do i  = UIS, UIE
       URBAN_SFLX_QTRC(i,j,iq) = 0.0_RP
    enddo
    enddo
    enddo
    !$acc end kernels

    select case ( URBAN_SFC_TYPE )
    case ( 'KUSAKA01' )

!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
       do k = UKS, UKE
          TRL(k,i,j) = URBAN_TRL(k,i,j)
          TBL(k,i,j) = URBAN_TBL(k,i,j)
          TGL(k,i,j) = URBAN_TGL(k,i,j)
       end do
       end do
       end do
       !$acc end kernels

!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
          TR(i,j) = URBAN_TR(i,j)
          TB(i,j) = URBAN_TB(i,j)
          TG(i,j) = URBAN_TG(i,j)
          TC(i,j) = URBAN_TC(i,j)
          QC(i,j) = URBAN_QC(i,j)
          UC(i,j) = URBAN_UC(i,j)
          RAINR(i,j) = URBAN_RAINR(i,j)
          RAINB(i,j) = URBAN_RAINB(i,j)
          RAING(i,j) = URBAN_RAING(i,j)
       end do
       end do
       !$acc end kernels


       ! universal time
       dsec = real( NOWDATE(5)*60.0_RP + NOWDATE(6), kind=RP ) / 3600.0_RP  ! [hour]
       tloc = NOWDATE(4)                                                    ! [hour]
       tloc = modulo(tloc-1,24)+1
       if ( tloc == 24 ) then
         tloc_next = 1
       else
         tloc_next = tloc + 1
       end if
       !--- Calculate AH at UTC
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
       if ( exists_urban(i,j) ) then
          URBAN_AH(i,j)  = ( 1.0_RP-dsec ) * AH_URB(i,j, tloc) &
                         + (        dsec ) * AH_URB(i,j, tloc_next)
          URBAN_AHL(i,j) = ( 1.0_RP-dsec ) * AHL_URB(i,j, tloc) &
                         + (        dsec ) * AHL_URB(i,j, tloc_next)
       end if
       enddo
       enddo
       !$acc end kernels

       call HYDROMETEOR_LHV( UIA, UIS, UIE, UJA, UJS, UJE, &
                             ATMOS_TEMP(:,:), LHV(:,:) )

       call URBAN_DYN_kusaka01( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                                ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                            & ! [IN]
                                ATMOS_U(:,:), ATMOS_V(:,:),                                  & ! [IN]
                                ATMOS_DENS(:,:), ATMOS_QV(:,:), LHV(:,:),                    & ! [IN]
                                REAL_Z1(:,:),                                                & ! [IN]
                                ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                    & ! [IN]
                                ATMOS_SFLX_LW(:,:,:), ATMOS_SFLX_SW(:,:,:),                  & ! [IN]
                                ATMOS_SFLX_water(:,:), ATMOS_SFLX_ENGI(:,:),                 & ! [IN]
                                URBAN_Z0M(:,:), URBAN_Z0H(:,:), URBAN_Z0E(:,:),              & ! [IN]
                                URBAN_ZD(:,:),                                               & ! [IN]
                                CDZ(:),                                                      & ! [IN]
                                TanSL_X(:,:), TanSL_Y(:,:),                                  & ! [IN]
                                LANDUSE_fact_urban(:,:),                                     & ! [IN]
                                dt,                                                          & ! [IN]
                                TRL(:,:,:), TBL(:,:,:), TGL(:,:,:),                          & ! [INOUT]
                                TR(:,:), TB(:,:), TG(:,:), TC(:,:), QC(:,:), UC(:,:),        & ! [INOUT]
                                RAINR(:,:), RAINB(:,:), RAING(:,:), URBAN_ROFF(:,:),         & ! [INOUT]
                                URBAN_SFC_TEMP(:,:),                                         & ! [OUT]
                                URBAN_SFC_albedo(:,:,:,:),                                   & ! [OUT]
                                URBAN_SFLX_MW(:,:), URBAN_SFLX_MU(:,:), URBAN_SFLX_MV(:,:),  & ! [OUT]
                                URBAN_SFLX_SH(:,:), URBAN_SFLX_LH(:,:), URBAN_SFLX_GH(:,:),  & ! [OUT]
                                URBAN_Ustar(:,:), URBAN_Tstar(:,:), URBAN_Qstar(:,:),        & ! [OUT]
                                URBAN_Wstar(:,:),                                            & ! [OUT]
                                URBAN_RLmo(:,:),                                             & ! [OUT]
                                URBAN_U10(:,:), URBAN_V10(:,:), URBAN_T2(:,:), URBAN_Q2(:,:) ) ! [OUT]

       !-----------------------------------------------------------
       ! anthropogenic heat fluxes
       !-----------------------------------------------------------
       !$omp parallel do
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
       if ( exists_urban(i,j) ) then
          !URBAN_SFLX_SHEX(i,j) = URBAN_AH (i,j) / LANDUSE_fact_urban(i,j) ! Sensible heat flux [W/m2]
          !URBAN_SFLX_LHEX(i,j) = URBAN_AHL(i,j) / LANDUSE_fact_urban(i,j) ! Latent heat flux   [W/m2]
          !URBAN_SFLX_SH  (i,j) = URBAN_SFLX_SH(i,j) + URBAN_SFLX_SHEX(i,j)
          !URBAN_SFLX_LH  (i,j) = URBAN_SFLX_LH(i,j) + URBAN_SFLX_LHEX(i,j)
          URBAN_SFLX_SHEX(i,j) = URBAN_AH (i,j)     ! Sensible anthropogenic heat flux [W/m2]
          URBAN_SFLX_LHEX(i,j) = URBAN_AHL(i,j)     ! Latent anthropogenic heat flux [W/m2]
       end if
       end do
       end do
       !$acc end kernels

!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
       if ( exists_urban(i,j) ) then
          do k = UKS, UKE
             URBAN_TRL_t(k,i,j) = ( TRL(k,i,j) - URBAN_TRL(k,i,j) ) / dt
             URBAN_TBL_t(k,i,j) = ( TBL(k,i,j) - URBAN_TBL(k,i,j) ) / dt
             URBAN_TGL_t(k,i,j) = ( TGL(k,i,j) - URBAN_TGL(k,i,j) ) / dt
          end do
       end if
       end do
       end do
       !$acc end kernels

!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
       if ( exists_urban(i,j) ) then
          URBAN_TR_t(i,j) = ( TR(i,j) - URBAN_TR(i,j) ) / dt
          URBAN_TB_t(i,j) = ( TB(i,j) - URBAN_TB(i,j) ) / dt
          URBAN_TG_t(i,j) = ( TG(i,j) - URBAN_TG(i,j) ) / dt
          URBAN_TC_t(i,j) = ( TC(i,j) - URBAN_TC(i,j) ) / dt
          URBAN_QC_t(i,j) = ( QC(i,j) - URBAN_QC(i,j) ) / dt
          URBAN_UC_t(i,j) = ( UC(i,j) - URBAN_UC(i,j) ) / dt
          URBAN_RAINR_t(i,j) = ( RAINR(i,j) - URBAN_RAINR(i,j) ) / dt
          URBAN_RAINB_t(i,j) = ( RAINB(i,j) - URBAN_RAINB(i,j) ) / dt
          URBAN_RAING_t(i,j) = ( RAING(i,j) - URBAN_RAING(i,j) ) / dt
       end if
       end do
       end do
       !$acc end kernels

       if ( .NOT. ATMOS_HYDROMETEOR_dry ) then
          !$omp parallel do
          !$acc kernels
          do j = UJS, UJE
          do i = UIS, UIE
          if ( exists_urban(i,j) ) then
             URBAN_SFLX_QTRC(i,j,I_QV) = URBAN_SFLX_LH  (i,j) / LHV(i,j)
             URBAN_SFLX_QVEX(i,j)      = URBAN_SFLX_LHEX(i,j) / LHV(i,j)
          end if
          enddo
          enddo
          !$acc end kernels
       endif

    end select

    ! Surface flux for chemical tracers
    if ( ATMOS_sw_phy_ch ) then
       call ATMOS_PHY_CH_driver_URBAN_flux( URBAN_SFLX_QTRC(:,:,:) ) ! [INOUT]
    endif

    call FILE_HISTORY_in( URBAN_TR_t(:,:), 'URBAN_TR_t', 'tendency of URBAN_TR', 'K/s',     dim_type='XY' )
    call FILE_HISTORY_in( URBAN_TB_t(:,:), 'URBAN_TB_t', 'tendency of URBAN_TB', 'K/s',     dim_type='XY' )
    call FILE_HISTORY_in( URBAN_TG_t(:,:), 'URBAN_TG_t', 'tendency of URBAN_TG', 'K/s',     dim_type='XY' )
    call FILE_HISTORY_in( URBAN_TC_t(:,:), 'URBAN_TC_t', 'tendency of URBAN_TC', 'K/s',     dim_type='XY' )
    call FILE_HISTORY_in( URBAN_QC_t(:,:), 'URBAN_QC_t', 'tendency of URBAN_QC', 'kg/kg/s', dim_type='XY' )
    call FILE_HISTORY_in( URBAN_UC_t(:,:), 'URBAN_UC_t', 'tendency of URBAN_UC', 'm/s2',    dim_type='XY' )

    call FILE_HISTORY_in( URBAN_TRL_t(:,:,:), 'URBAN_TRL_t', 'tendency of URBAN_TRL', 'K/s', dim_type='UXY' )
    call FILE_HISTORY_in( URBAN_TBL_t(:,:,:), 'URBAN_TBL_t', 'tendency of URBAN_TBL', 'K/s', dim_type='UXY' )
    call FILE_HISTORY_in( URBAN_TGL_t(:,:,:), 'URBAN_TGL_t', 'tendency of URBAN_TGL', 'K/s', dim_type='UXY' )

    call FILE_HISTORY_in( URBAN_RAINR_t(:,:), 'URBAN_RAINR_t', 'tendency of URBAN_RAINR', 'kg/m2/s', dim_type='XY' )
    call FILE_HISTORY_in( URBAN_RAINB_t(:,:), 'URBAN_RAINB_t', 'tendency of URBAN_RAINB', 'kg/m2/s', dim_type='XY' )
    call FILE_HISTORY_in( URBAN_RAING_t(:,:), 'URBAN_RAING_t', 'tendency of URBAN_RAING', 'kg/m2/s', dim_type='XY' )
    call FILE_HISTORY_in( URBAN_ROFF   (:,:), 'URBAN_ROFF',    'urban runoff water',      'kg/m2/s', dim_type='XY' )

    if ( STATISTICS_checktotal ) then

       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TRL_t (:,:,:), 'URBAN_TRL_t', &
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTVOL      )
       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TBL_t (:,:,:), 'URBAN_TBL_t', &
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTVOL      )
       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TGL_t (:,:,:), 'URBAN_TGL_t', &
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTVOL      )

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TR_t(:,:), 'URBAN_TR_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TB_t(:,:), 'URBAN_TB_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TG_t(:,:), 'URBAN_TG_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TC_t(:,:), 'URBAN_TC_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_QC_t(:,:), 'URBAN_QC_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_UC_t(:,:), 'URBAN_UC_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAINR_t(:,:), 'URBAN_RAINR_t', &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAINB_t(:,:), 'URBAN_RAINB_t', &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAING_t(:,:), 'URBAN_RAING_t', &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_ROFF(:,:),    'URBAN_ROFF',    &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
    endif


    !########## Set Surface Boundary to coupler ##########
    call URBAN_SURFACE_SET( countup=.true. )

    !$acc end data

    call PROF_rapend  ('URB_CalcTend', 1)

    return
  end subroutine URBAN_driver_calc_tendency

  !-----------------------------------------------------------------------------
  !> Urban step
  subroutine URBAN_driver_update
    use scale_time, only: &
       dt => TIME_DTSEC_URBAN
    use mod_urban_vars, only: &
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
       URBAN_vars_check
    use mod_urban_admin, only: &
       URBAN_DYN_TYPE
    use scale_landuse, only: &
       exists_urban => LANDUSE_exists_urban
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('URB_Update', 1)

    !########## Get Surface Boundary from coupler ##########
    call URBAN_SURFACE_GET

    !########## Dynamics / Update variables ##########
    select case ( URBAN_DYN_TYPE )
    case ( 'KUSAKA01' )

!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
       if ( exists_urban(i,j) ) then
          do k = UKS, UKE
             URBAN_TRL(k,i,j) = URBAN_TRL(k,i,j) + URBAN_TRL_t(k,i,j) * dt
             URBAN_TBL(k,i,j) = URBAN_TBL(k,i,j) + URBAN_TBL_t(k,i,j) * dt
             URBAN_TGL(k,i,j) = URBAN_TGL(k,i,j) + URBAN_TGL_t(k,i,j) * dt
          end do
       end if
       end do
       end do
       !$acc end kernels

!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = UJS, UJE
       do i = UIS, UIE
       if ( exists_urban(i,j) ) then
          URBAN_TR(i,j) = URBAN_TR(i,j) + URBAN_TR_t(i,j) * dt
          URBAN_TB(i,j) = URBAN_TB(i,j) + URBAN_TB_t(i,j) * dt
          URBAN_TG(i,j) = URBAN_TG(i,j) + URBAN_TG_t(i,j) * dt
          URBAN_TC(i,j) = URBAN_TC(i,j) + URBAN_TC_t(i,j) * dt
          URBAN_QC(i,j) = max( real(URBAN_QC(i,j) + URBAN_QC_t(i,j) * dt, RP), 0.0_RP )
          URBAN_UC(i,j) = max( real(URBAN_UC(i,j) + URBAN_UC_t(i,j) * dt, RP), 0.0_RP )
          URBAN_RAINR(i,j) = max( real(URBAN_RAINR(i,j) + URBAN_RAINR_t(i,j) * dt, RP), 0.0_RP )
          URBAN_RAINB(i,j) = max( real(URBAN_RAINB(i,j) + URBAN_RAINB_t(i,j) * dt, RP), 0.0_RP )
          URBAN_RAING(i,j) = max( real(URBAN_RAING(i,j) + URBAN_RAING_t(i,j) * dt, RP), 0.0_RP )
       end if
       end do
       end do
       !$acc end kernels

    end select

    call URBAN_vars_check

    call PROF_rapend  ('URB_Update', 1)

    return
  end subroutine URBAN_driver_update

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
       ATMOS_SFLX_water, &
       ATMOS_SFLX_ENGI
    use mod_cpl_vars, only: &
       CPL_getATM_URB
    implicit none

    real(RP) :: ATMOS_SFLX_rad_dn(UIA,UJA,N_RAD_DIR,N_RAD_RGN)

    integer  :: i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('URB_SfcExch', 3)

    !$acc data create(ATMOS_SFLX_rad_dn)

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
                            ATMOS_SFLX_water (:,:),     & ! [OUT]
                            ATMOS_SFLX_ENGI  (:,:)      ) ! [OUT]
    endif

!OCL XFILL
    !$omp parallel do
    !$acc kernels
    do j = UJS, UJE
    do i = UIS, UIE
       ATMOS_SFLX_LW(i,j,I_R_direct ) = ATMOS_SFLX_rad_dn(i,j,I_R_direct ,I_R_IR)    ! IR, direct
       ATMOS_SFLX_LW(i,j,I_R_diffuse) = ATMOS_SFLX_rad_dn(i,j,I_R_diffuse,I_R_IR)    ! IR, diffuse

       ATMOS_SFLX_SW(i,j,I_R_direct ) = ATMOS_SFLX_rad_dn(i,j,I_R_direct ,I_R_NIR) & ! NIR, direct
                                      + ATMOS_SFLX_rad_dn(i,j,I_R_direct ,I_R_VIS)   ! VIS, direct
       ATMOS_SFLX_SW(i,j,I_R_diffuse) = ATMOS_SFLX_rad_dn(i,j,I_R_diffuse,I_R_NIR) & ! NIR, diffuse
                                      + ATMOS_SFLX_rad_dn(i,j,I_R_diffuse,I_R_VIS)   ! VIS, diffuse
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    call PROF_rapend  ('URB_SfcExch', 3)

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
       URBAN_SFLX_SHEX,  &
       URBAN_SFLX_LHEX,  &
       URBAN_SFLX_QVEX,  &
       URBAN_SFLX_GH,    &
       URBAN_SFLX_QTRC,  &
       URBAN_Z0M,        &
       URBAN_Z0H,        &
       URBAN_Z0E,        &
       URBAN_U10,        &
       URBAN_V10,        &
       URBAN_T2,         &
       URBAN_Q2
    use mod_cpl_vars, only: &
       CPL_putURB
    use scale_landuse, only: &
       exists_urban => LANDUSE_exists_urban
    implicit none

    ! arguments
    logical, intent(in) :: countup
    !---------------------------------------------------------------------------

    call PROF_rapstart('URB_SfcExch', 3)

    if ( URBAN_do ) then
       call CPL_putURB( URBAN_SFC_TEMP  (:,:),     & ! [IN]
                        URBAN_SFC_albedo(:,:,:,:), & ! [IN]
                        URBAN_Z0M       (:,:),     & ! [IN]
                        URBAN_Z0H       (:,:),     & ! [IN]
                        URBAN_Z0E       (:,:),     & ! [IN]
                        URBAN_SFLX_MW   (:,:),     & ! [IN]
                        URBAN_SFLX_MU   (:,:),     & ! [IN]
                        URBAN_SFLX_MV   (:,:),     & ! [IN]
                        URBAN_SFLX_SH   (:,:),     & ! [IN]
                        URBAN_SFLX_LH   (:,:),     & ! [IN]
                        URBAN_SFLX_SHEX (:,:),     & ! [IN]
                        URBAN_SFLX_LHEX (:,:),     & ! [IN]
                        URBAN_SFLX_QVEX (:,:),     & ! [IN]
                        URBAN_SFLX_GH   (:,:),     & ! [IN]
                        URBAN_SFLX_QTRC (:,:,:),   & ! [IN]
                        URBAN_U10       (:,:),     & ! [IN]
                        URBAN_V10       (:,:),     & ! [IN]
                        URBAN_T2        (:,:),     & ! [IN]
                        URBAN_Q2        (:,:),     & ! [IN]
                        exists_urban    (:,:),     & ! [IN]
                        countup                    ) ! [IN]
    endif

    call PROF_rapend  ('URB_SfcExch', 3)

    return
  end subroutine URBAN_SURFACE_SET

end module mod_urban_driver
