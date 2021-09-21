!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom boundary of atmosphere (surface)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_sf_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_driver_setup
  public :: ATMOS_PHY_SF_driver_step

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
  subroutine ATMOS_PHY_SF_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_sf_bulk, only: &
       ATMOS_PHY_SF_bulk_setup
    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_sw_phy_sf
    use mod_atmos_phy_sf_vars, only: &
       SFC_Z0M   => ATMOS_PHY_SF_SFC_Z0M,   &
       SFC_Z0H   => ATMOS_PHY_SF_SFC_Z0H,   &
       SFC_Z0E   => ATMOS_PHY_SF_SFC_Z0E,   &
       SFLX_MW   => ATMOS_PHY_SF_SFLX_MW,   &
       SFLX_MU   => ATMOS_PHY_SF_SFLX_MU,   &
       SFLX_MV   => ATMOS_PHY_SF_SFLX_MV,   &
       SFLX_SH   => ATMOS_PHY_SF_SFLX_SH,   &
       SFLX_LH   => ATMOS_PHY_SF_SFLX_LH,   &
       SFLX_QTRC => ATMOS_PHY_SF_SFLX_QTRC
    use mod_cpl_admin, only: &
       CPL_sw
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_sf ) then

       if ( CPL_sw ) then
          LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'Coupler is enabled.'
       else
          ! setup library component
          select case( ATMOS_PHY_SF_TYPE )
          case ( 'BULK' )
             call ATMOS_PHY_SF_bulk_setup
          case ( 'CONST' )
             call ATMOS_PHY_SF_const_setup
          case default
             LOG_ERROR("ATMOS_PHY_SF_driver_setup",*) 'invalid Surface flux type(', trim(ATMOS_PHY_SF_TYPE), '). CHECK!'
             call PRC_abort
          end select
       endif

    else

       LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'this component is never called.'
       LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'surface fluxes are set to zero.'
       SFLX_MW(:,:,:) = 0.0_RP
       SFLX_MU(:,:,:) = 0.0_RP
       SFLX_MV(:,:,:) = 0.0_RP
       SFLX_SH(:,:,:) = 0.0_RP
       SFLX_LH(:,:,:) = 0.0_RP
       LOG_INFO("ATMOS_PHY_SF_driver_setup",*) 'SFC_TEMP, SFC_albedo is set in ATMOS_PHY_SF_vars.'

    endif

    SFLX_QTRC(:,:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_SF_driver_setup

  !-----------------------------------------------------------------------------
  !> time step
  subroutine ATMOS_PHY_SF_driver_step
    use scale_const, only: &
       UNDEF  => CONST_UNDEF,  &
       EPS    => CONST_EPS,    &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN, &
       CPdry  => CONST_CPdry,  &
       CVvap  => CONST_CVvap
    use scale_time, only: &
       dt_SF => TIME_DTSEC_ATMOS_PHY_SF
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       I_QV
    use scale_atmos_phy_sf_bulk, only: &
       ATMOS_PHY_SF_bulk_flux
    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_flux
    use mod_atmos_vars, only: &
       ATMOS_vars_calc_diagnostics, &
       DENS, &
       RHOU, &
       RHOV, &
       MOMZ, &
       RHOE, &
       RHOQ, &
       POTT, &
       TEMP, &
       PRES, &
       W,    &
       U,    &
       V,    &
       QV,   &
       CZ,   &
       FZ,   &
       Z1
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn, &
       SFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn
    use mod_atmos_phy_bl_vars, only: &
       PBL_Zi => ATMOS_PHY_BL_Zi
    use mod_atmos_phy_sf_vars, only: &
       SFC_DENS   => ATMOS_PHY_SF_SFC_DENS,   &
       SFC_PRES   => ATMOS_PHY_SF_SFC_PRES,   &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_Z0M    => ATMOS_PHY_SF_SFC_Z0M,    &
       SFC_Z0H    => ATMOS_PHY_SF_SFC_Z0H,    &
       SFC_Z0E    => ATMOS_PHY_SF_SFC_Z0E,    &
       SFLX_MW    => ATMOS_PHY_SF_SFLX_MW,    &
       SFLX_MU    => ATMOS_PHY_SF_SFLX_MU,    &
       SFLX_MV    => ATMOS_PHY_SF_SFLX_MV,    &
       SFLX_SH    => ATMOS_PHY_SF_SFLX_SH,    &
       SFLX_LH    => ATMOS_PHY_SF_SFLX_LH,    &
       SFLX_QTRC  => ATMOS_PHY_SF_SFLX_QTRC,  &
       U10        => ATMOS_PHY_SF_U10,        &
       V10        => ATMOS_PHY_SF_V10,        &
       T2         => ATMOS_PHY_SF_T2,         &
       Q2         => ATMOS_PHY_SF_Q2,         &
       Ustar      => ATMOS_PHY_SF_Ustar,      &
       Tstar      => ATMOS_PHY_SF_Tstar,      &
       Qstar      => ATMOS_PHY_SF_Qstar,      &
       Wstar      => ATMOS_PHY_SF_Wstar,      &
       RLmo       => ATMOS_PHY_SF_RLmo
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_atmos_admin, only: &
       ATMOS_PHY_SF_TYPE
    implicit none

    real(RP) :: ATM_W   (IA,JA,ADM_lall)
    real(RP) :: ATM_U   (IA,JA,ADM_lall)
    real(RP) :: ATM_V   (IA,JA,ADM_lall)
    real(RP) :: ATM_DENS(IA,JA,ADM_lall)
    real(RP) :: ATM_TEMP(IA,JA,ADM_lall)
    real(RP) :: ATM_PRES(IA,JA,ADM_lall)
    real(RP) :: ATM_QV  (IA,JA,ADM_lall)
    real(RP) :: SFLX_QV (IA,JA,ADM_lall)

    real(RP) :: work, dz

    integer :: i, j, iq, l
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall

       ! update surface density, surface pressure
       call BOTTOM_estimate( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                             DENS(:,:,:,l), PRES(:,:,:,l), QV(:,:,:,l), & ! [IN]
                             SFC_TEMP(:,:,l),                           & ! [IN]
                             FZ(:,:,:,l),                               & ! [IN]
                             SFC_DENS(:,:,l), SFC_PRES(:,:,l)           ) ! [OUT]

    end do

    if ( .NOT. CPL_sw ) then

       !$omp parallel do
       do l = 1, ADM_lall
       do j = JS, JE
       do i = IS, IE
          ATM_W   (i,j,l) = W   (KS,i,j,l)
          ATM_U   (i,j,l) = U   (KS,i,j,l)
          ATM_V   (i,j,l) = V   (KS,i,j,l)
          ATM_DENS(i,j,l) = DENS(KS,i,j,l)
          ATM_TEMP(i,j,l) = TEMP(KS,i,j,l)
          ATM_PRES(i,j,l) = PRES(KS,i,j,l)
          ATM_QV  (i,j,l) = QV  (KS,i,j,l)
       end do
       end do
       end do

       select case ( ATMOS_PHY_SF_TYPE )
       case ( 'BULK' )
          do l = 1, ADM_lall
             call ATMOS_PHY_SF_bulk_flux( &
                  IA, IS, IE, JA, JS, JE, &
                  ATM_W(:,:,l), ATM_U(:,:,l), ATM_V(:,:,l),          & ! [IN]
                  ATM_TEMP(:,:,l), ATM_PRES(:,:,l), ATM_QV(:,:,l),   & ! [IN]
                  SFC_DENS(:,:,l), SFC_TEMP(:,:,l), SFC_PRES(:,:,l), & ! [IN]
                  SFC_Z0M(:,:,l), SFC_Z0H(:,:,l), SFC_Z0E(:,:,l),    & ! [IN]
                  PBL_Zi(:,:,l), Z1(:,:,l),                          & ! [IN]
                  SFLX_MW(:,:,l), SFLX_MU(:,:,l), SFLX_MV(:,:,l),    & ! [OUT]
                  SFLX_SH(:,:,l), SFLX_LH(:,:,l), SFLX_QV(:,:,l),    & ! [OUT]
                  Ustar(:,:,l), Tstar(:,:,l), Qstar(:,:,l),          & ! [OUT]
                  Wstar(:,:,l), RLmo(:,:,l),                         & ! [OUT]
                  U10(:,:,l), V10(:,:,l), T2(:,:,l), Q2(:,:,l)       ) ! [OUT]
          end do

       case ( 'CONST' )

          do l = 1, ADM_lall
             call ATMOS_PHY_SF_const_flux( &
                  IA, IS, IE, JA, JS, JE, &
                  ATM_W(:,:,l), ATM_U(:,:,l), ATM_V(:,:,l), ATM_TEMP(:,:,l), & ! [IN]
                  Z1(:,:,l), SFC_DENS(:,:,l),                                & ! [IN]
                  SFLX_MW(:,:,l), SFLX_MU(:,:,l), SFLX_MV(:,:,l),            & ! [OUT]
                  SFLX_SH(:,:,l), SFLX_LH(:,:,l), SFLX_QV(:,:,l),            & ! [OUT]
                  U10(:,:,l), V10(:,:,l)                                     ) ! [OUT]
             Ustar(:,:,l) = UNDEF
             Tstar(:,:,l) = UNDEF
             Qstar(:,:,l) = UNDEF
             Wstar(:,:,l) = UNDEF
             RLmo (:,:,l) = UNDEF
             T2(:,:,l) = ATM_TEMP(:,:,l)
             Q2(:,:,l) = ATM_QV(:,:,l)

          end do

       end select

       if ( .not. ATMOS_HYDROMETEOR_dry ) then
          SFLX_QTRC(:,:,I_QV,:) = SFLX_QV(:,:,:)
       end if

    endif

    do l = 1, ADM_lall

       !$omp parallel do private(dz)
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          dz = CZ(KS+1,i,j,l) - CZ(KS,i,j,l)
          MOMZ(KS,i,j,l) = MOMZ(KS,i,j,l) + SFLX_MW(i,j,l) / dz * dt_SF

          dz = FZ(KS,i,j,l) - FZ(KS-1,i,j,l)
          RHOU(KS,i,j,l) = RHOU(KS,i,j,l) + SFLX_MU(i,j,l) / dz * dt_SF
          RHOV(KS,i,j,l) = RHOV(KS,i,j,l) + SFLX_MV(i,j,l) / dz * dt_SF
          RHOE(KS,i,j,l) = RHOE(KS,i,j,l) + SFLX_SH(i,j,l) / dz * dt_SF
       enddo
       enddo

       if ( .not. ATMOS_HYDROMETEOR_dry ) then
          !$omp parallel do private(work)
          do j = JS, JE
          do i = IS, IE
             work = SFLX_QTRC(i,j,I_QV,l) / ( FZ(KS,i,j,l) - FZ(KS-1,i,j,l) ) * dt_SF
             DENS(KS,i,j,l) = DENS(KS,i,j,l) + work
             RHOE(KS,i,j,l) = RHOE(KS,i,j,l) + work * CVvap * SFC_TEMP(i,j,l)
             RHOQ(KS,i,j,I_QV,l) = RHOQ(KS,i,j,I_QV,l) + work
          enddo
          enddo
       end if

    end do

    call history_output


    call ATMOS_vars_calc_diagnostics

    return
  end subroutine ATMOS_PHY_SF_driver_step

  subroutine history_output
    use scale_const, only: &
       UNDEF => CONST_UNDEF
   !  use scale_file_history, only: &
   !     FILE_HISTORY_in
    use mod_history, only: &
       history_in
    use scale_atmos_hydrostatic, only: &
       barometric_law_mslp => ATMOS_HYDROSTATIC_barometric_law_mslp
    use mod_atmos_vars, only: &
       TEMP, &
       PRES, &
       QV,   &
       CZ
    use mod_atmos_phy_sf_vars, only: &
       SFC_DENS   => ATMOS_PHY_SF_SFC_DENS,   &
       SFC_PRES   => ATMOS_PHY_SF_SFC_PRES,   &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo => ATMOS_PHY_SF_SFC_albedo, &
       SFC_Z0M    => ATMOS_PHY_SF_SFC_Z0M,    &
       SFC_Z0H    => ATMOS_PHY_SF_SFC_Z0H,    &
       SFC_Z0E    => ATMOS_PHY_SF_SFC_Z0E,    &
       SFLX_MW    => ATMOS_PHY_SF_SFLX_MW,    &
       SFLX_MU    => ATMOS_PHY_SF_SFLX_MU,    &
       SFLX_MV    => ATMOS_PHY_SF_SFLX_MV,    &
       SFLX_SH    => ATMOS_PHY_SF_SFLX_SH,    &
       SFLX_LH    => ATMOS_PHY_SF_SFLX_LH,    &
       SFLX_GH    => ATMOS_PHY_SF_SFLX_GH,    &
       U10        => ATMOS_PHY_SF_U10,        &
       V10        => ATMOS_PHY_SF_V10,        &
       T2         => ATMOS_PHY_SF_T2,         &
       Q2         => ATMOS_PHY_SF_Q2

    real(RP) :: MSLP  (IA,JA,ADM_lall) ! mean sea-level pressure [Pa]

    real(RP) :: Uabs10(IA,JA,ADM_lall) ! 10m absolute wind [m/s]

    integer :: i, j, l

    do l = 1, ADM_lall

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
          Uabs10(i,j,l) = sqrt( U10(i,j,l)**2 + V10(i,j,l)**2 )
       end do
       end do

       call barometric_law_mslp( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 PRES(:,:,:,l), TEMP(:,:,:,l), QV(:,:,:,l), & ! [IN]
                                 CZ(:,:,:,l),                               & ! [IN]
                                 MSLP(:,:,l)                                ) ! [OUT]
    end do


   !  call FILE_HISTORY_in( SFC_DENS  (:,:,:),                     'SFC_DENS',        'surface atmospheric density',         'kg/m3'   )
   !  call FILE_HISTORY_in( SFC_PRES  (:,:,:),                     'SFC_PRES',        'surface atmospheric pressure',        'Pa'      )
   !  call FILE_HISTORY_in( SFC_TEMP  (:,:,:),                     'SFC_TEMP',        'surface skin temperature (merged)',   'K'       )
   !  call FILE_HISTORY_in( SFC_albedo(:,:,I_R_direct ,I_R_IR ,:), 'SFC_ALB_IR_dir' , 'surface albedo (IR, direct, merged)', '1'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFC_albedo(:,:,I_R_diffuse,I_R_IR ,:), 'SFC_ALB_IR_dif' , 'surface albedo (IR, diffuse,merged)', '1'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFC_albedo(:,:,I_R_direct ,I_R_NIR,:), 'SFC_ALB_NIR_dir', 'surface albedo (NIR,direct, merged)', '1'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFC_albedo(:,:,I_R_diffuse,I_R_NIR,:), 'SFC_ALB_NIR_dif', 'surface albedo (NIR,diffuse,merged)', '1'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFC_albedo(:,:,I_R_direct ,I_R_VIS,:), 'SFC_ALB_VIS_dir', 'surface albedo (VIS,direct, merged)', '1'       , fill_halo= .true. )
   !  call FILE_HISTORY_in( SFC_albedo(:,:,I_R_diffuse,I_R_VIS,:), 'SFC_ALB_VIS_dif', 'surface albedo (VIS,diffuse,merged)', '1'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFC_Z0M   (:,:,:),                     'SFC_Z0M',         'roughness length (momentum)',         'm'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFC_Z0H   (:,:,:),                     'SFC_Z0H',         'roughness length (heat)',             'm'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFC_Z0E   (:,:,:),                     'SFC_Z0E',         'roughness length (vapor)',            'm'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_MW   (:,:,:),                     'MWFLX',           'w-momentum flux (merged)',            'kg/m/s2' )
   !  call FILE_HISTORY_in( SFLX_MU   (:,:,:),                     'MUFLX',           'u-momentum flux (merged)',            'kg/m/s2' )
   !  call FILE_HISTORY_in( SFLX_MV   (:,:,:),                     'MVFLX',           'v-momentum flux (merged)',            'kg/m/s2' )
   !  call FILE_HISTORY_in( SFLX_SH   (:,:,:),                     'SHFLX',           'sensible heat flux (merged)',         'W/m2'    , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_LH   (:,:,:),                     'LHFLX',           'latent heat flux (merged)',           'W/m2'    , fill_halo=.true. )
   !  call FILE_HISTORY_in( SFLX_GH   (:,:,:),                     'GHFLX',           'ground heat flux (merged)',           'W/m2'    , fill_halo=.true. )
   !  call FILE_HISTORY_in( Uabs10    (:,:,:),                     'Uabs10',          '10m absolute wind',                   'm/s'     , fill_halo=.true. )
   !  call FILE_HISTORY_in( U10       (:,:,:),                     'U10',             '10m x-wind',                          'm/s'     , fill_halo=.true. )
   !  call FILE_HISTORY_in( V10       (:,:,:),                     'V10',             '10m y-wind',                          'm/s'     , fill_halo=.true. )
   !  call FILE_HISTORY_in( T2        (:,:,:),                     'T2 ',             '2m air temperature',                  'K'       , fill_halo=.true. )
   !  call FILE_HISTORY_in( Q2        (:,:,:),                     'Q2 ',             '2m specific humidity',                'kg/kg'   , fill_halo=.true. )
   !  call FILE_HISTORY_in( MSLP      (:,:,:),                     'MSLP',            'mean sea-level pressure',             'Pa'      )

    do l = 1, ADM_lall

       call history_in( 'SFC_DENS',        var_g(  SFC_DENS  (:,:,l) ) )                       !  'surface atmospheric density',         'kg/m3'
       call history_in( 'SFC_PRES',        var_g(  SFC_PRES  (:,:,l) ) )                       !  'surface atmospheric pressure',        'Pa'
       call history_in( 'SFC_TEMP',        var_g(  SFC_TEMP  (:,:,l) ) )                       !  'surface skin temperature (merged)',   'K'
       call history_in( 'SFC_ALB_IR_dir' , var_g(  SFC_albedo(:,:,I_R_direct ,I_R_IR ,l) ) )   !  'surface albedo (IR, direct, merged)', '1'
       call history_in( 'SFC_ALB_IR_dif' , var_g(  SFC_albedo(:,:,I_R_diffuse,I_R_IR ,l) ) )   !  'surface albedo (IR, diffuse,merged)', '1'
       call history_in( 'SFC_ALB_NIR_dir', var_g(  SFC_albedo(:,:,I_R_direct ,I_R_NIR,l) ) )   !  'surface albedo (NIR,direct, merged)', '1'
       call history_in( 'SFC_ALB_NIR_dif', var_g(  SFC_albedo(:,:,I_R_diffuse,I_R_NIR,l) ) )   !  'surface albedo (NIR,diffuse,merged)', '1'
       call history_in( 'SFC_ALB_VIS_dir', var_g(  SFC_albedo(:,:,I_R_direct ,I_R_VIS,l) ) )   !  'surface albedo (VIS,direct, merged)', '1'
       call history_in( 'SFC_ALB_VIS_dif', var_g(  SFC_albedo(:,:,I_R_diffuse,I_R_VIS,l) ) )   !  'surface albedo (VIS,diffuse,merged)', '1'
       call history_in( 'SFC_Z0M',         var_g(  SFC_Z0M   (:,:,l) ) )                       !  'roughness length (momentum)',         'm'
       call history_in( 'SFC_Z0H',         var_g(  SFC_Z0H   (:,:,l) ) )                       !  'roughness length (heat)',             'm'
       call history_in( 'SFC_Z0E',         var_g(  SFC_Z0E   (:,:,l) ) )                       !  'roughness length (vapor)',            'm'
       call history_in( 'MWFLX',           var_g(  SFLX_MW   (:,:,l) ) )                       !  'w-momentum flux (merged)',            'kg/m/s2'
       call history_in( 'MUFLX',           var_g(  SFLX_MU   (:,:,l) ) )                       !  'u-momentum flux (merged)',            'kg/m/s2'
       call history_in( 'MVFLX',           var_g(  SFLX_MV   (:,:,l) ) )                       !  'v-momentum flux (merged)',            'kg/m/s2'
       call history_in( 'SHFLX',           var_g(  SFLX_SH   (:,:,l) ) )                       !  'sensible heat flux (merged)',         'W/m2'
       call history_in( 'LHFLX',           var_g(  SFLX_LH   (:,:,l) ) )                       !  'latent heat flux (merged)',           'W/m2'
       call history_in( 'GHFLX',           var_g(  SFLX_GH   (:,:,l) ) )                       !  'ground heat flux (merged)',           'W/m2'
       call history_in( 'Uabs10',          var_g(  Uabs10    (:,:,l) ) )                       !  '10m absolute wind',                   'm/s'
       call history_in( 'U10',             var_g(  U10       (:,:,l) ) )                       !  '10m x-wind',                          'm/s'
       call history_in( 'V10',             var_g(  V10       (:,:,l) ) )                       !  '10m y-wind',                          'm/s'
       call history_in( 'T2',              var_g(  T2        (:,:,l) ) )                       !  '2m air temperature',                  'K'
       call history_in( 'Q2',              var_g(  Q2        (:,:,l) ) )                       !  '2m specific humidity',                'kg/kg'
       call history_in( 'MSLP',            var_g(  MSLP      (:,:,l) ) )                       !  'mean sea-level pressure',             'Pa'

    enddo


    return
  end subroutine history_output

  ! the following function is conversion to use history_in(), which will be replaced by FILE_HISTOTY_in in future.
  function var_g(var_ij)
     implicit none
     real(RP) :: var_ij (IA,JA)
     real(RP) :: var_g  (ADM_gall_in, ADM_KNONE)
     integer  :: i, j, g
     do j = 1, JA
     do i = 1, IA
        g = i + ( j - 1 ) * ADM_imax
        var_g(g,ADM_KNONE) = var_ij(i,j)
     enddo
     enddo
  end function var_g


end module mod_atmos_phy_sf_driver
