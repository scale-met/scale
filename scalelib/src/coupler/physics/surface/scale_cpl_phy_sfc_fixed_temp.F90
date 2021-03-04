!-------------------------------------------------------------------------------
!> module coupler / surface fixed temp model
!!
!! @par Description
!!          Surface fixed temperature model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_cpl_phy_sfc_fixed_temp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_PHY_SFC_FIXED_TEMP_setup
  public :: CPL_PHY_SFC_FIXED_TEMP_finalize
  public :: CPL_PHY_SFC_FIXED_TEMP

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
  logical :: initialized = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_PHY_SFC_FIXED_TEMP_setup
    implicit none
    !---------------------------------------------------------------------------

    if ( initialized ) return

    LOG_NEWLINE
    LOG_INFO("CPL_PHY_SFC_FIXED_TEMP_setup",*) 'Setup'

    initialized = .true.

    return
  end subroutine CPL_PHY_SFC_FIXED_TEMP_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine CPL_PHY_SFC_FIXED_TEMP_finalize

    initialized = .false.

    return
  end subroutine CPL_PHY_SFC_FIXED_TEMP_finalize

  !-----------------------------------------------------------------------------
  subroutine CPL_PHY_SFC_fixed_temp( &
       IA, IS, IE,          &
       JA, JS, JE,          &
       TMPA, PRSA,          &
       WA, UA, VA,          &
       RHOA, QVA, LH,       &
       Z1, PBL,             &
       RHOS, PRSS,          &
       RFLXD,               &
       TMPS, WSTR, QVEF,    &
       ALBEDO,              &
       Rb, Z0M, Z0H, Z0E,   &
       calc_flag, dt,       &
       ZMFLX, XMFLX, YMFLX, &
       SHFLX, LHFLX, QVFLX, &
       GFLX,  &
       Ustar, Tstar, Qstar, &
       Wstar,               &
       RLmo,                &
       U10, V10, T2, Q2     )
    use scale_const, only: &
       EPS   => CONST_EPS, &
       UNDEF => CONST_UNDEF, &
       EPSvap => CONST_EPSvap, &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       Rvap  => CONST_Rvap,  &
       STB   => CONST_STB
    use scale_atmos_saturation, only: &
!       qsat => ATMOS_SATURATION_pres2qsat_all
       qsat => ATMOS_SATURATION_dens2qsat_all
    use scale_bulkflux, only: &
       BULKFLUX, &
       BULKFLUX_diagnose_surface
    implicit none

    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: TMPA     (IA,JA)                     ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in)  :: PRSA     (IA,JA)                     ! pressure    at the lowest atmospheric layer [Pa]
    real(RP), intent(in)  :: WA       (IA,JA)                     ! velocity w  at the lowest atmospheric layer [m/s]
    real(RP), intent(in)  :: UA       (IA,JA)                     ! velocity u  at the lowest atmospheric layer [m/s]
    real(RP), intent(in)  :: VA       (IA,JA)                     ! velocity v  at the lowest atmospheric layer [m/s]
    real(RP), intent(in)  :: RHOA     (IA,JA)                     ! density     at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in)  :: QVA      (IA,JA)                     ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in)  :: LH       (IA,JA)                     ! latent heat at the lowest atmospheric layer [J/kg]
    real(RP), intent(in)  :: Z1       (IA,JA)                     ! cell center height at the lowest atmospheric layer [m]
    real(RP), intent(in)  :: PBL      (IA,JA)                     ! the top of atmospheric mixing layer [m]
    real(RP), intent(in)  :: RHOS     (IA,JA)                     ! density  at the surface [kg/m3]
    real(RP), intent(in)  :: PRSS     (IA,JA)                     ! pressure at the surface [Pa]
    real(RP), intent(in)  :: RFLXD    (IA,JA,N_RAD_DIR,N_RAD_RGN) ! downward radiation flux at the surface (direct/diffuse,IR/near-IR/VIS) [J/m2/s]
    real(RP), intent(in)  :: TMPS     (IA,JA)                     ! surface temperature [K]
    real(RP), intent(in)  :: WSTR     (IA,JA)                     ! amount of the water storage [kg/m2]
    real(RP), intent(in)  :: QVEF     (IA,JA)                     ! efficiency of evaporation (0-1)
    real(RP), intent(in)  :: ALBEDO   (IA,JA,N_RAD_DIR,N_RAD_RGN) ! surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
    real(RP), intent(in)  :: Rb       (IA,JA)                     ! stomata resistance [1/s]
    real(RP), intent(in)  :: Z0M      (IA,JA)                     ! roughness length for momemtum [m]
    real(RP), intent(in)  :: Z0H      (IA,JA)                     ! roughness length for heat     [m]
    real(RP), intent(in)  :: Z0E      (IA,JA)                     ! roughness length for vapor    [m]
    logical,  intent(in)  :: calc_flag(IA,JA)                     ! to decide calculate or not
    real(DP), intent(in)  :: dt                                   ! delta time

    real(RP), intent(out) :: ZMFLX    (IA,JA)                     ! z-momentum      flux at the surface [kg/m/s2]
    real(RP), intent(out) :: XMFLX    (IA,JA)                     ! x-momentum      flux at the surface [kg/m/s2]
    real(RP), intent(out) :: YMFLX    (IA,JA)                     ! y-momentum      flux at the surface [kg/m/s2]
    real(RP), intent(out) :: SHFLX    (IA,JA)                     ! sensible heat   flux at the surface [J/m2/s]
    real(RP), intent(out) :: LHFLX    (IA,JA)                     ! latent heat     flux at the surface [J/m2/s]
    real(RP), intent(out) :: QVFLX    (IA,JA)                     ! water vapor     flux at the surface [kg/m2/s]
    real(RP), intent(out) :: GFLX     (IA,JA)                     ! subsurface heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: Ustar    (IA,JA)                     ! friction velocity         [m/s]
    real(RP), intent(out) :: Tstar    (IA,JA)                     ! temperature scale         [K]
    real(RP), intent(out) :: Qstar    (IA,JA)                     ! moisture scale            [kg/kg]
    real(RP), intent(out) :: Wstar    (IA,JA)                     ! convective velocity scale [m/s]
    real(RP), intent(out) :: RLmo     (IA,JA)                     ! inversed Obukhov length   [1/m]
    real(RP), intent(out) :: U10      (IA,JA)                     ! velocity u  at 10m [m/s]
    real(RP), intent(out) :: V10      (IA,JA)                     ! velocity v  at 10m [m/s]
    real(RP), intent(out) :: T2       (IA,JA)                     ! temperature at 2m  [K]
    real(RP), intent(out) :: Q2       (IA,JA)                     ! water vapor at 2m  [kg/kg]

    real(RP) :: emis ! surface longwave emission                 [J/m2/s]
    real(RP) :: LWD  ! surface downward longwave  radiation flux [J/m2/s]
    real(RP) :: LWU  ! surface upward   longwave  radiation flux [J/m2/s]
    real(RP) :: SWD  ! surface downward shortwave radiation flux [J/m2/s]
    real(RP) :: SWU  ! surface upward   shortwave radiation flux [J/m2/s]
    real(RP) :: res  ! residual

    real(RP) :: Uabs ! absolute velocity               [m/s]
    real(RP) :: Ra   ! Aerodynamic resistance (=1/Ce)  [1/s]

    real(RP) :: QVsat      ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS(IA,JA) ! water vapor mixing ratio at surface            [kg/kg]
    real(RP) :: Rtot       ! total gas constant
    real(RP) :: qdry       ! dry air mass ratio [kg/kg]

    real(RP) :: FracU10(IA,JA) ! calculation parameter for U10 [1]
    real(RP) :: FracT2 (IA,JA) ! calculation parameter for T2  [1]
    real(RP) :: FracQ2 (IA,JA) ! calculation parameter for Q2  [1]

    real(RP) :: MFLUX

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'coupler / physics / surface / FIXED-TEMP'

    ! calculate surface flux
    !$omp parallel do schedule(dynamic) collapse(2) &
#ifndef __GFORTRAN__
    !$omp default(none) &
    !$omp shared(IS,IE,JS,JE,EPS,UNDEF,Rdry,CPdry,bulkflux,dt, &
    !$omp        calc_flag,TMPA,QVA,QVS,LH,WA,UA,VA,Z1,PBL,PRSA,TMPS,WSTR,PRSS,RHOS,QVEF,Z0M,Z0H,Z0E,ALBEDO,RFLXD,Rb, &
    !$omp        FracU10,FracT2,FracQ2, &
    !$omp        SHFLX,LHFLX,QVFLX,GFLX,ZMFLX,XMFLX,YMFLX,Ustar,Tstar,Qstar,Wstar,RLmo,U10,V10,T2,Q2) &
#else
    !$omp default(shared) &
#endif
    !$omp private(qdry,Rtot,QVsat,Uabs,Ra,res,emis,LWD,LWU,SWD,SWU,MFLUX)
    do j = JS, JE
    do i = IS, IE
       if ( calc_flag(i,j) ) then

!          qdry = 1.0_RP - QVA(i,j)
!          Rtot = qdry * Rdry + QVA(i,j) * Rvap
!          call qsat( TMPS(i,j), PRSS(i,j), qdry, QVsat )
          call qsat( TMPS(i,j), RHOS(i,j), QVsat )

          QVS(i,j) = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
                   + (        QVEF(i,j) ) * QVsat

          Uabs = sqrt( WA(i,j)**2 + UA(i,j)**2 + VA(i,j)**2 )

          call BULKFLUX( TMPA(i,j), TMPS(i,j),                  & ! [IN]
                         PRSA(i,j), PRSS(i,j),                  & ! [IN]
                         QVA (i,j), QVS (i,j),                  & ! [IN]
                         Uabs, Z1(i,j), PBL(i,j),               & ! [IN]
                         Z0M(i,j), Z0H(i,j), Z0E(i,j),          & ! [IN]
                         Ustar(i,j), Tstar(i,j), Qstar(i,j),    & ! [OUT]
                         Wstar(i,j), RLmo(i,j), Ra,             & ! [OUT]
                         FracU10(i,j), FracT2(i,j), FracQ2(i,j) ) ! [OUT]

          if ( Uabs < EPS ) then
             ZMFLX(i,j) = 0.0_RP
             XMFLX(i,j) = 0.0_RP
             YMFLX(i,j) = 0.0_RP
          else
             MFLUX = - RHOS(i,j) * Ustar(i,j)**2
             ZMFLX(i,j) = MFLUX * WA(i,j) / Uabs
             XMFLX(i,j) = MFLUX * UA(i,j) / Uabs
             YMFLX(i,j) = MFLUX * VA(i,j) / Uabs
          end if
          SHFLX(i,j) = -RHOS(i,j) * Ustar(i,j) * Tstar(i,j) * CPdry
          QVFLX(i,j) = -RHOS(i,j) * Ustar(i,j) * Qstar(i,j) * Ra / ( Ra+Rb(i,j) )
          QVFLX(i,j) = min( QVFLX(i,j), WSTR(i,j) / real(dt,RP) )
          LHFLX(i,j) = QVFLX(i,j) * LH(i,j)

          emis = ( 1.0_RP-ALBEDO(i,j,I_R_diffuse,I_R_IR) ) * STB * TMPS(i,j)**4

          LWD  = RFLXD(i,j,I_R_diffuse,I_R_IR)
          LWU  = RFLXD(i,j,I_R_diffuse,I_R_IR)  * ALBEDO(i,j,I_R_diffuse,I_R_IR) + emis
          SWD  = RFLXD(i,j,I_R_direct ,I_R_NIR) &
               + RFLXD(i,j,I_R_diffuse,I_R_NIR) &
               + RFLXD(i,j,I_R_direct ,I_R_VIS) &
               + RFLXD(i,j,I_R_diffuse,I_R_VIS)
          SWU  = RFLXD(i,j,I_R_direct ,I_R_NIR) * ALBEDO(i,j,I_R_direct ,I_R_NIR) &
               + RFLXD(i,j,I_R_diffuse,I_R_NIR) * ALBEDO(i,j,I_R_diffuse,I_R_NIR) &
               + RFLXD(i,j,I_R_direct ,I_R_VIS) * ALBEDO(i,j,I_R_direct ,I_R_VIS) &
               + RFLXD(i,j,I_R_diffuse,I_R_VIS) * ALBEDO(i,j,I_R_diffuse,I_R_VIS)

          ! calculation for residual
          res = SWD - SWU + LWD - LWU - SHFLX(i,j) - QVFLX(i,j) * LH(i,j)

          ! put residual in ground heat flux (positive for downward)
          GFLX(i,j) = res

       else ! not calculate surface flux
          ZMFLX(i,j) = UNDEF
          XMFLX(i,j) = UNDEF
          YMFLX(i,j) = UNDEF
          SHFLX(i,j) = UNDEF
          LHFLX(i,j) = UNDEF
          QVFLX(i,j) = UNDEF
          GFLX (i,j) = UNDEF
          Ustar(i,j) = UNDEF
          Tstar(i,j) = UNDEF
          Qstar(i,j) = UNDEF
          Wstar(i,j) = UNDEF
          RLmo (i,j) = UNDEF
          U10  (i,j) = UNDEF
          V10  (i,j) = UNDEF
          T2   (i,j) = UNDEF
          Q2   (i,j) = UNDEF
       endif
    enddo
    enddo

    call BULKFLUX_diagnose_surface( IA, IS, IE, JA, JS, JE, &
                                    UA(:,:), VA(:,:),                      & ! (in)
                                    TMPA(:,:), QVA(:,:),                   & ! (in)
                                    TMPS(:,:), QVS(:,:),                   & ! (in)
                                    Z1(:,:), Z0M(:,:), Z0H(:,:), Z0E(:,:), & ! (in)
                                    U10(:,:), V10(:,:), T2(:,:), Q2(:,:),  & ! (out)
                                    mask = calc_flag(:,:),                 & ! (in)
                                    FracU10 = FracU10(:,:),                & ! (in)
                                    FracT2 = FracT2(:,:),                  & ! (in)
                                    FracQ2 = FracQ2(:,:)                   ) ! (in)

    return
  end subroutine CPL_PHY_SFC_fixed_temp

end module scale_cpl_phy_sfc_fixed_temp
