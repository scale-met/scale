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

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( initialized ) return

    LOG_NEWLINE
    LOG_INFO("CPL_PHY_SFC_FIXED_TEMP_setup",*) 'Setup'

    initialized = .true.

    return
  end subroutine CPL_PHY_SFC_FIXED_TEMP_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_PHY_SFC_fixed_temp( &
       IA, IS, IE,          &
       JA, JS, JE,          &
       TMPA, PRSA,          &
       WA, UA, VA,          &
       RHOA, QVA, LHV,      &
       Z1, PBL,             &
       RHOS, PRSS,          &
       RFLXD,               &
       TMPS, QVEF,          &
       ALBEDO,              &
       Rb, Z0M, Z0H, Z0E,   &
       fact_area, dt,       &
       ZMFLX, XMFLX, YMFLX, &
       SHFLX, LHFLX, GFLX,  &
       U10, V10, T2, Q2     )
    use scale_const, only: &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       Rvap  => CONST_Rvap,  &
       STB   => CONST_STB
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
       BULKFLUX
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
    real(RP), intent(in)  :: LHV      (IA,JA)                     ! latent heat of vaporization [J/kg]
    real(RP), intent(in)  :: Z1       (IA,JA)                     ! cell center height at the lowest atmospheric layer [m]
    real(RP), intent(in)  :: PBL      (IA,JA)                     ! the top of atmospheric mixing layer [m]
    real(RP), intent(in)  :: RHOS     (IA,JA)                     ! density  at the surface [kg/m3]
    real(RP), intent(in)  :: PRSS     (IA,JA)                     ! pressure at the surface [Pa]
    real(RP), intent(in)  :: RFLXD    (IA,JA,N_RAD_DIR,N_RAD_RGN) ! downward radiation flux at the surface (direct/diffuse,IR/near-IR/VIS) [J/m2/s]
    real(RP), intent(in)  :: TMPS     (IA,JA)                     ! surface temperature [K]
    real(RP), intent(in)  :: QVEF     (IA,JA)                     ! efficiency of evaporation (0-1)
    real(RP), intent(in)  :: ALBEDO   (IA,JA,N_RAD_DIR,N_RAD_RGN) ! surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
    real(RP), intent(in)  :: Rb       (IA,JA)                     ! stomata resistance [1/s]
    real(RP), intent(in)  :: Z0M      (IA,JA)                     ! roughness length for momemtum [m]
    real(RP), intent(in)  :: Z0H      (IA,JA)                     ! roughness length for heat     [m]
    real(RP), intent(in)  :: Z0E      (IA,JA)                     ! roughness length for vapor    [m]
    real(RP), intent(in)  :: fact_area(IA,JA)                     ! to decide calculate or not
    real(DP), intent(in)  :: dt                                   ! delta time
    real(RP), intent(out) :: ZMFLX    (IA,JA)                     ! z-momentum      flux at the surface [kg/m/s2]
    real(RP), intent(out) :: XMFLX    (IA,JA)                     ! x-momentum      flux at the surface [kg/m/s2]
    real(RP), intent(out) :: YMFLX    (IA,JA)                     ! y-momentum      flux at the surface [kg/m/s2]
    real(RP), intent(out) :: SHFLX    (IA,JA)                     ! sensible heat   flux at the surface [J/m2/s]
    real(RP), intent(out) :: LHFLX    (IA,JA)                     ! latent heat     flux at the surface [J/m2/s]
    real(RP), intent(out) :: GFLX     (IA,JA)                     ! subsurface heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: U10      (IA,JA)                     ! velocity u  at 10m [m/s]
    real(RP), intent(out) :: V10      (IA,JA)                     ! velocity v  at 10m [m/s]
    real(RP), intent(out) :: T2       (IA,JA)                     ! temperature at 2m  [K]
    real(RP), intent(out) :: Q2       (IA,JA)                     ! water vapor at 2m  [kg/kg]

    real(RP) :: emis    ! surface longwave emission                 [J/m2/s]
    real(RP) :: LWD     ! surface downward longwave  radiation flux [J/m2/s]
    real(RP) :: LWU     ! surface upward   longwave  radiation flux [J/m2/s]
    real(RP) :: SWD     ! surface downward shortwave radiation flux [J/m2/s]
    real(RP) :: SWU     ! surface upward   shortwave radiation flux [J/m2/s]
    real(RP) :: res     ! residual

    real(RP) :: Ustar   ! friction velocity               [m]
    real(RP) :: Tstar   ! friction potential temperature  [K]
    real(RP) :: Qstar   ! friction water vapor mass ratio [kg/kg]
    real(RP) :: Uabs    ! modified absolute velocity      [m/s]
    real(RP) :: Ra      ! Aerodynamic resistance (=1/Ce)  [1/s]

    real(RP) :: QVsat   ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS     ! water vapor mixing ratio at surface            [kg/kg]
    real(RP) :: Rtot    ! total gas constant
    real(RP) :: qdry    ! dry air mass ratio [kg/kg]

    real(RP) :: FracU10 ! calculation parameter for U10 [1]
    real(RP) :: FracT2  ! calculation parameter for T2  [1]
    real(RP) :: FracQ2  ! calculation parameter for Q2  [1]

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'coupler / physics / surface / FIXED-TEMP'

    ! calculate surface flux
    !$omp parallel do default(none) &
    !$omp private(qdry,Rtot,QVsat,QVS,Ustar,Tstar,Qstar,Uabs,Ra,FracU10,FracT2,FracQ2,res,emis,LWD,LWU,SWD,SWU) &
    !$omp shared(IS,IE,JS,JE,Rdry,CPdry,bulkflux, &
    !$omp        fact_area,TMPA,QVA,LHV,UA,VA,WA,Z1,PBL,PRSA,TMPS,PRSS,RHOS,QVEF,Z0M,Z0H,Z0E,ALBEDO,RFLXD,Rb, &
    !$omp        SHFLX,LHFLX,GFLX,ZMFLX,XMFLX,YMFLX,U10,V10,T2,Q2)
    do j = JS, JE
    do i = IS, IE
       if ( fact_area(i,j) > 0.0_RP ) then

          qdry = 1.0_RP - QVA(i,j)
          Rtot = qdry * Rdry + QVA(i,j) * Rvap

          call qsat( TMPS(i,j), PRSS(i,j), qdry, QVsat )

          QVS = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
              + (        QVEF(i,j) ) * QVsat

          call BULKFLUX( Ustar,     & ! [OUT]
                         Tstar,     & ! [OUT]
                         Qstar,     & ! [OUT]
                         Uabs,      & ! [OUT]
                         Ra,        & ! [OUT]
                         FracU10,   & ! [OUT]
                         FracT2,    & ! [OUT]
                         FracQ2,    & ! [OUT]
                         TMPA(i,j), & ! [IN]
                         TMPS(i,j), & ! [IN]
                         PRSA(i,j), & ! [IN]
                         PRSS(i,j), & ! [IN]
                         QVA (i,j), & ! [IN]
                         QVS,       & ! [IN]
                         UA  (i,j), & ! [IN]
                         VA  (i,j), & ! [IN]
                         Z1  (i,j), & ! [IN]
                         PBL (i,j), & ! [IN]
                         Z0M (i,j), & ! [IN]
                         Z0H (i,j), & ! [IN]
                         Z0E (i,j)  ) ! [IN]

          ZMFLX(i,j) = -RHOS(i,j) * Ustar * Ustar / Uabs * WA(i,j)
          XMFLX(i,j) = -RHOS(i,j) * Ustar * Ustar / Uabs * UA(i,j)
          YMFLX(i,j) = -RHOS(i,j) * Ustar * Ustar / Uabs * VA(i,j)
          SHFLX(i,j) = -RHOS(i,j) * Ustar * Tstar * CPdry
          LHFLX(i,j) = -RHOS(i,j) * Ustar * Qstar * LHV(i,j) * Ra / ( Ra+Rb(i,j) )

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
          res = SWD - SWU + LWD - LWU - SHFLX(i,j) - LHFLX(i,j)

          ! put residual in ground heat flux
          GFLX(i,j) = -res

          ! diagnostic variables considering unstable/stable state
          !U10(i,j) = FracU10 * UA(i,j)
          !V10(i,j) = FracU10 * VA(i,j)
          !T2 (i,j) = ( 1.0_RP - FracT2 ) * TMPS(i,j) + FracT2 * TMPA(i,j)
          !Q2 (i,j) = ( 1.0_RP - FracQ2 ) * QVS       + FracQ2 * QVA (i,j)

          ! diagnostic variables for neutral state
          U10(i,j) = UA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
          V10(i,j) = VA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
          T2 (i,j) = TMPS(i,j) + ( TMPA(i,j) - TMPS(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                           / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
          Q2 (i,j) = QVS       + (  QVA(i,j) - QVS       ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                           / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )

       else ! not calculate surface flux
          ZMFLX(i,j) = 0.0_RP
          XMFLX(i,j) = 0.0_RP
          YMFLX(i,j) = 0.0_RP
          SHFLX(i,j) = 0.0_RP
          LHFLX(i,j) = 0.0_RP
          GFLX (i,j) = 0.0_RP
          U10  (i,j) = 0.0_RP
          V10  (i,j) = 0.0_RP
          T2   (i,j) = 0.0_RP
          Q2   (i,j) = 0.0_RP
       endif
    enddo
    enddo

    return
  end subroutine CPL_PHY_SFC_fixed_temp

end module scale_cpl_phy_sfc_fixed_temp
