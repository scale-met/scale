!-------------------------------------------------------------------------------
!> module land / physics / snow / diagnostics
!!
!! @par Description
!!       Momentum fluxes and Diagnostic variables
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_land_phy_snow_diagnos
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_SNOW_DIAGS

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
  subroutine LAND_PHY_SNOW_DIAGS( &
       LIA, LIS, LIE, LJA, LJS, LJE, &
       SNOW_frac,           &
       TMPA, PRSA,          &
       WA, UA, VA,          &
       RHOA, QVA,           &
       Z1, PBL,             &
       RHOS, PRSS, LST1,    &
       QVEF,                &
       Z0M, Z0H, Z0E,       &
       ZMFLX, XMFLX, YMFLX, &
       Ustar, Tstar, Qstar, &
       Wstar,               &
       RLmo,                &
       U10, V10,            &
       T2, Q2               )
    use scale_const, only: &
      EPS   => CONST_EPS,   &
      UNDEF => CONST_UNDEF, &
      PRE00 => CONST_PRE00, &
      Rdry  => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      Rvap  => CONST_Rvap,  &
      STB   => CONST_STB
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_dens2qsat_all
!      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
      BULKFLUX, &
      BULKFLUX_diagnose_surface
    implicit none
    integer, intent(in) :: LIA, LIS, LIE
    integer, intent(in) :: LJA, LJS, LJE

    real(RP), intent(in) :: SNOW_frac(LIA,LJA) ! snow fraction
    real(RP), intent(in) :: TMPA(LIA,LJA) ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: PRSA(LIA,LJA) ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: WA  (LIA,LJA) ! velocity w at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: UA  (LIA,LJA) ! velocity u at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: VA  (LIA,LJA) ! velocity v at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: RHOA(LIA,LJA) ! density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: QVA (LIA,LJA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Z1  (LIA,LJA) ! cell center height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL (LIA,LJA) ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: RHOS(LIA,LJA) ! density  at the surface [kg/m3]
    real(RP), intent(in) :: PRSS(LIA,LJA) ! pressure at the surface [Pa]

    real(RP), intent(in) :: LST1  (LIA,LJA) ! land surface temperature [K]
    real(RP), intent(in) :: QVEF  (LIA,LJA) ! efficiency of evaporation (0-1)

    real(RP), intent(in) :: Z0M   (LIA,LJA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H   (LIA,LJA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (LIA,LJA) ! roughness length for vapor [m]

    real(RP), intent(out) :: ZMFLX(LIA,LJA) ! z-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: XMFLX(LIA,LJA) ! x-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: YMFLX(LIA,LJA) ! y-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: Ustar(LIA,LJA) ! friction velocity         [m/s]
    real(RP), intent(out) :: Tstar(LIA,LJA) ! temperature scale         [K]
    real(RP), intent(out) :: Qstar(LIA,LJA) ! moisture scale            [kg/kg]
    real(RP), intent(out) :: Wstar(LIA,LJA) ! convective velocity scale [m/s]
    real(RP), intent(out) :: RLmo (LIA,LJA) ! inversed Obukhov length   [1/m]
    real(RP), intent(out) :: U10  (LIA,LJA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (LIA,LJA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (LIA,LJA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (LIA,LJA) ! water vapor at 2m [kg/kg]


    ! works
    real(RP) :: Uabs ! modified absolute velocity [m/s]
    real(RP) :: Ra   ! Aerodynamic resistance (=1/Ce) [1/s]

    real(RP) :: QVsat        ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS(LIA,LJA) ! water vapor mixing ratio at surface [kg/kg]
!    real(RP) :: Rtot         ! total gas constant
!    real(RP) :: qdry         ! dry air mass ratio [kg/kg]

    real(RP) :: FracU10(LIA,LJA) ! calculation parameter for U10 [-]
    real(RP) :: FracT2 (LIA,LJA) ! calculation parameter for T2 [-]
    real(RP) :: FracQ2 (LIA,LJA) ! calculation parameter for Q2 [-]

    real(RP) :: MFLUX

    logical :: calc_flag(LIA,LJA)

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_INFO("LAND_PHY_SNOW_DIAGS",*) 'Snow surface diagnostic'

    ! calculate surface flux
    !$omp parallel do &
    !$omp private(QVsat,Uabs,MFLUX)
    do j = LJS, LJE
    do i = LIS, LIE

      if( SNOW_frac(i,j) > 0.0_RP ) then

         calc_flag(i,j) = .true.

!        qdry = 1.0_RP - QVA(i,j)
!        Rtot = qdry * Rdry + QVA(i,j) * Rvap

!        call qsat( LST1(i,j),  PRSS(i,j), qdry, & ! [IN]
!                   QVsat                        ) ! [OUT]
        call qsat( LST1(i,j),  RHOS(i,j), & ! [IN]
                   QVsat                  ) ! [OUT]

        QVS(i,j)  = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * QVsat

        Uabs = sqrt( WA(i,j)**2 + UA(i,j)**2 + VA(i,j)**2 )

        call BULKFLUX( TMPA(i,j), LST1(i,j),                  & ! [IN]
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

      else

         ! not calculate surface flux

         calc_flag(i,j) = .false.

         ZMFLX(i,j) = UNDEF
         XMFLX(i,j) = UNDEF
         YMFLX(i,j) = UNDEF
         Ustar(i,j) = UNDEF
         Tstar(i,j) = UNDEF
         Qstar(i,j) = UNDEF
         Wstar(i,j) = UNDEF
         RLmo (i,j) = UNDEF
         U10  (i,j) = UNDEF
         V10  (i,j) = UNDEF
         T2   (i,j) = UNDEF
         Q2   (i,j) = UNDEF

      end if

    end do
    end do

    call BULKFLUX_diagnose_surface( LIA, LIS, LIE, LJA, LJS, LJE, &
                                    UA(:,:), VA(:,:),                      & ! (in)
                                    TMPA(:,:), QVA(:,:),                   & ! (in)
                                    LST1(:,:), QVS(:,:),                   & ! (in)
                                    Z1(:,:), Z0M(:,:), Z0H(:,:), Z0E(:,:), & ! (in)
                                    U10(:,:), V10(:,:), T2(:,:), Q2(:,:),  & ! (out)
                                    mask = calc_flag(:,:),                 & ! (in)
                                    FracU10 = FracU10(:,:),                & ! (in)
                                    FracT2 = FracT2(:,:),                  & ! (in)
                                    FracQ2 = FracQ2(:,:)                   ) ! (in)

    return
  end subroutine LAND_PHY_SNOW_DIAGS

end module scale_land_phy_snow_diagnos
