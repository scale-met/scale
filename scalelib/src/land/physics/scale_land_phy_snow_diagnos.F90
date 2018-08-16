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
       U10, V10,            &
       T2, Q2               )
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

    real(RP), intent(in) :: LST1   (LIA,LJA) ! land surface temperature [K]
    real(RP), intent(in) :: QVEF  (LIA,LJA) ! efficiency of evaporation (0-1)

    real(RP), intent(in) :: Z0M   (LIA,LJA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H   (LIA,LJA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (LIA,LJA) ! roughness length for vapor [m]

    real(RP), intent(out) :: ZMFLX(LIA,LJA) ! z-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: XMFLX(LIA,LJA) ! x-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: YMFLX(LIA,LJA) ! y-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: U10  (LIA,LJA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (LIA,LJA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (LIA,LJA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (LIA,LJA) ! water vapor at 2m [kg/kg]


    ! works
    real(RP) :: Ustar  ! friction velocity [m]
    real(RP) :: Tstar  ! friction potential temperature [K]
    real(RP) :: Qstar  ! friction water vapor mass ratio [kg/kg]
    real(RP) :: Uabs   ! modified absolute velocity [m/s]
    real(RP) :: Ra     ! Aerodynamic resistance (=1/Ce) [1/s]

    real(RP) :: QVsat  ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS    ! water vapor mixing ratio at surface [kg/kg]
    real(RP) :: Rtot   ! total gas constant
    real(RP) :: qdry   ! dry air mass ratio [kg/kg]

    real(RP) :: FracU10 ! calculation parameter for U10 [-]
    real(RP) :: FracT2  ! calculation parameter for T2 [-]
    real(RP) :: FracQ2  ! calculation parameter for Q2 [-]

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_INFO("LAND_PHY_SNOW_DIAGS",*) 'Snow surface diagnostic'

    ! calculate surface flux
    do j = LJS, LJE
    do i = LIS, LIE

      if( SNOW_frac(i,j) > 0.0_RP ) then

        qdry = 1.0_RP - QVA(i,j)
        Rtot = qdry * Rdry + QVA(i,j) * Rvap

        call qsat( LST1(i,j),  PRSS(i,j), qdry, & ! [IN]
                   QVsat                        ) ! [OUT]

        QVS  = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * QVsat

        call BULKFLUX( &
            Ustar,     & ! [OUT]
            Tstar,     & ! [OUT]
            Qstar,     & ! [OUT]
            Uabs,      & ! [OUT]
            Ra,        & ! [OUT]
            FracU10,   & ! [OUT]
            FracT2,    & ! [OUT]
            FracQ2,    & ! [OUT]
            TMPA(i,j), & ! [IN]
            LST1(i,j), & ! [IN]
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

        ! diagnostic variables for neutral state
        U10(i,j) = UA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        V10(i,j) = VA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        T2 (i,j) = LST1(i,j) + ( TMPA(i,j) - LST1(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                         / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
        Q2 (i,j) = QVS       + (  QVA(i,j) - QVS       ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                         / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )

      else

        ! not calculate surface flux
        ZMFLX(i,j) = 0.0_RP
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
        U10  (i,j) = 0.0_RP
        V10  (i,j) = 0.0_RP
        T2   (i,j) = 0.0_RP
        Q2   (i,j) = 0.0_RP

      end if

    end do
    end do

    return
  end subroutine LAND_PHY_SNOW_DIAGS

end module scale_land_phy_snow_diagnos
