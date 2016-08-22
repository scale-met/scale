!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)  [new]
!! @li      2014-04-11 (A.Noda)       [mod] add the grayzone module
!! @li      2014-05-01 (Y.Sato)       [mod] move grayzone module to mod_user
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_sf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_setup

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
  abstract interface
     subroutine sf( &
          ATM_TEMP, ATM_PRES, ATM_W, ATM_U, ATM_V,     &
          ATM_DENS,                                    &
          ATM_QTRC,                                    &
          ATM_Z1, dt,                                  &
          SFC_DENS, SFC_PRES,                          &
          SFLX_LW_dn, SFLX_SW_dn,                      &
          SFC_TEMP, SFC_albedo,                        &
          SFC_Z0M, SFC_Z0H, SFC_Z0E,                   &
          SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_LH, &
          SFLX_QTRC,                                   &
          U10, V10, T2, Q2                             )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       implicit none

       real(RP), intent(in)    :: ATM_TEMP  (IA,JA)    ! temperature at the lowermost layer (cell center) [K]
       real(RP), intent(in)    :: ATM_PRES  (IA,JA)    ! pressure    at the lowermost layer (cell center) [Pa]
       real(RP), intent(in)    :: ATM_W     (IA,JA)    ! velocity w  at the lowermost layer (cell center) [m/s]
       real(RP), intent(in)    :: ATM_U     (IA,JA)    ! velocity u  at the lowermost layer (cell center) [m/s]
       real(RP), intent(in)    :: ATM_V     (IA,JA)    ! velocity v  at the lowermost layer (cell center) [m/s]
       real(RP), intent(in)    :: ATM_DENS  (IA,JA)    ! density     at the lowermost layer (cell center) [kg/m3]
       real(RP), intent(in)    :: ATM_QTRC  (IA,JA,QA) ! tracer      at the lowermost layer (cell center) [kg/kg]
       real(RP), intent(in)    :: ATM_Z1    (IA,JA)    ! height of the lowermost grid from surface (cell center) [m]
       real(DP), intent(in)    :: dt                   ! delta time
       real(RP), intent(in)    :: SFC_DENS  (IA,JA)    ! density     at the surface atmosphere [kg/m3]
       real(RP), intent(in)    :: SFC_PRES  (IA,JA)    ! pressure    at the surface atmosphere [Pa]
       real(RP), intent(in)    :: SFLX_LW_dn(IA,JA)    ! downward longwave  radiation flux at the surface [J/m2/s]
       real(RP), intent(in)    :: SFLX_SW_dn(IA,JA)    ! downward shortwave radiation flux at the surface [J/m2/s]
       real(RP), intent(in)    :: SFC_TEMP  (IA,JA)    ! temperature at the surface skin [K]
       real(RP), intent(in)    :: SFC_albedo(IA,JA,2)  ! surface albedo (LW/SW) [0-1]
       real(RP), intent(inout) :: SFC_Z0M   (IA,JA)    ! surface roughness length (momentum) [m]
       real(RP), intent(inout) :: SFC_Z0H   (IA,JA)    ! surface roughness length (heat) [m]
       real(RP), intent(inout) :: SFC_Z0E   (IA,JA)    ! surface roughness length (vapor) [m]
       real(RP), intent(out)   :: SFLX_MW   (IA,JA)    ! surface flux for z-momentum    (area center)   [m/s*kg/m2/s]
       real(RP), intent(out)   :: SFLX_MU   (IA,JA)    ! surface flux for x-momentum    (area center)   [m/s*kg/m2/s]
       real(RP), intent(out)   :: SFLX_MV   (IA,JA)    ! surface flux for y-momentum    (area center)   [m/s*kg/m2/s]
       real(RP), intent(out)   :: SFLX_SH   (IA,JA)    ! surface flux for sensible heat (area center)   [J/m2/s]
       real(RP), intent(out)   :: SFLX_LH   (IA,JA)    ! surface flux for latent   heat (area center)   [J/m2/s]
       real(RP), intent(out)   :: SFLX_QTRC (IA,JA,QA) ! surface flux for tracer mass   (area center)   [kg/m2/s]
       real(RP), intent(out)   :: U10       (IA,JA)    ! velocity u        at 10m height
       real(RP), intent(out)   :: V10       (IA,JA)    ! velocity v        at 10m height
       real(RP), intent(out)   :: T2        (IA,JA)    ! temperature t     at  2m height
       real(RP), intent(out)   :: Q2        (IA,JA)    ! water vapor q     at  2m height
     end subroutine sf
  end interface
  procedure(sf), pointer :: ATMOS_PHY_SF => NULL()
  public :: ATMOS_PHY_SF

contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_setup( SF_TYPE )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef SF
    use NAME(scale_atmos_phy_, SF,), only: &
       NAME(ATMOS_PHY_, SF, _setup), &
       NAME(ATMOS_PHY_, SF,)
#else
    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_setup, &
       ATMOS_PHY_SF_const
    use scale_atmos_phy_sf_bulk, only: &
       ATMOS_PHY_SF_bulk_setup, &
       ATMOS_PHY_SF_bulk
#endif
    implicit none

    character(len=*), intent(in) :: SF_TYPE
    !---------------------------------------------------------------------------

    select case( SF_TYPE )
    case('CONST')
       call ATMOS_PHY_SF_const_setup( SF_TYPE )
       ATMOS_PHY_SF => ATMOS_PHY_SF_const
    case('BULK')
       call ATMOS_PHY_SF_bulk_setup( SF_TYPE )
       ATMOS_PHY_SF => ATMOS_PHY_SF_bulk
    case('OFF', 'COUPLE')
       ! do nothing
    case default
       write(*,*) 'xxx ATMPS_PHY_SF_TYPE is invalid'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_PHY_SF_setup

end module scale_atmos_phy_sf
