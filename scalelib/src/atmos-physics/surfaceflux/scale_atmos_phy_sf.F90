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
!!          2014-04-11 (A.Noda)       [mod] add the grayzone module
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
          SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
          DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, SST,             & ! (in)
          CZ, ctime                                            ) ! (in)
       use scale_precision
       use scale_grid_index
       use scale_tracer
       implicit none

       real(RP), intent(out) :: SFLX_MOMZ(IA,JA)
       real(RP), intent(out) :: SFLX_MOMX(IA,JA)
       real(RP), intent(out) :: SFLX_MOMY(IA,JA)
       real(RP), intent(out) :: SFLX_POTT(IA,JA)
       real(RP), intent(out) :: SFLX_QV  (IA,JA)

       real(RP), intent(in)  :: DENS(KA,IA,JA)
       real(RP), intent(in)  :: MOMZ(KA,IA,JA)
       real(RP), intent(in)  :: MOMX(KA,IA,JA)
       real(RP), intent(in)  :: MOMY(KA,IA,JA)
       real(RP), intent(in)  :: RHOT(KA,IA,JA)
       real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
       real(RP), intent(in)  :: SST (1,IA,JA)

       real(RP), intent(in)  :: CZ(KA)
       real(DP), intent(in)  :: ctime
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
    use scale_atmos_phy_sf_louis, only: &
       ATMOS_PHY_SF_louis_setup, &
       ATMOS_PHY_SF_louis
#endif
    implicit none

    character(len=H_SHORT), intent(in) :: SF_TYPE
    !---------------------------------------------------------------------------

    select case( SF_TYPE )
    case ( 'CONST')
       call ATMOS_PHY_SF_const_setup( SF_TYPE )
       ATMOS_PHY_SF => ATMOS_PHY_SF_const
    case ( 'LOUIS')
       call ATMOS_PHY_SF_louis_setup( SF_TYPE )
       ATMOS_PHY_SF => ATMOS_PHY_SF_louis
    end select

    return
  end subroutine ATMOS_PHY_SF_setup

end module scale_atmos_phy_sf
