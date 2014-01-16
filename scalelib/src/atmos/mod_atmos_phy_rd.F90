!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process wrapper
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_setup

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
  abstract interface
     subroutine rd( &
          flux_rad, flux_top,  & ! [out]
          solins, cosSZA,      & ! [out]
          DENS, RHOT, QTRC,    & ! [in]
          temp_sfc, param_sfc, & ! [in]
          CZ, FZ, CDZ, RCDZ,   & ! [in]
          REAL_lon, REAL_lat,  & ! [in]
          TIME_NOWDATE         ) ! [in]
       use mod_precision
       use mod_grid_index
       use mod_tracer
       implicit none

       real(RP), intent(out) :: flux_rad(KA,IA,JA,2,2)
       real(RP), intent(out) :: flux_top(IA,JA,2)
       real(RP), intent(out) :: solins(IA,JA)
       real(RP), intent(out) :: cosSZA(IA,JA)
       real(RP), intent(in)  :: DENS(KA,IA,JA)
       real(RP), intent(in)  :: RHOT(KA,IA,JA)
       real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
       real(RP), intent(in)  :: temp_sfc(IA,JA)
       real(RP), intent(in)  :: param_sfc(5)
       real(RP), intent(in)  :: CZ(KA)
       real(RP), intent(in)  :: FZ(KA-1)
       real(RP), intent(in)  :: CDZ(KA)
       real(RP), intent(in)  :: RCDZ(KA)
       real(RP), intent(in)  :: REAL_lon(IA,JA)
       real(RP), intent(in)  :: REAL_lat(IA,JA)
       integer , intent(in)  :: TIME_NOWDATE(6)
     end subroutine rd
  end interface
  procedure(rd), pointer, public :: ATMOS_PHY_RD => NULL()

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_setup( RD_TYPE )
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    use mod_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef RD
    use NAME(mod_atmos_phy_rd_, RD,), only: &
       NAME(ATMOS_PHY_RD_, RD, _setup), &
       NAME(ATMOS_PHY_RD_, RD,)
#else
    use mod_atmos_phy_rd_mstrnX, only: &
       ATMOS_PHY_RD_mstrnX_setup, &
       ATMOS_PHY_RD_mstrnX
    use mod_atmos_phy_rd_dycoms2, only: &
       ATMOS_PHY_RD_dycoms2_setup, &
       ATMOS_PHY_RD_dycoms2
    use mod_atmos_phy_rd_dummy, only: &
       ATMOS_PHY_RD_dummy_setup, &
       ATMOS_PHY_RD_dummy
#endif
    implicit none
    character(len=IO_SYSCHR), intent(in) :: RD_TYPE

    !---------------------------------------------------------------------------

    select case ( RD_TYPE )
    case ( 'MSTRNX' )
       call ATMOS_PHY_RD_mstrnX_setup( RD_TYPE )
       ATMOS_PHY_RD => ATMOS_PHY_RD_mstrnX
    case ( 'DYCOMSII' )
       call ATMOS_PHY_RD_dycoms2_setup( RD_TYPE )
       ATMOS_PHY_RD => ATMOS_PHY_RD_dycoms2
    case default
       write(*,*) 'xxx invalid Radiation type(', trim(RD_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_PHY_RD_setup

end module mod_atmos_phy_rd
