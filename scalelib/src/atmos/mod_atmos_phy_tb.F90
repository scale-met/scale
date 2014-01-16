!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_tb
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
  public :: ATMOS_PHY_TB_setup

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
     subroutine tb( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_qtrc,                & ! (out)
       tke, nu_C, Ri, Pr,                           & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC           ) ! (in)
       use mod_precision
       use mod_grid_index
       use mod_tracer
       implicit none
       real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
       real(RP), intent(out) :: qflx_sgs_qtrc(KA,IA,JA,QA,3)

       real(RP), intent(out) :: tke (KA,IA,JA) ! TKE
       real(RP), intent(out) :: nu_C(KA,IA,JA) ! eddy viscosity (center)
       real(RP), intent(out) :: Pr  (KA,IA,JA) ! Prantle number
       real(RP), intent(out) :: Ri  (KA,IA,JA) ! Richardson number

       real(RP), intent(in)  :: MOMZ(KA,IA,JA)
       real(RP), intent(in)  :: MOMX(KA,IA,JA)
       real(RP), intent(in)  :: MOMY(KA,IA,JA)
       real(RP), intent(in)  :: RHOT(KA,IA,JA)
       real(RP), intent(in)  :: DENS(KA,IA,JA)
       real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
     end subroutine tb
  end interface
  procedure(tb), pointer, public :: ATMOS_PHY_TB => NULL()

  !-----------------------------------------------------------------------------
contains

  subroutine ATMOS_PHY_TB_setup( &
       TB_TYPE, &
       CDZ, CDX, CDY, &
       FDZ, FDX, FDY, &
       CZ, FZ         )
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    use mod_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef TB
    use NAME(mod_atmos_phy_tb_, TB,), only: &
       NAME(ATMOS_PHY_TB_, TB, _setup), &
       NAME(ATMOS_PHY_TB_, TB,)
#else
    use mod_atmos_phy_tb_smg, only: &
       ATMOS_PHY_TB_smg_setup, &
       ATMOS_PHY_TB_smg
    use mod_atmos_phy_tb_dummy, only: &
       ATMOS_PHY_TB_dummy_setup, &
       ATMOS_PHY_TB_dummy
#endif
    implicit none
    character(len=IO_SYSCHR), intent(in) :: TB_TYPE


    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: FDX(IA-1)
    real(RP), intent(in) :: FDY(JA-1)
    real(RP), intent(in) :: CZ(KA)
    real(RP), intent(in) :: FZ(0:KA)
    !---------------------------------------------------------------------------

    select case( TB_TYPE )
    case ( 'SMAGORINSKY' )
       call ATMOS_PHY_tb_smg_setup( &
            TB_TYPE, &
            CDZ, CDX, CDY, &
            FDZ, FDX, FDY, &
            CZ, FZ )
       ATMOS_PHY_TB => ATMOS_PHY_TB_smg
    end select
  end subroutine ATMOS_PHY_TB_setup

end module mod_atmos_phy_tb
