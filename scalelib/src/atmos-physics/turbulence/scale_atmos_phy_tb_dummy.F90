!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          dummy code
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-01-22 (S.Nishizawa)       [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_tb_dummy
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
  public :: ATMOS_PHY_TB_dummy_setup
  public :: ATMOS_PHY_TB_dummy

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
  subroutine ATMOS_PHY_TB_dummy_setup( &
       TYPE_TB, &
       CDZ, CDX, CDY,   &
       FDZ, FDX, FDY,   &
       CZ, FZ )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: TYPE_TB

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: FDX(IA-1)
    real(RP), intent(in) :: FDY(JA-1)
    real(RP), intent(in) :: CZ(KA)
    real(RP), intent(in) :: FZ(0:KA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Turbulence Dummy'

    if ( TYPE_TB /= 'DUMMY' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_TB_TYPE is not Dummy. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_TB_dummy_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_dummy( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_qtrc,                & ! (out)
       tke, nu_C, Ri, Pr,                           & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC           ) ! (in)
    implicit none

    ! SGS flux
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
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Turbulence(dummy)'

    ! do nothing

    return
  end subroutine ATMOS_PHY_TB_dummy

end module scale_atmos_phy_tb_dummy
