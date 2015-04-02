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
       TYPE_TB,       &
       CDZ, CDX, CDY, &
       CZ             )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: TYPE_TB

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA,IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Turbulence Dummy'

    if ( TYPE_TB /= 'DUMMY' ) then
       write(*,*) 'xxx ATMOS_PHY_TB_TYPE is not Dummy. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_TB_dummy_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_dummy( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       tke,                                         & ! (inout) diagnostic variables
       nu_C, Ri, Pr,                                & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          & ! (in)
       sflx_mw, sflx_mu, sflx_mv, sflx_sh, sflx_qv, & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)
    implicit none

    ! SGS flux
    real(RP), intent(out)   :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_rhoq(KA,IA,JA,QA,3)

    real(RP), intent(inout) :: tke (KA,IA,JA) ! TKE

    real(RP), intent(out)   :: nu_C(KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out)   :: Pr  (KA,IA,JA) ! Prantle number
    real(RP), intent(out)   :: Ri  (KA,IA,JA) ! Richardson number

    real(RP), intent(in)    :: MOMZ(KA,IA,JA)
    real(RP), intent(in)    :: MOMX(KA,IA,JA)
    real(RP), intent(in)    :: MOMY(KA,IA,JA)
    real(RP), intent(in)    :: RHOT(KA,IA,JA)
    real(RP), intent(in)    :: DENS(KA,IA,JA)
    real(RP), intent(in)    :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)    :: sflx_mw(IA,JA)
    real(RP), intent(in)    :: sflx_mu(IA,JA)
    real(RP), intent(in)    :: sflx_mv(IA,JA)
    real(RP), intent(in)    :: sflx_sh(IA,JA)
    real(RP), intent(in)    :: sflx_qv(IA,JA)

    real(RP), intent(in)    :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)    :: J13G (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J33G              !< (3,3) element of Jacobian matrix
    real(RP), intent(in)    :: MAPF (IA,JA,2,4)  !< map factorp
    real(RP), intent(in)    :: dt
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Turbulence(dummy)'

    ! do nothing
    qflx_sgs_momz(:,:,:,:)   = 0.0_RP
    qflx_sgs_momx(:,:,:,:)   = 0.0_RP
    qflx_sgs_momy(:,:,:,:)   = 0.0_RP
    qflx_sgs_rhot(:,:,:,:)   = 0.0_RP
    qflx_sgs_rhoq(:,:,:,:,:) = 0.0_RP

    tke (:,:,:) = 0.0_RP

    nu_C(:,:,:) = 0.0_RP
    Pr  (:,:,:) = 0.0_RP
    Ri  (:,:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_TB_dummy

end module scale_atmos_phy_tb_dummy
