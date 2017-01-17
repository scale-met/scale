!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cumulus Parameterization
!!
!! @par Description
!!         Cumulus convection parameterization
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-06-27 (S.Matsugishi) [new]
!!
!<
module scale_atmos_phy_cp
  !------------------------------------------------------------------------------
  !
  !+++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  use scale_atmos_phy_mp
  private

  abstract interface
     subroutine cp( &
          DENS,           &
          MOMZ,           &
          MOMX,           &
          MOMY,           &
          RHOT,           &
          QTRC,           &
          DENS_t_CP,      &
          MOMZ_t_CP,      &
          MOMX_t_CP,      &
          MOMY_t_CP,      &
          RHOT_t_CP,      &
          RHOQ_t_CP,      &
          MFLX_cloudbase, &
          SFLX_convrain,  &
          cloudtop,       &
          cloudbase,      &
          cldfrac_dp,     &
          cldfrac_sh,     &
          kf_nca,         &
          kf_w0avg        )
       use scale_precision
       use scale_stdio
       use scale_prof
       use scale_grid_index
       use scale_tracer
       use scale_atmos_phy_mp
       real(RP), intent(in)    :: DENS(KA,IA,JA)
       real(RP), intent(in)    :: MOMX(KA,IA,JA)
       real(RP), intent(in)    :: MOMY(KA,IA,JA)
       real(RP), intent(in)    :: MOMZ(KA,IA,JA)
       real(RP), intent(in)    :: RHOT(KA,IA,JA)
       real(RP), intent(in)    :: QTRC(KA,IA,JA,QA)
       real(RP), intent(inout) :: DENS_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMZ_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMX_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: MOMY_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: RHOT_t_CP(KA,IA,JA)
       real(RP), intent(inout) :: RHOQ_t_CP(KA,IA,JA,QA_MP)
       real(RP), intent(inout) :: MFLX_cloudbase(IA,JA)
       real(RP), intent(inout) :: SFLX_convrain(IA,JA)
       real(RP), intent(inout) :: cloudtop(IA,JA)
       real(RP), intent(inout) :: cloudbase(IA,JA)
       real(RP), intent(inout) :: cldfrac_dp(KA,IA,JA)
       real(RP), intent(inout) :: cldfrac_sh(KA,IA,JA)
       real(RP), intent(inout) :: kf_nca(IA,JA)
       real(RP), intent(inout) :: kf_W0avg(KA,IA,JA)
     end subroutine cp
  end interface
  procedure(cp), pointer :: ATMOS_PHY_CP => NULL()
  public :: ATMOS_PHY_CP
  public :: ATMOS_PHY_CP_setup

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup Cumulus parameterization
  !------------------------------------------------------------------------------
  subroutine ATMOS_PHY_CP_setup( CP_TYPE )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef CP
    use NAME(scale_atmos_phy_mp_, CP,), only: &
       NAME(ATMOS_PHY_CP_, CP, _setup), &
       NAME(ATMOS_PHY_CP_, CP,), &
#else
    use scale_atmos_phy_cp_kf, only: &
         ATMOS_PHY_CP_kf_setup, &
         ATMOS_PHY_CP_kf
#endif
    implicit none

    character(len=*), intent(in) :: CP_TYPE
    !------------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** => ', trim(CP_TYPE), ' is selected.'

#ifdef CP
    call NAME(ATMOS_PHY_CP_, CP, _setup)( CP_TYPE )
    ATMOS_PHY_CP => NAME(ATMOS_PHY_CP_, CP,)
#else
    select case( CP_TYPE )
    case('OFF')
       ! do nothing
    case('KF')
       call ATMOS_PHY_CP_kf_setup( CP_TYPE )
       ATMOS_PHY_CP => ATMOS_PHY_CP_kf
    case default
       write(*,*) 'xxx invalid Cumulus parameterization type(', trim(CP_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select
#endif

    return
  end subroutine ATMOS_PHY_CP_setup

end module scale_atmos_phy_cp
