!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          Aerosol Microphysics process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_ae
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
  public :: ATMOS_PHY_AE_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  abstract interface
     subroutine ae( &
          DENS, &
          MOMZ, &
          MOMX, &
          MOMY, &
          RHOT, &
          QTRC  )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       real(RP), intent(inout) :: DENS(KA,IA,JA)
       real(RP), intent(inout) :: MOMZ(KA,IA,JA)
       real(RP), intent(inout) :: MOMX(KA,IA,JA)
       real(RP), intent(inout) :: MOMY(KA,IA,JA)
       real(RP), intent(inout) :: RHOT(KA,IA,JA)
       real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
     end subroutine ae

     subroutine er( RE, QTRC, RH )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       real(RP), intent(out) :: Re  (KA,IA,JA,AE_QA) ! effective radius
       real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
       real(RP), intent(in)  :: RH  (KA,IA,JA)       ! relative humidity         [0-1]
     end subroutine er
  end interface

  procedure(ae), pointer :: ATMOS_PHY_AE => NULL()
  procedure(er), pointer :: ATMOS_PHY_AE_EffectiveRadius => NULL()
  public :: ATMOS_PHY_AE
  public :: ATMOS_PHY_AE_EffectiveRadius

  real(RP), public, pointer :: AE_DENS(:) ! aerosol density [kg/m3]=[g/L]

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
  !> Setup
  subroutine ATMOS_PHY_AE_setup( AE_TYPE )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef AE
    use NAME(scale_atmos_phy_ae_, AE,), only: &
       NAME(ATMOS_PHY_AE_, AE, _setup), &
       NAME(ATMOS_PHY_AE_, AE,)
       NAME(ATMOS_PHY_AE_, AE, _EffectiveRadius)
#else
    use scale_atmos_phy_ae_dummy, only: &
       ATMOS_PHY_AE_dummy_setup, &
       ATMOS_PHY_AE_dummy, &
       ATMOS_PHY_AE_dummy_EffectiveRadius
#endif
    implicit none

    character(len=H_SHORT), intent(in) :: AE_TYPE
    !---------------------------------------------------------------------------

    select case( AE_TYPE )
    case ( 'DUMMY', 'NONE' )
       call ATMOS_PHY_AE_dummy_setup( AE_TYPE )
       ATMOS_PHY_AE => ATMOS_PHY_AE_dummy
       ATMOS_PHY_AE_EffectiveRadius => ATMOS_PHY_AE_dummy_EffectiveRadius
    case default
       write(*,*) 'xxx invalid aerosol type(', AE_TYPE, '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine ATMOS_PHY_AE_setup

end module scale_atmos_phy_ae
