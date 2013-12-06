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
module mod_atmos_phy_ae
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_tracer
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
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
     subroutine ae
     end subroutine ae
     subroutine er( RE, QTRC, RH )
       use mod_precision
       use mod_index
       use mod_tracer
       real(RP), intent(out) :: Re  (KA,IA,JA,AE_QA) ! effective radius
       real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
       real(RP), intent(in)  :: RH  (KA,IA,JA)       ! relative humidity         [0-1]
     end subroutine er
  end interface
  procedure(ae), public, pointer :: ATMOS_PHY_AE => NULL()
  procedure(er), public, pointer :: ATMOS_PHY_AE_EffectiveRadius => NULL()

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
    use mod_stdio, only: &
       IO_SYSCHR
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef AE
    use NAME(mod_atmos_phy_ae_, AE,), only: &
       NAME(ATMOS_PHY_AE_, AE, _setup), &
       NAME(ATMOS_PHY_AE_, AE,)
       NAME(ATMOS_PHY_AE_, AE, _EffectiveRadius)
#else
    use mod_atmos_phy_ae_dummy, only: &
       ATMOS_PHY_AE_dummy_setup, &
       ATMOS_PHY_AE_dummy, &
       ATMOS_PHY_AE_dummy_EffectiveRadius
#endif
    implicit none
    character(len=IO_SYSCHR), intent(in) :: AE_TYPE
    !---------------------------------------------------------------------------

    select case( AE_TYPE )
    case ( 'DUMMY' )
       call ATMOS_PHY_AE_dummy_setup( AE_TYPE )
       ATMOS_PHY_AE => ATMOS_PHY_AE_dummy
       ATMOS_PHY_AE_EffectiveRadius => ATMOS_PHY_AE_dummy_EffectiveRadius
    end select

    return
  end subroutine ATMOS_PHY_AE_setup

end module mod_atmos_phy_ae 
