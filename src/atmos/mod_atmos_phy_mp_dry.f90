!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Dummy Cloud Microphysics for dry atmosphere
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-10-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
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
  public :: ATMOS_PHY_MP_setup
  public :: ATMOS_PHY_MP

  public :: ATMOS_PHY_MP_CloudFraction
  public :: ATMOS_PHY_MP_EffectiveRadius

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: MP_DENS(MP_QA) = 0.0_RP ! hydrometeor density [kg/m3]=[g/L]

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
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_MP
    implicit none

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Dry Atmosphere'

    if ( ATMOS_TYPE_PHY_MP /= 'DRY' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_MP is not DRY. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(dummy)'

    return
  end subroutine ATMOS_PHY_MP

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_CloudFraction( &
       cldfrac, &
       QTRC     )
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA)
    !---------------------------------------------------------------------------

    cldfrac(:,:,:) = 0.D0 ! dummy

    return
  end subroutine ATMOS_PHY_MP_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0  )
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QA) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
    !---------------------------------------------------------------------------

    Re(:,:,:,:) = 8.E-6_RP ! dummy

    return
  end subroutine ATMOS_PHY_MP_EffectiveRadius

end module mod_atmos_phy_mp
