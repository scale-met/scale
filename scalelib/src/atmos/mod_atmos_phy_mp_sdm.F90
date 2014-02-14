!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Dummy module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "macro_thermodyn.h"
module mod_atmos_phy_mp_sdm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  use mod_tracer_sdm

  use mod_process, only: &
     PRC_MPIstop
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  public :: ATMOS_PHY_MP_sdm_setup
  public :: ATMOS_PHY_MP_sdm
  public :: ATMOS_PHY_MP_sdm_CloudFraction
  public :: ATMOS_PHY_MP_sdm_EffectiveRadius
  public :: ATMOS_PHY_MP_sdm_Mixingratio
  public :: ATMOS_PHY_MP_sdm_restart_out
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]
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
  subroutine ATMOS_PHY_MP_sdm_setup( MP_TYPE, DENS, RHOT, QTRC )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: MP_TYPE

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'

    if ( IO_L ) write(IO_FID_LOG,*) 'xxx SDM module is not used by default. Please contact developers.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  subroutine ATMOS_PHY_MP_sdm( &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use mod_tracer, only: &
       QAD => QA
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)
    !---------------------------------------------------------------------------

    if ( IO_L ) write(IO_FID_LOG,*) 'xxx SDM module is not used by default. Please contact developers.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm

  !----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_sdm_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0  )
    use mod_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)        ! density                   [kg/m3]
    !---------------------------------------------------------------------------

    if ( IO_L ) write(IO_FID_LOG,*) 'xxx SDM module is not used by default. Please contact developers.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm_EffectiveRadius

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sdm_CloudFraction( &
       cldfrac, &
       QTRC     )
    use mod_tracer, only: &
       QAD => QA
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QAD)
    !---------------------------------------------------------------------------

    if ( IO_L ) write(IO_FID_LOG,*) 'xxx SDM module is not used by default. Please contact developers.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_Mixingratio( &
       Qe,    &
       QTRC0  )
    use mod_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,MP_QAD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    !---------------------------------------------------------------------------

    if ( IO_L ) write(IO_FID_LOG,*) 'xxx SDM module is not used by default. Please contact developers.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm_Mixingratio

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_in
    implicit none
    !---------------------------------------------------------------------------

    if ( IO_L ) write(IO_FID_LOG,*) 'xxx SDM module is not used by default. Please contact developers.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm_restart_in

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_out(otime)
    implicit none

    real(DP), intent(in) :: otime
    !---------------------------------------------------------------------------

    if ( IO_L ) write(IO_FID_LOG,*) 'xxx SDM module is not used by default. Please contact developers.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm_restart_out

end module mod_atmos_phy_mp_sdm
!-------------------------------------------------------------------------------
