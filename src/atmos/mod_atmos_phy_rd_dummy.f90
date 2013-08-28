!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          dummy code
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-03-21 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd
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
  public :: ATMOS_PHY_RD_setup
  public :: ATMOS_PHY_RD

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
  subroutine ATMOS_PHY_RD_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_RD
    use mod_atmos_vars_sf, only: &
       SWD, &
       LWD
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ dummy radiation process'

    if ( ATMOS_TYPE_PHY_RD /= 'NONE' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_RD is not NONE. Check!'
       call PRC_MPIstop
    endif

    SWD(:,:) = 0.0_RP
    LWD(:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_RD_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD
    use mod_const, only: &
       PI => CONST_PI
    use mod_time, only: &
       dt   => TIME_DTSEC_ATMOS_PHY_RD, &
       step => TIME_NOWSTEP
    use mod_atmos_vars_sf, only: &
       SWD, &
       LWD
    implicit none

    real(RP) :: factor

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Radiation(dummy)'

!    factor = ( 1.0_RP - cos ( dble(step)*dt/86400.0_RP * 2.0_RP*PI ) ) * 0.5_RP
!    SWD(:,:) = - 900.0_RP * factor
!    LWD(:,:) = - 400 - 10.0_RP * factor
    SWD(:,:) = 400.0_RP
    LWD(:,:) = 400.0_RP

    return
  end subroutine ATMOS_PHY_RD

end module mod_atmos_phy_rd
