!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Turbulence
!!
!! @par Description
!!          dummy code
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2013-01-22 (S.Nishizawa)       [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_tb
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
  public :: ATMOS_PHY_TB_setup
  public :: ATMOS_PHY_TB

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
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
contains

  subroutine ATMOS_PHY_TB_setup()
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_TB
    implicit none

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Turbulence]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Dummy'

    if ( ATMOS_TYPE_PHY_TB .ne. 'DUMMY' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_TB is not Dummy. Check!'
       call PRC_MPIstop
    end if



    return
  end subroutine ATMOS_PHY_TB_setup

  subroutine ATMOS_PHY_TB
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Turbulence(dummy)'

    return
  end subroutine ATMOS_PHY_TB


end module mod_atmos_phy_tb
