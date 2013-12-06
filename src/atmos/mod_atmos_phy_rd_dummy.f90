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
       ATMOS_PHY_RD_TYPE
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ dummy radiation process'

    if ( ATMOS_PHY_RD_TYPE /= 'NONE' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_RD_TYPE is not NONE. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_RD_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD( update_flag, history_flag )
    implicit none
    logical, intent(in) :: update_flag
    logical, intent(in), optional :: history_flag
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Radiation(dummy)'

    return
  end subroutine ATMOS_PHY_RD

end module mod_atmos_phy_rd
