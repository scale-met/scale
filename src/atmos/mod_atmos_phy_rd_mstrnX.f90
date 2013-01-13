!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          mstrnX
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-26 (H.Yashiro) [new]
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
!    use mod_stdio, only: &
!       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_RD
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ MstrnX radiation process'

    if ( trim(ATMOS_TYPE_PHY_RD) .ne. 'MSTRNX' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_RD is not MSTRNX. Check!'
       call PRC_MPIstop
    end if

    return
  end subroutine ATMOS_PHY_RD_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine ATMOS_PHY_RD

end module mod_atmos_phy_rd
