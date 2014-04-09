!-------------------------------------------------------------------------------
!> module grid index for land
!!
!! @par Description
!!          Grid Index module for land
!!
!! @author Team SCALE
!!
!! @par History
!!
!<
!-------------------------------------------------------------------------------
module scale_land_grid_index
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_GRID_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
#ifdef FIXED_INDEX
  include "inc_index.h"
#else
  integer, public :: LKMAX = 1 ! # of computational cells: z for land

  integer, public :: LKS ! start point of inner domain: z for land, local
  integer, public :: LKE ! end   point of inner domain: z for land, local
#endif

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
  subroutine LAND_GRID_INDEX_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

#ifndef FIXED_INDEX
    namelist / PARAM_LAND_INDEX / &
       LKMAX
#endif

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[LAND_INDEX]/Categ[COMMON]'

#ifndef FIXED_INDEX
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_INDEX,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_INDEX. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND_INDEX)

    LKS  = 1
    LKE  = LKMAX
#endif

  end subroutine LAND_GRID_INDEX_setup

end module scale_land_grid_index
