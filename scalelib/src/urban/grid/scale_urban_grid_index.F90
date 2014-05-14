!-------------------------------------------------------------------------------
!> module grid index for urban
!!
!! @par Description
!!          Grid Index module for urban
!!
!! @author Team SCALE
!!
!! @par History
!!
!<
!-------------------------------------------------------------------------------
module scale_urban_grid_index
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
  public :: URBAN_GRID_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
#ifdef FIXED_INDEX
  include "inc_index.h"
#else
  integer, public :: UKMAX = 1 ! # of computational cells: z for urban

  integer, public :: UKS ! start point of inner domain: z for urban, local
  integer, public :: UKE ! end   point of inner domain: z for urban, local
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
  subroutine URBAN_GRID_INDEX_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

#ifndef FIXED_INDEX
    namelist / PARAM_URBAN_INDEX / &
       UKMAX
#endif

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[URBAN_INDEX]/Categ[COMMON]'

#ifndef FIXED_INDEX
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_INDEX,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_INDEX. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_URBAN_INDEX)

    UKS  = 1
    UKE  = UKMAX
#endif

  end subroutine URBAN_GRID_INDEX_setup

end module scale_urban_grid_index
