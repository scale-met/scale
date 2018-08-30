#include "scalelib.h"
module scale_spnudge
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  logical, public :: SPNUDGE_uv = .false.
  logical, public :: SPNUDGE_uv_divfree = .true.
  integer, public :: SPNUDGE_uv_lm = 3
  integer, public :: SPNUDGE_uv_mm = 3 

  public :: SPNUDGE_setup
  
  contains
  
  subroutine SPNUDGE_setup
    implicit none
    
    namelist /PARAM_SPNUDGE/ &
      SPNUDGE_uv, &
      SPNUDGE_uv_divfree, &
      SPNUDGE_uv_lm, &
      SPNUDGE_uv_mm

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SPNUDGE_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SPNUDGE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SPNUDGE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SPNUDGE_setup",*) 'Not appropriate names in namelist PARAM_SPNUDGE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SPNUDGE)

  end subroutine SPNUDGE_setup
  
end module scale_spnudge
