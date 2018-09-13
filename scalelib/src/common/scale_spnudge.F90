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
  logical, public :: SPNUDGE_uv_divfree = .false.
  integer, public :: SPNUDGE_uv_lm = 3
  integer, public :: SPNUDGE_uv_mm = 3 
  real(RP), public :: SPNUDGE_uv_tau

  real(RP), public :: SPNUDGE_uv_alpha

  logical, public :: SPNUDGE_pt = .false.
  integer, public :: SPNUDGE_pt_lm = 3
  integer, public :: SPNUDGE_pt_mm = 3 
  real(RP), public :: SPNUDGE_pt_tau

  real(RP), public :: SPNUDGE_pt_alpha

  
  public :: SPNUDGE_setup
  
  contains
  
  subroutine SPNUDGE_setup
    implicit none
    
    namelist /PARAM_SPNUDGE/ &
      SPNUDGE_uv, &
      SPNUDGE_uv_divfree, &
      SPNUDGE_uv_lm, &
      SPNUDGE_uv_mm, &
      SPNUDGE_uv_tau, &
      SPNUDGE_pt, &
      SPNUDGE_pt_lm, &
      SPNUDGE_pt_mm, &
      SPNUDGE_pt_tau

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

    if( SPNUDGE_uv_tau <= 0.0_RP ) then
        SPNUDGE_uv_alpha = 0
    else
        SPNUDGE_uv_alpha = 1.0_RP / SPNUDGE_uv_tau
    endif

    if( SPNUDGE_pt_tau <= 0.0_RP ) then
        SPNUDGE_pt_alpha = 0
    else
        SPNUDGE_pt_alpha = 1.0_RP / SPNUDGE_pt_tau
    endif

    
  end subroutine SPNUDGE_setup
  
end module scale_spnudge
