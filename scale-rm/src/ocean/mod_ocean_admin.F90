!-------------------------------------------------------------------------------
!> module Ocean admin
!!
!! @par Description
!!          Ocean submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_ocean_admin
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_ADMIN_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: OCEAN_do   = .true. ! main switch for the model

  character(len=H_SHORT), public :: OCEAN_DYN_TYPE = 'NONE'
  character(len=H_SHORT), public :: OCEAN_SFC_TYPE = 'FIXED-TEMP'
  character(len=H_SHORT), public :: OCEAN_ALB_TYPE = 'NAKAJIMA00'
  character(len=H_SHORT), public :: OCEAN_RGN_TYPE = 'MOON07'

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
  subroutine OCEAN_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    NAMELIST / PARAM_OCEAN / &
       OCEAN_DYN_TYPE, &
       OCEAN_SFC_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[ADMIN] / Categ[OCEAN] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_ADMIN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_ADMIN_setup",*) 'Not appropriate names in namelist PARAM_OCEAN. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_OCEAN)

    !-----< module component check >-----

    LOG_NEWLINE
    LOG_INFO("OCEAN_ADMIN_setup",*) 'Ocean model components '

    if ( OCEAN_DYN_TYPE /= 'OFF' .AND. OCEAN_DYN_TYPE /= 'NONE' ) then
       LOG_INFO("OCEAN_ADMIN_setup",*) 'Ocean model     : ON, ', trim(OCEAN_DYN_TYPE)
       OCEAN_do = .true.
    else
       LOG_INFO("OCEAN_ADMIN_setup",*) 'Ocean model     : OFF'
       OCEAN_do = .false.
    endif

    if ( OCEAN_do ) then

       LOG_INFO_CONT(*) '+ Ocean surface   model : ', trim(OCEAN_SFC_TYPE)
       LOG_INFO_CONT(*) '+ Ocean albedo    model : ', trim(OCEAN_ALB_TYPE)
       LOG_INFO_CONT(*) '+ Ocean roughness model : ', trim(OCEAN_RGN_TYPE)

    end if

    return
  end subroutine OCEAN_ADMIN_setup

end module mod_ocean_admin
