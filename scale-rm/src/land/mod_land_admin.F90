!-------------------------------------------------------------------------------
!> module Land admin
!!
!! @par Description
!!          Land submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_land_admin
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_ADMIN_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: LAND_DYN_TYPE = 'NONE'
                                                  ! 'OFF'
                                                  ! 'BUCKET'
                                                  ! 'INIT'
  character(len=H_SHORT), public :: LAND_SFC_TYPE = 'SKIN'
                                                  ! 'FIXED-TEMP'
  character(len=H_SHORT), public :: SNOW_TYPE     = 'NONE'
                                                  ! 'OFF'
                                                  ! 'KY90'
  logical,                public :: LAND_do
  logical,                public :: SNOW_sw

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
  subroutine LAND_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    use scale_land_grid_cartesC_index, only: &
       LKMAX
    implicit none

    namelist / PARAM_LAND / &
       LAND_DYN_TYPE, &
       LAND_SFC_TYPE, &
       SNOW_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAND_ADMIN_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("LAND_ADMIN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("LAND_ADMIN_setup",*) 'Not appropriate names in namelist PARAM_LAND. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_LAND)

    !-----< module component check >-----

    LOG_NEWLINE
    LOG_INFO("LAND_ADMIN_setup",*) 'Land model components '

    if ( LAND_DYN_TYPE /= 'OFF' .AND. LAND_DYN_TYPE /= 'NONE' ) then
       if ( LKMAX < 0 ) then
          LOG_ERROR("LAND_ADMIN_setup",*) 'LAND_DYN_TYPE is set but LKMAX < 0'
          call PRC_abort
       end if
       LOG_INFO_CONT(*) 'Land model           : ON, ', trim(LAND_DYN_TYPE)
       LAND_do = .true.
    else
       LOG_INFO_CONT(*) 'Land model           : OFF'
       LAND_do = .false.
    endif

    if ( LAND_do ) then

       if ( SNOW_TYPE /= 'OFF' .AND. SNOW_TYPE /= 'NONE' ) then
          LOG_INFO_CONT(*) '+ Snow  physics      : ON, ', trim(SNOW_TYPE)
          SNOW_sw = .true.
       else
          LOG_INFO_CONT(*) '+ Snow  physics      : OFF'
          SNOW_sw = .false.
       endif

       LOG_INFO_CONT(*) '+ Land surface model : ', trim(LAND_SFC_TYPE)

    end if

    return
  end subroutine LAND_ADMIN_setup

end module mod_land_admin
