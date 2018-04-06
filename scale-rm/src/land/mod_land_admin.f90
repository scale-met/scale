!-------------------------------------------------------------------------------
!> module Land admin
!!
!! @par Description
!!          Land submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_land_admin
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
  public :: LAND_ADMIN_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: LAND_DYN_TYPE = 'NONE'
  character(len=H_SHORT), public :: LAND_SFC_TYPE = 'SKIN'
  character(len=H_SHORT), public :: SNOW_TYPE     = 'NONE'

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
    implicit none

    NAMELIST / PARAM_LAND / &
       LAND_DYN_TYPE, &
       LAND_SFC_TYPE, &
       SNOW_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[LAND] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Land model components ***'

    if ( LAND_DYN_TYPE /= 'OFF' .AND. LAND_DYN_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** + Land model : ON, ', trim(LAND_DYN_TYPE)
       LAND_do = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** + Land model : OFF'
       LAND_do = .false.
    endif

    if ( LAND_do ) then

       if ( SNOW_TYPE /= 'OFF' .AND. SNOW_TYPE /= 'NONE' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** + Snow  physics : ON, ', trim(SNOW_TYPE)
          SNOW_sw = .true.
       else
          if( IO_L ) write(IO_FID_LOG,*) '*** + Snow  physics : OFF'
          SNOW_sw = .false.
       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** + Land surface model : ', trim(LAND_SFC_TYPE)

    end if

    return
  end subroutine LAND_ADMIN_setup

end module mod_land_admin
