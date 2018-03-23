!-------------------------------------------------------------------------------
!> module Urban admin
!!
!! @par Description
!!          Urban submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_urban_admin
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
  public :: URBAN_ADMIN_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: URBAN_do   = .true. ! main switch for the model

  character(len=H_SHORT), public :: URBAN_DYN_TYPE = 'NONE'

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
  subroutine URBAN_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    NAMELIST / PARAM_URBAN / &
       URBAN_DYN_TYPE
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[URBAN] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_URBAN)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Urban model components ***'

    if ( URBAN_DYN_TYPE /= 'OFF' .AND. URBAN_DYN_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** + Urban model     : ON, ', trim(URBAN_DYN_TYPE)
       URBAN_do = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** + Urban model     : OFF'
       URBAN_do = .false.
    endif

    return
  end subroutine URBAN_ADMIN_setup

end module mod_urban_admin
