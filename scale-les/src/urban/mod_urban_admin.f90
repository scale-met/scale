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
  use scale_debug
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_ADMIN_setup
  public :: URBAN_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: URBAN_do   = .true. ! main switch for the model

  character(len=H_SHORT), public :: URBAN_TYPE = 'NONE'

  logical,                public :: URBAN_sw

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
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_URBAN / &
       URBAN_do,  &
       URBAN_TYPE
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[URBAN] / Origin[SCALE-LES]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_URBAN)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Urban model components ***'

    if ( URBAN_TYPE == 'OFF' .OR. URBAN_TYPE == 'NONE' ) then
       URBAN_do = .false. ! force off
    endif

    if ( URBAN_do ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Urban model : ON'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Urban model : OFF'
    endif

    if ( URBAN_TYPE /= 'OFF' .AND. URBAN_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Urban physics : ON, ', trim(URBAN_TYPE)
       URBAN_sw = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Urban physics : OFF'
       URBAN_sw = .false.
    endif

    return
  end subroutine URBAN_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine URBAN_ADMIN_getscheme( &
       scheme_name     )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(out) :: scheme_name
    !---------------------------------------------------------------------------

    scheme_name = URBAN_TYPE

    return
  end subroutine URBAN_ADMIN_getscheme

end module mod_urban_admin
