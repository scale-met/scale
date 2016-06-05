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
  use scale_debug
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_ADMIN_setup
  public :: LAND_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: LAND_do   = .true. ! main switch for the model

  character(len=H_SHORT), public :: LAND_TYPE = 'NONE'

  logical,                public :: LAND_sw

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
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_LAND / &
       LAND_do,  &
       LAND_TYPE

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
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Land model components ***'

    if ( LAND_TYPE == 'OFF' .OR. LAND_TYPE == 'NONE' ) then
       LAND_do = .false. ! force off
    endif

    if ( LAND_do ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Land  model     : ON'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Land  model     : OFF'
    endif

    if ( LAND_TYPE /= 'OFF' .AND. LAND_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** + Land  physics : ON, ', trim(LAND_TYPE)
       LAND_sw = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** + Land  physics : OFF'
       LAND_sw = .false.
    endif

    return
  end subroutine LAND_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine LAND_ADMIN_getscheme( &
       scheme_name     )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(out) :: scheme_name
    !---------------------------------------------------------------------------

    scheme_name = LAND_TYPE

    return
  end subroutine LAND_ADMIN_getscheme

end module mod_land_admin
