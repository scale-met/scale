!-------------------------------------------------------------------------------
!> module Ocean admin
!!
!! @par Description
!!          Ocean submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_ocean_admin
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
  public :: OCEAN_ADMIN_setup
  public :: OCEAN_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: OCEAN_do   = .true. ! main switch for the model

  character(len=H_SHORT), public :: OCEAN_TYPE = 'NONE'

  logical,                public :: OCEAN_sw

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
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_OCEAN / &
       OCEAN_do,  &
       OCEAN_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[OCEAN] / Origin[SCALE-LES]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean model components ***'

    if ( OCEAN_TYPE == 'OFF' .OR. OCEAN_TYPE == 'NONE' ) then
       OCEAN_do = .false. ! force off
    endif

    if ( OCEAN_do ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocean model : ON'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Ocean model : OFF'
    endif

    if ( OCEAN_TYPE /= 'OFF' .AND. OCEAN_TYPE /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** +Ocean physics : ON, ', trim(OCEAN_TYPE)
       OCEAN_sw = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** +Ocean physics : OFF'
       OCEAN_sw = .false.
    endif

    return
  end subroutine OCEAN_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine OCEAN_ADMIN_getscheme( &
       scheme_name     )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(out) :: scheme_name
    !---------------------------------------------------------------------------

    scheme_name = OCEAN_TYPE

    return
  end subroutine OCEAN_ADMIN_getscheme

end module mod_ocean_admin
