!-------------------------------------------------------------------------------
!> module Coupler admin
!!
!! @par Description
!!          Coupler submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_cpl_admin
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
  public :: CPL_ADMIN_setup
  public :: CPL_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: CPL_do   = .true. ! main switch for the model

  character(len=H_SHORT), public :: CPL_TYPE_AtmOcn = 'NONE'         !< atmos-ocean coupler type
  character(len=H_SHORT), public :: CPL_TYPE_AtmLnd = 'NONE'         !< atmos-land coupler type
  character(len=H_SHORT), public :: CPL_TYPE_AtmUrb = 'NONE'         !< atmos-urban coupler type
  logical,                public :: CPL_OCN_SFC_TEMP_UPDATE = .true. !< update Ocean Surface Temperature?
  logical,                public :: CPL_LND_SFC_TEMP_UPDATE = .true. !< update Land  Surface Temperature?
  logical,                public :: CPL_URB_SFC_TEMP_UPDATE = .true. !< update Urban Surface Temperature?

  logical,                public :: CPL_sw                    !< do coupler calculation?
  logical,                public :: CPL_sw_AtmOcn             !< do atmos-ocean coupler calculation?
  logical,                public :: CPL_sw_AtmLnd             !< do atmos-land  coupler calculation?
  logical,                public :: CPL_sw_AtmUrb             !< do atmos-urban coupler calculation?

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
  subroutine CPL_ADMIN_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_CPL / &
       CPL_do,                  &
       CPL_TYPE_AtmOcn,         &
       CPL_TYPE_AtmLnd,         &
       CPL_TYPE_AtmUrb,         &
       CPL_OCN_SFC_TEMP_UPDATE, &
       CPL_LND_SFC_TEMP_UPDATE, &
       CPL_URB_SFC_TEMP_UPDATE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[CPL] / Origin[SCALE-LES]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CPL)

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler components ***'

    ! Atoms-Ocean Switch
    if ( CPL_TYPE_AtmOcn /= 'OFF' .AND. CPL_TYPE_AtmOcn /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean Coupler : ON, ', trim(CPL_TYPE_AtmOcn)
       CPL_sw_AtmOcn = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Ocean Coupler : OFF'
       CPL_sw_AtmOcn = .false.
    endif

    ! Atoms-Land Switch
    if ( CPL_TYPE_AtmLnd /= 'OFF' .AND. CPL_TYPE_AtmLnd /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land  Coupler : ON, ', trim(CPL_TYPE_AtmLnd)
       CPL_sw_AtmLnd = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Land  Coupler : OFF'
       CPL_sw_AtmLnd = .false.
    endif

    ! Atoms-Urban Switch
    if ( CPL_TYPE_AtmUrb /= 'OFF' .AND. CPL_TYPE_AtmUrb /= 'NONE' ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Urban Coupler : ON, ', trim(CPL_TYPE_AtmUrb)
       CPL_sw_AtmUrb = .true.
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Atmos-Urban Coupler : OFF'
       CPL_sw_AtmUrb = .false.
    endif

    if ( CPL_sw_AtmLnd .OR. CPL_sw_AtmOcn .OR. CPL_sw_AtmUrb ) then
       CPL_sw = .true.
    else
       CPL_sw = .false.
       CPL_do = .false.
    endif

    if ( CPL_do ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Coupler : ON'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Coupler : OFF'
    endif

    return
  end subroutine CPL_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine CPL_ADMIN_getscheme( &
       component_name, &
       scheme_name     )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*),       intent(in)  :: component_name
    character(len=H_SHORT), intent(out) :: scheme_name
    !---------------------------------------------------------------------------

    select case (component_name)
    case ("AtmOcn")
       scheme_name = CPL_TYPE_AtmOcn
    case ("AtmLnd")
       scheme_name = CPL_TYPE_AtmLnd
    case ("AtmUrb")
       scheme_name = CPL_TYPE_AtmUrb
    case default
       write(*,*) 'xxx Unsupported component_name. Check!', trim(component_name)
       call PRC_MPIstop
    end select

    return
  end subroutine CPL_ADMIN_getscheme

end module mod_cpl_admin
