!-------------------------------------------------------------------------------
!> module CONVERT driver
!!
!! @par Description
!!          administrator of convert tools
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_convert
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CONVERT_setup
  public :: CONVERT

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, public            :: CONVERT_TYPE = -1
  integer, public, parameter :: I_IGNORE     =  0
  integer, public, parameter :: I_CNVTOPO    =  1
  integer, public, parameter :: I_CNVLANDUSE =  2

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CONVERT_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_cnvtopo, only: &
       CNVTOPO_setup
    use mod_cnvlanduse, only: &
       CNVLANDUSE_setup
    implicit none

    character(len=H_SHORT) :: CONVERT_name = 'NONE'

    NAMELIST / PARAM_CONVERT / &
       CONVERT_name

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CNV BOUNDARY]/Categ[preprocess]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CONVERT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CONVERT. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CONVERT)

    select case(CONVERT_name)
    case('NONE')
       CONVERT_TYPE = I_IGNORE

    case('CNVTOPO')
       CONVERT_TYPE = I_CNVTOPO
       call CNVTOPO_setup

    case('CNVLANDUSE')
       CONVERT_TYPE = I_CNVLANDUSE
       call CNVLANDUSE_setup

    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(CONVERT_name)
       call PRC_MPIstop
    endselect

    return
  end subroutine CONVERT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CONVERT
    use scale_process, only: &
       PRC_MPIstop
    use mod_cnvtopo, only: &
       CNVTOPO
    use mod_cnvlanduse, only: &
       CNVLANDUSE
    implicit none
    !---------------------------------------------------------------------------

    if ( CONVERT_TYPE == I_IGNORE ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  CONVERT BOUNDARY DATA ++++++'
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ START CONVERT BOUNDARY DATA ++++++'

       select case(CONVERT_TYPE)
       case(I_CNVTOPO)
          call CNVTOPO

       case(I_CNVLANDUSE)
          call CNVLANDUSE

       case default
          write(*,*) ' xxx Unsupported TYPE:', CONVERT_TYPE
          call PRC_MPIstop
       endselect

       if( IO_L ) write(IO_FID_LOG,*) '++++++ END   CONVERT BOUNDARY DATA ++++++'
    endif

    return
  end subroutine CONVERT

end module mod_convert
