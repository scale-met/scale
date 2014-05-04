!-------------------------------------------------------------------------------
!> module CNV boundary
!!
!! @par Description
!!          subroutines for preparing topography data
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_cnvboundary
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

  use scale_process, only: &
     PRC_MPIstop
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
  integer, public, save      :: CONVERT_TYPE = -1
  integer, public, parameter :: I_IGNORE    =  0

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CONVERT_setup
    implicit none

    character(len=H_SHORT) :: CONVERT_initname = 'OFF'

    NAMELIST / PARAM_CONVERT / &
       CONVERT_initname

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CONVERT)

    select case(trim(CONVERT_initname))
    case default
       write(*,*) ' xxx Unsupported TYPE:', trim(CONVERT_initname)
    endselect

    return
  end subroutine CONVERT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CONVERT
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( CONVERT_TYPE == I_IGNORE ) then
      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  CONVERT BOUNDARY DATA ++++++'
    else
      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '++++++ START CONVERT BOUNDARY DATA ++++++'

      select case(CONVERT_TYPE)
      case default
         write(*,*) ' xxx Unsupported TYPE:', CONVERT_TYPE
      endselect

      if( IO_L ) write(IO_FID_LOG,*) '++++++ END   CONVERT BOUNDARY DATA ++++++'
    endif

    return
  end subroutine CONVERT

end module mod_cnvboundary
