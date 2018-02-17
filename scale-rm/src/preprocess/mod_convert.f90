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
  use scale_atmos_grid_cartesC_index
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
  logical :: CONVERT_TOPO    = .false.
  logical :: CONVERT_LANDUSE = .false.
  logical :: CONVERT_2D      = .false.

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
    use mod_cnv2d, only: &
       CNV2D_setup
    implicit none

    NAMELIST / PARAM_CONVERT / &
       CONVERT_TOPO,    &
       CONVERT_LANDUSE, &
       CONVERT_2D

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[convert] / Categ[preprocess] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CONVERT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CONVERT. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_CONVERT)

    ! set up TOPO
    if( CONVERT_TOPO ) then
       call CNVTOPO_setup
    end if

    ! set up LANDUSE
    if( CONVERT_LANDUSE ) then
       call CNVLANDUSE_setup
    end if

    ! set up LANDUSE
    if( CONVERT_2D ) then
       call CNV2D_setup
    end if

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
    use mod_cnv2d, only: &
       CNV2D
    implicit none
    !---------------------------------------------------------------------------

    if ( CONVERT_TOPO .OR. CONVERT_LANDUSE .OR. CONVERT_2D ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ START CONVERT BOUNDARY DATA ++++++'

       if( CONVERT_TOPO ) then
          call CNVTOPO
       end if

       if( CONVERT_LANDUSE ) then
          call CNVLANDUSE
       end if

       if( CONVERT_2D ) then
          call CNV2D
       end if

       if( IO_L ) write(IO_FID_LOG,*) '++++++ END   CONVERT BOUNDARY DATA ++++++'
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '++++++ SKIP  CONVERT BOUNDARY DATA ++++++'
    endif

    return
  end subroutine CONVERT

end module mod_convert
