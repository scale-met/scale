!-------------------------------------------------------------------------------
!> module TOPOGRAPHY
!!
!! @par Description
!!          Topography module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_topography
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TOPO_setup
  public :: TOPO_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: TOPO_Zsfc(:,:)   !< absolute ground height [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TOPO_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_FILECHR), private :: TOPO_IN_BASENAME  = ''                     !< basename of the input  file
  character(len=IO_FILECHR), private :: TOPO_OUT_BASENAME = ''                     !< basename of the output file
  character(len=IO_SYSCHR),  private :: TOPO_OUT_TITLE    = 'SCALE-LES TOPOGRAPHY' !< title    of the output file
  character(len=IO_SYSCHR),  private :: TOPO_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TOPO_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_TOPO / &
       TOPO_IN_BASENAME,  &
       TOPO_OUT_BASENAME, &
       TOPO_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TOPOGRAPHY]/Categ[GRID]'

    allocate( TOPO_Zsfc(IA,JA) )
    TOPO_Zsfc(:,:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPO,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_TOPO. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_TOPO)

    ! read from file
    call TOPO_read

    return
  end subroutine TOPO_setup

  !-----------------------------------------------------------------------------
  !> Read topography
  subroutine TOPO_read
    use mod_fileio, only: &
       FILEIO_read
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input topography file ***'

    if ( TOPO_IN_BASENAME /= '' ) then

       call FILEIO_read( TOPO_Zsfc(:,:),                        & ! [OUT]
                         TOPO_IN_BASENAME, 'TOPO', 'XY', step=1 ) ! [IN]

       ! fill IHALO & JHALO
       call COMM_vars8( TOPO_Zsfc(:,:), 1 )
       call COMM_wait ( TOPO_Zsfc(:,:), 1 )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** topography file is not specified.'
    endif

    return
  end subroutine TOPO_read

  !-----------------------------------------------------------------------------
  !> Write topography
  subroutine TOPO_write
    use mod_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( TOPO_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output topography file ***'

       call FILEIO_write( TOPO_Zsfc(:,:), TOPO_OUT_BASENAME, TOPO_OUT_TITLE, & ! [IN]
                          'TOPO', 'Topography', 'm', 'XY', TOPO_OUT_DTYPE    ) ! [IN]

    endif

    return
  end subroutine TOPO_write

end module mod_topography
