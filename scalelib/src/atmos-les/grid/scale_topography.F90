!-------------------------------------------------------------------------------
!> module TOPOGRAPHY
!!
!! @par Description
!!          Topography module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_topography
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TOPO_setup
  public :: TOPO_fillhalo
  public :: TOPO_write

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: TOPO_exist = .false. !< topography exists?

  real(RP), public, allocatable :: TOPO_Zsfc(:,:) !< absolute ground height [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TOPO_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG), private :: TOPO_IN_BASENAME  = ''                     !< basename of the input  file
  character(len=H_LONG), private :: TOPO_OUT_BASENAME = ''                     !< basename of the output file
  character(len=H_MID),  private :: TOPO_OUT_TITLE    = 'SCALE-LES TOPOGRAPHY' !< title    of the output file
  character(len=H_MID),  private :: TOPO_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TOPO_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_TOPO / &
       TOPO_IN_BASENAME,  &
       TOPO_OUT_BASENAME, &
       TOPO_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[TOPOGRAPHY] / Categ[ATMOS-LES GRID] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_TOPO. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_TOPO)

    allocate( TOPO_Zsfc(IA,JA) )
    TOPO_Zsfc(:,:) = 0.0_RP

    ! read from file
    call TOPO_read

    return
  end subroutine TOPO_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine TOPO_fillhalo( Zsfc )
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    real(RP), intent(inout), optional :: Zsfc(IA,JA)
    !---------------------------------------------------------------------------

    if ( present(Zsfc) ) then
       call COMM_vars8( Zsfc(:,:), 1 )
       call COMM_wait ( Zsfc(:,:), 1 )
    else
       call COMM_vars8( TOPO_Zsfc(:,:), 1 )
       call COMM_wait ( TOPO_Zsfc(:,:), 1 )
    end if

    return
  end subroutine TOPO_fillhalo

  !-----------------------------------------------------------------------------
  !> Read topography
  subroutine TOPO_read
    use gtool_file, only: &
       FileRead
    use scale_fileio, only: &
       FILEIO_read
    use scale_process, only: &
       PRC_MPIstop, &
       PRC_myrank
    use scale_grid, only: &
       GRID_CX, &
       GRID_CY

    implicit none

    real(RP) :: tmp_CX(IA), tmp_CY(JA)
    real(RP) :: epsilon
    integer  :: i, j
    !---------------------------------------------------------------------------

    epsilon = 0.1_RP**(RP-1)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input topography file ***'

    if ( TOPO_IN_BASENAME /= '' ) then

       call FILEIO_read( TOPO_Zsfc(:,:),                        & ! [OUT]
                         TOPO_IN_BASENAME, 'TOPO', 'XY', step=1 ) ! [IN]

       call TOPO_fillhalo

       TOPO_exist = .true.

       call FileRead( tmp_CX(:),  TOPO_IN_BASENAME, 'CX', 1, PRC_myrank )
       call FileRead( tmp_CY(:),  TOPO_IN_BASENAME, 'CY', 1, PRC_myrank )

       do i = 1, IA
         if( abs(tmp_CX(i) - GRID_CX(i)) > epsilon ) then
           write( IO_FID_LOG,'(A)')  '*** X position in TOPO_IN_BASENAME is different from GRID_IN_BASENAME ***'
           write( IO_FID_LOG,* )  "I", i, tmp_CX(i), GRID_CX(i)
           call PRC_MPIstop
         endif
       enddo
       do j = 1, JA
         if( abs(tmp_CY(j) - GRID_CY(j)) > epsilon ) then
           write( IO_FID_LOG,'(A)')  '*** Y position in TOPO_IN_BASENAME is different from GRID_IN_BASENAME ***'
           write( IO_FID_LOG,* )  "J", j, tmp_CY(j), GRID_CY(j)
           call PRC_MPIstop
         endif
       enddo

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** topography file is not specified.'

       TOPO_exist = .false.
    endif

    return
  end subroutine TOPO_read

  !-----------------------------------------------------------------------------
  !> Write topography
  subroutine TOPO_write
    use scale_fileio, only: &
       FILEIO_write
    implicit none
    !---------------------------------------------------------------------------

    if ( TOPO_OUT_BASENAME /= '' ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Output topography file ***'

       call FILEIO_write( TOPO_Zsfc(:,:), TOPO_OUT_BASENAME, TOPO_OUT_TITLE, & ! [IN]
                          'TOPO', 'Topography', 'm', 'XY',   TOPO_OUT_DTYPE, & ! [IN]
                          nozcoord=.true. )

    endif

    return
  end subroutine TOPO_write

end module scale_topography
