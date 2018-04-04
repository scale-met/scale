!-------------------------------------------------------------------------------
!> module GRID (cartesian) for ocean
!!
!! @par Description
!!          Grid module for cartesian coordinate for ocean
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_grid_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_ocean_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_GRID_CARTESC_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: OCEAN_GRID_CARTESC_CZ (:) !< center coordinate [m]: z, local=global
  real(RP), public, allocatable :: OCEAN_GRID_CARTESC_FZ (:) !< face   coordinate [m]: z, local=global
  real(RP), public, allocatable :: OCEAN_GRID_CARTESC_CDZ(:) !< z-length of control volume [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: OCEAN_GRID_CARTESC_read
  private :: OCEAN_GRID_CARTESC_generate

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: FZ_MAX = 100
  real(RP), private :: FZ(FZ_MAX) ! face coordinate without surface (=0 m)

  character(len=H_LONG) :: OCEAN_GRID_CARTESC_IN_BASENAME  = ''
  logical               :: OCEAN_GRID_CARTESC_IN_AGGREGATE

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_GRID_CARTESC_setup
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_AGGREGATE
    implicit none

    namelist / PARAM_OCEAN_GRID_CARTESC / &
       OCEAN_GRID_CARTESC_IN_BASENAME,  &
       OCEAN_GRID_CARTESC_IN_AGGREGATE, &
       FZ

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CartesC] / Categ[OCEAN GRID] / Origin[SCALElib]'

    FZ(:) = 0.0_RP

    OCEAN_GRID_CARTESC_IN_AGGREGATE = FILE_AGGREGATE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_GRID_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_GRID_CARTESC. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_OCEAN_GRID_CARTESC)

    allocate( OCEAN_GRID_CARTESC_CZ (OKS  :OKE) )
    allocate( OCEAN_GRID_CARTESC_FZ (OKS-1:OKE) )
    allocate( OCEAN_GRID_CARTESC_CDZ(OKS  :OKE) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean grid information ***'

    if ( OCEAN_GRID_CARTESC_IN_BASENAME /= '' ) then
       call OCEAN_GRID_CARTESC_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found input grid file. Grid position is calculated.'

       call OCEAN_GRID_CARTESC_generate
    endif

    if ( OKE == OKS ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Single layer. ODZ = ', OCEAN_GRID_CARTESC_CDZ(1)
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|====== Vertical Coordinate ======|'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|   k       z      zh      dz   k |'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|         [m]     [m]     [m]     |'
       k = OKS-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',OCEAN_GRID_CARTESC_FZ(k),'        ',k,' | Atmosphere interface'
       do k = OKS, OKE-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,OCEAN_GRID_CARTESC_CZ(k),'        ',OCEAN_GRID_CARTESC_CDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',OCEAN_GRID_CARTESC_FZ(k),'       |',k,' | '
       enddo
       k = OKE
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,OCEAN_GRID_CARTESC_CZ(k),'        ',OCEAN_GRID_CARTESC_CDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',OCEAN_GRID_CARTESC_FZ(k),'        ',k,' | layer of no motion'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|=================================|'
    endif

    return
  end subroutine OCEAN_GRID_CARTESC_setup

  !-----------------------------------------------------------------------------
  !> Read ocean grid
  subroutine OCEAN_GRID_CARTESC_read
    use scale_file, only: &
       FILE_open, &
       FILE_read
    use scale_prc, only: &
       PRC_myrank
    implicit none

    integer  :: fid
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input ocean grid file ***'

    call FILE_open( OCEAN_GRID_CARTESC_IN_BASENAME, fid, rankid=PRC_myrank, aggregate=OCEAN_GRID_CARTESC_IN_AGGREGATE )

    call FILE_read( fid, 'OCZ',  OCEAN_GRID_CARTESC_CZ (:) )
    call FILE_read( fid, 'OCDZ', OCEAN_GRID_CARTESC_CDZ(:) )
    call FILE_read( fid, 'OFZ',  OCEAN_GRID_CARTESC_FZ (:) )

    return
  end subroutine OCEAN_GRID_CARTESC_read

  !-----------------------------------------------------------------------------
  !> Generate ocean grid
  ! It uses FZ, not DZ. Note, LAND_GRID_CARTESC_generate uses DZ
  subroutine OCEAN_GRID_CARTESC_generate
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    if ( OKA == 1 .and. FZ(1) == 0.0_RP ) then
       FZ(1) = 10.0_RP ! to avoid zero thickness (tentative)
    end if
    OCEAN_GRID_CARTESC_FZ(OKS-1) = 0.0_RP
    do k = OKS, OKE
       OCEAN_GRID_CARTESC_FZ(k) = FZ(k)
    enddo

    do k = OKS, OKE
       OCEAN_GRID_CARTESC_CDZ(k) = OCEAN_GRID_CARTESC_FZ(k) - OCEAN_GRID_CARTESC_FZ(k-1)
       OCEAN_GRID_CARTESC_CZ (k) = OCEAN_GRID_CARTESC_CDZ(k) / 2.0_RP + OCEAN_GRID_CARTESC_FZ(k-1)
    enddo

    return
  end subroutine OCEAN_GRID_CARTESC_generate

end module scale_ocean_grid_cartesC

