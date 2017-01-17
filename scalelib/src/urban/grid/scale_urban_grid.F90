!-------------------------------------------------------------------------------
!> module GRID (cartesian) for urban
!!
!! @par Description
!!          Grid module for cartesian coordinate for urban
!!
!! @author Team SCALE
!!
!! @par History
!!
!<
!-------------------------------------------------------------------------------
module scale_urban_grid
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_GRID_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: GRID_UCZ  (:)  !< center coordinate [m]: z, local=global
  real(RP), public, allocatable :: GRID_UFZ  (:)  !< face   coordinate [m]: z, local=global
  real(RP), public, allocatable :: GRID_UCDZ (:)  !< z-length of control volume [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: URBAN_GRID_read
  private :: URBAN_GRID_generate

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: UDZ(200)

  character(len=H_LONG), private :: URBAN_GRID_IN_BASENAME  = ''
  character(len=H_LONG), private :: URBAN_GRID_OUT_BASENAME = ''

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_GRID_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_URBAN_GRID / &
       URBAN_GRID_IN_BASENAME,  &
       URBAN_GRID_OUT_BASENAME, &
       UDZ

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID] / Categ[URBAN GRID] / Origin[SCALElib]'

    UDZ(:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_GRID,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_GRID. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_URBAN_GRID)

    allocate( GRID_UCZ (UKS  :UKE) )
    allocate( GRID_UFZ (UKS-1:UKE) )
    allocate( GRID_UCDZ(UKS  :UKE) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Urban grid information ***'

    if ( URBAN_GRID_IN_BASENAME /= '' ) then
       call URBAN_GRID_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found input grid file. Grid position is calculated.'

       call URBAN_GRID_generate
    endif

    if ( UKE == UKS ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Single layer. LDZ = ', UDZ(1)
    else
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|====== Vertical Coordinate ======|'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|   k       z      zh      dz   k |'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|         [m]     [m]     [m]     |'
       k = UKS-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_UFZ(k),'        ',k,' | Atmosphere interface'
       do k = UKS, UKE-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,GRID_UCZ(k),'        ',GRID_UCDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_UFZ(k),'       |',k,' | '
       enddo
       k = UKE
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,GRID_UCZ(k),'        ',GRID_UCDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_UFZ(k),'        ',k,' | bedrock'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|=================================|'
    endif

    return
  end subroutine URBAN_GRID_setup

  !-----------------------------------------------------------------------------
  !> Read urban grid
  subroutine URBAN_GRID_read
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    implicit none

    character(len=H_LONG) :: bname
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input urban grid file ***'

    write(bname,'(A,A,F15.3)') trim(URBAN_GRID_IN_BASENAME)

    call FileRead( GRID_UCZ(:),  bname, 'UCZ',  1, PRC_myrank )
    call FileRead( GRID_UCDZ(:), bname, 'UCDZ', 1, PRC_myrank )
    call FileRead( GRID_UFZ(:),  bname, 'UFZ',  1, PRC_myrank )

    return
  end subroutine URBAN_GRID_read

  !-----------------------------------------------------------------------------
  !> Generate urban grid
  subroutine URBAN_GRID_generate
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    do k = UKS, UKE
       GRID_UCDZ(k) = UDZ(k)
    enddo

    GRID_UFZ(UKS-1) = 0.0_RP

    do k = UKS, UKE
       GRID_UCZ(k) = GRID_UCDZ(k) / 2.0_RP + GRID_UFZ(k-1)
       GRID_UFZ(k) = GRID_UCDZ(k)          + GRID_UFZ(k-1)
    enddo

    return
  end subroutine URBAN_GRID_generate

end module scale_urban_grid
