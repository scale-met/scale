!-------------------------------------------------------------------------------
!> module GRID (cartesian) for land
!!
!! @par Description
!!          Grid module for cartesian coordinate for land
!!
!! @author Team SCALE
!!
!! @par History
!!
!<
!-------------------------------------------------------------------------------
module scale_land_grid
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_GRID_allocate
  public :: LAND_GRID_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: GRID_LCZ  (:)  !< center coordinate [m]: z, local=global
  real(RP), public, allocatable :: GRID_LFZ  (:)  !< face   coordinate [m]: z, local=global
  real(RP), public, allocatable :: GRID_LCDZ (:)  !< z-length of control volume [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: LAND_GRID_read
  private :: LAND_GRID_generate

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: LDZ(:)

  character(len=H_LONG), private :: LAND_GRID_IN_BASENAME  = ''
  character(len=H_LONG), private :: LAND_GRID_OUT_BASENAME = ''

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_GRID_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    namelist / PARAM_LAND_GRID / &
       LAND_GRID_IN_BASENAME,  &
       LAND_GRID_OUT_BASENAME, &
       LDZ

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CARTESIAN_PLANE]/Categ[LAND_GRID]'

    allocate( LDZ(LKS:LKE) )

    LDZ(:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_GRID,iostat=ierr)

    if( ierr < 0 ) then !--- missing
      if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_GRID. Check!'
      call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND_GRID)

    call LAND_GRID_allocate

    if ( LAND_GRID_IN_BASENAME /= '' ) then
      call LAND_GRID_read
    else
      if( IO_L ) write(IO_FID_LOG,*) '*** Not found input grid file. Generate!'

      call LAND_GRID_generate
    endif

    return
  end subroutine LAND_GRID_setup

  !-----------------------------------------------------------------------------
  !> Allocate arrays
  subroutine LAND_GRID_allocate
    implicit none
    !---------------------------------------------------------------------------

    ! local domain
    allocate( GRID_LCZ (LKS  :LKE) )
    allocate( GRID_LFZ (LKS-1:LKE) )
    allocate( GRID_LCDZ(LKS  :LKE) )

  end subroutine LAND_GRID_allocate

  !-----------------------------------------------------------------------------
  !> Read land grid
  subroutine LAND_GRID_read
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank
    implicit none

    character(len=H_LONG) :: bname
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input land grid file ***'

    write(bname,'(A,A,F15.3)') trim(LAND_GRID_IN_BASENAME)

    call FileRead( GRID_LCZ(:),  bname, 'LCZ',  1, PRC_myrank )
    call FileRead( GRID_LCDZ(:), bname, 'LCDZ', 1, PRC_myrank )
    call FileRead( GRID_LFZ(:),  bname, 'LFZ',  1, PRC_myrank )

    return
  end subroutine LAND_GRID_read

  !-----------------------------------------------------------------------------
  !> Generate land grid
  subroutine LAND_GRID_generate
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    do k = LKS, LKE
      GRID_LCDZ(k) = LDZ(k)
    enddo

    GRID_LFZ(LKS-1) = 0.0_RP

    do k = LKS, LKE
      GRID_LCZ(k) = GRID_LCDZ(k) / 2.0_RP + GRID_LFZ(k-1)
      GRID_LFZ(k) = GRID_LCDZ(k)          + GRID_LFZ(k-1)
    end do

    return
  end subroutine LAND_GRID_generate

end module scale_land_grid
