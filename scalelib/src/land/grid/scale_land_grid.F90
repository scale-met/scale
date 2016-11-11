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
  real(RP), private :: LDZ(100)

  character(len=H_LONG), private :: LAND_GRID_IN_BASENAME  = ''
  character(len=H_LONG), private :: LAND_GRID_OUT_BASENAME = ''

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_GRID_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_LAND_GRID / &
       LAND_GRID_IN_BASENAME,  &
       LAND_GRID_OUT_BASENAME, &
       LDZ

    integer :: ierr
    integer :: k

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID] / Categ[LAND GRID] / Origin[SCALElib]'

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
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_LAND_GRID)

    allocate( GRID_LCZ (LKS  :LKE) )
    allocate( GRID_LFZ (LKS-1:LKE) )
    allocate( GRID_LCDZ(LKS  :LKE) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Land grid information ***'

    if ( LAND_GRID_IN_BASENAME /= '' ) then
       call LAND_GRID_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found input grid file. Grid position is calculated.'

       call LAND_GRID_generate
    endif

    if ( LKE == LKS ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Single layer. LDZ = ', LDZ(1)
    else
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|====== Vertical Coordinate ======|'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|   k       z      zh      dz   k |'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|         [m]     [m]     [m]     |'
       k = LKS-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_LFZ(k),'        ',k,' | Atmosphere interface'
       do k = LKS, LKE-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,GRID_LCZ(k),'        ',GRID_LCDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_LFZ(k),'       |',k,' | '
       enddo
       k = LKE
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,GRID_LCZ(k),'        ',GRID_LCDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_LFZ(k),'        ',k,' | bedrock'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|=================================|'
    endif

    return
  end subroutine LAND_GRID_setup

  !-----------------------------------------------------------------------------
  !> Read land grid
  subroutine LAND_GRID_read
    use gtool_file, only: &
       FileRead
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_grid, only: &
       GRID_CBFZ, &
       GRID_CBFX, &
       GRID_CBFY
    implicit none

    character(len=H_LONG) :: bname
    real(RP) :: tmp_CBFZ(KA), tmp_CBFX(IA), tmp_CBFY(JA)
    integer :: i, j, k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input land grid file ***'

    write(bname,'(A,A,F15.3)') trim(LAND_GRID_IN_BASENAME)

    call FileRead( GRID_LCZ(:),  bname, 'LCZ',  1, PRC_myrank )
    call FileRead( GRID_LCDZ(:), bname, 'LCDZ', 1, PRC_myrank )
    call FileRead( GRID_LFZ(:),  bname, 'LFZ',  1, PRC_myrank )

    call FileRead( tmp_CBFZ(:),  bname, 'CBFZ', 1, PRC_myrank )
    call FileRead( tmp_CBFY(:),  bname, 'CBFX', 1, PRC_myrank )
    call FileRead( tmp_CBFY(:),  bname, 'CBFY', 1, PRC_myrank )

    do i = 1, IA
       if( tmp_CBFX(i) /= GRID_CBFX(i) ) then
          write(*,*) 'xxx Buffer layer in LAND_GRID_IN_BASENAME is different from GRID_IN_BASENAME'
          call PRC_MPIstop
       endif
    enddo

    do j = 1, JA
       if( tmp_CBFY(j) /= GRID_CBFY(j) ) then
          write(*,*) 'xxx Buffer layer in LAND_GRID_IN_BASENAME is different from GRID_IN_BASENAME'
          call PRC_MPIstop
       endif
    enddo

    do k = 1, KA
       if ( tmp_CBFZ(k) /= GRID_CBFZ(k) ) then
          write(*,*) 'xxx Buffer layer in LAND_GRID_IN_BASENAME is different from GRID_IN_BASENAME'
          call PRC_MPIstop
       endif
    enddo

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
    enddo

    return
  end subroutine LAND_GRID_generate

end module scale_land_grid
