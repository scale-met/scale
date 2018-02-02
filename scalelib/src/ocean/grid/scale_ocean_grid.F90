!-------------------------------------------------------------------------------
!> module GRID (cartesian) for ocean
!!
!! @par Description
!!          Grid module for cartesian coordinate for ocean
!!
!! @author Team SCALE
!!
!! @par History
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_grid
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_ocean_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_GRID_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: GRID_OCZ (:) !< center coordinate [m]: z, local=global
  real(RP), public, allocatable :: GRID_OFZ (:) !< face   coordinate [m]: z, local=global
  real(RP), public, allocatable :: GRID_OCDZ(:) !< z-length of control volume [m]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: OCEAN_GRID_read
  private :: OCEAN_GRID_generate

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: FZ(100) ! face coordinate without surface (=0 m)

  character(len=H_LONG), private :: OCEAN_GRID_IN_BASENAME  = ''
  character(len=H_LONG), private :: OCEAN_GRID_OUT_BASENAME = ''

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_GRID_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_OCEAN_GRID / &
       OCEAN_GRID_IN_BASENAME,  &
       OCEAN_GRID_OUT_BASENAME, &
       FZ

    integer :: ierr
    integer :: k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[GRID] / Categ[OCEAN GRID] / Origin[SCALElib]'

    FZ(:) = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_GRID,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_GRID. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_OCEAN_GRID)

    allocate( GRID_OCZ (OKS  :OKE) )
    allocate( GRID_OFZ (OKS-1:OKE) )
    allocate( GRID_OCDZ(OKS  :OKE) )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean grid information ***'

    if ( OCEAN_GRID_IN_BASENAME /= '' ) then
       call OCEAN_GRID_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found input grid file. Grid position is calculated.'

       call OCEAN_GRID_generate
    endif

    if ( OKE == OKS ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Single layer. ODZ = ', GRID_OCDZ(1)
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
       '|            ',GRID_OFZ(k),'        ',k,' | Atmosphere interface'
       do k = OKS, OKE-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,GRID_OCZ(k),'        ',GRID_OCDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_OFZ(k),'       |',k,' | '
       enddo
       k = OKE
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I4,F8.3,A,F8.3,A)') &
       '|',k,GRID_OCZ(k),'        ',GRID_OCDZ(k),'     | '
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F8.3,A,I4,A)') &
       '|            ',GRID_OFZ(k),'        ',k,' | layer of no motion'
       if( IO_L ) write(IO_FID_LOG,'(1x,A)') &
       '|=================================|'
    endif

    return
  end subroutine OCEAN_GRID_setup

  !-----------------------------------------------------------------------------
  !> Read ocean grid
  subroutine OCEAN_GRID_read
    use scale_file, only: &
       FILE_Read
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_grid, only: &
       GRID_CBFZ, &
       GRID_CBFX, &
       GRID_CBFY
    implicit none

    character(len=H_LONG) :: bname
    real(RP)              :: tmp_CBFZ(KA)
    real(RP)              :: tmp_CBFX(IA)
    real(RP)              :: tmp_CBFY(JA)

    integer  :: i, j, k
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input ocean grid file ***'

    write(bname,'(A,A,F15.3)') trim(OCEAN_GRID_IN_BASENAME)

    call FILE_Read( bname, 'OCZ',  GRID_OCZ (:) )
    call FILE_Read( bname, 'OCDZ', GRID_OCDZ(:) )
    call FILE_Read( bname, 'OFZ',  GRID_OFZ (:) )
                                                 
    call FILE_Read( bname, 'CBFZ', tmp_CBFZ (:) )
    call FILE_Read( bname, 'CBFX', tmp_CBFY (:) )
    call FILE_Read( bname, 'CBFY', tmp_CBFY (:) )

    do i = 1, IA
       if ( tmp_CBFX(i) /= GRID_CBFX(i) ) then
          write(*,*) 'xxx Buffer layer in OCEAN_GRID_IN_BASENAME is different from GRID_IN_BASENAME'
          call PRC_MPIstop
       endif
    enddo

    do j = 1, JA
       if ( tmp_CBFY(j) /= GRID_CBFY(j) ) then
          write(*,*) 'xxx Buffer layer in OCEAN_GRID_IN_BASENAME is different from GRID_IN_BASENAME'
          call PRC_MPIstop
       endif
    enddo

    do k = 1, KA
       if ( tmp_CBFZ(k) /= GRID_CBFZ(k) ) then
          write(*,*) 'xxx Buffer layer in OCEAN_GRID_IN_BASENAME is different from GRID_IN_BASENAME'
          call PRC_MPIstop
       endif
    enddo

    return
  end subroutine OCEAN_GRID_read

  !-----------------------------------------------------------------------------
  !> Generate ocean grid
  ! It uses OFZ, not ODZ. Note, LAND_GRID_generate uses LDZ
  subroutine OCEAN_GRID_generate
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    GRID_OFZ(OKS-1) = 0.0_RP
    do k = OKS, OKE
       GRID_OFZ(k) = FZ(k)
    enddo

    do k = OKS, OKE
       GRID_OCDZ(k) = GRID_OFZ(k) - GRID_OFZ(k-1)
       GRID_OCZ (k) = GRID_OCDZ(k) / 2.0_RP + GRID_OFZ(k-1)
    enddo

    return
  end subroutine OCEAN_GRID_generate

end module scale_ocean_grid

