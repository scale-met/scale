!-------------------------------------------------------------------------------
!> module GRID
!!
!! @par Description
!!          Grid module for plane cartesian coordinate
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!!
!<
!-------------------------------------------------------------------------------
module mod_grid
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GRID_setup
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public,         parameter :: GRID_IHALO = 2     ! # of halo cells: x
  integer, public,         parameter :: GRID_JHALO = 2     ! # of halo cells: y
  integer, public,         parameter :: GRID_KHALO = 2     ! # of halo cells: z

  integer, public,              save :: GRID_IMAX = 60     ! # of computational cells: x
  integer, public,              save :: GRID_JMAX = 60     ! # of computational cells: y
  integer, public,              save :: GRID_KMAX = 320    ! # of computational cells: z

  real(8), public,              save :: GRID_DX = 40.D0    ! center/face length [m]: x
  real(8), public,              save :: GRID_DY = 40.D0    ! center/face length [m]: y
  real(8), public,              save :: GRID_DZ = 40.D0    ! layer/interface length [m]: z
  real(8), public,              save :: GRID_RDX           ! inverse DX
  real(8), public,              save :: GRID_RDY           ! inverse DY
  real(8), public,              save :: GRID_RDZ           ! inverse DZ

  integer, public,              save :: GRID_IA            ! # of x whole cells (with HALO)
  integer, public,              save :: GRID_JA            ! # of y whole cells (with HALO)
  integer, public,              save :: GRID_KA            ! # of z whole cells (with HALO)

  integer, public,              save :: GRID_IS, GRID_IE   ! start/end of inner domain: x
  integer, public,              save :: GRID_JS, GRID_JE   ! start/end of inner domain: y
  integer, public,              save :: GRID_KS, GRID_KE   ! start/end of inner domain: z, layer
  integer, public,              save :: GRID_WS, GRID_WE   ! start/end of inner domain: z, interface

  integer, public,              save :: GRID_ISG, GRID_IEG ! start/end of inner domain: x, global
  integer, public,              save :: GRID_JSG, GRID_JEG ! start/end of inner domain: y, global

  real(8), public, allocatable, save :: GRID_CX(:)         ! center coordinate [m]: x
  real(8), public, allocatable, save :: GRID_FX(:)         ! face   coordinate [m]: x
  real(8), public, allocatable, save :: GRID_CY(:)         ! center coordinate [m]: y
  real(8), public, allocatable, save :: GRID_FY(:)         ! face   coordinate [m]: y
  real(8), public, allocatable, save :: GRID_CZ(:)         ! center coordinate [m]: z
  real(8), public, allocatable, save :: GRID_FZ(:)         ! face   coordinate [m]: z

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup horizontal&vertical grid
  !-----------------------------------------------------------------------------
  subroutine GRID_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_2Drank,  &
       PRC_nmax,    &
       PRC_NUM_X,   &
       PRC_NUM_Y
    implicit none

    NAMELIST / PARAM_GRID / &
       GRID_IMAX, &
       GRID_JMAX, &
       GRID_KMAX, &
       GRID_DX,   &
       GRID_DY,   &
       GRID_DZ

    integer :: i, j, k, ii, jj
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[CARTESIAN_PLANE]/Categ[GRID]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_GRID,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_GRID. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_GRID)

    ! domain setting for MPI (divided only i-direction)
    GRID_IA = GRID_IHALO + GRID_IMAX + GRID_IHALO
    GRID_JA = GRID_JHALO + GRID_JMAX + GRID_JHALO
    GRID_KA = GRID_KHALO + GRID_KMAX + GRID_KHALO

    ! horizontal index (local domain)
    GRID_IS = GRID_IHALO + 1
    GRID_IE = GRID_IHALO + GRID_IMAX
    GRID_JS = GRID_JHALO + 1
    GRID_JE = GRID_JHALO + GRID_JMAX

    ! horizontal index (global domain)
    GRID_ISG = GRID_IHALO + 1         + PRC_2Drank(PRC_myrank,1) * GRID_IMAX
    GRID_IEG = GRID_IHALO + GRID_IMAX + PRC_2Drank(PRC_myrank,1) * GRID_IMAX
    GRID_JSG = GRID_JHALO + 1         + PRC_2Drank(PRC_myrank,2) * GRID_JMAX
    GRID_JEG = GRID_JHALO + GRID_JMAX + PRC_2Drank(PRC_myrank,2) * GRID_JMAX

    ! vertical index
    GRID_KS = GRID_KHALO + 1
    GRID_KE = GRID_KHALO + GRID_KMAX
    GRID_WS = GRID_KS - 1
    GRID_WE = GRID_KE

    allocate( GRID_CX(  GRID_IA) )
    allocate( GRID_FX(  GRID_IA) )
    allocate( GRID_CY(  GRID_JA) )
    allocate( GRID_FY(  GRID_JA) )
    allocate( GRID_CZ(  GRID_KA) )
    allocate( GRID_FZ(0:GRID_KA) )

    ! horizontal coordinate: uniform interval
    do i = 1, GRID_IA
       ii = i + PRC_2Drank(PRC_myrank,1) * GRID_IMAX

       GRID_FX(i) = GRID_DX * dble( ii - GRID_IHALO )
       GRID_CX(i) = GRID_FX(i) - 0.5D0 * GRID_DX
    enddo

    do j = 1, GRID_JA
       jj = j + PRC_2Drank(PRC_myrank,2) * GRID_JMAX

       GRID_FY(j) = GRID_DY * dble( jj - GRID_JHALO )
       GRID_CY(j) = GRID_FY(j) - 0.5D0 * GRID_DY
    enddo

    ! vertical coordinate: uniform interval
    do k = 0, GRID_KA
       GRID_FZ(k) = GRID_DZ * dble(k - GRID_KHALO)
       if (k /= 0) GRID_CZ(k) = GRID_FZ(k) - 0.5D0 * GRID_DZ
    enddo

    GRID_RDX = 1.D0 / GRID_DX
    GRID_RDY = 1.D0 / GRID_DY
    GRID_RDZ = 1.D0 / GRID_DZ

    if( IO_L ) write(IO_FID_LOG,*) '*** GRID INFORMATION ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f6.0,A,f6.0,A,f6.0)') '*** delta X, Y, Z [m]                  :', &
                                                             GRID_DX," : ",GRID_DY," : ",GRID_DZ
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)')       '*** No. of Grid (including HALO,1node) :', &
                                                             GRID_IA," x ",GRID_JA," x ",GRID_KA
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6)')            '*** Global Grid Index (X)              :', &
                                                             GRID_ISG," - ",GRID_IEG
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6)')            '*** Global Grid Index (Y)              :', &
                                                             GRID_JSG," - ",GRID_JEG
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)')       '*** No. of Computational Grid (1 node) :', &
                                                             GRID_IMAX, " x ", &
                                                             GRID_JMAX, " x ", &
                                                             GRID_KMAX
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f6.1,A,f6.1,A,f6.1)') '*** No. of Domain size [km]   (1 node) :', &
                                                             GRID_IMAX*GRID_DX*1.D-3, " x ", &
                                                             GRID_JMAX*GRID_DY*1.D-3, " x ", &
                                                             GRID_KMAX*GRID_DZ*1.D-3
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)')       '*** No. of Computational Grid (global) :', &
                                                             GRID_IMAX*PRC_NUM_X, " x ", &
                                                             GRID_JMAX*PRC_NUM_Y, " x ", &
                                                             GRID_KMAX
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f6.1,A,f6.1,A,f6.1)') '*** No. of Domain size [km]   (global) :', &
                                                             GRID_IMAX*PRC_NUM_X*GRID_DX*1.D-3, " x ", &
                                                             GRID_JMAX*PRC_NUM_Y*GRID_DY*1.D-3, " x ", &
                                                             GRID_KMAX*GRID_DZ*1.D-3

!    if( IO_L ) write(IO_FID_LOG,*) '*** value X'
!    if( IO_L ) write(IO_FID_LOG,'(24E10.3)')  (GRID_CX(i),i=1,GRID_IA)
!    if( IO_L ) write(IO_FID_LOG,*) '*** value Y'
!    if( IO_L ) write(IO_FID_LOG,'(24E10.3)')  (GRID_CY(j),j=1,GRID_JA)
!    if( IO_L ) write(IO_FID_LOG,*) '*** value Z'
!    if( IO_L ) write(IO_FID_LOG,'(12E15.8)')  (GRID_CZ(k),k=1,GRID_KA)

    return
  end subroutine GRID_setup

end module mod_grid
