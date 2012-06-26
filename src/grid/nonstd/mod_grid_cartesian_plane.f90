!-------------------------------------------------------------------------------
!> module GRID
!!
!! @par Description
!!          Grid module for plane cartesian coordinate
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_grid
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FILECHR, &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GRID_setup
  public :: GRID_read
  public :: GRID_write
  public :: GRID_generate

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_index.h"

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, save :: ISG, IEG ! start/end of inner domain: x, global
  integer, public, save :: JSG, JEG ! start/end of inner domain: y, global

  integer, public, save :: GRID_IJA           ! # of x*y
  integer, public, save :: GRID_IJS, GRID_IJE ! start/end of inner domain

  real(8), public, allocatable, save :: GRID_CZ(:)      ! center coordinate [m]: z, local=global
  real(8), public, allocatable, save :: GRID_CX(:)      ! center coordinate [m]: x, local
  real(8), public, allocatable, save :: GRID_CY(:)      ! center coordinate [m]: y, local
  real(8), public, allocatable, save :: GRID_CDZ(:)     ! z-length of control volume [m]
  real(8), public, allocatable, save :: GRID_CDX(:)     ! x-length of control volume [m]
  real(8), public, allocatable, save :: GRID_CDY(:)     ! y-length of control volume [m]
  real(8), public, allocatable, save :: GRID_RCDZ(:)    ! reciprocal of center-dz
  real(8), public, allocatable, save :: GRID_RCDX(:)    ! reciprocal of center-dx
  real(8), public, allocatable, save :: GRID_RCDY(:)    ! reciprocal of center-dy

  real(8), public, allocatable, save :: GRID_FZ(:)      ! face   coordinate [m]: z, local=global
  real(8), public, allocatable, save :: GRID_FX(:)      ! face   coordinate [m]: x, local
  real(8), public, allocatable, save :: GRID_FY(:)      ! face   coordinate [m]: y, local
  real(8), public, allocatable, save :: GRID_FDZ(:)     ! z-length of grid(k)-to-grid(k-1) [m]
  real(8), public, allocatable, save :: GRID_FDX(:)     ! x-length of grid(i)-to-grid(i-1) [m]
  real(8), public, allocatable, save :: GRID_FDY(:)     ! y-length of grid(j)-to-grid(j-1) [m]
  real(8), public, allocatable, save :: GRID_RFDZ(:)    ! reciprocal of face-dz
  real(8), public, allocatable, save :: GRID_RFDX(:)    ! reciprocal of face-dx
  real(8), public, allocatable, save :: GRID_RFDY(:)    ! reciprocal of face-dy

  logical, public, allocatable, save :: GRID_CZ_mask(:) ! main/buffer region mask: z
  logical, public, allocatable, save :: GRID_CX_mask(:) ! main/buffer region mask: x
  logical, public, allocatable, save :: GRID_CY_mask(:) ! main/buffer region mask: y

  real(8), public, allocatable, save :: GRID_CBFZ(:)    ! center buffer factor [0-1]: z
  real(8), public, allocatable, save :: GRID_CBFX(:)    ! center buffer factor [0-1]: x
  real(8), public, allocatable, save :: GRID_CBFY(:)    ! center buffer factor [0-1]: y
  real(8), public, allocatable, save :: GRID_FBFZ(:)    ! face   buffer factor [0-1]: z
  real(8), public, allocatable, save :: GRID_FBFX(:)    ! face   buffer factor [0-1]: x
  real(8), public, allocatable, save :: GRID_FBFY(:)    ! face   buffer factor [0-1]: y

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=IO_FILECHR), private :: GRID_IN_BASENAME  = ''
  character(len=IO_FILECHR), private :: GRID_OUT_BASENAME = ''

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup horizontal&vertical grid
  !-----------------------------------------------------------------------------
  subroutine GRID_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_2Drank,  &
       PRC_NUM_X,   &
       PRC_NUM_Y
    use mod_fileio, only: &
       FIO_setgridinfo
    implicit none

    namelist / PARAM_GRID / &
       GRID_IN_BASENAME,  &
       GRID_OUT_BASENAME

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

    ! horizontal index (global domain)
    ISG = IHALO + 1    + PRC_2Drank(PRC_myrank,1) * IMAX
    IEG = IHALO + IMAX + PRC_2Drank(PRC_myrank,1) * IMAX
    JSG = JHALO + 1    + PRC_2Drank(PRC_myrank,2) * JMAX
    JEG = JHALO + JMAX + PRC_2Drank(PRC_myrank,2) * JMAX

    if( IO_L ) write(IO_FID_LOG,*) '*** GRID INFORMATION ***'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,f6.0)')         '*** delta Z, X, Y [m]                  :', &
                                                       DXYZ
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)') '*** No. of Computational Grid (global) :', &
                                                       KMAX,           ' x ',                      &
                                                       IMAX*PRC_NUM_X, ' x ',                      &
                                                       JMAX*PRC_NUM_Y

    allocate( GRID_FZ  (0:KA) )
    allocate( GRID_FX  (0:IA) )
    allocate( GRID_FY  (0:JA) )
    allocate( GRID_FDZ (KA-1) )
    allocate( GRID_FDX (IA-1) )
    allocate( GRID_FDY (JA-1) )
    allocate( GRID_RFDZ(KA-1) )
    allocate( GRID_RFDX(IA-1) )
    allocate( GRID_RFDY(JA-1) )

    allocate( GRID_CZ  (KA) )
    allocate( GRID_CX  (IA) )
    allocate( GRID_CY  (JA) )
    allocate( GRID_CDZ (KA) )
    allocate( GRID_CDX (IA) )
    allocate( GRID_CDY (JA) )
    allocate( GRID_RCDZ(KA) )
    allocate( GRID_RCDX(IA) )
    allocate( GRID_RCDY(JA) )

    allocate( GRID_CZ_mask(KA) )
    allocate( GRID_CX_mask(IA) )
    allocate( GRID_CY_mask(JA) )

    allocate( GRID_CBFZ(KA) )
    allocate( GRID_CBFX(IA) )
    allocate( GRID_CBFY(JA) )
    allocate( GRID_FBFZ(KA) )
    allocate( GRID_FBFX(IA) )
    allocate( GRID_FBFY(JA) )

    if ( GRID_IN_BASENAME /= '' ) then
       call GRID_read
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found input grid file. Generate!'
       
       call GRID_generate
    endif

    if ( GRID_OUT_BASENAME /= '' ) then
       call GRID_write
    endif

    call FIO_setgridinfo( DXYZ, IMAX, JMAX )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)') '*** No. of Grid (including HALO,1node) :', &
                                                       KA," x ",IA," x ",JA
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6)')      '*** Global Grid Index (X)              :', &
                                                       ISG," - ",IEG
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6)')      '*** Global Grid Index (Y)              :', &
                                                       JSG," - ",JEG
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I6,A,I6,A,I6)') '*** No. of Computational Grid (1 node) :', &
                                                       KMAX,' x ',IMAX,' x ',JMAX

    if( IO_L ) write(IO_FID_LOG,*) '*** Domain size [km] (this node) :'
    if( IO_L ) write(IO_FID_LOG,'(1x,6(A,f8.3))') '  X:',                    &
                  GRID_FX(0) *1.D-3, ' -HALO- ', GRID_FX(IS-1)*1.D-3, ' | ', &
                  GRID_CX(IS)*1.D-3, ' - ',      GRID_CX(IE)  *1.D-3, ' | ', &
                  GRID_FX(IE)*1.D-3, ' -HALO- ', GRID_FX(IA)  *1.D-3
    if( IO_L ) write(IO_FID_LOG,'(1x,6(A,f8.3))') '  Y:',                    &
                  GRID_FY(0) *1.D-3, ' -HALO- ', GRID_FY(JS-1)*1.D-3, ' | ', &
                  GRID_CY(JS)*1.D-3, ' - ',      GRID_CY(JE)  *1.D-3, ' | ', &
                  GRID_FY(JE)*1.D-3, ' -HALO- ', GRID_FY(JA)  *1.D-3

    return
  end subroutine GRID_setup

  !-----------------------------------------------------------------------------
  !> Read horizontal&vertical grid
  !-----------------------------------------------------------------------------
  subroutine GRID_read
    use mod_fileio, only: &
       FIO_input_1D
    implicit none

    character(len=IO_FILECHR) :: bname

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Input grid file ***'

    write(bname,'(A,A,F15.3)') trim(GRID_IN_BASENAME)

    call FIO_input_1D( GRID_CZ(:), bname, 'CZ', 'Z1D', 1, KA,   1, .false. )
    call FIO_input_1D( GRID_CX(:), bname, 'CX', 'Z1D', 1, IA,   1, .false. )
    call FIO_input_1D( GRID_CY(:), bname, 'CY', 'Z1D', 1, JA,   1, .false. )
    call FIO_input_1D( GRID_FZ(:), bname, 'FZ', 'Z1D', 1, KA+1, 1, .false. )
    call FIO_input_1D( GRID_FX(:), bname, 'FX', 'Z1D', 1, IA+1, 1, .false. )
    call FIO_input_1D( GRID_FY(:), bname, 'FY', 'Z1D', 1, JA+1, 1, .false. )

    call FIO_input_1D( GRID_CDZ(:), bname, 'CDZ', 'Z1D', 1, KA,   1, .false. )
    call FIO_input_1D( GRID_CDX(:), bname, 'CDX', 'Z1D', 1, IA,   1, .false. )
    call FIO_input_1D( GRID_CDY(:), bname, 'CDY', 'Z1D', 1, JA,   1, .false. )
    call FIO_input_1D( GRID_FDZ(:), bname, 'FDZ', 'Z1D', 1, KA-1, 1, .false. )
    call FIO_input_1D( GRID_FDX(:), bname, 'FDX', 'Z1D', 1, IA-1, 1, .false. )
    call FIO_input_1D( GRID_FDY(:), bname, 'FDY', 'Z1D', 1, JA-1, 1, .false. )

    GRID_RCDZ(:) = 1.D0 / GRID_CDZ(:)
    GRID_RCDX(:) = 1.D0 / GRID_CDX(:)
    GRID_RCDY(:) = 1.D0 / GRID_CDY(:)
    GRID_RFDZ(:) = 1.D0 / GRID_FDZ(:)
    GRID_RFDX(:) = 1.D0 / GRID_FDX(:)
    GRID_RFDY(:) = 1.D0 / GRID_FDY(:)

    call FIO_input_1D( GRID_CBFZ(:), bname, 'CBFZ', 'Z1D', 1, KA, 1, .false. )
    call FIO_input_1D( GRID_CBFX(:), bname, 'CBFX', 'Z1D', 1, IA, 1, .false. )
    call FIO_input_1D( GRID_CBFY(:), bname, 'CBFY', 'Z1D', 1, JA, 1, .false. )
    call FIO_input_1D( GRID_FBFZ(:), bname, 'FBFZ', 'Z1D', 1, KA, 1, .false. )
    call FIO_input_1D( GRID_FBFX(:), bname, 'FBFX', 'Z1D', 1, IA, 1, .false. )
    call FIO_input_1D( GRID_FBFY(:), bname, 'FBFY', 'Z1D', 1, JA, 1, .false. )

    do k = 1, KA
       if ( GRID_CBFZ(k) == 0.D0 ) then
          GRID_CZ_mask(k) = .true.
       else
          GRID_CZ_mask(k) = .false.
       endif
    enddo

    do i = 1, IA
       if ( GRID_CBFX(i) == 0.D0 ) then
          GRID_CX_mask(i) = .true.
       else
          GRID_CX_mask(i) = .false.
       endif
    enddo

    do j = 1, JA
       if ( GRID_CBFY(j) == 0.D0 ) then
          GRID_CY_mask(j) = .true.
       else
          GRID_CY_mask(j) = .false.
       endif
    enddo

    return
  end subroutine GRID_read

  !-----------------------------------------------------------------------------
  !> Write horizontal&vertical grid
  !-----------------------------------------------------------------------------
  subroutine GRID_write
    use mod_fileio_h, only: &
       FIO_HMID, &
       FIO_REAL8
    use mod_fileio, only: &
       FIO_output_1D
    implicit none

    character(len=IO_FILECHR) :: bname
    character(len=FIO_HMID)   :: desc
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Output grid file ***'

    write(bname,'(A,A,F15.3)') trim(GRID_OUT_BASENAME)
    desc  = 'SCALE3 GRID VARS.'

    call FIO_output_1D( GRID_CZ(:), bname, desc, '',                   &
                       'CZ', 'Grid Center Position Z', '', 'm',        &
                       FIO_REAL8, 'Z1D', 1, KA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_CX(:), bname, desc, '',                   &
                       'CX', 'Grid Center Position X', '', 'm',        &
                       FIO_REAL8, 'Z1D', 1, IA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_CY(:), bname, desc, '',                   &
                       'CY', 'Grid Center Position Y', '', 'm',        &
                       FIO_REAL8, 'Z1D', 1, JA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FZ(:), bname, desc, '',                     &
                       'FZ', 'Grid Face Position Z', '', 'm',            &
                       FIO_REAL8, 'Z1D', 1, KA+1, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FX(:), bname, desc, '',                     &
                       'FX', 'Grid Face Position X', '', 'm',            &
                       FIO_REAL8, 'Z1D', 1, IA+1, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FY(:), bname, desc, '',                     &
                       'FY', 'Grid Face Position Y', '', 'm',            &
                       FIO_REAL8, 'Z1D', 1, JA+1, 1, 0.D0, 0.D0, .false. )

    call FIO_output_1D( GRID_CDZ(:), bname, desc, '',                  &
                       'CDZ', 'Grid Cell length Z', '', 'm',           &
                       FIO_REAL8, 'Z1D', 1, KA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_CDX(:), bname, desc, '',                  &
                       'CDX', 'Grid Cell length X', '', 'm',           &
                       FIO_REAL8, 'Z1D', 1, IA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_CDY(:), bname, desc, '',                  &
                       'CDY', 'Grid Cell length Y', '', 'm',           &
                       FIO_REAL8, 'Z1D', 1, JA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FDZ(:), bname, desc, '',                    &
                       'FDZ', 'Grid distance Z', '', 'm',                &
                       FIO_REAL8, 'Z1D', 1, KA-1, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FDX(:), bname, desc, '',                    &
                       'FDX', 'Grid distance X', '', 'm',                &
                       FIO_REAL8, 'Z1D', 1, IA-1, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FDY(:), bname, desc, '',                    &
                       'FDY', 'Grid distance Y', '', 'm',                &
                       FIO_REAL8, 'Z1D', 1, JA-1, 1, 0.D0, 0.D0, .false. )

    call FIO_output_1D( GRID_CBFZ(:), bname, desc, '',                 &
                       'CBFZ', 'Boundary factor Center Z', '', '0-1',  &
                       FIO_REAL8, 'Z1D', 1, KA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_CBFX(:), bname, desc, '',                 &
                       'CBFX', 'Boundary factor Center X', '', '0-1',  &
                       FIO_REAL8, 'Z1D', 1, IA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_CBFY(:), bname, desc, '',                 &
                       'CBFY', 'Boundary factor Center Y', '', '0-1',  &
                       FIO_REAL8, 'Z1D', 1, JA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FBFZ(:), bname, desc, '',                 &
                       'FBFZ', 'Boundary factor Face Z', '', '0-1',    &
                       FIO_REAL8, 'Z1D', 1, KA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FBFX(:), bname, desc, '',                 &
                       'FBFX', 'Boundary factor Face X', '', '0-1',    &
                       FIO_REAL8, 'Z1D', 1, IA, 1, 0.D0, 0.D0, .false. )
    call FIO_output_1D( GRID_FBFY(:), bname, desc, '',                 &
                       'FBFY', 'Boundary factor Face Y', '', '0-1',    &
                       FIO_REAL8, 'Z1D', 1, JA, 1, 0.D0, 0.D0, .false. )

    return
  end subroutine GRID_write

  !-----------------------------------------------------------------------------
  !> Generate horizontal&vertical grid
  !-----------------------------------------------------------------------------
  subroutine GRID_generate
    use mod_process, only: &
       PRC_MPIstop, &
       PRC_myrank,  &
       PRC_2Drank,  &
       PRC_NUM_X,   &
       PRC_NUM_Y
    implicit none

    integer :: IAG ! # of x whole cells (global, with HALO)
    integer :: JAG ! # of y whole cells (global, with HALO)

    real(8), allocatable :: GRID_FXG(:)      ! face   coordinate [m]: x, global
    real(8), allocatable :: GRID_FYG(:)      ! face   coordinate [m]: y, global
    real(8), allocatable :: GRID_CXG(:)      ! center coordinate [m]: x, global
    real(8), allocatable :: GRID_CYG(:)      ! center coordinate [m]: y, global
    real(8), allocatable :: GRID_FBFXG(:)    ! face   buffer factor [0-1]: x, global
    real(8), allocatable :: GRID_FBFYG(:)    ! face   buffer factor [0-1]: y, global
    real(8), allocatable :: GRID_CBFXG(:)    ! center buffer factor [0-1]: x, global
    real(8), allocatable :: GRID_CBFYG(:)    ! center buffer factor [0-1]: y, global
    logical, allocatable :: GRID_CXG_mask(:)
    logical, allocatable :: GRID_CYG_mask(:)

    real(8), allocatable :: buffz(:), buffx(:), buffy(:)
    real(8)              :: bufftotz, bufftotx, bufftoty

    integer :: kbuff, ibuff, jbuff
    integer :: kmain, imain, jmain

    integer :: k, i, j, ii, jj
    !---------------------------------------------------------------------------

    !##### coordinate in global domain #####

    ! array size (global domain)
    IAG = IHALO + IMAX*PRC_NUM_X + IHALO
    JAG = JHALO + JMAX*PRC_NUM_Y + JHALO

    allocate( GRID_FXG     (0:IAG) )
    allocate( GRID_FYG     (0:JAG) )
    allocate( GRID_CXG     (  IAG) )
    allocate( GRID_CYG     (  JAG) )
    allocate( GRID_CBFXG   (  IAG) )
    allocate( GRID_CBFYG   (  JAG) )
    allocate( GRID_FBFXG   (  IAG) )
    allocate( GRID_FBFYG   (  JAG) )
    allocate( GRID_CXG_mask(  IAG) )
    allocate( GRID_CYG_mask(  JAG) )

    allocate( buffx(0:IAG) )
    allocate( buffy(0:JAG) )

    ! X-direction
    ! calculate buffer grid size
    buffx(0) = DXYZ
    bufftotx = 0.D0

    do i = 1, IAG
       if( bufftotx >= BUFFER_DX ) exit

       buffx(i) = buffx(i-1) * BUFFFACT
       bufftotx = bufftotx + buffx(i)
    enddo
    ibuff = i - 1
    imain = IAG - 2*ibuff - 2*IHALO

    if ( imain < 1 ) then
       write(*,*) 'xxx Not appropriate buffer length for global domain(X). Check!', BUFFER_DX
       call PRC_MPIstop
    endif

    ! horizontal coordinate (global domaim)
    GRID_FXG(IHALO) = 0.D0
    do i = IHALO-1, 0, -1
       GRID_FXG(i) = GRID_FXG(i+1) - buffx(ibuff)
    enddo

    do i = 1, IHALO
       GRID_CXG     (i) = 0.5D0 * ( GRID_FXG(i)+GRID_FXG(i-1) )
       GRID_CXG_mask(i) = .false.
    enddo

    if ( ibuff > 0 ) then
       do i = IHALO+1, IHALO+ibuff
          GRID_FXG     (i) = GRID_FXG(i-1) + buffx(ibuff+IHALO+1-i)
          GRID_CXG     (i) = 0.5D0 * ( GRID_FXG(i)+GRID_FXG(i-1) )
          GRID_CXG_mask(i) = .false.
       enddo
    endif

    do i = IHALO+ibuff+1, IHALO+ibuff+imain
       GRID_FXG     (i) = GRID_FXG(i-1) + DXYZ
       GRID_CXG     (i) = 0.5D0 * ( GRID_FXG(i)+GRID_FXG(i-1) )
       GRID_CXG_mask(i) = .true.
    enddo

    if ( ibuff > 0 ) then
       do i = IHALO+ibuff+imain+1, IHALO+ibuff+imain+ibuff
          GRID_FXG     (i) = GRID_FXG(i-1) + buffx(i-IHALO-ibuff-imain)
          GRID_CXG     (i) = 0.5D0 * ( GRID_FXG(i)+GRID_FXG(i-1) )
          GRID_CXG_mask(i) = .false.
       enddo
    endif

    do i = IHALO+ibuff+imain+ibuff+1, IHALO+ibuff+imain+ibuff+IHALO
       GRID_FXG     (i) = GRID_FXG(i-1) + buffx(ibuff)
       GRID_CXG     (i) = 0.5D0 * ( GRID_FXG(i)+GRID_FXG(i-1) )
       GRID_CXG_mask(i) = .false.
    enddo

    ! calc buffer factor (global domaim)
    GRID_CBFXG(:) = 0.D0
    GRID_FBFXG(:) = 0.D0
    do i = 1, IHALO
       GRID_CBFXG(i) = 1.D0
       GRID_FBFXG(i) = 1.D0
    enddo

    if ( ibuff > 0 ) then
       do i = IHALO+1, IHALO+ibuff
          GRID_CBFXG(i) = (bufftotx-GRID_CXG(i)) / bufftotx
          GRID_FBFXG(i) = (bufftotx-GRID_FXG(i)) / bufftotx
       enddo

       do i = IHALO+ibuff+imain+1, IHALO+ibuff+imain+ibuff
          GRID_CBFXG(i) = (bufftotx-GRID_FXG(IAG-IHALO)+GRID_CXG(i)) / bufftotx
          GRID_FBFXG(i) = (bufftotx-GRID_FXG(IAG-IHALO)+GRID_FXG(i)) / bufftotx
       enddo
    endif

    do i = IHALO+ibuff+imain+ibuff+1, IHALO+ibuff+imain+ibuff+IHALO
       GRID_CBFXG(i) = 1.D0
       GRID_FBFXG(i) = 1.D0
    enddo

    GRID_CBFXG(:) = max( min( GRID_CBFXG(:), 1.D0 ), 0.D0 )
    GRID_FBFXG(:) = max( min( GRID_FBFXG(:), 1.D0 ), 0.D0 )

    ! Y-direction
    ! calculate buffer grid size
    buffy(0) = DXYZ
    bufftoty = 0.D0

    do j = 1, JAG
       if( bufftoty >= BUFFER_DY ) exit

       buffy(j) = buffy(j-1) * BUFFFACT
       bufftoty = bufftoty + buffy(j)
    enddo
    jbuff = j - 1
    jmain = JAG - 2*jbuff - 2*JHALO

    if ( jmain < 1 ) then
       write(*,*) 'xxx Not appropriate buffer length for global domain(Y). Check!', BUFFER_DY
       call PRC_MPIstop
    endif

    ! horizontal coordinate (global domaim)
    GRID_FYG(JHALO) = 0.D0
    do j = JHALO-1, 0, -1
       GRID_FYG(j) = GRID_FYG(j+1) - buffy(jbuff)
    enddo

    do j = 1, JHALO
       GRID_CYG     (j) = 0.5D0 * ( GRID_FYG(j)+GRID_FYG(j-1) )
       GRID_CYG_mask(j) = .false.
    enddo

    if ( jbuff > 0 ) then
       do j = JHALO+1, JHALO+jbuff
          GRID_FYG     (j) = GRID_FYG(j-1) + buffy(jbuff+JHALO+1-j)
          GRID_CYG     (j) = 0.5D0 * ( GRID_FYG(j)+GRID_FYG(j-1) )
          GRID_CYG_mask(j) = .false.
       enddo
    endif

    do j = JHALO+jbuff+1, JHALO+jbuff+jmain
       GRID_FYG     (j) = GRID_FYG(j-1) + DXYZ
       GRID_CYG     (j) = 0.5D0 * ( GRID_FYG(j)+GRID_FYG(j-1) )
       GRID_CYG_mask(j) = .true.
    enddo

    if ( jbuff > 0 ) then
       do j = JHALO+jbuff+jmain+1, JHALO+jbuff+jmain+jbuff
          GRID_FYG     (j) = GRID_FYG(j-1) + buffy(j-JHALO-jbuff-jmain)
          GRID_CYG     (j) = 0.5D0 * ( GRID_FYG(j)+GRID_FYG(j-1) )
          GRID_CYG_mask(j) = .false.
       enddo
    endif

    do j = JHALO+jbuff+jmain+jbuff+1, JHALO+jbuff+jmain+jbuff+JHALO
       GRID_FYG     (j) = GRID_FYG(j-1) + buffy(jbuff)
       GRID_CYG     (j) = 0.5D0 * ( GRID_FYG(j)+GRID_FYG(j-1) )
       GRID_CYG_mask(j) = .false.
    enddo

    ! calc buffer factor (global domaim)
    GRID_CBFYG(:) = 0.D0
    GRID_FBFYG(:) = 0.D0
    do j = 1, JHALO
       GRID_CBFYG(j) = 1.D0
       GRID_FBFYG(j) = 1.D0
    enddo

    if ( jbuff > 0 ) then
       do j = JHALO+1, JHALO+jbuff
          GRID_CBFYG(j) = (bufftoty-GRID_CYG(j)) / bufftoty
          GRID_FBFYG(j) = (bufftoty-GRID_FYG(j)) / bufftoty
       enddo

       do j = JHALO+jbuff+jmain+1, JHALO+jbuff+jmain+jbuff
          GRID_CBFYG(j) = (bufftoty-GRID_FYG(JAG-JHALO)+GRID_CYG(j)) / bufftoty
          GRID_FBFYG(j) = (bufftoty-GRID_FYG(JAG-JHALO)+GRID_FYG(j)) / bufftoty
       enddo
    endif

    do j = JHALO+jbuff+jmain+jbuff+1, JHALO+jbuff+jmain+jbuff+JHALO
       GRID_CBFYG(j) = 1.D0
       GRID_FBFYG(j) = 1.D0
    enddo
    GRID_CBFYG(:) = max( min( GRID_CBFYG(:), 1.D0 ), 0.D0 )
    GRID_FBFYG(:) = max( min( GRID_FBFYG(:), 1.D0 ), 0.D0 )

    deallocate( buffx )
    deallocate( buffy )

    !##### coordinate in local domain #####

    allocate( buffz(0:KA) )

    ! Z-direction
    ! calculate buffer grid size
    buffz(0) = DXYZ
    bufftotz = 0.D0

    do k = 1, KA
       if( bufftotz >= BUFFER_DZ ) exit

       buffz(k) = buffz(k-1) * BUFFFACT
       bufftotz = bufftotz + buffz(k)
    enddo
    kbuff = k - 1
    kmain = KA - kbuff - 2*KHALO

    if ( kmain < 1 ) then
       write(*,*) 'xxx Not appropriate buffer length for global domain(Z). Check!', BUFFER_DZ
       call PRC_MPIstop
    endif

    ! vartical coordinate (local=global domaim)
    GRID_FZ(KHALO) = 0.D0
    do k = KHALO-1, 0, -1
       GRID_FZ(k) = GRID_FZ(k+1) - DXYZ
    enddo

    do k = 1, KHALO
       GRID_CZ     (k) = 0.5D0 * ( GRID_FZ(k)+GRID_FZ(k-1) )
       GRID_CZ_mask(k) = .true.
    enddo

    do k = KHALO+1, KHALO+kmain
       GRID_FZ     (k) = GRID_FZ(k-1) + DXYZ
       GRID_CZ     (k) = 0.5D0 * ( GRID_FZ(k)+GRID_FZ(k-1) )
       GRID_CZ_mask(k) = .true.
    enddo

    if ( kbuff > 0 ) then
       do k = KHALO+kmain+1, KHALO+kmain+kbuff
          GRID_FZ     (k) = GRID_FZ(k-1) + buffz(k-KHALO-kmain)
          GRID_CZ     (k) = 0.5D0 * ( GRID_FZ(k)+GRID_FZ(k-1) )
          GRID_CZ_mask(k) = .false.
       enddo
    endif

    do k = KHALO+kmain+kbuff+1, KHALO+kmain+kbuff+KHALO
       GRID_FZ     (k) = GRID_FZ(k-1) + buffz(kbuff)
       GRID_CZ     (k) = 0.5D0 * ( GRID_FZ(k)+GRID_FZ(k-1) )
       GRID_CZ_mask(k) = .false.
    enddo

    ! calc buffer factor (global domaim)
    GRID_CBFZ(:) = 0.D0
    GRID_FBFZ(:) = 0.D0
    if ( kbuff > 0 ) then
       do k = KHALO+kmain+1, KHALO+kmain+kbuff
          GRID_CBFZ(k) = (bufftotz-GRID_FZ(KE)+GRID_CZ(k)) / bufftotz
          GRID_FBFZ(k) = (bufftotz-GRID_FZ(KE)+GRID_FZ(k)) / bufftotz
       enddo
    endif

    do k = KHALO+kmain+kbuff+1, KHALO+kmain+kbuff+KHALO
       GRID_CBFZ(k) = 1.D0
       GRID_FBFZ(k) = 1.D0
    enddo
    GRID_CBFZ(:) = max( min( GRID_CBFZ(:), 1.D0 ), 0.D0 )
    GRID_FBFZ(:) = max( min( GRID_FBFZ(:), 1.D0 ), 0.D0 )

    deallocate( buffz )

    ! vartical coordinate (local domaim)
    do k = 1, KA
       GRID_CDZ (k) = GRID_FZ(k) - GRID_FZ(k-1)
       GRID_RCDZ(k) = 1.D0 / GRID_CDZ(k)
    enddo

    do k = 1, KA-1
       GRID_FDZ (k) = GRID_CZ(k+1)-GRID_CZ(k)
       GRID_RFDZ(k) = 1.D0 / GRID_FDZ(k)
    enddo

    ! X-direction
    ! horizontal coordinate (local domaim)
    do i = 0, IA
       ii = i + PRC_2Drank(PRC_myrank,1) * IMAX

       GRID_FX(i) = GRID_FXG(ii)
    enddo

    do i = 1, IA
       ii = i + PRC_2Drank(PRC_myrank,1) * IMAX

       GRID_CX     (i) = GRID_CXG     (ii)
       GRID_CX_mask(i) = GRID_CXG_mask(ii)

       GRID_CBFX(i) = GRID_CBFXG(ii)
       GRID_FBFX(i) = GRID_FBFXG(ii)

       GRID_CDX (i) = GRID_FX(i) - GRID_FX(i-1)
       GRID_RCDX(i) = 1.D0 / GRID_CDX(i)
    enddo

    do i = 1, IA-1
       GRID_FDX (i) = GRID_CX(i+1)-GRID_CX(i)
       GRID_RFDX(i) = 1.D0 / GRID_FDX(i)
    enddo

    ! Y-direction
    ! horizontal coordinate (local domaim)
    do j = 0, JA
       jj = j + PRC_2Drank(PRC_myrank,2) * JMAX

       GRID_FY(j) = GRID_FYG(jj)
    enddo

    do j = 1, JA
       jj = j + PRC_2Drank(PRC_myrank,2) * JMAX

       GRID_CY     (j) = GRID_CYG     (jj)
       GRID_CY_mask(j) = GRID_CYG_mask(jj)

       GRID_CBFY(j) = GRID_CBFYG(jj)
       GRID_FBFY(j) = GRID_FBFYG(jj)

       GRID_CDY (j) = GRID_FY(j) - GRID_FY(j-1)
       GRID_RCDY(j) = 1.D0 / GRID_CDY(j)
    enddo

    do j = 1, JA-1
       GRID_FDY (j) = GRID_CY(j+1)-GRID_CY(j)
       GRID_RFDY(j) = 1.D0 / GRID_FDY(j)
    enddo

    ! report
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Main/buffer Grid (global) :'
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))')  '  Z: buffer = ', kbuff,' x 1, main = ',kmain
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))')  '  X: buffer = ', ibuff,' x 2, main = ',imain
    if( IO_L ) write(IO_FID_LOG,'(1x,2(A,I6))')  '  Y: buffer = ', jbuff,' x 2, main = ',jmain

    if( IO_L ) write(IO_FID_LOG,*) '*** Domain size [km] (global) :'
    if( IO_L ) write(IO_FID_LOG,'(1x,7(A,f8.3))') '  Z:',                &
                  GRID_FZ(0)       *1.D-3, ' -HALO-                   ', &
                  GRID_FZ(KS-1)    *1.D-3, ' | ',                        &
                  GRID_CZ(KS)      *1.D-3, ' - ',                        &
                  GRID_CZ(KE-kbuff)*1.D-3, ' | ',                        &
                  GRID_FZ(KE-kbuff)*1.D-3, ' -buffer- ',                 &
                  GRID_FZ(KE)      *1.D-3, ' -HALO- ',                   &
                  GRID_FZ(KA)      *1.D-3
    if( IO_L ) write(IO_FID_LOG,'(1x,8(A,f8.3))') '  X:',        &
                  GRID_FXG(0)              *1.D-3, ' -HALO- ',   &
                  GRID_FXG(IHALO)          *1.D-3, ' -buffer- ', &
                  GRID_FXG(IHALO+ibuff)    *1.D-3, ' | ',        &
                  GRID_CXG(IHALO+ibuff+1)  *1.D-3, ' - ',        &
                  GRID_CXG(IAG-IHALO-ibuff)*1.D-3, ' | ',        &
                  GRID_FXG(IAG-IHALO-ibuff)*1.D-3, ' -buffer- ', &
                  GRID_FXG(IAG-IHALO)      *1.D-3, ' -HALO- ',   &
                  GRID_FXG(IAG)            *1.D-3
    if( IO_L ) write(IO_FID_LOG,'(1x,8(A,f8.3))') '  Y:',        &
                  GRID_FYG(0)              *1.D-3, ' -HALO- ',   &
                  GRID_FYG(JHALO)          *1.D-3, ' -buffer- ', &
                  GRID_FYG(JHALO+jbuff)    *1.D-3, ' | ',        &
                  GRID_CYG(JHALO+jbuff+1)  *1.D-3, ' - ',        &
                  GRID_CYG(JAG-JHALO-jbuff)*1.D-3, ' | ',        &
                  GRID_FYG(JAG-JHALO-jbuff)*1.D-3, ' -buffer- ', &
                  GRID_FYG(JAG-JHALO)      *1.D-3, ' -HALO- ',   &
                  GRID_FYG(JAG)            *1.D-3

! for debug
!    if( IO_L ) write(IO_FID_LOG,*)
!    if( IO_L ) write(IO_FID_LOG,*) ' ', 0, GRID_FZ(0)
!    do k = 1, KA-1
!       if( IO_L ) write(IO_FID_LOG,*) k, GRID_CZ(k), GRID_CBFZ(k), GRID_CDZ(k), GRID_CZ_mask(k)
!       if( IO_L ) write(IO_FID_LOG,*) ' ', k, GRID_FZ(k), GRID_FBFZ(k), GRID_FDZ(k)
!    enddo
!    k = KA
!    if( IO_L ) write(IO_FID_LOG,*) k, GRID_CZ(k), GRID_CBFZ(k), GRID_CDZ(k), GRID_CZ_mask(k)
!    if( IO_L ) write(IO_FID_LOG,*) ' ', k, GRID_FZ(k), GRID_FBFZ(k)

!    if( IO_L ) write(IO_FID_LOG,*)
!    if( IO_L ) write(IO_FID_LOG,*) ' ', 0, GRID_FX(0)
!    do i = 1, IA-1
!       if( IO_L ) write(IO_FID_LOG,*) i, GRID_CX(i), GRID_CBFX(i), GRID_CDX(i), GRID_CX_mask(i)
!       if( IO_L ) write(IO_FID_LOG,*) ' ', i, GRID_FX(i), GRID_FBFX(i), GRID_FDX(i)
!    enddo
!    i = IA
!    if( IO_L ) write(IO_FID_LOG,*) i, GRID_CX(i), GRID_CBFX(i), GRID_CDX(i), GRID_CX_mask(i)
!    if( IO_L ) write(IO_FID_LOG,*) ' ', i, GRID_FX(i), GRID_FBFX(i)

!    if( IO_L ) write(IO_FID_LOG,*)
!    if( IO_L ) write(IO_FID_LOG,*) ' ', 0, GRID_FY(0)
!    do j = 1, JA-1
!       if( IO_L ) write(IO_FID_LOG,*) j, GRID_CY(j), GRID_CBFY(j), GRID_CDY(j), GRID_CY_mask(j)
!       if( IO_L ) write(IO_FID_LOG,*) ' ', j, GRID_FY(j), GRID_FBFY(j), GRID_FDY(j)
!    enddo
!    j = JA
!    if( IO_L ) write(IO_FID_LOG,*) j, GRID_CY(j), GRID_CBFY(j), GRID_CDY(j), GRID_CY_mask(j)
!    if( IO_L ) write(IO_FID_LOG,*) ' ', j, GRID_FY(j), GRID_FBFY(j)

    deallocate( GRID_FXG      )
    deallocate( GRID_FYG      )
    deallocate( GRID_CXG      )
    deallocate( GRID_CYG      )
    deallocate( GRID_FBFXG    )
    deallocate( GRID_FBFYG    )
    deallocate( GRID_CBFXG    )
    deallocate( GRID_CBFYG    )
    deallocate( GRID_CXG_mask )
    deallocate( GRID_CYG_mask )

    return
  end subroutine GRID_generate

end module mod_grid
