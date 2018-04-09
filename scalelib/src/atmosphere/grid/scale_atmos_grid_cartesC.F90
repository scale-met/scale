!-------------------------------------------------------------------------------
!> module atmosphere / grid / cartesC
!!
!! @par Description
!!          Atmospheric grid module for the cartesianC coordinate
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_grid_cartesC
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_GRID_CARTESC_setup
  public :: ATMOS_GRID_CARTESC_allocate
  public :: ATMOS_GRID_CARTESC_generate

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=7), public, parameter :: ATMOS_GRID_CARTESC_NAME = 'cartesC'

  real(RP), public              :: DZ, DX, DY

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CZ   (:) !< center coordinate [m]: z, local
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FZ   (:) !< face   coordinate [m]: z, local
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CDZ  (:) !< z-length of control volume [m]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FDZ  (:) !< z-length of grid(i+1) to grid(i) [m]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_RCDZ (:) !< reciprocal of center-dz
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_RFDZ (:) !< reciprocal of face-dz
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CBFZ (:) !< center buffer factor (0-1): z
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FBFZ (:) !< face   buffer factor (0-1): z

  ! land
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_LCZ  (:) !< center coordinate [m]: z, local land
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_LFZ (:) !< face   coordinate [m]: z, local
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_LCDZ(:) !< z-length of control volume [m]

  ! horizontal
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CX   (:) !< center coordinate [m]: x, local
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CY   (:) !< center coordinate [m]: y, local
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FX   (:) !< face   coordinate [m]: x, local
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FY   (:) !< face   coordinate [m]: y, local

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CDX  (:) !< x-length of control volume [m]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CDY  (:) !< y-length of control volume [m]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FDX  (:) !< x-length of grid(i+1) to grid(i) [m]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FDY  (:) !< y-length of grid(j+1) to grid(j) [m]

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_RCDX (:) !< reciprocal of center-dx
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_RCDY (:) !< reciprocal of center-dy
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_RFDX (:) !< reciprocal of face-dx
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_RFDY (:) !< reciprocal of face-dy

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CBFX (:) !< center buffer factor (0-1): x
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CBFY (:) !< center buffer factor (0-1): y
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FBFX (:) !< face   buffer factor (0-1): x
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FBFY (:) !< face   buffer factor (0-1): y

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CXG  (:) !< center coordinate [m]: x, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CYG  (:) !< center coordinate [m]: y, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FXG  (:) !< face   coordinate [m]: x, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FYG  (:) !< face   coordinate [m]: y, global

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CDXG (:) !< center coordinate [m]: x, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CDYG (:) !< center coordinate [m]: y, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FDXG (:) !< center coordinate [m]: x, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FDYG (:) !< center coordinate [m]: y, global

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CBFXG(:) !< center buffer factor (0-1): x, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_CBFYG(:) !< center buffer factor (0-1): y, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FBFXG(:) !< face   buffer factor (0-1): x, global
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_FBFYG(:) !< face   buffer factor (0-1): y, global

  real(RP), public              :: ATMOS_GRID_CARTESC_DOMAIN_CENTER_X !< center position of global domain [m]: x
  real(RP), public              :: ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y !< center position of global domain [m]: y

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
  !> Setup
  subroutine ATMOS_GRID_CARTESC_setup( &
       basename, &
       aggregate )
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_AGGREGATE
    implicit none
    character(len=*), intent(in), optional :: basename
    logical,          intent(in), optional :: aggregate

    character(len=H_LONG) :: ATMOS_GRID_CARTESC_IN_BASENAME  = ''
    logical               :: ATMOS_GRID_CARTESC_IN_AGGREGATE

    real(RP)              :: OFFSET_X     = 0.0_RP
    real(RP)              :: OFFSET_Y     = 0.0_RP

    real(RP)              :: BUFFER_DZ  =  0.0_RP !< thickness of buffer region [m]: z
    real(RP)              :: BUFFER_DX  =  0.0_RP !< thickness of buffer region [m]: x
    real(RP)              :: BUFFER_DY  =  0.0_RP !< thickness of buffer region [m]: y
    real(RP)              :: BUFFFACT   =  1.0_RP !< default strech factor for dx/dy/dz of buffer region
    real(RP)              :: BUFFFACT_Z = -1.0_RP !< strech factor for dz of buffer region
    real(RP)              :: BUFFFACT_X = -1.0_RP !< strech factor for dx of buffer region
    real(RP)              :: BUFFFACT_Y = -1.0_RP !< strech factor for dy of buffer region

    integer               :: BUFFER_NZ = -1    !< thickness of buffer region by number of grids: z
    integer               :: BUFFER_NX = -1    !< thickness of buffer region by number of grids: x
    integer               :: BUFFER_NY = -1    !< thickness of buffer region by number of grids: y

    integer,    parameter :: FZ_MAX = 300 !< limit of index size for user defined z
    real(RP)              :: FZ(FZ_MAX)   !< user defined center coordinate [m]: z, local=global

    namelist / PARAM_ATMOS_GRID_CARTESC / &
       ATMOS_GRID_CARTESC_IN_BASENAME,  &
       ATMOS_GRID_CARTESC_IN_AGGREGATE, &
       DZ,         &
       DX,         &
       DY,         &
       BUFFER_DZ,  &
       BUFFER_DX,  &
       BUFFER_DY,  &
       BUFFER_NZ,  &
       BUFFER_NX,  &
       BUFFER_NY,  &
       BUFFFACT,   &
       BUFFFACT_Z, &
       BUFFFACT_X, &
       BUFFFACT_Y, &
       FZ,         &
       OFFSET_X, &
       OFFSET_Y

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[CartesC] / Categ[ATMOSPHER GRID] / Origin[SCALElib]'

    FZ(:) = -1.0_RP

    if ( present(basename)  ) ATMOS_GRID_CARTESC_IN_BASENAME = basename
    if ( present(aggregate) ) then
       ATMOS_GRID_CARTESC_IN_AGGREGATE = aggregate
    else
       ATMOS_GRID_CARTESC_IN_AGGREGATE = FILE_AGGREGATE
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_GRID_CARTESC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_GRID_CARTESC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_GRID_CARTESC_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_GRID_CARTESC. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_GRID_CARTESC)


    call ATMOS_GRID_CARTESC_allocate

    if ( ATMOS_GRID_CARTESC_IN_BASENAME /= '' ) then

       call ATMOS_GRID_CARTESC_read( ATMOS_GRID_CARTESC_IN_BASENAME, ATMOS_GRID_CARTESC_IN_AGGREGATE )

    else

       call ATMOS_GRID_CARTESC_generate( &
            DZ, DX, DY, FZ(:), FZ_MAX,         &
            OFFSET_X, OFFSET_Y,                &
            BUFFER_DZ, BUFFER_DX, BUFFER_DY,   &
            BUFFER_NZ, BUFFER_NX, BUFFER_NY,   &
            BUFFFACT,                          &
            BUFFFACT_Z, BUFFFACT_X, BUFFFACT_Y )

    end if

    call ATMOS_GRID_CARTESC_output_info

    return
  end subroutine ATMOS_GRID_CARTESC_setup

  !-----------------------------------------------------------------------------
  ! private
  !-----------------------------------------------------------------------------

  subroutine ATMOS_GRID_CARTESC_allocate
    implicit none
    !---------------------------------------------------------------------------

    ! local domain
    allocate( ATMOS_GRID_CARTESC_CZ  (  KA) )
    allocate( ATMOS_GRID_CARTESC_CX  (  IA) )
    allocate( ATMOS_GRID_CARTESC_CY  (  JA) )
    allocate( ATMOS_GRID_CARTESC_FZ  (0:KA) )
    allocate( ATMOS_GRID_CARTESC_FX  (0:IA) )
    allocate( ATMOS_GRID_CARTESC_FY  (0:JA) )

    allocate( ATMOS_GRID_CARTESC_CDZ (KA)   )
    allocate( ATMOS_GRID_CARTESC_CDX (IA)   )
    allocate( ATMOS_GRID_CARTESC_CDY (JA)   )
    allocate( ATMOS_GRID_CARTESC_FDZ (KA-1) )
    allocate( ATMOS_GRID_CARTESC_FDX (IA-1) )
    allocate( ATMOS_GRID_CARTESC_FDY (JA-1) )

    allocate( ATMOS_GRID_CARTESC_RCDZ(KA)   )
    allocate( ATMOS_GRID_CARTESC_RCDX(IA)   )
    allocate( ATMOS_GRID_CARTESC_RCDY(JA)   )
    allocate( ATMOS_GRID_CARTESC_RFDZ(KA-1) )
    allocate( ATMOS_GRID_CARTESC_RFDX(IA-1) )
    allocate( ATMOS_GRID_CARTESC_RFDY(JA-1) )

    allocate( ATMOS_GRID_CARTESC_CBFZ(  KA) )
    allocate( ATMOS_GRID_CARTESC_CBFX(  IA) )
    allocate( ATMOS_GRID_CARTESC_CBFY(  JA) )
    allocate( ATMOS_GRID_CARTESC_FBFZ(0:KA) )
    allocate( ATMOS_GRID_CARTESC_FBFX(0:IA) )
    allocate( ATMOS_GRID_CARTESC_FBFY(0:JA) )

    ! global domain
    allocate( ATMOS_GRID_CARTESC_CXG  (  IAG) )
    allocate( ATMOS_GRID_CARTESC_CYG  (  JAG) )
    allocate( ATMOS_GRID_CARTESC_FXG  (0:IAG) )
    allocate( ATMOS_GRID_CARTESC_FYG  (0:JAG) )

    allocate( ATMOS_GRID_CARTESC_CDXG (IAG)   )
    allocate( ATMOS_GRID_CARTESC_CDYG (JAG)   )
    allocate( ATMOS_GRID_CARTESC_FDXG (IAG-1) )
    allocate( ATMOS_GRID_CARTESC_FDYG (JAG-1) )

    allocate( ATMOS_GRID_CARTESC_CBFXG(  IAG) )
    allocate( ATMOS_GRID_CARTESC_CBFYG(  JAG) )
    allocate( ATMOS_GRID_CARTESC_FBFXG(0:IAG) )
    allocate( ATMOS_GRID_CARTESC_FBFYG(0:JAG) )

    return
  end subroutine ATMOS_GRID_CARTESC_allocate

  !-----------------------------------------------------------------------------
  !> Read horizontal&vertical grid
  subroutine ATMOS_GRID_CARTESC_read( &
       basename, aggregate )
    use scale_file, only: &
       FILE_open, &
       FILE_read
    use scale_prc, only: &
       PRC_myrank
    implicit none

    character(len=*), intent(in) :: basename
    logical, intent(in), optional :: aggregate

    integer :: fid

    real(RP) :: FDXG(0:IAG), FDYG(0:JAG)
    real(RP) :: FDX(0:IA), FDY(0:JA)
    !---------------------------------------------------------------------------


    call FILE_open( basename, fid, rankid=PRC_myrank, aggregate=aggregate )

    call FILE_read( fid, 'CZ', ATMOS_GRID_CARTESC_CZ(:) )
    call FILE_read( fid, 'CX', ATMOS_GRID_CARTESC_CX(:) )
    call FILE_read( fid, 'CY', ATMOS_GRID_CARTESC_CY(:) )

    call FILE_read( fid, 'FZ', ATMOS_GRID_CARTESC_FZ(:) )
    call FILE_read( fid, 'FX', ATMOS_GRID_CARTESC_FX(:) )
    call FILE_read( fid, 'FY', ATMOS_GRID_CARTESC_FY(:) )

    call FILE_read( fid, 'CDZ', ATMOS_GRID_CARTESC_CDZ(:) )
    call FILE_read( fid, 'CDX', ATMOS_GRID_CARTESC_CDX(:) )
    call FILE_read( fid, 'CDY', ATMOS_GRID_CARTESC_CDY(:) )

    call FILE_read( fid, 'FDZ', ATMOS_GRID_CARTESC_FDZ(:) )
    call FILE_read( fid, 'FDX',                    FDX(:) )
    call FILE_read( fid, 'FDY',                    FDY(:) )
    ATMOS_GRID_CARTESC_FDX(:) = FDX(1:IA-1)
    ATMOS_GRID_CARTESC_FDY(:) = FDY(1:JA-1)

    ATMOS_GRID_CARTESC_RCDZ(:) = 1.0_RP / ATMOS_GRID_CARTESC_CDZ(:)
    ATMOS_GRID_CARTESC_RCDX(:) = 1.0_RP / ATMOS_GRID_CARTESC_CDX(:)
    ATMOS_GRID_CARTESC_RCDY(:) = 1.0_RP / ATMOS_GRID_CARTESC_CDY(:)
    ATMOS_GRID_CARTESC_RFDZ(:) = 1.0_RP / ATMOS_GRID_CARTESC_FDZ(:)
    ATMOS_GRID_CARTESC_RFDX(:) = 1.0_RP / ATMOS_GRID_CARTESC_FDX(:)
    ATMOS_GRID_CARTESC_RFDY(:) = 1.0_RP / ATMOS_GRID_CARTESC_FDY(:)

    call FILE_read( fid, 'CBFZ', ATMOS_GRID_CARTESC_CBFZ(:) )
    call FILE_read( fid, 'CBFX', ATMOS_GRID_CARTESC_CBFX(:) )
    call FILE_read( fid, 'CBFY', ATMOS_GRID_CARTESC_CBFY(:) )
    call FILE_read( fid, 'FBFZ', ATMOS_GRID_CARTESC_FBFZ(:) )
    call FILE_read( fid, 'FBFX', ATMOS_GRID_CARTESC_FBFX(:) )
    call FILE_read( fid, 'FBFY', ATMOS_GRID_CARTESC_FBFY(:) )

    call FILE_read( fid, 'CXG', ATMOS_GRID_CARTESC_CXG(:) )
    call FILE_read( fid, 'CYG', ATMOS_GRID_CARTESC_CYG(:) )
    call FILE_read( fid, 'FXG', ATMOS_GRID_CARTESC_FXG(:) )
    call FILE_read( fid, 'FYG', ATMOS_GRID_CARTESC_FYG(:) )

    call FILE_read( fid, 'CDXG', ATMOS_GRID_CARTESC_CDXG(:) )
    call FILE_read( fid, 'CDYG', ATMOS_GRID_CARTESC_CDYG(:) )
    call FILE_read( fid, 'FDXG',                    FDXG(:) )
    call FILE_read( fid, 'FDYG',                    FDYG(:) )
    ATMOS_GRID_CARTESC_FDXG(:) = FDXG(1:IA-1)
    ATMOS_GRID_CARTESC_FDYG(:) = FDYG(1:JA-1)


    ATMOS_GRID_CARTESC_DOMAIN_CENTER_X = 0.5_RP * ( ATMOS_GRID_CARTESC_FXG(IHALO) + ATMOS_GRID_CARTESC_FXG(IAG-IHALO) )
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y = 0.5_RP * ( ATMOS_GRID_CARTESC_FYG(JHALO) + ATMOS_GRID_CARTESC_FYG(JAG-JHALO) )

    return
  end subroutine ATMOS_GRID_CARTESC_read

  !-----------------------------------------------------------------------------
  !> Generate horizontal&vertical grid
  subroutine ATMOS_GRID_CARTESC_generate( &
       DZ, DX, DY, FZ, FZ_MAX,            &
       OFFSET_X, OFFSET_Y,                &
       BUFFER_DZ, BUFFER_DX, BUFFER_DY,   &
       BUFFER_NZ, BUFFER_NX, BUFFER_NY,   &
       BUFFFACT,                          &
       BUFFFACT_Z, BUFFFACT_X, BUFFFACT_Y )
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    use scale_prc_cartesC, only: &
       PRC_2Drank,  &
       PRC_NUM_X,   &
       PRC_NUM_Y
    implicit none
    real(RP), intent(in) :: DZ, DX, DY
    real(RP), intent(in), optional :: FZ(:)
    integer,  intent(in), optional :: FZ_MAX
    real(RP), intent(in), optional :: OFFSET_X, OFFSET_Y
    real(RP), intent(in), optional :: BUFFER_DZ, BUFFER_DX, BUFFER_DY
    integer,  intent(in), optional :: BUFFER_NZ, BUFFER_NX, BUFFER_NY
    real(RP), intent(in), optional :: BUFFFACT
    real(RP), intent(in), optional :: BUFFFACT_Z, BUFFFACT_X, BUFFFACT_Y

    real(RP), allocatable :: buffz(:), buffx(:), buffy(:)
    real(RP)              :: bufftotz, bufftotx, bufftoty
    real(RP)              :: fact

    integer :: kbuff, ibuff, jbuff
    integer :: kmain, imain, jmain

    real(RP) :: dz_tmp

    logical :: use_user_input

    integer :: k, i, j, ii, jj
    !---------------------------------------------------------------------------

    !##### coordinate in global domain #####

    allocate( buffx(0:IAG) )
    allocate( buffy(0:JAG) )


    ! X-direction
    ! calculate buffer grid size

    fact = -1.0_RP
    if ( present(BUFFFACT_X) ) fact = BUFFFACT_X
    if ( fact < 0.0_RP .and. present(BUFFFACT) ) fact = BUFFFACT
    if ( fact < 0.0_RP ) fact = 1.0_RP

    buffx(0) = DX
    bufftotx = 0.0_RP
    ibuff = -1
    if ( present(BUFFER_NX) ) ibuff = BUFFER_NX
    if ( ibuff > 0 ) then
       if ( 2*ibuff > IMAXG ) then
          LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer grid size (', ibuff, &
                     'x2) must be smaller than global domain size (X). Use smaller BUFFER_NX!'
          call PRC_abort
       endif

       do i = 1, ibuff
          buffx(i) = buffx(i-1) * fact
          bufftotx = bufftotx + buffx(i)
       enddo
       imain = IMAXG - 2*ibuff
    else if ( present(BUFFER_DZ) ) then
       do i = 1, IAG
          if( bufftotx >= BUFFER_DX ) exit
          buffx(i) = buffx(i-1) * fact
          bufftotx = bufftotx + buffx(i)
       enddo
       ibuff = i - 1
       imain = IMAXG - 2*ibuff

       if ( imain < 0 ) then
          LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer length (', bufftotx, &
                     'x2[m]) must be smaller than global domain size (X). Use smaller BUFFER_DX!'
          call PRC_abort
       endif
    else
       ibuff = 0
       imain = IMAXG
    endif

    ! horizontal coordinate (global domain)
    if ( present(OFFSET_X) ) then
       ATMOS_GRID_CARTESC_FXG(IHALO) = OFFSET_X
    else
       ATMOS_GRID_CARTESC_FXG(IHALO) = 0.0_RP
    end if
    do i = IHALO-1, 0, -1
       ATMOS_GRID_CARTESC_FXG(i) = ATMOS_GRID_CARTESC_FXG(i+1) - buffx(ibuff)
    enddo

    do i = 1, IHALO
       ATMOS_GRID_CARTESC_CXG(i) = 0.5_RP * ( ATMOS_GRID_CARTESC_FXG(i)+ATMOS_GRID_CARTESC_FXG(i-1) )
    enddo

    if ( ibuff > 0 ) then
       do i = IHALO+1, IHALO+ibuff
          ATMOS_GRID_CARTESC_FXG(i) = ATMOS_GRID_CARTESC_FXG(i-1) + buffx(ibuff+IHALO+1-i)
          ATMOS_GRID_CARTESC_CXG(i) = 0.5_RP * ( ATMOS_GRID_CARTESC_FXG(i)+ATMOS_GRID_CARTESC_FXG(i-1) )
       enddo
    endif

    do i = IHALO+ibuff+1, IHALO+ibuff+imain
       ATMOS_GRID_CARTESC_FXG(i) = ATMOS_GRID_CARTESC_FXG(i-1) + DX
       ATMOS_GRID_CARTESC_CXG(i) = 0.5_RP * ( ATMOS_GRID_CARTESC_FXG(i)+ATMOS_GRID_CARTESC_FXG(i-1) )
    enddo

    if ( ibuff > 0 ) then
       do i = IHALO+ibuff+imain+1, IHALO+ibuff+imain+ibuff
          ATMOS_GRID_CARTESC_FXG(i) = ATMOS_GRID_CARTESC_FXG(i-1) + buffx(i-IHALO-ibuff-imain)
          ATMOS_GRID_CARTESC_CXG(i) = 0.5_RP * ( ATMOS_GRID_CARTESC_FXG(i)+ATMOS_GRID_CARTESC_FXG(i-1) )
       enddo
    endif

    do i = IHALO+ibuff+imain+ibuff+1, IHALO+ibuff+imain+ibuff+IHALO
       ATMOS_GRID_CARTESC_FXG(i) = ATMOS_GRID_CARTESC_FXG(i-1) + buffx(ibuff)
       ATMOS_GRID_CARTESC_CXG(i) = 0.5_RP * ( ATMOS_GRID_CARTESC_FXG(i)+ATMOS_GRID_CARTESC_FXG(i-1) )
    enddo

    do i = 1, IAG
       ATMOS_GRID_CARTESC_CDXG(i) = ATMOS_GRID_CARTESC_FXG(i) - ATMOS_GRID_CARTESC_FXG(i-1)
    end do
    do i = 1, IAG-1
       ATMOS_GRID_CARTESC_FDXG(i) = ATMOS_GRID_CARTESC_CXG(i+1)-ATMOS_GRID_CARTESC_CXG(i)
    end do

    ! calc buffer factor (global domain)
    ATMOS_GRID_CARTESC_CBFXG(:) = 0.0_RP
    ATMOS_GRID_CARTESC_FBFXG(:) = 0.0_RP
    do i = 1, IHALO
       ATMOS_GRID_CARTESC_CBFXG(i) = 1.0_RP
    enddo
    do i = 0, IHALO
       ATMOS_GRID_CARTESC_FBFXG(i) = 1.0_RP
    enddo

    if ( ibuff > 0 ) then
       do i = IHALO+1, IHALO+ibuff
          ATMOS_GRID_CARTESC_CBFXG(i) = (bufftotx+ATMOS_GRID_CARTESC_FXG(IHALO)-ATMOS_GRID_CARTESC_CXG(i)) / bufftotx
          ATMOS_GRID_CARTESC_FBFXG(i) = (bufftotx+ATMOS_GRID_CARTESC_FXG(IHALO)-ATMOS_GRID_CARTESC_FXG(i)) / bufftotx
       enddo

       do i = IHALO+ibuff+imain+1, IHALO+ibuff+imain+ibuff
          ATMOS_GRID_CARTESC_CBFXG(i) = (bufftotx-ATMOS_GRID_CARTESC_FXG(IAG-IHALO)+ATMOS_GRID_CARTESC_CXG(i)) / bufftotx
          ATMOS_GRID_CARTESC_FBFXG(i) = (bufftotx-ATMOS_GRID_CARTESC_FXG(IAG-IHALO)+ATMOS_GRID_CARTESC_FXG(i)) / bufftotx
       enddo
    endif

    do i = IHALO+ibuff+imain+ibuff+1, IHALO+ibuff+imain+ibuff+IHALO
       ATMOS_GRID_CARTESC_CBFXG(i) = 1.0_RP
       ATMOS_GRID_CARTESC_FBFXG(i) = 1.0_RP
    enddo

    ATMOS_GRID_CARTESC_CBFXG(:) = max( min( ATMOS_GRID_CARTESC_CBFXG(:), 1.0_RP ), 0.0_RP )
    ATMOS_GRID_CARTESC_FBFXG(:) = max( min( ATMOS_GRID_CARTESC_FBFXG(:), 1.0_RP ), 0.0_RP )

    ! Y-direction
    ! calculate buffer grid size

    fact = -1.0_RP
    if ( present(BUFFFACT_Y) ) fact = BUFFFACT_Y
    if ( fact < 0.0_RP .and. present(BUFFFACT) )fact = BUFFFACT
    if ( fact < 0.0_RP ) fact = 1.0_RP

    buffy(0) = DY
    bufftoty = 0.0_RP
    jbuff = -1
    if ( present(BUFFER_NY) ) jbuff = BUFFER_NY
    if ( jbuff > 0 ) then
       if ( 2*jbuff > JMAXG ) then
          LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer grid size (', jbuff, &
                     'x2) must be smaller than global domain size (Y). Use smaller BUFFER_NY!'
          call PRC_abort
       endif

       do j = 1, jbuff
          buffy(j) = buffy(j-1) * fact
          bufftoty = bufftoty + buffy(j)
       enddo
       jmain = JMAXG - 2*jbuff
    else if ( present(BUFFER_DY) ) then
       do j = 1, JAG
          if( bufftoty >= BUFFER_DY ) exit
          buffy(j) = buffy(j-1) * fact
          bufftoty = bufftoty + buffy(j)
       enddo
       jbuff = j - 1
       jmain = JMAXG - 2*jbuff

       if ( jmain < 0 ) then
          LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer length (', bufftoty, &
                     'x2[m]) must be smaller than global domain size (Y). Use smaller BUFFER_DY!'
          call PRC_abort
       endif
    else
       jbuff = 0
       jmain = JMAXG
    endif

    ! horizontal coordinate (global domain)
    if ( present(OFFSET_Y) ) then
       ATMOS_GRID_CARTESC_FYG(JHALO) = OFFSET_Y
    else
       ATMOS_GRID_CARTESC_FYG(JHALO) = 0.0_RP
    end if
    do j = JHALO-1, 0, -1
       ATMOS_GRID_CARTESC_FYG(j) = ATMOS_GRID_CARTESC_FYG(j+1) - buffy(jbuff)
    enddo

    do j = 1, JHALO
       ATMOS_GRID_CARTESC_CYG(j) = 0.5_RP * ( ATMOS_GRID_CARTESC_FYG(j)+ATMOS_GRID_CARTESC_FYG(j-1) )
    enddo

    if ( jbuff > 0 ) then
       do j = JHALO+1, JHALO+jbuff
          ATMOS_GRID_CARTESC_FYG(j) = ATMOS_GRID_CARTESC_FYG(j-1) + buffy(jbuff+JHALO+1-j)
          ATMOS_GRID_CARTESC_CYG(j) = 0.5_RP * ( ATMOS_GRID_CARTESC_FYG(j)+ATMOS_GRID_CARTESC_FYG(j-1) )
       enddo
    endif

    do j = JHALO+jbuff+1, JHALO+jbuff+jmain
       ATMOS_GRID_CARTESC_FYG(j) = ATMOS_GRID_CARTESC_FYG(j-1) + DY
       ATMOS_GRID_CARTESC_CYG(j) = 0.5_RP * ( ATMOS_GRID_CARTESC_FYG(j)+ATMOS_GRID_CARTESC_FYG(j-1) )
    enddo

    if ( jbuff > 0 ) then
       do j = JHALO+jbuff+jmain+1, JHALO+jbuff+jmain+jbuff
          ATMOS_GRID_CARTESC_FYG(j) = ATMOS_GRID_CARTESC_FYG(j-1) + buffy(j-JHALO-jbuff-jmain)
          ATMOS_GRID_CARTESC_CYG(j) = 0.5_RP * ( ATMOS_GRID_CARTESC_FYG(j)+ATMOS_GRID_CARTESC_FYG(j-1) )
       enddo
    endif

    do j = JHALO+jbuff+jmain+jbuff+1, JHALO+jbuff+jmain+jbuff+JHALO
       ATMOS_GRID_CARTESC_FYG(j) = ATMOS_GRID_CARTESC_FYG(j-1) + buffy(jbuff)
       ATMOS_GRID_CARTESC_CYG(j) = 0.5_RP * ( ATMOS_GRID_CARTESC_FYG(j)+ATMOS_GRID_CARTESC_FYG(j-1) )
    enddo

    do j = 1, JAG
       ATMOS_GRID_CARTESC_CDYG(j) = ATMOS_GRID_CARTESC_FYG(j) - ATMOS_GRID_CARTESC_FYG(j-1)
    end do
    do j = 1, JAG-1
       ATMOS_GRID_CARTESC_FDYG(j) = ATMOS_GRID_CARTESC_CYG(j+1)-ATMOS_GRID_CARTESC_CYG(j)
    end do

    ! calc buffer factor (global domain)
    ATMOS_GRID_CARTESC_CBFYG(:) = 0.0_RP
    ATMOS_GRID_CARTESC_FBFYG(:) = 0.0_RP
    do j = 1, JHALO
       ATMOS_GRID_CARTESC_CBFYG(j) = 1.0_RP
    enddo
    do j = 0, JHALO
       ATMOS_GRID_CARTESC_FBFYG(j) = 1.0_RP
    enddo

    if ( jbuff > 0 ) then
       do j = JHALO+1, JHALO+jbuff
          ATMOS_GRID_CARTESC_CBFYG(j) = (bufftoty+ATMOS_GRID_CARTESC_FYG(JHALO)-ATMOS_GRID_CARTESC_CYG(j)) / bufftoty
          ATMOS_GRID_CARTESC_FBFYG(j) = (bufftoty+ATMOS_GRID_CARTESC_FYG(JHALO)-ATMOS_GRID_CARTESC_FYG(j)) / bufftoty
       enddo

       do j = JHALO+jbuff+jmain+1, JHALO+jbuff+jmain+jbuff
          ATMOS_GRID_CARTESC_CBFYG(j) = (bufftoty-ATMOS_GRID_CARTESC_FYG(JAG-JHALO)+ATMOS_GRID_CARTESC_CYG(j)) / bufftoty
          ATMOS_GRID_CARTESC_FBFYG(j) = (bufftoty-ATMOS_GRID_CARTESC_FYG(JAG-JHALO)+ATMOS_GRID_CARTESC_FYG(j)) / bufftoty
       enddo
    endif

    do j = JHALO+jbuff+jmain+jbuff+1, JHALO+jbuff+jmain+jbuff+JHALO
       ATMOS_GRID_CARTESC_CBFYG(j) = 1.0_RP
       ATMOS_GRID_CARTESC_FBFYG(j) = 1.0_RP
    enddo
    ATMOS_GRID_CARTESC_CBFYG(:) = max( min( ATMOS_GRID_CARTESC_CBFYG(:), 1.0_RP ), 0.0_RP )
    ATMOS_GRID_CARTESC_FBFYG(:) = max( min( ATMOS_GRID_CARTESC_FBFYG(:), 1.0_RP ), 0.0_RP )

    deallocate( buffx )
    deallocate( buffy )

    !##### coordinate in local domain #####

    allocate( buffz(0:KA) )

    use_user_input = .false.
    if ( present(FZ) ) then
       if ( maxval(FZ(:)) > 0.0_RP ) then ! try to use input from namelist
          LOG_INFO("ATMOS_GRID_CARTESC_generate",*) 'Z coordinate is given from NAMELIST.'

          if ( KMAX < 2 ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'KMAX must be larger than 1. Check!', KMAX
             call PRC_abort
          endif

          if ( KMAX > FZ_MAX ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'KMAX must be smaller than ', FZ_MAX, '. Check!', KMAX
             call PRC_abort
          endif

          if ( minval(FZ(1:KMAX)) <= 0.0_RP ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'FZ must be positive. Check! minval(FZ(1:KMAX))=', minval(FZ(1:KMAX))
             call PRC_abort
          endif

          use_user_input = .true.
       endif
    end if

    if ( use_user_input ) then ! input from namelist

       ! Z-direction
       ! calculate buffer grid size

       kbuff = -1
       if ( present(BUFFER_NZ) ) kbuff = BUFFER_NZ
       if ( kbuff > 0 ) then
          if ( kbuff > KMAX ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer grid size (', kbuff, &
                        ') must be smaller than global domain size (Z). Use smaller BUFFER_NZ!'
             call PRC_abort
          endif

          bufftotz = 0.0_RP
          do k = KMAX, KMAX-kbuff+1, -1
             bufftotz = bufftotz + ( FZ(k) - FZ(k-1) )
          enddo
          kmain = KMAX - kbuff
       else if ( present(BUFFER_DZ) ) then
          if ( BUFFER_DZ > FZ(KMAX) ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer length (', BUFFER_DZ, &
                        '[m]) must be smaller than global domain size (Z). Use smaller BUFFER_DZ!'
             call PRC_abort
          endif

          bufftotz = 0.0_RP
          do k = KMAX, 2, -1
             if( bufftotz >= BUFFER_DZ ) exit
             bufftotz = bufftotz + ( FZ(k) - FZ(k-1) )
          enddo
          kbuff = KMAX - k
          kmain = k
       else
          bufftotz = 0.0_RP
          kbuff = 0
          kmain = KMAX
       endif

       ! vertical coordinate (local=global domain)
       ATMOS_GRID_CARTESC_FZ(KS-1) = 0.0_RP

       dz_tmp = FZ(1)
       do k = KS-2, 0, -1
          ATMOS_GRID_CARTESC_FZ(k) = ATMOS_GRID_CARTESC_FZ(k+1) - dz_tmp
       enddo

       do k = KS, KE
          ATMOS_GRID_CARTESC_FZ(k) = FZ(k-KS+1)
       enddo

       dz_tmp = FZ(KMAX) - FZ(KMAX-1)
       do k = KE+1, KA
          ATMOS_GRID_CARTESC_FZ(k) = ATMOS_GRID_CARTESC_FZ(k-1) + dz_tmp
       enddo

       do k = 1, KA
          ATMOS_GRID_CARTESC_CZ(k) = 0.5_RP * ( ATMOS_GRID_CARTESC_FZ(k)+ATMOS_GRID_CARTESC_FZ(k-1) )
       enddo

    else ! calc using DZ

       ! Z-direction
       ! calculate buffer grid size

       fact = -1.0_RP
       if ( present(BUFFFACT_Z) ) fact = BUFFFACT_Z
       if ( fact < 0.0_RP .and. present(BUFFFACT) ) fact = BUFFFACT
       if ( fact < 0.0_RP ) fact = 1.0_RP

       buffz(0) = DZ
       bufftotz = 0.0_RP
       kbuff = -1
       if ( present(BUFFER_NZ) ) kbuff = BUFFER_NZ
       if ( kbuff > 0 ) then
          if ( kbuff > KMAX ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer grid size (', kbuff, &
                        ') must be smaller than global domain size (Z). Use smaller BUFFER_NZ!'
             call PRC_abort
          endif

          do k = 1, kbuff
             buffz(k) = buffz(k-1) * fact
             bufftotz = bufftotz + buffz(k)
          enddo
          kmain = KMAX - kbuff
       else if ( present(BUFFER_DZ) ) then
          do k = 1, KA
             if( bufftotz >= BUFFER_DZ ) exit
             buffz(k) = buffz(k-1) * fact
             bufftotz = bufftotz + buffz(k)
          enddo
          kbuff = k - 1
          kmain = KMAX - kbuff

          if ( kmain < 0 ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_generate",*) 'Buffer length (', bufftotz, &
                        '[m]) must be smaller than global domain size (Z). Use smaller BUFFER_DZ!'
             call PRC_abort
          endif
       else
          kbuff = 0
          kmain = KMAX
       endif

       ! vertical coordinate (local=global domain)
       ATMOS_GRID_CARTESC_FZ(KS-1) = 0.0_RP
       do k = KS-2, 0, -1
          ATMOS_GRID_CARTESC_FZ(k) = ATMOS_GRID_CARTESC_FZ(k+1) - DZ
       enddo

       do k = 1, KS-1
          ATMOS_GRID_CARTESC_CZ(k) = 0.5_RP * ( ATMOS_GRID_CARTESC_FZ(k)+ATMOS_GRID_CARTESC_FZ(k-1) )
       enddo

       do k = KS, KS+kmain-1
          ATMOS_GRID_CARTESC_FZ(k) = ATMOS_GRID_CARTESC_FZ(k-1) + DZ
          ATMOS_GRID_CARTESC_CZ(k) = 0.5_RP * ( ATMOS_GRID_CARTESC_FZ(k)+ATMOS_GRID_CARTESC_FZ(k-1) )
       enddo

       if ( kbuff > 0 ) then
          do k = KS+kmain, KE
             ATMOS_GRID_CARTESC_FZ(k) = ATMOS_GRID_CARTESC_FZ(k-1) + buffz(k-KS-kmain+1)
             ATMOS_GRID_CARTESC_CZ(k) = 0.5_RP * ( ATMOS_GRID_CARTESC_FZ(k)+ATMOS_GRID_CARTESC_FZ(k-1) )
          enddo
       endif

       do k = KE+1, KA
          ATMOS_GRID_CARTESC_FZ(k) = ATMOS_GRID_CARTESC_FZ(k-1) + buffz(kbuff)
          ATMOS_GRID_CARTESC_CZ(k) = 0.5_RP * ( ATMOS_GRID_CARTESC_FZ(k)+ATMOS_GRID_CARTESC_FZ(k-1) )
       enddo

    endif

    ! calc buffer factor (global domain)
    ATMOS_GRID_CARTESC_CBFZ(:) = 0.0_RP
    ATMOS_GRID_CARTESC_FBFZ(:) = 0.0_RP
    if ( kbuff > 0 ) then
       do k = KS+kmain, KE
          ATMOS_GRID_CARTESC_CBFZ(k) = (bufftotz-ATMOS_GRID_CARTESC_FZ(KE)+ATMOS_GRID_CARTESC_CZ(k)) / bufftotz
          ATMOS_GRID_CARTESC_FBFZ(k) = (bufftotz-ATMOS_GRID_CARTESC_FZ(KE)+ATMOS_GRID_CARTESC_FZ(k)) / bufftotz
       enddo
    endif

    do k = KE+1, KA
       ATMOS_GRID_CARTESC_CBFZ(k) = 1.0_RP
       ATMOS_GRID_CARTESC_FBFZ(k) = 1.0_RP
    enddo
    ATMOS_GRID_CARTESC_CBFZ(:) = max( min( ATMOS_GRID_CARTESC_CBFZ(:), 1.0_RP ), 0.0_RP )
    ATMOS_GRID_CARTESC_FBFZ(:) = max( min( ATMOS_GRID_CARTESC_FBFZ(:), 1.0_RP ), 0.0_RP )

    deallocate( buffz )

    ! vertical coordinate (local domain)
    do k = 1, KA
       ATMOS_GRID_CARTESC_CDZ (k) = ATMOS_GRID_CARTESC_FZ(k) - ATMOS_GRID_CARTESC_FZ(k-1)
       ATMOS_GRID_CARTESC_RCDZ(k) = 1.0_RP / ATMOS_GRID_CARTESC_CDZ(k)
    enddo

    do k = 1, KA-1
       ATMOS_GRID_CARTESC_FDZ (k) = ATMOS_GRID_CARTESC_CZ(k+1)-ATMOS_GRID_CARTESC_CZ(k)
       ATMOS_GRID_CARTESC_RFDZ(k) = 1.0_RP / ATMOS_GRID_CARTESC_FDZ(k)
    enddo

    ! X-direction
    ! horizontal coordinate (local domain)
    do i = 0, IA
       ii = i + PRC_2Drank(PRC_myrank,1) * IMAX

       ATMOS_GRID_CARTESC_FX(i) = ATMOS_GRID_CARTESC_FXG(ii)
    enddo

    ii = PRC_2Drank(PRC_myrank,1) * IMAX
    ATMOS_GRID_CARTESC_FBFX(0) = ATMOS_GRID_CARTESC_FBFXG(ii)
    do i = 1, IA
       ii = i + PRC_2Drank(PRC_myrank,1) * IMAX

       ATMOS_GRID_CARTESC_CX  (i) = ATMOS_GRID_CARTESC_CXG  (ii)
       ATMOS_GRID_CARTESC_CBFX(i) = ATMOS_GRID_CARTESC_CBFXG(ii)
       ATMOS_GRID_CARTESC_FBFX(i) = ATMOS_GRID_CARTESC_FBFXG(ii)

       ATMOS_GRID_CARTESC_CDX (i) = ATMOS_GRID_CARTESC_FX(i) - ATMOS_GRID_CARTESC_FX(i-1)
       ATMOS_GRID_CARTESC_RCDX(i) = 1.0_RP / ATMOS_GRID_CARTESC_CDX(i)
    enddo

    do i = 1, IA-1
       ATMOS_GRID_CARTESC_FDX (i) = ATMOS_GRID_CARTESC_CX(i+1)-ATMOS_GRID_CARTESC_CX(i)
       ATMOS_GRID_CARTESC_RFDX(i) = 1.0_RP / ATMOS_GRID_CARTESC_FDX(i)
    enddo

    ! Y-direction
    ! horizontal coordinate (local domain)
    do j = 0, JA
       jj = j + PRC_2Drank(PRC_myrank,2) * JMAX

       ATMOS_GRID_CARTESC_FY(j) = ATMOS_GRID_CARTESC_FYG(jj)
    enddo

    jj = PRC_2Drank(PRC_myrank,2) * JMAX
    ATMOS_GRID_CARTESC_FBFY(0) = ATMOS_GRID_CARTESC_FBFYG(jj)
    do j = 1, JA
       jj = j + PRC_2Drank(PRC_myrank,2) * JMAX

       ATMOS_GRID_CARTESC_CY  (j) = ATMOS_GRID_CARTESC_CYG  (jj)
       ATMOS_GRID_CARTESC_CBFY(j) = ATMOS_GRID_CARTESC_CBFYG(jj)
       ATMOS_GRID_CARTESC_FBFY(j) = ATMOS_GRID_CARTESC_FBFYG(jj)

       ATMOS_GRID_CARTESC_CDY (j) = ATMOS_GRID_CARTESC_FY(j) - ATMOS_GRID_CARTESC_FY(j-1)
       ATMOS_GRID_CARTESC_RCDY(j) = 1.0_RP / ATMOS_GRID_CARTESC_CDY(j)
    enddo

    do j = 1, JA-1
       ATMOS_GRID_CARTESC_FDY (j) = ATMOS_GRID_CARTESC_CY(j+1)-ATMOS_GRID_CARTESC_CY(j)
       ATMOS_GRID_CARTESC_RFDY(j) = 1.0_RP / ATMOS_GRID_CARTESC_FDY(j)
    enddo

    ATMOS_GRID_CARTESC_DOMAIN_CENTER_X = 0.5_RP * ( ATMOS_GRID_CARTESC_FXG(IHALO) + ATMOS_GRID_CARTESC_FXG(IAG-IHALO) )
    ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y = 0.5_RP * ( ATMOS_GRID_CARTESC_FYG(JHALO) + ATMOS_GRID_CARTESC_FYG(JAG-JHALO) )

    ! report
    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_generate",*)                   'Grid information '
    LOG_INFO_CONT('(1x,A,3(1x,F9.3))') 'delta Z, X, Y [m]        :', DZ, DX, DY

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_generate",*)                'Main/buffer Grid (global) :'
    LOG_INFO_CONT('(1x,2(A,I6))')   'Z: buffer = ', kbuff,' x 1, main = ',kmain
    LOG_INFO_CONT('(1x,2(A,I6))')   'X: buffer = ', ibuff,' x 2, main = ',imain
    LOG_INFO_CONT('(1x,2(A,I6))')   'Y: buffer = ', jbuff,' x 2, main = ',jmain

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_generate",*)                'Domain size [km] (global) :'
    LOG_INFO_CONT('(1x,7(A,F9.3))') 'Z:', &
                                                  ATMOS_GRID_CARTESC_FZ(0)       *1.E-3_RP, ' -HALO-                    ',   &
                                                  ATMOS_GRID_CARTESC_FZ(KS-1)    *1.E-3_RP, ' | ',        &
                                                  ATMOS_GRID_CARTESC_CZ(KS)      *1.E-3_RP, ' - ',        &
                                                  ATMOS_GRID_CARTESC_CZ(KE-kbuff)*1.E-3_RP, ' | ',        &
                                                  ATMOS_GRID_CARTESC_FZ(KE-kbuff)*1.E-3_RP, ' -buffer- ', &
                                                  ATMOS_GRID_CARTESC_FZ(KE)      *1.E-3_RP, ' -HALO- ',   &
                                                  ATMOS_GRID_CARTESC_FZ(KA)      *1.E-3_RP
    LOG_INFO_CONT('(1x,8(A,F9.3))') 'X:', &
                                                  ATMOS_GRID_CARTESC_FXG(0)              *1.E-3_RP, ' -HALO- ',   &
                                                  ATMOS_GRID_CARTESC_FXG(IHALO)          *1.E-3_RP, ' -buffer- ', &
                                                  ATMOS_GRID_CARTESC_FXG(IHALO+ibuff)    *1.E-3_RP, ' | ',        &
                                                  ATMOS_GRID_CARTESC_CXG(IHALO+ibuff+1)  *1.E-3_RP, ' - ',        &
                                                  ATMOS_GRID_CARTESC_CXG(IAG-IHALO-ibuff)*1.E-3_RP, ' | ',        &
                                                  ATMOS_GRID_CARTESC_FXG(IAG-IHALO-ibuff)*1.E-3_RP, ' -buffer- ', &
                                                  ATMOS_GRID_CARTESC_FXG(IAG-IHALO)      *1.E-3_RP, ' -HALO- ',   &
                                                  ATMOS_GRID_CARTESC_FXG(IAG)            *1.E-3_RP
    LOG_INFO_CONT('(1x,8(A,F9.3))') 'Y:', &
                                                  ATMOS_GRID_CARTESC_FYG(0)              *1.E-3_RP, ' -HALO- ',   &
                                                  ATMOS_GRID_CARTESC_FYG(JHALO)          *1.E-3_RP, ' -buffer- ', &
                                                  ATMOS_GRID_CARTESC_FYG(JHALO+jbuff)    *1.E-3_RP, ' | ',        &
                                                  ATMOS_GRID_CARTESC_CYG(JHALO+jbuff+1)  *1.E-3_RP, ' - ',        &
                                                  ATMOS_GRID_CARTESC_CYG(JAG-JHALO-jbuff)*1.E-3_RP, ' | ',        &
                                                  ATMOS_GRID_CARTESC_FYG(JAG-JHALO-jbuff)*1.E-3_RP, ' -buffer- ', &
                                                  ATMOS_GRID_CARTESC_FYG(JAG-JHALO)      *1.E-3_RP, ' -HALO- ',   &
                                                  ATMOS_GRID_CARTESC_FYG(JAG)            *1.E-3_RP

    return
  end subroutine ATMOS_GRID_CARTESC_generate

  !-----------------------------------------------------------------------------
  !> Output information
  subroutine ATMOS_GRID_CARTESC_output_info
    integer :: k

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_output_info",*)                'Center Position of Grid (global) :'
    LOG_INFO_CONT('(1x,A,F12.3)')   'X: ', ATMOS_GRID_CARTESC_DOMAIN_CENTER_X
    LOG_INFO_CONT('(1x,A,F12.3)')   'Y: ', ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y


    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_output_info",*)                'Domain size [km] (local) :'
    LOG_INFO_CONT('(1x,6(A,F9.3))') 'X:',                                                          &
                                                  ATMOS_GRID_CARTESC_FX(0) *1.E-3_RP, ' -HALO- ', ATMOS_GRID_CARTESC_FX(IS-1)*1.E-3_RP, ' | ', &
                                                  ATMOS_GRID_CARTESC_CX(IS)*1.E-3_RP, ' - ',      ATMOS_GRID_CARTESC_CX(IE)  *1.E-3_RP, ' | ', &
                                                  ATMOS_GRID_CARTESC_FX(IE)*1.E-3_RP, ' -HALO- ', ATMOS_GRID_CARTESC_FX(IA)  *1.E-3_RP
    LOG_INFO_CONT('(1x,6(A,F9.3))') 'Y:',                    &
                                                  ATMOS_GRID_CARTESC_FY(0) *1.E-3_RP, ' -HALO- ', ATMOS_GRID_CARTESC_FY(JS-1)*1.E-3_RP, ' | ', &
                                                  ATMOS_GRID_CARTESC_CY(JS)*1.E-3_RP, ' - ',      ATMOS_GRID_CARTESC_CY(JE)  *1.E-3_RP, ' | ', &
                                                  ATMOS_GRID_CARTESC_FY(JE)*1.E-3_RP, ' -HALO- ', ATMOS_GRID_CARTESC_FY(JA)  *1.E-3_RP


    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_output_info",'(1x,A)') 'Vertical Coordinate'
    LOG_INFO_CONT('(1x,A)') '|===============================================|'
    LOG_INFO_CONT('(1x,A)') '|    k        z       zh       dz   buffer    k |'
    LOG_INFO_CONT('(1x,A)') '|           [m]      [m]      [m]   factor      |'

    do k = KA, KE+1, -1
    LOG_INFO_CONT('(1x,A,F9.2,A,F9.2,I5,A)')  '|              ',ATMOS_GRID_CARTESC_FZ(k),'         ', ATMOS_GRID_CARTESC_FBFZ(k),k,' |'
    LOG_INFO_CONT('(1x,A,I5,F9.2,A,2F9.2,A)') '|',k,ATMOS_GRID_CARTESC_CZ(k),'         ',ATMOS_GRID_CARTESC_CDZ(k), ATMOS_GRID_CARTESC_CBFZ(k),'      |'
    enddo

    k = KE
    LOG_INFO_CONT('(1x,A,F9.2,A,F9.2,I5,A)')  '|              ',ATMOS_GRID_CARTESC_FZ(k),'         ', ATMOS_GRID_CARTESC_FBFZ(k),k,' | KE = TOA'
    LOG_INFO_CONT('(1x,A,I5,F9.2,A,2F9.2,A)') '|',k,ATMOS_GRID_CARTESC_CZ(k),'         ',ATMOS_GRID_CARTESC_CDZ(k), ATMOS_GRID_CARTESC_CBFZ(k),'      |'

    do k = KE-1, KS, -1
    LOG_INFO_CONT('(1x,A,F9.2,A,F9.2,I5,A)')  '|              ',ATMOS_GRID_CARTESC_FZ(k),'         ', ATMOS_GRID_CARTESC_FBFZ(k),k,' |'
    LOG_INFO_CONT('(1x,A,I5,F9.2,A,2F9.2,A)') '|',k,ATMOS_GRID_CARTESC_CZ(k),'         ',ATMOS_GRID_CARTESC_CDZ(k), ATMOS_GRID_CARTESC_CBFZ(k),'      |'
    enddo

    k = KS-1
    LOG_INFO_CONT('(1x,A,F9.2,A,F9.2,I5,A)')  '|              ',ATMOS_GRID_CARTESC_FZ(k),'         ', ATMOS_GRID_CARTESC_FBFZ(k),k,' | KS-1 = surface'
    LOG_INFO_CONT('(1x,A,I5,F9.2,A,2F9.2,A)') '|',k,ATMOS_GRID_CARTESC_CZ(k),'         ',ATMOS_GRID_CARTESC_CDZ(k), ATMOS_GRID_CARTESC_CBFZ(k),'      |'

    do k = KS-2, 1, -1
    LOG_INFO_CONT('(1x,A,F9.2,A,F9.2,I5,A)')  '|              ',ATMOS_GRID_CARTESC_FZ(k),'         ', ATMOS_GRID_CARTESC_FBFZ(k),k,' |'
    LOG_INFO_CONT('(1x,A,I5,F9.2,A,2F9.2,A)') '|',k,ATMOS_GRID_CARTESC_CZ(k),'         ',ATMOS_GRID_CARTESC_CDZ(k), ATMOS_GRID_CARTESC_CBFZ(k),'      |'
    enddo

    k = 0
    LOG_INFO_CONT('(1x,A,F9.2,A,F9.2,I5,A)') '|              ',ATMOS_GRID_CARTESC_FZ(k),'         ', ATMOS_GRID_CARTESC_FBFZ(k),k,' |'

    LOG_INFO_CONT('(1x,A)') '|===============================================|'


    return
  end subroutine ATMOS_GRID_CARTESC_output_info

end module scale_atmos_grid_cartesC
