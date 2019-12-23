!-------------------------------------------------------------------------------
!> module TOPOGRAPHY
!!
!! @par Description
!!          Topography module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_topography
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TOPOGRAPHY_setup
  public :: TOPOGRAPHY_fillhalo
  public :: TOPOGRAPHY_write
  public :: TOPOGRAPHY_calc_tan_slope

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: TOPOGRAPHY_exist = .false. !< topography exists?

  real(RP), public, allocatable :: TOPOGRAPHY_Zsfc   (:,:) !< absolute ground height [m]
  real(RP), public, allocatable :: TOPOGRAPHY_TanSL_X(:,:) !< tan(slope_x)
  real(RP), public, allocatable :: TOPOGRAPHY_TanSL_Y(:,:) !< tan(slope_y)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: TOPOGRAPHY_read

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: TOPOGRAPHY_IN_BASENAME = ''                     !< basename of the input  file
  character(len=H_LONG),  private :: TOPOGRAPHY_IN_VARNAME  = 'topo'                 !< variable name of topo in the input  file
  logical,                private :: TOPOGRAPHY_IN_AGGREGATE                          !> switch to use aggregated file
  logical,                private :: TOPOGRAPHY_IN_CHECK_COORDINATES = .false.        !> switch for check of coordinates
  character(len=H_LONG),  private :: TOPOGRAPHY_OUT_BASENAME = ''                     !< basename of the output file
  logical,                private :: TOPOGRAPHY_OUT_AGGREGATE                         !> switch to use aggregated file
  character(len=H_MID),   private :: TOPOGRAPHY_OUT_TITLE    = 'SCALE-RM TOPOGRAPHY'  !< title    of the output file
  character(len=H_SHORT), private :: TOPOGRAPHY_OUT_DTYPE    = 'DEFAULT'              !< REAL4 or REAL8

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine TOPOGRAPHY_setup
    use scale_file, only: &
       FILE_AGGREGATE
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_TOPOGRAPHY / &
       TOPOGRAPHY_IN_BASENAME,          &
       TOPOGRAPHY_IN_VARNAME,           &
       TOPOGRAPHY_IN_AGGREGATE,         &
       TOPOGRAPHY_IN_CHECK_COORDINATES, &
       TOPOGRAPHY_OUT_BASENAME,         &
       TOPOGRAPHY_OUT_AGGREGATE,        &
       TOPOGRAPHY_OUT_DTYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("TOPOGRAPHY_setup",*) 'Setup'

    TOPOGRAPHY_IN_AGGREGATE  = FILE_AGGREGATE
    TOPOGRAPHY_OUT_AGGREGATE = FILE_AGGREGATE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPOGRAPHY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("TOPOGRAPHY_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("TOPOGRAPHY_setup",*) 'Not appropriate names in namelist PARAM_TOPOGRAPHY. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TOPOGRAPHY)

    allocate( TOPOGRAPHY_Zsfc   (IA,JA) )
    allocate( TOPOGRAPHY_TanSL_X(IA,JA) )
    allocate( TOPOGRAPHY_TanSL_Y(IA,JA) )
    TOPOGRAPHY_Zsfc(:,:) = 0.0_RP
    TOPOGRAPHY_TanSL_X(:,:) = 0.0_RP
    TOPOGRAPHY_TanSL_Y(:,:) = 0.0_RP

    ! read from file
    call TOPOGRAPHY_read

    return
  end subroutine TOPOGRAPHY_setup

  !-----------------------------------------------------------------------------
  !> HALO Communication
  subroutine TOPOGRAPHY_fillhalo( Zsfc, FILL_BND )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(inout), optional :: Zsfc(IA,JA)
    logical,  intent(in),    optional :: FILL_BND

    logical :: FILL_BND_
    !---------------------------------------------------------------------------

    FILL_BND_ = .false.
    if ( present(FILL_BND) ) FILL_BND_ = FILL_BND

    if ( present(Zsfc) ) then
       call COMM_vars8( Zsfc(:,:), 1 )
       call COMM_wait ( Zsfc(:,:), 1, FILL_BND_ )
    else
       call COMM_vars8( TOPOGRAPHY_Zsfc(:,:), 1 )
       call COMM_wait ( TOPOGRAPHY_Zsfc(:,:), 1, FILL_BND_ )
    end if

    return
  end subroutine TOPOGRAPHY_fillhalo

  !-----------------------------------------------------------------------------
  !> Read topography
  subroutine TOPOGRAPHY_read
    use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_read, &
       FILE_CARTESC_flush, &
       FILE_CARTESC_check_coordinates, &
       FILE_CARTESC_close
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer :: fid
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("TOPOGRAPHY_read",*) 'Input topography file '

    if ( TOPOGRAPHY_IN_BASENAME /= '' ) then

       call FILE_CARTESC_open( TOPOGRAPHY_IN_BASENAME, fid, aggregate=TOPOGRAPHY_IN_AGGREGATE )
       call FILE_CARTESC_read( fid, TOPOGRAPHY_IN_VARNAME, 'XY', TOPOGRAPHY_Zsfc(:,:) )

       call FILE_CARTESC_flush( fid )

       if ( TOPOGRAPHY_IN_CHECK_COORDINATES ) then
          call FILE_CARTESC_check_coordinates( fid )
       end if

       call FILE_CARTESC_close( fid )

       call TOPOGRAPHY_fillhalo( FILL_BND=.false. )

       TOPOGRAPHY_exist = .true.

    else
       LOG_INFO_CONT(*) 'topography file is not specified.'

       TOPOGRAPHY_exist = .false.
    endif

    return
  end subroutine TOPOGRAPHY_read

  !-----------------------------------------------------------------------------
  !> Write topography
  subroutine TOPOGRAPHY_write
    use scale_file_cartesC, only: &
       FILE_CARTESC_write
    implicit none
    !---------------------------------------------------------------------------

    if ( TOPOGRAPHY_OUT_BASENAME /= '' .and. TOPOGRAPHY_OUT_BASENAME /= TOPOGRAPHY_IN_BASENAME ) then

       LOG_NEWLINE
       LOG_INFO("TOPOGRAPHY_write",*) 'Output topography file '

       call TOPOGRAPHY_fillhalo( FILL_BND=.false. )

       call FILE_CARTESC_write( TOPOGRAPHY_Zsfc(:,:),                                    & ! [IN]
                                TOPOGRAPHY_OUT_BASENAME, TOPOGRAPHY_OUT_TITLE,           & ! [IN]
                                'topo', 'Topography', 'm', 'XY',   TOPOGRAPHY_OUT_DTYPE, & ! [IN]
                                standard_name="surface_altitude",                        & ! [IN]
                                haszcoord=.false., aggregate=TOPOGRAPHY_OUT_AGGREGATE    ) ! [IN]

    endif

    return
  end subroutine TOPOGRAPHY_write

  subroutine TOPOGRAPHY_calc_tan_slope( &
       IA, IS, IE, JA, JS, JE, &
       RCDX, RCDY, MAPF )
       integer,  intent(in) :: IA, IS, IE
       integer,  intent(in) :: JA, JS, JE
       real(RP), intent(in) :: RCDX(IA), RCDY(JA)
       real(RP), intent(in) :: MAPF(IA,JA,2)

       integer :: i, j

       do j = JS, JE
       do i = IS, IE
          TOPOGRAPHY_TanSL_X(i,j) = ( ( TOPOGRAPHY_Zsfc(i+1,j) + TOPOGRAPHY_Zsfc(i  ,j) ) * 0.5_RP &
                                    - ( TOPOGRAPHY_Zsfc(i  ,j) + TOPOGRAPHY_Zsfc(i-1,j) ) * 0.5_RP ) &
                                  * RCDX(i) * MAPF(i,j,1)
          TOPOGRAPHY_TanSL_Y(i,j) = ( ( TOPOGRAPHY_Zsfc(i,j+1) + TOPOGRAPHY_Zsfc(i,j  ) ) * 0.5_RP &
                                    - ( TOPOGRAPHY_Zsfc(i,j  ) + TOPOGRAPHY_Zsfc(i,j-1) ) * 0.5_RP ) &
                                  * RCDY(j) * MAPF(i,j,2)
       end do
       end do

       call TOPOGRAPHY_fillhalo( TOPOGRAPHY_TanSL_X(:,:), .true. )
       call TOPOGRAPHY_fillhalo( TOPOGRAPHY_TanSL_Y(:,:), .true. )

       return
     end subroutine TOPOGRAPHY_calc_tan_slope

end module scale_topography
