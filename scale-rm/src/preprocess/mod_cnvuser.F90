!-------------------------------------------------------------------------------
!> module Convert 2D user data
!!
!! @par Description
!!          subroutines for preparing 2D data
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cnvuser
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CNVUSER_setup
  public :: CNVUSER

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNVUSER_write

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: CNVUSER_FILE_TYPE           = ''                     ! '' : do nothing
                                                                                          ! 'TILE': tile data
                                                                                          ! 'GrADS': GrADS data

  character(len=H_LONG),  private :: CNVUSER_OUT_BASENAME        = ''                     ! basename of the output file
  character(len=H_MID),   private :: CNVUSER_OUT_TITLE           = 'SCALE-RM 2D Boundary' ! title    of the output file
  character(len=H_SHORT), private :: CNVUSER_OUT_VARNAME         = ''                     ! name  of the variable
  character(len=H_MID),   private :: CNVUSER_OUT_VARDESC         = ''                     ! title of the variable
  character(len=H_SHORT), private :: CNVUSER_OUT_VARUNIT         = ''                     ! units of the variable
  character(len=H_SHORT), private :: CNVUSER_OUT_DTYPE           = 'DEFAULT'              ! REAL4 or REAL8
  real(DP),               private :: CNVUSER_OUT_DT              = -1_DP                  ! sec
  logical,                private :: CNVUSER_OUT_AGGREGATE

  integer               , private :: CNVUSER_NSTEPS              = 1                      ! # of time steps

  ! TILE data
  character(len=H_SHORT), private :: CNVUSER_TILE_DTYPE          = 'real4'                ! data type in the tiled data
  real(RP),               private :: CNVUSER_TILE_DLAT
  real(RP),               private :: CNVUSER_TILE_DLON
  character(len=H_LONG),  private :: CNVUSER_TILE_DIR            = ''
  character(len=H_LONG),  private :: CNVUSER_TILE_CATALOGUE      = ''

  ! GrADS data

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNVUSER_setup
    use scale_prc, only: &
       PRC_abort
    use scale_file, only: &
       FILE_AGGREGATE
    use mod_cnv2d, only: &
       CNV2D_tile_init, &
       CNV2D_grads_init
    implicit none

    character(len=H_SHORT) :: CNVUSER_INTERP_TYPE  = 'LINEAR'
    integer                :: CNVUSER_INTERP_LEVEL = 5

    character(len=H_LONG)  :: CNVUSER_GrADS_FILENAME = ''
    character(len=H_SHORT) :: CNVUSER_GrADS_VARNAME  = ''
    character(len=H_SHORT) :: CNVUSER_GrADS_LATNAME  = 'lat'
    character(len=H_SHORT) :: CNVUSER_GrADS_LONNAME  = 'lon'

    namelist / PARAM_CNVUSER / &
       CNVUSER_FILE_TYPE,      &
       CNVUSER_NSTEPS,         &
       CNVUSER_INTERP_TYPE,    &
       CNVUSER_INTERP_LEVEL,   &
       CNVUSER_TILE_DTYPE,     &
       CNVUSER_TILE_DLAT,      &
       CNVUSER_TILE_DLON,      &
       CNVUSER_TILE_DIR,       &
       CNVUSER_TILE_CATALOGUE, &
       CNVUSER_GrADS_FILENAME, &
       CNVUSER_GrADS_VARNAME,  &
       CNVUSER_GrADS_LATNAME,  &
       CNVUSER_GrADS_LONNAME,  &
       CNVUSER_OUT_BASENAME,   &
       CNVUSER_OUT_TITLE,      &
       CNVUSER_OUT_VARNAME,    &
       CNVUSER_OUT_VARDESC,    &
       CNVUSER_OUT_VARUNIT,    &
       CNVUSER_OUT_DTYPE,      &
       CNVUSER_OUT_DT

    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVUSER_setup",*) 'Setup'

    CNVUSER_OUT_AGGREGATE = FILE_AGGREGATE

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVUSER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVUSER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVUSER_setup",*) 'Not appropriate names in namelist PARAM_CNVUSER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVUSER)

    select case ( CNVUSER_FILE_TYPE )
    case ( '' )
       ! do nothing
    case ( 'TILE' )
       call CNV2D_tile_init( CNVUSER_TILE_DTYPE,                   &
                             CNVUSER_TILE_DLAT, CNVUSER_TILE_DLON, &
                             CNVUSER_TILE_DIR,                     &
                             CNVUSER_TILE_CATALOGUE,               &
                             CNVUSER_INTERP_TYPE,                  &
                             interp_level = CNVUSER_INTERP_LEVEL   )
    case ( 'GrADS' )
       call CNV2D_grads_init( CNVUSER_GrADS_FILENAME,             &
                              CNVUSER_GrADS_VARNAME,              &
                              CNVUSER_GrADS_LATNAME,              &
                              CNVUSER_GrADS_LONNAME,              &
                              CNVUSER_INTERP_TYPE,                &
                              interp_level = CNVUSER_INTERP_LEVEL )
    case default
       if ( CNVUSER_FILE_TYPE .ne. '' .and. CNVUSER_OUT_BASENAME == '' ) then
          LOG_ERROR('CNVUSER_setup',*) 'CNVUSER_OUT_BASENAME is required'
          call PRC_abort
       end if
    end select

    return
  end subroutine CNVUSER_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVUSER
    use scale_file_cartesc, only: &
       FILE_CARTESC_create, &
       FILE_CARTESC_def_var, &
       FILE_CARTESC_enddef, &
       FILE_CARTESC_close
    use scale_time, only: &
       TIME_NOWDATE
    use mod_cnv2d, only: &
       CNV2D_exec
    implicit none

    real(RP) :: var(IA,JA)

    integer :: fid, vid
    integer :: step

    if ( CNVUSER_FILE_TYPE == '' ) return

    call FILE_CARTESC_create( CNVUSER_OUT_BASENAME,             & ! [IN]
                              CNVUSER_OUT_TITLE,                & ! [IN]
                              CNVUSER_OUT_DTYPE,                & ! [IN]
                              fid,                              & ! [OUT]
                              date = TIME_NOWDATE,              & ! [IN]
                              haszcoord = .false.,              & ! [IN]
                              aggregate = CNVUSER_OUT_AGGREGATE ) ! [IN]

    call FILE_CARTESC_def_var( fid,                       & ! [IN]
                               CNVUSER_OUT_VARNAME,       & ! [IN]
                               CNVUSER_OUT_VARDESC,       & ! [IN]
                               CNVUSER_OUT_VARUNIT,       & ! [IN]
                               'XYT',                     & ! [IN]
                               CNVUSER_OUT_DTYPE,         & ! [IN]
                               vid,                       & ! [OUT]
                               timeintv = CNVUSER_OUT_DT, & ! [IN]
                               nsteps = CNVUSER_NSTEPS    ) ! [IN]

    call FILE_CARTESC_enddef(fid)

    do step = 1, CNVUSER_NSTEPS

       LOG_PROGRESS(*) 'step = ', step

       call CNV2D_exec( var(:,:), step = step )

       call CNVUSER_write( fid, vid, var(:,:), CNVUSER_OUT_DT, step )

    end do

    call FILE_CARTESC_close( fid )

    return
  end subroutine CNVUSER

  ! private

  subroutine CNVUSER_write( &
       fid, vid,  &
       VAR,       &
       timeintv,  &
       istep      )
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
    implicit none
    integer,  intent(in) :: fid, vid
    real(RP), intent(in) :: VAR(IA,JA,1)
    real(DP), intent(in) :: timeintv
    integer,  intent(in) :: istep

    real(DP) :: timeofs
    !---------------------------------------------------------------------------

    timeofs = ( istep - 1 ) * timeintv
    call FILE_CARTESC_write_var( fid, vid,            & ! [IN]
                                 var(:,:,:),          & ! [IN]
                                 CNVUSER_OUT_VARNAME, & ! [IN]
                                 'XYT',               & ! [IN]
                                 timeintv,            & ! [IN]
                                 timetarg = 1,        & ! [IN]
                                 timeofs = timeofs    ) ! [IN]

    return
  end subroutine CNVUSER_write


end module mod_cnvuser
