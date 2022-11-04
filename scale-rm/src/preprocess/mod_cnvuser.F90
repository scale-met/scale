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
  use scale_prc, only: &
     PRC_abort
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
  integer, parameter     :: n_vars_max                              = 10                       ! Maximum number of variables
  integer                :: CNVUSER_NVARS                           = 1
  character(len=H_SHORT) :: CNVUSER_FILE_TYPE(n_vars_max)           = ''                       ! '' : do nothing
                                                                                               ! 'TILE': tile data
                                                                                               ! 'GrADS': GrADS data
  character(len=H_SHORT) :: CNVUSER_INTERP_TYPE(n_vars_max)         = 'LINEAR'
  integer                :: CNVUSER_INTERP_LEVEL(n_vars_max)        = 5

  character(len=H_LONG)  :: CNVUSER_OUT_BASENAME(n_vars_max)        = ''                       ! basename of the output file
  character(len=H_MID)   :: CNVUSER_OUT_TITLE(n_vars_max)           = 'SCALE-RM User Boundary' ! title    of the output file
  character(len=H_SHORT) :: CNVUSER_OUT_VARNAME(n_vars_max)         = ''                       ! name  of the variable
  character(len=H_MID)   :: CNVUSER_OUT_VARDESC(n_vars_max)         = ''                       ! title of the variable
  character(len=H_SHORT) :: CNVUSER_OUT_VARUNIT(n_vars_max)         = ''                       ! units of the variable
  character(len=H_SHORT) :: CNVUSER_OUT_DTYPE(n_vars_max)           = 'DEFAULT'                ! REAL4 or REAL8
  real(DP)               :: CNVUSER_OUT_DT(n_vars_max)              = -1_DP                    ! sec
  logical                :: CNVUSER_OUT_AGGREGATE

  integer                :: CNVUSER_NSTEPS(n_vars_max)              = 1                        ! # of time steps

  ! TILE data
  character(len=H_SHORT) :: CNVUSER_TILE_DTYPE(n_vars_max)          = 'real4'                  ! data type in the tiled data
  real(RP)               :: CNVUSER_TILE_DLAT(n_vars_max)
  real(RP)               :: CNVUSER_TILE_DLON(n_vars_max)
  character(len=H_LONG)  :: CNVUSER_TILE_DIR(n_vars_max)            = ''
  character(len=H_LONG)  :: CNVUSER_TILE_CATALOGUE(n_vars_max)      = ''

  ! GrADS data
  character(len=H_LONG)  :: CNVUSER_GrADS_FILENAME(n_vars_max)     = ''
  character(len=H_SHORT) :: CNVUSER_GrADS_VARNAME(n_vars_max)      = ''
  character(len=H_SHORT) :: CNVUSER_GrADS_LATNAME(n_vars_max)      = 'lat'
  character(len=H_SHORT) :: CNVUSER_GrADS_LONNAME(n_vars_max)      = 'lon'
  character(len=H_SHORT) :: CNVUSER_GrADS_LEVNAME(n_vars_max)      = ''
  character(len=H_SHORT) :: CNVUSER_GrADS_HEIGHT_PLEV(n_vars_max)  = 'HGT'

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup

  subroutine CNVUSER_setup
    use scale_file, only: &
       FILE_AGGREGATE

    namelist / PARAM_CNVUSER /    &
       CNVUSER_NVARS,             &
       CNVUSER_FILE_TYPE,         &
       CNVUSER_NSTEPS,            &
       CNVUSER_INTERP_TYPE,       &
       CNVUSER_INTERP_LEVEL,      &
       CNVUSER_TILE_DTYPE,        &
       CNVUSER_TILE_DLAT,         &
       CNVUSER_TILE_DLON,         &
       CNVUSER_TILE_DIR,          &
       CNVUSER_TILE_CATALOGUE,    &
       CNVUSER_GrADS_FILENAME,    &
       CNVUSER_GrADS_VARNAME,     &
       CNVUSER_GrADS_LATNAME,     &
       CNVUSER_GrADS_LONNAME,     &
       CNVUSER_GrADS_LEVNAME,     &
       CNVUSER_GrADS_HEIGHT_PLEV, &
       CNVUSER_OUT_BASENAME,      &
       CNVUSER_OUT_TITLE,         &
       CNVUSER_OUT_VARNAME,       &
       CNVUSER_OUT_VARDESC,       &
       CNVUSER_OUT_VARUNIT,       &
       CNVUSER_OUT_DTYPE,         &
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

  end subroutine CNVUSER_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVUSER
    integer :: i
    !---------------------------------------------------------------------------

    do i = 1, CNVUSER_NVARS
      if ( CNVUSER_FILE_TYPE(i) == '' ) cycle

      select case ( CNVUSER_FILE_TYPE(i) )
        case ( 'TILE' )
          call CNVUSER_prepare_TILE( &
            CNVUSER_TILE_DTYPE(i), &
            CNVUSER_TILE_DLAT(i), &
            CNVUSER_TILE_DLON(i), &
            CNVUSER_TILE_DIR(i), &
            CNVUSER_TILE_CATALOGUE(i), &
            CNVUSER_INTERP_TYPE(i), &
            CNVUSER_INTERP_LEVEL(i) &
            )
        case ( 'GrADS' )
          call CNVUSER_prepare_GrADS( &
            CNVUSER_GrADS_FILENAME(i), &
            CNVUSER_GrADS_VARNAME(i), &
            CNVUSER_GrADS_LATNAME(i), &
            CNVUSER_GrADS_LONNAME(i), &
            CNVUSER_INTERP_TYPE(i), &
            CNVUSER_INTERP_LEVEL(i), &
            CNVUSER_OUT_VARNAME(i) &
            )
        case ( 'GrADS-3D' )
          call CNVUSER_prepare_GrADS_3D( &
            CNVUSER_GrADS_VARNAME(i), &
            CNVUSER_OUT_VARNAME(i) &
            )
        case default
           LOG_ERROR('CNVUSER',*) 'CNVUSER_FILE_TYPE is invalid: ', CNVUSER_FILE_TYPE(i)
           LOG_ERROR_CONT(*)            'It must be "TILE" or "GrADS".'
           call PRC_abort
      end select

      if ( CNVUSER_OUT_BASENAME(i) == '' .or. CNVUSER_OUT_VARNAME(i) == '' ) then
         LOG_INFO("CNVUSER", *) CNVUSER_OUT_BASENAME(1:CNVUSER_NVARS), CNVUSER_OUT_VARNAME(1:CNVUSER_NVARS)
         LOG_ERROR('CNVUSER',*) 'CNVUSER_OUT_BASENAME and CNVUSER_OUT_VARNAME are required'
         call PRC_abort
      end if

      select case ( CNVUSER_FILE_TYPE(i) )
        case ( 'TILE' )
          call CNVUSER_execute_TILE_GrADS( &
            CNVUSER_NSTEPS(i), &
            CNVUSER_OUT_BASENAME(i), &
            CNVUSER_OUT_TITLE(i), &
            CNVUSER_OUT_VARNAME(i), &
            CNVUSER_OUT_VARDESC(i), &
            CNVUSER_OUT_VARUNIT(i), &
            CNVUSER_OUT_DTYPE(i), &
            CNVUSER_OUT_DT(i) &
            )
        case ( 'GrADS' )
          call CNVUSER_execute_TILE_GrADS( &
            CNVUSER_NSTEPS(i), &
            CNVUSER_OUT_BASENAME(i), &
            CNVUSER_OUT_TITLE(i), &
            CNVUSER_OUT_VARNAME(i), &
            CNVUSER_OUT_VARDESC(i), &
            CNVUSER_OUT_VARUNIT(i), &
            CNVUSER_OUT_DTYPE(i), &
            CNVUSER_OUT_DT(i) &
            )
        case ( 'GrADS-3D' )
          call CNVUSER_execute_GrADS_3D( &
            CNVUSER_NSTEPS(i), &
            CNVUSER_GrADS_FILENAME(i), &
            CNVUSER_GrADS_VARNAME(i), &
            CNVUSER_GrADS_LATNAME(i), &
            CNVUSER_GrADS_LONNAME(i), &
            CNVUSER_GrADS_LEVNAME(i), &
            CNVUSER_GrADS_HEIGHT_PLEV(i), &
            CNVUSER_INTERP_TYPE(i), &
            CNVUSER_INTERP_LEVEL(i), &
            CNVUSER_OUT_BASENAME(i), &
            CNVUSER_OUT_TITLE(i), &
            CNVUSER_OUT_VARNAME(i), &
            CNVUSER_OUT_VARDESC(i), &
            CNVUSER_OUT_VARUNIT(i), &
            CNVUSER_OUT_DTYPE(i), &
            CNVUSER_OUT_DT(i) &
            )
      end select
    end do

  end subroutine CNVUSER

  ! private

  subroutine CNVUSER_prepare_TILE( &
    CNVUSER_TILE_DTYPE, &
    CNVUSER_TILE_DLAT, &
    CNVUSER_TILE_DLON, &
    CNVUSER_TILE_DIR, &
    CNVUSER_TILE_CATALOGUE, &
    CNVUSER_INTERP_TYPE, &
    CNVUSER_INTERP_LEVEL &
    )

    use mod_cnv2d, only: &
       CNV2D_tile_init

    character(len=H_SHORT), intent(in) :: CNVUSER_TILE_DTYPE
    real(RP),               intent(in) :: CNVUSER_TILE_DLAT
    real(RP),               intent(in) :: CNVUSER_TILE_DLON
    character(len=H_LONG),  intent(in) :: CNVUSER_TILE_DIR
    character(len=H_LONG),  intent(in) :: CNVUSER_TILE_CATALOGUE
    character(len=H_SHORT), intent(in) :: CNVUSER_INTERP_TYPE
    integer,                intent(in) :: CNVUSER_INTERP_LEVEL

    call CNV2D_tile_init( CNVUSER_TILE_DTYPE,                   &
                          CNVUSER_TILE_DLAT, CNVUSER_TILE_DLON, &
                          CNVUSER_TILE_DIR,                     &
                          CNVUSER_TILE_CATALOGUE,               &
                          CNVUSER_INTERP_TYPE,                  &
                          interp_level = CNVUSER_INTERP_LEVEL   )
  end subroutine CNVUSER_prepare_TILE

  subroutine CNVUSER_prepare_GrADS( &
    CNVUSER_GrADS_FILENAME, &
    CNVUSER_GrADS_VARNAME, &
    CNVUSER_GrADS_LATNAME, &
    CNVUSER_GrADS_LONNAME, &
    CNVUSER_INTERP_TYPE, &
    CNVUSER_INTERP_LEVEL, &
    CNVUSER_OUT_VARNAME &
    )

    use mod_cnv2d, only: &
       CNV2D_grads_init

    character(len=H_LONG),  intent(in) :: CNVUSER_GrADS_FILENAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_VARNAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_LATNAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_LONNAME
    character(len=H_SHORT), intent(in) :: CNVUSER_INTERP_TYPE
    integer,                intent(in) :: CNVUSER_INTERP_LEVEL
    character(len=H_SHORT), intent(inout) :: CNVUSER_OUT_VARNAME

    call CNV2D_grads_init( CNVUSER_GrADS_FILENAME, &
                           CNVUSER_GrADS_VARNAME,  &
                           CNVUSER_GrADS_LATNAME,  &
                           CNVUSER_GrADS_LONNAME,  &
                           CNVUSER_INTERP_TYPE,    &
                           CNVUSER_INTERP_LEVEL    )
    if ( CNVUSER_OUT_VARNAME == '' ) CNVUSER_OUT_VARNAME = CNVUSER_GrADS_VARNAME
  end subroutine CNVUSER_prepare_GrADS

  subroutine CNVUSER_prepare_GrADS_3D( &
    CNVUSER_GrADS_VARNAME, &
    CNVUSER_OUT_VARNAME &
    )

    character(len=H_SHORT), intent(in)    :: CNVUSER_GrADS_VARNAME
    character(len=H_SHORT), intent(inout) :: CNVUSER_OUT_VARNAME

    if ( CNVUSER_OUT_VARNAME == '' ) CNVUSER_OUT_VARNAME = CNVUSER_GrADS_VARNAME
  end subroutine CNVUSER_prepare_GrADS_3D

  subroutine CNVUSER_execute_GrADS_3D( &
    CNVUSER_NSTEPS, &
    CNVUSER_GrADS_FILENAME, &
    CNVUSER_GrADS_VARNAME, &
    CNVUSER_GrADS_LATNAME, &
    CNVUSER_GrADS_LONNAME, &
    CNVUSER_GrADS_LEVNAME, &
    CNVUSER_GrADS_HEIGHT_PLEV, &
    CNVUSER_INTERP_TYPE, &
    CNVUSER_INTERP_LEVEL, &
    CNVUSER_OUT_BASENAME, &
    CNVUSER_OUT_TITLE, &
    CNVUSER_OUT_VARNAME, &
    CNVUSER_OUT_VARDESC, &
    CNVUSER_OUT_VARUNIT, &
    CNVUSER_OUT_DTYPE, &
    CNVUSER_OUT_DT &
    )

    use scale_const, only: &
       PI  => CONST_PI,  &
       D2R => CONST_D2R
    use scale_file_cartesc, only: &
       FILE_CARTESC_create, &
       FILE_CARTESC_def_var, &
       FILE_CARTESC_write_var, &
       FILE_CARTESC_enddef, &
       FILE_CARTESC_close
    use scale_file_grads, only: &
       FILE_GrADS_open,      &
       FILE_GrADS_get_shape, &
       FILE_GrADS_varid,     &
       FILE_GrADS_isOneD,    &
       FILE_GrADS_read
    use scale_time, only: &
       TIME_NOWDATE
    use scale_atmos_grid_cartesC, only: &
       CX => ATMOS_GRID_CARTESC_CX, &
       CY => ATMOS_GRID_CARTESC_CY
    use scale_atmos_grid_cartesC_real, only: &
       LAT => ATMOS_GRID_CARTESC_REAL_LAT, &
       LON => ATMOS_GRID_CARTESC_REAL_LON, &
       CZ  => ATMOS_GRID_CARTESC_REAL_CZ
    use scale_interp, only: &
       INTERP_factor3d_linear_xy,     &
       INTERP_factor3d_weight, &
       INTERP_interp3d
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy

    integer,                intent(in) :: CNVUSER_NSTEPS
    character(len=H_LONG),  intent(in) :: CNVUSER_GrADS_FILENAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_VARNAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_LATNAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_LONNAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_LEVNAME
    character(len=H_SHORT), intent(in) :: CNVUSER_GrADS_HEIGHT_PLEV
    character(len=H_SHORT), intent(in) :: CNVUSER_INTERP_TYPE
    integer,                intent(inout) :: CNVUSER_INTERP_LEVEL
    character(len=H_LONG),  intent(in) :: CNVUSER_OUT_BASENAME
    character(len=H_MID),   intent(in) :: CNVUSER_OUT_TITLE
    character(len=H_SHORT), intent(in) :: CNVUSER_OUT_VARNAME
    character(len=H_MID),   intent(in) :: CNVUSER_OUT_VARDESC
    character(len=H_SHORT), intent(in) :: CNVUSER_OUT_VARUNIT
    character(len=H_SHORT), intent(in) :: CNVUSER_OUT_DTYPE
    real(DP),               intent(in) :: CNVUSER_OUT_DT

    integer  :: i, j, k, step
    integer  :: file_id, var_id
    integer  :: data_shape(3)
    integer  :: fid, vid
    real(DP) :: timeofs
    real(RP) :: var_3D(KA,IA,JA,1)

    ! interpolation
    integer :: nLEV, nLON, nLAT
    logical :: zonal, pole

    integer,  allocatable :: idx_i     (:,:,:)
    integer,  allocatable :: idx_j     (:,:,:)
    real(RP), allocatable :: hfact     (:,:,:)
    integer,  allocatable :: idx_k     (:,:,:,:,:)
    real(RP), allocatable :: vfact     (:,:,:,:)
    real(RP), allocatable :: data_org  (:,:,:)
    real(RP), allocatable :: X_org     (:,:)
    real(RP), allocatable :: Y_org     (:,:)
    real(RP), allocatable :: LAT_org   (:,:)
    real(RP), allocatable :: LON_org   (:,:)
    real(RP), allocatable :: LEV_org   (:,:,:)
    real(RP), allocatable :: LAT_1d    (:)
    real(RP), allocatable :: LON_1d    (:)
    real(RP), allocatable :: LEV_1d    (:)

    call FILE_GrADS_open( CNVUSER_GrADS_FILENAME, & ! [IN]
                          file_id                 ) ! [OUT]

    call FILE_GrADS_get_shape( file_id,               & ! [IN]
                               CNVUSER_GrADS_VARNAME, & ! [IN]
                               data_shape(:)          ) ! [OUT]
    nLEV = data_shape(1)
    nLON = data_shape(2)
    nLAT = data_shape(3)

    ! interporation
    select case ( trim(CNVUSER_INTERP_TYPE) )
    case ( 'LINEAR' )
       CNVUSER_INTERP_LEVEL = 4
    case ( 'DIST-WEIGHT' )
       ! do nothing
    end select

    allocate( idx_i(     IA,JA,CNVUSER_INTERP_LEVEL) )
    allocate( idx_j(     IA,JA,CNVUSER_INTERP_LEVEL) )
    allocate( hfact(     IA,JA,CNVUSER_INTERP_LEVEL) )
    allocate( idx_k(KA,2,IA,JA,CNVUSER_INTERP_LEVEL) )
    allocate( vfact(KA,  IA,JA,CNVUSER_INTERP_LEVEL) )

    allocate( data_org( nLEV, nLON, nLAT ) )
    allocate( X_org   (       nLON, nLAT ) )
    allocate( Y_org   (       nLON, nLAT ) )
    allocate( LAT_org (       nLON, nLAT ) )
    allocate( LON_org (       nLON, nLAT ) )
    allocate( LEV_org ( nLEV, nLON, nLAT ) )
    allocate( LAT_1d  (             nLAT ) )
    allocate( LON_1d  (       nLON       ) )
    allocate( LEV_1d  ( nLEV             ) )

    ! lat
    call FILE_GrADS_varid( file_id,               & ! [IN]
                           CNVUSER_GrADS_LATNAME, & ! [IN]
                           var_id                 ) ! [OUT]

    if ( FILE_GrADS_isOneD( file_id, var_id ) ) then
       call FILE_GrADS_read( file_id,  & ! [IN]
                             var_id,   & ! [IN]
                             LAT_1d(:) ) ! [OUT]

       !$omp parallel do
       do j = 1, nLAT
       do i = 1, nLON
          LAT_org(i,j) = LAT_1d(j) * D2R
       end do
       end do
    else
       call FILE_GrADS_read( file_id,     & ! [IN]
                             var_id,      & ! [IN]
                             LAT_org(:,:) ) ! [OUT]

       !$omp parallel do
       do j = 1, nLAT
       do i = 1, nLON
          LAT_org(i,j) = LAT_org(i,j) * D2R
       end do
       end do
    end if

    ! lon
    call FILE_GrADS_varid( file_id,               & ! [IN]
                           CNVUSER_GrADS_LONNAME, & ! [IN]
                           var_id                 ) ! [OUT]

    if ( FILE_GrADS_isOneD( file_id, var_id ) ) then
       call FILE_GrADS_read( file_id,  & ! [IN]
                             var_id,   & ! [IN]
                             LON_1d(:) ) ! [OUT]

       !$omp parallel do
       do j = 1, nLAT
       do i = 1, nLON
          LON_org(i,j) = LON_1d(i) * D2R
       end do
       end do
    else
       call FILE_GrADS_read( file_id,     & ! [IN]
                             var_id,      & ! [IN]
                             LON_org(:,:) ) ! [OUT]

       !$omp parallel do
       do j = 1, nLAT
       do i = 1, nLON
          LON_org(i,j) = LON_org(i,j) * D2R
       end do
       end do
    end if

    ! lev
    select case ( trim(CNVUSER_GrADS_LEVNAME) )
    case ( 'zlev' )
       call FILE_GrADS_varid( file_id,               & ! [IN]
                              CNVUSER_GrADS_LEVNAME, & ! [IN]
                              var_id                 ) ! [OUT]

       call FILE_GrADS_read( file_id,  & ! [IN]
                             var_id,   & ! [IN]
                             LEV_1d(:) ) ! [OUT]

       !$omp parallel do collapse(2)
       do j = 1, nLAT
       do i = 1, nLON
       do k = 1, nLEV
          LEV_org(k,i,j) = LEV_1d(k)
       end do
       end do
       end do
    case ( 'plev' )
       call FILE_GrADS_varid( file_id,                   & ! [IN]
                              CNVUSER_GrADS_HEIGHT_PLEV, & ! [IN]
                              var_id                     ) ! [OUT]

       call FILE_GrADS_read( file_id,       & ! [IN]
                             var_id,        & ! [IN]
                             LEV_org(:,:,:) ) ! [OUT]
    case default
       LOG_ERROR("CNVUSER_setup",*) 'Invalid proparty in CNVUSER_GrADS_LEVNAME: ', trim(CNVUSER_GrADS_LEVNAME)
       call PRC_abort
    end select

    ! prepare to read target data
    call FILE_GrADS_varid( file_id,               & ! [IN]
                           CNVUSER_GrADS_VARNAME, & ! [IN]
                           var_id                 ) ! [OUT]

    select case ( CNVUSER_INTERP_TYPE )
    case ( 'LINEAR' )
       call MAPPROJECTION_lonlat2xy( nLON, 1, nLON, & ! [IN]
                                     nLAT, 1, nLAT, & ! [IN]
                                     LON_org(:,:),  & ! [IN]
                                     LAT_org(:,:),  & ! [IN]
                                     X_org  (:,:),  & ! [OUT]
                                     Y_org  (:,:)   ) ! [OUT]

       zonal = ( maxval(LON_org(:,:)) - minval(LAT_org(:,:)) ) > 2.0_RP * PI * 0.9_RP
       pole  = ( maxval(LAT_org(:,:)) > PI * 0.5_RP * 0.9_RP ) .or. ( minval(LAT_org(:,:)) < - PI * 0.5_RP * 0.9_RP )

       call INTERP_factor3d_linear_xy( nLEV, 1, nLEV,    & ! [IN]
                                       nLON, nLAT,       & ! [IN]
                                       KA, 1, KA,        & ! [IN]
                                       IA, JA,           & ! [IN]
                                       X_org(:,:),       & ! [IN]
                                       Y_org(:,:),       & ! [IN]
                                       LEV_org(:,:,:),   & ! [IN]
                                       CX(:),            & ! [IN]
                                       CY(:),            & ! [IN]
                                       CZ(:,:,:),        & ! [IN]
                                       idx_i(:,:,:),     & ! [OUT]
                                       idx_j(:,:,:),     & ! [OUT]
                                       hfact(:,:,:),     & ! [OUT]
                                       idx_k(:,:,:,:,:), & ! [OUT]
                                       vfact(:,:,:,:),   & ! [OUT]
                                       zonal = zonal,    & ! [IN]
                                       pole  = pole      ) ! [IN]
    case ( 'DIST-WEIGHT' )
       call INTERP_factor3d_weight( CNVUSER_INTERP_LEVEL, & ! [IN]
                                    nLEV, 1, nLEV,        & ! [IN]
                                    nLON, nLAT,           & ! [IN]
                                    KA, 1, KA,            & ! [IN]
                                    IA, JA,               & ! [IN]
                                    LON_org(:,:),         & ! [IN]
                                    LAT_org(:,:),         & ! [IN]
                                    LEV_org(:,:,:),       & ! [IN]
                                    LON(:,:),             & ! [IN]
                                    LAT(:,:),             & ! [IN]
                                    CZ(:,:,:),            & ! [IN]
                                    idx_i(:,:,:),         & ! [OUT]
                                    idx_j(:,:,:),         & ! [OUT]
                                    hfact(:,:,:),         & ! [OUT]
                                    idx_k(:,:,:,:,:),     & ! [OUT]
                                    vfact(:,:,:,:)        ) ! [OUT]
    end select

    call FILE_CARTESC_create( CNVUSER_OUT_BASENAME,             & ! [IN]
                              CNVUSER_OUT_TITLE,                & ! [IN]
                              CNVUSER_OUT_DTYPE,                & ! [IN]
                              fid,                              & ! [OUT]
                              date = TIME_NOWDATE,              & ! [IN]
                              haszcoord = .true.,               & ! [IN]
                              aggregate = CNVUSER_OUT_AGGREGATE ) ! [IN]

    call FILE_CARTESC_def_var( fid,                       & ! [IN]
                               CNVUSER_OUT_VARNAME,       & ! [IN]
                               CNVUSER_OUT_VARDESC,       & ! [IN]
                               CNVUSER_OUT_VARUNIT,       & ! [IN]
                               'ZXYT',                    & ! [IN]
                               CNVUSER_OUT_DTYPE,         & ! [IN]
                               vid,                       & ! [OUT]
                               timeintv = CNVUSER_OUT_DT, & ! [IN]
                               nsteps = CNVUSER_NSTEPS    ) ! [IN]

    call FILE_CARTESC_enddef(fid)

    call FILE_GrADS_open( CNVUSER_GrADS_FILENAME, & ! [IN]
                          file_id                 ) ! [OUT]

    call FILE_GrADS_varid( file_id,               & ! [IN]
                           CNVUSER_GrADS_VARNAME, & ! [IN]
                           var_id                 ) ! [OUT]

    do step = 1, CNVUSER_NSTEPS
       LOG_PROGRESS(*) 'step = ', step
       timeofs = ( step - 1 ) * CNVUSER_OUT_DT

       call FILE_GrADS_read( file_id, var_id, & ! [IN]
                             data_org(:,:,:), & ! [OUT]
                             step = step      ) ! [IN]

       call INTERP_interp3d( CNVUSER_INTERP_LEVEL, & ! [IN]
                             nLEV, 1, nLEV,        & ! [IN]
                             nLON, nLAT,           & ! [IN]
                             KA, 1, KA,            & ! [IN]
                             IA, JA,               & ! [IN]
                             idx_i(:,:,:),         & ! [IN]
                             idx_j(:,:,:),         & ! [IN]
                             hfact(:,:,:),         & ! [IN]
                             idx_k(:,:,:,:,:),     & ! [IN]
                             vfact(:,:,:,:),       & ! [IN]
                             LEV_org(:,:,:),       & ! [IN]
                             CZ(:,:,:),            & ! [IN]
                             data_org(:,:,:),      & ! [IN]
                             var_3D(:,:,:,1)       ) ! [OUT]

       call FILE_CARTESC_write_var( fid, vid,            & ! [IN]
                                    var_3D(:,:,:,:),     & ! [IN]
                                    CNVUSER_OUT_VARNAME, & ! [IN]
                                    'ZXYT',              & ! [IN]
                                    CNVUSER_OUT_DT,      & ! [IN]
                                    timetarg = 1,        & ! [IN]
                                    timeofs  = timeofs   ) ! [IN]
    end do

    call FILE_CARTESC_close( fid )

    deallocate( idx_i )
    deallocate( idx_j )
    deallocate( hfact )
    deallocate( idx_k )
    deallocate( vfact )

    deallocate( data_org )
    deallocate( X_org    )
    deallocate( Y_org    )
    deallocate( LAT_org  )
    deallocate( LON_org  )
    deallocate( LEV_org  )
    deallocate( LAT_1d   )
    deallocate( LON_1d   )
    deallocate( LEV_1d   )
  end subroutine CNVUSER_execute_GrADS_3D

  subroutine CNVUSER_execute_TILE_GrADS( &
    CNVUSER_NSTEPS, &
    CNVUSER_OUT_BASENAME, &
    CNVUSER_OUT_TITLE, &
    CNVUSER_OUT_VARNAME, &
    CNVUSER_OUT_VARDESC, &
    CNVUSER_OUT_VARUNIT, &
    CNVUSER_OUT_DTYPE, &
    CNVUSER_OUT_DT &
    )

    use scale_file_cartesc, only: &
       FILE_CARTESC_create, &
       FILE_CARTESC_def_var, &
       FILE_CARTESC_enddef, &
       FILE_CARTESC_close
    use scale_time, only: &
       TIME_NOWDATE
    use mod_cnv2d, only: &
       CNV2D_exec

    integer,                intent(in) :: CNVUSER_NSTEPS
    character(len=H_LONG),  intent(in) :: CNVUSER_OUT_BASENAME
    character(len=H_MID),   intent(in) :: CNVUSER_OUT_TITLE
    character(len=H_SHORT), intent(in) :: CNVUSER_OUT_VARNAME
    character(len=H_MID),   intent(in) :: CNVUSER_OUT_VARDESC
    character(len=H_SHORT), intent(in) :: CNVUSER_OUT_VARUNIT
    character(len=H_SHORT), intent(in) :: CNVUSER_OUT_DTYPE
    real(DP),               intent(in) :: CNVUSER_OUT_DT

    real(RP) :: var(IA, JA)

    integer :: fid, vid
    integer :: step

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

       call CNVUSER_write( CNVUSER_OUT_VARNAME, fid, vid, var(:,:), CNVUSER_OUT_DT, step )

    end do

    call FILE_CARTESC_close( fid )
  end subroutine CNVUSER_execute_TILE_GrADS

  subroutine CNVUSER_write( &
       CNVUSER_OUT_VARNAME, &
       fid, vid,  &
       VAR,       &
       timeintv,  &
       istep      )
    use scale_file_cartesC, only: &
       FILE_CARTESC_write_var

    character(len=H_SHORT), intent(in) :: CNVUSER_OUT_VARNAME
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
