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
  private :: CNVUSER_prepare_TILE
  private :: CNVUSER_prepare_GrADS
  private :: CNVUSER_prepare_GrADS_3D
  private :: CNVUSER_execute_GrADS_3D
  private :: CNVUSER_execute_TILE_GrADS
  private :: CNVUSER_write

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical                  :: CNVUSER_OUT_AGGREGATE

  type, abstract :: t_param
    character(len=H_SHORT) :: INTERP_TYPE
    integer                :: INTERP_LEVEL

    character(len=H_LONG)  :: OUT_BASENAME     ! basename of the output file
    character(len=H_MID)   :: OUT_TITLE        ! title    of the output file
    character(len=H_SHORT) :: OUT_VARNAME      ! name  of the variable
    character(len=H_MID)   :: OUT_VARDESC      ! title of the variable
    character(len=H_SHORT) :: OUT_VARUNIT      ! units of the variable
    character(len=H_SHORT) :: OUT_DTYPE        ! REAL4 or REAL8
    real(DP)               :: OUT_DT           ! sec

    integer                :: NSTEPS           ! # of time steps
  end type t_param

  ! TILE data
  type, extends(t_param) :: t_tile
    character(len=H_SHORT) :: TILE_DTYPE       ! data type in the tiled data
    real(RP)               :: TILE_DLAT
    real(RP)               :: TILE_DLON
    character(len=H_LONG)  :: TILE_DIR
    character(len=H_LONG)  :: TILE_CATALOGUE
  end type t_tile

  ! GrADS data
  type, extends(t_param) :: t_grads
    character(len=H_LONG)  :: GrADS_FILENAME
    character(len=H_SHORT) :: GrADS_VARNAME
    character(len=H_SHORT) :: GrADS_LATNAME
    character(len=H_SHORT) :: GrADS_LONNAME
  end type t_grads

  ! GrADS-3D data
  type, extends(t_param) :: t_grads_3d
    character(len=H_LONG)  :: GrADS_FILENAME
    character(len=H_SHORT) :: GrADS_VARNAME
    character(len=H_SHORT) :: GrADS_LATNAME
    character(len=H_SHORT) :: GrADS_LONNAME
    character(len=H_SHORT) :: GrADS_LEVNAME
    character(len=H_SHORT) :: GrADS_HEIGHT_PLEV
  end type t_grads_3d

  type t_param_wrapper
    class(t_param), allocatable :: param
  end type t_param_wrapper

  type(t_param_wrapper), allocatable :: params(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup

  subroutine CNVUSER_setup
    use scale_file, only: &
       FILE_AGGREGATE

    character(len=H_SHORT) :: CNVUSER_FILE_TYPE
    character(len=H_SHORT) :: CNVUSER_INTERP_TYPE
    integer                :: CNVUSER_INTERP_LEVEL

    character(len=H_LONG)  :: CNVUSER_OUT_BASENAME
    character(len=H_MID)   :: CNVUSER_OUT_TITLE
    character(len=H_SHORT) :: CNVUSER_OUT_VARNAME
    character(len=H_MID)   :: CNVUSER_OUT_VARDESC
    character(len=H_SHORT) :: CNVUSER_OUT_VARUNIT
    character(len=H_SHORT) :: CNVUSER_OUT_DTYPE
    real(DP)               :: CNVUSER_OUT_DT

    integer                :: CNVUSER_NSTEPS

    ! TILE data
    character(len=H_SHORT) :: CNVUSER_TILE_DTYPE
    real(RP)               :: CNVUSER_TILE_DLAT
    real(RP)               :: CNVUSER_TILE_DLON
    character(len=H_LONG)  :: CNVUSER_TILE_DIR
    character(len=H_LONG)  :: CNVUSER_TILE_CATALOGUE

    ! GrADS data
    character(len=H_LONG)  :: CNVUSER_GrADS_FILENAME
    character(len=H_SHORT) :: CNVUSER_GrADS_VARNAME
    character(len=H_SHORT) :: CNVUSER_GrADS_LATNAME
    character(len=H_SHORT) :: CNVUSER_GrADS_LONNAME
    character(len=H_SHORT) :: CNVUSER_GrADS_LEVNAME
    character(len=H_SHORT) :: CNVUSER_GrADS_HEIGHT_PLEV

    namelist / PARAM_CNVUSER /    &
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

    integer  :: ierr, n_vars
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVUSER_setup",*) 'Setup'

    CNVUSER_OUT_AGGREGATE = FILE_AGGREGATE

    !--- read namelist
    rewind(IO_FID_CONF)
    n_vars = 0
    do
      CNVUSER_FILE_TYPE = ''
      read(IO_FID_CONF,nml=PARAM_CNVUSER,iostat=ierr)
      if( ierr < 0 ) then !--- missing
         exit
      elseif( ierr > 0 ) then !--- fatal error
         LOG_ERROR("CNVUSER_setup",*) 'Not appropriate names in namelist PARAM_CNVUSER. Check!'
         call PRC_abort
      endif
      if (CNVUSER_FILE_TYPE == "") cycle
      n_vars = n_vars + 1
    end do

    allocate(params(n_vars))
    rewind(IO_FID_CONF)
    n_vars = 0
    do
      ! Default value
      CNVUSER_FILE_TYPE           = ''
      CNVUSER_INTERP_TYPE         = 'LINEAR'
      CNVUSER_INTERP_LEVEL        = 5

      CNVUSER_OUT_BASENAME        = ''
      CNVUSER_OUT_TITLE           = 'SCALE-RM User Boundary'
      CNVUSER_OUT_VARNAME         = ''
      CNVUSER_OUT_VARDESC         = ''
      CNVUSER_OUT_VARUNIT         = ''
      CNVUSER_OUT_DTYPE           = 'DEFAULT'
      CNVUSER_OUT_DT              = -1_DP

      CNVUSER_NSTEPS              = 1


      CNVUSER_TILE_DTYPE          = 'real4'
      CNVUSER_TILE_DLAT           = -1
      CNVUSER_TILE_DLON           = -1
      CNVUSER_TILE_DIR            = ''
      CNVUSER_TILE_CATALOGUE      = ''


      CNVUSER_GrADS_FILENAME     = ''
      CNVUSER_GrADS_VARNAME      = ''
      CNVUSER_GrADS_LATNAME      = 'lat'
      CNVUSER_GrADS_LONNAME      = 'lon'
      CNVUSER_GrADS_LEVNAME      = ''
      CNVUSER_GrADS_HEIGHT_PLEV  = 'HGT'

      read(IO_FID_CONF,nml=PARAM_CNVUSER,iostat=ierr)
      if( ierr /= 0 ) exit
      LOG_NML(PARAM_CNVUSER)
      if (CNVUSER_FILE_TYPE == "") cycle

      n_vars = n_vars + 1

      select case (CNVUSER_FILE_TYPE)
      case ("TILE")
        allocate(t_tile::params(n_vars)%param)
      case ("GrADS")
        allocate(t_grads::params(n_vars)%param)
      case ("GrADS-3D")
        allocate(t_grads_3d::params(n_vars)%param)
      case default
         LOG_ERROR('CNVUSER_setup',*) 'CNVUSER_FILE_TYPE is invalid: ', CNVUSER_FILE_TYPE
         LOG_ERROR_CONT(*)            'It must be "TILE" or "GrADS".'
         call PRC_abort
      end select

      associate (param => params(n_vars)%param)
        param%INTERP_TYPE = CNVUSER_INTERP_TYPE
        param%INTERP_LEVEL = CNVUSER_INTERP_LEVEL
        param%OUT_BASENAME = CNVUSER_OUT_BASENAME
        param%OUT_TITLE = CNVUSER_OUT_TITLE
        param%OUT_VARNAME = CNVUSER_OUT_VARNAME
        param%OUT_VARDESC = CNVUSER_OUT_VARDESC
        param%OUT_VARUNIT = CNVUSER_OUT_VARUNIT
        param%OUT_DTYPE = CNVUSER_OUT_DTYPE
        param%OUT_DT = CNVUSER_OUT_DT
        param%NSTEPS = CNVUSER_NSTEPS
      end associate

      select type (param => params(n_vars)%param)
      type is (t_tile)
        param%TILE_DTYPE      = CNVUSER_TILE_DTYPE
        param%TILE_DLAT       = CNVUSER_TILE_DLAT
        param%TILE_DLON       = CNVUSER_TILE_DLON
        param%TILE_DIR        = CNVUSER_TILE_DIR
        param%TILE_CATALOGUE  = CNVUSER_TILE_CATALOGUE
      type is (t_grads)
        param%GrADS_FILENAME = CNVUSER_GrADS_FILENAME
        param%GrADS_VARNAME  = CNVUSER_GrADS_VARNAME
        param%GrADS_LATNAME  = CNVUSER_GrADS_LATNAME
        param%GrADS_LONNAME  = CNVUSER_GrADS_LONNAME
        if ( param%OUT_VARNAME == '' ) param%OUT_VARNAME = param%GrADS_VARNAME
      type is (t_grads_3d)
        param%GrADS_FILENAME    = CNVUSER_GrADS_FILENAME
        param%GrADS_VARNAME     = CNVUSER_GrADS_VARNAME
        param%GrADS_LATNAME     = CNVUSER_GrADS_LATNAME
        param%GrADS_LONNAME     = CNVUSER_GrADS_LONNAME
        param%GrADS_LEVNAME     = CNVUSER_GrADS_LEVNAME
        param%GrADS_HEIGHT_PLEV = CNVUSER_GrADS_HEIGHT_PLEV
        if ( param%OUT_VARNAME == '' ) param%OUT_VARNAME = param%GrADS_VARNAME
      end select
    end do

  end subroutine CNVUSER_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVUSER
    integer :: i
    !---------------------------------------------------------------------------

    do i = 1, size(params)
    associate (param => params(i)%param)
      if ( param%OUT_BASENAME == '' .or. param%OUT_VARNAME == '' ) then
         LOG_ERROR('CNVUSER',*) 'CNVUSER_OUT_BASENAME and CNVUSER_OUT_VARNAME are required'
         call PRC_abort
      end if

      select type (param)
      type is (t_tile)
        call CNVUSER_prepare_TILE(param)
        call CNVUSER_execute_TILE_GrADS(param)
      type is (t_grads)
        call CNVUSER_prepare_GrADS(param)
        call CNVUSER_execute_TILE_GrADS(param)
      type is (t_grads_3d)
        call CNVUSER_prepare_GrADS_3D(param)
        call CNVUSER_execute_GrADS_3D(param)
      end select

    end associate
    end do

  end subroutine CNVUSER

  ! private

  subroutine CNVUSER_prepare_TILE(tile)
    use mod_cnv2d, only: &
       CNV2D_tile_init
    type(t_tile), intent(in) :: tile

    call CNV2D_tile_init( tile%TILE_DTYPE,                   &
                          tile%TILE_DLAT, tile%TILE_DLON,    &
                          tile%TILE_DIR,                     &
                          tile%TILE_CATALOGUE,               &
                          tile%INTERP_TYPE,                  &
                          interp_level = tile%INTERP_LEVEL   )
  end subroutine CNVUSER_prepare_TILE

  subroutine CNVUSER_prepare_GrADS(grads)
    use mod_cnv2d, only: &
       CNV2D_grads_init
    type(t_grads), intent(inout) :: grads

    call CNV2D_grads_init( grads%GrADS_FILENAME, &
                           grads%GrADS_VARNAME,  &
                           grads%GrADS_LATNAME,  &
                           grads%GrADS_LONNAME,  &
                           grads%INTERP_TYPE,    &
                           grads%INTERP_LEVEL    )
  end subroutine CNVUSER_prepare_GrADS

  subroutine CNVUSER_prepare_GrADS_3D(grads_3d)
    type(t_grads_3d), intent(inout) :: grads_3d

  end subroutine CNVUSER_prepare_GrADS_3D

  subroutine CNVUSER_execute_GrADS_3D(grads_3d)
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

    type(t_grads_3d), intent(inout) :: grads_3d

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

    call FILE_GrADS_open( grads_3d%GrADS_FILENAME, & ! [IN]
                          file_id                  ) ! [OUT]

    call FILE_GrADS_get_shape( file_id,                & ! [IN]
                               grads_3d%GrADS_VARNAME, & ! [IN]
                               data_shape(:)           ) ! [OUT]
    nLEV = data_shape(1)
    nLON = data_shape(2)
    nLAT = data_shape(3)

    ! interporation
    select case ( trim(grads_3d%INTERP_TYPE) )
    case ( 'LINEAR' )
       grads_3d%INTERP_LEVEL = 4
    case ( 'DIST-WEIGHT' )
       ! do nothing
    end select

    allocate( idx_i(     IA,JA,grads_3d%INTERP_LEVEL) )
    allocate( idx_j(     IA,JA,grads_3d%INTERP_LEVEL) )
    allocate( hfact(     IA,JA,grads_3d%INTERP_LEVEL) )
    allocate( idx_k(KA,2,IA,JA,grads_3d%INTERP_LEVEL) )
    allocate( vfact(KA,  IA,JA,grads_3d%INTERP_LEVEL) )

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
    call FILE_GrADS_varid( file_id,                & ! [IN]
                           grads_3d%GrADS_LATNAME, & ! [IN]
                           var_id                  ) ! [OUT]

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
    call FILE_GrADS_varid( file_id,                & ! [IN]
                           grads_3d%GrADS_LONNAME, & ! [IN]
                           var_id                  ) ! [OUT]

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
    select case ( trim(grads_3d%GrADS_LEVNAME) )
    case ( 'zlev' )
       call FILE_GrADS_varid( file_id,                & ! [IN]
                              grads_3d%GrADS_LEVNAME, & ! [IN]
                              var_id                  ) ! [OUT]

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
       call FILE_GrADS_varid( file_id,                    & ! [IN]
                              grads_3d%GrADS_HEIGHT_PLEV, & ! [IN]
                              var_id                      ) ! [OUT]

       call FILE_GrADS_read( file_id,       & ! [IN]
                             var_id,        & ! [IN]
                             LEV_org(:,:,:) ) ! [OUT]
    case default
       LOG_ERROR("CNVUSER_execute_GrADS_3D",*) 'Invalid property in grads_3d%GrADS_LEVNAME: ', trim(grads_3d%GrADS_LEVNAME), ' for ', trim(grads_3d%GrADS_VARNAME)
       call PRC_abort
    end select

    ! prepare to read target data
    call FILE_GrADS_varid( file_id,                & ! [IN]
                           grads_3d%GrADS_VARNAME, & ! [IN]
                           var_id                  ) ! [OUT]

    select case ( grads_3d%INTERP_TYPE )
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
       call INTERP_factor3d_weight( grads_3d%INTERP_LEVEL, & ! [IN]
                                    nLEV, 1, nLEV,         & ! [IN]
                                    nLON, nLAT,            & ! [IN]
                                    KA, 1, KA,             & ! [IN]
                                    IA, JA,                & ! [IN]
                                    LON_org(:,:),          & ! [IN]
                                    LAT_org(:,:),          & ! [IN]
                                    LEV_org(:,:,:),        & ! [IN]
                                    LON(:,:),              & ! [IN]
                                    LAT(:,:),              & ! [IN]
                                    CZ(:,:,:),             & ! [IN]
                                    idx_i(:,:,:),          & ! [OUT]
                                    idx_j(:,:,:),          & ! [OUT]
                                    hfact(:,:,:),          & ! [OUT]
                                    idx_k(:,:,:,:,:),      & ! [OUT]
                                    vfact(:,:,:,:)         ) ! [OUT]
    end select

    call FILE_CARTESC_create( grads_3d%OUT_BASENAME,            & ! [IN]
                              grads_3d%OUT_TITLE,               & ! [IN]
                              grads_3d%OUT_DTYPE,               & ! [IN]
                              fid,                              & ! [OUT]
                              date = TIME_NOWDATE,              & ! [IN]
                              haszcoord = .true.,               & ! [IN]
                              aggregate = CNVUSER_OUT_AGGREGATE ) ! [IN]

    call FILE_CARTESC_def_var( fid,                        & ! [IN]
                               grads_3d%OUT_VARNAME,       & ! [IN]
                               grads_3d%OUT_VARDESC,       & ! [IN]
                               grads_3d%OUT_VARUNIT,       & ! [IN]
                               'ZXYT',                     & ! [IN]
                               grads_3d%OUT_DTYPE,         & ! [IN]
                               vid,                        & ! [OUT]
                               timeintv = grads_3d%OUT_DT, & ! [IN]
                               nsteps = grads_3d%NSTEPS    ) ! [IN]

    call FILE_CARTESC_enddef(fid)

    call FILE_GrADS_open( grads_3d%GrADS_FILENAME, & ! [IN]
                          file_id                  ) ! [OUT]

    call FILE_GrADS_varid( file_id,                & ! [IN]
                           grads_3d%GrADS_VARNAME, & ! [IN]
                           var_id                  ) ! [OUT]

    do step = 1, grads_3d%NSTEPS
       LOG_PROGRESS(*) 'step = ', step
       timeofs = ( step - 1 ) * grads_3d%OUT_DT

       call FILE_GrADS_read( file_id, var_id, & ! [IN]
                             data_org(:,:,:), & ! [OUT]
                             step = step      ) ! [IN]

       call INTERP_interp3d( grads_3d%INTERP_LEVEL, & ! [IN]
                             nLEV, 1, nLEV,         & ! [IN]
                             nLON, nLAT,            & ! [IN]
                             KA, 1, KA,             & ! [IN]
                             IA, JA,                & ! [IN]
                             idx_i(:,:,:),          & ! [IN]
                             idx_j(:,:,:),          & ! [IN]
                             hfact(:,:,:),          & ! [IN]
                             idx_k(:,:,:,:,:),      & ! [IN]
                             vfact(:,:,:,:),        & ! [IN]
                             LEV_org(:,:,:),        & ! [IN]
                             CZ(:,:,:),             & ! [IN]
                             data_org(:,:,:),       & ! [IN]
                             var_3D(:,:,:,1)        ) ! [OUT]

       call FILE_CARTESC_write_var( fid, vid,             & ! [IN]
                                    var_3D(:,:,:,: ),     & ! [IN]
                                    grads_3d%OUT_VARNAME, & ! [IN]
                                    'ZXYT',               & ! [IN]
                                    grads_3d%OUT_DT,      & ! [IN]
                                    timetarg = 1,         & ! [IN]
                                    timeofs  = timeofs    ) ! [IN]
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

  subroutine CNVUSER_execute_TILE_GrADS(param)

    use scale_file_cartesc, only: &
       FILE_CARTESC_create, &
       FILE_CARTESC_def_var, &
       FILE_CARTESC_enddef, &
       FILE_CARTESC_close
    use scale_time, only: &
       TIME_NOWDATE
    use mod_cnv2d, only: &
       CNV2D_exec

    class(t_param), intent(in) :: param

    real(RP) :: var(IA, JA)

    integer :: fid, vid
    integer :: step

    call FILE_CARTESC_create( param%OUT_BASENAME,               & ! [IN]
                              param%OUT_TITLE,                  & ! [IN]
                              param%OUT_DTYPE,                  & ! [IN]
                              fid,                              & ! [OUT]
                              date = TIME_NOWDATE,              & ! [IN]
                              haszcoord = .false.,              & ! [IN]
                              aggregate = CNVUSER_OUT_AGGREGATE ) ! [IN]

    call FILE_CARTESC_def_var( fid,                     & ! [IN]
                               param%OUT_VARNAME,       & ! [IN]
                               param%OUT_VARDESC,       & ! [IN]
                               param%OUT_VARUNIT,       & ! [IN]
                               'XYT',                   & ! [IN]
                               param%OUT_DTYPE,         & ! [IN]
                               vid,                     & ! [OUT]
                               timeintv = param%OUT_DT, & ! [IN]
                               nsteps = param%NSTEPS    ) ! [IN]

    call FILE_CARTESC_enddef(fid)

    do step = 1, param%NSTEPS
       LOG_PROGRESS(*) 'step = ', step

       call CNV2D_exec( var(:,:), step = step )

       call CNVUSER_write( param%OUT_VARNAME, fid, vid, var(:,:), param%OUT_DT, step )

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
