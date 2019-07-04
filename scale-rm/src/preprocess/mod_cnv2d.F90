!-------------------------------------------------------------------------------
!> module Convert 2D data
!!
!! @par Description
!!          subroutines for preparing 2D data (convert from external file)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cnv2d
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
  public :: CNV2D_setup
  public :: CNV2D_tile_init
  public :: CNV2D_grads_init
  public :: CNV2D_exec

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNV2D_init

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: I_TILE  = 1
  integer, parameter :: I_GrADS = 2
  integer :: CNV2D_ftype

  real(RP) :: DOMAIN_LATS, DOMAIN_LATE
  real(RP) :: DOMAIN_LONS, DOMAIN_LONE
  integer  :: GLOBAL_IA

  integer               :: nLON, nLAT
  real(RP), allocatable :: data_org(:,:)
  real(RP), allocatable :: LAT_org (:,:)
  real(RP), allocatable :: LON_org (:,:)
  real(RP), allocatable :: LAT_1d  (:)
  real(RP), allocatable :: LON_1d  (:)

  ! interpolation
  integer               :: itp_lev
  integer,  allocatable :: idx_i(:,:,:)
  integer,  allocatable :: idx_j(:,:,:)
  real(RP), allocatable :: hfact(:,:,:)
  logical :: zonal, pole

  logical :: first = .true.

  ! TILE data
  character(len=H_SHORT) :: CNV2D_tile_dtype
  character(len=H_LONG)  :: CNV2D_tile_dir
  integer, parameter     :: TILE_nlim = 100
  integer                :: TILE_nmax
  character(len=H_LONG)  :: TILE_fname(TILE_nlim)
  logical                :: TILE_hit  (TILE_nlim)
  integer                :: TILE_JS   (TILE_nlim)
  integer                :: TILE_JE   (TILE_nlim)
  integer                :: TILE_IS   (TILE_nlim)
  integer                :: TILE_IE   (TILE_nlim)
  real(RP)               :: TILE_DLAT, TILE_DLON
  integer                :: dom_is, dom_ie, dom_js, dom_je


  ! GrADS data
  integer :: CNV2D_grads_fid
  integer :: CNV2D_grads_vid

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNV2D_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       D2R => CONST_D2R
    use scale_atmos_grid_cartesC_real, only: &
       LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       LONUY => ATMOS_GRID_CARTESC_REAL_LONUY
    implicit none

    integer                :: CNV2D_NSTEPS      = 1
    character(len=H_SHORT) :: CNV2D_INTERP_TYPE = 'LINEAR'

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNV2D_setup",*) 'Setup'

    !--- read namelist
!!$    rewind(IO_FID_CONF)
!!$    read(IO_FID_CONF,nml=PARAM_CNV2D,iostat=ierr)
!!$    if( ierr < 0 ) then !--- missing
!!$       LOG_INFO("CNV2D_setup",*) 'Not found namelist. Default used.'
!!$    elseif( ierr > 0 ) then !--- fatal error
!!$       LOG_ERROR("CNV2D_setup",*) 'Not appropriate names in namelist PARAM_CNV2D. Check!'
!!$       call PRC_abort
!!$    endif
!!$    LOG_NML(PARAM_CNV2D)


    DOMAIN_LATS = minval( LATXV(:,:) )
    DOMAIN_LATE = maxval( LATXV(:,:) )
    DOMAIN_LONS = minval( LONUY(:,:) )
    DOMAIN_LONE = maxval( LONUY(:,:) )

    LOG_INFO("CNV2D_setup",*) 'Domain Information'
    LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
    LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R


    return
  end subroutine CNV2D_setup

  subroutine CNV2D_tile_init( &
       dtype, &
       dlat, dlon, &
       dir, &
       catalogue, &
       interp_type, &
       interp_level  )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF  => CONST_UNDEF, &
       D2R    => CONST_D2R
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info,   &
       FILE_TILEDATA_get_latlon, &
       FILE_TILEDATA_get_data
    implicit none
    character(len=*), intent(in) :: dtype
    real(RP),         intent(in) :: dlat, dlon
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: catalogue
    character(len=*), intent(in) :: interp_type

    integer, intent(in), optional :: interp_level

    real(RP) :: CNV3D_TILE_DLAT ! degree
    real(RP) :: CNV3D_TILE_DLON ! degree

    character(len=H_LONG) :: fname

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    !--- read namelist
    TILE_DLAT = dlon * D2R
    TILE_DLON = dlat * D2R

    ! catalogue file
    fname = trim(DIR)//'/'//trim(CATALOGUE)

    call FILE_TILEDATA_get_info( TILE_nlim,                                          & ! [IN]
                                 TILE_DLAT, TILE_DLON,                               & ! [IN]
                                 DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                 fname,                                              & ! [IN]
                                 GLOBAL_IA,                                          & ! [OUT]
                                 TILE_nmax,                                          & ! [OUT]
                                 TILE_fname(:), TILE_hit(:),                         & ! [OUT]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                 nLAT, nLON, dom_js, dom_je, dom_is, dom_ie,         & ! [OUT]
                                 zonal, pole                                         ) ! [OUT]

    allocate( LAT_1d(nLAT) )
    allocate( LON_1d(nLON) )

    call FILE_TILEDATA_get_latlon( nLAT, nLON,                     & ! [IN]
                                   GLOBAL_IA,                      & ! [IN]
                                   TILE_nmax,                      & ! [IN]
                                   TILE_DLAT, TILE_DLON,           & ! [IN]
                                   TILE_hit(:),                    & ! [IN]
                                   TILE_JS(:), TILE_JE(:),         & ! [IN]
                                   TILE_IS(:), TILE_IE(:),         & ! [IN]
                                   dom_js, dom_je, dom_is, dom_ie, & ! [IN]
                                   LAT_1d(:), LON_1d(:)            ) ! [OUT]

    CNV2D_ftype = I_TILE
    CNV2D_tile_dtype = DTYPE
    CNV2D_tile_dir   = dir

    call CNV2D_init( interp_type, &
                     interp_level = interp_level, &
                     ll_struct = .true.           )

    return
  end subroutine CNV2D_tile_init


  subroutine CNV2D_grads_init( &
       FILE_NAME,   &
       VAR_NAME,    &
       interp_type, &
       interp_level )
    use scale_file_grads, only: &
       FILE_GrADS_open, &
       FILE_GrADS_get_shape, &
       FILE_GrADS_varid
    implicit none
    character(len=*), intent(in) :: FILE_NAME
    character(len=*), intent(in) :: VAR_NAME
    character(len=*), intent(in) :: INTERP_TYPE

    integer, intent(in), optional :: interp_level

    integer :: file_id, var_id
    integer :: shape(2)

    call FILE_GrADS_open( FILE_NAME, & ! [IN]
                          file_id    ) ! [OUT]

    call FILE_GrADS_get_shape( file_id, VAR_NAME, &
                               shape(:)           )
    nLON = shape(1)
    nLAT = shape(2)

    call FILE_GrADS_varid( file_id, VAR_NAME, & ! [IN]
                           var_id             ) ! [OUT]

    CNV2D_ftype = I_GRADS
    CNV2D_GrADS_fid = file_id
    CNV2D_GrADS_vid = var_id

    call CNV2D_init( interp_type, &
                     interp_level = interp_level, &
                     ll_struct = .false.          )

    return
  end subroutine CNV2D_grads_init

  subroutine CNV2D_exec( &
       var, &
       step )
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_data
    use scale_file_grads, only: &
       FILE_GrADS_read
    use scale_interp, only: &
       INTERP_interp2d
    implicit none
    real(RP), intent(out) :: var(IA,JA)
    integer,  intent(in), optional :: step

    select case ( CNV2D_ftype )
    case ( I_TILE )
       call FILE_TILEDATA_get_data( nLAT, nLON,                                     & ! [IN]
                                    CNV2D_TILE_DIR,                                 & ! [IN]
                                    GLOBAL_IA,                                      & ! [IN]
                                    TILE_nmax,                                      & ! [IN]
                                    TILE_DLAT, TILE_DLON,                           & ! [IN]
                                    TILE_fname(:), TILE_hit(:),                     & ! [IN]
                                    TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:), & ! [IN]
                                    dom_js, dom_je, dom_is, dom_ie,                 & ! [IN]
                                    CNV2D_TILE_DTYPE,                               & ! [IN]
                                    data_org(:,:),                                  & ! [OUT]
                                    min_value = -9000.0_RP, yrevers = .true.,       & ! [IN]
                                    step = step                                     ) ! [IN]
    case ( I_GrADS )
       call FILE_GrADS_read( CNV2D_grads_fid, CNV2D_grads_vid, & ! [IN]
                             data_org(:,:),                    & ! [OUT]
                             step = step                       ) ! [IN]
    end select

    call INTERP_interp2d( itp_lev,                    & ! [IN]
                          nLON,nLAT,                  & ! [IN]
                          IA, JA,                     & ! [IN]
                          idx_i(:,:,:), idx_j(:,:,:), & ! [IN]
                          hfact(:,:,:),               & ! [IN]
                          data_org(:,:),              & ! [IN]
                          var(:,:)                    ) ! [OUT]

    return
  end subroutine CNV2D_exec


  ! private

  subroutine CNV2D_init( &
       interp_type, &
       interp_level, &
       search_limit, &
       ll_struct     )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       PI => CONST_PI
    use scale_interp, only: &
       INTERP_factor2d_linear_latlon, &
       INTERP_factor2d_linear_xy,     &
       INTERP_factor2d_weight
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
    use scale_atmos_grid_cartesC, only: &
       CX => ATMOS_GRID_CARTESC_CX, &
       CY => ATMOS_GRID_CARTESC_CY
    use scale_atmos_grid_cartesC_real, only: &
       LAT => ATMOS_GRID_CARTESC_REAL_LAT, &
       LON => ATMOS_GRID_CARTESC_REAL_LON
    implicit none

    character(len=*), intent(in) :: INTERP_TYPE
    integer,          intent(in), optional :: INTERP_LEVEL
    real(RP),         intent(in), optional :: search_limit
    logical,          intent(in), optional :: ll_struct

    real(RP), allocatable :: X_org(:,:), Y_org(:,:)
    !---------------------------------------------------------------------------

    if ( .not. first ) deallocate( idx_i, idx_j, hfact, data_org )

    select case ( INTERP_TYPE )
    case ( 'LINEAR' )
       itp_lev = 4
    case ( 'DIST-WEIGHT' )
       if ( present(INTERP_LEVEL) ) then
          itp_lev = INTERP_LEVEL
       else
          LOG_ERROR('CNV2D_init',*) 'INTERP_LEVEL is required for the DIST-WEIGHT interpolation'
          call PRC_abort
       end if
    case default
       LOG_ERROR('CNV2D_init',*) 'INTERP_TYPE is invalid: ', trim(INTERP_TYPE)
       call PRC_abort
    end select

    ! interporation
    allocate( idx_i(IA,JA,itp_lev) )
    allocate( idx_j(IA,JA,itp_lev) )
    allocate( hfact(IA,JA,itp_lev) )

    select case ( INTERP_TYPE )
    case ( 'LINEAR' )
       if ( ll_struct ) then
          call INTERP_factor2d_linear_latlon( nLON, nLAT,                 & ! [IN]
                                              IA, JA,                     & ! [IN]
                                              LON_1d(:), LAT_1d(:),       & ! [IN]
                                              LON(:,:), LAT(:,:),         & ! [IN]
                                              idx_i(:,:,:), idx_j(:,:,:), & ! [OUT]
                                              hfact(:,:,:)                ) ! [OUT]
       else
          allocate( X_org(nLON,nLAT), Y_org(nLON,nLAT) )
          call MAPPROJECTION_lonlat2xy( nLON, 1, nLON, nLAT, 1, nLAT, & ! [IN]
                                        LON_org(:,:), LAT_org(:,:),   & ! [IN]
                                        X_org(:,:), Y_org(:,:)        ) ! [OUT]

          zonal = ( maxval(LON_org) - minval(LAT_org) ) > 2.0_RP * PI * 0.9_RP
          pole = ( maxval(LAT_org) > PI * 0.5_RP * 0.9_RP ) .or. ( minval(LAT_org) < - PI * 0.5_RP * 0.9_RP )
          call INTERP_factor2d_linear_xy( nLON, nLAT,                 & ! [IN]
                                          IA, JA,                     & ! [IN]
                                          X_org(:,:), Y_org(:,:),     & ! [IN]
                                          CX(:), CY(:),               & ! [IN]
                                          idx_i(:,:,:), idx_j(:,:,:), & ! [OUT]
                                          hfact(:,:,:),               & ! [OUT]
                                          zonal = zonal, pole = pole  ) ! [IN]
          deallocate( X_org, Y_org )
       end if
    case ( 'DIST-WEIGHT' )
       call INTERP_factor2d_weight( itp_lev,       & ! [IN]
                                    nLON, nLAT,                   & ! [IN]
                                    IA, JA,                       & ! [IN]
                                    LON_org(:,:), LAT_org(:,:),   & ! [IN]
                                    LON(:,:), LAT(:,:),           & ! [IN]
                                    idx_i(:,:,:), idx_j(:,:,:),   & ! [OUT]
                                    hfact(:,:,:),                 & ! [OUT]
                                    latlon_structure = ll_struct, & ! [IN]
                                    lon_1d = LON_1d(:),           & ! [IN]
                                    lat_1d = LAT_1d(:),           & ! [IN]
                                    search_limit = search_limit   ) ! [IN]
    end select

    allocate( data_org(nLON,nLAT) )

    first = .false.

    return
  end subroutine CNV2D_init

  ! private

end module mod_cnv2d
