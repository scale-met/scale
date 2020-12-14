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
module mod_gm_cnv2d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CNV2D_tile_init
  public :: CNV2D_GRADS_init
  public :: CNV2D_convert

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
  integer,  private, parameter   :: TILE_nlim = 1000

  integer,  private              :: CNV2D_ftype
  integer,  private, parameter   :: I_TILE  = 1
  integer,  private, parameter   :: I_GrADS = 2

  ! TILE data
  character(len=H_SHORT), private :: CNV2D_tile_dtype
  character(len=H_LONG),  private :: CNV2D_tile_dir
  real(RP), private               :: CNV2D_TILE_dlat
  real(RP), private               :: CNV2D_TILE_dlon
  integer,  private               :: CNV2D_TILE_nlat
  integer,  private               :: CNV2D_TILE_nlon

  integer,  private               :: TILEDATA_GLOBAL_IA
  integer,  private               :: TILEDATA_TILE_nmax
  character(len=H_LONG),  private :: TILEDATA_TILE_fname(TILE_nlim)
  logical,  private               :: TILEDATA_TILE_hit  (TILE_nlim)
  integer,  private               :: TILEDATA_TILE_JS   (TILE_nlim)
  integer,  private               :: TILEDATA_TILE_JE   (TILE_nlim)
  integer,  private               :: TILEDATA_TILE_IS   (TILE_nlim)
  integer,  private               :: TILEDATA_TILE_IE   (TILE_nlim)
  integer,  private               :: TILEDATA_DOMAIN_IS
  integer,  private               :: TILEDATA_DOMAIN_IE
  integer,  private               :: TILEDATA_DOMAIN_JS
  integer,  private               :: TILEDATA_DOMAIN_JE
  logical,  private               :: TILEDATA_zonal
  logical,  private               :: TILEDATA_pole

  ! GrADS data
  integer,  private               :: CNV2D_GRADS_fid
  integer,  private               :: CNV2D_GRADS_vid
  integer,  private               :: CNV2D_GRADS_nlat
  integer,  private               :: CNV2D_GRADS_nlon

  real(RP), private, allocatable :: DATA_org(:,:)
  real(RP), private, allocatable :: LAT_org (:,:)
  real(RP), private, allocatable :: LON_org (:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine CNV2D_tile_init( &
       dtype,      &
       dlat,       &
       dlon,       &
       dir,        &
       catalogue,  &
       interp_type )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       D2R   => CONST_D2R
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info,   &
       FILE_TILEDATA_get_latlon, &
       FILE_TILEDATA_get_data
    implicit none

    character(len=*), intent(in) :: dtype
    real(RP),         intent(in) :: dlat
    real(RP),         intent(in) :: dlon
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: catalogue
    character(len=*), intent(in) :: interp_type

    real(RP)              :: DOMAIN_LATS, DOMAIN_LATE
    real(RP)              :: DOMAIN_LONS, DOMAIN_LONE
    character(len=H_LONG) :: fname

    real(RP), allocatable :: LAT_1d(:)
    real(RP), allocatable :: LON_1d(:)

    integer :: i, j
    !---------------------------------------------------------------------------

    CNV2D_ftype      = I_TILE
    CNV2D_TILE_dtype = dtype
    CNV2D_TILE_dir   = dir
    CNV2D_TILE_dlat  = dlat * D2R
    CNV2D_TILE_dlon  = dlon * D2R

    ! read catalogue file
    DOMAIN_LATS =  -90.0_RP * D2R
    DOMAIN_LATE =   90.0_RP * D2R
    DOMAIN_LONS = -180.0_RP * D2R
    DOMAIN_LONE =  180.0_RP * D2R
    fname       = trim(dir)//'/'//trim(catalogue)

    call FILE_TILEDATA_get_info( TILE_nlim,                        & ! [IN]
                                 CNV2D_TILE_dlat, CNV2D_TILE_dlon, & ! [IN]
                                 DOMAIN_LATS, DOMAIN_LATE,         & ! [IN]
                                 DOMAIN_LONS, DOMAIN_LONE,         & ! [IN]
                                 fname,                            & ! [IN]
                                 TILEDATA_GLOBAL_IA,               & ! [OUT]
                                 TILEDATA_TILE_nmax,               & ! [OUT]
                                 TILEDATA_TILE_fname(:),           & ! [OUT]
                                 TILEDATA_TILE_hit  (:),           & ! [OUT]
                                 TILEDATA_TILE_JS   (:),           & ! [OUT]
                                 TILEDATA_TILE_JE   (:),           & ! [OUT]
                                 TILEDATA_TILE_IS   (:),           & ! [OUT]
                                 TILEDATA_TILE_IE   (:),           & ! [OUT]
                                 CNV2D_TILE_nlat, CNV2D_TILE_nlon, & ! [OUT]
                                 TILEDATA_DOMAIN_JS,               & ! [OUT]
                                 TILEDATA_DOMAIN_JE,               & ! [OUT]
                                 TILEDATA_DOMAIN_IS,               & ! [OUT]
                                 TILEDATA_DOMAIN_IE,               & ! [OUT]
                                 TILEDATA_zonal, TILEDATA_pole     ) ! [OUT]

    allocate( LAT_org (CNV2D_TILE_nlon,CNV2D_TILE_nlat) )
    allocate( LON_org (CNV2D_TILE_nlon,CNV2D_TILE_nlat) )
    allocate( DATA_org(CNV2D_TILE_nlon,CNV2D_TILE_nlat) )

    allocate( LAT_1d(CNV2D_TILE_nlat) )
    allocate( LON_1d(CNV2D_TILE_nlon) )

    ! read lat, lon
    call FILE_TILEDATA_get_latlon( CNV2D_TILE_nlat, CNV2D_TILE_nlon,       & ! [IN]
                                   TILEDATA_DOMAIN_JS, TILEDATA_DOMAIN_IS, & ! [IN]
                                   CNV2D_TILE_dlat, CNV2D_TILE_dlon,       & ! [IN]
                                   LAT_1d(:), LON_1d(:)                    ) ! [OUT]

    do j = 1, CNV2D_TILE_nlat
    do i = 1, CNV2D_TILE_nlon
       LAT_org(i,j) = LAT_1d(j)
       LON_org(i,j) = LON_1d(i)
    enddo
    enddo

    deallocate( LAT_1d )
    deallocate( LON_1d )

    call CNV2D_init( interp_type )

    return
  end subroutine CNV2D_tile_init

  !-----------------------------------------------------------------------------
  subroutine CNV2D_GRADS_init( &
       FILE_NAME,    &
       VAR_NAME,     &
       LAT_NAME,     &
       LON_NAME,     &
       interp_type,  &
       interp_level, &
       search_limit, &
       POSTFIX       )
    use scale_file_grads, only: &
       FILE_GrADS_open,      &
       FILE_GrADS_get_shape, &
       FILE_GrADS_varid,     &
       FILE_GrADS_isOneD,    &
       FILE_GrADS_read
    use scale_const, only: &
       D2R => CONST_D2R
    implicit none

    character(len=*), intent(in) :: FILE_NAME
    character(len=*), intent(in) :: VAR_NAME
    character(len=*), intent(in) :: LAT_NAME
    character(len=*), intent(in) :: LON_NAME
    character(len=*), intent(in) :: INTERP_TYPE
    integer,          intent(in), optional :: interp_level
    real(RP),         intent(in), optional :: search_limit
    character(len=*), intent(in), optional :: POSTFIX

    integer :: file_id, var_id
    integer :: shape(2)

    real(RP), allocatable :: LAT_1d(:)
    real(RP), allocatable :: LON_1d(:)

    integer :: i, j
    !---------------------------------------------------------------------------

    call FILE_GrADS_open( FILE_NAME, & ! [IN]
                          file_id    ) ! [OUT]

    call FILE_GrADS_varid( file_id, VAR_NAME, & ! [IN]
                           var_id             ) ! [OUT]

    call FILE_GrADS_get_shape( file_id, VAR_NAME, & ! [IN]
                               shape(:)           ) ! [OUT]

    CNV2D_ftype      = I_GRADS
    CNV2D_GRADS_fid  = file_id
    CNV2D_GRADS_vid  = var_id
    CNV2D_GRADS_nlon = shape(1)
    CNV2D_GRADS_nlat = shape(2)

    allocate( LAT_org (CNV2D_GRADS_nlon,CNV2D_GRADS_nlat) )
    allocate( LON_org (CNV2D_GRADS_nlon,CNV2D_GRADS_nlat) )
    allocate( DATA_org(CNV2D_GRADS_nlon,CNV2D_GRADS_nlat) )

    allocate( LAT_1d(CNV2D_GRADS_nlat) )
    allocate( LON_1d(CNV2D_GRADS_nlon) )

    ! read lat
    call FILE_GrADS_varid( file_id, LAT_NAME, & ! [IN]
                           var_id             ) ! [OUT]

    if ( FILE_GrADS_isOneD( file_id, var_id ) ) then
       call FILE_GrADS_read( file_id, var_id, & ! [IN]
                             LAT_1d(:)        ) ! [OUT]

       do j = 1, CNV2D_GRADS_nlat
       do i = 1, CNV2D_GRADS_nlon
          LAT_org(i,j) = LAT_1d(j) * D2R
       enddo
       enddo
    else
       call FILE_GrADS_read( file_id, var_id,  & ! [IN]
                             LAT_org(:,:),     & ! [OUT]
                             postfix = POSTFIX ) ! [IN]

       do j = 1, CNV2D_GRADS_nlat
       do i = 1, CNV2D_GRADS_nlon
          LAT_org(i,j) = LAT_org(i,j) * D2R
       enddo
       enddo
    endif

    ! read lon
    call FILE_GrADS_varid( file_id, LON_NAME, & ! [IN]
                           var_id             ) ! [OUT]

    if ( FILE_GrADS_isOneD( file_id, var_id ) ) then
       call FILE_GrADS_read( file_id, var_id, & ! [IN]
                             LON_1d(:)        ) ! [OUT]

       do j = 1, CNV2D_GRADS_nlat
       do i = 1, CNV2D_GRADS_nlon
          LON_org(i,j) = LON_1d(i) * D2R
       enddo
       enddo
    else
       call FILE_GrADS_read( file_id, var_id,  & ! [IN]
                             LON_org(:,:),     & ! [OUT]
                             postfix = POSTFIX ) ! [IN]

       do j = 1, CNV2D_GRADS_nlat
       do i = 1, CNV2D_GRADS_nlon
          LON_org(i,j) = LON_org(i,j) * D2R
       enddo
       enddo
    endif

    deallocate( LAT_1d )
    deallocate( LON_1d )

    call CNV2D_init( interp_type )

    return
  end subroutine CNV2D_GRADS_init

  !-----------------------------------------------------------------------------
  subroutine CNV2D_init( &
       interp_type )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in) :: INTERP_TYPE
    !---------------------------------------------------------------------------

    select case( INTERP_TYPE )
    case('LINEAR')
    case('DIST-WEIGHT')
    case default
       LOG_ERROR('CNV2D_init',*) 'INTERP_TYPE is invalid: ', trim(INTERP_TYPE)
       call PRC_abort
    end select

    return
  end subroutine CNV2D_init

  !-----------------------------------------------------------------------------
  subroutine CNV2D_convert( &
       var,       &
       var_pl,    &
       step,      &
       min_value, &
       yrevers    )
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_data
    use scale_file_grads, only: &
       FILE_GrADS_read
    implicit none

    real(RP), intent(out) :: var   (ADM_gall   ,ADM_lall   )
    real(RP), intent(out) :: var_pl(ADM_gall_pl,ADM_lall_pl)
    integer,  intent(in), optional :: step
    real(RP), intent(in), optional :: min_value
    logical,  intent(in), optional :: yrevers
    !---------------------------------------------------------------------------

    select case(CNV2D_ftype)
    case(I_TILE)

       call FILE_TILEDATA_get_data( CNV2D_TILE_nlat, CNV2D_TILE_nlon, & ! [IN]
                                    CNV2D_TILE_DIR,                   & ! [IN]
                                    TILEDATA_GLOBAL_IA,               & ! [IN]
                                    TILEDATA_TILE_nmax,               & ! [IN]
                                    CNV2D_TILE_dlat, CNV2D_TILE_dlon, & ! [IN]
                                    TILEDATA_TILE_fname(:),           & ! [IN]
                                    TILEDATA_TILE_hit  (:),           & ! [IN]
                                    TILEDATA_TILE_JS   (:),           & ! [IN]
                                    TILEDATA_TILE_JE   (:),           & ! [IN]
                                    TILEDATA_TILE_IS   (:),           & ! [IN]
                                    TILEDATA_TILE_IE   (:),           & ! [IN]
                                    TILEDATA_DOMAIN_JS,               & ! [IN]
                                    TILEDATA_DOMAIN_JE,               & ! [IN]
                                    TILEDATA_DOMAIN_IS,               & ! [IN]
                                    TILEDATA_DOMAIN_IE,               & ! [IN]
                                    CNV2D_TILE_DTYPE,                 & ! [IN]
                                    DATA_org(:,:),                    & ! [OUT]
                                    step      = step,                 & ! [IN]
                                    min_value = min_value,            & ! [IN]
                                    yrevers   = yrevers               ) ! [IN]

    case(I_GrADS)

       call FILE_GrADS_read( CNV2D_GRADS_fid, & ! [IN]
                             CNV2D_GRADS_vid, & ! [IN]
                             DATA_org(:,:),   & ! [OUT]
                             step = step      ) ! [IN]

    end select

    ! interpolation

    return
  end subroutine CNV2D_convert

end module mod_gm_cnv2d
