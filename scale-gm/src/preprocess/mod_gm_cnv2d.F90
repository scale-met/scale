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
  integer,  private, parameter    :: TILE_nlim = 1000

  integer,  private               :: CNV2D_ftype
  integer,  private, parameter    :: I_TILE  = 1
  integer,  private, parameter    :: I_GrADS = 2

  ! TILE data
  character(len=H_SHORT), private :: CNV2D_tile_dtype
  character(len=H_LONG),  private :: CNV2D_tile_dir
  real(RP), private               :: CNV2D_TILE_dlat
  real(RP), private               :: CNV2D_TILE_dlon

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

  real(RP), private, allocatable  :: DATA_org(:,:)
  real(RP), private, allocatable  :: LAT_org (:,:)
  real(RP), private, allocatable  :: LON_org (:,:)
  integer,  private               :: CNV2D_nlat
  integer,  private               :: CNV2D_nlon

  integer,  private               :: npoints
  integer,  private, allocatable  :: idx_i   (:,:,:)
  integer,  private, allocatable  :: idx_j   (:,:,:)
  real(RP), private, allocatable  :: hfact   (:,:,:)
  integer,  private, allocatable  :: idx_i_pl(:,:,:)
  integer,  private, allocatable  :: idx_j_pl(:,:,:)
  real(RP), private, allocatable  :: hfact_pl(:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine CNV2D_tile_init( &
       dtype,       &
       dlat,        &
       dlon,        &
       dir,         &
       catalogue,   &
       interp_type, &
       interp_level )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_have_pl
    use scale_const, only: &
       PI  => CONST_PI,  &
       D2R => CONST_D2R
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info,   &
       FILE_TILEDATA_get_latlon, &
       FILE_TILEDATA_get_data
    use mod_grd, only: &
       GRD_LAT,    &
       GRD_LAT_pl, &
       GRD_LON,    &
       GRD_LON_pl
    implicit none

    character(len=*), intent(in) :: dtype
    real(RP),         intent(in) :: dlat
    real(RP),         intent(in) :: dlon
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: catalogue
    character(len=*), intent(in) :: interp_type
    integer,          intent(in) :: interp_level

    real(RP)              :: DOMAIN_LATS, DOMAIN_LATE
    real(RP)              :: DOMAIN_LONS, DOMAIN_LONE
    real(RP)              :: lon_swap
    character(len=H_LONG) :: fname

    real(RP), allocatable :: LAT_1d(:)
    real(RP), allocatable :: LON_1d(:)

    integer :: ij, l, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNV2D_tile_init",*) 'Setup tiled data'

    CNV2D_ftype      = I_TILE
    CNV2D_TILE_dtype = dtype
    CNV2D_TILE_dir   = dir
    CNV2D_TILE_dlat  = dlat * D2R
    CNV2D_TILE_dlon  = dlon * D2R

    ! read catalogue file
    DOMAIN_LATS =   90.0_RP * D2R
    DOMAIN_LATE =  -90.0_RP * D2R
    DOMAIN_LONS =  180.0_RP * D2R
    DOMAIN_LONE = -180.0_RP * D2R

    do l = 1, ADM_lall
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       ij = suf(i,j)

       if ( GRD_LON(ij,l) < 0.0_RP ) then
          lon_swap = GRD_LON(ij,l) + 2.0_RP * PI
       else
          lon_swap = GRD_LON(ij,l)
       endif

       DOMAIN_LATS = min( DOMAIN_LATS, GRD_LAT(ij,l) )
       DOMAIN_LATE = max( DOMAIN_LATE, GRD_LAT(ij,l) )
       DOMAIN_LONS = min( DOMAIN_LONS, GRD_LON(ij,l) )
       DOMAIN_LONE = max( DOMAIN_LONE, GRD_LON(ij,l) )
    enddo
    enddo
    enddo

    if ( PRC_have_pl ) then
    ij = ADM_gslf_pl
    do l = 1, ADM_lall_pl

       if ( GRD_LON_pl(ij,l) < 0.0_RP ) then
          lon_swap = GRD_LON_pl(ij,l) + 2.0_RP * PI
       else
          lon_swap = GRD_LON_pl(ij,l)
       endif

       DOMAIN_LATS = min( DOMAIN_LATS, GRD_LAT_pl(ij,l) )
       DOMAIN_LATE = max( DOMAIN_LATE, GRD_LAT_pl(ij,l) )
       DOMAIN_LONS = min( DOMAIN_LONS, GRD_LON_pl(ij,l) )
       DOMAIN_LONE = max( DOMAIN_LONE, GRD_LON_pl(ij,l) )
    enddo
    endif

    fname = trim(dir)//'/'//trim(catalogue)

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
                                 CNV2D_nlat, CNV2D_nlon,           & ! [OUT]
                                 TILEDATA_DOMAIN_JS,               & ! [OUT]
                                 TILEDATA_DOMAIN_JE,               & ! [OUT]
                                 TILEDATA_DOMAIN_IS,               & ! [OUT]
                                 TILEDATA_DOMAIN_IE,               & ! [OUT]
                                 TILEDATA_zonal, TILEDATA_pole     ) ! [OUT]

    LOG_INFO("CNV2D_tile_init",*) 'CNV2D_nlat = ', CNV2D_nlat
    LOG_INFO("CNV2D_tile_init",*) 'CNV2D_nlon = ', CNV2D_nlon

    allocate( LAT_org (CNV2D_nlon,CNV2D_nlat) )
    allocate( LON_org (CNV2D_nlon,CNV2D_nlat) )
    allocate( DATA_org(CNV2D_nlon,CNV2D_nlat) )

    allocate( LAT_1d(CNV2D_nlat) )
    allocate( LON_1d(CNV2D_nlon) )

    ! read lat, lon
    call FILE_TILEDATA_get_latlon( CNV2D_nlat, CNV2D_nlon,                 & ! [IN]
                                   TILEDATA_DOMAIN_JS, TILEDATA_DOMAIN_IS, & ! [IN]
                                   CNV2D_TILE_dlat, CNV2D_TILE_dlon,       & ! [IN]
                                   LAT_1d(:), LON_1d(:)                    ) ! [OUT]

    do j = 1, CNV2D_nlat
    do i = 1, CNV2D_nlon
       LAT_org(i,j) = LAT_1d(j)
       LON_org(i,j) = LON_1d(i)
    enddo
    enddo

    deallocate( LAT_1d )
    deallocate( LON_1d )

    call CNV2D_init( interp_type, interp_level )

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
    character(len=*), intent(in) :: interp_type
    integer,          intent(in) :: interp_level
    real(RP),         intent(in), optional :: search_limit
    character(len=*), intent(in), optional :: POSTFIX

    integer :: file_id, var_id
    integer :: shape(2)

    real(RP), allocatable :: LAT_1d(:)
    real(RP), allocatable :: LON_1d(:)

    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNV2D_GRADS_init",*) 'Setup GrADS format data'

    call FILE_GrADS_open( FILE_NAME, & ! [IN]
                          file_id    ) ! [OUT]

    call FILE_GrADS_varid( file_id, VAR_NAME, & ! [IN]
                           var_id             ) ! [OUT]

    call FILE_GrADS_get_shape( file_id, VAR_NAME, & ! [IN]
                               shape(:)           ) ! [OUT]

    CNV2D_ftype     = I_GRADS
    CNV2D_GRADS_fid = file_id
    CNV2D_GRADS_vid = var_id
    CNV2D_nlon      = shape(1)
    CNV2D_nlat      = shape(2)

    LOG_INFO("CNV2D_GRADS_init",*) 'CNV2D_nlat = ', CNV2D_nlat
    LOG_INFO("CNV2D_GRADS_init",*) 'CNV2D_nlon = ', CNV2D_nlon

    allocate( LAT_org (CNV2D_nlon,CNV2D_nlat) )
    allocate( LON_org (CNV2D_nlon,CNV2D_nlat) )
    allocate( DATA_org(CNV2D_nlon,CNV2D_nlat) )

    allocate( LAT_1d(CNV2D_nlat) )
    allocate( LON_1d(CNV2D_nlon) )

    ! read lat
    call FILE_GrADS_varid( file_id, LAT_NAME, & ! [IN]
                           var_id             ) ! [OUT]

    if ( FILE_GrADS_isOneD( file_id, var_id ) ) then
       call FILE_GrADS_read( file_id, var_id, & ! [IN]
                             LAT_1d(:)        ) ! [OUT]

       do j = 1, CNV2D_nlat
       do i = 1, CNV2D_nlon
          LAT_org(i,j) = LAT_1d(j) * D2R
       enddo
       enddo
    else
       call FILE_GrADS_read( file_id, var_id,  & ! [IN]
                             LAT_org(:,:),     & ! [OUT]
                             postfix = POSTFIX ) ! [IN]

       do j = 1, CNV2D_nlat
       do i = 1, CNV2D_nlon
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

       do j = 1, CNV2D_nlat
       do i = 1, CNV2D_nlon
          LON_org(i,j) = LON_1d(i) * D2R
       enddo
       enddo
    else
       call FILE_GrADS_read( file_id, var_id,  & ! [IN]
                             LON_org(:,:),     & ! [OUT]
                             postfix = POSTFIX ) ! [IN]

       do j = 1, CNV2D_nlat
       do i = 1, CNV2D_nlon
          LON_org(i,j) = LON_org(i,j) * D2R
       enddo
       enddo
    endif

    deallocate( LAT_1d )
    deallocate( LON_1d )

    call CNV2D_init( interp_type, interp_level )

    return
  end subroutine CNV2D_GRADS_init

  !-----------------------------------------------------------------------------
  subroutine CNV2D_init( &
       interp_type, &
       interp_level )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_RGN_ndiamond, &
       PRC_have_pl
    use scale_const, only: &
       PI  => CONST_PI,  &
       EPS => CONST_EPS
    use scale_sort, only: &
       SORT_exec
    use mod_grd, only: &
       GRD_LAT,    &
       GRD_LAT_pl, &
       GRD_LON,    &
       GRD_LON_pl
    implicit none

    character(len=*), intent(in) :: interp_type
    integer,          intent(in) :: interp_level

    real(RP) :: lon_swap
    real(RP) :: a, c

    real(RP), allocatable :: dist_list(:)
    integer,  allocatable :: idxi_list(:)
    integer,  allocatable :: idxj_list(:)
    real(RP) :: dlat, dlon, wk
    real(RP) :: ref_dist, dist, sum
    integer  :: global_grid

    integer  :: ij, l, n, i, j, ii, jj
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNV2D_init",*) 'Calculate interpolation factor'

    select case( interp_type )
    case('LINEAR')
       LOG_INFO("CNV2D_init",*) 'interp_type is LINEAR'
       npoints = 4
    case ( 'DIST-WEIGHT' )
       LOG_INFO("CNV2D_init",*) 'interp_type is DIST-WEIGHT'
       npoints = interp_level

       global_grid = PRC_RGN_ndiamond * 4**ADM_glevel + 2
       ref_dist    = sqrt( 4.0_RP*PI / real(global_grid,kind=RP) )

    case default
       LOG_ERROR('CNV2D_init',*) 'interp_type is invalid: ', trim(interp_type)
       call PRC_abort
    end select

    allocate( idx_i   (ADM_gall   ,ADM_lall   ,npoints) )
    allocate( idx_j   (ADM_gall   ,ADM_lall   ,npoints) )
    allocate( hfact   (ADM_gall   ,ADM_lall   ,npoints) )
    allocate( idx_i_pl(ADM_gall_pl,ADM_lall_pl,npoints) )
    allocate( idx_j_pl(ADM_gall_pl,ADM_lall_pl,npoints) )
    allocate( hfact_pl(ADM_gall_pl,ADM_lall_pl,npoints) )
    idx_i   (:,:,:) = -1
    idx_j   (:,:,:) = -1
    hfact   (:,:,:) = 0.0_RP
    idx_i_pl(:,:,:) = -1
    idx_j_pl(:,:,:) = -1
    hfact_pl(:,:,:) = 0.0_RP

    select case( interp_type )
    case('LINEAR')

       do l = 1, ADM_lall
       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij = suf(i,j)

          if ( GRD_LON(ij,l) < 0.0_RP ) then
             lon_swap = GRD_LON(ij,l) + 2.0_RP * PI
          else
             lon_swap = GRD_LON(ij,l)
          endif

          if ( GRD_LAT(ij,l) <= LAT_org(1,1) ) then

             c = 0.D0
             idx_j(ij,l,1) = 1
             idx_j(ij,l,2) = 1
             idx_j(ij,l,3) = 1
             idx_j(ij,l,4) = 1

          elseif( GRD_LAT(ij,l) > LAT_org(1,CNV2D_nlat) ) then

             c = 1.D0
             idx_j(ij,l,1) = CNV2D_nlat
             idx_j(ij,l,2) = CNV2D_nlat
             idx_j(ij,l,3) = CNV2D_nlat
             idx_j(ij,l,4) = CNV2D_nlat

          else

             do jj = 1, CNV2D_nlat-1
                if (       GRD_LAT(ij,l) >  LAT_org(1,jj  ) &
                     .AND. GRD_LAT(ij,l) <= LAT_org(1,jj+1) ) then

                   c = ( GRD_LAT(ij,l) - LAT_org(1,jj) ) / ( LAT_org(1,jj+1) - LAT_org(1,jj) )
                   idx_j(ij,l,1) = jj+1
                   idx_j(ij,l,2) = jj+1
                   idx_j(ij,l,3) = jj
                   idx_j(ij,l,4) = jj

                   exit
                endif
             enddo

          endif

          if ( lon_swap <= LON_org(1,1) ) then

             a = 0.D0
             idx_i(ij,l,1) = 1
             idx_i(ij,l,2) = 1
             idx_i(ij,l,3) = 1
             idx_i(ij,l,4) = 1

          elseif( lon_swap >  LON_org(CNV2D_nlon,1) ) then

             a = 1.D0
             idx_i(ij,l,1) = CNV2D_nlon
             idx_i(ij,l,2) = CNV2D_nlon
             idx_i(ij,l,3) = CNV2D_nlon
             idx_i(ij,l,4) = CNV2D_nlon

          else
             do ii = 1, CNV2D_nlon-1
                if (       lon_swap >  LON_org(ii  ,1) &
                     .AND. lon_swap <= LON_org(ii+1,1) ) then

                   a = ( lon_swap - LON_org(ii,1) ) / ( LON_org(ii+1,1) - LON_org(ii,1) )
                   idx_i(ij,l,1) = ii+1
                   idx_i(ij,l,2) = ii
                   idx_i(ij,l,3) = ii+1
                   idx_i(ij,l,4) = ii

                   exit
                endif
             enddo
          endif

          hfact(ij,l,1) = (      c ) * (      a )
          hfact(ij,l,2) = (      c ) * ( 1.D0-a )
          hfact(ij,l,3) = ( 1.D0-c ) * (      a )
          hfact(ij,l,4) = ( 1.D0-c ) * ( 1.D0-a )
       enddo
       enddo
       enddo

       if ( PRC_have_pl ) then
       ij = ADM_gslf_pl
       do l = 1, ADM_lall_pl

          if ( GRD_LON_pl(ij,l) < 0.0_RP ) then
             lon_swap = GRD_LON_pl(ij,l) + 2.0_RP * PI
          else
             lon_swap = GRD_LON_pl(ij,l)
          endif

          if ( GRD_LAT_pl(ij,l) <= LAT_org(1,1) ) then

             c = 0.D0
             idx_j_pl(ij,l,1) = 1
             idx_j_pl(ij,l,2) = 1
             idx_j_pl(ij,l,3) = 1
             idx_j_pl(ij,l,4) = 1

          elseif( GRD_LAT_pl(ij,l) > LAT_org(1,CNV2D_nlat) ) then

             c = 1.D0
             idx_j_pl(ij,l,1) = CNV2D_nlat
             idx_j_pl(ij,l,2) = CNV2D_nlat
             idx_j_pl(ij,l,3) = CNV2D_nlat
             idx_j_pl(ij,l,4) = CNV2D_nlat

          else

             do jj = 1, CNV2D_nlat-1
                if (       GRD_LAT_pl(ij,l) >  LAT_org(1,jj  ) &
                     .AND. GRD_LAT_pl(ij,l) <= LAT_org(1,jj+1) ) then

                   c = ( GRD_LAT_pl(ij,l) - LAT_org(1,jj) ) / ( LAT_org(1,jj+1) - LAT_org(1,jj) )
                   idx_j_pl(ij,l,1) = jj+1
                   idx_j_pl(ij,l,2) = jj+1
                   idx_j_pl(ij,l,3) = jj
                   idx_j_pl(ij,l,4) = jj

                   exit
                endif
             enddo

          endif

          if ( lon_swap <= LON_org(1,1) ) then

             a = 0.D0
             idx_i_pl(ij,l,1) = 1
             idx_i_pl(ij,l,2) = 1
             idx_i_pl(ij,l,3) = 1
             idx_i_pl(ij,l,4) = 1

          elseif( lon_swap >  LON_org(CNV2D_nlon,1) ) then

             a = 1.D0
             idx_i_pl(ij,l,1) = CNV2D_nlon
             idx_i_pl(ij,l,2) = CNV2D_nlon
             idx_i_pl(ij,l,3) = CNV2D_nlon
             idx_i_pl(ij,l,4) = CNV2D_nlon

          else
             do ii = 1, CNV2D_nlon-1
                if (       lon_swap >  LON_org(ii  ,1) &
                     .AND. lon_swap <= LON_org(ii+1,1) ) then

                   a = ( lon_swap - LON_org(ii,1) ) / ( LON_org(ii+1,1) - LON_org(ii,1) )
                   idx_i_pl(ij,l,1) = ii+1
                   idx_i_pl(ij,l,2) = ii
                   idx_i_pl(ij,l,3) = ii+1
                   idx_i_pl(ij,l,4) = ii

                   exit
                endif
             enddo
          endif

          hfact_pl(ij,l,1) = (      c ) * (      a )
          hfact_pl(ij,l,2) = (      c ) * ( 1.D0-a )
          hfact_pl(ij,l,3) = ( 1.D0-c ) * (      a )
          hfact_pl(ij,l,4) = ( 1.D0-c ) * ( 1.D0-a )
       enddo
       endif

    case('DIST-WEIGHT')

       allocate( dist_list(npoints) )
       allocate( idxi_list(npoints) )
       allocate( idxj_list(npoints) )

       do l = 1, ADM_lall
       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij = suf(i,j)

          dist_list(:) = 9.999E+15_RP
          idxi_list(:) = -1
          idxj_list(:) = -1

          do jj = 1, CNV2D_nlat
          do ii = 1, CNV2D_nlon
             dlat = GRD_LAT(ij,l) - LAT_org(ii,jj)
             dlon = GRD_LON(ij,l) - LON_org(ii,jj) + 4.0_RP * PI
             if    ( dlon >= 4.0_RP * PI ) then
                dlon = dlon - 4.0_RP * PI
             elseif( dlon >= 2.0_RP * PI ) then
                dlon = dlon - 2.0_RP * PI
             endif

             if( abs( dlat ) > ref_dist ) cycle
             if( abs( dlon ) > ref_dist ) cycle

             wk = sin(0.5_RP*dlat)**2 + cos(GRD_LAT(ij,l)) * cos(LAT_org(ii,jj)) * sin(0.5_RP*dlon)**2

             dist = 2.0_RP * asin( min(sqrt(wk),1.0_RP) )

             if ( dist <= ref_dist ) then
                if ( dist <= dist_list(npoints) ) then ! replace last(=longest) value
                   dist_list(npoints) = dist
                   idxi_list(npoints) = ii
                   idxj_list(npoints) = jj

                   ! sort by ascending order
                   call SORT_exec( npoints,      & ! [IN]
                                   dist_list(:), & ! [INOUT]
                                   idxi_list(:), & ! [INOUT]
                                   idxj_list(:)  ) ! [INOUT]
                endif
             endif
          enddo
          enddo

          do n = 1, npoints
             idx_i(ij,l,n) = idxi_list(n)
             idx_j(ij,l,n) = idxj_list(n)
          enddo

          if ( abs(dist_list(1)) < EPS ) then
             hfact(ij,l,:) = 0.0_RP
             hfact(ij,l,1) = 1.0_RP
          else
             sum = 0.0_RP
             do n = 1, npoints
                hfact(ij,l,n) = 1.0_RP / dist_list(n)

                sum = sum + hfact(ij,l,n)
             enddo

             if ( sum > 0.0_RP ) then ! normalize
                do n = 1, npoints
                   hfact(ij,l,n) = hfact(ij,l,n) / sum
                enddo
             endif
          endif
       enddo
       enddo
       enddo

       if ( PRC_have_pl ) then
       ij = ADM_gslf_pl
       do l = 1, ADM_lall_pl

          dist_list(:) = 9.999E+15_RP
          idxi_list(:) = -1
          idxj_list(:) = -1

          do jj = 1, CNV2D_nlat
          do ii = 1, CNV2D_nlon
             dlat = GRD_LAT_pl(ij,l) - LAT_org(ii,jj)
             dlon = GRD_LON_pl(ij,l) - LON_org(ii,jj) + 4.0_RP * PI
             if    ( dlon >= 4.0_RP * PI ) then
                dlon = dlon - 4.0_RP * PI
             elseif( dlon >= 2.0_RP * PI ) then
                dlon = dlon - 2.0_RP * PI
             endif

             if( abs( dlat ) > ref_dist ) cycle
             if( abs( dlon ) > ref_dist ) cycle

             wk = sin(0.5_RP*dlat)**2 + cos(GRD_LAT_pl(ij,l)) * cos(LAT_org(ii,jj)) * sin(0.5_RP*dlon)**2

             dist = 2.0_RP * asin( min(sqrt(wk),1.0_RP) )

             if ( dist <= ref_dist ) then
                if ( dist <= dist_list(npoints) ) then ! replace last(=longest) value
                   dist_list(npoints) = dist
                   idxi_list(npoints) = ii
                   idxj_list(npoints) = jj

                   ! sort by ascending order
                   call SORT_exec( npoints,      & ! [IN]
                                   dist_list(:), & ! [INOUT]
                                   idxi_list(:), & ! [INOUT]
                                   idxj_list(:)  ) ! [INOUT]
                endif
             endif
          enddo
          enddo

          do n = 1, npoints
             idx_i_pl(ij,l,n) = idxi_list(n)
             idx_j_pl(ij,l,n) = idxj_list(n)
          enddo

          if ( abs(dist_list(1)) < EPS ) then
             hfact_pl(ij,l,:) = 0.0_RP
             hfact_pl(ij,l,1) = 1.0_RP
          else
             sum = 0.0_RP
             do n = 1, npoints
                hfact_pl(ij,l,n) = 1.0_RP / dist_list(n)

                sum = sum + hfact_pl(ij,l,n)
             enddo

             if ( sum > 0.0_RP ) then ! normalize
                do n = 1, npoints
                   hfact_pl(ij,l,n) = hfact_pl(ij,l,n) / sum
                enddo
             endif
          endif
       enddo
       endif

       deallocate( dist_list )
       deallocate( idxi_list )
       deallocate( idxj_list )

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
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
    use scale_prc_icoA, only: &
       PRC_have_pl
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

    real(RP) :: wgt, val
    real(RP) :: wgt_sum, val_sum
    integer  :: i, j, l, n, ij, i2, j2
    !---------------------------------------------------------------------------

    select case(CNV2D_ftype)
    case(I_TILE)

       call FILE_TILEDATA_get_data( CNV2D_nlat, CNV2D_nlon, & ! [IN]
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
    do l = 1, ADM_lall
    do j = ADM_gmin, ADM_gmax
    do i = ADM_gmin, ADM_gmax
       ij = suf(i,j)

       wgt_sum   = 0.0_RP
       val_sum   = 0.0_RP
       var(ij,l) = 0.0_RP

       do n = 1, npoints
          i2  = idx_i(ij,l,n)
          j2  = idx_j(ij,l,n)
          wgt = hfact(ij,l,n)

          if ( wgt > EPS ) then
             val = DATA_org(i2,j2)

             if ( val-UNDEF > EPS ) then
                wgt_sum = wgt_sum + wgt
                val_sum = val_sum + wgt * val
             endif
          endif
       enddo

       if ( wgt_sum > EPS ) then
          var(ij,l) = val_sum / wgt_sum
       endif
    enddo
    enddo
    enddo

    if ( PRC_have_pl ) then
    ij = ADM_gslf_pl
    do l = 1, ADM_lall_pl

       wgt_sum      = 0.0_RP
       val_sum      = 0.0_RP
       var_pl(ij,l) = 0.0_RP

       do n = 1, npoints
          i2  = idx_i_pl(ij,l,n)
          j2  = idx_j_pl(ij,l,n)
          wgt = hfact_pl(ij,l,n)

          if ( wgt > EPS ) then
             val = DATA_org(i2,j2)

             if ( val-UNDEF > EPS ) then
                wgt_sum = wgt_sum + wgt
                val_sum = val_sum + wgt * val
             endif
          endif
       enddo

       if ( wgt_sum > EPS ) then
          var_pl(ij,l) = val_sum / wgt_sum
       endif
    enddo
    endif

    return
  end subroutine CNV2D_convert

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_gm_cnv2d
