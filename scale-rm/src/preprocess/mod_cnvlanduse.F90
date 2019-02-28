!-------------------------------------------------------------------------------
!> module Convert LandUseIndex
!!
!! @par Description
!!          subroutines for preparing landuse index data (convert from external file)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_cnvlanduse
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
  public :: CNVLANDUSE_setup
  public :: CNVLANDUSE

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: CNVLANDUSE_DoNothing
  logical, public :: CNVLANDUSE_UseGLCCv2 = .false.
  logical, public :: CNVLANDUSE_UseLU100M = .false.
  logical, public :: CNVLANDUSE_UseJIBIS  = .false.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CNVLANDUSE_GLCCv2
  private :: CNVLANDUSE_LU100M
  private :: CNVLANDUSE_JIBIS

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private :: CNVLANDUSE_limit_urban_fraction = 1.0_RP !< fraction limiter for urban area

  real(RP), private :: DOMAIN_LATS, DOMAIN_LATE
  real(RP), private :: DOMAIN_LONS, DOMAIN_LONE
  real(RP), private :: DOMAIN_DLAT
  real(RP), private, allocatable :: DOMAIN_DXY(:,:)

  real(RP), private, parameter :: d_large = 1e20_RP
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CNVLANDUSE_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=H_SHORT) :: CNVLANDUSE_name = 'NONE' ! keep backward compatibility

    namelist / PARAM_CNVLANDUSE / &
       CNVLANDUSE_name,                &
       CNVLANDUSE_UseGLCCv2,           &
       CNVLANDUSE_UseLU100M,           &
!       CNVLANDUSE_UseJIBIS
       CNVLANDUSE_limit_urban_fraction

    integer  :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CNVLANDUSE_setup",*) 'Setup'
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVLANDUSE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVLANDUSE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVLANDUSE_setup",*) 'Not appropriate names in namelist PARAM_CNVLANDUSE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVLANDUSE)

    select case(CNVLANDUSE_name)
    case('NONE')
       ! do nothing
    case('GLCCv2')
       CNVLANDUSE_UseGLCCv2 = .true.
       CNVLANDUSE_UseLU100M = .false.
       CNVLANDUSE_UseJIBIS  = .false.
    case('LU100M')
       CNVLANDUSE_UseGLCCv2 = .false.
       CNVLANDUSE_UseLU100M = .true.
       CNVLANDUSE_UseJIBIS  = .false.
    case('COMBINE')
       CNVLANDUSE_UseGLCCv2 = .true.
       CNVLANDUSE_UseLU100M = .true.
       CNVLANDUSE_UseJIBIS  = .true.
    case default
       LOG_ERROR("CNVLANDUSE_setup",*) 'Unsupported TYPE: ', trim(CNVLANDUSE_name)
       call PRC_abort
    endselect

    CNVLANDUSE_DoNothing = .true.

    if ( CNVLANDUSE_UseGLCCv2 ) then
       CNVLANDUSE_DoNothing = .false.
       LOG_INFO("CNVLANDUSE_setup",*) 'Use GLCC ver.2, global 30 arcsec. data'
       if ( CNVLANDUSE_UseLU100M ) then
          LOG_INFO("CNVLANDUSE_setup",*) 'Use KSJ landuse 100m data for Japan region'
          LOG_INFO("CNVLANDUSE_setup",*) 'Overwrite Japan region'
          if ( CNVLANDUSE_UseJIBIS ) then
             LOG_INFO("CNVLANDUSE_setup",*) 'Use J-IBIS map 100m data for Japan region'
             LOG_INFO("CNVLANDUSE_setup",*) 'Overwrite Japan region (PFT only)'
          endif
       endif
    elseif ( CNVLANDUSE_UseLU100M ) then
       CNVLANDUSE_DoNothing = .false.
       LOG_INFO("CNVLANDUSE_setup",*) 'Use KSJ landuse 100m data, Japan region only'
       if ( CNVLANDUSE_UseJIBIS ) then
          LOG_INFO("CNVLANDUSE_setup",*) 'Use J-IBIS map 100m data for Japan region'
          LOG_INFO("CNVLANDUSE_setup",*) 'Overwrite Japan region (PFT only)'
       endif
    endif

    if ( CNVLANDUSE_DoNothing ) then
       LOG_INFO("CNVLANDUSE_setup",*) 'Do nothing for landuse index'
    endif

    return
  end subroutine CNVLANDUSE_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CNVLANDUSE
    use scale_const, only: &
       EPS => CONST_EPS, &
       D2R => CONST_D2R
    use scale_prc, only: &
       PRC_abort
    use scale_sort, only: &
       SORT_exec
    use scale_landuse, only: &
       LANDUSE_frac_land,  &
       LANDUSE_frac_lake,  &
       LANDUSE_frac_urban, &
       LANDUSE_PFT_mosaic, &
       LANDUSE_PFT_nmax,   &
       LANDUSE_frac_PFT,   &
       LANDUSE_index_PFT,  &
       LANDUSE_calc_fact,  &
       LANDUSE_fillhalo
    use scale_atmos_grid_cartesC_real, only: &
       LATXV => ATMOS_GRID_CARTESC_REAL_LATXV, &
       LONUY => ATMOS_GRID_CARTESC_REAL_LONUY, &
       DLAT  => ATMOS_GRID_CARTESC_REAL_DLAT,  &
       AREA  => ATMOS_GRID_CARTESC_REAL_AREA
    implicit none

    real(RP) :: PFT_weight(-2:LANDUSE_PFT_nmax,IA,JA)
    integer  :: PFT_idx(LANDUSE_PFT_nmax)
    real(RP) :: lake_wgt, ocean_wgt, urban_wgt, land_wgt
    real(RP) :: allsum
    real(RP) :: zerosw

    integer :: i, j, p
    !---------------------------------------------------------------------------

    if ( CNVLANDUSE_DoNothing ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'skip  convert landuse data'
    else
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start convert landuse data'

       DOMAIN_LATS = minval( LATXV(:,:) )
       DOMAIN_LATE = maxval( LATXV(:,:) )
       DOMAIN_LONS = minval( LONUY(:,:) )
       DOMAIN_LONE = maxval( LONUY(:,:) )
       DOMAIN_DLAT = maxval( DLAT(:,:) )

       allocate( DOMAIN_DXY(IA,JA) )
       do j = 1, JA
       do i = 1, IA
          DOMAIN_DXY(i,j) = sqrt( AREA(i,j) )
       end do
       end do

       LOG_INFO("CNVLANDUSE",*) 'Domain Information'
       LOG_INFO_CONT(*) 'Domain (LAT)    :', DOMAIN_LATS/D2R, DOMAIN_LATE/D2R
       LOG_INFO_CONT(*) '       (LON)    :', DOMAIN_LONS/D2R, DOMAIN_LONE/D2R


       !$omp parallel do
!OCL XFILL
       do j = 1, JA
       do i = 1, IA
          PFT_weight(:,i,j) = 0.0_RP
       end do
       end do

       if ( CNVLANDUSE_UseGLCCv2 ) then
          call CNVLANDUSE_GLCCv2( PFT_weight(:,:,:) ) ! [INOUT]
       endif

       if ( CNVLANDUSE_UseLU100M ) then
          call CNVLANDUSE_LU100M( PFT_weight(:,:,:) ) ! [INOUT]
       endif

       if ( CNVLANDUSE_UseJIBIS ) then
          call CNVLANDUSE_JIBIS( PFT_weight(:,:,:) ) ! [INOUT]
       endif

       deallocate( DOMAIN_DXY )

       !$omp parallel do &
       !$omp private(lake_wgt,ocean_wgt,urban_wgt,land_wgt,allsum,zerosw,PFT_idx)
       do j = JS, JE
       do i = IS, IE

          lake_wgt  = PFT_weight(-2,i,j)
          ocean_wgt = PFT_weight(-1,i,j)
          urban_wgt = PFT_weight( 0,i,j)
          land_wgt  = sum( PFT_weight(1:,i,j) )

          do p = 1, LANDUSE_PFT_nmax
             PFT_idx(p) = p
          end do
          call SORT_exec( LANDUSE_PFT_nmax, PFT_weight(1:,i,j), PFT_idx(:) )


          ! land fraction : 1 - ocean / total
          allsum = lake_wgt + ocean_wgt + urban_wgt + land_wgt
          zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
          LANDUSE_frac_land (i,j) = ( allsum-ocean_wgt ) * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

          ! lake fraction : lake / ( total - ocean )
          allsum = lake_wgt + urban_wgt + land_wgt
          zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
          LANDUSE_frac_lake (i,j) = lake_wgt * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

          ! urban fraction : urban / ( total - ocean - lake )
          allsum = urban_wgt + land_wgt
          zerosw = 0.5_RP - sign( 0.5_RP, allsum-EPS )
          LANDUSE_frac_urban(i,j) = urban_wgt * ( 1.0_RP-zerosw ) / ( allsum-zerosw )

          ! PFT fraction : PFT / sum( PFT(1:mosaic) )
          allsum = sum( PFT_weight(LANDUSE_PFT_nmax-LANDUSE_PFT_mosaic+1:,i,j) )
          if ( allsum > EPS ) then
             do p = 1, LANDUSE_PFT_mosaic
                LANDUSE_frac_PFT (i,j,p) = PFT_weight(LANDUSE_PFT_nmax-p+1,i,j) / allsum
                LANDUSE_index_PFT(i,j,p) = PFT_idx(LANDUSE_PFT_nmax-p+1)
             enddo
             ! if no second PFT, set to same as PFT1
             if ( abs(LANDUSE_frac_PFT(i,j,1)-1.0_RP) <= EPS ) then
                LANDUSE_frac_PFT (i,j,:) = 0.0_RP
                LANDUSE_frac_PFT (i,j,1) = 1.0_RP
                LANDUSE_index_PFT(i,j,:) = PFT_idx(LANDUSE_PFT_nmax)
             endif
          else ! if no PFT, set to bare ground
             LANDUSE_frac_PFT (i,j,:) = 0.0_RP
             LANDUSE_frac_PFT (i,j,1) = 1.0_RP
             LANDUSE_index_PFT(i,j,:) = 1
          endif

       enddo
       enddo

       if ( CNVLANDUSE_limit_urban_fraction < 1.0_RP ) then
          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
             if ( LANDUSE_frac_urban(i,j) == 1.0_RP ) then ! if no PFT, set to grassland
                LANDUSE_frac_PFT (i,j,:) = 0.0_RP
                LANDUSE_frac_PFT (i,j,1) = 1.0_RP
                LANDUSE_index_PFT(i,j,:) = 2     ! Grassland
             endif
             LANDUSE_frac_urban(i,j) = min( LANDUSE_frac_urban(i,j), CNVLANDUSE_limit_urban_fraction )
          enddo
          enddo
       endif

       ! calculate landuse factors
       call LANDUSE_fillhalo( FILL_BND=.true. )
       call LANDUSE_calc_fact

       LOG_PROGRESS(*) 'end   convert landuse data'

    endif

    return
  end subroutine CNVLANDUSE

  !-----------------------------------------------------------------------------
  !> Convert from GLCCv2
  subroutine CNVLANDUSE_GLCCv2( PFT_weight )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       RADIUS => CONST_RADIUS, &
       D2R    => CONST_D2R
    use scale_atmos_grid_cartesC, only: &
       FX => ATMOS_GRID_CARTESC_FX, &
       FY => ATMOS_GRID_CARTESC_FY, &
       CX => ATMOS_GRID_CARTESC_CX, &
       CY => ATMOS_GRID_CARTESC_CY
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info, &
       FILE_TILEDATA_get_data
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
    use scale_landuse, only: &
       LANDUSE_PFT_nmax
    implicit none

    real(RP), intent(inout) :: PFT_weight(-2:LANDUSE_PFT_nmax,IA,JA)

    character(len=H_LONG) :: GLCCv2_IN_DIR        = '.'    !< directory contains GLCCv2 files (GrADS format)
    character(len=H_LONG) :: GLCCv2_IN_CATALOGUE  = ''     !< metadata files for GLCCv2

    namelist / PARAM_CNVLANDUSE_GLCCv2 / &
       GLCCv2_IN_DIR,       &
       GLCCv2_IN_CATALOGUE

    ! GLCCv2 data
    integer,  parameter :: categ_nmax = 25
    real(RP), parameter :: GLCCv2_DLAT = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.
    real(RP), parameter :: GLCCv2_DLON = 30.0_RP / 60.0_RP / 60.0_RP ! 30 arc sec.

    integer  :: lookuptable(1:categ_nmax)
    data lookuptable /  0, & !   1 Urban and Built-Up Land        ->  0 urban
                        7, & !   2 Dryland Cropland and Pasture   ->  7 Dryland Cropland and Pasture
                        8, & !   3 Irrigated Cropland and Pasture ->  8 Irrigated Cropland and Pasture
                        9, & !   4 Mixed Cropland and Pasture     ->  9 Mixed Cropland and Pasture
                        5, & !   5 Cropland/Grassland Mosaic      ->  5 Cropland/Grassland Mosaic
                        6, & !   6 Cropland/Woodland Mosaic       ->  6 Cropland/Woodland Mosaic
                        2, & !   7 Grassland                      ->  2 Grassland
                        3, & !   8 Shrubland                      ->  3 Shrubland
                        4, & !   9 Mixed Shrubland/Grassland      ->  4 Mixed Shrubland/Grassland
                        4, & !  10 Savanna                        ->  4 Mixed Shrubland/Grassland
                       11, & !  11 Deciduous Broadleaf Forest     -> 11 Deciduous Broadleaf Forest
                       12, & !  12 Deciduous Needleleaf Forest    -> 12 Deciduous Needleleaf Forest
                       13, & !  13 Evergreen Broadleaf Forest     -> 13 Deciduous Broadleaf Forest
                       14, & !  14 Evergreen Needleleaf Forest    -> 14 Deciduous Needleleaf Forest
                       15, & !  15 Mixed Forest                   -> 15 Mixed Forest
                       -2, & !  16 Water Bodies                   -> -2 Lake/River
                       10, & !  17 Herbaceous Wetland             -> 10 Paddy
                       10, & !  18 Wooded Wetland                 -> 10 Paddy
                        1, & !  19 Barren or Sparsely Vegetated   ->  1 Dessert
                       16, & !  20 Herbaceous Tundra              -> 16 Tundra
                       16, & !  21 Wooded Tundra                  -> 16 Tundra
                       16, & !  22 Mixed Tundra                   -> 16 Tundra
                       16, & !  23 Bare Ground Tundra             -> 16 Tundra
                       17, & !  24 Snow or Ice                    -> 17 Gracier
                       -1  / !  25+ Sea Surface                   -> -1 Sea Surface

    !---------------------------------------------------------------------------

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 100
    integer               :: TILE_nmax
    character(len=H_LONG) :: TILE_fname(TILE_nlim)
    logical               :: TILE_hit  (TILE_nlim)
    integer               :: TILE_JS   (TILE_nlim)
    integer               :: TILE_JE   (TILE_nlim)
    integer               :: TILE_IS   (TILE_nlim)
    integer               :: TILE_IE   (TILE_nlim)
    real(RP)              :: TILE_DLAT, TILE_DLON

    integer,  allocatable :: LANDUSE(:,:)
    real(RP), allocatable :: LATH(:,:), LONH(:,:)
    real(RP), allocatable :: XH(:,:), YH(:,:)
    integer               :: nLONH, nLATH

    integer  :: GLOBAL_IA

    character(len=H_LONG) :: fname

    logical  :: hit, no_hit_x
    real(RP) :: dmin, d
    integer  :: min_i, min_j
    real(RP) :: limit

    integer :: ish, ieh, jsh, jeh
    logical :: zonal, pole
    integer :: lu
    integer :: ierr
    integer :: i, j, ii, jj, p
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVLANDUSE_GLCCv2,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVLANDUSE_GLCCv2",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVLANDUSE_GLCCv2",*) 'Not appropriate names in namelist PARAM_CNVLANDUSE_GLCCv2. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVLANDUSE_GLCCv2)


    TILE_DLAT = GLCCv2_DLAT * D2R
    TILE_DLON = GLCCv2_DLON * D2R

    ! catalogue file
    fname = trim(GLCCv2_IN_DIR)//'/'//trim(GLCCv2_IN_CATALOGUE)

    call FILE_TILEDATA_get_info( TILE_nlim,                                          & ! [IN]
                                 TILE_DLAT, TILE_DLON,                               & ! [IN]
                                 DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                 fname,                                              & ! [IN]
                                 GLOBAL_IA,                                          & ! [OUT]
                                 TILE_nmax,                                          & ! [OUT]
                                 TILE_fname(:), TILE_hit(:),                         & ! [OUT]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                 nLATH, nLONH, jsh, jeh, ish, ieh, zonal, pole       ) ! [OUT]

    allocate( LANDUSE(nLONH,nLATH) )
    allocate( LATH   (nLONH,nLATH) )
    allocate( LONH   (nLONH,nLATH) )
    allocate( YH     (nLONH,nLATH) )
    allocate( XH     (nLONH,nLATH) )

    call FILE_TILEDATA_get_data( nLATH, nLONH,                                   & ! [IN]
                                 GLCCv2_IN_DIR,                                  & ! [IN]
                                 GLOBAL_IA,                                      & ! [IN]
                                 TILE_nmax,                                      & ! [IN]
                                 TILE_DLAT, TILE_DLON,                           & ! [IN]
                                 TILE_fname(:), TILE_hit(:),                     & ! [IN]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:), & ! [IN]
                                 jsh, jeh, ish, ieh,                             & ! [IN]
                                 "INT1",                                         & ! [IN]
                                 LANDUSE(:,:), LATH(:,:), LONH(:,:)              ) ! [OUT]

    call MAPPROJECTION_lonlat2xy( nLONH, 1, nLONH, nLATH, 1, nLATH, &
                                  LONH(:,:), LATH(:,:), & ! [IN]
                                  XH(:,:), YH(:,:)      ) ! [OUT]

    limit = ( TILE_DLAT * RADIUS * 1.5_RP )**2
    !$omp parallel do collapse(2) &
    !$omp private(hit,no_hit_x,lu,p,dmin,min_i,min_j,d)
    do j = 1, JA
    do i = 1, IA
       hit = .false.
       dmin = d_large
       min_i = -1
       min_j = -1
       do jj = 1, nLATH
          no_hit_x = .true.
          do ii = 1, nLONH
             if (       FX(i-1) <= XH(ii,jj) .and. XH(ii,jj) < FX(i) &
                  .and. FY(j-1) <= YH(ii,jj) .and. YH(ii,jj) < FY(j) ) then
                no_hit_x = .false.
                lu = min( LANDUSE(ii,jj), categ_nmax )
                if ( 1 <= lu ) then
                   p = lookuptable(lu)
                   PFT_weight(p,i,j) = PFT_weight(p,i,j) + cos(LATH(ii,jj)) ! area weight
                   hit = .true.
                   cycle
                end if
             end if
             d = ( XH(ii,jj)-CX(i) )**2 + ( YH(ii,jj)-CY(j) )**2
             lu = min( LANDUSE(ii,jj), categ_nmax )
             if ( d < dmin .and. 1 <= lu ) then
                dmin = d
                min_i = ii
                min_j = jj
             end if
          end do
          if ( hit .and. no_hit_x ) exit
       end do
       if ( ( .not. hit ) .and. dmin < limit ) then
          lu = min( LANDUSE(min_i,min_j), categ_nmax )
          p = lookuptable(lu)
          PFT_weight(p,i,j) = 1.0_RP
       end if
    end do
    end do

    deallocate( LANDUSE, LATH, LONH, YH, XH )

    return
  end subroutine CNVLANDUSE_GLCCv2

  !-----------------------------------------------------------------------------
  !> Convert from KSJ landuse 100m mesh
  subroutine CNVLANDUSE_LU100M( PFT_weight )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF2 => CONST_UNDEF2, &
       RADIUS => CONST_RADIUS, &
       D2R    => CONST_D2R
    use scale_atmos_grid_cartesC, only: &
       FX => ATMOS_GRID_CARTESC_FX, &
       FY => ATMOS_GRID_CARTESC_FY, &
       CX => ATMOS_GRID_CARTESC_CX, &
       CY => ATMOS_GRID_CARTESC_CY
    use scale_file_tiledata, only: &
       FILE_TILEDATA_get_info, &
       FILE_TILEDATA_get_data
    use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
    use scale_landuse, only: &
       LANDUSE_PFT_nmax
    implicit none

    real(RP), intent(inout) :: PFT_weight(-2:LANDUSE_PFT_nmax,IA,JA)

    character(len=H_LONG) :: LU100M_IN_DIR       = '.'     !< directory contains LU100M files (GrADS format)
    character(len=H_LONG) :: LU100M_IN_CATALOGUE = ''      !< metadata files for LU100M

    namelist / PARAM_CNVLANDUSE_LU100M / &
       LU100M_IN_DIR,       &
       LU100M_IN_CATALOGUE

    ! LU100M data
    integer,  parameter :: categ_nmax = 16
    real(RP), parameter :: LU100M_DLAT = 5.0_RP / 60.0_RP / 100.0_RP
    real(RP), parameter :: LU100M_DLON = 7.5_RP / 60.0_RP / 100.0_RP

    integer  :: lookuptable(0:16)
    data lookuptable / -1, & ! -999 missing      -> -1 Sea Surface
                       10, & !  1 paddy          -> 10 Paddy
                        9, & !  2 cropland       ->  9 Mixed Cropland and Pasture
                        1, & !  3 UNDEF          ->  1 Dessert
                        1, & !  4 UNDEF          ->  1 Dessert
                       11, & !  5 forest         -> 11 Deciduous Broadleaf Forest
                        1, & !  6 bareground     ->  1 Dessert
                        0, & !  7 urban building ->  0 Urban and Built-up Land
                        1, & !  8 UNDEF          ->  1 Dessert
                        0, & !  9 motorway       ->  0 Urban and Built-up Land
                        0, & ! 10 urban ground   ->  0 Urban and Built-up Land
                       -2, & ! 11 lake,river     -> -2 Lake/River
                        1, & ! 12 UNDEF          ->  1 Dessert
                        1, & ! 13 UNDEF          ->  1 Dessert
                        1, & ! 14 seashore       ->  1 Dessert
                       -1, & ! 15 ocean          -> -1 Sea Surface
                        2  / ! 16 golf course    ->  2 Grassland

    !---------------------------------------------------------------------------

    ! data catalogue list
    integer, parameter    :: TILE_nlim = 1000
    integer               :: TILE_nmax
    character(len=H_LONG) :: TILE_fname(TILE_nlim)
    logical               :: TILE_hit  (TILE_nlim)
    integer               :: TILE_JS   (TILE_nlim)
    integer               :: TILE_JE   (TILE_nlim)
    integer               :: TILE_IS   (TILE_nlim)
    integer               :: TILE_IE   (TILE_nlim)
    real(RP)              :: TILE_DLAT, TILE_DLON

    integer,  allocatable :: LANDUSE(:,:)
    real(RP), allocatable :: LATH   (:,:)
    real(RP), allocatable :: LONH   (:,:)
    real(RP), allocatable :: XH(:,:), YH(:,:)
    integer               :: nLONH, nLATH

    integer  :: GLOBAL_IA

    character(len=H_LONG) :: fname

    logical  :: hit, no_hit_x
    real(RP) :: dmin, d
    integer  :: min_i, min_j
    real(RP) :: limit

    integer :: ish, ieh, jsh, jeh
    logical :: zonal, pole
    integer :: lu
    integer :: ierr
    integer :: i, j, ii, jj, p
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CNVLANDUSE_LU100M,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CNVLANDUSE_LU100M",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CNVLANDUSE_LU100M",*) 'Not appropriate names in namelist PARAM_CNVLANDUSE_LU100M. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CNVLANDUSE_LU100M)


    TILE_DLAT = LU100M_DLAT * D2R
    TILE_DLON = LU100M_DLON * D2R

    ! catalogue file
    fname = trim(LU100M_IN_DIR)//'/'//trim(LU100M_IN_CATALOGUE)

    call FILE_TILEDATA_get_info( TILE_nlim,                                          & ! [IN]
                                 TILE_DLAT, TILE_DLON,                               & ! [IN]
                                 DOMAIN_LATS, DOMAIN_LATE, DOMAIN_LONS, DOMAIN_LONE, & ! [IN]
                                 fname,                                              & ! [IN]
                                 GLOBAL_IA,                                          & ! [OUT]
                                 TILE_nmax,                                          & ! [OUT]
                                 TILE_fname(:), TILE_hit(:),                         & ! [OUT]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:),     & ! [OUT]
                                 nLATH, nLONH, jsh, jeh, ish, ieh, zonal, pole       ) ! [OUT]

    if ( .not. any(TILE_hit(1:TILE_nmax) ) ) return

    allocate( LANDUSE(nLONH,nLATH) )
    allocate( LATH   (nLONH,nLATH) )
    allocate( LONH   (nLONH,nLATH) )
    allocate( YH     (nLONH,nLATH) )
    allocate( XH     (nLONH,nLATH) )

    call FILE_TILEDATA_get_data( nLATH, nLONH,                                   & ! [IN]
                                 LU100M_IN_DIR,                                  & ! [IN]
                                 GLOBAL_IA,                                      & ! [IN]
                                 TILE_nmax,                                      & ! [IN]
                                 TILE_DLAT, TILE_DLON,                           & ! [IN]
                                 TILE_fname(:), TILE_hit(:),                     & ! [IN]
                                 TILE_JS(:), TILE_JE(:), TILE_IS(:), TILE_IE(:), & ! [IN]
                                 jsh, jeh, ish, ieh,                             & ! [IN]
                                 "REAL4",                                        & ! [IN]
                                 LANDUSE(:,:), LATH(:,:), LONH(:,:),             & ! [OUT]
                                 min_value = -999                                ) ! [IN]

    call MAPPROJECTION_lonlat2xy( nLONH, 1, nLONH, nLATH, 1, nLATH, &
                                  LONH(:,:), LATH(:,:), & ! [IN]
                                  XH(:,:), YH(:,:)      ) ! [OUT]

    limit = ( TILE_DLAT * RADIUS * 1.5_RP )**2
    !$omp parallel do collapse(2) &
    !$omp private(hit,no_hit_x,lu,p,dmin,min_i,min_j,d)
    do j = 1, JA
    do i = 1, IA
       hit = .false.
       dmin = d_large
       min_i = -1
       min_j = -1
       do jj = 1, nLATH
          no_hit_x = .true.
          do ii = 1, nLONH
             if (       FX(i-1) <= XH(ii,jj) .and. XH(ii,jj) < FX(i) &
                  .and. FY(j-1) <= YH(ii,jj) .and. YH(ii,jj) < FY(j) ) then
                no_hit_x = .false.
                lu = LANDUSE(ii,jj)
                if ( lu /= UNDEF2 ) then
                   p = lookuptable( max(0,lu) ) ! -999 to 0
                   PFT_weight(p,i,j) = PFT_weight(p,i,j) + 1.0_RP
                   hit = .true.
                   cycle
                end if
             end if
             d = ( XH(ii,jj)-CX(i) )**2 + ( YH(ii,jj)-CY(j) )**2
             lu = LANDUSE(ii,jj)
             if ( d < dmin .and. lu /= UNDEF2 ) then
                dmin = d
                min_i = ii
                min_j = jj
             end if
          end do
          if ( hit .and. no_hit_x ) exit
       end do
       if ( ( .not. hit ) .and. dmin < limit ) then
          lu = LANDUSE(min_i,min_j)
          if ( lu /= UNDEF2 ) then
             p = lookuptable( max(0,lu) ) ! -999 to 0
             PFT_weight(p,i,j) = 1.0_RP
          end if
       end if
    end do
    end do

    deallocate( LANDUSE, LATH, LONH, YH, XH )

    return
  end subroutine CNVLANDUSE_LU100M

  !-----------------------------------------------------------------------------
  !> Convert from Japan Integrated Biodiversity Information System database 100m mesh
  subroutine CNVLANDUSE_JIBIS( PFT_weight )
    use scale_landuse, only: &
       LANDUSE_PFT_nmax
    implicit none

    real(RP), intent(inout) :: PFT_weight(-2:LANDUSE_PFT_nmax,IA,JA)
    !---------------------------------------------------------------------------

    return
  end subroutine CNVLANDUSE_JIBIS

end module mod_cnvlanduse
