!-------------------------------------------------------------------------------
!> module Atmosphere GRID CartesC Real(real space)
!!
!! @par Description
!!          Grid module for orthogonal curvelinear, terrain-following coordinate
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_grid_cartesC_real
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
  public :: ATMOS_GRID_CARTESC_REAL_setup
  public :: ATMOS_GRID_CARTESC_REAL_calc_Z
  public :: ATMOS_GRID_CARTESC_REAL_calc_areavol
  public :: ATMOS_GRID_CARTESC_REAL_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON  !< position of base point in real world [rad,0-2pi]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT  !< position of base point in real world [rad,-pi,pi]

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_CZ  (:,:,:)    !< geopotential height [m] (zxy)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_CZUY(:,:,:)    !< geopotential height [m] (zuy)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_CZXV(:,:,:)    !< geopotential height [m] (zxv)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_CZUV(:,:,:)    !< geopotential height [m] (zuv)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_FZ  (:,:,:)    !< geopotential height [m] (wxy)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_FZUY(:,:,:)    !< geopotential height [m] (wuy)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_FZXV(:,:,:)    !< geopotential height [m] (wxv)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_FZUV(:,:,:)    !< geopotential height [m] (wuv)
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_F2H (:,:,:,:)  !< coefficient for interpolation from full to half levels

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LON  (:,:)     !< longitude [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LONUY(:,:)     !< longitude at staggered point (uy) [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LONXV(:,:)     !< longitude at staggered point (xv) [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LONUV(:,:)     !< longitude at staggered point (uv) [rad,0-2pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LAT  (:,:)     !< latitude  [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LATUY(:,:)     !< latitude  at staggered point (uy) [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LATXV(:,:)     !< latitude  at staggered point (xv) [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_LATUV(:,:)     !< latitude  at staggered point (uv) [rad,-pi,pi]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_DLON (:,:)     !< delta longitude
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_DLAT (:,:)     !< delta latitude

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_Z1  (:,:)      !< Height of the lowermost grid from surface (cell center) [m]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_ASPECT_MAX     !< maximum aspect ratio of the grid cell
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_ASPECT_MIN     !< minimum aspect ratio of the grid cell

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_PHI (:,:,:)    !< geopotential [m2/s2] (cell center)

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREA     (:,:)   !< horizontal area ( xy, normal z) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAZUY_X(:,:,:) !< virtical   area (zuy, normal x) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAZXV_Y(:,:,:) !< virtical   area (zxv, normal y) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAWUY_X(:,:,:) !< virtical   area (wuy, normal x) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAWXV_Y(:,:,:) !< virtical   area (wxv, normal y) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAUY   (:,:)   !< horizontal area ( uy, normal z) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAZXY_X(:,:,:) !< virtical   area (zxy, normal x) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAZUV_Y(:,:,:) !< virtical   area (zuv, normal y) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAXV   (:,:)   !< horizontal area ( xv, normal z) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAZUV_X(:,:,:) !< virtical   area (zuv, normal x) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_AREAZXY_Y(:,:,:) !< virtical   area (zxy, normal y) [m2]

  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTAREA         !< total area (xy, local) [m2]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTAREAUY       !< total area (uy, local) [m2]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTAREAXV       !< total area (xv, local) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(:) !< total area (zuy, normal x) [m2]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y(:) !< total area (zxv, normal y) [m2]

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_VOL   (:,:,:)  !< control volume (zxy) [m3]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_VOLWXY(:,:,:)  !< control volume (wxy) [m3]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_VOLZUY(:,:,:)  !< control volume (zuy) [m3]
  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_VOLZXV(:,:,:)  !< control volume (zxv) [m3]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTVOL         !< total volume (zxy, local) [m3]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTVOLWXY      !< total volume (wxy, local) [m3]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTVOLZUY      !< total volume (zuy, local) [m3]
  real(RP), public              :: ATMOS_GRID_CARTESC_REAL_TOTVOLZXV      !< total volume (zxv, local) [m3]

  real(RP), public, allocatable :: ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(:,:,:) !< domain latlon catalogue [rad]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_GRID_CARTESC_REAL_calc_latlon

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_GRID_CARTESC_REAL_setup
    use scale_prc, only: &
       PRC_nprocs,  &
       PRC_abort
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, &
       ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
       ATMOS_GRID_CARTESC_CZ,              &
       ATMOS_GRID_CARTESC_FZ
    use scale_topography, only: &
       TOPOGRAPHY_exist
    use scale_mapprojection, only: &
       MAPPROJECTION_setup
    use scale_interp_vert, only: &
       INTERP_VERT_setcoef
    implicit none

    character(len=H_LONG) :: DOMAIN_CATALOGUE_FNAME  = 'latlon_domain_catalogue.txt' !< metadata files for lat-lon domain for all processes
    logical               :: DOMAIN_CATALOGUE_OUTPUT = .false.

    namelist / PARAM_DOMAIN_CATALOGUE / &
       DOMAIN_CATALOGUE_FNAME,  &
       DOMAIN_CATALOGUE_OUTPUT

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_REAL_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_DOMAIN_CATALOGUE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_GRID_CARTESC_REAL_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_GRID_CARTESC_REAL_setup",*) 'Not appropriate names in namelist PARAM_DOMAIN_CATALOGUE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_DOMAIN_CATALOGUE)

    allocate( ATMOS_GRID_CARTESC_REAL_LON  (  IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LAT  (  IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LONUY(0:IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LONXV(  IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LONUV(0:IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LATUY(0:IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LATXV(  IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_LATUV(0:IA,0:JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_DLON (  IA,  JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_DLAT (  IA,  JA) )
    !$acc enter data create(ATMOS_GRID_CARTESC_REAL_LON, ATMOS_GRID_CARTESC_REAL_LAT, ATMOS_GRID_CARTESC_REAL_LONUY, ATMOS_GRID_CARTESC_REAL_LONXV, ATMOS_GRID_CARTESC_REAL_LONUV, ATMOS_GRID_CARTESC_REAL_LATUY, ATMOS_GRID_CARTESC_REAL_LATXV, ATMOS_GRID_CARTESC_REAL_LATUV, ATMOS_GRID_CARTESC_REAL_DLON, ATMOS_GRID_CARTESC_REAL_DLAT)

    allocate( ATMOS_GRID_CARTESC_REAL_CZ  (  KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_CZUY(  KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_CZXV(  KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_CZUV(  KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_FZ  (0:KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_FZUY(0:KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_FZXV(0:KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_FZUV(0:KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_F2H (KA,2,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_Z1 (     IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_PHI(  KA,IA,JA) )
    !$acc enter data create(ATMOS_GRID_CARTESC_REAL_CZ, ATMOS_GRID_CARTESC_REAL_CZUY, ATMOS_GRID_CARTESC_REAL_CZXV, ATMOS_GRID_CARTESC_REAL_CZUV, ATMOS_GRID_CARTESC_REAL_FZ, ATMOS_GRID_CARTESC_REAL_FZUY, ATMOS_GRID_CARTESC_REAL_FZXV, ATMOS_GRID_CARTESC_REAL_FZUV, ATMOS_GRID_CARTESC_REAL_F2H, ATMOS_GRID_CARTESC_REAL_Z1, ATMOS_GRID_CARTESC_REAL_PHI)

    allocate( ATMOS_GRID_CARTESC_REAL_AREA     (     IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAZUY_X(KA  ,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAZXV_Y(KA  ,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAWUY_X(KA+1,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAWXV_Y(KA+1,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAUY   (     IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAZXY_X(KA,  IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAZUV_Y(KA,  IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAXV   (     IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAZUV_X(KA,  IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_AREAZXY_Y(KA,  IA,JA) )
    !$acc enter data create(ATMOS_GRID_CARTESC_REAL_AREA, ATMOS_GRID_CARTESC_REAL_AREAZUY_X, ATMOS_GRID_CARTESC_REAL_AREAZXV_Y, ATMOS_GRID_CARTESC_REAL_AREAWUY_X, ATMOS_GRID_CARTESC_REAL_AREAWXV_Y, ATMOS_GRID_CARTESC_REAL_AREAUY, ATMOS_GRID_CARTESC_REAL_AREAZXY_X, ATMOS_GRID_CARTESC_REAL_AREAZUV_Y, ATMOS_GRID_CARTESC_REAL_AREAXV, ATMOS_GRID_CARTESC_REAL_AREAZUV_X, ATMOS_GRID_CARTESC_REAL_AREAZXY_Y)

    allocate( ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(IA) )
    allocate( ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y(JA) )
    !$acc enter data create(ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X, ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y)

    allocate( ATMOS_GRID_CARTESC_REAL_VOL   (  KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_VOLWXY(0:KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_VOLZUY(  KA,IA,JA) )
    allocate( ATMOS_GRID_CARTESC_REAL_VOLZXV(  KA,IA,JA) )
    !$acc enter data create(ATMOS_GRID_CARTESC_REAL_VOL, ATMOS_GRID_CARTESC_REAL_VOLWXY, ATMOS_GRID_CARTESC_REAL_VOLZUY, ATMOS_GRID_CARTESC_REAL_VOLZXV)

    allocate( ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(PRC_nprocs,2,2) )

    ! setup map projection
    call MAPPROJECTION_setup( ATMOS_GRID_CARTESC_DOMAIN_CENTER_X, ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y )

    ! calc longitude & latitude
    call ATMOS_GRID_CARTESC_REAL_calc_latlon( DOMAIN_CATALOGUE_FNAME, DOMAIN_CATALOGUE_OUTPUT )

    ! calc real height
    call ATMOS_GRID_CARTESC_REAL_calc_Z

    call INTERP_VERT_setcoef( KA, KS, KE, IA, 1,  IA, JA, 1,  JA, & ! [IN]
                              TOPOGRAPHY_exist,                  & ! [IN]
                              ATMOS_GRID_CARTESC_CZ(:),          & ! [IN]
                              ATMOS_GRID_CARTESC_FZ(:),          & ! [IN]
                              ATMOS_GRID_CARTESC_REAL_CZ(:,:,:), & ! [IN]
                              ATMOS_GRID_CARTESC_REAL_FZ(:,:,:)  ) ! [IN]

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_setup

  !-----------------------------------------------------------------------------
  !> Calc longitude & latitude
  subroutine ATMOS_GRID_CARTESC_REAL_calc_latlon( &
       catalogue_fname, &
       catalogue_output )
    use scale_prc, only: &
       PRC_abort, &
       PRC_nprocs,  &
       PRC_IsMaster
    use scale_const, only: &
       PI  => CONST_PI, &
       D2R => CONST_D2R
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CX, &
       ATMOS_GRID_CARTESC_CY, &
       ATMOS_GRID_CARTESC_FX, &
       ATMOS_GRID_CARTESC_FY
    use scale_comm_cartesC, only: &
       COMM_gather,  &
       COMM_bcast
    use scale_mapprojection, only: &
       MAPPROJECTION_basepoint_lon, &
       MAPPROJECTION_basepoint_lat, &
       MAPPROJECTION_xy2lonlat
    implicit none

    character(len=*), intent(in) :: catalogue_fname  !< metadata files for lat-lon domain for all processes
    logical,          intent(in) :: catalogue_output

    integer, parameter :: I_MIN  = 1
    integer, parameter :: I_MAX  = 2
    integer, parameter :: I_LON  = 1
    integer, parameter :: I_LAT  = 2

    real(RP) :: mine (2,2)            !< send    buffer of lon-lat [deg]
    real(RP) :: whole(2,2,PRC_nprocs) !< recieve buffer of lon-lat [deg]

    character(len=H_LONG) :: fname

    integer  :: i, j
    integer  :: fid, ierr
    !---------------------------------------------------------------------------

    ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON = MAPPROJECTION_basepoint_lon * D2R
    ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT = MAPPROJECTION_basepoint_lat * D2R

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_REAL_calc_latlon",*) 'Base position in the global domain (lat,lon)'
    LOG_INFO_CONT(*) '-> (',ATMOS_GRID_CARTESC_REAL_BASEPOINT_LON/D2R,',',ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT/D2R,')'

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
       call MAPPROJECTION_xy2lonlat( ATMOS_GRID_CARTESC_CX(i), ATMOS_GRID_CARTESC_CY(j), ATMOS_GRID_CARTESC_REAL_LON  (i,j), ATMOS_GRID_CARTESC_REAL_LAT  (i,j) )
    enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_LON, ATMOS_GRID_CARTESC_REAL_LAT) async

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 0, IA
       call MAPPROJECTION_xy2lonlat( ATMOS_GRID_CARTESC_FX(i), ATMOS_GRID_CARTESC_CY(j), ATMOS_GRID_CARTESC_REAL_LONUY(i,j), ATMOS_GRID_CARTESC_REAL_LATUY(i,j) )
    enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_LONUY, ATMOS_GRID_CARTESC_REAL_LATUY) async

    !$omp parallel do collapse(2)
    do j = 0, JA
    do i = 1, IA
       call MAPPROJECTION_xy2lonlat( ATMOS_GRID_CARTESC_CX(i), ATMOS_GRID_CARTESC_FY(j), ATMOS_GRID_CARTESC_REAL_LONXV(i,j), ATMOS_GRID_CARTESC_REAL_LATXV(i,j) )
    enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_LONXV, ATMOS_GRID_CARTESC_REAL_LATXV) async

    !$omp parallel do collapse(2)
    do j = 0, JA
    do i = 0, IA
       call MAPPROJECTION_xy2lonlat( ATMOS_GRID_CARTESC_FX(i), ATMOS_GRID_CARTESC_FY(j), ATMOS_GRID_CARTESC_REAL_LONUV(i,j), ATMOS_GRID_CARTESC_REAL_LATUV(i,j) )
    enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_LONUV, ATMOS_GRID_CARTESC_REAL_LATUV) async

    !$omp workshare
!OCL ZFILL
    ATMOS_GRID_CARTESC_REAL_DLON(:,:) = 0.0_RP
!OCL ZFILL
    ATMOS_GRID_CARTESC_REAL_DLAT(:,:) = 0.0_RP
    !$omp end workshare

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       ATMOS_GRID_CARTESC_REAL_DLON(i,j) = abs( ATMOS_GRID_CARTESC_REAL_LONUY(i,j) - ATMOS_GRID_CARTESC_REAL_LONUY(i-1,j) )
       ATMOS_GRID_CARTESC_REAL_DLAT(i,j) = abs( ATMOS_GRID_CARTESC_REAL_LATXV(i,j) - ATMOS_GRID_CARTESC_REAL_LATXV(i,j-1) )
       if( ATMOS_GRID_CARTESC_REAL_DLON(i,j) > 2.0_RP*PI - ATMOS_GRID_CARTESC_REAL_DLON(i,j)  ) ATMOS_GRID_CARTESC_REAL_DLON(i,j) = 2.0_RP*PI - ATMOS_GRID_CARTESC_REAL_DLON(i,j)

       if (      ATMOS_GRID_CARTESC_REAL_DLON(i,j) == 0.0_RP &
            .OR. ATMOS_GRID_CARTESC_REAL_DLAT(i,j) == 0.0_RP ) then
          LOG_ERROR("ATMOS_GRID_CARTESC_REAL_calc_latlon",*) 'Invalid grid distance in lat-lon! i,j=', i,j
          LOG_ERROR_CONT(*) 'Lon(i-1),Lon(i),dlon=', ATMOS_GRID_CARTESC_REAL_LONUY(i-1,j)/D2R,ATMOS_GRID_CARTESC_REAL_LONUY(i,j)/D2R,ATMOS_GRID_CARTESC_REAL_DLON(i,j)/D2R
          LOG_ERROR_CONT(*) 'Lat(j-1),Lat(j),dlat=', ATMOS_GRID_CARTESC_REAL_LATXV(i,j-1)/D2R,ATMOS_GRID_CARTESC_REAL_LATXV(i,j)/D2R,ATMOS_GRID_CARTESC_REAL_DLAT(i,j)/D2R
          call PRC_abort
       endif
    enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_DLON, ATMOS_GRID_CARTESC_REAL_DLAT) async

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_REAL_calc_latlon",*) 'Position on the earth (Local)'
    LOG_INFO_CONT('(1x,A,F10.5,A,F9.5,A,A,F10.5,A,F9.5,A)') &
                               'NW(',ATMOS_GRID_CARTESC_REAL_LON(IS,JE)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IS,JE)/D2R,')', &
                            ' - NE(',ATMOS_GRID_CARTESC_REAL_LON(IE,JE)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IE,JE)/D2R,')'

    LOG_INFO_CONT('(1x,A)') '             |                          |'
    LOG_INFO_CONT('(1x,A,F10.5,A,F9.5,A,A,F10.5,A,F9.5,A)') &
                               'SW(',ATMOS_GRID_CARTESC_REAL_LON(IS,JS)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IS,JS)/D2R,')', &
                            ' - SE(',ATMOS_GRID_CARTESC_REAL_LON(IE,JS)/D2R,',',ATMOS_GRID_CARTESC_REAL_LAT(IE,JS)/D2R,')'

    mine(I_MIN,I_LON) = minval(ATMOS_GRID_CARTESC_REAL_LONUV(IS-1:IE,JS-1:JE)) / D2R
    mine(I_MAX,I_LON) = maxval(ATMOS_GRID_CARTESC_REAL_LONUV(IS-1:IE,JS-1:JE)) / D2R
    mine(I_MIN,I_LAT) = minval(ATMOS_GRID_CARTESC_REAL_LATUV(IS-1:IE,JS-1:JE)) / D2R
    mine(I_MAX,I_LAT) = maxval(ATMOS_GRID_CARTESC_REAL_LATUV(IS-1:IE,JS-1:JE)) / D2R

    call COMM_gather( 2, 2, mine(:,:), whole(:,:,:)  ) ! everytime do for online nesting

    if ( PRC_IsMaster ) then
       if ( catalogue_output ) then

          fid = IO_get_available_fid()
          call IO_get_fname(fname, catalogue_fname)
          open( fid,                  &
                file   = fname,       &
                form   = 'formatted', &
                status = 'replace',   &
                iostat = ierr         )

          if ( ierr /= 0 ) then
             LOG_ERROR("ATMOS_GRID_CARTESC_REAL_calc_latlon",*) 'cannot create latlon-catalogue file!'
             call PRC_abort
          endif

          do i = 1, PRC_nprocs ! for offline nesting
             write(fid,'(I8,8F32.24)',iostat=ierr) i, whole(I_MIN,I_LON,i), whole(I_MAX,I_LON,i), & ! LON: MIN, MAX
                                                      whole(I_MIN,I_LAT,i), whole(I_MAX,I_LAT,i)    ! LAT: MIN, MAX
             if ( ierr /= 0 ) exit
          enddo

          close(fid)

       endif

       do i = 1, PRC_nprocs ! for online nesting
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_MIN,I_LON) = whole(I_MIN,I_LON,i)
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_MAX,I_LON) = whole(I_MAX,I_LON,i)
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_MIN,I_LAT) = whole(I_MIN,I_LAT,i)
          ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(i,I_MAX,I_LAT) = whole(I_MAX,I_LAT,i)
       enddo
    endif

    call COMM_bcast( PRC_nprocs, 2, 2, ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE(:,:,:) )

    !$acc wait

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_calc_latlon

  !-----------------------------------------------------------------------------
  !> Convert Xi to Z coordinate
  subroutine ATMOS_GRID_CARTESC_REAL_calc_Z
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CZ,  &
       ATMOS_GRID_CARTESC_FZ,  &
       ATMOS_GRID_CARTESC_CDX, &
       ATMOS_GRID_CARTESC_CDY
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_coordinates_atmos
    use scale_topography, only: &
       Zsfc => TOPOGRAPHY_Zsfc
    use scale_landuse, only: &
       LANDUSE_frac_land
    implicit none

    real(DP) :: Htop
    real(RP) :: Zs
    real(RP) :: DFZ

    real(RP) :: dz1, dz2

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    Htop = ATMOS_GRID_CARTESC_FZ(KE) - ATMOS_GRID_CARTESC_FZ(KS-1)

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA
    do i = 1, IA
       Zs = Zsfc(i,j)
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZ(k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_CZ) async

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA
    do i = 1, IA-1
       Zs = ( Zsfc(i,j) + Zsfc(i+1,j) ) * 0.5_RP
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZUY(k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    enddo
    !$omp parallel do private(zs)
    do j = 1, JA
       Zs = Zsfc(IA,j)
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZUY(k,IA,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_CZUY) async

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA-1
    do i = 1, IA
       Zs = ( Zsfc(i,j) + Zsfc(i,j+1) ) * 0.5_RP
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZXV(k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    enddo
    do i = 1, IA
       Zs = Zsfc(i,JA)
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZXV(k,i,JA) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_CZXV) async

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA-1
    do i = 1, IA-1
       Zs = ( Zsfc(i,j) + Zsfc(i+1,j) + Zsfc(i,j+1) + Zsfc(i+1,j+1) ) * 0.25_RP
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZUV(k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    enddo
    !$omp parallel do private(zs)
    do j = 1, JA-1
       Zs = ( Zsfc(IA,j) + Zsfc(IA,j+1) ) * 0.5_RP
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZUV(k,IA,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    do i = 1, IA-1
       Zs = ( Zsfc(i,JA) + Zsfc(i+1,JA) ) * 0.5_RP
       do k = 1, KA
          ATMOS_GRID_CARTESC_REAL_CZUV(k,i,JA) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
       enddo
    enddo
    Zs = Zsfc(IA,JA)
    do k = 1, KA
       ATMOS_GRID_CARTESC_REAL_CZUV(k,IA,JA) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_CZ(k) + Zs
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_CZUV) async

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA
    do i = 1, IA
       Zs = Zsfc(i,j)
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZ (k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       end do
    end do
    end do
    !$acc update device(ATMOS_GRID_CARTESC_REAL_FZ) async

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA
    do i = 1, IA-1
       Zs = ( Zsfc(i,j) + Zsfc(i+1,j) ) * 0.5_RP
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZUY(k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       end do
    end do
    end do
    !$omp parallel do private(zs)
    do j = 1, JA
       Zs = Zsfc(IA,j)
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZUY(k,IA,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       end do
    end do
    !$acc update device(ATMOS_GRID_CARTESC_REAL_FZUY) async

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA-1
    do i = 1, IA
       Zs = ( Zsfc(i,j) + Zsfc(i,j+1) ) * 0.5_RP
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZXV(k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       enddo
    enddo
    enddo
    !$omp parallel do private(zs)
    do i = 1, IA
       Zs = Zsfc(i,JA)
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZXV(k,i,JA) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       enddo
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_FZXV) async

    !$omp parallel do private(zs) collapse(2)
    do j = 1, JA-1
    do i = 1, IA-1
       Zs = ( Zsfc(i,j) + Zsfc(i+1,j) + Zsfc(i,j+1) + Zsfc(i+1,j+1) ) * 0.25_RP
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZUV(k,i,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       enddo
    enddo
    enddo
    !$omp parallel do private(zs)
    do j = 1, JA-1
       Zs = ( Zsfc(IA,j) + Zsfc(IA,j+1) ) * 0.5_RP
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZUV(k,IA,j) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       enddo
    enddo
    !$omp parallel do private(zs)
    do i = 1, IA-1
       Zs = ( Zsfc(i,JA) + Zsfc(i+1,JA) ) * 0.5_RP
       do k = 0, KA
          ATMOS_GRID_CARTESC_REAL_FZUV(k,i,JA) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
       enddo
    enddo
    Zs = Zsfc(IA,JA)
    do k = 0, KA
       ATMOS_GRID_CARTESC_REAL_FZUV(k,IA,JA) = ( Htop - Zs ) / Htop * ATMOS_GRID_CARTESC_FZ(k) + Zs
    enddo
    !$acc update device(ATMOS_GRID_CARTESC_REAL_FZUV) async

    !$omp parallel do private(dz1,dz2) collapse(2)
    do j = 1, JA
    do i = 1, IA
       do k = KS, KE-1
          dz1 = ATMOS_GRID_CARTESC_REAL_FZ(k+1,i,j) - ATMOS_GRID_CARTESC_REAL_FZ(k  ,i,j)
          dz2 = ATMOS_GRID_CARTESC_REAL_FZ(k  ,i,j) - ATMOS_GRID_CARTESC_REAL_FZ(k-1,i,j)
          ATMOS_GRID_CARTESC_REAL_F2H(k,1,i,j) = dz2 / ( dz1 + dz2 )
          ATMOS_GRID_CARTESC_REAL_F2H(k,2,i,j) = dz1 / ( dz1 + dz2 )
       end do
       ATMOS_GRID_CARTESC_REAL_F2H(1:KS-1,1,i,j) = 0.5_RP
       ATMOS_GRID_CARTESC_REAL_F2H(1:KS-1,2,i,j) = 0.5_RP
       ATMOS_GRID_CARTESC_REAL_F2H(KE:KA ,1,i,j) = 0.5_RP
       ATMOS_GRID_CARTESC_REAL_F2H(KE:KA ,2,i,j) = 0.5_RP
    end do
    end do
    !$acc update device(ATMOS_GRID_CARTESC_REAL_F2H) async

    !$omp workshare
!OCL ZFILL
    ATMOS_GRID_CARTESC_REAL_Z1(:,:) = ATMOS_GRID_CARTESC_REAL_CZ(KS,:,:) - Zsfc(:,:)
    !$acc update device(ATMOS_GRID_CARTESC_REAL_Z1) async

!OCL ZFILL
    ATMOS_GRID_CARTESC_REAL_PHI(:,:,:) = GRAV * ATMOS_GRID_CARTESC_REAL_CZ(:,:,:)
    !$omp end workshare
    !$acc update device(ATMOS_GRID_CARTESC_REAL_Z1, ATMOS_GRID_CARTESC_REAL_PHI) async

    ATMOS_GRID_CARTESC_REAL_ASPECT_MAX = -1.E+30_RP
    ATMOS_GRID_CARTESC_REAL_ASPECT_MIN =  1.E+30_RP

    !$omp parallel do private(dfz) collapse(2) &
    !$omp reduction(max:ATMOS_GRID_CARTESC_REAL_ASPECT_MAX) &
    !$omp reduction(min:ATMOS_GRID_CARTESC_REAL_ASPECT_MIN)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DFZ = ATMOS_GRID_CARTESC_REAL_FZ(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZ(k-1,i,j)
       ATMOS_GRID_CARTESC_REAL_ASPECT_MAX = max( ATMOS_GRID_CARTESC_REAL_ASPECT_MAX, ATMOS_GRID_CARTESC_CDX(i) / DFZ, ATMOS_GRID_CARTESC_CDY(j) / DFZ )
       ATMOS_GRID_CARTESC_REAL_ASPECT_MIN = min( ATMOS_GRID_CARTESC_REAL_ASPECT_MIN, ATMOS_GRID_CARTESC_CDX(i) / DFZ, ATMOS_GRID_CARTESC_CDY(j) / DFZ )
    enddo
    enddo
    enddo

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_CARTESC_REAL_calc_Z",*) 'Minimum & maximum lowermost CZ'
    LOG_INFO_CONT(*) '-> (',minval( ATMOS_GRID_CARTESC_REAL_CZ(KS,:,:) ),',',maxval( ATMOS_GRID_CARTESC_REAL_CZ(KS,:,:) ),')'
    LOG_INFO("ATMOS_GRID_CARTESC_REAL_calc_Z",*) 'Minimum & maximum aspect ratio'
    LOG_INFO_CONT(*) '-> (',ATMOS_GRID_CARTESC_REAL_ASPECT_MIN,',',ATMOS_GRID_CARTESC_REAL_ASPECT_MAX,')'

    ! set latlon and z to fileio module
    call FILE_CARTESC_set_coordinates_atmos( ATMOS_GRID_CARTESC_REAL_CZ,  ATMOS_GRID_CARTESC_REAL_FZ,                                                                      & ! [IN]
                                             ATMOS_GRID_CARTESC_REAL_LON, ATMOS_GRID_CARTESC_REAL_LONUY, ATMOS_GRID_CARTESC_REAL_LONXV, ATMOS_GRID_CARTESC_REAL_LONUV,     & ! [IN]
                                             ATMOS_GRID_CARTESC_REAL_LAT, ATMOS_GRID_CARTESC_REAL_LATUY, ATMOS_GRID_CARTESC_REAL_LATXV, ATMOS_GRID_CARTESC_REAL_LATUV,     & ! [IN]
                                             Zsfc, LANDUSE_frac_land                                                                                                      ) ! [IN]

    !$acc wait

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_calc_Z

  !-----------------------------------------------------------------------------
  !> Calc control area/volume
  subroutine ATMOS_GRID_CARTESC_REAL_calc_areavol( &
       MAPF )
    use scale_prc_cartesC, only: &
       PRC_TwoD
    use scale_const, &
       UNDEF => CONST_UNDEF
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CDX, &
       ATMOS_GRID_CARTESC_FDX, &
       ATMOS_GRID_CARTESC_CDY, &
       ATMOS_GRID_CARTESC_FDY
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_file_cartesC, only: &
       FILE_CARTESC_set_areavol_atmos
    use scale_topography, only: &
       TOPOGRAPHY_Zsfc
    use scale_landuse, only: &
       LANDUSE_frac_land
    implicit none

    real(RP), intent(in) :: MAPF(IA,JA,2,4)

    real(RP) :: AREAUV(IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ATMOS_GRID_CARTESC_REAL_AREA     (:,:)   = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAZUY_X(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAZXV_Y(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAWUY_X(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAWXV_Y(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAUY   (:,:)   = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAZXY_X(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAZUV_Y(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAXV   (:,:)   = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAZUV_Y(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_AREAZXY_Y(:,:,:) = 0.0_RP

    ATMOS_GRID_CARTESC_REAL_TOTAREA         = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTAREAUY       = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTAREAXV       = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y(:) = 0.0_RP

    ATMOS_GRID_CARTESC_REAL_VOL   (:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_VOLWXY(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_VOLZUY(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_VOLZXV(:,:,:) = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTVOL    = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTVOLWXY = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTVOLZUY = 0.0_RP
    ATMOS_GRID_CARTESC_REAL_TOTVOLZXV = 0.0_RP

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       ATMOS_GRID_CARTESC_REAL_AREA  (i,j) = ATMOS_GRID_CARTESC_CDX(i) * ATMOS_GRID_CARTESC_CDY(j) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
       ATMOS_GRID_CARTESC_REAL_AREAXV(i,j) = ATMOS_GRID_CARTESC_CDX(i) * ATMOS_GRID_CARTESC_FDY(j) / ( MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) )
    end do
    end do
    if ( PRC_TwoD ) then
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
          ATMOS_GRID_CARTESC_REAL_AREAUY(i,j) = ATMOS_GRID_CARTESC_REAL_AREA  (i,j)
                                  AREAUV(i,j) = ATMOS_GRID_CARTESC_REAL_AREAXV(i,j)

       end do
       end do
    else
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
          ATMOS_GRID_CARTESC_REAL_AREAUY(i,j) = ATMOS_GRID_CARTESC_FDX(i) * ATMOS_GRID_CARTESC_CDY(j) / ( MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) )
                                  AREAUV(i,j) = ATMOS_GRID_CARTESC_FDX(i) * ATMOS_GRID_CARTESC_FDY(j) / ( MAPF(i,j,1,I_UV) * MAPF(i,j,2,I_UV) )

       end do
       end do
    end if

#ifdef QUICKDEBUG
    ATMOS_GRID_CARTESC_REAL_AREA  (1:IS-1,:)  = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREA  (IE+1:IA,:) = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREA  (:,1:JS-1)  = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREA  (:,JE+1:JA) = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAXV(1:IS-1,:)  = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAXV(IE+1:IA,:) = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAXV(:,1:JS-1)  = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAXV(:,JE+1:JA) = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAUY(1:IS-1,:)  = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAUY(IE+1:IA,:) = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAUY(:,1:JS-1)  = UNDEF
    ATMOS_GRID_CARTESC_REAL_AREAUY(:,JE+1:JA) = UNDEF
                            AREAUV(1:IS-1,:)  = UNDEF
                            AREAUV(IE+1:IA,:) = UNDEF
                            AREAUV(:,1:JS-1)  = UNDEF
                            AREAUV(:,JE+1:JA) = UNDEF
#endif
    call COMM_vars8( ATMOS_GRID_CARTESC_REAL_AREA  (:,:), 1 )
    call COMM_vars8( ATMOS_GRID_CARTESC_REAL_AREAXV(:,:), 2 )
    call COMM_vars8( ATMOS_GRID_CARTESC_REAL_AREAUY(:,:), 3 )
    call COMM_vars8(                         AREAUV(:,:), 4 )

    !$omp parallel do &
    !$omp reduction(+:ATMOS_GRID_CARTESC_REAL_TOTAREA,ATMOS_GRID_CARTESC_REAL_TOTAREAXV,ATMOS_GRID_CARTESC_REAL_TOTAREAUY)
    do j = JS, JE
    do i = IS, IE
       ATMOS_GRID_CARTESC_REAL_TOTAREA   = ATMOS_GRID_CARTESC_REAL_TOTAREA   + ATMOS_GRID_CARTESC_REAL_AREA  (i,j)
       ATMOS_GRID_CARTESC_REAL_TOTAREAXV = ATMOS_GRID_CARTESC_REAL_TOTAREAXV + ATMOS_GRID_CARTESC_REAL_AREAXV(i,j)
       ATMOS_GRID_CARTESC_REAL_TOTAREAUY = ATMOS_GRID_CARTESC_REAL_TOTAREAUY + ATMOS_GRID_CARTESC_REAL_AREAUY(i,j)
    end do
    end do

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
       do k = KS, KE
          ATMOS_GRID_CARTESC_REAL_AREAZXV_Y(k,i,j) = ATMOS_GRID_CARTESC_CDX(i) / MAPF(i,j,1,I_XV) * ( ATMOS_GRID_CARTESC_REAL_FZXV(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZXV(k-1,i,j) )
       end do
       do k = KS-1, KE
          ATMOS_GRID_CARTESC_REAL_AREAWXV_Y(k,i,j) = ATMOS_GRID_CARTESC_CDX(i) / MAPF(i,j,1,I_XV) * ( ATMOS_GRID_CARTESC_REAL_CZXV(k+1,i,j) - ATMOS_GRID_CARTESC_REAL_CZXV(k,i,j) )
       end do
       do k = KS, KE
          ATMOS_GRID_CARTESC_REAL_AREAZXY_X(k,i,j) = ATMOS_GRID_CARTESC_CDY(j) / MAPF(i,j,2,I_XY) * ( ATMOS_GRID_CARTESC_REAL_FZ  (k,i,j) - ATMOS_GRID_CARTESC_REAL_FZ  (k-1,i,j) )
          ATMOS_GRID_CARTESC_REAL_AREAZXY_Y(k,i,j) = ATMOS_GRID_CARTESC_CDX(i) / MAPF(i,j,1,I_XY) * ( ATMOS_GRID_CARTESC_REAL_FZ  (k,i,j) - ATMOS_GRID_CARTESC_REAL_FZ  (k-1,i,j) )
       end do
    end do
    end do
    !$acc update device(ATMOS_GRID_CARTESC_REAL_AREAZXV_Y, ATMOS_GRID_CARTESC_REAL_AREAWXV_Y, ATMOS_GRID_CARTESC_REAL_AREAZXY_X, ATMOS_GRID_CARTESC_REAL_AREAZXY_Y) async

    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
       do k = KS, KE
          ATMOS_GRID_CARTESC_REAL_AREAZUY_X(k,i,j) = ATMOS_GRID_CARTESC_CDY(j) / MAPF(i,j,2,I_UY) * ( ATMOS_GRID_CARTESC_REAL_FZUY(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUY(k-1,i,j) )
       end do
       do k = KS-1, KE
          ATMOS_GRID_CARTESC_REAL_AREAWUY_X(k,i,j) = ATMOS_GRID_CARTESC_CDY(j) / MAPF(i,j,2,I_UY) * ( ATMOS_GRID_CARTESC_REAL_CZUY(k+1,i,j) - ATMOS_GRID_CARTESC_REAL_CZUY(k,i,j) )
       end do
       do k = KS, KE
          ATMOS_GRID_CARTESC_REAL_AREAZUV_Y(k,i,j) = ATMOS_GRID_CARTESC_CDX(i) / MAPF(i,j,1,I_UV) * ( ATMOS_GRID_CARTESC_REAL_FZUV(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUV(k-1,i,j) )
          ATMOS_GRID_CARTESC_REAL_AREAZUV_X(k,i,j) = ATMOS_GRID_CARTESC_CDY(j) / MAPF(i,j,2,I_UV) * ( ATMOS_GRID_CARTESC_REAL_FZUV(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUV(k-1,i,j) )
       end do
    end do
    end do
    !$acc update device(ATMOS_GRID_CARTESC_REAL_AREAZUY_X, ATMOS_GRID_CARTESC_REAL_AREAWUY_X, ATMOS_GRID_CARTESC_REAL_AREAZUV_Y, ATMOS_GRID_CARTESC_REAL_AREAZUV_X) async

    call COMM_wait( ATMOS_GRID_CARTESC_REAL_AREA  (:,:), 1 )
    call COMM_wait( ATMOS_GRID_CARTESC_REAL_AREAXV(:,:), 2 )
    call COMM_wait( ATMOS_GRID_CARTESC_REAL_AREAUY(:,:), 3 )
    call COMM_wait(                         AREAUV(:,:), 4 )
    !$acc update device(ATMOS_GRID_CARTESC_REAL_AREA, ATMOS_GRID_CARTESC_REAL_AREAXV, ATMOS_GRID_CARTESC_REAL_AREAUY) async


    !$omp parallel do collapse(2)
    do j = 1,  JA
    do i = IS, IE
    do k = KS, KE
       ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y(j) = ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y(j) + ATMOS_GRID_CARTESC_REAL_AREAZXV_Y(k,i,j)
    end do
    end do
    end do
    !$omp parallel do collapse(2)
    do j = JS, JE
    do i = 1,  IA
    do k = KS, KE
       ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(i) = ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X(i) + ATMOS_GRID_CARTESC_REAL_AREAZUY_X(k,i,j)
    end do
    end do
    end do
    !$acc update device(ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y, ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X) async


    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
       do k = KS, KE
          ATMOS_GRID_CARTESC_REAL_VOL   (k,i,j) = ( ATMOS_GRID_CARTESC_REAL_FZ  (k,i,j) - ATMOS_GRID_CARTESC_REAL_FZ  (k-1,i,j) ) * ATMOS_GRID_CARTESC_REAL_AREA (i,j)
          ATMOS_GRID_CARTESC_REAL_VOLZXV(k,i,j) = ( ATMOS_GRID_CARTESC_REAL_FZXV(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZXV(k-1,i,j) ) * ATMOS_GRID_CARTESC_REAL_AREAXV(i,j)
       end do
       do k = KS-1, KE
          ATMOS_GRID_CARTESC_REAL_VOLWXY(k,i,j) = ( ATMOS_GRID_CARTESC_REAL_CZ(k+1,i,j) - ATMOS_GRID_CARTESC_REAL_CZ(k,i,j) ) * ATMOS_GRID_CARTESC_REAL_AREA(i,j)
       end do
    end do
    end do
    if ( PRC_TwoD ) then
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          ATMOS_GRID_CARTESC_REAL_VOLZUY(k,i,j) = ATMOS_GRID_CARTESC_REAL_VOL(k,i,j)
       end do
       end do
       end do
    else
       !$omp parallel do collapse(2)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE
          ATMOS_GRID_CARTESC_REAL_VOLZUY(k,i,j) = ( ATMOS_GRID_CARTESC_REAL_FZUY(k,i,j) - ATMOS_GRID_CARTESC_REAL_FZUY(k-1,i,j) ) * ATMOS_GRID_CARTESC_REAL_AREAUY(i,j)
       end do
       end do
       end do
    end if
    !$acc update device(ATMOS_GRID_CARTESC_REAL_VOL, ATMOS_GRID_CARTESC_REAL_VOLZXV, ATMOS_GRID_CARTESC_REAL_VOLWXY, ATMOS_GRID_CARTESC_REAL_VOLZUY) async

    !$omp parallel do &
    !$omp reduction(+:ATMOS_GRID_CARTESC_REAL_TOTVOL,ATMOS_GRID_CARTESC_REAL_TOTVOLZXV,ATMOS_GRID_CARTESC_REAL_TOTVOLWXY,ATMOS_GRID_CARTESC_REAL_TOTVOLZUY)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ATMOS_GRID_CARTESC_REAL_TOTVOL    = ATMOS_GRID_CARTESC_REAL_TOTVOL    + ATMOS_GRID_CARTESC_REAL_VOL   (k,i,j)
       ATMOS_GRID_CARTESC_REAL_TOTVOLZXV = ATMOS_GRID_CARTESC_REAL_TOTVOLZXV + ATMOS_GRID_CARTESC_REAL_VOLZXV(k,i,j)
       ATMOS_GRID_CARTESC_REAL_TOTVOLWXY = ATMOS_GRID_CARTESC_REAL_TOTVOLWXY + ATMOS_GRID_CARTESC_REAL_VOLWXY(k,i,j)
       ATMOS_GRID_CARTESC_REAL_TOTVOLZUY = ATMOS_GRID_CARTESC_REAL_TOTVOLZUY + ATMOS_GRID_CARTESC_REAL_VOLZUY(k,i,j)
    enddo
    enddo
    enddo

    ! set latlon and z to fileio module
    call FILE_CARTESC_set_areavol_atmos( ATMOS_GRID_CARTESC_REAL_AREA,   ATMOS_GRID_CARTESC_REAL_AREAZUY_X, ATMOS_GRID_CARTESC_REAL_AREAZXV_Y,                         & ! [IN]
                                                                         ATMOS_GRID_CARTESC_REAL_AREAWUY_X, ATMOS_GRID_CARTESC_REAL_AREAWXV_Y,                         & ! [IN]
                                         ATMOS_GRID_CARTESC_REAL_AREAUY, ATMOS_GRID_CARTESC_REAL_AREAZXY_X, ATMOS_GRID_CARTESC_REAL_AREAZUV_Y,                         & ! [IN]
                                         ATMOS_GRID_CARTESC_REAL_AREAXV, ATMOS_GRID_CARTESC_REAL_AREAZUV_X, ATMOS_GRID_CARTESC_REAL_AREAZXY_Y,                         & ! [IN]
                                         ATMOS_GRID_CARTESC_REAL_VOL, ATMOS_GRID_CARTESC_REAL_VOLWXY, ATMOS_GRID_CARTESC_REAL_VOLZUY, ATMOS_GRID_CARTESC_REAL_VOLZXV   ) ! [IN]

    !$acc wait

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_calc_areavol

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine ATMOS_GRID_CARTESC_REAL_finalize
    use scale_interp_vert, only: &
       INTERP_VERT_finalize
    implicit none
    !---------------------------------------------------------------------------

    !$acc exit data delete(ATMOS_GRID_CARTESC_REAL_LON, ATMOS_GRID_CARTESC_REAL_LAT, ATMOS_GRID_CARTESC_REAL_LONUY, ATMOS_GRID_CARTESC_REAL_LONXV, ATMOS_GRID_CARTESC_REAL_LONUV, ATMOS_GRID_CARTESC_REAL_LATUY, ATMOS_GRID_CARTESC_REAL_LATXV, ATMOS_GRID_CARTESC_REAL_LATUV, ATMOS_GRID_CARTESC_REAL_DLON, ATMOS_GRID_CARTESC_REAL_DLAT)
    deallocate( ATMOS_GRID_CARTESC_REAL_LON   )
    deallocate( ATMOS_GRID_CARTESC_REAL_LAT   )
    deallocate( ATMOS_GRID_CARTESC_REAL_LONUY )
    deallocate( ATMOS_GRID_CARTESC_REAL_LONXV )
    deallocate( ATMOS_GRID_CARTESC_REAL_LONUV )
    deallocate( ATMOS_GRID_CARTESC_REAL_LATUY )
    deallocate( ATMOS_GRID_CARTESC_REAL_LATXV )
    deallocate( ATMOS_GRID_CARTESC_REAL_LATUV )
    deallocate( ATMOS_GRID_CARTESC_REAL_DLON )
    deallocate( ATMOS_GRID_CARTESC_REAL_DLAT )

    !$acc exit data delete(ATMOS_GRID_CARTESC_REAL_CZ, ATMOS_GRID_CARTESC_REAL_CZUY, ATMOS_GRID_CARTESC_REAL_CZXV, ATMOS_GRID_CARTESC_REAL_CZUV, ATMOS_GRID_CARTESC_REAL_FZ, ATMOS_GRID_CARTESC_REAL_FZUY, ATMOS_GRID_CARTESC_REAL_FZXV, ATMOS_GRID_CARTESC_REAL_FZUV, ATMOS_GRID_CARTESC_REAL_F2H, ATMOS_GRID_CARTESC_REAL_Z1, ATMOS_GRID_CARTESC_REAL_PHI)
    deallocate( ATMOS_GRID_CARTESC_REAL_CZ   )
    deallocate( ATMOS_GRID_CARTESC_REAL_CZUY )
    deallocate( ATMOS_GRID_CARTESC_REAL_CZXV )
    deallocate( ATMOS_GRID_CARTESC_REAL_CZUV )
    deallocate( ATMOS_GRID_CARTESC_REAL_FZ   )
    deallocate( ATMOS_GRID_CARTESC_REAL_FZUY )
    deallocate( ATMOS_GRID_CARTESC_REAL_FZXV )
    deallocate( ATMOS_GRID_CARTESC_REAL_FZUV )
    deallocate( ATMOS_GRID_CARTESC_REAL_F2H  )
    deallocate( ATMOS_GRID_CARTESC_REAL_Z1  )
    deallocate( ATMOS_GRID_CARTESC_REAL_PHI )

    !acc exit data delete(ATMOS_GRID_CARTESC_REAL_AREA, ATMOS_GRID_CARTESC_REAL_AREAZUY_X, ATMOS_GRID_CARTESC_REAL_AREAZXV_Y, ATMOS_GRID_CARTESC_REAL_AREAWUY_X, ATMOS_GRID_CARTESC_REAL_AREAWXV_Y, ATMOS_GRID_CARTESC_REAL_AREAUY, ATMOS_GRID_CARTESC_REAL_AREAZXY_X, ATMOS_GRID_CARTESC_REAL_AREAZUV_Y, ATMOS_GRID_CARTESC_REAL_AREAXV, ATMOS_GRID_CARTESC_REAL_AREAZUV_X, ATMOS_GRID_CARTESC_REAL_ZXY_Y)
    deallocate( ATMOS_GRID_CARTESC_REAL_AREA      )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAZUY_X )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAZXV_Y )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAWUY_X )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAWXV_Y )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAUY    )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAZXY_X )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAZUV_Y )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAXV    )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAZUV_X )
    deallocate( ATMOS_GRID_CARTESC_REAL_AREAZXY_Y )

    !$acc exit data delete(ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X, ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y)
    deallocate( ATMOS_GRID_CARTESC_REAL_TOTAREAZUY_X )
    deallocate( ATMOS_GRID_CARTESC_REAL_TOTAREAZXV_Y )

    !$acc exit data delete(ATMOS_GRID_CARTESC_REAL_VOL, ATMOS_GRID_CARTESC_REAL_VOLWXY, ATMOS_GRID_CARTESC_REAL_VOLZUY, ATMOS_GRID_CARTESC_REAL_VOLZXV)
    deallocate( ATMOS_GRID_CARTESC_REAL_VOL    )
    deallocate( ATMOS_GRID_CARTESC_REAL_VOLWXY )
    deallocate( ATMOS_GRID_CARTESC_REAL_VOLZUY )
    deallocate( ATMOS_GRID_CARTESC_REAL_VOLZXV )

    deallocate( ATMOS_GRID_CARTESC_REAL_DOMAIN_CATALOGUE )

    call INTERP_VERT_finalize

    return
  end subroutine ATMOS_GRID_CARTESC_REAL_finalize

end module scale_atmos_grid_cartesC_real
