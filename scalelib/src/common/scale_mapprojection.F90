!-------------------------------------------------------------------------------
!> module Map projection
!!
!! @par Description
!!          Map projection module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_mapprojection
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MAPPROJECTION_setup

  public :: MAPPROJECTION_get_param
  public :: MAPPROJECTION_get_param_None
  public :: MAPPROJECTION_get_param_LambertConformal
  public :: MAPPROJECTION_get_param_PolarStereographic
  public :: MAPPROJECTION_get_param_Mercator
  public :: MAPPROJECTION_get_param_EquidistantCylindrical

  public :: MAPPROJECTION_xy2lonlat
  public :: MAPPROJECTION_xy2lonlat_None
  public :: MAPPROJECTION_xy2lonlat_LambertConformal
  public :: MAPPROJECTION_xy2lonlat_PolarStereographic
  public :: MAPPROJECTION_xy2lonlat_Mercator
  public :: MAPPROJECTION_xy2lonlat_EquidistantCylindrical

  public :: MAPPROJECTION_lonlat2xy
  public :: MAPPROJECTION_lonlat2xy_None
  public :: MAPPROJECTION_lonlat2xy_LambertConformal
  public :: MAPPROJECTION_lonlat2xy_PolarStereographic
  public :: MAPPROJECTION_lonlat2xy_Mercator
  public :: MAPPROJECTION_lonlat2xy_EquidistantCylindrical

  public :: MAPPROJECTION_mapfactor
  public :: MAPPROJECTION_mapfactor_None
  public :: MAPPROJECTION_mapfactor_LambertConformal
  public :: MAPPROJECTION_mapfactor_PolarStereographic
  public :: MAPPROJECTION_mapfactor_Mercator
  public :: MAPPROJECTION_mapfactor_EquidistantCylindrical

  public :: MAPPROJECTION_rotcoef
  public :: MAPPROJECTION_rotcoef_None
  public :: MAPPROJECTION_rotcoef_LambertConformal
  public :: MAPPROJECTION_rotcoef_PolarStereographic
  public :: MAPPROJECTION_rotcoef_Mercator
  public :: MAPPROJECTION_rotcoef_EquidistantCylindrical

  interface MAPPROJECTION_xy2lonlat
     module procedure MAPPROJECTION_xy2lonlat_0D
     module procedure MAPPROJECTION_xy2lonlat_2D
  end interface
  interface MAPPROJECTION_lonlat2xy
     module procedure MAPPROJECTION_lonlat2xy_0D
     module procedure MAPPROJECTION_lonlat2xy_2D
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public :: MAPPROJECTION_basepoint_lon = 135.221_RP ! position of base point (domain center) in real world [deg]
  real(RP), public :: MAPPROJECTION_basepoint_lat =  34.653_RP ! position of base point (domain center) in real world [deg]

  type, public :: mappinginfo
     character(len=H_SHORT) :: mapping_name
     real(DP)               :: false_easting
     real(DP)               :: false_northing
     real(DP)               :: longitude_of_central_meridian
     real(DP)               :: longitude_of_projection_origin
     real(DP)               :: latitude_of_projection_origin
     real(DP)               :: straight_vertical_longitude_from_pole
     real(DP)               :: standard_parallel(2)
     real(DP)               :: rotation
  end type mappinginfo
  type(mappinginfo), public :: MAPPROJECTION_mappinginfo

  type, public :: mappingparam
     real(DP) :: basepoint_x
     real(DP) :: basepoint_y
     real(DP) :: basepoint_lon
     real(DP) :: basepoint_lat
     real(DP) :: rotation
     real(DP) :: rot_fact_sin
     real(DP) :: rot_fact_cos
     real(DP) :: hemisphere ! hemisphere flag: 1=north, -1=south
     real(DP) :: x          ! position of the pole or equator in the model [m]
     real(DP) :: y          !
     real(DP) :: fact       ! pre-calc factor
     real(DP) :: c          ! conformal factor
  end type mappingparam

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: MAPPROJECTION_type = 'NONE' !< map projection type
                                               ! 'NONE'
                                               ! 'LC'
                                               ! 'PS'
                                               ! 'MER'
                                               ! 'EC'

  real(DP), private :: MAPPROJECTION_basepoint_x                ! position of base point in the model [m]
  real(DP), private :: MAPPROJECTION_basepoint_y                ! position of base point in the model [m]
  real(DP), private :: MAPPROJECTION_basepoint_lonrad
  real(DP), private :: MAPPROJECTION_basepoint_latrad
  real(DP), private :: MAPPROJECTION_rotation =  0.0_DP         ! rotation factor

  real(DP), private :: MAPPROJECTION_LC_lat1  = 30.0_DP ! standard latitude1 for L.C. projection [deg]
  real(DP), private :: MAPPROJECTION_LC_lat2  = 60.0_DP ! standard latitude2 for L.C. projection [deg]
  real(DP), private :: MAPPROJECTION_PS_lat             ! standard latitude1 for P.S. projection [deg]
  real(DP), private :: MAPPROJECTION_M_lat    =  0.0_DP ! standard latitude1 for Mer. projection [deg]
  real(DP), private :: MAPPROJECTION_EC_lat   =  0.0_DP ! standard latitude1 for E.C. projection [deg]

  type(mappingparam), private :: MAPPROJECTION_mappingparam

  real(DP), private :: PI
  real(DP), private :: D2R
  real(DP), private :: RADIUS

  interface
     subroutine xy2lonlat_s( x, y, param, lon, lat )
       use scale_precision
       import mappingparam
       real(RP),           intent(in)  :: x, y
       type(mappingparam), intent(in)  :: param
       real(RP),           intent(out) :: lon, lat
     end subroutine xy2lonlat_s
     subroutine lonlat2xy_s( lon, lat, param, x, y )
       use scale_precision
       import mappingparam
       real(RP),           intent(in)  :: lon, lat
       type(mappingparam), intent(in)  :: param
       real(RP),           intent(out) :: x, y
     end subroutine lonlat2xy_s
     subroutine mapfactor_s( lat, param, m1, m2 )
       use scale_precision
       import mappingparam
       real(RP),           intent(in)  :: lat
       type(mappingparam), intent(in)  :: param
       real(RP),           intent(out) :: m1, m2
     end subroutine mapfactor_s
     subroutine rotcoef_s( lon, lat, param, rotc_cos, rotc_sin )
       use scale_precision
       import mappingparam
       real(RP),           intent(in)  :: lon, lat
       type(mappingparam), intent(in)  :: param
       real(RP),           intent(out) :: rotc_cos, rotc_sin
     end subroutine rotcoef_s
  end interface
  procedure(xy2lonlat_s), pointer :: xy2lonlat => NULL()
  procedure(lonlat2xy_s), pointer :: lonlat2xy => NULL()
  procedure(mapfactor_s), pointer :: mapfactor => NULL()
  procedure(rotcoef_s),   pointer :: rotcoef   => NULL()

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MAPPROJECTION_setup( DOMAIN_CENTER_X, DOMAIN_CENTER_Y )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF     => CONST_UNDEF,  &
       PI_RP     => CONST_PI,     &
       D2R_RP    => CONST_D2R,    &
       RADIUS_RP => CONST_RADIUS
    implicit none

    real(RP), intent(in) :: DOMAIN_CENTER_X !< center position of global domain [m]: x
    real(RP), intent(in) :: DOMAIN_CENTER_Y !< center position of global domain [m]: y

    namelist / PARAM_MAPPROJECTION / &
       MAPPROJECTION_basepoint_lon, &
       MAPPROJECTION_basepoint_lat, &
       MAPPROJECTION_basepoint_x,   &
       MAPPROJECTION_basepoint_y,   &
       MAPPROJECTION_type,          &
       MAPPROJECTION_rotation,      &
       MAPPROJECTION_LC_lat1,       &
       MAPPROJECTION_LC_lat2,       &
       MAPPROJECTION_PS_lat,        &
       MAPPROJECTION_M_lat,         &
       MAPPROJECTION_EC_lat

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_setup",*) 'Setup'

    PI     = real(PI_RP,     kind=DP)
    D2R    = real(D2R_RP,    kind=DP)
    RADIUS = real(RADIUS_RP, kind=DP)

    MAPPROJECTION_basepoint_x = DOMAIN_CENTER_X
    MAPPROJECTION_basepoint_y = DOMAIN_CENTER_Y
    MAPPROJECTION_PS_lat      = UNDEF


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MAPPROJECTION,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MAPPROJECTION_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MAPPROJECTION_setup",*) 'Not appropriate names in namelist PARAM_MAPPROJECTION. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MAPPROJECTION)

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_setup",*) 'Map projection information'
    LOG_INFO_CONT('(1x,A,F15.3)') 'Basepoint(x)    [m] : ', MAPPROJECTION_basepoint_x
    LOG_INFO_CONT('(1x,A,F15.3)') 'Basepoint(y)    [m] : ', MAPPROJECTION_basepoint_y
    LOG_INFO_CONT(*)              'Map projection type : ', trim(MAPPROJECTION_type)

    MAPPROJECTION_mappinginfo%longitude_of_central_meridian         = UNDEF
    MAPPROJECTION_mappinginfo%straight_vertical_longitude_from_pole = UNDEF
    MAPPROJECTION_mappinginfo%standard_parallel(:)                  = UNDEF

    MAPPROJECTION_mappinginfo%longitude_of_projection_origin = MAPPROJECTION_basepoint_lon
    MAPPROJECTION_mappinginfo%latitude_of_projection_origin  = MAPPROJECTION_basepoint_lat
    MAPPROJECTION_mappinginfo%false_easting                  = MAPPROJECTION_basepoint_x
    MAPPROJECTION_mappinginfo%false_northing                 = MAPPROJECTION_basepoint_y
    MAPPROJECTION_mappinginfo%rotation                       = MAPPROJECTION_rotation


    select case(trim(MAPPROJECTION_type))
    case('NONE')
       LOG_INFO_CONT(*) '=> NO map projection'
       MAPPROJECTION_mappinginfo%mapping_name = ""

       xy2lonlat => MAPPROJECTION_xy2lonlat_None
       lonlat2xy => MAPPROJECTION_lonlat2xy_None
       mapfactor => MAPPROJECTION_mapfactor_None
       rotcoef   => MAPPROJECTION_rotcoef_None
    case('LC')
       LOG_INFO_CONT(*) '=> Lambert Conformal projection'
       MAPPROJECTION_mappinginfo%mapping_name = "lambert_conformal_conic"
       MAPPROJECTION_mappinginfo%standard_parallel(:) = (/ MAPPROJECTION_LC_lat1, MAPPROJECTION_LC_lat2 /)
       MAPPROJECTION_mappinginfo%longitude_of_central_meridian = MAPPROJECTION_basepoint_lon

       xy2lonlat => MAPPROJECTION_xy2lonlat_LambertConformal
       lonlat2xy => MAPPROJECTION_lonlat2xy_LambertConformal
       mapfactor => MAPPROJECTION_mapfactor_LambertConformal
       rotcoef   => MAPPROJECTION_rotcoef_LambertConformal
    case('PS')
       LOG_INFO_CONT(*) '=> Polar Stereographic projection'
       if( MAPPROJECTION_PS_lat == UNDEF ) MAPPROJECTION_PS_lat = MAPPROJECTION_basepoint_lat
       MAPPROJECTION_mappinginfo%mapping_name = "polar_stereographic"
       MAPPROJECTION_mappinginfo%straight_vertical_longitude_from_pole = MAPPROJECTION_basepoint_lon
       MAPPROJECTION_mappinginfo%standard_parallel(1) = MAPPROJECTION_PS_lat

       xy2lonlat => MAPPROJECTION_xy2lonlat_PolarStereographic
       lonlat2xy => MAPPROJECTION_lonlat2xy_PolarStereographic
       mapfactor => MAPPROJECTION_mapfactor_PolarStereographic
       rotcoef   => MAPPROJECTION_rotcoef_PolarStereographic
    case('MER')
       LOG_INFO_CONT(*) '=> Mercator projection'
       MAPPROJECTION_mappinginfo%mapping_name = "mercator"
       MAPPROJECTION_mappinginfo%standard_parallel(1) = MAPPROJECTION_M_lat

       xy2lonlat => MAPPROJECTION_xy2lonlat_Mercator
       lonlat2xy => MAPPROJECTION_lonlat2xy_Mercator
       mapfactor => MAPPROJECTION_mapfactor_Mercator
       rotcoef   => MAPPROJECTION_rotcoef_Mercator
    case('EC')
       LOG_INFO_CONT(*) '=> Equidistant Cylindrical projection'
       MAPPROJECTION_mappinginfo%mapping_name = "equirectangular"
       MAPPROJECTION_mappinginfo%standard_parallel(1) = MAPPROJECTION_EC_lat
       MAPPROJECTION_mappinginfo%longitude_of_central_meridian = MAPPROJECTION_basepoint_lon

       xy2lonlat => MAPPROJECTION_xy2lonlat_EquidistantCylindrical
       lonlat2xy => MAPPROJECTION_lonlat2xy_EquidistantCylindrical
       mapfactor => MAPPROJECTION_mapfactor_EquidistantCylindrical
       rotcoef   => MAPPROJECTION_rotcoef_EquidistantCylindrical
    case default
       LOG_ERROR("MAPPROJECTION_setup",*) 'Unsupported MAPPROJECTION_type. STOP'
       call PRC_abort
    endselect

    call MAPPROJECTION_get_param( MAPPROJECTION_mappinginfo, & ! (in)
                                  MAPPROJECTION_mappingparam ) ! (out)

    return
  end subroutine MAPPROJECTION_setup

  subroutine MAPPROJECTION_get_param( &
       info, &
       param )
    use scale_prc, only: &
       PRC_abort
    implicit none
    type(mappinginfo),  intent(in)  :: info
    type(mappingparam), intent(out) :: param

    select case( info%mapping_name )
    case( "" )
       call MAPPROJECTION_get_param_None
    case( "lambert_conformal_conic" )
       call MAPPROJECTION_get_param_LambertConformal( info, & ! (in)
                                                      param ) ! (out)
       param%basepoint_lon = info%longitude_of_central_meridian * D2R
    case( "polar_stereographic" )
       call MAPPROJECTION_get_param_PolarStereographic( info, & ! (in)
                                                        param ) ! (out)
       param%basepoint_lon = info%straight_vertical_longitude_from_pole * D2R
    case( "mercator" )
       call MAPPROJECTION_get_param_Mercator( info, & ! (in)
                                              param ) ! (out)
       param%basepoint_lon = info%longitude_of_projection_origin * D2R
    case( "equirectangular" )
       call MAPPROJECTION_get_param_EquidistantCylindrical( info, & ! (in)
                                                            param ) ! (out)
       param%basepoint_lon = info%longitude_of_central_meridian * D2R
    case default
       LOG_ERROR("MAPPROJECTION_set_param",*) 'Unsupported mapping type. STOP'
       call PRC_abort
    endselect

    param%basepoint_x = info%false_easting
    param%basepoint_y = info%false_northing

    param%basepoint_lat = info%latitude_of_projection_origin * D2R

    param%rotation     = info%rotation * D2R
    param%rot_fact_sin = sin(MAPPROJECTION_mappingparam%rotation)
    param%rot_fact_cos = cos(MAPPROJECTION_mappingparam%rotation)

    return
  end subroutine MAPPROJECTION_get_param

  !-----------------------------------------------------------------------------
  !> (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_xy2lonlat_0D( &
       x, y,    &
       lon, lat )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    real(RP), intent(in)  :: x, y
    real(RP), intent(out) :: lon, lat ! [rad]

    !---------------------------------------------------------------------------

    if ( x == UNDEF .or. y == UNDEF ) then
       lon = UNDEF
       lat = UNDEF
       return
    end if

    call xy2lonlat( x, y,                       & ! (in)
                    MAPPROJECTION_mappingparam, & ! (in)
                    lon, lat                    ) ! (out)

    return
  end subroutine MAPPROJECTION_xy2lonlat_0D

  subroutine MAPPROJECTION_xy2lonlat_2D( &
       IA, IS, IE, JA, JS, JE, &
       x, y,    &
       lon, lat )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: x(IA,JA)
    real(RP), intent(in)  :: y(IA,JA)
    real(RP), intent(out) :: lon(IA,JA) ! [rad]
    real(RP), intent(out) :: lat(IA,JA) ! [rad]

    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       call MAPPROJECTION_xy2lonlat_0D( x(i,j), y(i,j), lon(i,j), lat(i,j) )
    end do
    end do

    return
  end subroutine MAPPROJECTION_xy2lonlat_2D

  !-----------------------------------------------------------------------------
  !> (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_lonlat2xy_0D( &
       lon, lat, &
       x, y      )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    real(RP), intent(in)  :: lon, lat ! [rad]
    real(RP), intent(out) :: x, y

    !---------------------------------------------------------------------------

    if ( lon == UNDEF .or. lat == UNDEF ) then
       x = UNDEF
       y = UNDEF

       return
    end if

    call lonlat2xy( lon, lat,                   & ! (in)
                    MAPPROJECTION_mappingparam, & ! (in)
                    x, y                       ) ! (out)

    return
  end subroutine MAPPROJECTION_lonlat2xy_0D

  subroutine MAPPROJECTION_lonlat2xy_2D( &
       IA, IS, IE, JA, JS, JE, &
       lon, lat, &
       x, y      )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: lon(IA,JA) ! [rad]
    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: x(IA,JA)
    real(RP), intent(out) :: y(IA,JA)

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       call MAPPROJECTION_lonlat2xy_0D( lon(i,j), lat(i,j), x(i,j), y(i,j) )
    end do
    end do

    return
  end subroutine MAPPROJECTION_lonlat2xy_2D

  !-----------------------------------------------------------------------------
  !> (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_mapfactor( &
       IA, IS, IE, JA, JS, JE, &
       lat,   &
       m1, m2 )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)

    real(RP) :: mm1, mm2
    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do &
    !$omp private(mm1,mm2)
    do j = JS, JE
    do i = IS, IE
       call mapfactor( lat(i,j),                   & ! (in)
                       MAPPROJECTION_mappingparam, & ! (in)
                       m1(i,j), m2(i,j)            ) ! (out)
    end do
    end do

    return
  end subroutine MAPPROJECTION_mapfactor

  !-----------------------------------------------------------------------------
  !> u(lat,lon) = cos u(x,y) - sin v(x,y)
  !> v(lat,lon) = sin u(x,y) + cos v(x,y)
  subroutine MAPPROJECTION_rotcoef( &
       IA, IS, IE, JA, JS, JE, &
       lon, lat, &
       rotc      )
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: lon(IA,JA)    ! [rad]
    real(RP), intent(in)  :: lat(IA,JA)    ! [rad]
    real(RP), intent(out) :: rotc(IA,JA,2) !< rotc(:,:,1)->cos, rotc(:,:,2)->sin

    integer :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       call rotcoef( lon(i,j), lat(i,j),         & ! (in)
                     MAPPROJECTION_mappingparam, & ! (in)
                     rotc(i,j,1), rotc(i,j,2)    ) ! (out)
    end do
    end do

    return
  end subroutine MAPPROJECTION_rotcoef

  !-----------------------------------------------------------------------------
  !> No projection
  subroutine MAPPROJECTION_get_param_None
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine MAPPROJECTION_get_param_None

  !-----------------------------------------------------------------------------
  !> No projection, lon,lat are determined by gnomonic projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_xy2lonlat_None( &
       x, y,    &
       param,   &
       lon, lat )
    implicit none
    real(RP),           intent(in)  :: x, y
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: lon, lat ! [rad]

    real(DP) :: xx, yy
    real(DP) :: gno1, gno2
    real(DP) :: rho, gmm
    real(DP) :: gmmc, gmms
    real(DP) :: latc, lats
    !---------------------------------------------------------------------------

    call xy_rotate( x, y,  & ! (in)
                    param, & ! (in)
                    xx, yy ) ! (out)

    gno1 = ( x - param%basepoint_x ) / RADIUS
    gno2 = ( y - param%basepoint_y ) / RADIUS

    rho = sqrt( gno1 * gno1 + gno2 * gno2 )
    gmm = atan( rho )

    if ( rho == 0.0_DP ) then
       lon = param%basepoint_lon
       lat = param%basepoint_lat
    else
       gmmc = cos(gmm)
       gmms = sin(gmm)
       latc = cos(param%basepoint_lat)
       lats = sin(param%basepoint_lat)
       lon = param%basepoint_lon &
           + atan( gno1*gmms / ( rho * latc * gmmc &
                              - gno2 * lats * gmms ) )
       lat = asin(      lats * gmmc       &
                 + gno2*latc * gmms / rho )
    endif

    return
  end subroutine MAPPROJECTION_xy2lonlat_None

  subroutine MAPPROJECTION_lonlat2xy_None( &
       lon, lat, &
       param,    &
       x, y      )
    implicit none
    real(RP),           intent(in)  :: lon, lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: x, y

    real(DP) :: xx, yy
    real(DP) :: lat_d, lat0_d, dlon
    real(DP) :: lat_d_c, lat_d_s, lat0_d_c, lat0_d_s, dlon_c, dlon_s
    real(DP) :: gno1, gno2
    real(DP) :: cos_gmm
    !---------------------------------------------------------------------------
    ! http://mathworld.wolfram.com/GnomonicProjection.html

    lat_d  = real(lat,kind=DP)
    lat0_d = param%basepoint_lat

    dlon = lon - param%basepoint_lon

    lat_d_c = cos(lat_d)
    lat_d_s = sin(lat_d)
    lat0_d_c = cos(lat0_d)
    lat0_d_s = sin(lat0_d)
    dlon_c = cos(dlon)
    dlon_s = sin(dlon)

    cos_gmm = lat0_d_s * lat_d_s &
            + lat0_d_c * lat_d_c * dlon_c

    gno1 = ( lat_d_c  * dlon_s  ) / cos_gmm
    gno2 = ( lat0_d_c * lat_d_s &
           - lat0_d_s * lat_d_c * dlon_c ) / cos_gmm

    xx = gno1 * RADIUS + param%basepoint_x
    yy = gno2 * RADIUS + param%basepoint_y

    call xy_unrotate( xx, yy, & ! (in)
                      param,  & ! (in)
                      x, y    ) ! (out)

    return
  end subroutine MAPPROJECTION_lonlat2xy_None

  !-----------------------------------------------------------------------------
  !> No projection: m1=m2=1
  subroutine MAPPROJECTION_mapfactor_None( &
       lat,   &
       param, &
       m1, m2 )
    implicit none
    real(RP),           intent(in)  :: lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: m1, m2
    !---------------------------------------------------------------------------

    m1 = 1.0_RP
    m2 = m1

    return
  end subroutine MAPPROJECTION_mapfactor_None

  !-----------------------------------------------------------------------------
  !> No projection:
  subroutine MAPPROJECTION_rotcoef_None( &
       lon, lat,          &
       param,             &
       rotc_cos, rotc_sin )
    implicit none
    real(RP),           intent(in)  :: lon, lat
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: rotc_cos, rotc_sin
    !---------------------------------------------------------------------------

    rotc_cos = param%rot_fact_cos
    rotc_sin = param%rot_fact_sin

    return
  end subroutine MAPPROJECTION_rotcoef_None

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection
  subroutine MAPPROJECTION_get_param_LambertConformal( &
       info, &
       param )
    use scale_prc, only: &
       PRC_abort
    implicit none
    type(mappinginfo),  intent(in)  :: info
    type(mappingparam), intent(out) :: param

    real(DP) :: LC_lat1, LC_lat2
    real(DP) :: basepoint_lat
    real(DP) :: basepoint_x, basepoint_y

    real(DP) :: lat1rot, lat2rot
    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    LC_lat1 = info%standard_parallel(1)
    LC_lat2 = info%standard_parallel(2)

    basepoint_lat = info%latitude_of_projection_origin
    basepoint_x   = info%false_easting
    basepoint_y   = info%false_northing

    if ( LC_lat1 >= LC_lat2 ) then
       LOG_ERROR("MAPPROJECTION_get_param_LambertConformal",*) 'Please set LC_lat1 < LC_lat2 in degree. STOP'
       call PRC_abort
    endif

    ! check hemisphere: 1=north, -1=south
    param%hemisphere = sign(1.0_DP,LC_lat1+LC_lat2)

    lat1rot = 0.5_DP*PI - param%hemisphere * LC_lat1 * D2R
    lat2rot = 0.5_DP*PI - param%hemisphere * LC_lat2 * D2R

    ! calc conformal factor c
    param%c = ( log( sin(lat1rot)        ) - log( sin(lat2rot)        ) ) &
            / ( log( tan(0.5_DP*lat1rot) ) - log( tan(0.5_DP*lat2rot) ) )

    ! pre-calc factor
    param%fact = sin(lat1rot) / param%c / tan(0.5_DP*lat1rot)**param%c

    ! calc (x,y) at pole point
    dlon = 0.0_DP

    latrot = 0.5_DP*PI - param%hemisphere * basepoint_lat * D2R

    dist = param%fact * RADIUS * tan(0.5_DP*latrot)**param%c

    param%x = basepoint_x -                    dist * sin(param%c*dlon)
    param%y = basepoint_y + param%hemisphere * dist * cos(param%c*dlon)

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_get_param_LambertConformal",*) 'Input parameters'
    LOG_INFO_CONT(*) 'LC_lat1    = ', LC_lat1
    LOG_INFO_CONT(*) 'LC_lat2    = ', LC_lat2
    LOG_INFO_CONT(*) 'hemisphere = ', param%hemisphere
    LOG_INFO_CONT(*) 'LC_c       = ', param%c
    LOG_INFO_CONT(*) 'LC_fact    = ', param%fact
    LOG_INFO_CONT(*) 'pole_x     = ', param%x
    LOG_INFO_CONT(*) 'pole_y     = ', param%y

    return
  end subroutine MAPPROJECTION_get_param_LambertConformal

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_xy2lonlat_LambertConformal( &
       x, y,    &
       param,   &
       lon, lat )
    implicit none
    real(RP),           intent(in)  :: x, y
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: lon, lat ! [rad]

    real(DP) :: xx, yy, dist
    !---------------------------------------------------------------------------

    call xy_rotate( x, y,  & ! (in)
                    param, & ! (in)
                    xx, yy ) ! (out)
    xx =                      ( xx - param%x ) / RADIUS / param%fact
    yy = - param%hemisphere * ( yy - param%y ) / RADIUS / param%fact

    dist = sqrt( xx*xx + yy*yy )

    lon = param%basepoint_lon + atan2(xx,yy) / param%c

    ! check hemisphere: 1=north, -1=south
    lat = param%hemisphere * ( 0.5_DP*PI - 2.0_DP*atan( dist**(1.0_DP/param%c) ) )

    return
  end subroutine MAPPROJECTION_xy2lonlat_LambertConformal

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_lonlat2xy_LambertConformal( &
       lon, lat, &
       param,    &
       x, y      )
    implicit none
    real(RP),           intent(in)  :: lon, lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: x, y

    real(DP) :: xx, yy
    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - param%basepoint_lon
    if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

    latrot = 0.5_DP*PI - param%hemisphere * lat

    dist = param%fact * RADIUS * tan(0.5_DP*latrot)**param%c

    xx = param%x +                    dist * sin(param%c*dlon)
    yy = param%y - param%hemisphere * dist * cos(param%c*dlon)

    call xy_unrotate( xx, yy, & ! (in)
                      param,  & ! (in)
                      x, y    ) ! (out)

    return
  end subroutine MAPPROJECTION_lonlat2xy_LambertConformal

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (lon,lat) -> (m1=m2)
  subroutine MAPPROJECTION_mapfactor_LambertConformal( &
       lat,   &
       param, &
       m1, m2 )
    implicit none
    real(RP),           intent(in)  :: lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: m1, m2

    real(DP) :: latrot
    !---------------------------------------------------------------------------

    latrot = 0.5_DP*PI - param%hemisphere * lat

    m1 = param%fact / sin(latrot) * param%c * tan(0.5_DP*latrot)**param%c
    m2 = m1

    return
  end subroutine MAPPROJECTION_mapfactor_LambertConformal

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_rotcoef_LambertConformal( &
       lon, lat,          &
       param,             &
       rotc_cos, rotc_sin )
    implicit none
    real(RP),           intent(in)  :: lon, lat
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: rotc_cos, rotc_sin

    real(DP) :: dlon
    real(DP) :: alpha
    !---------------------------------------------------------------------------

    dlon = lon - param%basepoint_lon
    if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

    alpha = - param%c * dlon * param%hemisphere &
            + param%rotation

    rotc_cos = cos( alpha )
    rotc_sin = sin( alpha )

    return
  end subroutine MAPPROJECTION_rotcoef_LambertConformal

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection
  subroutine MAPPROJECTION_get_param_PolarStereographic( &
       info, &
       param )
    implicit none
    type(mappinginfo),  intent(in)  :: info
    type(mappingparam), intent(out) :: param

    real(DP) :: PS_lat
    real(DP) :: basepoint_lat
    real(DP) :: basepoint_x, basepoint_y

    real(DP) :: lat0
    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    PS_lat = info%standard_parallel(1)

    basepoint_lat = info%latitude_of_projection_origin
    basepoint_x   = info%false_easting
    basepoint_y   = info%false_northing

    ! check hemisphere: 1=north, -1=south
    param%hemisphere = sign(1.0_DP,PS_lat)

    lat0 = param%hemisphere * PS_lat * D2R

    ! pre-calc factor
    param%fact = 1.0_DP + sin(lat0)

    ! calc (x,y) at pole point
    dlon = 0.0_DP

    latrot = 0.5_DP*PI - param%hemisphere * basepoint_lat * D2R

    dist = param%fact * RADIUS * tan(0.5_DP*latrot)

    param%x = basepoint_x -                    dist * sin(dlon)
    param%y = basepoint_y + param%hemisphere * dist * cos(dlon)

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_get_param_PolarStereographic",*) 'PS_lat1    = ', PS_lat
    LOG_INFO("MAPPROJECTION_get_param_PolarStereographic",*) 'hemisphere = ', param%hemisphere
    LOG_INFO("MAPPROJECTION_get_param_PolarStereographic",*) 'PS_fact    = ', param%fact
    LOG_INFO("MAPPROJECTION_get_param_PolarStereographic",*) 'pole_x     = ', param%x
    LOG_INFO("MAPPROJECTION_get_param_PolarStereographic",*) 'pole_y     = ', param%y

    return
  end subroutine MAPPROJECTION_get_param_PolarStereographic

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_xy2lonlat_PolarStereographic( &
       x, y,    &
       param,   &
       lon, lat )
    implicit none
    real(RP),           intent(in)  :: x, y
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: lon, lat ! [rad]

    real(DP) :: xx, yy, dist
    !---------------------------------------------------------------------------

    call xy_rotate( x, y,  & ! (in)
                    param, & ! (in)
                    xx, yy ) ! (out)

    xx =                      ( xx - param%x ) / RADIUS / param%fact
    yy = - param%hemisphere * ( yy - param%y ) / RADIUS / param%fact

    dist = sqrt( xx*xx + yy*yy )

    lon = param%basepoint_lon + atan2(xx,yy)
    lat = param%hemisphere * ( 0.5_DP*PI - 2.0_DP*atan(dist) )

    return
  end subroutine MAPPROJECTION_xy2lonlat_PolarStereographic

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_lonlat2xy_PolarStereographic( &
       lon, lat, &
       param,    &
       x, y      )
    implicit none
    real(RP),           intent(in)  :: lon, lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: x, y

    real(DP) :: dlon, latrot, dist
    real(DP) :: xx, yy
    !---------------------------------------------------------------------------

    dlon = lon - param%basepoint_lon

    latrot = 0.5_DP*PI - param%hemisphere * lat

    dist = param%fact * RADIUS * tan(0.5_DP*latrot)

    xx = param%x +                    dist * sin(dlon)
    yy = param%y - param%hemisphere * dist * cos(dlon)

    call xy_unrotate( xx, yy, & ! (in)
                      param,  & ! (in)
                      x, y    ) ! (out)

    return
  end subroutine MAPPROJECTION_lonlat2xy_PolarStereographic

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (lon,lat) -> (m1=m2)
  subroutine MAPPROJECTION_mapfactor_PolarStereographic( &
       lat,   &
       param, &
       m1, m2 )
    implicit none
    real(RP),           intent(in)  :: lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: m1, m2
    !---------------------------------------------------------------------------

    m1 = param%fact / ( 1.0_DP + sin(param%hemisphere*lat) )
    m2 = m1

    return
  end subroutine MAPPROJECTION_mapfactor_PolarStereographic

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_rotcoef_PolarStereographic( &
       lon, lat,          &
       param,             &
       rotc_cos, rotc_sin )
    implicit none
    real(RP),           intent(in)  :: lon, lat
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: rotc_cos, rotc_sin

    real(DP) :: dlon
    real(DP) :: alpha
    !---------------------------------------------------------------------------

    dlon = lon - param%basepoint_lon
    if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

    alpha = - dlon * param%hemisphere &
            + param%rotation

    rotc_cos = cos( alpha )
    rotc_sin = sin( alpha )

    return
  end subroutine MAPPROJECTION_rotcoef_PolarStereographic

  !-----------------------------------------------------------------------------
  !> Mercator projection
  subroutine MAPPROJECTION_get_param_Mercator( &
       info, &
       param )
    use scale_prc, only: &
       PRC_abort
    implicit none
    type(mappinginfo),  intent(in)  :: info
    type(mappingparam), intent(out) :: param

    real(DP) :: M_lat
    real(DP) :: basepoint_lat
    real(DP) :: basepoint_x, basepoint_y

    real(DP) :: lat0
    real(DP) :: latrot, dist
    !---------------------------------------------------------------------------

    M_lat = info%standard_parallel(1)

    basepoint_lat = info%latitude_of_projection_origin
    basepoint_x   = info%false_easting
    basepoint_y   = info%false_northing

    lat0 = M_lat * D2R

    ! pre-calc factor
    param%fact = cos(lat0)

    if ( param%fact == 0.0_DP ) then
       LOG_ERROR("MAPPROJECTION_get_param_Mercator",*) 'M_lat cannot be set to pole point! value=', M_lat
       call PRC_abort
    endif

    ! calc (x,y) at (lon,lat) = (base,0)
    latrot = 0.5_DP*PI - basepoint_lat * D2R

    dist = 1.0_DP / tan(0.5_DP*latrot)

    param%x = basepoint_x
    param%y = basepoint_y - RADIUS * param%fact * log(dist)

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_get_param_Mercator",*) 'M_lat  = ', M_lat
    LOG_INFO("MAPPROJECTION_get_param_Mercator",*) 'M_fact = ', param%fact
    LOG_INFO("MAPPROJECTION_get_param_Mercator",*) 'eq_x   = ', param%x
    LOG_INFO("MAPPROJECTION_get_param_Mercator",*) 'eq_y   = ', param%y

    return
  end subroutine MAPPROJECTION_get_param_Mercator

  !-----------------------------------------------------------------------------
  !> Mercator projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_xy2lonlat_Mercator( &
       x, y,    &
       param,   &
       lon, lat )
    implicit none
    real(RP),           intent(in)  :: x, y
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: lon, lat ! [rad]

    real(DP) :: xx, yy
    !---------------------------------------------------------------------------

    call xy_rotate( x, y,  & ! (in)
                    param, & ! (in)
                    xx, yy ) ! (out)

    xx = ( xx - param%x ) / RADIUS / param%fact
    yy = ( yy - param%y ) / RADIUS / param%fact

    lon = xx + param%basepoint_lon

    lat = 0.5_DP*PI - 2.0_DP*atan( 1.0_DP/exp(yy) )

    return
  end subroutine MAPPROJECTION_xy2lonlat_Mercator

  !-----------------------------------------------------------------------------
  !> Mercator projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_lonlat2xy_Mercator( &
       lon, lat, &
       param,    &
       x, y      )
    implicit none
    real(RP),           intent(in)  :: lon, lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: x, y

    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - param%basepoint_lon

    latrot = 0.5_DP*PI - lat

    dist = 1.0_DP / tan(0.5_DP*latrot)

    x = param%x + RADIUS * param%fact * dlon
    y = param%y + RADIUS * param%fact * log(dist)

    return
  end subroutine MAPPROJECTION_lonlat2xy_Mercator

  !-----------------------------------------------------------------------------
  !> Mercator projection: (lon,lat) -> (m1=m2)
  subroutine MAPPROJECTION_mapfactor_Mercator( &
       lat,   &
       param, &
       m1, m2 )
    implicit none
    real(RP),           intent(in)  :: lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: m1, m2
    !---------------------------------------------------------------------------

    m1 = param%fact / cos(real(lat,kind=DP))
    m2 = m1

    return
  end subroutine MAPPROJECTION_mapfactor_Mercator

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_rotcoef_Mercator( &
       lon, lat,          &
       param,             &
       rotc_cos, rotc_sin )
    implicit none
    real(RP),           intent(in)  :: lon, lat
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: rotc_cos, rotc_sin
    !---------------------------------------------------------------------------

    rotc_cos = param%rot_fact_cos
    rotc_sin = param%rot_fact_sin

    return
  end subroutine MAPPROJECTION_rotcoef_Mercator

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection
  subroutine MAPPROJECTION_get_param_EquidistantCylindrical( &
       info, &
       param )
    use scale_prc, only: &
       PRC_abort
    implicit none
    type(mappinginfo),  intent(in)  :: info
    type(mappingparam), intent(out) :: param

    real(DP) :: EC_lat
    real(DP) :: basepoint_lat
    real(DP) :: basepoint_x, basepoint_y

    real(DP) :: lat0
    !---------------------------------------------------------------------------

    EC_lat = info%standard_parallel(1)

    basepoint_lat = info%latitude_of_projection_origin
    basepoint_x   = info%false_easting
    basepoint_y   = info%false_northing

    lat0 = EC_lat * D2R

    ! pre-calc factor
    param%fact = cos(lat0)

    if ( param%fact == 0.0_DP ) then
       LOG_ERROR("MAPPROJECTION_get_param_EquidistantCylindrical",*) 'EC_lat cannot be set to pole point! value=', EC_lat
       call PRC_abort
    endif

    param%x = basepoint_x
    param%y = basepoint_y - RADIUS * basepoint_lat * D2R

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_get_param_EquidistantCylindrical",*) 'EC_lat  = ', EC_lat
    LOG_INFO("MAPPROJECTION_get_param_EquidistantCylindrical",*) 'EC_fact = ', param%fact
    LOG_INFO("MAPPROJECTION_get_param_EquidistantCylindrical",*) 'eq_x    = ', param%x
    LOG_INFO("MAPPROJECTION_get_param_EquidistantCylindrical",*) 'eq_y    = ', param%y

    return
  end subroutine MAPPROJECTION_get_param_EquidistantCylindrical

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_xy2lonlat_EquidistantCylindrical( &
       x, y,    &
       param,   &
       lon, lat )
    use scale_prc, only: &
       PRC_abort
    implicit none
    real(RP),           intent(in)  :: x, y
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: lon, lat ! [rad]

    real(DP) :: xx, yy
    !---------------------------------------------------------------------------

    call xy_rotate( x, y,  & ! (in)
                    param, & ! (in)
                    xx, yy ) ! (out)

    xx = ( xx - param%x ) / RADIUS / param%fact
    yy = ( yy - param%y ) / RADIUS

    lon = xx + param%basepoint_lon
    lat = yy

    if ( abs(lat) >  0.5_DP*PI ) then
       LOG_ERROR("MAPPROJECTION_EquidistantCylindrical_xy2lonlat",*) 'Invalid latitude range! value=', lat
       call PRC_abort
    endif

    return
  end subroutine MAPPROJECTION_xy2lonlat_EquidistantCylindrical

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_lonlat2xy_EquidistantCylindrical( &
       lon, lat, &
       param,    &
       x, y      )
    implicit none
    real(RP),           intent(in)  :: lon, lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: x, y

    real(DP) :: dlon
    real(DP) :: xx, yy
    !---------------------------------------------------------------------------

    dlon = lon - param%basepoint_lon

    xx = param%x + RADIUS * param%fact * dlon
    yy = param%y + RADIUS * lat

    call xy_unrotate( xx, yy, & ! (in)
                      param,  & ! (in)
                      x, y    ) ! (out)

    return
  end subroutine MAPPROJECTION_lonlat2xy_EquidistantCylindrical

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (lon,lat) -> (m1,m2)
  subroutine MAPPROJECTION_mapfactor_EquidistantCylindrical( &
       lat,   &
       param, &
       m1, m2 )
    implicit none
    real(RP), intent(in)  :: lat ! [rad]
    type(mappingparam), intent(in)  :: param
    real(RP), intent(out) :: m1, m2

    real(DP) :: mm1, mm2
    !---------------------------------------------------------------------------

    mm1 = param%fact / cos(real(lat,kind=DP))
    mm2 = 1.0_DP

    m1 = sqrt( (param%rot_fact_cos * mm1)**2 + (param%rot_fact_sin * mm2)**2 )
    m2 = sqrt( (param%rot_fact_cos * mm2)**2 + (param%rot_fact_sin * mm1)**2 )

    return
  end subroutine MAPPROJECTION_mapfactor_EquidistantCylindrical

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_rotcoef_EquidistantCylindrical( &
       lon, lat,          &
       param,             &
       rotc_cos, rotc_sin )
    implicit none
    real(RP),           intent(in)  :: lon, lat
    type(mappingparam), intent(in)  :: param
    real(RP),           intent(out) :: rotc_cos, rotc_sin
    !---------------------------------------------------------------------------

    rotc_cos = param%rot_fact_cos
    rotc_sin = param%rot_fact_sin

    return
  end subroutine MAPPROJECTION_rotcoef_EquidistantCylindrical

! private

  subroutine xy_rotate( &
       x, y,  &
       param, &
       xx, yy )
    implicit none
    real(RP),           intent(in)  :: x, y
    type(mappingparam), intent(in)  :: param
    real(DP),           intent(out) :: xx, yy

    real(DP) :: xd, yd

    xd = x - param%basepoint_x
    yd = y - param%basepoint_y

    xx = param%basepoint_x + xd * param%rot_fact_cos &
                           - yd * param%rot_fact_sin
    yy = param%basepoint_y + yd * param%rot_fact_cos &
                           + xd * param%rot_fact_sin

    return
  end subroutine xy_rotate

  subroutine xy_unrotate( &
       xx, yy, &
       param,  &
       x, y    )
    implicit none
    real(DP),           intent(in)  :: xx, yy
    type(mappingparam), intent(in) :: param
    real(RP),           intent(out) :: x, y

    real(DP) :: xxd, yyd

    xxd = xx - param%basepoint_x
    yyd = yy - param%basepoint_y

    x = param%basepoint_x + ( xxd * param%rot_fact_cos &
                            + yyd * param%rot_fact_sin )
    y = param%basepoint_y + ( yyd * param%rot_fact_cos &
                            - xxd * param%rot_fact_sin )

    return
  end subroutine xy_unrotate

end module scale_mapprojection
