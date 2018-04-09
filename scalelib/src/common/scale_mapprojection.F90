!-------------------------------------------------------------------------------
!> module Map projection
!!
!! @par Description
!!          Map projection module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-10-24 (H.Yashiro)  [new]
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
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MAPPROJECTION_setup
  public :: MAPPROJECTION_xy2lonlat
  public :: MAPPROJECTION_lonlat2xy
  public :: MAPPROJECTION_mapfactor
  public :: MAPPROJECTION_rotcoef
  public :: MAPPROJECTION_get_attributes

  interface MAPPROJECTION_rotcoef
     module procedure MAPPROJECTION_rotcoef_0D
     module procedure MAPPROJECTION_rotcoef_2D
  end interface MAPPROJECTION_rotcoef

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public :: MAPPROJECTION_basepoint_lon = 135.221_RP ! position of base point (domain center) in real world [deg]
  real(RP), public :: MAPPROJECTION_basepoint_lat =  34.653_RP ! position of base point (domain center) in real world [deg]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MAPPROJECTION_None_setup
  private :: MAPPROJECTION_LambertConformal_setup
  private :: MAPPROJECTION_PolarStereographic_setup
  private :: MAPPROJECTION_Mercator_setup
  private :: MAPPROJECTION_EquidistantCylindrical_setup

  private :: MAPPROJECTION_None_xy2lonlat
  private :: MAPPROJECTION_LambertConformal_xy2lonlat
  private :: MAPPROJECTION_PolarStereographic_xy2lonlat
  private :: MAPPROJECTION_Mercator_xy2lonlat
  private :: MAPPROJECTION_EquidistantCylindrical_xy2lonlat

  private :: MAPPROJECTION_None_lonlat2xy
  private :: MAPPROJECTION_LambertConformal_lonlat2xy
  private :: MAPPROJECTION_PolarStereographic_lonlat2xy
  private :: MAPPROJECTION_Mercator_lonlat2xy
  private :: MAPPROJECTION_EquidistantCylindrical_lonlat2xy

  private :: MAPPROJECTION_None_mapfactor
  private :: MAPPROJECTION_LambertConformal_mapfactor
  private :: MAPPROJECTION_PolarStereographic_mapfactor
  private :: MAPPROJECTION_Mercator_mapfactor
  private :: MAPPROJECTION_EquidistantCylindrical_mapfactor

  private :: MAPPROJECTION_None_rotcoef_2D
  private :: MAPPROJECTION_LambertConformal_rotcoef_2D
  private :: MAPPROJECTION_PolarStereographic_rotcoef_2D
  private :: MAPPROJECTION_Mercator_rotcoef_2D
  private :: MAPPROJECTION_EquidistantCylindrical_rotcoef_2D

  private :: MAPPROJECTION_None_rotcoef_0D
  private :: MAPPROJECTION_LambertConformal_rotcoef_0D
  private :: MAPPROJECTION_PolarStereographic_rotcoef_0D
  private :: MAPPROJECTION_Mercator_rotcoef_0D
  private :: MAPPROJECTION_EquidistantCylindrical_rotcoef_0D

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

  real(DP), private :: MAPPROJECTION_hemisphere         ! hemisphere flag: 1=north, -1=south

  real(DP), private :: MAPPROJECTION_basepoint_x        ! position of base point in the model [m]
  real(DP), private :: MAPPROJECTION_basepoint_y        ! position of base point in the model [m]

  real(DP), private :: MAPPROJECTION_pole_x             ! position of north/south pole in the model [m]
  real(DP), private :: MAPPROJECTION_pole_y             ! position of north/south pole in the model [m]
  real(DP), private :: MAPPROJECTION_eq_x               ! position of equator at the base lon. in the model [m]
  real(DP), private :: MAPPROJECTION_eq_y               ! position of equator at the base lon. in the model [m]

  real(DP), private :: MAPPROJECTION_rotation =  0.0_DP ! rotation factor (only for 'NONE' type)

  real(DP), private :: MAPPROJECTION_LC_lat1  = 30.0_DP ! standard latitude1 for L.C. projection [deg]
  real(DP), private :: MAPPROJECTION_LC_lat2  = 60.0_DP ! standard latitude2 for L.C. projection [deg]
  real(DP), private :: MAPPROJECTION_LC_c               ! conformal factor
  real(DP), private :: MAPPROJECTION_LC_fact            ! pre-calc factor

  real(DP), private :: MAPPROJECTION_PS_lat             ! standard latitude1 for P.S. projection [deg]
  real(DP), private :: MAPPROJECTION_PS_fact            ! pre-calc factor

  real(DP), private :: MAPPROJECTION_M_lat    =  0.0_DP ! standard latitude1 for Mer. projection [deg]
  real(DP), private :: MAPPROJECTION_M_fact             ! pre-calc factor

  real(DP), private :: MAPPROJECTION_EC_lat   =  0.0_DP ! standard latitude1 for E.C. projection [deg]
  real(DP), private :: MAPPROJECTION_EC_fact            ! pre-calc factor

  real(DP), private :: PI
  real(DP), private :: D2R
  real(DP), private :: RADIUS

  character(len=H_SHORT) :: MAPPROJECTION_mapping
  real(DP), private :: MAPPROJECTION_false_easting
  real(DP), private :: MAPPROJECTION_false_northing
  real(DP), private :: MAPPROJECTION_longitude_of_central_meridian
  real(DP), private :: MAPPROJECTION_longitude_of_projection_origin
  real(DP), private :: MAPPROJECTION_latitude_of_projection_origin
  real(DP), private :: MAPPROJECTION_straight_vertical_longitude_from_pole
  real(DP), private :: MAPPROJECTION_standard_parallel(2)

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

    MAPPROJECTION_basepoint_x = UNDEF
    MAPPROJECTION_basepoint_y = UNDEF
    MAPPROJECTION_PS_lat      = UNDEF

    MAPPROJECTION_basepoint_x = DOMAIN_CENTER_X
    MAPPROJECTION_basepoint_y = DOMAIN_CENTER_Y

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
    LOG_INFO_CONT('(1x,A,F15.3)') 'Basepoint(x)       : ', MAPPROJECTION_basepoint_x
    LOG_INFO_CONT('(1x,A,F15.3)') 'Basepoint(y)       : ', MAPPROJECTION_basepoint_y
    LOG_INFO_CONT(*)              'Map projection type: ', trim(MAPPROJECTION_type)

    MAPPROJECTION_mapping = ""
    MAPPROJECTION_false_easting = UNDEF
    MAPPROJECTION_false_northing = UNDEF
    MAPPROJECTION_longitude_of_central_meridian = UNDEF
    MAPPROJECTION_longitude_of_projection_origin = UNDEF
    MAPPROJECTION_latitude_of_projection_origin = UNDEF
    MAPPROJECTION_straight_vertical_longitude_from_pole = UNDEF
    MAPPROJECTION_standard_parallel(:) = UNDEF

    if( MAPPROJECTION_PS_lat == UNDEF ) MAPPROJECTION_PS_lat = MAPPROJECTION_basepoint_lat

    select case(trim(MAPPROJECTION_type))
    case('NONE')
       LOG_INFO_CONT(*) '=> NO map projection'
       call MAPPROJECTION_None_setup
    case('LC')
       LOG_INFO_CONT(*) '=> Lambert Conformal projection'
       call MAPPROJECTION_LambertConformal_setup
    case('PS')
       LOG_INFO_CONT(*) '=> Polar Stereographic projection'
       call MAPPROJECTION_PolarStereographic_setup
    case('MER')
       LOG_INFO_CONT(*) '=> Mercator projection'
       call MAPPROJECTION_Mercator_setup
    case('EC')
       LOG_INFO_CONT(*) '=> Equidistant Cylindrical projection'
       call MAPPROJECTION_EquidistantCylindrical_setup
    case default
       LOG_ERROR("MAPPROJECTION_setup",*) 'Unsupported MAPPROJECTION_type. STOP'
       call PRC_abort
    endselect

    return
  end subroutine MAPPROJECTION_setup

  !-----------------------------------------------------------------------------
  !> (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_xy2lonlat( &
       x,   &
       y,   &
       lon, &
       lat  )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in)  :: x
    real(RP), intent(in)  :: y
    real(RP), intent(out) :: lon ! [rad]
    real(RP), intent(out) :: lat ! [rad]
    !---------------------------------------------------------------------------

    select case(MAPPROJECTION_type)
    case('NONE')
       call MAPPROJECTION_None_xy2lonlat( x, y, lon, lat )
    case('LC')
       call MAPPROJECTION_LambertConformal_xy2lonlat( x, y, lon, lat )
    case('PS')
       call MAPPROJECTION_PolarStereographic_xy2lonlat( x, y, lon, lat )
    case('MER')
       call MAPPROJECTION_Mercator_xy2lonlat( x, y, lon, lat )
    case('EC')
       call MAPPROJECTION_EquidistantCylindrical_xy2lonlat( x, y, lon, lat )
    case default
       LOG_ERROR("MAPPROJECTION_xy2lonlat",*) 'Unsupported MAPPROJECTION_type. STOP'
       call PRC_abort
    endselect

    return
  end subroutine MAPPROJECTION_xy2lonlat

  !-----------------------------------------------------------------------------
  !> (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y
    !---------------------------------------------------------------------------

    select case(MAPPROJECTION_type)
    case('NONE')
       call MAPPROJECTION_None_lonlat2xy( lon, lat, x, y )
    case('LC')
       call MAPPROJECTION_LambertConformal_lonlat2xy( lon, lat, x, y )
    case('PS')
       call MAPPROJECTION_PolarStereographic_lonlat2xy( lon, lat, x, y )
    case('MER')
       call MAPPROJECTION_Mercator_lonlat2xy( lon, lat, x, y )
    case('EC')
       call MAPPROJECTION_EquidistantCylindrical_lonlat2xy( lon, lat, x, y )
    case default
       LOG_ERROR("MAPPROJECTION_lonlat2xy",*) 'Unsupported MAPPROJECTION_type. STOP'
       call PRC_abort
    endselect

    return
  end subroutine MAPPROJECTION_lonlat2xy

  !-----------------------------------------------------------------------------
  !> (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_mapfactor( &
       lat, &
       m1,  &
       m2   )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)
    !---------------------------------------------------------------------------

    select case(MAPPROJECTION_type)
    case('NONE')
       call MAPPROJECTION_None_mapfactor( lat, m1, m2 )
    case('LC')
       call MAPPROJECTION_LambertConformal_mapfactor( lat, m1, m2 )
    case('PS')
       call MAPPROJECTION_PolarStereographic_mapfactor( lat, m1, m2 )
    case('MER')
       call MAPPROJECTION_Mercator_mapfactor( lat, m1, m2 )
    case('EC')
       call MAPPROJECTION_EquidistantCylindrical_mapfactor( lat, m1, m2 )
    case default
       LOG_ERROR("MAPPROJECTION_mapfactor",*) 'Unsupported MAPPROJECTION_type. STOP'
       call PRC_abort
    endselect

    return
  end subroutine MAPPROJECTION_mapfactor

  !-----------------------------------------------------------------------------
  !> u(lat,lon) = cos u(x,y) - sin v(x,y)
  !> v(lat,lon) = sin u(x,y) + cos v(x,y)
  subroutine MAPPROJECTION_rotcoef_0D( &
       rotc, &
       lon,  &
       lat   )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(out) :: rotc(2) !< rotc(:,:,1)->cos, rotc(:,:,2)->sin
    real(RP), intent(in)  :: lon   ! [rad]
    real(RP), intent(in)  :: lat   ! [rad]
    !---------------------------------------------------------------------------

    select case(MAPPROJECTION_type)
    case('NONE')
       call MAPPROJECTION_None_rotcoef_0D( rotc )
    case('LC')
       call MAPPROJECTION_LambertConformal_rotcoef_0D( rotc, lon, lat )
    case('PS')
       call MAPPROJECTION_PolarStereographic_rotcoef_0D( rotc, lon, lat )
    case('MER')
       call MAPPROJECTION_Mercator_rotcoef_0D( rotc )
    case('EC')
       call MAPPROJECTION_EquidistantCylindrical_rotcoef_0D( rotc )
    case default
       LOG_ERROR("MAPPROJECTION_rotcoef_0D",*) 'Unsupported MAPPROJECTION_type. STOP'
       call PRC_abort
    endselect

    return
  end subroutine MAPPROJECTION_rotcoef_0D

  !-----------------------------------------------------------------------------
  !> u(lat,lon) = cos u(x,y) - sin v(x,y)
  !> v(lat,lon) = sin u(x,y) + cos v(x,y)
  subroutine MAPPROJECTION_rotcoef_2D( &
       rotc, &
       lon,  &
       lat   )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2) !< rotc(:,:,1)->cos, rotc(:,:,2)->sin
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]
    !---------------------------------------------------------------------------

    select case(MAPPROJECTION_type)
    case('NONE')
       call MAPPROJECTION_None_rotcoef_2D( rotc )
    case('LC')
       call MAPPROJECTION_LambertConformal_rotcoef_2D( rotc, lon, lat )
    case('PS')
       call MAPPROJECTION_PolarStereographic_rotcoef_2D( rotc, lon, lat )
    case('MER')
       call MAPPROJECTION_Mercator_rotcoef_2D( rotc )
    case('EC')
       call MAPPROJECTION_EquidistantCylindrical_rotcoef_2D( rotc )
    case default
       LOG_ERROR("MAPPROJECTION_rotcoef_2D",*) 'Unsupported MAPPROJECTION_type. STOP'
       call PRC_abort
    endselect

    return
  end subroutine MAPPROJECTION_rotcoef_2D


  !-----------------------------------------------------------------------------
  !> Get mapping attributes
  subroutine MAPPROJECTION_get_attributes( &
       mapping,                               &
       false_easting,                         &
       false_northing,                        &
       longitude_of_central_meridian,         &
       longitude_of_projection_origin,        &
       latitude_of_projection_origin,         &
       straight_vertical_longitude_from_pole, &
       standard_parallel                      )
    implicit none

    character(len=*), intent(out) :: mapping

    real(DP), intent(out), optional :: false_easting
    real(DP), intent(out), optional :: false_northing
    real(DP), intent(out), optional :: longitude_of_central_meridian
    real(DP), intent(out), optional :: longitude_of_projection_origin
    real(DP), intent(out), optional :: latitude_of_projection_origin
    real(DP), intent(out), optional :: straight_vertical_longitude_from_pole
    real(DP), intent(out), optional :: standard_parallel(2)
    !---------------------------------------------------------------------------

    mapping = MAPPROJECTION_mapping

    if ( present(false_easting)                         ) &
                 false_easting                          = MAPPROJECTION_false_easting
    if ( present(false_northing)                        ) &
                 false_northing                         = MAPPROJECTION_false_northing
    if ( present(longitude_of_central_meridian)         ) &
                 longitude_of_central_meridian          = MAPPROJECTION_longitude_of_central_meridian
    if ( present(longitude_of_projection_origin)        ) &
                 longitude_of_projection_origin         = MAPPROJECTION_longitude_of_projection_origin
    if ( present(latitude_of_projection_origin)         ) &
                 latitude_of_projection_origin          = MAPPROJECTION_latitude_of_projection_origin
    if ( present(straight_vertical_longitude_from_pole) ) &
                 straight_vertical_longitude_from_pole  = MAPPROJECTION_straight_vertical_longitude_from_pole
    if ( present(standard_parallel)                     ) &
                 standard_parallel(:)                   = MAPPROJECTION_standard_parallel(:)

    return
  end subroutine MAPPROJECTION_get_attributes

  !-----------------------------------------------------------------------------
  !> No projection
  subroutine MAPPROJECTION_None_setup
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_None_setup",*) 'MAPPROJECTION_rotation   = ', MAPPROJECTION_rotation

    MAPPROJECTION_mapping = ""

    return
  end subroutine MAPPROJECTION_None_setup

  !-----------------------------------------------------------------------------
  !> No projection, lon,lat are determined by gnomonic projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_None_xy2lonlat( &
       x,   &
       y,   &
       lon, &
       lat  )
    implicit none

    real(RP), intent(in)  :: x
    real(RP), intent(in)  :: y
    real(RP), intent(out) :: lon ! [rad]
    real(RP), intent(out) :: lat ! [rad]

    real(DP) :: gno(2)
    real(DP) :: rho, gmm
    !---------------------------------------------------------------------------

    gno(1) = ( (y-MAPPROJECTION_basepoint_y) * sin(MAPPROJECTION_rotation*D2R) &
             + (x-MAPPROJECTION_basepoint_x) * cos(MAPPROJECTION_rotation*D2R) ) / RADIUS
    gno(2) = ( (y-MAPPROJECTION_basepoint_y) * cos(MAPPROJECTION_rotation*D2R) &
             - (x-MAPPROJECTION_basepoint_x) * sin(MAPPROJECTION_rotation*D2R) ) / RADIUS

    rho = sqrt( gno(1) * gno(1) + gno(2) * gno(2) )
    gmm = atan( rho )

    if ( rho == 0.0_DP ) then
       lon = MAPPROJECTION_basepoint_lon * D2R
       lat = MAPPROJECTION_basepoint_lat * D2R
    else
       lon = MAPPROJECTION_basepoint_lon * D2R &
           + atan( gno(1)*sin(gmm) / ( rho   *cos(MAPPROJECTION_basepoint_lat*D2R)*cos(gmm) &
                                     - gno(2)*sin(MAPPROJECTION_basepoint_lat*D2R)*sin(gmm) ) )
       lat = asin(        sin(MAPPROJECTION_basepoint_lat*D2R)*cos(gmm)       &
                 + gno(2)*cos(MAPPROJECTION_basepoint_lat*D2R)*sin(gmm) / rho )
    endif

    return
  end subroutine MAPPROJECTION_None_xy2lonlat

  !-----------------------------------------------------------------------------
  !> No projection
  subroutine MAPPROJECTION_None_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y

    real(DP) :: lat_d, lat0_d, dlon
    real(DP) :: gno(2)
    real(DP) :: cos_gmm
    !---------------------------------------------------------------------------
    ! http://mathworld.wolfram.com/GnomonicProjection.html

    lat_d  = real(lat,kind=DP)
    lat0_d = real(MAPPROJECTION_basepoint_lat*D2R,kind=DP)

    dlon = lon - MAPPROJECTION_basepoint_lon * D2R

    cos_gmm = sin(lat0_d) * sin(lat_d) &
            + cos(lat0_d) * cos(lat_d) * cos(dlon)

    gno(1) = ( cos(lat_d)  * sin(dlon)  ) / cos_gmm
    gno(2) = ( cos(lat0_d) * sin(lat_d) &
             - sin(lat0_d) * cos(lat_d) * cos(dlon) ) / cos_gmm

    x = MAPPROJECTION_basepoint_x + ( gno(1) * cos(MAPPROJECTION_rotation*D2R) &
                           - gno(2) * sin(MAPPROJECTION_rotation*D2R) ) * RADIUS
    y = MAPPROJECTION_basepoint_y + ( gno(1) * sin(MAPPROJECTION_rotation*D2R) &
                           + gno(2) * cos(MAPPROJECTION_rotation*D2R) ) * RADIUS

    return
  end subroutine MAPPROJECTION_None_lonlat2xy

  !-----------------------------------------------------------------------------
  !> No projection: m1=m2=1
  subroutine MAPPROJECTION_None_mapfactor( &
       lat, &
       m1,  &
       m2   )
    implicit none

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)
    !---------------------------------------------------------------------------

    m1 = 1.0_RP
    m2 = m1

    return
  end subroutine MAPPROJECTION_None_mapfactor

  !-----------------------------------------------------------------------------
  !> No projection:
  subroutine MAPPROJECTION_None_rotcoef_0D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(2)
    !---------------------------------------------------------------------------

    rotc(1) = 1.0_RP
    rotc(2) = 0.0_RP

    return
  end subroutine MAPPROJECTION_None_rotcoef_0D

  !-----------------------------------------------------------------------------
  !> No projection:
  subroutine MAPPROJECTION_None_rotcoef_2D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MAPPROJECTION_None_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection
  subroutine MAPPROJECTION_LambertConformal_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP) :: lat1rot, lat2rot
    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    if ( MAPPROJECTION_LC_lat1 >= MAPPROJECTION_LC_lat2 ) then
       LOG_ERROR("MAPPROJECTION_LambertConformal_setup",*) 'Please set MAPPROJECTION_LC_lat1 < MAPPROJECTION_LC_lat2 in degree. STOP'
       call PRC_abort
    endif

    ! check hemisphere: 1=north, -1=south
    MAPPROJECTION_hemisphere = sign(1.0_DP,MAPPROJECTION_LC_lat1+MAPPROJECTION_LC_lat2)

    lat1rot = 0.5_DP*PI - MAPPROJECTION_hemisphere * MAPPROJECTION_LC_lat1 * D2R
    lat2rot = 0.5_DP*PI - MAPPROJECTION_hemisphere * MAPPROJECTION_LC_lat2 * D2R

    ! calc conformal factor c
    MAPPROJECTION_LC_c = ( log( sin(lat1rot)        ) - log( sin(lat2rot)        ) ) &
              / ( log( tan(0.5_DP*lat1rot) ) - log( tan(0.5_DP*lat2rot) ) )

    ! pre-calc factor
    MAPPROJECTION_LC_fact = sin(lat1rot) / MAPPROJECTION_LC_c / tan(0.5_DP*lat1rot)**MAPPROJECTION_LC_c

    ! calc (x,y) at pole point
    dlon = 0.0_DP

    latrot = 0.5_DP*PI - MAPPROJECTION_hemisphere * MAPPROJECTION_basepoint_lat * D2R

    dist = MAPPROJECTION_LC_fact * RADIUS * tan(0.5_DP*latrot)**MAPPROJECTION_LC_c

    MAPPROJECTION_pole_x = MAPPROJECTION_basepoint_x -                   dist * sin(MAPPROJECTION_LC_c*dlon)
    MAPPROJECTION_pole_y = MAPPROJECTION_basepoint_y + MAPPROJECTION_hemisphere * dist * cos(MAPPROJECTION_LC_c*dlon)

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_LambertConformal_setup",*) 'MAPPROJECTION_LC_lat1    = ', MAPPROJECTION_LC_lat1
    LOG_INFO("MAPPROJECTION_LambertConformal_setup",*) 'MAPPROJECTION_LC_lat2    = ', MAPPROJECTION_LC_lat2
    LOG_INFO("MAPPROJECTION_LambertConformal_setup",*) 'MAPPROJECTION_hemisphere = ', MAPPROJECTION_hemisphere
    LOG_INFO("MAPPROJECTION_LambertConformal_setup",*) 'MAPPROJECTION_LC_c       = ', MAPPROJECTION_LC_c
    LOG_INFO("MAPPROJECTION_LambertConformal_setup",*) 'MAPPROJECTION_LC_fact    = ', MAPPROJECTION_LC_fact
    LOG_INFO("MAPPROJECTION_LambertConformal_setup",*) 'MAPPROJECTION_pole_x     = ', MAPPROJECTION_pole_x
    LOG_INFO("MAPPROJECTION_LambertConformal_setup",*) 'MAPPROJECTION_pole_y     = ', MAPPROJECTION_pole_y

    MAPPROJECTION_mapping = "lambert_conformal_conic"
    MAPPROJECTION_standard_parallel(:) = (/ MAPPROJECTION_LC_lat1, MAPPROJECTION_LC_lat2 /)
    MAPPROJECTION_longitude_of_central_meridian = MAPPROJECTION_basepoint_lon
    MAPPROJECTION_latitude_of_projection_origin = MAPPROJECTION_basepoint_lat
    MAPPROJECTION_false_easting  = MAPPROJECTION_basepoint_x
    MAPPROJECTION_false_northing = MAPPROJECTION_basepoint_y

    return
  end subroutine MAPPROJECTION_LambertConformal_setup

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_LambertConformal_xy2lonlat( &
       x,   &
       y,   &
       lon, &
       lat  )
    implicit none

    real(RP), intent(in)  :: x
    real(RP), intent(in)  :: y
    real(RP), intent(out) :: lon ! [rad]
    real(RP), intent(out) :: lat ! [rad]

    real(DP) :: xx, yy, dist
    !---------------------------------------------------------------------------

    xx =                    ( x - MAPPROJECTION_pole_x ) / RADIUS / MAPPROJECTION_LC_fact
    yy = -MAPPROJECTION_hemisphere * ( y - MAPPROJECTION_pole_y ) / RADIUS / MAPPROJECTION_LC_fact

    dist = sqrt( xx*xx + yy*yy )

    lon = MAPPROJECTION_basepoint_lon * D2R + atan2(xx,yy) / MAPPROJECTION_LC_c

    ! check hemisphere: 1=north, -1=south
    lat = MAPPROJECTION_hemisphere * ( 0.5_DP*PI - 2.0_DP*atan( dist**(1.0_DP/MAPPROJECTION_LC_c) ) )

    return
  end subroutine MAPPROJECTION_LambertConformal_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_LambertConformal_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y

    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - MAPPROJECTION_basepoint_lon * D2R
    if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

    latrot = 0.5_DP*PI - MAPPROJECTION_hemisphere * lat

    dist = MAPPROJECTION_LC_fact * RADIUS * tan(0.5_DP*latrot)**MAPPROJECTION_LC_c

    x = MAPPROJECTION_pole_x +                   dist * sin(MAPPROJECTION_LC_c*dlon)
    y = MAPPROJECTION_pole_y - MAPPROJECTION_hemisphere * dist * cos(MAPPROJECTION_LC_c*dlon)

    return
  end subroutine MAPPROJECTION_LambertConformal_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (lon,lat) -> (m1=m2)
  subroutine MAPPROJECTION_LambertConformal_mapfactor( &
       lat, &
       m1,  &
       m2   )
    implicit none

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)

    real(DP) :: latrot
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       latrot = 0.5_DP*PI - MAPPROJECTION_hemisphere * lat(i,j)

       m1(i,j) = MAPPROJECTION_LC_fact / sin(latrot) * MAPPROJECTION_LC_c * tan(0.5_DP*latrot)**MAPPROJECTION_LC_c
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MAPPROJECTION_LambertConformal_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_LambertConformal_rotcoef_0D( &
       rotc, &
       lon,  &
       lat   )
    implicit none

    real(RP), intent(out) :: rotc(2)
    real(RP), intent(in)  :: lon   ! [rad]
    real(RP), intent(in)  :: lat   ! [rad]

    real(DP) :: dlon
    real(DP) :: alpha
    !---------------------------------------------------------------------------

    dlon = lon - MAPPROJECTION_basepoint_lon * D2R
    if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

    alpha = - MAPPROJECTION_LC_c * dlon * MAPPROJECTION_hemisphere

    rotc(1) = cos( alpha )
    rotc(2) = sin( alpha )

    return
  end subroutine MAPPROJECTION_LambertConformal_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_LambertConformal_rotcoef_2D( &
       rotc, &
       lon,  &
       lat   )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]

    real(DP) :: dlon
    real(DP) :: alpha

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       dlon = lon(i,j) - MAPPROJECTION_basepoint_lon * D2R
       if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
       if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

       alpha = - MAPPROJECTION_LC_c * dlon * MAPPROJECTION_hemisphere

       rotc(i,j,1) = cos( alpha )
       rotc(i,j,2) = sin( alpha )
    enddo
    enddo

    return
  end subroutine MAPPROJECTION_LambertConformal_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection
  subroutine MAPPROJECTION_PolarStereographic_setup
    implicit none

    real(DP) :: lat0
    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    ! check hemisphere: 1=north, -1=south
    MAPPROJECTION_hemisphere = sign(1.0_DP,MAPPROJECTION_PS_lat)

    lat0 = MAPPROJECTION_hemisphere * MAPPROJECTION_PS_lat * D2R

    ! pre-calc factor
    MAPPROJECTION_PS_fact = 1.0_DP + sin(lat0)

    ! calc (x,y) at pole point
    dlon = 0.0_DP

    latrot = 0.5_DP*PI - MAPPROJECTION_hemisphere * MAPPROJECTION_basepoint_lat * D2R

    dist = MAPPROJECTION_PS_fact * RADIUS * tan(0.5_DP*latrot)

    MAPPROJECTION_pole_x = MAPPROJECTION_basepoint_x -                   dist * sin(dlon)
    MAPPROJECTION_pole_y = MAPPROJECTION_basepoint_y + MAPPROJECTION_hemisphere * dist * cos(dlon)

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_PolarStereographic_setup",*) 'MAPPROJECTION_PS_lat1    = ', MAPPROJECTION_PS_lat
    LOG_INFO("MAPPROJECTION_PolarStereographic_setup",*) 'MAPPROJECTION_hemisphere = ', MAPPROJECTION_hemisphere
    LOG_INFO("MAPPROJECTION_PolarStereographic_setup",*) 'MAPPROJECTION_PS_fact    = ', MAPPROJECTION_PS_fact
    LOG_INFO("MAPPROJECTION_PolarStereographic_setup",*) 'MAPPROJECTION_pole_x     = ', MAPPROJECTION_pole_x
    LOG_INFO("MAPPROJECTION_PolarStereographic_setup",*) 'MAPPROJECTION_pole_y     = ', MAPPROJECTION_pole_y

    MAPPROJECTION_mapping = "polar_stereographic"
    MAPPROJECTION_straight_vertical_longitude_from_pole = MAPPROJECTION_basepoint_lon
    MAPPROJECTION_latitude_of_projection_origin = MAPPROJECTION_basepoint_lat
    MAPPROJECTION_standard_parallel(1) = MAPPROJECTION_PS_lat
    MAPPROJECTION_false_easting = MAPPROJECTION_basepoint_x
    MAPPROJECTION_false_northing = MAPPROJECTION_basepoint_y

    return
  end subroutine MAPPROJECTION_PolarStereographic_setup

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_PolarStereographic_xy2lonlat( &
       x,   &
       y,   &
       lon, &
       lat  )
    implicit none

    real(RP), intent(in)  :: x
    real(RP), intent(in)  :: y
    real(RP), intent(out) :: lon ! [rad]
    real(RP), intent(out) :: lat ! [rad]

    real(DP) :: xx, yy, dist
    !---------------------------------------------------------------------------

    xx =                    ( x - MAPPROJECTION_pole_x ) / RADIUS / MAPPROJECTION_PS_fact
    yy = -MAPPROJECTION_hemisphere * ( y - MAPPROJECTION_pole_y ) / RADIUS / MAPPROJECTION_PS_fact

    dist = sqrt( xx*xx + yy*yy )

    lon = MAPPROJECTION_basepoint_lon * D2R + atan2(xx,yy)

    lat = MAPPROJECTION_hemisphere * ( 0.5_DP*PI - 2.0_DP*atan(dist) )

    return
  end subroutine MAPPROJECTION_PolarStereographic_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_PolarStereographic_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y

    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - MAPPROJECTION_basepoint_lon * D2R

    latrot = 0.5_DP*PI - MAPPROJECTION_hemisphere * lat

    dist = MAPPROJECTION_PS_fact * RADIUS * tan(0.5_DP*latrot)

    x = MAPPROJECTION_pole_x +                   dist * sin(dlon)
    y = MAPPROJECTION_pole_y - MAPPROJECTION_hemisphere * dist * cos(dlon)

    return
  end subroutine MAPPROJECTION_PolarStereographic_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (lon,lat) -> (m1=m2)
  subroutine MAPPROJECTION_PolarStereographic_mapfactor( &
       lat, &
       m1,  &
       m2   )
    implicit none

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       m1(i,j) = MAPPROJECTION_PS_fact / ( 1.0_DP + sin(MAPPROJECTION_hemisphere*lat(i,j)) )
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MAPPROJECTION_PolarStereographic_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_PolarStereographic_rotcoef_0D( &
       rotc, &
       lon,  &
       lat   )
    implicit none

    real(RP), intent(out) :: rotc(2)
    real(RP), intent(in)  :: lon   ! [rad]
    real(RP), intent(in)  :: lat   ! [rad]

    real(DP) :: dlon
    real(DP) :: alpha
    !---------------------------------------------------------------------------

    dlon = lon - MAPPROJECTION_basepoint_lon * D2R
    if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

    alpha = - dlon * MAPPROJECTION_hemisphere

    rotc(1) = cos( alpha )
    rotc(2) = sin( alpha )

    return
  end subroutine MAPPROJECTION_PolarStereographic_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_PolarStereographic_rotcoef_2D( &
       rotc, &
       lon,  &
       lat   )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]

    real(DP) :: dlon
    real(DP) :: alpha

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       dlon = lon(i,j) - MAPPROJECTION_basepoint_lon * D2R
       if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
       if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

       alpha = - dlon * MAPPROJECTION_hemisphere

       rotc(i,j,1) = cos( alpha )
       rotc(i,j,2) = sin( alpha )
    enddo
    enddo

    return
  end subroutine MAPPROJECTION_PolarStereographic_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Mercator projection
  subroutine MAPPROJECTION_Mercator_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP) :: lat0
    real(DP) :: latrot, dist
    !---------------------------------------------------------------------------

    lat0 = MAPPROJECTION_M_lat * D2R

    ! pre-calc factor
    MAPPROJECTION_M_fact = cos(lat0)

    if ( MAPPROJECTION_M_fact == 0.0_DP ) then
       LOG_ERROR("MAPPROJECTION_Mercator_setup",*) 'MAPPROJECTION_M_lat cannot be set to pole point! value=', MAPPROJECTION_M_lat
       call PRC_abort
    endif

    ! calc (x,y) at (lon,lat) = (base,0)
    latrot = 0.5_DP*PI - MAPPROJECTION_basepoint_lat * D2R

    dist = 1.0_DP / tan(0.5_DP*latrot)

    MAPPROJECTION_eq_x = MAPPROJECTION_basepoint_x
    MAPPROJECTION_eq_y = MAPPROJECTION_basepoint_y - RADIUS * MAPPROJECTION_M_fact * log(dist)

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_Mercator_setup",*) 'MAPPROJECTION_M_lat      = ', MAPPROJECTION_M_lat
    LOG_INFO("MAPPROJECTION_Mercator_setup",*) 'MAPPROJECTION_M_fact     = ', MAPPROJECTION_M_fact
    LOG_INFO("MAPPROJECTION_Mercator_setup",*) 'MAPPROJECTION_eq_x       = ', MAPPROJECTION_eq_x
    LOG_INFO("MAPPROJECTION_Mercator_setup",*) 'MAPPROJECTION_eq_y       = ', MAPPROJECTION_eq_y

    MAPPROJECTION_mapping = "mercator"
    MAPPROJECTION_longitude_of_projection_origin = MAPPROJECTION_basepoint_lon
    MAPPROJECTION_standard_parallel(1) = MAPPROJECTION_M_lat
    MAPPROJECTION_false_easting = MAPPROJECTION_basepoint_x
    MAPPROJECTION_false_northing = MAPPROJECTION_basepoint_y

    return
  end subroutine MAPPROJECTION_Mercator_setup

  !-----------------------------------------------------------------------------
  !> Mercator projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_Mercator_xy2lonlat( &
       x,   &
       y,   &
       lon, &
       lat  )
    implicit none

    real(RP), intent(in)  :: x
    real(RP), intent(in)  :: y
    real(RP), intent(out) :: lon ! [rad]
    real(RP), intent(out) :: lat ! [rad]

    real(DP) :: xx, yy
    !---------------------------------------------------------------------------

    xx = ( x - MAPPROJECTION_eq_x ) / RADIUS / MAPPROJECTION_M_fact
    yy = ( y - MAPPROJECTION_eq_y ) / RADIUS / MAPPROJECTION_M_fact

    lon = xx + MAPPROJECTION_basepoint_lon * D2R

    lat = 0.5_DP*PI - 2.0_DP*atan( 1.0_DP/exp(yy) )

    return
  end subroutine MAPPROJECTION_Mercator_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Mercator projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_Mercator_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y

    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - MAPPROJECTION_basepoint_lon * D2R

    latrot = 0.5_DP*PI - lat

    dist = 1.0_DP / tan(0.5_DP*latrot)

    x = MAPPROJECTION_eq_x + RADIUS * MAPPROJECTION_M_fact * dlon
    y = MAPPROJECTION_eq_y + RADIUS * MAPPROJECTION_M_fact * log(dist)

    return
  end subroutine MAPPROJECTION_Mercator_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Mercator projection: (lon,lat) -> (m1=m2)
  subroutine MAPPROJECTION_Mercator_mapfactor( &
       lat, &
       m1,  &
       m2   )
    implicit none

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       m1(i,j) = MAPPROJECTION_M_fact / cos(real(lat(i,j),kind=DP))
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MAPPROJECTION_Mercator_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_Mercator_rotcoef_0D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(2)
    !---------------------------------------------------------------------------

    rotc(1) = 1.0_RP
    rotc(2) = 0.0_RP

    return
  end subroutine MAPPROJECTION_Mercator_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_Mercator_rotcoef_2D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MAPPROJECTION_Mercator_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection
  subroutine MAPPROJECTION_EquidistantCylindrical_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP) :: lat0
    !---------------------------------------------------------------------------

    lat0 = MAPPROJECTION_EC_lat * D2R

    ! pre-calc factor
    MAPPROJECTION_EC_fact = cos(lat0)

    if ( MAPPROJECTION_EC_fact == 0.0_DP ) then
       LOG_ERROR("MAPPROJECTION_EquidistantCylindrical_setup",*) 'MAPPROJECTION_EC_lat cannot be set to pole point! value=', MAPPROJECTION_EC_lat
       call PRC_abort
    endif

    MAPPROJECTION_eq_x = MAPPROJECTION_basepoint_x
    MAPPROJECTION_eq_y = MAPPROJECTION_basepoint_y - RADIUS * MAPPROJECTION_basepoint_lat * D2R

    LOG_NEWLINE
    LOG_INFO("MAPPROJECTION_EquidistantCylindrical_setup",*) 'MAPPROJECTION_EC_lat     = ', MAPPROJECTION_EC_lat
    LOG_INFO("MAPPROJECTION_EquidistantCylindrical_setup",*) 'MAPPROJECTION_EC_fact    = ', MAPPROJECTION_EC_fact
    LOG_INFO("MAPPROJECTION_EquidistantCylindrical_setup",*) 'MAPPROJECTION_eq_x       = ', MAPPROJECTION_eq_x
    LOG_INFO("MAPPROJECTION_EquidistantCylindrical_setup",*) 'MAPPROJECTION_eq_y       = ', MAPPROJECTION_eq_y

    MAPPROJECTION_mapping = "equirectangular"
    MAPPROJECTION_standard_parallel(1) = MAPPROJECTION_EC_lat
    MAPPROJECTION_longitude_of_central_meridian = MAPPROJECTION_basepoint_lon
    MAPPROJECTION_false_easting  = MAPPROJECTION_basepoint_x
    MAPPROJECTION_false_northing = MAPPROJECTION_basepoint_y

    return
  end subroutine MAPPROJECTION_EquidistantCylindrical_setup

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (x,y) -> (lon,lat)
  subroutine MAPPROJECTION_EquidistantCylindrical_xy2lonlat( &
       x,   &
       y,   &
       lon, &
       lat  )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in)  :: x
    real(RP), intent(in)  :: y
    real(RP), intent(out) :: lon ! [rad]
    real(RP), intent(out) :: lat ! [rad]

    real(DP) :: xx, yy
    !---------------------------------------------------------------------------

    xx = ( x - MAPPROJECTION_eq_x ) / RADIUS / MAPPROJECTION_EC_fact
    yy = ( y - MAPPROJECTION_eq_y ) / RADIUS

    lon = xx + MAPPROJECTION_basepoint_lon * D2R

    lat = yy

    if ( abs(lat) >  0.5_DP*PI ) then
       LOG_ERROR("MAPPROJECTION_EquidistantCylindrical_xy2lonlat",*) 'Invalid latitude range! value=', lat
       call PRC_abort
    endif

    return
  end subroutine MAPPROJECTION_EquidistantCylindrical_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (lon,lat) -> (x,y)
  subroutine MAPPROJECTION_EquidistantCylindrical_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y

    real(DP) :: dlon
    !---------------------------------------------------------------------------

    dlon = lon - MAPPROJECTION_basepoint_lon * D2R

    x = MAPPROJECTION_eq_x + RADIUS * MAPPROJECTION_EC_fact * dlon
    y = MAPPROJECTION_eq_y + RADIUS * lat

    return
  end subroutine MAPPROJECTION_EquidistantCylindrical_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (lon,lat) -> (m1,m2)
  subroutine MAPPROJECTION_EquidistantCylindrical_mapfactor( &
       lat, &
       m1,  &
       m2   )
    implicit none

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       m1(i,j) = MAPPROJECTION_EC_fact / cos(real(lat(i,j),kind=DP))
       m2(i,j) = 1.0_RP
    enddo
    enddo

    return
  end subroutine MAPPROJECTION_EquidistantCylindrical_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_EquidistantCylindrical_rotcoef_0D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(2)
    !---------------------------------------------------------------------------

    rotc(1) = 1.0_RP
    rotc(2) = 0.0_RP

    return
  end subroutine MAPPROJECTION_EquidistantCylindrical_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MAPPROJECTION_EquidistantCylindrical_rotcoef_2D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MAPPROJECTION_EquidistantCylindrical_rotcoef_2D

end module scale_mapprojection
