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
module scale_mapproj
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_const, only: &
     UNDEF  => CONST_UNDEF
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MPRJ_setup
  public :: MPRJ_xy2lonlat
  public :: MPRJ_lonlat2xy
  public :: MPRJ_mapfactor
  public :: MPRJ_rotcoef

  interface MPRJ_rotcoef
     module procedure MPRJ_rotcoef_0D
     module procedure MPRJ_rotcoef_2D
  end interface MPRJ_rotcoef

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public :: MPRJ_basepoint_lon = 135.221_RP ! position of base point (domain center) in real world [deg]
  real(RP), public :: MPRJ_basepoint_lat =  34.653_RP ! position of base point (domain center) in real world [deg]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MPRJ_None_setup
  private :: MPRJ_LambertConformal_setup
  private :: MPRJ_PolarStereographic_setup
  private :: MPRJ_Mercator_setup
  private :: MPRJ_EquidistantCylindrical_setup

  private :: MPRJ_None_xy2lonlat
  private :: MPRJ_LambertConformal_xy2lonlat
  private :: MPRJ_PolarStereographic_xy2lonlat
  private :: MPRJ_Mercator_xy2lonlat
  private :: MPRJ_EquidistantCylindrical_xy2lonlat

  private :: MPRJ_None_lonlat2xy
  private :: MPRJ_LambertConformal_lonlat2xy
  private :: MPRJ_PolarStereographic_lonlat2xy
  private :: MPRJ_Mercator_lonlat2xy
  private :: MPRJ_EquidistantCylindrical_lonlat2xy

  private :: MPRJ_None_mapfactor
  private :: MPRJ_LambertConformal_mapfactor
  private :: MPRJ_PolarStereographic_mapfactor
  private :: MPRJ_Mercator_mapfactor
  private :: MPRJ_EquidistantCylindrical_mapfactor

  private :: MPRJ_None_rotcoef_2D
  private :: MPRJ_LambertConformal_rotcoef_2D
  private :: MPRJ_PolarStereographic_rotcoef_2D
  private :: MPRJ_Mercator_rotcoef_2D
  private :: MPRJ_EquidistantCylindrical_rotcoef_2D

  private :: MPRJ_None_rotcoef_0D
  private :: MPRJ_LambertConformal_rotcoef_0D
  private :: MPRJ_PolarStereographic_rotcoef_0D
  private :: MPRJ_Mercator_rotcoef_0D
  private :: MPRJ_EquidistantCylindrical_rotcoef_0D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: MPRJ_type = 'NONE' !< map projection type
                                               ! 'NONE'
                                               ! 'LC'
                                               ! 'PS'
                                               ! 'MER'
                                               ! 'EC'

  real(DP), private :: MPRJ_hemisphere         ! hemisphere flag: 1=north, -1=south

  real(DP), private :: MPRJ_basepoint_x        ! position of base point in the model [m]
  real(DP), private :: MPRJ_basepoint_y        ! position of base point in the model [m]

  real(DP), private :: MPRJ_pole_x             ! position of north/south pole in the model [m]
  real(DP), private :: MPRJ_pole_y             ! position of north/south pole in the model [m]
  real(DP), private :: MPRJ_eq_x               ! position of equator at the base lon. in the model [m]
  real(DP), private :: MPRJ_eq_y               ! position of equator at the base lon. in the model [m]

  real(DP), private :: MPRJ_rotation =  0.0_DP ! rotation factor (only for 'NONE' type)

  real(DP), private :: MPRJ_LC_lat1  = 30.0_DP ! standard latitude1 for L.C. projection [deg]
  real(DP), private :: MPRJ_LC_lat2  = 60.0_DP ! standard latitude2 for L.C. projection [deg]
  real(DP), private :: MPRJ_LC_c               ! conformal factor
  real(DP), private :: MPRJ_LC_fact            ! pre-calc factor

  real(DP), private :: MPRJ_PS_lat             ! standard latitude1 for P.S. projection [deg]
  real(DP), private :: MPRJ_PS_fact            ! pre-calc factor

  real(DP), private :: MPRJ_M_lat              ! standard latitude1 for Mer. projection [deg]
  real(DP), private :: MPRJ_M_fact             ! pre-calc factor

  real(DP), private :: MPRJ_EC_lat             ! standard latitude1 for E.C. projection [deg]
  real(DP), private :: MPRJ_EC_fact            ! pre-calc factor

  real(DP), private :: PI
  real(DP), private :: D2R
  real(DP), private :: RADIUS

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MPRJ_setup( DOMAIN_CENTER_X, DOMAIN_CENTER_Y )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
     PI_RP     => CONST_PI,     &
     D2R_RP    => CONST_D2R,    &
     RADIUS_RP => CONST_RADIUS
    implicit none

    real(RP), intent(in) :: DOMAIN_CENTER_X !< center position of global domain [m]: x
    real(RP), intent(in) :: DOMAIN_CENTER_Y !< center position of global domain [m]: y

    namelist / PARAM_MAPPROJ / &
       MPRJ_basepoint_lon, &
       MPRJ_basepoint_lat, &
       MPRJ_basepoint_x,   &
       MPRJ_basepoint_y,   &
       MPRJ_type,          &
       MPRJ_rotation,      &
       MPRJ_LC_lat1,       &
       MPRJ_LC_lat2,       &
       MPRJ_PS_lat,        &
       MPRJ_M_lat,         &
       MPRJ_EC_lat

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[MAPPROJ]/Categ[GRID]'

    PI     = real(PI_RP,     kind=DP)
    D2R    = real(D2R_RP,    kind=DP)
    RADIUS = real(RADIUS_RP, kind=DP)

    MPRJ_basepoint_x = UNDEF
    MPRJ_basepoint_y = UNDEF
    MPRJ_PS_lat      = UNDEF
    MPRJ_M_lat       = UNDEF
    MPRJ_EC_lat      = UNDEF

    MPRJ_basepoint_x = DOMAIN_CENTER_X
    MPRJ_basepoint_y = DOMAIN_CENTER_Y

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MAPPROJ,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MAPPROJ. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MAPPROJ)


    if( MPRJ_PS_lat      == UNDEF ) MPRJ_PS_lat      = MPRJ_basepoint_lat
    if( MPRJ_M_lat       == UNDEF ) MPRJ_M_lat       = MPRJ_basepoint_lat
    if( MPRJ_EC_lat      == UNDEF ) MPRJ_EC_lat      = MPRJ_basepoint_lat

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_basepoint_lon:', MPRJ_basepoint_lon
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_basepoint_lat:', MPRJ_basepoint_lat
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_basepoint_x  :', MPRJ_basepoint_x
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_basepoint_y  :', MPRJ_basepoint_y

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Projection type:', trim(MPRJ_type)
    select case(trim(MPRJ_type))
    case('NONE')
       if( IO_L ) write(IO_FID_LOG,*) '=> NO map projection'
       call MPRJ_None_setup
    case('LC')
       if( IO_L ) write(IO_FID_LOG,*) '=> Lambert Conformal projection'
       call MPRJ_LambertConformal_setup
    case('PS')
       if( IO_L ) write(IO_FID_LOG,*) '=> Polar Stereographic projection'
       call MPRJ_PolarStereographic_setup
    case('MER')
       if( IO_L ) write(IO_FID_LOG,*) '=> Mercator projection'
       call MPRJ_Mercator_setup
    case('EC')
       if( IO_L ) write(IO_FID_LOG,*) '=> Equidistant Cylindrical projection'
       call MPRJ_EquidistantCylindrical_setup
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_setup

  !-----------------------------------------------------------------------------
  !> (x,y) -> (lon,lat)
  subroutine MPRJ_xy2lonlat( &
       x,   &
       y,   &
       lon, &
       lat  )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(in)  :: x
    real(RP), intent(in)  :: y
    real(RP), intent(out) :: lon ! [rad]
    real(RP), intent(out) :: lat ! [rad]
    !---------------------------------------------------------------------------

    select case(MPRJ_type)
    case('NONE')
       call MPRJ_None_xy2lonlat( x, y, lon, lat )
    case('LC')
       call MPRJ_LambertConformal_xy2lonlat( x, y, lon, lat )
    case('PS')
       call MPRJ_PolarStereographic_xy2lonlat( x, y, lon, lat )
    case('MER')
       call MPRJ_Mercator_xy2lonlat( x, y, lon, lat )
    case('EC')
       call MPRJ_EquidistantCylindrical_xy2lonlat( x, y, lon, lat )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_xy2lonlat

  !-----------------------------------------------------------------------------
  !> (lon,lat) -> (x,y)
  subroutine MPRJ_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y
    !---------------------------------------------------------------------------

    select case(MPRJ_type)
    case('NONE')
       call MPRJ_None_lonlat2xy( lon, lat, x, y )
    case('LC')
       call MPRJ_LambertConformal_lonlat2xy( lon, lat, x, y )
    case('PS')
       call MPRJ_PolarStereographic_lonlat2xy( lon, lat, x, y )
    case('MER')
       call MPRJ_Mercator_lonlat2xy( lon, lat, x, y )
    case('EC')
       call MPRJ_EquidistantCylindrical_lonlat2xy( lon, lat, x, y )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_lonlat2xy

  !-----------------------------------------------------------------------------
  !> (x,y) -> (lon,lat)
  subroutine MPRJ_mapfactor( &
       lat, &
       m1,  &
       m2   )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(in)  :: lat(IA,JA) ! [rad]
    real(RP), intent(out) :: m1 (IA,JA)
    real(RP), intent(out) :: m2 (IA,JA)
    !---------------------------------------------------------------------------

    select case(MPRJ_type)
    case('NONE')
       call MPRJ_None_mapfactor( lat, m1, m2 )
    case('LC')
       call MPRJ_LambertConformal_mapfactor( lat, m1, m2 )
    case('PS')
       call MPRJ_PolarStereographic_mapfactor( lat, m1, m2 )
    case('MER')
       call MPRJ_Mercator_mapfactor( lat, m1, m2 )
    case('EC')
       call MPRJ_EquidistantCylindrical_mapfactor( lat, m1, m2 )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_mapfactor

  !-----------------------------------------------------------------------------
  !> u(lat,lon) = cos u(x,y) - sin v(x,y)
  !> v(lat,lon) = sin u(x,y) + cos v(x,y)
  subroutine MPRJ_rotcoef_0D( &
       rotc, &
       lon,  &
       lat   )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: rotc(2) !< rotc(:,:,1)->cos, rotc(:,:,2)->sin
    real(RP), intent(in)  :: lon   ! [rad]
    real(RP), intent(in)  :: lat   ! [rad]
    !---------------------------------------------------------------------------

    select case(MPRJ_type)
    case('NONE')
       call MPRJ_None_rotcoef_0D( rotc )
    case('LC')
       call MPRJ_LambertConformal_rotcoef_0D( rotc, lon, lat )
    case('PS')
       call MPRJ_PolarStereographic_rotcoef_0D( rotc, lon, lat )
    case('MER')
       call MPRJ_Mercator_rotcoef_0D( rotc )
    case('EC')
       call MPRJ_EquidistantCylindrical_rotcoef_0D( rotc )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_rotcoef_0D

  !-----------------------------------------------------------------------------
  !> u(lat,lon) = cos u(x,y) - sin v(x,y)
  !> v(lat,lon) = sin u(x,y) + cos v(x,y)
  subroutine MPRJ_rotcoef_2D( &
       rotc, &
       lon,  &
       lat   )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2) !< rotc(:,:,1)->cos, rotc(:,:,2)->sin
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]
    !---------------------------------------------------------------------------

    select case(MPRJ_type)
    case('NONE')
       call MPRJ_None_rotcoef_2D( rotc )
    case('LC')
       call MPRJ_LambertConformal_rotcoef_2D( rotc, lon, lat )
    case('PS')
       call MPRJ_PolarStereographic_rotcoef_2D( rotc, lon, lat )
    case('MER')
       call MPRJ_Mercator_rotcoef_2D( rotc )
    case('EC')
       call MPRJ_EquidistantCylindrical_rotcoef_2D( rotc )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> No projection
  subroutine MPRJ_None_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_rotation:', MPRJ_rotation

    return
  end subroutine MPRJ_None_setup

  !-----------------------------------------------------------------------------
  !> No projection, lon,lat are determined by gnomonic projection: (x,y) -> (lon,lat)
  subroutine MPRJ_None_xy2lonlat( &
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

    gno(1) = ( (y-MPRJ_basepoint_y) * sin(MPRJ_rotation*D2R) &
             + (x-MPRJ_basepoint_x) * cos(MPRJ_rotation*D2R) ) / RADIUS
    gno(2) = ( (y-MPRJ_basepoint_y) * cos(MPRJ_rotation*D2R) &
             - (x-MPRJ_basepoint_x) * sin(MPRJ_rotation*D2R) ) / RADIUS

    rho = sqrt( gno(1) * gno(1) + gno(2) * gno(2) )
    gmm = atan( rho )

    if ( rho == 0.0_DP ) then
       lon = MPRJ_basepoint_lon * D2R
       lat = MPRJ_basepoint_lat * D2R
    else
       lon = MPRJ_basepoint_lon * D2R &
           + atan( gno(1)*sin(gmm) / ( rho   *cos(MPRJ_basepoint_lat*D2R)*cos(gmm) &
                                     - gno(2)*sin(MPRJ_basepoint_lat*D2R)*sin(gmm) ) )
       lat = asin(        sin(MPRJ_basepoint_lat*D2R)*cos(gmm)       &
                 + gno(2)*cos(MPRJ_basepoint_lat*D2R)*sin(gmm) / rho )
    endif

    return
  end subroutine MPRJ_None_xy2lonlat

  !-----------------------------------------------------------------------------
  !> No projection
  subroutine MPRJ_None_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y

    real(DP) :: lat_d
    real(DP) :: gno(2)
    real(DP) :: cos_gmm
    !---------------------------------------------------------------------------
    ! http://mathworld.wolfram.com/GnomonicProjection.html

    lat_d = real(lat,kind=DP)
    cos_gmm = sin(MPRJ_basepoint_lat*D2R) * sin(lat_d) &
            + cos(MPRJ_basepoint_lat*D2R) * cos(lat_d) * cos(lon - MPRJ_basepoint_lon*D2R)

    gno(1) = (cos(lat_d) * sin(lon - MPRJ_basepoint_lon*D2R)) / cos_gmm
    gno(2) = (cos(MPRJ_basepoint_lat*D2R) * sin(lat_d) &
            - sin(MPRJ_basepoint_lat*D2R) * cos(lat_d) * cos(lon - MPRJ_basepoint_lon*D2R)) / cos_gmm

    x = MPRJ_basepoint_x + (gno(1) * cos(MPRJ_rotation*D2R) &
                          - gno(2) * sin(MPRJ_rotation*D2R)) * RADIUS
    y = MPRJ_basepoint_y + (gno(1) * sin(MPRJ_rotation*D2R) &
                          + gno(2) * cos(MPRJ_rotation*D2R)) * RADIUS

    return
  end subroutine MPRJ_None_lonlat2xy

  !-----------------------------------------------------------------------------
  !> No projection: m1=m2=1
  subroutine MPRJ_None_mapfactor( &
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
  end subroutine MPRJ_None_mapfactor

  !-----------------------------------------------------------------------------
  !> No projection:
  subroutine MPRJ_None_rotcoef_0D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(2)
    !---------------------------------------------------------------------------

    rotc(1) = 1.0_RP
    rotc(2) = 0.0_RP

    return
  end subroutine MPRJ_None_rotcoef_0D

  !-----------------------------------------------------------------------------
  !> No projection:
  subroutine MPRJ_None_rotcoef_2D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MPRJ_None_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection
  subroutine MPRJ_LambertConformal_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(DP) :: lat1rot, lat2rot
    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    if ( MPRJ_LC_lat1 >= MPRJ_LC_lat2 ) then
       write(*,*) ' xxx Please set MPRJ_LC_lat1 < MPRJ_LC_lat2 in degree. STOP'
       call PRC_MPIstop
    endif

    ! check hemisphere: 1=north, -1=south
    MPRJ_hemisphere = sign(1.0_DP,MPRJ_LC_lat1+MPRJ_LC_lat2)

    lat1rot = 0.5_DP*PI - MPRJ_hemisphere * MPRJ_LC_lat1 * D2R
    lat2rot = 0.5_DP*PI - MPRJ_hemisphere * MPRJ_LC_lat2 * D2R

    ! calc conformal factor c
    MPRJ_LC_c = ( log( sin(lat1rot) ) - log( sin(lat2rot) ) ) &
              / ( log( tan(0.5_DP*lat1rot) ) - log( tan(0.5_DP*lat2rot) ) )

    ! pre-calc factor
    MPRJ_LC_fact = sin(lat1rot) / MPRJ_LC_c / tan(0.5_DP*lat1rot)**MPRJ_LC_c

    ! calc (x,y) at pole point
    dlon = 0.0_DP

    latrot = 0.5_DP*PI - MPRJ_hemisphere * MPRJ_basepoint_lat * D2R

    dist = MPRJ_LC_fact * RADIUS * tan(0.5_DP*latrot)**MPRJ_LC_c

    MPRJ_pole_x = MPRJ_basepoint_x -                   dist * sin(MPRJ_LC_c*dlon)
    MPRJ_pole_y = MPRJ_basepoint_y + MPRJ_hemisphere * dist * cos(MPRJ_LC_c*dlon)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_LC_lat1   :', MPRJ_LC_lat1
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_LC_lat2   :', MPRJ_LC_lat2
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_hemisphere:', MPRJ_hemisphere
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_LC_c      :', MPRJ_LC_c
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_LC_fact   :', MPRJ_LC_fact
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_pole_x    :', MPRJ_pole_x
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_pole_y    :', MPRJ_pole_y

    return
  end subroutine MPRJ_LambertConformal_setup

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (x,y) -> (lon,lat)
  subroutine MPRJ_LambertConformal_xy2lonlat( &
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

    xx =                    ( x - MPRJ_pole_x ) / RADIUS / MPRJ_LC_fact
    yy = -MPRJ_hemisphere * ( y - MPRJ_pole_y ) / RADIUS / MPRJ_LC_fact

    dist = sqrt( xx*xx + yy*yy )

    lon = MPRJ_basepoint_lon * d2r + atan2(xx,yy) / MPRJ_LC_c
    lon = mod( lon+2.0_DP*PI, 2.0_DP*PI )

    ! check hemisphere: 1=north, -1=south
    lat = MPRJ_hemisphere * ( 0.5_DP*PI - 2.0_DP*atan( dist**(1.0_DP/MPRJ_LC_c) ) )

    return
  end subroutine MPRJ_LambertConformal_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (lon,lat) -> (x,y)
  subroutine MPRJ_LambertConformal_lonlat2xy( &
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

    dlon = lon - MPRJ_basepoint_lon * D2R
    if ( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if ( dlon < -PI ) dlon = dlon + PI*2.0_DP

    latrot = 0.5_DP*PI - lat

    dist = MPRJ_LC_fact * RADIUS * tan(0.5_DP*latrot)**MPRJ_LC_c

    x = MPRJ_pole_x +                   dist * sin(MPRJ_LC_c*dlon)
    y = MPRJ_pole_y - MPRJ_hemisphere * dist * cos(MPRJ_LC_c*dlon)

    return
  end subroutine MPRJ_LambertConformal_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection: (lon,lat) -> (m1=m2)
  subroutine MPRJ_LambertConformal_mapfactor( &
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
       latrot = 0.5_DP*PI - MPRJ_hemisphere * lat(i,j)

       m1(i,j) = MPRJ_LC_fact / sin(latrot) * MPRJ_LC_c * tan(0.5_DP*latrot)**MPRJ_LC_c
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MPRJ_LambertConformal_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_LambertConformal_rotcoef_0D( &
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

    dlon = lon - MPRJ_basepoint_lon * D2R
    if( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if( dlon < -PI ) dlon = dlon + PI*2.0_DP
    alpha = - MPRJ_LC_c * dlon * MPRJ_hemisphere
    rotc(1) = cos( alpha )
    rotc(2) = sin( alpha )

    return
  end subroutine MPRJ_LambertConformal_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MPRJ_LambertConformal_rotcoef_2D( &
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
       dlon = lon(i,j) - MPRJ_basepoint_lon * D2R
       if( dlon >  PI ) dlon = dlon - PI*2.0_DP
       if( dlon < -PI ) dlon = dlon + PI*2.0_DP
       alpha = - MPRJ_LC_c * dlon * MPRJ_hemisphere
       rotc(i,j,1) = cos( alpha )
       rotc(i,j,2) = sin( alpha )
    enddo
    enddo

    return
  end subroutine MPRJ_LambertConformal_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection
  subroutine MPRJ_PolarStereographic_setup
    implicit none

    real(DP) :: lat0
    real(DP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    ! check hemisphere: 1=north, -1=south
    MPRJ_hemisphere = sign(1.0_DP,MPRJ_PS_lat)

    lat0 = MPRJ_hemisphere * MPRJ_PS_lat * D2R

    ! pre-calc factor
    MPRJ_PS_fact = 1.0_DP + sin(lat0)

    ! calc (x,y) at pole point
    dlon = 0.0_DP

    latrot = 0.5_DP*PI - MPRJ_hemisphere * MPRJ_basepoint_lat * D2R

    dist = MPRJ_PS_fact * RADIUS * tan(0.5_DP*latrot)

    MPRJ_pole_x = MPRJ_basepoint_x -                   dist * sin(dlon)
    MPRJ_pole_y = MPRJ_basepoint_y + MPRJ_hemisphere * dist * cos(dlon)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_PS_lat1   :', MPRJ_PS_lat
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_hemisphere:', MPRJ_hemisphere
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_PS_fact   :', MPRJ_PS_fact
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_pole_x    :', MPRJ_pole_x
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_pole_y    :', MPRJ_pole_y

    return
  end subroutine MPRJ_PolarStereographic_setup

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (x,y) -> (lon,lat)
  subroutine MPRJ_PolarStereographic_xy2lonlat( &
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

    xx =                    ( x - MPRJ_pole_x ) / RADIUS / MPRJ_PS_fact
    yy = -MPRJ_hemisphere * ( y - MPRJ_pole_y ) / RADIUS / MPRJ_PS_fact

    dist = sqrt( xx*xx + yy*yy )

    lon = MPRJ_basepoint_lon * D2R + atan2(xx,yy)
    lon = mod( lon+2.0_DP*PI, 2.0_DP*PI )
    lat = MPRJ_hemisphere * ( 0.5_DP*PI - 2.0_DP*atan(dist) )

    return
  end subroutine MPRJ_PolarStereographic_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (lon,lat) -> (x,y)
  subroutine MPRJ_PolarStereographic_lonlat2xy( &
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

    dlon = lon - MPRJ_basepoint_lon * D2R

    latrot = 0.5_DP*PI - lat

    dist = MPRJ_PS_fact * RADIUS * tan(0.5_DP*latrot)

    x = MPRJ_pole_x +                   dist * sin(dlon)
    y = MPRJ_pole_y - MPRJ_hemisphere * dist * cos(dlon)

    return
  end subroutine MPRJ_PolarStereographic_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (lon,lat) -> (m1=m2)
  subroutine MPRJ_PolarStereographic_mapfactor( &
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
       m1(i,j) = MPRJ_PS_fact / ( 1.0_DP + sin(MPRJ_hemisphere*lat(i,j)) )
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MPRJ_PolarStereographic_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_PolarStereographic_rotcoef_0D( &
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

    dlon = lon - MPRJ_basepoint_lon * D2R
    if( dlon >  PI ) dlon = dlon - PI*2.0_DP
    if( dlon < -PI ) dlon = dlon + PI*2.0_DP
    alpha = - dlon * MPRJ_hemisphere
    rotc(1) = cos( alpha )
    rotc(2) = sin( alpha )

    return
  end subroutine MPRJ_PolarStereographic_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MPRJ_PolarStereographic_rotcoef_2D( &
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
       dlon = lon(i,j) - MPRJ_basepoint_lon * D2R
       if( dlon >  PI ) dlon = dlon - PI*2.0_DP
       if( dlon < -PI ) dlon = dlon + PI*2.0_DP
       alpha = - dlon * MPRJ_hemisphere
       rotc(i,j,1) = cos( alpha )
       rotc(i,j,2) = sin( alpha )
    enddo
    enddo

    return
  end subroutine MPRJ_PolarStereographic_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Mercator projection
  subroutine MPRJ_Mercator_setup
    implicit none

    real(DP) :: lat0
    real(DP) :: latrot, dist
    !---------------------------------------------------------------------------

    lat0 = MPRJ_M_lat * D2R

    ! pre-calc factor
    MPRJ_M_fact = cos(lat0)

    ! calc (x,y) at (lon,lat) = (base,0)
    latrot = 0.5_DP*PI - MPRJ_basepoint_lat * D2R

    dist = 1.0_DP / tan(0.5_DP*latrot)

    MPRJ_eq_x = MPRJ_basepoint_x
    MPRJ_eq_y = MPRJ_basepoint_y - RADIUS * MPRJ_M_fact * log(dist)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_M_lat :', MPRJ_M_lat
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_M_fact:', MPRJ_M_fact
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_eq_x  :', MPRJ_eq_x
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_eq_y  :', MPRJ_eq_y

    return
  end subroutine MPRJ_Mercator_setup

  !-----------------------------------------------------------------------------
  !> Mercator projection: (x,y) -> (lon,lat)
  subroutine MPRJ_Mercator_xy2lonlat( &
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

    xx = ( x - MPRJ_eq_x ) / RADIUS / MPRJ_M_fact
    yy = ( y - MPRJ_eq_y ) / RADIUS / MPRJ_M_fact

    lon = xx + MPRJ_basepoint_lon * D2R
    lon = mod( lon+2.0_DP*PI, 2.0_DP*PI )
    lat = 0.5_DP*PI - 2.0_DP*atan( 1.0_DP/exp(yy) )

    return
  end subroutine MPRJ_Mercator_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Mercator projection: (lon,lat) -> (x,y)
  subroutine MPRJ_Mercator_lonlat2xy( &
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

    dlon = lon - MPRJ_basepoint_lon * D2R

    latrot = 0.5_DP*PI - lat

    dist = 1.0_DP / tan(0.5_DP*latrot)

    x = MPRJ_eq_x + RADIUS * MPRJ_M_fact * dlon
    y = MPRJ_eq_y + RADIUS * MPRJ_M_fact * log(dist)

    return
  end subroutine MPRJ_Mercator_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Mercator projection: (lon,lat) -> (m1=m2)
  subroutine MPRJ_Mercator_mapfactor( &
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
       m1(i,j) = MPRJ_M_fact / cos(real(lat(i,j),kind=DP))
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MPRJ_Mercator_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_Mercator_rotcoef_0D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(2)
    !---------------------------------------------------------------------------

    rotc(1) = 1.0_RP
    rotc(2) = 0.0_RP

    return
  end subroutine MPRJ_Mercator_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MPRJ_Mercator_rotcoef_2D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MPRJ_Mercator_rotcoef_2D

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection
  subroutine MPRJ_EquidistantCylindrical_setup
    implicit none

    real(DP) :: lat0
    !---------------------------------------------------------------------------

    lat0 = MPRJ_EC_lat * D2R

    ! pre-calc factor
    MPRJ_EC_fact = cos(lat0)

    MPRJ_eq_x = MPRJ_basepoint_x
    MPRJ_eq_y = MPRJ_basepoint_y - RADIUS * MPRJ_basepoint_lat * D2R

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_EC_lat :', MPRJ_EC_lat
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_EC_fact:', MPRJ_EC_fact
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_eq_x   :', MPRJ_eq_x
    if( IO_L ) write(IO_FID_LOG,*) '+++ MPRJ_eq_y   :', MPRJ_eq_y

    return
  end subroutine MPRJ_EquidistantCylindrical_setup

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (x,y) -> (lon,lat)
  subroutine MPRJ_EquidistantCylindrical_xy2lonlat( &
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

    xx = ( x - MPRJ_eq_x ) / RADIUS / MPRJ_EC_fact
    yy = ( y - MPRJ_eq_y ) / RADIUS

    lon = xx + MPRJ_basepoint_lon * D2R
    lon = mod( lon+2.0_DP*PI, 2.0_DP*PI )
    lat = yy

    return
  end subroutine MPRJ_EquidistantCylindrical_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (lon,lat) -> (x,y)
  subroutine MPRJ_EquidistantCylindrical_lonlat2xy( &
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

    dlon = lon - MPRJ_basepoint_lon * D2R

    x = MPRJ_eq_x + RADIUS * MPRJ_EC_fact * dlon
    y = MPRJ_eq_y + RADIUS * lat

    return
  end subroutine MPRJ_EquidistantCylindrical_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection: (lon,lat) -> (m1,m2)
  subroutine MPRJ_EquidistantCylindrical_mapfactor( &
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
       m1(i,j) = MPRJ_EC_fact / cos(real(lat(i,j),kind=DP))
       m2(i,j) = 1.0_RP
    enddo
    enddo

    return
  end subroutine MPRJ_EquidistantCylindrical_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_EquidistantCylindrical_rotcoef_0D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(2)
    !---------------------------------------------------------------------------

    rotc(1) = 1.0_RP
    rotc(2) = 0.0_RP

    return
  end subroutine MPRJ_EquidistantCylindrical_rotcoef_0D

  !-----------------------------------------------------------------------------
  subroutine MPRJ_EquidistantCylindrical_rotcoef_2D( &
       rotc )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MPRJ_EquidistantCylindrical_rotcoef_2D

end module scale_mapproj
