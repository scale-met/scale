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
     PI     => CONST_PI,     &
     D2R    => CONST_D2R,    &
     RADIUS => CONST_RADIUS, &
     CONST_UNDEF
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

  private :: MPRJ_None_rotcoef
  private :: MPRJ_LambertConformal_rotcoef
  private :: MPRJ_PolarStereographic_rotcoef
  private :: MPRJ_Mercator_rotcoef
  private :: MPRJ_EquidistantCylindrical_rotcoef

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

  real(RP), private :: MPRJ_hemisphere            ! hemisphere flag: 1=north, -1=south

  real(RP), private :: MPRJ_basepoint_x =  0.0_RP ! position of base point in the model  [m]
  real(RP), private :: MPRJ_basepoint_y =  0.0_RP ! position of base point in the model  [m]

  real(RP), private :: MPRJ_pole_x                ! position of north/south pole in the model [m]
  real(RP), private :: MPRJ_pole_y                ! position of north/south pole in the model [m]
  real(RP), private :: MPRJ_eq_x                  ! position of equator at the base lon. in the model [m]
  real(RP), private :: MPRJ_eq_y                  ! position of equator at the base lon. in the model [m]

  real(RP), private :: MPRJ_rotation    =  0.0_RP ! rotation factor (only for 'NONE' type)

  real(RP), private :: MPRJ_LC_lat1     = 30.0_RP ! standard latitude1 for LC projection
  real(RP), private :: MPRJ_LC_lat2     = 60.0_RP ! standard latitude2 for LC projection
  real(RP), private :: MPRJ_LC_c                  ! conformal factor
  real(RP), private :: MPRJ_LC_fact               ! pre-calc factor

  real(RP), private :: MPRJ_PS_lat      =  0.0_RP ! standard latitude1 for PS projection
  real(RP), private :: MPRJ_PS_fact               ! pre-calc factor

  real(RP), private :: MPRJ_M_lat       =  0.0_RP ! standard latitude1 for Mer. projection
  real(RP), private :: MPRJ_M_fact                ! pre-calc factor

  real(RP), private :: MPRJ_EC_lat      =  0.0_RP ! standard latitude1 for Mer. projection
  real(RP), private :: MPRJ_EC_fact               ! pre-calc factor

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MPRJ_setup( DOMAIN_CENTER_X, DOMAIN_CENTER_Y )
    use scale_process, only: &
       PRC_MPIstop
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

    MPRJ_basepoint_x = DOMAIN_CENTER_X
    MPRJ_basepoint_y = DOMAIN_CENTER_Y

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
  subroutine MPRJ_rotcoef( &
       rotc, &
       lon, lat )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2) !< rotc(:,:,1)->cos, rotc(:,:,2)->sin
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]
    !---------------------------------------------------------------------------

    select case(MPRJ_type)
    case('NONE')
       call MPRJ_None_rotcoef( rotc, lon, lat )
    case('LC')
       call MPRJ_LambertConformal_rotcoef( rotc, lon, lat )
    case('PS')
       call MPRJ_PolarStereographic_rotcoef( rotc, lon, lat )
    case('MER')
       call MPRJ_Mercator_rotcoef( rotc, lon, lat )
    case('EC')
       call MPRJ_EquidistantCylindrical_rotcoef( rotc, lon, lat )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_rotcoef

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

    real(RP) :: gno(2)
    real(RP) :: rho, gmm
    !---------------------------------------------------------------------------

    gno(1) = ( (y-MPRJ_basepoint_y) * sin(MPRJ_rotation*D2R) &
             + (x-MPRJ_basepoint_x) * cos(MPRJ_rotation*D2R) ) / RADIUS
    gno(2) = ( (y-MPRJ_basepoint_y) * cos(MPRJ_rotation*D2R) &
             - (x-MPRJ_basepoint_x) * sin(MPRJ_rotation*D2R) ) / RADIUS

    rho = sqrt( gno(1) * gno(1) + gno(2) * gno(2) )
    gmm = atan( rho )

    if ( rho == 0.0_RP ) then
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
    !---------------------------------------------------------------------------

    ! Todo: must use inverse gnomonic projection
    write(*,*) ' xxx inverse transformation is not implemented for MPRJ_type=NONE. STOP'
    call PRC_MPIstop

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
  subroutine MPRJ_None_rotcoef( &
       rotc, &
       lon, lat )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MPRJ_None_rotcoef

  !-----------------------------------------------------------------------------
  !> Lambert Conformal projection
  subroutine MPRJ_LambertConformal_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP) :: lat1rot, lat2rot
    real(RP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    if ( MPRJ_LC_lat1 >= MPRJ_LC_lat2 ) then
       write(*,*) ' xxx Please set MPRJ_LC_lat1 < MPRJ_LC_lat2 in degree. STOP'
       call PRC_MPIstop
    endif

    ! check hemisphere: 1=north, -1=south
    MPRJ_hemisphere = sign(1.0_RP,MPRJ_LC_lat1+MPRJ_LC_lat2)

    lat1rot = 0.5_RP*PI - MPRJ_LC_lat1 * D2R
    lat2rot = 0.5_RP*PI - MPRJ_LC_lat2 * D2R

    ! calc conformal factor c
    MPRJ_LC_c = ( log( sin(lat1rot) ) - log( sin(lat2rot) ) ) &
              / ( log( tan(0.5_RP*lat1rot) ) - log( tan(0.5_RP*lat2rot) ) )

    ! pre-calc factor
    MPRJ_LC_fact = sin(lat1rot) / MPRJ_LC_c / tan(0.5_RP*lat1rot)**MPRJ_LC_c

    ! calc (x,y) at pole point
    dlon = 0.0_RP

    latrot = 0.5_RP*PI - MPRJ_basepoint_lat * D2R

    dist = MPRJ_LC_fact * tan(0.5_RP*latrot)**MPRJ_LC_c

    MPRJ_pole_x = MPRJ_basepoint_x - MPRJ_hemisphere * RADIUS * dist * sin(MPRJ_LC_c*dlon)
    MPRJ_pole_y = MPRJ_basepoint_y +                   RADIUS * dist * cos(MPRJ_LC_c*dlon)

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

    real(RP) :: xx, yy, dist
    !---------------------------------------------------------------------------

    xx =  MPRJ_hemisphere * x - MPRJ_pole_x
    yy = -MPRJ_hemisphere * y + MPRJ_pole_y

    dist = sqrt( xx*xx + yy*yy ) / RADIUS

    lon = MPRJ_basepoint_lon * d2r + atan2(MPRJ_hemisphere*xx,yy) / MPRJ_LC_c
    lon = mod( lon+2.0_RP*PI, 2.0_RP*PI )

    ! check hemisphere: 1=north, -1=south
    lat = MPRJ_hemisphere * ( 0.5_RP*PI - 2.0_RP*ATAN( (dist/MPRJ_LC_fact)**(1.0_RP/MPRJ_LC_c) ) )

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

    real(RP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - MPRJ_basepoint_lon * D2R
    dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )

    latrot = 0.5_RP*PI - lat

    dist = MPRJ_LC_fact * tan(0.5_RP*latrot)**MPRJ_LC_c

    x = MPRJ_pole_x + MPRJ_hemisphere * RADIUS * dist * sin(MPRJ_LC_c*dlon)
    y = MPRJ_pole_y -                   RADIUS * dist * cos(MPRJ_LC_c*dlon)

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

    real(RP) :: latrot
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       latrot = 0.5_RP*PI - lat(i,j)

       m1(i,j) = MPRJ_LC_fact / sin(latrot) * MPRJ_LC_c * tan(0.5_RP*latrot)**MPRJ_LC_c
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MPRJ_LambertConformal_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_LambertConformal_rotcoef( &
       rotc, &
       lon, lat )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]

    real(RP) :: dlon
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       dlon = lon(i,j) - MPRJ_basepoint_lon * D2R
       dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )
       rotc(i,j,1) = cos( MPRJ_LC_c * dlon ) * MPRJ_hemisphere
       rotc(i,j,2) = sin( MPRJ_LC_c * dlon )
    enddo
    enddo

    return
  end subroutine MPRJ_LambertConformal_rotcoef

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection
  subroutine MPRJ_PolarStereographic_setup
    implicit none

    real(RP) :: lat0
    real(RP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    ! check hemisphere: 1=north, -1=south
    MPRJ_hemisphere = sign(1.0_RP,MPRJ_PS_lat)

    lat0 = MPRJ_PS_lat * D2R

    ! pre-calc factor
    MPRJ_PS_fact = 1.0_RP + sin(lat0)

    ! calc (x,y) at pole point
    dlon = 0.0_RP

    latrot = 0.5_RP*PI - MPRJ_basepoint_lat * D2R

    dist = MPRJ_PS_fact * tan(0.5_RP*latrot)

    MPRJ_pole_x = MPRJ_basepoint_x -                   RADIUS * dist * sin(dlon)
    MPRJ_pole_y = MPRJ_basepoint_y + MPRJ_hemisphere * RADIUS * dist * cos(dlon)

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

    real(RP) :: xx, yy, dist
    !---------------------------------------------------------------------------

    xx =  MPRJ_hemisphere * x - MPRJ_pole_x
    yy = -MPRJ_hemisphere * y + MPRJ_pole_y

    dist = sqrt( xx*xx + yy*yy ) / RADIUS

    lon = MPRJ_basepoint_lon * d2r + atan2(MPRJ_hemisphere*xx,yy)
    lon = mod( lon+2.0_RP*PI, 2.0_RP*PI )

    ! check hemisphere: 1=north, -1=south
    lat = MPRJ_hemisphere * ( 0.5_RP*PI - 2.0_RP*atan(dist/MPRJ_PS_fact) )

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

    real(RP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - MPRJ_basepoint_lon * D2R
    dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )

    latrot = 0.5_RP*PI - lat

    dist = MPRJ_PS_fact * tan(0.5_RP*latrot)

    x = MPRJ_pole_x + MPRJ_hemisphere * RADIUS * dist * sin(dlon)
    y = MPRJ_pole_y -                   RADIUS * dist * cos(dlon)

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
       m1(i,j) = MPRJ_LC_fact / ( 1.0_RP + sin(lat(i,j)) )
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MPRJ_PolarStereographic_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_PolarStereographic_rotcoef( &
       rotc, &
       lon, lat )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]

    real(RP) :: dlon
    integer :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       dlon = lon(i,j) - MPRJ_basepoint_lon * D2R
       dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )
       rotc(i,j,1) = cos( dlon ) * MPRJ_hemisphere
       rotc(i,j,2) = sin( dlon )
    enddo
    enddo

    return
  end subroutine MPRJ_PolarStereographic_rotcoef

  !-----------------------------------------------------------------------------
  !> Mercator projection
  subroutine MPRJ_Mercator_setup
    implicit none

    real(RP) :: lat0
    real(RP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    lat0 = MPRJ_M_lat * D2R

    ! pre-calc factor
    MPRJ_M_fact = cos(lat0)

    dlon = MPRJ_basepoint_lon * D2R
    dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )

    latrot = 0.5_RP*PI - MPRJ_basepoint_lat * D2R

    dist = 1.0_RP / tan(0.5_RP*latrot)

    MPRJ_eq_x = MPRJ_basepoint_x - RADIUS * MPRJ_M_fact * dlon
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

    real(RP) :: xx, yy
    !---------------------------------------------------------------------------

    xx = ( x - MPRJ_eq_x ) / RADIUS / MPRJ_M_fact
    yy = ( y - MPRJ_eq_y ) / RADIUS / MPRJ_M_fact

    lon = xx
    lat = 0.5_RP*PI - 2.0_RP*atan( 1.0_RP/exp(yy) )

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

    real(RP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    dlon = lon - MPRJ_basepoint_lon * D2R
    dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )

    latrot = 0.5_RP*PI - lat

    dist = 1.0_RP / tan(0.5_RP*latrot)

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
       m1(i,j) = MPRJ_M_fact / cos(lat(i,j))
       m2(i,j) = m1(i,j)
    enddo
    enddo

    return
  end subroutine MPRJ_Mercator_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_Mercator_rotcoef( &
       rotc, &
       lon, lat )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MPRJ_Mercator_rotcoef

  !-----------------------------------------------------------------------------
  !> Equidistant Cylindrical projection
  subroutine MPRJ_EquidistantCylindrical_setup
    implicit none

    real(RP) :: lat0
    real(RP) :: dlon
    !---------------------------------------------------------------------------

    lat0 = MPRJ_EC_lat * D2R

    ! pre-calc factor
    MPRJ_EC_fact = cos(lat0)

    dlon = MPRJ_basepoint_lon * D2R
    dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )

    MPRJ_eq_x = MPRJ_basepoint_x - RADIUS * MPRJ_EC_fact * dlon
    MPRJ_eq_y = MPRJ_basepoint_y - RADIUS * MPRJ_basepoint_lat

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

    real(RP) :: xx, yy, dist
    !---------------------------------------------------------------------------

    xx = ( x - MPRJ_eq_x ) / RADIUS / MPRJ_EC_fact
    yy = ( y - MPRJ_eq_y ) / RADIUS

    lon = xx
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

    real(RP) :: dlon
    !---------------------------------------------------------------------------

    dlon = lon - MPRJ_basepoint_lon * D2R
    dlon = mod( dlon+2.0_RP*PI, 2.0_RP*PI )

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
       m1(i,j) = MPRJ_EC_fact / cos(lat(i,j))
       m2(i,j) = 1.0_RP
    enddo
    enddo

    return
  end subroutine MPRJ_EquidistantCylindrical_mapfactor

  !-----------------------------------------------------------------------------
  subroutine MPRJ_EquidistantCylindrical_rotcoef( &
       rotc, &
       lon, lat )
    implicit none

    real(RP), intent(out) :: rotc(IA,JA,2)
    real(RP), intent(in)  :: lon (IA,JA) ! [rad]
    real(RP), intent(in)  :: lat (IA,JA) ! [rad]
    !---------------------------------------------------------------------------

    rotc(:,:,1) = 1.0_RP
    rotc(:,:,2) = 0.0_RP

    return
  end subroutine MPRJ_EquidistantCylindrical_rotcoef

end module scale_mapproj
