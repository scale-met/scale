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

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MPRJ_None_setup
  private :: MPRJ_LambertConformal_setup
  private :: MPRJ_PolarStereo_setup
!  private :: MPRJ_Mercator_setup

  private :: MPRJ_None_xy2lonlat
  private :: MPRJ_LambertConformal_xy2lonlat
  private :: MPRJ_PolarStereo_xy2lonlat
!  private :: MPRJ_Mercator_xy2lonlat

  private :: MPRJ_None_lonlat2xy
  private :: MPRJ_LambertConformal_lonlat2xy
  private :: MPRJ_PolarStereo_lonlat2xy
!  private :: MPRJ_Mercator_lonlat2xy

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: MPRJ_type = 'NONE' !< map projection type
                                                 ! 'NONE'
                                                 ! 'LC'
                                                 ! 'PS'
                                                 ! 'MER'

  real(RP), private :: MPRJ_hemisphere               ! hemisphere flag: 1=north, -1=south

  real(RP), private :: MPRJ_pole_x                   ! position of north/south pole in the model [m]
  real(RP), private :: MPRJ_pole_y                   ! position of north/south pole in the model [m]

  real(RP), private :: MPRJ_basepoint_lon = 135.2_RP ! position of base point in real world [deg]
  real(RP), private :: MPRJ_basepoint_lat =  34.7_RP ! position of base point in real world [deg]
  real(RP), private :: MPRJ_basepoint_x   =   0.0_RP ! position of base point in the model  [m]
  real(RP), private :: MPRJ_basepoint_y   =   0.0_RP ! position of base point in the model  [m]

  real(RP), private :: MPRJ_rotation      =   0.0_RP ! rotation factor (only for 'NONE' type)

  real(RP), private :: MPRJ_LC_lon                   ! standard longitude for LC projection
  real(RP), private :: MPRJ_LC_lat1 = 30.0_RP        ! standard latitude1 for LC projection
  real(RP), private :: MPRJ_LC_lat2 = 60.0_RP        ! standard latitude2 for LC projection
  real(RP), private :: MPRJ_LC_c                     ! conformal factor
  real(RP), private :: MPRJ_LC_fact                  ! pre-calc factor

  real(RP), private :: MPRJ_PS_lon                   ! standard longitude for PS projection
  real(RP), private :: MPRJ_PS_lat  =  0.0_RP        ! standard latitude1 for PS projection
  real(RP), private :: MPRJ_PS_fact                  ! pre-calc factor

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MPRJ_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_MAPPROJ / &
       MPRJ_type,          &
       MPRJ_rotation,      &
       MPRJ_basepoint_lon, &
       MPRJ_basepoint_lat, &
       MPRJ_basepoint_x,   &
       MPRJ_basepoint_y,   &
       MPRJ_LC_lon,        &
       MPRJ_LC_lat1,       &
       MPRJ_LC_lat2

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[MAPPROJ]/Categ[GRID]'

    MPRJ_LC_lon = CONST_UNDEF
    MPRJ_PS_lon = CONST_UNDEF

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
       call MPRJ_PolarStereo_setup
    case('MER')
       if( IO_L ) write(IO_FID_LOG,*) '=> Mercator projection'
!       call MPRJ_Mercator_setup
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_setup

  !-----------------------------------------------------------------------------
  !> No projection
  subroutine MPRJ_None_setup
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine MPRJ_None_setup

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

    if ( MPRJ_LC_lon == CONST_UNDEF ) then
       MPRJ_LC_lon = MPRJ_basepoint_lon
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
    dlon = ( MPRJ_basepoint_lon - MPRJ_LC_lon ) * D2R

    latrot = 0.5_RP*PI - MPRJ_basepoint_lat * D2R

    dist = RADIUS * MPRJ_LC_fact * tan(0.5_RP*latrot)**MPRJ_LC_c

    MPRJ_pole_x = MPRJ_basepoint_x - MPRJ_hemisphere * dist * sin(MPRJ_LC_c*dlon)
    MPRJ_pole_y = MPRJ_basepoint_y +                   dist * cos(MPRJ_LC_c*dlon)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_hemisphere   :', MPRJ_hemisphere
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_pole_x       :', MPRJ_pole_x
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_pole_y       :', MPRJ_pole_y
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_lon:', MPRJ_basepoint_lon
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_lat:', MPRJ_basepoint_lat
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_x  :', MPRJ_basepoint_x
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_y  :', MPRJ_basepoint_y
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_LC_lon       :', MPRJ_LC_lon
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_LC_lat1      :', MPRJ_LC_lat1
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_LC_lat2      :', MPRJ_LC_lat2
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_LC_c         :', MPRJ_LC_c
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_LC_fact      :', MPRJ_LC_fact

    return
  end subroutine MPRJ_LambertConformal_setup

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection
  subroutine MPRJ_PolarStereo_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP) :: lat1rot
    real(RP) :: dlon, latrot, dist
    !---------------------------------------------------------------------------

    if ( MPRJ_PS_lon == CONST_UNDEF ) then
       MPRJ_PS_lon = MPRJ_basepoint_lon
    endif

    ! check hemisphere: 1=north, -1=south
    MPRJ_hemisphere = sign(1.0_RP,MPRJ_PS_lat)

    lat1rot = 0.5_RP*PI - MPRJ_PS_lat * D2R

    ! pre-calc factor
    MPRJ_PS_fact = 1.0_RP + cos(lat1rot)

    ! calc (x,y) at pole point
    dlon = ( MPRJ_basepoint_lon - MPRJ_PS_lon ) * D2R

    latrot = 0.5_RP*PI - MPRJ_basepoint_lat * D2R

    dist = RADIUS * MPRJ_PS_fact * tan(0.5_RP*latrot)

    MPRJ_pole_x = MPRJ_basepoint_x -                   dist * sin(dlon)
    MPRJ_pole_y = MPRJ_basepoint_y + MPRJ_hemisphere * dist * cos(dlon)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_hemisphere   :', MPRJ_hemisphere
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_pole_x       :', MPRJ_pole_x
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_pole_y       :', MPRJ_pole_y
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_lon:', MPRJ_basepoint_lon
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_lat:', MPRJ_basepoint_lat
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_x  :', MPRJ_basepoint_x
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_basepoint_y  :', MPRJ_basepoint_y
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_PS_lon       :', MPRJ_PS_lon
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_PS_lat1      :', MPRJ_PS_lat
    if( IO_L ) write(IO_FID_LOG,*) 'MPRJ_PS_fact      :', MPRJ_PS_fact

    return
  end subroutine MPRJ_PolarStereo_setup

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
       call MPRJ_PolarStereo_xy2lonlat( x, y, lon, lat )
    case('MER')
!       call MPRJ_Mercator_xy2lonlat( x, y, lon, lat )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_xy2lonlat

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

    gno(1) = ( y * sin(MPRJ_rotation*D2R) + x * cos(MPRJ_rotation*D2R) ) / RADIUS
    gno(2) = ( y * cos(MPRJ_rotation*D2R) - x * sin(MPRJ_rotation*D2R) ) / RADIUS

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

    lon = MPRJ_LC_lon * d2r + atan2(MPRJ_hemisphere*xx,yy) / MPRJ_LC_c
    lon = mod(lon+2.0_RP*PI,2.0_RP*PI)

    ! check hemisphere: 1=north, -1=south
    lat = MPRJ_hemisphere * ( 0.5_RP*PI - 2.0*ATAN( (dist/MPRJ_LC_fact)**(1.0_RP/MPRJ_LC_c) ) )

    return
  end subroutine MPRJ_LambertConformal_xy2lonlat

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (x,y) -> (lon,lat)
  subroutine MPRJ_PolarStereo_xy2lonlat( &
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

    lon = MPRJ_PS_lon * d2r + atan2(MPRJ_hemisphere*xx,yy)
    lon = mod(lon+2.0_RP*PI,2.0_RP*PI)

    ! check hemisphere: 1=north, -1=south
    lat = MPRJ_hemisphere * ( 0.5_RP*PI - 2.0*ATAN(dist/MPRJ_PS_fact) )

    return
  end subroutine MPRJ_PolarStereo_xy2lonlat

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
       call MPRJ_PolarStereo_lonlat2xy( lon, lat, x, y )
    case('MER')
!       call MPRJ_Mercator_lonlat2xy( lon, lat, x, y )
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
    endselect

    return
  end subroutine MPRJ_lonlat2xy

  !-----------------------------------------------------------------------------
  !> No projection
  subroutine MPRJ_None_lonlat2xy( &
       lon, &
       lat, &
       x,   &
       y    )
    implicit none

    real(RP), intent(in)  :: lon ! [rad]
    real(RP), intent(in)  :: lat ! [rad]
    real(RP), intent(out) :: x
    real(RP), intent(out) :: y
    !---------------------------------------------------------------------------

    ! dummy
    x = 0.0_RP
    y = 0.0_RP

    return
  end subroutine MPRJ_None_lonlat2xy

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

    dlon = lon - MPRJ_LC_lon * D2R
    dlon = mod(dlon+2.0_RP*PI,2.0_RP*PI)

    latrot = 0.5_RP*PI - lat

    dist = MPRJ_LC_fact * tan(0.5*latrot)**MPRJ_LC_c

    x = MPRJ_pole_x + MPRJ_hemisphere * dist * sin(MPRJ_LC_c*dlon)
    y = MPRJ_pole_y -                   dist * cos(MPRJ_LC_c*dlon)

    return
  end subroutine MPRJ_LambertConformal_lonlat2xy

  !-----------------------------------------------------------------------------
  !> Polar Stereographic projection: (lon,lat) -> (x,y)
  subroutine MPRJ_PolarStereo_lonlat2xy( &
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

    dlon = lon - MPRJ_LC_lon * D2R
    dlon = mod(dlon+2.0_RP*PI,2.0_RP*PI)

    latrot = 0.5_RP*PI - lat

    dist = MPRJ_LC_fact * tan(0.5*latrot)

    x = MPRJ_pole_x + MPRJ_hemisphere * dist * sin(dlon)
    y = MPRJ_pole_y -                   dist * cos(dlon)

    return
  end subroutine MPRJ_PolarStereo_lonlat2xy

end module scale_mapproj
