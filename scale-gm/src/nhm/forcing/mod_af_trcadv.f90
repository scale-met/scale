!-------------------------------------------------------------------------------
!
!+  Module of Updating UVW for Tracer Advection Test
!
!-------------------------------------------------------------------------------
module mod_af_trcadv
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module is for the Dyn Core Test Initialization.
  !
  !
  !++ Current Corresponding Author : R.Yoshida
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      12-10-19   artificial updating  R.Yoshida
  !
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: test11_velocity
  public :: test12_velocity

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: Sp_Unit_East
  private :: Sp_Unit_North

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: a  = 6371220.0_RP        ! Earth's Radius [m]
  real(RP), private, parameter :: Rd = 287.0_RP            ! Ideal gas const dry air [J/kg*K]
  real(RP), private, parameter :: g  = 9.80616_RP         ! Gravity [m/s2]
  real(RP), private, parameter :: cp = 1004.5_RP          ! Specific heat capacity [J/kg*K]
  real(RP), private, parameter :: pi = 3.141592653589793238D0 ! pi

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> original by H.Miura, 20130612 R.Yoshida
  subroutine test11_velocity( &
       time, &
       lon,  &
       lat,  &
       zf,   &
       zh,   &
       vx,   &
       vy,   &
       vz,   &
       w     )
    implicit none

    real(RP),intent(in)  :: time
    real(RP),intent(in)  :: lon
    real(RP),intent(in)  :: lat
    real(RP),intent(in)  :: zf
    real(RP),intent(in)  :: zh
    real(RP),intent(out) :: vx
    real(RP),intent(out) :: vy
    real(RP),intent(out) :: vz
    real(RP),intent(out) :: w

    real(RP), parameter :: tau     = 12.0_RP * 86400.0_RP ! period of motion 12 days
    real(RP), parameter :: u0      = 2.0_RP*pi*a/tau    ! 2 pi a / 12 days
    real(RP), parameter :: k0      = 10.0_RP*a/tau      ! Velocity Magnitude
    real(RP), parameter :: omega0  = 23000.0_RP*pi/tau  ! Velocity Magnitude
    real(RP), parameter :: T0      = 300.0_RP           ! temperature
    real(RP), parameter :: H       = Rd * T0 / g      ! scale height
    real(RP), parameter :: p0      = 100000.0_RP        ! reference pressure (Pa)

    real(RP) :: u ! Zonal wind      [m/s]
    real(RP) :: v ! Meridional wind [m/s]
    real(RP) :: dlon, lonp, bs, height, p, ptop, s, ud
    real(RP) :: east(3), nrth(3)
    !---------------------------------------------------------------------------

    dlon = 2.0_RP * pi * time / tau
    lonp = lon - dlon
    bs   = 0.2_RP

    ! full level
    height = zf
    p      = p0 * exp(-zf/H)
    ptop   = p0 * exp(-12000.0_RP/H)

    s = 1.0_RP + exp( (ptop-p0) / (bs*ptop) ) &
             - exp( (p-p0)    / (bs*ptop) ) &
             - exp( (ptop-p)  / (bs*ptop) )

    ud = (omega0*a) / (bs*ptop) * cos(lonp) * cos(lat)**2 * cos(dlon) &
       * ( -exp( (p-p0)/(bs*ptop) ) + exp( (ptop-p)/(bs*ptop) ) )

    u = k0 * sin(2.0_RP*lat)  * sin(lonp)**2 * cos(0.5_RP*dlon) + u0 * cos(lat) + ud
    v = k0 * sin(2.0_RP*lonp) * cos(lat)     * cos(0.5_RP*dlon)

    east = Sp_Unit_East (lon)
    nrth = Sp_Unit_North(lon,lat)

    vx = east(1) * u + nrth(1) * v
    vy = east(2) * u + nrth(2) * v
    vz = east(3) * u + nrth(3) * v

    height = zh
    p      = p0 * exp(-zh/H)
    ptop   = p0 * exp(-12000.0_RP/H)

    s = 1.0 + exp( (ptop-p0) / (bs*ptop) ) &
            - exp( (p-p0)    / (bs*ptop) ) &
            - exp( (ptop-p)  / (bs*ptop) )

    w = -(Rd*T0) / (g*p) * omega0 * sin(lonp) * cos(lat) * cos(dlon) * s

    return
  end subroutine test11_velocity

  !-----------------------------------------------------------------------------
  !> original by H.Miura, 20130612 R.Yoshida
  subroutine test12_velocity( &
       time, &
       lon,  &
       lat,  &
       zf,   &
       zh,   &
       vx,   &
       vy,   &
       vz,   &
       w     )
    implicit none

    real(RP),intent(in)  :: time
    real(RP),intent(in)  :: lon
    real(RP),intent(in)  :: lat
    real(RP),intent(in)  :: zf
    real(RP),intent(in)  :: zh
    real(RP),intent(out) :: vx
    real(RP),intent(out) :: vy
    real(RP),intent(out) :: vz
    real(RP),intent(out) :: w

    real(RP), parameter :: tau  = 1.0_RP * 86400.0_RP ! period of motion 1 day (in s)
    real(RP), parameter :: u0   = 40.0_RP           ! Zonal velocity magnitude (m/s)
    real(RP), parameter :: w0   = 0.15_RP          ! Vertical velocity magnitude (m/s), changed in v5
    real(RP), parameter :: T0   = 300.0_RP          ! temperature
    real(RP), parameter :: H    = Rd * T0 / g     ! scale height
    real(RP), parameter :: K    = 5.0_RP            ! number of Hadley-like cells
    real(RP), parameter :: ztop = 12000.0_RP        ! model top (m)
    real(RP), parameter :: p0   = 100000.0_RP       ! reference pressure (Pa)

    real(RP) :: u ! Zonal wind      [m/s]
    real(RP) :: v ! Meridional wind [m/s]
    real(RP) :: height, p, rho, t, rho0
    real(RP) :: east(3), nrth(3)
    !---------------------------------------------------------------------------

    t    = T0
    rho0 = p0 / (Rd*t)

    height = zf
    p      = p0 * exp(-zf/H)
    rho    = p / (Rd*t)

    u = u0 * cos(lat)
    v = -(rho0/rho) * a*w0*pi / (K*ztop) * cos(lat) * sin(K*lat) &
      * cos(pi*height/ztop) * cos(pi*time/tau)

    east = Sp_Unit_East (lon)
    nrth = Sp_Unit_North(lon,lat)

    vx = east(1) * u + nrth(1) * v
    vy = east(2) * u + nrth(2) * v
    vz = east(3) * u + nrth(3) * v

    height = zh
    p      = p0 * exp(-zh/H)
    rho    = p / (Rd*t)

    w = (rho0/rho) * (w0/K) * (-2.0_RP*sin(K*lat)*sin(lat)+K*cos(lat)*cos(K*lat)) &
      * sin(pi*height/ztop) * cos(pi*time/tau)

    return
  end subroutine test12_velocity

  !-----------------------------------------------------------------------------
  function Sp_Unit_East( lon ) result( unit_east )
    implicit none

    real(RP), intent(in) :: lon ! [rad]
    real(RP)             :: unit_east(3)
    !---------------------------------------------------------------------------

    unit_east(1) = -sin(lon) ! x-direction
    unit_east(2) =  cos(lon) ! y-direction
    unit_east(3) = 0.0_RP      ! z-direction

    return
  end function Sp_Unit_East

  !-----------------------------------------------------------------------------
  function Sp_Unit_North( lon, lat ) result( unit_north )
    implicit none

    real(RP), intent(in) :: lon, lat ! [rad]
    real(RP)             :: unit_north(3)
    !---------------------------------------------------------------------------

    unit_north(1) = -sin(lat) * cos(lon) ! x-direction
    unit_north(2) = -sin(lat) * sin(lon) ! y-direction
    unit_north(3) =  cos(lat)            ! z-direction

    return
  end function Sp_Unit_North

end module mod_af_trcadv
