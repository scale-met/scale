!-------------------------------------------------------------------------------
!> Module Held-Suarez forcing
!!
!! @par Description
!!          This module contains subroutines for forcing term of Held-Suarez test
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_af_heldsuarez
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: af_heldsuarez_init
  public :: af_heldsuarez

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: sigma_b = 0.7_RP
  real(RP), private, parameter :: Kf      = 1.0_RP / ( 1.0_RP * 86400.0_RP )

  real(RP), private            :: T_eq0   = 315.0_RP ! equatorial maximum temperature [K]
  real(RP), private            :: DT_y    =  60.0_RP ! meridional Equator–pole temperature difference [K]
  real(RP), private, parameter :: Dth_z   =  10.0_RP ! [K]
  real(RP), private, parameter :: Ka      = 1.0_RP / (40.0_RP * 86400.0_RP )
  real(RP), private, parameter :: Ks      = 1.0_RP / ( 4.0_RP * 86400.0_RP )

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_heldsuarez_init( moist_case )
    implicit none

    logical, intent(in) :: moist_case
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[af_heldsuarez]/Category[nhm forcing]'

    if ( moist_case ) then
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Moist H-S testcase by Thatcher and Jablonowski (2016)'

       T_eq0 = 294.0_RP
       DT_y  =  65.0_RP
    else ! original HS parameter
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** Held and Suarez (1994) testcase'

       T_eq0 = 315.0_RP
       DT_y  =  60.0_RP
    endif

    return
  end subroutine af_heldsuarez_init

  !-----------------------------------------------------------------------------
  subroutine af_heldsuarez( &
       ijdim, &
       lat,   &
       pre,   &
       tem,   &
       vx,    &
       vy,    &
       vz,    &
       fvx,   &
       fvy,   &
       fvz,   &
       fe     )
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin, &
       kmax => ADM_kmax
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: lat(ijdim)
    real(RP), intent(in)  :: pre(ijdim,kdim)
    real(RP), intent(in)  :: tem(ijdim,kdim)
    real(RP), intent(in)  :: vx (ijdim,kdim)
    real(RP), intent(in)  :: vy (ijdim,kdim)
    real(RP), intent(in)  :: vz (ijdim,kdim)
    real(RP), intent(out) :: fvx(ijdim,kdim)
    real(RP), intent(out) :: fvy(ijdim,kdim)
    real(RP), intent(out) :: fvz(ijdim,kdim)
    real(RP), intent(out) :: fe (ijdim,kdim)

    real(RP) :: sigma, factor
    real(RP) :: T_eq, coslat, sinlat, ap0, Kt

    integer  :: ij, k
    !---------------------------------------------------------------------------

    do k  = kmin, kmax
    do ij = 1,    ijdim
       ! Rayleigh damping of low-level winds as the boundary-layer scheme for the horizontal velocity
       sigma  = pre(ij,k) / ( 0.5_RP * ( pre(ij,kmin) + pre(ij,kmin-1) ) )
       factor = max( (sigma-sigma_b) / (1.0_RP-sigma_b), 0.0_RP )

       fvx(ij,k) = -Kf * factor * vx(ij,k)
       fvy(ij,k) = -Kf * factor * vy(ij,k)
       fvz(ij,k) = -Kf * factor * vz(ij,k)

       ! Newtonian temperature relaxation as the idealized radiation
       sinlat = abs( sin(lat(ij)) )
       coslat = abs( cos(lat(ij)) )
       ap0    = abs( pre(ij,k) / PRE00 )

       Kt     = Ka + ( Ks - Ka ) * factor * coslat**4

       T_eq   = max( 200.0_RP , ( T_eq0 - DT_y*sinlat*sinlat - Dth_z*log(ap0)*coslat*coslat ) * ap0**(Rdry/CPdry) )

       fe(ij,k) = -Kt * ( tem(ij,k) - T_eq ) * CVdry
    enddo
    enddo

    fvx(:,kmin-1) = 0.0_RP
    fvy(:,kmin-1) = 0.0_RP
    fvz(:,kmin-1) = 0.0_RP
    fe (:,kmin-1) = 0.0_RP
    fvx(:,kmax+1) = 0.0_RP
    fvy(:,kmax+1) = 0.0_RP
    fvz(:,kmax+1) = 0.0_RP
    fe (:,kmax+1) = 0.0_RP

    return
  end subroutine af_heldsuarez

end module mod_af_heldsuarez
