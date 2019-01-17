!-------------------------------------------------------------------------------
!> module Surface bulk flux
!!
!! @par Description
!!          calculation of surface bulk flux
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_bulkflux
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
  public :: BULKFLUX_setup

  abstract interface
     subroutine bc( &
          Ustar,   & ! (out)
          Tstar,   & ! (out)
          Qstar,   & ! (out)
          Wstar,   & ! (out)
          Ra,      & ! (out)
          FracU10, & ! (out)
          FracT2,  & ! (out)
          FracQ2,  & ! (out)
          T1,      & ! (in)
          T0,      & ! (in)
          P1,      & ! (in)
          P0,      & ! (in)
          Q1,      & ! (in)
          Q0,      & ! (in)
          Uabs,    & ! (in)
          Z1,      & ! (in)
          PBL,     & ! (in)
          Z0M,     & ! (in)
          Z0H,     & ! (in)
          Z0E      ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(out) :: Ustar   ! friction velocity [m/s]
       real(RP), intent(out) :: Tstar   ! friction temperature [K]
       real(RP), intent(out) :: Qstar   ! friction mixing rate [kg/kg]
       real(RP), intent(out) :: Wstar   ! free convection velocity scale [m/s]
       real(RP), intent(out) :: Ra      ! Aerodynamic resistance (=1/Ce)
       real(RP), intent(out) :: FracU10 ! calculation parameter for U10 [-]
       real(RP), intent(out) :: FracT2  ! calculation parameter for T2 [-]
       real(RP), intent(out) :: FracQ2  ! calculation parameter for Q2 [-]

       real(RP), intent(in) :: T1  ! tempearature at the lowest atmospheric layer [K]
       real(RP), intent(in) :: T0  ! skin temperature [K]
       real(RP), intent(in) :: P1  ! pressure at the lowest atmospheric layer [Pa]
       real(RP), intent(in) :: P0  ! surface pressure [Pa]
       real(RP), intent(in) :: Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
       real(RP), intent(in) :: Q0  ! surface mixing ratio [kg/kg]
       real(RP), intent(in) :: Uabs! absolute velocity at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: Z1  ! height at the lowest atmospheric layer [m]
       real(RP), intent(in) :: PBL ! the top of atmospheric mixing layer [m]
       real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
       real(RP), intent(in) :: Z0H ! roughness length of heat [m]
       real(RP), intent(in) :: Z0E ! roughness length of moisture [m]
     end subroutine bc
  end interface

  procedure(bc), pointer :: BULKFLUX => NULL()
  public :: BULKFLUX

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: BULKFLUX_U95
  private :: BULKFLUX_B91W01
  private :: fm_unstable
  private :: fh_unstable
  private :: fm_stable
  private :: fh_stable

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: BULKFLUX_type = 'B91W01' ! 'U95', 'B94', and 'B91W01'

  logical,  private :: BULKFLUX_NK2018 = .true. !> Nishizawa and Kitamura (2018)

  integer,  private :: BULKFLUX_itr_sa_max = 5  ! maximum iteration number for successive approximation
  integer,  private :: BULKFLUX_itr_nr_max = 10 ! maximum iteration number for Newton-Raphson method

  real(RP), private :: BULKFLUX_err_min = 1.0E-4_RP ! minimum value of error

  real(RP), private :: BULKFLUX_WSCF ! empirical scaling factor of Wstar (Beljaars 1994)

  ! limiter
  real(RP), private :: BULKFLUX_Uabs_min  = 1.0E-2_RP ! minimum of Uabs [m/s]
  real(RP), private :: BULKFLUX_Wstar_min = 1.0E-4_RP ! minimum of W* [m/s]

  logical,  private :: flag_W01

contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine BULKFLUX_setup( dx )
    use scale_prc, only: &
       PRC_abort
    implicit none
    real(RP), intent(in) :: dx

    namelist / PARAM_BULKFLUX / &
       BULKFLUX_type,       &
       BULKFLUX_NK2018,     &
       BULKFLUX_itr_sa_max, &
       BULKFLUX_itr_nr_max, &
       BULKFLUX_err_min,    &
       BULKFLUX_WSCF,       &
       BULKFLUX_Uabs_min,   &
       BULKFLUX_Wstar_min

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("BULKFLUX_setup",*) 'Setup'

    ! WSCF = 1.2 for dx > 1 km (Beljaars 1994)
    ! lim_{dx->0} WSCF = 0  for LES (Fig. 6 Kitamura and Ito 2016 BLM)
    BULKFLUX_WSCF = 1.2_RP * min(dx*1.e-3_RP, 1.0_RP)


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BULKFLUX,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("BULKFLUX_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("BULKFLUX_setup",*) 'Not appropriate names in namelist PARAM_BULKFLUX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_BULKFLUX)

    LOG_NEWLINE
    LOG_INFO("BULKFLUX_setup",*) 'Scheme for surface bulk flux : ', trim(BULKFLUX_type)
    select case(BULKFLUX_type)
    case('U95')
       LOG_INFO_CONT(*) '=> Uno et al.(1995)'
       BULKFLUX => BULKFLUX_U95
    case('B91W01')
       LOG_INFO_CONT(*) '=> Beljaars (1991) and Wilson (2001)'
       FLAG_W01 = .true.
       BULKFLUX => BULKFLUX_B91W01
    case('B94')
       LOG_INFO_CONT(*) '=> Beljaars (1994)'
       FLAG_W01 = .false.
       BULKFLUX => BULKFLUX_B91W01
    case default
       LOG_ERROR("BULKFLUX_setup",*) 'Unsupported BULKFLUX_type. STOP'
       call PRC_abort
    end select

    return
  end subroutine BULKFLUX_setup

  !-----------------------------------------------------------------------------
  ! ref. Uno et al. (1995)
  !-----------------------------------------------------------------------------
  subroutine BULKFLUX_U95( &
      Ustar,   & ! (out)
      Tstar,   & ! (out)
      Qstar,   & ! (out)
      Wstar,   & ! (out)
      Ra,      & ! (out)
      FracU10, & ! (out)
      FracT2,  & ! (out)
      FracQ2,  & ! (out)
      T1,      & ! (in)
      T0,      & ! (in)
      P1,      & ! (in)
      P0,      & ! (in)
      Q1,      & ! (in)
      Q0,      & ! (in)
      Uabs,    & ! (in)
      Z1,      & ! (in)
      PBL,     & ! (in)
      Z0M,     & ! (in)
      Z0H,     & ! (in)
      Z0E      ) ! (in)
    use scale_const, only: &
      GRAV   => CONST_GRAV,   &
      KARMAN => CONST_KARMAN, &
      Rdry   => CONST_Rdry,   &
      CPdry  => CONST_CPdry,  &
      PRE00  => CONST_PRE00
    implicit none

    ! parameter
    real(RP), parameter :: tPrn = 0.74_RP    ! turbulent Prandtl number (Businger et al. 1971)
    real(RP), parameter :: LFb  = 9.4_RP     ! Louis factor b (Louis 1979)
    real(RP), parameter :: LFbp = 4.7_RP     ! Louis factor b' (Louis 1979)
    real(RP), parameter :: LFdm = 7.4_RP     ! Louis factor d for momemtum (Louis 1979)
    real(RP), parameter :: LFdh = 5.3_RP     ! Louis factor d for heat (Louis 1979)

    ! argument
    real(RP), intent(out) :: Ustar   ! friction velocity [m/s]
    real(RP), intent(out) :: Tstar   ! friction temperature [K]
    real(RP), intent(out) :: Qstar   ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: Wstar   ! free convection velocity scale [m/s]
    real(RP), intent(out) :: Ra      ! Aerodynamic resistance (=1/Ce)
    real(RP), intent(out) :: FracU10 ! calculation parameter for U10 [-]
    real(RP), intent(out) :: FracT2  ! calculation parameter for T2 [-]
    real(RP), intent(out) :: FracQ2  ! calculation parameter for Q2 [-]

    real(RP), intent(in) :: T1  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: T0  ! skin temperature [K]
    real(RP), intent(in) :: P1  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: P0  ! surface pressure [Pa]
    real(RP), intent(in) :: Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Q0  ! surface mixing ratio [kg/kg]
    real(RP), intent(in) :: Uabs! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Z1  ! height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(in) :: Z0H ! roughness length of heat [m]
    real(RP), intent(in) :: Z0E ! roughness length of moisture [m]

    ! work
    real(RP) :: UabsW
    real(RP) :: RiB0, RiB ! bulk Richardson number [-]
    real(RP) :: C0Z1, C010, C002 ! initial drag coefficient [-]
    real(RP) :: CmZ1, ChZ1, CqZ1, fmZ1, fhZ1, t0thZ1, q0qeZ1
    real(RP) :: Cm10, fm10
    real(RP) :: Cm02, Ch02, Cq02, fm02, fh02, t0th02, q0qe02
    real(RP) :: TH1, TH0
    real(RP) :: logZ1Z0M, log10Z0m, log02Z0m
    real(RP) :: logZ0MZ0E
    real(RP) :: logZ0MZ0H
    !---------------------------------------------------------------------------

    logZ1Z0m = log( Z1      / Z0M )
    log10Z0m = log( 10.0_RP / Z0M )
    log02Z0m = log( 2.0_RP  / Z0M )

    logZ0MZ0E = max( log( Z0M/Z0E ), 1.0_RP )
    logZ0MZ0H = max( log( Z0M/Z0H ), 1.0_RP )

    UabsW = max( Uabs, BULKFLUX_Uabs_min )
    TH1  = T1 * ( P0 / P1 )**( Rdry / CPdry )
    TH0  = T0

    ! bulk Richardson number
    RiB0 = GRAV * Z1 * ( TH1 - TH0 ) / ( TH1 * UabsW**2 )

    C0Z1 = ( KARMAN / logZ1Z0M )**2
    C010 = ( KARMAN / log10Z0M )**2
    C002 = ( KARMAN / log02Z0M )**2

    RiB = RiB0

    if( RiB0 > 0.0_RP ) then
      ! stable condition
      fmZ1 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fhZ1 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fm10 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fm02 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh02 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
    else
      ! unstable condition
      fmZ1 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0Z1 * sqrt( Z1      / Z0M ) * sqrt( abs( RiB ) ) )
      fhZ1 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0Z1 * sqrt( Z1      / Z0M ) * sqrt( abs( RiB ) ) )
      fm10 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C010 * sqrt( 10.0_RP / Z0M ) * sqrt( abs( RiB ) ) )
      fm02 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C002 * sqrt( 2.0_RP  / Z0M ) * sqrt( abs( RiB ) ) )
      fh02 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C002 * sqrt( 2.0_RP  / Z0M ) * sqrt( abs( RiB ) ) )
    end if

    t0thZ1 = 1.0_RP / ( 1.0_RP + logZ0MZ0H / logZ1Z0M / sqrt( fmZ1 ) * fhZ1 )
    q0qeZ1 = 1.0_RP / ( 1.0_RP + logZ0MZ0E / logZ1Z0M / sqrt( fmZ1 ) * fhZ1 )
    t0th02 = 1.0_RP / ( 1.0_RP + logZ0MZ0H / log02Z0M / sqrt( fm02 ) * fh02 )
    q0qe02 = 1.0_RP / ( 1.0_RP + logZ0MZ0E / log02Z0M / sqrt( fm02 ) * fh02 )

    RiB = RiB * t0thZ1

    if( RiB0 > 0.0_RP ) then
      ! stable condition
      fmZ1 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fhZ1 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fm10 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fm02 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh02 = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
    else
      ! unstable condition
      fmZ1 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0Z1 * sqrt( Z1      / Z0M ) * sqrt( abs( RiB ) ) )
      fhZ1 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0Z1 * sqrt( Z1      / Z0M ) * sqrt( abs( RiB ) ) )
      fm10 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C010 * sqrt( 10.0_RP / Z0M ) * sqrt( abs( RiB ) ) )
      fm02 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C002 * sqrt( 2.0_RP  / Z0M ) * sqrt( abs( RiB ) ) )
      fh02 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C002 * sqrt( 2.0_RP  / Z0M ) * sqrt( abs( RiB ) ) )
    end if

    t0thZ1 = 1.0_RP / ( 1.0_RP + logZ0MZ0H / logZ1Z0M / sqrt( fmZ1 ) * fhZ1 )
    q0qeZ1 = 1.0_RP / ( 1.0_RP + logZ0MZ0E / logZ1Z0M / sqrt( fmZ1 ) * fhZ1 )
    t0th02 = 1.0_RP / ( 1.0_RP + logZ0MZ0H / log02Z0M / sqrt( fm02 ) * fh02 )
    q0qe02 = 1.0_RP / ( 1.0_RP + logZ0MZ0E / log02Z0M / sqrt( fm02 ) * fh02 )

    CmZ1 = C0Z1 * fmZ1
    ChZ1 = C0Z1 * fhZ1 * t0thZ1 / tPrn
    CqZ1 = C0Z1 * fhZ1 * q0qeZ1 / tPrn
    Cm10 = C010 * fm10
    Cm02 = C002 * fm02
    Ch02 = C002 * fh02 * t0th02 / tPrn
    Cq02 = C002 * fh02 * q0qe02 / tPrn

    Ustar = sqrt( CmZ1 ) * UabsW
    Tstar = ChZ1 * UabsW / Ustar * ( TH1 - TH0 )
    Qstar = CqZ1 * UabsW / Ustar * ( Q1  - Q0  )
    Wstar = 0.0_RP

    FracU10 = sqrt( CmZ1 / Cm10 )
    FracT2  = ChZ1 / Ch02 * sqrt( Cm02 / CmZ1 )
    FracQ2  = CqZ1 / Cq02 * sqrt( Cm02 / CmZ1 )

    Ra = 1.0_RP / ( CqZ1 * UabsW )

    return
  end subroutine BULKFLUX_U95

  !-----------------------------------------------------------------------------
  !
  ! refs. Beljaars (1991) and Wilson (2001)
  !
  ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
  ! you should fix the stability functions (fm_unstable, fh_unstable, fm_stable, and fh_stable).
  !
  ! Iteration method: refs. JMA-NHM Description Note II, Mar 2008
  !
  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine BULKFLUX_B91W01( &
      Ustar,   & ! (out)
      Tstar,   & ! (out)
      Qstar,   & ! (out)
      Wstar,   & ! (out)
      Ra,      & ! (out)
      FracU10, & ! (out)
      FracT2,  & ! (out)
      FracQ2,  & ! (out)
      T1,      & ! (in)
      T0,      & ! (in)
      P1,      & ! (in)
      P0,      & ! (in)
      Q1,      & ! (in)
      Q0,      & ! (in)
      Uabs,    & ! (in)
      Z1,      & ! (in)
      PBL,     & ! (in)
      Z0M,     & ! (in)
      Z0H,     & ! (in)
      Z0E      ) ! (in)
    use scale_const, only: &
      GRAV    => CONST_GRAV,    &
      KARMAN  => CONST_KARMAN,  &
      Rdry    => CONST_Rdry,    &
      Rvap    => CONST_Rvap,    &
      CPdry   => CONST_CPdry,   &
      CPvap   => CONST_CPvap,   &
      EPSTvap => CONST_EPSTvap, &
      EPS     => CONST_EPS,     &
      PRE00   => CONST_PRE00
    implicit none

    ! parameter
    real(DP), parameter :: dIL = 1.0E-6_DP ! delta [1/m]

    ! argument
    real(RP), intent(out) :: Ustar   ! friction velocity [m/s]
    real(RP), intent(out) :: Tstar   ! friction temperature [K]
    real(RP), intent(out) :: Qstar   ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: Wstar   ! free convection velocity scale [m/s]
    real(RP), intent(out) :: Ra      ! Aerodynamic resistance (=1/Ce)
    real(RP), intent(out) :: FracU10 ! calculation parameter for U10 [-]
    real(RP), intent(out) :: FracT2  ! calculation parameter for T2 [-]
    real(RP), intent(out) :: FracQ2  ! calculation parameter for Q2 [-]

    real(RP), intent(in) :: T1  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: T0  ! skin temperature [K]
    real(RP), intent(in) :: P1  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: P0  ! surface pressure [Pa]
    real(RP), intent(in) :: Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Q0  ! mixing ratio at surface [kg/kg]
    real(RP), intent(in) :: Uabs! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Z1  ! height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(in) :: Z0H ! roughness length of heat [m]
    real(RP), intent(in) :: Z0E ! roughness length of moisture [m]

    ! work
    integer :: n

    real(DP) :: IL   ! inversed Obukhov length [1/m]
    real(DP) :: res  ! residual
    real(DP) :: dres ! d(residual)/dIL

    real(DP) :: RiB0 ! bulk Richardson number [no unit]

    real(DP) :: UabsC
    real(DP) :: UstarC, dUstarC
    real(DP) :: TstarC, dTstarC
    real(DP) :: QstarC, dQstarC
    real(DP) :: WstarC, dWstar

    real(DP) :: FracU10US, FracU10S, FracU10C
    real(DP) :: FracT2US,  FracT2S,  FracT2C
    real(DP) :: FracQ2US,  FracQ2S,  FracQ2C

    real(DP) :: Rtot, CPtot
    real(DP) :: TH1, TH0
    real(DP) :: TV1, TV0, TVM
    real(DP) :: sw

    real(DP) :: BFLX, dBFLX

    real(DP) :: DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E
    real(DP) :: log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E
    real(DP) :: log_10ovZ0M, log_02ovZ0H, log_02ovZ0E

    real(DP) :: RzM, RzH, RzE
    !---------------------------------------------------------------------------

    ! convert to DP
    if ( BULKFLUX_NK2018 ) then
       DP_Z1  = real( Z1*2.0_RP,  kind=DP )
    else
       DP_Z1  = real( Z1,  kind=DP )
    end if
    DP_Z0M = real( Z0M, kind=DP )
    DP_Z0H = real( Z0H, kind=DP )
    DP_Z0E = real( Z0E, kind=DP )

    UabsC = max( Uabs, BULKFLUX_Uabs_min )

    Rtot  = Rdry  * ( 1.0_DP - Q1 ) + Rvap  * Q1
    CPtot = CPdry * ( 1.0_DP - Q1 ) + CPvap * Q1
    TH1 = T1 * ( P0 / P1 )**( Rtot / CPtot )
    TH0 = T0
    TV1 = TH1 * ( 1.0_DP + EPSTvap * Q1 )
    TV0 = TH0 * ( 1.0_DP + EPSTvap * Q0 )
    TVM = ( TV1 + TV0 ) * 0.5_RP

    RzM = 1.0_DP - DP_Z0M / DP_Z1
    RzH = 1.0_DP - DP_Z0H / DP_Z1
    RzE = 1.0_DP - DP_Z0E / DP_Z1

    ! make log constant
    log_Z1ovZ0M = log( DP_Z1 / DP_Z0M )
    log_Z1ovZ0H = log( DP_Z1 / DP_Z0H )
    log_Z1ovZ0E = log( DP_Z1 / DP_Z0E )

    log_10ovZ0M = log( 10.0_DP / DP_Z0M )
    log_02ovZ0H = log(  2.0_DP / DP_Z0H )
    log_02ovZ0E = log(  2.0_DP / DP_Z0E )

    ! bulk Richardson number at initial step
    RiB0 = GRAV * DP_Z1 * ( TV1 - TV0 ) / ( TVM * UabsC**2 )

    ! inversed Obukhov length at initial step
    IL = RiB0 / DP_Z1 * log_Z1ovZ0M**2 / log_Z1ovZ0H

    ! free convection velocity scale at initial step
    WstarC = BULKFLUX_Wstar_min
    dWstar = BULKFLUX_Wstar_min

    ! Successive approximation
    do n = 1, BULKFLUX_itr_sa_max

       call calc_scales_B91W01( &
            IL, UabsC, TH1, TH0, Q1, Q0, PBL,      & ! (in)
            log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E, & ! (in)
            DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E,         & ! (in)
            RzM, RzH, RzE,                         & ! (in)
            WstarC,                                & ! (inout)
            UstarC, TstarC, QstarC, BFLX           ) ! (out)

       ! estimate the inversed Obukhov length
       IL = - KARMAN * GRAV * BFLX / ( UstarC**3 * TH0 )
    end do


    ! Newton-Raphson method
    do n = 1, BULKFLUX_itr_nr_max

       dWstar = WstarC

       call calc_scales_B91W01( &
            IL, UabsC, TH1, TH0, Q1, Q0, PBL,      & ! (in)
            log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E, & ! (in)
            DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E,         & ! (in)
            RzM, RzH, RzE,                         & ! (in)
            WstarC,                                & ! (inout)
            UstarC, TstarC, QstarC, BFLX           ) ! (out)

       call calc_scales_B91W01( &
            IL+dIL, UabsC, TH1, TH0, Q1, Q0, PBL,  & ! (in)
            log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E, & ! (in)
            DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E,         & ! (in)
            RzM, RzH, RzE,                         & ! (in)
            dWstar,                                & ! (inout)
            dUstarC, dTstarC, dQstarC, dBFLX       ) ! (out)

       res = IL + KARMAN * GRAV * BFLX / ( UstarC**3 * TH0 )

       ! calculate d(residual)/dIL
       dres = 1.0_DP + KARMAN * GRAV / ( TH0 * dIL ) * ( dBFLX / dUstarC**3 - BFLX / UstarC**3 )

       ! stop iteration to prevent numerical error
       if( abs( dres ) < EPS ) exit

       ! convergence test with error levels
       if( abs( res/dres ) < BULKFLUX_err_min ) exit

       ! avoid sign changing
       if( IL * ( IL - res / dres ) < 0.0_RP ) then
          if ( abs(IL) <= EPS ) exit
          IL = sign(EPS, IL)
       end if

       ! update the inversed Obukhov length
       IL = IL - res / dres
    end do

    ! Successive approximation after Newton-Raphson method
    if( .NOT. abs( res/dres ) < BULKFLUX_err_min ) then

       call calc_scales_B91W01( &
            IL, UabsC, TH1, TH0, Q1, Q0, PBL,      & ! (in)
            log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E, & ! (in)
            DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E,         & ! (in)
            RzM, RzH, RzE,                         & ! (in)
            WstarC,                                & ! (inout)
            UstarC, TstarC, QstarC, BFLX           ) ! (out)

       ! estimate the inversed Obukhov length
       IL = - KARMAN * GRAV * BFLX / ( UstarC**3 * TH0 )
    end if


    ! calculate Ustar, Tstar, and Qstar based on IL

    call calc_scales_B91W01( &
         IL, UabsC, TH1, TH0, Q1, Q0, PBL,      & ! (in)
         log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E, & ! (in)
         DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E,         & ! (in)
         RzM, RzH, RzE,                         & ! (in)
         WstarC,                                & ! (inout)
         UstarC, TstarC, QstarC, BFLX           ) ! (out)

    FracU10US = ( log_10ovZ0M - fm_unstable(10.0_DP,IL) + fm_unstable(DP_Z0M,IL) ) &
              / ( log_Z1ovZ0M - fm_unstable(  DP_Z1,IL) + fm_unstable(DP_Z0M,IL) )
    FracT2US  = ( log_02ovZ0H - fm_unstable( 2.0_DP,IL) + fm_unstable(DP_Z0H,IL) ) &
              / ( log_Z1ovZ0H - fm_unstable(  DP_Z1,IL) + fm_unstable(DP_Z0H,IL) )
    FracQ2US  = ( log_02ovZ0E - fm_unstable( 2.0_DP,IL) + fm_unstable(DP_Z0E,IL) ) &
              / ( log_Z1ovZ0E - fm_unstable(  DP_Z1,IL) + fm_unstable(DP_Z0E,IL) )

    FracU10S = ( log_10ovZ0M - fm_stable(10.0_DP,IL) + fm_stable(DP_Z0M,IL) ) &
             / ( log_Z1ovZ0M - fm_stable(  DP_Z1,IL) + fm_stable(DP_Z0M,IL) )
    FracT2S  = ( log_02ovZ0H - fm_stable( 2.0_DP,IL) + fm_stable(DP_Z0H,IL) ) &
             / ( log_Z1ovZ0H - fm_stable(  DP_Z1,IL) + fm_stable(DP_Z0H,IL) )
    FracQ2S  = ( log_02ovZ0E - fm_stable( 2.0_DP,IL) + fm_stable(DP_Z0E,IL) ) &
             / ( log_Z1ovZ0E - fm_stable(  DP_Z1,IL) + fm_stable(DP_Z0E,IL) )

    sw = 0.5_DP - sign( 0.5_DP, IL ) ! if unstable, sw = 1

    FracU10C = ( sw ) * FracU10US + ( 1.0_DP-sw ) * FracU10S
    FracT2C  = ( sw ) * FracT2US  + ( 1.0_DP-sw ) * FracT2S
    FracQ2C  = ( sw ) * FracQ2US  + ( 1.0_DP-sw ) * FracQ2S

    ! revert to RP
    Ustar   = real( UstarC,   kind=RP )
    Tstar   = real( TstarC,   kind=RP )
    Qstar   = real( QstarC,   kind=RP )
    Wstar   = real( WstarC,   kind=RP )
    FracU10 = real( FracU10C, kind=RP )
    FracT2  = real( FracT2C,  kind=RP )
    FracQ2  = real( FracQ2C,  kind=RP )

    Ra = max( ( Q1 - Q0 ) / real(UstarC * QstarC + EPS,kind=RP), EPS )

    return
  end subroutine BULKFLUX_B91W01

  subroutine calc_scales_B91W01( &
       IL, Uabs, TH1, TH0, Q1, Q0, PBL,       &
       log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E, &
       DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E,         &
       RzM, RzH, RzE,                         &
       Wstar,                                 &
       UstarC, TstarC, QstarC, BFLX           )
    use scale_const, only: &
       GRAV    => CONST_GRAV,    &
       KARMAN  => CONST_KARMAN,  &
       EPSTvap => CONST_EPSTvap, &
       EPS     => CONST_EPS
    implicit none
    real(DP), intent(in) :: IL
    real(DP), intent(in) :: Uabs
    real(DP), intent(in) :: TH1
    real(DP), intent(in) :: TH0
    real(RP), intent(in) :: Q1
    real(RP), intent(in) :: Q0
    real(RP), intent(in) :: PBL
    real(DP), intent(in) :: log_Z1ovZ0M
    real(DP), intent(in) :: log_Z1ovZ0H
    real(DP), intent(in) :: log_Z1ovZ0E
    real(DP), intent(in) :: DP_Z1
    real(DP), intent(in) :: DP_Z0M
    real(DP), intent(in) :: DP_Z0H
    real(DP), intent(in) :: DP_Z0E
    real(DP), intent(in) :: RzM
    real(DP), intent(in) :: RzH
    real(DP), intent(in) :: RzE

    real(DP), intent(inout) :: Wstar

    real(DP), intent(out) :: UstarC
    real(DP), intent(out) :: TstarC
    real(DP), intent(out) :: QstarC
    real(DP), intent(out) :: BFLX

    real(DP) :: UabsUS, UabsS
    real(DP) :: UstarUS, UstarS
    real(DP) :: TstarUS, TstarS
    real(DP) :: QstarUS, QstarS
    real(DP) :: denoM, denoH, denoE

    real(DP) :: tmp, sw

    ! unstable condition
    if ( BULKFLUX_NK2018 ) then
       denoM = log_Z1ovZ0M &
             - fmm_unstable(DP_Z1,IL) + fmm_unstable(DP_Z0M,IL) * DP_Z0M / DP_Z1 &
             + RzM * ( fm_unstable(DP_Z0M,IL) - 1.0_DP )
       tmp = fhm_unstable(DP_Z1,IL)
       denoH = log_Z1ovZ0H &
             - tmp + fhm_unstable(DP_Z0H,IL) * DP_Z0H / DP_Z1 &
             + RzH * ( fh_unstable(DP_Z0H,IL) - 1.0_DP )
       denoE = log_Z1ovZ0E &
             - tmp + fhm_unstable(DP_Z0E,IL) * DP_Z0E / DP_Z1 &
             + RzE * ( fh_unstable(DP_Z0E,IL) - 1.0_DP )
    else
       denoM = log_Z1ovZ0M - fm_unstable(DP_Z1,IL) + fm_unstable(DP_Z0M,IL)
       tmp = fh_unstable(DP_Z1,IL)
       denoH = log_Z1ovZ0H - tmp + fh_unstable(DP_Z0H,IL)
       denoE = log_Z1ovZ0E - tmp + fh_unstable(DP_Z0E,IL)
    end if
    UabsUS  = max( sqrt( Uabs**2 + (BULKFLUX_WSCF*Wstar)**2 ), real(BULKFLUX_Uabs_min, kind=DP) )
    UstarUS = KARMAN / denoM * UabsUS
    TstarUS = KARMAN / denoH * ( TH1 - TH0 )
    QstarUS = KARMAN / denoE * ( Q1  - Q0  )

    ! stable condition
    if ( BULKFLUX_NK2018 ) then
       denoM = log_Z1ovZ0M &
             - fmm_stable(DP_Z1,IL) + fmm_stable(DP_Z0M,IL) * DP_Z0M / DP_Z1 &
             + RzM * ( fm_stable(DP_Z0M,IL) - 1.0_DP )
       tmp = fhm_stable(DP_Z1,IL)
       denoH = log_Z1ovZ0H &
             - tmp + fhm_stable(DP_Z0H,IL) * DP_Z0H / DP_Z1 &
             + RzH * ( fh_stable(DP_Z0H,IL) - 1.0_DP )
       denoE = log_Z1ovZ0E &
             - tmp + fhm_stable(DP_Z0E,IL) * DP_Z0E / DP_Z1 &
             + RzE * ( fh_stable(DP_Z0E,IL) - 1.0_DP )
    else
       denoM = log_Z1ovZ0M - fm_stable(DP_Z1,IL) + fm_stable(DP_Z0M,IL)
       tmp = fh_stable(DP_Z1,IL)
       denoH = log_Z1ovZ0H - tmp + fh_stable(DP_Z0H,IL)
       denoE = log_Z1ovZ0E - tmp + fh_stable(DP_Z0E,IL)
    end if
    UabsS  = max( Uabs, real(BULKFLUX_Uabs_min, kind=DP) )
    UstarS = KARMAN / denoM * UabsS
    TstarS = KARMAN / denoH * ( TH1 - TH0 )
    QstarS = KARMAN / denoE * ( Q1  - Q0  )

    sw = 0.5_DP - sign( 0.5_DP, IL ) ! if unstable, sw = 1

    UstarC = ( sw ) * UstarUS + ( 1.0_DP-sw ) * UstarS
    TstarC = ( sw ) * TstarUS + ( 1.0_DP-sw ) * TstarS
    QstarC = ( sw ) * QstarUS + ( 1.0_DP-sw ) * QstarS

    ! estimate buoyancy flux
    BFLX = - UstarC * TstarC * ( 1.0_RP + EPSTvap * Q0 ) - EPSTvap * UstarC * QstarC * TH0

    ! update free convection velocity scale
    tmp = PBL * GRAV / TH0 * BFLX
    sw  = 0.5_DP + sign( 0.5_DP, tmp ) ! if tmp is plus, sw = 1
    Wstar = ( tmp * sw )**( 1.0_DP / 3.0_DP )

    ! avoid zero division with UstarC = 0
    UstarC = sign( max( abs(UstarC), EPS ), UstarC )

    return
  end subroutine calc_scales_B91W01

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in unstable condition
  function fm_unstable( Z, IL )
    use scale_const, only: &
         PI => CONST_PI
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fm_unstable

    ! Wilson (2001)
    real(DP), parameter :: gamma = 3.6_DP

    ! works
    real(DP) :: R
    real(DP) :: r4R
    !---------------------------------------------------------------------------

    R = min( Z * IL, 0.0_DP )

    if ( flag_W01 ) then
       ! Wilson (2001)
       fm_unstable = 3.0_DP * log( ( 1.0_DP + sqrt( 1.0_DP + gamma * (-R)**(2.0_DP/3.0_DP) ) ) * 0.5_DP )
    else
       ! Beljaars and Holtslag (1994), originally Paulson (1974) and Dyer (1974)
       r4R = sqrt( sqrt( 1.0_DP - 16.0_DP * R ) )
       fm_unstable = log( ( 1.0_DP + r4R )**2 * ( 1.0_DP + r4R * r4R ) * 0.125_DP ) - 2.0_DP * atan( r4R ) + PI * 0.5_DP
    end if

    return
  end function fm_unstable
  function fmm_unstable( Z, IL )
    use scale_const, only: &
      PI  => CONST_PI, &
      EPS => CONST_EPS
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fmm_unstable

    ! Wilson (2001)
    real(DP), parameter :: gamma = 3.6_DP

    ! works
    real(DP) :: R
    real(DP) :: f, r3
    real(DP) :: r4R, r2R
    !---------------------------------------------------------------------------

    R = min( Z * IL, 0.0_DP )

    if ( flag_W01 ) then
       ! Wilson (2001)
       r3 = (-R)**(1.0_DP/3.0_DP)
       if ( R > -EPS ) then
          fmm_unstable = 9.0_DP / 20.0_DP * gamma * r3**2
       else
          f = sqrt( 1.0_DP + gamma * r3**2 )
          fmm_unstable = 3.0_DP * log( ( 1.0_DP + f ) * 0.5_DP ) &
                       + 1.5_DP / ( sqrt(gamma)**3 * R) * asinh( sqrt(gamma) * r3 ) &
                       + 1.5_DP * f / ( gamma * r3**2 ) &
                       - 1.0_DP
       end if
    else
       ! Beljaars and Holtslag (1994), originally Paulson (1974) and Dyer (1974)
       if ( R > -EPS ) then ! |R|<EPS, now R < 0
          fmm_unstable = - 2.0_DP * R
       else
          r2R = sqrt( 1.0_DP - 16.0_DP * R )
          r4R = sqrt( r2R )
          fmm_unstable = log( ( 1.0_DP + r4R )**2 * ( 1.0_DP + r2R ) * 0.125_DP ) &
                        - 2.0_DP * atan( r4R ) &
                        + ( 1.0_DP - r4R*r2R ) / ( 12.0_DP * R ) &
                        + PI * 0.5_DP - 1.0_DP
       end if
    end if

    return
  end function fmm_unstable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in unstable condition
  function fh_unstable( Z, IL )
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fh_unstable

    ! Wilson (2001)
    real(DP), parameter :: Pt = 0.95_DP ! turbulent Prandtl number
    real(DP), parameter :: gamma = 7.9_DP

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = min( Z * IL, 0.0_DP )

    if ( flag_W01 ) then
       ! Wilson (2001)
       fh_unstable = 3.0_DP * log( ( 1.0_DP + sqrt( 1.0_DP + gamma * (-R)**(2.0_DP/3.0_DP) ) ) * 0.5_DP ) * PT &
                   + log( Z ) * ( 1.0_DP - Pt )
    else
       ! Beljaars and Holtslag (1994), originally Paulson (1974); Dyer (1974)
       fh_unstable = 2.0_DP * log( ( 1.0_DP + sqrt( 1.0_DP - 16.0_DP * R ) ) * 0.5_DP )
    end if

    return
  end function fh_unstable
  function fhm_unstable( Z, IL )
    use scale_const, only: &
      EPS => CONST_EPS
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fhm_unstable

    ! Wilson (2001)
    real(DP), parameter :: Pt = 0.95_DP ! turbulent Prandtl number
    real(DP), parameter :: gamma = 7.9_DP

    ! works
    real(DP) :: R
    real(DP) :: f, r3
    real(DP) :: r2R
    !---------------------------------------------------------------------------

    R = min( Z * IL, 0.0_DP )

    if ( flag_W01 ) then
       ! Wilson (2001)
       r3 = (-R)**(1.0_DP/3.0_DP)
       if ( R > -EPS ) then
          fhm_unstable = 9.0_DP / 20.0_DP * gamma * r3**2
       else
          f = sqrt( 1.0_DP + gamma * r3**2 )
          fhm_unstable = 3.0_DP * log( ( 1.0_DP + f ) * 0.5_DP ) &
                       + 1.5_DP / ( sqrt(gamma)**3 * R) * asinh( sqrt(gamma) * r3 ) &
                       + 1.5_DP * f / ( gamma * r3**2 ) &
                       - 1.0_DP
          fhm_unstable = fhm_unstable * Pt + ( log( Z ) - 1.0_DP ) * ( 1.0_DP - Pt )
       end if
    else
       ! Beljaars and Holtslag (1994), originally Paulson (1974); Dyer (1974)
       if ( R > -EPS ) then ! |R| < EPS, now R < 0
          fhm_unstable = - 4.0_DP * R
       else
          r2R = sqrt( 1.0_DP - 16.0_DP * R )
          fhm_unstable = 2.0_DP * log( ( 1.0_DP + r2R ) * 0.5_DP ) &
                       + ( 1.0_DP - r2R ) / ( 8.0_DP * R ) &
                       - 1.0_DP
       end if
    end if

    return
  end function fhm_unstable

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in stable condition
  function fm_stable( Z, IL )
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fm_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(DP), parameter :: a = 1.0_DP
    real(DP), parameter :: b = 0.667_DP
    real(DP), parameter :: c = 5.0_DP
    real(DP), parameter :: d = 0.35_DP

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = max( Z * IL, 0.0_DP )

    ! Holtslag and DeBruin (1988)
#if defined(PGI) || defined(SX)
    fm_stable = - a*R - b*( R - c/d )*exp( -min( d*R, 1.E+3_RP ) ) - b*c/d ! apply exp limiter
#else
    fm_stable = - a*R - b*( R - c/d )*exp( -d*R ) - b*c/d
#endif

    return
  end function fm_stable
  function fmm_stable( Z, IL )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fmm_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(DP), parameter :: a = 1.0_DP
    real(DP), parameter :: b = 0.667_DP
    real(DP), parameter :: c = 5.0_DP
    real(DP), parameter :: d = 0.35_DP

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = max( Z * IL, 0.0_DP )

    ! Holtslag and DeBruin (1988)
    if ( R < EPS ) then
       fmm_stable = - 0.5_DP * ( a + b * c + d ) * R
    else
       fmm_stable = b * ( d*R - c + 1.0_DP ) / ( d**2 * R ) &
#if defined(__PGI) || defined(__ES2)
    ! apply exp limiter
                    * exp( -min( d*R, 1.E+3_DP ) ) &
#else
                    * exp( -d*R ) &
#endif
                  - a * R * 0.5_DP &
                  - b * ( c*d*R - c + 1.0_DP ) / ( d**2 * R )
    end if

    return
  end function fmm_stable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in stable condition
  function fh_stable( Z, IL )
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fh_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(DP), parameter :: a = 1.0_DP
    real(DP), parameter :: b = 0.667_DP
    real(DP), parameter :: c = 5.0_DP
    real(DP), parameter :: d = 0.35_DP

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = max( Z * IL, 0.0_DP )

    ! Beljaars and Holtslag (1991)
#if defined(PGI) || defined(SX)
    fh_stable = 1.0_DP - ( 1.0_DP + 2.0_DP/3.0_DP * a*R )**1.5_DP - b*( R - c/d )*exp( -min( d*R, 1.E+3_RP ) ) - b*c/d ! apply exp limiter
#else
    fh_stable = 1.0_DP - ( 1.0_DP + 2.0_DP/3.0_DP * a*R )**1.5_DP - b*( R - c/d )*exp( -d*R ) - b*c/d
#endif

    return
  end function fh_stable
  function fhm_stable( Z, IL )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fhm_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(DP), parameter :: a = 1.0_DP
    real(DP), parameter :: b = 0.667_DP
    real(DP), parameter :: c = 5.0_DP
    real(DP), parameter :: d = 0.35_DP

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = max( Z * IL, 0.0_DP )

    ! Beljaars and Holtslag (1991)
    if ( R < EPS ) then
       fhm_stable = - 0.5_DP * ( a + b*c + b ) * R
    else
       fhm_stable = b * ( d*R - c + 1.0_DP ) / ( d**2 * R ) &
#if defined(__PGI) || defined(__ES2)
                    * exp( -min( d*R, 1.E+3_DP) ) &
#else
                    * exp( -d*R ) &
#endif
                  - 3.0_DP * sqrt( 1.0_DP + 2.0_DP*a*R/3.0_DP )**5 / ( 5.0_DP * a * R ) &
                  - b * ( c*d*R - c + 1.0_DP ) / ( d**2 * R ) &
                  + 0.6_RP / ( a * R ) + 1.0_DP
    end if

    return
  end function fhm_stable

end module scale_bulkflux
