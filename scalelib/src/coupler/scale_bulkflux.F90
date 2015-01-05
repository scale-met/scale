!-------------------------------------------------------------------------------
!> module Surface bulk flux
!!
!! @par Description
!!          calculation of surface bulk flux
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
module scale_bulkflux
  !-----------------------------------------------------------------------------
  !
  !++ used modules
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
  public :: BULKFLUX_setup

  abstract interface
     subroutine bc( &
          Ustar,   & ! (out)
          Tstar,   & ! (out)
          Qstar,   & ! (out)
          Uabs,    & ! (out)
          Ta,      & ! (in)
          Ts,      & ! (in)
          Pa,      & ! (in)
          Ps,      & ! (in)
          Qa,      & ! (in)
          Qs,      & ! (in)
          Ua,      & ! (in)
          Va,      & ! (in)
          Za,      & ! (in)
          PBL,     & ! (in)
          Z0M,     & ! (in)
          Z0H,     & ! (in)
          Z0E      ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(out) :: Ustar ! friction velocity [m/s]
       real(RP), intent(out) :: Tstar ! friction temperature [K]
       real(RP), intent(out) :: Qstar ! friction mixing rate [kg/kg]
       real(RP), intent(out) :: Uabs  ! modified absolute velocity [m/s]

       real(RP), intent(in) :: Ta  ! tempearature at the lowest atmospheric layer [K]
       real(RP), intent(in) :: Ts  ! skin temperature [K]
       real(RP), intent(in) :: Pa  ! pressure at the lowest atmospheric layer [Pa]
       real(RP), intent(in) :: Ps  ! surface pressure [Pa]
       real(RP), intent(in) :: Qa  ! mixing ratio at the lowest atmospheric layer [kg/kg]
       real(RP), intent(in) :: Qs  ! surface mixing ratio [kg/kg]
       real(RP), intent(in) :: Ua  ! zonal wind at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: Va  ! meridional wind at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: Za  ! height at the lowest atmospheric layer [m]
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
  character(len=H_SHORT), private :: BULKFLUX_TYPE = 'B91W01'

  real(RP), private :: BULKFLUX_WSCF = 1.2E+0_RP ! empirical scaling factor of Wstar (Beljaars 1994)

  ! limiter
  real(RP), private :: BULKFLUX_Uabs_min  = 1.0E-2_RP ! minimum of Uabs [m/s]
  real(RP), private :: BULKFLUX_RiB_min   = 1.0E-4_RP ! minimum of RiB [no unit]
  real(RP), private :: BULKFLUX_Wstar_min = 1.0E-4_RP ! minimum of W* [m/s]

contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine BULKFLUX_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_BULKFLUX / &
       BULKFLUX_TYPE,      &
       BULKFLUX_WSCF,      &
       BULKFLUX_Uabs_min,  &
       BULKFLUX_RiB_min,   &
       BULKFLUX_Wstar_min

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '*** Bulk coefficient parameter'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BULKFLUX,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_BULKFLUX. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_BULKFLUX)

    select case( BULKFLUX_TYPE )
    case ( 'U95' )
       BULKFLUX => BULKFLUX_U95
    case ( 'B91W01' )
       BULKFLUX => BULKFLUX_B91W01
    case default
       write(*,*) ' xxx Unsupported TYPE. STOP'
       call PRC_MPIstop
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
      Uabs,    & ! (out)
      Ta,      & ! (in)
      Ts,      & ! (in)
      Pa,      & ! (in)
      Ps,      & ! (in)
      Qa,      & ! (in)
      Qs,      & ! (in)
      Ua,      & ! (in)
      Va,      & ! (in)
      Za,      & ! (in)
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

    ! argument
    real(RP), intent(out) :: Ustar ! friction velocity [m/s]
    real(RP), intent(out) :: Tstar ! friction temperature [K]
    real(RP), intent(out) :: Qstar ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: Uabs  ! modified absolute velocity [m/s]

    real(RP), intent(in) :: Ta  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: Ts  ! skin temperature [K]
    real(RP), intent(in) :: Pa  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: Ps  ! surface pressure [Pa]
    real(RP), intent(in) :: Qa  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Qs  ! surface mixing ratio [kg/kg]
    real(RP), intent(in) :: Ua  ! zonal wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Va  ! meridional wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Za  ! height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(in) :: Z0H ! roughness length of heat [m]
    real(RP), intent(in) :: Z0E ! roughness length of moisture [m]

    ! constant
    real(RP), parameter :: tPrn = 0.74_RP    ! turbulent Prandtl number (Businger et al. 1971)
    real(RP), parameter :: LFb  = 9.4_RP     ! Louis factor b (Louis 1979)
    real(RP), parameter :: LFbp = 4.7_RP     ! Louis factor b' (Louis 1979)
    real(RP), parameter :: LFdm = 7.4_RP     ! Louis factor d for momemtum (Louis 1979)
    real(RP), parameter :: LFdh = 5.3_RP     ! Louis factor d for heat (Louis 1979)

    ! work
    real(RP) :: RiB0, RiB ! bulk Richardson number [no unit]
    real(RP) :: C0 ! initial drag coefficient [no unit]
    real(RP) :: fm, fh, t0th, q0qe
    real(RP) :: Tha, Ths
    !---------------------------------------------------------------------------

    Uabs = max( sqrt( Ua**2 + Va**2 ), BULKFLUX_Uabs_min )
    Tha  = Ta * ( PRE00 / Pa )**( Rdry / CPdry )
    Ths  = Ts * ( PRE00 / Ps )**( Rdry / CPdry )

    RiB0 = GRAV * Za * ( Tha - Ths ) / ( Tha * Uabs**2 )
    if( abs( RiB0 ) < BULKFLUX_RiB_min ) then
      RiB0 = sign( BULKFLUX_RiB_min, RiB0 )
    end if

    C0  = ( KARMAN / log( Za/Z0M ) )**2
    RiB = RiB0

    if( RiB0 >= 0.0_RP ) then
      ! stable condition
      fm = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh = fm
    else
      ! unstable condition
      fm = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0 * sqrt( Za/Z0M ) * sqrt( abs( RiB ) ) )
      fh = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0 * sqrt( Za/Z0M ) * sqrt( abs( RiB ) ) )
    end if

    t0th = 1.0_RP / ( 1.0_RP + log( Z0M/Z0H ) / log( Za/Z0M ) / sqrt( fm ) * fh )
    q0qe = 1.0_RP / ( 1.0_RP + log( Z0M/Z0E ) / log( Za/Z0M ) / sqrt( fm ) * fh )
    RiB  = RiB * t0th

    if( RiB0 >= 0.0_RP ) then
      ! stable condition
      fm = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh = fm
    else
      ! unstable condition
      fm = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0 * sqrt( Za/Z0M ) * sqrt( abs( RiB ) ) )
      fh = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0 * sqrt( Za/Z0M ) * sqrt( abs( RiB ) ) )
    end if

    t0th = 1.0_RP / ( 1.0_RP + log( Z0M/Z0H ) / log( Za/Z0M ) / sqrt( fm ) * fh )
    q0qe = 1.0_RP / ( 1.0_RP + log( Z0M/Z0E ) / log( Za/Z0M ) / sqrt( fm ) * fh )

    Ustar = sqrt( C0 * fm ) * Uabs
    Tstar = C0 * fh * t0th / tPrn * Uabs / Ustar * ( Tha - Ths )
    Qstar = C0 * fh * q0qe / tPrn * Uabs / Ustar * ( Qa  - Qs  )

    return
  end subroutine BULKFLUX_U95

  !-----------------------------------------------------------------------------
  !
  ! refs. Beljaars (1991) and Wilson (2001)
  !
  ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
  ! you should fix the stability functions (fm_unstable, fh_unstable, fm_stable, and fh_stable).
  !
  !-----------------------------------------------------------------------------
  subroutine BULKFLUX_B91W01( &
       Ustar,   & ! (out)
       Tstar,   & ! (out)
       Qstar,   & ! (out)
       Uabs,    & ! (out)
       Ta,      & ! (in)
       Ts,      & ! (in)
       Pa,      & ! (in)
       Ps,      & ! (in)
       Qa,      & ! (in)
       Qs,      & ! (in)
       Ua,      & ! (in)
       Va,      & ! (in)
       Za,      & ! (in)
       PBL,     & ! (in)
       Z0M,     & ! (in)
       Z0H,     & ! (in)
       Z0E      ) ! (in)
    use scale_const, only: &
       GRAV    => CONST_GRAV,    &
       KARMAN  => CONST_KARMAN,  &
       Rdry    => CONST_Rdry,    &
       CPdry   => CONST_CPdry,   &
       EPS     => CONST_EPS,     &
       EPSTvap => CONST_EPSTvap, &
       PRE00   => CONST_PRE00
    implicit none

    ! parameter
    real(RP), parameter :: Pt = 0.95_RP ! turbulent Prandtl number

    ! argument
    real(RP), intent(out) :: Ustar ! friction velocity [m/s]
    real(RP), intent(out) :: Tstar ! friction temperature [K]
    real(RP), intent(out) :: Qstar ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: Uabs  ! modified absolute velocity [m/s]

    real(RP), intent(in) :: Ta  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: Ts  ! skin temperature [K]
    real(RP), intent(in) :: Pa  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: Ps  ! surface pressure [Pa]
    real(RP), intent(in) :: Qa  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Qs  ! surface mixing ratio [kg/kg]
    real(RP), intent(in) :: Ua  ! zonal wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Va  ! meridional wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Za  ! height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(in) :: Z0H ! roughness length of heat [m]
    real(RP), intent(in) :: Z0E ! roughness length of moisture [m]

    ! constant
    integer,  parameter :: nmax    = 100        ! maximum iteration number

    real(RP), parameter :: res_min = 1.0E-4_RP
    real(RP), parameter :: dL      = 1.0E-8_RP  ! delta Obukhov length [m]

    ! variables
    integer :: n

    real(RP) :: res    ! residual
    real(RP) :: dres   ! d(residual)/dL

    real(RP) :: L ! Obukhov length [m]
    real(RP) :: RiB0 ! bulk Richardson number [no unit]
    real(RP) :: Wstar, dWstar ! free convection velocity scale [m/s]

    real(RP) :: UabsUS, UabsS, dUabsUS, dUabsS
    real(RP) :: UstarUS, UstarS, dUstar, dUstarUS, dUstarS
    real(RP) :: TstarUS, TstarS, dTstar, dTstarUS, dTstarS
    real(RP) :: QstarUS, QstarS, dQstar, dQstarUS, dQstarS

    real(RP) :: Tha, Ths
    real(RP) :: sw
    !---------------------------------------------------------------------------

    Uabs = max( sqrt( Ua**2 + Va**2 ), BULKFLUX_Uabs_min )
    Tha  = Ta * ( PRE00 / Pa )**( Rdry / CPdry )
    Ths  = Ts * ( PRE00 / Ps )**( Rdry / CPdry )

    ! initial bulk Richardson number
    RiB0 = GRAV * Za * ( Tha - Ths ) / ( Tha * Uabs**2 )
    if( abs( RiB0 ) < BULKFLUX_RiB_min ) then
      RiB0 = sign( BULKFLUX_RiB_min, RiB0 )
    end if

    ! initial Obukhov length assumed by neutral condition
    L = Za / RiB0 * log(Za/Z0H) / log(Za/Z0M)**2

    ! initial free convection velocity scale
    Wstar  = BULKFLUX_Wstar_min
    dWstar = BULKFLUX_Wstar_min

    do n = 1, nmax
      ! unstable condition
      UabsUS  = max( sqrt( Ua**2 + Va**2 + (BULKFLUX_WSCF*Wstar)**2 ), BULKFLUX_Uabs_min )
      UstarUS = KARMAN / ( log(Za/Z0M) - fm_unstable(Za,L) + fm_unstable(Z0M,L) ) * UabsUS
      TstarUS = KARMAN / ( log(Za/Z0H) - fh_unstable(Za,L) + fh_unstable(Z0H,L) ) / Pt * ( Tha - Ths )
      QstarUS = KARMAN / ( log(Za/Z0E) - fh_unstable(Za,L) + fh_unstable(Z0E,L) ) / Pt * ( Qa  - Qs  )

      ! stable condition
      UabsS  = max( sqrt( Ua**2 + Va**2 ), BULKFLUX_Uabs_min )
      UstarS = KARMAN / ( log(Za/Z0M) - fm_stable(Za,L) + fm_stable(Z0M,L) ) * UabsS
      TstarS = KARMAN / ( log(Za/Z0H) - fh_stable(Za,L) + fh_stable(Z0H,L) ) / Pt * ( Tha - Ths )
      QstarS = KARMAN / ( log(Za/Z0E) - fh_stable(Za,L) + fh_stable(Z0E,L) ) / Pt * ( Qa  - Qs  )

      sw = 0.5_RP - sign( 0.5_RP, L ) ! if unstable, sw = 1

      Uabs  = ( sw ) * UabsUS  + ( 1.0_RP-sw ) * UabsS
      Ustar = ( sw ) * UstarUS + ( 1.0_RP-sw ) * UstarS
      Tstar = ( sw ) * TstarUS + ( 1.0_RP-sw ) * TstarS
      Qstar = ( sw ) * QstarUS + ( 1.0_RP-sw ) * QstarS

      ! update free convection velocity scale (unstable condition only)
      Wstar = ( -PBL * GRAV / Ta * Ustar * Tstar * sw )**(1.0_RP/3.0_RP)

      ! calculate residual
      res = L - Ustar**2 * Ta / ( KARMAN * GRAV * Tstar )

      ! unstable condition
      dUabsUS  = max( sqrt( Ua**2 + Va**2 + (BULKFLUX_WSCF*dWstar)**2 ), BULKFLUX_Uabs_min )
      dUstarUS = KARMAN / ( log(Za/Z0M) - fm_unstable(Za,L+dL) + fm_unstable(Z0M,L+dL) ) * dUabsUS
      dTstarUS = KARMAN / ( log(Za/Z0H) - fh_unstable(Za,L+dL) + fh_unstable(Z0H,L+dL) ) / Pt * ( Tha - Ths )
      dQstarUS = KARMAN / ( log(Za/Z0E) - fh_unstable(Za,L+dL) + fh_unstable(Z0E,L+dL) ) / Pt * ( Qa  - Qs  )
      ! stable condition
      dUabsS  = max( sqrt( Ua**2 + Va**2 ), BULKFLUX_Uabs_min )
      dUstarS = KARMAN / ( log(Za/Z0M) - fm_stable(Za,L+dL) + fm_stable(Z0M,L+dL) ) * dUabsS
      dTstarS = KARMAN / ( log(Za/Z0H) - fh_stable(Za,L+dL) + fh_stable(Z0H,L+dL) ) / Pt * ( Tha - Ths )
      dQstarS = KARMAN / ( log(Za/Z0E) - fh_stable(Za,L+dL) + fh_stable(Z0E,L+dL) ) / Pt * ( Qa  - Qs  )

      sw = 0.5_RP - sign( 0.5_RP, L+dL ) ! if unstable, sw = 1

      dUstar = ( sw ) * dUstarUS + ( 1.0_RP-sw ) * dUstarS
      dTstar = ( sw ) * dTstarUS + ( 1.0_RP-sw ) * dTstarS
      dQstar = ( sw ) * dQstarUS + ( 1.0_RP-sw ) * dQstarS

      ! update d(free convection velocity scale) (unstable condition only)
      dWstar = ( -PBL * GRAV / Ta * dUstar * dTstar * sw )**(1.0_RP/3.0_RP)

      ! calculate d(residual)/dL
      dres = ( (L+dL) - dUstar**2 * Ta / ( KARMAN * GRAV * dTstar ) - res ) / dL

      if( abs( res ) < res_min .or. dres < EPS ) then
        ! finish iteration
        exit
      end if

      ! update Obukhov length
      L = L - res / dres

    end do

    if( n > nmax ) then
      if( IO_L ) write(IO_FID_LOG,*) 'Warning: reach maximum iteration in the function of BULKFLUX_B91W01.'

      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Residual                            [m]     :', res
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- delta Residual                      [m]     :', dres
      if( IO_L ) write(IO_FID_LOG,*) ''
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- air tempearature                    [K]     :', Ta
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- surface temperature                 [K]     :', Ts
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- pressure                            [Pa]    :', Pa
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- surface pressure                    [Pa]    :', Ps
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- water vapor mass ratio              [kg/kg] :', Qa
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- surface water vapor mass ratio      [kg/kg] :', Qs
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- zonal wind                          [m/s]   :', Ua
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- meridional wind                     [m/s]   :', Va
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- cell center height                  [m]     :', Za
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- the top of atmospheric mixing layer [m]     :', PBL
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- roughness length of momentum        [m]     :', Z0M
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- roughness length of heat            [m]     :', Z0H
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- roughness length of moisture        [m]     :', Z0E
      if( IO_L ) write(IO_FID_LOG,*) ''
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- friction velocity                   [m]     :', Ustar
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- friction potential temperature      [K]     :', Tstar
      if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- friction water vapor mass ratio     [kg/kg] :', Qstar
    end if

    return
  end subroutine BULKFLUX_B91W01

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in unstable condition
  function fm_unstable( Z, L )
!    use scale_const, only: &
!      PI  => CONST_PI
    implicit none

    ! argument
    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L

    ! function
    real(RP) :: fm_unstable

    ! works
    real(RP) :: R
!    real(RP) :: r4R
    !---------------------------------------------------------------------------

    R = min( Z/L, 0.0_RP )

    ! Wilson (2001)
    fm_unstable = 3.0_RP * log( ( 1.0_RP + sqrt( 1.0_RP + 3.6_RP * (-R)**(2.0_RP/3.0_RP) ) ) * 0.5_RP )

    ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
    ! you should comment out the above line (Wilson 2001) and uncomment the below lines (Paulson 1974; Dyer 1974).
    !
    !! Paulson (1974); Dyer (1974)
    !r4R = ( 1.0_RP - 16.0_RP * R )**0.25_RP
    !fm_unstable = log( ( 1.0_RP + r4R )**2 * ( 1.0_RP + r4R * r4R ) * 0.125_RP ) - 2.0_RP * atan( r4R ) + PI * 0.5_RP

    return
  end function fm_unstable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in unstable condition
  function fh_unstable( Z, L )
    implicit none

    ! argument
    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L

    ! function
    real(RP) :: fh_unstable

    ! works
    real(RP) :: R
    !---------------------------------------------------------------------------

    R = min( Z/L, 0.0_RP )

    ! Wilson (2001)
    fh_unstable = 3.0_RP * log( ( 1.0_RP + sqrt( 1.0_RP + 7.9_RP * (-R)**(2.0_RP/3.0_RP) ) ) * 0.5_RP )

    ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
    ! you should comment out the above line (Wilson 2001) and uncomment the below lines (Paulson 1974; Dyer 1974).
    !
    !! Paulson (1974); Dyer (1974)
    !fh_unstable = 2.0_RP * log( ( 1.0_RP + sqrt( 1.0_RP - 16.0_RP * R ) ) * 0.5_RP )

    return
  end function fh_unstable

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in stable condition
  function fm_stable( Z, L )
    implicit none

    ! argument
    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L

    ! function
    real(RP) :: fm_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(RP), parameter :: a = 1.0_RP
    real(RP), parameter :: b = 0.667_RP
    real(RP), parameter :: c = 5.0_RP
    real(RP), parameter :: d = 0.35_RP

    ! works
    real(RP) :: R
    !---------------------------------------------------------------------------

    R = max( Z/L, 0.0_RP )

    ! Holtslag and DeBruin (1988)
    fm_stable = - a*R - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fm_stable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in stable condition
  function fh_stable( Z, L )
    implicit none

    ! argument
    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L

    ! function
    real(RP) :: fh_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(RP), parameter :: a = 1.0_RP
    real(RP), parameter :: b = 0.667_RP
    real(RP), parameter :: c = 5.0_RP
    real(RP), parameter :: d = 0.35_RP

    ! works
    real(RP) :: R
    !---------------------------------------------------------------------------

    R = max( Z/L, 0.0_RP )

    ! Beljaars and Holtslag (1991)
    fh_stable = 1.0_RP - ( 1.0_RP + 2.0_RP/3.0_RP * a*R )**1.5_RP - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fh_stable

end module scale_bulkflux
