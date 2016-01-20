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
          T1,      & ! (in)
          T0,      & ! (in)
          P1,      & ! (in)
          P0,      & ! (in)
          Q1,      & ! (in)
          Q0,      & ! (in)
          U1,      & ! (in)
          V1,      & ! (in)
          Z1,      & ! (in)
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

       real(RP), intent(in) :: T1  ! tempearature at the lowest atmospheric layer [K]
       real(RP), intent(in) :: T0  ! skin temperature [K]
       real(RP), intent(in) :: P1  ! pressure at the lowest atmospheric layer [Pa]
       real(RP), intent(in) :: P0  ! surface pressure [Pa]
       real(RP), intent(in) :: Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
       real(RP), intent(in) :: Q0  ! surface mixing ratio [kg/kg]
       real(RP), intent(in) :: U1  ! zonal wind at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: V1  ! meridional wind at the lowest atmospheric layer [m/s]
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
  character(len=H_SHORT), private :: BULKFLUX_TYPE = 'B91W01'

  integer,  private :: BULKFLUX_itr_max = 100 ! maximum iteration number

  real(RP), private :: BULKFLUX_res_min = 1.0E+0_RP ! minimum value of residual
  real(RP), private :: BULKFLUX_err_min = 1.0E-2_RP ! minimum value of error

  real(RP), private :: BULKFLUX_WSCF = 1.2_RP ! empirical scaling factor of Wstar (Beljaars 1994)

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
       BULKFLUX_itr_max,   &
       BULKFLUX_res_min,   &
       BULKFLUX_err_min,   &
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
      T1,      & ! (in)
      T0,      & ! (in)
      P1,      & ! (in)
      P0,      & ! (in)
      Q1,      & ! (in)
      Q0,      & ! (in)
      U1,      & ! (in)
      V1,      & ! (in)
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
    real(RP), intent(out) :: Ustar ! friction velocity [m/s]
    real(RP), intent(out) :: Tstar ! friction temperature [K]
    real(RP), intent(out) :: Qstar ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: Uabs  ! modified absolute velocity [m/s]

    real(RP), intent(in) :: T1  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: T0  ! skin temperature [K]
    real(RP), intent(in) :: P1  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: P0  ! surface pressure [Pa]
    real(RP), intent(in) :: Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Q0  ! surface mixing ratio [kg/kg]
    real(RP), intent(in) :: U1  ! zonal wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: V1  ! meridional wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Z1  ! height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(in) :: Z0H ! roughness length of heat [m]
    real(RP), intent(in) :: Z0E ! roughness length of moisture [m]

    ! work
    real(RP) :: RiB0, RiB ! bulk Richardson number [no unit]
    real(RP) :: C0 ! initial drag coefficient [no unit]
    real(RP) :: fm, fh, t0th, q0qe
    real(RP) :: TH1, TH0
    real(RP) :: logZ1Z0M
    real(RP) :: logZ0MZ0E
    real(RP) :: logZ0MZ0H
    !---------------------------------------------------------------------------

    logZ1Z0m = log( Z1/Z0M )
    logZ0MZ0E = max( log( Z0M/Z0E ), 1.0_RP )
    logZ0MZ0H = max( log( Z0M/Z0H ), 1.0_RP )

    Uabs = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
    TH1  = T1 * ( PRE00 / P1 )**( Rdry / CPdry )
    TH0  = T0 * ( PRE00 / P0 )**( Rdry / CPdry )

    RiB0 = GRAV * Z1 * ( TH1 - TH0 ) / ( TH1 * Uabs**2 )
    if( abs( RiB0 ) < BULKFLUX_RiB_min ) then
      RiB0 = sign( BULKFLUX_RiB_min, RiB0 )
    end if

    C0  = ( KARMAN / logZ1Z0M )**2
    RiB = RiB0

    if( RiB0 > 0.0_RP ) then
      ! stable condition
      fm = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh = fm
    else
      ! unstable condition
      fm = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0 * sqrt( Z1/Z0M ) * sqrt( abs( RiB ) ) )
      fh = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0 * sqrt( Z1/Z0M ) * sqrt( abs( RiB ) ) )
    end if

    t0th = 1.0_RP / ( 1.0_RP + logZ0MZ0H / logZ1Z0M / sqrt( fm ) * fh )
    q0qe = 1.0_RP / ( 1.0_RP + logZ0MZ0E / logZ1Z0M / sqrt( fm ) * fh )
    RiB  = RiB * t0th

    if( RiB0 > 0.0_RP ) then
      ! stable condition
      fm = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh = fm
    else
      ! unstable condition
      fm = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0 * sqrt( Z1/Z0M ) * sqrt( abs( RiB ) ) )
      fh = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0 * sqrt( Z1/Z0M ) * sqrt( abs( RiB ) ) )
    end if

    t0th = 1.0_RP / ( 1.0_RP + logZ0MZ0H / logZ1Z0M / sqrt( fm ) * fh )
    q0qe = 1.0_RP / ( 1.0_RP + logZ0MZ0E / logZ1Z0M / sqrt( fm ) * fh )

    Ustar = sqrt( C0 * fm ) * Uabs
    Tstar = C0 * fh * t0th / tPrn * Uabs / Ustar * ( TH1 - TH0 )
    Qstar = C0 * fh * q0qe / tPrn * Uabs / Ustar * ( Q1  - Q0  )

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
       RP_Ustar,   & ! (out)
       RP_Tstar,   & ! (out)
       RP_Qstar,   & ! (out)
       RP_Uabs,    & ! (out)
       RP_T1,      & ! (in)
       RP_T0,      & ! (in)
       RP_P1,      & ! (in)
       RP_P0,      & ! (in)
       RP_Q1,      & ! (in)
       RP_Q0,      & ! (in)
       RP_U1,      & ! (in)
       RP_V1,      & ! (in)
       RP_Z1,      & ! (in)
       RP_PBL,     & ! (in)
       RP_Z0M,     & ! (in)
       RP_Z0H,     & ! (in)
       RP_Z0E      ) ! (in)
    use scale_const, only: &
       GRAV    => CONST_GRAV,    &
       KARMAN  => CONST_KARMAN,  &
       Rdry    => CONST_Rdry,    &
       CPdry   => CONST_CPdry,   &
       EPS     => CONST_EPS,     &
       PRE00   => CONST_PRE00
    implicit none

    ! parameter
    real(DP), parameter :: dL = 1.0D-6 ! delta Obukhov length [m]

    real(DP), parameter :: Pt = 0.95D+0 ! turbulent Prandtl number

    ! argument
    real(RP), intent(out) :: RP_Ustar ! friction velocity [m/s]
    real(RP), intent(out) :: RP_Tstar ! friction temperature [K]
    real(RP), intent(out) :: RP_Qstar ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: RP_Uabs  ! modified absolute velocity [m/s]

    real(RP), intent(in) :: RP_T1  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: RP_T0  ! skin temperature [K]
    real(RP), intent(in) :: RP_P1  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: RP_P0  ! surface pressure [Pa]
    real(RP), intent(in) :: RP_Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: RP_Q0  ! surface mixing ratio [kg/kg]
    real(RP), intent(in) :: RP_U1  ! zonal wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: RP_V1  ! meridional wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: RP_Z1  ! height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: RP_PBL ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: RP_Z0M ! roughness length of momentum [m]
    real(RP), intent(in) :: RP_Z0H ! roughness length of heat [m]
    real(RP), intent(in) :: RP_Z0E ! roughness length of moisture [m]

    ! work
    integer :: n

    real(DP) :: Ustar ! friction velocity [m/s]
    real(DP) :: Tstar ! friction temperature [K]
    real(DP) :: Qstar ! friction mixing rate [kg/kg]
    real(DP) :: Uabs  ! modified absolute velocity [m/s]

    real(DP) :: T1  ! tempearature at the lowest atmospheric layer [K]
    real(DP) :: T0  ! skin temperature [K]
    real(DP) :: P1  ! pressure at the lowest atmospheric layer [Pa]
    real(DP) :: P0  ! surface pressure [Pa]
    real(DP) :: Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(DP) :: Q0  ! surface mixing ratio [kg/kg]
    real(DP) :: U1  ! zonal wind at the lowest atmospheric layer [m/s]
    real(DP) :: V1  ! meridional wind at the lowest atmospheric layer [m/s]
    real(DP) :: Z1  ! height at the lowest atmospheric layer [m]
    real(DP) :: PBL ! the top of atmospheric mixing layer [m]
    real(DP) :: Z0M ! roughness length of momentum [m]
    real(DP) :: Z0H ! roughness length of heat [m]
    real(DP) :: Z0E ! roughness length of moisture [m]

    real(DP) :: L      ! Obukhov length [m]
    real(DP) :: res    ! residual
    real(DP) :: dres   ! d(residual)/dL

    real(DP) :: RiB0 ! bulk Richardson number [no unit]
    real(DP) :: Wstar, dWstar ! free convection velocity scale [m/s]

    real(DP) :: UabsUS, UabsS, dUabsUS, dUabsS
    real(DP) :: UstarUS, UstarS, dUstar, dUstarUS, dUstarS
    real(DP) :: TstarUS, TstarS, dTstar, dTstarUS, dTstarS
    real(DP) :: QstarUS, QstarS, dQstar, dQstarUS, dQstarS

    real(DP) :: TH1, TH0
    real(DP) :: sw, tmp

    real(DP) :: log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E
    !---------------------------------------------------------------------------

    ! convert to DP
    T1  = real( RP_T1,  kind=DP )
    T0  = real( RP_T0,  kind=DP )
    P1  = real( RP_P1,  kind=DP )
    P0  = real( RP_P0,  kind=DP )
    Q1  = real( RP_Q1,  kind=DP )
    Q0  = real( RP_Q0,  kind=DP )
    U1  = real( RP_U1,  kind=DP )
    V1  = real( RP_V1,  kind=DP )
    Z1  = real( RP_Z1,  kind=DP )
    PBL = real( RP_PBL, kind=DP )
    Z0M = real( RP_Z0M, kind=DP )
    Z0H = real( RP_Z0H, kind=DP )
    Z0E = real( RP_Z0E, kind=DP )

    Uabs = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
    TH1  = T1 * ( PRE00 / P1 )**( Rdry / CPdry )
    TH0  = T0 * ( PRE00 / P0 )**( Rdry / CPdry )

    ! make log constant
    log_Z1ovZ0M = log( Z1 / Z0M )
    log_Z1ovZ0H = log( Z1 / Z0H )
    log_Z1ovZ0E = log( Z1 / Z0E )

    ! initial bulk Richardson number
    RiB0 = GRAV * Z1 * ( TH1 - TH0 ) / ( TH1 * Uabs**2 )
    if( abs( RiB0 ) < BULKFLUX_RiB_min ) then
      RiB0 = sign( BULKFLUX_RiB_min, RiB0 )
    end if

    ! initial Obukhov length assumed by neutral condition
    L = Z1 / RiB0 * log_Z1ovZ0H / log_Z1ovZ0M**2

    ! initial free convection velocity scale
    Wstar  = BULKFLUX_Wstar_min
    dWstar = BULKFLUX_Wstar_min

    do n = 1, BULKFLUX_itr_max
      ! unstable condition
      UabsUS  = max( sqrt( U1**2 + V1**2 + (BULKFLUX_WSCF*Wstar)**2 ), BULKFLUX_Uabs_min )
      UstarUS = KARMAN / ( log_Z1ovZ0M - fm_unstable(Z1,L) + fm_unstable(Z0M,L) ) * UabsUS
      TstarUS = KARMAN / ( log_Z1ovZ0H - fh_unstable(Z1,L) + fh_unstable(Z0H,L) ) / Pt * ( TH1 - TH0 )
      QstarUS = KARMAN / ( log_Z1ovZ0E - fh_unstable(Z1,L) + fh_unstable(Z0E,L) ) / Pt * ( Q1  - Q0  )

      ! stable condition
      UabsS  = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
      UstarS = KARMAN / ( log_Z1ovZ0M - fm_stable(Z1,L) + fm_stable(Z0M,L) ) * UabsS
      TstarS = KARMAN / ( log_Z1ovZ0H - fh_stable(Z1,L) + fh_stable(Z0H,L) ) / Pt * ( TH1 - TH0 )
      QstarS = KARMAN / ( log_Z1ovZ0E - fh_stable(Z1,L) + fh_stable(Z0E,L) ) / Pt * ( Q1  - Q0  )

      sw = 0.5D+0 - sign( 0.5D+0, L ) ! if unstable, sw = 1

      Uabs  = ( sw ) * UabsUS  + ( 1.0D+0-sw ) * UabsS
      Ustar = ( sw ) * UstarUS + ( 1.0D+0-sw ) * UstarS
      Tstar = ( sw ) * TstarUS + ( 1.0D+0-sw ) * TstarS
      Qstar = ( sw ) * QstarUS + ( 1.0D+0-sw ) * QstarS

      ! avoid zero division with Tstar = 0
      sw    = 0.5D+0 + sign( 0.5D+0, abs(Tstar) - EPS )
      Tstar = ( sw ) * Tstar + ( 1.0D+0-sw ) * EPS

      ! update free convection velocity scale
      tmp   = -PBL * GRAV / T1 * Ustar * Tstar
      sw    = 0.5D+0 + sign( 0.5D+0, tmp ) ! if tmp is plus, sw = 1
      Wstar = ( tmp * sw )**( 1.0D+0 / 3.0D+0 )

      ! calculate residual
      res = L - Ustar**2 * T1 / ( KARMAN * GRAV * Tstar )

      ! unstable condition
      dUabsUS  = max( sqrt( U1**2 + V1**2 + (BULKFLUX_WSCF*dWstar)**2 ), BULKFLUX_Uabs_min )
      dUstarUS = KARMAN / ( log_Z1ovZ0M - fm_unstable(Z1,L+dL) + fm_unstable(Z0M,L+dL) ) * dUabsUS
      dTstarUS = KARMAN / ( log_Z1ovZ0H - fh_unstable(Z1,L+dL) + fh_unstable(Z0H,L+dL) ) / Pt * ( TH1 - TH0 )
      dQstarUS = KARMAN / ( log_Z1ovZ0E - fh_unstable(Z1,L+dL) + fh_unstable(Z0E,L+dL) ) / Pt * ( Q1  - Q0  )
      ! stable condition
      dUabsS  = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
      dUstarS = KARMAN / ( log_Z1ovZ0M - fm_stable(Z1,L+dL) + fm_stable(Z0M,L+dL) ) * dUabsS
      dTstarS = KARMAN / ( log_Z1ovZ0H - fh_stable(Z1,L+dL) + fh_stable(Z0H,L+dL) ) / Pt * ( TH1 - TH0 )
      dQstarS = KARMAN / ( log_Z1ovZ0E - fh_stable(Z1,L+dL) + fh_stable(Z0E,L+dL) ) / Pt * ( Q1  - Q0  )

      sw = 0.5D+0 - sign( 0.5D+0, L+dL ) ! if unstable, sw = 1

      dUstar = ( sw ) * dUstarUS + ( 1.0D+0-sw ) * dUstarS
      dTstar = ( sw ) * dTstarUS + ( 1.0D+0-sw ) * dTstarS
      dQstar = ( sw ) * dQstarUS + ( 1.0D+0-sw ) * dQstarS

      ! avoid zero division with dTstar = 0
      sw     = 0.5D+0 + sign( 0.5D+0, abs(dTstar) - EPS )
      dTstar = ( sw ) * dTstar + ( 1.0D+0-sw ) * EPS

      ! update d(free convection velocity scale)
      tmp    = -PBL * GRAV / T1 * dUstar * dTstar
      sw     = 0.5D+0 + sign( 0.5D+0, tmp ) ! if tmp is plus, sw = 1
      dWstar = ( tmp * sw )**( 1.0D+0 / 3.0D+0 )

      ! calculate d(residual)/dL
      dres = 1.0D+0 - T1 / ( KARMAN * GRAV * dL ) * ( dUstar**2 / dTstar - Ustar**2 / Tstar )

      ! convergence test with residual and error levels
      if( abs( res      ) < BULKFLUX_res_min .or. &
          abs( res/dres ) < BULKFLUX_err_min      ) then
        exit
      end if

      ! stop iteration to prevent numerical error
      if( abs( dres ) < EPS ) then
        exit
      end if

      ! update Obukhov length
      L = L - res / dres

    end do

    ! revert to RP
    RP_Ustar = real( Ustar, kind=RP )
    RP_Tstar = real( Tstar, kind=RP )
    RP_Qstar = real( Qstar, kind=RP )
    RP_Uabs  = real( Uabs,  kind=RP )

    if( n > BULKFLUX_itr_max ) then
      if( IO_L ) write(IO_FID_LOG,'(A)'       ) 'Warning: reach maximum iteration in the function of BULKFLUX_B91W01.'
      if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
      if( IO_L ) write(IO_FID_LOG,'(A,I32)'   ) 'DEBUG --- LOOP number                     [no unit] :', n
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- Residual                        [m]       :', res
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- delta Residual                  [m]       :', dres
      if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- air tempearature                [K]       :', T1
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- surface temperature             [K]       :', T0
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- pressure                        [Pa]      :', P1
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- surface pressure                [Pa]      :', P0
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- water vapor mass ratio          [kg/kg]   :', Q1
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- surface water vapor mass ratio  [kg/kg]   :', Q0
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- zonal wind                      [m/s]     :', U1
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- meridional wind                 [m/s]     :', V1
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- cell center height              [m]       :', Z1
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- atmospheric mixing layer height [m]       :', PBL
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- roughness length of momentum    [m]       :', Z0M
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- roughness length of heat        [m]       :', Z0H
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- roughness length of moisture    [m]       :', Z0E
      if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- friction velocity               [m]       :', Ustar
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- friction potential temperature  [K]       :', Tstar
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- friction water vapor mass ratio [kg/kg]   :', Qstar
      if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- Obukhov length                  [m]       :', L
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- bulk Richardson number          [no unit] :', RiB0
      if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- free convection velocity scale  [m/s]     :', Wstar
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
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: L

    ! function
    real(DP) :: fm_unstable

    ! works
    real(DP) :: R
!    real(DP) :: r4R
    !---------------------------------------------------------------------------

    R = min( Z/L, 0.0D+0 )

    ! Wilson (2001)
    fm_unstable = 3.0D+0 * log( ( 1.0D+0 + sqrt( 1.0D+0 + 3.6D+0 * (-R)**(2.0D+0/3.0D+0) ) ) * 0.5D+0 )

    ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
    ! you should comment out the above line (Wilson 2001) and uncomment the below lines (Paulson 1974; Dyer 1974).
    !
    !! Paulson (1974); Dyer (1974)
    !r4R = ( 1.0D+0 - 16.0D+0 * R )**0.25D+0
    !fm_unstable = log( ( 1.0D+0 + r4R )**2 * ( 1.0D+0 + r4R * r4R ) * 0.125D+0 ) - 2.0D+0 * atan( r4R ) + PI * 0.5D+0

    return
  end function fm_unstable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in unstable condition
  function fh_unstable( Z, L )
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: L

    ! function
    real(DP) :: fh_unstable

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = min( Z/L, 0.0D+0 )

    ! Wilson (2001)
    fh_unstable = 3.0D+0 * log( ( 1.0D+0 + sqrt( 1.0D+0 + 7.9D+0 * (-R)**(2.0D+0/3.0D+0) ) ) * 0.5D+0 )

    ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
    ! you should comment out the above line (Wilson 2001) and uncomment the below lines (Paulson 1974; Dyer 1974).
    !
    !! Paulson (1974); Dyer (1974)
    !fh_unstable = 2.0D+0 * log( ( 1.0D+0 + sqrt( 1.0D+0 - 16.0D+0 * R ) ) * 0.5D+0 )

    return
  end function fh_unstable

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in stable condition
  function fm_stable( Z, L )
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: L

    ! function
    real(DP) :: fm_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(DP), parameter :: a = 1.0D+0
    real(DP), parameter :: b = 0.667D+0
    real(DP), parameter :: c = 5.0D+0
    real(DP), parameter :: d = 0.35D+0

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = max( Z/L, 0.0D+0 )

    ! Holtslag and DeBruin (1988)
    fm_stable = - a*R - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fm_stable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in stable condition
  function fh_stable( Z, L )
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: L

    ! function
    real(DP) :: fh_stable

    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(DP), parameter :: a = 1.0D+0
    real(DP), parameter :: b = 0.667D+0
    real(DP), parameter :: c = 5.0D+0
    real(DP), parameter :: d = 0.35D+0

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = max( Z/L, 0.0D+0 )

    ! Beljaars and Holtslag (1991)
    fh_stable = 1.0D+0 - ( 1.0D+0 + 2.0D+0/3.0D+0 * a*R )**1.5D+0 - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fh_stable

end module scale_bulkflux
