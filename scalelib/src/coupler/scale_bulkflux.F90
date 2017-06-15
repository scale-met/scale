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
          FracU10, & ! (out)
          FracT2,  & ! (out)
          FracQ2,  & ! (out)
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

       real(RP), intent(out) :: Ustar   ! friction velocity [m/s]
       real(RP), intent(out) :: Tstar   ! friction temperature [K]
       real(RP), intent(out) :: Qstar   ! friction mixing rate [kg/kg]
       real(RP), intent(out) :: Uabs    ! modified absolute velocity [m/s]
       real(RP), intent(out) :: FracU10 ! calculation parameter for U10 [-]
       real(RP), intent(out) :: FracT2  ! calculation parameter for T2 [-]
       real(RP), intent(out) :: FracQ2  ! calculation parameter for Q2 [-]

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
  character(len=H_SHORT), private :: BULKFLUX_type = 'B91W01'

  integer,  private :: BULKFLUX_itr_sa_max = 5  ! maximum iteration number for successive approximation
  integer,  private :: BULKFLUX_itr_nr_max = 10 ! maximum iteration number for Newton-Raphson method

  real(RP), private :: BULKFLUX_err_min = 1.0E-3_RP ! minimum value of error

  real(RP), private :: BULKFLUX_WSCF = 1.2_RP ! empirical scaling factor of Wstar (Beljaars 1994)

  ! limiter
  real(RP), private :: BULKFLUX_Uabs_min  = 1.0E-2_RP ! minimum of Uabs [m/s]
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
       BULKFLUX_type,       &
       BULKFLUX_itr_sa_max, &
       BULKFLUX_itr_nr_max, &
       BULKFLUX_err_min,    &
       BULKFLUX_WSCF,       &
       BULKFLUX_Uabs_min,   &
       BULKFLUX_Wstar_min

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[BULKFLUX] / Categ[COUPLER] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BULKFLUX,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_BULKFLUX. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_BULKFLUX)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Scheme for surface bulk flux : ', trim(BULKFLUX_type)
    select case(BULKFLUX_type)
    case('U95')
       if( IO_L ) write(IO_FID_LOG,*) '*** => Uno et al.(1995)'
       BULKFLUX => BULKFLUX_U95
    case('B91W01')
       if( IO_L ) write(IO_FID_LOG,*) '*** => Beljaars (1991) and Wilson (2001)'
       BULKFLUX => BULKFLUX_B91W01
    case default
       write(*,*) 'xxx Unsupported BULKFLUX_type. STOP'
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
      FracU10, & ! (out)
      FracT2,  & ! (out)
      FracQ2,  & ! (out)
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
    real(RP), intent(out) :: Ustar   ! friction velocity [m/s]
    real(RP), intent(out) :: Tstar   ! friction temperature [K]
    real(RP), intent(out) :: Qstar   ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: Uabs    ! modified absolute velocity [m/s]
    real(RP), intent(out) :: FracU10 ! calculation parameter for U10 [-]
    real(RP), intent(out) :: FracT2  ! calculation parameter for T2 [-]
    real(RP), intent(out) :: FracQ2  ! calculation parameter for Q2 [-]

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
    real(RP) :: RiB0, RiB ! bulk Richardson number [-]
    real(RP) :: C0Z1, C010, C002 ! initial drag coefficient [-]
    real(RP) :: CmZ1, ChZ1, CqZ1, fmZ1, fhZ1, t0thZ1, q0qeZ1
    real(RP) :: Cm10, Ch10, Cq10, fm10
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

    Uabs = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
    TH1  = T1 * ( PRE00 / P1 )**( Rdry / CPdry )
    TH0  = T0 * ( PRE00 / P0 )**( Rdry / CPdry )

    ! bulk Richardson number
    RiB0 = GRAV * Z1 * ( TH1 - TH0 ) / ( TH1 * Uabs**2 )

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

    Ustar = sqrt( CmZ1 ) * Uabs
    Tstar = ChZ1 * Uabs / Ustar * ( TH1 - TH0 )
    Qstar = CqZ1 * Uabs / Ustar * ( Q1  - Q0  )

    FracU10 = sqrt( CmZ1 / Cm10 )
    FracT2  = ChZ1 / Ch02 * sqrt( Cm02 / CmZ1 )
    FracQ2  = CqZ1 / Cq02 * sqrt( Cm02 / CmZ1 )

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
  subroutine BULKFLUX_B91W01( &
      Ustar,   & ! (out)
      Tstar,   & ! (out)
      Qstar,   & ! (out)
      Uabs,    & ! (out)
      FracU10, & ! (out)
      FracT2,  & ! (out)
      FracQ2,  & ! (out)
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
      GRAV    => CONST_GRAV,    &
      KARMAN  => CONST_KARMAN,  &
      Rdry    => CONST_Rdry,    &
      CPdry   => CONST_CPdry,   &
      EPSTvap => CONST_EPSTvap, &
      EPS     => CONST_EPS,     &
      PRE00   => CONST_PRE00
    implicit none

    ! parameter
    real(DP), parameter :: dIL = 1.0E-6_DP ! delta [1/m]

    real(DP), parameter :: Pt = 0.95_DP ! turbulent Prandtl number

    ! argument
    real(RP), intent(out) :: Ustar   ! friction velocity [m/s]
    real(RP), intent(out) :: Tstar   ! friction temperature [K]
    real(RP), intent(out) :: Qstar   ! friction mixing rate [kg/kg]
    real(RP), intent(out) :: Uabs    ! modified absolute velocity [m/s]
    real(RP), intent(out) :: FracU10 ! calculation parameter for U10 [-]
    real(RP), intent(out) :: FracT2  ! calculation parameter for T2 [-]
    real(RP), intent(out) :: FracQ2  ! calculation parameter for Q2 [-]

    real(RP), intent(in) :: T1  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: T0  ! skin temperature [K]
    real(RP), intent(in) :: P1  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: P0  ! surface pressure [Pa]
    real(RP), intent(in) :: Q1  ! mixing ratio at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Q0  ! mixing ratio at surface [kg/kg]
    real(RP), intent(in) :: U1  ! zonal wind at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: V1  ! meridional wind at the lowest atmospheric layer [m/s]
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
    real(DP) :: Wstar, dWstar ! free convection velocity scale [m/s]

    real(DP) :: UabsUS, UabsS, UabsC
    real(DP) :: dUabsUS, dUabsS

    real(DP) :: UstarUS, UstarS, UstarC
    real(DP) :: TstarUS, TstarS, TstarC
    real(DP) :: QstarUS, QstarS, QstarC

    real(DP) :: dUstarUS, dUstarS, dUstarC
    real(DP) :: dTstarUS, dTstarS, dTstarC
    real(DP) :: dQstarUS, dQstarS, dQstarC

    real(DP) :: FracU10US, FracU10S, FracU10C
    real(DP) :: FracT2US,  FracT2S,  FracT2C
    real(DP) :: FracQ2US,  FracQ2S,  FracQ2C

    real(DP) :: TH1, TH0, THM
    real(DP) :: TV1, TV0, TVM
    real(DP) :: QM
    real(DP) :: sw, tmp

    real(DP) :: BFLX, dBFLX

    real(DP) :: DP_Z1, DP_Z0M, DP_Z0H, DP_Z0E
    real(DP) :: log_Z1ovZ0M, log_Z1ovZ0H, log_Z1ovZ0E
    real(DP) :: log_10ovZ0M, log_02ovZ0H, log_02ovZ0E
    !---------------------------------------------------------------------------

    ! convert to DP
    DP_Z1  = real( Z1,  kind=DP )
    DP_Z0M = real( Z0M, kind=DP )
    DP_Z0H = real( Z0H, kind=DP )
    DP_Z0E = real( Z0E, kind=DP )

    UabsC = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )

    TH1 = T1 * ( PRE00 / P1 )**( Rdry / CPdry )
    TH0 = T0 * ( PRE00 / P0 )**( Rdry / CPdry )
    THM = ( TH1 + TH0 ) * 0.5_RP
    QM  = ( Q1  + Q0  ) * 0.5_RP
    TV1 = TH1 * ( 1.0_DP + EPSTvap * Q1 )
    TV0 = TH0 * ( 1.0_DP + EPSTvap * Q0 )
    TVM = ( TV1 + TV0 ) * 0.5_RP

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
    Wstar  = BULKFLUX_Wstar_min
    dWstar = BULKFLUX_Wstar_min

    ! Successive approximation
    do n = 1, BULKFLUX_itr_sa_max
      ! unstable condition
      UabsUS  = max( sqrt( U1**2 + V1**2 + (BULKFLUX_WSCF*Wstar)**2 ), real( BULKFLUX_Uabs_min, kind=DP ) )
      UstarUS = KARMAN / ( log_Z1ovZ0M - fm_unstable(DP_Z1,IL) + fm_unstable(DP_Z0M,IL) ) * UabsUS
      TstarUS = KARMAN / ( log_Z1ovZ0H - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
      QstarUS = KARMAN / ( log_Z1ovZ0E - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

      ! stable condition
      UabsS  = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
      UstarS = KARMAN / ( log_Z1ovZ0M - fm_stable(DP_Z1,IL) + fm_stable(DP_Z0M,IL) ) * UabsS
      TstarS = KARMAN / ( log_Z1ovZ0H - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
      QstarS = KARMAN / ( log_Z1ovZ0E - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

      sw = 0.5_DP - sign( 0.5_DP, IL ) ! if unstable, sw = 1

      UstarC = ( sw ) * UstarUS + ( 1.0_DP-sw ) * UstarS
      TstarC = ( sw ) * TstarUS + ( 1.0_DP-sw ) * TstarS
      QstarC = ( sw ) * QstarUS + ( 1.0_DP-sw ) * QstarS

      ! estimate buoyancy flux
      BFLX = - UstarC * TstarC * ( 1.0_RP + EPSTvap * QM ) - EPSTvap * UstarC * QstarC * THM

      ! update free convection velocity scale
      tmp = PBL * GRAV / T1 * BFLX
      sw  = 0.5_DP + sign( 0.5_DP, tmp ) ! if tmp is plus, sw = 1

      Wstar = ( tmp * sw )**( 1.0_DP / 3.0_DP )

      ! avoid zero division with UstarC = 0
      sw = 0.5_DP + sign( 0.5_DP, abs( UstarC ) - EPS )

      UstarC = ( sw ) * UstarC + ( 1.0_DP-sw ) * EPS

      ! estimate the inversed Obukhov length
      IL = - KARMAN * GRAV * BFLX / ( UstarC**3 * THM )
    end do

    ! Newton-Raphson method
    do n = 1, BULKFLUX_itr_nr_max
      ! unstable condition
      UabsUS  = max( sqrt( U1**2 + V1**2 + (BULKFLUX_WSCF*Wstar)**2 ), real( BULKFLUX_Uabs_min, kind=DP ) )
      UstarUS = KARMAN / ( log_Z1ovZ0M - fm_unstable(DP_Z1,IL) + fm_unstable(DP_Z0M,IL) ) * UabsUS
      TstarUS = KARMAN / ( log_Z1ovZ0H - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
      QstarUS = KARMAN / ( log_Z1ovZ0E - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

      ! stable condition
      UabsS  = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
      UstarS = KARMAN / ( log_Z1ovZ0M - fm_stable(DP_Z1,IL) + fm_stable(DP_Z0M,IL) ) * UabsS
      TstarS = KARMAN / ( log_Z1ovZ0H - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
      QstarS = KARMAN / ( log_Z1ovZ0E - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

      sw = 0.5_DP - sign( 0.5_DP, IL ) ! if unstable, sw = 1

      UstarC = ( sw ) * UstarUS + ( 1.0_DP-sw ) * UstarS
      TstarC = ( sw ) * TstarUS + ( 1.0_DP-sw ) * TstarS
      QstarC = ( sw ) * QstarUS + ( 1.0_DP-sw ) * QstarS

      ! estimate buoyancy flux
      BFLX = - UstarC * TstarC * ( 1.0_RP + EPSTvap * QM ) - EPSTvap * UstarC * QstarC * THM

      ! update free convection velocity scale
      tmp = PBL * GRAV / T1 * BFLX
      sw  = 0.5_DP + sign( 0.5_DP, tmp ) ! if tmp is plus, sw = 1

      Wstar = ( tmp * sw )**( 1.0_DP / 3.0_DP )

      ! avoid zero division with UstarC = 0
      sw = 0.5_DP + sign( 0.5_DP, abs( UstarC ) - EPS )

      UstarC = ( sw ) * UstarC + ( 1.0_DP-sw ) * EPS

      ! calculate residual
      res = IL + KARMAN * GRAV * BFLX / ( UstarC**3 * THM )

      ! unstable condition
      dUabsUS  = max( sqrt( U1**2 + V1**2 + (BULKFLUX_WSCF*dWstar)**2 ), real( BULKFLUX_Uabs_min, kind=DP ) )
      dUstarUS = KARMAN / ( log_Z1ovZ0M - fm_unstable(DP_Z1,IL+dIL) + fm_unstable(DP_Z0M,IL+dIL) ) * dUabsUS
      dTstarUS = KARMAN / ( log_Z1ovZ0H - fh_unstable(DP_Z1,IL+dIL) + fh_unstable(DP_Z0H,IL+dIL) ) / Pt * ( TH1 - TH0 )
      dQstarUS = KARMAN / ( log_Z1ovZ0E - fh_unstable(DP_Z1,IL+dIL) + fh_unstable(DP_Z0E,IL+dIL) ) / Pt * ( Q1  - Q0  )

      ! stable condition
      dUabsS  = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
      dUstarS = KARMAN / ( log_Z1ovZ0M - fm_stable(DP_Z1,IL+dIL) + fm_stable(DP_Z0M,IL+dIL) ) * dUabsS
      dTstarS = KARMAN / ( log_Z1ovZ0H - fh_stable(DP_Z1,IL+dIL) + fh_stable(DP_Z0H,IL+dIL) ) / Pt * ( TH1 - TH0 )
      dQstarS = KARMAN / ( log_Z1ovZ0E - fh_stable(DP_Z1,IL+dIL) + fh_stable(DP_Z0E,IL+dIL) ) / Pt * ( Q1  - Q0  )

      sw = 0.5_DP - sign( 0.5_DP, IL+dIL ) ! if unstable, sw = 1

      dUstarC = ( sw ) * dUstarUS + ( 1.0_DP-sw ) * dUstarS
      dTstarC = ( sw ) * dTstarUS + ( 1.0_DP-sw ) * dTstarS
      dQstarC = ( sw ) * dQstarUS + ( 1.0_DP-sw ) * dQstarS

      ! estimate buoyancy flux
      dBFLX = - dUstarC * dTstarC * ( 1.0_RP + EPSTvap * QM ) - EPSTvap * dUstarC * dQstarC * THM

      ! update d(free convection velocity scale)
      tmp = -PBL * GRAV / T1 * dBFLX
      sw  = 0.5_DP + sign( 0.5_DP, tmp ) ! if tmp is plus, sw = 1

      dWstar = ( tmp * sw )**( 1.0_DP / 3.0_DP )

      ! avoid zero division with dUstarC = 0
      sw = 0.5_DP + sign( 0.5_DP, abs( dUstarC ) - EPS )

      dUstarC = ( sw ) * dUstarC + ( 1.0_DP-sw ) * EPS

      ! calculate d(residual)/dIL
      dres = 1.0_DP + KARMAN * GRAV / ( THM * dIL ) * ( dBFLX / dUstarC**3 - BFLX / UstarC**3 )

      ! stop iteration to prevent numerical error
      if( abs( dres ) < EPS ) exit

      ! convergence test with error levels
      if( abs( res/dres ) < BULKFLUX_err_min ) exit

      ! avoid sign changing
      if( IL * ( IL - res / dres ) < 0.0_RP ) exit

      ! update the inversed Obukhov length
      IL = IL - res / dres
    end do

    ! Successive approximation after Newton-Raphson method
    if( .NOT. abs( res/dres ) < BULKFLUX_err_min ) then
      ! unstable condition
      UabsUS  = max( sqrt( U1**2 + V1**2 + (BULKFLUX_WSCF*Wstar)**2 ), real( BULKFLUX_Uabs_min, kind=DP ) )
      UstarUS = KARMAN / ( log_Z1ovZ0M - fm_unstable(DP_Z1,IL) + fm_unstable(DP_Z0M,IL) ) * UabsUS
      TstarUS = KARMAN / ( log_Z1ovZ0H - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
      QstarUS = KARMAN / ( log_Z1ovZ0E - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

      ! stable condition
      UabsS  = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
      UstarS = KARMAN / ( log_Z1ovZ0M - fm_stable(DP_Z1,IL) + fm_stable(DP_Z0M,IL) ) * UabsS
      TstarS = KARMAN / ( log_Z1ovZ0H - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
      QstarS = KARMAN / ( log_Z1ovZ0E - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

      sw = 0.5_DP - sign( 0.5_DP, IL ) ! if unstable, sw = 1

      UstarC = ( sw ) * UstarUS + ( 1.0_DP-sw ) * UstarS
      TstarC = ( sw ) * TstarUS + ( 1.0_DP-sw ) * TstarS
      QstarC = ( sw ) * QstarUS + ( 1.0_DP-sw ) * QstarS

      ! estimate buoyancy flux
      BFLX = - UstarC * TstarC * ( 1.0_RP + EPSTvap * QM ) - EPSTvap * UstarC * QstarC * THM

      ! update free convection velocity scale
      tmp = PBL * GRAV / T1 * BFLX
      sw  = 0.5_DP + sign( 0.5_DP, tmp ) ! if tmp is plus, sw = 1

      Wstar = ( tmp * sw )**( 1.0_DP / 3.0_DP )

      ! avoid zero division with UstarC = 0
      sw = 0.5_DP + sign( 0.5_DP, abs( UstarC ) - EPS )

      UstarC = ( sw ) * UstarC + ( 1.0_DP-sw ) * EPS

      ! estimate the inversed Obukhov length
      IL = - KARMAN * GRAV * BFLX / ( UstarC**3 * THM )
    end if

    ! calculate Ustar, Tstar, and Qstar based on IL

    ! unstable condition
    UabsUS  = max( sqrt( U1**2 + V1**2 + (BULKFLUX_WSCF*Wstar)**2 ), real( BULKFLUX_Uabs_min, kind=DP ) )
    UstarUS = KARMAN / ( log_Z1ovZ0M - fm_unstable(DP_Z1,IL) + fm_unstable(DP_Z0M,IL) ) * UabsUS
    TstarUS = KARMAN / ( log_Z1ovZ0H - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
    QstarUS = KARMAN / ( log_Z1ovZ0E - fh_unstable(DP_Z1,IL) + fh_unstable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

    FracU10US = ( log_10ovZ0M - fm_unstable(10.0_DP,IL) + fm_unstable(DP_Z0M,IL) ) &
              / ( log_Z1ovZ0M - fm_unstable(  DP_Z1,IL) + fm_unstable(DP_Z0M,IL) )
    FracT2US  = ( log_02ovZ0H - fm_unstable( 2.0_DP,IL) + fm_unstable(DP_Z0H,IL) ) &
              / ( log_Z1ovZ0H - fm_unstable(  DP_Z1,IL) + fm_unstable(DP_Z0H,IL) )
    FracQ2US  = ( log_02ovZ0E - fm_unstable( 2.0_DP,IL) + fm_unstable(DP_Z0E,IL) ) &
              / ( log_Z1ovZ0E - fm_unstable(  DP_Z1,IL) + fm_unstable(DP_Z0E,IL) )

    ! stable condition
    UabsS  = max( sqrt( U1**2 + V1**2 ), BULKFLUX_Uabs_min )
    UstarS = KARMAN / ( log_Z1ovZ0M - fm_stable(DP_Z1,IL) + fm_stable(DP_Z0M,IL) ) * UabsS
    TstarS = KARMAN / ( log_Z1ovZ0H - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0H,IL) ) / Pt * ( TH1 - TH0 )
    QstarS = KARMAN / ( log_Z1ovZ0E - fh_stable(DP_Z1,IL) + fh_stable(DP_Z0E,IL) ) / Pt * ( Q1  - Q0  )

    FracU10S = ( log_10ovZ0M - fm_stable(10.0_DP,IL) + fm_stable(DP_Z0M,IL) ) &
             / ( log_Z1ovZ0M - fm_stable(  DP_Z1,IL) + fm_stable(DP_Z0M,IL) )
    FracT2S  = ( log_02ovZ0H - fm_stable( 2.0_DP,IL) + fm_stable(DP_Z0H,IL) ) &
             / ( log_Z1ovZ0H - fm_stable(  DP_Z1,IL) + fm_stable(DP_Z0H,IL) )
    FracQ2S  = ( log_02ovZ0E - fm_stable( 2.0_DP,IL) + fm_stable(DP_Z0E,IL) ) &
             / ( log_Z1ovZ0E - fm_stable(  DP_Z1,IL) + fm_stable(DP_Z0E,IL) )

    sw = 0.5_DP - sign( 0.5_DP, IL ) ! if unstable, sw = 1

    UstarC   = ( sw ) * UstarUS   + ( 1.0_DP-sw ) * UstarS
    TstarC   = ( sw ) * TstarUS   + ( 1.0_DP-sw ) * TstarS
    QstarC   = ( sw ) * QstarUS   + ( 1.0_DP-sw ) * QstarS
    UabsC    = ( sw ) * UabsUS    + ( 1.0_DP-sw ) * UabsS
    FracU10C = ( sw ) * FracU10US + ( 1.0_DP-sw ) * FracU10S
    FracT2C  = ( sw ) * FracT2US  + ( 1.0_DP-sw ) * FracT2S
    FracQ2C  = ( sw ) * FracQ2US  + ( 1.0_DP-sw ) * FracQ2S

    ! revert to RP
    Ustar   = real( UstarC,   kind=RP )
    Tstar   = real( TstarC,   kind=RP )
    Qstar   = real( QstarC,   kind=RP )
    Uabs    = real( UabsC,    kind=RP )
    FracU10 = real( FracU10C, kind=RP )
    FracT2  = real( FracT2C,  kind=RP )
    FracQ2  = real( FracQ2C,  kind=RP )

    return
  end subroutine BULKFLUX_B91W01

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

    ! works
    real(DP) :: R
    real(DP) :: r4R
    !---------------------------------------------------------------------------

    R = min( Z * IL, 0.0_DP )

    ! Wilson (2001)
    fm_unstable = 3.0_DP * log( ( 1.0_DP + sqrt( 1.0_DP + 3.6_DP * (-R)**(2.0_DP/3.0_DP) ) ) * 0.5_DP )

    ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
    ! you should comment out the above line (Wilson 2001) and uncomment the below lines (Paulson 1974; Dyer 1974).
    !
    !! Paulson (1974); Dyer (1974)
    !r4R = ( 1.0_DP - 16.0_DP * R )**0.25_DP
    !fm_unstable = log( ( 1.0_DP + r4R )**2 * ( 1.0_DP + r4R * r4R ) * 0.125_DP ) - 2.0_DP * atan( r4R ) + PI * 0.5_DP

    return
  end function fm_unstable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in unstable condition
  function fh_unstable( Z, IL )
    implicit none

    ! argument
    real(DP), intent(in) :: Z
    real(DP), intent(in) :: IL

    ! function
    real(DP) :: fh_unstable

    ! works
    real(DP) :: R
    !---------------------------------------------------------------------------

    R = min( Z * IL, 0.0_DP )

    ! Wilson (2001)
    fh_unstable = 3.0_DP * log( ( 1.0_DP + sqrt( 1.0_DP + 7.9_DP * (-R)**(2.0_DP/3.0_DP) ) ) * 0.5_DP )

    ! If you want to run with the original Beljaars scheme (Beljaars and Holtslag 1994),
    ! you should comment out the above line (Wilson 2001) and uncomment the below lines (Paulson 1974; Dyer 1974).
    !
    !! Paulson (1974); Dyer (1974)
    !fh_unstable = 2.0_DP * log( ( 1.0_DP + sqrt( 1.0_DP - 16.0_DP * R ) ) * 0.5_DP )

    return
  end function fh_unstable

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
#if defined(__PGI) || defined(__ES2)
    fm_stable = - a*R - b*( R - c/d )*exp( -min( d*R, 1.E+3_RP ) ) - b*c/d ! apply exp limiter
#else
    fm_stable = - a*R - b*( R - c/d )*exp( -d*R ) - b*c/d
#endif

    return
  end function fm_stable

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
#if defined(__PGI) || defined(__ES2)
    fh_stable = 1.0_DP - ( 1.0_DP + 2.0_DP/3.0_DP * a*R )**1.5_DP - b*( R - c/d )*exp( -min( d*R, 1.E+3_RP ) ) - b*c/d ! apply exp limiter
#else
    fh_stable = 1.0_DP - ( 1.0_DP + 2.0_DP/3.0_DP * a*R )**1.5_DP - b*( R - c/d )*exp( -d*R ) - b*c/d
#endif

    return
  end function fh_stable

end module scale_bulkflux
