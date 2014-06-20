!-------------------------------------------------------------------------------
!> module COUPLER / Surface Bulk coefficient
!!
!! @par Description
!!          calculation of Bulk coefficient at the surface
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
module scale_cpl_bulkcoef
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
  public :: CPL_bulkcoef_setup

  abstract interface
     subroutine bc( &
          Cm,   & ! (out)
          Ch,   & ! (out)
          Ce,   & ! (out)
          Ta,   & ! (in)
          Ts,   & ! (in)
          Pa,   & ! (in)
          Ps,   & ! (in)
          Ua,   & ! (in)
          Za,   & ! (in)
          Z0M,  & ! (in)
          Z0H,  & ! (in)
          Z0E   ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(out) :: Cm   ! momentum bulk coefficient [no unit]
       real(RP), intent(out) :: Ch   ! heat bulk coefficient [no unit]
       real(RP), intent(out) :: Ce   ! moisture bulk coefficient [no unit]

       real(RP), intent(in) :: Ta  ! tempearature at the lowest atmospheric layer [K]
       real(RP), intent(in) :: Ts  ! skin temperature [K]
       real(RP), intent(in) :: Pa  ! pressure at the lowest atmospheric layer [Pa]
       real(RP), intent(in) :: Ps  ! surface pressure [Pa]
       real(RP), intent(in) :: Ua  ! wind speed at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: Za  ! height at the lowest atmospheric layer [m]
       real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
       real(RP), intent(in) :: Z0H ! roughness length of heat [m]
       real(RP), intent(in) :: Z0E ! roughness length of moisture [m]
     end subroutine bc
  end interface

  procedure(bc), pointer :: CPL_bulkcoef => NULL()
  public :: CPL_bulkcoef

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: CPL_bulkcoef_uno
  private :: CPL_bulkcoef_beljaars
  private :: fm_unstable
  private :: fh_unstable
  private :: fm_stable
  private :: fh_stable

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: bulkcoef_TYPE = 'BH91'

  ! limiter
  real(RP), private :: U_min  = 1.0E-2_RP ! minimum absolute velocity
  real(RP), private :: Cm_min = 1.0E-8_RP ! minimum bulk coef. of u,v,w
  real(RP), private :: Ch_min = 1.0E-8_RP !                       T
  real(RP), private :: Ce_min = 1.0E-8_RP !                       q
  real(RP), private :: Cm_max = 1.0_RP    ! maximum bulk coef. of u,v,w
  real(RP), private :: Ch_max = 1.0_RP    !                       T
  real(RP), private :: Ce_max = 1.0_RP    !                       q

contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine CPL_bulkcoef_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT) :: CPL_bulkcoef_TYPE

    real(RP) :: CPL_bulkcoef_Cm_min
    real(RP) :: CPL_bulkcoef_Ch_min
    real(RP) :: CPL_bulkcoef_Ce_min
    real(RP) :: CPL_bulkcoef_Cm_max
    real(RP) :: CPL_bulkcoef_Ch_max
    real(RP) :: CPL_bulkcoef_Ce_max

    NAMELIST / PARAM_CPL_BULKCOEF / &
       CPL_bulkcoef_TYPE,    &
       CPL_bulkcoef_Cm_min,  &
       CPL_bulkcoef_Ch_min,  &
       CPL_bulkcoef_Ce_min,  &
       CPL_bulkcoef_Cm_max,  &
       CPL_bulkcoef_Ch_max,  &
       CPL_bulkcoef_Ce_max

    integer :: ierr
    !---------------------------------------------------------------------------

    CPL_bulkcoef_TYPE   = bulkcoef_TYPE
    CPL_bulkcoef_Cm_min = Cm_min
    CPL_bulkcoef_Ch_min = Ch_min
    CPL_bulkcoef_Ce_min = Ce_min
    CPL_bulkcoef_Cm_max = Cm_max
    CPL_bulkcoef_Ch_max = Ch_max
    CPL_bulkcoef_Ce_max = Ce_max

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Bulk coefficient parameter'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_BULKCOEF,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_BULKCOEF. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CPL_BULKCOEF)

    bulkcoef_TYPE = CPL_bulkcoef_TYPE
    Cm_min        = CPL_bulkcoef_Cm_min
    Ch_min        = CPL_bulkcoef_Ch_min
    Ce_min        = CPL_bulkcoef_Ce_min
    Cm_max        = CPL_bulkcoef_Cm_max
    Ch_max        = CPL_bulkcoef_Ch_max
    Ce_max        = CPL_bulkcoef_Ce_max

    select case( bulkcoef_TYPE )
    case ( 'U95' )
       CPL_bulkcoef => CPL_bulkcoef_uno
    case ( 'BH91' )
       CPL_bulkcoef => CPL_bulkcoef_beljaars
    case default
       write(*,*) 'xxx invalid bulk scheme (', trim(bulkcoef_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine CPL_bulkcoef_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_bulkcoef_uno( &
      Cm,   & ! (out)
      Ch,   & ! (out)
      Ce,   & ! (out)
      Ta,   & ! (in)
      Ts,   & ! (in)
      Pa,   & ! (in)
      Ps,   & ! (in)
      Ua,   & ! (in)
      Za,   & ! (in)
      Z0M,  & ! (in)
      Z0H,  & ! (in)
      Z0E   ) ! (in)
    use scale_const, only: &
      GRAV   => CONST_GRAV,   &
      KARMAN => CONST_KARMAN, &
      RovCP  => CONST_RovCP
    implicit none

    ! argument
    real(RP), intent(out) :: Cm   ! momentum bulk coefficient [no unit]
    real(RP), intent(out) :: Ch   ! heat bulk coefficient [no unit]
    real(RP), intent(out) :: Ce   ! moisture bulk coefficient [no unit]

    real(RP), intent(in) :: Ta  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: Ts  ! skin temperature [K]
    real(RP), intent(in) :: Pa  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: Ps  ! surface pressure [Pa]
    real(RP), intent(in) :: Ua  ! wind speed at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Za  ! height at the lowest atmospheric layer [m]
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
    real(RP) :: RiB, RiBT ! bulk Richardson number [no unit]
    real(RP) :: C0 ! initial drag coefficient [no unit]
    real(RP) :: fm, fh, t0th, q0qe

    C0   = ( KARMAN / log(      Za/Z0M ) )**2
    RiBT = GRAV * Za * ( Ta*(Ps/Pa)**RovCP - Ts ) / ( Ta*(Ps/Pa)**RovCP * max(Ua,U_min)**2 )
    RiB  = RiBT

    if( RiBT >= 0.0_RP ) then
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

    if( RiBT >= 0.0_RP ) then
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

    Cm = C0 * fm
    Ch = C0 * fh * t0th / tPrn
    Ce = C0 * fh * q0qe / tPrn

    Cm = min( max( Cm, Cm_min ), Cm_max )
    Ch = min( max( Ch, Ch_min ), Ch_max )
    Ce = min( max( Ce, Ce_min ), Ce_max )

    return
  end subroutine CPL_bulkcoef_uno

  !-----------------------------------------------------------------------------
  subroutine CPL_bulkcoef_beljaars( &
      Cm,   & ! (out)
      Ch,   & ! (out)
      Ce,   & ! (out)
      Ta,   & ! (in)
      Ts,   & ! (in)
      Pa,   & ! (in)
      Ps,   & ! (in)
      Ua,   & ! (in)
      Za,   & ! (in)
      Z0M,  & ! (in)
      Z0H,  & ! (in)
      Z0E   ) ! (in)
    use scale_const, only: &
      GRAV   => CONST_GRAV,   &
      KARMAN => CONST_KARMAN, &
      RovCP  => CONST_RovCP
    implicit none

    ! argument
    real(RP), intent(out) :: Cm   ! momentum bulk coefficient [no unit]
    real(RP), intent(out) :: Ch   ! heat bulk coefficient [no unit]
    real(RP), intent(out) :: Ce   ! moisture bulk coefficient [no unit]

    real(RP), intent(in) :: Ta  ! tempearature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: Ts  ! skin temperature [K]
    real(RP), intent(in) :: Pa  ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: Ps  ! surface pressure [Pa]
    real(RP), intent(in) :: Ua  ! wind speed at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: Za  ! height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(in) :: Z0H ! roughness length of heat [m]
    real(RP), intent(in) :: Z0E ! roughness length of moisture [m]

    ! constant
    integer,  parameter :: nmax    = 100        ! maximum iteration number

    real(RP), parameter :: RiB_min = 1.0E-4_RP
    real(RP), parameter :: res_min = 1.0E-6_RP
    real(RP), parameter :: dL      = 1.0E-8_RP  ! delta Obukhov length [m]

    ! variables
    integer :: n
    real(RP) :: RiB0, RiB ! bulk Richardson number [no unit]
    real(RP) :: L ! Obukhov length [m]
    real(RP) :: res, dres
    real(RP) :: CmUS, ChUS, CeUS
    real(RP) :: CmS, ChS, CeS
    real(RP) :: dCm, dCh, dCe
    real(RP) :: dCmUS, dChUS, dCeUS
    real(RP) :: dCmS, dChS, dCeS
    real(RP) :: sw
    !---------------------------------------------------------------------------

    RiB0 = GRAV * Za * ( Ta*(Ps/Pa)**RovCP - Ts ) / ( Ta*(Ps/Pa)**RovCP * max(Ua,U_min)**2 )
    if( abs( RiB0 ) < RiB_min ) then
      RiB0 = RiB_min
    end if

    ! The initial Obukhov length is assumed by bulk Richardson number.
    L = Za / RiB0

    do n = 1, nmax
      ! unstable condition
      CmUS = KARMAN**2 / ( log(Za/Z0M) - fm_unstable(Za,L) + fm_unstable(Z0M,L) )**2
      ChUS = KARMAN**2 / ( log(Za/Z0M) - fm_unstable(Za,L) + fm_unstable(Z0M,L) ) &
                       / ( log(Za/Z0H) - fh_unstable(Za,L) + fh_unstable(Z0H,L) )
      CeUS = KARMAN**2 / ( log(Za/Z0M) - fm_unstable(Za,L) + fm_unstable(Z0M,L) ) &
                       / ( log(Za/Z0E) - fh_unstable(Za,L) + fh_unstable(Z0E,L) )
      ! stable condition
      CmS = KARMAN**2 / ( log(Za/Z0M) - fm_stable(Za,L) + fm_stable(Z0M,L) )**2
      ChS = KARMAN**2 / ( log(Za/Z0M) - fm_stable(Za,L) + fm_stable(Z0M,L) ) &
                      / ( log(Za/Z0H) - fh_stable(Za,L) + fh_stable(Z0H,L) )
      CeS = KARMAN**2 / ( log(Za/Z0M) - fm_stable(Za,L) + fm_stable(Z0M,L) ) &
                      / ( log(Za/Z0E) - fh_stable(Za,L) + fh_stable(Z0E,L) )

      sw = 0.5_RP - sign( 0.5_RP, L ) ! if unstable, sw = 1

      Cm = ( sw ) * CmUS + ( 1.0_RP-sw ) * CmS 
      Ch = ( sw ) * ChUS + ( 1.0_RP-sw ) * ChS 
      Ce = ( sw ) * CeUS + ( 1.0_RP-sw ) * CeS 

      ! calculate residual
      res = L - Za * Cm**1.5_RP / ( KARMAN * Ch * RiB0 )

      ! unstable condition
      dCmUS = KARMAN**2 / ( log(Za/Z0M) - fm_unstable(Za,L+dL) + fm_unstable(Z0M,L+dL) )**2
      dChUS = KARMAN**2 / ( log(Za/Z0M) - fm_unstable(Za,L+dL) + fm_unstable(Z0M,L+dL) ) &
                        / ( log(Za/Z0H) - fh_unstable(Za,L+dL) + fh_unstable(Z0H,L+dL) )
      dCeUS = KARMAN**2 / ( log(Za/Z0M) - fm_unstable(Za,L+dL) + fm_unstable(Z0M,L+dL) ) &
                        / ( log(Za/Z0E) - fh_unstable(Za,L+dL) + fh_unstable(Z0E,L+dL) )
      ! stable condition
      dCmS = KARMAN**2 / ( log(Za/Z0M) - fm_stable(Za,L+dL) + fm_stable(Z0M,L+dL) )**2
      dChS = KARMAN**2 / ( log(Za/Z0M) - fm_stable(Za,L+dL) + fm_stable(Z0M,L+dL) ) &
                       / ( log(Za/Z0H) - fh_stable(Za,L+dL) + fh_stable(Z0H,L+dL) )
      dCeS = KARMAN**2 / ( log(Za/Z0M) - fm_stable(Za,L+dL) + fm_stable(Z0M,L+dL) ) &
                       / ( log(Za/Z0E) - fh_stable(Za,L+dL) + fh_stable(Z0E,L+dL) )

      sw = 0.5_RP - sign( 0.5_RP, L+dL ) ! if unstable, sw = 1

      dCm = ( sw ) * dCmUS + ( 1.0_RP-sw ) * dCmS 
      dCh = ( sw ) * dChUS + ( 1.0_RP-sw ) * dChS 
      dCe = ( sw ) * dCeUS + ( 1.0_RP-sw ) * dCeS 

      ! calculate d(residual)/dL
      dres = ( (L+dL) - Za * dCm**1.5_RP / ( KARMAN * dCh * RiB0 ) - res ) / dL

      ! update Obukhov length
      L = L - res / dres

      if( abs( res ) < res_min ) then
        ! finish iteration
        exit
      end if

    end do

    if( n > nmax ) then
      if( IO_L ) write(IO_FID_LOG,*) 'Warning: reach maximum iteration in the function of CPL_bulkcoef_beljaars.'
    end if

    Cm = min( max( Cm, Cm_min ), Cm_max )
    Ch = min( max( Ch, Ch_min ), Ch_max )
    Ce = min( max( Ce, Ce_min ), Ce_max )

    ! update bulk Richardson nubmer
    RiB = Za * Cm**1.5_RP / ( L * KARMAN * Ch )

    return
  end subroutine CPL_bulkcoef_beljaars

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in unstable condition
  function fm_unstable( Z, L )
    use scale_const, only: &
      PI  => CONST_PI,  &
      EPS => CONST_EPS
    implicit none

    ! argument
    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L

    ! function
    real(RP) :: fm_unstable

    ! works
    real(RP) :: R
    real(RP) :: r4R
    !---------------------------------------------------------------------------

    R = Z / min(L,-EPS) ! should be negative
    r4R = ( 1.0_RP - 16.0_RP * R )**0.25_RP

    ! Paulson (1974); Dyer (1974)
    fm_unstable = log( ( 1.0_RP + r4R )**2 * ( 1.0_RP + r4R * r4R ) * 0.125_RP ) &
                - 2.0_RP * atan( r4R ) + PI * 0.5_RP

    return
  end function fm_unstable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in unstable condition
  function fh_unstable( Z, L )
    use scale_const, only: &
      EPS => CONST_EPS
    implicit none

    ! argument
    real(RP), intent(in) :: Z
    real(RP), intent(in) :: L

    ! function
    real(RP) :: fh_unstable

    ! works
    real(RP) :: R
    !---------------------------------------------------------------------------

    R = Z / min(L,-EPS) ! should be negative

    ! Paulson (1974); Dyer (1974)
    fh_unstable = 2.0_RP * log( ( 1.0_RP + sqrt( 1.0_RP - 16.0_RP * R ) ) * 0.5_RP )

    return
  end function fh_unstable

  !-----------------------------------------------------------------------------
  ! stability function for momemtum in stable condition
  function fm_stable( Z, L )
    use scale_const, only: &
      EPS => CONST_EPS
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

    R = Z / max(L,EPS) ! should be positive

    ! Holtslag and DeBruin (1988)
    fm_stable = - a*R - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fm_stable

  !-----------------------------------------------------------------------------
  ! stability function for heat/vapor in stable condition
  function fh_stable( Z, L )
    use scale_const, only: &
      EPS => CONST_EPS
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

    R = Z / max(L,EPS) ! should be positive

    ! Beljaars and Holtslag (1991)
    fh_stable = 1.0_RP - ( 1.0_RP + 2.0_RP/3.0_RP * a*R )**1.5_RP - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fh_stable

end module scale_cpl_bulkcoef
