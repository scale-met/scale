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
          R10m, & ! (out)
          R02h, & ! (out)
          R02e, & ! (out)
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
       real(RP), intent(out) :: R10m ! lapse rate for 10m momentum [0-1]
       real(RP), intent(out) :: R02h ! lapse rate for 2m heat [0-1]
       real(RP), intent(out) :: R02e ! lapse rate for 2m moisture [0-1]

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
  private :: fmUS
  private :: fhUS
  private :: fmS
  private :: fhS

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private, save :: bulkcoef_TYPE = 'BH91'

  ! limiter
  real(RP), private, save :: Cm_min = 1.0E-8_RP ! minimum bulk coef. of u,v,w
  real(RP), private, save :: Ch_min = 1.0E-8_RP !                       T
  real(RP), private, save :: Ce_min = 1.0E-8_RP !                       q
  real(RP), private, save :: Cm_max = 1.0_RP    ! maximum bulk coef. of u,v,w
  real(RP), private, save :: Ch_max = 1.0_RP    !                       T
  real(RP), private, save :: Ce_max = 1.0_RP    !                       q

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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL_BULKCOEF)

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

  subroutine CPL_bulkcoef_uno( &
      Cm,   & ! (out)
      Ch,   & ! (out)
      Ce,   & ! (out)
      R10m, & ! (out)
      R02h, & ! (out)
      R02e, & ! (out)
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
    real(RP), intent(out) :: R10m ! lapse rate for 10m momentum [0-1]
    real(RP), intent(out) :: R02h ! lapse rate for 2m heat [0-1]
    real(RP), intent(out) :: R02e ! lapse rate for 2m moisture [0-1]

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
    real(RP) :: C0, C0_10, C0_02 ! initial drag coefficient [no unit]
    real(RP) :: fm, fh, t0th, q0qe
    real(RP) :: fm10, fm02, fh02

    C0    = ( KARMAN / log(      Za/Z0M ) )**2
    C0_10 = ( KARMAN / log( 10.0_RP/Z0M ) )**2
    C0_02 = ( KARMAN / log(  2.0_RP/Z0M ) )**2

    RiBT  = GRAV * Za * ( Ta - Ts*(Pa/Ps)**RovCP ) / ( Ta * Ua**2 )
    RiB   = RiBT

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
      fm   = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh   = fm
      fm10 = fm
      fm02 = fm
      fh02 = fm
    else
      ! unstable condition
      fm   = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0    * sqrt(      Za/Z0M ) * sqrt( abs( RiB ) ) )
      fh   = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0    * sqrt(      Za/Z0M ) * sqrt( abs( RiB ) ) )
      fm10 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0_10 * sqrt( 10.0_RP/Z0M ) * sqrt( abs( RiB ) ) )
      fm02 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0_02 * sqrt(  2.0_RP/Z0M ) * sqrt( abs( RiB ) ) )
      fh02 = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0_02 * sqrt(  2.0_RP/Z0M ) * sqrt( abs( RiB ) ) )
    end if

    t0th = 1.0_RP / ( 1.0_RP + log( Z0M/Z0H ) / log( Za/Z0M ) / sqrt( fm ) * fh )
    q0qe = 1.0_RP / ( 1.0_RP + log( Z0M/Z0E ) / log( Za/Z0M ) / sqrt( fm ) * fh )

    Cm = C0 * fm
    Ch = C0 * fh * t0th / tPrn
    Ce = C0 * fh * q0qe / tPrn

    R10m = log( 10.0_RP/Z0M ) / log( Za/Z0M ) / sqrt( fm10 / fm )
    R02h = ( log( Z0M/Z0H ) / log( Za/Z0M ) + sqrt( fm02 ) / fh02 * log( 2.0_RP/Z0M ) ) &
         / ( log( Z0M/Z0H ) / log( Za/Z0M ) + sqrt( fm   ) / fh   * log(     Za/Z0M ) )
    R02e = ( log( Z0M/Z0E ) / log( Za/Z0M ) + sqrt( fm02 ) / fh02 * log( 2.0_RP/Z0M ) ) &
         / ( log( Z0M/Z0E ) / log( Za/Z0M ) + sqrt( fm   ) / fh   * log(     Za/Z0M ) )

    Cm = min( max( Cm, Cm_min ), Cm_max )
    Ch = min( max( Ch, Ch_min ), Ch_max )
    Ce = min( max( Ce, Ce_min ), Ce_max )

    R10m = min( max( R10m, 0.0_RP ), 1.0_RP )
    R02h = min( max( R02h, 0.0_RP ), 1.0_RP )
    R02e = min( max( R02e, 0.0_RP ), 1.0_RP )

    return
  end subroutine CPL_bulkcoef_uno

  subroutine CPL_bulkcoef_beljaars( &
      Cm,   & ! (out)
      Ch,   & ! (out)
      Ce,   & ! (out)
      R10m, & ! (out)
      R02h, & ! (out)
      R02e, & ! (out)
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
    real(RP), intent(out) :: R10m ! lapse rate for 10m momentum [0-1]
    real(RP), intent(out) :: R02h ! lapse rate for 2m heat [0-1]
    real(RP), intent(out) :: R02e ! lapse rate for 2m moisture [0-1]

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
    real(RP), parameter :: dL      = 1.0E-10_RP ! delta Obukhov length [m]

    ! variables
    integer :: n
    real(RP) :: RiB0, RiB ! bulk Richardson number [no unit]
    real(RP) :: L ! Obukhov length [m]
    real(RP) :: res, dres
    real(RP) :: dCm, dCh, dCe
    real(RP) :: tmp
    real(RP) :: fm, fh
    real(RP) :: fm10, fm02, fh02

    RiB0 = GRAV * Za * ( Ta - Ts*(Pa/Ps)**RovCP ) / ( Ta * Ua**2 )
    if( abs( RiB0 ) < RiB_min ) then
      RiB0 = RiB_min
    end if

    ! The initial Obukhov length is assumed by bulk Richardson number.
    L = Za / RiB0

    do n = 1, nmax
      if( L < 0.0_RP ) then
        ! unstable condition
        Cm = KARMAN**2 / ( log(Za/Z0M) - fmUS(Za/L) + fmUS(Z0M/L) )**2
        Ch = KARMAN**2 / ( log(Za/Z0M) - fmUS(Za/L) + fmUS(Z0M/L) ) &
                       / ( log(Za/Z0H) - fhUS(Za/L) + fhUS(Z0H/L) )
        Ce = KARMAN**2 / ( log(Za/Z0M) - fmUS(Za/L) + fmUS(Z0M/L) ) &
                       / ( log(Za/Z0E) - fhUS(Za/L) + fhUS(Z0E/L) )
      else
        ! stable condition
        Cm = KARMAN**2 / ( log(Za/Z0M) - fmS(Za/L) + fmS(Z0M/L) )**2
        Ch = KARMAN**2 / ( log(Za/Z0M) - fmS(Za/L) + fmS(Z0M/L) ) &
                       / ( log(Za/Z0H) - fhS(Za/L) + fhS(Z0H/L) )
        Ce = KARMAN**2 / ( log(Za/Z0M) - fmS(Za/L) + fmS(Z0M/L) ) &
                       / ( log(Za/Z0E) - fhS(Za/L) + fhS(Z0E/L) )
      end if
      ! calculate residual
      res = L - Za * Cm**1.5_RP / ( KARMAN * Ch * RiB0 )

      if( L+dL < 0.0_RP ) then
        ! unstable condition
        dCm = KARMAN**2 / ( log(Za/Z0M) - fmUS(Za/(L+dL)) + fmUS(Z0M/(L+dL)) )**2
        dCh = KARMAN**2 / ( log(Za/Z0M) - fmUS(Za/(L+dL)) + fmUS(Z0M/(L+dL)) ) &
                        / ( log(Za/Z0H) - fhUS(Za/(L+dL)) + fhUS(Z0H/(L+dL)) )
        dCe = KARMAN**2 / ( log(Za/Z0M) - fmUS(Za/(L+dL)) + fmUS(Z0M/(L+dL)) ) &
                        / ( log(Za/Z0E) - fhUS(Za/(L+dL)) + fhUS(Z0E/(L+dL)) )
      else
        ! stable condition
        dCm = KARMAN**2 / ( log(Za/Z0M) - fmS(Za/(L+dL)) + fmS(Z0M/(L+dL)) )**2
        dCh = KARMAN**2 / ( log(Za/Z0M) - fmS(Za/(L+dL)) + fmS(Z0M/(L+dL)) ) &
                        / ( log(Za/Z0H) - fhS(Za/(L+dL)) + fhS(Z0H/(L+dL)) )
        dCe = KARMAN**2 / ( log(Za/Z0M) - fmS(Za/(L+dL)) + fmS(Z0M/(L+dL)) ) &
                        / ( log(Za/Z0E) - fhS(Za/(L+dL)) + fhS(Z0E/(L+dL)) )
      end if
      ! calculate d(residual)
      dres = (L+dL) - Za * dCm**1.5_RP / ( KARMAN * dCh * RiB0 )

      ! update Obukhov length
      L = L - res / ( dres - res ) * dL

      if( abs( res ) < res_min ) then
        ! finish iteration
        exit
      end if

    end do

    if( n > nmax ) then
      if( IO_L ) write(IO_FID_LOG,*) 'Error: reach maximum iteration in the function of CPL_bulkcoef_beljaars.'
    end if

    ! update bulk Richardson nubmer
    !RiB = Za * Cm**1.5_RP / ( L * KARMAN * Ch )

    if( L < 0.0_RP ) then
      ! unstable condition
      fm   = fmUS(      Za/L )
      fh   = fhUS(      Za/L )
      fm10 = fmUS( 10.0_RP/L )
      fm02 = fmUS(  2.0_RP/L )
      fh02 = fhUS(  2.0_RP/L )
    else
      ! stable condition
      fm   = fmS(      Za/L )
      fh   = fhS(      Za/L )
      fm10 = fmS( 10.0_RP/L )
      fm02 = fmS(  2.0_RP/L )
      fh02 = fhS(  2.0_RP/L )
    end if

    R10m = log( 10.0_RP/Z0M ) / log( Za/Z0M ) / sqrt( fm10 / fm )
    R02h = ( log( Z0M/Z0H ) / log( Za/Z0M ) + sqrt( fm02 ) / fh02 * log( 2.0_RP/Z0M ) ) &
         / ( log( Z0M/Z0H ) / log( Za/Z0M ) + sqrt( fm   ) / fh   * log(     Za/Z0M ) )
    R02e = ( log( Z0M/Z0E ) / log( Za/Z0M ) + sqrt( fm02 ) / fh02 * log( 2.0_RP/Z0M ) ) &
         / ( log( Z0M/Z0E ) / log( Za/Z0M ) + sqrt( fm   ) / fh   * log(     Za/Z0M ) )

    Cm = min( max( Cm, Cm_min ), Cm_max )
    Ch = min( max( Ch, Ch_min ), Ch_max )
    Ce = min( max( Ce, Ce_min ), Ce_max )

    R10m = min( max( R10m, 0.0_RP ), 1.0_RP )
    R02h = min( max( R02h, 0.0_RP ), 1.0_RP )
    R02e = min( max( R02e, 0.0_RP ), 1.0_RP )

    return
  end subroutine CPL_bulkcoef_beljaars

  ! stability function for momemtum in unstable condition
  function fmUS( R )
    use scale_const, only: &
      PI => CONST_PI
    implicit none

    ! argument
    real(RP), intent(in) :: R

    ! function
    real(RP) :: fmUS
    
    ! Paulson (1974); Dyer (1974)
    fmUS = log( ( 1.0_RP + ( 1.0_RP - 16.0_RP * R )**0.25_RP )**2 * ( 1.0_RP + sqrt( 1.0_RP - 16.0_RP * R ) ) * 0.125_RP ) &
         - 2.0_RP * atan( ( 1.0_RP - 16.0_RP * R )**0.25_RP ) + PI * 0.5_RP

    return
  end function fmUS

  ! stability function for heat/vapor in unstable condition
  function fhUS( R )
    implicit none

    ! argument
    real(RP), intent(in) :: R

    ! function
    real(RP) :: fhUS
    
    ! Paulson (1974); Dyer (1974)
    fhUS = 2.0_RP * log( ( 1.0_RP + sqrt( 1.0_RP - 16.0_RP * R ) ) * 0.5_RP )

    return
  end function fhUS

  ! stability function for momemtum in stable condition
  function fmS( R )
    implicit none

    ! argument
    real(RP), intent(in) :: R

    ! function
    real(RP) :: fmS
    
    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(RP), parameter :: a = 1.0_RP
    real(RP), parameter :: b = 0.667_RP
    real(RP), parameter :: c = 5.0_RP
    real(RP), parameter :: d = 0.35_RP

    ! Holtslag and DeBruin (1988)
    fmS = - a*R - b*( R - c/d )*exp( -d*R ) - b*c/d

    return
  end function fmS

  ! stability function for heat/vapor in stable condition
  function fhS( R )
    implicit none

    ! argument
    real(RP), intent(in) :: R

    ! function
    real(RP) :: fhS
    
    ! parameters of stability functions (Beljaars and Holtslag 1991)
    real(RP), parameter :: a = 1.0_RP
    real(RP), parameter :: b = 0.667_RP
    real(RP), parameter :: c = 5.0_RP
    real(RP), parameter :: d = 0.35_RP

    ! Beljaars and Holtslag (1991)
    fhS = 1.0_RP - ( 1.0_RP + 2.0_RP/3.0_RP * a*R )**1.5_RP - b*( R - c/d )*exp( -d*R ) - b*c/d 

    return
  end function fhS

end module scale_cpl_bulkcoef
