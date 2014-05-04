!-------------------------------------------------------------------------------
!> module OCEAN / Surface roughness length
!!
!! @par Description
!!          roughness length of sea surface
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_ocean_roughness
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
  public :: OCEAN_roughness_setup

  abstract interface
     subroutine rl( &
          Z0W,           & ! (inout)
          Z0M, Z0H, Z0E, & ! (out)
          Uabs, ZA       ) ! (in)
       use scale_precision
       implicit none

       real(RP), intent(inout) :: Z0W ! roughness length of sea surface [m]

       real(RP), intent(out) :: Z0M ! roughness length of momentum [m]
       real(RP), intent(out) :: Z0H ! roughness length of heat [m]
       real(RP), intent(out) :: Z0E ! roughness length of vapor [m]

       real(RP), intent(in) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: ZA   ! height at 1st atm. layer [m]
     end subroutine rl
  end interface

  procedure(rl), pointer :: OCEAN_roughness => NULL()
  public :: OCEAN_roughness

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: OCEAN_roughness_miller92
  private :: OCEAN_roughness_moon07

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private, save :: ROUGHNESS_TYPE = 'MILLER' ! sea roughness length scheme

  real(RP), private, save :: CM0   = 1.0E-3_RP ! bulk coef. for U*
  real(RP), private, save :: visck = 1.5E-5_RP ! kinematic viscosity

  real(RP), private, save :: Ustar_min = 1.0E-3_RP ! minimum fiction velocity
  real(RP), private, save :: Z0M_min   = 1.0E-5_RP ! minimum roughness length for u,v,w
  real(RP), private, save :: Z0H_min   = 1.0E-5_RP !                              T
  real(RP), private, save :: Z0E_min   = 1.0E-5_RP !                              q

  real(RP), private, save :: Z0MI = 0.0E-0_RP ! base roughness rength for u,v,w
  real(RP), private, save :: Z0MR = 1.8E-2_RP ! rough factor for u,v,w
  real(RP), private, save :: Z0MS = 1.1E-1_RP ! smooth factor for u,v,w
  real(RP), private, save :: Z0HI = 1.4E-5_RP ! base roughness rength for T
  real(RP), private, save :: Z0HR = 0.0E-0_RP ! rough factor for T
  real(RP), private, save :: Z0HS = 4.0E-1_RP ! smooth factor for T
  real(RP), private, save :: Z0EI = 1.3E-4_RP ! base roughness rength for q
  real(RP), private, save :: Z0ER = 0.0E-0_RP ! rough factor for q
  real(RP), private, save :: Z0ES = 6.2E-1_RP ! smooth factor for q

contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine OCEAN_roughness_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT) :: OCEAN_roughness_TYPE

    real(RP) :: OCEAN_roughness_CM0
    real(RP) :: OCEAN_roughness_visck
    real(RP) :: OCEAN_roughness_Z0M_min
    real(RP) :: OCEAN_roughness_Z0H_min
    real(RP) :: OCEAN_roughness_Z0E_min
    real(RP) :: OCEAN_roughness_Z0MI
    real(RP) :: OCEAN_roughness_Z0MR
    real(RP) :: OCEAN_roughness_Z0MS
    real(RP) :: OCEAN_roughness_Z0HI
    real(RP) :: OCEAN_roughness_Z0HR
    real(RP) :: OCEAN_roughness_Z0HS
    real(RP) :: OCEAN_roughness_Z0EI
    real(RP) :: OCEAN_roughness_Z0ER
    real(RP) :: OCEAN_roughness_Z0ES

    NAMELIST / PARAM_OCEAN_ROUGHNESS / &
       OCEAN_roughness_TYPE,    &
       OCEAN_roughness_CM0,     &
       OCEAN_roughness_visck,   &
       OCEAN_roughness_Z0M_min, &
       OCEAN_roughness_Z0H_min, &
       OCEAN_roughness_Z0E_min, &
       OCEAN_roughness_Z0MI,    &
       OCEAN_roughness_Z0MR,    &
       OCEAN_roughness_Z0MS,    &
       OCEAN_roughness_Z0HI,    &
       OCEAN_roughness_Z0HR,    &
       OCEAN_roughness_Z0HS,    &
       OCEAN_roughness_Z0EI,    &
       OCEAN_roughness_Z0ER,    &
       OCEAN_roughness_Z0ES

    integer :: ierr
    !---------------------------------------------------------------------------


    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean roughness length parameter'

    OCEAN_roughness_TYPE    = ROUGHNESS_TYPE
    OCEAN_roughness_CM0     = CM0
    OCEAN_roughness_visck   = visck
    OCEAN_roughness_Z0M_min = Z0M_min
    OCEAN_roughness_Z0H_min = Z0H_min
    OCEAN_roughness_Z0E_min = Z0E_min
    OCEAN_roughness_Z0MI    = Z0MI
    OCEAN_roughness_Z0MR    = Z0MR
    OCEAN_roughness_Z0MS    = Z0MS
    OCEAN_roughness_Z0HI    = Z0HI
    OCEAN_roughness_Z0HR    = Z0HR
    OCEAN_roughness_Z0HS    = Z0HS
    OCEAN_roughness_Z0EI    = Z0EI
    OCEAN_roughness_Z0ER    = Z0ER
    OCEAN_roughness_Z0ES    = Z0ES

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_ROUGHNESS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_ROUGHNESS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_OCEAN_ROUGHNESS)

    ROUGHNESS_TYPE = OCEAN_roughness_TYPE
    CM0            = OCEAN_roughness_CM0
    visck          = OCEAN_roughness_visck
    Z0M_min        = OCEAN_roughness_Z0M_min
    Z0H_min        = OCEAN_roughness_Z0H_min
    Z0E_min        = OCEAN_roughness_Z0E_min
    Z0MI           = OCEAN_roughness_Z0MI
    Z0MR           = OCEAN_roughness_Z0MR
    Z0MS           = OCEAN_roughness_Z0MS
    Z0HI           = OCEAN_roughness_Z0HI
    Z0HR           = OCEAN_roughness_Z0HR
    Z0HS           = OCEAN_roughness_Z0HS
    Z0EI           = OCEAN_roughness_Z0EI
    Z0ER           = OCEAN_roughness_Z0ER
    Z0ES           = OCEAN_roughness_Z0ES

    select case( ROUGHNESS_TYPE )
    case ( 'MILLER92' )
       OCEAN_roughness => OCEAN_roughness_miller92
    case ( 'MOON07' )
       OCEAN_roughness => OCEAN_roughness_moon07
    case default
       write(*,*) 'xxx invalid sea roughness length scheme (', trim(ROUGHNESS_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine OCEAN_roughness_setup

  subroutine OCEAN_roughness_miller92( &
      Z0W,           & ! (inout)
      Z0M, Z0H, Z0E, & ! (out)
      Uabs, ZA       ) ! (in)
    use scale_const, only: &
      GRAV => CONST_GRAV
    implicit none
    
    ! argument
    real(RP), intent(inout) :: Z0W ! roughness length of sea surface [m]

    real(RP), intent(out) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(out) :: Z0H ! roughness length of heat [m]
    real(RP), intent(out) :: Z0E ! roughness length of vapor [m]

    real(RP), intent(in) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: ZA   ! height at 1st atm. layer [m]

    ! work
    real(RP) :: Ustar
    !---------------------------------------------------------------------------

    !--- friction velocity at u, v, and w points
    Ustar = max ( sqrt ( CM0 ) * Uabs , Ustar_min )

    !--- roughness lengths at u, v, and w points
    Z0M = max( Z0MI + Z0MR/GRAV * Ustar*Ustar + Z0MS*visck / Ustar, Z0M_min )
    Z0H = max( Z0HI + Z0HR/GRAV * Ustar*Ustar + Z0HS*visck / Ustar, Z0H_min )
    Z0E = max( Z0EI + Z0ER/GRAV * Ustar*Ustar + Z0ES*visck / Ustar, Z0E_min )

    return
  end subroutine OCEAN_roughness_miller92

  subroutine OCEAN_roughness_moon07( &
      Z0W,           & ! (inout)
      Z0M, Z0H, Z0E, & ! (out)
      Uabs, ZA       ) ! (in)
    use scale_const, only: &
      GRAV   => CONST_GRAV,   &
      KARMAN => CONST_KARMAN
    implicit none

    ! argument
    real(RP), intent(inout) :: Z0W ! roughness length of sea surface [m]

    real(RP), intent(out) :: Z0M ! roughness length of momentum [m]
    real(RP), intent(out) :: Z0H ! roughness length of heat [m]
    real(RP), intent(out) :: Z0E ! roughness length of vapor [m]

    real(RP), intent(in) :: Uabs ! absolute velocity at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: ZA   ! height at 1st atm. layer [m]

    ! parameter
    integer, parameter :: nmax = 10 ! maximum iteration number

    ! work
    integer :: n

    real(RP) :: Ustar, U10M
    !---------------------------------------------------------------------------

    ! IL-JU MOON, ISAAC GINIS, TETSU HARA, AND BIJU THOMAS, 2007:
    ! A Physics-Based Parameterization of Air-Sea Momentum Flux at High Wind Speeds
    ! and Its Impact on Hurricane Intensity Predictions, Mon. Wea. Rev., 135, 2869-2878

    Z0M = Z0W
    do n = 1, nmax
      Ustar = KARMAN * Uabs / log( ZA / Z0M )
      U10M  = Ustar / KARMAN * log( 10.0_RP / Z0M )

      if( U10M <= 12.5_RP ) then
        Z0M = 0.0185_RP * Ustar**2 / GRAV
      else
        Z0M = 1.0E-3_RP * ( 0.085_RP * ( -0.56_RP * Ustar**2 + 20.255_RP * Ustar + 2.458_RP) - 0.58_RP )
      end if
    end do

    !  Fairall et al. TOGA V3.0
    !  Fairall et al. (2003) JCLI, vol. 16, 571-591. Eq. (28)
    Z0H = min( 5.5E-5_RP * ( Z0M * Ustar / visck )**(-0.6_RP), 1.1E-4_RP )
    Z0E = min( 5.5E-5_RP * ( Z0M * Ustar / visck )**(-0.6_RP), 1.1E-4_RP ) ! Z0H = Z0E

    ! update Z0W
    Z0W = Z0M

    return
  end subroutine OCEAN_roughness_moon07

end module scale_ocean_roughness
