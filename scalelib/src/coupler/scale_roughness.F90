!-------------------------------------------------------------------------------
!> module Surface roughness length
!!
!! @par Description
!!          surface roughness length
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_roughness
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ROUGHNESS_setup

  abstract interface
     subroutine rl( &
          Z0M_t, &
          Z0H_t, &
          Z0E_t, &
          Z0M,   &
          Z0H,   &
          Z0E,   &
          UA,    &
          VA,    &
          Z1,    &
          dt     )
       use scale_precision
       use scale_grid_index
       implicit none

       real(RP), intent(out) :: Z0M_t(IA,JA) ! tendency of roughness length for momentum [m]
       real(RP), intent(out) :: Z0H_t(IA,JA) ! tendency of roughness length for heat [m]
       real(RP), intent(out) :: Z0E_t(IA,JA) ! tendency of roughness length for vapor [m]

       real(RP), intent(in) :: Z0M(IA,JA) ! roughness length for momentum [m]
       real(RP), intent(in) :: Z0H(IA,JA) ! roughness length for heat [m]
       real(RP), intent(in) :: Z0E(IA,JA) ! roughness length for vapor [m]
       real(RP), intent(in) :: UA (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: VA (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
       real(RP), intent(in) :: Z1 (IA,JA) ! cell center height at the lowest atmospheric layer [m]
       real(DP), intent(in) :: dt         ! delta time
     end subroutine rl
  end interface

  procedure(rl), pointer :: ROUGHNESS => NULL()
  public :: ROUGHNESS

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ROUGHNESS_miller92_setup
  private :: ROUGHNESS_moon07_setup
  private :: ROUGHNESS_const_setup
  private :: ROUGHNESS_miller92
  private :: ROUGHNESS_moon07
  private :: ROUGHNESS_const

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: ROUGHNESS_type = 'MOON07' ! surface roughness length scheme

  real(RP), private :: ROUGHNESS_visck          = 1.5E-5_RP ! kinematic viscosity
  real(RP), private :: ROUGHNESS_Ustar_min      = 1.0E-3_RP ! minimum fiction velocity
  real(RP), private :: ROUGHNESS_Z0M_min        = 1.0E-5_RP ! minimum roughness length for momentum [m]
  real(RP), private :: ROUGHNESS_Z0H_min        = 1.0E-5_RP ! minimum roughness length for heat     [m]
  real(RP), private :: ROUGHNESS_Z0E_min        = 1.0E-5_RP ! minimum roughness length for moisture [m]

  real(RP), private :: ROUGHNESS_miller92_CM0   = 1.0E-3_RP ! bulk coef. for U*
  real(RP), private :: ROUGHNESS_miller92_Z0MI  = 0.0E-0_RP ! base roughness length for momentum
  real(RP), private :: ROUGHNESS_miller92_Z0MR  = 1.8E-2_RP ! rough factor          for momentum
  real(RP), private :: ROUGHNESS_miller92_Z0MS  = 1.1E-1_RP ! smooth factor         for momentum
  real(RP), private :: ROUGHNESS_miller92_Z0HI  = 1.4E-5_RP ! base roughness length for heat
  real(RP), private :: ROUGHNESS_miller92_Z0HR  = 0.0E-0_RP ! rough factor          for heat
  real(RP), private :: ROUGHNESS_miller92_Z0HS  = 4.0E-1_RP ! smooth factor         for heat
  real(RP), private :: ROUGHNESS_miller92_Z0EI  = 1.3E-4_RP ! base roughness length for moisture
  real(RP), private :: ROUGHNESS_miller92_Z0ER  = 0.0E-0_RP ! rough factor          for moisture
  real(RP), private :: ROUGHNESS_miller92_Z0ES  = 6.2E-1_RP ! smooth factor         for moisture

  integer,  private :: ROUGHNESS_moon07_itelim  = 10        ! maximum iteration number

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ROUGHNESS_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ROUGHNESS / &
       ROUGHNESS_type,      &
       ROUGHNESS_visck,     &
       ROUGHNESS_Ustar_min, &
       ROUGHNESS_Z0M_min,   &
       ROUGHNESS_Z0H_min,   &
       ROUGHNESS_Z0E_min

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ROUGHNESS] / Categ[COUPLER] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ROUGHNESS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ROUGHNESS. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ROUGHNESS)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Scheme for ocean roughness length : ', trim(ROUGHNESS_type)
    select case(ROUGHNESS_type)
    case('MILLER92')
       if( IO_L ) write(IO_FID_LOG,*) '*** => Miller (1992)'
       ROUGHNESS => ROUGHNESS_miller92
       call ROUGHNESS_miller92_setup
    case('MOON07')
       if( IO_L ) write(IO_FID_LOG,*) '*** => Moon et al. (2007)'
       ROUGHNESS => ROUGHNESS_moon07
       call ROUGHNESS_moon07_setup
    case('CONST')
       if( IO_L ) write(IO_FID_LOG,*) '*** => Constant.'
       ROUGHNESS => ROUGHNESS_const
       call ROUGHNESS_const_setup
    case default
       write(*,*) 'xxx Unsupported BULKFLUX_type. STOP'
       call PRC_MPIstop
    end select

    return
  end subroutine ROUGHNESS_setup

  !-----------------------------------------------------------------------------
  subroutine ROUGHNESS_miller92_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ROUGHNESS_MILLER92 / &
       ROUGHNESS_miller92_CM0,    &
       ROUGHNESS_miller92_Z0MI,   &
       ROUGHNESS_miller92_Z0MR,   &
       ROUGHNESS_miller92_Z0MS,   &
       ROUGHNESS_miller92_Z0HI,   &
       ROUGHNESS_miller92_Z0HR,   &
       ROUGHNESS_miller92_Z0HS,   &
       ROUGHNESS_miller92_Z0EI,   &
       ROUGHNESS_miller92_Z0ER,   &
       ROUGHNESS_miller92_Z0ES

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ROUGHNESS_MILLER92,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ROUGHNESS_MILLER92. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ROUGHNESS_MILLER92)

    return
  end subroutine ROUGHNESS_miller92_setup

  !-----------------------------------------------------------------------------
  subroutine ROUGHNESS_moon07_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ROUGHNESS_MOON07 / &
       ROUGHNESS_moon07_itelim

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ROUGHNESS_MOON07,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ROUGHNESS_MOON07. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ROUGHNESS_MOON07)

    return
  end subroutine ROUGHNESS_moon07_setup

  !-----------------------------------------------------------------------------
  subroutine ROUGHNESS_const_setup
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine ROUGHNESS_const_setup

  !-----------------------------------------------------------------------------
  subroutine ROUGHNESS_miller92( &
       Z0M_t, & ! [OUT]
       Z0H_t, & ! [OUT]
       Z0E_t, & ! [OUT]
       Z0M,   & ! [IN]
       Z0H,   & ! [IN]
       Z0E,   & ! [IN]
       UA,    & ! [IN]
       VA,    & ! [IN]
       Z1,    & ! [IN]
       dt     ) ! [IN]
    use scale_const, only: &
       GRAV => CONST_GRAV
    implicit none

    ! arguments
    real(RP), intent(out) :: Z0M_t(IA,JA) ! tendency of roughness length for momentum [m]
    real(RP), intent(out) :: Z0H_t(IA,JA) ! tendency of roughness length for heat [m]
    real(RP), intent(out) :: Z0E_t(IA,JA) ! tendency of roughness length for vapor [m]

    real(RP), intent(in) :: Z0M(IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in) :: Z0H(IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E(IA,JA) ! roughness length for vapor [m]
    real(RP), intent(in) :: UA (IA,JA) ! velocity u at the lowest atomspheric layer [m/s]
    real(RP), intent(in) :: VA (IA,JA) ! velocity v at the lowest atomspheric layer [m/s]
    real(RP), intent(in) :: Z1 (IA,JA) ! cell center height at the lowest atmospheric layer [m]
    real(DP), intent(in) :: dt         ! delta time

    ! works
    real(RP) :: Z0M1(IA,JA)
    real(RP) :: Z0H1(IA,JA)
    real(RP) :: Z0E1(IA,JA)

    real(RP) :: Uabs, Ustar

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE

       Uabs  = sqrt( UA(i,j)**2 + VA(i,j)**2 )
       Ustar = max( sqrt( ROUGHNESS_miller92_CM0 ) * Uabs, ROUGHNESS_Ustar_min )

       Z0M1(i,j) = max( ROUGHNESS_miller92_Z0MI &
                      + ROUGHNESS_miller92_Z0MR / GRAV * Ustar * Ustar &
                      + ROUGHNESS_miller92_Z0MS * ROUGHNESS_visck / Ustar, &
                        ROUGHNESS_Z0M_min )
       Z0H1(i,j) = max( ROUGHNESS_miller92_Z0HI &
                      + ROUGHNESS_miller92_Z0HR / GRAV * Ustar * Ustar &
                      + ROUGHNESS_miller92_Z0HS * ROUGHNESS_visck / Ustar, &
                        ROUGHNESS_Z0H_min )
       Z0E1(i,j) = max( ROUGHNESS_miller92_Z0EI &
                      + ROUGHNESS_miller92_Z0ER / GRAV * Ustar * Ustar &
                      + ROUGHNESS_miller92_Z0ES * ROUGHNESS_visck / Ustar, &
                        ROUGHNESS_Z0E_min )

       ! calculate tendency
       Z0M_t(i,j) = ( Z0M1(i,j) - Z0M(i,j) ) / dt
       Z0H_t(i,j) = ( Z0H1(i,j) - Z0H(i,j) ) / dt
       Z0E_t(i,j) = ( Z0E1(i,j) - Z0E(i,j) ) / dt

    enddo
    enddo

    return
  end subroutine ROUGHNESS_miller92

  !-----------------------------------------------------------------------------
  !> IL-JU MOON, ISAAC GINIS, TETSU HARA, AND BIJU THOMAS, 2007:
  !> A Physics-Based Parameterization of Air-Sea Momentum Flux at High Wind Speeds
  !> and Its Impact on Hurricane Intensity Predictions, Mon. Wea. Rev., 135, 2869-2878
  subroutine ROUGHNESS_moon07( &
       Z0M_t, & ! [OUT]
       Z0H_t, & ! [OUT]
       Z0E_t, & ! [OUT]
       Z0M,   & ! [IN]
       Z0H,   & ! [IN]
       Z0E,   & ! [IN]
       UA,    & ! [IN]
       VA,    & ! [IN]
       Z1,    & ! [IN]
       dt     ) ! [IN]
    use scale_const, only: &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN
    implicit none

    ! arguments
    real(RP), intent(out) :: Z0M_t(IA,JA) ! tendency of roughness length for momentum [m]
    real(RP), intent(out) :: Z0H_t(IA,JA) ! tendency of roughness length for heat [m]
    real(RP), intent(out) :: Z0E_t(IA,JA) ! tendency of roughness length for vapor [m]

    real(RP), intent(in) :: Z0M(IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in) :: Z0H(IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E(IA,JA) ! roughness length for vapor [m]
    real(RP), intent(in) :: UA (IA,JA) ! velocity u at the lowest atomspheric layer [m/s]
    real(RP), intent(in) :: VA (IA,JA) ! velocity v at the lowest atomspheric layer [m/s]
    real(RP), intent(in) :: Z1 (IA,JA) ! cell center height at the lowest atmospheric layer [m]
    real(DP), intent(in) :: dt         ! delta time

    ! works
    real(RP) :: Z0M1(IA,JA)
    real(RP) :: Z0H1(IA,JA)
    real(RP) :: Z0E1(IA,JA)

    real(RP) :: Ustar(IA,JA)
    real(RP) :: Uabs (IA,JA)
    real(RP) :: U10M

    integer  :: ite
    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       Z0M1(i,j) = max( Z0M(i,j), ROUGHNESS_Z0M_min )
       Uabs(i,j) = sqrt( UA(i,j)**2 + VA(i,j)**2 )
    enddo
    enddo

    do ite = 1, ROUGHNESS_moon07_itelim
       do j = JS, JE
       do i = IS, IE
          Ustar(i,j) = max( KARMAN * Uabs(i,j) / log( Z1(i,j)/Z0M1(i,j) ), ROUGHNESS_Ustar_min )
          U10M = Ustar(i,j) / KARMAN * log( 10.0_RP/Z0M1(i,j) )

          if ( U10M <= 12.5_RP ) then
             Z0M1(i,j) = max( 0.0185_RP * Ustar(i,j)**2 / GRAV, ROUGHNESS_Z0M_min )
          else
             Z0M1(i,j) = 1.0E-3_RP * ( 0.085_RP * (  -0.56_RP*Ustar(i,j)**2 &
                                                  + 20.255_RP*Ustar(i,j)    &
                                                  +  2.458_RP               ) - 0.58_RP )
          endif
       enddo
       enddo
    enddo

    !  Fairall et al. TOGA V3.0
    !  Fairall et al. (2003) JCLI, vol. 16, 571-591. Eq. (28)
    do j = JS, JE
    do i = IS, IE
       Z0H1(i,j) = min( 5.5E-5_RP / ( Z0M1(i,j) * Ustar(i,j) / ROUGHNESS_visck )**0.6_RP, 1.1E-4_RP )
       Z0E1(i,j) = Z0H1(i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       ! limiter
       Z0M1(i,j) = max( Z0M1(i,j), ROUGHNESS_Z0M_min )
       Z0H1(i,j) = max( Z0H1(i,j), ROUGHNESS_Z0H_min )
       Z0E1(i,j) = max( Z0E1(i,j), ROUGHNESS_Z0E_min )

       ! calculate tendency
       Z0M_t(i,j) = ( Z0M1(i,j) - Z0M(i,j) ) / dt
       Z0H_t(i,j) = ( Z0H1(i,j) - Z0H(i,j) ) / dt
       Z0E_t(i,j) = ( Z0E1(i,j) - Z0E(i,j) ) / dt
    enddo
    enddo

    return
  end subroutine ROUGHNESS_moon07

  !-----------------------------------------------------------------------------
  subroutine ROUGHNESS_const( &
       Z0M_t, & ! [OUT]
       Z0H_t, & ! [OUT]
       Z0E_t, & ! [OUT]
       Z0M,   & ! [IN]
       Z0H,   & ! [IN]
       Z0E,   & ! [IN]
       UA,    & ! [IN]
       VA,    & ! [IN]
       Z1,    & ! [IN]
       dt     ) ! [IN]
    implicit none

    ! arguments
    real(RP), intent(out) :: Z0M_t(IA,JA) ! tendency of roughness length for momentum [m]
    real(RP), intent(out) :: Z0H_t(IA,JA) ! tendency of roughness length for heat [m]
    real(RP), intent(out) :: Z0E_t(IA,JA) ! tendency of roughness length for vapor [m]

    real(RP), intent(in) :: Z0M(IA,JA) ! roughness length for momentum [m]
    real(RP), intent(in) :: Z0H(IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E(IA,JA) ! roughness length for vapor [m]
    real(RP), intent(in) :: UA (IA,JA) ! velocity u at the lowest atomspheric layer [m/s]
    real(RP), intent(in) :: VA (IA,JA) ! velocity v at the lowest atomspheric layer [m/s]
    real(RP), intent(in) :: Z1 (IA,JA) ! cell center height at the lowest atmospheric layer [m]
    real(DP), intent(in) :: dt         ! delta time

    Z0M_t(:,:) = 0.0_RP
    Z0H_t(:,:) = 0.0_RP
    Z0E_t(:,:) = 0.0_RP

    return
  end subroutine ROUGHNESS_const

end module scale_roughness
