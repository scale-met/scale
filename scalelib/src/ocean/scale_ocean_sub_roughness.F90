!-------------------------------------------------------------------------------
!> module OCEAN / Surface roughness length
!!
!! @par Description
!!          roughness length of sea surface
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-28 (T.Yamaura)   [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_ocean_roughness
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
  public :: OCEAN_roughness_setup

  abstract interface
     subroutine rl( &
          Z0,  &
          Z0M, &
          Z0H, &
          Z0E, &
          UA,  &
          VA,  &
          WA   )
       use scale_precision
       use scale_grid_index
       implicit none

       real(RP), intent(inout) :: Z0 (IA,JA) ! roughness length for save     [m]
       real(RP), intent(out)   :: Z0M(IA,JA) ! roughness length for momentum [m]
       real(RP), intent(out)   :: Z0H(IA,JA) ! roughness length for heat     [m]
       real(RP), intent(out)   :: Z0E(IA,JA) ! roughness length for moisture [m]
       real(RP), intent(in)    :: UA (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
       real(RP), intent(in)    :: VA (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
       real(RP), intent(in)    :: WA (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
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
  private :: OCEAN_roughness_const_setup
  private :: OCEAN_roughness_miller92_setup
  private :: OCEAN_roughness_moon07_setup
  private :: OCEAN_roughness_const
  private :: OCEAN_roughness_miller92
  private :: OCEAN_roughness_moon07

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_SHORT), private :: OCEAN_roughness_TYPE = 'MILLER92' ! sea roughness length scheme

  real(RP), private :: OCEAN_roughness_visck          = 1.5E-5_RP ! kinematic viscosity
  real(RP), private :: OCEAN_roughness_Ustar_min      = 1.0E-3_RP ! minimum fiction velocity
  real(RP), private :: OCEAN_roughness_Z0M_min        = 1.0E-5_RP ! minimum roughness length for momentum [m]
  real(RP), private :: OCEAN_roughness_Z0H_min        = 1.0E-5_RP ! minimum roughness length for heat     [m]
  real(RP), private :: OCEAN_roughness_Z0E_min        = 1.0E-5_RP ! minimum roughness length for moisture [m]

  real(RP), private :: OCEAN_roughness_const_Z0M      = 1.0E-5_RP ! constant roughness length for momentum [m]
  real(RP), private :: OCEAN_roughness_const_Z0H      = 1.0E-5_RP ! constant roughness length for momentum [m]
  real(RP), private :: OCEAN_roughness_const_Z0E      = 1.0E-5_RP ! constant roughness length for momentum [m]

  real(RP), private :: OCEAN_roughness_miller92_CM0   = 1.0E-3_RP ! bulk coef. for U*
  real(RP), private :: OCEAN_roughness_miller92_Z0MI  = 0.0E-0_RP ! base roughness length for momentum
  real(RP), private :: OCEAN_roughness_miller92_Z0MR  = 1.8E-2_RP ! rough factor          for momentum
  real(RP), private :: OCEAN_roughness_miller92_Z0MS  = 1.1E-1_RP ! smooth factor         for momentum
  real(RP), private :: OCEAN_roughness_miller92_Z0HI  = 1.4E-5_RP ! base roughness length for heat
  real(RP), private :: OCEAN_roughness_miller92_Z0HR  = 0.0E-0_RP ! rough factor          for heat
  real(RP), private :: OCEAN_roughness_miller92_Z0HS  = 4.0E-1_RP ! smooth factor         for heat
  real(RP), private :: OCEAN_roughness_miller92_Z0EI  = 1.3E-4_RP ! base roughness length for moisture
  real(RP), private :: OCEAN_roughness_miller92_Z0ER  = 0.0E-0_RP ! rough factor          for moisture
  real(RP), private :: OCEAN_roughness_miller92_Z0ES  = 6.2E-1_RP ! smooth factor         for moisture

  integer,  private :: OCEAN_roughness_moon07_itelim  = 10        ! maximum iteration number

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_roughness_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_OCEAN_ROUGHNESS / &
       OCEAN_roughness_TYPE,      &
       OCEAN_roughness_visck,     &
       OCEAN_roughness_Ustar_min, &
       OCEAN_roughness_Z0M_min,   &
       OCEAN_roughness_Z0H_min,   &
       OCEAN_roughness_Z0E_min

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Ocean roughness length parameter'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_ROUGHNESS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_ROUGHNESS. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_ROUGHNESS)

    select case( OCEAN_roughness_TYPE )
    case ('CONST')
       OCEAN_roughness => OCEAN_roughness_const
       call OCEAN_roughness_const_setup
    case ('MILLER92')
       OCEAN_roughness => OCEAN_roughness_miller92
       call OCEAN_roughness_miller92_setup
    case ('MOON07')
       OCEAN_roughness => OCEAN_roughness_moon07
       call OCEAN_roughness_moon07_setup
    case default
       write(*,*) 'xxx invalid sea roughness length scheme (', trim(OCEAN_roughness_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

    return
  end subroutine OCEAN_roughness_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_roughness_const_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_OCEAN_ROUGHNESS_CONST / &
       OCEAN_roughness_const_Z0M,   &
       OCEAN_roughness_const_Z0H,   &
       OCEAN_roughness_const_Z0E

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_ROUGHNESS_CONST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_ROUGHNESS_CONST. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_ROUGHNESS_CONST)

    return
  end subroutine OCEAN_roughness_const_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_roughness_miller92_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_OCEAN_ROUGHNESS_MILLER92 / &
       OCEAN_roughness_miller92_CM0,    &
       OCEAN_roughness_miller92_Z0MI,   &
       OCEAN_roughness_miller92_Z0MR,   &
       OCEAN_roughness_miller92_Z0MS,   &
       OCEAN_roughness_miller92_Z0HI,   &
       OCEAN_roughness_miller92_Z0HR,   &
       OCEAN_roughness_miller92_Z0HS,   &
       OCEAN_roughness_miller92_Z0EI,   &
       OCEAN_roughness_miller92_Z0ER,   &
       OCEAN_roughness_miller92_Z0ES

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_ROUGHNESS_MILLER92,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_ROUGHNESS_MILLER92. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_ROUGHNESS_MILLER92)

    return
  end subroutine OCEAN_roughness_miller92_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_roughness_moon07_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_OCEAN_ROUGHNESS_MOON07 / &
       OCEAN_roughness_moon07_itelim

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_ROUGHNESS_MOON07,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_OCEAN_ROUGHNESS_MOON07. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_OCEAN_ROUGHNESS_MOON07)

    return
  end subroutine OCEAN_roughness_moon07_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_roughness_const( &
       Z0,  &
       Z0M, &
       Z0H, &
       Z0E, &
       UA,  &
       VA,  &
       WA   )
    implicit none

    real(RP), intent(inout) :: Z0 (IA,JA) ! roughness length for save     [m]
    real(RP), intent(out)   :: Z0M(IA,JA) ! roughness length for momentum [m]
    real(RP), intent(out)   :: Z0H(IA,JA) ! roughness length for heat     [m]
    real(RP), intent(out)   :: Z0E(IA,JA) ! roughness length for moisture [m]
    real(RP), intent(in)    :: UA (IA,JA) ! velocity u at the lowest atomspheric layer [m/s]
    real(RP), intent(in)    :: VA (IA,JA) ! velocity v at the lowest atomspheric layer [m/s]
    real(RP), intent(in)    :: WA (IA,JA) ! velocity w at the lowest atomspheric layer [m/s]

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = 1, JA
    do i = 1, IA
       Z0M(i,j) = OCEAN_roughness_const_Z0M
       Z0H(i,j) = OCEAN_roughness_const_Z0H
       Z0E(i,j) = OCEAN_roughness_const_Z0E
    end do
    end do

    return
  end subroutine OCEAN_roughness_const

  !-----------------------------------------------------------------------------
  subroutine OCEAN_roughness_miller92( &
       Z0,  &
       Z0M, &
       Z0H, &
       Z0E, &
       UA,  &
       VA,  &
       WA   )
    use scale_const, only: &
       GRAV => CONST_GRAV
    implicit none

    real(RP), intent(inout) :: Z0 (IA,JA) ! roughness length for save     [m]
    real(RP), intent(out)   :: Z0M(IA,JA) ! roughness length for momentum [m]
    real(RP), intent(out)   :: Z0H(IA,JA) ! roughness length for heat     [m]
    real(RP), intent(out)   :: Z0E(IA,JA) ! roughness length for moisture [m]
    real(RP), intent(in)    :: UA (IA,JA) ! velocity u at the lowest atomspheric layer [m/s]
    real(RP), intent(in)    :: VA (IA,JA) ! velocity v at the lowest atomspheric layer [m/s]
    real(RP), intent(in)    :: WA (IA,JA) ! velocity w at the lowest atomspheric layer [m/s]

    real(RP) :: Uabs, Ustar
    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       Uabs  = sqrt( UA(i,j)**2 + VA(i,j)**2 )
       Ustar = max( sqrt( OCEAN_roughness_miller92_CM0 ) * Uabs, OCEAN_roughness_Ustar_min )

       Z0M(i,j) = max( OCEAN_roughness_miller92_Z0MI &
                     + OCEAN_roughness_miller92_Z0MR / GRAV * Ustar * Ustar &
                     + OCEAN_roughness_miller92_Z0MS * OCEAN_roughness_visck / Ustar, &
                       OCEAN_roughness_Z0M_min )
       Z0H(i,j) = max( OCEAN_roughness_miller92_Z0HI &
                     + OCEAN_roughness_miller92_Z0HR / GRAV * Ustar * Ustar &
                     + OCEAN_roughness_miller92_Z0HS * OCEAN_roughness_visck / Ustar, &
                       OCEAN_roughness_Z0H_min )
       Z0E(i,j) = max( OCEAN_roughness_miller92_Z0EI &
                     + OCEAN_roughness_miller92_Z0ER / GRAV * Ustar * Ustar &
                     + OCEAN_roughness_miller92_Z0ES * OCEAN_roughness_visck / Ustar, &
                       OCEAN_roughness_Z0E_min )
    enddo
    enddo

    return
  end subroutine OCEAN_roughness_miller92

  !-----------------------------------------------------------------------------
  !> IL-JU MOON, ISAAC GINIS, TETSU HARA, AND BIJU THOMAS, 2007:
  !> A Physics-Based Parameterization of Air-Sea Momentum Flux at High Wind Speeds
  !> and Its Impact on Hurricane Intensity Predictions, Mon. Wea. Rev., 135, 2869-2878
  subroutine OCEAN_roughness_moon07( &
       Z0,  &
       Z0M, &
       Z0H, &
       Z0E, &
       UA,  &
       VA,  &
       WA   )
    use scale_const, only: &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN
    use scale_grid_real, only: &
       Z1 => REAL_Z1
    implicit none

    real(RP), intent(inout) :: Z0 (IA,JA) ! roughness length for save     [m]
    real(RP), intent(out)   :: Z0M(IA,JA) ! roughness length for momentum [m]
    real(RP), intent(out)   :: Z0H(IA,JA) ! roughness length for heat     [m]
    real(RP), intent(out)   :: Z0E(IA,JA) ! roughness length for moisture [m]
    real(RP), intent(in)    :: UA (IA,JA) ! velocity u at the lowest atomspheric layer [m/s]
    real(RP), intent(in)    :: VA (IA,JA) ! velocity v at the lowest atomspheric layer [m/s]
    real(RP), intent(in)    :: WA (IA,JA) ! velocity w at the lowest atomspheric layer [m/s]

    real(RP) :: Ustar(IA,JA)
    real(RP) :: Uabs (IA,JA)
    real(RP) :: U10M

    integer  :: ite
    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       Z0M(i,j)  = max( Z0(i,j), OCEAN_roughness_Z0M_min )
       Uabs(i,j) = sqrt( UA(i,j)**2 + VA(i,j)**2 )
    enddo
    enddo

    do ite = 1, OCEAN_roughness_moon07_itelim
       do j = JS, JE
       do i = IS, IE
          Ustar(i,j) = max( KARMAN * Uabs(i,j) / log( Z1(i,j)/Z0M(i,j) ), OCEAN_roughness_Ustar_min )
          U10M = Ustar(i,j) / KARMAN * log( 10.0_RP/Z0M(i,j) )

          if ( U10M <= 12.5_RP ) then
             Z0M(i,j) = max( 0.0185_RP * Ustar(i,j)**2 / GRAV, OCEAN_roughness_Z0M_min )
          else
             Z0M(i,j) = 1.0E-3_RP * ( 0.085_RP * (  -0.56_RP*Ustar(i,j)**2 &
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
       Z0H(i,j) = min( 5.5E-5_RP / ( Z0M(i,j) * Ustar(i,j) / OCEAN_roughness_visck )**0.6_RP, 1.1E-4_RP )
       Z0E(i,j) = Z0H(i,j)
    enddo
    enddo

    ! limiter
    do j = JS, JE
    do i = IS, IE
       Z0M(i,j) = max( Z0M(i,j), OCEAN_roughness_Z0M_min )
       Z0H(i,j) = max( Z0H(i,j), OCEAN_roughness_Z0H_min )
       Z0E(i,j) = max( Z0E(i,j), OCEAN_roughness_Z0E_min )
       ! update and save Z0
       Z0 (i,j) = Z0M(i,j)
    enddo
    enddo

    return
  end subroutine OCEAN_roughness_moon07

end module scale_ocean_roughness
