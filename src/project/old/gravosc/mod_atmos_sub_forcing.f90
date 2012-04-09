!-------------------------------------------------------------------------------
!> module Atmosphere / adiitional forcing
!!
!! @par Description
!!          Adiitional forcing
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-01-26 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_forcing
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_FORCING_setup
  public :: ATMOS_FORCING
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(8), private, save :: OSC_AMP_POTT = 15.D0
  real(8), private, save :: OSC_FREQ_MIN = 10.D0
  real(8), private, save :: ENV_THETA    = 300.D0
  real(8), private, save :: ZC_BBL = 3.D3      ! center location [m]: z
  real(8), private, save :: XC_BBL = 18.D3     ! center location [m]: x
  real(8), private, save :: YC_BBL = 18.D3     ! center location [m]: y
  real(8), private, save :: ZR_BBL = 2.D3      ! bubble radius   [m]: z
  real(8), private, save :: XR_BBL = 4.D3      ! bubble radius   [m]: x
  real(8), private, save :: YR_BBL = 4.D3      ! bubble radius   [m]: y

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Initialize Reference state
  !-----------------------------------------------------------------------------
  subroutine ATMOS_FORCING_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_ATMOS_FORCING / &
       OSC_AMP_POTT, &
       OSC_FREQ_MIN, &
       ENV_THETA, &
       ZC_BBL,    &
       XC_BBL,    &
       YC_BBL,    &
       ZR_BBL,    &
       XR_BBL,    &
       YR_BBL

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[FORCING]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_FORCING,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_FORCING. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_FORCING)

    return
  end subroutine ATMOS_FORCING_setup

  !-----------------------------------------------------------------------------
  !> FORCING
  !-----------------------------------------------------------------------------
  subroutine ATMOS_FORCING
    use mod_const, only : &
       PI => CONST_PI
    use mod_time, only: &
       TIME_NOWSEC
    use mod_grid, only : &
       KA => GRID_KA, &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KS => GRID_KS, &
       KE => GRID_KE, &
       IS => GRID_IS, &
       IE => GRID_IE, &
       JS => GRID_JS, &
       JE => GRID_JE, &
       GRID_CZ,       &
       GRID_CX,       &
       GRID_CY
    use mod_atmos_vars, only: &
       var => atmos_var, &
       I_DENS,      &
       I_RHOT
    implicit none

    real(8) :: dist, force

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dist = ( (GRID_CZ(k)-ZC_BBL)/ZR_BBL )**2 &
            + ( (GRID_CX(i)-XC_BBL)/XR_BBL )**2 &
            + ( (GRID_CY(j)-YC_BBL)/YR_BBL )**2

       if ( dist <= 1.D0 ) then
          force = OSC_AMP_POTT * dsin( 2.D0 * PI * TIME_NOWSEC / (OSC_FREQ_MIN*60.D0) )

          var(k,i,j,I_RHOT) = ( ENV_THETA + force * dcos( 0.5D0*PI*dsqrt(dist) )**2 ) &
                            * var(k,i,j,I_DENS)
       endif
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_FORCING

end module mod_atmos_forcing
