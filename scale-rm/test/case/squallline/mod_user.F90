!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Cold pool forcing
!!          Redelsperger et al. (2000) Q.J.R.Meteorol.Soc., 126, pp.823-863
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  real(DP), private, save :: TIME0
  real(RP), private, save :: pi2
  integer,  private, save :: Ktop

  logical,  private, save :: USER_do  = .true.
  real(DP), private, save :: FORCE_DURATION = 1200.D0
  real(RP), private, save :: SHIFT_X = 12.0E0_RP
  real(RP), private, save :: SHIFT_Y = -2.0E0_RP
  real(RP), private, save :: DT_MAX  = -6.7e-3_RP
  real(RP), private, save :: DQ_MAX  = -1.675e-6_RP
  real(RP), private, save :: POOL_TOP  = 2.5e3_RP
  real(RP), private, save :: POOL_CX   = 100.e3_RP
  real(RP), private, save :: POOL_CY0  = 100.e3_RP
  real(RP), private, save :: POOL_RX   = 7.e3_RP
  real(RP), private, save :: POOL_RY   = 6.e3_RP
  real(RP), private, save :: POOL_DIST = 15.e3_RP
  integer,  private, save :: POOL_NUM  = 4


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine USER_tracer_setup

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    use scale_time, only: &
       NOWSEC => TIME_NOWSEC
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ
    implicit none

    namelist / PARAM_USER / &
       USER_do,        &
       FORCE_DURATION, &
       DT_MAX,         &
       DQ_MAX,         &
       SHIFT_X,        &
       SHIFT_Y,        &
       POOL_CX,        &
       POOL_CY0,       &
       POOL_TOP,       &
       POOL_RX,        &
       POOL_RY,        &
       POOL_DIST,      &
       POOL_NUM

    integer :: k
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    if ( USER_do ) then
       LOG_INFO("USER_setup",*) 'Enable cold pool forcing'
       TIME0 = NOWSEC

       pi2 = atan(1.0) * 2.0_RP

       Ktop = KE
       do k = KS, KE
          if ( CZ(k) > POOL_TOP ) then
             Ktop = k-1
             exit
          endif
       enddo
    endif

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calc tendency
  subroutine USER_calc_tendency
    use scale_time, only: &
       NOWSEC => TIME_NOWSEC, &
       DTSEC  => TIME_DTSEC
    use mod_atmos_vars, only: &
       DENS, &
       QTRC, &
       DENS_tp, &
       RHOT_tp, &
       RHOQ_tp
    use scale_atmos_grid_cartesC, only: &
       CX => ATMOS_GRID_CARTESC_CX, &
       CY => ATMOS_GRID_CARTESC_CY
    use scale_atmos_hydrometeor, only: &
       I_QV
    implicit none

    real(RP) :: dt, dq
    real(RP) :: time
    real(RP) :: fact, dist
    real(RP) :: POOL_CY

    integer :: k, i, j, n
    !---------------------------------------------------------------------------

    if ( USER_do ) then

       time = NOWSEC - TIME0

       if ( time <= FORCE_DURATION ) then
          do j = JS, JE
          do i = IS, IE
             dt = 0.0_RP
             dq = 0.0_RP

             do n = 1, POOL_NUM
                POOL_CY = POOL_CY0 - real(2*n-POOL_NUM-1) * 0.5 * POOL_DIST
                dist = ( (CX(i)-POOL_CX+SHIFT_X*time)/POOL_RX )**2 &
                     + ( (CY(j)-POOL_CY+SHIFT_Y*time)/POOL_RY )**2

                if ( dist < 1.0_RP ) then
                   fact = cos( pi2*sqrt(dist) )**2
                   dt = dt + DT_MAX * fact
                   dq = dq + DQ_MAX * fact
                endif
             enddo

             do k = KS, Ktop
                RHOT_tp(k,i,j)      = RHOT_tp(k,i,j)      + dt * DENS(k,i,j)

                dq = max( dq, -QTRC(k,i,j,I_QV)/DTSEC )
                RHOQ_tp(k,i,j,I_QV) = RHOQ_tp(k,i,j,I_QV) + dq * DENS(k,i,j)
                DENS_tp(k,i,j)      = DENS_tp(k,i,j)      + dq * DENS(k,i,j)
             enddo

          enddo
          enddo

       endif

    endif

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_update
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_update

end module mod_user
