!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User module for the inertial oscillation test case
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
  ! geostrophic wind
  real(RP), private :: Ug   = 10.0_RP
  real(RP), private :: Vg   =  0.0_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine USER_tracer_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
         Ug, &
         Vg

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_tracer_setup",*) 'Setup'
    LOG_INFO("USER_tracer_setup",*) 'Inertial oscillation experiment'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_tracer_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_tracer_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    implicit none

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
  !> Calculate tendency
  subroutine USER_calc_tendency
    use scale_atmos_dyn, only: &
       CORIOLIS
    use mod_atmos_vars, only: &
       DENS,    &
       RHOU_tp, &
       RHOV_tp
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE

       ! geostrophic forcing
       do k = KS, KE
          RHOU_tp(k,i,j) = RHOU_tp(k,i,j) - CORIOLIS(i,j) * Vg * DENS(k,i,j)
          RHOV_tp(k,i,j) = RHOV_tp(k,i,j) + CORIOLIS(i,j) * Ug * DENS(k,i,j)
       end do

    end do
    end do

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_update
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_update

end module mod_user
