!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          calc perturbation
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
  logical,  private, save :: USER_do = .false. !< do user step?

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
    implicit none

    namelist / PARAM_USER / &
       USER_do

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
    LOG_NML(PARAM_USER)

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
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_update
    use scale_atmos_grid_cartesC_real, only : &
       CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_time, only : &
       NOWSEC => TIME_NOWSEC
    use mod_atmos_vars, only: &
       QTRC
    implicit none

    integer  :: modsec
    real(RP) :: dist

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       LOG_INFO("USER_update",*) 'Add rain.'
       USER_do = .false.

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          dist = ( ( CZ(k,i,j) - 3000.0_RP ) / 50.0_RP )**2

          QTRC(k,i,j,I_QR) = 1.E-3 / ( 1.0_RP + dist )
!          if (       FZ(k,  i,j) >= 3000.0_RP &
!               .AND. FZ(k-1,i,j) <  3000.0_RP ) then
!
!             LOG_INFO("USER_update",*) k,i,j,CZ(k,i,j),dist,1.E-3/( 1.0_RP + dist )
!             QTRC(k,i,j,I_QR) = 1.E-3
!          endif
       enddo
       enddo
       enddo
    endif

    return
  end subroutine USER_update

end module mod_user
