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
  logical,  private :: USER_do   = .false. !< do user step?

  real(RP), private :: SFC_THETA = 300.0_RP ! surface potential temperature [K]
  real(RP), private :: ENV_BVF   =  0.01_RP ! Brunt Vaisala frequencies of environment [1/s]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine USER_tracer_setup
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
       USER_do,   &
       SFC_THETA, &
       ENV_BVF

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
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_atmos_grid_cartesC, only : &
       CZ => ATMOS_GRID_CARTESC_CZ
    use mod_atmos_vars, only: &
       DENS, &
       RHOT
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none

    real(RP) :: PT_diff(KA,IA,JA)

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PT_diff(k,i,j) = RHOT(k,i,j)/DENS(k,i,j) - SFC_THETA * exp( ENV_BVF*ENV_BVF / GRAV * CZ(k) )
       enddo
       enddo
       enddo

       call FILE_HISTORY_in( PT_diff(:,:,:), 'PT_diff', 'PT perturbation', 'K' )
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
