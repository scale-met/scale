!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Calculate the perturbations from background fields in mountain wave test case.  
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

  use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       RHOT
  

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

  real(RP), private, allocatable :: init_U(:,:,:)
  real(RP), private, allocatable :: init_PT(:,:,:)
  
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
    implicit none

    !---------------------------------------------------------------------------

    
    ! calculate diagnostic value and input to history buffer

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_update
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_atmos_grid_cartesC, only : &
       CZ => ATMOS_GRID_CARTESC_CZ
    use scale_time, only: &
       NOWSEC  => TIME_NOWSEC, &
       DTSEC => TIME_DTSEC
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none

    real(RP) :: PT_diff(KA,IA,JA)
    real(RP) :: U_diff(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       if (.not. allocated( init_U )) then
          allocate( init_U(KA,IA,JA) )
          allocate( init_PT(KA,IA,JA) )
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE    
             init_U(k,i,j) = MOMX(k,i,j) * 2.0_RP / (DENS(k,i,j) + DENS(k,i+1,j))
             init_PT(k,i,j) = RHOT(k,i,j)/DENS(k,i,j)
          enddo
          enddo
          enddo
       end if
       
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PT_diff(k,i,j) = RHOT(k,i,j)/DENS(k,i,j) - init_PT(k,i,j)
          U_diff(k,i,j) = MOMX(k,i,j) * 2.0_RP / (DENS(k,i,j) + DENS(k,i+1,j)) - init_U(k,i,j)
       enddo
       enddo
       enddo

       call FILE_HISTORY_in( PT_diff(:,:,:), 'PT_diff', 'PT perturbation', 'K' )
       call FILE_HISTORY_in( U_diff(:,:,:), 'U_diff', 'U perturbation', 'm/s', dim_type='ZXHY' )
    endif

    return
  end subroutine USER_update

end module mod_user
