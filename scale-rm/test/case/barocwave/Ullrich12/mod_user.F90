!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Set coriolis parameter and boundary conditions for baroclinic wave in a channel.
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
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer

  use scale_time, only: &
       DTSEC => TIME_DTSEC, &
       NOWTSEC => TIME_NOWSEC

  use scale_const, only: &
       PI  => CONST_PI,    &
       OHM => CONST_OHM,   &
       RPlanet => CONST_RADIUS, &
       Rdry => CONST_Rdry

  use scale_atmos_refstate, only: &
       ATMOS_REFSTATE_pres, &
       ATMOS_REFSTATE_temp, &
       ATMOS_REFSTATE_dens, &
       ATMOS_REFSTATE_pott, &
       ATMOS_REFSTATE_qv,   &
       ATMOS_REFSTATE_write

  use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       PRES

  use scale_prc

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

  real(RP), private, allocatable :: RHOT_bc(:,:,:)
  real(RP), private, allocatable :: DENS_bc(:,:,:)

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
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'
    LOG_INFO("USER_setup",*) 'User procedure in test/case/barocwave/Ullrich12'

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

    !
    allocate( RHOT_bc(KA,IA,2) )
    allocate( DENS_bc(KA,IA,2) )


    ! Save some information of inital fields to set boundary conditions.
    !
    RHOT_bc(:,:,1) = RHOT(:,:,JS) - 0.5_RP*(RHOT(:,:,JS+1) - RHOT(:,:,JS))
    RHOT_bc(:,:,2) = RHOT(:,:,JE-1) + 1.5_RP*(RHOT(:,:,JE) - RHOT(:,:,JE-1))

    DENS_bc(:,:,1) = DENS(:,:,JS) - 0.5_RP*(DENS(:,:,JS+1) - DENS(:,:,JS))
    DENS_bc(:,:,2) = DENS(:,:,JE-1) + 1.5_RP*(DENS(:,:,JE) - DENS(:,:,JE-1))


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
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_prc_cartesC, only: &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_atmos_grid_cartesC, only : &
       CX => ATMOS_GRID_CARTESC_CX, &
       CZ => ATMOS_GRID_CARTESC_CZ
    use scale_time, only: &
       DTSEC => TIME_DTSEC
    use scale_file_history, only: &
       FILE_HISTORY_in
    implicit none


    integer :: k, i, j

    !---------------------------------------------------------------------------

    if ( .not. USER_do ) return

    ! Apply the boundary condition at y=+Ly and y=-Ly

    if ( .NOT. PRC_HAS_N ) then
       MOMY(:,:,JE)   = 0.0_RP
       do j = 1, JHALO
          MOMY(:,:,JE+j  ) = - MOMY(:,:,JE-j  )
          DENS(:,:,JE+j) =  2.0_RP*DENS_bc(:,:,2) - DENS(:,:,JE-j+1)
          MOMX(:,:,JE+j) = - MOMX(:,:,JE-j+1)
          MOMZ(:,:,JE+j) = - MOMZ(:,:,JE-j+1)
          RHOT(:,:,JE+j) = 2.0_RP*RHOT_bc(:,:,2) - RHOT(:,:,JE-j+1)
       enddo
    end if

    if ( .NOT. PRC_HAS_S ) then
       MOMY(:,:,JS-1) = 0.0_RP
       do j = 1, JHALO
          if ( j < JHALO ) MOMY(:,:,JS-j-1) = - MOMY(:,:,JS+j-1)
          DENS(:,:,JS-j) = 2.0_RP*DENS_bc(:,:,1) - DENS(:,:,JS+j-1)
          MOMX(:,:,JS-j) = - MOMX(:,:,JS+j-1)
          MOMZ(:,:,JS-j) = - MOMZ(:,:,JS+j-1)
          RHOT(:,:,JS-j) = 2.0_RP*RHOT_bc(:,:,1) - RHOT(:,:,JS+j-1)
       enddo
    end if

    return
  end subroutine USER_update

end module mod_user
