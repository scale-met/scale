!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User module for the Rossby wave experiment
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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_finalize
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
  real(RP), private :: JET_DRHO =  1.0E-4_RP
  real(RP), private :: JET_CZ =    5.0E3_RP
  real(RP), private :: JET_RZ =    5.0E3_RP
  real(RP), private :: JET_RY = 2000.0E3_RP

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
         JET_DRHO, &
         JET_CZ, &
         JET_RZ, &
         JET_RY

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'
    LOG_INFO("USER_setup",*) 'Rosby wave experiment'

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
  !> Finalization
  subroutine USER_finalize
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_finalize

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    use scale_const, only: &
       PI   => CONST_PI, &
       GRAV => CONST_GRAV, &
       P00  => CONST_PRE00, &
       R    => CONST_Rdry, &
       CP   => CONST_CPdry, &
       CV   => CONST_CVdry
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
       FXG  => ATMOS_GRID_CARTESC_FXG, &
       FZ   => ATMOS_GRID_CARTESC_FZ, &
       CZ   => ATMOS_GRID_CARTESC_CZ, &
       CX   => ATMOS_GRID_CARTESC_CX, &
       CY   => ATMOS_GRID_CARTESC_CY, &
       RCDX => ATMOS_GRID_CARTESC_RCDX, &
       RCDY => ATMOS_GRID_CARTESC_RCDY
    use scale_coriolis, only: &
       CORIOLIS_f, &
       CORIOLIS_beta
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       RHOT
    implicit none

    real(RP) :: pres(KA,IA,JA), pres_toa, temp
    real(RP) :: Lx, dy
    real(RP) :: RovCP
    real(RP) :: wn_k, wn_l !> wave number
    real(RP) :: c          !> phase speed

    integer  :: k, i, j

    !---------------------------------------------------------------------------

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_vars8( RHOT(:,:,:), 2 )
    call COMM_wait ( DENS(:,:,:), 1 )
    call COMM_wait ( RHOT(:,:,:), 2 )

    RovCP = R / CP

    ! pressure at TOA
    pres_toa = P00 * ( R * RHOT(KE,ISB,JSB) / P00 )**( CP/CV ) & ! at k=KE
             - GRAV * DENS(KE,1,1) * ( FZ(KE) - CZ(KE) )

    ! wave numbers
    Lx = FXG(IAG-IHALO) - FXG(IHALO)
    wn_k = 2.0_RP * PI / Lx
    wn_l = PI / JET_RY
    c = CORIOLIS_beta / ( wn_k**2 + wn_l**2 )

    if( IO_L ) then
       LOG_NEWLINE
       LOG_INFO("USER_mkinit",*) 'Rossby wave experiment on a beta-plane'
       LOG_INFO("USER_mkinit",*) 'Wave number [x-direction]: ', wn_k
       LOG_INFO("USER_mkinit",*) 'Wave number [y-direction]: ', wn_l
       LOG_INFO("USER_mkinit",*) 'Phase speed [x-direction]: ', c
    end if

    ! initial profile
    do j = 1, JA
       dy = CY(j) - ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y
       if ( abs(dy) <= JET_RY ) then
          do i = 1, IA
          do k = KS, KE
             DENS(k,i,j) = DENS(k,i,j) &
                         * ( 1.0_RP + JET_DRHO * cos( 0.5_RP * wn_l * dy )**2 &
                                    * sin( wn_k * CX(i) ) &
                                    * exp( - ( ( CZ(k) - JET_CZ ) / JET_CZ )**2 ) )
          end do
          end do
       end if
    end do

    ! hydrostatic balance
    do j = 1, JA
    do i = 1, IA
       pres(KE,i,j) = pres_toa + GRAV * DENS(KE,i,j) * ( FZ(KE) - CZ(KE) )
       do k = KE-1, KS, -1
          pres(k,i,j) = pres(k+1,i,j) + GRAV * ( DENS(k,i,j) + DENS(k+1,i,j) ) * ( CZ(k+1) - CZ(k) ) * 0.5_RP
          temp = pres(k,i,j) / ( R * DENS(k,i,j) )
          RHOT(k,i,j) = DENS(k,i,j) * temp * ( P00 / pres(k,i,j) )**RovCP
       end do
    end do
    end do

    ! geostrophic balance
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX(k,i,j) = - ( ( pres(k,i,j+1) + pres(k,i+1,j+1) ) &
                       - ( pres(k,i,j-1) + pres(k,i+1,j-1) ) ) * 0.25_RP * RCDY(j) &
                     / CORIOLIS_f(i,j)
       MOMY(k,i,j) =   ( ( pres(k,i+1,j) + pres(k,i+1,j+1) ) &
                       - ( pres(k,i-1,j) + pres(k,i-1,j+1) ) ) * 0.25_RP * RCDX(i) &
                     / CORIOLIS_f(i,j)
    end do
    end do
    end do

    call COMM_vars8( MOMX(:,:,:), 1 )
    call COMM_vars8( MOMY(:,:,:), 2 )
    call COMM_wait ( MOMX(:,:,:), 1 )
    call COMM_wait ( MOMY(:,:,:), 2 )

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
  !> User step
  subroutine USER_update
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_update

end module mod_user
