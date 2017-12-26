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
module mod_user
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
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  !> Config
  subroutine USER_config
    use scale_process, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
         JET_DRHO, &
         JET_CZ, &
         JET_RZ, &
         JET_RY

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Rosby wave experiment'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    implicit none

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    use scale_const, only: &
       PI   => CONST_PI, &
       GRAV => CONST_GRAV, &
       P00  => CONST_PRE00, &
       R    => CONST_Rdry, &
       CP   => CONST_CPdry, &
       CV   => CONST_CVdry
    use scale_grid, only: &
       GRID_DOMAIN_CENTER_Y, &
       FXG  => GRID_FXG, &
       FZ   => GRID_FZ, &
       CZ   => GRID_CZ, &
       CX   => GRID_CX, &
       CY   => GRID_CY, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use scale_atmos_dyn, only: &
       CORIOLIS
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       RHOT
    use mod_atmos_dyn_driver, only: &
       ATMOS_DYN_coriolis_beta
    implicit none

    real(RP) :: pres(KA,IA,JA), pres_toa, temp
    real(RP) :: Lx, dy
    real(RP) :: RovCP
    real(RP) :: wn_k, wn_l !> wave number
    real(RP) :: c          !> phase speed
    integer :: k, i, j

    RovCP = R / CP

    ! pressure at TOA
    pres_toa = P00 * ( R * RHOT(KE,1,1) / P00 )**( CP/CV ) & ! at k=KE
             - GRAV * DENS(KE,1,1) * ( FZ(KE) - CZ(KE) )

    ! wave numbers
    Lx = FXG(IAG-IHALO) - FXG(IHALO)
    wn_k = 2.0_RP * PI / Lx
    wn_l = PI / JET_RY
    c = ATMOS_DYN_coriolis_beta / ( wn_k**2 + wn_l**2 )

    if( IO_L ) then
       write(IO_FID_LOG,*) '+++ Wave number [x-direction]: ', wn_k
       write(IO_FID_LOG,*) '+++ Wave number [y-direction]: ', wn_l
       write(IO_FID_LOG,*) '+++ Phase speed [x-direction]: ', c
    end if

    ! initial profile
    do j = 1, JA
       dy = CY(j) - GRID_DOMAIN_CENTER_Y
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
                     / CORIOLIS(i,j)
       MOMY(k,i,j) =   ( ( pres(k,i+1,j) + pres(k,i+1,j+1) ) &
                       - ( pres(k,i-1,j) + pres(k,i-1,j+1) ) ) * 0.25_RP * RCDX(i) &
                     / CORIOLIS(i,j)
    end do
    end do
    end do

    call COMM_vars8(MOMX, 1)
    call COMM_vars8(MOMY, 2)
    call COMM_wait (MOMX, 1)
    call COMM_wait (MOMY, 2)

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_step

end module mod_user
