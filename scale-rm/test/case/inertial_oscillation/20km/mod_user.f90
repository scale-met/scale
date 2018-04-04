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
  real(RP), private :: Ug   = 10.0_RP
  real(RP), private :: Vg   =  0.0_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
         Ug, &
         Vg

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Inertial oscillation experiment'

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
    implicit none

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
  end subroutine USER_step

end module mod_user
