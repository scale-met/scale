!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          calc perturbation
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-09-05 (H.Yashiro)   [new]
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
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_setup0
  public :: USER_setup
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
  logical,  private, save :: USER_do = .false. !< do user step?

  real(RP), private, save :: SFC_THETA = 300.0_RP ! surface potential temperature [K]
  real(RP), private, save :: ENV_BVF   =  0.01_RP ! Brunt Vaisala frequencies of environment [1/s]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup0
  subroutine USER_setup0
  end subroutine USER_setup0

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_USER / &
       USER_do,   &
       SFC_THETA, &
       ENV_BVF

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

    ! run once (only for the diagnostic value)
    call USER_step

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       GRAV  => CONST_GRAV
    use scale_grid, only : &
       CZ => GRID_CZ
    use scale_time, only: &
       DTSEC => TIME_DTSEC
    use mod_atmos_vars, only: &
       DENS, &
       RHOT
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP) :: PT_diff(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PT_diff(k,i,j) = RHOT(k,i,j)/DENS(k,i,j) - SFC_THETA * exp( ENV_BVF*ENV_BVF / GRAV * CZ(k) )
       enddo
       enddo
       enddo

       call HIST_in( PT_diff(:,:,:), 'PT_diff', 'PT perturbation', 'K' )
    endif

    return
  end subroutine USER_step

end module mod_user
