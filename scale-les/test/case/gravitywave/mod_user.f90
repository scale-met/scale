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
  use mod_precision
  use mod_index
  use mod_stdio
  use mod_prof
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
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
  !> Setup
  subroutine USER_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only: &
       GRAV  => CONST_GRAV
    use mod_grid, only : &
       CZ => GRID_CZ
    use mod_time, only: &
       DTSEC => TIME_DTSEC
    use mod_atmos_vars, only: &
       DENS, &
       RHOT
    use mod_history, only: &
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

       call HIST_in( PT_diff(:,:,:), 'PT_diff', 'PT perturbation', 'K', DTSEC)
    endif

    return
  end subroutine USER_step

end module mod_user
