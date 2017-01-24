!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          Calculate the perturbations from background fields in mountain wave test case.  
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-09-30 (Y.Kawai)     [modify]
!! @li      2016-??-?? (Y.Kawai)     [copied from test case of gravity wave]
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
  logical,  private, save :: USER_do = .false. !< do user step?

  real(RP), private, allocatable :: init_U(:,:,:)
  real(RP), private, allocatable :: init_PT(:,:,:)
  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_USER / &
       USER_do


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
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none

    !---------------------------------------------------------------------------

    
    ! calculate diagnostic value and input to history buffer
    call USER_step


    return
  end subroutine USER_resume

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
       NOWSEC  => TIME_NOWSEC, &
       DTSEC => TIME_DTSEC
    use scale_history, only: &
       HIST_in
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

       call HIST_in( PT_diff(:,:,:), 'PT_diff', 'PT perturbation', 'K' )
       call HIST_in( U_diff(:,:,:), 'U_diff', 'U perturbation', 'm/s', xdim='half' )
    endif

    return
  end subroutine USER_step

end module mod_user
