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

  !-----------------------------------------------------------------------------
contains
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
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    call USER_step

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
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid_real, only : &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
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
       if( IO_L ) write(IO_FID_LOG,*) '*** Add rain.'
       USER_do = .false.

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          dist = ( ( CZ(k,i,j) - 3000.0_RP ) / 50.0_RP )**2

          QTRC(k,i,j,I_QR) = 1.E-3 / ( 1.0_RP + dist )
!          if (       FZ(k,  i,j) >= 3000.0_RP &
!               .AND. FZ(k-1,i,j) <  3000.0_RP ) then
!
!             if( IO_L ) write(IO_FID_LOG,*) k,i,j,CZ(k,i,j),dist,1.E-3/( 1.0_RP + dist )
!             QTRC(k,i,j,I_QR) = 1.E-3
!          endif
       enddo
       enddo
       enddo
    endif

    return
  end subroutine USER_step

end module mod_user
