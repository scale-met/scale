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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

    call USER_step

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  subroutine USER_step
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid_real, only : &
       CZ => REAL_CZ
    use scale_time, only : &
       NOWSEC => TIME_NOWSEC
    use mod_atmos_vars, only: &
       QTRC
    implicit none

    integer :: modsec

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Add rain.'

       modsec = mod(int(NOWSEC),300)

       if ( modsec < 1 ) then
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if (       CZ(k,i,j) >= 2950.0_RP &
                  .AND. CZ(k,i,j) <  3050.0_RP ) then

                QTRC(k,i,j,I_QR) = 1.E-3
             endif
          enddo
          enddo
          enddo
       endif
    endif

    return
  end subroutine USER_step

end module mod_user
