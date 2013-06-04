!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          User defined module
!!
!! @author H.Tomita and SCALE developers
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
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
  !++ included parameters
  !
  include "inc_precision.h"
  include "inc_index.h"
  include 'inc_tracer.h'

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
  logical, private, save :: USER_do = .false.

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  !-----------------------------------------------------------------------------
  subroutine USER_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
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
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Step
  !-----------------------------------------------------------------------------
  subroutine USER_step
    use mod_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use mod_topography, only: &
       TOPO_CZ, &
       TOPO_FZ, &
       TOPO_CGSQRT, &
       TOPO_FGSQRT, &
       TOPO_Gvec
    use mod_time, only: &
       TIME_DTSEC
    use mod_history, only: &
       HIST_in
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: TEST(IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call HIST_in( TOPO_CZ    (:,:,:), 'TOPO_CZ',     'TOPO_CZ',     '', TIME_DTSEC )
    call HIST_in( TOPO_FZ    (:,:,:), 'TOPO_FZ',     'TOPO_FZ',     '', TIME_DTSEC )
    call HIST_in( TOPO_CGSQRT(:,:,:), 'TOPO_CGSQRT', 'TOPO_CGSQRT', '', TIME_DTSEC )
    call HIST_in( TOPO_FGSQRT(:,:,:), 'TOPO_FGSQRT', 'TOPO_FGSQRT', '', TIME_DTSEC )
    call HIST_in( TOPO_Gvec(:,:,:,1), 'TOPO_GX',     'TOPO_GX',     '', TIME_DTSEC )
    call HIST_in( TOPO_Gvec(:,:,:,2), 'TOPO_GY',     'TOPO_GY',     '', TIME_DTSEC )

!    do j = 1, JA
!    do i = 1, IA
!       TEST(i,j) = PRC_myrank * 10000.0_RP &
!                 + j * 100.0_RP            &
!                 + i
!    enddo
!    enddo
!
!    call COMM_vars8( TEST(:,:), 1 )
!    call COMM_wait ( TEST(:,:), 1 )
!
!    do j = 1, JA
!       if( IO_L ) write(IO_FID_LOG,'(44I5)') ( int(TEST(i,j)),i=1,IA )
!    enddo
!
!
!    call HIST_in( TEST(:,:), 'TEST', 'TEST', '', TIME_DTSEC )

    if ( USER_do ) then
       call PRC_MPIstop
    endif

    return
  end subroutine USER_step

end module mod_user
