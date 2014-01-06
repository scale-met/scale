!-------------------------------------------------------------------------------
!> module PAPI
!!
!! @par Description
!!          PAPI toolbox
!!
!! @author Team SCALE
!!
!<
module mod_papi
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use dc_types, only: &
     DP, SP
  use mod_precision
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_SYSCHR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PAPI_rapstart
  public :: PAPI_rapstop
  public :: PAPI_rapreport

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
  integer(DP),private :: papi_flpops    = 0   !> total floating point operations since the first call
  real(SP),   private :: papi_real_time = 0.0 !> total realtime since the first PAPI_flops() call
  real(SP),   private :: papi_proc_time = 0.0 !> total process time since the first PAPI_flops() call
  real(SP),   private :: papi_mflops    = 0.0 !> Mflop/s achieved since the previous call
  integer,    private :: papi_check

  integer,                  private, parameter :: PAPI_rapnlimit = 100
  integer,                  private            :: PAPI_rapnmax   = 0
  character(len=IO_SYSCHR), private            :: PAPI_rapname(PAPI_rapnlimit)
  real(DP),                 private            :: PAPI_raptstr(PAPI_rapnlimit)
  real(DP),                 private            :: PAPI_rapttot(PAPI_rapnlimit)
  integer,                  private            :: PAPI_rapnstr(PAPI_rapnlimit)
  integer,                  private            :: PAPI_rapnend(PAPI_rapnlimit)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Start flop counter
  subroutine PAPI_rapstart
    implicit none
    !---------------------------------------------------------------------------

    call PAPIF_flops( PAPI_real_time, PAPI_proc_time, PAPI_flpops, PAPI_mflops, PAPI_check )

    return
  end subroutine PAPI_rapstart

  !-----------------------------------------------------------------------------
  !> Stop flop counter
  subroutine PAPI_rapstop
    implicit none
    !---------------------------------------------------------------------------

    call PAPIF_flops( PAPI_real_time, PAPI_proc_time, PAPI_flpops, PAPI_mflops, PAPI_check )

    return
  end subroutine PAPI_rapstop

  !-----------------------------------------------------------------------------
  !> Report flop
  subroutine PAPI_rapreport
    use mod_stdio, only: &
       IO_LOG_SUPPRESS, &
       IO_LOG_ALLNODE
    use mod_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_nmax,   &
       PRC_MPItimestat
    implicit none

    real(DP) :: avgvar(PAPI_rapnlimit)
    real(DP) :: maxvar(PAPI_rapnlimit)
    real(DP) :: minvar(PAPI_rapnlimit)
    integer  :: maxidx(PAPI_rapnlimit)
    integer  :: minidx(PAPI_rapnlimit)

    real(DP) :: papi_gflop
    real(DP) :: stats(PAPI_rapnlimit)
    !---------------------------------------------------------------------------

    papi_gflop = real(papi_flpops,kind=8) / 1024.0_DP**3

    if ( IO_LOG_ALLNODE ) then ! report for each node

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** PAPI Report [Local PE information]'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** Real time          [sec] : ', papi_real_time
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** CPU  time          [sec] : ', papi_proc_time
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOP             [GFLOP] : ', papi_gflop
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOPS by PAPI    [GFLOPS] : ', papi_mflops/1024.0_DP
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOP / CPU Time [GFLOPS] : ', papi_gflop/papi_proc_time

    else
       stats(1) = real(papi_real_time,kind=8)
       stats(2) = real(papi_proc_time,kind=8)
       stats(3) = papi_gflop

       call PRC_MPItimestat( avgvar(1:3), &
                             maxvar(1:3), &
                             minvar(1:3), &
                             maxidx(1:3), &
                             minidx(1:3), &
                             stats (1:3)  )

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** PAPI Report'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** Real time [sec]',' T(avg)=',avgvar(1), &
                  ', T(max)=',maxvar(1),'[',maxidx(1),']',', T(min)=',minvar(1),'[',minidx(1),']'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** CPU  time [sec]',' T(avg)=',avgvar(2), &
                  ', T(max)=',maxvar(2),'[',maxidx(2),']',', T(min)=',minvar(2),'[',minidx(2),']'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** FLOP    [GFLOP]',' N(avg)=',avgvar(3), &
                  ', N(max)=',maxvar(3),'[',maxidx(3),']',', N(min)=',minvar(3),'[',minidx(3),']'
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3,A,I6,A)') &
                  '*** TOTAL FLOP    [GFLOP] : ', avgvar(3)*PRC_nmax, '(',PRC_nmax,' PEs)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') &
                  '*** FLOPS        [GFLOPS] : ', avgvar(3)*PRC_nmax/maxvar(2)
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') &
                  '*** FLOPS per PE [GFLOPS] : ', avgvar(3)/maxvar(2)
       if( IO_L ) write(IO_FID_LOG,*)

       if ( IO_LOG_SUPPRESS ) then ! report to STDOUT
          if ( PRC_myrank == PRC_master ) then ! master node
             write(*,*)
             write(*,*) '*** PAPI Report'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** Real time [sec]',' T(avg)=',avgvar(1), &
                  ', T(max)=',maxvar(1),'[',maxidx(1),']',', T(min)=',minvar(1),'[',minidx(1),']'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** CPU  time [sec]',' T(avg)=',avgvar(2), &
                  ', T(max)=',maxvar(2),'[',maxidx(2),']',', T(min)=',minvar(2),'[',minidx(2),']'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                  '*** FLOP    [GFLOP]',' N(avg)=',avgvar(3), &
                  ', N(max)=',maxvar(3),'[',maxidx(3),']',', N(min)=',minvar(3),'[',minidx(3),']'
             write(*,*)
             write(*,'(1x,A,F15.3,A,I6,A)') &
                  '*** TOTAL FLOP    [GFLOP] : ', avgvar(3)*PRC_nmax, '(',PRC_nmax,' PEs)'
             write(*,'(1x,A,F15.3)') &
                  '*** FLOPS        [GFLOPS] : ', avgvar(3)*PRC_nmax/maxvar(2)
             write(*,'(1x,A,F15.3)') &
                  '*** FLOPS per PE [GFLOPS] : ', avgvar(3)/maxvar(2)
          endif
       endif
    endif

    return
  end subroutine PAPI_rapreport

end module mod_papi
!-------------------------------------------------------------------------------
