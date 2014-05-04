!-------------------------------------------------------------------------------
!> module profiler
!!
!! @par Description
!!          Time counter & FLOP counter(PAPI) toolbox
!!
!! @author Team SCALE
!!
!<
module scale_prof
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PROF_rapstart
  public :: PROF_rapend
  public :: PROF_rapreport

#ifdef _PAPI_
  public :: PROF_PAPI_rapstart
  public :: PROF_PAPI_rapstop
  public :: PROF_PAPI_rapreport
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: PROF_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: PROF_rapnlimit = 100
  integer,                private,      save :: PROF_rapnmax   = 0
  character(len=H_SHORT), private,      save :: PROF_rapname(PROF_rapnlimit)
  real(DP),               private,      save :: PROF_raptstr(PROF_rapnlimit)
  real(DP),               private,      save :: PROF_rapttot(PROF_rapnlimit)
  integer,                private,      save :: PROF_rapnstr(PROF_rapnlimit)
  integer,                private,      save :: PROF_rapnend(PROF_rapnlimit)

#ifdef _PAPI_
  integer(DP),private :: PROF_PAPI_flops     = 0   !> total floating point operations since the first call
  real(SP),   private :: PROF_PAPI_real_time = 0.0 !> total realtime since the first PROF_PAPI_flops() call
  real(SP),   private :: PROF_PAPI_proc_time = 0.0 !> total process time since the first PROF_PAPI_flops() call
  real(SP),   private :: PROF_PAPI_mflops    = 0.0 !> Mflop/s achieved since the previous call
  integer,    private :: PROF_PAPI_check
#endif

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Get item ID or register item
  function PROF_rapid( rapname ) result(id)
    implicit none

    character(len=*), intent(in) :: rapname !< name of item

    integer :: id
    character (len=H_SHORT) :: trapname
    !---------------------------------------------------------------------------

    trapname = trim(rapname)

    do id = 1, PROF_rapnmax
       if( trapname == PROF_rapname(id) ) return
    enddo

    PROF_rapnmax     = PROF_rapnmax + 1
    id               = PROF_rapnmax
    PROF_rapname(id) = trapname

    PROF_rapnstr(id) = 0
    PROF_rapnend(id) = 0
    PROF_rapttot(id) = 0.0_DP

  end function PROF_rapid

  !-----------------------------------------------------------------------------
  !> Start raptime
  subroutine PROF_rapstart( rapname )
    use scale_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname !< name of item

    integer :: id
    !---------------------------------------------------------------------------

    id = PROF_rapid( rapname )

    PROF_raptstr(id) = PRC_MPItime()
    PROF_rapnstr(id) = PROF_rapnstr(id) + 1

    !if( IO_L ) write(IO_FID_LOG,*) rapname, PROF_rapnstr(id)

#ifdef _FAPP_
call START_COLLECTION( rapname )
#endif

    return
  end subroutine PROF_rapstart

  !-----------------------------------------------------------------------------
  !> Save raptime
  subroutine PROF_rapend( rapname )
    use scale_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname !< name of item

    integer :: id
    !---------------------------------------------------------------------------

    id = PROF_rapid( rapname )

    PROF_rapttot(id) = PROF_rapttot(id) + ( PRC_MPItime()-PROF_raptstr(id) )
    PROF_rapnend(id) = PROF_rapnend(id) + 1

#ifdef _FAPP_
call STOP_COLLECTION( rapname )
#endif

    return
  end subroutine PROF_rapend

  !-----------------------------------------------------------------------------
  !> Report raptime
  subroutine PROF_rapreport
    use scale_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_MPItimestat
    implicit none

    real(DP) :: avgvar(PROF_rapnlimit)
    real(DP) :: maxvar(PROF_rapnlimit)
    real(DP) :: minvar(PROF_rapnlimit)
    integer  :: maxidx(PROF_rapnlimit)
    integer  :: minidx(PROF_rapnlimit)

    integer :: id
    !---------------------------------------------------------------------------

    do id = 1, PROF_rapnmax
       if ( PROF_rapnstr(id) /= PROF_rapnend(id) ) then
           write(*,*) '*** Mismatch Report',id,PROF_rapname(id),PROF_rapnstr(id),PROF_rapnend(id)
       endif
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Computational Time Report'

    if ( IO_LOG_ALLNODE ) then ! report for each node

       do id = 1, PROF_rapnmax
          if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,I7)') &
                     '*** ID=',id,' : ',PROF_rapname(id),' T=',PROF_rapttot(id),' N=',PROF_rapnstr(id)
       enddo

    else

       call PRC_MPItimestat( avgvar      (1:PROF_rapnmax), &
                             maxvar      (1:PROF_rapnmax), &
                             minvar      (1:PROF_rapnmax), &
                             maxidx      (1:PROF_rapnmax), &
                             minidx      (1:PROF_rapnmax), &
                             PROF_rapttot(1:PROF_rapnmax)  )

       do id = 1, PROF_rapnmax
          if( IO_L ) write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                     '*** ID=',id,' : ',PROF_rapname(id), &
                     ' T(avg)=',avgvar(id), &
                     ', T(max)=',maxvar(id),'[',maxidx(id),']', &
                     ', T(min)=',minvar(id),'[',minidx(id),']', &
                     ' N=',PROF_rapnstr(id)
       enddo

       if ( IO_LOG_SUPPRESS ) then ! report to STDOUT
          if ( PRC_myrank == PRC_master ) then ! master node
             write(*,*) '*** Computational Time Report'
             do id = 1, PROF_rapnmax
                write(*,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                     '*** ID=',id,' : ',PROF_rapname(id), &
                     ' T(avg)=',avgvar(id), &
                     ', T(max)=',maxvar(id),'[',maxidx(id),']', &
                     ', T(min)=',minvar(id),'[',minidx(id),']', &
                     ' N=',PROF_rapnstr(id)
             enddo
          endif
       endif
    endif

    return
  end subroutine PROF_rapreport

#ifdef _PAPI_
  !-----------------------------------------------------------------------------
  !> Start flop counter
  subroutine PROF_PAPI_rapstart
    implicit none
    !---------------------------------------------------------------------------

    call PAPIF_flops( PROF_PAPI_real_time, PROF_PAPI_proc_time, PROF_PAPI_flops, PROF_PAPI_mflops, PROF_PAPI_check )

    return
  end subroutine PROF_PAPI_rapstart

  !-----------------------------------------------------------------------------
  !> Stop flop counter
  subroutine PROF_PAPI_rapstop
    implicit none
    !---------------------------------------------------------------------------

    call PAPIF_flops( PROF_PAPI_real_time, PROF_PAPI_proc_time, PROF_PAPI_flops, PROF_PAPI_mflops, PROF_PAPI_check )

    return
  end subroutine PROF_PAPI_rapstop

  !-----------------------------------------------------------------------------
  !> Report flop
  subroutine PROF_PAPI_rapreport
    use scale_process, only: &
       PRC_master, &
       PRC_myrank, &
       PRC_nmax,   &
       PRC_MPItimestat
    implicit none

    real(DP) :: avgvar(3)
    real(DP) :: maxvar(3)
    real(DP) :: minvar(3)
    integer  :: maxidx(3)
    integer  :: minidx(3)

    real(DP) :: PROF_PAPI_gflop
    real(DP) :: stats(3)
    !---------------------------------------------------------------------------

    PROF_PAPI_gflop = real(PROF_PAPI_flops,kind=8) / 1024.0_DP**3

    if ( IO_LOG_ALLNODE ) then ! report for each node

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** PAPI Report [Local PE information]'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** Real time          [sec] : ', PROF_PAPI_real_time
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** CPU  time          [sec] : ', PROF_PAPI_proc_time
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOP             [GFLOP] : ', PROF_PAPI_gflop
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOPS by PAPI   [GFLOPS] : ', PROF_PAPI_mflops/1024.0_DP
       if( IO_L ) write(IO_FID_LOG,'(1x,A,F15.3)') '*** FLOP / CPU Time [GFLOPS] : ', PROF_PAPI_gflop/PROF_PAPI_proc_time

    else
       stats(1) = real(PROF_PAPI_real_time,kind=8)
       stats(2) = real(PROF_PAPI_proc_time,kind=8)
       stats(3) = PROF_PAPI_gflop

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
  end subroutine PROF_PAPI_rapreport
#endif

end module scale_prof
!-------------------------------------------------------------------------------
