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
  public :: PROF_setup
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
  private :: get_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                private, parameter :: PROF_rapnlimit = 100
  integer,                private,      save :: PROF_rapnmax   = 0
  character(len=H_SHORT), private,      save :: PROF_rapname(PROF_rapnlimit)
  integer,                private,      save :: PROF_grpnmax   = 0
  character(len=H_SHORT), private,      save :: PROF_grpname(PROF_rapnlimit)
  integer,                private,      save :: PROF_grpid  (PROF_rapnlimit)
  real(DP),               private,      save :: PROF_raptstr(PROF_rapnlimit)
  real(DP),               private,      save :: PROF_rapttot(PROF_rapnlimit)
  integer,                private,      save :: PROF_rapnstr(PROF_rapnlimit)
  integer,                private,      save :: PROF_rapnend(PROF_rapnlimit)
  integer,                private,      save :: PROF_raplevel(PROF_rapnlimit)

  integer,                private, parameter :: PROF_default_rap_level = 2
  integer,                private,      save :: PROF_rap_level = 2

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
  subroutine PROF_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_PROF / &
         PROF_rap_level

    integer :: ierr

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[PROF] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PROF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_PROF. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_PROF)

    return
  end subroutine PROF_setup

  !-----------------------------------------------------------------------------
  !> Start raptime
  subroutine PROF_rapstart( rapname, level )
    use scale_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname !< name of item
    integer, intent(in), optional :: level  !< level of item (default is 2)

    integer :: id
    integer :: level_
    !---------------------------------------------------------------------------

    if ( present(level) ) then
       level_ = level
    else
       level_ = PROF_default_rap_level
    end if

    if ( level_ > PROF_rap_level ) return

    id = get_rapid( rapname, level_ )

    PROF_raptstr(id) = PRC_MPItime()
    PROF_rapnstr(id) = PROF_rapnstr(id) + 1

    !if( IO_L ) write(IO_FID_LOG,*) rapname, PROF_rapnstr(id)

#ifdef _FAPP_
    call FAPP_START( trim(PROF_grpname(get_grpid(rapname))), id, level_ )
#endif
#ifdef _FINEPA_
    call START_COLLECTION( rapname )
#endif

    return
  end subroutine PROF_rapstart

  !-----------------------------------------------------------------------------
  !> Save raptime
  subroutine PROF_rapend( rapname, level )
    use scale_process, only: &
       PRC_MPItime
    implicit none

    character(len=*), intent(in) :: rapname !< name of item
    integer, intent(in), optional :: level  !< level of item

    integer :: id
    integer :: level_
    !---------------------------------------------------------------------------

    if ( present(level) ) then
       if ( level > PROF_rap_level ) return
    end if

    id = get_rapid( rapname, level_ )

    if ( level_ > PROF_rap_level ) return

    PROF_rapttot(id) = PROF_rapttot(id) + ( PRC_MPItime()-PROF_raptstr(id) )
    PROF_rapnend(id) = PROF_rapnend(id) + 1

#ifdef _FINEPA_
    call STOP_COLLECTION( rapname )
#endif
#ifdef _FAPP_
    call FAPP_STOP( trim(PROF_grpname(PROF_grpid(id))), id, level_ )
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

    integer :: id, gid
    integer :: fid
    !---------------------------------------------------------------------------

    do id = 1, PROF_rapnmax
       if ( PROF_rapnstr(id) /= PROF_rapnend(id) ) then
           write(*,*) '*** Mismatch Report',id,PROF_rapname(id),PROF_rapnstr(id),PROF_rapnend(id)
       endif
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Computational Time Report'
    if( IO_L ) write(IO_FID_LOG,*) '*** Rap level is ', PROF_rap_level

    if ( IO_LOG_ALLNODE ) then ! report for each node

       do gid = 1, PROF_rapnmax
          do id = 1, PROF_rapnmax
             if ( PROF_raplevel(id) .le. PROF_rap_level .and. &
                  PROF_grpid(id) == gid .and. &
                  IO_L ) then
                write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,I7)') &
                     '*** ID=',id,' : ',PROF_rapname(id),' T=',PROF_rapttot(id),' N=',PROF_rapnstr(id)
             end if
          enddo
       enddo

    else

       call PRC_MPItimestat( avgvar      (1:PROF_rapnmax), &
                             maxvar      (1:PROF_rapnmax), &
                             minvar      (1:PROF_rapnmax), &
                             maxidx      (1:PROF_rapnmax), &
                             minidx      (1:PROF_rapnmax), &
                             PROF_rapttot(1:PROF_rapnmax)  )

       fid = -1
       if ( IO_LOG_SUPPRESS ) then ! report to STDOUT
          if ( PRC_myrank == PRC_master ) then
             write(*,*) '*** Computational Time Report'
             fid = 6 ! master node
          end if
       else
          if ( IO_L ) fid = IO_FID_LOG
       end if

       do gid = 1, PROF_rapnmax
          do id = 1, PROF_rapnmax
             if ( PROF_raplevel(id) .le. PROF_rap_level .and. &
                  PROF_grpid(id) == gid .and. &
                  fid > 0 ) then
                write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I7)') &
                     '*** ID=',id,' : ',PROF_rapname(id), &
                     ' T(avg)=',avgvar(id), &
                     ', T(max)=',maxvar(id),'[',maxidx(id),']', &
                     ', T(min)=',minvar(id),'[',minidx(id),']', &
                     ' N=',PROF_rapnstr(id)
             end if
          enddo
       enddo

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
    real(DP) :: statistics(3)
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
       statistics(1) = real(PROF_PAPI_real_time,kind=8)
       statistics(2) = real(PROF_PAPI_proc_time,kind=8)
       statistics(3) = PROF_PAPI_gflop

       call PRC_MPItimestat( avgvar(1:3), &
                             maxvar(1:3), &
                             minvar(1:3), &
                             maxidx(1:3), &
                             minidx(1:3), &
                             statistics (1:3)  )

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

  !-----------------------------------------------------------------------------
  !> Get item ID or register item
  function get_rapid( rapname, level ) result(id)
    implicit none

    character(len=*), intent(in)    :: rapname !< name of item
    integer,          intent(inout) :: level   !< level of item

    integer :: id
    character (len=H_SHORT) :: trapname
    !---------------------------------------------------------------------------

    trapname = trim(rapname)

    do id = 1, PROF_rapnmax
       if( trapname == PROF_rapname(id) )  then
          level = PROF_raplevel(id)
          return
       end if
    enddo

    PROF_rapnmax     = PROF_rapnmax + 1
    id               = PROF_rapnmax
    PROF_rapname(id) = trapname
    
    PROF_grpid(id) = get_grpid(trapname)

    PROF_rapnstr(id) = 0
    PROF_rapnend(id) = 0
    PROF_rapttot(id) = 0.0_DP

    PROF_raplevel(id) = level

    return
  end function get_rapid

  !-----------------------------------------------------------------------------
  !> Get group ID
  function get_grpid( rapname ) result(gid)
    implicit none

    character(len=*), intent(in)    :: rapname !< name of item

    integer :: gid
    integer :: idx
    character(len=H_SHORT) :: grpname
    !---------------------------------------------------------------------------

    idx = index(rapname," ")
    if ( idx > 1 ) then
       grpname = rapname(1:idx-1)
    else
       grpname = rapname
    end if

    do gid = 1, PROF_grpnmax
       if( grpname == PROF_grpname(gid) ) return
    enddo

    PROF_grpnmax      = PROF_grpnmax + 1
    gid               = PROF_grpnmax
    PROF_grpname(gid) = grpname

    return
  end function get_grpid

end module scale_prof
!-------------------------------------------------------------------------------
