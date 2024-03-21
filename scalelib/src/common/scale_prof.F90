!-------------------------------------------------------------------------------
!> module profiler
!!
!! @par Description
!!          Time counter & FLOP counter(PAPI) toolbox
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_prof
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: PROF_setup
  public :: PROF_finalize
  public :: PROF_setprefx
  public :: PROF_rapstart
  public :: PROF_rapend
  public :: PROF_rapreport

#ifdef PAPI
  public :: PROF_PAPI_rapstart
  public :: PROF_PAPI_rapstop
  public :: PROF_PAPI_rapreport
#endif

  public :: PROF_valcheck

  interface PROF_valcheck
     module procedure PROF_valcheck_SP_1D
     module procedure PROF_valcheck_SP_2D
     module procedure PROF_valcheck_SP_3D
     module procedure PROF_valcheck_SP_4D
     module procedure PROF_valcheck_SP_5D
     module procedure PROF_valcheck_SP_6D
     module procedure PROF_valcheck_DP_1D
     module procedure PROF_valcheck_DP_2D
     module procedure PROF_valcheck_DP_3D
     module procedure PROF_valcheck_DP_4D
     module procedure PROF_valcheck_DP_5D
     module procedure PROF_valcheck_DP_6D
  end interface PROF_valcheck

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
  integer,                private, parameter :: PROF_rapnlimit = 300
  character(len=H_SHORT), private            :: PROF_prefix    = ''
  integer,                private            :: PROF_rapnmax   = 0
  character(len=H_SHORT), private            :: PROF_rapname(PROF_rapnlimit)
  integer,                private            :: PROF_grpnmax   = 0
  character(len=H_SHORT), private            :: PROF_grpname(PROF_rapnlimit)
  integer,                private            :: PROF_grpid  (PROF_rapnlimit)
  real(DP),               private            :: PROF_raptstr(PROF_rapnlimit)
  real(DP),               private            :: PROF_rapttot(PROF_rapnlimit)
  integer,                private            :: PROF_rapnstr(PROF_rapnlimit)
  integer,                private            :: PROF_rapcnt (PROF_rapnlimit)
  integer,                private            :: PROF_rapnend(PROF_rapnlimit)
  integer,                private            :: PROF_raplevel(PROF_rapnlimit)

  integer,                private, parameter :: PROF_default_rap_level = 2
  integer,                private            :: PROF_rap_level         = 2
  logical,                private            :: PROF_mpi_barrier       = .false.

#ifdef PAPI
  integer(DP),            private            :: PROF_PAPI_flops     = 0   !> total floating point operations since the first call
  real(SP),               private            :: PROF_PAPI_real_time = 0.0 !> total realtime since the first PROF_PAPI_flops() call
  real(SP),               private            :: PROF_PAPI_proc_time = 0.0 !> total process time since the first PROF_PAPI_flops() call
  real(SP),               private            :: PROF_PAPI_mflops    = 0.0 !> Mflop/s achieved since the previous call
  integer,                private            :: PROF_PAPI_check
#endif

  character(len=7),       private            :: PROF_header
  character(len=16),      private            :: PROF_item
  real(DP),               private            :: PROF_max
  real(DP),               private            :: PROF_min
  real(DP),               private            :: PROF_sum

  logical,                private            :: PROF_barrier_flag

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine PROF_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_PROF / &
         PROF_rap_level, &
         PROF_mpi_barrier

    integer :: ierr

    LOG_NEWLINE
    LOG_INFO("PROF_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PROF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("PROF_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("PROF_setup",*) 'Not appropriate names in namelist PARAM_PROF. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_PROF)

    LOG_NEWLINE
    LOG_INFO("PROF_setup",*) 'Rap output level              = ', PROF_rap_level
    LOG_INFO("PROF_setup",*) 'Add MPI_barrier in every rap? = ', PROF_mpi_barrier

    PROF_prefix = ''

    return
  end subroutine PROF_setup

  subroutine PROF_finalize

    PROF_rap_level   = 2
    PROF_mpi_barrier = .false.
#ifdef PAPI
    PROF_PAPI_flops     = 0
    PROF_PAPI_real_time = 0.0
    PROF_PAPI_proc_time = 0.0
    PROF_PAPI_mflops    = 0.0
#endif

    PROF_rapnmax = 0
    PROF_grpnmax = 0

    return
  end subroutine PROF_finalize

  !-----------------------------------------------------------------------------
  subroutine PROF_setprefx( &
       prefxname )
    implicit none

    character(len=*), intent(in) :: prefxname !< prefix

    !---------------------------------------------------------------------------

    PROF_prefix = prefxname

    return
  end subroutine PROF_setprefx

  !-----------------------------------------------------------------------------
  !> Start raptime
  subroutine PROF_rapstart( rapname_base, level, disable_barrier )
    use scale_prc, only: &
       PRC_MPIbarrier, &
       PRC_MPItime
    implicit none

    character(len=*), intent(in)           :: rapname_base    !< name  of item
    integer,          intent(in), optional :: level           !< level of item
    logical,          intent(in), optional :: disable_barrier !< disable barrier if .true.

    character(len=H_SHORT) :: rapname !< name of item with prefix

    integer :: id
    integer :: level_
    integer :: tn
    integer :: i
    logical :: disable_barrier_
    !$ integer :: omp_get_thread_num
    !---------------------------------------------------------------------------

    tn = 0
    !$ tn = omp_get_thread_num()
    if ( tn > 0 ) return

    if ( present(level) ) then
       level_ = level
    else
       level_ = PROF_default_rap_level
    endif

    if ( present(disable_barrier) ) then
       disable_barrier_ = disable_barrier
    else
       disable_barrier_ = .false.
    endif

    if( level_ > PROF_rap_level ) return

    if ( len_trim(PROF_prefix) > 0 ) then
       rapname = trim(PROF_prefix)//" "//trim(rapname_base)
    else
       rapname = rapname_base
    end if

    id = get_rapid( rapname, level_ )

    PROF_rapcnt(id) = PROF_rapcnt(id) + 1

    if ( PROF_rapcnt(id) > 1 ) return

    if ( ( .not. disable_barrier_ ) .and. PROF_mpi_barrier ) call PRC_MPIbarrier

    PROF_raptstr(id) = PRC_MPItime()
    PROF_rapnstr(id) = PROF_rapnstr(id) + 1

    !LOG_INFO("PROF_rapstart",'(1x,A,I8)') rapname, PROF_rapnstr(id)
    !call flush(IO_FID_LOG)

#ifdef FAPP
    i = index(rapname," ")
    if ( i == 0 .or. i > len_trim(rapname)) then
       call FAPP_START( rapname, id, level_ )
    else
       call FAPP_START( rapname(1:i-1)//"_"//trim(rapname(i+1:)), id, level_ )
    end if
#endif

    return
  end subroutine PROF_rapstart

  !-----------------------------------------------------------------------------
  !> Save raptime
  subroutine PROF_rapend( rapname_base, level, disable_barrier )
    use scale_prc, only: &
       PRC_MPIbarrier, &
       PRC_MPItime
    implicit none

    character(len=*), intent(in)           :: rapname_base    !< name  of item
    integer,          intent(in), optional :: level           !< level of item
    logical,          intent(in), optional :: disable_barrier !< disable barrier if .true.

    character(len=H_SHORT) :: rapname !< name of item with prefix

    integer :: id
    integer :: level_
    integer :: tn
    integer :: i
    logical :: disable_barrier_
    !$ integer :: omp_get_thread_num
    !---------------------------------------------------------------------------

    tn = 0
    !$ tn = omp_get_thread_num()
    if ( tn > 0 ) return

    if ( present(level) ) then
       if( level > PROF_rap_level ) return
    endif

    if ( len_trim(PROF_prefix) > 0 ) then
       rapname = trim(PROF_prefix)//" "//trim(rapname_base)
    else
       rapname = rapname_base
    end if

    if ( present(disable_barrier) ) then
       disable_barrier_ = disable_barrier
    else
       disable_barrier_ = .false.
    endif


    level_ = -1
    id = get_rapid( rapname, level_ )

    if( level_ > PROF_rap_level ) return

    PROF_rapcnt(id) = PROF_rapcnt(id) - 1

    if ( PROF_rapcnt(id) > 0 ) return

#ifdef FAPP
    i = index(rapname," ")
    if ( i == 0 .or. i > len_trim(rapname)) then
       call FAPP_STOP( rapname, id, level_ )
    else
       call FAPP_STOP( rapname(1:i-1)//"_"//trim(rapname(i+1:)), id, level_ )
    end if
#endif

    PROF_rapttot(id) = PROF_rapttot(id) + ( PRC_MPItime()-PROF_raptstr(id) )
    PROF_rapnend(id) = PROF_rapnend(id) + 1

    if ( ( .not. disable_barrier_ ) .and. PROF_mpi_barrier ) call PRC_MPIbarrier

    return
  end subroutine PROF_rapend

  !-----------------------------------------------------------------------------
  !> Report raptime
  subroutine PROF_rapreport
    use scale_prc, only: &
       PRC_MPItimestat, &
       PRC_IsMaster
    implicit none

    real(DP) :: avgvar(PROF_rapnlimit)
    real(DP) :: maxvar(PROF_rapnlimit)
    real(DP) :: minvar(PROF_rapnlimit)
    integer  :: maxidx(PROF_rapnlimit)
    integer  :: minidx(PROF_rapnlimit)

    integer  :: idx(PROF_rapnlimit)
    integer  :: cnt

    integer :: id, gid
    integer :: fid
    integer :: i, j
    !---------------------------------------------------------------------------

    do id = 1, PROF_rapnmax
       if ( PROF_rapnstr(id) /= PROF_rapnend(id) ) then
           LOG_WARN("PROF_rapreport",*) 'Mismatch Report',id,PROF_rapname(id),PROF_rapnstr(id),PROF_rapnend(id)
       endif
    enddo

    LOG_NEWLINE
    LOG_INFO("PROF_rapreport",'(1x,A,I2,A)') 'Computational Time Report (Rap level = ', PROF_rap_level, ')'

    if ( IO_LOG_ALLNODE ) then ! report for each node

       do gid = 1, PROF_rapnmax
          cnt = 0
          do id = 1, PROF_rapnmax
             if (       PROF_raplevel(id) <= PROF_rap_level &
                  .AND. PROF_grpid   (id) == gid            ) then
                cnt = cnt + 1
                idx(cnt) = id
             end if
          end do
          do j = 1, cnt-1
             do i = j+1, cnt
                if ( PROF_raplevel(idx(i)) < PROF_raplevel(idx(j)) .or. &
                     ( PROF_raplevel(idx(i)) == PROF_raplevel(idx(j)) &
                     .and. PROF_rapname(idx(i)) < PROF_rapname(idx(j)) ) ) then
                   id = idx(i)
                   idx(i) = idx(j)
                   idx(j) = id
                end if
             end do
          end do
          do i = 1, cnt
             id = idx(i)
             LOG_INFO_CONT('(1x,2A,I2,A,F10.3,A,I9)') &
                  PROF_rapname(id), ' lev=', PROF_raplevel(id), &
                  ': T=',PROF_rapttot(id),' N=',PROF_rapnstr(id)
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
          if ( PRC_IsMaster ) then
             write(*,*) 'INFO  [PROF_rapreport] Computational Time Report'
             fid = IO_FID_STDOUT ! master node
          endif
       else
          if ( IO_L ) fid = IO_FID_LOG
       endif

       if ( fid > 0 ) then
          do gid = 1, PROF_rapnmax
             cnt = 0
             do id = 1, PROF_rapnmax
                if (       PROF_raplevel(id) <= PROF_rap_level &
                     .AND. PROF_grpid   (id) == gid            ) then
                   cnt = cnt + 1
                   idx(cnt) = id
                end if
             end do
             do j = 1, cnt-1
                do i = j+1, cnt
                   if ( PROF_raplevel(idx(i)) < PROF_raplevel(idx(j)) .or. &
                        ( PROF_raplevel(idx(i)) == PROF_raplevel(idx(j)) &
                        .and. PROF_rapname(idx(i)) < PROF_rapname(idx(j)) ) ) then
                      id = idx(i)
                      idx(i) = idx(j)
                      idx(j) = id
                   end if
                end do
             end do
             do i = 1, cnt
                id = idx(i)
                write(fid,'(1x,2A,I2,A,F10.3,2(A,F10.3,A,I6,A),A,I9)') &
                     PROF_rapname(id), ' lev=', PROF_raplevel(id), &
                     ': T(avg)=',avgvar(id), &
                     ', T(max)=',maxvar(id),'[',maxidx(id),']', &
                     ', T(min)=',minvar(id),'[',minidx(id),']', &
                     ', N=',PROF_rapnstr(id)
             end do
          enddo
       end if

    endif

    return
  end subroutine PROF_rapreport

#ifdef PAPI
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
    use scale_prc, only: &
       PRC_MPItimestat, &
       PRC_nprocs,      &
       PRC_IsMaster
    implicit none

    real(DP) :: avgvar(3)
    real(DP) :: maxvar(3)
    real(DP) :: minvar(3)
    integer  :: maxidx(3)
    integer  :: minidx(3)

    real(DP) :: zerosw
    real(DP) :: PROF_PAPI_gflop
    real(DP) :: statistics(3)
    !---------------------------------------------------------------------------

    PROF_PAPI_gflop = real(PROF_PAPI_flops,kind=8) / 1024.0_DP**3

    if ( IO_LOG_ALLNODE ) then ! report for each node

       LOG_NEWLINE
       LOG_INFO("PROF_PAPI_rapreport",*) 'PAPI Report [Local PE information]'
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3)') 'Real time          [sec] : ', PROF_PAPI_real_time
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3)') 'CPU  time          [sec] : ', PROF_PAPI_proc_time
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3)') 'FLOP             [GFLOP] : ', PROF_PAPI_gflop
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3)') 'FLOPS by PAPI   [GFLOPS] : ', PROF_PAPI_mflops/1024.0_DP
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3)') 'FLOP / CPU Time [GFLOPS] : ', PROF_PAPI_gflop/PROF_PAPI_proc_time

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

       zerosw = 0.5_DP - sign(0.5_DP,maxvar(2)-1.D-12) ! if maxvar(2) = 0 then zerosw = 1

       LOG_NEWLINE
       LOG_INFO("PROF_PAPI_rapreport",*) 'PAPI Report'
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,A,F10.3,A,F10.3,A,I6,A,A,F10.3,A,I6,A)') &
                  'Real time [sec]',' T(avg)=',avgvar(1), &
                  ', T(max)=',maxvar(1),'[',maxidx(1),']',', T(min)=',minvar(1),'[',minidx(1),']'
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,A,F10.3,A,F10.3,A,I6,A,A,F10.3,A,I6,A)') &
                  'CPU  time [sec]',' T(avg)=',avgvar(2), &
                  ', T(max)=',maxvar(2),'[',maxidx(2),']',', T(min)=',minvar(2),'[',minidx(2),']'
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,A,F10.3,A,F10.3,A,I6,A,A,F10.3,A,I6,A)') &
                  'FLOP    [GFLOP]',' N(avg)=',avgvar(3), &
                  ', N(max)=',maxvar(3),'[',maxidx(3),']',', N(min)=',minvar(3),'[',minidx(3),']'
       LOG_NEWLINE
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3,A,I6,A)') &
                  'TOTAL FLOP    [GFLOP] : ', avgvar(3)*PRC_nprocs, '(',PRC_nprocs,' PEs)'
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3)') &
                  'FLOPS        [GFLOPS] : ', avgvar(3)*PRC_nprocs * ( 1.0_DP-zerosw ) / ( maxvar(2)+zerosw )
       LOG_INFO("PROF_PAPI_rapreport",'(1x,A,F15.3)') &
                  'FLOPS per PE [GFLOPS] : ', avgvar(3)            * ( 1.0_DP-zerosw ) / ( maxvar(2)+zerosw )
       LOG_NEWLINE

       if ( IO_LOG_SUPPRESS ) then ! report to STDOUT
          if ( PRC_IsMaster ) then ! master node
             write(*,*)
             write(*,*) '*** PAPI Report'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I6,A,A,F10.3,A,I6,A)') &
                  '*** Real time [sec]',' T(avg)=',avgvar(1), &
                  ', T(max)=',maxvar(1),'[',maxidx(1),']',', T(min)=',minvar(1),'[',minidx(1),']'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I6,A,A,F10.3,A,I6,A)') &
                  '*** CPU  time [sec]',' T(avg)=',avgvar(2), &
                  ', T(max)=',maxvar(2),'[',maxidx(2),']',', T(min)=',minvar(2),'[',minidx(2),']'
             write(*,'(1x,A,A,F10.3,A,F10.3,A,I6,A,A,F10.3,A,I6,A)') &
                  '*** FLOP    [GFLOP]',' N(avg)=',avgvar(3), &
                  ', N(max)=',maxvar(3),'[',maxidx(3),']',', N(min)=',minvar(3),'[',minidx(3),']'
             write(*,*)
             write(*,'(1x,A,F15.3,A,I6,A)') &
                  '*** TOTAL FLOP    [GFLOP] : ', avgvar(3)*PRC_nprocs, '(',PRC_nprocs,' PEs)'
             write(*,'(1x,A,F15.3)') &
                  '*** FLOPS        [GFLOPS] : ', avgvar(3)*PRC_nprocs * ( 1.0_DP-zerosw ) / ( maxvar(2)+zerosw )
             write(*,'(1x,A,F15.3)') &
                  '*** FLOPS per PE [GFLOPS] : ', avgvar(3)            * ( 1.0_DP-zerosw ) / ( maxvar(2)+zerosw )
          endif
       endif
    endif

    return
  end subroutine PROF_PAPI_rapreport
#endif

  !-----------------------------------------------------------------------------
  !> Get item ID or register item
  function get_rapid( rapname, level ) result(id)
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in)    :: rapname !< name of item
    integer,          intent(inout) :: level   !< level of item

    character (len=H_SHORT) :: trapname

    integer :: lev
    integer :: id
    !---------------------------------------------------------------------------

    trapname  = trim(rapname)

    do id = 1, PROF_rapnmax
       if ( trapname == PROF_rapname(id) ) then
          lev = PROF_raplevel(id)
#ifdef QUICKDEBUG
          if ( level > 0 .and. lev .ne. level ) then
             LOG_ERROR("PROF_get_rapid",*) 'level is different ', trim(rapname), lev, level
             call PRC_abort
          end if
#endif
          level = lev
          return
       endif
    enddo

    PROF_rapnmax     = PROF_rapnmax + 1
    id               = PROF_rapnmax
    PROF_rapname(id) = trapname

    PROF_rapnstr(id) = 0
    PROF_rapnend(id) = 0
    PROF_rapttot(id) = 0.0_DP
    PROF_rapcnt (id) = 0

    PROF_grpid   (id) = get_grpid(trapname)
    PROF_raplevel(id) = level

    return
  end function get_rapid

  !-----------------------------------------------------------------------------
  !> Get group ID
  function get_grpid( rapname ) result(gid)
    implicit none

    character(len=*), intent(in)    :: rapname !< name of item

    character(len=H_SHORT) :: grpname

    integer :: gid
    integer :: idx
    !---------------------------------------------------------------------------

    idx = index(rapname," ")
    if ( idx > 1 ) then
       grpname = rapname(1:idx-1)
    else
       grpname = rapname
    endif

    do gid = 1, PROF_grpnmax
       if( grpname == PROF_grpname(gid) ) return
    enddo

    PROF_grpnmax      = PROF_grpnmax + 1
    gid               = PROF_grpnmax
    PROF_grpname(gid) = grpname

    return
  end function get_grpid

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_1D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_SP_1D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_1D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_2D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_SP_2D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_2D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_3D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_SP_3D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_3D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_4D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_SP_4D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_4D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_5D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_SP_5D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_5D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_6D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(SP),         intent(in) :: var(:,:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_SP_6D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_6D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_1D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_DP_1D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_1D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_2D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_DP_2D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_2D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_3D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_DP_3D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_3D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_4D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_DP_4D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_4D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_5D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_DP_5D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_5D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_6D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*), intent(in) :: header
    character(len=*), intent(in) :: varname
    real(DP),         intent(in) :: var(:,:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    LOG_INFO("PROF_valcheck_DP_6D",'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_6D

end module scale_prof
