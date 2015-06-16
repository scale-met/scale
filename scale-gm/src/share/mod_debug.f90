!-------------------------------------------------------------------------------
!> Module debug
!!
!! @par Description
!!         Debug utility module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_debug
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use mod_precision
  use mod_adm, only: &
     ADM_LOG_FID, &
     ADM_NSYS,    &
     ADM_MAXFNAME
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: DEBUG_dampdata
  public :: DEBUG_dampascii4D
  public :: DEBUG_dampascii3D
  public :: DEBUG_rapstart
  public :: DEBUG_rapend
  public :: DEBUG_rapreport

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: DEBUG_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                 private, parameter :: DEBUG_rapnlimit = 100
  integer,                 private            :: DEBUG_rapnmax   = 0
  character(len=ADM_NSYS), private            :: DEBUG_rapname(DEBUG_rapnlimit)
  real(RP),                 private            :: DEBUG_raptstr(DEBUG_rapnlimit)
  real(RP),                 private            :: DEBUG_rapttot(DEBUG_rapnlimit)
  integer,                 private            :: DEBUG_rapnstr(DEBUG_rapnlimit)
  integer,                 private            :: DEBUG_rapnend(DEBUG_rapnlimit)

#ifdef PAPI_OPS
  integer(8),public :: papi_flpops      ! total floating point operations since the first call
  real(4),   public :: papi_real_time_o ! total realtime since the first PAPI_flops() call
  real(4),   public :: papi_proc_time_o ! total process time since the first PAPI_flops() call
  real(4),   public :: papi_mflops      ! Mflop/s achieved since the previous call
  integer,   public :: papi_check
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Damp all data
  subroutine DEBUG_dampdata( &
      basename, & !--- [IN]
      var,      & !--- [IN]
      var_pl    ) !--- [IN]
    use mod_misc, only: &
       MISC_make_idstr, &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_prc_me
    implicit none

    character(len=*), intent(in) :: basename
    real(RP),          intent(in) :: var   (:,:,:,:)
    real(RP),          intent(in) :: var_pl(:,:,:,:)

    integer :: shp(4)

    character(LEN=ADM_MAXFNAME) :: fname

    integer :: fid
    !---------------------------------------------------------------------------

    shp(:) = shape(var)

    call MISC_make_idstr(fname,trim(basename),'pe',ADM_prc_me)
    fid = MISC_get_available_fid()
    open( unit   = fid,                           &
          file   = trim(fname),                   &
          form   = 'unformatted',                 &
          access = 'direct',                      &
          recl   = shp(1)*shp(2)*shp(3)*shp(4)*8, &
          status = 'unknown'                      )

       write(fid,rec=1) var

    close(fid)

    if ( ADM_have_pl ) then
       shp(:) = shape(var_pl)

       fname = trim(basename)//'.pl'
       fid = MISC_get_available_fid()
       open( unit   = fid,                           &
             file   = trim(fname),                   &
             form   = 'unformatted',                 &
             access = 'direct',                      &
             recl   = shp(1)*shp(2)*shp(3)*shp(4)*8, &
             status = 'unknown'                      )

          write(fid,rec=1) var_pl

       close(fid)

    endif

  end subroutine DEBUG_dampdata

  !-----------------------------------------------------------------------------
  !> Damp all data
  subroutine DEBUG_dampascii4D( &
      basename, & !--- [IN]
      var,      & !--- [IN]
      var_pl    ) !--- [IN]
    use mod_misc, only: &
       MISC_make_idstr, &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_prc_me
    implicit none

    character(len=*), intent(in) :: basename
    real(RP),          intent(in) :: var   (:,:,:,:)
    real(RP),          intent(in) :: var_pl(:,:,:,:)

    integer :: shp(4)

    character(LEN=ADM_MAXFNAME) :: fname

    integer :: fid
    integer :: i1,i2,i3,i4
    !---------------------------------------------------------------------------

    shp(:) = shape(var)

    call MISC_make_idstr(fname,trim(basename),'txt',ADM_prc_me)
    fid = MISC_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'unknown'    )

       do i4 = 1, shp(4)
       do i3 = 1, shp(3)
       do i2 = 1, shp(2)
       do i1 = 1, shp(1)
          write(fid,*) "(",i1,",",i2,",",i3,",",i4,")=",var(i1,i2,i3,i4)
       enddo
       enddo
       enddo
       enddo

    close(fid)

    if ( ADM_have_pl ) then
       shp(:) = shape(var_pl)

       fname = trim(basename)//'.txtpl'
       fid = MISC_get_available_fid()
       open( unit   = fid,         &
             file   = trim(fname), &
             form   = 'formatted', &
             status = 'unknown'    )

          do i4 = 1, shp(4)
          do i3 = 1, shp(3)
          do i2 = 1, shp(2)
          do i1 = 1, shp(1)
             write(fid,*) "(",i1,",",i2,",",i3,",",i4,")=",var_pl(i1,i2,i3,i4)
          enddo
          enddo
          enddo
          enddo

       close(fid)

    endif

  end subroutine DEBUG_dampascii4D

  !-----------------------------------------------------------------------------
  !> Damp all data
  subroutine DEBUG_dampascii3D( &
      basename, & !--- [IN]
      var,      & !--- [IN]
      var_pl    ) !--- [IN]
    use mod_misc, only: &
       MISC_make_idstr, &
       MISC_get_available_fid
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_prc_me
    implicit none

    character(len=*), intent(in) :: basename
    real(RP),          intent(in) :: var   (:,:,:)
    real(RP),          intent(in) :: var_pl(:,:,:)

    integer :: shp(3)

    character(LEN=ADM_MAXFNAME) :: fname

    integer :: fid
    integer :: i1,i2,i3
    !---------------------------------------------------------------------------

    shp(:) = shape(var)

    call MISC_make_idstr(fname,trim(basename),'txt',ADM_prc_me)
    fid = MISC_get_available_fid()
    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'unknown'    )

       do i3 = 1, shp(3)
       do i2 = 1, shp(2)
       do i1 = 1, shp(1)
          write(fid,*) "(",i1,",",i2,",",i3,")=",var(i1,i2,i3)
       enddo
       enddo
       enddo

    close(fid)

    if ( ADM_have_pl ) then
       shp(:) = shape(var_pl)

       fname = trim(basename)//'.txtpl'
       fid = MISC_get_available_fid()
       open( unit   = fid,         &
             file   = trim(fname), &
             form   = 'formatted', &
             status = 'unknown'    )

          do i3 = 1, shp(3)
          do i2 = 1, shp(2)
          do i1 = 1, shp(1)
             write(fid,*) "(",i1,",",i2,",",i3,")=",var_pl(i1,i2,i3)
          enddo
          enddo
          enddo

       close(fid)

    endif

  end subroutine DEBUG_dampascii3D

  !-----------------------------------------------------------------------------
  function DEBUG_rapid( rapname ) result(id)
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    if ( DEBUG_rapnmax >= 1 ) then
       do id = 1, DEBUG_rapnmax
          if( trim(rapname) == trim(DEBUG_rapname(id)) ) return
       enddo
    endif

    DEBUG_rapnmax     = DEBUG_rapnmax + 1
    id                = DEBUG_rapnmax
    DEBUG_rapname(id) = trim(rapname)
    DEBUG_raptstr(id) = 0.0_RP
    DEBUG_rapttot(id) = 0.0_RP
    DEBUG_rapnstr(id) = 0
    DEBUG_rapnend(id) = 0

  end function DEBUG_rapid

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapstart( rapname )
    use mod_adm, only: &
       ADM_MPItime
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    DEBUG_raptstr(id) = ADM_MPItime()
    DEBUG_rapnstr(id) = DEBUG_rapnstr(id) + 1

    !write(ADM_LOG_FID,*) rapname, DEBUG_rapnstr(id)

#ifdef _FAPP_
    call fapp_start( rapname, id, 1 )
#endif
#ifdef _FINEPA_
    call START_COLLECTION( rapname )
#endif

    return
  end subroutine DEBUG_rapstart

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapend( rapname )
    use mod_adm, only: &
       ADM_MPItime
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    DEBUG_rapttot(id) = DEBUG_rapttot(id) + ( ADM_MPItime()-DEBUG_raptstr(id) )
    DEBUG_rapnend(id) = DEBUG_rapnend(id) + 1

#ifdef _FAPP_
    call fapp_stop( rapname, id, 1 )
#endif
#ifdef _FINEPA_
    call STOP_COLLECTION( rapname )
#endif

    return
  end subroutine DEBUG_rapend

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapreport
    use mod_adm, only: &
       ADM_COMM_WORLD, &
       ADM_prc_all,    &
       ADM_proc_stop
    implicit none

    real(RP) :: sendbuf(1)
    real(RP) :: recvbuf(ADM_prc_all)

    real(RP) :: globalavg, globalmax, globalmin
#ifdef PAPI_OPS
    real(RP) :: globalsum, total_flops
#endif

    integer :: datatype

    integer :: ierr
    integer :: id
    !---------------------------------------------------------------------------

    if ( RP == DP ) then
       datatype = MPI_DOUBLE_PRECISION
    elseif( RP == SP ) then
       datatype = MPI_REAL
    else
       write(*,*) 'xxx precision is not supportd'
       call ADM_proc_stop
    endif

    if ( DEBUG_rapnmax >= 1 ) then

       do id = 1, DEBUG_rapnmax
          if ( DEBUG_rapnstr(id) /= DEBUG_rapnend(id) ) then
              write(*,*) '*** Mismatch Report',id,DEBUG_rapname(id),DEBUG_rapnstr(id),DEBUG_rapnend(id)
          endif
       enddo

       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '*** Computational Time Report'

       !do id = 1, DEBUG_rapnmax
       !   write(ADM_LOG_FID,'(1x,A,I3.3,A,A,A,F10.3,A,I7)') &
       !   '*** ID=',id,' : ',DEBUG_rapname(id),' T=',DEBUG_rapttot(id),' N=',DEBUG_rapnstr(id)
       !enddo

       do id = 1, DEBUG_rapnmax
          sendbuf(1) = DEBUG_rapttot(id)
          call MPI_Allgather( sendbuf,        &
                              1,              &
                              datatype,       &
                              recvbuf,        &
                              1,              &
                              datatype,       &
                              ADM_COMM_WORLD, &
                              ierr            )

          globalavg = sum( recvbuf(:) ) / real(ADM_prc_all,kind=RP)
          globalmax = maxval( recvbuf(:) )
          globalmin = minval( recvbuf(:) )

          write(ADM_LOG_FID,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,F10.3,A,I7)') &
                            '*** ID=',   id,                &
                            ' : ',       DEBUG_rapname(id), &
                            '  T(avg)=', globalavg,         &
                            ', T(max)=', globalmax,         &
                            ', T(min)=', globalmin,         &
                            ', N=',      DEBUG_rapnstr(id)
       enddo
    else
       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '*** Computational Time Report: NO item.'
    endif

#ifdef PAPI_OPS
    ! [add] PAPI R.Yoshida 20121022
    !write(ADM_LOG_FID,*) ' *** Type: Instructions'
    !write(ADM_LOG_FID,*) ' --- Real Time:',papi_real_time_i*2.0_RP,' Proc. Time:',papi_proc_time_i*2.0_RP
    !write(ADM_LOG_FID,*) ' --- flop inst:',papi_flpins*2,'  Gflins/s:',papi_mflins*2.0_RP/1.0d3  !GIGA
    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '********* PAPI report *********'
    write(ADM_LOG_FID,*) '*** Type: Operations'
    write(ADM_LOG_FID,*) '--- Wall clock Time      [sec] (this PE):', papi_real_time_o
    write(ADM_LOG_FID,*) '--- Processor Time       [sec] (this PE):', papi_proc_time_o
    write(ADM_LOG_FID,*) '--- Floating Operations [FLOP] (this PE):', papi_flpops
    write(ADM_LOG_FID,*) '--- FLOPS by PAPI     [MFLOPS] (this PE):', papi_mflops
    write(ADM_LOG_FID,*) '--- FLOP / Time       [MFLOPS] (this PE):', papi_flpops / papi_proc_time_o / 1024.0_RP**2 !GIGA
    write(ADM_LOG_FID,*)

    sendbuf(1) = real(papi_proc_time_o,kind=RP)
    call MPI_Allgather( sendbuf,        &
                        1,              &
                        datatype,       &
                        recvbuf,        &
                        1,              &
                        datatype,       &
                        ADM_COMM_WORLD, &
                        ierr            )

    globalavg = sum( recvbuf(:) ) / real(ADM_prc_all,kind=RP)
    globalmax = maxval( recvbuf(:) )
    globalmin = minval( recvbuf(:) )

    call COMM_Stat_avg( real(papi_proc_time_o,kind=RP), globalavg )
    call COMM_Stat_max( real(papi_proc_time_o,kind=RP), globalmax )
    call COMM_Stat_min( real(papi_proc_time_o,kind=RP), globalmin )

    write(ADM_LOG_FID,'(1x,A,F10.3,A,F10.3,A,F10.3)') &
                      '--- Processor Time        [sec] (avg)=', globalavg, &
                                                    ', (max)=', globalmax, &
                                                    ', (min)=', globalmin

    sendbuf(1) = real(papi_flpops,kind=RP)
    call MPI_Allgather( sendbuf,        &
                        1,              &
                        datatype,       &
                        recvbuf,        &
                        1,              &
                        datatype,       &
                        ADM_COMM_WORLD, &
                        ierr            )

    globalsum = sum( recvbuf(:) )
    globalavg = globalsum / real(ADM_prc_all,kind=RP)
    globalmax = maxval( recvbuf(:) )
    globalmin = minval( recvbuf(:) )

    total_flops = globalsum / globalmax / 1024.0_RP**3

    write(ADM_LOG_FID,'(1x,A,F10.3,A,F10.3,A,F10.3)') &
                      '--- Floating Operations [GFLOP] (avg)=', globalavg / 1024.0_RP**3, &
                                                    ', (max)=', globalmax / 1024.0_RP**3, &
                                                    ', (min)=', globalmin / 1024.0_RP**3
    write(ADM_LOG_FID,'(1x,A,F10.3)') &
                      '--- Total Flops [GFLOPS] (all PE):',total_flops

    call PAPIF_shutdown
#endif

    return
  end subroutine DEBUG_rapreport

end module mod_debug
