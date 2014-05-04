!-------------------------------------------------------------------------------
!> module MONITOR
!!
!! @par Description
!!          Monitor output module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-03-22 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_monitor
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: MONIT_setup
  public :: MONIT_reg
  public :: MONIT_put
  public :: MONIT_in
  public :: MONIT_write
  public :: MONIT_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: MONIT_writeheader

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer                :: MONIT_FID = -1                    !< fileID for monitor output file

  character(len=H_LONG)  :: MONITOR_OUT_BASENAME  = 'monitor' !< filename of monitor output
  logical                :: MONITOR_USEDEVATION   = .true.    !< use deviation from first step?
  integer                :: MONITOR_STEP_INTERVAL = 1         !< step interval

  integer, parameter     :: MONIT_req_limit = 1000            !< number limit for item request
  integer                :: MONIT_req_nmax = 0                !< number of requested item
  character(len=H_SHORT) :: MONIT_req_item(MONIT_req_limit)   !< name of requested monitor item

  integer                             :: MONIT_id_count = 0 !< number of item to output
  character(len=H_SHORT), allocatable :: MONIT_item (:)     !< name                of the item
  character(len=H_MID)  , allocatable :: MONIT_desc (:)     !< description         of the item
  character(len=H_SHORT), allocatable :: MONIT_unit (:)     !< unit                of the item
  character(len=H_SHORT), allocatable :: MONIT_ktype(:)     !< vertical layer type of the item
  integer,                allocatable :: MONIT_kmax (:)     !< # of vertical grid  of the item
  real(RP),               allocatable :: MONIT_var  (:)     !< value               of the item
  logical,                allocatable :: MONIT_first(:)     !< first time?         of the item
  real(RP),               allocatable :: MONIT_var0 (:)     !< value at first time of the item

  real(RP), parameter :: eps = 1.E-10_RP !< epsilon for timesec

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MONIT_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_MONITOR / &
       MONITOR_OUT_BASENAME, &
       MONITOR_USEDEVATION,  &
       MONITOR_STEP_INTERVAL

    character(len=H_SHORT) :: ITEM  !> name of monitor item

    NAMELIST / MONITITEM / &
       ITEM

    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[MONITOR]/Categ[IO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MONITOR,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MONITOR. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MONITOR)

    ! listup monitor request
    rewind(IO_FID_CONF)
    do n = 1, MONIT_req_limit
       read(IO_FID_CONF,nml=MONITITEM,iostat=ierr)
       if( ierr /= 0 ) exit
    enddo
    MONIT_req_nmax = n - 1

    if    ( MONIT_req_nmax > MONIT_req_limit ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** request of monitor file is exceed! n >', MONIT_req_limit
    elseif( MONIT_req_nmax == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** No monitor file specified.'
       return
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Number of requested monitor item: ', MONIT_req_nmax
       if( IO_L ) write(IO_FID_LOG,*) '*** Monitor output interval [step]  : ', MONITOR_STEP_INTERVAL
    endif

    allocate( MONIT_item (MONIT_req_nmax) )
    allocate( MONIT_desc (MONIT_req_nmax) )
    allocate( MONIT_unit (MONIT_req_nmax) )
    allocate( MONIT_ktype(MONIT_req_nmax) )
    allocate( MONIT_kmax (MONIT_req_nmax) )
    allocate( MONIT_var  (MONIT_req_nmax) )
    allocate( MONIT_var0 (MONIT_req_nmax) )
    allocate( MONIT_first(MONIT_req_nmax) )

    rewind(IO_FID_CONF)
    do n = 1, MONIT_req_nmax
       ! set default
       ITEM  = 'unknown'

       read(IO_FID_CONF,nml=MONITITEM,iostat=ierr)
       if( ierr /= 0 ) exit

       MONIT_req_item(n) = ITEM
    enddo

    return
  end subroutine MONIT_setup

  !-----------------------------------------------------------------------------
  !> Search existing item, or matching check between requested and registered item
  subroutine MONIT_reg( &
      itemid, &
      item,   &
      desc,   &
      unit,   &
      ndim    )
    implicit none

    integer,          intent(out) :: itemid !< index number of the item
    character(len=*), intent(in)  :: item   !< name         of the item
    character(len=*), intent(in)  :: desc   !< description  of the item
    character(len=*), intent(in)  :: unit   !< unit         of the item
    integer,          intent(in)  :: ndim   !< dimension    of the item

    character(len=8) :: lname

    integer :: n, nmax, reqid
    !---------------------------------------------------------------------------

    !--- search existing item
    itemid = -1
    nmax = min( MONIT_id_count, MONIT_req_nmax )
    do n = 1, nmax
       if ( item == MONIT_item(n) ) then ! match existing item
          itemid = n
          return
       endif
    enddo

    if ( itemid < 0 ) then ! request-register matching check
       do n = 1, MONIT_req_nmax
          if ( item == MONIT_req_item(n) ) then
             MONIT_id_count = MONIT_id_count + 1
             itemid = MONIT_id_count
             reqid  = n

             ! new file registration
             MONIT_item(itemid) = trim(item)
             MONIT_desc(itemid) = trim(desc)
             MONIT_unit(itemid) = trim(unit)
             if    ( ndim == 2 ) then
                MONIT_ktype(itemid) = 'ZSFC'
                MONIT_kmax (itemid) = 1
             elseif( ndim == 3 ) then
                write(lname,'(A,I4.4)') 'ZDEF', KMAX
                MONIT_ktype(itemid) = lname
                MONIT_kmax (itemid) = KMAX
             endif
             MONIT_var  (itemid) = 0.0_RP
             MONIT_var0 (itemid) = 0.0_RP
             MONIT_first(itemid) = .true.

             if( IO_L ) write(IO_FID_LOG,*) ' *** [MONIT] Item registration No.= ', itemid
             if( IO_L ) write(IO_FID_LOG,*) ' ] Name           : ', trim(MONIT_item (itemid))
             if( IO_L ) write(IO_FID_LOG,*) ' ] Description    : ', trim(MONIT_desc (itemid))
             if( IO_L ) write(IO_FID_LOG,*) ' ] Unit           : ', trim(MONIT_unit (itemid))
             if( IO_L ) write(IO_FID_LOG,*) ' ] Vert. type     : ', trim(MONIT_ktype(itemid))
             if( IO_L ) write(IO_FID_LOG,*) ' ] # of layer     : ', MONIT_kmax(itemid)
          endif
       enddo
    endif

    return
  end subroutine MONIT_reg

  !-----------------------------------------------------------------------------
  !> Put total value to the monitor buffer
  subroutine MONIT_put( &
      itemid, &
      var     )
    use scale_stats, only: &
       STAT_total
    implicit none

    integer,  intent(in) :: itemid     !< index number of the item
    real(RP), intent(in) :: var(:,:,:) !< value

    real(RP) :: total
    !---------------------------------------------------------------------------

    if( itemid <= 0 ) return

    call STAT_total( total, var(:,:,:), MONIT_item(itemid) )

    MONIT_var(itemid) = total ! overwrite by last put

    if ( MONITOR_USEDEVATION ) then
       if ( MONIT_first(itemid) ) then
          MONIT_var  (itemid) = 0.0_RP
          MONIT_var0 (itemid) = total
          MONIT_first(itemid) = .false.
       else
          MONIT_var  (itemid) = total - MONIT_var0(itemid) ! overwrite by last put
       endif
    else
       MONIT_var(itemid) = total ! overwrite by last put
    endif

    return
  end subroutine MONIT_put

  !-----------------------------------------------------------------------------
  !> Wrapper routine of MONIT_reg+MONIT_put
  subroutine MONIT_in( &
      var,  &
      item, &
      desc, &
      unit, &
      ndim  )
    implicit none

    real(RP),         intent(in) :: var(:,:,:) !< value
    character(len=*), intent(in) :: item       !< name        of the item
    character(len=*), intent(in) :: desc       !< description of the item
    character(len=*), intent(in) :: unit       !< unit        of the item
    integer,          intent(in) :: ndim       !< dimension   of the item

    integer :: itemid
    !---------------------------------------------------------------------------

    call MONIT_reg( itemid, item, desc, unit, ndim )
    call MONIT_put( itemid, var(:,:,:) )

    return
  end subroutine MONIT_in

  !-----------------------------------------------------------------------------
  !> Flush monitor buffer to formatted file
  subroutine MONIT_write( memo )
    use scale_time, only: &
       NOWSTEP => TIME_NOWSTEP
    implicit none

    character(len=4), intent(in) :: memo !< note

    logical, save :: firsttime = .true.

    integer :: n
    !---------------------------------------------------------------------------

    if( MONIT_id_count == 0 ) return

    call PROF_rapstart('FILE O ASCII')

    if (firsttime) then
       firsttime = .false.
       call MONIT_writeheader
    endif

    if ( MONIT_FID > 0 ) then

       if ( mod(NOWSTEP-1,MONITOR_STEP_INTERVAL) == 0 ) then
          write(MONIT_FID,'(A,i7,A,A,A)',advance='no') 'STEP=',NOWSTEP,' (',memo,')'
          do n = 1, MONIT_id_count
             write(MONIT_FID,'(A,E15.8)',advance='no') ' ',MONIT_var(n)
          enddo
          write(MONIT_FID,*)

          if( IO_L ) write(IO_FID_LOG,*) '*** Write monitor'
       endif

    endif

    call PROF_rapend  ('FILE O ASCII')

    return
  end subroutine MONIT_write

  !-----------------------------------------------------------------------------
  !> Open file and write header at the first time
  subroutine MONIT_writeheader
    use scale_process, only: &
       PRC_myrank, &
       PRC_master, &
       PRC_MPIstop
    implicit none

    character(len=H_LONG) :: fname !< name of monitor file for each process

    logical :: MONIT_L
    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [MONITOR] Output item list '
    if( IO_L ) write(IO_FID_LOG,*) '*** Number of monitor item :', MONIT_req_nmax
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)') 'NAME           :description                                     ', &
                                            '               :UNIT           :Layername'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)') '=====================================================', &
                                            '====================================================='
    do n = 1, MONIT_id_count
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,A)') MONIT_item(n), MONIT_desc(n), MONIT_unit(n), MONIT_ktype(n)
    enddo
    if( IO_L ) write(IO_FID_LOG,'(1x,A,A)') '=====================================================', &
                                            '====================================================='

    if ( PRC_myrank == PRC_master ) then ! master node
       MONIT_L = .true.
    else
       MONIT_L = IO_LOG_ALLNODE
    endif

    if ( MONIT_L ) then

       !--- Open logfile
       MONIT_FID = IO_get_available_fid()
       call IO_make_idstr(fname,trim(MONITOR_OUT_BASENAME),'pe',PRC_myrank)
       open( unit   = MONIT_FID,  &
             file   = trim(fname),  &
             form   = 'formatted',  &
             iostat = ierr          )
       if ( ierr /= 0 ) then
          write(*,*) 'xxx File open error! :', trim(fname)
          call PRC_MPIstop
       endif

       write(MONIT_FID,'(A)',advance='no') '                   '
       do n = 1, MONIT_id_count
          write(MONIT_FID,'(A,A16)',advance='no') MONIT_item(n)
       enddo
       write(MONIT_FID,*)

    endif

    return
  end subroutine MONIT_writeheader

  !-----------------------------------------------------------------------------
  !> Close file
  subroutine MONIT_finalize
    use scale_process, only: &
       PRC_myrank
    implicit none

    character(len=H_LONG) :: fname !< name of monitor file for each process
    !---------------------------------------------------------------------------

    if ( MONIT_FID > 0 ) then
       close(MONIT_FID)

       call IO_make_idstr(fname,trim(MONITOR_OUT_BASENAME),'pe',PRC_myrank)
       if( IO_L ) write(IO_FID_LOG,*) '*** [MONITOR] File Close'
       if( IO_L ) write(IO_FID_LOG,*) '*** closed filename: ', fname
    endif

    return
  end subroutine MONIT_finalize

end module scale_monitor
!-------------------------------------------------------------------------------
