!-------------------------------------------------------------------------------
!> module Monitor
!!
!! @par Description
!!          Monitor output module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-22 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_monitor
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_FILECHR
  use mod_fileio_h, only: &
     FIO_HSHORT, &
     FIO_HMID
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
  !++ included parameters
  !
  include "inc_index.h"

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
  integer,                   private,      save :: MONIT_FID = -1

  character(len=IO_FILECHR), private,      save :: MONITOR_OUT_BASENAME = 'monitor'

  integer,                   private, parameter :: MONIT_req_limit = 1000 !> number limit for monitor item request
  character(len=FIO_HSHORT), private,      save :: MONIT_req_item   (MONIT_req_limit)

  integer,                   private,              save :: MONIT_req_nmax = 0 !> number of requested item
  character(len=FIO_HSHORT), private, allocatable, save :: MONIT_item   (:)
  character(len=FIO_HMID),   private, allocatable, save :: MONIT_desc   (:)
  character(len=FIO_HSHORT), private, allocatable, save :: MONIT_unit   (:)
  character(len=FIO_HSHORT), private, allocatable, save :: MONIT_ktype  (:)
  integer,                   private, allocatable, save :: MONIT_kmax   (:)
  real(8),                   private, allocatable, save :: MONIT_var    (:)

  integer,                   private,              save :: MONIT_id_count = 1 !> number of registered item

  real(8), private, parameter :: eps = 1.D-10 !> epsilon for timesec

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  subroutine MONIT_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    NAMELIST / PARAM_MONITOR / &
       MONITOR_OUT_BASENAME

    character(len=FIO_HSHORT) :: ITEM  !> name of monitor item

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
    endif

    allocate( MONIT_item   (MONIT_req_nmax) )
    allocate( MONIT_desc   (MONIT_req_nmax) )
    allocate( MONIT_unit   (MONIT_req_nmax) )
    allocate( MONIT_ktype  (MONIT_req_nmax) )
    allocate( MONIT_kmax   (MONIT_req_nmax) )
    allocate( MONIT_var    (MONIT_req_nmax) )

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
  subroutine MONIT_reg( &
      itemid, &
      item,   &
      desc,   &
      unit,   &
      ktype   )
    implicit none

    integer,          intent(out) :: itemid
    character(len=*), intent( in) :: item
    character(len=*), intent( in) :: desc
    character(len=*), intent( in) :: unit
    character(len=*), intent( in) :: ktype

    character(len=8) :: lname

    integer :: n, nmax, reqid
    !---------------------------------------------------------------------------

    !--- search existing item
    itemid = -1
    nmax = min( MONIT_id_count, MONIT_req_nmax )
    do n = 1, nmax
       if ( trim(item) == trim(MONIT_item(n)) ) then ! match existing item
          itemid = n
          return
       endif
    enddo

    if ( itemid < 0 ) then ! request-register matching check
       do n = 1, MONIT_req_nmax
          if ( trim(item) == MONIT_req_item(n) ) then
             itemid = MONIT_id_count
             reqid  = n
             MONIT_id_count = MONIT_id_count + 1

             ! new file registration
             MONIT_item(itemid) = trim(item)
             MONIT_desc(itemid) = trim(desc)
             MONIT_unit(itemid) = trim(unit)
             if    ( trim(ktype) == '2D' ) then
                MONIT_ktype(itemid) = 'ZSFC'
                MONIT_kmax (itemid) = 1
             elseif( trim(ktype) == '3D' ) then
                write(lname,'(A,I4.4)') 'ZDEF', KMAX
                MONIT_ktype(itemid) = lname
                MONIT_kmax (itemid) = KMAX
             endif
             MONIT_var(itemid) = 0.D0

             if( IO_L ) write(IO_FID_LOG,*) '*** [MONIT] Item registration No.= ', itemid
             if( IO_L ) write(IO_FID_LOG,*) '] Name           : ', trim(MONIT_item (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Description    : ', trim(MONIT_desc (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Unit           : ', trim(MONIT_unit (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Vert. type     : ', trim(MONIT_ktype(itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] # of layer     : ', MONIT_kmax(itemid)
          endif
       enddo
    endif

    return
  end subroutine MONIT_reg

  !-----------------------------------------------------------------------------
  subroutine MONIT_put( &
      itemid, &
      var     )
    use mod_geometrics, only: &
       area => GEOMETRICS_area, &
       vol  => GEOMETRICS_vol
    implicit none

    integer, intent(in) :: itemid
    real(8), intent(in) :: var(:,:,:)

    real(8) :: total

    integer :: k, i, j
    !---------------------------------------------------------------------------

    total = 0.D0
    if ( MONIT_kmax(itemid) == KMAX ) then ! 3D
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          total = total + var(k,i,j) * vol(k,i,j)
       enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          total = total + var(1,i,j) * area(1,i,j)
       enddo
       enddo
    endif

    MONIT_var(itemid) = total ! overwrite by last put

    return
  end subroutine MONIT_put

  !-----------------------------------------------------------------------------
  subroutine MONIT_in( &
      var,   &
      item,  &
      desc,  &
      unit,  &
      ktype  )
    implicit none

    real(8),          intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: ktype

    integer :: itemid
    !---------------------------------------------------------------------------

    call MONIT_reg( itemid, item, desc, unit, ktype )
    
    if ( itemid > 0 ) then
       call MONIT_put( itemid, var(:,:,:) )
    endif

    return
  end subroutine MONIT_in

  !-----------------------------------------------------------------------------
  subroutine MONIT_write( memo )
    use mod_time, only: &
       NOWSTEP => TIME_NOWSTEP
    implicit none

    character(len=4), intent(in) :: memo

    logical, save :: firsttime = .true.

    integer :: n
    !---------------------------------------------------------------------------

    if( MONIT_id_count == 1 ) return

    if (firsttime) then
       firsttime = .false.
       call MONIT_writeheader
    endif

    write(MONIT_FID,'(A,i7,A,A,A)',advance='no') 'STEP=',NOWSTEP,' (',memo,')'
    do n = 1, MONIT_id_count-1
       write(MONIT_FID,'(A,E15.8)',advance='no') ' ',MONIT_var(n)
    enddo
    write(MONIT_FID,*)

    return
  end subroutine MONIT_write

  !-----------------------------------------------------------------------------
  subroutine MONIT_writeheader
    use mod_stdio, only : &
       IO_get_available_fid, &
       IO_make_idstr,        &
       IO_FILECHR
    use mod_process, only : &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    character(len=IO_FILECHR) :: fname !< name of monitor file for each process

    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [MONITOR] Output item list '
    if( IO_L ) write(IO_FID_LOG,*) '*** Number of monitor item :', MONIT_req_nmax
    if( IO_L ) write(IO_FID_LOG,*) 'NAME           :description                                     ', &
                                   '               :UNIT           :Layername'
    if( IO_L ) write(IO_FID_LOG,*) '=====================================================', &
                                   '====================================================='
    do n = 1, MONIT_id_count-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,A)') MONIT_item(n), MONIT_desc(n), MONIT_unit(n), MONIT_ktype(n)
    enddo
    if( IO_L ) write(IO_FID_LOG,*) '=====================================================', &
                                   '====================================================='

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
    do n = 1, MONIT_id_count-1
       write(MONIT_FID,'(A,A16)',advance='no') MONIT_item(n)
    enddo
    write(MONIT_FID,*)

    return
  end subroutine MONIT_writeheader

  !-----------------------------------------------------------------------------
  subroutine MONIT_finalize
    use mod_stdio, only : &
       IO_make_idstr
    use mod_process, only : &
       PRC_myrank
    implicit none

    character(len=IO_FILECHR) :: fname !< name of monitor file for each process
    !---------------------------------------------------------------------------

    if ( MONIT_FID > 0 ) then
       close(MONIT_FID)

       call IO_make_idstr(fname,trim(MONITOR_OUT_BASENAME),'pe',PRC_myrank)
       if( IO_L ) write(IO_FID_LOG,*) '*** [MONITOR] File Close'
       if( IO_L ) write(IO_FID_LOG,*) '*** closed filename: ', fname
    endif

    return
  end subroutine MONIT_finalize

end module mod_monitor
!-------------------------------------------------------------------------------
