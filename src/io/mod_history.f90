!-------------------------------------------------------------------------------
!> module History
!!
!! @par Description
!!          History output module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-12-05 (H.Yashiro)  [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!!
!<
!-------------------------------------------------------------------------------
module mod_history
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
  public :: HIST_setup
  public :: HIST_reg
  public :: HIST_put
  public :: HIST_in
  public :: HIST_write
  public :: HIST_outputlist

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
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  character(len=IO_FILECHR), private,      save :: HISTORY_OUT_BASENAME = 'history'
  integer,                   private,      save :: HISTORY_DTYPE

  integer,                   private, parameter :: HIST_req_limit = 1000 !> number limit for history item request
  character(len=FIO_HSHORT), private,      save :: HIST_req_item   (HIST_req_limit)
  real(8),                   private,      save :: HIST_req_tintsec(HIST_req_limit)
  logical,                   private,      save :: HIST_req_tavg   (HIST_req_limit)

  integer,                   private,              save :: HIST_req_nmax = 0 !> number of requested item
  character(len=FIO_HSHORT), private, allocatable, save :: HIST_item   (:)
  character(len=FIO_HMID),   private, allocatable, save :: HIST_desc   (:)
  character(len=FIO_HSHORT), private, allocatable, save :: HIST_unit   (:)
  character(len=FIO_HSHORT), private, allocatable, save :: HIST_ktype  (:)
  integer,                   private, allocatable, save :: HIST_kmax   (:)
  real(8),                   private, allocatable, save :: HIST_tintsec(:)
  logical,                   private, allocatable, save :: HIST_tavg   (:)

  real(8),                   private, allocatable, save :: HIST_varsum (:,:,:,:)
  integer,                   private, allocatable, save :: HIST_step   (:)
  real(8),                   private, allocatable, save :: HIST_tstrsec(:)
  real(8),                   private, allocatable, save :: HIST_tsumsec(:)

  integer,                   private,              save :: HIST_id_count = 1 !> number of registered item

  real(8), private, parameter :: eps = 1.D-10 !> epsilon for timesec

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  subroutine HIST_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only: &
       TIME_ymdhms2sec
    use mod_fileio_h, only: &
       FIO_REAL4, &
       FIO_REAL8, &
       FIO_preclist
    implicit none

    real(8)                   :: HISTORY_DEFAULT_TINTERVAL = 1.D0
    character(len=FIO_HSHORT) :: HISTORY_DEFAULT_TUNIT     = "MIN"
    logical                   :: HISTORY_DEFAULT_AVERAGE   = .false.
    character(len=FIO_HSHORT) :: HISTORY_DATATYPE          = "REAL4"

    NAMELIST / PARAM_HISTORY / &
       HISTORY_OUT_BASENAME,      &
       HISTORY_DEFAULT_TINTERVAL, &
       HISTORY_DEFAULT_TUNIT,     &
       HISTORY_DEFAULT_AVERAGE,   &
       HISTORY_DATATYPE

    character(len=FIO_HSHORT) :: ITEM  !> name of history item
    real(8)                   :: TINT  !> time interval to output
    character(len=FIO_HSHORT) :: TUNIT !> time unit
    logical                   :: TAVG  !> time average  to output

    NAMELIST / HISTITEM / &
       ITEM,  &
       TINT,  &
       TUNIT, &
       TAVG

    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HISTORY]/Categ[IO]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_HISTORY,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_HISTORY. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_HISTORY)

    if    ( trim(HISTORY_DATATYPE) == 'REAL4' ) then
       HISTORY_DTYPE = FIO_REAL4
    elseif( trim(HISTORY_DATATYPE) == 'REAL8' ) then
       HISTORY_DTYPE = FIO_REAL8
    else
       write(*,*) 'xxx Not appropriate DATATYPE. Check!', HISTORY_DATATYPE
       call PRC_MPIstop
    endif

    ! listup history request
    rewind(IO_FID_CONF)
    do n = 1, HIST_req_limit
       read(IO_FID_CONF,nml=HISTITEM,iostat=ierr)
       if( ierr /= 0 ) exit
    enddo
    HIST_req_nmax = n - 1

    if    ( HIST_req_nmax > HIST_req_limit ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** request of history file is exceed! n >', HIST_req_limit
    elseif( HIST_req_nmax == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** No history file specified.'
       return
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Number of requested history item             : ', HIST_req_nmax
       if( IO_L ) write(IO_FID_LOG,*) '*** Output data type                             : ', HISTORY_DATATYPE
       if( IO_L ) write(IO_FID_LOG,*) '*** Memory usage for history data buffer [Mbyte] : ', &
                                      IMAX*JMAX*KMAX*HIST_req_nmax*FIO_preclist(HISTORY_DTYPE)/1024/1024
    endif

    allocate( HIST_item   (HIST_req_nmax) )
    allocate( HIST_desc   (HIST_req_nmax) )
    allocate( HIST_unit   (HIST_req_nmax) )
    allocate( HIST_ktype  (HIST_req_nmax) )
    allocate( HIST_kmax   (HIST_req_nmax) )
    allocate( HIST_tintsec(HIST_req_nmax) )
    allocate( HIST_tavg   (HIST_req_nmax) )

    allocate( HIST_varsum (KMAX,IMAX,JMAX,HIST_req_nmax) )
    allocate( HIST_step   (HIST_req_nmax) )
    allocate( HIST_tstrsec(HIST_req_nmax) )
    allocate( HIST_tsumsec(HIST_req_nmax) )

    rewind(IO_FID_CONF)
    do n = 1, HIST_req_nmax
       ! set default
       ITEM  = 'unknown'
       TINT  = HISTORY_DEFAULT_TINTERVAL
       TAVG  = HISTORY_DEFAULT_AVERAGE
       TUNIT = HISTORY_DEFAULT_TUNIT

       read(IO_FID_CONF,nml=HISTITEM,iostat=ierr)
       if( ierr /= 0 ) exit

       HIST_req_item(n) = ITEM
       call TIME_ymdhms2sec( HIST_req_tintsec(n), TINT, TUNIT )
       HIST_req_tavg(n) = TAVG

       if ( HIST_req_tintsec(n) <= 0.D0 ) then
          write(*,*) 'xxx Not appropriate time interval. Check!', ITEM, TINT
          call PRC_MPIstop
       endif
    enddo

    return
  end subroutine HIST_setup

  !-----------------------------------------------------------------------------
  subroutine HIST_reg( &
      itemid, &
      item,   &
      desc,   &
      unit,   &
      ktype   )
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
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
    nmax = min( HIST_id_count, HIST_req_nmax )
    do n = 1, nmax
       if ( trim(item) == trim(HIST_item(n)) ) then ! match existing item
          itemid = n
          return
       endif
    enddo

    if ( itemid < 0 ) then ! request-register matching check
       do n = 1, HIST_req_nmax
          if ( trim(item) == HIST_req_item(n) ) then
             itemid = HIST_id_count
             reqid  = n
             HIST_id_count = HIST_id_count + 1

             ! new file registration
             HIST_item(itemid) = trim(item)
             HIST_desc(itemid) = trim(desc)
             HIST_unit(itemid) = trim(unit)
             if    ( trim(ktype) == '2D' ) then
                HIST_ktype(itemid) = 'ZSFC'
                HIST_kmax (itemid) = 1
             elseif( trim(ktype) == '3D' ) then
                write(lname,'(A,I4.4)') 'ZDEF', KMAX
                HIST_ktype(itemid) = lname
                HIST_kmax (itemid) = KMAX
             endif
             HIST_tintsec(itemid) = HIST_req_tintsec(reqid)
             HIST_tavg   (itemid) = HIST_req_tavg(reqid)

             HIST_varsum(:,:,:,itemid) = 0.D0
             HIST_step        (itemid) = 1
             HIST_tstrsec     (itemid) = NOWSEC
             HIST_tsumsec     (itemid) = 0.D0

             if( IO_L ) write(IO_FID_LOG,*) '*** [HIST] Item registration No.= ', itemid
             if( IO_L ) write(IO_FID_LOG,*) '] Name           : ', trim(HIST_item (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Description    : ', trim(HIST_desc (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Unit           : ', trim(HIST_unit (itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] Vert. type     : ', trim(HIST_ktype(itemid))
             if( IO_L ) write(IO_FID_LOG,*) '] # of layer     : ', HIST_kmax(itemid)
             if( IO_L ) write(IO_FID_LOG,*) '] Interval [sec] : ', HIST_tintsec(itemid)
             if( IO_L ) write(IO_FID_LOG,*) '] Average?       : ', HIST_tavg   (itemid)
          endif
       enddo
    endif

    return
  end subroutine HIST_reg

  !-----------------------------------------------------------------------------
  subroutine HIST_put( &
      itemid, &
      var,    &
      dt      )
    implicit none

    integer, intent(in) :: itemid
    real(8), intent(in) :: var(:,:,:)
    real(8), intent(in) :: dt

    integer :: ksize, kstr, kend
    !---------------------------------------------------------------------------

    if ( HIST_kmax(itemid) == KMAX ) then ! 3D
       ksize = KMAX
       kstr  = KS
       kend  = KE
    else
       ksize = 1
       kstr  = 1
       kend  = 1
    endif

    if ( HIST_tavg(itemid) ) then
       HIST_varsum(1:ksize,1:IMAX,1:JMAX,itemid) = HIST_varsum(1:ksize,1:IMAX,1:JMAX,itemid) &
                                                 + var(kstr:kend,IS:IE,JS:JE) * dt
    else
       HIST_varsum(1:ksize,1:IMAX,1:JMAX,itemid) = var(kstr:kend,IS:IE,JS:JE)
    endif
    HIST_tsumsec(itemid) = HIST_tsumsec(itemid) + dt

    return
  end subroutine HIST_put

  !-----------------------------------------------------------------------------
  subroutine HIST_in( &
      var,   &
      item,  &
      desc,  &
      unit,  &
      ktype, &
      dt     )
    implicit none

    real(8),          intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: item
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: ktype
    real(8),          intent(in) :: dt

    integer :: itemid
    !---------------------------------------------------------------------------

    call HIST_reg( itemid, item, desc, unit, ktype )
    
    if ( itemid > 0 ) then
       call HIST_put( itemid, var(:,:,:), dt )
    endif

    return
  end subroutine HIST_in

  !-----------------------------------------------------------------------------
  subroutine HIST_write
    use mod_time, only: &
       NOWSEC => TIME_NOWSEC
    use mod_fileio, only: &
#ifdef CONFIG_HDF5
       FIO_output, &
       FIO_getfid
#else
       FIO_output
#endif
    implicit none

#ifdef CONFIG_HDF5
    integer :: fid
!    character(LEN=64) :: gname
#endif

    real(8) :: var(KMAX,IMAX,JMAX)

    logical, save :: firsttime = .true.

    integer :: n, ksize

    !---------------------------------------------------------------------------

    if( HIST_id_count == 1 ) return

    if (firsttime) then
       firsttime = .false.
       call HIST_outputlist
    endif

#ifdef CONFIG_HDF5
    call FIO_getfid(fid, trim(HISTORY_OUT_BASENAME), 1, 'SCALE3 HISTORY OUTPUT', '')
    call fio_new_group(fid, NOWSEC, 6)
#endif

    do n = 1, HIST_id_count-1
       
       if ( HIST_tsumsec(n) - HIST_tintsec(n) > -eps ) then
          ksize = HIST_kmax(n)

          if ( HIST_tavg(n) ) then
             var(1:HIST_kmax(n),:,:) = HIST_varsum(1:HIST_kmax(n),:,:,n) / HIST_tsumsec(n)
          else
             var(1:HIST_kmax(n),:,:) = HIST_varsum(1:HIST_kmax(n),:,:,n)
          endif

          call FIO_output( var(:,:,:),                     & ! data
                           trim(HISTORY_OUT_BASENAME),     & ! package name
                           'SCALE3 HISTORY OUTPUT',        & ! package desc
                           '',                             & ! package note
                           HIST_item(n),                   & ! data item name
                           HIST_desc(n),                   & ! data desc
                           '',                             & ! data note
                           HIST_unit(n),                   & ! data unit
                           HISTORY_DTYPE,                  & ! datatype
                           HIST_ktype(n),                  & ! vertical layer name
                           1,                              & ! start index of k
                           HIST_kmax(n),                   & ! end   index pf k
                           HIST_step(n),                   & ! step number
                           HIST_tstrsec(n),                & ! package name
                           HIST_tstrsec(n)+HIST_tsumsec(n) ) ! package name

          HIST_varsum(:,:,:,n) = 0.D0
          HIST_step(n)         = HIST_step(n) + 1
          HIST_tstrsec(n)      = NOWSEC
          HIST_tsumsec(n)      = 0.D0
       endif

    enddo

#ifdef CONFIG_HDF5
!    call fio_get_group(fid, gid)
    call fio_close_group(fid)
#endif

    return
  end subroutine HIST_write

  !-----------------------------------------------------------------------------
  subroutine HIST_outputlist
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [HIST] Output item list '
    if( IO_L ) write(IO_FID_LOG,*) '*** Number of history item :', HIST_req_nmax
    if( IO_L ) write(IO_FID_LOG,*) 'NAME           :UNIT           :Layername         :interval[sec]:avg'
    if( IO_L ) write(IO_FID_LOG,*) '============================================================================'
 
    do n = 1, HIST_id_count-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,1x,f13.3,1x,L)') HIST_item(n), HIST_unit(n), HIST_ktype(n), HIST_tintsec(n), HIST_tavg(n)
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '============================================================================'

    return
  end subroutine HIST_outputlist

end module mod_history
!-------------------------------------------------------------------------------
