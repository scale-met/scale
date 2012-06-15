!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!> module Gtool_History
!!
!! @par Description
!!          module library for history output 
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-06-11 (S.Nishizawa)  [new] imported from SCALE-LES
!!
!<
!-------------------------------------------------------------------------------
module gtool_history
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_process, only: &
     PRC_MPIstop
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L,       &
     IO_FILECHR
  use gtool_file_h, only: &
     File_HSHORT, &
     File_HMID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: HistoryInit
  public :: HistoryAddVariable
  public :: HistoryPutAxis
  public :: HistoryPutAdditionalAxis
  public :: HistoryPut
  public :: HistoryWriteAll
  public :: HistoryGet
  public :: HistoryOutputList
  public :: HistoryFinalize

  interface HistoryPut
     module procedure HistoryPut1D
     module procedure HistoryPut2D
     module procedure HistoryPut3D
  end interface HistoryPut

  interface HistoryGet
     module procedure HistoryGet1D
     module procedure HistoryGet2D
     module procedure HistoryGet3D
  end interface HistoryGet
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_index.h"
  include "inc_precision.h"

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
  character(len=File_HMID),   private,      save :: HISTORY_TITLE
  character(len=File_HMID),   private,      save :: HISTORY_SOURCE
  character(len=File_HMID),   private,      save :: HISTORY_INSTITUTION
  character(len=File_HMID),   private,      save :: HISTORY_TIME_UNITS

  integer,                    private, parameter :: History_req_limit = 1000 !> number limit for history item request
  character(len=IO_FILECHR),  private,      save :: History_req_basename(History_req_limit)
  character(len=File_HSHORT), private,      save :: History_req_item    (History_req_limit)
  real(8),                    private,      save :: History_req_tintsec (History_req_limit)
  logical,                    private,      save :: History_req_tavg    (History_req_limit)
  integer,                    private,      save :: History_req_dtype   (History_req_limit)

  integer,                    private,              save :: History_req_nmax = 0 !> number of requested item
  character(len=File_HSHORT), private, allocatable, save :: History_item   (:)
  integer,                    private, allocatable, save :: History_fid    (:)
  integer,                    private, allocatable, save :: History_vid    (:)
  real(8),                    private, allocatable, save :: History_tintsec(:)
  logical,                    private, allocatable, save :: History_tavg   (:)

  real(DP),                   private, allocatable, save :: History_varsum (:,:)
  logical,                    private, allocatable, save :: History_size   (:)
  real(8),                    private, allocatable, save :: History_tstrsec(:)
  real(8),                    private, allocatable, save :: History_tsumsec(:)

  integer,                    private,              save :: History_id_count = 1 !> number of registered item

  character(len=File_HSHORT), private, allocatable, save :: History_dim_name(:)
  integer,                    private, allocatable, save :: History_dim_size(:)
  character(len=File_HMID),   private, allocatable, save :: History_dim_desc(:)
  character(len=File_HSHORT), private, allocatable, save :: History_dim_units(:)
  integer,                    private, allocatable, save :: History_dim_type(:)

  real(8), private, parameter :: eps = 1.D-10 !> epsilon for timesec

contains
  !-----------------------------------------------------------------------------
  subroutine HistoryInit( &
       title, source, institution,                         & ! (in)
       dim_name, dim_size, dim_desc, dim_units,            & ! (in)
       dim_type,                                           & ! (in) optional
       default_basename,                                   & ! (in) optional
       default_tinterval, default_tunit, default_taverage, & ! (in) optional
       default_datatype                                    & ! (in) optional
       )
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_time, only: &
       TIME_ymdhms2sec
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8, &
       File_preclist
    implicit none

    character(len=*), intent(in)           :: title
    character(len=*), intent(in)           :: source
    character(len=*), intent(in)           :: institution
    character(len=*), intent(in)           :: dim_name(:)
    integer,          intent(in)           :: dim_size(:)
    character(len=*), intent(in)           :: dim_desc(:)
    character(len=*), intent(in)           :: dim_units(:)
    character(len=*), intent(in), optional :: dim_type(:)
    character(len=*), intent(in), optional :: default_basename
    real(DP)        , intent(in), optional :: default_tinterval
    character(len=*), intent(in), optional :: default_tunit
    logical,          intent(in), optional :: default_taverage
    character(len=*), intent(in), optional :: default_datatype

    integer :: ndims
    character(len=IO_FILECHR)  :: HISTORY_DEFAULT_BASENAME  = 'history'
    real(DP)                   :: HISTORY_DEFAULT_TINTERVAL = 1.0D0
    character(len=File_HSHORT) :: HISTORY_DEFAULT_TUNIT     = 'sec'
    logical                    :: HISTORY_DEFAULT_TAVERAGE  = .false.
    character(len=File_HSHORT) :: HISTORY_DEFAULT_DATATYPE  = 'REAL4'

    NAMELIST / PARAM_HISTORY / &
         HISTORY_TITLE,             &
         HISTORY_SOURCE,            &
         HISTORY_INSTITUTION,       &
         HISTORY_TIME_UNITS,        &
         HISTORY_DEFAULT_BASENAME,  &
         HISTORY_DEFAULT_TINTERVAL, &
         HISTORY_DEFAULT_TUNIT,     &
         HISTORY_DEFAULT_TAVERAGE,  &
         HISTORY_DEFAULT_DATATYPE

    character(len=IO_FILECHR)  :: BASENAME  !> file base name
    character(len=File_HSHORT) :: ITEM      !> name of history item
    real(8)                    :: TINTERVAL !> time interval to output
    character(len=File_HSHORT) :: TUNIT     !> time unit
    logical                    :: TAVERAGE  !> time average to output
    character(len=File_HSHORT) :: DATATYPE  !> data type

    NAMELIST / HISTITEM / &
       BASENAME,  &
       ITEM,      &
       TINTERVAL, &
       TUNIT,     &
       TAVERAGE,  &
       DATATYPE

    integer :: ierr
    integer :: n
    integer :: arysize, memsize
    intrinsic size
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[HISTORY]/Categ[IO]'

    ndims = size(dim_name)
    if ( size(dim_size)  /= ndims .or. &
         size(dim_desc)  /= ndims .or. &
         size(dim_units) /= ndims      &
         ) then
       write(*,*) "xxx size of dimensions are mismatch"
       call PRC_MPIstop
    end if

    !--- read namelist
    HISTORY_TITLE       = title
    HISTORY_SOURCE      = source
    HISTORY_INSTITUTION = institution
    HISTORY_TIME_UNITS  = 'sec'
    if ( present(default_basename) ) then
       HISTORY_DEFAULT_BASENAME = default_basename
    end if
    if ( present(default_tinterval) ) then
       HISTORY_DEFAULT_TINTERVAL = default_tinterval
       if ( present(default_tunit) ) then
          HISTORY_DEFAULT_TUNIT = default_tunit
       end if
    end if
    if ( present(default_taverage) ) then
       HISTORY_DEFAULT_TAVERAGE = default_taverage
    end if
    if ( present(default_datatype) ) then
       HISTORY_DEFAULT_DATATYPE = default_datatype
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_HISTORY,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_HISTORY. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_HISTORY)

    allocate(History_dim_name (ndims))
    allocate(History_dim_size (ndims))
    allocate(History_dim_desc (ndims))
    allocate(History_dim_units(ndims))
    allocate(History_dim_type (ndims))
    do n = 1, ndims
      History_dim_name(n)  = dim_name(n)
      History_dim_size(n)  = dim_size(n)
      History_dim_desc(n)  = dim_desc(n)
      History_dim_units(n) = dim_units(n)
      History_dim_type(n)  = File_REAL4
      if ( present(dim_type) ) then
         if ( size(dim_type) >= n ) then
            if    ( trim(dim_type(n)) == 'REAL4' ) then
               History_dim_type(n) = File_REAL4
            elseif( trim(dim_type(n)) == 'REAL8' ) then
               History_dim_type(n) = File_REAL8
            else
               write(*,*) 'xxx Not appropriate dim_type. Check!', dim_type(n), n
               call PRC_MPIstop
            endif
         end if
      end if
    end do

    ! listup history request
    rewind(IO_FID_CONF)
    do n = 1, History_req_limit
       read(IO_FID_CONF,nml=HISTITEM,iostat=ierr)
       if( ierr /= 0 ) exit
    enddo
    History_req_nmax = n - 1

    if    ( History_req_nmax > History_req_limit ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** request of history file is exceed! n >', History_req_limit
    elseif( History_req_nmax == 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** No history file specified.'
       return
    endif

    arysize = 1
    do n = 1, ndims
       arysize = arysize * dim_size(n)
    end do

    allocate( History_item   (History_req_nmax) ); History_item(:) = ''
    allocate( History_fid    (History_req_nmax) )
    allocate( History_vid    (History_req_nmax) )
    allocate( History_tintsec(History_req_nmax) )
    allocate( History_tavg   (History_req_nmax) )

    allocate( History_varsum (arysize,History_req_nmax) )
    allocate( History_size   (History_req_nmax) )
    allocate( History_tstrsec(History_req_nmax) )
    allocate( History_tsumsec(History_req_nmax) )

    rewind(IO_FID_CONF)
    memsize = 0
    do n = 1, History_req_nmax
       ! set default
       BASENAME  = HISTORY_DEFAULT_BASENAME
       ITEM      = 'unknown'
       TINTERVAL = HISTORY_DEFAULT_TINTERVAL
       TUNIT     = HISTORY_DEFAULT_TUNIT
       TAVERAGE  = HISTORY_DEFAULT_TAVERAGE
       DATATYPE  = HISTORY_DEFAULT_DATATYPE

       read(IO_FID_CONF,nml=HISTITEM,iostat=ierr)
       if( ierr /= 0 ) exit

       History_req_item(n) = ITEM
       History_req_basename(n) = BASENAME
       call TIME_ymdhms2sec( History_req_tintsec(n), TINTERVAL, TUNIT )
       History_req_tavg(n) = TAVERAGE

       if ( History_req_tintsec(n) <= 0.D0 ) then
          write(*,*) 'xxx Not appropriate time interval. Check!', ITEM, TINTERVAL
          call PRC_MPIstop
       endif

       if    ( trim(DATATYPE) == 'REAL4' ) then
          History_req_dtype(n) = File_REAL4
       elseif( trim(DATATYPE) == 'REAL8' ) then
          History_req_dtype(n) = File_REAL8
       else
          write(*,*) 'xxx Not appropriate DATATYPE. Check!', DATATYPE
          call PRC_MPIstop
       endif

       memsize = memsize + arysize*File_preclist(History_req_dtype(n))
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '*** Number of requested history item             : ', History_req_nmax
    if( IO_L ) write(IO_FID_LOG,*) '*** Output default data type                             : ', HISTORY_DEFAULT_DATATYPE
    if( IO_L ) write(IO_FID_LOG,*) '*** Memory usage for history data buffer [Mbyte] : ', memsize/1024/1024


    return
  end subroutine HistoryInit

  !-----------------------------------------------------------------------------
  subroutine HistoryAddVariable( &
      varname, &
      dims,    &
      desc,    &
      units,   &
      itemid   )
    use gtool_file, only : &
         FileCreate, &
         FileAddVariable
    implicit none

    character(len=*), intent( in) :: varname
    character(len=*), intent( in) :: dims(:)
    character(len=*), intent( in) :: desc
    character(len=*), intent( in) :: units
    integer,          intent(out), optional :: itemid

    integer :: id
    integer :: ary_size
    integer :: nmax, reqid
    integer :: n, m, l
    intrinsic size
    !---------------------------------------------------------------------------

    !--- search existing item
    id = -1
    nmax = min( History_id_count, History_req_nmax )
    do n = 1, nmax
       if ( trim(varname) == trim(History_item(n)) ) then ! match existing item
          if ( present(itemid) ) itemid = n
          return
       endif
    enddo

    if ( id < 0 ) then ! request-register matching check
       do n = 1, History_req_nmax
          if ( trim(varname) == History_req_item(n) ) then
             id = History_id_count
             reqid  = n
             History_id_count = History_id_count + 1

             ! new file registration
             call FileCreate(History_fid(id),                           & ! (out)
                  trim(History_req_basename(reqid)),                    & ! (in)
                  HISTORY_TITLE, HISTORY_SOURCE, HISTORY_INSTITUTION,   & ! (in)
                  History_dim_name, History_dim_size, History_dim_desc, & ! (in)
                  History_dim_units, History_dim_type,                  & ! (in)
                  time_units = HISTORY_TIME_UNITS                       & ! (in)
                  )

             call FileAddVariable(History_vid(id), & ! (out)
                  History_fid(id),                 & ! (in)
                  varname, desc, units, dims,      & ! (in)
                  History_req_dtype  (reqid),      & ! (in)
                  History_req_tintsec(reqid),      & ! (in)
                  History_req_tavg   (reqid)      & ! (in)
                  )

             ary_size = 1
             do m = 1, size(dims)
                do l = 1, size(History_dim_name)
                   if ( trim(dims(m)) == trim(History_dim_name(l)) ) then
                      if ( History_dim_size(l) > 0 ) then
                         ary_size = ary_size * History_dim_size(l)
                         exit
                      end if
                   end if
                end do
             end do

             History_tintsec(id) = History_req_tintsec(reqid)
             History_tavg   (id) = History_req_tavg   (reqid)

             History_varsum(:,id) =  0.D0
             History_size    (id) = ary_size
             History_tstrsec (id) = -1.D0
             History_tsumsec (id) =  0.D0

             if( IO_L ) write(IO_FID_LOG,*) '*** [HIST] Item registration No.= ', id
             if( IO_L ) write(IO_FID_LOG,*) '] Name           : ', trim(History_item(id))
             if( IO_L ) write(IO_FID_LOG,*) '] Description    : ', trim(desc)
             if( IO_L ) write(IO_FID_LOG,*) '] Unit           : ', trim(units)
             if( IO_L ) write(IO_FID_LOG,*) '] size           : ', ary_size
             if( IO_L ) write(IO_FID_LOG,*) '] Interval [sec] : ', History_tintsec(id)
             if( IO_L ) write(IO_FID_LOG,*) '] Average?       : ', History_tavg   (id)
          endif
       enddo
    endif

    if ( present(itemid) ) itemid = id

    return
  end subroutine HistoryAddVariable

  !-----------------------------------------------------------------------------
  subroutine HistoryPutAxis( &
       dim, & ! (in)
       val  & ! (in)
       )
    use gtool_file, only : &
         FilePutAxis
    implicit none

    character(len=*), intent(in) :: dim
    real(DP),         intent(in) :: val(:)

    integer :: m
    logical :: flag = .false.
    intrinsic size
    !---------------------------------------------------------------------------

    do m = 1, size(History_dim_name)
       if ( trim(History_dim_name(m)) == trim(dim) ) then ! dimension is found
          flag = .true.
          exit
       end if
    end do

    if ( .not. flag ) then ! dimension was not found
       write(*,*) "xxx dimension name is invalid: ", dim
       call PRC_MPIstop
    end if

    do m = 1, History_req_nmax
       call FilePutAxis( History_fid(m), dim, val )
    end do

    return
  end subroutine HistoryPutAxis

  !-----------------------------------------------------------------------------
  subroutine HistoryPutAdditionalAxis( &
       name,  & ! (in)
       desc,  & ! (in)
       units, & ! (in)
       dim,   & ! (in)
       var,   & ! (in)
       dtype  & ! (in) optional
       )
    use gtool_file_h, only: &
       File_REAL4, &
       File_REAL8
    use gtool_file, only : &
       FilePutAdditionalAxis
    implicit none

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: dim
    real(DP),         intent(in) :: var(:)
    character(len=*), intent(in), optional :: dtype

    integer :: type
    integer :: m
    intrinsic size
    !---------------------------------------------------------------------------

    if ( present(dtype) ) then
       if ( trim(dtype) == 'REAL4' ) then
          type = File_REAL4
       else if ( trim(dtype) == 'REAL8' ) then
          type = File_REAL8
       else
          write(*,*) 'xxx Not appropriate dtype. Check!', dtype
          call PRC_MPIstop
       end if
    else
       if ( DP == 4 ) then
          type = File_REAL4
       else if ( DP == 8 ) then
          type = File_REAL8
       end if
    end if

    do m = 1, size(History_dim_name)
       if ( trim(History_dim_name(m)) == trim(dim) ) then
          if ( History_dim_size(m) /= size(var) ) then
             write(*,*) 'xxx size of var is not match to dim. Check!', size(var), History_dim_size(m), '(', dim, ')'
             call PRC_MPIstop
          end if
          exit
       end if
    end do

    do m = 1, History_req_nmax
       call FilePutAdditionalAxis( History_fid(m),     & ! (in)
            name, desc, units, dim, type, var, size(var) ) ! (in)
    end do

    return
  end subroutine HistoryPutAdditionalAxis

  !-----------------------------------------------------------------------------
  ! interface HistoryPut
  !-----------------------------------------------------------------------------
  subroutine HistoryPut1D( &
       varname, &
       var,     &
       time,    &
       dt,      &
       itemid,  &
       force    )
    implicit none

    character(len=*), intent(in)           :: varname
    real(DP),         intent(in)           :: var(:)
    real(DP),         intent(in)           :: time
    real(DP),         intent(in)           :: dt
    integer,          intent(in), optional :: itemid
    logical,          intent(in), optional :: force

    integer :: id, n
    !---------------------------------------------------------------------------

    if ( present(itemid) ) then
       id = itemid
    else
       ! search item id
       id = -1
       do n = 1, History_req_nmax
          if ( trim(varname) == trim(History_item(n)) ) then
             id = n
             exit
          end if
       end do
    end if

    if ( id < 0 ) return

    if ( History_tavg(id) ) then
       History_varsum(1:History_size(id),id) = &
            History_varsum(1:History_size(id),id) &
            + var(:) * dt
    else
       History_varsum(1:History_size(id),id) = var(:)
    endif
    History_tsumsec(id) = History_tsumsec(id) + dt

    if ( ( present(force) .and. force ) .or. &
         History_tsumsec(id) - History_tintsec(id) > -eps ) then
       call HistoryWrite(id, time)
    end if

    return
  end subroutine HistoryPut1D
  subroutine HistoryPut2D( &
       varname, &
       var,     &
       time,    &
       dt,      &
       itemid,  &
       force    )
    implicit none

    character(len=*), intent(in)           :: varname
    real(DP),         intent(in)           :: var(:,:)
    real(DP),         intent(in)           :: time
    real(DP),         intent(in)           :: dt
    integer,          intent(in), optional :: itemid
    logical,          intent(in), optional :: force

    intrinsic size
    !---------------------------------------------------------------------------

    call HistoryPut1D(varname, &
         reshape(var(:,:),(/size(var)/)), &
         time, dt, itemid, force )

    return
  end subroutine HistoryPut2D
  subroutine HistoryPut3D( &
       varname, &
       var,     &
       time,    &
       dt,      &
       itemid,  &
       force    )
    implicit none

    character(len=*), intent(in)           :: varname
    real(DP),         intent(in)           :: var(:,:,:)
    real(DP),         intent(in)           :: time
    real(DP),         intent(in)           :: dt
    integer,          intent(in), optional :: itemid
    logical,          intent(in), optional :: force

    intrinsic size
    !---------------------------------------------------------------------------

    call HistoryPut1D(varname, &
         reshape(var(:,:,:),(/size(var)/)), &
         time, dt, itemid, force )

    return
  end subroutine HistoryPut3D

  !-----------------------------------------------------------------------------
  subroutine HistoryWriteAll( &
       time & ! (in)
       )
    implicit none

    real(DP), intent(in) :: time

    integer :: n

    do n = 1, History_req_nmax
       call HistoryWrite( n, time )
    end do

    return
  end subroutine HistoryWriteAll

  !-----------------------------------------------------------------------------
  ! interface HistoryGet
  !-----------------------------------------------------------------------------
  subroutine HistoryGet1D( &
       var, &
       basename, &
       varname, &
       step, &
       single )
    use gtool_file, only : &
         FileRead
    implicit none

    real(DP),         intent(out) :: var(:)
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: varname
    integer,          intent( in) :: step
    logical,          intent( in), optional :: single

    call FileRead(var,                   & ! (out)
         basename, varname, step, single & ! (in)
         )

  end subroutine HistoryGet1D
  subroutine HistoryGet2D( &
       var, &
       basename, &
       varname, &
       step, &
       allow_missing )
    use gtool_file, only : &
         FileRead
    implicit none

    real(DP),         intent(out) :: var(:,:)
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: varname
    integer,          intent( in) :: step
    logical,          intent( in), optional :: allow_missing
    !---------------------------------------------------------------------------

    call FileRead(var,                                   & ! (out)
         basename, varname, step, allow_missing, .false. & ! (in)
         )

    return
  end subroutine HistoryGet2D
  subroutine HistoryGet3D( &
       var, &
       basename, &
       varname, &
       step, &
       allow_missing )
    use gtool_file, only : &
         FileRead
    implicit none

    real(DP),         intent(out) :: var(:,:,:)
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: varname
    integer,          intent( in) :: step
    logical,          intent( in), optional :: allow_missing
    !---------------------------------------------------------------------------

    call FileRead(var,                                   & ! (out)
         basename, varname, step, allow_missing, .false. & ! (in)
         )

    return
  end subroutine HistoryGet3D

  !-----------------------------------------------------------------------------
  subroutine HistoryOutputList
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** [HIST] Output item list '
    if( IO_L ) write(IO_FID_LOG,*) '*** Number of history item :', History_req_nmax
    if( IO_L ) write(IO_FID_LOG,*) 'NAME           :size         :interval[sec]:avg'
    if( IO_L ) write(IO_FID_LOG,*) '============================================================================'
 
    do n = 1, History_id_count-1
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I10,1x,f13.3,1x,L)') History_item(n), History_size(n), History_tintsec(n), History_tavg(n)
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '============================================================================'

    return
  end subroutine HistoryOutputList

  !-----------------------------------------------------------------------------
  subroutine HistoryFinalize
    use gtool_file, only : &
       FileClose
    implicit none

    integer :: n
    !---------------------------------------------------------------------------

    do n = 1, History_req_nmax
       call FileClose( History_fid(n) )
    end do

    return
  end subroutine HistoryFinalize

  !-----------------------------------------------------------------------------
  ! private
  !-----------------------------------------------------------------------------
  subroutine HistoryWrite( &
       itemid, &
       time )
    use mod_time, only: &
       TIME_sec2ymdhms
    use gtool_file, only: &
       FileWrite
    implicit none

    integer, intent(in) :: itemid
    real(8), intent(in) :: time

    real(8) :: stime, etime
    logical, save :: firsttime = .true.
    !---------------------------------------------------------------------------

    if( History_id_count == 1 ) return

    if (firsttime) then
       firsttime = .false.
       call HistoryOutputList
    endif

    if ( History_tsumsec(itemid) > eps ) then
       if ( History_tavg(itemid) ) then
          History_varsum(1:History_size(itemid),itemid) = &
               History_varsum(1:History_size(itemid),itemid) / History_tsumsec(itemid)
       end if

       if ( History_tstrsec(itemid) < 0.D0 ) then ! first time
          stime = time - History_tsumsec(itemid)
       else
          stime = History_tstrsec(itemid)
       end if
       etime = History_tstrsec(itemid) + History_tsumsec(itemid) ! neary equal to time

       ! convert time units
       call TIME_sec2ymdhms( stime, stime, HISTORY_TIME_UNITS )
       call TIME_sec2ymdhms( etime, etime, HISTORY_TIME_UNITS )
       
       call FileWrite( History_vid(itemid),                & ! vid
            History_varsum(1:History_size(itemid),itemid), & ! data
            stime, etime                                   ) ! start & end time

    endif

    History_varsum(:,itemid) = 0.0_DP
    History_tstrsec (itemid) = time
    History_tsumsec (itemid) = 0.D0

    return
  end subroutine HistoryWrite

end module gtool_history
!-------------------------------------------------------------------------------
