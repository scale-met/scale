!-------------------------------------------------------------------------------
!> Module file I/O
!!
!! @par Description
!!         In this module, the geometrics of the icosahedral grid such as area are calculated
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_fio
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_prof
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FIO_setup
  public :: FIO_input
  public :: FIO_seek
  public :: FIO_output
  public :: FIO_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- character length
  integer, public, parameter :: FIO_HSHORT =  16
  integer, public, parameter :: FIO_HMID   =  64
  integer, public, parameter :: FIO_HLONG  = 256

  !--- data type
  integer, public, parameter :: FIO_REAL4    = 0
  integer, public, parameter :: FIO_REAL8    = 1
  integer, public, parameter :: FIO_INTEGER4 = 2
  integer, public, parameter :: FIO_INTEGER8 = 3

  !--- data endian
  integer, public, parameter :: FIO_UNKNOWN_ENDIAN = 0
  integer, public, parameter :: FIO_LITTLE_ENDIAN  = 1
  integer, public, parameter :: FIO_BIG_ENDIAN     = 2

  !--- topology
  integer, public, parameter :: FIO_ICOSAHEDRON = 0
  integer, public, parameter :: FIO_IGA_LCP     = 1
  integer, public, parameter :: FIO_IGA_MLCP    = 2

  !--- file mode (partial or complete)
  integer, public, parameter :: FIO_SPLIT_FILE = 0
  integer, public, parameter :: FIO_INTEG_FILE = 1

  !--- proccessor type
  integer, public, parameter :: FIO_SINGLE_PROC = 0
  integer, public, parameter :: FIO_MULTI_PROC  = 1

  !--- action type
  integer, public, parameter :: FIO_FREAD   = 0
  integer, public, parameter :: FIO_FWRITE  = 1
  integer, public, parameter :: FIO_FAPPEND = 2

  !--- data dump type
  integer, public, parameter :: FIO_DUMP_OFF      = 0
  integer, public, parameter :: FIO_DUMP_HEADER   = 1
  integer, public, parameter :: FIO_DUMP_ALL      = 2
  integer, public, parameter :: FIO_DUMP_ALL_MORE = 3

  !--- struct for package infomation
  type, public :: headerinfo
     character(LEN=FIO_HLONG) :: fname
     character(LEN=FIO_HMID)  :: description
     character(LEN=FIO_HLONG) :: note
     integer                  :: num_of_data
     integer                  :: fmode
     integer                  :: endiantype
     integer                  :: grid_topology
     integer                  :: glevel
     integer                  :: rlevel
     integer                  :: num_of_rgn
     integer, pointer         :: rgnid(:)
  endtype headerinfo

  !--- struct for data infomation
  type, public :: datainfo
     character(LEN=FIO_HSHORT) :: varname
     character(LEN=FIO_HMID)   :: description
     character(LEN=FIO_HSHORT) :: unit
     character(LEN=FIO_HSHORT) :: layername
     character(LEN=FIO_HLONG)  :: note
     integer(8)                :: datasize
     integer                   :: datatype
     integer                   :: num_of_layer
     integer                   :: step
     integer(8)                :: time_start
     integer(8)                :: time_end
  endtype datainfo

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                  private, parameter :: FIO_nmaxfile = 64
  character(LEN=FIO_HLONG), private            :: FIO_fname_list(FIO_nmaxfile)
  integer,                  private            :: FIO_fid_list  (FIO_nmaxfile)
  integer,                  private            :: FIO_fid_count = 1

  type(headerinfo), private :: hinfo
  type(datainfo),   private :: dinfo

  integer, private, parameter :: max_num_of_data = 2500 !--- max time step num
  integer, private, parameter :: preclist(0:3) = (/ 4, 8, 4, 8 /)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine FIO_setup
    use mod_adm, only : &
       ADM_prc_tab, &
       ADM_prc_me,  &
       ADM_glevel,  &
       ADM_rlevel,  &
       ADM_lall
    implicit none

    integer, allocatable :: prc_tab(:)
    !---------------------------------------------------------------------------

    ! dummy call
    call PROF_rapstart('FILEIO_in')
    call PROF_rapend  ('FILEIO_in')
    call PROF_rapstart('FILEIO_out')
    call PROF_rapend  ('FILEIO_out')

    allocate( prc_tab(ADM_lall) )
    prc_tab(1:ADM_lall) = ADM_prc_tab(1:ADM_lall,ADM_prc_me)-1

    call fio_syscheck()
    call fio_put_commoninfo( FIO_SPLIT_FILE,  &
                             FIO_BIG_ENDIAN,  &
                             FIO_ICOSAHEDRON, &
                             ADM_glevel,      &
                             ADM_rlevel,      &
                             ADM_lall,        &
                             prc_tab          )

    deallocate(prc_tab)

    allocate( hinfo%rgnid(ADM_lall) )

    return
  end subroutine FIO_setup

  !-----------------------------------------------------------------------------
  subroutine FIO_getfid( &
       fid,      &
       basename, &
       rwtype,   &
       pkg_desc, &
       pkg_note  )
    use mod_adm, only : &
       ADM_prc_me
    implicit none

    integer,          intent(out) :: fid
    character(LEN=*), intent( in) :: basename
    integer,          intent( in) :: rwtype
    character(LEN=*), intent( in) :: pkg_desc
    character(LEN=*), intent( in) :: pkg_note

    character(LEN=FIO_HSHORT) :: rwname(0:2)
    data rwname / 'READ','WRITE','APPEND' /

    character(LEN=FIO_HLONG) :: fname
    integer                  :: n
    !---------------------------------------------------------------------------

    !--- search existing file
    fid = -1
    do n = 1, FIO_fid_count
       if ( trim(basename)==trim(FIO_fname_list(n)) ) fid = FIO_fid_list(n)
    enddo

    if ( fid < 0 ) then ! file registration
       !--- register new file and open
       call fio_mk_fname(fname,trim(basename),'pe',ADM_prc_me-1,6)
       call fio_register_file(n,fname)

       if ( rwtype == FIO_FREAD ) then

!          call fio_dump_finfo(n,FIO_BIG_ENDIAN,FIO_DUMP_HEADER) ! dump to stdout(check)
          call fio_fopen(n,FIO_FREAD)
          call fio_read_allinfo(n)

       elseif( rwtype == FIO_FWRITE ) then

          call fio_fopen(n,FIO_FWRITE)
          call fio_put_write_pkginfo(n,pkg_desc,pkg_note)

       endif

       write(ADM_LOG_FID,*) '*** [FIO] File registration : ',trim(rwname(rwtype)),'-', n
       write(ADM_LOG_FID,*) '*** filename: ', trim(fname)

       FIO_fname_list(FIO_fid_count) = trim(basename)
       FIO_fid_list  (FIO_fid_count) = n
       FIO_fid_count = FIO_fid_count + 1
       fid = n
    endif

    return
  end subroutine FIO_getfid

  !-----------------------------------------------------------------------------
  subroutine FIO_input( &
       var,           &
       basename,      &
       varname,       &
       layername,     &
       k_start,       &
       k_end,         &
       step,          &
       allow_missingq ) !--- optional
    use mod_adm, only : &
       ADM_proc_stop, &
       ADM_gall,      &
       ADM_lall
    implicit none

    real(RP),         intent(out) :: var(:,:,:)
    character(LEN=*), intent( in) :: basename
    character(LEN=*), intent( in) :: varname
    character(LEN=*), intent( in) :: layername
    integer,          intent( in) :: k_start, k_end
    integer,          intent( in) :: step

    logical, intent(in), optional :: allow_missingq !--- if data is missing, set value to zero

    real(SP) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in')

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FREAD, "", "" )

    !--- seek data ID and get information
    call fio_seek_datainfo(did,fid,varname,step)
    call fio_get_datainfo(fid,did,dinfo)

!    write(ADM_LOG_FID,*) dinfo%varname
!    write(ADM_LOG_FID,*) dinfo%description
!    write(ADM_LOG_FID,*) dinfo%unit
!    write(ADM_LOG_FID,*) dinfo%layername
!    write(ADM_LOG_FID,*) dinfo%note
!    write(ADM_LOG_FID,*) dinfo%datasize
!    write(ADM_LOG_FID,*) dinfo%datatype
!    write(ADM_LOG_FID,*) dinfo%num_of_layer
!    write(ADM_LOG_FID,*) dinfo%step
!    write(ADM_LOG_FID,*) dinfo%time_start
!    write(ADM_LOG_FID,*) dinfo%time_end

    !--- verify
    if ( did == -1 ) then
       if ( present(allow_missingq) ) then
          if ( allow_missingq ) then
             write(ADM_LOG_FID,*) '*** [INPUT]/[FIO] data not found! : ', &
                                  'varname= ',trim(varname),', step=',step
             write(ADM_LOG_FID,*) '*** [INPUT]/[FIO] Q Value is set to 0.'

             var(:,k_start:k_end,:) = 0.0_RP

             call PROF_rapend('FILEIO_in')
             return
          endif
       endif

       write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                            'varname= ',trim(varname),', step=',step
       call ADM_proc_stop
    endif

    if ( trim(dinfo%layername) /= trim(layername) ) then
       write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                            "[",trim(dinfo%layername),":",trim(layername),"]"
       call ADM_proc_stop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch! ', &
                            dinfo%num_of_layer,k_end-k_start+1
       call ADM_proc_stop
    endif

    !--- read data
    if ( dinfo%datatype == FIO_REAL4 ) then

       call fio_read_data(fid,did,var4(:,:,:))
       var(:,k_start:k_end,:) = real(var4(:,1:dinfo%num_of_layer,:),kind=RP)

    elseif( dinfo%datatype == FIO_REAL8 ) then

       call fio_read_data(fid,did,var8(:,:,:))
       var(:,k_start:k_end,:) = real(var8(:,1:dinfo%num_of_layer,:),kind=RP)

    endif

    call PROF_rapend('FILEIO_in')

    return
  end subroutine FIO_input

  !-----------------------------------------------------------------------------
  subroutine FIO_seek( &
       start_step,       &
       num_of_step,      &
       data_date,        &
       prec,             &
       basename,         &
       varname,          &
       layername,        &
       k_start,          &
       k_end,            &
       ctime,            &
       cdate,            &
       opt_periodic_year )
    use mod_adm, only: &
       ADM_proc_stop
    use scale_calendar, only: &
       CALENDAR_daysec2date,    &
       CALENDAR_date2daysec,    &
       CALENDAR_adjust_daysec,  &
       CALENDAR_combine_daysec
    implicit none

    integer,          intent(inout) :: start_step
    integer,          intent(inout) :: num_of_step
    integer,          intent(inout) :: data_date(6,max_num_of_data)
    integer,          intent(inout) :: prec
    character(len=*), intent(in)    :: basename
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: layername
    integer,          intent(in)    :: k_start, k_end
    real(RP),         intent(in)    :: ctime
    integer,          intent(in)    :: cdate(6)
    logical,          intent(in)    :: opt_periodic_year

    real(RP) :: midtime !--- [sec]

    integer  :: midday
    real(DP) :: midsec
    real(DP) :: midms
    integer  :: offset_year

    logical :: startflag
    integer :: did, fid
    integer :: i
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in')

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FREAD, "", "" )

    startflag = .false.

    do i = 1, max_num_of_data
       !--- seek data ID and get information
       call fio_seek_datainfo(did,fid,varname,i)
       call fio_get_datainfo(fid,did,dinfo)

!       write(ADM_LOG_FID,*) dinfo%varname
!       write(ADM_LOG_FID,*) dinfo%description
!       write(ADM_LOG_FID,*) dinfo%unit
!       write(ADM_LOG_FID,*) dinfo%layername
!       write(ADM_LOG_FID,*) dinfo%note
!       write(ADM_LOG_FID,*) dinfo%datasize
!       write(ADM_LOG_FID,*) dinfo%datatype
!       write(ADM_LOG_FID,*) dinfo%num_of_layer
!       write(ADM_LOG_FID,*) dinfo%step
!       write(ADM_LOG_FID,*) dinfo%time_start
!       write(ADM_LOG_FID,*) dinfo%time_end

       if ( did == -1 ) then
          num_of_step = i - 1
          exit
       endif

       !--- verify
       if ( trim(dinfo%layername) /= trim(layername) ) then
          write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                               "[",trim(dinfo%layername),":",trim(layername),"]"
          call ADM_proc_stop
       elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
          write(ADM_LOG_FID,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch!', &
                               dinfo%num_of_layer,k_end-k_start+1
          call ADM_proc_stop
       endif

       ! [fix] H.Yashiro 20111011 : specify int kind=8
       midtime = real( int( (dinfo%time_start+dinfo%time_end)*0.5_RP+1.0_RP,kind=8 ),kind=RP )

       midday      = 0
       midsec      = midtime
       offset_year = 0

       call CALENDAR_adjust_daysec( midday, midsec ) ! [INOUT]

       call CALENDAR_daysec2date( data_date(:,i), & ! [OUT]
                                  midms,          & ! [OUT]
                                  midday,         & ! [IN]
                                  midsec,         & ! [IN]
                                  offset_year     ) ! [IN]

       if ( opt_periodic_year ) then
          data_date(1,i) = cdate(1)

          call CALENDAR_date2daysec( midday,         & ! [OUT]
                                     midsec,         & ! [OUT]
                                     data_date(:,i), & ! [IN]
                                     midms,          & ! [IN]
                                     offset_year     ) ! [IN]

          midtime = CALENDAR_combine_daysec( midday, midsec )
       endif

       if (       ( .not. startflag ) &
            .AND. ( ctime < midtime ) ) then
          startflag  = .true.
          start_step = i
          prec       = preclist(dinfo%datatype)
       endif
    enddo

    call PROF_rapend('FILEIO_in')

    return
  end subroutine FIO_seek

  !-----------------------------------------------------------------------------
  subroutine FIO_output( &
       var,       &
       basename,  &
       pkg_desc,  &
       pkg_note,  &
       varname,   &
       data_desc, &
       data_note, &
       unit,      &
       dtype,     &
       layername, &
       k_start,   &
       k_end,     &
       step,      &
       t_start,   &
       t_end      )
    use mod_adm, only: &
       ADM_proc_stop, &
       ADM_gall,      &
       ADM_lall
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4
    implicit none

    real(RP),          intent(in) :: var(:,:,:)
    character(LEN=*), intent(in) :: basename
    character(LEN=*), intent(in) :: pkg_desc
    character(LEN=*), intent(in) :: pkg_note
    character(LEN=*), intent(in) :: varname
    character(LEN=*), intent(in) :: data_desc
    character(LEN=*), intent(in) :: data_note
    character(LEN=*), intent(in) :: unit
    integer,          intent(in) :: dtype
    character(LEN=*), intent(in) :: layername
    integer,          intent(in) :: k_start, k_end
    integer,          intent(in) :: step
    real(RP),          intent(in) :: t_start, t_end

    real(SP) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_out')

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    dinfo%varname      = varname
    dinfo%description  = data_desc
    dinfo%unit         = unit
    dinfo%layername    = layername
    dinfo%note         = data_note
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * preclist(dtype), kind=8 )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start, kind=8 )
    dinfo%time_end     = int( t_end,   kind=8 )

!    write(ADM_LOG_FID,*) dinfo%varname
!    write(ADM_LOG_FID,*) dinfo%description
!    write(ADM_LOG_FID,*) dinfo%unit
!    write(ADM_LOG_FID,*) dinfo%layername
!    write(ADM_LOG_FID,*) dinfo%note
!    write(ADM_LOG_FID,*) dinfo%datasize
!    write(ADM_LOG_FID,*) dinfo%datatype
!    write(ADM_LOG_FID,*) dinfo%num_of_layer
!    write(ADM_LOG_FID,*) dinfo%step
!    write(ADM_LOG_FID,*) dinfo%time_start
!    write(ADM_LOG_FID,*) dinfo%time_end

    if ( dtype == FIO_REAL4 ) then

       var4(:,k_start:k_end,:)=real(var(:,k_start:k_end,:),kind=SP)
       where( var4(:,:,:) < (UNDEF4+1.0) )
          var4(:,:,:) = UNDEF4
       endwhere

       call fio_put_write_datainfo_data(did,fid,dinfo,var4(:,:,:))

    elseif( dtype == FIO_REAL8 ) then

       var8(:,k_start:k_end,:)=real(var(:,k_start:k_end,:),kind=DP)

       call fio_put_write_datainfo_data(did,fid,dinfo,var8(:,:,:))
    else
       write(ADM_LOG_FID,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call ADM_proc_stop
    endif

    call PROF_rapend('FILEIO_out')

    return
  end subroutine FIO_output

  !-----------------------------------------------------------------------------
  subroutine FIO_finalize
    use mod_adm, only : &
       ADM_prc_me
    implicit none

    character(LEN=FIO_HLONG) :: fname
    integer                  :: n
    !---------------------------------------------------------------------------

    do n = 1, FIO_fid_count
       call fio_fclose(FIO_fid_list(n))

       write(ADM_LOG_FID,*) '*** [FIO] File Close : NO.', n
       call fio_mk_fname(fname,trim(FIO_fname_list(n)),'pe',ADM_prc_me-1,6)
       write(ADM_LOG_FID,*) '*** closed filename: ', trim(fname)
    enddo

    return
  end subroutine FIO_finalize

end module mod_fio
