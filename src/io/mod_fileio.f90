!-------------------------------------------------------------------------------
!> module FILE I/O
!!
!! @par Description
!!          Fortran-C interface for file I/O
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-FIO
!!
!<
!-------------------------------------------------------------------------------
module mod_fileio
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_fileio_h
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
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,             parameter, private :: FIO_nmaxfile = 64
  character(LEN=FIO_HLONG), save, private :: FIO_fname_list(FIO_nmaxfile)
  integer,                  save, private :: FIO_fid_list  (FIO_nmaxfile)
  integer,                  save, private :: FIO_fid_count = 1

  type(headerinfo), private :: hinfo 
  type(datainfo),   private :: dinfo 

  integer, private, parameter :: max_num_of_data = 2500 !--- max time step num
  integer, parameter, private :: preclist(0:3) = (/ 4, 8, 4, 8 /)

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  subroutine FIO_setup
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_myrank
    use mod_grid, only: &
       GRID_DX,  &
       GRID_IMAX
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[FILEIO]/Categ[IO]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Maximum limit for file registration   : ', FIO_nmaxfile
    if( IO_L ) write(IO_FID_LOG,*) '*** Maximum limit for timestep in one file: ', max_num_of_data

    call fio_syscheck()
    call fio_put_commoninfo( FIO_MPIIO_NOUSE, &
                             FIO_SPLIT_FILE,  &
                             FIO_BIG_ENDIAN,  &
                             FIO_CARTESIAN,   &
                             int(GRID_DX),    &
                             GRID_IMAX,       &
                             1,               &
                             PRC_myrank       )

    allocate( hinfo%rgnid(1) )

    return
  end subroutine FIO_setup

  !-----------------------------------------------------------------------------
  subroutine FIO_getfid( &
      fid,      &
      basename, &
      rwtype,   &
      pkg_desc, &
      pkg_note  )
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_myrank
    implicit none

    integer,          intent(out) :: fid
    character(LEN=*), intent( in) :: basename
    integer,          intent( in) :: rwtype
    character(LEN=*), intent( in) :: pkg_desc
    character(LEN=*), intent( in) :: pkg_note

    character(LEN=FIO_HSHORT) :: rwname(0:2)
    data rwname / 'READ','WRITE','APPEND' / ! [fix] H.Yashiro 20110912

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
       call fio_mk_fname(fname,trim(basename),'pe',PRC_myrank,6)
       call fio_register_file(n,fname)

       if ( rwtype == FIO_FREAD ) then

!          call fio_dump_finfo(n,FIO_BIG_ENDIAN,FIO_DUMP_HEADER) ! dump to stdout(check)
          call fio_fopen(n,FIO_FREAD)
          call fio_read_allinfo(n)

       elseif( rwtype == FIO_FWRITE ) then

          call fio_fopen(n,FIO_FWRITE)
          call fio_put_write_pkginfo(n,pkg_desc,pkg_note)

       endif

       if( IO_L ) write(IO_FID_LOG,*) '*** [FIO] File registration : ',trim(rwname(rwtype)),' -',n
       if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(fname)

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
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid, only: &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX
    implicit none

    real(8),          intent(out) :: var(:,:,:)
    character(LEN=*), intent( in) :: basename
    character(LEN=*), intent( in) :: varname
    character(LEN=*), intent( in) :: layername
    integer,          intent( in) :: k_start, k_end
    integer,          intent( in) :: step

    logical, intent(in), optional :: allow_missingq !--- if data is missing, set value to zero

    real(4) :: var4(IMAX,JMAX,k_start:k_end)
    real(8) :: var8(IMAX,JMAX,k_start:k_end)

    integer :: did, fid
    !---------------------------------------------------------------------------

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FREAD, "", "" )

    !--- seek data ID and get information
    call fio_seek_datainfo(did,fid,varname,step)
    call fio_get_datainfo(fid,did,dinfo)

!    if( IO_L ) write(IO_FID_LOG,*) dinfo%varname
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%description
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%unit
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%layername
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%note
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%datasize
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%datatype
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%num_of_layer
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%step
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%time_start
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%time_end

    !--- verify
    if ( did == -1 ) then
       if ( present(allow_missingq) ) then ! [bugfix] H.Yashiro 20110912
          if ( allow_missingq ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                                            'varname= ',trim(varname),', step=',step
             if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] Q Value is set to 0.'

             var(:,:,k_start:k_end) = 0.D0

             return
          endif
       else
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                                         'varname= ',trim(varname),', step=',step
          call PRC_MPIstop
       endif
    endif

    if ( trim(dinfo%layername) /= trim(layername) ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                                      "[",trim(dinfo%layername),":",trim(layername),"]"
       call PRC_MPIstop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch! ', &
                                      dinfo%num_of_layer,k_end-k_start+1
       call PRC_MPIstop
    endif

    !--- read data
    if ( dinfo%datatype == FIO_REAL4 ) then

       call fio_read_data(fid,did,var4(:,:,:))
       var(:,:,k_start:k_end) = real(var4(:,:,1:dinfo%num_of_layer),kind=8)

    elseif( dinfo%datatype == FIO_REAL8 ) then

       call fio_read_data(fid,did,var8(:,:,:))
       var(:,:,k_start:k_end) = var8(:,:,1:dinfo%num_of_layer)

    endif

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
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only :&
      TIME_date2sec, &
      TIME_sec2date
    implicit none

    integer,          intent(inout) :: start_step
    integer,          intent(inout) :: num_of_step
    integer,          intent(inout) :: data_date(6,max_num_of_data)
    integer,          intent(inout) :: prec
    character(LEN=*), intent(   in) :: basename
    character(LEN=*), intent(   in) :: varname
    character(LEN=*), intent(   in) :: layername
    integer,          intent(   in) :: k_start, k_end
    real(8),          intent(   in) :: ctime
    integer,          intent(   in) :: cdate(6)
    logical,          intent(   in) :: opt_periodic_year

    real(8) :: midtime  !--- [sec]
    real(8) :: microsec !--- [sec]
    logical :: startflag
    integer :: did, fid
    integer :: i
    !---------------------------------------------------------------------------

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FREAD, "", "" )

    startflag = .false.

    do i = 1, max_num_of_data
       !--- seek data ID and get information
       call fio_seek_datainfo(did,fid,varname,i)
       call fio_get_datainfo(fid,did,dinfo)

!       if( IO_L ) write(IO_FID_LOG,*) dinfo%varname
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%description
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%unit
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%layername
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%note
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%datasize
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%datatype
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%num_of_layer
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%step
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%time_start
!       if( IO_L ) write(IO_FID_LOG,*) dinfo%time_end

       if ( did == -1 ) then
          num_of_step = i - 1
          exit
       endif

       !--- verify
       if ( trim(dinfo%layername) /= trim(layername) ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                               "[",trim(dinfo%layername),":",trim(layername),"]"
          call PRC_MPIstop
       elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch!', &
                               dinfo%num_of_layer,k_end-k_start+1
          call PRC_MPIstop
       endif

       ! [fix] H.Yashiro 20111011 : specify int kind=8
       midtime = dble( int( (dinfo%time_start+dinfo%time_end)*0.5D0+1.D0,kind=8 ) ) * 1.D-3
       call TIME_sec2date( data_date(:,i), microsec, midtime )

       if ( opt_periodic_year ) then
          data_date(1,i) = cdate(1)
          call TIME_date2sec( midtime, data_date(:,i), microsec )
       endif

       if (       ( .not. startflag ) &
            .AND. ( ctime < midtime ) ) then
          startflag  = .true.
          start_step = i
          prec       = preclist(dinfo%datatype)
       endif
    enddo

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
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF4
    use mod_grid, only: &
       IMAX => GRID_IMAX, &
       JMAX => GRID_JMAX
    implicit none

    real(8),          intent(in) :: var(:,:,:)
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
    real(8),          intent(in) :: t_start, t_end

    real(4) :: var4(IMAX,JMAX,k_start:k_end)
    real(8) :: var8(IMAX,JMAX,k_start:k_end)

    integer :: did, fid
    !---------------------------------------------------------------------------

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    dinfo%varname      = varname
    dinfo%description  = data_desc
    dinfo%unit         = unit
    dinfo%layername    = layername
    dinfo%note         = data_note
    dinfo%datasize     = IMAX * JMAX * (k_end-k_start+1) * preclist(dtype)
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start * 1.D3, kind=8 )
    dinfo%time_end     = int( t_end   * 1.D3, kind=8 )

!    if( IO_L ) write(IO_FID_LOG,*) dinfo%varname
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%description
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%unit
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%layername
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%note
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%datasize
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%datatype
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%num_of_layer
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%step
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%time_start
!    if( IO_L ) write(IO_FID_LOG,*) dinfo%time_end

    if ( dtype == FIO_REAL4 ) then

       var4(:,k_start:k_end,:)=real(var(:,:,k_start:k_end),kind=4)
       where( var4(:,:,:) < (CONST_UNDEF4+1.0) )
          var4(:,:,:) = CONST_UNDEF4
       endwhere

       call fio_put_write_datainfo_data(did,fid,dinfo,var4(:,:,:))

    elseif( dtype == FIO_REAL8 ) then

       var8(:,:,k_start:k_end)=var(:,:,k_start:k_end)

       call fio_put_write_datainfo_data(did,fid,dinfo,var8(:,:,:))
    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call PRC_MPIstop
    endif

    return
  end subroutine FIO_output

  !-----------------------------------------------------------------------------
  subroutine FIO_finalize
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_myrank
    implicit none

    character(LEN=FIO_HLONG) :: fname
    integer                  :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    do n = 1, FIO_fid_count-1
       call fio_fclose(FIO_fid_list(n))

       if( IO_L ) write(IO_FID_LOG,*) '*** [FIO] File Close : NO.', n
       call fio_mk_fname(fname,trim(FIO_fname_list(n)),'pe',PRC_myrank,6)
       if( IO_L ) write(IO_FID_LOG,*) '*** closed filename: ', trim(fname)
    enddo

    return
  end subroutine FIO_finalize

end module mod_fileio
!-------------------------------------------------------------------------------
