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
  use mod_stdio, only: &
     IO_FILECHR,  &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  use mod_fileio_h
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FIO_setup
  public :: FIO_setgridinfo
  public :: FIO_getfid
  public :: FIO_input
  public :: FIO_seek
  public :: FIO_output
  public :: FIO_finalize
  public :: FIO_input_1D
  public :: FIO_output_1D
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
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
#ifdef CONFIG_HDF5
  character(len=IO_FILECHR), private,      save :: GROUP_OUT_BASENAME = 'time_step'
#endif
  integer,                   private, parameter :: FIO_nfile_max = 64 ! number limit of file step
  integer,                   private, parameter :: FIO_nstep_max = 2500 ! number limit of time step

  character(LEN=FIO_HLONG),  private,      save :: FIO_fname_list(FIO_nfile_max)
  integer,                   private,      save :: FIO_fid_list  (FIO_nfile_max)
  integer,                   private,      save :: FIO_fid_count = 1

  type(headerinfo), private :: hinfo 
  type(datainfo),   private :: dinfo 

  integer, private, save :: FIO_IMAX
  integer, private, save :: FIO_JMAX

  !-----------------------------------------------------------------------------  
contains
  !-----------------------------------------------------------------------------
  subroutine FIO_setup
    use mod_process, only: &
       PRC_myrank
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[FILEIO]/Categ[IO]'
#ifdef CONFIG_HDF5
    if( IO_L ) write(IO_FID_LOG,*) '+++ HDF5 format / SZIP compression'
#endif
    if( IO_L ) write(IO_FID_LOG,*) '*** Maximum limit for file registration    : ', FIO_nfile_max
    if( IO_L ) write(IO_FID_LOG,*) '*** Maximum limit for timestep in one file : ', FIO_nstep_max

    ! only for register
    call TIME_rapstart('FILE I')
    call TIME_rapend  ('FILE I')
    call TIME_rapstart('FILE O')
    call TIME_rapend  ('FILE O')

    call fio_syscheck

    call fio_put_commoninfo( FIO_MPIIO_NOUSE,     &
                             FIO_SPLIT_FILE,      &
                             FIO_BIG_ENDIAN,      &
                             FIO_CARTESIAN,       &
                             0,                   & !--- will update after setgridinfo
                             1,                   & !--- will update after setgridinfo
#ifdef CONFIG_HDF5
                             1,                   & !--- will update after setgridinfo
                             1,                   & !--- will update after setgridinfo
#endif
                             1,                   &
                             PRC_myrank           )

    allocate( hinfo%rgnid(1) )

    return
  end subroutine FIO_setup

  !-----------------------------------------------------------------------------
  subroutine FIO_setgridinfo( &
      resolution, &
      imax,       &
      jmax        )
    use mod_process, only: &
       PRC_myrank
    implicit none

    real(RP), intent(in) :: resolution
    integer, intent(in) :: imax
    integer, intent(in) :: jmax

    !---------------------------------------------------------------------------

    FIO_IMAX = imax
    FIO_JMAX = jmax

    !--- reset comoninfo
    call fio_put_commoninfo( FIO_MPIIO_NOUSE, &
                             FIO_SPLIT_FILE,  &
                             FIO_BIG_ENDIAN,  &
                             FIO_CARTESIAN,   &
                             int(resolution), &
                             imax*jmax,       &
#ifdef CONFIG_HDF5
                             imax,            &
                             jmax,            &
#endif
                             1,               &
                             PRC_myrank       )

    return
  end subroutine FIO_setgridinfo

  !-----------------------------------------------------------------------------
  subroutine FIO_getfid( &
      fid,      &
      basename, &
      rwtype,   &
      pkg_desc, &
      pkg_note  )
    use mod_process, only: &
       PRC_myrank
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
       call fio_mk_fname(fname,trim(basename),'pe',PRC_myrank,6)
       call fio_register_file(n,fname)

       if ( rwtype == FIO_FREAD ) then

          call fio_fopen(n,FIO_FREAD)
          call fio_read_allinfo(n)

       elseif( rwtype == FIO_FWRITE ) then

          call fio_fopen(n,FIO_FWRITE)
          call fio_put_write_pkginfo(n,pkg_desc,pkg_note)

       endif

       if( IO_L ) write(IO_FID_LOG,'(1x,3(A),i3)') '*** [FIO] File registration : ',trim(rwname(rwtype)),' -',n
       if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(fname)

       FIO_fname_list(FIO_fid_count) = trim(basename)
       FIO_fid_list  (FIO_fid_count) = n
       FIO_fid_count                 = FIO_fid_count + 1

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
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP),          intent(out) :: var(:,:,:)
    character(LEN=*), intent( in) :: basename
    character(LEN=*), intent( in) :: varname
    character(LEN=*), intent( in) :: layername
    integer,          intent( in) :: k_start, k_end
    integer,          intent( in) :: step

    logical, intent(in), optional :: allow_missingq !--- if data is missing, set value to zero

    real(4) :: var4(k_start:k_end,FIO_IMAX,FIO_JMAX)
    real(8) :: var8(k_start:k_end,FIO_IMAX,FIO_JMAX)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE I')

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
       if ( present(allow_missingq) ) then
          if ( allow_missingq ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                                            'varname= ',trim(varname),', step=',step
             if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO] Q Value is set to 0.'

             var(k_start:k_end,:,:) = 0.D0

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
       var(k_start:k_end,:,:) = real(var4(1:dinfo%num_of_layer,:,:),kind=RP)

    elseif( dinfo%datatype == FIO_REAL8 ) then

       call fio_read_data(fid,did,var8(:,:,:))
       var(k_start:k_end,:,:) = real(var8(1:dinfo%num_of_layer,:,:),kind=RP)

    endif

    call TIME_rapend  ('FILE I')

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
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only :&
      TIME_date2sec, &
      TIME_sec2date
    implicit none

    integer,          intent(inout) :: start_step
    integer,          intent(inout) :: num_of_step
    integer,          intent(inout) :: data_date(6,FIO_nstep_max)
    integer,          intent(inout) :: prec
    character(LEN=*), intent(   in) :: basename
    character(LEN=*), intent(   in) :: varname
    character(LEN=*), intent(   in) :: layername
    integer,          intent(   in) :: k_start, k_end
    real(RP),          intent(   in) :: ctime
    integer,          intent(   in) :: cdate(6)
    logical,          intent(   in) :: opt_periodic_year

    real(RP) :: midtime  !--- [sec]
    logical :: startflag
    integer :: did, fid
    integer :: i
    !---------------------------------------------------------------------------

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FREAD, "", "" )

    startflag = .false.

    do i = 1, FIO_nstep_max
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
       call TIME_sec2date( data_date(:,i), midtime )

       if ( opt_periodic_year ) then
          data_date(1,i) = cdate(1)
          call TIME_date2sec( midtime, data_date(:,i) )
       endif

       if (       ( .not. startflag ) &
            .AND. ( ctime < midtime ) ) then
          startflag  = .true.
          start_step = i
          prec       = FIO_preclist(dinfo%datatype)
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
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF4, &
       CONST_UNDEF8
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

    real(4) :: var4(k_start:k_end,FIO_IMAX,FIO_JMAX)
    real(8) :: var8(k_start:k_end,FIO_IMAX,FIO_JMAX)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O')

    !--- search/register file
    call FIO_getfid( fid, basename, FIO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    dinfo%varname      = varname
    dinfo%description  = data_desc
    dinfo%unit         = unit
    dinfo%layername    = layername
    dinfo%note         = data_note
    dinfo%datasize     = FIO_IMAX * FIO_JMAX * (k_end-k_start+1) * FIO_preclist(dtype)
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

       var4(k_start:k_end,:,:) = real(var(k_start:k_end,:,:),kind=4)
       where( var4(:,:,:) < (CONST_UNDEF4+1.0) )
          var4(:,:,:) = CONST_UNDEF4
       endwhere

       call fio_put_write_datainfo_data(did,fid,dinfo,var4(:,:,:))

    elseif( dtype == FIO_REAL8 ) then

       var8(k_start:k_end,:,:) = real(var(k_start:k_end,:,:),kind=8)
       where( var8(:,:,:) < (CONST_UNDEF8+1.0) )
          var8(:,:,:) = CONST_UNDEF8
       endwhere
       call fio_put_write_datainfo_data(did,fid,dinfo,var8(:,:,:))
    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call PRC_MPIstop
    endif

    call TIME_rapend  ('FILE O')

    return
  end subroutine FIO_output

  !-----------------------------------------------------------------------------
  subroutine FIO_finalize
    use mod_process, only: &
       PRC_myrank
    implicit none

    character(LEN=FIO_HLONG) :: fname
    integer                  :: n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    do n = 1, FIO_fid_count-1
       call fio_fclose(FIO_fid_list(n))

       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3)') '*** [FIO] File Close : NO.', n
       call fio_mk_fname(fname,trim(FIO_fname_list(n)),'pe',PRC_myrank,6)
       if( IO_L ) write(IO_FID_LOG,*) '*** closed filename: ', trim(fname)
    enddo

    return
  end subroutine FIO_finalize

  !-----------------------------------------------------------------------------
  subroutine FIO_getfid_1D( &
      fid,      &
      basename, &
      rwtype,   &
      pkg_desc, &
      pkg_note, &
      single    )
    use mod_process, only: &
       PRC_master, &
       PRC_myrank
    implicit none

    integer,          intent(out) :: fid
    character(LEN=*), intent( in) :: basename
    integer,          intent( in) :: rwtype
    character(LEN=*), intent( in) :: pkg_desc
    character(LEN=*), intent( in) :: pkg_note
    logical,          intent( in) :: single

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
       if ( single ) then
          fname = trim(basename)//'.peall'
       else
          call fio_mk_fname(fname,trim(basename),'pe',PRC_myrank,6)
       endif
       call fio_register_file(n,fname)

       if ( rwtype == FIO_FREAD ) then

          call fio_fopen(n,FIO_FREAD)
          call fio_read_allinfo_novalid(n)

       elseif( rwtype == FIO_FWRITE ) then

          call fio_fopen(n,FIO_FWRITE)

          !--- append package header to the file
          if ( single ) then
             hinfo%fmode        = FIO_INTEG_FILE
             hinfo%num_of_rgn   = 1
             hinfo%rgnid(1)     = PRC_master
          else
             hinfo%fmode        = FIO_SPLIT_FILE
             hinfo%num_of_rgn   = 1
             hinfo%rgnid(1)     = PRC_myrank
          endif

          hinfo%description   = pkg_desc
          hinfo%note          = pkg_note
          hinfo%endiantype    = FIO_BIG_ENDIAN
          hinfo%grid_topology = FIO_NONE
          hinfo%glevel        = 0
          hinfo%rlevel        = 1

          call fio_put_pkginfo  (n,hinfo)
          call fio_write_pkginfo(n)

       endif

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '*** [FIO1D] File registration : ',trim(rwname(rwtype)),' -',n
       if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(fname)

       FIO_fname_list(FIO_fid_count) = trim(basename)
       FIO_fid_list  (FIO_fid_count) = n
       FIO_fid_count                   = FIO_fid_count + 1

       fid = n
    endif

    return
  end subroutine FIO_getfid_1D

  !-----------------------------------------------------------------------------
  subroutine FIO_input_1D( &
      var,       &
      basename,  &
      varname,   &
      layername, &
      k_start,   &
      k_end,     &
      step,      &
      single     )
    use mod_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    implicit none

    real(RP),          intent(out) :: var(:)
    character(LEN=*), intent( in) :: basename
    character(LEN=*), intent( in) :: varname
    character(LEN=*), intent( in) :: layername
    integer,          intent( in) :: k_start, k_end
    integer,          intent( in) :: step
    logical,          intent( in) :: single

    real(4) :: var4(k_start:k_end)
    real(8) :: var8(k_start:k_end)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE I')

    !--- search/register file
    call FIO_getfid_1D( fid, basename, FIO_FREAD, '', '', single )

    !--- seek data ID and get information
    call fio_seek_datainfo(did,fid,varname,step)
    call fio_get_datainfo (fid,did,dinfo)

    !--- verify
    if ( did == -1 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO1D] data not found! : ', &
                                      'varname= ',trim(varname),', step=',step
       call PRC_MPIstop
    endif

    if ( trim(dinfo%layername) /= trim(layername) ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO1D] layername mismatch! ', &
                                      "[",trim(dinfo%layername),":",trim(layername),"]"
       call PRC_MPIstop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[FIO1D] num_of_layer mismatch! ', &
                                      dinfo%num_of_layer,k_end-k_start+1
       call PRC_MPIstop
    endif

    !--- read data
    if ( dinfo%datatype == FIO_REAL4 ) then

       call fio_read_data(fid,did,var4(:))
       var(k_start:k_end) = real(var4(1:dinfo%num_of_layer),kind=RP)

    elseif( dinfo%datatype == FIO_REAL8 ) then

       call fio_read_data(fid,did,var8(:))
       var(k_start:k_end) = real(var8(1:dinfo%num_of_layer),kind=RP)

    endif

    call TIME_rapend  ('FILE I')

    return
  end subroutine FIO_input_1D

  !-----------------------------------------------------------------------------
  subroutine FIO_output_1D( &
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
      t_end,     &
      single     )
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       CONST_UNDEF4, &
       CONST_UNDEF8
    implicit none

    real(RP),          intent(in) :: var(:)
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
    logical,          intent(in) :: single

    real(4) :: var4(k_start:k_end)
    real(8) :: var8(k_start:k_end)

    integer :: did, fid
#ifdef CONFIG_HDF5
    integer :: gid
#endif

    !---------------------------------------------------------------------------

    call TIME_rapstart('FILE O')

    !--- search/register file
    call FIO_getfid_1D( fid, basename, FIO_FWRITE, pkg_desc, pkg_note, single )

    !--- append data to the file
    dinfo%varname      = varname
    dinfo%description  = data_desc
    dinfo%unit         = unit
    dinfo%layername    = layername
    dinfo%note         = data_note
    dinfo%datasize     = (k_end-k_start+1) * FIO_preclist(dtype)
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start * 1.D3, kind=8 )
    dinfo%time_end     = int( t_end   * 1.D3, kind=8 )

    if ( dtype == FIO_REAL4 ) then

       var4(k_start:k_end)=real(var(k_start:k_end),kind=4)
       where( var4(:) < (CONST_UNDEF4+1.0) )
          var4(:) = CONST_UNDEF4
       endwhere

       call fio_put_write_datainfo_data(did,fid,dinfo,var4(:))

    elseif( dtype == FIO_REAL8 ) then

       var8(k_start:k_end)=real(var(k_start:k_end),kind=8)
       where( var8(:) < (CONST_UNDEF8+1.0) )
          var8(:) = CONST_UNDEF8
       endwhere

       call fio_put_write_datainfo_data(did,fid,dinfo,var8(:))
    else

       if( IO_L ) write(IO_FID_LOG,*) 'xxx [OUTPUT]/[FIO1D] Unsupported datatype!', dtype
       call PRC_MPIstop

    endif

    call TIME_rapend  ('FILE O')

    return
  end subroutine FIO_output_1D

end module mod_fileio
!-------------------------------------------------------------------------------
