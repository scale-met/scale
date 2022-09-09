!-------------------------------------------------------------------------------
!> Module file I/O
!!
!! @par Description
!!         This module is container for file I/O (PaNDa format)
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
  use mod_io_param
  use scale_io
  use scale_prof
  use scale_atmos_grid_icoA_index
  use iso_c_binding
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
  public :: FIO_close
  public :: FIO_finalize

  interface FIO_input
     module procedure FIO_input_SP
     module procedure FIO_input_DP
  end interface FIO_input

  interface FIO_output
     module procedure FIO_output_SP
     module procedure FIO_output_DP
  end interface FIO_output

  public :: cstr
  public :: cstr2
  public :: fstr
  interface fstr
     module procedure fstr1
     module procedure fstr2
  end interface fstr

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  include 'fio_c.inc'
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,               private, parameter :: FIO_nmaxfile = 64
  character(len=H_LONG), private            :: FIO_fname_list(FIO_nmaxfile)
  integer,               private            :: FIO_fid_list  (FIO_nmaxfile)
  integer,               private            :: FIO_fid_count = 1

  type(datainfo),   private :: dinfo

  integer, private, parameter :: max_num_of_data = 2500 !--- max time step num
  integer, private, parameter :: preclist(0:3) = (/ 4, 8, 4, 8 /)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine FIO_setup
    use scale_prc_icoA, only: &
       PRC_RGN_level, &
       PRC_RGN_local, &
       PRC_RGN_l2r
    implicit none

    integer, allocatable :: prc_tab(:)
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[fio]/Category[common share]'

    allocate( prc_tab(PRC_RGN_local) )
    prc_tab(1:PRC_RGN_local) = PRC_RGN_l2r(1:PRC_RGN_local)-1

    ierr = fio_syscheck()
    ierr = fio_put_commoninfo( IO_SPLIT_FILE,  & ! [IN]
                               IO_BIG_ENDIAN,  & ! [IN]
                               IO_ICOSAHEDRON, & ! [IN]
                               ADM_glevel,     & ! [IN]
                               PRC_RGN_level,  & ! [IN]
                               PRC_RGN_local,  & ! [IN]
                               prc_tab         ) ! [IN]

    deallocate(prc_tab)

    return
  end subroutine FIO_setup

  !-----------------------------------------------------------------------------
  !> Get file ID of given basename.
  subroutine FIO_getfid( &
       fid,      &
       basename, &
       rwtype,   &
       pkg_desc, &
       pkg_note  )
    use scale_prc, only: &
       PRC_myrank
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of file
    integer,          intent(in)  :: rwtype   !< file access type
    character(len=*), intent(in)  :: pkg_desc !< package(file) description
    character(len=*), intent(in)  :: pkg_note !< package(file) note

    character(len=H_SHORT) :: rwname(0:2)
    data rwname / 'READ','WRITE','APPEND' /

    character(len=H_LONG) :: fname
    integer               :: ierr
    integer               :: n
    !---------------------------------------------------------------------------

    !--- search existing file
    fid = -1
    do n = 1, FIO_fid_count
       if ( basename == FIO_fname_list(n) ) fid = FIO_fid_list(n)
    enddo

    if ( fid < 0 ) then ! file registration
       !--- register new file and open
       call fio_mk_fname(fname,cstr(basename),cstr('pe'),PRC_myrank,6)
       call fstr(fname)
       fid = fio_register_file(fname)

       if ( rwtype == IO_FREAD ) then

!          call fio_dump_finfo(n,FIO_BIG_ENDIAN,FIO_DUMP_HEADER) ! dump to stdout(check)
          ierr = fio_fopen(fid,rwtype)
          ierr = fio_read_allinfo(fid)

       elseif( rwtype == IO_FWRITE ) then

          ierr = fio_fopen(fid,rwtype)
          ierr = fio_put_write_pkginfo(fid,cstr(pkg_desc),cstr(pkg_note))

       elseif( rwtype == IO_FAPPEND ) then

          ierr = fio_fopen(fid,rwtype)
          ierr = fio_read_pkginfo(fid)
          ierr = fio_write_pkginfo(fid)

       endif

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I3)') '*** [FIO] File registration (ADVANCED) : ', &
                            trim(rwname(rwtype)),' - ', FIO_fid_count
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A)') '*** fid= ', fid, ', name: ', trim(fname(1:IO_HLONG))

       FIO_fname_list(FIO_fid_count) = trim(basename)
       FIO_fid_list  (FIO_fid_count) = fid
       FIO_fid_count = FIO_fid_count + 1
    endif

    return
  end subroutine FIO_getfid

  !-----------------------------------------------------------------------------
  !> Input(read) one variable at one step.
  subroutine FIO_input_SP( &
       var,           &
       basename,      &
       varname,       &
       layername,     &
       k_start,       &
       k_end,         &
       step,          &
       allow_missingq ) !--- optional
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(SP),         intent(out) :: var(:,:,:) !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename   !< basename of file
    character(len=*), intent(in)  :: varname    !< variable name
    character(len=*), intent(in)  :: layername  !< layer name
    integer,          intent(in)  :: k_start    !< start index of vertical level
    integer,          intent(in)  :: k_end      !< end   index of vertical level
    integer,          intent(in)  :: step       !< step to be read

    logical, intent(in), optional :: allow_missingq !< if data is missing, set value to zero, else execution stops.

    real(SP), target :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP), target :: var8(ADM_gall,k_start:k_end,ADM_lall)

    character(len=IO_HSHORT) :: lname

    integer :: did, fid
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_CARTESC_in',2)

    !--- search/register file
    call FIO_getfid( fid, basename, IO_FREAD, '', '' )

    !--- seek data ID and get information
    did = fio_seek_datainfo(fid,cstr(varname),step)
    ierr = fio_get_datainfo(dinfo,fid,did)
    call fstr( lname, dinfo%layername )

    !--- verify
    if ( did == -1 ) then
       if ( present(allow_missingq) ) then
          if ( allow_missingq ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[FIO] data not found! : ', &
                                  'varname= ',trim(varname),', step=',step
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[FIO] Q Value is set to 0.'

             var(:,k_start:k_end,:) = 0.0_SP

             call PROF_rapend('FILE_CARTESC_in',2)
             return
          endif
       else
          write(*,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                               'varname= ',trim(varname),', step=',step
          call PRC_abort
       endif
    endif

    if ( lname /= layername ) then
       write(*,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                            '[',trim(lname),':',trim(layername),']'
       call PRC_abort
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       write(*,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch! ', &
                            dinfo%num_of_layer,k_end-k_start+1
       call PRC_abort
    endif

    !--- read data
    if ( dinfo%datatype == IO_REAL4 ) then

       ierr = fio_read_data(fid,did,c_loc(var4))
       var(:,k_start:k_end,:) = real( var4(:,1:dinfo%num_of_layer,:), kind=SP )

    elseif( dinfo%datatype == IO_REAL8 ) then

       ierr = fio_read_data(fid,did,c_loc(var8))
       var(:,k_start:k_end,:) = real( var8(:,1:dinfo%num_of_layer,:), kind=SP )

    endif

    call PROF_rapend('FILE_CARTESC_in',2)

    return
  end subroutine FIO_input_SP

  !-----------------------------------------------------------------------------
  !> Input(read) one variable at one step.
  subroutine FIO_input_DP( &
       var,           &
       basename,      &
       varname,       &
       layername,     &
       k_start,       &
       k_end,         &
       step,          &
       allow_missingq ) !--- optional
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(DP),         intent(out) :: var(:,:,:) !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename   !< basename of file
    character(len=*), intent(in)  :: varname    !< variable name
    character(len=*), intent(in)  :: layername  !< layer name
    integer,          intent(in)  :: k_start    !< start index of vertical level
    integer,          intent(in)  :: k_end      !< end   index of vertical level
    integer,          intent(in)  :: step       !< step to be read

    logical, intent(in), optional :: allow_missingq !< if data is missing, set value to zero, else execution stops.

    real(SP), target :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP), target :: var8(ADM_gall,k_start:k_end,ADM_lall)

    character(len=IO_HSHORT) :: lname

    integer :: did, fid
    integer :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_CARTESC_in',2)

    !--- search/register file
    call FIO_getfid( fid, basename, IO_FREAD, '', '' )

    !--- seek data ID and get information
    did = fio_seek_datainfo(fid,cstr(varname),step)
    ierr = fio_get_datainfo(dinfo,fid,did)
    call fstr(lname, dinfo%layername)

    !--- verify
    if ( did == -1 ) then
       if ( present(allow_missingq) ) then
          if ( allow_missingq ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[FIO] data not found! : ', &
                                  'varname= ',trim(varname),', step=',step
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[FIO] Q Value is set to 0.'

             var(:,k_start:k_end,:) = 0.0_DP

             call PROF_rapend('FILE_CARTESC_in',2)
             return
          endif
       else
          write(*,*) 'xxx [INPUT]/[FIO] data not found! : ', &
                               'varname= ',trim(varname),', step=',step
          call PRC_abort
       endif
    endif

    if ( lname /= layername ) then
       write(*,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                            '[',trim(lname),':',trim(layername),']'
       call PRC_abort
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       write(*,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch! ', &
                            dinfo%num_of_layer,k_end-k_start+1
       call PRC_abort
    endif

    !--- read data
    if ( dinfo%datatype == IO_REAL4 ) then

       ierr = fio_read_data(fid,did,c_loc(var4))
       var(:,k_start:k_end,:) = real( var4(:,1:dinfo%num_of_layer,:), kind=DP )

    elseif( dinfo%datatype == IO_REAL8 ) then

       ierr = fio_read_data(fid,did,c_loc(var8))
       var(:,k_start:k_end,:) = real( var8(:,1:dinfo%num_of_layer,:), kind=DP )

    endif

    call PROF_rapend('FILE_CARTESC_in',2)

    return
  end subroutine FIO_input_DP

  !-----------------------------------------------------------------------------
  !> Read in all steps of given `varname`, returns total data
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
    use scale_prc, only: &
       PRC_abort
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
    character(len=*), intent(in)    :: layername      ! for verification only
    integer,          intent(in)    :: k_start, k_end ! for verification only
    real(DP),         intent(in)    :: ctime
    integer,          intent(in)    :: cdate(6)       ! cdate(1) is only used only when opt_periodic_year is T.
    logical,          intent(in)    :: opt_periodic_year

    real(DP) :: midtime ! [sec]
    integer  :: midday, offset_year
    real(DP) :: midsec, midms

    character(len=IO_HSHORT) :: lname
    logical  :: startflag
    integer  :: did, fid
    integer  :: ierr
    integer  :: i
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_CARTESC_in',2)

    !--- search/register file
    call FIO_getfid( fid, basename, IO_FREAD, '', '' )

    startflag = .false.

    do i = 1, max_num_of_data
       !--- seek data ID and get information
       did = fio_seek_datainfo(fid,cstr(varname),i)
       ierr = fio_get_datainfo(dinfo,fid,did)
       call fstr(lname, dinfo%layername)

       if ( did == -1 ) then
          num_of_step = i - 1
          exit
       endif

       !--- verify
       if ( lname /= layername ) then
          write(*,*) 'xxx [INPUT]/[FIO] layername mismatch! ', &
                               '[',trim(lname),':',trim(layername),']'
          call PRC_abort
       elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
          write(*,*) 'xxx [INPUT]/[FIO] num_of_layer mismatch!', &
                               dinfo%num_of_layer,k_end-k_start+1
          call PRC_abort
       endif

       midtime = real( int( (dinfo%time_start+dinfo%time_end)*0.5_DP+1.0_DP, kind=DP ), kind=DP )

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

    call PROF_rapend('FILE_CARTESC_in',2)

    return
  end subroutine FIO_seek

  !-----------------------------------------------------------------------------
  !> Append data with data header
  subroutine FIO_output_SP( &
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
       append     )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4
    implicit none

    real(SP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: basename
    character(len=*), intent(in) :: pkg_desc
    character(len=*), intent(in) :: pkg_note
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: data_desc
    character(len=*), intent(in) :: data_note
    character(len=*), intent(in) :: unit
    integer,          intent(in) :: dtype
    character(len=*), intent(in) :: layername
    integer,          intent(in) :: k_start, k_end
    integer,          intent(in) :: step
    real(DP),         intent(in) :: t_start, t_end

    logical,intent(in), optional :: append

    real(SP), target :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP), target :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_CARTESC_out',2)

    !--- search/register file
    call FIO_getfid( fid, basename, IO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    call cstr2(dinfo%varname,     varname)
    call cstr2(dinfo%description, data_desc)
    call cstr2(dinfo%unit,        unit)
    call cstr2(dinfo%layername,   layername)
    call cstr2(dinfo%note,        data_note)
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * preclist(dtype), kind=DP )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start, kind=DP )
    dinfo%time_end     = int( t_end,   kind=DP )

    if ( dtype == IO_REAL4 ) then

       var4(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=SP )
       where( var4(:,:,:) < (UNDEF4+1.0_SP) )
          var4(:,:,:) = UNDEF4
       endwhere

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var4))

    elseif( dtype == IO_REAL8 ) then

       var8(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=DP )

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var8))

    else
       write(*,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call PRC_abort
    endif

    call PROF_rapend('FILE_CARTESC_out',2)

    return
  end subroutine FIO_output_SP

  !-----------------------------------------------------------------------------
  !> Append data with data header
  subroutine FIO_output_DP( &
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
       append     )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF4 => CONST_UNDEF4
    implicit none

    real(DP),         intent(in) :: var(:,:,:)
    character(len=*), intent(in) :: basename
    character(len=*), intent(in) :: pkg_desc
    character(len=*), intent(in) :: pkg_note
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: data_desc
    character(len=*), intent(in) :: data_note
    character(len=*), intent(in) :: unit
    integer,          intent(in) :: dtype
    character(len=*), intent(in) :: layername
    integer,          intent(in) :: k_start, k_end
    integer,          intent(in) :: step
    real(DP),         intent(in) :: t_start, t_end

    logical,intent(in), optional :: append

    real(SP), target :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP), target :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILE_CARTESC_out',2)

    !--- search/register file
    call FIO_getfid( fid, basename, IO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    call cstr2(dinfo%varname,     varname)
    call cstr2(dinfo%description, data_desc)
    call cstr2(dinfo%unit,        unit)
    call cstr2(dinfo%layername,   layername)
    call cstr2(dinfo%note,        data_note)
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * preclist(dtype), kind=DP )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1
    dinfo%step         = step
    dinfo%time_start   = int( t_start, kind=DP )
    dinfo%time_end     = int( t_end,   kind=DP )

    if ( dtype == IO_REAL4 ) then

       var4(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=SP )
       where( var4(:,:,:) < (UNDEF4+1.0_SP) )
          var4(:,:,:) = UNDEF4
       endwhere

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var4))

    elseif( dtype == IO_REAL8 ) then

       var8(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=DP )

       did = fio_put_write_datainfo_data(fid,dinfo,c_loc(var8))

    else
       write(*,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call PRC_abort
    endif

    call PROF_rapend('FILE_CARTESC_out',2)

    return
  end subroutine FIO_output_DP

  !-----------------------------------------------------------------------------
  subroutine FIO_close( &
       basename )
    use scale_prc, only: &
       PRC_myrank
    implicit none

    character(len=*), intent(in) :: basename

    character(len=H_LONG) :: fname

    integer :: fid
    integer :: ierr
    integer :: n
    !---------------------------------------------------------------------------

    !--- search/register file
    do n = 1, FIO_fid_count
       if ( basename == FIO_fname_list(n) ) then
          fid = FIO_fid_list(n)

          ierr = fio_fclose(fid)
          call fio_mk_fname(fname,cstr(FIO_fname_list(n)),cstr('pe'),PRC_myrank,6)
          call fstr(fname)

          if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A)') &
          '*** [FIO] File close (ADVANCED) fid= ', fid, ', name: ', trim(fname)

          ! remove closed file info from the list
          FIO_fname_list(n) = ''
          FIO_fid_list  (n) = -1
       endif
    enddo

    return
  end subroutine FIO_close

  !-----------------------------------------------------------------------------
  subroutine FIO_finalize
    use scale_prc, only: &
       PRC_myrank
    implicit none

    character(len=H_LONG) :: fname
    integer               :: n, fid
    integer               :: ierr
    !---------------------------------------------------------------------------

    do n = 1, FIO_fid_count
       fid = FIO_fid_list(n)

       ierr = fio_fclose(fid)
       call fio_mk_fname(fname,cstr(FIO_fname_list(n)),cstr('pe'),PRC_myrank,6)
       call fstr(fname)

       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A)') &
       '*** [FIO] File close (ADVANCED) fid= ', fid, ', name: ', trim(fname)
    enddo

    return
  end subroutine FIO_finalize

  function cstr(str)
    character(*), intent(in) :: str
    character(:,c_char), allocatable, target :: cstr
    cstr = trim(str) // c_null_char
  end function cstr

  subroutine cstr2(cstr, fstr)
    character(c_char), intent(out) :: cstr(:)
    character(len=*),  intent(in)  :: fstr
    integer :: i, j
    integer :: l
    l = min( len(fstr), size(cstr)-1 )
    cstr(l+1) = c_null_char
    do i = l, 1, -1
       if ( fstr(i:i) == " " ) then
          cstr(i) = c_null_char
       else
          exit
       end if
    end do
    do j = i, 1, -1
       cstr(j) = fstr(j:j)
    end do
    return
  end subroutine cstr2

  subroutine fstr1(str)
    character(len=*), intent(inout) :: str
    integer :: i, j
    do i = 1, len(str)
       if ( str(i:i) == c_null_char ) exit
    end do
    do j = i, len(str)
       str(j:j) = " "
    end do
    return
  end subroutine fstr1

  subroutine fstr2(fstr, cstr)
    character(len=*),  intent(out) :: fstr
    character(c_char), intent(in)  :: cstr(:)
    integer :: i, j
    integer :: l
    l = min( len(fstr), size(cstr) )
    do i = 1, l
       if ( cstr(i) == c_null_char ) exit
       fstr(i:i) = cstr(i)
    end do
    do j = i, len(fstr)
       fstr(j:j) = " "
    end do
    return
  end subroutine fstr2

end module mod_fio
