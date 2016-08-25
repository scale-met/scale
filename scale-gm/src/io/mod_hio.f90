!-------------------------------------------------------------------------------
!> Module file I/O
!!
!! @par Description
!!         This module is container for file I/O (HDF5 format)
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_hio
  !
  !++ Used modules
  !
  use scale_precision
  use mod_io_param
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: HIO_setup
  public :: HIO_input
  public :: HIO_seek
  public :: HIO_output
  public :: HIO_close
  public :: HIO_finalize

  interface HIO_input
     module procedure HIO_input_SP
     module procedure HIO_input_DP
  end interface HIO_input

  interface HIO_output
     module procedure HIO_output_SP
     module procedure HIO_output_DP
  end interface HIO_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !> struct for package infomation
  type, public :: headerinfo
     character(len=H_LONG)  :: fname         !< file name
     character(len=H_MID)   :: description   !< file description
     character(len=H_LONG)  :: note          !< longer note of file
     integer                :: num_of_var    !< number of data
     integer                :: fmode         !< file mode(0,1,2)
     integer                :: endiantype    !< endian type(0,1,2)
     integer                :: grid_topology !< grid topology(0,1,2)
     integer                :: glevel        !< glevel
     integer                :: rlevel        !< rlevel
     integer                :: num_of_rgn    !< number of region
     integer, pointer       :: rgnid(:)      !< array of region id
  endtype headerinfo

  !> struct for data infomation
  type, public :: datainfo
     character(len=H_SHORT) :: varname      !< variable name
     character(len=H_MID)   :: description  !< variable description
     character(len=H_SHORT) :: unit         !< unit of variable
     character(len=H_SHORT) :: layername    !< layer name
     character(len=H_LONG)  :: note         !< longer note of variable
     integer(DP)            :: datasize     !< data size
     integer                :: datatype     !< data type(0,1,2,3)
     integer                :: num_of_layer !< number of layer
     integer                :: num_of_step  !< number of step
  endtype datainfo

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,               private, parameter :: HIO_nmaxfile = 64
  character(len=H_LONG), private            :: HIO_fname_list(HIO_nmaxfile)
  integer,               private            :: HIO_fid_list  (HIO_nmaxfile)
  integer,               private            :: HIO_fid_count = 1

  type(headerinfo), private :: hinfo
  type(datainfo),   private :: dinfo

  integer, private, parameter :: max_num_of_data = 2500 !--- max time step num
  integer, private, parameter :: preclist(0:3) = (/ 4, 8, 4, 8 /)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup this module.
  !!
  !! Must be called first.
  !!
  subroutine HIO_setup
    use mod_adm, only: &
       ADM_prc_me,  &
       ADM_prc_tab, &
       ADM_glevel,  &
       ADM_rlevel,  &
       ADM_lall
    implicit none

    integer, allocatable :: prc_tab(:)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[hio]/Category[common share]'

    allocate( prc_tab(ADM_lall) )

    prc_tab(1:ADM_lall) = ADM_prc_tab(1:ADM_lall,ADM_prc_me)-1

    call hio_syscheck()
    call hio_put_commoninfo( IO_SPLIT_FILE,  & ! [IN]
                             IO_BIG_ENDIAN,  & ! [IN]
                             IO_ICOSAHEDRON, & ! [IN]
                             ADM_glevel,     & ! [IN]
                             ADM_rlevel,     & ! [IN]
                             ADM_lall,       & ! [IN]
                             prc_tab         ) ! [IN]

    deallocate(prc_tab)

    allocate( hinfo%rgnid(ADM_lall) )

    return
  end subroutine HIO_setup

  !-----------------------------------------------------------------------------
  !> Get file ID of given basename.
  !!
  !! Open it if not opened yet.
  !!
  subroutine HIO_getfid( &
       fid,      &
       basename, &
       rwtype,   &
       pkg_desc, &
       pkg_note  )
    use mod_adm, only: &
       ADM_prc_me
    implicit none

    integer,          intent(out) :: fid      !< file ID
    character(len=*), intent(in)  :: basename !< basename of file
    integer,          intent(in)  :: rwtype   !< file access type
    character(len=*), intent(in)  :: pkg_desc !< package(file) description
    character(len=*), intent(in)  :: pkg_note !< package(file) note

    character(len=H_SHORT) :: rwname(0:2)
    data rwname / 'READ','WRITE','APPEND' /

    character(len=H_LONG) :: fname
    integer               :: n
    !---------------------------------------------------------------------------

    !--- search existing file
    fid = -1
    do n = 1, HIO_fid_count
       if ( basename == HIO_fname_list(n) ) fid = HIO_fid_list(n)
    enddo

    if ( fid < 0 ) then ! file registration
       !--- register new file and open
       call hio_mk_fname(fname,trim(basename),'pe',ADM_prc_me-1,6)
       call hio_register_file(fid,fname)

       if ( rwtype == IO_FREAD ) then

!          call hio_dump_finfo(n,HIO_BIG_ENDIAN,HIO_DUMP_HEADER) ! dump to stdout(check)
          call hio_fopen(fid,IO_FREAD)
          call hio_read_allinfo(fid)

       elseif( rwtype == IO_FWRITE ) then

          call hio_fopen(fid,IO_FWRITE)
          call hio_put_write_pkginfo(fid,pkg_desc,pkg_note)

       elseif( rwtype == IO_FAPPEND ) then

          call hio_fopen(fid,IO_FAPPEND)
          call hio_read_pkginfo(fid)
          call hio_write_pkginfo(fid)

       endif

       if( IO_L ) write(IO_FID_LOG,'(1x,A,A,A,I3)') '*** [HIO] File registration (ADVANCED) : ', &
                            trim(rwname(rwtype)),' - ', HIO_fid_count
       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A)') '*** fid= ', fid, ', name: ', trim(fname)

       HIO_fname_list(HIO_fid_count) = trim(basename)
       HIO_fid_list  (HIO_fid_count) = fid
       HIO_fid_count = HIO_fid_count + 1
    endif

    return
  end subroutine HIO_getfid

  !-----------------------------------------------------------------------------
  !> Input(read) one variable at one step.
  subroutine HIO_input_SP( &
       var,           &
       basename,      &
       varname,       &
       layername,     &
       k_start,       &
       k_end,         &
       step,          &
       allow_missingq ) !--- optional
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
    implicit none

    real(SP),         intent(out) :: var(:,:,:) !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename   !< basename of file
    character(len=*), intent(in)  :: varname    !< variable name
    character(len=*), intent(in)  :: layername  !< layer name
    integer,          intent(in)  :: k_start    !< start index of vertical level
    integer,          intent(in)  :: k_end      !< end   index of vertical level
    integer,          intent(in)  :: step       !< step to be read

    logical, intent(in), optional :: allow_missingq !< if data is missing, set value to zero, else execution stops.

    real(SP) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer(DP) :: ts !! time_start
    integer(DP) :: te !! time_end

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in',2)

    !--- search/register file
    call HIO_getfid( fid, basename, IO_FREAD, '', '' )

    !--- seek data ID and get information
    call hio_seek_datainfo(did,fid,varname,step)
    call hio_get_datainfo(fid,did,dinfo)

    !--- verify
    if ( did == -1 ) then
       if ( present(allow_missingq) ) then
          if ( allow_missingq ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[HIO] data not found! : ', &
                                  'varname= ',trim(varname),', step=',step
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[HIO] Q Value is set to 0.'

             var(:,k_start:k_end,:) = 0.0_SP

             call PROF_rapend('FILEIO_in',2)
             return
          endif
       else
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] data not found! : ', &
                               'varname= ',trim(varname),', step=',step
          call PRC_MPIstop
       endif
    endif

    if ( dinfo%layername /= layername ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] layername mismatch! ', &
                            '[',trim(dinfo%layername),':',trim(layername),']'
       call PRC_MPIstop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] num_of_layer mismatch! ', &
                            dinfo%num_of_layer,k_end-k_start+1
       call PRC_MPIstop
    endif

    !--- read data
    if ( dinfo%datatype == IO_REAL4 ) then

       call hio_read_data(fid,did,step,ts,te,var4(:,:,:))
       var(:,k_start:k_end,:) = real( var4(:,1:dinfo%num_of_layer,:), kind=SP )

    elseif( dinfo%datatype == IO_REAL8 ) then

       call hio_read_data(fid,did,step,ts,te,var8(:,:,:))
       var(:,k_start:k_end,:) = real( var8(:,1:dinfo%num_of_layer,:), kind=SP )

    endif

    call PROF_rapend('FILEIO_in',2)

    return
  end subroutine HIO_input_SP

  !-----------------------------------------------------------------------------
  !> Input(read) one variable at one step.
  subroutine HIO_input_DP( &
       var,           &
       basename,      &
       varname,       &
       layername,     &
       k_start,       &
       k_end,         &
       step,          &
       allow_missingq ) !--- optional
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
    implicit none

    real(DP),         intent(out) :: var(:,:,:) !< variable(ij,k,l)
    character(len=*), intent(in)  :: basename   !< basename of file
    character(len=*), intent(in)  :: varname    !< variable name
    character(len=*), intent(in)  :: layername  !< layer name
    integer,          intent(in)  :: k_start    !< start index of vertical level
    integer,          intent(in)  :: k_end      !< end   index of vertical level
    integer,          intent(in)  :: step       !< step to be read

    logical, intent(in), optional :: allow_missingq !< if data is missing, set value to zero, else execution stops.

    real(SP) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer(DP) :: ts !! time_start
    integer(DP) :: te !! time_end

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in',2)

    !--- search/register file
    call HIO_getfid( fid, basename, IO_FREAD, '', '' )

    !--- seek data ID and get information
    call hio_seek_datainfo(did,fid,varname,step)
    call hio_get_datainfo(fid,did,dinfo)

    !--- verify
    if ( did == -1 ) then
       if ( present(allow_missingq) ) then
          if ( allow_missingq ) then
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[HIO] data not found! : ', &
                                  'varname= ',trim(varname),', step=',step
             if( IO_L ) write(IO_FID_LOG,*) '*** [INPUT]/[HIO] Q Value is set to 0.'

             var(:,k_start:k_end,:) = 0.0_DP

             call PROF_rapend('FILEIO_in',2)
             return
          endif
       else
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] data not found! : ', &
                               'varname= ',trim(varname),', step=',step
          call PRC_MPIstop
       endif
    endif

    if ( dinfo%layername /= layername ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] layername mismatch! ', &
                            '[',trim(dinfo%layername),':',trim(layername),']'
       call PRC_MPIstop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] num_of_layer mismatch! ', &
                            dinfo%num_of_layer,k_end-k_start+1
       call PRC_MPIstop
    endif

    !--- read data
    if ( dinfo%datatype == IO_REAL4 ) then

       call hio_read_data(fid,did,step,ts,te,var4(:,:,:))
       var(:,k_start:k_end,:) = real( var4(:,1:dinfo%num_of_layer,:), kind=DP )

    elseif( dinfo%datatype == IO_REAL8 ) then

       call hio_read_data(fid,did,step,ts,te,var8(:,:,:))
       var(:,k_start:k_end,:) = real( var8(:,1:dinfo%num_of_layer,:), kind=DP )

    endif

    call PROF_rapend('FILEIO_in',2)

    return
  end subroutine HIO_input_DP

  !-----------------------------------------------------------------------------
  !> Read in all steps of given `varname`, returns total data
  !! size(`num_of_step`) and mid of (time_start+time_end) of each
  !! step(`data_date`).
  !!
  !! `start_step` is maximum step where `ctime < 0.5*(ts(step)+te(step))` is true.
  !!
  !! `prec` is presicion, 4 or 8.
  !!
  !! If `opt_periodic_year` is T, data_date(:,1) is set as cdate(1) on
  !! return, else cdate(:) is neglected.
  !!
  subroutine HIO_seek( &
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
    use scale_process, only: &
       PRC_MPIstop
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

    integer(DP), allocatable :: ts(:)
    integer(DP), allocatable :: te(:)
    integer                  :: num_of_var

    real(DP) :: midtime ! [sec]
    integer  :: midday, offset_year
    real(DP) :: midsec, midms

    logical  :: startflag
    integer  :: did, fid
    integer  :: i
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_in',2)

    !--- search/register file
    call HIO_getfid( fid, basename, IO_FREAD, '', '' )

    startflag = .false.

    call hio_get_num_of_var(fid,num_of_var)

    !--- seek data ID and get information
    call hio_seek_datainfo(did,fid,varname,i) ! i is meaningless, for compatibility
    call hio_get_datainfo (fid,did,dinfo)     ! here dinfo must contain actual ts(:) and te(:).
    num_of_step = dinfo%num_of_step

    allocate( ts(dinfo%num_of_step) )
    allocate( te(dinfo%num_of_step) )
    call hio_get_timeinfo(fid,did,ts,te)

    !--- verify
    if ( dinfo%layername /= layername ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] layername mismatch! ', &
                            '[',trim(dinfo%layername),':',trim(layername),']'
       call PRC_MPIstop
    elseif( dinfo%num_of_layer /= k_end-k_start+1 ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[HIO] num_of_layer mismatch!', &
                            dinfo%num_of_layer,k_end-k_start+1
       call PRC_MPIstop
    endif

    do i = 1, num_of_step
       midtime = real( int( (ts(i)+te(i))*0.5_DP+1.0_DP, kind=DP ), kind=DP )

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

    deallocate( ts )
    deallocate( te )

    call PROF_rapend('FILEIO_in',2)

    return
  end subroutine HIO_seek

  !-----------------------------------------------------------------------------
  !> Append data with data header
  subroutine HIO_output_SP( &
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
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
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

    real(SP) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer(DP) :: ts, te

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_out',2)

    !--- search/register file
    call HIO_getfid( fid, basename, IO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    dinfo%varname      = varname
    dinfo%description  = data_desc
    dinfo%unit         = unit
    dinfo%layername    = layername
    dinfo%note         = data_note
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * preclist(dtype), kind=DP )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1

    ts                 = int( t_start, kind=DP )
    te                 = int( t_end,   kind=DP )

    if ( dtype == IO_REAL4 ) then

       var4(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=SP )
       where( var4(:,:,:) < (UNDEF4+1.0_SP) )
          var4(:,:,:) = UNDEF4
       endwhere

       call hio_put_write_datainfo_data(did,fid,step,ts,te,dinfo,var4(:,:,:))

    elseif( dtype == IO_REAL8 ) then

       var8(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=DP )

       call hio_put_write_datainfo_data(did,fid,step,ts,te,dinfo,var8(:,:,:))

    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call PRC_MPIstop
    endif

    call PROF_rapend('FILEIO_out',2)

    return
  end subroutine HIO_output_SP

  !-----------------------------------------------------------------------------
  !> Append data with data header
  subroutine HIO_output_DP( &
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
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_gall, &
       ADM_lall
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

    real(SP) :: var4(ADM_gall,k_start:k_end,ADM_lall)
    real(DP) :: var8(ADM_gall,k_start:k_end,ADM_lall)

    integer(DP) :: ts, te

    integer :: did, fid
    !---------------------------------------------------------------------------

    call PROF_rapstart('FILEIO_out',2)

    !--- search/register file
    call HIO_getfid( fid, basename, IO_FWRITE, pkg_desc, pkg_note )

    !--- append data to the file
    dinfo%varname      = varname
    dinfo%description  = data_desc
    dinfo%unit         = unit
    dinfo%layername    = layername
    dinfo%note         = data_note
    dinfo%datasize     = int( ADM_gall * ADM_lall * (k_end-k_start+1) * preclist(dtype), kind=DP )
    dinfo%datatype     = dtype
    dinfo%num_of_layer = k_end-k_start+1

    ts                 = int( t_start, kind=DP )
    te                 = int( t_end,   kind=DP )

    if ( dtype == IO_REAL4 ) then

       var4(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=SP )
       where( var4(:,:,:) < (UNDEF4+1.0_SP) )
          var4(:,:,:) = UNDEF4
       endwhere

       call hio_put_write_datainfo_data(did,fid,step,ts,te,dinfo,var4(:,:,:))

    elseif( dtype == IO_REAL8 ) then

       var8(:,k_start:k_end,:) = real( var(:,k_start:k_end,:), kind=DP )

       call hio_put_write_datainfo_data(did,fid,step,ts,te,dinfo,var8(:,:,:))

    else
       if( IO_L ) write(IO_FID_LOG,*) 'xxx [OUTPUT]/[FIO] Unsupported datatype!', dtype
       call PRC_MPIstop
    endif

    call PROF_rapend('FILEIO_out',2)

    return
  end subroutine HIO_output_DP

  !-----------------------------------------------------------------------------
  subroutine HIO_close( &
       basename )
    use mod_adm, only: &
       ADM_prc_me
    implicit none

    character(len=*), intent(in) :: basename

    character(len=H_LONG) :: fname

    integer :: fid
    integer :: n
    !---------------------------------------------------------------------------

    !--- search/register file
    do n = 1, HIO_fid_count
       if ( basename == HIO_fname_list(n) ) then
          fid = HIO_fid_list(n)

          call hio_fclose(fid)
          call hio_mk_fname(fname,trim(HIO_fname_list(n)),'pe',ADM_prc_me-1,6)

          if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A)') &
          '*** [HIO] File close (ADVANCED) fid= ', fid, ', name: ', trim(fname)

          ! remove closed file info from the list
          HIO_fname_list(n) = ''
          HIO_fid_list  (n) = -1
       endif
    enddo

    return
  end subroutine HIO_close

  !-----------------------------------------------------------------------------
  subroutine HIO_finalize
    use mod_adm, only: &
       ADM_prc_me
    implicit none

    character(len=H_LONG) :: fname
    integer               :: n, fid
    !---------------------------------------------------------------------------

    do n = 1, HIO_fid_count
       fid = HIO_fid_list(n)

       call hio_fclose(fid)
       call hio_mk_fname(fname,trim(HIO_fname_list(n)),'pe',ADM_prc_me-1,6)

       if( IO_L ) write(IO_FID_LOG,'(1x,A,I3,A,A)') &
       '*** [HIO] File close (poh5) fid= ', fid, ', name: ', trim(fname)
    enddo

    return
  end subroutine HIO_finalize

end module mod_hio
