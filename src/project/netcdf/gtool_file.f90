!-------------------------------------------------------------------------------
!> module FILE
!!
!! @par Description
!!          unified hundring of various kinds of files
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-06-12 (S.Nishizawa) [new] Imported from SCALE-LES
!!
!<
!-------------------------------------------------------------------------------
module gtool_file
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_process, only: &
     PRC_MPIstop
  use mod_stdio, only: &
     IO_FILECHR,  &
     IO_FID_LOG,  &
     IO_L
  use gtool_file_h
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: FileOpen
  public :: FileCreate
  public :: FileAddVariable
  public :: FilePutAxis
  public :: FilePutAdditionalAxis
  public :: FileRead
  public :: FileWrite
  public :: FileClose
  public :: FileCloseAll

  interface FileRead
     module procedure FileRead1D
     module procedure FileRead2D
     module procedure FileRead3D
  end interface FileRead
  interface FileWrite
     module procedure FileWrite1D
     module procedure FileWrite2D
     module procedure FileWrite3D
  end interface FileWrite

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
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
  integer,                   private, parameter :: File_nfile_max = 64 ! number limit of file
  integer,                   private, parameter :: File_nvar_max  = 128 ! number limit of variables
  integer,                   private, parameter :: File_nstep_max = 2500 ! number limit of time step

  character(LEN=File_HLONG), private,      save :: File_bname_list(File_nfile_max)
  integer,                   private,      save :: File_fid_list  (File_nfile_max)
  integer,                   private,      save :: File_fid_count = 1
  character(LEN=File_HLONG), private,      save :: File_vname_list  (File_nvar_max)
  integer,                   private,      save :: File_vid_fid_list(File_nvar_max)
  integer,                   private,      save :: File_vid_list    (File_nvar_max)
  integer,                   private,      save :: File_vid_count = 1

contains
  !-----------------------------------------------------------------------------
  subroutine FileCreate( &
       fid,         & ! (out)
       basename,    & ! (in)
       title,       & ! (in)
       source,      & ! (in)
       institution, & ! (in)
       dim_name,    & ! (in)
       dim_size,    & ! (in)
       dim_desc,    & ! (in)
       dim_units,   & ! (in)
       dim_type,    & ! (in)
       single,      & ! (in) optional
       time_units   & ! (in) optional
       )
    use mod_process, only: &
       PRC_master, &
       PRC_myrank
    implicit none

    integer,          intent(out)           :: fid
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: title
    character(LEN=*), intent( in)           :: source
    character(LEN=*), intent( in)           :: institution
    character(LEN=*), intent( in)           :: dim_name(:)
    integer,          intent( in)           :: dim_size(:)
    character(LEN=*), intent( in)           :: dim_desc(:)
    character(LEN=*), intent( in)           :: dim_units(:)
    integer,          intent( in)           :: dim_type(:)
    character(LEN=*), intent( in), optional :: time_units
    logical,          intent( in), optional :: single

    character(len=File_HSHORT) :: time_units_
    integer :: nodeid
    logical :: single_ = .false.
    logical :: existed
    integer :: error

    if ( present(time_units) ) then
       time_units_ = time_units
    else
       time_units_ = 'sec'
    end if

    if ( present(single) ) single_ = single

    call FileGetfid(  &
         fid,         & ! (out)
         existed,     & ! (out)
         basename,    & ! (in)
         File_FWRITE, & ! (in)
         single_      & ! (in)
         )
         
    if ( existed ) return

    !--- append package header to the file
    if ( single_ ) then
       nodeid        = PRC_master
    else
       nodeid        = PRC_myrank
    endif

    call file_set_global_attributes( fid,                 & ! (in)
         title, source, institution, time_units_, nodeid, & ! (in)
         error                                            ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx failed to set global attributes'
       call PRC_MPIstop
    end if
    call file_set_dim_info( fid,                            & ! (in)
         size(dim_name),                                    & ! (in)
         dim_name, dim_size, dim_desc, dim_units, dim_type, & ! (in)
         error                                              ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx failed to set dimension information'
       call PRC_MPIstop
    end if

    return
  end subroutine FileCreate

  !-----------------------------------------------------------------------------
  subroutine FileOpen( &
      fid,       & ! (out)
      basename,  & ! (in)
      mode,      & ! (in)
      single     & ! (in) optional
      )
    implicit none

    integer,          intent(out) :: fid
    character(LEN=*), intent( in) :: basename
    integer,          intent( in) :: mode
    logical,          intent( in), optional :: single

    logical :: existed
    logical :: single_ = .false.

    if ( present(single) ) single_ = single

    call FileGetfid( fid,        & ! (out)
         existed,                & ! (out)
         basename, mode, single_ ) ! (in)
         
    return
  end subroutine FileOpen

  !-----------------------------------------------------------------------------
  subroutine FilePutAxis( &
       fid,      & ! (in)
       dim_name, & ! (in)
       val       & ! (in)
       )
    integer,          intent(in) :: fid
    character(len=*), intent(in) :: dim_name
    real(DP),         intent(in) :: val(:)

    integer error

    call file_put_axis( fid, & ! (in)
         dim_name, val, DP,  & ! (in)
         error               ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx failed to put axis value'
       call PRC_MPIstop
    end if
    
    return
  end subroutine FilePutAxis

  !-----------------------------------------------------------------------------
  subroutine FilePutAdditionalAxis( &
       fid,      & ! (in)
       name,     & ! (in)
       desc,     & ! (in)
       units,    & ! (in)
       dim_name, & ! (in)
       dtype,    & ! (in)
       val,      & ! (in)
       size      ) ! (in)
    integer,          intent(in) :: fid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: desc
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: dim_name
    integer,          intent(in) :: dtype
    real(DP),         intent(in) :: val(:)
    integer,          intent(in) :: size

    integer error

    call file_put_additional_axis( fid,                     & ! (in)
         name, desc, units, dim_name, dtype, val, size, DP, & ! (in)
         error                                              ) ! (out)
    if ( error /= SUCCESS_CODE .and. error /= ALREADY_EXISTED_CODE ) then
       write(*,*) 'xxx failed to put additional axis'
       call PRC_MPIstop
    end if

    return
  end subroutine FilePutAdditionalAxis
  !-----------------------------------------------------------------------------
  subroutine FileAddVariable( &
       vid,     & ! (out)
       fid,     & ! (in)
       varname, & ! (in)
       desc,    & ! (in)
       units,   & ! (in)
       dims,    & ! (in)
       dtype,   & ! (in)
       tint,    & ! (in) optional
       tavg     & ! (in) optional
       )
    integer,          intent(out) :: vid
    integer,          intent( in) :: fid
    character(len=*), intent( in) :: varname
    character(len=*), intent( in) :: desc
    character(len=*), intent( in) :: units
    character(len=*), intent( in) :: dims(:)
    integer,          intent( in) :: dtype
    real(DP),         intent( in), optional :: tint
    logical,          intent( in), optional :: tavg

    real(8) :: tint8
    integer :: itavg
    integer :: error
    integer :: n

    intrinsic size
    !---------------------------------------------------------------------------

    vid = -1
    do n = 1, File_vid_count
       if ( File_vid_fid_list(n) == fid .and. &
            trim(varname) == trim(File_vname_list(n)) ) then
          vid = File_vid_list(n)
       end if
    enddo

    if ( vid < 0 ) then ! variable registration
       !--- register new variable
       if( IO_L ) write(IO_FID_LOG,*) '*** [File] Var registration'
       if( IO_L ) write(IO_FID_LOG,*) '*** variable name: ', trim(varname)

       if ( present(tint) ) then
          tint8 = tint
       else
          tint8 = -1.0D0
       end if
       if ( present(tavg) .and. tavg ) then
          itavg = 1
       else
          itavg = 0
       end if
       call file_add_variable( vid,                             & ! (out)
            fid, varname, desc, units, dims, size(dims), dtype, & ! (in)
            tint8, itavg,                                       & ! (in)
            error                                               ) ! (out)
       if ( error /= SUCCESS_CODE ) then
          write(*,*) 'xxx failed to add variable'
          call PRC_MPIstop
       end if

       File_vname_list  (File_vid_count) = trim(varname)
       File_vid_list    (File_vid_count) = vid
       File_vid_fid_list(File_vid_count) = fid
       File_vid_count                    = File_vid_count + 1
    endif

    return
  end subroutine FileAddVariable

  !-----------------------------------------------------------------------------
  ! interface File_read
  !-----------------------------------------------------------------------------
  subroutine FileRead1D( &
      var,      & ! (out)
      basename, & ! (in)
      varname,  & ! (in)
      step,     & ! (in)
      single    & ! (in) optional
      )
    implicit none

    real(DP),         intent(out) :: var(:)
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: varname
    integer,          intent( in) :: step
    logical,          intent( in), optional :: single

    integer :: fid
    type(datainfo) :: dinfo 
    integer :: dim_size(1)
    integer :: error
    integer :: n

    intrinsic shape
    !---------------------------------------------------------------------------

    call FileOpen( fid,              & ! (out)
         basename, File_FREAD, single ) ! (in)

    !--- get data information
    call file_get_datainfo( dinfo, & ! (out)
         fid, varname, step,       & ! (in)
         error                     ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to get data information'
       call PRC_MPIstop
    end if

    if ( dinfo%rank /= 1 ) then
       write(*,*) 'xxx rank is not 1', dinfo%rank
       call PRC_MPIstop
    end if
    dim_size(:) = shape(var)
    do n = 1, 1
       if ( dinfo%dim_size(n) /= dim_size(n) ) then
          write(*,*) 'xxx shape is different: ', varname, n, dinfo%dim_size(n), dim_size(n)
          call PRC_MPIstop
       end if
    end do

    call file_read_data( var(:), & ! (out)
         dinfo, DP,              & ! (in)
         error                   ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to get data value'
       call PRC_MPIstop
    end if

    return
  end subroutine FileRead1D
  subroutine FileRead2D( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      step,          & ! (in)
      allow_missing, & ! (in) optional
      single         & ! (in) optional
      )
    implicit none

    real(DP),         intent(out)           :: var(:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)           :: step
    logical,          intent( in), optional :: allow_missing !--- if data is missing, set value to zero
    logical,          intent( in), optional :: single

    integer :: fid
    type(datainfo) :: dinfo 
    integer :: dim_size(2)
    integer :: error
    integer :: n

    intrinsic shape
    !---------------------------------------------------------------------------

    !--- search/register file
    call FileOpen( fid,              & ! (out)
         basename, File_FREAD, single ) ! (in)

    !--- get data information
    call file_get_datainfo( dinfo, & ! (out)
         fid, varname, step,       & ! (in)
         error                     ) ! (out)

    !--- verify
    if ( error /= SUCCESS_CODE ) then
       if ( present(allow_missing) .and. allow_missing) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[File] data not found! : ', &
               'varname= ',trim(varname),', step=',step
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[File] Value is set to 0.'
          var(:,:) = 0.D0
       else
          write(*,*) 'xxx faild to get data information'
          call PRC_MPIstop
       end if
    end if

    if ( dinfo%rank /= 2 ) then
       write(*,*) 'xxx rank is not 2', dinfo%rank
       call PRC_MPIstop
    end if
    dim_size(:) = shape(var)
    do n = 1, 2
       if ( dinfo%dim_size(n) /= dim_size(n) ) then
          write(*,*) 'xxx shape is different: ', varname, n, dinfo%dim_size(n), dim_size(n)
          call PRC_MPIstop
       end if
    end do

    call file_read_data( var(:,:), & ! (out)
         dinfo, DP,                & ! (in)
         error                     ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to get data value'
       call PRC_MPIstop
    end if

    return
  end subroutine FileRead2D
  subroutine FileRead3D( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      step,          & ! (in)
      allow_missing, & ! (in) optional
      single         & ! (in) optional
      )
    implicit none

    real(DP),         intent(out)           :: var(:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)           :: step
    logical,          intent( in), optional :: allow_missing !--- if data is missing, set value to zero
    logical,          intent( in), optional :: single

    integer :: fid
    type(datainfo) :: dinfo 
    integer :: dim_size(3)
    integer :: error
    integer :: n

    intrinsic shape
    !---------------------------------------------------------------------------

    !--- search/register file
    call FileOpen( fid,               & ! (out)
         basename, File_FREAD, single ) ! (in)

    !--- get data information
    call file_get_datainfo( dinfo, & ! (out)
         fid, varname, step,       & ! (in)
         error                     ) ! (out)

    !--- verify
    if ( error /= SUCCESS_CODE ) then
       if ( present(allow_missing) .and. allow_missing) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[File] data not found! : ', &
               'varname= ',trim(varname),', step=',step
          if( IO_L ) write(IO_FID_LOG,*) 'xxx [INPUT]/[File] Value is set to 0.'
          var(:,:,:) = 0.D0
       else
          write(*,*) 'xxx faild to get data information'
          call PRC_MPIstop
       end if
    end if

    if ( dinfo%rank /= 3 ) then
       write(*,*) 'xxx rank is not 3', dinfo%rank
       call PRC_MPIstop
    end if
    dim_size(:) = shape(var)
    do n = 1, 3
       if ( dinfo%dim_size(n) /= dim_size(n) ) then
          write(*,*) 'xxx shape is different: ', varname, n, dinfo%dim_size(n), dim_size(n)
          call PRC_MPIstop
       end if
    end do

    call file_read_data( var(:,:,:), & ! (out)
         dinfo, DP,                  & ! (in)
         error                       ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to get data value'
       call PRC_MPIstop
    end if

    return
  end subroutine FileRead3D

  !-----------------------------------------------------------------------------
  ! interface FileWrite
  !-----------------------------------------------------------------------------
  subroutine FileWrite1D( &
      vid,     & ! (in)
      var,     & ! (in)
      t_start, & ! (in)
      t_end    & ! (in)
      )
    implicit none

    real(DP), intent(in) :: var(:)
    integer,  intent(in) :: vid
    real(DP), intent(in) :: t_start
    real(DP), intent(in) :: t_end

    real(8) :: ts, te

    integer :: error
    !---------------------------------------------------------------------------

    ts = t_start
    te = t_end
    call file_write_data( vid, var(:), ts, te, DP, & ! (in)
         error                                     ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to write data'
       call PRC_MPIstop
    end if

    return
  end subroutine FileWrite1D
  subroutine FileWrite2D( &
      vid,     & ! (in)
      var,     & ! (in)
      t_start, & ! (in)
      t_end    & ! (in)
      )
    implicit none

    real(DP), intent(in) :: var(:,:)
    integer,  intent(in) :: vid
    real(DP), intent(in) :: t_start
    real(DP), intent(in) :: t_end

    real(8) :: ts, te
    integer :: error
    !---------------------------------------------------------------------------

    ts = t_start
    te = t_end
    call file_write_data( vid, var(:,:), ts, te, DP, & ! (in)
         error                                       ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to write data'
       call PRC_MPIstop
    end if

    return
  end subroutine FileWrite2D
  subroutine FileWrite3D( &
      vid,     & ! (in)
      var,     & ! (in)
      t_start, & ! (in)
      t_end    & ! (in)
      )
    implicit none

    real(DP), intent(in) :: var(:,:,:)
    integer,  intent(in) :: vid
    real(DP), intent(in) :: t_start
    real(DP), intent(in) :: t_end

    real(8) :: ts, te
    integer :: error
    !---------------------------------------------------------------------------

    ts = t_start
    te = t_end
    call file_write_data( vid, var(:,:,:), ts, te, DP, & ! (in)
         error                                         ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to write data'
       call PRC_MPIstop
    end if

    return
  end subroutine FileWrite3D

  !-----------------------------------------------------------------------------
  subroutine FileClose( &
       fid & ! (in)
       )
    use mod_process, only: &
       PRC_myrank
    implicit none

    integer, intent(in) :: fid

    character(LEN=File_HLONG) :: fname
    integer                   :: error
    integer                   :: n
    !---------------------------------------------------------------------------

    if ( fid < 0 ) return

    do n = 1, File_fid_count-1
       if ( File_fid_list(n) == fid ) exit
    end do
    if ( fid /= File_fid_list(n) ) then
       write(*,*) 'xxx invalid fid' , fid
       call PRC_MPIstop
    end if
    call file_close( fid , & ! (in)
         error             ) ! (out)
    if ( error == SUCCESS_CODE ) then
       if( IO_L ) write(IO_FID_LOG,'(1x,A,i3)') '*** [File] File Close : NO.', n
       call FileMakeFname(fname,trim(File_bname_list(n)),'pe',PRC_myrank,6)
       if( IO_L ) write(IO_FID_LOG,*) '*** closed filename: ', trim(fname)
    else if ( error /= ALREADY_CLOSED_CODE ) then
       write(*,*) 'xxx failed to close file'
       call PRC_MPIstop
    end if

    do n = 1, File_fid_count-1
       if ( File_fid_list(n) == fid ) then
          File_fid_list(n) = -1
          File_bname_list(n) = ''
       end if
    end do

    return
  end subroutine FileClose
  !-----------------------------------------------------------------------------
  subroutine FileCloseAll
    implicit none

    integer n
    !---------------------------------------------------------------------------

    do n = 1, File_fid_count-1
       call FileClose( File_fid_list(n) )
    enddo

    return
  end subroutine FileCloseAll

  !-----------------------------------------------------------------------------
  ! private
  !-----------------------------------------------------------------------------
  subroutine FileMakeFname( &
       fname,    & ! (out)
       basename, & ! (in)
       prefix,   & ! (in)
       myrank,   & ! (in)
       len       ) ! (in)
    character(len=*), intent(out) :: fname
    character(len=*), intent( in) :: basename
    character(len=*), intent( in) :: prefix
    integer,          intent( in) :: myrank
    integer,          intent( in) :: len

    !                           12345678901234567
    character(len=17) :: fmt = "(A, '.', A, I*.*)"
    integer :: n
    !---------------------------------------------------------------------------

    if ( len < 1 .or. len > 9 ) then
       write(*,*) 'xxx len is invalid'
       call PRC_MPIstop
    end if

    write(fmt(14:14),'(I1)') len
    write(fmt(16:16),'(I1)') len
    write(fname, fmt) trim(basename), trim(prefix), myrank

    return
  end subroutine FileMakeFname
  !-----------------------------------------------------------------------------
  subroutine FileGetfid( &
      fid,        &
      existed,    &
      basename,   &
      mode,       &
      single      )
    use mod_process, only: &
       PRC_master, &
       PRC_myrank
    implicit none

    integer,          intent(out) :: fid
    logical,          intent(out) :: existed
    character(LEN=*), intent( in) :: basename
    integer,          intent( in) :: mode
    logical,          intent( in) :: single


    character(LEN=File_HSHORT) :: rwname(0:2)
    data rwname / 'READ','WRITE','APPEND' /

    character(LEN=File_HLONG) :: fname
    integer                  :: n

    integer :: error
    !---------------------------------------------------------------------------

    !--- search existing file
    fid = -1
    do n = 1, File_fid_count-1
       if ( trim(basename)==trim(File_bname_list(n)) ) fid = File_fid_list(n)
    enddo

    if ( fid >= 0 ) then
       existed = .true.
       return
    end if

    !--- register new file and open
    if ( single ) then
       fname = trim(basename)//'.peall'
    else
       call FileMakeFname(fname,trim(basename),'pe',PRC_myrank,6)
    endif
    
    call file_open( fid, & ! (out)
         fname, mode,    & ! (in)
         error           ) ! (out)
    if ( error /= SUCCESS_CODE ) then
       write(*,*) 'xxx faild to open file'
       call PRC_MPIstop
    end if

    if( IO_L ) write(IO_FID_LOG,*) '*** [File1D] File registration : ',trim(rwname(mode)),' -', fid
    if( IO_L ) write(IO_FID_LOG,*) '*** filename: ', trim(fname)
    
    File_bname_list(File_fid_count) = trim(basename)
    File_fid_list  (File_fid_count) = fid
    File_fid_count                  = File_fid_count + 1

    existed = .false.

    return
  end subroutine FileGetfid

end module gtool_file
!-------------------------------------------------------------------------------
