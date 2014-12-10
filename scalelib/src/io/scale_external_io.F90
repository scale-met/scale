!-------------------------------------------------------------------------------
!> module FILE I/O (netcdf)
!!
!! @par Description
!!          general file I/O module
!!          frontend interface of netCDF I/O routine
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-01 (R.Yoshida)   [new]
!!
!<
module scale_external_io
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_file, only: &
     FileMakeFname
  use scale_precision
  use scale_stdio
  use scale_grid_index
  use scale_process, only: &
     PRC_MPIstop
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ExternalFileGetShape
  public :: ExternalFileGetGlobalAttV
  public :: ExternalFileGetGlobalAttC
  public :: ExternalFileVarExistence
  public :: ExternalFileRead
  public :: ExternalFileReadOffset

  interface ExternalFileGetGlobalAttV
     module procedure ExternalFileGetGlobalAttVInteger
     module procedure ExternalFileGetGlobalAttVRealSP
     module procedure ExternalFileGetGlobalAttVRealDP
  end interface ExternalFileGetGlobalAttV

  interface ExternalFileRead
     module procedure ExternalFileRead2DRealSP
     module procedure ExternalFileRead2DRealDP
     module procedure ExternalFileRead3DRealSP
     module procedure ExternalFileRead3DRealDP
     module procedure ExternalFileRead4DRealSP
     module procedure ExternalFileRead4DRealDP
  end interface ExternalFileRead

  interface ExternalFileReadOffset
     !module procedure ExternalFileRead2DRealSP
     !module procedure ExternalFileRead2DRealDP
     module procedure ExternalFileReadOffset3DRealSP
     module procedure ExternalFileReadOffset3DRealDP
     module procedure ExternalFileReadOffset4DRealSP
     module procedure ExternalFileReadOffset4DRealDP
  end interface ExternalFileReadOffset

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ExternalFileMakeFName
  private :: ExternalTakeDimension
  private :: ConvertArrayOrder
  private :: handle_err

  interface ConvertArrayOrder
     module procedure ConvertArrayOrderWRF2DSP
     module procedure ConvertArrayOrderWRF2DDP
     module procedure ConvertArrayOrderWRF3DSP
     module procedure ConvertArrayOrderWRF3DDP
     module procedure ConvertArrayOrderWRF4DSP
     module procedure ConvertArrayOrderWRF4DDP
  end interface ConvertArrayOrder

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: iSCALE  = 1  ! use gtool, coz it's not external
  integer, private, parameter :: iWRFARW = 2
  integer, private, parameter :: iNICAM  = 3
  integer, private, parameter :: iJMAMSM = 4


  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  ! ExternalFileGetShape
  !-----------------------------------------------------------------------------
  subroutine ExternalFileGetShape( &
      dims,          & ! (out)
      timelen,       & ! (out)
      mdlid,         & ! (in)
      basename,      & ! (in)
      myrank,        & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    integer,          intent(out)           :: dims(:)
    integer,          intent(out)           :: timelen
    integer,          intent( in)           :: mdlid
    character(LEN=*), intent( in)           :: basename
    integer,          intent( in)           :: myrank
    logical,          intent( in), optional :: single

    integer :: status
    integer :: ncid, unlimid
    integer :: dims_org(7)
    character(len=NF90_MAX_NAME) :: tname
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire( ncid, unlimitedDimID=unlimid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_dimension( ncid, unlimid, tname, timelen )
    if (status .ne. nf90_noerr) call handle_err(status)

    if( trim(tname)=='time' .or. trim(tname)=='Time' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'Time Dimension Name: '//trim(tname)
    else
       write(*,*) 'xxx Not appropriate time dimension is used in the external file. Check!'
       call PRC_MPIstop
    endif

    call ExternalTakeDimension( dims_org(:),ncid,mdlid )

    ! convert dimension order for return to scale-system
    if( mdlid == iWRFARW )then   !MODEL ID: WRF-ARW
       dims(1) = dims_org(3)
       dims(2) = dims_org(1)
       dims(3) = dims_org(2)
       dims(4) = dims_org(6)
       dims(5) = dims_org(4)
       dims(6) = dims_org(5)
       dims(7) = dims_org(7)
    else
       dims(:) = dims_org(:)
    endif

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    return
  end subroutine ExternalFileGetShape

  !-----------------------------------------------------------------------------
  ! ExternalFileGet Global Attribute (value, real, single precision)
  !-----------------------------------------------------------------------------
  subroutine ExternalFileGetGlobalAttVInteger( &
      var,           & ! (out)
      mdlid,         & ! (in)
      basename,      & ! (in)
      attname,       & ! (in)
      myrank,        & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    integer,          intent(out)           :: var(:)
    integer,          intent( in)           :: mdlid
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: attname
    integer,          intent( in)           :: myrank
    logical,          intent( in), optional :: single

    integer, allocatable :: work(:)

    integer :: status
    integer :: i, ncid, length
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_attribute(ncid, nf90_global, trim(attname), len=length)
    if (status .ne. nf90_noerr) call handle_err(status)

    allocate( work(length) )

    status = nf90_get_att(ncid, nf90_global, trim(attname), work)
    if (status .ne. nf90_noerr) call handle_err(status)

    do i = 1, length
       var(i) = work(i)
    enddo

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)
    deallocate( work )

    return
  end subroutine ExternalFileGetGlobalAttVInteger

  !-----------------------------------------------------------------------------
  ! ExternalFileGet Global Attribute (value, real, single precision)
  !-----------------------------------------------------------------------------
  subroutine ExternalFileGetGlobalAttVRealSP( &
      var,           & ! (out)
      mdlid,         & ! (in)
      basename,      & ! (in)
      attname,       & ! (in)
      myrank,        & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(SP),         intent(out)           :: var(:)
    integer,          intent( in)           :: mdlid
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: attname
    integer,          intent( in)           :: myrank
    logical,          intent( in), optional :: single

    real(SP), allocatable :: work(:)

    integer :: status
    integer :: i, ncid, length
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_attribute(ncid, nf90_global, trim(attname), len=length)
    if (status .ne. nf90_noerr) call handle_err(status)

    allocate( work(length) )

    status = nf90_get_att(ncid, nf90_global, trim(attname), work)
    if (status .ne. nf90_noerr) call handle_err(status)

    do i = 1, length
       var(i) = work(i)
    enddo

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)
    deallocate( work )

    return
  end subroutine ExternalFileGetGlobalAttVRealSP

  !-----------------------------------------------------------------------------
  ! ExternalFileGet Global Attribute (value, real, single precision)
  !-----------------------------------------------------------------------------
  subroutine ExternalFileGetGlobalAttVRealDP( &
      var,           & ! (out)
      mdlid,         & ! (in)
      basename,      & ! (in)
      attname,       & ! (in)
      myrank,        & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(DP),         intent(out)           :: var(:)
    integer,          intent( in)           :: mdlid
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: attname
    integer,          intent( in)           :: myrank
    logical,          intent( in), optional :: single

    real(DP), allocatable :: work(:)

    integer :: status
    integer :: i, ncid, length
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_attribute(ncid, nf90_global, trim(attname), len=length)
    if (status .ne. nf90_noerr) call handle_err(status)

    allocate( work(length) )

    status = nf90_get_att(ncid, nf90_global, trim(attname), work)
    if (status .ne. nf90_noerr) call handle_err(status)

    do i = 1, length
       var(i) = work(i)
    enddo

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)
    deallocate( work )

    return
  end subroutine ExternalFileGetGlobalAttVRealDP

  !-----------------------------------------------------------------------------
  ! ExternalFileGet Global Attribute (character)
  !-----------------------------------------------------------------------------
  subroutine ExternalFileGetGlobalAttC( &
      chr,           & ! (out)
      mdlid,         & ! (in)
      basename,      & ! (in)
      attname,       & ! (in)
      myrank,        & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    character(LEN=*), intent(out)           :: chr(:)
    integer,          intent( in)           :: mdlid
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: attname
    integer,          intent( in)           :: myrank
    logical,          intent( in), optional :: single

    integer :: status
    integer :: ncid, length
    character(len=H_LONG) :: fname = ''
    character(len=80)     :: work
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_attribute(ncid, nf90_global, trim(attname), len=length)
    if (status .ne. nf90_noerr) call handle_err(status)

    if( len(work) < length ) then
       write(*,*) 'xxx Not enough space to put attribute values. [externalio/scalelib]'
       call PRC_MPIstop
    endif

    status = nf90_get_att(ncid, nf90_global, trim(attname), work)
    if (status .ne. nf90_noerr) call handle_err(status)

    chr = trim(work)

    return
  end subroutine ExternalFileGetGlobalAttC

  !-----------------------------------------------------------------------------
  !> Check Existence of a Variable
  !-----------------------------------------------------------------------------
  subroutine ExternalFileVarExistence( &
      existence,     & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    logical,          intent(out)           :: existence
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)           :: myrank
    integer,          intent( in)           :: mdlid
    logical,          intent( in), optional :: single

    integer :: ncid, varid
    integer :: status

    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .eq. nf90_noerr) then
       existence = .true.
    else
       existence = .false.
       if( IO_L ) write(IO_FID_LOG,*) '+++ not exist variable: ', trim(varname)
    endif

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    return
  end subroutine ExternalFileVarExistence

  !-----------------------------------------------------------------------------
  !> File Read
  !-----------------------------------------------------------------------------
  subroutine ExternalFileRead2DRealSP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,          & ! (in)
      te,        & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      nx,            & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(SP),         intent(out)            :: var(:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    integer,          intent( in)            :: nx
    logical,          intent( in), optional :: single

    real(SP), allocatable :: var_org(:,:)
    integer :: ncid, varid
    integer :: status
    integer :: precis

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! based on specified dimension size
    allocate( var_org(nx,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_FLOAT) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileRead2DSP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, var_org(:,:), start = (/ 1,ts /), &
                            count = (/ nx,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nx )

    deallocate( var_org )

    return
  end subroutine ExternalFileRead2DRealSP
  subroutine ExternalFileRead2DRealDP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,          & ! (in)
      te,        & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      nx,            & ! (in)
      single         & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(DP),         intent(out)            :: var(:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    integer,          intent( in)            :: nx
    logical,          intent( in), optional :: single

    real(DP), allocatable :: var_org(:,:)
    integer :: ncid, varid
    integer :: status
    integer :: precis

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! based on specified dimension size
    allocate( var_org(nx,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_DOUBLE) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileRead2DDP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, var_org(:,:), start = (/ 1,ts /), &
                            count = (/ nx,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nx )

    deallocate( var_org )

    return
  end subroutine ExternalFileRead2DRealDP
  subroutine ExternalFileRead3DRealSP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,          & ! (in)
      te,        & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag          & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(SP),         intent(out)            :: var(:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag

    real(SP), allocatable :: var_org(:,:,:)
    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    allocate( var_org(nx,ny,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_FLOAT) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileRead3DSP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, var_org(:,:,:), start = (/ 1,1,ts /), &
                            count = (/ nx,ny,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nx,ny )

    deallocate( var_org )

    return
  end subroutine ExternalFileRead3DRealSP
  subroutine ExternalFileRead3DRealDP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,          & ! (in)
      te,        & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag          & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(DP),         intent(out)            :: var(:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag

    real(DP), allocatable :: var_org(:,:,:)
    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    allocate( var_org(nx,ny,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_DOUBLE) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileRead3DDP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, var_org(:,:,:), start = (/ 1,1,ts /), &
                            count = (/ nx,ny,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nx,ny )

    deallocate( var_org )

    return
  end subroutine ExternalFileRead3DRealDP
  subroutine ExternalFileRead4DRealSP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,          & ! (in)
      te,        & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag,         & ! (in) optional
      zstag,         & ! (in) optional
      landgrid       & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(SP),         intent(out)            :: var(:,:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag
    logical,          intent( in), optional :: zstag
    logical,          intent( in), optional :: landgrid

    real(SP), allocatable :: var_org(:,:,:,:)
    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny, nz
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    nz = dims(3)
    if ( present(zstag) .and. zstag ) then
       nz = dims(6)
    endif
    if ( present(landgrid) .and. landgrid ) then
       nz = dims(7)
    endif
    allocate( var_org(nx,ny,nz,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_FLOAT) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileRead4DSP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, var_org(:,:,:,:), start = (/ 1,1,1,ts /), &
                            count = (/ nx,ny,nz,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nz,nx,ny )

    deallocate( var_org )

    return
  end subroutine ExternalFileRead4DRealSP
  subroutine ExternalFileRead4DRealDP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,          & ! (in)
      te,        & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag,         & ! (in) optional
      zstag,         & ! (in) optional
      landgrid       & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(DP),         intent(out)            :: var(:,:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag
    logical,          intent( in), optional :: zstag
    logical,          intent( in), optional :: landgrid

    real(DP), allocatable :: var_org(:,:,:,:)
    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny, nz
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    nz = dims(3)
    if ( present(zstag) .and. zstag ) then
       nz = dims(6)
    endif
    if ( present(landgrid) .and. landgrid ) then
       nz = dims(7)
    endif
    allocate( var_org(nx,ny,nz,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_DOUBLE) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileRead4DDP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, var_org(:,:,:,:), start = (/ 1,1,1,ts /), &
                            count = (/ nx,ny,nz,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nz,nx,ny )

    deallocate( var_org )

    return
  end subroutine ExternalFileRead4DRealDP

  subroutine ExternalFileReadOffset3DRealSP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,            & ! (in)
      te,            & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag          & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(SP),         intent(out)            :: var(:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag

    real(SP),   allocatable :: var_org(:,:,:)
    integer(2), allocatable :: short(:,:,:)

    real(4) :: scale_factor, add_offset

    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    allocate( var_org(nx,ny,tcount) )
    allocate( short  (nx,ny,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "scale_factor", scale_factor)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "add_offset", add_offset)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_SHORT) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileReadOffset4DSP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, short(:,:,:), start = (/ 1,1,ts /), &
                            count = (/ nx,ny,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    var_org(:,:,:) = real( short(:,:,:) )*scale_factor + add_offset

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nx,ny )

    deallocate( var_org )
    deallocate( short )

    return
  end subroutine ExternalFileReadOffset3DRealSP
  subroutine ExternalFileReadOffset3DRealDP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,            & ! (in)
      te,            & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag          & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(DP),         intent(out)            :: var(:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag

    real(DP),   allocatable :: var_org(:,:,:)
    integer(2), allocatable :: short(:,:,:)

    real(4) :: scale_factor, add_offset

    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    allocate( var_org(nx,ny,tcount) )
    allocate( short  (nx,ny,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "scale_factor", scale_factor)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "add_offset", add_offset)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_SHORT) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileReadOffset4DSP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, short(:,:,:), start = (/ 1,1,ts /), &
                            count = (/ nx,ny,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    var_org(:,:,:) =  real( real(short(:,:,:),kind=SP)*scale_factor + add_offset, kind=DP )

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nx,ny )

    deallocate( var_org )
    deallocate( short )

    return
  end subroutine ExternalFileReadOffset3DRealDP
  subroutine ExternalFileReadOffset4DRealSP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,            & ! (in)
      te,            & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag,         & ! (in) optional
      zstag,         & ! (in) optional
      landgrid       & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(SP),         intent(out)            :: var(:,:,:,:)
    character(LEN=*), intent( in)           :: basename
    character(LEN=*), intent( in)           :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag
    logical,          intent( in), optional :: zstag
    logical,          intent( in), optional :: landgrid

    real(SP),   allocatable :: var_org(:,:,:,:)
    integer(2), allocatable :: short(:,:,:,:)

    real(4) :: scale_factor, add_offset

    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny, nz
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    nz = dims(3)
    if ( present(zstag) .and. zstag ) then
       nz = dims(6)
    endif
    if ( present(landgrid) .and. landgrid ) then
       nz = dims(7)
    endif
    allocate( var_org(nx,ny,nz,tcount) )
    allocate( short  (nx,ny,nz,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "scale_factor", scale_factor)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "add_offset", add_offset)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_SHORT) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileReadOffset4DSP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, short(:,:,:,:), start = (/ 1,1,1,ts /), &
                            count = (/ nx,ny,nz,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    var_org(:,:,:,:) = real( short(:,:,:,:) )*scale_factor + add_offset

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nz,nx,ny )

    deallocate( var_org )
    deallocate( short )

    return
  end subroutine ExternalFileReadOffset4DRealSP
  subroutine ExternalFileReadOffset4DRealDP( &
      var,           & ! (out)
      basename,      & ! (in)
      varname,       & ! (in)
      ts,            & ! (in)
      te,            & ! (in)
      myrank,        & ! (in)
      mdlid,         & ! (in)
      single,        & ! (in) optional
      xstag,         & ! (in) optional
      ystag,         & ! (in) optional
      zstag,         & ! (in) optional
      landgrid       & ! (in) optional
      )
    use netcdf  ![external lib]
    implicit none

    real(DP),         intent(out)            :: var(:,:,:,:)
    character(LEN=*), intent( in)            :: basename
    character(LEN=*), intent( in)            :: varname
    integer,          intent( in)            :: ts
    integer,          intent( in)            :: te
    integer,          intent( in)            :: myrank
    integer,          intent( in)            :: mdlid
    logical,          intent( in), optional :: single
    logical,          intent( in), optional :: xstag
    logical,          intent( in), optional :: ystag
    logical,          intent( in), optional :: zstag
    logical,          intent( in), optional :: landgrid

    real(DP),   allocatable :: var_org(:,:,:,:)
    integer(2), allocatable :: short(:,:,:,:)

    real(4) :: scale_factor, add_offset

    integer :: ncid, varid
    integer :: status
    integer :: precis
    integer :: nx, ny, nz
    integer :: dims(7)

    integer :: tcount
    character(len=H_LONG) :: fname = ''
    logical :: single_ = .false.

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    tcount = te - ts + 1

    if ( present(single) ) then
       single_ = single
    else
       single_ = .false.
    endif

    call ExternalFileMakeFname( fname,mdlid,basename,myrank,single_ )

    status = nf90_open( trim(fname), nf90_nowrite, ncid )
    if (status .ne. nf90_noerr) call handle_err(status)

    ! retrieve dimension size in data original order
    call ExternalTakeDimension( dims(:),ncid,mdlid )
    nx = dims(1)
    if ( present(xstag) .and. xstag ) then
       nx = dims(4)
    endif
    ny = dims(2)
    if ( present(ystag) .and. ystag ) then
       ny = dims(5)
    endif
    nz = dims(3)
    if ( present(zstag) .and. zstag ) then
       nz = dims(6)
    endif
    if ( present(landgrid) .and. landgrid ) then
       nz = dims(7)
    endif
    allocate( var_org(nx,ny,nz,tcount) )
    allocate( short  (nx,ny,nz,tcount) )

    status = nf90_inq_varid( ncid, trim(varname), varid )
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "scale_factor", scale_factor)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_get_att(ncid, varid, "add_offset", add_offset)
    if (status .ne. nf90_noerr) call handle_err(status)

    status = nf90_inquire_variable( ncid, varid, xtype=precis )
    if(status /= nf90_NoErr) call handle_err(status)
    if(precis /= NF90_SHORT) then
       write(*,*) 'xxx Internal Error: [scale_external_io]/[ExternalFileReadOffset4DSP]'
       call PRC_MPIstop
    endif

    status = nf90_get_var( ncid, varid, short(:,:,:,:), start = (/ 1,1,1,ts /), &
                            count = (/ nx,ny,nz,tcount /) )
    if (status .ne. nf90_noerr) call handle_err(status)

    var_org(:,:,:,:) = real( real(short(:,:,:,:),kind=SP)*scale_factor + add_offset, kind=DP )

    status = nf90_close(ncid)
    if (status .ne. nf90_noerr) call handle_err(status)

    call ConvertArrayOrder( var,var_org,tcount,nz,nx,ny )

    deallocate( var_org )
    deallocate( short )

    return
  end subroutine ExternalFileReadOffset4DRealDP

  !-----------------------------------------------------------------------------
  ! ExternalMakeFName
  !-----------------------------------------------------------------------------
  subroutine ExternalFileMakeFName( &
      fname,         & ! (out)
      mdlid,         & ! (in)
      basename,      & ! (in)
      myrank,        & ! (in)
      single         & ! (in) optional
      )
    implicit none

    character(LEN=*), intent(out) :: fname
    integer,          intent( in)  :: mdlid
    character(LEN=*), intent( in) :: basename
    integer,          intent( in)  :: myrank
    logical,          intent( in)  :: single

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if( mdlid == iWRFARW )then          !TYPE: WRF-ARW
       if ( single ) then
          fname = trim(basename)
       else
          call FileMakeFname(fname,trim(basename),'_',myrank,4)
       endif
    elseif( mdlid == iNICAM )then      !TYPE: NICAM-NETCDF
       if ( single ) then
          fname = trim(basename)//'.peall.nc'
       else
          call FileMakeFname(fname,trim(basename),'anl.pe',myrank,6)
       endif
    !elseif( mdlid == iJMAMSM )then      !TYPE: JMA-MSM
    !   if ( single ) then
    !      fname = trim(basename)//'.anl'
    !   else
    !      call FileMakeFname(fname,trim(basename),'anl.pe',myrank,6)
    !   endif
    else
       write(*,*) 'xxx failed, wrong filetype: [scale_external_io]/[ExternalFileMakeFName]'
       call PRC_MPIstop
    endif

    return
  end subroutine ExternalFileMakeFName

  !-----------------------------------------------------------------------------
  ! ExternalMakeFName
  !-----------------------------------------------------------------------------
  subroutine ExternalTakeDimension( &
      dims,          & ! (out)
      ncid,          & ! (in)
      mdlid          & ! (in)
      )
    use netcdf  ![external lib]
    implicit none

    integer,          intent(out)  :: dims(:)
    integer,          intent( in)  :: ncid
    integer,          intent( in)  :: mdlid

    integer :: dimid
    integer :: status

    intrinsic size
    intrinsic shape
    !---------------------------------------------------------------------------

    if( mdlid == iWRFARW )then   !MODEL ID: WRF-ARW
       status = nf90_inq_dimid( ncid, "west_east", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(1) )
       status = nf90_inq_dimid( ncid, "south_north", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(2) )
       status = nf90_inq_dimid( ncid, "bottom_top", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(3) )
       status = nf90_inq_dimid( ncid, "west_east_stag", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(4) )
       status = nf90_inq_dimid( ncid, "south_north_stag", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(5) )
       status = nf90_inq_dimid( ncid, "bottom_top_stag", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(6) )
       status = nf90_inq_dimid( ncid, "soil_layers_stag", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(7) )

    elseif( mdlid == iNICAM )then   !MODEL ID: NICAM-NETCDF
       status = nf90_inq_dimid( ncid, "lon", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(1) )
       status = nf90_inq_dimid( ncid, "lat", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(2) )
       status = nf90_inq_dimid( ncid, "lev", dimid )
       status = nf90_inquire_dimension( ncid, dimid,len=dims(3) )

    else
       write(*,*) 'xxx This external file format is not supported, Sorry.'
       call PRC_MPIstop
    endif

    return
  end subroutine ExternalTakeDimension

  !-----------------------------------------------------------------------------
  subroutine ConvertArrayOrderWRF2DSP( &
      var,           & ! (out)
      var_org,       & ! (in)
      tcount,        & ! (in)
      nx             & ! (in)
      )
    implicit none

    real(SP),         intent(out)  :: var(:,:)
    real(SP),         intent( in)  :: var_org(:,:)
    integer,          intent( in)  :: tcount
    integer,          intent( in)  :: nx
    integer :: n, i
    intrinsic shape

    do n = 1, tcount
    do i = 1, nx
       var(i,n) = var_org(i,n)
    end do
    end do

    return
  end subroutine ConvertArrayOrderWRF2DSP
  subroutine ConvertArrayOrderWRF2DDP( &
      var,           & ! (out)
      var_org,       & ! (in)
      tcount,        & ! (in)
      nx             & ! (in)
      )
    implicit none

    real(DP),         intent(out)  :: var(:,:)
    real(DP),         intent( in)  :: var_org(:,:)
    integer,          intent( in)  :: tcount
    integer,          intent( in)  :: nx
    integer :: n, i
    intrinsic shape

    do n = 1, tcount
    do i = 1, nx
       var(i,n) = var_org(i,n)
    end do
    end do

    return
  end subroutine ConvertArrayOrderWRF2DDP
  subroutine ConvertArrayOrderWRF3DSP( &
      var,           & ! (out)
      var_org,       & ! (in)
      tcount,        & ! (in)
      nx,            & ! (in)
      ny             & ! (in)
      )
    implicit none

    real(SP),         intent(out)  :: var(:,:,:)
    real(SP),         intent( in)  :: var_org(:,:,:)
    integer,          intent( in)  :: tcount
    integer,          intent( in)  :: nx
    integer,          intent( in)  :: ny
    integer :: n, i, j
    intrinsic shape

    do n = 1, tcount
    do j = 1, ny
    do i = 1, nx
       var(i,j,n) = var_org(i,j,n)
    end do
    end do
    end do

    return
  end subroutine ConvertArrayOrderWRF3DSP
  subroutine ConvertArrayOrderWRF3DDP( &
      var,           & ! (out)
      var_org,       & ! (in)
      tcount,        & ! (in)
      nx,            & ! (in)
      ny             & ! (in)
      )
    implicit none

    real(DP),         intent(out)  :: var(:,:,:)
    real(DP),         intent( in)  :: var_org(:,:,:)
    integer,          intent( in)  :: tcount
    integer,          intent( in)  :: nx
    integer,          intent( in)  :: ny
    integer :: n, i, j
    intrinsic shape

    do n = 1, tcount
    do j = 1, ny
    do i = 1, nx
       var(i,j,n) = var_org(i,j,n)
    end do
    end do
    end do

    return
  end subroutine ConvertArrayOrderWRF3DDP
  subroutine ConvertArrayOrderWRF4DSP( &
      var,           & ! (out)
      var_org,       & ! (in)
      tcount,        & ! (in)
      nz,            & ! (in)
      nx,            & ! (in)
      ny             & ! (in)
      )
    implicit none

    real(SP),         intent(out)  :: var(:,:,:,:)
    real(SP),         intent( in)  :: var_org(:,:,:,:)
    integer,          intent( in)  :: tcount
    integer,          intent( in)  :: nz
    integer,          intent( in)  :: nx
    integer,          intent( in)  :: ny
    integer :: n, k, i, j
    intrinsic shape

    do n = 1, tcount
    do j = 1, ny
    do i = 1, nx
    do k = 1, nz
       var(k,i,j,n) = var_org(i,j,k,n)
    end do
    end do
    end do
    end do

    return
  end subroutine ConvertArrayOrderWRF4DSP
  subroutine ConvertArrayOrderWRF4DDP( &
      var,           & ! (out)
      var_org,       & ! (in)
      tcount,        & ! (in)
      nz,            & ! (in)
      nx,            & ! (in)
      ny             & ! (in)
      )
    implicit none

    real(DP),         intent(out)  :: var(:,:,:,:)
    real(DP),         intent( in)  :: var_org(:,:,:,:)
    integer,          intent( in)  :: tcount
    integer,          intent( in)  :: nz
    integer,          intent( in)  :: nx
    integer,          intent( in)  :: ny
    integer :: n, k, i, j
    intrinsic shape

    do n = 1, tcount
    do j = 1, ny
    do i = 1, nx
    do k = 1, nz
       var(k,i,j,n) = var_org(i,j,k,n)
    end do
    end do
    end do
    end do

    return
  end subroutine ConvertArrayOrderWRF4DDP

  !-----------------------------------------------------------------------------
  subroutine handle_err(status)
    use netcdf  ![external lib]
    implicit none
    integer, intent(in) :: status

    write(*,*) nf90_strerror(status)
    stop

    return
  end subroutine handle_err


end module scale_external_io
!-------------------------------------------------------------------------------
