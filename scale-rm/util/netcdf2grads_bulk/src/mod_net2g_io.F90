!-------------------------------------------------------------------------------------------
!> module NET2G io
!!
!! @par Description
!!          data output module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-07-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_io
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_net2g_vars
  use mod_net2g_error

  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: io_log_init
  public :: io_close_all
  public :: io_create_ctl
  public :: io_write_vars

  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: io_set_fname
  private :: io_control_output_fid

  !-----------------------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer :: nfile = 0
  logical, private :: open_file = .false.

  !-----------------------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------------------

  !> open logfile
  !-----------------------------------------------------------------------------------------
  subroutine io_log_init( &
      irank,   & ! [in   ]
      LOUT     ) ! [inout]
    implicit none

    integer, intent(in)    :: irank
    logical, intent(inout) :: LOUT

    character(6) :: num
    !---------------------------------------------------------------------------

    LOUT = .false.
    if ( LOG_ALL_OUTPUT ) then
       LOUT = .true.
    else
       if ( irank == master ) LOUT = .true.
    endif

    if ( trim(LOG_BASENAME) .eq. "STDOUT" ) then
       FID_LOG = FID_STD
       open_file = .false.
    else
       FID_LOG = FID_LOGF
       open_file = .true.
    endif

    if ( LOG_LEVEL > 0 ) then
       LOG_DBUG = .false.
    else
       LOG_DBUG = .true.
    endif

    if ( open_file .and. LOUT ) then
       write( num,'(I6.6)' ) irank
       open ( FID_LOG, file=trim(LOG_BASENAME)//".pe"//num, &
              status='replace', form='formatted' )
    endif

    return
  end subroutine io_log_init


  !> close all files for binary data output
  !-----------------------------------------------------------------------------------------
  subroutine io_close_all()
    implicit none

    integer :: i, fid
    !---------------------------------------------------------------------------

    do i=1, nfile
       fid = FID_DAT + i
       close( fid )
    enddo

    if ( open_file .and. LOUT ) close ( FID_LOG )

    return
  end subroutine io_close_all


  !> create control file
  !-----------------------------------------------------------------------------------------
  subroutine io_create_ctl( &
      varname,    & ! [in]
      atype,      & ! [in]
      ctype,      & ! [in]
      vtype,      & ! [in]
      idom,       & ! [in]
      nx, ny,     & ! [in]
      zz,         & ! [in]
      nt,         & ! [in]
      cx, cy,     & ! [in]
      vgrid,      & ! [in]
      zlev        ) ! [in]
    implicit none

    character(CMID), intent(in) :: varname
    integer,         intent(in) :: atype, ctype, vtype
    integer,         intent(in) :: idom
    integer,         intent(in) :: nx, ny, zz, nt
    real(SP), intent(in) :: cx(:), cy(:)
    real(SP), intent(in) :: vgrid(:), zlev(:)

    character(2)    :: cdom
    character(3)    :: clev
    character(CLNG) :: fname
    character(CLNG) :: fname2
    !---------------------------------------------------------------------------

     write(cdom,'(i2.2)') idom
     select case( atype )
     case ( a_slice, a_conv )
        if ( Z_MERGE_OUT ) then
           clev = "-3d"
        else
           write(clev,'(i3.3)') zz
        endif
     case ( a_max )
        clev = "max"
     case ( a_min )
        clev = "min"
     case ( a_sum )
        clev = "sum"
     case ( a_ave )
        clev = "ave"
     end select
     if ( vtype == vt_2d ) clev = "-2d"
    if ( T_MERGE_OUT ) then
        fname = trim(ODIR)//'/'//trim(varname)//'_d'//cdom//'z'//clev
        fname2 = trim(varname)//'_d'//cdom//'z'//clev
     else
        fname = trim(ODIR)//'/'//trim(varname)//'_d'//cdom//'z'//clev//'_'//trim(FTIME)
        fname2 = trim(varname)//'_d'//cdom//'z'//clev//'_'//trim(FTIME)
     endif
     if ( LOUT .and. LOG_DBUG ) write( FID_LOG, '(1X,A,A)') "Create ctl file: ", trim(fname)
     open ( FID_CTL, file=trim(fname)//".ctl", form="formatted", access="sequential" )

     write( FID_CTL, '(a,1x,a)') "DSET", "^"//trim(fname2)//".grd"
     write( FID_CTL, '(a)') "TITLE SCALE3 data output"
     write( FID_CTL, '(a)') "OPTIONS BIG_ENDIAN"
     write( FID_CTL, '(a,1x,ES15.7)') "UNDEF", -9.9999001E+30
     write( FID_CTL, '(a,3x,i7,1x,a)') "XDEF", nx, "LEVELS"
     write( FID_CTL, '(5(1x,ES15.7))') cx(1:nx)*1.d-3
     write( FID_CTL, '(a,3x,i7,1x,a)') "YDEF", ny, "LEVELS"
     write( FID_CTL, '(5(1x,ES15.7))') cy(1:ny)*1.d-3

     select case( vtype )
     case ( vt_urban, vt_land, vt_tpmsk )
        write( FID_CTL, '(a,3x,i7,1x,a,1x,a)') "ZDEF", 1, "linear", "1 1"
     case ( vt_2d )
        write( FID_CTL, '(a,3x,i7,1x,a,1x,ES15.7)') "ZDEF", 1, "LEVELS", zlev(zz)*1.d-3
     case ( vt_3d, vt_height )
        if ( Z_MERGE_OUT ) then
           write( FID_CTL, '(a,3x,i7,1x,a)') "ZDEF", ZCOUNT, "LEVELS"
           if ( atype == a_conv .and. ctype == c_pres ) then
              write( FID_CTL, '(5(1x,ES15.7))') vgrid(1:ZCOUNT)
           else
              write( FID_CTL, '(5(1x,ES15.7))') vgrid(1:ZCOUNT)*1.d-3
           endif
        else
           write( FID_CTL, '(a,3x,i7,1x,a,1x,ES15.7)') "ZDEF", 1, "LEVELS", zlev(zz)*1.d-3
        endif

     end select

     write( FID_CTL, '(a,3x,i5,1x,a,1x,a,3x,a)') "TDEF", nt, "LINEAR", trim(STIME), trim(DELT)
     write( FID_CTL, '(a,3x,i2)') "VARS", 1
     if ( Z_MERGE_OUT ) then
        write( FID_CTL, '(a,1x,i7,1x,i2,1x,a)') trim(varname), ZCOUNT, 99, "NONE"
     else
        write( FID_CTL, '(a,1x,i7,1x,i2,1x,a)') trim(varname), 0, 99, "NONE"
     endif
     write( FID_CTL, '(a)') "ENDVARS"

     close( FID_CTL )
    return
  end subroutine io_create_ctl


  !> write data file
  !-----------------------------------------------------------------------------------------
  subroutine io_write_vars( &
      var_2d,      & ! [in]
      varname,     & ! [in]
      atype,       & ! [in]
      vtype,       & ! [in]
      idom,        & ! [in]
      nx, ny, zz,  & ! [in]
      irec         ) ! [in]
    implicit none

    real(SP),        intent(in) :: var_2d(:,:)
    character(CMID), intent(in) :: varname
    integer,         intent(in) :: atype, vtype
    integer,         intent(in) :: idom
    integer,         intent(in) :: nx, ny, zz
    integer,         intent(in) :: irec

    integer         :: fid
    character(CLNG) :: fname
    !---------------------------------------------------------------------------

    call io_set_fname( varname, idom, zz, atype, vtype, fname )
    call io_control_output_fid( fname, nx, ny, fid )
    if ( LOUT ) write( FID_LOG, '(1X,A,A)') "+++ Output data file: ", trim(fname)
    if ( LOUT .and. LOG_DBUG ) write( FID_LOG, '(1X,A,I5,A,I4)') "+++ Data record num:", irec, "  FID:", fid
    if ( LOUT .and. LOG_DBUG ) write( FID_LOG, *) "+++ Check data range: ", maxval(var_2d), minval(var_2d)

    write( fid, rec=irec ) var_2d(:,:)

    return
  end subroutine io_write_vars


  !> set file name for binary data output
  !-----------------------------------------------------------------------------------------
  subroutine io_set_fname( &
      varname,  & ! [in]
      idom,     & ! [in]
      zz,       & ! [in]
      atype,    & ! [in]
      vtype,    & ! [in]
      fname     ) ! [out]
    implicit none

    character(CMID), intent(in)  :: varname
    integer,         intent(in)  :: idom, zz
    integer,         intent(in)  :: atype, vtype
    character(CLNG), intent(out) :: fname

    character(2) :: cdom
    character(3) :: clev
    !---------------------------------------------------------------------------

    write(cdom,'(i2.2)') idom
    select case( atype )
    case ( a_slice, a_conv )
        if ( Z_MERGE_OUT ) then
           clev = "-3d"
        else
           write(clev,'(i3.3)') zz
        endif
    case ( a_max )
       clev = "max"
    case ( a_min )
       clev = "min"
    case ( a_sum )
       clev = "sum"
    case ( a_ave )
       clev = "ave"
    end select
    if ( vtype == vt_2d ) clev = "-2d"

    if ( T_MERGE_OUT ) then
       fname = trim(ODIR)//'/'//trim(varname)//'_d'//cdom//'z'//clev//".grd"
    else
       fname = trim(ODIR)//'/'//trim(varname)//'_d'//cdom//'z'//clev//'_'//trim(FTIME)//".grd"
    endif

    return
  end subroutine io_set_fname


  !> control fid for binary data output
  !-----------------------------------------------------------------------------------------
  subroutine io_control_output_fid( &
      fname,    & ! [in]
      nx, ny,   & ! [in]
      fid       ) ! [out]
    implicit none

    character(CLNG), intent(in)  :: fname
    integer,         intent(in)  :: nx, ny
    integer,         intent(out) :: fid

    integer    :: i
    integer(8) :: irecl
    logical    :: existence
    !---------------------------------------------------------------------------

    existence = .false.
    if ( nfile > 0 ) then
       do i=1, nfile
          if ( trim(fname_bank(i)) .eq. trim(fname) ) then
             existence = .true.
             exit
          endif
       enddo
    else
       i = 1
    endif

    if ( existence ) then
       fid = FID_DAT + i
    else
       fid = FID_DAT + i
       fname_bank(i) = trim(fname)
       nfile = nfile + 1

       irecl = int(nx,kind=8) * int(ny,kind=8) * 4_8
       open( fid, file=trim(fname), form="unformatted", access="direct", recl=irecl)
    endif

    return
  end subroutine io_control_output_fid

end module mod_net2g_io
