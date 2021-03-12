!-------------------------------------------------------------------------------
!> module file_h
!!
!! @par Description
!!          header information for FILE module
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_file_h
  use scale_precision, only: &
     DP
  use iso_c_binding
  implicit none
  private

#include "scale_file_const.h"

  !--- character length
  integer, public, parameter :: FILE_HSHORT = File_HSHORT
  integer, public, parameter :: FILE_HMID   = File_HMID 
  integer, public, parameter :: FILE_HLONG  = File_HLONG

  !--- data type
  integer, public, parameter :: FILE_REAL4    = File_REAL4
  integer, public, parameter :: FILE_REAL8    = File_REAL8
  integer, public, parameter :: FILE_INTEGER2 = File_INTEGER2
  integer, public, parameter :: FILE_INTEGER4 = File_INTEGER4
  integer, public, parameter :: FILE_INTEGER8 = File_INTEGER8
  integer, public, parameter :: FILE_TEXT     = File_TEXT

  !--- action type
  integer, public, parameter :: FILE_FREAD   = File_FREAD
  integer, public, parameter :: FILE_FWRITE  = File_FWRITE
  integer, public, parameter :: FILE_FAPPEND = File_FAPPEND

  !--- return codes
  integer, public, parameter :: FILE_ERROR_CODE           = ERROR_CODE
  integer, public, parameter :: FILE_SUCCESS_CODE         = SUCCESS_CODE
  integer, public, parameter :: FILE_ALREADY_CLOSED_CODE  = ALREADY_CLOSED_CODE
  integer, public, parameter :: FILE_ALREADY_EXISTED_CODE = ALREADY_EXISTED_CODE

  !--- max
  integer, public, parameter :: FILE_FILE_MAX = FILE_MAX
  integer, public, parameter :: FILE_VAR_MAX  = VAR_MAX
  integer, public, parameter :: FILE_RANK_MAX = RANK_MAX

  !--- missing value
#define DBL(v) v
  real(DP), public, parameter :: FILE_RMISS = DBL(RMISS)_DP

  !--- struct for data infomation
  type, public, Bind(C) :: datainfo
     character(c_char) :: varname(FILE_HSHORT)
     character(c_char) :: description(FILE_HMID)
     character(c_char) :: units(FILE_HSHORT)
     character(c_char) :: standard_name(FILE_HMID)
     integer           :: datatype
     integer           :: rank
     character(c_char) :: dim_name(FILE_HSHORT,RANK_MAX)
     integer           :: dim_size(RANK_MAX)
     integer           :: step
     real(DP)          :: time_start
     real(DP)          :: time_end
     character(c_char) :: time_units(FILE_HMID)
     character(c_char) :: calendar(FILE_HSHORT)
     integer           :: natts
     character(c_char) :: att_name(FILE_HSHORT,ATT_MAX)
     integer           :: att_type(ATT_MAX)
     integer           :: att_len (ATT_MAX)
     integer           :: fid
  endtype datainfo

  integer, public, parameter :: FILE_preclist(0:3) = (/ 4, 8, 4, 8 /)

  character(len=FILE_HSHORT), public :: FILE_dtypelist(0:4)

  data FILE_dtypelist / "REAL4", "REAL8", "INTEGER2", "INTEGER4", "INTEGER8" /

  public :: cstr
  public :: fstr

  interface fstr
     module procedure fstr1
     module procedure fstr2
  end interface fstr

contains
  

  function cstr(str)
    character(*), intent(in) :: str
    character(:,c_char), allocatable, target :: cstr
    cstr = trim(str)//C_null_char
  end function cstr

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
    do i = 1, len(fstr)
       if ( cstr(i) == c_null_char ) exit
       fstr(i:i) = cstr(i)
    end do
    do j = i, len(fstr)
       fstr(j:j) = " "
    end do
    return
  end subroutine fstr2


end module scale_file_h
