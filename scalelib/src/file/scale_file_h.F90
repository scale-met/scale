!-------------------------------------------------------------------------------
!> module file_h
!!
!! @par Description
!!          header information for FILE module
!!
!! @author Team SCALE
!!
!<
module scale_file_h
  use scale_precision, only: &
     DP
  implicit none

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

  !--- action type
  integer, public, parameter :: FILE_FREAD   = File_FREAD
  integer, public, parameter :: FILE_FWRITE  = File_FWRITE
  integer, public, parameter :: FILE_FAPPEND = File_FAPPEND

  !--- return codes
  integer, public, parameter :: FILE_ERROR_CODE           = ERROR_CODE
  integer, public, parameter :: FILE_SUCCESS_CODE         = SUCCESS_CODE
  integer, public, parameter :: FILE_ALREADY_CLOSED_CODE  = ALREADY_CLOSED_CODE
  integer, public, parameter :: FILE_ALREADY_EXISTED_CODE = ALREADY_EXISTED_CODE

  integer, public, parameter :: FILE_FILE_MAX = FILE_MAX
  integer, public, parameter :: FILE_VAR_MAX  = VAR_MAX
  integer, public, parameter :: FILE_RANK_MAX = RANK_MAX

  !--- missing value
  real(DP), public, parameter :: FILE_RMISS = RMISS

  !--- struct for data infomation
  type, public :: datainfo
     character(len=FILE_HSHORT) :: varname
     character(len=FILE_HMID)   :: description
     character(len=FILE_HSHORT) :: units
     integer                    :: datatype
     integer                    :: rank
     character(len=FILE_HSHORT) :: dim_name(RANK_MAX)
     integer                    :: dim_size(RANK_MAX)
     integer                    :: step
     real(DP)                   :: time_start
     real(DP)                   :: time_end
     character(len=FILE_HMID)   :: time_units
     character(len=FILE_HSHORT) :: calendar
     integer                    :: fid
  endtype datainfo

  integer, public, parameter :: FILE_preclist(0:3) = (/ 4, 8, 4, 8 /)

  character(len=FILE_HSHORT), public :: FILE_dtypelist(0:4)

  data FILE_dtypelist / "REAL4", "REAL8", "INTEGER2", "INTEGER4", "INTEGER8" /

end module scale_file_h
