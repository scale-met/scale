!-------------------------------------------------------------------------------
!> module FILE I/O HEADER
!!
!! @par Description
!!          File I/O module (Parameter Container)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-06-13 (S.Nishizawa) [new] Imported from SCALE-LES
!!
!<
!-------------------------------------------------------------------------------
module gtool_file_h
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  use dc_types, only: &
     DP
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- character length 
  integer, public, parameter :: File_HSHORT =  16
  integer, public, parameter :: File_HMID   =  64
  integer, public, parameter :: File_HLONG  = 256

  !--- data type 
  integer, public, parameter :: File_REAL4    = 0
  integer, public, parameter :: File_REAL8    = 1
  integer, public, parameter :: File_INTEGER4 = 2
  integer, public, parameter :: File_INTEGER8 = 3

  !--- action type 
  integer, public, parameter :: File_FREAD   = 0
  integer, public, parameter :: File_FWRITE  = 1
  integer, public, parameter :: File_FAPPEND = 2

  !--- return codes
  integer, public, parameter :: ERROR_CODE          = -1
  integer, public, parameter :: SUCCESS_CODE        =  0
  integer, public, parameter :: ALREADY_CLOSED_CODE =  1
  integer, public, parameter :: ALREADY_EXISTED_CODE = 2

  integer, public, parameter :: MAX_RANK = 10

  !--- struct for data infomation
  type, public :: datainfo
     character(LEN=File_HSHORT) :: varname
     character(LEN=File_HMID)   :: description
     character(LEN=File_HSHORT) :: units
     integer                    :: datatype
     integer                    :: rank
     character(LEN=File_HSHORT) :: dim_name(MAX_RANK)
     integer                    :: dim_size(MAX_RANK)
     integer                    :: step
     real(DP)                    :: time_start
     real(DP)                    :: time_end
     integer                    :: fid
  endtype datainfo

  integer, public, parameter :: File_preclist(0:3) = (/ 4, 8, 4, 8 /)

end module gtool_file_h
!-------------------------------------------------------------------------------
