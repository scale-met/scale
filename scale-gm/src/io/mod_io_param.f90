!-------------------------------------------------------------------------------
!> Module file I/O parameter
!!
!! @par Description
!!         This module is parameter list for file I/O (common)
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_io_param
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !--- character length
  integer, public, parameter :: IO_HSHORT =   32 !< character length for short var.
  integer, public, parameter :: IO_HMID   =  128 !< character length for middle var.
  integer, public, parameter :: IO_HLONG  = 1024 !< character length for long var.

  !--- data type
  integer, public, parameter :: IO_REAL4    = 0 !< ID for 4byte real
  integer, public, parameter :: IO_REAL8    = 1 !< ID for 8byte real
  integer, public, parameter :: IO_INTEGER4 = 2 !< ID for 4byte int
  integer, public, parameter :: IO_INTEGER8 = 3 !< ID for 8byte int

  !--- data endian
  integer, public, parameter :: IO_UNKNOWN_ENDIAN = 0 !< ID for unknown endian
  integer, public, parameter :: IO_LITTLE_ENDIAN  = 1 !< ID for little endian
  integer, public, parameter :: IO_BIG_ENDIAN     = 2 !< ID for big endian

  !--- topology
  integer, public, parameter :: IO_ICOSAHEDRON = 0 !< ID for ico grid
  integer, public, parameter :: IO_IGA_LCP     = 1 !< ID for LCP grid
  integer, public, parameter :: IO_IGA_MLCP    = 2 !< ID for MLCP grid

  !--- file mode (partial or complete)
  integer, public, parameter :: IO_SPLIT_FILE = 0 !< ID for split(partical) file
  integer, public, parameter :: IO_INTEG_FILE = 1 !< ID for integrated(complete) file

  !--- proccessor type
  integer, public, parameter :: IO_SINGLE_PROC = 0 !< ID for single processor
  integer, public, parameter :: IO_MULTI_PROC  = 1 !< ID for multi processor

  !--- action type
  integer, public, parameter :: IO_FREAD   = 0 !< ID for read file
  integer, public, parameter :: IO_FWRITE  = 1 !< ID for write file
  integer, public, parameter :: IO_FAPPEND = 2 !< ID for append file

  !--- data dump type
  integer, public, parameter :: IO_DUMP_OFF      = 0 !< Dumping off
  integer, public, parameter :: IO_DUMP_HEADER   = 1 !< Dump header only
  integer, public, parameter :: IO_DUMP_ALL      = 2 !< Dump all
  integer, public, parameter :: IO_DUMP_ALL_MORE = 3 !< Dump all and more

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
end module mod_io_param
