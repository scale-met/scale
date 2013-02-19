!-------------------------------------------------------------------------------
!> module FILE I/O HEADER
!!
!! @par Description
!!          File I/O module (Parameter Container)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-FIO
!!
!<
!-------------------------------------------------------------------------------
module mod_fileio_c_h
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
  integer, public, parameter :: FIO_HSHORT =  16
  integer, public, parameter :: FIO_HMID   =  64
  integer, public, parameter :: FIO_HLONG  = 256

  !--- data type 
  integer, public, parameter :: FIO_REAL4    = 0
  integer, public, parameter :: FIO_REAL8    = 1
  integer, public, parameter :: FIO_INTEGER4 = 2
  integer, public, parameter :: FIO_INTEGER8 = 3

  !--- data endian 
  integer, public, parameter :: FIO_UNKNOWN_ENDIAN = 0
  integer, public, parameter :: FIO_LITTLE_ENDIAN  = 1
  integer, public, parameter :: FIO_BIG_ENDIAN     = 2

  !--- topology 
  integer, public, parameter :: FIO_ICOSAHEDRON = 0
  integer, public, parameter :: FIO_IGA_LCP     = 1
  integer, public, parameter :: FIO_IGA_MLCP    = 2
  integer, public, parameter :: FIO_CARTESIAN   = 3
  integer, public, parameter :: FIO_NONE        = 4

  !--- file mode (partial or complete) 
  integer, public, parameter :: FIO_SPLIT_FILE = 0
  integer, public, parameter :: FIO_INTEG_FILE = 1

  !--- proccessor type 
  integer, public, parameter :: FIO_SINGLE_PROC = 0
  integer, public, parameter :: FIO_MULTI_PROC  = 1

  !--- MPI-IO 
  integer, public, parameter :: FIO_MPIIO_NOUSE = 0
  integer, public, parameter :: FIO_MPIIO_USE   = 1

  !--- action type 
  integer, public, parameter :: FIO_FREAD   = 0
  integer, public, parameter :: FIO_FWRITE  = 1
  integer, public, parameter :: FIO_FAPPEND = 2

  !--- data dump type 
  integer, public, parameter :: FIO_DUMP_OFF    = 0
  integer, public, parameter :: FIO_DUMP_HEADER = 1
  integer, public, parameter :: FIO_DUMP_ALL    = 2

  !--- struct for package infomation
  type, public :: headerinfo
     character(LEN=FIO_HLONG) :: fname
     character(LEN=FIO_HMID)  :: description
     character(LEN=FIO_HLONG) :: note
     integer                  :: num_of_data
     integer                  :: fmode
     integer                  :: endiantype
     integer                  :: grid_topology
     integer                  :: glevel
     integer                  :: rlevel
     integer                  :: num_of_rgn
     integer, pointer         :: rgnid(:)
  endtype headerinfo

  !--- struct for data infomation
  type, public :: datainfo
     character(LEN=FIO_HSHORT) :: varname
     character(LEN=FIO_HMID)   :: description
     character(LEN=FIO_HSHORT) :: unit
     character(LEN=FIO_HSHORT) :: layername
     character(LEN=FIO_HLONG)  :: note
     integer(8)                :: datasize
     integer                   :: datatype
     integer                   :: num_of_layer
     integer                   :: step
     integer(8)                :: time_start
     integer(8)                :: time_end
  endtype datainfo

  integer, public, parameter :: FIO_preclist(0:3) = (/ 4, 8, 4, 8 /)

end module mod_fileio_c_h
!-------------------------------------------------------------------------------
