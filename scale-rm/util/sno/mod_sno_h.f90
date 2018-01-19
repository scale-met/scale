!-------------------------------------------------------------------------------
!> Module SNO (RM)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE NetCDF Operator (SNO)
!!          header module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_sno_h
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_MID), public, parameter :: TOOLNAME = "SCALE NetCDF Operator"

  integer, public, parameter :: item_limit = 1000  !< limit of item
  integer, public, parameter :: step_limit = 10000 !< limit of steps          for each item
  integer, public, parameter :: dim_limit  = 3     !< limit of dimension rank for each item

  integer, public, parameter :: I_map_p = 1
  integer, public, parameter :: I_map_i = 2
  integer, public, parameter :: I_map_j = 3

  ! struct for common infomation
  type, public :: commoninfo
    character(len=H_MID)   :: title               ! title of the file
    character(len=H_MID)   :: source              ! for file header
    character(len=H_MID)   :: institute           ! for file header
    character(len=5)       :: periodic(3)         ! periodic condition?            (z:y)
    integer                :: gridsize(3)         ! total grid size     in global  (z:y), always including halo
    integer                :: halosize(3)         ! halo  grid size     in global  (z:y), always existing
    character(len=H_MID)   :: time_units
    real(DP)               :: time_start(1)
    integer                :: xatt_size_global(1)
    integer                :: xatt_halo_global(2)
    integer                :: yatt_size_global(1)
    integer                :: yatt_halo_global(2)

    character(len=H_SHORT) :: minfo_mapping_name
    real(DP)               :: minfo_false_easting                        (1)
    real(DP)               :: minfo_false_northing                       (1)
    real(DP)               :: minfo_longitude_of_central_meridian        (1)
    real(DP)               :: minfo_longitude_of_projection_origin       (1)
    real(DP)               :: minfo_latitude_of_projection_origin        (1)
    real(DP)               :: minfo_straight_vertical_longitude_from_pole(1)
    real(DP)               :: minfo_standard_parallel                    (2)
  end type commoninfo

  ! struct for axis infomation
  type, public :: axisinfo
    character(len=H_SHORT) :: varname
    character(len=H_MID)   :: description
    character(len=H_SHORT) :: units
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name(dim_limit)
    integer                :: dim_size(dim_limit)
    logical                :: transpose
    logical                :: regrid
    real(RP), allocatable  :: AXIS_1d(:)
    real(RP), allocatable  :: AXIS_2d(:,:)
    real(RP), allocatable  :: AXIS_3d(:,:,:)
  end type axisinfo

  ! struct for item infomation
  type, public :: iteminfo
    character(len=H_SHORT) :: varname
    character(len=H_MID)   :: description
    character(len=H_SHORT) :: units
    integer                :: datatype
    integer                :: dim_rank
    character(len=H_SHORT) :: dim_name  (dim_limit)
    integer                :: dim_size  (dim_limit)
    logical                :: transpose
    integer                :: step_nmax
    real(DP)               :: time_start(step_limit)
    real(DP)               :: time_end  (step_limit)
    real(DP)               :: dt
    character(len=H_MID)   :: time_units
    real(RP), allocatable  :: VAR_1d(:)
    real(RP), allocatable  :: VAR_2d(:,:)
    real(RP), allocatable  :: VAR_3d(:,:,:)
  end type iteminfo

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
end module mod_sno_h
