!-------------------------------------------------------------------------------
!> Pre-defined grid index
!!
!! @par Description
!!          Include file of grid index
!!          If you set environment "SCALE_USE_FIXEDINDEX=T", this file is used
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ grid parameters (500m res., 17km isotropic, 22km model top)
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: KMAX   =  128 ! # of computational cells: z
  integer, public, parameter :: IMAX   =   32 ! # of computational cells: x
  integer, public, parameter :: JMAX   =   32 ! # of computational cells: y

  integer, public, parameter :: IBLOCK =    8 ! block size for cache blocking: x
  integer, public, parameter :: JBLOCK =    8 ! block size for cache blocking: y
