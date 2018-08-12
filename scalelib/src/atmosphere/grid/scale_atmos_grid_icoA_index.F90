!-------------------------------------------------------------------------------
!> module atmosphere / grid / icoA index
!!
!! @par Description
!!          Atmospheric grid Index module for the icosaphdralA grid
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_grid_icoA_index
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_GRID_ICOA_INDEX_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  ! Identifier of triangle element (i-axis-side or j-axis side)
  integer, public, parameter :: ADM_TI = 1
  integer, public, parameter :: ADM_TJ = 2

  ! Identifier of arc element (i-axis-side, ij-axis side, or j-axis side)
  integer, public, parameter :: ADM_AI  = 1
  integer, public, parameter :: ADM_AIJ = 2
  integer, public, parameter :: ADM_AJ  = 3

  ! Identifier of 1 variable
  integer, public, parameter :: ADM_KNONE = 1

  ! dimension of the spacial vector
  integer, public, parameter :: ADM_nxyz = 3

  ! number of HALOs
  integer, public, parameter :: KHALO  = 1  !< # of halo cells: z
  integer, public            :: IHALO  = 1  !< # of halo cells: x
  integer, public            :: JHALO  = 1  !< # of halo cells: y

  ! main parameter
  integer, public            :: ADM_glevel           ! grid   division level
  integer, public            :: ADM_vlayer           ! number of vertical layer

  ! region
  integer, public            :: ADM_lall             ! number of regular region per process
  integer, public, parameter :: ADM_lall_pl     =  2 ! number of pole    region per process

  ! horizontal grid
  integer, public            :: ADM_gall             ! number of horizontal grid per regular region
  integer, public            :: ADM_gall_in          ! number of horizontal grid (inner part)
  integer, public            :: ADM_gall_1d          ! number of horizontal grid (1D)
  integer, public            :: ADM_gmin             ! start index of 1D horizontal grid
  integer, public            :: ADM_gmax             ! end   index of 1D horizontal grid

  integer, public            :: ADM_iall             ! number of horizontal grid per regular region (i-axis)
  integer, public            :: ADM_imin             ! start index of 1D horizontal grid
  integer, public            :: ADM_imax             ! end   index of 1D horizontal grid
  integer, public            :: ADM_jall             ! number of horizontal grid per regular region (j-axis)
  integer, public            :: ADM_jmin             ! start index of 1D horizontal grid
  integer, public            :: ADM_jmax             ! end   index of 1D horizontal grid

  integer, public            :: ADM_vlink       = -1 ! maximum number of vertex linkage, ICO:5, PSP:6, LCP, MLCP:k
  integer, public            :: ADM_gall_pl          ! number of horizontal grid for pole region
  integer, public, parameter :: ADM_gslf_pl     =  1 ! index for pole point
  integer, public, parameter :: ADM_gmin_pl     =  2 ! start index of grid around the pole point
  integer, public            :: ADM_gmax_pl          ! end   index of grid around the pole point

  ! vertical grid
  integer, public            :: ADM_kall             ! number of vertical grid
  integer, public            :: ADM_kmin             ! start index of vertical grid
  integer, public            :: ADM_kmax             ! end   index of vertical grid

  ! for physics grid
  integer, public             :: KS
  integer, public             :: KE
  integer, public             :: KA
  integer, public, parameter  :: IS = 1
  integer, public             :: IE
  integer, public             :: IA
  integer, public, parameter  :: JS = 1
  integer, public             :: JE
  integer, public             :: JA

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup index
  subroutine ATMOS_GRID_ICOA_INDEX_setup

    return
  end subroutine ATMOS_GRID_ICOA_INDEX_setup


end module scale_atmos_grid_icoA_index
