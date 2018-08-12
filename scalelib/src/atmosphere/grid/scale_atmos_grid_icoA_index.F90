!-------------------------------------------------------------------------------
!> module atmosphere / grid / icosahedralA / index
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
  integer, public             :: IMAX = -1
  integer, public, parameter  :: IS = 1
  integer, public             :: IE
  integer, public             :: IA
  integer, public             :: JMAX = -1
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
    use scale_prc, only: &
       PRC_abort
    use scale_prc_icoA, only: &
       PRC_RGN_level, &
       PRC_RGN_local, &
       PRC_RGN_vlink
    implicit none

    integer :: GRID_LEVEL = -1 !> grid division level
    integer :: KMAX       = -1 !< # of computational cells: z, local

    namelist / PARAM_ATMOS_GRID_ICOA_INDEX / &
       GRID_LEVEL, &
       KMAX

    integer :: nmax

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_GRID_ICOA_INDEX_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_GRID_ICOA_INDEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_GRID_ICOA_INDEX_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_GRID_ICOA_INDEX_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_GRID_ICOA_INDEX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_GRID_ICOA_INDEX)

    ! Error if glevel & rlevel are not defined
    if ( GRID_LEVEL < 1 ) then
       LOG_ERROR("ATMOS_GRID_ICOA_INDEX_setup",*) 'GRID_LEVEL must be >= 1! ', GRID_LEVEL
       call PRC_abort
    endif

    if ( KMAX < 2 ) then
       LOG_ERROR("ATMOS_GRID_ICOA_INDEX_setup",*) 'KMAX must be >= 2! ', KMAX
       call PRC_abort
    end if

    ADM_glevel  = GRID_LEVEL
    ADM_vlayer  = KMAX

    ADM_lall    = PRC_RGN_local

    nmax        = 2**(ADM_glevel-PRC_RGN_level)
    ADM_gall_1d = 1 + nmax + 1
    ADM_gmin    = 1 + 1
    ADM_gmax    = 1 + nmax

    ADM_gall    = ( 1+nmax+1 ) * ( 1+nmax+1 )
    ADM_gall_in = (   nmax+1 ) * (   nmax+1 )

    ADM_imin    = 1 + 1
    ADM_imax    = 1 + nmax
    ADM_iall    = 1 + nmax + 1
    ADM_jmin    = 1 + 1
    ADM_jmax    = 1 + nmax
    ADM_jall    = 1 + nmax + 1

    ADM_vlink   = PRC_RGN_vlink
    ADM_gall_pl = PRC_RGN_vlink + 1
    ADM_gmax_pl = PRC_RGN_vlink + 1

    if ( ADM_vlayer == 1 ) then
       ADM_kall = 1
       ADM_kmin = 1
       ADM_kmax = 1
    else
       ADM_kall = 1 + ADM_vlayer + 1
       ADM_kmin = 1 + 1
       ADM_kmax = 1 + ADM_vlayer
    endif

    ! for physics grid
    KS = ADM_kmin
    KE = ADM_kmax
    KA = ADM_kall

    IMAX = nmax
    ! IS = 1
    IE = nmax + 1
    IA = IE

    JMAX = nmax
    ! JS = 1
    JE = nmax + 1
    JA = JE

    return
  end subroutine ATMOS_GRID_ICOA_INDEX_setup


end module scale_atmos_grid_icoA_index
