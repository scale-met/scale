!-------------------------------------------------------------------------------
!> Pre-defined grid index
!!
!! @par Description
!!          Include file of grid index
!!          If you set environment "ENABLE_FIXEDINDEX=T", this file is used
!!          If you set this file, execute following:
!!          "make fixedindex glevel=xxx rlevel=xxx nmpi=xxx zlayer=xxx diamond=xxx".
!<
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Total inner grid (2D) = 163842
  ! Total inner grid (3D) = 15401148
  ! Total grid       (2D) = 16900
  ! Total grid       (3D) = 16224000
  !-----------------------------------------------------------------------------

  ! main parameter
  integer, public, parameter :: ADM_glevel      = 7
  integer, public, parameter :: ADM_rlevel      = 0
  integer, public, parameter :: ADM_vlayer      = 94
  integer, public, parameter :: ADM_DMD         = 10
  integer, public            :: ADM_prc_all     = 10

  ! region
  integer, public, parameter :: ADM_rgn_nmax    = 10
  integer, public, parameter :: ADM_lall        = 1
  integer, public, parameter :: ADM_rgn_nmax_pl = 2 ! number of pole    region
  integer, public, parameter :: ADM_lall_pl     = 2 ! number of pole    region per process

  ! horizontal grid
  integer, public, parameter :: ADM_gall        = 16900
  integer, public, parameter :: ADM_gall_in     = 16641
  integer, public, parameter :: ADM_gall_1d     = 130
  integer, public, parameter :: ADM_gmin        = 2
  integer, public, parameter :: ADM_gmax        = 129

  integer, public, parameter :: ADM_gall_pl     = 6
  integer, public, parameter :: ADM_gslf_pl     = 1
  integer, public, parameter :: ADM_gmin_pl     = 2
  integer, public, parameter :: ADM_gmax_pl     = 6

  ! vertical grid
  integer, public, parameter :: ADM_kall        = 96
  integer, public, parameter :: ADM_kmin        = 2
  integer, public, parameter :: ADM_kmax        = 95

  ! List vectors
  integer, public, parameter :: ADM_IooJoo_nmax = 16384
  integer, public, parameter :: ADM_IooJmo_nmax = 16512
  integer, public, parameter :: ADM_IooJop_nmax = 16512
  integer, public, parameter :: ADM_IooJmp_nmax = 16640
  integer, public, parameter :: ADM_ImoJoo_nmax = 16512
  integer, public, parameter :: ADM_ImoJmo_nmax = 16641
  integer, public, parameter :: ADM_ImoJop_nmax = 16641
  integer, public, parameter :: ADM_ImoJmp_nmax = 16770
  integer, public, parameter :: ADM_IopJoo_nmax = 16512
  integer, public, parameter :: ADM_IopJmo_nmax = 16641
  integer, public, parameter :: ADM_IopJop_nmax = 16641
  integer, public, parameter :: ADM_IopJmp_nmax = 16770
  integer, public, parameter :: ADM_ImpJoo_nmax = 16640
  integer, public, parameter :: ADM_ImpJmo_nmax = 16770
  integer, public, parameter :: ADM_ImpJop_nmax = 16770
  integer, public, parameter :: ADM_ImpJmp_nmax = 16900

!-------------------------------------------------------------------------------
