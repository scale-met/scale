#! ruby -Ks

require 'erb'

glevel  = ARGV[0].to_i
rlevel  = ARGV[1].to_i
nmpi    = ARGV[2].to_i
zlayer  = ARGV[3].to_i
diamond = ARGV[4].to_i

regionall = 2**rlevel * 2**rlevel * diamond
grid1D    = 2**(glevel-rlevel)

contents = <<EOS
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
  ! Total inner grid (2D) = <%= regionall * grid1D * grid1D + 2 %>
  ! Total inner grid (3D) = <%= ( regionall * grid1D * grid1D + 2 ) * zlayer %>
  ! Total grid       (2D) = <%= (1+grid1D+1) * (1+grid1D+1) %>
  ! Total grid       (3D) = <%= ( regionall * (1+grid1D+1) * (1+grid1D+1) ) * (1+zlayer+1) %>
  !-----------------------------------------------------------------------------

  ! main parameter
  integer, public, parameter :: ADM_glevel      = <%= glevel  %>
  integer, public, parameter :: ADM_rlevel      = <%= rlevel  %>
  integer, public, parameter :: ADM_vlayer      = <%= zlayer  %>
  integer, public, parameter :: ADM_DMD         = <%= diamond %>
  integer, public            :: ADM_prc_all     = <%= nmpi    %>

  ! region
  integer, public, parameter :: ADM_rgn_nmax    = <%= regionall %>
  integer, public, parameter :: ADM_lall        = <%= regionall / nmpi  %>
  integer, public, parameter :: ADM_rgn_nmax_pl = 2 ! number of pole    region
  integer, public, parameter :: ADM_lall_pl     = 2 ! number of pole    region per process

  ! horizontal grid
  integer, public, parameter :: ADM_gall        = <%= (1+grid1D+1) * (1+grid1D+1) %>
  integer, public, parameter :: ADM_gall_in     = <%= (  grid1D+1) * (  grid1D+1) %>
  integer, public, parameter :: ADM_gall_1d     = <%= 1 + grid1D + 1 %>
  integer, public, parameter :: ADM_gmin        = <%= 1 + 1          %>
  integer, public, parameter :: ADM_gmax        = <%= 1 + grid1D     %>

  integer, public, parameter :: ADM_gall_pl     = 6
  integer, public, parameter :: ADM_gslf_pl     = 1
  integer, public, parameter :: ADM_gmin_pl     = 2
  integer, public, parameter :: ADM_gmax_pl     = 6

  ! vertical grid
  integer, public, parameter :: ADM_kall        = <%= 1 + zlayer + 1 %>
  integer, public, parameter :: ADM_kmin        = <%= 1 + 1          %>
  integer, public, parameter :: ADM_kmax        = <%= 1 + zlayer     %>

  ! List vectors
  integer, public, parameter :: ADM_IooJoo_nmax = <%= ( grid1D   ) * ( grid1D   ) %>
  integer, public, parameter :: ADM_IooJmo_nmax = <%= ( grid1D   ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_IooJop_nmax = <%= ( grid1D   ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_IooJmp_nmax = <%= ( grid1D   ) * ( grid1D+2 ) %>
  integer, public, parameter :: ADM_ImoJoo_nmax = <%= ( grid1D+1 ) * ( grid1D   ) %>
  integer, public, parameter :: ADM_ImoJmo_nmax = <%= ( grid1D+1 ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_ImoJop_nmax = <%= ( grid1D+1 ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_ImoJmp_nmax = <%= ( grid1D+1 ) * ( grid1D+2 ) %>
  integer, public, parameter :: ADM_IopJoo_nmax = <%= ( grid1D+1 ) * ( grid1D   ) %>
  integer, public, parameter :: ADM_IopJmo_nmax = <%= ( grid1D+1 ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_IopJop_nmax = <%= ( grid1D+1 ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_IopJmp_nmax = <%= ( grid1D+1 ) * ( grid1D+2 ) %>
  integer, public, parameter :: ADM_ImpJoo_nmax = <%= ( grid1D+2 ) * ( grid1D   ) %>
  integer, public, parameter :: ADM_ImpJmo_nmax = <%= ( grid1D+2 ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_ImpJop_nmax = <%= ( grid1D+2 ) * ( grid1D+1 ) %>
  integer, public, parameter :: ADM_ImpJmp_nmax = <%= ( grid1D+2 ) * ( grid1D+2 ) %>

!-------------------------------------------------------------------------------
EOS

erb = ERB.new(contents).result
open("./inc_index.h", "w") {|f| f.write erb}
