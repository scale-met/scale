!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM common variables
!!
!! @par Description
!!          Common variables used in the Super Droplet Method (SDM)
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-06-12 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!!
!<
!-------------------------------------------------------------------------------
module sdm_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use rng_uniform_mt
  !-----------------------------------------------------------------------------
  implicit none
  public
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
!  character(len=H_LONG), save :: SD_IN_BASENAME = ''
!  character(len=H_LONG), save :: SD_OUT_BASENAME = ''
!  character(len=H_LONG), save :: RANDOM_IN_BASENAME = ''
!  character(len=H_LONG), save :: RANDOM_OUT_BASENAME = ''
  type(c_rng_uniform_mt), save :: rng_s2c
  integer, save :: fid_sd_i, fid_sd_o
  integer, save :: fid_random_i, fid_random_o
  logical, save :: sd_rest_flg_in = .false. ! restart flg of Super Droplet
  logical, save :: sd_first = .true. 
  !
  !++ Basic variable for SDM
  !
  !------------------------------------------------------------------------------
  integer(DP), allocatable :: sdn_s2c(:)   ! multipilicity

end module sdm_common
