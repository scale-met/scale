!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!          common subroutines for Atmospheric dynamical core used in finite diffrence scheme.
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-12-14 (Y.Kawai) [new] Add this module
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_numeric_fdm_def
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !**
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------

  !* Centeral difference 

  integer, parameter, public :: FLXEVALTYPE_CD2 = 101
  integer, parameter, public :: FLXEVALTYPE_CD4 = 102
  integer, parameter, public :: FLXEVALTYPE_CD6 = 103

  character(len=H_SHORT), parameter, public :: FLXEVALTYPENAME_CD2 = 'FDM_CD2'
  character(len=H_SHORT), parameter, public :: FLXEVALTYPENAME_CD4 = 'FDM_CD4'
  character(len=H_SHORT), parameter, public :: FLXEVALTYPENAME_CD6 = 'FDM_CD6'

  !* Upwind difference 
  
  integer, parameter, public :: FLXEVALTYPE_UD1 = 201  
  integer, parameter, public :: FLXEVALTYPE_UD3 = 202
  integer, parameter, public :: FLXEVALTYPE_UD5 = 203

  character(len=H_SHORT), parameter, public :: FLXEVALTYPENAME_UD1 = 'FDM_UD1'
  character(len=H_SHORT), parameter, public :: FLXEVALTYPENAME_UD3 = 'FDM_UD3'
  character(len=H_SHORT), parameter, public :: FLXEVALTYPENAME_UD5 = 'FDM_UD5'
  
  !* ID for location definted variable
  
  integer, parameter, public :: VL_XY        = 011
  integer, parameter, public :: VL_UY        =  21
  integer, parameter, public :: VL_XV        =  12  
  integer, parameter, public :: VL_ZXY       = 111
  integer, parameter, public :: VL_ZUY       = 121
  integer, parameter, public :: VL_ZXV       = 112
  integer, parameter, public :: VL_WXY       = 211
  
end module scale_atmos_numeric_fdm_def
