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
module scale_atmos_numeric_fdm_cd
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_index
  use scale_tracer

  use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY,  &
       I_UY,  &
       I_XV,  &
       I_UV
  
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

  public :: FDM_CD_setup

  !**
  public :: FDM_EvalFlxCD4_VarZXY
  public :: FDM_EvalFlxCD4_VarZUY
  public :: FDM_EvalFlxCD4_VarZXV
  public :: FDM_EvalFlxCD4_VarWXY
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------

  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

  integer :: IFS_OFF
  integer :: JFS_OFF
  real(RP), parameter :: FACT_N =   7.0_RP/12.0_RP
  real(RP), parameter :: FACT_F = - 1.0_RP/12.0_RP
  
contains
  subroutine FDM_CD_setup(IFS_OFF_, JFS_OFF_)
    integer, intent(in) :: IFS_OFF_, JFS_OFF_

    IFS_OFF = IFS_OFF_; JFS_OFF = JFS_OFF_    
  end subroutine FDM_CD_setup

  subroutine FDM_EvalFlxCD4_VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,          &  ! (inout)
       & Var_ZXY, VelX_ZUY, VelY_ZXV, VelZ_WXY,                            &  ! (in)
       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                          &  ! (in)
       & IIS, IIE, JJS, JJE, KS, KE )                                         ! (in)

    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZXY
    real(RP), intent(in), dimension(KA,IA,JA)     :: VelX_ZUY, VelY_ZXV, VelZ_WXY
    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
    real(RP), intent(in) :: J33G
    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
    
  end subroutine FDM_EvalFlxCD4_VarZXY

  subroutine FDM_EvalFlxCD4_VarZUY( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,          &  ! (inout)
       & Var_ZUY, VelX_ZUY, VelY_ZXV, VelZ_WXY,                         &  ! (in)
       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
       & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)

    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZUY
    real(RP), intent(in), dimension(KA,IA,JA)     :: VelX_ZUY, VelY_ZXV, VelZ_WXY
    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
    real(RP), intent(in) :: J33G
    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
    
  end subroutine FDM_EvalFlxCD4_VarZUY
  
  subroutine FDM_EvalFlxCD4_VarZXV( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,          &  ! (inout)
       & Var_ZXV, VelX_ZUY, VelY_ZXV, VelZ_WXY,                         &  ! (in)
       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
       & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)

    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZXV
    real(RP), intent(in), dimension(KA,IA,JA)     :: VelX_ZUY, VelY_ZXV, VelZ_WXY
    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
    real(RP), intent(in) :: J33G
    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
    
  end subroutine FDM_EvalFlxCD4_VarZXV

  subroutine FDM_EvalFlxCD4_VarWXY( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,          &  ! (inout)
       & Var_WXY, VelX_ZUY, VelY_ZXV, VelZ_WXY,                         &  ! (in)
       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
       & IIS, IIE, JJS, JJE, KS, KE )                         ! (in)

    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_WXY
    real(RP), intent(in), dimension(KA,IA,JA)     :: VelX_ZUY, VelY_ZXV, VelZ_WXY
    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
    real(RP), intent(in) :: J33G
    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
    
  end subroutine FDM_EvalFlxCD4_VarWXY

end module scale_atmos_numeric_fdm_cd
