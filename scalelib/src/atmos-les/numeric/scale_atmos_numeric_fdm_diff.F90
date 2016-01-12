!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!          common subroutines for Atmospheric dynamical core used in finite diffrence scheme.
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-01-06 (Y.Kawai) [new] Add this module
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_numeric_fdm_diff
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

#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif

  use scale_atmos_numeric_fdm_def, only: &
       & VL_ZXY, VL_ZUY, VL_ZXV, VL_WXY
       
  use scale_atmos_numeric_fdm_util, only: &
       & FDM_RhoVar2Var
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: FDM_EvalDiffFlxCD2_VarZXY
  public :: FDM_EvalDiffFlxCD2_VarZUY
  public :: FDM_EvalDiffFlxCD2_VarZXV
  public :: FDM_EvalDiffFlxCD2_VarWXY

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine FDM_EvalDiffFlxCD2_VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,         &  ! (inout)
    & RhoVar_ZXY, Dens_ZXY, DiffCoef, RFDX, RFDY, RFDZ,                       &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)    :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: RhoVar_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY
     real(RP), intent(in)                          :: DiffCoef
     real(RP), intent(in)                          :: RFDX(IA), RFDY(JA), RFDZ(KA)
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer,  intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: k, i, j
     real(RP) :: Var_ZXY(KA,IA,JA)

     do j = JJS-JHALO, JJE+JHALO
     do i = IIS-IHALO, IIE+IHALO
       do k=KS, KE
         Var_ZXY(k,i,j) = RhoVar_ZXY(k,i,j) / DENS_ZXY(k,i,j)
       enddo
     end do
     end do
     
     do j = JJS, JJE
     do i = IIS, IIE
       do k=KS, KE-1
         FlxZ_WXY(k,i,j) =   DiffCoef * RFDZ(k) * 0.5_RP * sum(Dens_ZXY(k:k+1,i,j)) &
              &            * (Var_ZXY(k,i,j) - Var_ZXY(k+1,i,j))
       enddo
       FlxZ_WXY(KS-1,i,j) = 0.0_RP
       FlxZ_WXY(  KE,i,j) = 0.0_RP
     enddo
     enddo

     do j = JJS, JJE
     do i = IIS-1, min(IIE,IEH)
       do k=KS, KE
         FlxX_ZUY(k,i,j) =   DiffCoef * RFDX(i) * 0.5_RP * sum(Dens_ZXY(k,i:i+1,j)) &
              &            * (Var_ZXY(k,i,j) - Var_ZXY(k,i+1,j))
       enddo
     enddo
     enddo

     do j = JJS-1, min(JJE,JEH)
     do i = IIS, IIE
       do k=KS, KE
         FlxY_ZXV(k,i,j) =   DiffCoef * RFDY(j) * 0.5_RP * sum(Dens_ZXY(k,i,j:j+1)) &
              &            * (Var_ZXY(k,i,j) - Var_ZXY(k,i,j+1))
       enddo
     enddo
     enddo
     
  end subroutine FDM_EvalDiffFlxCD2_VarZXY

  subroutine FDM_EvalDiffFlxCD2_VarZUY( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,         &  ! (inout)
    & RhoVar_ZUY, Dens_ZXY, DiffCoef, RCDX, RFDY, RFDZ,                       &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)    :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: RhoVar_ZUY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY
     real(RP), intent(in)                          :: DiffCoef
     real(RP), intent(in)                          :: RCDX(IA), RFDY(JA), RFDZ(KA)
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer,  intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: k, i, j
     real(RP) :: Var_ZUY(KA,IA,JA)

     do j = JJS-JHALO+1, JJE+JHALO
     do i = IIS-IHALO, IIE+IHALO-1
       do k=KS, KE
         Var_ZUY(k,i,j) = 2.0_RP * RhoVar_ZUY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i+1,j))
       enddo
     end do
     end do

     
     do j = JJS, JJE
     do i = IIS, IIE
       do k=KS, KE-1
         FlxZ_WUY(k,i,j) =   DiffCoef * RFDZ(k) * 0.25_RP * sum(Dens_ZXY(k:k+1,i:i+1,j)) &
           &               * (Var_ZUY(k,i,j) - Var_ZUY(k+1,i,j))
       enddo
       FlxZ_WUY(KS-1,i,j) = 0.0_RP
       FlxZ_WUY(  KE,i,j) = 0.0_RP
     enddo
     enddo

     do j = JJS, JJE
     do i = IIS, IIE+1
       do k=KS, KE
         FlxX_ZXY(k,i-1,j) =   DiffCoef * RCDX(i) * Dens_ZXY(k,i,j) &
           &                 * (Var_ZUY(k,i-1,j) - Var_ZUY(k,i,j))
       enddo
     enddo
     enddo

     do j = JJS-1, JJE
     do i = IIS, IIE
       do k=KS, KE
         FlxY_ZUV(k,i,j) =   DiffCoef * RFDY(j) * 0.25_RP * sum(Dens_ZXY(k,i:i+1,j:j+1)) &
           &               * (Var_ZUY(k,i,j) - Var_ZUY(k,i,j+1))
       enddo
     enddo
     enddo
     
  end subroutine FDM_EvalDiffFlxCD2_VarZUY

  subroutine FDM_EvalDiffFlxCD2_VarZXV( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,        &  ! (inout)
    & RhoVar_ZXV, Dens_ZXY, DiffCoef, RFDX, RCDY, RFDZ,                      &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                       &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE )                                              ! (in)

     real(RP), intent(out), dimension(KA,IA,JA)    :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
     real(RP), intent(in),  dimension(KA,IA,JA)    :: RhoVar_ZXV
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY
     real(RP), intent(in)                          :: DiffCoef
     real(RP), intent(in)                          :: RFDX(IA), RCDY(JA), RFDZ(KA)
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer,  intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: k, i, j
     real(RP) :: Var_ZXV(KA,IA,JA)

     do j = JJS-JHALO, JJE+JHALO-1
     do i = IIS-IHALO+1, IIE+IHALO
       do k=KS, KE-1
         Var_ZXV(k,i,j) = 2.0_RP * RhoVar_ZXV(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i,j+1))
       enddo
     end do
     end do

     
     do j = JJS, JJE
     do i = IIS, IIE
       do k=KS, KE-1
         FlxZ_WXV(k,i,j) =   DiffCoef * RFDZ(k) * 0.25_RP * sum(Dens_ZXY(k:k+1,i:i+1,j)) &
           &               * (Var_ZXV(k,i,j) - Var_ZXV(k+1,i,j))
       enddo
       FlxZ_WXV(KS-1,i,j) = 0.0_RP
       FlxZ_WXV(KE,i,j) = 0.0_RP
     enddo
     enddo

     
     do j = JJS, JJE
     do i = IIS-1, IIE
       do k=KS, KE
         FlxX_ZUV(k,i,j) =   DiffCoef * RFDX(i) * 0.25_RP * sum(Dens_ZXY(k,i:i+1,j:j+1)) &
           &                 * (Var_ZXV(k,i,j) - Var_ZXV(k,i+1,j))
       enddo
     enddo
     enddo

     do j = JJS, JJE+1
     do i = IIS, IIE
       do k=KS, KE
         FlxY_ZXY(k,i,j-1) =   DiffCoef * RCDY(j) * Dens_ZXY(k,i,j) &
           &               * (Var_ZXV(k,i,j-1) - Var_ZXV(k,i,j))
       enddo
     enddo
     enddo

  end subroutine FDM_EvalDiffFlxCD2_VarZXV

  subroutine FDM_EvalDiffFlxCD2_VarWXY( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,       &  ! (inout)
    & RhoVar_WXY, Dens_ZXY, DiffCoef, RFDX, RFDY, RCDZ,                     &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                      &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE )                                             ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)    :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: RhoVar_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY
     real(RP), intent(in)                          :: DiffCoef
     real(RP), intent(in)                          :: RFDX(IA), RFDY(JA), RCDZ(KA)
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer,  intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: k, i, j
     real(RP) :: Var_WXY(KA,IA,JA)

     do j = JJS-JHALO, JJE+JHALO
     do i = IIS-IHALO, IIE+IHALO
       do k=KS, KE-1
         Var_WXY(k,i,j) = 2.0_RP * RhoVar_WXY(k,i,j) / (DENS_ZXY(k,i,j) + DENS_ZXY(k+1,i,j))
       enddo
       Var_WXY(KS-1,i,j) = 0.0_RP!- Var_WXY(KS+1,i,j)
       Var_WXY(  KE,i,j) = 0.0_RP!- Var_WXY(KE-2,i,j)
     end do
     end do

     do j = JJS, JJE
     do i = IIS, IIE
       do k=KS, KE
         FlxZ_ZXY(k-1,i,j) =   DiffCoef * RCDZ(k) * Dens_ZXY(k,i,j) &
              &            * (Var_WXY(k-1,i,j) - Var_WXY(k,i,j))
       enddo
!       FlxZ_ZXY(KS-1,i,j) = 0.0_RP
!       FlxZ_ZXY(KE-1,i,j) = 0.0_RP
     enddo
     enddo

     do j = JJS, JJE
     do i = IIS-1, IIE
       do k=KS, KE-1
         FlxX_WUY(k,i,j) =   DiffCoef * RFDX(i) * 0.25_RP * sum(Dens_ZXY(k:k+1,i:i+1,j)) &
              &            * (Var_WXY(k,i,j) - Var_WXY(k,i+1,j))
       enddo
       FlxX_WUY(KE,i,j) = 0.0_RP
     enddo
     enddo

     do j = JJS-1, JJE
        
     do i = IIS, IIE
       do k=KS, KE-1
         FlxY_WXV(k,i,j) =   DiffCoef * RFDY(j) * 0.25_RP * sum(Dens_ZXY(k:k+1,i,j:j+1)) &
              &            * (Var_WXY(k,i,j) - Var_WXY(k,i,j+1))
       enddo
       FlxY_WXV(KE,i,j) = 0.0_RP
     enddo
     enddo
     
  end subroutine FDM_EvalDiffFlxCD2_VarWXY
  
end module scale_atmos_numeric_fdm_diff
