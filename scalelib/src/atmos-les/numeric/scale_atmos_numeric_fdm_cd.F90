!-------------------------------------------------------------------------------
!> module Atmosphere / Numeric FDM utility (Central difference)
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

  ! Weights used in 4th order spacial interpolation.
  real(RP), parameter :: FACT_N =   7.0_RP/12.0_RP
  real(RP), parameter :: FACT_F = - 1.0_RP/12.0_RP

contains
  subroutine FDM_CD_setup(IFS_OFF_, JFS_OFF_)
    integer, intent(in) :: IFS_OFF_, JFS_OFF_

    IFS_OFF = IFS_OFF_; JFS_OFF = JFS_OFF_    
  end subroutine FDM_CD_setup

  !** --  2nd order -- *********************************************************************


  
  !** --  4th order -- *********************************************************************
  
  subroutine FDM_EvalFlxCD4_VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,             &  ! (inout)
    & Var_ZXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: i, j, k
     real(RP) :: VelZ_W(KA)
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_W) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z). 
       do k = KS, KE-1
         VelZ_W(k) = 2.0_RP * MomZ_WXY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k+1,i,j))
       enddo
       !- FluxZ at z surface within interior region      
       do k=KS+1, KE-2
          FlxZ_WXY(k,i,j) = VelZ_W(k) * (   FACT_N * (Var_ZXY(  k,i,j) + Var_ZXY(k+1,i,j))   & !
            &                             + FACT_F * (Var_ZXY(k-1,i,j) + Var_ZXY(k+2,i,j)) )   ! [ZXY -> WXY] (4th order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_WXY(KS-1,i,j) = 0.0_RP
       !- FluxZ at z surface near bottom boundary
       FlxZ_WXY(KS,i,j)   = VelZ_W(KS)   * 0.5_RP*(Var_ZXY(KS,i,j) + Var_ZXY(KS+1,i,j))
       !- FluxZ at z surface near top boundary
       FlxZ_WXY(KE-1,i,j) = VelZ_W(KE-1) * 0.5_RP*(Var_ZXY(KE-1,i,j) + Var_ZXY(KE,i,j))
       !- FluxZ at top boundary
       FlxZ_WXY(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-IFS_OFF, min(IIE,IEH)
      do k = KS, KE       
         FlxX_ZUY(k,i,j) = &
              & 2.0_RP * MomX_ZUY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i+1,j))  & ! Half level velocity(x)
              &   * (   FACT_N * (Var_ZXY(k,  i,j) + Var_ZXY(k,i+1,j))            & ! [ZXY -> ZUY] (4th order)
              &       + FACT_F * (Var_ZXY(k,i-1,j) + Var_ZXY(k,i+2,j)) )            !
      enddo
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS-JFS_OFF, min(JJE,JEH)
    do i = IIS, IIE
      do k = KS, KE
         FlxY_ZXV(k,i,j) = &
              & 2.0_RP * MomY_ZXV(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i,j+1))  & ! Half level velocity(y)
              &   * (   FACT_N * (Var_ZXY(k,i,  j) + Var_ZXY(k,i,j+1))            & ! [ZXY -> ZXV] (4th order)
              &       + FACT_F * (Var_ZXY(k,i,j-1) + Var_ZXY(k,i,j+2)) )            !
      enddo
    enddo
    enddo
   
  end subroutine FDM_EvalFlxCD4_VarZXY
  
  subroutine FDM_EvalFlxCD4_VarZUY( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,             &  ! (inout)
    & Var_ZUY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: i, j, k
     real(RP) :: VelZ_W(KA)
     
     !** Calculate z-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelZ_W) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z) and [WXY -> WUY]
       do k = KS, KE-1
         VelZ_W(k) = sum( MomZ_WXY(k,i:i+1,j) / (Dens_ZXY(k,i:i+1,j) + Dens_ZXY(k+1,i:i+1,j)) )
       enddo
       !- FluxZ at z surface within interior region      
       do k=KS+1, KE-2
          FlxZ_WUY(k,i,j) = VelZ_W(k) * (   FACT_N * (Var_ZUY(  k,i,j) + Var_ZUY(k+1,i,j))   & !
            &                             + FACT_F * (Var_ZUY(k-1,i,j) + Var_ZUY(k+2,i,j)) )   ! [ZXY -> WXY] (4th order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_WUY(KS-1,i,j) = 0.0_RP
       !- FluxZ at z surface near bottom boundary
       FlxZ_WUY(KS,i,j)   = VelZ_W(KS)   * 0.5_RP*(Var_ZUY(KS,i,j) + Var_ZUY(KS+1,i,j))
       !- FluxZ at z surface near top boundary
       FlxZ_WUY(KE-1,i,j) = VelZ_W(KE-1) * 0.5_RP*(Var_ZUY(KE-1,i,j) + Var_ZUY(KE,i,j))
       !- FluxZ at top boundary
       FlxZ_WUY(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE+1
      do k = KS, KE       
         FlxX_ZXY(k,i-1,j) = &
              &  (    MomX_ZUY(k,i-1,j) / (Dens_ZXY(k,i-1,j) + Dens_ZXY(k,  i,j))   & ! Half level velocity(x) and [ ZUY -> ZXY ]
              &     + MomX_ZUY(k,  i,j) / (Dens_ZXY(k,  i,j) + Dens_ZXY(k,i+1,j)) ) & ! 
              &  *(   FACT_N * (Var_ZUY(k,i-1,j) + Var_ZUY(k,  i,j))                & ! [ZUY -> ZXY] (4th order)
              &     + FACT_F * (Var_ZUY(k,i-2,j) + Var_ZUY(k,i+1,j)) )                !
      enddo
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      do k = KS, KE
         FlxY_ZUV(k,i,j) = &
              &  (    MomY_ZXV(k,  i,j) / (Dens_ZXY(k,  i,j) + Dens_ZXY(k,  i,j+1))   & ! Half level velocity(x) and [ ZUY -> ZXY ]
              &     + MomY_ZXV(k,i+1,j) / (Dens_ZXY(k,i+1,j) + Dens_ZXY(k,i+1,j+1)) ) & ! 
              &  *(   FACT_N * (Var_ZUY(k,i,  j) + Var_ZUY(k,i,j+1))                  & ! [ZUY -> ZXY] (4th order)
              &     + FACT_F * (Var_ZUY(k,i,j-1) + Var_ZUY(k,i,j+2)) )                  !
      enddo
    enddo
    enddo
     
  end subroutine FDM_EvalFlxCD4_VarZUY
  
  subroutine FDM_EvalFlxCD4_VarZXV( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,          &  ! (inout)
    & Var_ZXV, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: i, j, k
     real(RP) :: VelZ_WXV(KA)
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_WXV) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z) and [WXY -> WXV]
       do k = KS, KE-1
         VelZ_WXV(k) = sum( MomZ_WXY(k,i,j:j+1) / (Dens_ZXY(k,i,j:j+1) + Dens_ZXY(k+1,i,j:j+1)) )
       enddo
       !- FluxZ at z surface within interior region      
       do k=KS+1, KE-2
          FlxZ_WXV(k,i,j) = VelZ_WXV(k) * (   FACT_N * (Var_ZXV(  k,i,j) + Var_ZXV(k+1,i,j))   & !
            &                               + FACT_F * (Var_ZXV(k-1,i,j) + Var_ZXV(k+2,i,j)) )   ! [ZXV -> WXV] (4th order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_WXV(KS-1,i,j) = 0.0_RP
       !- FluxZ at z surface near bottom boundary
       FlxZ_WXV(KS,i,j)   = VelZ_WXV(KS)   * 0.5_RP*(Var_ZXV(KS,i,j) + Var_ZXV(KS+1,i,j))
       !- FluxZ at z surface near top boundary
       FlxZ_WXV(KE-1,i,j) = VelZ_WXV(KE-1) * 0.5_RP*(Var_ZXV(KE-1,i,j) + Var_ZXV(KE,i,j))
       !- FluxZ at top boundary
       FlxZ_WXV(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      do k = KS, KE       
         FlxX_ZUV(k,i,j) = &
              &  (    MomX_ZUY(k,i,j-1) / (Dens_ZXY(k,i,j-1) + Dens_ZXY(k,i+1,j-1))   & ! Half level velocity(x) and [ ZUY -> ZUV ]
              &     + MomX_ZUY(k,i,  j) / (Dens_ZXY(k,i,j  ) + Dens_ZXY(k,i+1,j))   ) & ! 
              &  *(   FACT_N * (Var_ZXV(k,i  ,j) + Var_ZXV(k,i+1,j))                  & ! [ZXV -> ZUV] (4th order)
              &     + FACT_F * (Var_ZXV(k,i-1,j) + Var_ZXV(k,i+2,j)) )                  !
      enddo
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE+1
    do i = IIS, IIE
      do k = KS, KE
         FlxY_ZXY(k,i,j-1) = &
              &  (    MomY_ZXV(k,i,j-1) / (Dens_ZXY(k,i,j-1) + Dens_ZXY(k,i,j))       & ! Half level velocity(x) and [ ZXV -> ZXY ]
              &     + MomY_ZXV(k,i,j  ) / (Dens_ZXY(k,i,j-1) + Dens_ZXY(k,i,j)) )     & ! 
              &  *(   FACT_N * (Var_ZXV(k,i,j-1) + Var_ZXV(k,i,  j))                  & ! [ZXV -> ZXY] (4th order)
              &     + FACT_F * (Var_ZXV(k,i,j-2) + Var_ZXv(k,i,j+1)) )                  !
      enddo
    enddo
    enddo
           
  end subroutine FDM_EvalFlxCD4_VarZXV

  subroutine FDM_EvalFlxCD4_VarWXY( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,          &  ! (inout)
    & Var_WXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE    

     integer :: i, j, k
     real(RP) :: VelZ_ZXY(KA), sw
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_ZXY) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z) and [WXY -> ZXY]
       do k = KS+1, KE-1
         VelZ_ZXY(k-1) = 0.5_RP * sum(MomZ_WXY(k-1:k,i,j)) / Dens_ZXY(k,i,j)
       enddo
       !- FluxZ at z surface within interior region      
       do k = KS+2, KE-2
          FlxZ_ZXY(k-1,i,j) = VelZ_ZXY(k) * (   FACT_N * (Var_WXY(k-1,i,j) + Var_WXY(  k,i,j))   & !
            &                                 + FACT_F * (Var_WXY(k-2,i,j) + Var_WXY(k+1,i,j)) )   ! [WXY -> ZXV] (4th order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_ZXY(KS-1,i,j) = 0.0_RP
       !- FluxZ at z surface near bottom boundary
       FlxZ_ZXY(KS,i,j)   = VelZ_ZXY(KS)   * 0.5_RP*(Var_WXY(KS,i,j)   + Var_WXY(KS+1,i,j))
       !- FluxZ at z surface near top boundary
       FlxZ_ZXY(KE-2,i,j) = VelZ_ZXY(KE-2) * 0.5_RP*(Var_WXY(KE-2,i,j) + Var_WXY(KE-1,i,j))
       !- FluxZ at top boundary
       FlxZ_ZXY(KE-1,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      do k = KS, KE-1
         FlxX_WUY(k,i,j) = &
              &  (    MomX_ZUY(  k,i,j) / (Dens_ZXY(  k,i,j) + Dens_ZXY(  k,i+1,j))     & ! Half level velocity(x) and [ ZUY -> WUY ]
              &     + MomX_ZUY(k+1,i,j) / (Dens_ZXY(k+1,i,j) + Dens_ZXY(k+1,i+1,j))   ) & ! 
              &  *(   FACT_N * (Var_WXY(k,i  ,j) + Var_WXY(k,i+1,j))                    & ! [ZXV -> ZUV] (4th order)
              &     + FACT_F * (Var_WXY(k,i-1,j) + Var_WXY(k,i+2,j)) )                    !
      enddo
      FlxX_WUY(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      do k = KS, KE-1
         FlxY_WXV(k,i,j) = &
              &  (    MomY_ZXV(  k,i,j) / (Dens_ZXY(  k,i,j) + Dens_ZXY(  k,i,j+1))       & ! Half level velocity(y) and [ ZYV -> WYV ]
              &     + MomY_ZXV(k+1,i,j) / (Dens_ZXY(k+1,i,j) + Dens_ZXY(k+1,i,j+1)) )     & ! 
              &  *(   FACT_N * (Var_WXY(k,i,  j) + Var_WXY(k,i,j+1))                      & ! [WXY -> WXV] (4th order)
              &     + FACT_F * (Var_WXY(k,i,j-1) + Var_WXY(k,i,j+2)) )                      !
      enddo
      FlxY_WXV(KE,i,j) = 0.0_RP      
    enddo
    enddo

  end subroutine FDM_EvalFlxCD4_VarWXY

end module scale_atmos_numeric_fdm_cd
