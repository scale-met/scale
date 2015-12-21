!-------------------------------------------------------------------------------
!> module Atmosphere / Numeric FDM utility (Upwind difference)
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
module scale_atmos_numeric_fdm_ud
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

  public :: FDM_UD_setup

  !**
  public :: FDM_EvalFlxUD1_VarZXY
  public :: FDM_EvalFlxUD1_VarZUY
  public :: FDM_EvalFlxUD1_VarZXV
  public :: FDM_EvalFlxUD1_VarWXY
  
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
  subroutine FDM_UD_setup(IFS_OFF_, JFS_OFF_)
    integer, intent(in) :: IFS_OFF_, JFS_OFF_

    IFS_OFF = IFS_OFF_; JFS_OFF = JFS_OFF_
    
  end subroutine FDM_UD_setup

  !** --  1st order -- *********************************************************************
  
  subroutine FDM_EvalFlxUD1_VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,             &  ! (inout)
    & Var_ZXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: i, j, k
     real(RP) :: VelZ_W(KA), VelX_U(KA), VelY_V(KA)
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_W) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z). 
       do k = KS, KE-1
         VelZ_W(k) = 2.0_RP * MomZ_WXY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k+1,i,j))
       enddo
       !- FluxZ at z surface within interior region      
       do k=KS, KE-1
          FlxZ_WXY(k,i,j) = 0.5_RP * (       VelZ_W(k)  * (Var_ZXY(k+1,i,j) + Var_ZXY(k,i,j))   & !
            &                          - abs(VelZ_W(k)) * (Var_ZXY(k+1,i,j) - Var_ZXY(k,i,j)) )   ! [ZXY -> WXY] (1st order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_WXY(KS-1,i,j) = 0.0_RP
       !- FluxZ at top boundary
       FlxZ_WXY(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX_U) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-IFS_OFF, min(IIE,IEH)
      ! Half level velocity(x)
      do k = KS, KE
         VelX_U(k) = 2.0_RP * MomX_ZUY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i+1,j))  
      enddo
      do k = KS, KE
        FlxX_ZUY(k,i,j) =  0.5_RP * (       VelX_U(k)  * (Var_ZXY(k,i+1,j) + Var_ZXY(k,i,j))   & !
           &                          - abs(VelX_U(k)) * (Var_ZXY(k,i+1,j) - Var_ZXY(k,i,j)) )   ! [ZXY -> ZUY] (1st order)      
      enddo
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY_V) OMP_SCHEDULE_ collapse(2)
    do j = JJS-JFS_OFF, min(JJE,JEH)
    do i = IIS, IIE
      ! Half level velocity(x)
      do k = KS, KE
         VelY_V(k) = 2.0_RP * MomY_ZXV(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i,j+1))  
      enddo
      do k = KS, KE
        FlxY_ZXV(k,i,j) =  0.5_RP * (       VelY_V(k)  * (Var_ZXY(k,i,j+1) + Var_ZXY(k,i,j))   & !
           &                          - abs(VelY_V(k)) * (Var_ZXY(k,i,j+1) - Var_ZXY(k,i,j)) )   ! [ZXY -> ZXV] (1st order)      
      enddo
    enddo
    enddo
   
  end subroutine FDM_EvalFlxUD1_VarZXY
  
  subroutine FDM_EvalFlxUD1_VarZUY( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,             &  ! (inout)
    & Var_ZUY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: i, j, k
     real(RP) :: VelZ_W(KA), VelX_X(KA), VelY_V(KA)
     
     !** Calculate z-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelZ_W) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z) and [WXY -> WUY]
       do k = KS, KE-1
         VelZ_W(k) = sum( MomZ_WXY(k,i:i+1,j) / (Dens_ZXY(k,i:i+1,j) + Dens_ZXY(k+1,i:i+1,j)) )
       enddo
       !- FluxZ at z surface within interior region      
       do k=KS, KE-1
          FlxZ_WUY(k,i,j) = 0.5_RP * (       VelZ_W(k)  * (Var_ZUY(k+1,i,j) + Var_ZUY(k,i,j))   & !
            &                          - abs(VelZ_W(k)) * (Var_ZUY(k+1,i,j) - Var_ZUY(k,i,j)) )   ! [ZXY -> WXY] (1st order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_WUY(KS-1,i,j) = 0.0_RP
       !- FluxZ at top boundary
       FlxZ_WUY(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX_X) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE+1
      ! Half level velocity(x)
      do k = KS, KE
         VelX_X(k) = (    MomX_ZUY(k,i-1,j) / (Dens_ZXY(k,i-1,j) + Dens_ZXY(k,  i,j))   & ! Half level velocity(x) and [ ZUY -> ZXY ]
                  &     + MomX_ZUY(k,  i,j) / (Dens_ZXY(k,  i,j) + Dens_ZXY(k,i+1,j)) )   !
      enddo
      do k = KS, KE
        FlxX_ZXY(k,i-1,j) =  0.5_RP * (       VelX_X(k)  * (Var_ZUY(k,i,j) + Var_ZUY(k,i-1,j))   & !
           &                            - abs(VelX_X(k)) * (Var_ZUY(k,i,j) - Var_ZUY(k,i-1,j)) )   ! [ZUY -> ZXY] (1st order)      
      enddo
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY_V) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      ! Half level velocity(y)
      do k = KS, KE
         VelY_V(k) =  (    MomY_ZXV(k,  i,j) / (Dens_ZXY(k,  i,j) + Dens_ZXY(k,  i,j+1))   & ! Half level velocity(x) and [ ZXV -> ZUV ]
                   &     + MomY_ZXV(k,i+1,j) / (Dens_ZXY(k,i+1,j) + Dens_ZXY(k,i+1,j+1)) )   ! 
      enddo
      do k = KS, KE
        FlxY_ZUV(k,i,j) =  0.5_RP * (       VelY_V(k)  * (Var_ZUY(k,i,j+1) + Var_ZUY(k,i,j))    &  !
         &                            - abs(VelY_V(k)) * (Var_ZUY(k,i,j+1) - Var_ZUY(k,i,j)) )     ! [ZUY -> ZUV] (1st order)      
      enddo
    enddo
    enddo
     
  end subroutine FDM_EvalFlxUD1_VarZUY
  
  subroutine FDM_EvalFlxUD1_VarZXV( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,          &  ! (inout)
    & Var_ZXV, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE

     integer :: i, j, k
     real(RP) :: VelZ_WXV(KA), VelX_ZUV(KA), VelY_ZXY(KA)
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_WXV) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z) and [WXY -> WXV]
       do k = KS, KE-1
         VelZ_WXV(k) = sum( MomZ_WXY(k,i,j:j+1) / (Dens_ZXY(k,i,j:j+1) + Dens_ZXY(k+1,i,j:j+1)) )
       enddo
       !- FluxZ at z surface within interior region      
       do k=KS, KE-1
          FlxZ_WXV(k,i,j) = 0.5_RP * (       VelZ_WXV(k)  * (Var_ZXV(k+1,i,j) + Var_ZXV(k,i,j))   & !
            &                          - abs(VelZ_WXV(k)) * (Var_ZXV(k+1,i,j) - Var_ZXV(k,i,j)) )   ! [ZXY -> WXY] (1st order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_WXV(KS-1,i,j) = 0.0_RP
       !- FluxZ at top boundary
       FlxZ_WXV(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX_ZUV) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      ! Half level velocity(x)
      do k = KS, KE
         VelX_ZUV(k) = (    MomX_ZUY(k,i,j-1) / (Dens_ZXY(k,i,j-1) + Dens_ZXY(k,i+1,j-1))   & ! Half level velocity(x) and [ ZUY -> ZUV ]
                    &     + MomX_ZUY(k,i,  j) / (Dens_ZXY(k,i,j  ) + Dens_ZXY(k,i+1,j))   )   !          
      enddo
      do k = KS, KE
        FlxX_ZUV(k,i,j) =  0.5_RP * (       VelX_ZUV(k)  * (Var_ZXV(k,i+1,j) + Var_ZXV(k,i,j))   & !
           &                          - abs(VelX_ZUV(k)) * (Var_ZXV(k,i+1,j) - Var_ZXV(k,i,j)) )   ! [ZUY -> ZXY] (1st order)      
      enddo
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY_ZXY) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE+1
    do i = IIS, IIE
       do k = KS, KE
         VelY_ZXY(k) = (    MomY_ZXV(k,i,j-1) / (Dens_ZXY(k,i,j-1) + Dens_ZXY(k,i,j))       & ! Half level velocity(y) and [ ZXV -> ZXY ]
                    &     + MomY_ZXV(k,i,j  ) / (Dens_ZXY(k,i,j-1) + Dens_ZXY(k,i,j)) )       ! 
       enddo
       do k = KS, KE
         FlxY_ZXY(k,i,j-1) =  0.5_RP * (       VelY_ZXY(k)  * (Var_ZXV(k,i,j) + Var_ZXV(k,i,j-1))    &  !
            &                            - abs(VelY_ZXY(k)) * (Var_ZXV(k,i,j) - Var_ZXV(k,i,j-1)) )     ! [ZUY -> ZUV] (1st order)      
       enddo
    enddo
    enddo
           
  end subroutine FDM_EvalFlxUD1_VarZXV

  subroutine FDM_EvalFlxUD1_VarWXY( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,          &  ! (inout)
    & Var_WXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & IIS, IIE, JJS, JJE, KS, KE )                                               ! (in)

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE    

     integer :: i, j, k
     real(RP) :: VelZ_ZXY(KA), VelX_WUY(KA), VelY_WXV(KA), sw
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_ZXY) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       !- Half-level velocity(z) and [WXY -> ZXY]
       do k = KS+1, KE-1
         VelZ_ZXY(k-1) = 0.5_RP * sum(MomZ_WXY(k-1:k,i,j)) / Dens_ZXY(k,i,j)
       enddo
       !- FluxZ at z surface within interior region      
       do k=KS+1, KE-2
          FlxZ_ZXY(k-1,i,j) = 0.5_RP * (       VelZ_ZXY(k)  * (Var_WXY(k,i,j) + Var_WXY(k-1,i,j))   & !
            &                            - abs(VelZ_ZXY(k)) * (Var_WXY(k,i,j) - Var_WXY(k-1,i,j)) )   ! [ZXY -> WXY] (1st order)
       enddo
       !- FluxZ at bottom boundary
       FlxZ_ZXY(KS-1,i,j) = 0.0_RP
       !- FluxZ at top boundary
       FlxZ_ZXY(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,Vel_WUY) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      ! Half level velocity(x)
      do k = KS, KE-1
         VelX_WUY(k) = (    MomX_ZUY(  k,i,j) / (Dens_ZXY(  k,i,j) + Dens_ZXY(k,i+1,j))    & ! Half level velocity(x) and [ ZUY -> WUY ]
                    &     + MomX_ZUY(k+1,i,j) / (Dens_ZXY(k+1,i,j) + Dens_ZXY(k,i+1,j))   )  !          
      enddo
      do k = KS, KE-1
        FlxX_WUY(k,i,j) =  0.5_RP * (       VelX_WUY(k)  * (Var_WXY(k,i+1,j) + Var_WXY(k,i,j))   & !
           &                          - abs(VelX_WUY(k)) * (Var_WXY(k,i+1,j) - Var_WXY(k,i,j)) )   ! [WXY -> WUY] (1st order)      
      enddo
      FlxX_WUY(KE,i,j) = 0.0_RP
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY_WXV) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      ! Half level velocity(x)
      do k = KS, KE-1
         VelY_WXV(k) = (    MomY_ZXV(  k,i,j) / (Dens_ZXY(  k,i,j) + Dens_ZXY(  k,i,j+1))    & ! Half level velocity(x) and [ ZXV -> WXV ]
                    &     + MomY_ZXV(k+1,i,j) / (Dens_ZXY(k+1,i,j) + Dens_ZXY(k+1,i,j+1))   )  !          
      enddo
      do k = KS, KE-1
        FlxY_WXV(k,i,j) =  0.5_RP * (        VelY_WXV(k)  * (Var_WXY(k,i+1,j) + Var_WXY(k,i,j))   & !
           &                           - abs(VelY_WXV(k)) * (Var_WXY(k,i+1,j) - Var_WXY(k,i,j)) )   ! [WXY -> WUY] (1st order)      
      enddo
      FlxY_WXV(KE,i,j) = 0.0_RP      
    enddo
    enddo

  end subroutine FDM_EvalFlxUD1_VarWXY

end module scale_atmos_numeric_fdm_ud
