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

  use scale_atmos_numeric_fdm_def, only: &
       & VL_XY, VL_UY, VL_XV, &
       & VL_ZXY, VL_ZUY, VL_ZXV, VL_WXY

  use scale_atmos_numeric_fdm_util, only: &
       & FDM_AddDiffFlxX, &
       & FDM_AddDiffFlxY, &
       & FDM_AddDiffFlxZ, &       
       & FDM_MomX2VelXt, &
       & FDM_MomY2VelYt, &
       & FDM_MomZ2VelZt
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  public :: FDM_CD_setup

  !** 2nd order
  public :: FDM_EvalFlxCD2_VarZXY
  public :: FDM_EvalFlxCD2_VarZUY
  public :: FDM_EvalFlxCD2_VarZXV
  public :: FDM_EvalFlxCD2_VarWXY
  
  !** 4th order
  public :: FDM_EvalFlxCD4_VarZXY
  public :: FDM_EvalFlxCD4_VarZUY
  public :: FDM_EvalFlxCD4_VarZXV
  public :: FDM_EvalFlxCD4_VarWXY

  !** 6th order
  public :: FDM_EvalFlxCD6_VarZXY
  public :: FDM_EvalFlxCD6_VarZUY
  public :: FDM_EvalFlxCD6_VarZXV
  public :: FDM_EvalFlxCD6_VarWXY
  
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

  ! F_1/2 = (u_1/2*FAC1UD?) * slice( phi )
  !
  real(RP), parameter :: FAC1CD2(2) &
       & = (/ 1.0_RP,  1.0_RP /) / 2.0_RP
  
  real(RP), parameter :: FAC1CD4(4) &
       & = (/ -1.0_RP, 7.0_RP,  7.0_RP, -1.0_RP /) / 12.0_RP

  real(RP), parameter :: FAC1CD6(6) &
       & = (/ 1.0_RP, -8.0_RP, 37.0_RP,  37.0_RP, -8.0_RP,   1.0_RP /) / 60.0_RP
  
contains
  subroutine FDM_CD_setup(IFS_OFF_, JFS_OFF_)
    integer, intent(in) :: IFS_OFF_, JFS_OFF_

    IFS_OFF = IFS_OFF_; JFS_OFF = JFS_OFF_    
  end subroutine FDM_CD_setup

  !** --  2nd order -- *********************************************************************

  subroutine EvalFlxX_CD2(FlxX, VelX, Var, IOF, i, j, KS_, KE_)
    real(RP), intent(out)   :: FlxX(KA)
    real(RP), intent(in)    :: VelX(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: IOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxX(k) =   sum( VelX(k)*FAC1CD2(:)*Var(k,i+IOF:i+1+IOF,j) )           
    enddo
  end subroutine EvalFlxX_CD2

  subroutine EvalFlxY_CD2(FlxY, VelY, Var, JOF, i, j, KS_, KE_)
    real(RP), intent(out)    :: FlxY(KA)
    real(RP), intent(in)    :: VelY(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: JOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxY(k) =   sum( VelY(k)*FAC1CD2(:)*Var(k,i,j+JOF:j+1+JOF) )           
    enddo    
  end subroutine EvalFlxY_CD2

  subroutine EvalFlxZ_CD2(FlxZ, VelZ, Var, KOF, i, j, KS_, KE_)
    real(RP), intent(out)   :: FlxZ(KA)
    real(RP), intent(in)    :: VelZ(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: KOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
      FlxZ(k+KOF) = sum( VelZ(k+KOF)*FAC1CD2(:)*Var(k+KOF:k+1+KOF,i,j) )
    enddo
  end subroutine EvalFlxZ_CD2

  subroutine FDM_EvalFlxCD2_VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,             &  ! (inout)
    & Var_ZXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZUY, DiffFlxY_ZXV, DiffFlxZ_WXY )                                 ! (in)   

    
     real(RP), intent(out), dimension(KA,IA,JA)    :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Var_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer,  intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZUY, DiffFlxY_ZXV, DiffFlxZ_WXY

     real(RP), dimension(KA) :: VelZ, VelX, VelY
     integer :: i, j, k
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_W) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WXY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j, KS, KE-1)              
       FlxZ_WXY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WXY)) then
         call FDM_AddDiffFlxZ(FlxZ_WXY, DiffFlxZ_WXY, GSQRT(:,:,:,I_XYW), MAPF(:,:,:,I_XY), 0, i, j, KS, KE-1)
       end if
     enddo
     enddo

     !** Calculate x-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS-IFS_OFF, min(IIE,IEH)           
       call FDM_MomX2VelXt( VelX,                                        & ! (out)
         & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
         & GSQRT, MAPF, GeneralCoordFlag )

       call EvalFlxX_CD2(FlxX_ZUY(:,i,j), VelX, Var_ZXY, 0, i, j, KS, KE)              

       if (present(DiffFlxX_ZUY)) then
         call FDM_AddDiffFlxX(FlxX_ZUY, DiffFlxX_ZUY, GSQRT(:,:,:,I_UYZ), MAPF(:,:,:,I_UY), i, j, KS, KE)
       end if
     enddo
     enddo

     !** Calculate y-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
     do j = JJS-JFS_OFF, min(JJE,JEH)
     do i = IIS, IIE
       call FDM_MomY2VelYt( VelY,                                        & ! (out)
         & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
         & GSQRT, MAPF, GeneralCoordFlag )

       call EvalFlxY_CD2(FlxY_ZXV(:,i,j), VelY, Var_ZXY, 0, i, j, KS, KE)              

       if (present(DiffFlxY_ZXV)) then
         call FDM_AddDiffFlxY(FlxY_ZXV, DiffFlxY_ZXV, GSQRT(:,:,:,I_XVZ), MAPF(:,:,:,I_XV), i, j, KS, KE)
       end if
     enddo
     enddo
   
  end subroutine FDM_EvalFlxCD2_VarZXY
  
  subroutine FDM_EvalFlxCD2_VarZUY( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,             &  ! (inout)
    & Var_ZUY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZXY, DiffFlxY_ZUV, DiffFlxZ_WUY )                                 ! (in)   

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZXY, DiffFlxY_ZUV, DiffFlxZ_WUY     

     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY
     
     !** Calculate z-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WUY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j, KS, KE-1)              
       FlxZ_WUY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WUY)) then
          call FDM_AddDiffFlxZ(FlxZ_WUY, DiffFlxZ_WUY, GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), 0, i, j, KS, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE+1
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD2(FlxX_ZXY(:,i-1,j), VelX, Var_ZUY, -1, i, j, KS, KE)              

      if (present(DiffFlxX_ZXY)) then
        call FDM_AddDiffFlxX(FlxX_ZXY, DiffFlxX_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), i, j, KS, KE)
      end if
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD2(FlxY_ZUV(:,i,j), VelY, Var_ZUY, 0, i, j, KS, KE)              

      if (present(DiffFlxY_ZUV)) then
        call FDM_AddDiffFlxY(FlxY_ZUV, DiffFlxY_ZUV, GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), i, j, KS, KE)
      end if
    enddo
    enddo

  end subroutine FDM_EvalFlxCD2_VarZUY
  
  subroutine FDM_EvalFlxCD2_VarZXV( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,          &  ! (inout)
    & Var_ZXV, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZUV, DiffFlxY_ZXY, DiffFlxZ_WXV )                                 ! (in)       
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag     
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZUV, DiffFlxY_ZXY, DiffFlxZ_WXV
     
     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WXV(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j, KS, KE-1)              
       FlxZ_WXV(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WXV)) then
          call FDM_AddDiffFlxZ(FlxZ_WXV, DiffFlxZ_WXV, GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), 0, i, j, KS, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD2(FlxX_ZUV(:,i,j), VelX, Var_ZXV, 0, i, j, KS, KE)              

      if (present(DiffFlxX_ZUV)) then
        call FDM_AddDiffFlxX(FlxX_ZUV, DiffFlxX_ZUV, GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), i, j, KS, KE)
      end if
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE+1
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD2(FlxY_ZXY(:,i,j-1), VelY, Var_ZXV, -1, i, j, KS, KE)              

      if (present(DiffFlxY_ZXY)) then
        call FDM_AddDiffFlxY(FlxY_ZXY, DiffFlxY_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), i, j, KS, KE)
      end if
    enddo
    enddo
           
  end subroutine FDM_EvalFlxCD2_VarZXV

  subroutine FDM_EvalFlxCD2_VarWXY( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,             &  ! (inout)
    & Var_WXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)          
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_WUY, DiffFlxY_WXV, DiffFlxZ_ZXY )                                 ! (in)   
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE    
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_WUY, DiffFlxY_WXV, DiffFlxZ_ZXY

     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY

     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_ZXY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KS+1, KE-1)              
       FlxZ_ZXY(KE-1,i,j) = 0.0_RP
!       FlxZ_ZXY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_ZXY)) then
         call FDM_AddDiffFlxZ(FlxZ_ZXY, DiffFlxZ_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), -1, i, j, KS+1, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD2(FlxX_WUY(:,i,j), VelX, Var_WXY, 0, i, j, KS, KE-1)              
      FlxX_WUY(KE,i,j) = 0.0_RP

      if (present(DiffFlxX_WUY)) then
        call FDM_AddDiffFlxX(FlxX_WUY, DiffFlxX_WUY, GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), i, j, KS, KE-1)
      end if              
    enddo
    enddo
    
    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD2(FlxY_WXV(:,i,j), VelY, Var_WXY, 0, i, j, KS, KE-1)              
      FlxY_WXV(KE,i,j) = 0.0_RP
      
      if (present(DiffFlxY_WXV)) then
        call FDM_AddDiffFlxY(FlxY_WXV, DiffFlxY_WXV, GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), i, j, KS, KE-1)
      end if      
    enddo
    enddo

  end subroutine FDM_EvalFlxCD2_VarWXY  
  
  !** --  4th order -- *********************************************************************

  subroutine EvalFlxX_CD4(FlxX, VelX, Var, IOF, i, j, KS_, KE_)
    real(RP), intent(out) :: FlxX(KA)
    real(RP), intent(in)    :: VelX(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: IOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxX(k) =   sum( VelX(k)*FAC1CD4(:)*Var(k,i-1+IOF:i+2+IOF,j) )           
    enddo
  end subroutine EvalFlxX_CD4

  subroutine EvalFlxY_CD4(FlxY, VelY, Var, JOF, i, j, KS_, KE_)
    real(RP), intent(out) :: FlxY(KA)
    real(RP), intent(in)    :: VelY(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: JOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxY(k) =   sum( VelY(k)*FAC1CD4(:)*Var(k,i,j-1+JOF:j+2+JOF) )           
    enddo    
  end subroutine EvalFlxY_CD4

  subroutine EvalFlxZ_CD4(FlxZ, VelZ, Var, KOF, i, j, KS_, KE_)
    real(RP), intent(out) :: FlxZ(KA)
    real(RP), intent(in)    :: VelZ(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: KOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
      FlxZ(k+KOF) = sum( VelZ(k+KOF)*FAC1CD4(:)*Var(k-1+KOF:k+2+KOF,i,j) )
    enddo
  end subroutine EvalFlxZ_CD4
  
  subroutine FDM_EvalFlxCD4_VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,             &  ! (inout)
    & Var_ZXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZUY, DiffFlxY_ZXV, DiffFlxZ_WXY )                                 ! (in)   

    
     real(RP), intent(out), dimension(KA,IA,JA)    :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Var_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer,  intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZUY, DiffFlxY_ZXV, DiffFlxZ_WXY

     real(RP), dimension(KA) :: VelZ, VelX, VelY
     integer :: i, j, k
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_W) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WXY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j,   KS,   KS)
       call EvalFlxZ_CD4(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j, KS+1, KE-2)              
       call EvalFlxZ_CD2(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j, KE-1, KE-1)
       FlxZ_WXY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WXY)) then
         call FDM_AddDiffFlxZ(FlxZ_WXY, DiffFlxZ_WXY, GSQRT(:,:,:,I_XYW), MAPF(:,:,:,I_XY), 0, i, j, KS, KE-1)
       end if
     enddo
     enddo

     !** Calculate x-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS-IFS_OFF, min(IIE,IEH)           
       call FDM_MomX2VelXt( VelX,                                        & ! (out)
         & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
         & GSQRT, MAPF, GeneralCoordFlag )

       call EvalFlxX_CD4(FlxX_ZUY(:,i,j), VelX, Var_ZXY, 0, i, j, KS, KE)              

       if (present(DiffFlxX_ZUY)) then
         call FDM_AddDiffFlxX(FlxX_ZUY, DiffFlxX_ZUY, GSQRT(:,:,:,I_UYZ), MAPF(:,:,:,I_UY), i, j, KS, KE)
       end if
     enddo
     enddo

     !** Calculate y-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
     do j = JJS-JFS_OFF, min(JJE,JEH)
     do i = IIS, IIE
       call FDM_MomY2VelYt( VelY,                                        & ! (out)
         & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
         & GSQRT, MAPF, GeneralCoordFlag )

       call EvalFlxY_CD4(FlxY_ZXV(:,i,j), VelY, Var_ZXY, 0, i, j, KS, KE)              

       if (present(DiffFlxY_ZXV)) then
         call FDM_AddDiffFlxY(FlxY_ZXV, DiffFlxY_ZXV, GSQRT(:,:,:,I_XVZ), MAPF(:,:,:,I_XV), i, j, KS, KE)
       end if
     enddo
     enddo
   
  end subroutine FDM_EvalFlxCD4_VarZXY
  
  subroutine FDM_EvalFlxCD4_VarZUY( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,             &  ! (inout)
    & Var_ZUY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZXY, DiffFlxY_ZUV, DiffFlxZ_WUY )                                 ! (in)   

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZXY, DiffFlxY_ZUV, DiffFlxZ_WUY     

     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY
     
     !** Calculate z-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WUY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j,   KS,   KS)
       call EvalFlxZ_CD4(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j, KS+1, KE-2)              
       call EvalFlxZ_CD2(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j, KE-1, KE-1)
       FlxZ_WUY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WUY)) then
          call FDM_AddDiffFlxZ(FlxZ_WUY, DiffFlxZ_WUY, GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), 0, i, j, KS, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE+1
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD4(FlxX_ZXY(:,i-1,j), VelX, Var_ZUY, -1, i, j, KS, KE)              

      if (present(DiffFlxX_ZXY)) then
        call FDM_AddDiffFlxX(FlxX_ZXY, DiffFlxX_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), i, j, KS, KE)
      end if
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD4(FlxY_ZUV(:,i,j), VelY, Var_ZUY, 0, i, j, KS, KE)              

      if (present(DiffFlxY_ZUV)) then
        call FDM_AddDiffFlxY(FlxY_ZUV, DiffFlxY_ZUV, GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), i, j, KS, KE)
      end if
    enddo
    enddo

  end subroutine FDM_EvalFlxCD4_VarZUY
  
  subroutine FDM_EvalFlxCD4_VarZXV( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,          &  ! (inout)
    & Var_ZXV, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZUV, DiffFlxY_ZXY, DiffFlxZ_WXV )                                 ! (in)       
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag     
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZUV, DiffFlxY_ZXY, DiffFlxZ_WXV
     
     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WXV(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j,   KS,   KS)
       call EvalFlxZ_CD4(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j, KS+1, KE-2)              
       call EvalFlxZ_CD2(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j, KE-1, KE-1)
       FlxZ_WXV(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WXV)) then
          call FDM_AddDiffFlxZ(FlxZ_WXV, DiffFlxZ_WXV, GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), 0, i, j, KS, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD4(FlxX_ZUV(:,i,j), VelX, Var_ZXV, 0, i, j, KS, KE)              

      if (present(DiffFlxX_ZUV)) then
        call FDM_AddDiffFlxX(FlxX_ZUV, DiffFlxX_ZUV, GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), i, j, KS, KE)
      end if
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE+1
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD4(FlxY_ZXY(:,i,j-1), VelY, Var_ZXV, -1, i, j, KS, KE)              

      if (present(DiffFlxY_ZXY)) then
        call FDM_AddDiffFlxY(FlxY_ZXY, DiffFlxY_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), i, j, KS, KE)
      end if
    enddo
    enddo
           
  end subroutine FDM_EvalFlxCD4_VarZXV

  subroutine FDM_EvalFlxCD4_VarWXY( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,             &  ! (inout)
    & Var_WXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)          
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_WUY, DiffFlxY_WXV, DiffFlxZ_ZXY )                                 ! (in)   
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE    
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_WUY, DiffFlxY_WXV, DiffFlxZ_ZXY

     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY

     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_ZXY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KS+1, KS+1)
       call EvalFlxZ_CD4(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KS+2, KE-2)              
       call EvalFlxZ_CD2(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KE-1, KE-1)
       FlxZ_ZXY(KE-1,i,j) = 0.0_RP
!       FlxZ_ZXY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_ZXY)) then
         call FDM_AddDiffFlxZ(FlxZ_ZXY, DiffFlxZ_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), -1, i, j, KS+1, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD4(FlxX_WUY(:,i,j), VelX, Var_WXY, 0, i, j, KS, KE-1)              
      FlxX_WUY(KE,i,j) = 0.0_RP

      if (present(DiffFlxX_WUY)) then
        call FDM_AddDiffFlxX(FlxX_WUY, DiffFlxX_WUY, GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), i, j, KS, KE-1)
      end if              
    enddo
    enddo
    
    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD4(FlxY_WXV(:,i,j), VelY, Var_WXY, 0, i, j, KS, KE-1)              
      FlxY_WXV(KE,i,j) = 0.0_RP
      
      if (present(DiffFlxY_WXV)) then
        call FDM_AddDiffFlxY(FlxY_WXV, DiffFlxY_WXV, GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), i, j, KS, KE-1)
      end if      
    enddo
    enddo

  end subroutine FDM_EvalFlxCD4_VarWXY

  !** --  6th order -- *********************************************************************

  subroutine EvalFlxX_CD6(FlxX, VelX, Var, IOF, i, j, KS_, KE_)
    real(RP), intent(out) :: FlxX(KA)
    real(RP), intent(in)    :: VelX(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: IOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxX(k) =   sum( VelX(k)*FAC1CD6(:)*Var(k,i-2+IOF:i+3+IOF,j) )           
    enddo
  end subroutine EvalFlxX_CD6

  subroutine EvalFlxY_CD6(FlxY, VelY, Var, JOF, i, j, KS_, KE_)
    real(RP), intent(out) :: FlxY(KA)
    real(RP), intent(in)    :: VelY(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: JOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxY(k) =   sum( VelY(k)*FAC1CD6(:)*Var(k,i,j-2+JOF:j+3+JOF) )           
    enddo    
  end subroutine EvalFlxY_CD6

  subroutine EvalFlxZ_CD6(FlxZ, VelZ, Var, KOF, i, j, KS_, KE_)
    real(RP), intent(out) :: FlxZ(KA)
    real(RP), intent(in)    :: VelZ(KA), Var(KA,IA,JA)
    integer,  intent(in)    :: KOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
      FlxZ(k+KOF) = sum( VelZ(k+KOF)*FAC1CD6(:)*Var(k-2+KOF:k+3+KOF,i,j) )
    enddo
  end subroutine EvalFlxZ_CD6
  
  subroutine FDM_EvalFlxCD6_VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,             &  ! (inout)
    & Var_ZXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZUY, DiffFlxY_ZXV, DiffFlxZ_WXY )                                 ! (in)   

    
     real(RP), intent(out), dimension(KA,IA,JA)    :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Var_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer,  intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZUY, DiffFlxY_ZXV, DiffFlxZ_WXY

     real(RP), dimension(KA) :: VelZ, VelX, VelY
     integer :: i, j, k
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ_W) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WXY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j,   KS,   KS)
       call EvalFlxZ_CD4(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j, KS+1, KS+1)
       call EvalFlxZ_CD6(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j, KS+2, KE-3)
       call EvalFlxZ_CD4(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j, KS-2, KS-2)
       call EvalFlxZ_CD2(FlxZ_WXY(:,i,j), VelZ, Var_ZXY, 0, i, j, KE-1, KE-1)
       FlxZ_WXY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WXY)) then
         call FDM_AddDiffFlxZ(FlxZ_WXY, DiffFlxZ_WXY, GSQRT(:,:,:,I_XYW), MAPF(:,:,:,I_XY), 0, i, j, KS, KE-1)
       end if
     enddo
     enddo

     !** Calculate x-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS-IFS_OFF, min(IIE,IEH)           
       call FDM_MomX2VelXt( VelX,                                        & ! (out)
         & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
         & GSQRT, MAPF, GeneralCoordFlag )

       call EvalFlxX_CD6(FlxX_ZUY(:,i,j), VelX, Var_ZXY, 0, i, j, KS, KE)              

       if (present(DiffFlxX_ZUY)) then
         call FDM_AddDiffFlxX(FlxX_ZUY, DiffFlxX_ZUY, GSQRT(:,:,:,I_UYZ), MAPF(:,:,:,I_UY), i, j, KS, KE)
       end if
     enddo
     enddo

     !** Calculate y-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
     do j = JJS-JFS_OFF, min(JJE,JEH)
     do i = IIS, IIE
       call FDM_MomY2VelYt( VelY,                                        & ! (out)
         & VL_ZXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
         & GSQRT, MAPF, GeneralCoordFlag )

       call EvalFlxY_CD6(FlxY_ZXV(:,i,j), VelY, Var_ZXY, 0, i, j, KS, KE)              

       if (present(DiffFlxY_ZXV)) then
         call FDM_AddDiffFlxY(FlxY_ZXV, DiffFlxY_ZXV, GSQRT(:,:,:,I_XVZ), MAPF(:,:,:,I_XV), i, j, KS, KE)
       end if
     enddo
     enddo
   
  end subroutine FDM_EvalFlxCD6_VarZXY
  
  subroutine FDM_EvalFlxCD6_VarZUY( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,             &  ! (inout)
    & Var_ZUY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZXY, DiffFlxY_ZUV, DiffFlxZ_WUY )                                 ! (in)   

    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZUY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZXY, DiffFlxY_ZUV, DiffFlxZ_WUY     

     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY
     
     !** Calculate z-direction flux ********************************************************

     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WUY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j,   KS,   KS)
       call EvalFlxZ_CD4(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j, KS+1, KS+1)
       call EvalFlxZ_CD6(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j, KS+2, KE-3)
       call EvalFlxZ_CD4(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j, KE-2, KE-2)
       call EvalFlxZ_CD2(FlxZ_WUY(:,i,j), VelZ, Var_ZUY, 0, i, j, KE-1, KE-1)
       FlxZ_WUY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WUY)) then
          call FDM_AddDiffFlxZ(FlxZ_WUY, DiffFlxZ_WUY, GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), 0, i, j, KS, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE+1
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD6(FlxX_ZXY(:,i-1,j), VelX, Var_ZUY, -1, i, j, KS, KE)              

      if (present(DiffFlxX_ZXY)) then
        call FDM_AddDiffFlxX(FlxX_ZXY, DiffFlxX_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), i, j, KS, KE)
      end if
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_ZUY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD6(FlxY_ZUV(:,i,j), VelY, Var_ZUY, 0, i, j, KS, KE)              

      if (present(DiffFlxY_ZUV)) then
        call FDM_AddDiffFlxY(FlxY_ZUV, DiffFlxY_ZUV, GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), i, j, KS, KE)
      end if
    enddo
    enddo

  end subroutine FDM_EvalFlxCD6_VarZUY
  
  subroutine FDM_EvalFlxCD6_VarZXV( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,          &  ! (inout)
    & Var_ZXV, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)   
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_ZUV, DiffFlxY_ZXY, DiffFlxZ_WXV )                                 ! (in)       
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_ZXV
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag     
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_ZUV, DiffFlxY_ZXY, DiffFlxZ_WXV
     
     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY
     
     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_WXV(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j,   KS,   KS)
       call EvalFlxZ_CD4(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j, KS+1, KS+1)
       call EvalFlxZ_CD6(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j, KS+2, KE-3)
       call EvalFlxZ_CD4(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j, KE-2, KE-2)       
       call EvalFlxZ_CD2(FlxZ_WXV(:,i,j), VelZ, Var_ZXV, 0, i, j, KE-1, KE-1)
       FlxZ_WXV(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_WXV)) then
          call FDM_AddDiffFlxZ(FlxZ_WXV, DiffFlxZ_WXV, GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), 0, i, j, KS, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD6(FlxX_ZUV(:,i,j), VelX, Var_ZXV, 0, i, j, KS, KE)              

      if (present(DiffFlxX_ZUV)) then
        call FDM_AddDiffFlxX(FlxX_ZUV, DiffFlxX_ZUV, GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), i, j, KS, KE)
      end if
    enddo
    enddo

    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE+1
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_ZXV, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD6(FlxY_ZXY(:,i,j-1), VelY, Var_ZXV, -1, i, j, KS, KE)              

      if (present(DiffFlxY_ZXY)) then
        call FDM_AddDiffFlxY(FlxY_ZXY, DiffFlxY_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), i, j, KS, KE)
      end if
    enddo
    enddo
           
  end subroutine FDM_EvalFlxCD6_VarZXV

  subroutine FDM_EvalFlxCD6_VarWXY( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,             &  ! (inout)
    & Var_WXY, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                        &  ! (in)
    & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                        &  ! (in)          
    & IIS, IIE, JJS, JJE, KS, KE,                                             &  ! (in)
    & DiffFlxX_WUY, DiffFlxY_WXV, DiffFlxZ_ZXY )                                 ! (in)   
    
     real(RP), intent(out), dimension(KA,IA,JA)  :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Var_WXY
     real(RP), intent(in),  dimension(KA,IA,JA)  :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
     real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
     real(RP), intent(in)                          :: J33G
     real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
     logical,  intent(in)                          :: GeneralCoordFlag
     integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE    
     real(RP), intent(in),  dimension(KA,IA,JA), optional :: DiffFlxX_WUY, DiffFlxY_WXV, DiffFlxZ_ZXY

     integer :: i, j, k
     real(RP), dimension(KA) :: VelZ, VelX, VelY

     !** Calculate z-direction flux ********************************************************
    
     !$omp parallel do private(i,j,k,VelZ) OMP_SCHEDULE_ collapse(2)
     do j = JJS, JJE
     do i = IIS, IIE
       call FDM_MomZ2VelZt( VelZ,                                             & ! (out)
            & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,   & ! (in)
            & J13G, J23G, J33G, MAPF, GeneralCoordFlag )                        ! (in)

       FlxZ_ZXY(KS-1,i,j) = 0.0_RP
       call EvalFlxZ_CD2(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KS+1, KS+1)
       call EvalFlxZ_CD4(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KS+2, KS+2)
       call EvalFlxZ_CD6(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KS+3, KE-3)
       call EvalFlxZ_CD4(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KE-2, KE-2)       
       call EvalFlxZ_CD2(FlxZ_ZXY(:,i,j), VelZ, Var_WXY, -1, i, j, KE-1, KE-1)
       FlxZ_ZXY(KE-1,i,j) = 0.0_RP
!       FlxZ_ZXY(  KE,i,j) = 0.0_RP

       if (present(DiffFlxZ_ZXY)) then
         call FDM_AddDiffFlxZ(FlxZ_ZXY, DiffFlxZ_ZXY, GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), -1, i, j, KS+1, KE-1)
       end if
    enddo
    enddo

    !** Calculate x-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelX) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS-1, IIE
      call FDM_MomX2VelXt( VelX,                                        & ! (out)
        & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxX_CD6(FlxX_WUY(:,i,j), VelX, Var_WXY, 0, i, j, KS, KE-1)              
      FlxX_WUY(KE,i,j) = 0.0_RP

      if (present(DiffFlxX_WUY)) then
        call FDM_AddDiffFlxX(FlxX_WUY, DiffFlxX_WUY, GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), i, j, KS, KE-1)
      end if              
    enddo
    enddo
    
    !** Calculate y-direction flux ********************************************************

    !$omp parallel do private(i,j,k,VelY) OMP_SCHEDULE_ collapse(2)
    do j = JJS-1, JJE
    do i = IIS, IIE
      call FDM_MomY2VelYt( VelY,                                        & ! (out)
        & VL_WXY, i, j, KS, KE, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, & ! (in)
        & GSQRT, MAPF, GeneralCoordFlag )

      call EvalFlxY_CD6(FlxY_WXV(:,i,j), VelY, Var_WXY, 0, i, j, KS, KE-1)              
      FlxY_WXV(KE,i,j) = 0.0_RP
      
      if (present(DiffFlxY_WXV)) then
        call FDM_AddDiffFlxY(FlxY_WXV, DiffFlxY_WXV, GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), i, j, KS, KE-1)
      end if      
    enddo
    enddo

  end subroutine FDM_EvalFlxCD6_VarWXY
  
end module scale_atmos_numeric_fdm_cd
