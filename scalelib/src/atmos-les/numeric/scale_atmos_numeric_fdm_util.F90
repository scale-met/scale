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
module scale_atmos_numeric_fdm_util
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
       VL_ZXY, VL_ZUY, VL_ZXV, VL_WXY
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  interface FDM_EvolveVar
     module procedure EvolveVar
     module procedure EvolveVar_withVart
  end interface FDM_EvolveVar
  
  public :: FDM_util_setup
  public :: FDM_EvolveVar

  public :: FDM_AddDiffFlxX
  public :: FDM_AddDiffFlxY
  public :: FDM_AddDiffFlxZ

  public :: FDM_RhoVar2Var
  public :: FDM_MomX2VelXt
  public :: FDM_MomY2VelYt
  public :: FDM_MomZ2VelZt
  
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

  integer, save :: IFS_OFF
  integer, save :: JFS_OFF
  real(RP), parameter :: FACT_N =   7.0_RP/12.0_RP
  real(RP), parameter :: FACT_F = - 1.0_RP/12.0_RP
  
contains
  subroutine fdm_util_setup(IFS_OFF_, JFS_OFF_)
    integer, intent(in) :: IFS_OFF_, JFS_OFF_

    IFS_OFF = IFS_OFF_; JFS_OFF = JFS_OFF_
    
  end subroutine fdm_util_setup

  subroutine EvolveVar( VarA,                                 & ! (out)
       & Var0, FlxX, FlxY, FlxZ,                              & ! (in)
       & GSQRT, MAPF, dt, RDX, RDY, RDZ,                      & ! (in)
       & IIS, IIE, JJS, JJE, KS_, KE_,                        & ! (in)
       & lhist, advch_t, advcv_t )
    
    real(RP), intent(out), dimension(KA,IA,JA)   :: VarA
    real(RP), intent(in), dimension(KA,IA,JA)    :: Var0
    real(RP), intent(in), dimension(KA,IA,JA)    :: FlxX, FlxY, FlxZ
    real(RP), intent(in), dimension(KA,IA,JA)    :: GSQRT
    real(RP), intent(in), dimension(IA,JA,2)     :: MAPF
    real(RP), intent(in) :: dt, RDX(:), RDY(:), RDZ(:)
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS_, KE_
    logical, intent(in), optional :: lhist
    real(RP), intent(inout), optional, dimension(KA,IA,JA) :: advch_t, advcv_t

    integer :: i, j, k
    real(RP) :: advch, advcv

    !$omp  parallel do private(i,j,k,advch,advcv) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
      do k = KS_, KE_
        advcv = - ( FlxZ(k,i,j) - FlxZ(k-1,i,  j)   ) * RDZ(k)
        advch = - ( FlxX(k,i,j) - FlxX(k  ,i-1,j)   ) * RDX(i) &
             &  - ( FlxY(k,i,j) - FlxY(k  ,i,  j-1) ) * RDY(j)
        VarA(k,i,j) = Var0(k,i,j) &
             & + dt * (  ( advcv + advch ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j) )
          
#ifdef HIST_TEND
        if ( lhist ) then
          advcv_t(k,i,j) = advcv * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
          advch_t(k,i,j) = advch * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
       end if
#endif
      enddo
    enddo
    enddo      

  end subroutine EvolveVar

  subroutine EvolveVar_withVart( VarA,                 & ! (out)
       & Var0, FlxX, FlxY, FlxZ,                       & ! (in)
       & GSQRT, MAPF, dt, RDX, RDY, RDZ,               & ! (in)
       & IIS, IIE, JJS, JJE, KS_, KE_,                 & ! (in)
       & Var_t,                                        & ! (in)
       & lhist, advch_t, advcv_t )
    
    real(RP), intent(out), dimension(KA,IA,JA)   :: VarA
    real(RP), intent(in), dimension(KA,IA,JA)    :: Var0
    real(RP), intent(in), dimension(KA,IA,JA)    :: FlxX, FlxY, FlxZ
    real(RP), intent(in), dimension(KA,IA,JA)    :: GSQRT
    real(RP), intent(in), dimension(IA,JA,2)     :: MAPF
    real(RP), intent(in) :: dt, RDX(:), RDY(:), RDZ(:)
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS_, KE_
    real(RP), intent(in), dimension(KA,IA,JA) :: Var_t
    logical, intent(in), optional :: lhist
    real(RP), intent(inout), optional, dimension(KA,IA,JA) :: advch_t, advcv_t

    integer :: i, j, k
    real(RP) :: advch, advcv

    !$omp  parallel do private(i,j,k,advch,advcv) OMP_SCHEDULE_ collapse(2)
    do j = JJS, JJE
    do i = IIS, IIE
      do k = KS_, KE_
        advcv = - ( FlxZ(k,i,j) - FlxZ(k-1,i,  j)   ) * RDZ(k)
        advch = - ( FlxX(k,i,j) - FlxX(k  ,i-1,j)   ) * RDX(i) &
             &  - ( FlxY(k,i,j) - FlxY(k  ,i,  j-1) ) * RDY(j)
        VarA(k,i,j) = Var0(k,i,j) &
             & + dt * (  ( advcv + advch ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j) &
             &          + Var_t(k,i,j) )
#ifdef HIST_TEND
        if ( lhist ) then
          advcv_t(k,i,j) = advcv * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
          advch_t(k,i,j) = advch * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
        end if
#endif
      enddo
    enddo
    enddo
 
  end subroutine EvolveVar_WithVart

  !-- common routines -----------------------------------------------------------
  !
  
  subroutine FDM_AddDiffFlxZ(FlxZ, DiffFlxZ, GSQRT, MAPF, KOF, i, j, KS_, KE_)
    real(RP), intent(inout) :: FlxZ(KA,IA,JA)
    real(RP), intent(in) :: DiffFlxZ(KA,IA,JA)
    real(RP), intent(in) :: GSQRT(KA,IA,JA), MAPF(IA,JA,2)
    integer,  intent(in) :: KOF, i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxZ(k+KOF,i,j) = FlxZ(k+KOF,i,j) + GSQRT(k,i,j) * DiffFlxZ(k,i,j)  &
            &                        / (MAPF(i,j,1) * MAPF(i,j,2))
    enddo
    
  end subroutine FDM_AddDiffFlxZ

  subroutine FDM_AddDiffFlxX(FlxX, DiffFlxX, GSQRT, MAPF, i, j, KS_, KE_)
    real(RP), intent(inout) :: FlxX(KA,IA,JA)
    real(RP), intent(in) :: DiffFlxX(KA,IA,JA)
    real(RP), intent(in) :: GSQRT(KA,IA,JA), MAPF(IA,JA,2)
    integer,  intent(in) :: i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxX(k,i,j) = FlxX(k,i,j) + GSQRT(k,i,j) * DiffFlxX(k,i,j) / MAPF(i,j,2)
    enddo
    
  end subroutine FDM_AddDiffFlxX

  subroutine FDM_AddDiffFlxY(FlxY, DiffFlxY, GSQRT, MAPF, i, j, KS_, KE_)
    real(RP), intent(inout) :: FlxY(KA,IA,JA)
    real(RP), intent(in) :: DiffFlxY(KA,IA,JA)
    real(RP), intent(in) :: GSQRT(KA,IA,JA), MAPF(IA,JA,2)
    integer,  intent(in) :: i, j, KS_, KE_

    integer :: k

    do k=KS_, KE_
       FlxY(k,i,j) = FlxY(k,i,j) + GSQRT(k,i,j) * DiffFlxY(k,i,j) / MAPF(i,j,1)
    enddo
    
  end subroutine FDM_AddDiffFlxY
  
  ! Calculation of velocity with momentum and density
  !
  
  subroutine FDM_RhoVar2VarCore(  Var,                     & ! (out)
       & RhoVar, DENS_ZXY,                                 & ! (in)
       & IOF, JOF, KOF, IIS_, IIE_, JJS_, JJE_, KS_, KE_ )   ! (in)
    real(RP), intent(out)   :: Var(KA,IA,JA)
    real(RP), intent(in)    :: RhoVar(KA,IA,JA)
    real(RP), intent(in)    :: DENS_ZXY(KA,IA,JA)
    integer,  intent(in)    :: IOF, JOF, KOF
    integer,  intent(in)    :: IIS_, IIE_, JJS_, JJE_, KS_, KE_

    integer :: k, i, j

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JJS_, JJE_
    do i = IIS_, IIE_
      do k = KS_, KE_
        Var(k,i,j) = 2.0_RP * RhoVar(k,i,j) / (DENS_ZXY(k,i,j) + DENS_ZXY(k+KOF,i+IOF,j+JOF))
      enddo     
    enddo
    enddo
  end subroutine FDM_RhoVar2VarCore

  subroutine FDM_RhoVar2Var(  Var,                                        & ! (out)
       & RhoVar, VarLoc, DENS_ZXY, IIS_, IIE_, JJS_, JJE_, KS_, KE_ )       ! (in)
    real(RP), intent(out)   :: Var(KA,IA,JA)
    real(RP), intent(in)    :: RhoVar(KA,IA,JA)
    integer,  intent(in)    :: VarLoc
    real(RP), intent(in)    :: DENS_ZXY(KA,IA,JA)
    integer,  intent(in)    :: IIS_, IIE_, JJS_, JJE_, KS_, KE_

    integer :: IOF, JOF, KOF
    select case(VarLoc)
    case(VL_ZXY)
       IOF = 0; JOF = 0; KOF = 0
    case(VL_ZUY)
       IOF = 1; JOF = 0; KOF = 0
    case(VL_ZXV)
       IOF = 0; JOF = 1; KOF = 0
    case(VL_WXY)
       IOF = 0; JOF = 0; KOF = 1
    end select
    call FDM_RhoVar2VarCore(  Var,                         & ! (out)
       & RhoVar, DENS_ZXY,                                 & ! (in)
       & IOF, JOF, JOF, IIS_, IIE_, JJS_, JJE_, KS_, KE_ )   ! (in)    
  end subroutine FDM_RhoVar2Var
  
  subroutine FDM_MomX2VelXt( VelX, &
       & VarLocID, i, j, KS, KE, &
       & Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, &
       & GSQRT, MAPF, GeneralCoordFlag )
    real(RP), intent(out), dimension(KA) :: VelX
    integer,  intent(in) :: VarLocID, i, j, KS, KE
    real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
    real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT
    real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
    logical,  intent(in)                          :: GeneralCoordFlag

    integer :: k

    select case (VarLocID)
    case(VL_ZXY) ! VelX_ZUY
      if( GeneralCoordFlag ) then
        do k=KS, KE
          VelX(k) =   2.0_RP * MomX_ZUY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i+1,j)) &
                &   * GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY)
        enddo
      else
        do k=KS, KE
          VelX(k) =   2.0_RP * MomX_ZUY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i+1,j))
        enddo
      end if
    case(VL_ZUY) ! VelX_ZXV
      if( GeneralCoordFlag ) then
        do k=KS, KE
          VelX(k) =   sum( MomX_ZUY(k,i-1:i,j) / (Dens_ZXY(k,i-1:i,j) + Dens_ZXY(k,i:i+1,j)) ) &
                &   * GSQRT(k,i,j,I_XVZ) / MAPF(i,j,2,I_XV)
        enddo
      else
        do k=KS, KE
          VelX(k) = sum( MomX_ZUY(k,i-1:i,j) / (Dens_ZXY(k,i-1:i,j) + Dens_ZXY(k,i:i+1,j)) )
        enddo
      end if
    case(VL_ZXV) ! VelX_ZUV
      if( GeneralCoordFlag ) then
        do k=KS, KE
          VelX(k) =   sum( MomX_ZUY(k,i,j-1:j) / (Dens_ZXY(k,i,j-1:j) + Dens_ZXY(k,i+1,j-1:j)) ) &
                &   * GSQRT(k,i,j,I_UVZ) / MAPF(i,j,2,I_UV)
        enddo
      else
        do k=KS, KE
          VelX(k) = sum( MomX_ZUY(k,i,j-1:j) / (Dens_ZXY(k,i,j-1:j) + Dens_ZXY(k,i+1,j-1:j)) )
        enddo
      end if
    case(VL_WXY) ! VelX_WUY
     if( GeneralCoordFlag ) then
        do k=KS, KE-1
          VelX(k) =   sum( MomX_ZUY(k:k+1,i,j) / (Dens_ZXY(k:k+1,i,j) + Dens_ZXY(k:k+1,i+1,j)) ) &
                &   * GSQRT(k,i,j,I_UYW) / MAPF(i,j,2,I_UY)
        enddo
      else
        do k=KS, KE-1
          VelX(k) = sum( MomX_ZUY(k:k+1,i,j) / (Dens_ZXY(k:k+1,i,j) + Dens_ZXY(k:k+1,i+1,j)) )
        enddo
      end if

    end select
    
  end subroutine FDM_MomX2VelXt

  subroutine FDM_MomY2VelYt( VelY,               &
       & VarLocID, i, j, KS, KE,                 &
       & Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, &
       & GSQRT, MAPF, GeneralCoordFlag )

    real(RP), intent(out), dimension(KA) :: VelY
    integer,  intent(in) :: VarLocID, i, j, KS, KE
    real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
    real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT
    real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
    logical,  intent(in)                          :: GeneralCoordFlag
    
    integer :: k

    select case (VarLocID)
    case(VL_ZXY) ! VelY_ZXV
      if( GeneralCoordFlag ) then
        do k=KS, KE
          VelY(k) =   2.0_RP * MomY_ZXV(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i,j+1)) &
                &   * GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV)
        enddo
      else
        do k=KS, KE
          VelY(k) =   2.0_RP * MomY_ZXV(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k,i,j+1))
        enddo
      end if
    case(VL_ZUY) ! VelY_ZUV
      if( GeneralCoordFlag ) then
        do k=KS, KE
          VelY(k) =   sum( MomY_ZXV(k,i:i+1,j) / (Dens_ZXY(k,i:i+1,j) + Dens_ZXY(k,i:i+1,j+1)) ) &
                &   * GSQRT(k,i,j,I_UVZ) / MAPF(i,j,1,I_UV)
        enddo
      else
        do k=KS, KE
          VelY(k) =   sum( MomY_ZXV(k,i:i+1,j) / (Dens_ZXY(k,i:i+1,j) + Dens_ZXY(k,i:i+1,j+1)) )
        enddo
      end if
    case(VL_ZXV) ! VelY_ZXY
      if( GeneralCoordFlag ) then
        do k=KS, KE
          VelY(k) =   sum( MomY_ZXV(k,i,j-1:j) / (Dens_ZXY(k,i,j-1:j) + Dens_ZXY(k,i,j:j+1)) ) &
                &   * GSQRT(k,i,j,I_XYZ) / MAPF(i,j,2,I_XY)
        enddo
      else
        do k=KS, KE
          VelY(k) = sum( MomY_ZXV(k,i,j-1:j) / (Dens_ZXY(k,i,j-1:j) + Dens_ZXY(k,i,j:j+1)) )
        enddo
      end if
    case(VL_WXY) ! VelX_WXV
     if( GeneralCoordFlag ) then
        do k=KS, KE-1
          VelY(k) =   sum( MomY_ZXV(k:k+1,i,j) / (Dens_ZXY(k:k+1,i,j) + Dens_ZXY(k:k+1,i,j+1)) ) &
                &   * GSQRT(k,i,j,I_XVW) / MAPF(i,j,2,I_XV)
        enddo
      else
        do k=KS, KE-1
          VelY(k) = sum( MomY_ZXV(k:k+1,i,j) / (Dens_ZXY(k:k+1,i,j) + Dens_ZXY(k:k+1,i,j+1)) )
        enddo
      end if
    end select
    
  end subroutine FDM_MomY2VelYt
  
  subroutine FDM_MomZ2VelZt( VelZ, &
       & VarLocID, i, j, KS, KE, &
       & Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY, &
       & J13G, J23G, J33G, MAPF, GeneralCoordFlag )
    
    real(RP), intent(out), dimension(KA) :: VelZ
    integer,  intent(in) :: VarLocID, i, j, KS, KE
    real(RP), intent(in),  dimension(KA,IA,JA)    :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
    real(RP), intent(in),  dimension(KA,IA,JA,7)  :: J13G, J23G
    real(RP), intent(in)                          :: J33G
    real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF
    logical,  intent(in)                          :: GeneralCoordFlag

    integer :: k
    
    select case (VarLocID)
    case(VL_ZXY) ! VelZ_WXY
      if ( GeneralCoordFlag ) then
        do k = KS, KE-1
          VelZ(k) = &
            &   J33G / (MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY))                                      & !
            &     * 2.0_RP * MomZ_WXY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k+1,i,j))              & !
            & + J13G(k,i,j,I_XYW) / MAPF(i,j,2,I_XY) * 0.5_RP * sum(                              & ! Half level velocity(x) and [ ZUY -> WXY ]
            &     MomX_ZUY(k:k+1,i-1:i,j) / (Dens_ZXY(k:k+1,i-1:i,j) + Dens_ZXY(k:k+1,i:i+1,j)) ) & ! 
            & + J23G(k,i,j,I_XYW) / MAPF(i,j,1,I_XY) * 0.5_RP * sum(                              & ! Half level velocity(y) and [ ZXV -> WXY ]
            &     MomY_ZXV(k:k+1,i,j-1:j) / (Dens_ZXY(k:k+1,i,j-1:j) + Dens_ZXY(k:k+1,i,j:j+1)) )   ! 
        enddo
      else    
        do k = KS, KE-1
          VelZ(k) = 2.0_RP * MomZ_WXY(k,i,j) / (Dens_ZXY(k,i,j) + Dens_ZXY(k+1,i,j))
        enddo
      end if
    case(VL_ZUY) ! VelZ_WUY
      if ( GeneralCoordFlag ) then
        do k = KS, KE-1
          VelZ(k) = & 
            &   J33G / (MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY))                                             & !
            &     * sum( MomZ_WXY(k,i:i+1,j) / (Dens_ZXY(k,i:i+1,j) + Dens_ZXY(k+1,i:i+1,j)) )           & !
            & + J13G(k,i,j,I_UYW) / MAPF(i,j,2,I_UY) * sum(                                              & ! Half level velocity(x) and [ ZUY -> WUY ]
            &     MomX_ZUY(k:k+1,i,j) / (Dens_ZXY(k:k+1,i,j) + Dens_ZXY(k:k+1,i+1,j)) )                  & ! 
            & + J23G(k,i,j,I_UYW) / MAPF(i,j,1,I_UY) * 0.25_RP * sum(                                    & ! Half level velocity(y) and [ ZXV -> WUY ]
            &     MomY_ZXV(k:k+1,i:i+1,j-1:j) / (Dens_ZXY(k:k+1,i:i+1,j-1:j) + Dens_ZXY(k:k+1,i:i+1,j:j+1)) )
        enddo           
      else
        do k = KS, KE-1
          VelZ(k) = sum( MomZ_WXY(k,i:i+1,j) / (Dens_ZXY(k,i:i+1,j) + Dens_ZXY(k+1,i:i+1,j)) )
        enddo
      end if
    case(VL_ZXV) ! VelZ_WXV
      if ( GeneralCoordFlag ) then
        do k = KS, KE-1
          VelZ(k) = &
            &   J33G / (MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV))                                                   & !
            &     * sum( MomZ_WXY(k,i,j:j+1) / (Dens_ZXY(k,i,j:j+1) + Dens_ZXY(k+1,i,j:j+1)) )                 & !
            & + J13G(k,i,j,I_XVW) / MAPF(i,j,2,I_XV) * 0.25_RP * sum(                                          & ! Half level velocity(y) and [ ZXV -> WUY ]
            &     MomX_ZUY(k:k+1,i-1:i,j:j+1) / (Dens_ZXY(k:k+1,i-1:i,j:j+1) + Dens_ZXY(k:k+1,i:i+1,j:j+1)) )  & ! 
            & + J23G(k,i,j,I_XVW) / MAPF(i,j,1,I_XV) * sum(                                                    & ! Half level velocity(x) and [ ZUY -> WUY ]
            &       MomY_ZXV(k:k+1,i,j) / (Dens_ZXY(k:k+1,i,j) + Dens_ZXY(k:k+1,i,j+1)) )                        ! 
        enddo 
      else
        !- Half-level velocity(z) and [WXY -> WXV]
        do k = KS, KE-1
          VelZ(k) = sum( MomZ_WXY(k,i,j:j+1) / (Dens_ZXY(k,i,j:j+1) + Dens_ZXY(k+1,i,j:j+1)) )
        enddo
      end if
    case(VL_WXY) ! VelZ_ZXY
      if ( GeneralCoordFlag ) then
        do k = KS+1, KE-1
          VelZ(k-1) = &
            &   J33G / (MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY))                              & !
            &     * 0.5_RP * sum(MomZ_WXY(k-1:k,i,j)) / Dens_ZXY(k,i,j)                   & !
            & + J13G(k,i,j,I_XYZ) / MAPF(i,j,2,I_XY)                                      & ! Half level velocity(y) and [ ZXV -> WUY ]
            &     * 0.5_RP * sum(MomX_ZUY(k,i-1:i,j)) / Dens_ZXY(k,i,j)                   & !
            & + J23G(k,i,j,I_XYZ) / MAPF(i,j,1,I_XY)                                      & ! Half level velocity(x) and [ ZUY -> WUY ]
            &     * 0.5_RP * sum(MomY_ZXV(k,i,j-1:j)) / Dens_ZXY(k,i,j)                     !
        enddo
      else
        !- Half-level velocity(z) and [WXY -> ZXY]
        do k = KS+1, KE-1
          VelZ(k-1) = 0.5_RP * sum(MomZ_WXY(k-1:k,i,j)) / Dens_ZXY(k,i,j)
        enddo
      end if      
    end select
    
  end subroutine FDM_MomZ2VelZt
  
end module scale_atmos_numeric_fdm_util
