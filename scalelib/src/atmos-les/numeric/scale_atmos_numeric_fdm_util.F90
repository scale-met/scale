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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  public :: fdm_util_setup
  public :: FDM_EvolveVar
  
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
  subroutine fdm_util_setup(IFS_OFF_, JFS_OFF_)
    integer, intent(in) :: IFS_OFF_, JFS_OFF_

    IFS_OFF = IFS_OFF_; JFS_OFF = JFS_OFF_
    
  end subroutine fdm_util_setup


  subroutine FDM_EvolveVar( VarA,                             & ! (out)
       & Var0, FlxX, FlxY, FlxZ,                              & ! (in)
       & GSQRT, MAPF, dt, RDX, RDY, RDZ,                      & ! (in)
       & IIS, IIE, JJS, JJE, KS_, KE_,                        & ! (in)
       & Var_t,                                               & ! (in)
       & lhist, advch_t, advcv_t )
    
    real(RP), intent(out), dimension(KA,IA,JA)   :: VarA
    real(RP), intent(in), dimension(KA,IA,JA)    :: Var0
    real(RP), intent(in), dimension(KA,IA,JA)    :: FlxX, FlxY, FlxZ
    real(RP), intent(in), dimension(KA,IA,JA)    :: GSQRT
    real(RP), intent(in), dimension(IA,JA,2)     :: MAPF
    real(RP), intent(in) :: dt, RDX(:), RDY(:), RDZ(:)
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS_, KE_
    real(RP), intent(in), dimension(KA,IA,JA), optional :: Var_t
    logical, intent(in), optional :: lhist
    real(RP), intent(inout), optional, dimension(KA,IA,JA) :: advch_t, advcv_t

    integer :: i, j, k
    real(RP) :: advch, advcv

    if (present(Var_t)) then
      !$omp  parallel do private(i,j,k,advch,advcv) OMP_SCHEDULE_ collapse(2)
      do j = JJS, JJE
      do i = IIS, IIE
        do k = KS_, KE_
          advcv = - ( FlxZ(k,i,j) - FlxZ(k-1,i,  j)   ) * RDZ(k)
          advch = - ( FlxX(k,i,j) - FlxX(k  ,i-1,j)   ) * RDX(i) &
               &  - ( FlxY(k,i,j) - FlxY(k  ,i,  j-1) ) * RDY(j)
          VarA(k,i,j) = Var0(k,i,j) &
               & + dt * (  ( advcv + advch ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j) &
                         + Var_t(k,i,j) )
#ifdef HIST_TEND
          if ( lhist ) then
            advcv_t(k,i,j) = advcv * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
            advch_t(k,i,j) = advch * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
          end if
#endif
        enddo
      enddo
      enddo
    else
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
    endif  
  end subroutine FDM_EvolveVar

!!$  subroutine ATMOS_DYN_Flx4VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,        &  ! (inout)
!!$       & Var_ZXY, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$       & FlxEvalType, IIS, IIE, JJS, JJE, KS, KE )                        ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G
!!$    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
!!$    integer, intent(in) :: FlxEvalType
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$    select case(FlxEvalType)
!!$    case(FLXEVALTYPE_CD4)
!!$       call ATMOS_DYN_Flx4VarZXY_CD4( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,          &  ! (inout)
!!$            & Var_ZXY, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$            & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$            & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)       
!!$    case(FLXEVALTYPE_UD1)
!!$       call ATMOS_DYN_Flx4VarZXY_UD1( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,          &  ! (inout)
!!$            & Var_ZXY, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$            & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$            & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)
!!$     case default
!!$     end select
!!$     
!!$  end subroutine ATMOS_DYN_Flx4VarZXY
!!$
!!$  subroutine ATMOS_DYN_MomFlx4VarZXY( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,    &  ! (inout)
!!$       & MOMX_ZUY, MOMY_ZXV, MOMZ_WXY,                                 &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                      &  ! (in)
!!$       & FlxEvalType, IIS, IIE, JJS, JJE, KS, KE )                        ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MOMX_ZUY, MOMY_ZXV, MOMZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G    
!!$    real(RP), intent(in), dimension(KA,IA,JA,3)   :: num_diff
!!$    integer, intent(in) :: FlxEvalType
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$    select case(FlxEvalType)
!!$    case(FLXEVALTYPE_CD4)
!!$       call ATMOS_DYN_MomFlx4VarZXY_CD4( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY, &  ! (inout)
!!$       & MOMX_ZUY, MOMY_ZXV, MOMZ_WXY,                                 &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                      &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                     ! (in)       
!!$    case(FLXEVALTYPE_UD1)
!!$       call ATMOS_DYN_MomFlx4VarZXY_UD1( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,  &  ! (inout)
!!$            & MOMX_ZUY, MOMY_ZXV, MOMZ_WXY,                             &  ! (in)
!!$            & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                  &  ! (in)
!!$            & IIS, IIE, JJS, JJE, KS, KE )                                 ! (in)       
!!$     case default
!!$     end select
!!$    
!!$   end subroutine ATMOS_DYN_MomFlx4VarZXY
!!$  
!!$  !+ ------------------------------------------------------------------------------------------
!!$  
!!$  subroutine ATMOS_DYN_Flx4VarZXY_CD4( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,    &  ! (inout)
!!$       & Var_ZXY, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G
!!$    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$    integer :: i, j, k
!!$
!!$    !* Z Dir -----------------------------------------------------
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS,   JJE
!!$       do i = IIS,   IIE
!!$          do k = KS+1, KE-2
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, FlxZ_WXY(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k-1,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k  ,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k+1,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k+2,i,j) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,ZDIR) )
!!$#endif
!!$             FlxZ_WXY(k,i,j) = & 
!!$                  &   MomFlxZ_WXY(k,i,j)*(  FACT_N * ( Var_ZXY(k+1,i,j) + Var_ZXY(k  ,i,j) )   &
!!$                  &                       + FACT_F * ( Var_ZXY(k+2,i,j) + Var_ZXY(k-1,i,j) ) ) & ! [{x,y,z->x,y,w}]
!!$                  & + GSQRT(k,i,j,I_XYW) * num_diff(k,i,j,ZDIR)
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$       k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$       do j = JJS, JJE
!!$          do i = IIS, IIE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, FlxZ_WXY(KS,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(KS+1,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(KS  ,i,j) )
!!$             call CHECK( __LINE__, num_diff(KS,i,j,ZDIR) )
!!$             call CHECK( __LINE__, FlxZ_WXY(KE-1,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(KE-1,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(KE  ,i,j) )
!!$             call CHECK( __LINE__, num_diff(KE-1,i,j,I_RHOT,ZDIR) )
!!$#endif
!!$             FlxZ_WXY(KS-1,i,j) = 0.0_RP
!!$             FlxZ_WXY(KS  ,i,j) =   MomFlxZ_WXY(KS  ,i,j) * 0.5_RP * ( Var_ZXY(KS+1,i,j) + Var_ZXY(KS,i,j)   ) &
!!$                  &               + GSQRT(KS  ,i,j,I_XYW) * num_diff(KS  ,i,j,ZDIR)
!!$             FlxZ_WXY(KE-1,i,j) =   MomFlxZ_WXY(KE-1,i,j) * 0.5_RP * ( Var_ZXY(KE,i,j)   + Var_ZXY(KE-1,i,j) ) &
!!$                  &               + GSQRT(KE-1,i,j,I_XYW) * num_diff(KE-1,i,j,ZDIR)
!!$             FlxZ_WXY(KE  ,i,j) = 0.0_RP
!!$          enddo
!!$       enddo
!!$#ifdef DEBUG
!!$       k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$    
!!$    !* X Dir -----------------------------------------------------
!!$    
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS,   JJE
!!$       do i = IIS-1, IIE
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxX_ZUY(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i-1,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i  ,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i+1,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i+1,j) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,XDIR) )
!!$#endif
!!$             FlxX_ZUY(k,i,j) = MomFlxX_ZUY(k,i,j) &
!!$                  *(   FACT_N * ( Var_ZXY(k,i+1,j) + Var_ZXY(k,i  ,j) )   &
!!$                     + FACT_F * ( Var_ZXY(k,i+2,j) + Var_ZXY(k,i-1,j) ) ) & ! [{x,y,z->u,y,z}]
!!$                  + GSQRT(k,i,j,I_UYZ) * num_diff(k,i,j,XDIR)
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    !* Y Dir -----------------------------------------------------
!!$    
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-1, JJE
!!$       do i = IIS,   IIE
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k,i,j,YDIR) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i,j-1) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i,j  ) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i,j+1) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i,j+2) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,YDIR) )
!!$#endif
!!$             FlxY_ZXV(k,i,j) = MomFlxY_ZXV(k,i,j) &
!!$                              * ( FACT_N * ( Var_ZXY(k,i,j+1) + Var_ZXY(k,i,j  ) )   &
!!$                                + FACT_F * ( Var_ZXY(k,i,j+2) + Var_ZXY(k,i,j-1) ) ) & ! [{x,y,z->x,v,z}]
!!$                              + GSQRT(k,i,j,I_XVZ) * num_diff(k,i,j,YDIR)
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$  end subroutine ATMOS_DYN_Flx4VarZXY_CD4  
!!$
!!$  subroutine ATMOS_DYN_Flx4VarZXY_UD1( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,    &  ! (inout)
!!$       & Var_ZXY, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G
!!$    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$    integer :: i, j, k
!!$
!!$    !* Z Dir -----------------------------------------------------
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-1, IIE+1
!!$          do k = KS, KE-1
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxZ_WXY(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k  ,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k+1,i,j) )
!!$#endif
!!$             FlxZ_WXY(k,i,j) = 0.5_RP * (      MomFlxZ_WXY(k,i,j)  * ( Var_ZXY(k+1,i,j) + Var_ZXY(k,i,j)  ) &
!!$                                         - abs(MomFlxZ_WXY(k,i,j)) * ( Var_ZXY(k+1,i,j) - Var_ZXY(k,i,j) ) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-1, IIE+1
!!$          FlxZ_WXY(KS-1,i,j) = 0.0_RP
!!$          FlxZ_WXY(KE  ,i,j) = 0.0_RP
!!$       enddo
!!$    enddo
!!$
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    
!!$    !* X Dir -----------------------------------------------------
!!$          
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-2, IIE+1
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxX_ZUY(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i  ,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i+1,j) )
!!$#endif
!!$             FlxX_ZUY(k,i,j) = 0.5_RP * (       MomFlxX_ZUY(k,i,j)  * ( Var_ZXY(k,i+1,j) + Var_ZXY(k,i,j) ) &
!!$                                          - abs(MomFlxX_ZUY(k,i,j)) * ( Var_ZXY(k,i+1,j) - Var_ZXY(k,i,j) ) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    !* Y Dir -----------------------------------------------------
!!$    
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-2, JJE+1
!!$       do i = IIS-1,   IIE+1
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k,i,j,YDIR) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i,j  ) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i,j+1) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,YDIR) )
!!$#endif
!!$             FlxY_ZXV(k,i,j) = 0.5_RP * (       MomFlxY_ZXV(k,i,j)  * ( Var_ZXY(k,i,j+1) + Var_ZXY(k,i,j) ) &
!!$                                          - abs(MomFlxY_ZXV(k,i,j)) * ( Var_ZXY(k,i,j+1) - Var_ZXY(k,i,j) ) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$  end subroutine ATMOS_DYN_Flx4VarZXY_UD1
!!$
!!$  subroutine ATMOS_DYN_Flx4VarZUY_UD1( FlxX_ZXY, FlxY_ZUV, FlxZ_WUY,    &  ! (inout)
!!$       & Var_ZUY, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZXY, FlxY_ZUV, FlxZ_WUY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZUY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G
!!$    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$    integer :: i, j, k
!!$
!!$    !* Z Dir -----------------------------------------------------
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-1, IIE+1
!!$          do k = KS, KE-1
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxZ_WUY(k,i,j) )
!!$             call CHECK( __LINE__, MomFlxZ_WUY(k,i+1,j) )
!!$             call CHECK( __LINE__, Var_ZUY(k  ,i,j) )
!!$             call CHECK( __LINE__, Var_ZUY(k+1,i,j) )
!!$#endif
!!$             FlxZ_WUY(k,i,j) = 0.25_RP *(      (MomFlxZ_WXY(k,i+1,j) + MomFlxZ_WXY(k,i,j)) * ( Var_ZUY(k+1,i,j) + Var_ZUY(k,i,j) )   &
!!$                  &                       - abs(MomFlxZ_WXY(k,i+1,j) + MomFlxZ_WXY(k,i,j)) * ( Var_ZUY(k+1,i,j) - Var_ZUY(k,i,j) ) ) &
!!$                  & / (MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY))
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-1, IIE+1
!!$          FlxZ_WUY(KS-1,i,j) = 0.0_RP
!!$          FlxZ_WUY(KE  ,i,j) = 0.0_RP
!!$       enddo
!!$    enddo
!!$
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    
!!$    !* X Dir -----------------------------------------------------
!!$    ! note that x-index is added by -1
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-1, IIE+2
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxX_ZUY(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZUY(k,i  ,j) )
!!$             call CHECK( __LINE__, Var_ZUY(k,i-1,j) )
!!$#endif
!!$             FlxX_ZXY(k,i,j) = 0.25_RP * (      (MomFlxX_ZUY(k,i,j) + MomFlxX_ZUY(k,i-1,j)) * ( Var_ZUY(k,i,j) + Var_ZUY(k,i-1,j) ) &
!!$                                           - abs(MomFlxX_ZUY(k,i,j) + MomFlxX_ZUY(k,i-1,j)) * ( Var_ZUY(k,i,j) - Var_ZUY(k,i-1,j) ) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    !* Y Dir -----------------------------------------------------
!!$    
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-2, JJE+1
!!$       do i = IIS-1,   IIE+1
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k,i+1,j) )
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZUY(k,i,j  ) )
!!$             call CHECK( __LINE__, Var_ZUY(k,i,j+1) )
!!$#endif
!!$             FlxY_ZUV(k,i,j) = 0.25_RP * (      (MomFlxY_ZXV(k,i+1,j) + MomFlxY_ZXV(k,i,j)) * ( Var_ZUY(k,i,j+1) + Var_ZUY(k,i,j) ) &
!!$                                           - abs(MomFlxY_ZXV(k,i+1,j) + MomFlxY_ZXV(k,i,j)) * ( Var_ZUY(k,i,j+1) - Var_ZUY(k,i,j) ) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$  end subroutine ATMOS_DYN_Flx4VarZUY_UD1
!!$
!!$  subroutine ATMOS_DYN_Flx4VarZXV_UD1( FlxX_ZUV, FlxY_ZXY, FlxZ_WXV,    &  ! (inout)
!!$       & Var_ZXV, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUV, FlxY_ZXY, FlxZ_WXV
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_ZXV
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G
!!$    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$    integer :: i, j, k
!!$
!!$    !* Z Dir -----------------------------------------------------
!!$    
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)    
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-1, IIE+1
!!$          do k = KS, KE-1
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxZ_WXV(k,i,j) )
!!$             call CHECK( __LINE__, MomFlxZ_WXV(k,i,j+1) )
!!$             call CHECK( __LINE__, Var_ZXV(k  ,i,j) )
!!$             call CHECK( __LINE__, Var_ZXV(k+1,i,j) )
!!$#endif
!!$             FlxZ_WXV(k,i,j) = 0.25_RP *(      (MomFlxZ_WXY(k,i,j+1) + MomFlxZ_WXY(k,i,j)) * ( Var_ZXV(k+1,i,j) + Var_ZXV(k,i,j) )   &
!!$                  &                       - abs(MomFlxZ_WXY(k,i,j+1) + MomFlxZ_WXY(k,i,j)) * ( Var_ZXV(k+1,i,j) - Var_ZXV(k,i,j) ) ) &
!!$                  & / (MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV))
!!$          enddo
!!$          FlxZ_WXV(KS-1,i,j) = 0.0_RP
!!$          FlxZ_WXV(KE  ,i,j) = 0.0_RP          
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    
!!$    !* X Dir -----------------------------------------------------
!!$    ! note that x-index is added by -1
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-2, IIE+1
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxX_ZUY(k,i,j) )
!!$             call CHECK( __LINE__, MomFlxX_ZUY(k,i,j+1) )
!!$             call CHECK( __LINE__, Var_ZXV(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZXV(k,i+1,j) )
!!$#endif
!!$             FlxX_ZUV(k,i,j) = 0.25_RP * (      (MomFlxX_ZUY(k,i,j+1) + MomFlxX_ZUY(k,i,j)) * ( Var_ZXV(k,i+1,j) + Var_ZXV(k,i,j) ) &
!!$                                           - abs(MomFlxX_ZUY(k,i,j+1) + MomFlxX_ZUY(k,i,j)) * ( Var_ZXV(k,i+1,j) - Var_ZXV(k,i,j) ) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    !* Y Dir -----------------------------------------------------
!!$    
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-2, JJE+1
!!$       do i = IIS-1,   IIE+1
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k,i,j-1) )
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZXV(k,i,j  ) )
!!$             call CHECK( __LINE__, Var_ZXV(k,i,j-1) )
!!$#endif
!!$             FlxY_ZXY(k,i,j) = 0.25_RP * (      (MomFlxY_ZXV(k,i,j) + MomFlxY_ZXV(k,i,j-1)) * ( Var_ZXV(k,i,j) + Var_ZXV(k,i,j-1) ) &
!!$                                           - abs(MomFlxY_ZXV(k,i,j) + MomFlxY_ZXV(k,i,j-1)) * ( Var_ZXV(k,i,j) - Var_ZXV(k,i,j-1) ) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$  end subroutine ATMOS_DYN_Flx4VarZXV_UD1
!!$
!!$  subroutine ATMOS_DYN_Flx4VarWXY_UD1( FlxX_WUY, FlxY_WXV, FlxZ_ZXY,    &  ! (inout)
!!$       & Var_WXY, MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY,                &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                       &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                      ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_WUY, FlxY_WXV, FlxZ_ZXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: Var_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MomFlxX_ZUY, MomFlxY_ZXV, MomFlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G
!!$    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$    integer :: i, j, k
!!$
!!$    !* Z Dir -----------------------------------------------------
!!$    ! note that z-index is added by -1
!!$    
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)    
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-1, IIE+1
!!$          do k = KS+1, KE-1
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxZ_WXY(k-1,i,j) )
!!$             call CHECK( __LINE__, MomFlxZ_WXY(k,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k-1  ,i,j) )
!!$             call CHECK( __LINE__, Var_ZXY(k,i,j) )
!!$#endif
!!$             FlxZ_ZXY(k,i,j) = 0.25_RP * (    (MomFlxZ_WXY(k,i,j) + MomFlxZ_WXY(k-1,i,j)) * ( Var_WXY(k,i,j) + Var_WXY(k,i,j)   )   &
!!$                  &                      - abs(MomFlxZ_WXY(k,i,j) + MomFlxZ_WXY(k-1,i,j)) * ( Var_WXY(k,i,j) - Var_WXY(k-1,i,j) ) ) &
!!$                  &            / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
!!$                  
!!$          enddo
!!$          FlxZ_ZXY(KS-1,i,j) = 0.0_RP
!!$          FlxZ_ZXY(KE-1,i,j) = 0.0_RP          
!!$          FlxZ_ZXY(KE  ,i,j) = 0.0_RP          
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    
!!$    !* X Dir -----------------------------------------------------
!!$          
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-1, JJE+1
!!$       do i = IIS-2, IIE+1
!!$          do k = KS, KE-1
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxX_ZUY(k,  i,j) )
!!$             call CHECK( __LINE__, MomFlxX_ZUY(k+1,i,j) )
!!$             call CHECK( __LINE__, Var_WXY(k,i,  j) )
!!$             call CHECK( __LINE__, Var_WXY(k,i+1,j) )
!!$#endif
!!$             FlxX_WUY(k,i,j) = 0.25_RP * (    (MomFlxX_ZUY(k+1,i,j) + MomFlxX_ZUY(k,i,j)) * ( Var_WXY(k,i+1,j) + Var_WXY(k,i,j)   )   &
!!$                  &                      - abs(MomFlxX_ZUY(k+1,i,j) + MomFlxX_ZUY(k,i,j)) * ( Var_WXY(k,i+1,j) - Var_WXY(k,i,j) ) )
!!$          end do
!!$          FlxX_WUY(KE,i,j) = 0.0_RP          
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    !* Y Dir -----------------------------------------------------
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS-2, JJE+1
!!$       do i = IIS-1, IIE+1
!!$          do k = KS, KE-1
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k,  i,j) )
!!$             call CHECK( __LINE__, MomFlxY_ZXV(k+1,i,j) )
!!$             call CHECK( __LINE__, Var_WXY(k,i,  j) )
!!$             call CHECK( __LINE__, Var_WXY(k,i,j+1) )
!!$#endif
!!$             FlxY_WXV(k,i,j) = 0.25_RP * (    (MomFlxY_ZXV(k+1,i,j) + MomFlxY_ZXV(k,i,j)) * ( Var_WXY(k,i,j+1) + Var_WXY(k,i,j)   )   &
!!$                  &                      - abs(MomFlxY_ZXV(k+1,i,j) + MomFlxY_ZXV(k,i,j)) * ( Var_WXY(k,i,j+1) - Var_WXY(k,i,j) ) )
!!$          end do
!!$          FlxY_WXV(KE,i,j) = 0.0_RP
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$  end subroutine ATMOS_DYN_Flx4VarWXY_UD1
!!$  
!!$  subroutine ATMOS_DYN_MomFlx4VarZXY_CD4( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,    &  ! (inout)
!!$       & MOMX_ZUY, MOMY_ZXV, MOMZ_WXY,                                 &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                      &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                     ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MOMX_ZUY, MOMY_ZXV, MOMZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G    
!!$    real(RP), intent(in), dimension(KA,IA,JA,3)   :: num_diff
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$
!!$    integer :: i, j, k
!!$    real(DP) :: MOMX_WXY, MOMY_WXY
!!$
!!$    !** Flux of Z direction at (x, y, w)
!!$    
!!$    !$omp  parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    !$omp& private(MOMX_WXY, MOMY_WXY)
!!$    do j=JJS, JJE
!!$       do i=IIS, IIE
!!$          do k=KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MOMZ_WXY(k+1,i,j) )
!!$             call CHECK( __LINE__, MOMZ_WXY(k  ,i,j) )
!!$             call CHECK( __LINE__, MOMZ_WXY(k-1,i,j) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,ZDIR) )
!!$#endif             
!!$             MOMX_WXY = 0.25_RP*(  MOMX_ZUY(k+1,i,j) + MOMX_ZUY(k+1,i-1,j)     &
!!$                  &              + MOMX_ZUY(k  ,i,j) + MOMX_ZUY(k  ,i-1,j) )   ! [{u,y,z->x,y,w}]
!!$             MOMY_WXY = 0.25_RP*(  MOMY_ZXV(k+1,i,j) + MOMY_ZXV(k+1,i,j-1)     &
!!$                                 + MOMY_ZXV(k  ,i,j) + MOMY_ZXV(k  ,i,j-1) )   ! [{x,v,z->x,y,w}]
!!$             FlxZ_WXY(k,i,j) = &
!!$                  &   J33G * MOMZ_WXY(k,i,j) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) ) &
!!$                  & + J13G(k,i,j,I_XYW) * MOMX_WXY / MAPF(i,j,2,I_XY)                & 
!!$                  & + J23G(k,i,j,I_XYW) * MOMY_WXY / MAPF(i,j,1,I_XY)                & 
!!$                  & + GSQRT(k,i,j,I_XYW) * num_diff(k,i,j,ZDIR)                      &
!!$                  &   / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2, I_XY) )
!!$          end do
!!$       end do
!!$    end do
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    !** Flux of X direction at (x, y, w)
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS  , JJE
!!$       do i = IIS-IFS_OFF, min(IIE,IEH)
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MOMX_ZUY(k,i+1,j) )
!!$             call CHECK( __LINE__, MOMX_ZUY(k,i  ,j) )
!!$             call CHECK( __LINE__, MOMX_ZUY(k,i-1,j) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,XDIR) )
!!$#endif
!!$             FlxX_ZUY(k,i,j) =  GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY) &
!!$                  &           * ( MOMX_ZUY(k,i,j) + num_diff(k,i,j,XDIR) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$    !$omp  parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j=JJS-JFS_OFF, min(JJE,JEH)
!!$       do i=IIS, IIE
!!$          do k=KS, KE
!!$#ifdef DEBUG
!!$          call CHECK( __LINE__, MOMY_ZXV(k,i,j+1) )
!!$          call CHECK( __LINE__, MOMY_ZXV(k,i,j  ) )
!!$          call CHECK( __LINE__, MOMY_ZXV(k,i,j-1) )
!!$          call CHECK( __LINE__, num_diff(k,i,j,YDIR) )
!!$#endif             
!!$             FlxY_ZXV(k,i,j) = GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV) &
!!$                  &            * ( MOMY_ZXV(k,i,j) + num_diff(k,i,j,YDIR) )
!!$          end do
!!$       end do
!!$    end do
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$  end subroutine ATMOS_DYN_MomFlx4VarZXY_CD4
!!$
!!$  subroutine ATMOS_DYN_MomFlx4VarZXY_UD1( FlxX_ZUY, FlxY_ZXV, FlxZ_WXY,    &  ! (inout)
!!$       & MOMX_ZUY, MOMY_ZXV, MOMZ_WXY,                                 &  ! (in)
!!$       & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                      &  ! (in)
!!$       & IIS, IIE, JJS, JJE, KS, KE )                                     ! (in)
!!$
!!$    real(RP), intent(inout), dimension(KA,IA,JA)  :: FlxX_ZUY, FlxY_ZXV, FlxZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA)     :: MOMX_ZUY, MOMY_ZXV, MOMZ_WXY
!!$    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
!!$    real(RP), intent(in) :: J33G    
!!$    real(RP), intent(in), dimension(KA,IA,JA,3)   :: num_diff
!!$    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
!!$
!!$
!!$    integer :: i, j, k
!!$    real(DP) :: MOMX_WXY, MOMY_WXY
!!$
!!$    !** Flux of Z direction at (x, y, w)
!!$    
!!$    !$omp  parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    !$omp& private(MOMX_WXY, MOMY_WXY)
!!$    do j=JJS, JJE
!!$       do i=IIS, IIE
!!$          do k=KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MOMZ_WXY(k+1,i,j) )
!!$             call CHECK( __LINE__, MOMZ_WXY(k  ,i,j) )
!!$             call CHECK( __LINE__, MOMZ_WXY(k-1,i,j) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,ZDIR) )
!!$#endif             
!!$             MOMX_WXY = 0.25_RP*(  MOMX_ZUY(k+1,i,j) + MOMX_ZUY(k+1,i-1,j)     &
!!$                  &              + MOMX_ZUY(k  ,i,j) + MOMX_ZUY(k  ,i-1,j) )   ! [{u,y,z->x,y,w}]
!!$             MOMY_WXY = 0.25_RP*(  MOMY_ZXV(k+1,i,j) + MOMY_ZXV(k+1,i,j-1)     &
!!$                                 + MOMY_ZXV(k  ,i,j) + MOMY_ZXV(k  ,i,j-1) )   ! [{x,v,z->x,y,w}]
!!$             FlxZ_WXY(k,i,j) = &
!!$                  &   J33G * MOMZ_WXY(k,i,j) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) ) &
!!$                  & + J13G(k,i,j,I_XYW) * MOMX_WXY / MAPF(i,j,2,I_XY)                & 
!!$                  & + J23G(k,i,j,I_XYW) * MOMY_WXY / MAPF(i,j,1,I_XY)                & 
!!$                  & + GSQRT(k,i,j,I_XYW) * num_diff(k,i,j,ZDIR)                      &
!!$                  &   / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2, I_XY) )
!!$          end do
!!$       end do
!!$    end do
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$
!!$    !** Flux of X direction at (x, y, w)
!!$
!!$    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j = JJS  , JJE
!!$       do i = IIS-IFS_OFF, min(IIE,IEH)
!!$          do k = KS, KE
!!$#ifdef DEBUG
!!$             call CHECK( __LINE__, MOMX_ZUY(k,i+1,j) )
!!$             call CHECK( __LINE__, MOMX_ZUY(k,i  ,j) )
!!$             call CHECK( __LINE__, MOMX_ZUY(k,i-1,j) )
!!$             call CHECK( __LINE__, num_diff(k,i,j,XDIR) )
!!$#endif
!!$             FlxX_ZUY(k,i,j) =  GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY) &
!!$                  &           * ( MOMX_ZUY(k,i,j) + num_diff(k,i,j,XDIR) )
!!$          enddo
!!$       enddo
!!$    enddo
!!$#ifdef DEBUG
!!$    k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$    !$omp  parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!!$    do j=JJS-JFS_OFF, min(JJE,JEH)
!!$       do i=IIS, IIE
!!$          do k=KS, KE
!!$#ifdef DEBUG
!!$          call CHECK( __LINE__, MOMY_ZXV(k,i,j+1) )
!!$          call CHECK( __LINE__, MOMY_ZXV(k,i,j  ) )
!!$          call CHECK( __LINE__, MOMY_ZXV(k,i,j-1) )
!!$          call CHECK( __LINE__, num_diff(k,i,j,YDIR) )
!!$#endif             
!!$             FlxY_ZXV(k,i,j) = GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV) &
!!$                  &            * ( MOMY_ZXV(k,i,j) + num_diff(k,i,j,YDIR) )
!!$          end do
!!$       end do
!!$    end do
!!$#ifdef DEBUG
!!$       k = IUNDEF; i = IUNDEF; j = IUNDEF
!!$#endif
!!$       
!!$     end subroutine ATMOS_DYN_MomFlx4VarZXY_UD1
!!$     
end module scale_atmos_numeric_fdm_util
