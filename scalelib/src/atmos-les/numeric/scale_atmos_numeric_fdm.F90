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
module scale_atmos_numeric_fdm
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

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  public :: ATMOS_NUMERIC_FDM_setup
  public :: ATMOS_NUMERIC_FDM_evalFlux

  !* Cascade
  public :: ATMOS_NUMERIC_FDM_evolveVar
  
  !**
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  
  integer, parameter, public :: FLXEVALTYPE_CD2 = 101  
  integer, parameter, public :: FLXEVALTYPE_CD4 = 102
  integer, parameter, public :: FLXEVALTYPE_CD6 = 103

  integer, parameter, public :: FLXEVALTYPE_UD1 = 201  
  integer, parameter, public :: FLXEVALTYPE_UD3 = 202
  integer, parameter, public :: FLXEVALTYPE_UD5 = 203

  integer, parameter, public :: VL_XY  = 000
  integer, parameter, public :: VL_UY  = 010
  integer, parameter, public :: VL_XV  = 001  
  integer, parameter, public :: VL_ZXY = 100
  integer, parameter, public :: VL_ZUY = 110
  integer, parameter, public :: VL_ZXV = 101
  integer, parameter, public :: VL_WXY = 200

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

  abstract interface
     subroutine FDM_EvalFlux(FluxX, FluxY, FluxZ,                          &  ! (inout)
          & Var, VelX_ZUY, VelY_ZXV, VelZ_WXY,                              &  ! (in)
          & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                        &  ! (in)
          & IIS, IIE, JJS, JJE, KS_, KE_ )                                      ! (in)
       use scale_precision,  only: RP
       use scale_grid_index, only: IA, JA, KA
       real(RP), intent(inout), dimension(KA,IA,JA)  :: FluxX, FluxY, FluxZ
       real(RP), intent(in), dimension(KA,IA,JA)     :: Var
       real(RP), intent(in), dimension(KA,IA,JA)     :: VelX_ZUY, VelY_ZXV, VelZ_WXY
       real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
       real(RP), intent(in) :: J33G
       real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
       integer, intent(in) :: IIS, IIE, JJS, JJE, KS_, KE_
     end subroutine FDM_EvalFlux
  end interface
  
contains
  subroutine ATMOS_NUMERIC_FDM_setup(IFS_OFF_, JFS_OFF_)
    use scale_atmos_numeric_fdm_util, only: FDM_util_setup
    use scale_atmos_numeric_fdm_cd, only: FDM_CD_setup
    use scale_atmos_numeric_fdm_ud, only: FDM_UD_setup
    
    integer, intent(in) :: IFS_OFF_, JFS_OFF_
    
    IFS_OFF = IFS_OFF_; JFS_OFF = JFS_OFF_

    call FDM_util_setup(IFS_OFF, JFS_OFF)    
    call FDM_CD_setup(IFS_OFF, JFS_OFF)
    call FDM_UD_setup(IFS_OFF, JFS_OFF)
    
  end subroutine ATMOS_NUMERIC_FDM_setup

  subroutine ATMOS_NUMERIC_FDM_EvolveVar( VarA,                       & ! (out)
       & Var0, VarLoc, Flx,                                           & ! (in)
       & GSQRT, MAPF, dt, RDX, RDY, RDZ,                      & ! (in)
       & IIS, IIE, JJS, JJE, KS, KE,                          & ! (in)
       & Var_t,                                               & ! (in)
       & lhist, advch_t, advcv_t )

    use scale_atmos_numeric_fdm_util, only: FDM_EvolveVar
    
    real(RP), intent(out), dimension(KA,IA,JA)   :: VarA
    real(RP), intent(in), dimension(KA,IA,JA)    :: Var0
    integer, intent(in) :: VarLoc
    real(RP), intent(in), dimension(KA,IA,JA,3)  :: Flx
    real(RP), intent(in), dimension(KA,IA,JA,3)  :: GSQRT
    real(RP), intent(in), dimension(IA,JA,2)     :: MAPF
    real(RP), intent(in) :: dt, RDX(IA), RDY(JA), RDZ(KA)
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
    real(RP), intent(in), dimension(KA,IA,JA), optional :: Var_t
    logical, intent(in), optional :: lhist
    real(RP), intent(inout), optional, dimension(KA,IA,JA) :: advch_t, advcv_t

    integer :: IVL, IVL2D

    IVL = mapping_gridtransID(VarLoc)
    IVL2D = mapping_gridtransID(VarLoc, .true.)
    call FDM_EvolveVar( VarA,                                               & ! (out)
         & Var0, Flx(:,:,:,XDIR), Flx(:,:,:,YDIR), Flx(:,:,:,ZDIR),         & ! (in)
         & GSQRT(:,:,:,IVL), MAPF(:,:,IVL2D), dt, RDX, RDY, RDZ,            & ! (in)
         & IIS, IIE, JJS, JJE, KS, KE,                                      & ! (in)
         & Var_t,                                                           & ! (in)
         & lhist, advch_t, advcv_t )
        
  end subroutine ATMOS_NUMERIC_FDM_EvolveVar
  
  subroutine ATMOS_NUMERIC_FDM_EvalFlux(Flux,                                 &  ! (inout)
            & FlxEvalType, VarLoc, Var, VelX_ZUY, VelY_ZXV, VelZ_WXY,         &  ! (in)
            & J13G, J23G, J33G, GSQRT, MAPF, num_diff,                        &  ! (in)
            & IIS, IIE, JJS, JJE, KS, KE )                                       ! (in)

    use scale_atmos_numeric_fdm_cd, only: &
         & FDM_EvalFlxCD4_VarZXY, FDM_EvalFlxCD4_VarZUY, FDM_EvalFlxCD4_VarZXV, FDM_EvalFlxCD4_VarWXY 
    use scale_atmos_numeric_fdm_ud, only: &
         & FDM_EvalFlxUD1_VarZXY, FDM_EvalFlxUD1_VarZUY, FDM_EvalFlxUD1_VarZXV, FDM_EvalFlxUD1_VarWXY 
    
    real(RP), intent(inout), dimension(KA,IA,JA,3)  :: Flux
    integer, intent(in) :: FlxEvalType
    integer, intent(in) :: VarLoc    
    real(RP), intent(in), dimension(KA,IA,JA)     :: Var
    real(RP), intent(in), dimension(KA,IA,JA)     :: VelX_ZUY, VelY_ZXV, VelZ_WXY
    real(RP), intent(in), dimension(KA,IA,JA,7)   :: J13G, J23G, GSQRT, MAPF
    real(RP), intent(in) :: J33G
    real(RP), intent(in), dimension(KA,IA,JA,3) :: num_diff
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE

    procedure(FDM_EvalFlux), pointer :: Invoke_EvalFlux => NULL()

    !
    !

    ! Dynamical dispatching
    select case(FlxEvalType)
    case (FLXEVALTYPE_CD4)
       select case(VarLoc)
       case (VL_ZXY)
          Invoke_EvalFlux => FDM_EvalFlxCD4_VarZXY
       case (VL_ZUY)
          Invoke_EvalFlux => FDM_EvalFlxCD4_VarZUY          
       case (VL_ZXV)
          Invoke_EvalFlux => FDM_EvalFlxCD4_VarZXV          
       case (VL_WXY)
          Invoke_EvalFlux => FDM_EvalFlxCD4_VarWXY          
       case default
          write(*,*) "VarLoc", VarLoc, " in CD4 is not supported."
          stop          
       end select
    case (FLXEVALTYPE_UD1)
       select case(VarLoc)
       case (VL_ZXY)
          Invoke_EvalFlux => FDM_EvalFlxUD1_VarZXY
       case (VL_ZUY)
          Invoke_EvalFlux => FDM_EvalFlxUD1_VarZUY          
       case (VL_ZXV)
          Invoke_EvalFlux => FDM_EvalFlxUD1_VarZXV          
       case (VL_WXY)
          Invoke_EvalFlux => FDM_EvalFlxUD1_VarWXY          
       case default
          write(*,*) "VarLoc", VarLoc, " in UD1 is not supported."
          stop
       end select
    case default
       write(*,*) "FlxEvalType=", FlxEvalType, "is not supported."
       stop
    end select
    
    call Invoke_EvalFlux(Flux(:,:,:,XDIR), Flux(:,:,:,YDIR), Flux(:,:,:,ZDIR),          &  ! (inout)
         & Var, VelX_ZUY, VelY_ZXV, VelZ_WXY, J13G, J23G, J33G, GSQRT, MAPF, num_diff,  &  ! (in)
         & IIS, IIE, JJS, JJE, KS, KE )          

  end subroutine ATMOS_NUMERIC_FDM_EvalFlux

  !------------------------------

  function mapping_gridtransID(VarLoc, to2DFlag) result(Id)
    use scale_gridtrans, only: I_XYZ, I_UYZ, I_XVZ, I_XYW, I_XY, I_UY, I_XV
    integer, intent(in) :: VarLoc
    logical, intent(in), optional :: to2DFlag
    integer :: Id
    select case(VarLoc)
    case (VL_ZXY)
       Id = I_XYZ
       if(present(to2DFlag) .and. to2DFlag) Id = I_XY
    case (VL_ZUY)
       Id = I_UYZ
       if(present(to2DFlag) .and. to2DFlag) Id = I_UY       
    case (VL_ZXV)
       Id = I_XVZ
       if(present(to2DFlag) .and. to2DFlag) Id = I_XV       
    case (VL_WXY)
       Id = I_XYW
       if(present(to2DFlag) .and. to2DFlag) Id = I_XY       
    case (VL_XY)
       Id = I_XY
    case (VL_UY)
       Id = I_UY
    case (VL_XV)
       Id = I_XV
    case default
       write(*,*) "Fail to mapping gridtrans ID."
       stop
    end select
  end function mapping_gridtransID
    
end module scale_atmos_numeric_fdm

