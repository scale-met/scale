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

  use scale_atmos_numeric_fdm_def, only: &
       & FLXEVALTYPE_CD2, FLXEVALTYPENAME_CD2, &
       & FLXEVALTYPE_CD4, FLXEVALTYPENAME_CD4, &
       & FLXEVALTYPE_CD6, FLXEVALTYPENAME_CD6, &
       & FLXEVALTYPE_UD1, FLXEVALTYPENAME_UD1, &
       & FLXEVALTYPE_UD3, FLXEVALTYPENAME_UD3, &
       & FLXEVALTYPE_UD5, FLXEVALTYPENAME_UD5, &
       & VL_XY, VL_UY, VL_XV, &
       & VL_ZXY, VL_ZUY, VL_ZXV, VL_WXY
       
  use scale_atmos_numeric_fdm_util, only: &
       & ATMOS_NUMERIC_FDM_RhoVar2Var => FDM_RhoVar2Var, &
       & FDM_EvolveVar
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  public :: ATMOS_NUMERIC_FDM_Setup
  public :: ATMOS_NUMERIC_FDM_SetCoordMapInfo
  public :: FlxEvalTypeName2ID

  public :: ATMOS_NUMERIC_FDM_RhoVar2Var
  public :: ATMOS_NUMERIC_FDM_EvalFlux
  public :: ATMOS_NUMERIC_FDM_EvolveVar
    
  
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

  type :: CoordinateMapInfo
     real(RP), pointer :: GSQRT(:,:,:,:) => null()
     real(RP), pointer :: J13G(:,:,:,:)  => null()
     real(RP), pointer :: J23G(:,:,:,:)  => null()
     real(RP)          :: J33G
     real(RP), pointer :: MAPF(:,:,:,:)  => null()
  end type CoordinateMapInfo
  type(CoordinateMapInfo), save :: coordMapInfo
  
  abstract interface
     subroutine FDM_EvalFlux(FluxX, FluxY, FluxZ,                            &  ! (inout)
          & Var, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                     &  ! (in)
          & GSQRT, J13G, J23G, J33G, MAPF, GeneralCoordFlag,                 &  ! (in)                
          & IIS, IIE, JJS, JJE, KS_, KE_,                                    &  ! (in)
          & DiffFlxX, DiffFlxY, DiffFlxZ )                                      ! (in)

       use scale_precision,  only: RP
       use scale_grid_index, only: IA, JA, KA

       real(RP), intent(inout), dimension(KA,IA,JA)  :: FluxX, FluxY, FluxZ
       real(RP), intent(in), dimension(KA,IA,JA)     :: Var
       real(RP), intent(in), dimension(KA,IA,JA)     :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
       real(RP), intent(in),  dimension(KA,IA,JA,7)  :: GSQRT, J13G, J23G
       real(RP), intent(in)                          :: J33G
       real(RP), intent(in),  dimension(IA,JA,2,4)   :: MAPF       
       logical, intent(in)                           :: GeneralCoordFlag
       integer,  intent(in) :: IIS, IIE, JJS, JJE, KS_, KE_
       real(RP), intent(in), dimension(KA,IA,JA), optional :: DiffFlxX, DiffFlxY, DiffFlxZ       
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

  subroutine ATMOS_NUMERIC_FDM_SetCoordMapInfo(  &
       & GSQRT, J13G, J23G, J33G, MAPF           & ! (in)
       & )

    real(RP), intent(in), target :: GSQRT(KA,IA,JA,7)
    real(RP), intent(in), target :: J13G (KA,IA,JA,7)
    real(RP), intent(in), target :: J23G (KA,IA,JA,7)
    real(RP), intent(in)         :: J33G
    real(RP), intent(in), target :: MAPF (IA,JA,2,4)

    coordMapInfo%GSQRT => GSQRT
    coordMapInfo%J13G  => J13G
    coordMapInfo%J23G  => J23G
    coordMapInfo%J33G  =  J33G
    coordMapInfo%MAPF  => MAPF
    
  end subroutine ATMOS_NUMERIC_FDM_SetCoordMapInfo

  subroutine ATMOS_NUMERIC_FDM_EvolveVar( VarA,               & ! (out)
       & Var0, VarLoc, Flx,                                   & ! (in)
       & dt, RDX, RDY, RDZ,                                   & ! (in)
       & IIS, IIE, JJS, JJE, KS_, KE_,                        & ! (in)
       & Var_t,                                               & ! (in)
       & lhist, advch_t, advcv_t )

    use scale_process
    
    real(RP), intent(out), dimension(KA,IA,JA)   :: VarA
    real(RP), intent(in), dimension(KA,IA,JA)    :: Var0
    integer, intent(in) :: VarLoc
    real(RP), intent(in), dimension(KA,IA,JA,3)  :: Flx
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: RDX(IA), RDY(JA), RDZ(KA)
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS_, KE_
    real(RP), intent(in), dimension(KA,IA,JA), optional :: Var_t
    logical, intent(in), optional :: lhist
    real(RP), intent(inout), optional, dimension(KA,IA,JA) :: advch_t, advcv_t

    integer :: IVL, IVL2D

    IVL = mapping_gridtransID(VarLoc)
    IVL2D = mapping_gridtransID(VarLoc, .true.)
    call FDM_EvolveVar( VarA,                                               & ! (out)
         & Var0, Flx(:,:,:,XDIR), Flx(:,:,:,YDIR), Flx(:,:,:,ZDIR),         & ! (in)
         & coordMapInfo%GSQRT(:,:,:,IVL), coordMapInfo%MAPF(:,:,:,IVL2D),   & ! (in)
         & dt, RDX, RDY, RDZ,                                               & ! (in)
         & IIS, IIE, JJS, JJE, KS_, KE_,                                    & ! (in)
         & Var_t,                                                           & ! (in)
         & lhist, advch_t, advcv_t )
        
  end subroutine ATMOS_NUMERIC_FDM_EvolveVar

  subroutine ATMOS_NUMERIC_FDM_EvalFlux(Flux,                        &  ! (inout)
       & FlxEvalType, VarLoc, Var,                                   &  ! (in)
       & Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                     &  ! (in)
       & GeneralCoordFlag,                                           &  ! (in)
       & IIS, IIE, JJS, JJE, KS, KE,                                 &  ! (in)
       & DiffFlx )                                                      ! (in)

    use scale_atmos_numeric_fdm_cd, only: &
         & FDM_EvalFlxCD2_VarZXY, FDM_EvalFlxCD2_VarZUY, FDM_EvalFlxCD2_VarZXV, FDM_EvalFlxCD2_VarWXY, &
         & FDM_EvalFlxCD4_VarZXY, FDM_EvalFlxCD4_VarZUY, FDM_EvalFlxCD4_VarZXV, FDM_EvalFlxCD4_VarWXY, &
         & FDM_EvalFlxCD6_VarZXY, FDM_EvalFlxCD6_VarZUY, FDM_EvalFlxCD6_VarZXV, FDM_EvalFlxCD6_VarWXY 

    use scale_atmos_numeric_fdm_ud, only: &
         & FDM_EvalFlxUD1_VarZXY, FDM_EvalFlxUD1_VarZUY, FDM_EvalFlxUD1_VarZXV, FDM_EvalFlxUD1_VarWXY, &
         & FDM_EvalFlxUD3_VarZXY, FDM_EvalFlxUD3_VarZUY, FDM_EvalFlxUD3_VarZXV, FDM_EvalFlxUD3_VarWXY, &
         & FDM_EvalFlxUD5_VarZXY, FDM_EvalFlxUD5_VarZUY, FDM_EvalFlxUD5_VarZXV, FDM_EvalFlxUD5_VarWXY         
    
    real(RP), intent(inout), dimension(KA,IA,JA,3)  :: Flux
    integer, intent(in) :: FlxEvalType
    integer, intent(in) :: VarLoc    
    real(RP), intent(in), dimension(KA,IA,JA)     :: Var
    real(RP), intent(in), dimension(KA,IA,JA)     :: Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY
    logical, intent(in)                           :: GeneralCoordFlag
    integer, intent(in) :: IIS, IIE, JJS, JJE, KS, KE
    real(RP), intent(in), dimension(KA,IA,JA,3), optional   :: DiffFlx    
    
    procedure(FDM_EvalFlux), pointer :: Invoke_EvalFlux => NULL()

    !
    !

    ! Dynamical dispatching
    select case(FlxEvalType)
    case (FLXEVALTYPE_CD2)
       select case(VarLoc)
       case (VL_ZXY)
          Invoke_EvalFlux => FDM_EvalFlxCD2_VarZXY
       case (VL_ZUY)
          Invoke_EvalFlux => FDM_EvalFlxCD2_VarZUY
       case (VL_ZXV)
          Invoke_EvalFlux => FDM_EvalFlxCD2_VarZXV
       case (VL_WXY)
          Invoke_EvalFlux => FDM_EvalFlxCD2_VarWXY
       case default
          write(*,*) "VarLoc", VarLoc, " in CD2 is not supported."
          stop          
       end select
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
    case (FLXEVALTYPE_CD6)
       select case(VarLoc)
       case (VL_ZXY)
          Invoke_EvalFlux => FDM_EvalFlxCD6_VarZXY
       case (VL_ZUY)
          Invoke_EvalFlux => FDM_EvalFlxCD6_VarZUY
       case (VL_ZXV)
          Invoke_EvalFlux => FDM_EvalFlxCD6_VarZXV
       case (VL_WXY)
          Invoke_EvalFlux => FDM_EvalFlxCD6_VarWXY
       case default
          write(*,*) "VarLoc", VarLoc, " in CD6 is not supported."
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
          write(*,*) "VarLoc", VarLoc, " in FDM UD1 is not supported."
          stop
       end select
    case (FLXEVALTYPE_UD3)
       select case(VarLoc)
       case (VL_ZXY)
          Invoke_EvalFlux => FDM_EvalFlxUD3_VarZXY
       case (VL_ZUY)
          Invoke_EvalFlux => FDM_EvalFlxUD3_VarZUY          
       case (VL_ZXV)
          Invoke_EvalFlux => FDM_EvalFlxUD3_VarZXV          
       case (VL_WXY)
          Invoke_EvalFlux => FDM_EvalFlxUD3_VarWXY          
       case default
          write(*,*) "VarLoc", VarLoc, " in FDM UD3 is not supported."
          stop
       end select
    case (FLXEVALTYPE_UD5)
       select case(VarLoc)
       case (VL_ZXY)
          Invoke_EvalFlux => FDM_EvalFlxUD5_VarZXY
       case (VL_ZUY)
          Invoke_EvalFlux => FDM_EvalFlxUD5_VarZUY          
       case (VL_ZXV)
          Invoke_EvalFlux => FDM_EvalFlxUD5_VarZXV          
       case (VL_WXY)
          Invoke_EvalFlux => FDM_EvalFlxUD5_VarWXY          
       case default
          write(*,*) "VarLoc", VarLoc, " in FDM UD5 is not supported."
          stop
       end select
    case default
       write(*,*) "FlxEvalType=", FlxEvalType, "is not supported."
       stop
    end select

    if( present(DiffFlx) ) then
       call Invoke_EvalFlux( &
            & Flux(:,:,:,XDIR), Flux(:,:,:,YDIR), Flux(:,:,:,ZDIR),                         & ! (inout)
            & Var, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                                  & ! (in)
            & coordMapInfo%GSQRT, coordMapInfo%J13G, coordMapInfo%J23G, coordMapInfo%J33G,  & ! (in)
            & coordMapInfo%MAPF, GeneralCoordFlag,                                          & ! (in)  
            & IIS, IIE, JJS, JJE, KS, KE,                                                   & ! (in)
            & DiffFlx(:,:,:,XDIR), DiffFlx(:,:,:,YDIR), DiffFlx(:,:,:,ZDIR) )                 ! (in)
    else
       call Invoke_EvalFlux( &
            & Flux(:,:,:,XDIR), Flux(:,:,:,YDIR), Flux(:,:,:,ZDIR),                         & ! (inout)
            & Var, Dens_ZXY, MomX_ZUY, MomY_ZXV, MomZ_WXY,                                  & ! (in)
            & coordMapInfo%GSQRT, coordMapInfo%J13G, coordMapInfo%J23G, coordMapInfo%J33G,  & ! (in)
            & coordMapInfo%MAPF, GeneralCoordFlag,                                          & ! (in)  
            & IIS, IIE, JJS, JJE, KS, KE )                                                    ! (in)
    end if
  end subroutine ATMOS_NUMERIC_FDM_EvalFlux

  
  function FlxEvalTypeName2ID(name) result(id)
    character(H_SHORT), intent(in) :: name
    integer :: id

    select case(name)
    case(FLXEVALTYPENAME_CD2)
       id = FLXEVALTYPE_CD2
    case(FLXEVALTYPENAME_CD4)
       id = FLXEVALTYPE_CD4       
    case(FLXEVALTYPENAME_CD6)
       id = FLXEVALTYPE_CD6
    case(FLXEVALTYPENAME_UD1)
       id = FLXEVALTYPE_UD1
    case(FLXEVALTYPENAME_UD3)
       id = FLXEVALTYPE_UD3
    case(FLXEVALTYPENAME_UD5)
       id = FLXEVALTYPE_UD5
    case default
       write(*,*) 'FlxEvalTypeName=', trim(name), ' is invalid.'
       stop
    end select
    
  end function FlxEvalTypeName2ID
  
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

