#include "inc_openmp.h"
module scale_atmos_dyn_advtest
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
  use scale_process
#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif

  use scale_atmos_numeric_fdm_def, only: &
       & FLXEVALTYPE_CD2, FLXEVALTYPE_UD1,  &
       & FLXEVALTYPE_CD4, FLXEVALTYPE_UD3,  &
       & FLXEVALTYPE_CD6, FLXEVALTYPE_UD5,  &
       & VL_ZXY, VL_ZUY, VL_ZXV, VL_WXY
  
  use scale_atmos_numeric_fdm, only: &
       & ATMOS_NUMERIC_FDM_setup,           &
       & ATMOS_NUMERIC_FDM_SetCoordMapInfo, &
       & ATMOS_NUMERIC_FDM_RhoVar2Var,      &
       & ATMOS_NUMERIC_FDM_EvalFlux,        &
       & ATMOS_NUMERIC_FDM_EvolveVar

  !-----------------------------------------------------------------------------
  implicit none
  private

  integer :: FlxEvalTypeID
  
  public :: ATMOS_DYN_RK3_advtest
  public :: ATMOS_DYN_RK3_advtest_SetFluxEvalType
contains
  !-----------------------------------------------------------------------------
  !> Specify type of finite difference scheme
  subroutine ATMOS_DYN_RK3_advtest_SetFluxEvalType(FlxEvalTypeName)

    use scale_atmos_numeric_fdm, only: &
         & FlxEvalTypeName2ID
    
    character(len=H_SHORT), intent(in) :: FlxEvalTypeName

    FlxEvalTypeID = FlxEvalTypeName2ID(FlxEvalTypeName)
    
  end subroutine ATMOS_DYN_RK3_advtest_SetFluxEvalType

  subroutine ATMOS_DYN_RK3_advtest(PROG,  &
    & PROG0, DENS, mflx_hi, dt, RCDX, RCDY, RCDZ )
    real(RP), intent(inout) :: PROG(KA,IA,JA)
    real(RP), intent(in) :: PROG0(KA,IA,JA)
    real(RP), intent(in) :: mflx_hi(KA,IA,JA,3), DENS(KA,IA,JA,3)
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: RCDX(IA), RCDY(JA), RCDZ(KA)
    
    integer :: IIS, IIE, JJS, JJE

    
    real(RP) :: RHS_RKTMP(KA,IA,JA,4)
    real(RP) :: PROG_RKTMP(KA,IA,JA)

    call calc_RHS(RHS_RKTMP(:,:,:,1), PROG0, &
         & DENS, mflx_hi, RCDX, RCDY, RCDZ)

    call calc_RHS(RHS_RKTMP(:,:,:,2), PROG0 + 0.5_RP*dt*RHS_RKTMP(:,:,:,1), &
         & DENS, mflx_hi, RCDX, RCDY, RCDZ)
    
    call calc_RHS(RHS_RKTMP(:,:,:,3), PROG0 + 0.5_RP*dt*RHS_RKTMP(:,:,:,2), &
         & DENS, mflx_hi, RCDX, RCDY, RCDZ)

    call calc_RHS(RHS_RKTMP(:,:,:,4), PROG0 + dt*RHS_RKTMP(:,:,:,3), &
         & DENS, mflx_hi, RCDX, RCDY, RCDZ)

!!$
    PROG = PROG0 + dt*(RHS_RKTMP(:,:,:,1) + 2.0_RP*RHS_RKTMP(:,:,:,2) + 2.0_RP*RHS_RKTMP(:,:,:,3)  + RHS_RKTMP(:,:,:,4))/6.0_RP
!!$    PROG = PROG  + dt * RHS_RKTMP(:,:,:,1)
!!$    write(*,*) "--", PROG(KS,IS:IE,JS)
!!$    write(*,*) "**", RHS_RKTMP(KS,IS:IE,JS,1)

  end subroutine ATMOS_DYN_RK3_advtest

  subroutine calc_RHS(RHS, PROG, DENS, mflx_hi, RCDX, RCDY, RCDZ)
    real(RP), intent(inout) :: RHS(KA,IA,JA)
    real(RP), intent(in) :: PROG(KA,IA,JA), DENS(KA,IA,JA)
    real(RP), intent(in) :: mflx_hi(KA,IA,JA,3)
    real(RP), intent(in) :: RCDX(IA), RCDY(JA), RCDZ(KA)
    
    integer :: IIS, IIE, JJS, JJE
    real(RP) :: VARTMP(KA, IA, JA)
    real(RP) :: qflx_hi(KA, IA, JA, 3)
    real(RP) :: one(KA, IA, JA)
    real(RP) :: zero(KA, IA, JA)


    
    !####################################################################################
    ! PROG0 (advection test)
    !####################################################################################

    one = 1.0_RP
    zero = 0.0_RP
    
    ! RHOQ --> Q
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
      call ATMOS_NUMERIC_FDM_RhoVar2Var( VARTMP,                        &  ! (out)
         & PROG, VL_ZXY, one,                                           &  ! (in)
         & IIS-IHALO, IIE+IHALO, JJS-JHALO, JJE+JHALO, KS, KE )            ! (in)
    enddo
    enddo


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
    
      !-----< high order flux >-----
      ! * Input (Phi * Vec_i / Fact)
      !   Phi = POTT, Vec_i = MomFlx{X,Y,Z}, Fact = 1.0
      ! * Ouput
      !   [ POTT*MOMFlxX, POTT*MOMFlxY, POTT*MOMFlxZ ]
       call ATMOS_NUMERIC_FDM_EvalFlux( qflx_hi,                                                &  ! (inout)
        & FlxEvalTypeID, VL_ZXY,                                                               &  ! (in)
        & VARTMP, DENS, mflx_hi(:,:,:,XDIR), mflx_hi(:,:,:,YDIR), mflx_hi(:,:,:,ZDIR), .false., &  ! (in)
        & IIS, IIE, JJS, JJE, KS, KE  )                          ! (in)

    enddo
    enddo

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
    
      !-----< update rho*theta >-----
      call ATMOS_NUMERIC_FDM_EvolveVar( RHS,                      & ! (out)
         & zero, VL_ZXY, qflx_hi, 1.0_RP, RCDX, RCDY, RCDZ,            & ! (in)
         & IIS, IIE, JJS, JJE, KS, KE                                          & ! (in)
         & )

    enddo
    enddo
    
  end subroutine calc_RHS
end module scale_atmos_dyn_advtest
