!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          Runge-Kutta for Atmospheric dynamical process
!!          HEVE, no terrain + tracer FCT limiter
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-18 (S.Nishizawa) [new] splited from dynamical core
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tracer
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
       & FLXEVALTYPE_CD4, FLXEVALTYPE_UD1,  &
       & VL_ZXY, VL_ZUY, VL_ZXV, VL_WXY
  
  use scale_atmos_numeric_fdm, only: &
       & ATMOS_NUMERIC_FDM_setup,           &
       & ATMOS_NUMERIC_FDM_SetCoordMapInfo, &
       & ATMOS_NUMERIC_FDM_EvalFlux,        &
       & ATMOS_NUMERIC_FDM_EvolveVar

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_tracer_setup
  public :: ATMOS_DYN_tracer_update

  public :: ATMOS_DYN_tracer_SetFluxEvalType

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer :: IFS_OFF
  integer :: JFS_OFF

  integer :: FlxEvalTypeID
  
  !-----------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_tracer_setup( &
       BND_W, &
       BND_E, &
       BND_S, &
       BND_N  )
    implicit none

    logical, intent(in) :: BND_W
    logical, intent(in) :: BND_E
    logical, intent(in) :: BND_S
    logical, intent(in) :: BND_N
    !---------------------------------------------------------------------------

    IFS_OFF = 1
    JFS_OFF = 1
    if ( BND_W ) IFS_OFF = 0
    if ( BND_S ) JFS_OFF = 0

    !
    FlxEvalTypeID = FLXEVALTYPE_CD4
    call ATMOS_NUMERIC_FDM_setup(IFS_OFF, JFS_OFF)
    
    return
  end subroutine ATMOS_DYN_tracer_setup

  !-----------------------------------------------------------------------------
  !> Specify type of finite difference scheme
  subroutine ATMOS_DYN_tracer_SetFluxEvalType(FlxEvalTypeName)

    use scale_atmos_numeric_fdm, only: &
         & FlxEvalTypeName2ID
    
    character(len=H_SHORT), intent(in) :: FlxEvalTypeName

    FlxEvalTypeID = FlxEvalTypeName2ID(FlxEvalTypeName)
    
  end subroutine ATMOS_DYN_tracer_SetFluxEvalType

  !-----------------------------------------------------------------------------
  !> Update tracer
  subroutine ATMOS_DYN_tracer_update( &
       & QTRC,                                            &  ! (inout)
       & DENS, DENS00, mflx_hi, RHOQ_t, num_diff_q,       &  ! (in)
       & RCDX, RCDY, RCDZ, GSQRT, J13G, J23G, J33G, MAPF, &  ! (in)
       & dt, FLAG_FCT_ALONG_STREAM                        &  ! (in)
       & )

    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_fct
    use scale_gridtrans, only: &
         I_XYZ, &
         I_XY
    
    real(RP), intent(inout) :: QTRC(KA,IA,JA)
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: DENS00(KA,IA,JA)
    real(RP), intent(in) :: mflx_hi(KA,IA,JA,3)
    real(RP), intent(in) :: RHOQ_t(KA,IA,JA)
    real(RP), intent(in) :: num_diff_q(KA,IA,JA,3)
    real(RP), intent(in) :: RCDX(IA), RCDY(JA), RCDZ(KA)
    real(RP), intent(in), dimension(KA,IA,JA,7) :: GSQRT, J13G, J23G
    real(RP), intent(in) :: J33G
    real(RP), intent(in) :: MAPF(IA,JA,2,4)
    real(RP), intent(in) :: dt
    logical,  intent(in)  :: FLAG_FCT_ALONG_STREAM
    
    real(RP) :: qflx_hi(KA,IA,JA,3)    
    real(RP) :: qflx_lo(KA,IA,JA,3)
    real(RP) :: qflx_anti(KA,IA,JA,3)
    
    real(RP) :: VARTMP(KA,IA,JA)
    integer :: i, j, k
    integer :: IIS, IIE, JJS, JJE

    VARTMP = 1.0_RP

    call ATMOS_NUMERIC_FDM_SetCoordMapInfo(      &
       & GSQRT, J13G, J23G, J33G, MAPF           & ! (in)
       & )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
       
      !-----< high order flux >-----
      ! * Input (Phi * Vec_i / Fact)
      !   Phi = QTRC, Vec_i = MomFlx{X,Y,Z}, Fact = 1.0
      ! * Ouput
      !   [ QTRC*MOMFlxX, QTRC*MOMFlxY, QTRC*MOMFlxZ ]
      call ATMOS_NUMERIC_FDM_EvalFlux( qflx_hi,                                                 &  ! (inout)
        & FlxEvalTypeID, VL_ZXY,                                                                &  ! (in)
        & QTRC, VARTMP, mflx_hi(:,:,:,XDIR), mflx_hi(:,:,:,YDIR), mflx_hi(:,:,:,ZDIR), .false., &  ! (in)
        & IIS, IIE, JJS, JJE, KS, KE, num_diff_q )                                                 ! (in)

      !-----< low order flux >-----
      ! * Input (Phi * Vec_i / Fact)
      !   Phi = QTRC, Vec_i = MomFlx{X,Y,Z}, Fact = 1.0
      ! * Ouput
      !   [ QTRC*MOMFlxX, QTRC*MOMFlxY, QTRC*MOMFlxZ ]
      ! * Note
      !   In order to obtain interpolated values at upstream position in FCT scheme,
      !   fluxes at outer region are calculated. 
      call ATMOS_NUMERIC_FDM_EvalFlux( qflx_lo,                                                 &  ! (inout)
        & FLXEVALTYPE_UD1, VL_ZXY,                                                              &  ! (in)
        & QTRC, VARTMP, mflx_hi(:,:,:,XDIR), mflx_hi(:,:,:,YDIR), mflx_hi(:,:,:,ZDIR), .false., &  ! (in)
        & IIS-1, IIE+1, JJS-1, JJE+1, KS, KE )                                                     ! (in)

    enddo
    enddo


#ifdef HORDCFDDYN_FCT_Q_ON    
    call ATMOS_DYN_fct( qflx_anti,                    & ! (out)
                        QTRC, DENS00, DENS,           & ! (in)
                        qflx_hi, qflx_lo,             & ! (in)
                        mflx_hi,                      & ! (in)
                        RCDZ, RCDX, RCDY,             & ! (in)
                        GSQRT(:,:,:,I_XYZ),           & ! (in)
                        MAPF(:,:,:,I_XY), dt,         & ! (in)
                        FLAG_FCT_ALONG_STREAM         ) ! (in)
    write(*,*) "FCT called"
#else
    qflx_anti = 0.0_RP
#endif
    
    !--- < update rho*q >
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

      !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
      do j = JJS, JJE
      do i = IIS, IIE
        do k = KS, KE
          QTRC(k,i,j) = (  QTRC(k,i,j) * DENS00(k,i,j) &
                            + dt * ( - ( ( qflx_hi(k  ,i  ,j  ,ZDIR) - qflx_anti(k  ,i  ,j  ,ZDIR) &
                                         - qflx_hi(k-1,i  ,j  ,ZDIR) + qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k  ,i  ,j  ,XDIR) - qflx_anti(k  ,i  ,j  ,XDIR) &
                                         - qflx_hi(k  ,i-1,j  ,XDIR) + qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                       + ( qflx_hi(k  ,i  ,j  ,YDIR) - qflx_anti(k  ,i  ,j  ,YDIR) &
                                         - qflx_hi(k  ,i  ,j-1,YDIR) + qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) &
                                       ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) &
                                     + RHOQ_t(k,i,j) ) ) / DENS(k,i,j)
        enddo
      enddo
      enddo
              
    enddo
    enddo
      
    
  end subroutine ATMOS_DYN_tracer_update
  
end module scale_atmos_dyn_tracer
