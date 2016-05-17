!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          HEVE FVM scheme for tracer advection in Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-05-17 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tstep_tracer_fvm_heve
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
  public :: ATMOS_DYN_Tstep_tracer_fvm_heve_setup
  public :: ATMOS_DYN_Tstep_tracer_fvm_heve

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
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_tracer_fvm_heve_setup( type )
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    character(len=*), intent(in) :: type

    if ( type .ne. 'FVM-HEVE' ) then
       write(*,*) 'xxx Tstep_tracer_type is not "FVM-HEVE"!'
       call PRC_MPIstop
    end if

    return
  end subroutine ATMOS_DYN_Tstep_tracer_fvm_heve_setup

  subroutine ATMOS_DYN_Tstep_tracer_fvm_heve( &
       QTRCo, & ! (out)
       QTRC, QTRC0, RHOQ_t, &! (in)
       DENS0, DENS, & ! (in)
       mflx_hi, num_diff, & ! (in)
       GSQRT, MAPF, & ! (in)
       CDZ, RCDZ, RCDX, RCDY, & ! (in)
       dtl, & ! (in)
       FLAG_FCT_TRACER, & ! (in)
       FLAG_FCT_ALONG_STREAM ) ! (in)
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYZ, &
       I_XVZ, &
       I_XY,  &
       I_UY,  &
       I_XV
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_fct
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ_tracer, &
       ATMOS_DYN_FVM_fluxX_XYZ_tracer, &
       ATMOS_DYN_FVM_fluxY_XYZ_tracer
    use scale_atmos_dyn_fvm_flux_ud1, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxX_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxY_XYZ_ud1
    implicit none
    real(RP), intent(out) :: QTRCo   (KA,IA,JA)
    real(RP), intent(in)  :: QTRC    (KA,IA,JA)
    real(RP), intent(in)  :: QTRC0   (KA,IA,JA)
    real(RP), intent(in)  :: RHOQ_t  (KA,IA,JA)
    real(RP), intent(in)  :: DENS0   (KA,IA,JA)
    real(RP), intent(in)  :: DENS    (KA,IA,JA)
    real(RP), intent(in)  :: mflx_hi (KA,IA,JA,3)
    real(RP), intent(in)  :: num_diff(KA,IA,JA,3)
    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7)
    real(RP), intent(in)  :: MAPF    (IA,JA,2)
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: RCDX(IA)
    real(RP), intent(in)  :: RCDY(JA)
    real(RP), intent(in)  :: dtl
    logical,  intent(in)  :: FLAG_FCT_TRACER
    logical,  intent(in)  :: FLAG_FCT_ALONG_STREAM


    ! For tracer advection
    real(RP) :: qflx_hi  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(RP) :: qflx_lo  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3)  ! anti-diffusive flux

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: i, j, k
    !---------------------------------------------------------------------------

#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
    qflx_lo(:,:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! at (x, y, w)
       call ATMOS_DYN_FVM_fluxZ_XYZ_tracer( qflx_hi(:,:,:,ZDIR), & ! (out)
            mflx_hi(:,:,:,ZDIR), QTRC, GSQRT(:,:,:,I_XYW), & ! (in)
            num_diff(:,:,:,ZDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, y, z)
       call ATMOS_DYN_FVM_fluxX_XYZ_tracer( qflx_hi(:,:,:,XDIR), & ! (out)
            mflx_hi(:,:,:,XDIR), QTRC, GSQRT(:,:,:,I_UYZ), & ! (in)
            num_diff(:,:,:,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, v, z)
       call ATMOS_DYN_FVM_fluxY_XYZ_tracer( qflx_hi(:,:,:,YDIR), & ! (out)
            mflx_hi(:,:,:,YDIR), QTRC, GSQRT(:,:,:,I_XVZ), & ! (in)
            num_diff(:,:,:,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       
       if ( FLAG_FCT_TRACER ) then

          call ATMOS_DYN_FVM_fluxZ_XYZ_ud1( qflx_lo(:,:,:,ZDIR), & ! (out)
               mflx_hi(:,:,:,ZDIR), QTRC, GSQRT(:,:,:,I_XYZ), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          call ATMOS_DYN_FVM_fluxX_XYZ_ud1( qflx_lo(:,:,:,XDIR), & ! (out)
               mflx_hi(:,:,:,XDIR), QTRC, GSQRT(:,:,:,I_UYZ), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          call ATMOS_DYN_FVM_fluxY_XYZ_ud1( qflx_lo(:,:,:,YDIR), & ! (out)
               mflx_hi(:,:,:,YDIR), QTRC, GSQRT(:,:,:,I_XVZ), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)
       end if

    enddo
    enddo

    if ( FLAG_FCT_TRACER ) then

       call ATMOS_DYN_fct( qflx_anti,            & ! (out)
                           QTRC, DENS0, DENS,    & ! (in)
                           qflx_hi, qflx_lo,     & ! (in)
                           mflx_hi,              & ! (in)
                           RCDZ, RCDX, RCDY,     & ! (in)
                           GSQRT(:,:,:,I_XYZ),   & ! (in)
                           MAPF, dtl,            & ! (in)
                           FLAG_FCT_ALONG_STREAM ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             QTRCo(k,i,j) = ( QTRC0(k,i,j) * DENS0(k,i,j) &
                            + dtl * ( - ( ( qflx_hi(k  ,i  ,j  ,ZDIR) - qflx_anti(k  ,i  ,j  ,ZDIR) &
                                          - qflx_hi(k-1,i  ,j  ,ZDIR) + qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                        + ( qflx_hi(k  ,i  ,j  ,XDIR) - qflx_anti(k  ,i  ,j  ,XDIR) &
                                          - qflx_hi(k  ,i-1,j  ,XDIR) + qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                        + ( qflx_hi(k  ,i  ,j  ,YDIR) - qflx_anti(k  ,i  ,j  ,YDIR) &
                                          - qflx_hi(k  ,i  ,j-1,YDIR) + qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) &
                                      ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j,I_XYZ) &
                               + RHOQ_t(k,i,j) ) ) / DENS(k,i,j)
          enddo
          enddo
          enddo

       enddo
       enddo

    else ! skip FCT

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             QTRCo(k,i,j) = ( QTRC0(k,i,j) * DENS0(k,i,j) &
                            + dtl * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR)  ) * RCDZ(k) &
                                        + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR)  ) * RCDX(i) &
                                        + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR)  ) * RCDY(j) &
                                        ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j,I_XYZ) &
                            + RHOQ_t(k,i,j) ) ) / DENS(k,i,j)
          enddo
          enddo
          enddo

       enddo
       enddo

    end if

    return
  end subroutine ATMOS_DYN_Tstep_tracer_fvm_heve

end module scale_atmos_dyn_tstep_tracer_fvm_heve
