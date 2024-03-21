!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          HEVE FVM scheme for tracer advection in Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tstep_tracer_fvm_heve
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
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
    use scale_prc, only: &
       PRC_abort
    implicit none
    character(len=*), intent(in) :: type

    if ( type /= 'FVM-HEVE' ) then
       LOG_ERROR("ATMOS_DYN_Tstep_tracer_fvm_heve_setup",*) 'Tstep_tracer_type is not "FVM-HEVE"!'
       call PRC_abort
    end if

    return
  end subroutine ATMOS_DYN_Tstep_tracer_fvm_heve_setup

  subroutine ATMOS_DYN_Tstep_tracer_fvm_heve( &
       QTRCo, & ! (out)
       qflx_hi,  & ! (out)
       QTRC, QTRC0, RHOQ_t, &! (in)
       DENS0, DENS, & ! (in)
       mflx_hi, num_diff, & ! (in)
       GSQRT, MAPF, & ! (in)
       CDZ, RCDZ, RCDX, RCDY, & ! (in)
       TwoD, & ! (in)
       dtl, & ! (in)
       FLAG_FCT_TRACER, & ! (in)
       FLAG_FCT_ALONG_STREAM ) ! (in)
    use scale_atmos_dyn_fvm_fct, only: &
       ATMOS_DYN_FVM_fct
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ_tracer, &
       ATMOS_DYN_FVM_fluxX_XYZ_tracer, &
       ATMOS_DYN_FVM_fluxY_XYZ_tracer
    use scale_atmos_dyn_fvm_flux_ud1, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxX_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxY_XYZ_ud1
    implicit none
    real(RP), intent(inout) :: QTRCo   (KA,IA,JA) ! could be identical to QTRC0
    real(RP), intent(out) :: qflx_hi (KA,IA,JA,3) ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
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
    logical,  intent(in)  :: TwoD
    real(RP), intent(in)  :: dtl
    logical,  intent(in)  :: FLAG_FCT_TRACER
    logical,  intent(in)  :: FLAG_FCT_ALONG_STREAM


    ! For tracer advection
    real(RP) :: qflx_lo  (KA,IA,JA,3)  ! rho * vel(x,y,z) * phi,  monotone flux
    real(RP) :: qflx_anti(KA,IA,JA,3)  ! anti-diffusive flux

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: IIS0, JJS0
    integer  :: i, j, k
    !---------------------------------------------------------------------------

    !$acc data copy(QTRCo) &
    !$acc      copyout(qflx_hi) &
    !$acc      copyin(QTRC, QTRC0, RHOQ_t, DENS0, DENS, mflx_hi, num_diff, &
    !$acc             GSQRT, MAPF, CDZ, RCDZ, RCDX, RCDY) &
    !$acc      create(qflx_lo, qflx_anti)

#ifdef DEBUG
    !$acc kernels
    qflx_hi(:,:,:,:) = UNDEF
    qflx_lo(:,:,:,:) = UNDEF
    !$acc end kernels
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
       if ( .not. TwoD ) &
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
               mflx_hi(:,:,:,ZDIR), QTRC0, GSQRT(:,:,:,I_XYW), & ! (in)
               num_diff(:,:,:,ZDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          if ( .not. TwoD ) &
          call ATMOS_DYN_FVM_fluxX_XYZ_ud1( qflx_lo(:,:,:,XDIR), & ! (out)
               mflx_hi(:,:,:,XDIR), QTRC0, GSQRT(:,:,:,I_UYZ), & ! (in)
               num_diff(:,:,:,XDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          call ATMOS_DYN_FVM_fluxY_XYZ_ud1( qflx_lo(:,:,:,YDIR), & ! (out)
               mflx_hi(:,:,:,YDIR), QTRC0, GSQRT(:,:,:,I_XVZ), & ! (in)
               num_diff(:,:,:,YDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)
       end if

    enddo
    enddo

    if ( FLAG_FCT_TRACER ) then

       call ATMOS_DYN_FVM_fct( qflx_anti,        & ! (out)
                           QTRC0, DENS0, DENS,   & ! (in)
                           qflx_hi, qflx_lo,     & ! (in)
                           mflx_hi,              & ! (in)
                           RCDZ, RCDX, RCDY,     & ! (in)
                           GSQRT(:,:,:,I_XYZ),   & ! (in)
                           MAPF, TwoD, dtl,      & ! (in)
                           FLAG_FCT_ALONG_STREAM ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS-1, KE
             qflx_hi(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) - qflx_anti(k,i,j,ZDIR)
          end do
          end do
          end do

          if ( .not. TwoD ) then
             if ( IIS == IS ) then
                IIS0 = IIS-1
             else
                IIS0 = IIS
             end if
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS0, IIE
             do k = KS, KE
                qflx_hi(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) - qflx_anti(k,i,j,XDIR)
             end do
             end do
             end do
          end if

          if ( JJS == JS ) then
             JJS0 = JJS-1
          else
             JJS0 = JJS
          end if
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS0, JJE
          do i = IIS, IIE
          do k = KS, KE
             qflx_hi(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) - qflx_anti(k,i,j,YDIR)
          end do
          end do
          end do


          if ( TwoD ) then
             !$omp parallel do private(i,j,k) OMP_SCHEDULE_
             do j = JJS, JJE
             do k = KS, KE
                QTRCo(k,IS,j) = ( QTRC0(k,IS,j) * DENS0(k,IS,j) &
                                + dtl * ( - ( ( qflx_hi(k,IS,j,ZDIR) - qflx_hi(k-1,IS,j  ,ZDIR)  ) * RCDZ(k) &
                                            + ( qflx_hi(k,IS,j,YDIR) - qflx_hi(k  ,IS,j-1,YDIR)  ) * RCDY(j) &
                                          ) * MAPF(IS,j,2) / GSQRT(k,IS,j,I_XYZ) &
                                + RHOQ_t(k,IS,j) ) ) / DENS(k,IS,j)
             enddo
             enddo
          else
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
          end if

       enddo
       enddo

    else ! skip FCT

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          if ( TwoD ) then
             !$omp parallel do default(none) private(j,k) OMP_SCHEDULE_ &
             !$omp shared(JJS,JJE,IS,KS,KE,QTRCo,QTRC0,DENS0,dtl,qflx_hi,RCDZ,RCDY,MAPF) &
             !$omp shared(GSQRT,RHOQ_t,DENS,I_XYZ) 
             !$acc kernels
             do j = JJS, JJE
             do k = KS, KE
                QTRCo(k,IS,j) = ( QTRC0(k,IS,j) * DENS0(k,IS,j) &
                                + dtl * ( - ( ( qflx_hi(k,IS,j,ZDIR) - qflx_hi(k-1,IS,j  ,ZDIR)  ) * RCDZ(k) &
                                            + ( qflx_hi(k,IS,j,YDIR) - qflx_hi(k  ,IS,j-1,YDIR)  ) * RCDY(j) &
                                           ) * MAPF(IS,j,2) / GSQRT(k,IS,j,I_XYZ) &
                               + RHOQ_t(k,IS,j) ) ) / DENS(k,IS,j)
             enddo
             enddo
             !$acc end kernels
          else
             !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
             !$omp shared(JJS,JJE,IIS,IIE,KS,KE,QTRCo,QTRC0,DENS0,dtl,qflx_hi,RCDZ,RCDX,RCDY,MAPF) &
             !$omp shared(GSQRT,RHOQ_t,DENS,I_XYZ) 
             !$acc kernels
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
             !$acc end kernels
          end if

       enddo
       enddo

    end if
    !$acc end data

    return
  end subroutine ATMOS_DYN_Tstep_tracer_fvm_heve

end module scale_atmos_dyn_tstep_tracer_fvm_heve
