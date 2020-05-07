!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!          Flux-corrected transport (FCT) scheme for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_fvm_fct
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
  public :: ATMOS_DYN_fvm_fct

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: get_fact_fct

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------------
  !> Setup
 
  !-----------------------------------------------------------------------------
  !> Flux Correction Transport Limiter
  subroutine ATMOS_DYN_fvm_fct( &
       qflx_anti,           &
       phi_in, DENS0, DENS, &
       qflx_hi, qflx_lo,    &
       mflx_hi,             &
       rdz, rdx, rdy,       &
       GSQRT, MAPF,         &
       TwoD, dt,            &
       flag_vect )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       IUNDEF => CONST_UNDEF2, &
       EPSILON => CONST_EPS
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: qflx_anti(KA,IA,JA,3)

    real(RP), intent(in) :: phi_in(KA,IA,JA) ! physical quantity
    real(RP), intent(in) :: DENS0(KA,IA,JA)
    real(RP), intent(in) :: DENS (KA,IA,JA)

    real(RP), intent(in) :: qflx_hi(KA,IA,JA,3)
    real(RP), intent(in) :: qflx_lo(KA,IA,JA,3)
    real(RP), intent(in) :: mflx_hi(KA,IA,JA,3)

    real(RP), intent(in) :: RDZ(:)
    real(RP), intent(in) :: RDX(:)
    real(RP), intent(in) :: RDY(:)

    real(RP), intent(in) :: GSQRT(KA,IA,JA) !< vertical metrics {G}^1/2
    real(RP), intent(in) :: MAPF(IA,JA,2)   !< map factor

    logical,  intent(in) :: TwoD
    real(RP), intent(in) :: dt

    logical, intent(in) :: flag_vect

    ! work for FCT
    real(RP) :: phi_lo(KA,IA,JA)
    real(RP) :: pjpls(KA,IA,JA)
    real(RP) :: pjmns(KA,IA,JA)
    real(RP) :: qjpls(KA,IA,JA)
    real(RP) :: qjmns(KA,IA,JA)
    real(RP) :: rjpls(KA,IA,JA)
    real(RP) :: rjmns(KA,IA,JA)

    real(RP) :: qmin, qmax
    real(RP) :: zerosw, dirsw

    real(RP) :: fact(0:1,-1:1,-1:1)
    real(RP) :: rw, ru, rv
    real(RP) :: qa_in, qb_in
    real(RP) :: qa_lo, qb_lo

    integer :: k, i, j, ijs
    integer :: IIS, IIE, JJS, JJE
    !---------------------------------------------------------------------------

#ifdef DEBUG
    qflx_anti(:,:,:,:) = UNDEF

    pjpls(:,:,:) = UNDEF
    pjmns(:,:,:) = UNDEF
    qjpls(:,:,:) = UNDEF
    qjmns(:,:,:) = UNDEF
    rjpls(:,:,:) = UNDEF
    rjmns(:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,ZDIR) )
#endif
          qflx_anti(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) - qflx_lo(k,i,j,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_anti(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_anti(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       if ( .not. TwoD ) then
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS  , JJE
          do i = IIS-1, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_hi(k,i,j,XDIR) )
             call CHECK( __LINE__, qflx_lo(k,i,j,XDIR) )
#endif
             qflx_anti(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) - qflx_lo(k,i,j,XDIR)
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,YDIR) )
#endif
          qflx_anti(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) - qflx_lo(k,i,j,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update monotone scheme
       if ( TwoD ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j = JJS-1, JJE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, phi_in(k,IS,j) )
             call CHECK( __LINE__, qflx_lo(k  ,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_lo(k-1,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_lo(k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, qflx_lo(k  ,IS,j-1,YDIR) )
#endif
             phi_lo(k,IS,j) = ( phi_in(k,IS,j) * DENS0(k,IS,j) &
                            + dt * ( - ( ( qflx_lo(k,IS,j,ZDIR)-qflx_lo(k-1,IS,j  ,ZDIR) ) * RDZ(k) &
                                       + ( qflx_lo(k,IS,j,YDIR)-qflx_lo(k  ,IS,j-1,YDIR) ) * RDY(j) &
                                       ) * MAPF(IS,j,2) / GSQRT(k,IS,j)                  ) &
                             ) / DENS(k,IS,j)
          enddo
          enddo
       else
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, phi_in(k,i,j) )
             call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_lo(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_lo(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_lo(k  ,i  ,j-1,YDIR) )
#endif
             phi_lo(k,i,j) = ( phi_in(k,i,j) * DENS0(k,i,j) &
                             + dt * ( - ( ( qflx_lo(k,i,j,ZDIR)-qflx_lo(k-1,i  ,j  ,ZDIR) ) * RDZ(k) &
                                        + ( qflx_lo(k,i,j,XDIR)-qflx_lo(k  ,i-1,j  ,XDIR) ) * RDX(i) &
                                        + ( qflx_lo(k,i,j,YDIR)-qflx_lo(k  ,i  ,j-1,YDIR) ) * RDY(j) &
                                        ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)                 ) &
                             ) / DENS(k,i,j)
          enddo
          enddo
          enddo
       end if
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net incoming quantity change by antidiffusive flux
       if ( TwoD ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j = JJS, JJE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_anti(k  ,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,IS,j-1,YDIR) )
#endif
             pjpls(k,IS,j) = dt * ( ( max(0.0_RP,qflx_anti(k-1,IS ,j  ,ZDIR)) - min(0.0_RP,qflx_anti(k,IS,j,ZDIR)) ) * RDZ(k) &
                                  + ( max(0.0_RP,qflx_anti(k  ,IS,j-1,YDIR)) - min(0.0_RP,qflx_anti(k,IS,j,YDIR)) ) * RDY(j) &
                                 ) * MAPF(IS,j,2) / GSQRT(k,IS,j)
          enddo
          enddo
       else
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             pjpls(k,i,j) = dt * ( ( max(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) - min(0.0_RP,qflx_anti(k,i,j,ZDIR)) ) * RDZ(k) &
                                 + ( max(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) - min(0.0_RP,qflx_anti(k,i,j,XDIR)) ) * RDX(i) &
                                 + ( max(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) - min(0.0_RP,qflx_anti(k,i,j,YDIR)) ) * RDY(j) &
                                 ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
          enddo
          enddo
          enddo
       end if
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net outgoing quantity change by antidiffusive flux
       if ( TwoD ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_
          do j = JJS, JJE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_anti(k  ,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,IS,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,IS,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,IS,j-1,YDIR) )
#endif
             pjmns(k,IS,j) = dt * ( ( max(0.0_RP,qflx_anti(k,IS,j,ZDIR)) - min(0.0_RP,qflx_anti(k-1,IS,j  ,ZDIR)) ) * RDZ(k) &
                                  + ( max(0.0_RP,qflx_anti(k,IS,j,YDIR)) - min(0.0_RP,qflx_anti(k  ,IS,j-1,YDIR)) ) * RDY(j) &
                                  ) * MAPF(IS,j,2) / GSQRT(k,IS,j)
          enddo
          enddo
       else
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             pjmns(k,i,j) = dt * ( ( max(0.0_RP,qflx_anti(k,i,j,ZDIR)) - min(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) ) * RDZ(k) &
                                 + ( max(0.0_RP,qflx_anti(k,i,j,XDIR)) - min(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) ) * RDX(i) &
                                 + ( max(0.0_RP,qflx_anti(k,i,j,YDIR)) - min(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) ) * RDY(j) &
                                 ) * MAPF(i,j,1) * MAPF(i,j,2) / GSQRT(k,i,j)
          enddo
          enddo
          enddo
       end if
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc allowable range or quantity change by antidiffusive flux

       if (flag_vect) then

          if ( TwoD ) then
             i = IS
             !$omp parallel do private(j,k,rw,ru,rv,fact,qa_in,qb_in,qa_lo,qb_lo,qmax,qmin) OMP_SCHEDULE_
             do j = JJS, JJE
             do k = KS+1, KE-1
                rw = (mflx_hi(k,i,j,ZDIR)+mflx_hi(k-1,i  ,j  ,ZDIR)) * RDZ(k) ! 2 * rho * w / dz
                ru = 0.0_RP
                rv = (mflx_hi(k,i,j,YDIR)+mflx_hi(k  ,i  ,j-1,YDIR)) * RDY(j) ! 2 * rho * v / dy

                call get_fact_fct( fact, & ! (out)
                                   rw, ru, rv ) ! (in)

                qa_in = fact(1, 0, 1) * phi_in(k+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_in(k  ,i  ,j+1) &
                      + fact(1, 0, 0) * phi_in(k+1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_in(k+1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_in(k  ,i  ,j  )
                qb_in = fact(1, 0, 1) * phi_in(k-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_in(k  ,i  ,j-1) &
                      + fact(1, 0, 0) * phi_in(k-1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_in(k-1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_in(k  ,i  ,j  )
                qa_lo = fact(1, 0, 1) * phi_lo(k+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_lo(k  ,i  ,j+1) &
                      + fact(1, 0, 0) * phi_lo(k+1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_lo(k+1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_lo(k  ,i  ,j  )
                qb_lo = fact(1, 0, 1) * phi_lo(k-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_lo(k  ,i  ,j-1) &
                      + fact(1, 0, 0) * phi_lo(k-1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_lo(k-1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_lo(k  ,i  ,j  )

                qmax = max( &
                     phi_in(k,i,j), qa_in, qb_in, &
                     phi_lo(k,i,j), qa_lo, qb_lo  )
                qmin = min( &
                     phi_in(k,i,j), qa_in, qb_in, &
                     phi_lo(k,i,j), qa_lo, qb_lo  )
                qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
                qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
             end do
             end do

             ! k = KS
             !$omp parallel do private(j,rw,ru,rv,fact,qa_in,qb_in,qa_lo,qb_lo,qmax,qmin) OMP_SCHEDULE_
             do j = JJS, JJE
                rw = (mflx_hi(KS,i,j,ZDIR)                           ) * RDZ(KS)! 2 * rho * w / dz
                ru = 0.0_RP
                rv = (mflx_hi(KS,i,j,YDIR)+mflx_hi(KS  ,i  ,j-1,YDIR)) * RDY(j) ! 2 * rho * v / dy

                call get_fact_fct( fact, & ! (out)
                                   rw, ru, rv ) ! (in)

                qa_in = fact(1, 0, 1) * phi_in(KS+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_in(KS  ,i  ,j+1) &
                      + fact(1, 0, 0) * phi_in(KS+1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_in(KS+1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_in(KS  ,i  ,j  )
                qb_in = fact(1, 0, 1) * phi_in(KS  ,i  ,j-1) &
                      + fact(0, 0, 1) * phi_in(KS  ,i  ,j-1) &
                      + fact(1, 0, 0) * phi_in(KS  ,i  ,j  ) &
                      + fact(1, 0,-1) * phi_in(KS  ,i  ,j-1) &
                      + fact(0, 0, 0) * phi_in(KS  ,i  ,j  )
                qa_lo = fact(1, 0, 1) * phi_lo(KS+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_lo(KS  ,i  ,j+1) &
                      + fact(1, 0, 0) * phi_lo(KS+1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_lo(KS+1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_lo(KS  ,i  ,j  )
                qb_lo = fact(1, 0, 1) * phi_lo(KS  ,i  ,j-1) &
                      + fact(0, 0, 1) * phi_lo(KS  ,i  ,j-1) &
                      + fact(1, 0, 0) * phi_lo(KS  ,i  ,j  ) &
                      + fact(1, 0,-1) * phi_lo(KS  ,i  ,j-1) &
                      + fact(0, 0, 0) * phi_lo(KS  ,i  ,j  )

                qmax = max( &
                     phi_in(KS,i,j), qa_in, qb_in, &
                     phi_lo(KS,i,j), qa_lo, qb_lo  )
                qmin = min( &
                     phi_in(KS,i,j), qa_in, qb_in, &
                     phi_lo(KS,i,j), qa_lo, qb_lo  )
                qjpls(KS,i,j) = ( qmax - phi_lo(KS,i,j) ) * DENS(KS,i,j)
                qjmns(KS,i,j) = ( phi_lo(KS,i,j) - qmin ) * DENS(KS,i,j)
             end do

             ! k = KE
             !$omp parallel do private(j,rw,ru,rv,fact,qa_in,qb_in,qa_lo,qb_lo,qmax,qmin) OMP_SCHEDULE_
             do j = JJS, JJE
                rw = (                     mflx_hi(KE-1,i  ,j  ,ZDIR)) * RDZ(KE)! 2 * rho * w / dz
                ru = 0.0_RP
                rv = (mflx_hi(KE,i,j,YDIR)+mflx_hi(KE  ,i  ,j-1,YDIR)) * RDY(j) ! 2 * rho * v / dy

                call get_fact_fct( fact, & ! (out)
                                   rw, ru, rv ) ! (in)

                qa_in = fact(1, 0, 1) * phi_in(KE  ,i  ,j+1) &
                      + fact(0, 0, 1) * phi_in(KE  ,i  ,j+1) &
                      + fact(1, 0, 0) * phi_in(KE  ,i  ,j  ) &
                      + fact(1, 0,-1) * phi_in(KE  ,i  ,j-1) &
                      + fact(0, 0, 0) * phi_in(KE  ,i  ,j  )
                qb_in = fact(1, 0, 1) * phi_in(KE-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_in(KE  ,i  ,j-1) &
                      + fact(1, 0, 0) * phi_in(KE-1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_in(KE-1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_in(KE  ,i  ,j  )
                qa_lo = fact(1, 0, 1) * phi_lo(KE  ,i  ,j+1) &
                      + fact(0, 0, 1) * phi_lo(KE  ,i  ,j+1) &
                      + fact(1, 0, 0) * phi_lo(KE  ,i  ,j  ) &
                      + fact(1, 0,-1) * phi_lo(KE  ,i  ,j-1) &
                      + fact(0, 0, 0) * phi_lo(KE  ,i  ,j  )
                qb_lo = fact(1, 0, 1) * phi_lo(KE-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_lo(KE  ,i  ,j-1) &
                      + fact(1, 0, 0) * phi_lo(KE-1,i  ,j  ) &
                      + fact(1, 0,-1) * phi_lo(KE-1,i  ,j-1) &
                      + fact(0, 0, 0) * phi_lo(KE  ,i  ,j  )

                qmax = max( &
                     phi_in(KE,i,j), qa_in, qb_in, &
                     phi_lo(KE,i,j), qa_lo, qb_lo  )
                qmin = min( &
                     phi_in(KE,i,j), qa_in, qb_in, &
                     phi_lo(KE,i,j), qa_lo, qb_lo  )
                qjpls(KE,i,j) = ( qmax - phi_lo(KE,i,j) ) * DENS(KE,i,j)
                qjmns(KE,i,j) = ( phi_lo(KE,i,j) - qmin ) * DENS(KE,i,j)
             end do

          else
             !$omp parallel do private(i,j,k,rw,ru,rv,fact,qa_in,qb_in,qa_lo,qb_lo,qmax,qmin) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS+1, KE-1
                rw = (mflx_hi(k,i,j,ZDIR)+mflx_hi(k-1,i  ,j  ,ZDIR)) * RDZ(k) ! 2 * rho * w / dz
                ru = (mflx_hi(k,i,j,XDIR)+mflx_hi(k  ,i-1,j  ,XDIR)) * RDX(i) ! 2 * rho * u / dx
                rv = (mflx_hi(k,i,j,YDIR)+mflx_hi(k  ,i  ,j-1,YDIR)) * RDY(j) ! 2 * rho * v / dy

                call get_fact_fct( fact, & ! (out)
                                   rw, ru, rv ) ! (in)

                qa_in = fact(1, 1, 1) * phi_in(k+1,i+1,j+1) &
                      + fact(0, 1, 1) * phi_in(k  ,i+1,j+1) &
                      + fact(1, 0, 1) * phi_in(k+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_in(k  ,i  ,j+1) &
                      + fact(1,-1, 1) * phi_in(k+1,i-1,j+1) &
                      + fact(1, 1, 0) * phi_in(k+1,i+1,j  ) &
                      + fact(0, 1, 0) * phi_in(k  ,i+1,j  ) &
                      + fact(1, 0, 0) * phi_in(k+1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_in(k+1,i-1,j  ) &
                      + fact(1, 1,-1) * phi_in(k+1,i+1,j-1) &
                      + fact(0, 1,-1) * phi_in(k  ,i+1,j-1) &
                      + fact(1, 0,-1) * phi_in(k+1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_in(k+1,i-1,j-1) &
                      + fact(0, 0, 0) * phi_in(k  ,i  ,j  )
                qb_in = fact(1, 1, 1) * phi_in(k-1,i-1,j-1) &
                      + fact(0, 1, 1) * phi_in(k  ,i-1,j-1) &
                      + fact(1, 0, 1) * phi_in(k-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_in(k  ,i  ,j-1) &
                      + fact(1,-1, 1) * phi_in(k-1,i+1,j-1) &
                      + fact(1, 1, 0) * phi_in(k-1,i-1,j  ) &
                      + fact(0, 1, 0) * phi_in(k  ,i-1,j  ) &
                      + fact(1, 0, 0) * phi_in(k-1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_in(k-1,i+1,j  ) &
                      + fact(1, 1,-1) * phi_in(k-1,i-1,j+1) &
                      + fact(0, 1,-1) * phi_in(k  ,i-1,j-1) &
                      + fact(1, 0,-1) * phi_in(k-1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_in(k-1,i+1,j+1) &
                      + fact(0, 0, 0) * phi_in(k  ,i  ,j  )
                qa_lo = fact(1, 1, 1) * phi_lo(k+1,i+1,j+1) &
                      + fact(0, 1, 1) * phi_lo(k  ,i+1,j+1) &
                      + fact(1, 0, 1) * phi_lo(k+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_lo(k  ,i  ,j+1) &
                      + fact(1,-1, 1) * phi_lo(k+1,i-1,j+1) &
                      + fact(1, 1, 0) * phi_lo(k+1,i+1,j  ) &
                      + fact(0, 1, 0) * phi_lo(k  ,i+1,j  ) &
                      + fact(1, 0, 0) * phi_lo(k+1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_lo(k+1,i-1,j  ) &
                      + fact(1, 1,-1) * phi_lo(k+1,i+1,j-1) &
                      + fact(0, 1,-1) * phi_lo(k  ,i+1,j-1) &
                      + fact(1, 0,-1) * phi_lo(k+1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_lo(k+1,i-1,j-1) &
                      + fact(0, 0, 0) * phi_lo(k  ,i  ,j  )
                qb_lo = fact(1, 1, 1) * phi_lo(k-1,i-1,j-1) &
                      + fact(0, 1, 1) * phi_lo(k  ,i-1,j-1) &
                      + fact(1, 0, 1) * phi_lo(k-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_lo(k  ,i  ,j-1) &
                      + fact(1,-1, 1) * phi_lo(k-1,i+1,j-1) &
                      + fact(1, 1, 0) * phi_lo(k-1,i-1,j  ) &
                      + fact(0, 1, 0) * phi_lo(k  ,i-1,j  ) &
                      + fact(1, 0, 0) * phi_lo(k-1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_lo(k-1,i+1,j  ) &
                      + fact(1, 1,-1) * phi_lo(k-1,i-1,j+1) &
                      + fact(0, 1,-1) * phi_lo(k  ,i-1,j-1) &
                      + fact(1, 0,-1) * phi_lo(k-1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_lo(k-1,i+1,j+1) &
                      + fact(0, 0, 0) * phi_lo(k  ,i  ,j  )

                qmax = max( &
                     phi_in(k,i,j), qa_in, qb_in, &
                     phi_lo(k,i,j), qa_lo, qb_lo  )
                qmin = min( &
                     phi_in(k,i,j), qa_in, qb_in, &
                     phi_lo(k,i,j), qa_lo, qb_lo  )
                qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
                qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
             end do
             end do
             end do

             ! k = KS
             !$omp parallel do private(i,j,rw,ru,rv,fact,qa_in,qb_in,qa_lo,qb_lo,qmax,qmin) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
                rw = (mflx_hi(KS,i,j,ZDIR)                           ) * RDZ(KS)! 2 * rho * w / dz
                ru = (mflx_hi(KS,i,j,XDIR)+mflx_hi(KS  ,i-1,j  ,XDIR)) * RDX(i) ! 2 * rho * u / dx
                rv = (mflx_hi(KS,i,j,YDIR)+mflx_hi(KS  ,i  ,j-1,YDIR)) * RDY(j) ! 2 * rho * v / dy

                call get_fact_fct( fact, & ! (out)
                                   rw, ru, rv ) ! (in)

                qa_in = fact(1, 1, 1) * phi_in(KS+1,i+1,j+1) &
                      + fact(0, 1, 1) * phi_in(KS  ,i+1,j+1) &
                      + fact(1, 0, 1) * phi_in(KS+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_in(KS  ,i  ,j+1) &
                      + fact(1,-1, 1) * phi_in(KS+1,i-1,j+1) &
                      + fact(1, 1, 0) * phi_in(KS+1,i+1,j  ) &
                      + fact(0, 1, 0) * phi_in(KS  ,i+1,j  ) &
                      + fact(1, 0, 0) * phi_in(KS+1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_in(KS+1,i-1,j  ) &
                      + fact(1, 1,-1) * phi_in(KS+1,i+1,j-1) &
                      + fact(0, 1,-1) * phi_in(KS  ,i+1,j-1) &
                      + fact(1, 0,-1) * phi_in(KS+1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_in(KS+1,i-1,j-1) &
                      + fact(0, 0, 0) * phi_in(KS  ,i  ,j  )
                qb_in = fact(1, 1, 1) * phi_in(KS  ,i-1,j-1) &
                      + fact(0, 1, 1) * phi_in(KS  ,i-1,j-1) &
                      + fact(1, 0, 1) * phi_in(KS  ,i  ,j-1) &
                      + fact(0, 0, 1) * phi_in(KS  ,i  ,j-1) &
                      + fact(1,-1, 1) * phi_in(KS  ,i+1,j-1) &
                      + fact(1, 1, 0) * phi_in(KS  ,i-1,j  ) &
                      + fact(0, 1, 0) * phi_in(KS  ,i-1,j  ) &
                      + fact(1, 0, 0) * phi_in(KS  ,i  ,j  ) &
                      + fact(1,-1, 0) * phi_in(KS  ,i+1,j  ) &
                      + fact(1, 1,-1) * phi_in(KS  ,i-1,j+1) &
                      + fact(0, 1,-1) * phi_in(KS  ,i-1,j-1) &
                      + fact(1, 0,-1) * phi_in(KS  ,i  ,j-1) &
                      + fact(1,-1,-1) * phi_in(KS  ,i+1,j+1) &
                      + fact(0, 0, 0) * phi_in(KS  ,i  ,j  )
                qa_lo = fact(1, 1, 1) * phi_lo(KS+1,i+1,j+1) &
                      + fact(0, 1, 1) * phi_lo(KS  ,i+1,j+1) &
                      + fact(1, 0, 1) * phi_lo(KS+1,i  ,j+1) &
                      + fact(0, 0, 1) * phi_lo(KS  ,i  ,j+1) &
                      + fact(1,-1, 1) * phi_lo(KS+1,i-1,j+1) &
                      + fact(1, 1, 0) * phi_lo(KS+1,i+1,j  ) &
                      + fact(0, 1, 0) * phi_lo(KS  ,i+1,j  ) &
                      + fact(1, 0, 0) * phi_lo(KS+1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_lo(KS+1,i-1,j  ) &
                      + fact(1, 1,-1) * phi_lo(KS+1,i+1,j-1) &
                      + fact(0, 1,-1) * phi_lo(KS  ,i+1,j-1) &
                      + fact(1, 0,-1) * phi_lo(KS+1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_lo(KS+1,i-1,j-1) &
                      + fact(0, 0, 0) * phi_lo(KS  ,i  ,j  )
                qb_lo = fact(1, 1, 1) * phi_lo(KS  ,i-1,j-1) &
                      + fact(0, 1, 1) * phi_lo(KS  ,i-1,j-1) &
                      + fact(1, 0, 1) * phi_lo(KS  ,i  ,j-1) &
                      + fact(0, 0, 1) * phi_lo(KS  ,i  ,j-1) &
                      + fact(1,-1, 1) * phi_lo(KS  ,i+1,j-1) &
                      + fact(1, 1, 0) * phi_lo(KS  ,i-1,j  ) &
                      + fact(0, 1, 0) * phi_lo(KS  ,i-1,j  ) &
                      + fact(1, 0, 0) * phi_lo(KS  ,i  ,j  ) &
                      + fact(1,-1, 0) * phi_lo(KS  ,i+1,j  ) &
                      + fact(1, 1,-1) * phi_lo(KS  ,i-1,j+1) &
                      + fact(0, 1,-1) * phi_lo(KS  ,i-1,j-1) &
                      + fact(1, 0,-1) * phi_lo(KS  ,i  ,j-1) &
                      + fact(1,-1,-1) * phi_lo(KS  ,i+1,j+1) &
                      + fact(0, 0, 0) * phi_lo(KS  ,i  ,j  )

                qmax = max( &
                     phi_in(KS,i,j), qa_in, qb_in, &
                     phi_lo(KS,i,j), qa_lo, qb_lo  )
                qmin = min( &
                     phi_in(KS,i,j), qa_in, qb_in, &
                     phi_lo(KS,i,j), qa_lo, qb_lo  )
                qjpls(KS,i,j) = ( qmax - phi_lo(KS,i,j) ) * DENS(KS,i,j)
                qjmns(KS,i,j) = ( phi_lo(KS,i,j) - qmin ) * DENS(KS,i,j)
             end do
             end do

             ! k = KE
             !$omp parallel do private(i,j,rw,ru,rv,fact,qa_in,qb_in,qa_lo,qb_lo,qmax,qmin) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
                rw = (                     mflx_hi(KE-1,i  ,j  ,ZDIR)) * RDZ(KE)! 2 * rho * w / dz
                ru = (mflx_hi(KE,i,j,XDIR)+mflx_hi(KE  ,i-1,j  ,XDIR)) * RDX(i) ! 2 * rho * u / dx
                rv = (mflx_hi(KE,i,j,YDIR)+mflx_hi(KE  ,i  ,j-1,YDIR)) * RDY(j) ! 2 * rho * v / dy

                call get_fact_fct( fact, & ! (out)
                                   rw, ru, rv ) ! (in)

                qa_in = fact(1, 1, 1) * phi_in(KE  ,i+1,j+1) &
                      + fact(0, 1, 1) * phi_in(KE  ,i+1,j+1) &
                      + fact(1, 0, 1) * phi_in(KE  ,i  ,j+1) &
                      + fact(0, 0, 1) * phi_in(KE  ,i  ,j+1) &
                      + fact(1,-1, 1) * phi_in(KE  ,i-1,j+1) &
                      + fact(1, 1, 0) * phi_in(KE  ,i+1,j  ) &
                      + fact(0, 1, 0) * phi_in(KE  ,i+1,j  ) &
                      + fact(1, 0, 0) * phi_in(KE  ,i  ,j  ) &
                      + fact(1,-1, 0) * phi_in(KE  ,i-1,j  ) &
                      + fact(1, 1,-1) * phi_in(KE  ,i+1,j-1) &
                      + fact(0, 1,-1) * phi_in(KE  ,i+1,j-1) &
                      + fact(1, 0,-1) * phi_in(KE  ,i  ,j-1) &
                      + fact(1,-1,-1) * phi_in(KE  ,i-1,j-1) &
                      + fact(0, 0, 0) * phi_in(KE  ,i  ,j  )
                qb_in = fact(1, 1, 1) * phi_in(KE-1,i-1,j-1) &
                      + fact(0, 1, 1) * phi_in(KE  ,i-1,j-1) &
                      + fact(1, 0, 1) * phi_in(KE-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_in(KE  ,i  ,j-1) &
                      + fact(1,-1, 1) * phi_in(KE-1,i+1,j-1) &
                      + fact(1, 1, 0) * phi_in(KE-1,i-1,j  ) &
                      + fact(0, 1, 0) * phi_in(KE  ,i-1,j  ) &
                      + fact(1, 0, 0) * phi_in(KE-1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_in(KE-1,i+1,j  ) &
                      + fact(1, 1,-1) * phi_in(KE-1,i-1,j+1) &
                      + fact(0, 1,-1) * phi_in(KE  ,i-1,j-1) &
                      + fact(1, 0,-1) * phi_in(KE-1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_in(KE-1,i+1,j+1) &
                      + fact(0, 0, 0) * phi_in(KE  ,i  ,j  )
                qa_lo = fact(1, 1, 1) * phi_lo(KE  ,i+1,j+1) &
                      + fact(0, 1, 1) * phi_lo(KE  ,i+1,j+1) &
                      + fact(1, 0, 1) * phi_lo(KE  ,i  ,j+1) &
                      + fact(0, 0, 1) * phi_lo(KE  ,i  ,j+1) &
                      + fact(1,-1, 1) * phi_lo(KE  ,i-1,j+1) &
                      + fact(1, 1, 0) * phi_lo(KE  ,i+1,j  ) &
                      + fact(0, 1, 0) * phi_lo(KE  ,i+1,j  ) &
                      + fact(1, 0, 0) * phi_lo(KE  ,i  ,j  ) &
                      + fact(1,-1, 0) * phi_lo(KE  ,i-1,j  ) &
                      + fact(1, 1,-1) * phi_lo(KE  ,i+1,j-1) &
                      + fact(0, 1,-1) * phi_lo(KE  ,i+1,j-1) &
                      + fact(1, 0,-1) * phi_lo(KE  ,i  ,j-1) &
                      + fact(1,-1,-1) * phi_lo(KE  ,i-1,j-1) &
                      + fact(0, 0, 0) * phi_lo(KE  ,i  ,j  )
                qb_lo = fact(1, 1, 1) * phi_lo(KE-1,i-1,j-1) &
                      + fact(0, 1, 1) * phi_lo(KE  ,i-1,j-1) &
                      + fact(1, 0, 1) * phi_lo(KE-1,i  ,j-1) &
                      + fact(0, 0, 1) * phi_lo(KE  ,i  ,j-1) &
                      + fact(1,-1, 1) * phi_lo(KE-1,i+1,j-1) &
                      + fact(1, 1, 0) * phi_lo(KE-1,i-1,j  ) &
                      + fact(0, 1, 0) * phi_lo(KE  ,i-1,j  ) &
                      + fact(1, 0, 0) * phi_lo(KE-1,i  ,j  ) &
                      + fact(1,-1, 0) * phi_lo(KE-1,i+1,j  ) &
                      + fact(1, 1,-1) * phi_lo(KE-1,i-1,j+1) &
                      + fact(0, 1,-1) * phi_lo(KE  ,i-1,j-1) &
                      + fact(1, 0,-1) * phi_lo(KE-1,i  ,j-1) &
                      + fact(1,-1,-1) * phi_lo(KE-1,i+1,j+1) &
                      + fact(0, 0, 0) * phi_lo(KE  ,i  ,j  )

                qmax = max( &
                     phi_in(KE,i,j), qa_in, qb_in, &
                     phi_lo(KE,i,j), qa_lo, qb_lo  )
                qmin = min( &
                     phi_in(KE,i,j), qa_in, qb_in, &
                     phi_lo(KE,i,j), qa_lo, qb_lo  )
                qjpls(KE,i,j) = ( qmax - phi_lo(KE,i,j) ) * DENS(KE,i,j)
                qjmns(KE,i,j) = ( phi_lo(KE,i,j) - qmin ) * DENS(KE,i,j)
             end do
             end do
          end if

       else

          if ( twoD ) then
             i = IS
             !$omp parallel do private(j,k,qmax,qmin) OMP_SCHEDULE_
             do j = JJS, JJE
             do k = KS+1, KE-1
#ifdef DEBUG
                call CHECK( __LINE__, phi_in(k  ,i  ,j  ) )
                call CHECK( __LINE__, phi_in(k-1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(k+1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(k  ,i  ,j+1) )
                call CHECK( __LINE__, phi_in(k  ,i  ,j-1) )
                call CHECK( __LINE__, phi_lo(k  ,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(k-1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(k+1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(k  ,i  ,j+1) )
                call CHECK( __LINE__, phi_lo(k  ,i  ,j-1) )
#endif
                qmax = max( phi_in(k  ,i  ,j  ), &
                            phi_in(k+1,i  ,j  ), &
                            phi_in(k-1,i  ,j  ), &
                            phi_in(k  ,i  ,j+1), &
                            phi_in(k  ,i  ,j-1), &
                            phi_lo(k  ,i  ,j  ), &
                            phi_lo(k+1,i  ,j  ), &
                            phi_lo(k-1,i  ,j  ), &
                            phi_lo(k  ,i  ,j+1), &
                            phi_lo(k  ,i  ,j-1) )
                qmin = min( phi_in(k  ,i  ,j  ), &
                            phi_in(k+1,i  ,j  ), &
                            phi_in(k-1,i  ,j  ), &
                            phi_in(k  ,i  ,j+1), &
                            phi_in(k  ,i  ,j-1), &
                            phi_lo(k  ,i  ,j  ), &
                            phi_lo(k+1,i  ,j  ), &
                            phi_lo(k-1,i  ,j  ), &
                            phi_lo(k  ,i  ,j+1), &
                            phi_lo(k  ,i  ,j-1) )
                qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
                qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
             enddo
             enddo
#ifdef DEBUG
             k = IUNDEF; j = IUNDEF
#endif
             !$omp parallel do private(i,j,qmax,qmin) OMP_SCHEDULE_
             do j = JJS, JJE
#ifdef DEBUG
                call CHECK( __LINE__, phi_in(KS  ,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KS+1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KS  ,i  ,j+1) )
                call CHECK( __LINE__, phi_in(KS  ,i  ,j-1) )
                call CHECK( __LINE__, phi_lo(KS  ,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KS+1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KS  ,i  ,j+1) )
                call CHECK( __LINE__, phi_lo(KS  ,i  ,j-1) )
                call CHECK( __LINE__, phi_in(KE  ,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KE-1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KE  ,i  ,j+1) )
                call CHECK( __LINE__, phi_in(KE  ,i  ,j-1) )
                call CHECK( __LINE__, phi_lo(KE  ,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KE-1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KE  ,i  ,j+1) )
                call CHECK( __LINE__, phi_lo(KE  ,i  ,j-1) )
#endif
                qmax = max( phi_in(KS  ,i  ,j  ), &
                            phi_in(KS+1,i  ,j  ), &
                            phi_in(KS  ,i  ,j+1), &
                            phi_in(KS  ,i  ,j-1), &
                            phi_lo(KS  ,i  ,j  ), &
                            phi_lo(KS+1,i  ,j  ), &
                            phi_lo(KS  ,i  ,j+1), &
                            phi_lo(KS  ,i  ,j-1) )
                qmin = min( phi_in(KS  ,i  ,j  ), &
                            phi_in(KS+1,i  ,j  ), &
                            phi_in(KS  ,i  ,j+1), &
                            phi_in(KS  ,i  ,j-1), &
                            phi_lo(KS  ,i  ,j  ), &
                            phi_lo(KS+1,i  ,j  ), &
                            phi_lo(KS  ,i  ,j+1), &
                            phi_lo(KS  ,i  ,j-1) )
                qjmns(KS,i,j) = ( phi_lo(KS,i,j) - qmin ) * DENS(KS,i,j)
                qjpls(KS,i,j) = ( qmax - phi_lo(KS,i,j) ) * DENS(KS,i,j)

                qmax = max( phi_in(KE  ,i  ,j  ), &
                            phi_in(KE-1,i  ,j  ), &
                            phi_in(KE  ,i  ,j+1), &
                            phi_in(KE  ,i  ,j-1), &
                            phi_lo(KE  ,i  ,j  ), &
                            phi_lo(KE-1,i  ,j  ), &
                            phi_lo(KE  ,i  ,j+1), &
                            phi_lo(KE  ,i  ,j-1) )
                qmin = min( phi_in(KE  ,i  ,j  ), &
                            phi_in(KE-1,i  ,j  ), &
                            phi_in(KE  ,i  ,j+1), &
                            phi_in(KE  ,i  ,j-1), &
                            phi_lo(KE  ,i  ,j  ), &
                            phi_lo(KE-1,i  ,j  ), &
                            phi_lo(KE  ,i  ,j+1), &
                            phi_lo(KE  ,i  ,j-1) )
                qjpls(KE,i,j) = ( qmax - phi_lo(KE,i,j) ) * DENS(KE,i,j)
                qjmns(KE,i,j) = ( phi_lo(KE,i,j) - qmin ) * DENS(KE,i,j)
             enddo
#ifdef DEBUG
             k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          else
             !$omp parallel do private(i,j,k,qmax,qmin) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS+1, KE-1
#ifdef DEBUG
                call CHECK( __LINE__, phi_in(k  ,i  ,j  ) )
                call CHECK( __LINE__, phi_in(k-1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(k+1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(k  ,i-1,j  ) )
                call CHECK( __LINE__, phi_in(k  ,i+1,j  ) )
                call CHECK( __LINE__, phi_in(k  ,i  ,j+1) )
                call CHECK( __LINE__, phi_in(k  ,i  ,j-1) )
                call CHECK( __LINE__, phi_lo(k  ,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(k-1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(k+1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(k  ,i-1,j  ) )
                call CHECK( __LINE__, phi_lo(k  ,i+1,j  ) )
                call CHECK( __LINE__, phi_lo(k  ,i  ,j+1) )
                call CHECK( __LINE__, phi_lo(k  ,i  ,j-1) )
#endif
                qmax = max( phi_in(k  ,i  ,j  ), &
                            phi_in(k+1,i  ,j  ), &
                            phi_in(k-1,i  ,j  ), &
                            phi_in(k  ,i+1,j  ), &
                            phi_in(k  ,i-1,j  ), &
                            phi_in(k  ,i  ,j+1), &
                            phi_in(k  ,i  ,j-1), &
                            phi_lo(k  ,i  ,j  ), &
                            phi_lo(k+1,i  ,j  ), &
                            phi_lo(k-1,i  ,j  ), &
                            phi_lo(k  ,i+1,j  ), &
                            phi_lo(k  ,i-1,j  ), &
                            phi_lo(k  ,i  ,j+1), &
                            phi_lo(k  ,i  ,j-1) )
                qmin = min( phi_in(k  ,i  ,j  ), &
                            phi_in(k+1,i  ,j  ), &
                            phi_in(k-1,i  ,j  ), &
                            phi_in(k  ,i-1,j  ), &
                            phi_in(k  ,i+1,j  ), &
                            phi_in(k  ,i  ,j+1), &
                            phi_in(k  ,i  ,j-1), &
                            phi_lo(k  ,i  ,j  ), &
                            phi_lo(k+1,i  ,j  ), &
                            phi_lo(k-1,i  ,j  ), &
                            phi_lo(k  ,i-1,j  ), &
                            phi_lo(k  ,i+1,j  ), &
                            phi_lo(k  ,i  ,j+1), &
                            phi_lo(k  ,i  ,j-1) )
                qjpls(k,i,j) = ( qmax - phi_lo(k,i,j) ) * DENS(k,i,j)
                qjmns(k,i,j) = ( phi_lo(k,i,j) - qmin ) * DENS(k,i,j)
             enddo
             enddo
             enddo
#ifdef DEBUG
             k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
             !$omp parallel do private(i,j,qmax,qmin) OMP_SCHEDULE_ collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
#ifdef DEBUG
                call CHECK( __LINE__, phi_in(KS  ,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KS+1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KS  ,i-1,j  ) )
                call CHECK( __LINE__, phi_in(KS  ,i+1,j  ) )
                call CHECK( __LINE__, phi_in(KS  ,i  ,j+1) )
                call CHECK( __LINE__, phi_in(KS  ,i  ,j-1) )
                call CHECK( __LINE__, phi_lo(KS  ,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KS+1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KS  ,i-1,j  ) )
                call CHECK( __LINE__, phi_lo(KS  ,i+1,j  ) )
                call CHECK( __LINE__, phi_lo(KS  ,i  ,j+1) )
                call CHECK( __LINE__, phi_lo(KS  ,i  ,j-1) )
                call CHECK( __LINE__, phi_in(KE  ,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KE-1,i  ,j  ) )
                call CHECK( __LINE__, phi_in(KE  ,i-1,j  ) )
                call CHECK( __LINE__, phi_in(KE  ,i+1,j  ) )
                call CHECK( __LINE__, phi_in(KE  ,i  ,j+1) )
                call CHECK( __LINE__, phi_in(KE  ,i  ,j-1) )
                call CHECK( __LINE__, phi_lo(KE  ,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KE-1,i  ,j  ) )
                call CHECK( __LINE__, phi_lo(KE  ,i-1,j  ) )
                call CHECK( __LINE__, phi_lo(KE  ,i+1,j  ) )
                call CHECK( __LINE__, phi_lo(KE  ,i  ,j+1) )
                call CHECK( __LINE__, phi_lo(KE  ,i  ,j-1) )
#endif
                qmax = max( phi_in(KS  ,i  ,j  ), &
                            phi_in(KS+1,i  ,j  ), &
                            phi_in(KS  ,i+1,j  ), &
                            phi_in(KS  ,i-1,j  ), &
                            phi_in(KS  ,i  ,j+1), &
                            phi_in(KS  ,i  ,j-1), &
                            phi_lo(KS  ,i  ,j  ), &
                            phi_lo(KS+1,i  ,j  ), &
                            phi_lo(KS  ,i+1,j  ), &
                            phi_lo(KS  ,i-1,j  ), &
                            phi_lo(KS  ,i  ,j+1), &
                            phi_lo(KS  ,i  ,j-1) )
                qmin = min( phi_in(KS  ,i  ,j  ), &
                            phi_in(KS+1,i  ,j  ), &
                            phi_in(KS  ,i+1,j  ), &
                            phi_in(KS  ,i-1,j  ), &
                            phi_in(KS  ,i  ,j+1), &
                            phi_in(KS  ,i  ,j-1), &
                            phi_lo(KS  ,i  ,j  ), &
                            phi_lo(KS+1,i  ,j  ), &
                            phi_lo(KS  ,i+1,j  ), &
                            phi_lo(KS  ,i-1,j  ), &
                            phi_lo(KS  ,i  ,j+1), &
                            phi_lo(KS  ,i  ,j-1) )
                qjmns(KS,i,j) = ( phi_lo(KS,i,j) - qmin ) * DENS(KS,i,j)
                qjpls(KS,i,j) = ( qmax - phi_lo(KS,i,j) ) * DENS(KS,i,j)

                qmax = max( phi_in(KE  ,i  ,j  ), &
                            phi_in(KE-1,i  ,j  ), &
                            phi_in(KE  ,i+1,j  ), &
                            phi_in(KE  ,i-1,j  ), &
                            phi_in(KE  ,i  ,j+1), &
                            phi_in(KE  ,i  ,j-1), &
                            phi_lo(KE  ,i  ,j  ), &
                            phi_lo(KE-1,i  ,j  ), &
                            phi_lo(KE  ,i+1,j  ), &
                            phi_lo(KE  ,i-1,j  ), &
                            phi_lo(KE  ,i  ,j+1), &
                            phi_lo(KE  ,i  ,j-1) )
                qmin = min( phi_in(KE  ,i  ,j  ), &
                            phi_in(KE-1,i  ,j  ), &
                            phi_in(KE  ,i-1,j  ), &
                            phi_in(KE  ,i+1,j  ), &
                            phi_in(KE  ,i  ,j+1), &
                            phi_in(KE  ,i  ,j-1), &
                            phi_lo(KE  ,i  ,j  ), &
                            phi_lo(KE-1,i  ,j  ), &
                            phi_lo(KE  ,i-1,j  ), &
                            phi_lo(KE  ,i+1,j  ), &
                            phi_lo(KE  ,i  ,j+1), &
                            phi_lo(KE  ,i  ,j-1) )
                qjpls(KE,i,j) = ( qmax - phi_lo(KE,i,j) ) * DENS(KE,i,j)
                qjmns(KE,i,j) = ( phi_lo(KE,i,j) - qmin ) * DENS(KE,i,j)
             enddo
             enddo
#ifdef DEBUG
             k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          end if

       end if

       !--- incoming flux limitation factor (0-1)
       !$omp parallel do private(i,j,k,zerosw) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, pjpls(k,i,j) )
          call CHECK( __LINE__, qjpls(k,i,j) )
#endif
          ! if pjpls == 0, zerosw = 1 and rjpls = 0
          zerosw = 0.5_RP - sign( 0.5_RP, pjpls(k,i,j)-EPSILON )
          rjpls(k,i,j) = min( 1.0_RP, qjpls(k,i,j) * ( 1.0_RP-zerosw ) / ( pjpls(k,i,j)-zerosw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- outgoing flux limitation factor (0-1)
       !$omp parallel do private(i,j,k,zerosw) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, pjmns(k,i,j) )
          call CHECK( __LINE__, qjmns(k,i,j) )
#endif
          ! if pjmns == 0, zerosw = 1 and rjmns = 0
          zerosw = 0.5_RP - sign( 0.5_RP, pjmns(k,i,j)-EPSILON )
          rjmns(k,i,j) = min( 1.0_RP, qjmns(k,i,j) * ( 1.0_RP-zerosw ) / ( pjmns(k,i,j)-zerosw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    call COMM_vars8( rjpls(:,:,:), 1 )
    call COMM_vars8( rjmns(:,:,:), 2 )
    call COMM_wait ( rjpls(:,:,:), 1 )
    call COMM_wait ( rjmns(:,:,:), 2 )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !--- update high order flux with antidiffusive flux
       !$omp parallel do private(i,j,k,dirsw) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS , KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,ZDIR) )
          call CHECK( __LINE__, rjpls(k  ,i,j) )
          call CHECK( __LINE__, rjpls(k+1,i,j) )
          call CHECK( __LINE__, rjmns(k  ,i,j) )
          call CHECK( __LINE__, rjmns(k+1,i,j) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,ZDIR) )
          qflx_anti(k,i,j,ZDIR) = qflx_anti(k,i,j,ZDIR) &
                 * ( 1.0_RP &
                   - min( rjpls(k+1,i,j),rjmns(k  ,i,j) ) * (          dirsw ) &
                   - min( rjpls(k  ,i,j),rjmns(k+1,i,j) ) * ( 1.0_RP - dirsw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(KE,i,j,ZDIR) )
          call CHECK( __LINE__, rjpls(KE  ,i,j) )
          call CHECK( __LINE__, rjmns(KE  ,i,j) )
#endif
          qflx_anti(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_anti(KE  ,i,j,ZDIR) = 0.0_RP ! top    boundary
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       if ( .not. TwoD ) then
          if ( IIS == IS ) then
             ijs = IIS-1
          else
             ijs = IIS
          end if

          !$omp parallel do private(i,j,k,dirsw) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = ijs, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, qflx_anti(k,i,j,XDIR) )
             call CHECK( __LINE__, rjpls(k,i  ,j) )
             call CHECK( __LINE__, rjpls(k,i+1,j) )
             call CHECK( __LINE__, rjmns(k,i  ,j) )
             call CHECK( __LINE__, rjmns(k,i+1,j) )
#endif
             ! if qflx_anti > 0, dirsw = 1
             dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,XDIR) )
             qflx_anti(k,i,j,XDIR) = qflx_anti(k,i,j,XDIR) &
                 * ( 1.0_RP &
                   - min( rjpls(k,i+1,j),rjmns(k,i  ,j) ) * (          dirsw ) &
                   - min( rjpls(k,i  ,j),rjmns(k,i+1,j) ) * ( 1.0_RP - dirsw ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if

       if ( JJS == JS ) then
          ijs = JJS-1
       else
          ijs = JJS
       end if
       !$omp parallel do private(i,j,k,dirsw) OMP_SCHEDULE_ collapse(2)
       do j = ijs, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,YDIR) )
          call CHECK( __LINE__, rjpls(k,i,j+1) )
          call CHECK( __LINE__, rjpls(k,i,j  ) )
          call CHECK( __LINE__, rjmns(k,i,j  ) )
          call CHECK( __LINE__, rjmns(k,i,j+1) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,YDIR) )
          qflx_anti(k,i,j,YDIR) = qflx_anti(k,i,j,YDIR) &
                 * ( 1.0_RP &
                   - min( rjpls(k,i,j+1),rjmns(k,i,j  ) ) * (          dirsw ) &
                   - min( rjpls(k,i,j  ),rjmns(k,i,j+1) ) * ( 1.0_RP - dirsw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    return
  end subroutine ATMOS_DYN_fvm_fct

  !-----------------------------------------------------------------------------
  ! private procedure
  ! get factor for FCT
  subroutine get_fact_fct( &
       fact, &
       rw, ru, rv )
    use scale_const, only: &
       EPSILON => CONST_EPS
    implicit none
    real(RP), intent(out) :: fact(0:1,-1:1,-1:1)
    real(RP), intent(in)  :: rw, ru, rv

    real(RP) :: sign_uv, sign_uw, sign_vw ! uv>=0, uw>=0, vw>=0
    real(RP) :: ugev, ugew, vgew          ! u>=v, u>=w, u>=w
    real(RP) :: umax, vmax, wmax
    real(RP) :: vu, wu, uv, wv, uw, vw    ! |v/u|, |w/u|, ....
    real(RP) :: uzero, vzero, wzero
    !---------------------------------------------------------------------------

    ugev = sign(0.5_RP, abs(ru)-abs(rv)) + 0.5_RP ! u >= v
    ugew = sign(0.5_RP, abs(ru)-abs(rw)) + 0.5_RP ! u >= w
    vgew = sign(0.5_RP, abs(rv)-abs(rw)) + 0.5_RP ! v >= w

    uzero = sign(0.5_RP,abs(ru)-EPSILON) - 0.5_RP
    vzero = sign(0.5_RP,abs(rv)-EPSILON) - 0.5_RP
    wzero = sign(0.5_RP,abs(rw)-EPSILON) - 0.5_RP

    sign_uv = sign(0.5_RP, ru*rv) + 0.5_RP ! uv >= 0
    sign_uw = sign(0.5_RP, ru*rw) + 0.5_RP ! uw >= 0
    sign_vw = sign(0.5_RP, rv*rw) + 0.5_RP ! vw >= 0

    wu = abs( rw / ( ru+uzero ) * ( 1.0_RP+uzero ) )
    vu = abs( rv / ( ru+uzero ) * ( 1.0_RP+uzero ) )
    uv = abs( ru / ( rv+vzero ) * ( 1.0_RP+vzero ) )
    wv = abs( rw / ( rv+vzero ) * ( 1.0_RP+vzero ) )
    uw = abs( ru / ( rw+wzero ) * ( 1.0_RP+wzero ) )
    vw = abs( rv / ( rw+wzero ) * ( 1.0_RP+wzero ) )

    umax  = ugev * ugew * ( 1.0_RP+uzero ) ! u == max(u,v,w)
    vmax  = (1.0_RP-ugev) * vgew           ! v == max(u,v,w)
    wmax  = 1.0_RP - ugev * ugew - vmax    ! w == max(u,v,w)

    fact(0, 0, 0) = - ugev * ugew * uzero  ! 1.0 if max(u,v,w) < epsilon

    fact(1, 0, 0) = wmax * (1.0_RP-uw) * (1.0_RP-vw)
    fact(0, 1, 0) = umax * (1.0_RP-vu) * (1.0_RP-wu)
    fact(0, 0, 1) = vmax * (1.0_RP-uv) * (1.0_RP-wv)

    fact(1, 1, 1) =         sign_uv  *         sign_uw  * ( umax * vu*wu + vmax * uv*wv + wmax * uw*vw )
    fact(1,-1, 1) = (1.0_RP-sign_uv) * (1.0_RP-sign_uw) * ( umax * vu*wu + vmax * uv*wv + wmax * uw*vw )
    fact(1, 1,-1) = (1.0_RP-sign_uv) *         sign_uw  * ( umax * vu*wu + vmax * uv*wv + wmax * uw*vw )
    fact(1,-1,-1) =         sign_uv  * (1.0_RP-sign_uw) * ( umax * vu*wu + vmax * uv*wv + wmax * uw*vw )

    fact(1, 1, 0) =         sign_uw  * (1.0_RP-vmax) * ( ugew * wu * (1.0_RP-vu) + (1.0_RP-ugew) * uw * (1.0_RP-vw) )
    fact(1,-1, 0) = (1.0_RP-sign_uw) * (1.0_RP-vmax) * ( ugew * wu * (1.0_RP-vu) + (1.0_RP-ugew) * uw * (1.0_RP-vw) )
    fact(1, 0, 1) =         sign_vw  * (1.0_RP-umax) * ( vgew * wv * (1.0_RP-uv) + (1.0_RP-vgew) * vw * (1.0_RP-uw) )
    fact(1, 0,-1) = (1.0_RP-sign_vw) * (1.0_RP-umax) * ( vgew * wv * (1.0_RP-uv) + (1.0_RP-vgew) * vw * (1.0_RP-uw) )
    fact(0, 1, 1) =         sign_uv  * (1.0_RP-wmax) * ( ugev * vu * (1.0_RP-wu) + (1.0_RP-ugev) * uv * (1.0_RP-wv) )
    fact(0, 1,-1) = (1.0_RP-sign_uv) * (1.0_RP-wmax) * ( ugev * vu * (1.0_RP-wu) + (1.0_RP-ugev) * uv * (1.0_RP-wv) )

    return
  end subroutine get_fact_fct

end module scale_atmos_dyn_fvm_fct
