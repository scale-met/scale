!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!          common subroutines for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_common
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
  public :: ATMOS_DYN_wdamp_setup
  public :: ATMOS_DYN_fill_halo
  public :: ATMOS_DYN_Copy_boundary
  public :: ATMOS_DYN_Copy_boundary_tracer
  public :: ATMOS_DYN_divergence
  public :: ATMOS_DYN_prep_pres_linearization


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
  subroutine ATMOS_DYN_wdamp_setup( &
       wdamp_coef,              &
       wdamp_tau, wdamp_height, &
       FZ                       )
    use scale_const, only: &
       PI  => CONST_PI, &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(out) :: wdamp_coef(KA)
    real(RP), intent(in)  :: wdamp_tau
    real(RP), intent(in)  :: wdamp_height
    real(RP), intent(in)  :: FZ(0:KA)

    real(RP) :: alpha, sw

    integer :: k
    !---------------------------------------------------------------------------

    !$acc data copyout(wdamp_coef) copyin(FZ)

    if ( wdamp_height < 0.0_RP ) then
       !$acc kernels
       wdamp_coef(:) = 0.0_RP
       !$acc end kernels
    elseif( FZ(KE)-wdamp_height < EPS ) then
       !$acc kernels
       wdamp_coef(:) = 0.0_RP
       !$acc end kernels
    else
       alpha = 1.0_RP / wdamp_tau

       !$acc kernels
       do k = KS, KE
          sw = 0.5_RP + sign( 0.5_RP, FZ(k)-wdamp_height )

          wdamp_coef(k) = alpha * sw &
                        * 0.5_RP * ( 1.0_RP - cos( PI * (FZ(k)-wdamp_height) / (FZ(KE)-wdamp_height)) )
       enddo
       !$acc end kernels
       !$acc kernels
       wdamp_coef(   1:KS-1) = wdamp_coef(KS)
       !$acc end kernels
       !$acc kernels
       wdamp_coef(KE+1:KA  ) = wdamp_coef(KE)
       !$acc end kernels

       !$acc update host(wdamp_coef)
       LOG_NEWLINE
       LOG_INFO("ATMOS_DYN_wdamp_setup",*)                          'Setup Rayleigh damping coefficient'
       LOG_INFO_CONT('(1x,A)')                   '|=== Rayleigh Damping Coef ===|'
       LOG_INFO_CONT('(1x,A)')                   '|     k     zh[m]    coef[/s] |'
       do k = KA, KE+1, -1
       LOG_INFO_CONT('(1x,A,I5,F10.2,ES12.4,A)') '| ',k, FZ(k), wdamp_coef(k),' |'
       enddo
       k = KE
       LOG_INFO_CONT('(1x,A,I5,F10.2,ES12.4,A)') '| ',k, FZ(k), wdamp_coef(k),' | KE = TOA'
       do k = KE-1, KS, -1
       LOG_INFO_CONT('(1x,A,I5,F10.2,ES12.4,A)') '| ',k, FZ(k), wdamp_coef(k),' |'
       enddo
       k = KS-1
       LOG_INFO_CONT('(1x,A,I5,F10.2,ES12.4,A)') '| ',k, FZ(k), wdamp_coef(k),' | KS-1 = surface'
       do k = KS-2, 1, -1
       LOG_INFO_CONT('(1x,A,I5,F10.2,ES12.4,A)') '| ',k, FZ(k), wdamp_coef(k),' |'
       enddo
       k = 0
       LOG_INFO_CONT('(1x,A,I5,F10.2,12x,A)')    '| ',k, FZ(k),               ' |'
       LOG_INFO_CONT('(1x,A)')                   '|=============================|'
    endif

    !$acc end data

    return
  end subroutine ATMOS_DYN_wdamp_setup

  !-----------------------------------------------------------------------------

  subroutine ATMOS_DYN_fill_halo( var,             &
      fill_constval, lateral_halo, top_bottom_halo )
   implicit none

   real(RP), intent(inout) :: var(KA,IA,JA)
   real(RP), intent(in) :: fill_constval
   logical, intent(in) :: lateral_halo
   logical, intent(in) :: top_bottom_halo

   integer :: i, j, k
   !----------------------------

   !$acc data copy(var)

   if (lateral_halo) then
      !$acc kernels
!OCL XFILL
      do j = 1, JA
      do i = 1, ISB-1
      do k = 1, KA
         var(k,i,j) = fill_constval
      enddo
      enddo
      do i = IEB+1, IA
      do k = 1, KA
         var(k,i,j) = fill_constval
      enddo
      enddo
      enddo
      !$acc end kernels
      !$acc kernels
!OCL XFILL
      do j = 1, JSB-1
      do i = 1, IA
      do k = 1, KA
         var(k,i,j) = fill_constval
      enddo
      enddo
      enddo
      !$acc end kernels
      !$acc kernels
!OCL XFILL
      do j = JEB+1, JA
      do i = 1, IA
      do k = 1, KA
         var(k,i,j) = fill_constval
      enddo
      enddo
      enddo
      !$acc end kernels
   end if

   if (top_bottom_halo) then
      !$acc kernels
!OCL XFILL
      do j = JS, JE
      do i = IS, IE
         var(   1:KS-1,i,j) = fill_constval
         var(KE+1:KA  ,i,j) = fill_constval
      enddo
      enddo
      !$acc end kernels
   end if

   !$acc end data

   return
  end subroutine ATMOS_DYN_fill_halo

  subroutine ATMOS_DYN_Copy_boundary( &
       DENS,  MOMZ,  MOMX,  MOMY,  RHOT,  PROG, &
       DENS0, MOMZ0, MOMX0, MOMY0, RHOT0, PROG0, &
       BND_W, BND_E, BND_S, BND_N, TwoD )
    implicit none
    real(RP), intent(inout) :: DENS (KA,IA,JA)
    real(RP), intent(inout) :: MOMZ (KA,IA,JA)
    real(RP), intent(inout) :: MOMX (KA,IA,JA)
    real(RP), intent(inout) :: MOMY (KA,IA,JA)
    real(RP), intent(inout) :: RHOT (KA,IA,JA)
    real(RP), intent(inout) :: PROG (KA,IA,JA,VA)
    real(RP), intent(in)    :: DENS0(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ0(KA,IA,JA)
    real(RP), intent(in)    :: MOMX0(KA,IA,JA)
    real(RP), intent(in)    :: MOMY0(KA,IA,JA)
    real(RP), intent(in)    :: RHOT0(KA,IA,JA)
    real(RP), intent(in)    :: PROG0(KA,IA,JA,VA)
    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N
    logical,  intent(in)    :: TwoD

    integer :: k, i, j, iv

    !$acc data copy(DENS, MOMZ, MOMX, MOMY, RHOT, PROG) &
    !$acc      copyin(DENS0, MOMZ0, MOMX0, MOMY0, RHOT0, PROG0)

    if ( BND_W .and. (.not. TwoD) ) then
       !$omp parallel do default(none) private(j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,iv) &
       !$omp shared(JA,IS,KS,KE,DENS,DENS0,MOMZ,MOMZ0,MOMX,MOMX0,MOMY,MOMY0,RHOT,RHOT0,VA,PROG,PROG0)
       !$acc kernels
!OCL XFILL
       do j = 1, JA
       do i = 1, IS-1
       do k = KS, KE
          DENS(k,i,j) = DENS0(k,i,j)
          MOMZ(k,i,j) = MOMZ0(k,i,j)
          MOMX(k,i,j) = MOMX0(k,i,j)
          MOMY(k,i,j) = MOMY0(k,i,j)
          RHOT(k,i,j) = RHOT0(k,i,j)
          !$acc loop seq
          do iv = 1, VA
             PROG(k,i,j,iv) = PROG0(k,i,j,iv)
          end do
       enddo
       enddo
       enddo
       !$acc end kernels
    end if
    if ( BND_E .and. (.not. TwoD) ) then
       !$omp parallel do default(none) private(j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,iv) &
       !$omp shared(JA,IE,IA,KS,KE,DENS,DENS0,MOMZ,MOMZ0,MOMX,MOMX0,MOMY,MOMY0,RHOT,RHOT0,VA,PROG,PROG0)
       !$acc kernels
!OCL XFILL
       do j = 1, JA
       do i = IE+1, IA
       do k = KS, KE
          DENS(k,i,j) = DENS0(k,i,j)
          MOMZ(k,i,j) = MOMZ0(k,i,j)
          MOMX(k,i,j) = MOMX0(k,i,j)
          MOMY(k,i,j) = MOMY0(k,i,j)
          RHOT(k,i,j) = RHOT0(k,i,j)
          !$acc loop seq
          do iv = 1, VA
             PROG(k,i,j,iv) = PROG0(k,i,j,iv)
          end do
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
       !$acc kernels
!OCL XFILL
       do j = 1, JA
       do k = KS, KE
          MOMX(k,IE,j) = MOMX0(k,IE,j)
       enddo
       enddo
       !$acc end kernels
    end if
    if ( BND_S ) then
       !$omp parallel do default(none) private(j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,iv) &
       !$omp shared(JS,IA,KS,KE,DENS,DENS0,MOMZ,MOMZ0,MOMX,MOMX0,MOMY,MOMY0,RHOT,RHOT0,VA,PROG,PROG0)
       !$acc kernels
!OCL XFILL
       do j = 1, JS-1
       do i = 1, IA
       do k = KS, KE
          DENS(k,i,j) = DENS0(k,i,j)
          MOMZ(k,i,j) = MOMZ0(k,i,j)
          MOMX(k,i,j) = MOMX0(k,i,j)
          MOMY(k,i,j) = MOMY0(k,i,j)
          RHOT(k,i,j) = RHOT0(k,i,j)
          !$acc loop seq
          do iv = 1, VA
             PROG(k,i,j,iv) = PROG0(k,i,j,iv)
          end do
       enddo
       enddo
       enddo
       !$acc end kernels
    end if
    if ( BND_N ) then
       !$omp parallel do default(none) private(j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,iv) &
       !$omp shared(JA,JE,IA,KS,KE,DENS,DENS0,MOMZ,MOMZ0,MOMX,MOMX0,MOMY,MOMY0,RHOT,RHOT0,VA,PROG,PROG0)
       !$acc kernels
!OCL XFILL
       do j = JE+1, JA
       do i = 1, IA
       do k = KS, KE
          DENS(k,i,j) = DENS0(k,i,j)
          MOMZ(k,i,j) = MOMZ0(k,i,j)
          MOMX(k,i,j) = MOMX0(k,i,j)
          MOMY(k,i,j) = MOMY0(k,i,j)
          RHOT(k,i,j) = RHOT0(k,i,j)
          !$acc loop seq
          do iv = 1, VA
             PROG(k,i,j,iv) = PROG0(k,i,j,iv)
          end do
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
       !$acc kernels
!OCL XFILL
       do i = 1, IA
       do k = KS, KE
          MOMY(k,i,JE) = MOMY0(k,i,JE)
       enddo
       enddo
       !$acc end kernels
    end if

    !$acc end data

    return
  end subroutine ATMOS_DYN_Copy_boundary

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_Copy_boundary_tracer( &
       QTRC, QTRC0, &
       BND_W, BND_E, BND_S, BND_N, TwoD )
    implicit none
    real(RP), intent(inout) :: QTRC (KA,IA,JA)
    real(RP), intent(in)    :: QTRC0(KA,IA,JA)
    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N
    logical,  intent(in)    :: TwoD

    integer :: k, i, j

    !$acc data copy(QTRC) copyin(QTRC0)

    if ( BND_W .and. (.not. TwoD) ) then
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JA,IS,KS,KE,QTRC,QTRC0)
       !$acc kernels
!OCL XFILL
       do j = 1, JA
       do i = 1, IS-1
       do k = KS, KE
          QTRC(k,i,j) = QTRC0(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
    end if
    if ( BND_E .and. (.not. TwoD) ) then
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JA,IE,IA,KS,KE,QTRC,QTRC0)
       !$acc kernels
!OCL XFILL
       do j = 1, JA
       do i = IE+1, IA
       do k = KS, KE
          QTRC(k,i,j) = QTRC0(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
    end if
    if ( BND_S ) then
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,IA,KS,KE,QTRC,QTRC0)
       !$acc kernels
!OCL XFILL
       do j = 1, JS-1
       do i = 1, IA
       do k = KS, KE
          QTRC(k,i,j) = QTRC0(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
    end if
    if ( BND_N ) then
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JA,JE,IA,KS,KE,QTRC,QTRC0)
       !$acc kernels
!OCL XFILL
       do j = JE+1, JA
       do i = 1, IA
       do k = KS, KE
          QTRC(k,i,j) = QTRC0(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
    end if

    !$acc end data

    return
  end subroutine ATMOS_DYN_Copy_boundary_tracer

  !-----------------------------------------------------------------------------

  subroutine ATMOS_DYN_divergence( &
       DDIV, &
       MOMZ, MOMX, MOMY, &
       GSQRT, J13G, J23G, J33G, MAPF, &
       TwoD, &
       RCDZ, RCDX, RCDY, RFDZ, FDZ )
    implicit none
    real(RP), intent(out) :: DDIV(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: GSQRT(KA,IA,JA,7)
    real(RP), intent(in)  :: J13G(KA,IA,JA,7)
    real(RP), intent(in)  :: J23G(KA,IA,JA,7)
    real(RP), intent(in)  :: J33G
    real(RP), intent(in)  :: MAPF(IA,JA,2,7)
    logical,  intent(in)  :: TwoD
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: RCDX(IA)
    real(RP), intent(in)  :: RCDY(JA)
    real(RP), intent(in)  :: RFDZ(KA-1)
    real(RP), intent(in)  :: FDZ(KA-1)

    integer :: k, i, j

    call PROF_rapstart("DYN_divercence", 2)

    !$acc data copyout(DDIV) copyin(MOMZ, MOMX, MOMY, GSQRT, J13G, J23G, MAPF, RCDZ, RCDX, RCDY, RFDZ, FDZ)

    ! 3D divergence

    if ( TwoD ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_
       !$acc kernels
       do j = JS, JE+1
       do k = KS-1, KE+1
          DDIV(k,IS,j) = J33G * ( MOMZ(k,IS,j) - MOMZ(k-1,IS,j) ) * RCDZ(k) &
                              + ( ( MOMY(k+1,IS,j) + MOMY(k+1,IS,j-1) ) * J23G(k+1,IS,j,I_XYW) &
                                - ( MOMY(k-1,IS,j) + MOMY(k-1,IS,j-1) ) * J23G(k-1,IS,j,I_XYW) ) / ( FDZ(k)+FDZ(k-1) ) &
                      + MAPF(IS,j,2,I_XY) &
                      * ( MOMY(k,IS,j  ) * GSQRT(k,IS,j  ,I_XVZ) / MAPF(IS,j  ,1,I_XV) &
                        - MOMY(k,IS,j-1) * GSQRT(k,IS,j-1,I_XVZ) / MAPF(IS,j-1,1,I_XV) ) * RCDY(j)
       enddo
       enddo
       !$acc end kernels
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(j) OMP_SCHEDULE_
       !$acc kernels
       do j = JS, JE+1
          DDIV(KS,IS,j) = J33G * ( MOMZ(KS,IS,j) ) * RCDZ(KS) &
                         + ( ( MOMY(KS+1,IS,j) + MOMY(KS+1,IS,j-1) ) * J23G(KS+1,IS,j,I_XYW) &
                           - ( MOMY(KS  ,IS,j) + MOMY(KS  ,IS,j-1) ) * J23G(KS  ,IS,j,I_XYW) ) * RFDZ(KS) &
                       + MAPF(IS,j,2,I_XY) &
                       * ( MOMY(KS,IS,j  ) * GSQRT(KS,IS,j  ,I_XVZ) / MAPF(IS,j  ,1,I_XV) &
                         - MOMY(KS,IS,j-1) * GSQRT(KS,IS,j-1,I_XVZ) / MAPF(IS,j-1,1,I_XV) ) * RCDY(j)
          DDIV(KE,IS,j) = J33G * ( - MOMZ(KE-1,IS,j  ) ) * RCDZ(KE) &
                       + ( ( MOMY(KE  ,IS,j) + MOMY(KE  ,IS,j-1) ) * J23G(KE  ,IS,j,I_XYW) &
                         - ( MOMY(KE-1,IS,j) + MOMY(KE-1,IS,j-1) ) * J23G(KE-1,IS,j,I_XYW) ) * RFDZ(KE) &
                       + MAPF(IS,j,2,I_XY) &
                       * ( MOMY(KE,IS,j  ) * GSQRT(KE,IS,j  ,I_XVZ) / MAPF(IS,j  ,1,I_XV) &
                         - MOMY(KE,IS,j-1) * GSQRT(KE,IS,j-1,I_XVZ) / MAPF(IS,j-1,1,I_XV) ) * RCDY(j)
       enddo
       !$acc end kernels
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
    else
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       !$acc kernels
       do j = JS, JE+1
       do i = IS, IE+1
       do k = KS-1, KE+1
          DDIV(k,i,j) = J33G * ( MOMZ(k,i,j) - MOMZ(k-1,i  ,j  ) ) * RCDZ(k) &
                      + ( ( MOMX(k+1,i,j) + MOMX(k+1,i-1,j  ) ) * J13G(k+1,i,j,I_XYW) &
                        - ( MOMX(k-1,i,j) + MOMX(k-1,i-1,j  ) ) * J13G(k-1,i,j,I_XYW) &
                        + ( MOMY(k+1,i,j) + MOMY(k+1,i  ,j-1) ) * J23G(k+1,i,j,I_XYW) &
                        - ( MOMY(k-1,i,j) + MOMY(k-1,i  ,j-1) ) * J23G(k-1,i,j,I_XYW) ) / ( FDZ(k)+FDZ(k-1) ) &
                      + MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                      * ( ( MOMX(k,i  ,j  ) * GSQRT(k,i  ,j  ,I_UYZ) / MAPF(i  ,j  ,2,I_UY) &
                          - MOMX(k,i-1,j  ) * GSQRT(k,i-1,j  ,I_UYZ) / MAPF(i-1,j  ,2,I_UY) ) * RCDX(i) &
                        + ( MOMY(k,i  ,j  ) * GSQRT(k,i  ,j  ,I_XVZ) / MAPF(i  ,j  ,1,I_XV) &
                          - MOMY(k,i,  j-1) * GSQRT(k,i  ,j-1,I_XVZ) / MAPF(i  ,j-1,1,I_XV) ) * RCDY(j) )
       enddo
       enddo
       enddo
       !$acc end kernels
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       !$acc kernels
       do j = JS, JE+1
       do i = IS, IE+1
          DDIV(KS,i,j) = J33G * ( MOMZ(KS,i,j) ) * RCDZ(KS) &
                       + ( ( MOMX(KS+1,i,j) + MOMX(KS+1,i-1,j  ) ) * J13G(KS+1,i,j,I_XYW) &
                         - ( MOMX(KS-1,i,j) + MOMX(KS  ,i-1,j  ) ) * J13G(KS  ,i,j,I_XYW) &
                         + ( MOMY(KS+1,i,j) + MOMY(KS+1,i  ,j-1) ) * J23G(KS+1,i,j,I_XYW) &
                         - ( MOMY(KS  ,i,j) + MOMY(KS  ,i  ,j-1) ) * J23G(KS  ,i,j,I_XYW) ) * RFDZ(KS) &
                       + MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                       * ( ( MOMX(KS,i  ,j  ) * GSQRT(KS,i  ,j  ,I_UYZ) / MAPF(i  ,j  ,2,I_UY) &
                           - MOMX(KS,i-1,j  ) * GSQRT(KS,i-1,j  ,I_UYZ) / MAPF(i-1,j  ,2,I_UY) ) * RCDX(i) &
                         + ( MOMY(KS,i  ,j  ) * GSQRT(KS,i  ,j  ,I_XVZ) / MAPF(i  ,j  ,1,I_XV) &
                           - MOMY(KS,i,  j-1) * GSQRT(KS,i  ,j-1,I_XVZ) / MAPF(i  ,j-1,1,I_XV) ) * RCDY(j) )
          DDIV(KE,i,j) = J33G * ( - MOMZ(KE-1,i  ,j  ) ) * RCDZ(KE) &
                       + ( ( MOMX(KE  ,i,j) + MOMX(KE  ,i-1,j  ) ) * J13G(KE  ,i,j,I_XYW) &
                         - ( MOMX(KE-1,i,j) + MOMX(KE-1,i-1,j  ) ) * J13G(KE-1,i,j,I_XYW) &
                         + ( MOMY(KE  ,i,j) + MOMY(KE  ,i  ,j-1) ) * J23G(KE  ,i,j,I_XYW) &
                         - ( MOMY(KE-1,i,j) + MOMY(KE-1,i  ,j-1) ) * J23G(KE-1,i,j,I_XYW) ) * RFDZ(KE) &
                       + MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                       * ( ( MOMX(KE,i  ,j  ) * GSQRT(KE,i  ,j  ,I_UYZ) / MAPF(i  ,j  ,2,I_UY) &
                           - MOMX(KE,i-1,j  ) * GSQRT(KE,i-1,j  ,I_UYZ) / MAPF(i-1,j  ,2,I_UY) ) * RCDX(i) &
                         + ( MOMY(KE,i  ,j  ) * GSQRT(KE,i  ,j  ,I_XVZ) / MAPF(i  ,j  ,1,I_XV) &
                           - MOMY(KE,i,  j-1) * GSQRT(KE,i  ,j-1,I_XVZ) / MAPF(i  ,j-1,1,I_XV) ) * RCDY(j) )
       enddo
       enddo
       !$acc end kernels
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
    end if

    !$acc end data

    call PROF_rapend  ("DYN_divercence", 2)

    return
  end subroutine ATMOS_DYN_divergence


   !------------------------------------------------------------------------
   ! prepare thermodynamical data
   !   specific heat
   !   pressure data ( linearization )
   !
   ! pres = P0 * ( R * rhot / P0 )**(CP/CV)
   ! d pres / d rhot ~ CP*R/CV * ( R * rhot / P0 )**(R/CV)
   !                 = CP*R/CV * ( pres / P0 )**(R/CP)
   !                 = CP*R/CV * temp / pott
   !                 = CP/CV * pres / rhot
   ! pres ~ P0 * ( R * rhot0 / P0 ) ** (CP/CV) + CV*R/CP * ( pres / P0 )**(R/CP) * rhot'
   !------------------------------------------------------------------------

  subroutine ATMOS_DYN_prep_pres_linearization( &
   DPRES, RT2P, REF_rhot,                            & ! (out)
   RHOT, QTRC, REF_pres, AQ_R, AQ_CV, AQ_CP, AQ_MASS ) ! (in)
   use scale_const, only: &
      P0     => CONST_PRE00, &
      Rdry   => CONST_Rdry,  &
      CVdry  => CONST_CVdry, &
      CPdry  => CONST_CPdry    
   implicit none

   real(RP), intent(out) :: DPRES(KA,IA,JA)
   real(RP), intent(out) :: RT2P(KA,IA,JA)
   real(RP), intent(out) :: REF_rhot(KA,IA,JA)
   real(RP), intent(in) :: RHOT(KA,IA,JA)
   real(RP), intent(in) :: QTRC(KA,IA,JA,QA)
   real(RP), intent(in) :: REF_pres(KA,IA,JA)
   real(RP), intent(in) :: AQ_R(QA)
   real(RP), intent(in) :: AQ_CV(QA)
   real(RP), intent(in) :: AQ_CP(QA)
   real(RP), intent(in) :: AQ_MASS(QA)

   integer :: i, j, k
   integer :: iq    
   real(RP) :: QDRY (KA) ! dry air
   real(RP) :: Rtot (KA) ! total R
   real(RP) :: CVtot(KA) ! total CV
   real(RP) :: CPtot(KA) ! total CP
   real(RP) :: PRES      ! pressure

#ifdef DRY
   real(RP) :: CPovCV
#endif

   !--------------------------------------

#ifdef DRY
   CPovCV = CPdry / CVdry
#endif

!OCL XFILL
   !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
   !$omp shared(JA,IA,KS,KE) &
   !$omp shared(P0,Rdry,RHOT,AQ_R,AQ_CV,AQ_CP,QTRC,AQ_MASS,REF_rhot,REF_pres,CPdry,CVdry,QA,RT2P,DPRES) &
#ifdef DRY
   !$omp shared(CPovCV) &
#endif
   !$omp private(i,j,k,iq) &
   !$omp private(PRES,Rtot,CVtot,CPtot,QDRY)
   !$acc kernels copyout(DPRES, RT2P, REF_rhot) copyin(RHOT, QTRC, REF_pres, AQ_R, AQ_CV, AQ_CP, AQ_MASS)
   do j = 1, JA
   !$acc loop private(Rtot, CVtot, CPtot, QDRY)
   do i = 1, IA
      do k = KS, KE
         Rtot (k) = 0.0_RP
         CVtot(k) = 0.0_RP
         CPtot(k) = 0.0_RP
         QDRY (k) = 1.0_RP
      end do
      !$acc loop seq
      do iq = 1, QA
         do k = KS, KE
            Rtot (k) = Rtot (k) + AQ_R (iq) * QTRC(k,i,j,iq)
            CVtot(k) = CVtot(k) + AQ_CV(iq) * QTRC(k,i,j,iq)
            CPtot(k) = CPtot(k) + AQ_CP(iq) * QTRC(k,i,j,iq)
            QDRY (k) = QDRY (k) - QTRC(k,i,j,iq) * AQ_MASS(iq)
         enddo
      end do
      do k = KS, KE
         Rtot (k) = Rtot (k) + Rdry  * QDRY(k)
         CVtot(k) = CVtot(k) + CVdry * QDRY(k)
         CPtot(k) = CPtot(k) + CPdry * QDRY(k)
      end do
      do k = KS, KE
        PRES = P0 * ( Rtot(k) * RHOT(k,i,j) / P0 )**( CPtot(k) / CVtot(k) )
        RT2P(k,i,j) = CPtot(k) / CVtot(k) * PRES / RHOT(k,i,j)
        DPRES(k,i,j) = PRES - REF_pres(k,i,j)
        REF_rhot(k,i,j) = RHOT(k,i,j)
      end do
      DPRES(KS-1,i,j) = DPRES(KS+1,i,j) + ( REF_pres(KS+1,i,j) - REF_pres(KS-1,i,j) )
      DPRES(KE+1,i,j) = DPRES(KE-1,i,j) + ( REF_pres(KE-1,i,j) - REF_pres(KE+1,i,j) )
   end do
   end do
   !$acc end kernels

   return
 end subroutine ATMOS_DYN_prep_pres_linearization

end module scale_atmos_dyn_common
