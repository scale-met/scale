!-------------------------------------------------------------------------------
!> module atmosphere / physics / pbl / common
!!
!! @par Description
!!          Common routines for boundary layer turbulence
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_bl_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

#if defined DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
#endif
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_BL_tendency_tracer

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
  !> ATMOS_PHY_BL_tendency_tracer
  !! calculate tendency of tracers by the eddy viscosity
  !<
  subroutine ATMOS_PHY_BL_tendency_tracer( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, QTRC, SFLX_Q, &
       Kh, MASS,           &
       CZ, FZ, F2H, DDT,   &
       TRACER_NAME,        &
       RHOQ_t              )
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_file_history, only: &
       FILE_HISTORY_in
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP),         intent(in) :: DENS  (KA,IA,JA) !> density
    real(RP),         intent(in) :: QTRC  (KA,IA,JA) !> tracers
    real(RP),         intent(in) :: SFLX_Q(   IA,JA) !> surface flux
    real(RP),         intent(in) :: Kh    (KA,IA,JA) !> eddy diffusion coefficient @ half-level
    real(RP),         intent(in) :: MASS             !> mass
    real(RP),         intent(in) :: CZ (  KA,IA,JA)  !> z at the full level
    real(RP),         intent(in) :: FZ (0:KA,IA,JA)  !> z at the half level
    real(RP),         intent(in) :: F2H(KA,2,IA,JA)  !> coefficients to convert value at the full to half level
    real(DP),         intent(in) :: DDT              !> time step
    character(len=*), intent(in) :: TRACER_NAME      !> name of tracer (for history output)

    real(RP), intent(out) :: RHOQ_t(KA,IA,JA) !> tendency of tracers

    real(RP) :: QTRC_n(KA) !> value at the next time step
    real(RP) :: RHO   (KA)
    real(RP) :: RHOKh (KA)
    real(RP) :: a(KA)
    real(RP) :: b(KA)
    real(RP) :: c(KA)
    real(RP) :: d(KA)
    real(RP) :: rho_h
    real(RP) :: ap
    real(RP) :: sf_t

    real(RP) :: flx(KA,IA,JA)

    real(RP) :: CDZ(KA)
    real(RP) :: FDZ(KA)

    real(RP) :: dt

    integer :: KE_PBL
    integer :: k, i, j

    dt = real( DDT, kind=RP )

!OCL INDEPENDENT
    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE) &
    !$omp shared(RHOQ_t,DENS,QTRC,SFLX_Q,Kh,MASS,CZ,FZ,F2H,DT,flx) &
    !$omp private(QTRC_n,RHO,RHOKh,rho_h,a,b,c,d,ap,sf_t,CDZ,FDZ) &
    !$omp private(KE_PBL,k,i,j)
    do j = JS, JE
    do i = IS, IE

       KE_PBL = KE-1
       do k = KE-2, KS+1, -1
          if ( Kh(k,i,j) > 0.0_RP ) then
             KE_PBL = k + 1
          else
             exit
          end if
       end do

       do k = KS, KE_PBL
          CDZ(k) = FZ(k  ,i,j) - FZ(k-1,i,j)
          FDZ(k) = CZ(k+1,i,j) - CZ(k  ,i,j)
       end do

       sf_t = SFLX_Q(i,j) / CDZ(KS)

       RHO(KS) = DENS(KS,i,j) + dt * sf_t * MASS
       do k = KS+1, KE_PBL
          RHO(k) = DENS(k,i,j)
       end do

       ! dens * coefficient at the half level
       do k = KS, KE_PBL-1
          rho_h = F2H(k,1,i,j) * DENS(k+1,i,j) + F2H(k,2,i,j) * DENS(k,i,j)
          RHOKh(k) = rho_h * Kh(k,i,j)
       end do

       d(KS) = ( QTRC(KS,i,j) * DENS(KS,i,j) + dt * sf_t ) / RHO(KS)
       do k = KS+1, KE_PBL
          d(k) = QTRC(k,i,j)
       end do

       c(KS) = 0.0_RP
       do k = KS, KE_PBL-1
          ap = - dt * RHOKh(k) / FDZ(k)
          a(k) = ap / ( RHO(k) * CDZ(k) )
          b(k) = - a(k) - c(k) + 1.0_RP
          c(k+1) = ap / ( RHO(k+1) * CDZ(k+1) )
       end do
       a(KE_PBL) = 0.0_RP
       b(KE_PBL) = - c(KE_PBL) + 1.0_RP

       call MATRIX_SOLVER_tridiagonal( &
               KA, KS, KE_PBL, &
               a(:), b(:), c(:), d(:), & ! (in)
               QTRC_n(:)               ) ! (out)

       RHOQ_t(KS,i,j) = ( QTRC_n(KS) * RHO(KS) - QTRC(KS,i,j) * DENS(KS,i,j) ) / dt - sf_t
       do k = KS+1, KE_PBL
          RHOQ_t(k,i,j) = ( QTRC_n(k) - QTRC(k,i,j) ) * RHO(k) / dt
       end do
       do k = KE_PBL+1, KE
          RHOQ_t(k,i,j) = 0.0_RP
       end do

       do k = KS, KE_PBL-1
          flx(k,i,j) = - RHOKh(k) * ( QTRC_n(k+1) - QTRC_n(k) ) / FDZ(k)
       end do

    end do
    end do

    call FILE_HISTORY_in(flx(:,:,:), 'ZFLX_'//trim(TRACER_NAME)//'_BL', 'Z FLUX of DENS * '//trim(TRACER_NAME)//' (PBL)', 'kg/m2/s', fill_halo=.true.)

    return
  end subroutine ATMOS_PHY_BL_tendency_tracer

end module scale_atmos_phy_bl_common
