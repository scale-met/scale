#include "inc_openmp.h"
module scale_atmos_dyn_bc
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
#ifdef CHECK_MASS
  use mpi
#endif
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
  public :: ATMOS_DYN_bc_apply

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

  subroutine ATMOS_DYN_bc_apply( &
     & DENS, MOMX, MOMY, MOMZ, RHOT, PROG,      &
     & DENS0, MOMX0, MOMY0, MOMZ0, RHOT0, PROG0 &
     & )

  use scale_les_process, only: &
       PRC_HAS_E, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_S
  
  real(RP), intent(inout), dimension(KA,IA,JA) :: DENS, MOMX, MOMY, MOMZ, RHOT
  real(RP), intent(inout) :: PROG(KA,IA,JA,VA)
  real(RP), intent(in), dimension(KA,IA,JA) :: DENS0, MOMX0, MOMY0, MOMZ0, RHOT0
  real(RP), intent(in) :: PROG0(KA,IA,JA,VA)

  integer :: i, j, k, iv
  real(RP), parameter :: Ulid = 10.0_RP!3.47212902986050
  
  if( .not. PRC_HAS_W ) then
    !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = 1, JA
    do i = 1, IS-1
    do k = KS, KE
      DENS(k,i,j) = DENS0(k,i,j)
      MOMZ(k,i,j) = MOMZ0(k,i,j)
      MOMX(k,i,j) = MOMX0(k,i,j)
      MOMY(k,i,j) = MOMY0(k,i,j)
      RHOT(k,i,j) = RHOT0(k,i,j)
      do iv = 1, VA
        PROG(k,i,j,iv) = PROG0(k,i,j,iv)
      enddo
    enddo
    enddo
    enddo
    
 !
#ifdef DYNBC_W_RIGID
!    write(*,*) "DYN_BC_W_RIGID"
    do j = 1, JA
    do i = 1, IHALO
    do k = KS, KE
      DENS(k,IS-i,j) =   DENS(k,IS+i-1,j)
      MOMZ(k,IS-i,j) = - MOMZ(k,IS+i-1,j)
      MOMY(k,IS-i,j) = - MOMY(k,IS+i-1,j)
      RHOT(k,IS-i,j) =   RHOT(k,IS+i-1,j)
    enddo
    enddo
    enddo

    MOMX(KS:KE,IS-1,1:JA) = 0.0_RP
    do j = 1, JA
    do i = 1, IHALO-1
    do k = KS, KE
      MOMX(k,IS-i-1,j) = - MOMX(k,IS+i-1,j)
    enddo
    enddo
    enddo
#endif
  endif

  if( .not. PRC_HAS_E ) then
    !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
 !OCL XFILL
    do j = 1, JA
    do i = IE+1, IA
    do k = KS, KE
      DENS(k,i,j) = DENS0(k,i,j)
      MOMZ(k,i,j) = MOMZ0(k,i,j)
      MOMX(k,i,j) = MOMX0(k,i,j)
      MOMY(k,i,j) = MOMY0(k,i,j)
      RHOT(k,i,j) = RHOT0(k,i,j)
      do iv = 1, VA
        PROG(k,i,j,iv) = PROG0(k,i,j,iv)
      end do
    enddo
    enddo
    enddo
    !$omp parallel do private(j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do j = 1, JA
    do k = KS, KE
      MOMX(k,IE,j) = MOMX0(k,IE,j)
    enddo
    enddo
 
#ifdef DYNBC_E_RIGID
!    write(*,*) "DYN_BC_E_RIGID"
    do j = 1, JA
    do i = 1, IHALO
    do k = KS, KE
      DENS(k,IE+i,j) =   DENS(k,IE-i+1,j)
      MOMZ(k,IE+i,j) = - MOMZ(k,IE-i+1,j)
      MOMY(k,IE+i,j) = - MOMY(k,IE-i+1,j)
      RHOT(k,IE+i,j) =   RHOT(k,IE-i+1,j)
    enddo
    enddo
    enddo

    MOMX(KS:KE,IE,1:JA) = 0.0_RP
    do j = 1, JA
    do i = 1, IHALO-1
    do k = KS, KE
      MOMX(k,IE+i,j) = - MOMX(k,IE-i,j)
    enddo
    enddo
    enddo
#endif
  endif

  if( .not. PRC_HAS_N ) then
!OCL XFILL
    do j = JE+1, JA
    do i = 1, IA
    do k = KS, KE
      DENS(k,i,j) = DENS0(k,i,j)
      MOMZ(k,i,j) = MOMZ0(k,i,j)
      MOMX(k,i,j) = MOMX0(k,i,j)
      MOMY(k,i,j) = MOMY0(k,i,j)
      RHOT(k,i,j) = RHOT0(k,i,j)
      do iv = 1, VA
        PROG(k,i,j,iv) = PROG0(k,i,j,iv)
      end do
    enddo
    enddo
    enddo
    !$omp parallel do private(i,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
    do i = 1, IA
    do k = KS, KE
       MOMY(k,i,JE) = MOMY0(k,i,JE)
    enddo
    enddo

  endif

  if( .not. PRC_HAS_S ) then
!OCL XFILL
    do j = 1, JS-1
    do i = 1, IA
    do k = KS, KE
      DENS(k,i,j) = DENS0(k,i,j)
      MOMZ(k,i,j) = MOMZ0(k,i,j)
      MOMX(k,i,j) = MOMX0(k,i,j)
      MOMY(k,i,j) = MOMY0(k,i,j)
      RHOT(k,i,j) = RHOT0(k,i,j)
      do iv = 1, VA
        PROG(k,i,j,iv) = PROG0(k,i,j,iv)
      end do
    enddo
    enddo
    enddo
  endif

#ifdef DYNBC_B_RIGID
!    write(*,*) "DYN_BC_B_RIGID"
    do j = 1, JA
    do i = 1, IA
    do k = 1, KHALO
      DENS(KS-k,i,j) =   DENS(KS+k-1,i,j)
      MOMX(KS-k,i,j) = - MOMX(KS+k-1,i,j)
      MOMY(KS-k,i,j) = - MOMY(KS+k-1,i,j)
      RHOT(KS-k,i,j) =   RHOT(KS+k-1,i,j)
    enddo
    enddo
    enddo

    MOMZ(KS-1,1:IA,1:JA) = 0.0_RP
    do j = 1, JA
    do i = 1, IA
    do k = 1, KHALO-1
      MOMZ(KS-k-1,i,j) = - MOMZ(KS+k-1,i,j)
    enddo
    enddo
    enddo
    
#endif

#ifdef DYNBC_T_RIGID
!    write(*,*) "DYN_BC_T_RIGID,CAVITY"
    do j = 1, JA
    do i = 1, IA
    do k = 1, KHALO
      DENS(KE+k,i,j) =   DENS(KE-k+1,i,j)
      MOMX(KE+k,i,j) = - MOMX(KE-k+1,i,j)
      MOMY(KE+k,i,j) = - MOMY(KE-k+1,i,j)
      RHOT(KE+k,i,j) =   RHOT(KE-k+1,i,j)
    enddo
    enddo
    enddo
    
    MOMZ(KE,1:IA,1:JA) = 0.0_RP
    MOMX(KE+1,1:IA,1:JA) = DENS(KE+1,1:IA,1:JA) * ( &
         & 2.0_RP * Ulid - MOMX(KE,1:IA,1:JA)/DENS(KE,1:IA,1:JA) )
    do j = 1, JA
    do i = 1, IA
    do k = 1, KHALO-1
      MOMZ(KE+k,i,j) = - MOMZ(KE-k,i,j)
    enddo
    enddo
    enddo
#endif
    
end subroutine ATMOS_DYN_bc_apply


end module scale_atmos_dyn_bc
