!-------------------------------------------------------------------------------
!> module MATRIX
!!
!! @par Description
!!          solve matrix module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_matrix
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MATRIX_SOLVER_tridiagonal

  interface MATRIX_SOLVER_tridiagonal
     module procedure MATRIX_SOLVER_tridiagonal_1D
     module procedure MATRIX_SOLVER_tridiagonal_3D
  end interface MATRIX_SOLVER_tridiagonal

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
  !> solve tridiagonal matrix with Thomas's algorithm
!OCL SERIAL
  subroutine MATRIX_SOLVER_tridiagonal_1D( &
       KA, KS, KE, &
       ud, md, ld, &
       iv,         &
       ov          )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: ud(KA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA) ! input  vector

    real(RP), intent(out) :: ov(KA) ! output vector

    real(RP) :: c(KA)
    real(RP) :: d(KA)

    integer :: k
    !---------------------------------------------------------------------------

    ! foward reduction
    c(KS) = ud(KS) / md(KS)
    d(KS) = iv(KS) / md(KS)
    do k = KS+1, KE-1
       c(k) =           ud(k)            / ( md(k) - ld(k) * c(k-1) )
       d(k) = ( iv(k) - ld(k) * d(k-1) ) / ( md(k) - ld(k) * c(k-1) )
    enddo
    d(KE) = ( iv(KE) - ld(KE) * d(KE-1) ) / ( md(KE) - ld(KE) * c(KE-1) )

    ! backward substitution
    ov(KE) = d(KE)
    do k = KE-1, KS, -1
       ov(k) = d(k) - c(k) * ov(k+1)
    enddo

    return
  end subroutine MATRIX_SOLVER_tridiagonal_1D

  !-----------------------------------------------------------------------------
  !> solve tridiagonal matrix with Thomas's algorithm
  subroutine MATRIX_SOLVER_tridiagonal_3D( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       ud,         &
       md,         &
       ld,         &
       iv,         &
       ov,         &
       mask        )
    implicit none

    integer,  intent(in)  :: KA, KS, KE   ! array size
    integer,  intent(in)  :: IA, IS, IE   ! array size
    integer,  intent(in)  :: JA, JS, JE   ! array size
    real(RP), intent(in)  :: ud(KA,IA,JA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA,IA,JA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA,IA,JA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA,IA,JA) ! input  vector
    real(RP), intent(out) :: ov(KA,IA,JA) ! output vector

    logical,  intent(in), optional :: mask(IA,JA)

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( present(mask) ) then
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,KA,KS,KE,ud,md,ld,iv,ov,mask) &
       !$omp private(i,j)
       do j = JS, JE
       do i = IS, IE
          if ( mask(i,j) ) then
             call MATRIX_SOLVER_tridiagonal_1D( KA, KS, KE, &
                                                ud(:,i,j), md(:,i,j), ld(:,i,j), & ! (in)
                                                iv(:,i,j),                       & ! (in)
                                                ov(:,i,j)                        ) ! (out)
          end if
       enddo
       enddo
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,KA,KS,KE,ud,md,ld,iv,ov) &
       !$omp private(i,j)
       do j = JS, JE
       do i = IS, IE
          call MATRIX_SOLVER_tridiagonal_1D( KA, KS, KE, &
                                             ud(:,i,j), md(:,i,j), ld(:,i,j), & ! (in)
                                             iv(:,i,j),                       & ! (in)
                                             ov(:,i,j)                        ) ! (out)
       enddo
       enddo
    end if

    return
  end subroutine MATRIX_SOLVER_tridiagonal_3D

end module scale_matrix
