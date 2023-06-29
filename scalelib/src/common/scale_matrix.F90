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
     module procedure MATRIX_SOLVER_tridiagonal_2D
     module procedure MATRIX_SOLVER_tridiagonal_2D_trans
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
    real(RP) :: rdenom

    integer :: k
    !---------------------------------------------------------------------------

    ! foward reduction
    c(KS) = ud(KS) / md(KS)
    d(KS) = iv(KS) / md(KS)
    do k = KS+1, KE-1
       rdenom = 1.0_RP / ( md(k) - ld(k) * c(k-1) )
       c(k) =           ud(k)            * rdenom
       d(k) = ( iv(k) - ld(k) * d(k-1) ) * rdenom
    enddo

    ! backward substitution
    ov(KE) = ( iv(KE) - ld(KE) * d(KE-1) ) / ( md(KE) - ld(KE) * c(KE-1) )
    do k = KE-1, KS, -1
       ov(k) = d(k) - c(k) * ov(k+1)
    enddo

    return
  end subroutine MATRIX_SOLVER_tridiagonal_1D

  subroutine MATRIX_SOLVER_tridiagonal_2D( &
       KA, KS, KE, &
       IA, IS, IE, &
       ud, md, ld, &
       iv,         &
       ov          )
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE

    real(RP), intent(in)  :: ud(KA,IA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA,IA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA,IA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA,IA) ! input  vector

    real(RP), intent(out) :: ov(KA,IA) ! output vector

    real(RP) :: c(LSIZE,KS:KE)
    real(RP) :: d(LSIZE,KS:KE)
    real(RP) :: w(LSIZE,KS:KE)
    real(RP) :: rdenom

    integer :: k, i, ii, l
    !---------------------------------------------------------------------------

    do ii = IS, IE, LSIZE

       ! foward reduction
       do l = 1, LSIZE
          i = ii + l - 1
          if ( i <= IE ) then
             c(l,KS) = ud(KS,i) / md(KS,i)
             d(l,KS) = iv(KS,i) / md(KS,i)
          end if
       end do
       do k = KS+1, KE-1
          do l = 1, LSIZE
             i = ii + l - 1
             if ( i <= IE ) then
                rdenom = 1.0_RP / ( md(k,i) - ld(k,i) * c(l,k-1) )
                c(l,k) =             ud(k,i)              * rdenom
                d(l,k) = ( iv(k,i) - ld(k,i) * d(l,k-1) ) * rdenom
             end if
          end do
       end do

       ! backward substitution
       do l = 1, LSIZE
          i = ii + l - 1
          if ( i <= IE ) then
             w(l,KE) = ( iv(KE,i) - ld(KE,i) * d(l,KE-1) ) / ( md(KE,i) - ld(KE,i) * c(l,KE-1) )
          end if
       end do
       do k = KE-1, KS, -1
          do l = 1, LSIZE
             i = ii + l - 1
             if ( i <= IE ) then
                w(l,k) = d(l,k) - c(l,k) * w(l,k+1)
             end if
          end do
       enddo

       do l = 1, LSIZE
          i = ii + l - 1
          if ( i <= IE ) then
             do k = KS, KE
                ov(k,i) = w(l,k)
             end do
          end if
       end do

    end do

    return
  end subroutine MATRIX_SOLVER_tridiagonal_2D

  subroutine MATRIX_SOLVER_tridiagonal_2D_trans( &
       KA, KS, KE, &
       ud, md, ld, &
       iv,         &
       ov          )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: ud(LSIZE,KA) ! upper  diagonal
    real(RP), intent(in)  :: md(LSIZE,KA) ! middle diagonal
    real(RP), intent(in)  :: ld(LSIZE,KA) ! lower  diagonal
    real(RP), intent(in)  :: iv(LSIZE,KA) ! input  vector

    real(RP), intent(out) :: ov(LSIZE,KA) ! output vector

    real(RP) :: c(LSIZE,KS:KE)
    real(RP) :: d(LSIZE,KS:KE)
    real(RP) :: rdenom

    integer :: k, l
    !---------------------------------------------------------------------------

    ! foward reduction
    do l = 1, LSIZE
       c(l,KS) = ud(l,KS) / md(l,KS)
       d(l,KS) = iv(l,KS) / md(l,KS)
    end do
    do k = KS+1, KE-1
       do l = 1, LSIZE
          rdenom = 1.0_RP / ( md(l,k) - ld(l,k) * c(l,k-1) )
          c(l,k) =             ud(l,k)              * rdenom
          d(l,k) = ( iv(l,k) - ld(l,k) * d(l,k-1) ) * rdenom
       end do
    end do

    ! backward substitution
    do l = 1, LSIZE
       ov(l,KE) = ( iv(l,KE) - ld(l,KE) * d(l,KE-1) ) / ( md(l,KE) - ld(l,KE) * c(l,KE-1) )
    end do
    do k = KE-1, KS, -1
       do l = 1, LSIZE
          ov(l,k) = d(l,k) - c(l,k) * ov(l,k+1)
       end do
    end do

    return
  end subroutine MATRIX_SOLVER_tridiagonal_2D_trans

  !-----------------------------------------------------------------------------
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

    real(RP) :: udl(LSIZE,KA)
    real(RP) :: mdl(LSIZE,KA)
    real(RP) :: ldl(LSIZE,KA)
    real(RP) :: ivl(LSIZE,KA)
    real(RP) :: ovl(LSIZE,KA)
    integer  :: idx(LSIZE)
    integer  :: len

    integer :: i, j, k, l
    !---------------------------------------------------------------------------

    if ( present(mask) ) then
       !$omp parallel do schedule(dynamic) &
       !$omp private(len,udl,mdl,ldl,ivl,ovl,idx)
       do j = JS, JE
          len = 0
          do i = IS, IE
             if ( mask(i,j) ) then
                len = len + 1
                idx(len) = i
                if ( len == LSIZE ) then
                   do k = KS, KE
                   do l = 1, LSIZE
                      udl(l,k) = ud(k,idx(l),j)
                      mdl(l,k) = md(k,idx(l),j)
                      ldl(l,k) = ld(k,idx(l),j)
                      ivl(l,k) = iv(k,idx(l),j)
                   end do
                   end do
                   call MATRIX_SOLVER_tridiagonal_2D_trans( KA, KS, KE, &
                                                            udl(:,:), mdl(:,:), ldl(:,:),   & ! (in)
                                                            ivl(:,:),                       & ! (in)
                                                            ovl(:,:)                        ) ! (out)
                   do l = 1, LSIZE
                   do k = KS, KE
                      ov(k,idx(l),j) = ovl(l,k)
                   end do
                   end do
                   len = 0
                end if
             end if
          end do
          if ( len > 0 ) then
             do k = KS, KE
             do l = 1, len
                udl(l,k) = ud(k,idx(l),j)
                mdl(l,k) = md(k,idx(l),j)
                ldl(l,k) = ld(k,idx(l),j)
                ivl(l,k) = iv(k,idx(l),j)
             end do
#if defined DEBUG || defined QUICKDEBUG
             do l = len+1, LSIZE
                udl(l,k) = 0.0_RP
                mdl(l,k) = 1.0_RP
                ldl(l,k) = 0.0_RP
                ivl(l,k) = 0.0_RP
             end do
#endif
             end do
             call MATRIX_SOLVER_tridiagonal_2D_trans( KA, KS, KE, &
                                                      udl(:,:), mdl(:,:), ldl(:,:),   & ! (in)
                                                      ivl(:,:),                       & ! (in)
                                                      ovl(:,:)                        ) ! (out)
             do l = 1, len
             do k = KS, KE
                ov(k,idx(l),j) = ovl(l,k)
             end do
             end do
          end if
       end do
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ &
       !$omp shared(JS,JE,IA,IS,IE,KA,KS,KE,ud,md,ld,iv,ov) &
       !$omp private(j)
       do j = JS, JE
          call MATRIX_SOLVER_tridiagonal_2D( KA, KS, KE, IA, IS, IE, &
                                             ud(:,:,j), md(:,:,j), ld(:,:,j), & ! (in)
                                             iv(:,:,j),                       & ! (in)
                                             ov(:,:,j)                        ) ! (out)
       enddo
    end if

    return
  end subroutine MATRIX_SOLVER_tridiagonal_3D

end module scale_matrix
