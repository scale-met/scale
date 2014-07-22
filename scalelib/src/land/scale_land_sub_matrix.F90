!-------------------------------------------------------------------------------
!> module LAND / Physics
!!
!! @par Description
!!          solve matrix module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module scale_land_sub_matrix
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SOLVER_tridiagonal_matrix

  interface SOLVER_tridiagonal_matrix
     module procedure SOLVER_tridiagonal_matrix_1D
     module procedure SOLVER_tridiagonal_matrix_3D
  end interface SOLVER_tridiagonal_matrix

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
  subroutine SOLVER_tridiagonal_matrix_1D( &
       KA, &
       ud, &
       md, &
       ld, &
       iv, &
       ov  )
    implicit none

    integer,  intent(in)  :: KA     ! array size
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
    c(1) = ud(1) / md(1)
    d(1) = iv(1) / md(1)
    do k = 2, KA
       c(k) =           ud(k)            / ( md(k) - ld(k) * c(k-1) )
       d(k) = ( iv(k) - ld(k) * d(k-1) ) / ( md(k) - ld(k) * c(k-1) )
    enddo

    ! backward substitution
    ov(KA) = d(KA)
    do k = KA-1, 1, -1
       ov(k) = d(k) - c(k) * ov(k+1)
    enddo

    return
  end subroutine SOLVER_tridiagonal_matrix_1D

  !-----------------------------------------------------------------------------
  !> solve tridiagonal matrix with Thomas's algorithm
  subroutine SOLVER_tridiagonal_matrix_3D( &
       KA,         &
       IA, IS, IE, &
       JA, JS, JE, &
       ud,         &
       md,         &
       ld,         &
       iv,         &
       ov          )
    implicit none

    integer,  intent(in)  :: KA           ! array size
    integer,  intent(in)  :: IA, IS, IE   ! array size
    integer,  intent(in)  :: JA, JS, JE   ! array size
    real(RP), intent(in)  :: ud(KA,IA,JA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA,IA,JA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA,IA,JA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA,IA,JA) ! input  vector
    real(RP), intent(out) :: ov(KA,IA,JA) ! output vector

    real(RP) :: c(KA,IA,JA)
    real(RP) :: d(KA,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! foward reduction
    do j = JS, JE
    do i = IS, IE
       c(1,i,j) = ud(1,i,j) / md(1,i,j)
       d(1,i,j) = iv(1,i,j) / md(1,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 2,  KA
       c(k,i,j) =               ud(k,i,j)                &
                / ( md(k,i,j) - ld(k,i,j) * c(k-1,i,j) )
       d(k,i,j) = ( iv(k,i,j) - ld(k,i,j) * d(k-1,i,j) ) &
                / ( md(k,i,j) - ld(k,i,j) * c(k-1,i,j) )
    enddo
    enddo
    enddo

    ! backward substitution
    do j = JS, JE
    do i = IS, IE
       ov(KA,i,j) = d(KA,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KA-1, 1, -1
       ov(k,i,j) = d(k,i,j) - c(k,i,j) * ov(k+1,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine SOLVER_tridiagonal_matrix_3D

end module scale_land_sub_matrix
