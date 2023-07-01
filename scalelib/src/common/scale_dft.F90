#include "scalelib.h"
! special 2D-DFT routines for low wavenumber spectral transform
module scale_dft
  use mpi
  use scale_precision
  use scale_prc, only: &
    PRC_LOCAL_COMM_WORLD, &
    PRC_myrank
  use scale_prc_cartesC, only: &
    PRC_NUM_X, &
    PRC_NUM_Y, &
    PRC_2Drank
  use scale_comm_cartesC, only: &
    COMM_Datatype
  use scale_const, only: &
    PI => CONST_PI
  implicit none

  public :: DFT_setup
  public :: DFT_finalize
  public :: DFT_g2g
  public :: DFT_g2g_divfree

  integer, private :: IMAX, JMAX ! number of total grids for spectral transform
  integer, private :: IGS, JGS ! global index
  integer, private :: LMM, MMM ! maximum truncation wavenumber

  ! table for trigonometric function
  real(RP), private, allocatable :: table_x(:,:), table_y(:,:)
  real(RP), private, allocatable :: table_l(:), table_m(:)
  ! work array
  real(RP), private, allocatable :: work(:,:,:)

contains

  subroutine DFT_setup(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM)
    implicit none
    integer,intent(in) :: KA, IA, JA
    integer,intent(in) :: KS, KE, IS, IE, JS, JE
    integer,intent(in) :: LM, MM
    real(RP) :: x, y
    integer :: i, j, k, l, m

    IMAX = (IE-IS+1)*PRC_NUM_X
    JMAX = (JE-JS+1)*PRC_NUM_Y

    IGS = (IE-IS+1)*PRC_2Drank(PRC_myrank,1)+1
    JGS = (JE-JS+1)*PRC_2Drank(PRC_myrank,2)+1

    LMM = LM
    MMM = MM

    allocate( table_x(IS:IE,0:2*LM), table_y(JS:JE,0:2*MM) )
    allocate( table_l(0:2*LM), table_m(0:2*MM) )
    allocate( work(KA,0:2*LM,JA) )

    do i = IS, IE
      x = 2*PI/IMAX*(i-IS+IGS)
      table_x(i,0) = 1
    enddo

    do l = 1, LM
      do i = IS, IE
        x = 2*PI/IMAX*(i-IS+IGS)
        table_x(i,2*l-1) =  cos(l*x)
        table_x(i,2*l)   = -sin(l*x)
      enddo
    enddo

    do j = JS, JE
      y = 2*PI/JMAX*(j-JS+JGS)
      table_y(j,0) = 1
    enddo

    do m = 1, MM
      do j = JS, JE
        y = 2*PI/JMAX*(j-JS+JGS)
        table_y(j,2*m-1) =  cos(m*y)
        table_y(j,2*m)   = -sin(m*y)
      enddo
    enddo

    table_l(0) = 1
    do l = 1, lm
      table_l(2*l-1) = cos(PI*l/IMAX)
      table_l(2*l)   = sin(PI*l/IMAX)
    enddo

    table_m(0) = 1
    do m = 1, mm
      table_m(2*m-1) = cos(PI*m/JMAX)
      table_m(2*m)   = sin(PI*m/JMAX)
    enddo

    !$acc enter data copyin(table_x, table_y)
    !$acc enter data copyin(table_l, table_m)
    !$acc enter data create(work)

  end subroutine DFT_setup

  subroutine DFT_finalize

    !$acc exit data delete(table_x, table_y)
    !$acc exit data delete(table_l, table_m)
    !$acc exit data delete(work)
    deallocate( table_x, table_y )
    deallocate( table_l, table_m )
    deallocate( work )

    return
  end subroutine DFT_finalize

  subroutine DFT_g2s(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,f,s)
    ! Grid to spectral transformation
    integer,intent(in) :: KA, IA, JA
    integer,intent(in) :: KS, KE, IS, IE, JS, JE
    integer,intent(in) :: LM, MM
    real(RP), intent(in) :: f(KA, IA, JA) ! x: full level, y: full level
    real(RP), intent(out) :: s(KA,0:2*LM,0:2*MM)
    real(RP) ::  work_s(KA,0:2*LM,0:2*MM)
    real(RP) :: c, tb
    integer :: i, j, k, l, m
    integer :: ierr

    !$acc data copyin(f) copyout(s) create(work_s)

    c = 1.0_RP/IMAX
    !$acc kernels
    do j = JS, JE
      do l = 0, 2*LM
        do k = KS, KE
          work(k,l,j) = 0
        enddo
        !$acc loop seq
        do i = IS, IE
          tb = table_x(i,l)*c
          do k = KS, KE
            work(k,l,j) = work(k,l,j) + f(k,i,j)*tb
          enddo
        enddo
      enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do m = 0, 2*MM
      do l = 0, 2*LM
        do k = KS, KE
          work_s(k,l,m) = 0
        enddo
      enddo
    enddo
    !$acc end kernels

    c = 1.0_RP/JMAX
    !$acc kernels
    do m = 0, 2*MM
      do j = JS, JE
        tb = table_y(j,m)*c
        do l = 0, 2*LM
          do k = KS, KE
             !$acc atomic
            work_s(k,l,m) = work_s(k,l,m) + work(k,l,j)*tb
          enddo
        enddo
      enddo
    enddo
    !$acc end kernels

    call MPI_Allreduce(work_s, s, KA*(2*LM+1)*(2*MM+1), COMM_Datatype, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr)

    !$acc end data

  end subroutine DFT_g2s

  subroutine DFT_s2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,s,f)
    ! Grid to spectral transformation
    integer,intent(in) :: KA, IA, JA
    integer,intent(in) :: KS, KE, IS, IE, JS, JE
    integer,intent(in) :: LM, MM
    real(RP), intent(in) :: s(KA,0:2*LM,0:2*MM)
    real(RP), intent(out) :: f(KA, IA, JA) ! x: full level, y: full level
    real(RP) :: c, tb
    integer :: i, j, k, l, m

    !$acc data copyin(s) copyout(f)

    !$acc kernels
    do j = JS, JE
      do i = IS, IE
        do k = KS, KE
          f(k,i,j) = 0
        enddo
      enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JS, JE
      do l = 0, 2*LM
        do k = KS, KE
          work(k,l,j) = 0
        enddo
      enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do m = 0, 2*MM
      if( m == 0 ) then
        c = 1
      else
        c = 2
      endif
      do j = JS, JE
        tb = table_y(j,m)*c
        do l = 0, 2*LM
          do k = KS, KE
             !$acc atomic
             work(k,l,j) = work(k,l,j) + s(k,l,m)*tb
          enddo
        enddo
      enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JS, JE
      do i = IS, IE
        do k = KS, KE
          f(k,i,j) = 0
        enddo
      enddo
      do l = 0, 2*LM
        if( l == 0 ) then
          c = 1
        else
          c = 2
        endif
        do i = IS, IE
          tb = table_x(i,l)*c
          do k = KS, KE
             !$acc atomic
             f(k,i,j) = f(k,i,j) + work(k,l,j)*tb
          enddo
        enddo
      enddo
    enddo
    !$acc end kernels

    !$acc end data

  end subroutine DFT_s2g

  ! Grid to Grid transformation
  subroutine DFT_g2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,f)
    integer,intent(in) :: KA, IA, JA
    integer,intent(in) :: KS, KE, IS, IE, JS, JE
    integer,intent(in) :: LM, MM
    real(RP), intent(inout) :: f(KA, IA, JA) ! x,y: full or half level
    real(RP) :: s(KA,0:2*LM,0:2*MM)

    call DFT_g2s(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,f,s)
    call DFT_s2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,s,f)

  end subroutine DFT_g2g

  ! Grid to Grid transformation of vector field with horizontal divergence filter
  subroutine DFT_g2g_divfree(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,u,v)
    integer,intent(in) :: KA, IA, JA
    integer,intent(in) :: KS, KE, IS, IE, JS, JE
    integer,intent(in) :: LM, MM
    real(RP), intent(inout) :: u(KA, IA, JA) ! x: half level, y: full level
    real(RP), intent(inout) :: v(KA, IA, JA) ! x: full level, y: half level
    integer :: k, l, m
    real(RP) :: s1(KA,0:2*LM,0:2*MM)
    real(RP) :: s2(KA,0:2*LM,0:2*MM)
    real(RP) :: s3(KA,0:2*LM,0:2*MM)
    real(RP) :: a, b, fac

    !$acc data copy(u,v) create(s1,s2,s3)

    call DFT_g2s(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,u,s1)

    ! phase shift
    !$acc kernels
    do m = 0, 2*mm
      do l = 1, lm
        do k = KS, KE
          a = s1(k,2*l-1,m)
          b = s1(k,2*l,m)
          s1(k,2*l-1,m) = a*table_l(2*l-1) + b*table_l(2*l)
          s1(k,2*l,m)   = b*table_l(2*l-1) - a*table_l(2*l)
        enddo
      enddo
    enddo
    !$acc end kernels

    call DFT_g2s(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,v,s2)

    ! phase shift
    !$acc kernels
    do m = 1, mm
      do l = 0, 2*lm
        do k = KS, KE
          a = s2(k,l,2*m-1)
          b = s2(k,l,2*m)
          s2(k,l,2*m-1) = a*table_m(2*m-1) + b*table_m(2*m)
          s2(k,l,2*m)   = b*table_m(2*m-1) - a*table_m(2*m)
        enddo
      enddo
    enddo
    !$acc end kernels

    ! rotation
    !$acc kernels
    do m = 0, 2*mm
      do l = 0, 2*lm
        do k = KS, KE
          s3(k,l,m) = 0
        enddo
      enddo
    enddo
    !$acc end kernels

    ! ∂v/∂x
    !$acc kernels
    do m = 0, 2*mm
      do l = 1, lm
        do k = KS, KE
          s3(k,2*l-1,m) = -l*s2(k,2*l,m)
          s3(k,2*l,m) = l*s2(k,2*l-1,m)
        enddo
      enddo
    enddo
    !$acc end kernels

    ! + (- ∂u/∂y)
    !$acc kernels
    do m = 1, mm
      do l = 0, 2*lm
        do k = KS, KE
          s3(k,l,2*m-1) = s3(k,l,2*m-1) + m*s1(k,l,2*m)
          s3(k,l,2*m)   = s3(k,l,2*m)   - m*s1(k,l,2*m-1)
        enddo
      enddo
    enddo
    !$acc end kernels

    ! minus inverse laplacian ( stream function on model plane )
    !$acc kernels
    do k = KS, KE
      s3(k,0,0) = 0
    enddo
    !$acc end kernels

    !$acc kernels
    do l = 1, lm
      fac = 1.0_RP/ (l*l)
      do k = KS, KE
        s3(k,2*l-1,0) = s3(k,2*l-1,0)*fac
        s3(k,2*l,0)   = s3(k,2*l,0)*fac
      enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do m = 1, mm
      fac = 1.0_RP/ (m*m)
      do k = KS, KE
        s3(k,0,2*m-1) = s3(k,0,2*m-1)*fac
        s3(k,0,2*m)   = s3(k,0,2*m)*fac
      enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do m = 1, mm
      do l = 1, lm
        fac = 1.0_RP/ (l*l+m*m)
        do k = KS, KE
          s3(k,2*l-1,2*m-1) = s3(k,2*l-1,2*m-1)*fac
          s3(k,2*l-1,2*m)   = s3(k,2*l-1,2*m)*fac
          s3(k,2*l,2*m-1)   = s3(k,2*l,2*m-1)*fac
          s3(k,2*l,2*m)     = s3(k,2*l,2*m)*fac
        enddo
      enddo
    enddo
    !$acc end kernels

    ! divergence free 2D vector

    ! ∂ψ/∂y
    !$acc kernels
    do l = 1, 2*lm
      do k = KS, KE
        s1(k,l,0) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    do m = 1, mm
      do l = 0, 2*lm
        do k = KS, KE
          s1(k,l,2*m-1) = -m*s3(k,l,2*m)
          s1(k,l,2*m)   =  m*s3(k,l,2*m-1)
        enddo
      enddo
    enddo
    !$acc end kernels

    ! -∂ψ/∂x
    !$acc kernels
    do m = 1, 2*mm
      do k = KS, KE
        s2(k,0,m) = 0
      enddo
    enddo
    !$acc end kernels
    !$acc kernels
    do m = 0, 2*mm
      do l = 1, lm
        do k = KS, KE
          s2(k,2*l-1,m) =  l*s3(k,2*l,m)
          s2(k,2*l,m)   = -l*s3(k,2*l-1,m)
        enddo
      enddo
    enddo
    !$acc end kernels

    ! phase shift
    !$acc kernels
    do m = 0, 2*mm
      do l = 1, lm
        do k = KS, KE
          a = s1(k,2*l-1,m)
          b = s1(k,2*l,m)
          s1(k,2*l-1,m) = a*table_l(2*l-1) - b*table_l(2*l)
          s1(k,2*l,m)   = b*table_l(2*l-1) + a*table_l(2*l)
        enddo
      enddo
    enddo
    !$acc end kernels

    call DFT_s2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,s1,u)

    ! phase shift
    !$acc kernels
    do m = 1, mm
      do l = 0, 2*lm
        do k = KS, KE
          a = s2(k,l,2*m-1)
          b = s2(k,l,2*m)
          s2(k,l,2*m-1) = a*table_m(2*m-1) - b*table_m(2*m)
          s2(k,l,2*m) = b*table_m(2*m-1) + a*table_m(2*m)
        enddo
      enddo
    enddo
    !$acc end kernels

    call DFT_s2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,LM,MM,s2,v)

    !$acc end data

  end subroutine DFT_g2g_divfree

end module scale_dft
