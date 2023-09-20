!-------------------------------------------------------------------------------
!> module MATRIX
!!
!! @par Description
!!          solve matrix module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#ifndef _OPENACC
#ifdef USE_CUDALIB
#undef USE_CUDALIB
#endif
#endif

#include "scalelib.h"
module scale_matrix
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
#ifdef USE_CUDALIB
  use cusparse
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MATRIX_setup
  public :: MATRIX_finalize

  public :: MATRIX_SOLVER_tridiagonal
  public :: MATRIX_SOLVER_tridiagonal_1D_TA
  public :: MATRIX_SOLVER_tridiagonal_1D_CR
  public :: MATRIX_SOLVER_tridiagonal_1D_PCR

  interface MATRIX_SOLVER_tridiagonal
#ifdef _OPENACC
     module procedure MATRIX_SOLVER_tridiagonal_1D_CR
#else
     module procedure MATRIX_SOLVER_tridiagonal_1D_TA
#endif
     module procedure MATRIX_SOLVER_tridiagonal_2D
     module procedure MATRIX_SOLVER_tridiagonal_2D_block
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
#ifdef USE_CUDALIB
  type(cusparseHandle) :: handle
  real(RP), allocatable :: pbuffer(:)
  !$acc declare device_resident(pbuffer)
  integer(8) :: bufsize
  integer :: status
#endif
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MATRIX_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

!    namelist /PARAM_MATRIX/ &

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MATRIX_setup",*) 'Setup'

    !--- read namelist
!!$    rewind(IO_FID_CONF)
!!$    read(IO_FID_CONF,nml=PARAM_MATRIX,iostat=ierr)
!!$    if( ierr < 0 ) then !--- missing
!!$       LOG_INFO("MATRIX_setup",*) 'Not found namelist. Default used.'
!!$    elseif( ierr > 0 ) then !--- fatal error
!!$       LOG_ERROR("MATRIX_setup",*) 'Not appropriate names in namelist PARAM_MATRIX. Check!'
!!$       call PRC_abort
!!$    endif
!!$    LOG_NML(PARAM_MATRIX)

#ifdef USE_CUDALIB
    status = cusparseCreate(handle)
    if ( status /= CUSPARSE_STATUS_SUCCESS ) then
       LOG_ERROR("MATRIX_setup",*) "cusparseCreate failed: ", status
       call PRC_abort
    end if
    bufsize = -1
#endif
    return
  end subroutine MATRIX_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine MATRIX_finalize

#ifdef USE_CUDALIB
    status = cusparseDestroy(handle)
    if ( allocated(pbuffer) ) deallocate(pbuffer)
    bufsize = -1
#endif

    return
  end subroutine MATRIX_finalize

  !-----------------------------------------------------------------------------
  !> solve tridiagonal matrix with Thomas's algorithm
!OCL SERIAL
  subroutine MATRIX_SOLVER_tridiagonal_1D_TA( &
       KA, KS, KE, &
#ifdef _OPENACC
       work,       &
#endif
       ud, md, ld, &
       iv,         &
       ov          )
    !$acc routine seq
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: ud(KA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA) ! input  vector

    real(RP), intent(out) :: ov(KA) ! output vector

#ifdef _OPENACC
    real(RP), intent(out) :: work(KS:KE,2)
#define c_ta(k) work(k,1)
#define d_ta(k) work(k,2)
#else
    real(RP) :: c_ta(KS:KE)
    real(RP) :: d_ta(KS:KE)
#endif
    real(RP) :: rdenom

    integer :: k
    !---------------------------------------------------------------------------

    ! foward reduction
    c_ta(KS) = ud(KS) / md(KS)
    d_ta(KS) = iv(KS) / md(KS)
    do k = KS+1, KE-1
       rdenom = 1.0_RP / ( md(k) - ld(k) * c_ta(k-1) )
       c_ta(k) =           ud(k)               * rdenom
       d_ta(k) = ( iv(k) - ld(k) * d_ta(k-1) ) * rdenom
    enddo

    ! backward substitution
    ov(KE) = ( iv(KE) - ld(KE) * d_ta(KE-1) ) / ( md(KE) - ld(KE) * c_ta(KE-1) )
    do k = KE-1, KS, -1
       ov(k) = d_ta(k) - c_ta(k) * ov(k+1)
    enddo

    return
  end subroutine MATRIX_SOLVER_tridiagonal_1D_TA

  !-----------------------------------------------------------------------------
  !> solve tridiagonal matrix with Cyclic Reduction method
!OCL SERIAL
  subroutine MATRIX_SOLVER_tridiagonal_1D_CR( &
       KA, KS, KE, &
#ifdef _OPENACC
       work,       &
#endif
       ud, md, ld, &
       iv,         &
       ov          )
    !$acc routine vector
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: ud(KA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA) ! input  vector

    real(RP), intent(out) :: ov(KA) ! output vector

#ifdef _OPENACC
    real(RP), intent(out) :: work(KE-KS+1,4)
#define a1_cr(k) work(k,1)
#define b1_cr(k) work(k,2)
#define c1_cr(k) work(k,3)
#define x1_cr(k) work(k,4)
#else
    real(RP) :: a1_cr(KE-KS+1)
    real(RP) :: b1_cr(KE-KS+1)
    real(RP) :: c1_cr(KE-KS+1)
    real(RP) :: x1_cr(KE-KS+1)
#endif
    real(RP) :: f1, f2
    integer :: st
    integer :: lmax, kmax
    integer :: k, k1, k2, l

    kmax = KE - KS + 1
    lmax = floor( log(real(kmax,RP)) / log(2.0_RP) ) - 1

    a1_cr(:) = ld(KS:KE)
    b1_cr(:) = md(KS:KE)
    c1_cr(:) = ud(KS:KE)
    x1_cr(:) = iv(KS:KE)

    st = 1
    do l = 1, lmax
       !$omp parallel do private(k1,k2,f1,f2)
       !$acc loop private(k1,k2,f1,f2)
       do k = st*2, kmax, st*2
          k1 = k - st
          k2 = k + st
          f1 = a1_cr(k) / b1_cr(k1)
          if ( k2 > kmax ) then
             k2 = k1 ! dummy
             f2 = 0.0_RP
          else
             f2 = c1_cr(k) / b1_cr(k2)
          end if
          a1_cr(k) = - a1_cr(k1) * f1
          c1_cr(k) = - c1_cr(k2) * f2
          b1_cr(k) = b1_cr(k) - c1_cr(k1) * f1 - a1_cr(k2) * f2
          x1_cr(k) = x1_cr(k) - x1_cr(k1) * f1 - x1_cr(k2) * f2
       end do
       st = st * 2
    end do
    if ( kmax / st == 2 ) then
       ov(KS+st*2-1) = ( a1_cr(st*2) * x1_cr(st) - b1_cr(st) * x1_cr(st*2) ) &
               / ( a1_cr(st*2) * c1_cr(st) - b1_cr(st) * b1_cr(st*2) )
       ov(KS+st-1) = ( x1_cr(st) - c1_cr(st) * ov(KS+st*2-1) ) / b1_cr(st)
    else if ( kmax / st == 3 ) then
       k = st * 2
       k1 = st * 2 - st
       k2 = st * 2 + st
       f2 = c1_cr(k1) / b1_cr(k)
       c1_cr(k1) = - c1_cr(k) * f2
       b1_cr(k1) = b1_cr(k1) - a1_cr(k) * f2
       x1_cr(k1) = x1_cr(k1) - x1_cr(k) * f2

       f1 = a1_cr(k2) / b1_cr(k)
       a1_cr(k2) = - a1_cr(k) * f1
       b1_cr(k2) = b1_cr(k2) - c1_cr(k) * f1
       x1_cr(k2) = x1_cr(k2) - x1_cr(k) * f1

       ov(KS+k2-1) = ( a1_cr(k2) * x1_cr(k1) - b1_cr(k1) * x1_cr(k2) ) &
              / ( a1_cr(k2) * c1_cr(k1) - b1_cr(k1) * b1_cr(k2) )
       ov(KS+k1-1) = ( x1_cr(k1) - c1_cr(k1) * ov(KS+k2-1) ) / b1_cr(k1)
       ov(KS+k-1) = ( x1_cr(k) - a1_cr(k) * ov(KS+k1-1) - c1_cr(k) * ov(KS+k2-1) ) / b1_cr(k)
    end if

    do l = 1, lmax
       st = st / 2
       !$omp parallel do
       !$acc loop independent
       do k = st, kmax, st*2
          if ( k-st < 1 ) then
             ov(KS+k-1) = ( x1_cr(k) - c1_cr(k) * ov(KS+k+st-1) ) / b1_cr(k)
          elseif ( k+st <= kmax ) then
             ov(KS+k-1) = ( x1_cr(k) - a1_cr(k) * ov(KS+k-st-1) - c1_cr(k) * ov(KS+k+st-1) ) / b1_cr(k)
          else
             ov(KS+k-1) = ( x1_cr(k) - a1_cr(k) * ov(KS+k-st-1) ) / b1_cr(k)
          end if
       end do
    end do

    return
  end subroutine MATRIX_SOLVER_tridiagonal_1D_CR

  !-----------------------------------------------------------------------------
  !> solve tridiagonal matrix with Parallel Cyclic Reduction method
!OCL SERIAL
  subroutine MATRIX_SOLVER_tridiagonal_1D_PCR( &
       KA, KS, KE, &
#ifdef _OPENACC
       work,       &
#endif
       ud, md, ld, &
       iv,         &
       ov          )
    !$acc routine vector
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: ud(KA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA) ! input  vector

    real(RP), intent(out) :: ov(KA) ! output vector

#ifdef _OPENACC
    real(RP), intent(out) :: work(KE-KS+1,2,4)
#define a1_pcr(k,n) work(k,n,1)
#define b1_pcr(k,n) work(k,n,2)
#define c1_pcr(k,n) work(k,n,3)
#define x1_pcr(k,n) work(k,n,4)
#else
    real(RP) :: a1_pcr(KE-KS+1,2)
    real(RP) :: b1_pcr(KE-KS+1,2)
    real(RP) :: c1_pcr(KE-KS+1,2)
    real(RP) :: x1_pcr(KE-KS+1,2)
#endif
    real(RP) :: f1, f2
    integer :: st
    integer :: lmax, kmax
    integer :: iw1, iw2, iws
    integer :: k, k1, k2, l

    kmax = KE - KS + 1
    lmax = ceiling( log(real(kmax,RP)) / log(2.0_RP) )

    a1_pcr(:,1) = ld(KS:KE)
    b1_pcr(:,1) = md(KS:KE)
    c1_pcr(:,1) = ud(KS:KE)
    x1_pcr(:,1) = iv(KS:KE)

    st = 1
    iw1 = 1
    iw2 = 2
    do l = 1, lmax
       !$omp parallel do private(k1,k2,f1,f2)
       !$acc loop private(k1,k2,f1,f2)
       do k = 1, kmax
          k1 = k - st
          k2 = k + st
          if ( k1 < 1 ) then
             k1 = k ! dummy
             f1 = 0.0_RP
          else
             f1 = a1_pcr(k,iw1) / b1_pcr(k1,iw1)
          end if
          if ( k2 > kmax ) then
             k2 = k ! dummy
             f2 = 0.0_RP
          else
             f2 = c1_pcr(k,iw1) / b1_pcr(k2,iw1)
          end if
          a1_pcr(k,iw2) = - a1_pcr(k1,iw1) * f1
          c1_pcr(k,iw2) = - c1_pcr(k2,iw1) * f2
          b1_pcr(k,iw2) = b1_pcr(k,iw1) - c1_pcr(k1,iw1) * f1 - a1_pcr(k2,iw1) * f2
          x1_pcr(k,iw2) = x1_pcr(k,iw1) - x1_pcr(k1,iw1) * f1 - x1_pcr(k2,iw1) * f2
       end do
       st = st * 2
       iws = iw2
       iw2 = iw1
       iw1 = iws
    end do

    !$omp parallel do
    do k = 1, kmax
       ov(KS+k-1) = x1_pcr(k,iw1) / b1_pcr(k,iw1)
    end do

    return
  end subroutine MATRIX_SOLVER_tridiagonal_1D_PCR

!OCL SERIAL
  subroutine MATRIX_SOLVER_tridiagonal_2D( &
       KA, KS, KE, &
       IA, IS, IE, &
       ud, md, ld, &
       iv,         &
       ov          )
    use scale_prc, only: &
       PRC_abort
    implicit none
    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE

    real(RP), intent(in)  :: ud(KA,IA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA,IA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA,IA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA,IA) ! input  vector

    real(RP), intent(out) :: ov(KA,IA) ! output vector

#ifdef USE_CUDALIB
    integer(8) :: bsize
#elif defined(_OPENACC)
    real(RP) :: work(KS:KE,4) ! for CR
#else
    real(RP) :: c(LSIZE,KS:KE)
    real(RP) :: d(LSIZE,KS:KE)
    real(RP) :: w(LSIZE,KS:KE)
    real(RP) :: rdenom
#endif

    integer :: k, i, ii, l
    !---------------------------------------------------------------------------

    !$acc data copyin(ud,md,ld,iv) copyout(ov)

#ifdef USE_CUDALIB
    !$acc host_data use_device(ud,md,ld,iv)
#ifdef SINGLE
    status = cusparseSgtsv2StridedBatch_bufferSizeExt( &
#else
    status = cusparseDgtsv2StridedBatch_bufferSizeExt( &
#endif
         handle, &
         KE-KS+1, &
         ld(KS,IS), md(KS,IS), ud(KS,IS), &
         iv(KS,IS), &
         IE-IS+1, KA, &
         bsize )
    !$acc end host_data
    if ( status /= CUSPARSE_STATUS_SUCCESS ) then
       LOG_ERROR("MATRIX_SOLVER_tridiagonal_2D",*) "cusparseDgtsv2StridedBatch_bufferSizeExt failed: ", status
       call PRC_abort
    end if
    if ( bsize > bufsize ) then
       if ( allocated(pbuffer) ) deallocate( pbuffer )
       allocate( pbuffer(bsize/RP) )
       bufsize = bsize
    end if

    !$acc kernels
    ov(:,IS:IE) = iv(:,IS:IE)
    !$acc end kernels

    !$acc host_data use_device(ud,md,ld,ov,pbuffer)
#ifdef SINGLE
    status = cusparseSgtsv2StridedBatch( &
#else
    status = cusparseDgtsv2StridedBatch( &
#endif
         handle, &
         KE-KS+1, &
         ld(KS,IS), md(KS,IS), ud(KS,IS), & ! (in)
         ov(KS,IS), & ! (in,out)
         IE-IS+1, KA, &
         pbuffer )
    !$acc end host_data
    if ( status /= CUSPARSE_STATUS_SUCCESS ) then
       LOG_ERROR("MATRIX_SOLVER_tridiagonal_2D",*) "cusparseDgtsv2StridedBatch failed: ", status
       call PRC_abort
    end if

#elif defined(_OPENACC)

    !$acc kernels
    !$acc loop independent private(work)
    do i = IS, IE
       call MATRIX_SOLVER_tridiagonal_1D_CR( KA, KS, KE, &
                                             work(:,:), &
                                             ud(:,i), md(:,i), ld(:,i), &
                                             iv(:,i), &
                                             ov(:,i) )
    end do
    !$acc end kernels

#else

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

#endif

    !$acc end data

    return
  end subroutine MATRIX_SOLVER_tridiagonal_2D

!OCL SERIAL
  subroutine MATRIX_SOLVER_tridiagonal_2D_block( &
       KA, KS, KE, &
       ud, md, ld, &
       iv,         &
       ov          )
    !$acc routine vector
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: ud(KA,LSIZE) ! upper  diagonal
    real(RP), intent(in)  :: md(KA,LSIZE) ! middle diagonal
    real(RP), intent(in)  :: ld(KA,LSIZE) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA,LSIZE) ! input  vector

    real(RP), intent(out) :: ov(KA,LSIZE) ! output vector

    real(RP) :: c(LSIZE,KS:KE)
    real(RP) :: d(LSIZE,KS:KE)
    real(RP) :: work(LSIZE,KS:KE)
    real(RP) :: rdenom

    integer :: k, l
    !---------------------------------------------------------------------------

    ! foward reduction
    do l = 1, LSIZE
       c(l,KS) = ud(KS,l) / md(KS,l)
       d(l,KS) = iv(KS,l) / md(KS,l)
    end do
    do k = KS+1, KE-1
!OCL NOFULLUNROLL_PRE_SIMD
       do l = 1, LSIZE
          rdenom = 1.0_RP / ( md(k,l) - ld(k,l) * c(l,k-1) )
          c(l,k) =            ud(k,l)              * rdenom
          d(l,k) = ( iv(k,l) - ld(k,l) * d(l,k-1) ) * rdenom
       end do
    end do

    ! backward substitution
    do l = 1, LSIZE
       work(l,KE) = ( iv(KE,l) - ld(KE,l) * d(l,KE-1) ) / ( md(KE,l) - ld(KE,l) * c(l,KE-1) )
    end do
    do k = KE-1, KS, -1
!OCL NOFULLUNROLL_PRE_SIMD
       do l = 1, LSIZE
          work(l,k) = d(l,k) - c(l,k) * work(l,k+1)
       end do
    end do

    do l = 1, LSIZE
    do k = KS, KE
       ov(k,l) = work(l,k)
    end do
    end do

    return
  end subroutine MATRIX_SOLVER_tridiagonal_2D_block

!OCL SERIAL
  subroutine MATRIX_SOLVER_tridiagonal_2D_trans( &
       KA, KS, KE, &
       ud, md, ld, &
       iv,         &
       ov          )
    !$acc routine vector
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
!OCL NOFULLUNROLL_PRE_SIMD
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
!OCL NOFULLUNROLL_PRE_SIMD
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
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer,  intent(in)  :: KA, KS, KE   ! array size
    integer,  intent(in)  :: IA, IS, IE   ! array size
    integer,  intent(in)  :: JA, JS, JE   ! array size
    real(RP), intent(in)  :: ud(KA,IA,JA) ! upper  diagonal
    real(RP), intent(in)  :: md(KA,IA,JA) ! middle diagonal
    real(RP), intent(in)  :: ld(KA,IA,JA) ! lower  diagonal
    real(RP), intent(in)  :: iv(KA,IA,JA) ! input  vector
    real(RP), intent(out), target :: ov(KA,IA,JA) ! output vector

    logical,  intent(in), optional :: mask(IA,JA)

#ifdef USE_CUDALIB
    real(RP), pointer :: ovl(:,:,:)
    real(RP), target :: buf(KA,IA,JA)
    integer(8) :: bsize
#elif defined(_OPENACC)
    real(RP) :: work(KS:KE,4) ! for CR
#else
    real(RP) :: udl(LSIZE,KA)
    real(RP) :: mdl(LSIZE,KA)
    real(RP) :: ldl(LSIZE,KA)
    real(RP) :: ivl(LSIZE,KA)
    real(RP) :: ovl(LSIZE,KA)
    integer  :: idx(LSIZE)
    integer  :: len
#endif

    integer :: i, j, k, l
    !---------------------------------------------------------------------------


    !$acc data copyin(ud,md,ld,iv) copyout(ov)
    !$acc data copyin(mask) if( present(mask) )

#ifdef USE_CUDALIB
    !$acc host_data use_device(ud,md,ld,iv)
#ifdef SINGLE
    status = cusparseSgtsv2StridedBatch_bufferSizeExt( &
#else
    status = cusparseDgtsv2StridedBatch_bufferSizeExt( &
#endif
         handle, &
         KE-KS+1, &
         ld(KS,1,1), md(KS,1,1), ud(KS,1,1), &
         iv(KS,1,1), &
         IA*JA, KA, &
         bsize )
    !$acc end host_data
    if ( status /= CUSPARSE_STATUS_SUCCESS ) then
       LOG_ERROR("MATRIX_SOLVER_tridiagonal_3D",*) "cusparseDgtsv2StridedBatch_bufferSizeExt failed: ", status
       call PRC_abort
    end if
    if ( bsize > bufsize ) then
       if ( bufsize > 0 ) deallocate( pbuffer )
       allocate( pbuffer(bsize/RP) )
       bufsize = bsize
    end if

    if ( present(mask) ) then
       !$acc enter data create(buf)
       ovl => buf
    else
       ovl => ov
    end if
    !$acc kernels
    ovl(:,:,:) = iv(:,:,:)
    !$acc end kernels
    !$acc host_data use_device(ud,md,ld,ovl,pbuffer)
#ifdef SINGLE
    status = cusparseSgtsv2StridedBatch( &
#else
    status = cusparseDgtsv2StridedBatch( &
#endif
         handle, &
         KE-KS+1, &
         ld(KS,1,1), md(KS,1,1), ud(KS,1,1), & ! (in)
         ovl(KS,1,1), & ! (in,out)
         IA*JA, KA, &
         pbuffer )
    !$acc end host_data
    if ( status /= CUSPARSE_STATUS_SUCCESS ) then
       LOG_ERROR("MATRIX_SOLVER_tridiagonal_3D",*) "cusparseDgtsv2StridedBatch failed: ", status
       call PRC_abort
    end if

    if ( present(mask) ) then
       !$acc kernels
       !$acc loop independent collapse(3)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( mask(i,j) ) then
             ov(k,i,j) = ovl(k,i,j)
          end if
       end do
       end do
       end do
       !$acc end kernels
       !$acc exit data delete(buf)
    end if

#elif defined(_OPENACC)

    if ( present(mask) ) then
       !$acc kernels
       !$acc loop independent
       do j = JS, JE
       !$acc loop independent private(work)
       do i = IS, IE
          if ( mask(i,j) ) then
             call MATRIX_SOLVER_tridiagonal_1D_CR( KA, KS, KE, &
                                                   work(:,:), &
                                                   ud(:,i,j), md(:,i,j), ld(:,i,j), &
                                                   iv(:,i,j), &
                                                   ov(:,i,j) )
          end if
       end do
       end do
       !$acc end kernels
    else
       !$acc kernels
       !$acc loop independent
       do j = JS, JE
       !$acc loop independent private(work)
       do i = IS, IE
          call MATRIX_SOLVER_tridiagonal_1D_CR( KA, KS, KE, &
                                                work(:,:), &
                                                ud(:,i,j), md(:,i,j), ld(:,i,j), &
                                                iv(:,i,j), &
                                                ov(:,i,j) )
       end do
       end do
       !$acc end kernels
    end if

#else

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
!OCL NOFULLUNROLL_PRE_SIMD
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
!OCL NORECURRENCE
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
!OCL NOFULLUNROLL_PRE_SIMD
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
!OCL NORECURRENCE
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
#endif

    !$acc end data
    !$acc end data

    return
  end subroutine MATRIX_SOLVER_tridiagonal_3D

end module scale_matrix
