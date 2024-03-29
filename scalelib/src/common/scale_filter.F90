!-------------------------------------------------------------------------------
!> module FILTER
!!
!! @par Description
!!          Horizontal filter
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_filter
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
  public :: FILTER_hyperdiff

  interface FILTER_hyperdiff
     module procedure FILTER_hyperdiff_2D
     module procedure FILTER_hyperdiff_3D
  end interface FILTER_hyperdiff

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
  !> Hyper diffusion filter 2D
  subroutine FILTER_hyperdiff_2D( &
       IA, IS, IE, JA, JS, JE, &
       data, order, nite, &
       limiter_sign )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_prc, only: &
       PRC_abort
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(inout) :: data(IA,JA)

    integer,  intent(in)    :: order
    integer,  intent(in)    :: nite

    real(RP), intent(in), optional :: limiter_sign(IA,JA)

    real(RP), pointer :: p1(:,:)
    real(RP), pointer :: p2(:,:)
    real(RP), target  :: work1(IA,JA)
    real(RP), target  :: work2(IA,JA)

    logical :: limiter

    integer :: i, j
    integer :: ite, n

    call PROF_rapstart('FILTER', 3)

    if ( mod(order,2) .ne. 0 ) then
       LOG_ERROR("FILTER_hyperdiff_2D", *) "order must be even"
       call PRC_abort
    end if

    limiter = present( limiter_sign )

    !$acc data copy(data) create(work1, work2)
    !$acc data copyin(limiter_sign) if(limiter)

    ! reduce grid-scale variation
    do ite = 1, nite

       call COMM_vars8( data(:,:), 1 )
       call COMM_wait ( data(:,:), 1, .false. )

       !$omp parallel do
       !$acc kernels
       do j = 1, JA
       do i = 1, IA
          work2(i,j) = data(i,j)
       end do
       end do
       !$acc end kernels

       p1 => work2
       p2 => work1
       do n = 1, order
          !$omp parallel do
          !$acc kernels
          !$acc loop collapse(2) independent
          do j = max(JS,2), min(JE,JA-1)
          do i = max(IS,2), min(IE,IA-1)
             p2(i,j) = ( - p1(i+1,j) + p1(i,j)*2.0_RP - p1(i-1,j) &
                         - p1(i,j+1) + p1(i,j)*2.0_RP - p1(i,j-1) ) / 8.0_RP
          end do
          end do
          !$acc end kernels
          if ( JS == 1 ) then
             !$acc kernels
             !$acc loop independent
             do i = max(IS,2), min(IE,IA-1)
                p2(i,JS) = ( - p1(i+1,JS) + p1(i,JS)*2.0_RP - p1(i-1,JS) &
                             - p1(i,JS+1) + p1(i,JS)                     ) / 6.0_RP
             end do
             !$acc end kernels
          else
             !$acc kernels
             !$acc loop independent
             do i = max(IS,2), min(IE,IA-1)
                p2(i,JS-1) = p2(i,JS)
             end do
             !$acc end kernels
          end if
          if ( JE == JA ) then
             !$acc kernels
             !$acc loop independent
             do i = max(IS,2), min(IE,IA-1)
                p2(i,JE) = ( - p1(i+1,JE) + p1(i,JE)*2.0_RP - p1(i-1,JE) &
                                          + p1(i,JE)        - p1(i,JE-1) ) / 6.0_RP
             end do
             !$acc end kernels
          else
             !$acc kernels
             !$acc loop independent
             do i = max(IS,2), min(IE,IA-1)
                p2(i,JE+1) = p2(i,JE)
             end do
             !$acc end kernels
          end if
          if ( IS == 1 ) then
             !$acc kernels
             !$acc loop independent
             do j = max(JS,2), min(JE,JA-1)
                p2(IS,j) = ( - p1(IS+1,j) + p1(IS,j)                     &
                             - p1(IS,j+1) + p1(IS,j)*2.0_RP - p1(IS,j-1) ) / 6.0_RP
             end do
             !$acc end kernels
             if ( JS == 1 ) then
                !$acc kernels
                p2(IS,JS) = ( - p1(IS+1,JS) + p1(IS,JS)               &
                              - p1(IS,JS+1) + p1(IS,JS)               ) / 4.0_RP
                !$acc end kernels
             end if
             if ( JE == JA ) then
                !$acc kernels
                p2(IS,JE) = ( - p1(IS+1,JE) + p1(IS,JE)               &
                                            + p1(IS,JE) - p1(IS,JE-1) ) / 4.0_RP
                !$acc end kernels
             end if
          else
             !$acc kernels
             do j = max(JS-1,1), min(JE+1,JA)
                p2(IS-1,j) = p2(IS,j)
             end do
             !$acc end kernels
          end if
          if ( IE == IA ) then
             !$acc kernels
             !$acc loop independent
             do j = max(JS,2), min(JE,JA-1)
                p2(IE,j) = (              + p1(IE,j)        - p1(IE-1,j) &
                             - p1(IE,j+1) + p1(IE,j)*2.0_RP - p1(IE,j-1) ) / 6.0_RP
             end do
             !$acc end kernels
             if ( JS == 1 ) then
                !$acc kernels
                p2(IE,JS) = (               + p1(IE,JS) - p1(IE-1,JS) &
                              - p1(IE,JS+1) + p1(IE,JS)               ) / 4.0_RP
                !$acc end kernels
             end if
             if ( JE == JA ) then
                !$acc kernels
                p2(IE,JE) = (               + p1(IE,JE) - p1(IE-1,JE) &
                                            + p1(IE,JE) - p1(IE,JE-1) ) / 4.0_RP
                !$acc end kernels
             end if
          else
             !$acc kernels
             do j = max(JS-1,1), min(JE+1,JA)
                p2(IE+1,j) = p2(IE,j)
             end do
             !$acc end kernels
          end if

          call COMM_vars8( p2(:,:), 1 )
          call COMM_wait ( p2(:,:), 1, .false. )

          if ( mod(n,2) == 0 ) then
             p1 => work2
             p2 => work1
          else
             p1 => work1
             p2 => work2
          end if
       end do

       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          data(i,j) = data(i,j) - p1(i,j)
       end do
       end do
       !$acc end kernels

       if ( limiter ) then
          !$omp parallel do
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
             data(i,j) = sign( max( data(i,j) * limiter_sign(i,j), 0.0_RP ), limiter_sign(i,j) )
          end do
          end do
          !$acc end kernels
       end if
    end do

    !$acc end data
    !$acc end data

    call PROF_rapend('FILTER', 3)

    return
  end subroutine FILTER_hyperdiff_2D

  !-----------------------------------------------------------------------------
  !> Hyper diffusion filter 3D
  subroutine FILTER_hyperdiff_3D( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       data, order, nite, &
       limiter_sign )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_prc, only: &
       PRC_abort
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(inout) :: data(KA,IA,JA)

    integer,  intent(in)    :: order
    integer,  intent(in)    :: nite

    real(RP), intent(in), optional :: limiter_sign(KA,IA,JA)

    real(RP) :: work_data(IA,JA)
    real(RP) :: work_sign(IA,JA)

    logical :: flag

    integer :: k, i, j

    if ( present(limiter_sign) ) then
       flag = .true.
    else
       flag = .false.
    end if

    !$acc data copy(data) create(work_data)
    !$acc data copyin(limiter_sign) create(work_sign) if(flag)

    do k = KS, KE

       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          work_data(i,j) = data(k,i,j)
       end do
       end do
       !$acc end kernels
       if ( flag ) then
          !$omp parallel do
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
             work_sign(i,j) = limiter_sign(k,i,j)
          end do
          end do
          !$acc end kernels

          call FILTER_hyperdiff_2D( IA, IS, IE, JA, JS, JE, &
                                    work_data(:,:), order, nite, &
                                    limiter_sign = work_sign )

       else

          call FILTER_hyperdiff_2D( IA, IS, IE, JA, JS, JE, &
                                    work_data(:,:), order, nite )
       end if


       !$omp parallel do
       !$acc kernels
       do j = 1, JA
       do i = 1, IA
          data(k,i,j) = work_data(i,j)
       end do
       end do
       !$acc end kernels

    end do

    !$acc end data
    !$acc end data

    return
  end subroutine FILTER_hyperdiff_3D

end module scale_filter
