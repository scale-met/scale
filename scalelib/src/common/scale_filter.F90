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

    ! reduce grid-scale variation
    do ite = 1, nite

       call COMM_vars8( data(:,:), 1 )
       call COMM_wait ( data(:,:), 1, .true. )

       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          work2(i,j) = data(i,j)
       end do
       end do

       p1 => work2
       p2 => work1
       do n = 1, order
          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
             p2(i,j) = ( - p1(i+1,j) + p1(i,j)*2.0_RP - p1(i-1,j) &
                         - p1(i,j+1) + p1(i,j)*2.0_RP - p1(i,j-1) ) / 8.0_RP
          end do
          end do

          call COMM_vars8( p2(:,:), 1 )
          call COMM_wait ( p2(:,:), 1, .true. )

          if ( mod(n,2) == 0 ) then
             p1 => work2
             p2 => work1
          else
             p1 => work1
             p2 => work2
          end if
       end do

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
          data(i,j) = data(i,j) - p1(i,j)
       end do
       end do

       if ( limiter ) then
          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
             data(i,j) = sign( max( data(i,j) * limiter_sign(i,j), 0.0_RP ), limiter_sign(i,j) )
          end do
          end do
       end if
    end do

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
    real(RP), target :: work_sign(IA,JA)
    real(RP), pointer :: limiter(:,:)

    integer :: k, i, j

    if ( mod(order,2) .ne. 0 ) then
       LOG_ERROR("FILTER_hyperdiff_3D", *) "order must be even"
       call PRC_abort
    end if

    if ( present(limiter_sign) ) then
       limiter => work_sign
    else
       limiter => NULL()
    end if

    do k = KS, KE

       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
          work_data(i,j) = data(k,i,j)
       end do
       end do
       if ( present(limiter_sign) ) then
          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
             limiter(i,j) = limiter_sign(k,i,j)
          end do
          end do
       end if

       call FILTER_hyperdiff_2D( IA, IS, IE, JA, JS, JE, &
                                 work_data(:,:), order, nite, &
                                 limiter_sign = limiter )

       !$omp parallel do
       do j = 1, JA
       do i = 1, IA
          data(k,i,j) = work_data(i,j)
       end do
       end do

    end do


    return
  end subroutine FILTER_hyperdiff_3D

end module scale_filter
