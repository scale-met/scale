!-------------------------------------------------------------------------------
!> module SORT
!!
!! @par Description
!!          Sort data
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_sort
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SORT_exec
  public :: SORT_uniq_int_sorted
  public :: SORT_quicksort
  public :: SORT_heapsort
  public :: SORT_quickselect
  public :: SORT_quickselect_arg
  public :: SORT_quickselect_desc
  public :: SORT_quickselect_desc_arg

  interface SORT_exec
     module procedure SORT_exec_without_idx
     module procedure SORT_exec_with_idxs
     module procedure SORT_exec_with_idx
  end interface SORT_exec

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
  private :: swap
  private :: swap_i
  private :: partition
  private :: partition_arg
  private :: partition_desc
  private :: partition_desc_arg
  private :: median_of_three
  private :: median_of_three_arg
  private :: sample_second_min
  private :: sample_second_min_arg
  private :: sample_second_max
  private :: sample_second_max_arg

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> bubble sort
!OCL SERIAL
  subroutine SORT_exec_with_idxs( &
      npoints,      &
      val,          &
      idx_i, idx_j, &
      reverse       )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: npoints                ! number of interpolation points
    real(RP), intent(inout) :: val  (npoints)         ! value to sort
    integer,  intent(inout) :: idx_i(npoints)         ! i-index
    integer,  intent(inout) :: idx_j(npoints)         ! j-index

    logical,  intent(in), optional :: reverse

    real(RP) :: sig
    integer  :: itmp
    integer  :: jtmp
    real(RP) :: vtmp

    integer  :: n1, n2
    !---------------------------------------------------------------------------

    sig = 1.0_RP
    if ( present(reverse) ) then
       if ( reverse ) sig = -1.0_RP
    end if

    do n1 = 1, npoints-1
    do n2 = n1+1, npoints
       if ( val(n1) * sig > val(n2) * sig ) then
          itmp      = idx_i(n1)
          jtmp      = idx_j(n1)
          vtmp      = val  (n1)

          idx_i(n1) = idx_i(n2)
          idx_j(n1) = idx_j(n2)
          val  (n1) = val  (n2)

          idx_i(n2) = itmp
          idx_j(n2) = jtmp
          val  (n2) = vtmp
       endif
    enddo
    enddo

    return
  end subroutine SORT_exec_with_idxs

!OCL SERIAL
  subroutine SORT_exec_with_idx( &
      npoints,    &
      val, index, &
      reverse     )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: npoints                ! number of interpolation points
    real(RP), intent(inout) :: val  (npoints)         ! value to sort
    integer,  intent(inout) :: index(npoints)         ! index

    logical,  intent(in), optional :: reverse

    real(RP) :: sig
    integer  :: itmp
    real(RP) :: vtmp

    integer  :: n1, n2
    !---------------------------------------------------------------------------

    sig = 1.0_RP
    if ( present(reverse) ) then
       if ( reverse ) sig = -1.0_RP
    end if

    do n1 = 1, npoints-1
    do n2 = n1+1, npoints
       if ( val(n1) * sig > val(n2) * sig ) then
          itmp    = index(n1)
          vtmp    = val  (n1)

          index(n1) = index(n2)
          val  (n1) = val  (n2)

          index(n2) = itmp
          val  (n2) = vtmp
       endif
    enddo
    enddo

    return
  end subroutine SORT_exec_with_idx

!OCL SERIAL
  subroutine SORT_exec_without_idx( &
      npoints, &
      val,     &
      reverse  )
    !$acc routine seq
    implicit none
    integer,  intent(in)    :: npoints                ! number of interpolation points
    real(RP), intent(inout) :: val  (npoints)         ! value to sort

    logical,  intent(in), optional :: reverse

    real(RP) :: sig
    real(RP) :: vtmp

    integer  :: n1, n2
    !---------------------------------------------------------------------------

    sig = 1.0_RP
    if ( present(reverse) ) then
       if ( reverse ) sig = -1.0_RP
    end if

    do n1 = 1, npoints-1
    do n2 = n1+1, npoints
       if ( val(n1) * sig > val(n2) * sig ) then
          vtmp      = val  (n1)

          val  (n1) = val  (n2)

          val  (n2) = vtmp
       endif
    enddo
    enddo

    return
  end subroutine SORT_exec_without_idx

  subroutine SORT_uniq_int_sorted(n, ary, c)
    integer(DP), intent(in)  :: n
    integer(DP), intent(inout) :: ary(n)
    integer(DP), intent(out) :: c
    integer(DP) :: i, j

    c = 1
    j = ary(1)
    do i = 2, n
       if(j .ne. ary(i)) then
          j = ary(i)
          c = c + 1
          ary(c) = ary(i)
       end if
    end do

    return
  end subroutine SORT_uniq_int_sorted

  recursive subroutine SORT_quicksort(n, array)
    implicit none
    integer(DP), intent(in) :: n
    integer(DP), intent(inout) :: array(n)
    integer(DP) :: i, j, k, k1, kn2, kn, tmp
    integer, parameter :: threshold = 4

    if(n <= 1) return

    k1  = array(1)
    kn2 = array(n / 2)
    kn  = array(n)
    if(k1 < kn2) then
       if(kn2 < kn) then
          k = kn2 ! 1, n / 2, n
       else ! array(n) <= array(n / 2)
          if(k1 < kn) then
             k = kn ! 1, n, n / 2
          else
             k = k1 ! n, 1, n / 2
          end if
       end if
    else ! array(n / 2) <= array(1)
       if(k1 < kn) then
          k = k1 ! n / 2, 1, n
       else ! array(n) <= array(1)
          if(kn2 < kn) then
             k = kn ! n / 2, n, 1
          else ! array(n) <= array(n / 2)
             k = kn2 ! n, n / 2, 1
          end if
       end if
    end if

    i = 1
    j = n
    do
       do i = i, n, 1
          if(array(i) .ge. k) exit
       end do
       do j = j, 1, -1
          if(array(j) .le. k) exit
       end do
       if(i >= j) exit
       tmp = array(i)
       array(i) = array(j)
       array(j) = tmp
       i = i + 1
       j = j - 1
    end do

    if(i - 1 > threshold) then
       call SORT_quicksort(i - 1, array(1:(i - 1)))
    else if(i - 1 > 1) then
       call SORT_heapsort(i - 1, array(1:(i - 1)))
    end if
    if(n - j > threshold) then
       call SORT_quicksort(n - j, array((j + 1):n))
    else if(n - j > 1) then
       call SORT_heapsort(n - j, array((j + 1):n))
    end if

    return
  end subroutine SORT_quicksort

  subroutine SORT_heapsort(n, array)
    implicit none
    integer(DP), intent(in) :: n
    integer(DP), intent(inout) :: array(n)

    integer(DP) :: i, j, k, l
    integer(DP) :: t

    l = n / 2 + 1
    k = n
    do while(k .ne. 1)
       if(l .gt. 1) then
          l = l - 1
          t = array(l)
       else
          t = array(k)
          array(k) = array(1)
          k = k - 1
          if(k .eq. 1) then
             array(1) = t
             exit
          end if
       end if
       i = l
       j = l + l
       do while(j .le. k)
          if(j .lt. k)then
             if(array(j) .lt. array(j + 1)) j = j + 1
          end if
          if(t .lt. array(j))then
             array(i) = array(j)
             i = j
             j = j + j
          else
             j = k + 1
          endif
       enddo
       array(i) = t
    end do

    return
  end subroutine SORT_heapsort

  recursive subroutine SORT_quickselect(A, left, right, K)
    implicit none
    real(RP), intent(inout) :: A(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer :: middle, pivot

    if (left < right) then
      if ((right-left)/K >= 2) then
        call sample_second_min(A, left, right, K, pivot)
      else
        middle = (left + right) / 2
        call median_of_three(A, left, middle, right, pivot)
      end if
      call partition(A, left, right, pivot)
      if (K < pivot) then
        call SORT_quickselect(A, left, pivot-1, K)
      else if (K > pivot) then
        call SORT_quickselect(A, pivot+1, right, K)
      end if
    end if

    return
  end subroutine SORT_quickselect

  recursive subroutine SORT_quickselect_arg(A, X, left, right, K)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(inout) :: X(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer :: middle, pivot

    if (left < right) then
      if ((right-left)/K >= 2) then
        call sample_second_min_arg(A, X, left, right, K, pivot)
      else
        middle = (left + right) / 2
        call median_of_three_arg(A, X, left, middle, right, pivot)
      end if
      call partition_arg(A, X, left, right, pivot)
      if (K < pivot) then
        call SORT_quickselect_arg(A, X, left, pivot-1, K)
      else if (K > pivot) then
        call SORT_quickselect_arg(A, X, pivot+1, right, K)
      end if
    end if

    return
  end subroutine SORT_quickselect_arg

  recursive subroutine SORT_quickselect_desc(A, left, right, K)
    implicit none
    real(RP), intent(inout) :: A(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer :: middle, pivot

    if (left < right) then
      if ((right-left)/K >= 2) then
        call sample_second_max(A, left, right, K, pivot)
      else
        middle = (left + right) / 2
        call median_of_three(A, left, middle, right, pivot)
      end if
      call partition_desc(A, left, right, pivot)
      if (K < pivot) then
        call SORT_quickselect_desc(A, left, pivot-1, K)
      else if (K > pivot) then
        call SORT_quickselect_desc(A, pivot+1, right, K)
      end if
    end if

    return
  end subroutine SORT_quickselect_desc

  recursive subroutine SORT_quickselect_desc_arg(A, X, left, right, K)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(inout) :: X(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer :: middle, pivot

    if (left < right) then
      if ((right-left)/K >= 2) then
        call sample_second_max_arg(A, X, left, right, K, pivot)
      else
        middle = (left + right) / 2
        call median_of_three_arg(A, X, left, middle, right, pivot)
      end if
      call partition_desc_arg(A, X, left, right, pivot)
      if (K < pivot) then
        call SORT_quickselect_desc_arg(A, X, left, pivot-1, K)
      else if (K > pivot) then
        call SORT_quickselect_desc_arg(A, X, pivot+1, right, K)
      end if
    end if

    return
  end subroutine SORT_quickselect_desc_arg

  ! private

  subroutine swap(x, y)
    implicit none
    real(RP), intent(inout) :: x, y
    real(RP) :: tmp
    tmp = x
    x = y
    y = tmp

    return
  end subroutine swap

  subroutine swap_i(x, y)
    implicit none
    integer, intent(inout) :: x, y
    integer :: tmp
    tmp = x
    x = y
    y = tmp

    return
  end subroutine swap_i

  subroutine partition(A, left, right, pivot)
    implicit none
    real(RP), intent(inout) :: A(:)
    integer, intent(in) :: left, right
    integer, intent(inout) :: pivot
    real(RP) :: A_pivot
    integer :: idx, store_idx

    A_pivot = A(pivot)
    call swap(A(pivot), A(right))

    store_idx = left
    do idx = left, right-1
      if (A(idx) < A_pivot) then
        call swap(A(store_idx), A(idx))
        store_idx = store_idx + 1
      end if
    end do

    call swap(A(right), A(store_idx))
    pivot = store_idx

    return
  end subroutine partition

  subroutine partition_arg(A, X, left, right, pivot)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(inout) :: X(:)
    integer, intent(in) :: left, right
    integer, intent(inout) :: pivot
    real(RP) :: A_pivot
    integer :: idx, store_idx

    A_pivot = A(X(pivot))
    call swap_i(X(pivot), X(right))

    store_idx = left
    do idx = left, right-1
      if (A(X(idx)) < A_pivot) then
        call swap_i(X(store_idx), X(idx))
        store_idx = store_idx + 1
      end if
    end do

    call swap_i(X(right), X(store_idx))
    pivot = store_idx

    return
  end subroutine partition_arg

  subroutine partition_desc(A, left, right, pivot)
    implicit none
    real(RP), intent(inout) :: A(:)
    integer, intent(in) :: left, right
    integer, intent(inout) :: pivot
    real(RP) :: A_pivot
    integer :: idx, store_idx

    A_pivot = A(pivot)
    call swap(A(pivot), A(right))

    store_idx = left
    do idx = left, right-1
      if (A(idx) > A_pivot) then
        call swap(A(store_idx), A(idx))
        store_idx = store_idx + 1
      end if
    end do

    call swap(A(right), A(store_idx))
    pivot = store_idx

    return
  end subroutine partition_desc

  subroutine partition_desc_arg(A, X, left, right, pivot)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(inout) :: X(:)
    integer, intent(in) :: left, right
    integer, intent(inout) :: pivot
    real(RP) :: A_pivot
    integer :: idx, store_idx

    A_pivot = A(X(pivot))
    call swap_i(X(pivot), X(right))

    store_idx = left
    do idx = left, right-1
      if (A(X(idx)) > A_pivot) then
        call swap_i(X(store_idx), X(idx))
        store_idx = store_idx + 1
      end if
    end do

    call swap_i(X(right), X(store_idx))
    pivot = store_idx

    return
  end subroutine partition_desc_arg

  subroutine median_of_three(A, i1, i2, i3, i)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(in) :: i1, i2, i3
    integer, intent(out) :: i

    if (A(i1) < A(i2)) then
      if (A(i2) < A(i3)) then
        i = i2
      else if (A(i1) < A(i3)) then
        i = i3
      else
        i = i1
      end if
    else
      if (A(i1) < A(i3)) then
        i = i1
      else if (A(i2) < A(i3)) then
        i = i3
      else
        i = i2
      end if
    end if

    return
  end subroutine median_of_three

  subroutine median_of_three_arg(A, X, i1, i2, i3, i)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(in) :: X(:)
    integer, intent(in) :: i1, i2, i3
    integer, intent(out) :: i

    if (A(X(i1)) < A(X(i2))) then
      if (A(X(i2)) < A(X(i3))) then
        i = i2
      else if (A(X(i1)) < A(X(i3))) then
        i = i3
      else
        i = i1
      end if
    else
      if (A(X(i1)) < A(X(i3))) then
        i = i1
      else if (A(X(i2)) < A(X(i3))) then
        i = i3
      else
        i = i2
      end if
    end if

    return
  end subroutine median_of_three_arg

  subroutine sample_second_min(A, left, right, K, i_2)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer, intent(out) :: i_2
    real(RP) :: A_min, A_min_2
    integer :: i, j

    i = left
    i_2 = left
    A_min = huge(A)
    A_min_2 = huge(A)
    do j = left, right, K
      if (A(j) < A_min) then
        A_min_2 = A_min
        A_min = A(j)
        i_2 = i
        i = j
      else if (A(j) < A_min_2) then
        A_min_2 = A(j)
        i_2 = j
      end if
    end do

    return
  end subroutine sample_second_min

  subroutine sample_second_min_arg(A, X, left, right, K, i_2)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(in) :: X(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer, intent(out) :: i_2
    real(RP) :: A_min, A_min_2
    integer :: i, j

    i = left
    i_2 = left
    A_min = huge(A)
    A_min_2 = huge(A)
    do j = left, right, K
      if (A(X(j)) < A_min) then
        A_min_2 = A_min
        A_min = A(X(j))
        i_2 = i
        i = j
      else if (A(X(j)) < A_min_2) then
        A_min_2 = A(X(j))
        i_2 = j
      end if
    end do

    return
  end subroutine sample_second_min_arg

  subroutine sample_second_max(A, left, right, K, i_2)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer, intent(out) :: i_2
    real(RP) :: A_max, A_max_2
    integer :: i, j

    i = left
    i_2 = left
    A_max = 0.0d0 - huge(A)
    A_max_2 = 0.0d0 - huge(A)
    do j = left, right, K
      if (A(j) > A_max) then
        A_max_2 = A_max
        A_max = A(j)
        i_2 = i
        i = j
      else if (A(j) > A_max_2) then
        A_max_2 = A(j)
        i_2 = j
      end if
    end do

    return
  end subroutine sample_second_max

  subroutine sample_second_max_arg(A, X, left, right, K, i_2)
    implicit none
    real(RP), intent(in) :: A(:)
    integer, intent(in) :: X(:)
    integer, intent(in) :: left, right
    integer, intent(in) :: K
    integer, intent(out) :: i_2
    real(RP) :: A_max, A_max_2
    integer :: i, j

    i = left
    i_2 = left
    A_max = 0.0d0 - huge(A)
    A_max_2 = 0.0d0 - huge(A)
    do j = left, right, K
      if (A(X(j)) > A_max) then
        A_max_2 = A_max
        A_max = A(X(j))
        i_2 = i
        i = j
      else if (A(X(j)) > A_max_2) then
        A_max_2 = A(X(j))
        i_2 = j
      end if
    end do

    return
  end subroutine sample_second_max_arg

end module scale_sort
