! Auther:: Akio Kawano
! Copyright 2006 Akio Kawano. All rights reserved.

! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer. 
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! 3. All advertising materials mentioning features or use of this software must
! display the following acknowledgement:
!
! The names of its contributors may be used to endorse or promote products
! derived from this software without specific prior written permission. 
!
! THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
! DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "gadg_config.inc"
#include "gadg_macros.inc"

!module m_gadg_algorithm
!  use m_gadg_base
module gadg_algorithm
  use gadg_base
  implicit none
  private
  
  interface gadg_accumulate
    module procedure  gadg_acc_i32
    module procedure  gadg_acc_i64
    module procedure  gadg_acc_sp
    module procedure  gadg_acc_dp
  end interface
  
  public :: gadg_count_sort
  public :: gadg_vmat_height
  public :: gadg_accumulate
  
contains

  !###########################################################################
  subroutine gadg_count_sort(key, min_key, max_key, freq, tag, sorted, &
                             local_work, stat)
    ! counting sorting
    
    integer, intent(in)  :: key(1:)         ! sort keys
    integer, intent(in)  :: min_key         ! min(key)
    integer, intent(in)  :: max_key         ! max(key)
    integer, intent(out) :: freq(min_key:)
      ! freq(i) is count(key(:) == i)
      ! ubound(freq,1) must be at least max_key.
    integer, intent(out) :: tag(min_key:)
      ! tag(i) is the smallest of k such that key(sorted(k)) == i.
      ! ubound(tag,1) must be at least max_key.
    integer, intent(out) :: sorted(1:)
      ! Sorted indeces. key(sorted(i)) <= key(sorted(j))
      ! for all 1 <= i <= j <= size(key).
      ! size(sorted) is at least size(key).
    logical, intent(in), optional :: local_work
      ! If .true. then allocatable local array is used as the workarea.
      ! If .false. or not present, GADG_WORK_I is used as the workarea.
    integer, intent(out), optional :: stat
      ! State indicator for memory allocation/deallocation.
      ! Returns a non-zero value if memory allocation/deallocation fails.
      ! If not specified and an allocation/deallocation error condition
      ! occurs, program execution terminates.
    
    !=------------------------------------------------------------------------
    
    !#=verify [min_key] <= [max_key]
    /*||*/  if (.not.(min_key <= max_key)) then
    /*||*/    print *, 'min_key = ', min_key
    /*||*/    print *, 'max_key = ', max_key
    /*||*/    print *, 'Required: min_key <= max_key'
    /*||*/    call gadg_abort('m_gadg_algorithm::gadg_count_sort')
    /*||*/  end if
    
    !#=verify [ubound(freq, 1)] >= [max_key]
    /*||*/  if (.not.(ubound(freq, 1) >= max_key)) then
    /*||*/    print *, 'ubound(freq, 1) = ', ubound(freq, 1)
    /*||*/    print *, 'max_key = ', max_key
    /*||*/    print *, 'Required: ubound(freq, 1) >= max_key'
    /*||*/    call gadg_abort('m_gadg_algorithm::gadg_count_sort')
    /*||*/  end if
    
    !#=verify [ubound(tag, 1)] >= [max_key]
    /*||*/  if (.not.(ubound(tag, 1) >= max_key)) then
    /*||*/    print *, 'ubound(tag, 1) = ', ubound(tag, 1)
    /*||*/    print *, 'max_key = ', max_key
    /*||*/    print *, 'Required: ubound(tag, 1) >= max_key'
    /*||*/    call gadg_abort('m_gadg_algorithm::gadg_count_sort')
    /*||*/  end if
    
    !#=verify [size(sorted)] >= [size(key)]
    /*||*/  if (.not.(size(sorted) >= size(key))) then
    /*||*/    print *, 'size(sorted) = ', size(sorted)
    /*||*/    print *, 'size(key) = ', size(key)
    /*||*/    print *, 'Required: size(sorted) >= size(key)'
    /*||*/    call gadg_abort('m_gadg_algorithm::gadg_count_sort')
    /*||*/  end if
    
    if (GADG_CONF_VECTORIZE) then
      call i_csort_v(key, min_key, max_key, freq, tag, sorted, &
      &              local_work, stat)
    else
      call i_csort_s(key, min_key, max_key, freq, tag, sorted)
      if (present(stat)) stat = 0
    end if
  end subroutine gadg_count_sort
  
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine i_csort_s(key, min_key, max_key, freq, tag, sorted)
    ! distribution counting sorting for scalar processors
    
    integer, intent(in)  :: key(1:)
    integer, intent(in)  :: min_key
    integer, intent(in)  :: max_key
    integer, intent(out) :: freq(min_key:)
    integer, intent(out) :: tag(min_key:)
    integer, intent(out) :: sorted(1:)
    
    !=------------------------------------------------------------------------
    integer :: i, k, t, n_element
    
    n_element = size(key)
    
    !*** Count the occurrences of keys
    freq(min_key : max_key) = 0
    do i = 1, n_element
      k = key(i)
      freq(k) = freq(k) + 1
    end do
    
    !*** Get cumulative total
    freq(min_key) = freq(min_key) + 1   ! add 1 for all tag(:)
    call gadg_accumulate(freq, tag)
    freq(min_key) = freq(min_key) - 1
    
    !*** Write the sorted indeces
    do i = n_element, 1, -1
      k = key(i)
      t = tag(k)
      t = t - 1
      tag(k) = t
      sorted(t) = i
    end do
  end subroutine i_csort_s
  
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine i_csort_v(key, min_key, max_key, freq, tag, sorted, &
  &                    local_work, stat)
    ! distribution counting sorting for vector processors
    
    integer, intent(in)  :: key(1:)
    integer, intent(in)  :: min_key
    integer, intent(in)  :: max_key
    integer, intent(out) :: freq(min_key:)
    integer, intent(out) :: tag(min_key:)
    integer, intent(out) :: sorted(1:)
    logical, intent(in), optional :: local_work
    integer, intent(out), optional :: stat
    
    !=------------------------------------------------------------------------
    integer :: n_key_kind, worksiz, n_element
    integer :: st
    integer :: i, j, k, h
    integer, parameter :: w = GADG_CONF_VEC_REG_LEN
    integer, parameter :: IFACTOR = GADG_CONF_VEC_REG_LEN + 1
    integer, pointer :: work(:)
    logical :: local_work_act
      ! add 1 to avoid bank conflict
    
    n_key_kind = max_key - min_key + 1
    n_element = size(key)
    worksiz = IFACTOR * n_key_kind
    
    h = gadg_vmat_height(n_element)
    if (h < 5) then
      call i_csort_s(key, min_key, max_key, freq, tag, sorted)
      return
    end if
    
    local_work_act = .false.
    if (present(local_work)) local_work_act = local_work
    
    !*** Reserve GADG_WORK_I
    if (local_work_act) then
      allocate(work(worksiz), stat=st)
      if (st /= 0) then
        if (present(stat)) stat = st
        return
      end if
    else
      call gadg_reserve_work_i(worksiz, st)
      if (st /= 0) then
        call gadg_dealloc_all_works
        call gadg_reserve_work_i(worksiz, st)
        if (st /= 0) then
          if (present(stat)) stat = st
          return
        end if
      end if
      work => GADG_WORK_I(1 : worksiz)
    end if
    
    !*** Count the occurrences of extended keys
    work(1 : worksiz) = 0
    do j = 1, h
!CDIR NODEP, SHORTLOOP
      do i = 1, w
        k = key((i - 1) * h + j)
        k = (k - min_key) * IFACTOR + i
        work(k) = work(k) + 1
      end do
    end do
    do j = w * h + 1, n_element
      k = key(j)
      k = (k - min_key + 1) * IFACTOR
      work(k) = work(k) + 1
    end do
    
    !*** Get cumulative total
    work(1) = work(1) + 1   ! add 1 for all work(:)
    call gadg_accumulate(work(1 : worksiz))
    
    !*** Get the sorted indeces
    do j = n_element, w * h + 1, -1
      k = key(j)
      k = (k - min_key + 1) * IFACTOR
      work(k) = work(k) - 1
      sorted(work(k)) = j
    end do
    do j = h, 1, -1
!CDIR NODEP, SHORTLOOP
      do i = 1, w
        k = key((i - 1) * h + j)
        k = (k - min_key) * IFACTOR + i
        work(k) = work(k) - 1
        sorted(work(k)) = (i - 1) * h + j
      end do
    end do

    ! Write tag(:)
    tag(min_key : max_key) = work(1 : worksiz : IFACTOR)
    
    ! Write freq(:)
!CDIR NODEP
    do i = min_key, max_key - 1
      freq(i) = tag(i + 1) - tag(i)
    end do
    freq(max_key) = n_element + 1 - tag(max_key)
    
    if (local_work_act) then
      deallocate(work)
    end if
    
  end subroutine i_csort_v

  !###########################################################################
  function gadg_vmat_height(size_a) result(h)
    ! return the following h;
    !
    ! x(1) x(h+1) .. x((w-1)*h+1) x(w*h+1)
    ! x(2) x(h+2) .. x((w-1)*h+2)  :
    !  :    :         :           x(size_a)
    ! x(h) x(h+h) .. x((w-1)*h+h)            ,
    !
    ! where w = GADG_CONF_VEC_REG_LEN and h must be an odd integer or zero.

    integer, intent(in) :: size_a
    integer :: h
    
    !=------------------------------------------------------------------------
    
    ! 1st. trial
    h = size_a / GADG_CONF_VEC_REG_LEN
    if (h > 0 .and. GADG_IS_EVEN(h)) h = h - 1
    
  end function gadg_vmat_height
  
#define GTYPE integer(int32)
#define GADG_ACCUMULATE_  gadg_acc_i32
#define I_ACCUMULATE_SCALAR i_acc_s_i32
#define I_ACCUMULATE_VECTOR i_acc_v_i32
#include "gadg_algorithm_accum.inc"

#define GTYPE integer(int64)
#define GADG_ACCUMULATE_  gadg_acc_i64
#define I_ACCUMULATE_SCALAR i_acc_s_i64
#define I_ACCUMULATE_VECTOR i_acc_v_i64
#include "gadg_algorithm_accum.inc"

#define GTYPE real(sp)
#define GADG_ACCUMULATE_  gadg_acc_sp
#define I_ACCUMULATE_SCALAR i_acc_s_sp
#define I_ACCUMULATE_VECTOR i_acc_v_sp
#include "gadg_algorithm_accum.inc"

#define GTYPE real(dp)
#define GADG_ACCUMULATE_  gadg_acc_dp
#define I_ACCUMULATE_SCALAR i_acc_s_dp
#define I_ACCUMULATE_VECTOR i_acc_v_dp
#include "gadg_algorithm_accum.inc"
!end module m_gadg_algorithm
end module gadg_algorithm
