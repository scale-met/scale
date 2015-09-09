/*
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
*/
#ifndef GTYPE
  COMPILE_ERROR "'GTYPE' is not defined"
#endif
#ifndef GADG_ACCUMULATE_
  COMPILE_ERROR "'GADG_ACCUMULATE_' is not defined"
#endif
#ifndef I_ACCUMULATE_SCALAR
  COMPILE_ERROR "'I_ACCUMULATE_SCALAR' is not defined"
#endif
#ifndef I_ACCUMULATE_VECTOR
  COMPILE_ERROR "'I_ACCUMULATE_VECTOR' is not defined"
#endif

  !###########################################################################
  subroutine GADG_ACCUMULATE_(a, rslt)
    ! Accumulation, rslt(i) = a(1) + a(2) + ... + a(i).
    
    GTYPE, intent(inout) :: a(1:)
      ! Input array. If rslt is omitted,
      ! A is overwritten by the result on exit.
      ! Otherwise A is not changed.
    GTYPE, intent(out), optional :: rslt(1:)
      ! Result array. size(rslt) must be at least size(a).
      ! When rslt is omitted, the results are stored in a(:)
      
    !=------------------------------------------------------------------------
    if (present(rslt)) then
      !#=confirm(argument error) [size(a)] <= [size(rslt)]
        !#=begin
          if (.not.(size(a) <= size(rslt))) then
            print *, 'size(a) = ', size(a)
            print *, 'size(rslt) = ', size(rslt)
            print *, 'Required: size(a) <= size(rslt)'
            call gadg_abort('gadg_accumulate_', 'argument error')
          end if
        !#=end
    end if
    if (size(a) == 0) return
    
    if (GADG_CONF_VECTORIZE) then
      call I_ACCUMULATE_VECTOR(a, rslt)
    else
      call I_ACCUMULATE_SCALAR(a, rslt)
    end if
  end subroutine GADG_ACCUMULATE_
  
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine I_ACCUMULATE_SCALAR(a, rslt)
    ! Accumulation for scalar machine
    
    GTYPE, intent(inout) :: a(1:)
    GTYPE, intent(out), optional :: rslt(1:)
    !=------------------------------------------------------------------------
    integer :: i
    GTYPE :: tmp
    
    if (GADG_COMPILER_DEPENDENT_A) then
      if (present(rslt)) then
        rslt(1) = a(1)
        do i = 2, size(a)
          rslt(i) = rslt(i-1) + a(i)
        end do
      else
        do i = 2, size(a)
          a(i) = a(i) + a(i-1)
        end do
      end if
    else
      tmp = a(1)
      if (present(rslt)) then
        rslt(1) = tmp
        do i = 2, size(a)
          tmp = tmp + a(i)
          rslt(i) = tmp
        end do
      else
        do i = 2, size(a)
          tmp = tmp + a(i)
          a(i) = tmp
        end do
      end if
    end if
  end subroutine I_ACCUMULATE_SCALAR
  
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine I_ACCUMULATE_VECTOR(a, rslt)
    ! Accumulation for vector machine
    
    GTYPE, intent(inout) :: a(1:)
    GTYPE, intent(out), optional :: rslt(1:)
    
    !=------------------------------------------------------------------------
    integer :: i, j, size_a
    integer, parameter :: w = GADG_CONF_VEC_REG_LEN ! matrix width
    integer :: h  ! matrix height
    GTYPE :: s(GADG_CONF_VEC_REG_LEN)
    GTYPE, pointer :: r(:)
    logical :: present_rslt
    target :: a, rslt
!CDIR VREG(s)
    
    size_a = size(a)
    h = gadg_vmat_height(size_a)
    if (h <= 5) then
      call I_ACCUMULATE_SCALAR(a, rslt)
    else
      present_rslt = present(rslt)
      
      s(1 : w) = a(1 : (w - 1) * h + 1 : h)
      
      r => a(:)
      if (present_rslt) then
        r => rslt(:)
        rslt(1 : (w - 1) * h + 1 : h) = s(1 : w)
      end if
      do i = 2, h
        do j = 1, w
          s(j) = s(j) + a((j - 1) * h + i)
        end do
        do j = 1, w
          r((j - 1) * h + i) = s(j)
        end do
      end do
      
      !*** Accumulate the bottom row
      do j = 2, w
        s(j) = s(j) + s(j-1)
      end do
        
      !*** Consummation 
      do i = 1, h
        do j = 1, w - 1
          r(j * h + i) = r(j * h + i) + s(j)
        end do
      end do
      if (present_rslt) then
        do i = w * h + 1, size_a
          rslt(i) = a(i) + rslt(i - 1)
        end do
      else
        do i = w * h + 1, size_a
          a(i) = a(i) + a(i - 1)
        end do
      end if
    end if
  end subroutine I_ACCUMULATE_VECTOR
  
#undef GTYPE
#undef GADG_ACCUMULATE_
#undef I_ACCUMULATE_SCALAR
#undef I_ACCUMULATE_VECTOR
