! Auther:: Akio Kawano
! generator :: gadg_reduction.F90_fpl
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

!module m_gadg_reduction
!  use m_gadg_base
module gadg_reduction
  use gadg_base
  implicit none
  private
  interface gadg_reduction_minval
    module procedure rdct_minval_1, rdct_minval_2, rdct_minval_3, rdct_minval_4
  end interface
  
  interface gadg_reduction_maxval
    module procedure rdct_maxval_1, rdct_maxval_2, rdct_maxval_3, rdct_maxval_4
  end interface
  
  interface gadg_reduction_sum
    module procedure rdct_sum_1, rdct_sum_2, rdct_sum_3, rdct_sum_4
  end interface
  
  interface gadg_reduction_prod
    module procedure rdct_prod_1, rdct_prod_2, rdct_prod_3, rdct_prod_4
  end interface
  
  interface gadg_reduction_iand
    module procedure rdct_iand_1, rdct_iand_2
  end interface
  
  interface gadg_reduction_ior
    module procedure rdct_ior_1, rdct_ior_2
  end interface
  
  interface gadg_reduction_ieor
    module procedure rdct_ieor_1, rdct_ieor_2
  end interface
  
  interface gadg_reduction_minloc
    module procedure rdct_minloc_1, rdct_minloc_2, rdct_minloc_3, rdct_minloc_4
  end interface
  
  interface gadg_reduction_maxloc
    module procedure rdct_maxloc_1, rdct_maxloc_2, rdct_maxloc_3, rdct_maxloc_4
  end interface
  
  public :: gadg_reduction_minval
  public :: gadg_reduction_maxval
  public :: gadg_reduction_sum
  public :: gadg_reduction_prod
  public :: gadg_reduction_iand
  public :: gadg_reduction_ior
  public :: gadg_reduction_ieor
  public :: gadg_reduction_minloc
  public :: gadg_reduction_maxloc

contains
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_minval_1(a, rslt, dim_reduced, flush)
    ! rslt(i) = minval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = minval(rslt(i), minval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = minval(rslt(i), minval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minval_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minval_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_minval_1s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_minval_1v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_minval_1s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_minval_1v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minval_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minval_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minval_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minval_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minval_', 'internal error')
    end select
  end subroutine rdct_minval_1
  
#define GM_REDUCE_OP(R,A,Q) if (R > A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_1s1n
#define GM_SUBNAME_SPLIT    irdct_minval_1v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_1s2n
#define GM_SUBNAME_SPLIT    irdct_minval_1v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_1s1f
#define GM_SUBNAME_SPLIT    irdct_minval_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_1s2f
#define GM_SUBNAME_SPLIT    irdct_minval_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_minval_2(a, rslt, dim_reduced, flush)
    ! rslt(i) = minval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = minval(rslt(i), minval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = minval(rslt(i), minval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minval_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minval_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_minval_2s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_minval_2v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_minval_2s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_minval_2v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minval_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minval_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minval_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minval_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minval_', 'internal error')
    end select
  end subroutine rdct_minval_2
  
#define GM_REDUCE_OP(R,A,Q) if (R > A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_2s1n
#define GM_SUBNAME_SPLIT    irdct_minval_2v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_2s2n
#define GM_SUBNAME_SPLIT    irdct_minval_2v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_2s1f
#define GM_SUBNAME_SPLIT    irdct_minval_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_2s2f
#define GM_SUBNAME_SPLIT    irdct_minval_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(sp)
#define GTYPE2  real(sp)
  !###########################################################################
  subroutine rdct_minval_3(a, rslt, dim_reduced, flush)
    ! rslt(i) = minval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = minval(rslt(i), minval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = minval(rslt(i), minval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minval_3')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minval_3')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_minval_3s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_minval_3v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_minval_3s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_minval_3v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minval_3s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minval_3v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minval_3s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minval_3v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minval_', 'internal error')
    end select
  end subroutine rdct_minval_3
  
#define GM_REDUCE_OP(R,A,Q) if (R > A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_3s1n
#define GM_SUBNAME_SPLIT    irdct_minval_3v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_3s2n
#define GM_SUBNAME_SPLIT    irdct_minval_3v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_3s1f
#define GM_SUBNAME_SPLIT    irdct_minval_3v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_3s2f
#define GM_SUBNAME_SPLIT    irdct_minval_3v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(dp)
#define GTYPE2  real(dp)
  !###########################################################################
  subroutine rdct_minval_4(a, rslt, dim_reduced, flush)
    ! rslt(i) = minval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = minval(rslt(i), minval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = minval(rslt(i), minval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minval_4')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minval_4')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_minval_4s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_minval_4v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_minval_4s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_minval_4v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minval_4s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minval_4v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minval_4s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minval_4v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minval_', 'internal error')
    end select
  end subroutine rdct_minval_4
  
#define GM_REDUCE_OP(R,A,Q) if (R > A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_4s1n
#define GM_SUBNAME_SPLIT    irdct_minval_4v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_4s2n
#define GM_SUBNAME_SPLIT    irdct_minval_4v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minval_4s1f
#define GM_SUBNAME_SPLIT    irdct_minval_4v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minval_4s2f
#define GM_SUBNAME_SPLIT    irdct_minval_4v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_maxval_1(a, rslt, dim_reduced, flush)
    ! rslt(i) = maxval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = maxval(rslt(i), maxval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = maxval(rslt(i), maxval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxval_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxval_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_maxval_1s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_maxval_1v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_maxval_1s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_maxval_1v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxval_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxval_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxval_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxval_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxval_', 'internal error')
    end select
  end subroutine rdct_maxval_1
  
#define GM_REDUCE_OP(R,A,Q) if (R < A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_1s1n
#define GM_SUBNAME_SPLIT    irdct_maxval_1v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_1s2n
#define GM_SUBNAME_SPLIT    irdct_maxval_1v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_1s1f
#define GM_SUBNAME_SPLIT    irdct_maxval_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_1s2f
#define GM_SUBNAME_SPLIT    irdct_maxval_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_maxval_2(a, rslt, dim_reduced, flush)
    ! rslt(i) = maxval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = maxval(rslt(i), maxval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = maxval(rslt(i), maxval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxval_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxval_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_maxval_2s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_maxval_2v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_maxval_2s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_maxval_2v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxval_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxval_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxval_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxval_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxval_', 'internal error')
    end select
  end subroutine rdct_maxval_2
  
#define GM_REDUCE_OP(R,A,Q) if (R < A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_2s1n
#define GM_SUBNAME_SPLIT    irdct_maxval_2v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_2s2n
#define GM_SUBNAME_SPLIT    irdct_maxval_2v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_2s1f
#define GM_SUBNAME_SPLIT    irdct_maxval_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_2s2f
#define GM_SUBNAME_SPLIT    irdct_maxval_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(sp)
#define GTYPE2  real(sp)
  !###########################################################################
  subroutine rdct_maxval_3(a, rslt, dim_reduced, flush)
    ! rslt(i) = maxval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = maxval(rslt(i), maxval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = maxval(rslt(i), maxval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxval_3')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxval_3')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_maxval_3s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_maxval_3v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_maxval_3s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_maxval_3v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxval_3s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxval_3v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxval_3s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxval_3v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxval_', 'internal error')
    end select
  end subroutine rdct_maxval_3
  
#define GM_REDUCE_OP(R,A,Q) if (R < A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_3s1n
#define GM_SUBNAME_SPLIT    irdct_maxval_3v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_3s2n
#define GM_SUBNAME_SPLIT    irdct_maxval_3v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_3s1f
#define GM_SUBNAME_SPLIT    irdct_maxval_3v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_3s2f
#define GM_SUBNAME_SPLIT    irdct_maxval_3v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(dp)
#define GTYPE2  real(dp)
  !###########################################################################
  subroutine rdct_maxval_4(a, rslt, dim_reduced, flush)
    ! rslt(i) = maxval(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxval(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = maxval(rslt(i), maxval(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = maxval(rslt(i), maxval(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxval_4')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxval_4')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_maxval_4s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_maxval_4v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_maxval_4s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_maxval_4v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxval_4s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxval_4v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxval_4s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxval_4v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxval_', 'internal error')
    end select
  end subroutine rdct_maxval_4
  
#define GM_REDUCE_OP(R,A,Q) if (R < A) R = A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_4s1n
#define GM_SUBNAME_SPLIT    irdct_maxval_4v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_4s2n
#define GM_SUBNAME_SPLIT    irdct_maxval_4v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxval_4s1f
#define GM_SUBNAME_SPLIT    irdct_maxval_4v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxval_4s2f
#define GM_SUBNAME_SPLIT    irdct_maxval_4v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_sum_1(a, rslt, dim_reduced, flush)
    ! rslt(i) = sum(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = sum(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = sum(rslt(i), sum(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = sum(rslt(i), sum(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_sum_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_sum_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_sum_1s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_sum_1v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_sum_1s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_sum_1v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_sum_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_sum_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_sum_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_sum_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_sum_', 'internal error')
    end select
  end subroutine rdct_sum_1
  
#define GM_REDUCE_OP(R,A,Q) R = R + A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_1s1n
#define GM_SUBNAME_SPLIT    irdct_sum_1v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_1s2n
#define GM_SUBNAME_SPLIT    irdct_sum_1v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_1s1f
#define GM_SUBNAME_SPLIT    irdct_sum_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_1s2f
#define GM_SUBNAME_SPLIT    irdct_sum_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_sum_2(a, rslt, dim_reduced, flush)
    ! rslt(i) = sum(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = sum(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = sum(rslt(i), sum(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = sum(rslt(i), sum(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_sum_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_sum_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_sum_2s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_sum_2v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_sum_2s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_sum_2v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_sum_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_sum_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_sum_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_sum_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_sum_', 'internal error')
    end select
  end subroutine rdct_sum_2
  
#define GM_REDUCE_OP(R,A,Q) R = R + A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_2s1n
#define GM_SUBNAME_SPLIT    irdct_sum_2v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_2s2n
#define GM_SUBNAME_SPLIT    irdct_sum_2v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_2s1f
#define GM_SUBNAME_SPLIT    irdct_sum_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_2s2f
#define GM_SUBNAME_SPLIT    irdct_sum_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(sp)
#define GTYPE2  real(sp)
  !###########################################################################
  subroutine rdct_sum_3(a, rslt, dim_reduced, flush)
    ! rslt(i) = sum(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = sum(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = sum(rslt(i), sum(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = sum(rslt(i), sum(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_sum_3')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_sum_3')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_sum_3s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_sum_3v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_sum_3s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_sum_3v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_sum_3s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_sum_3v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_sum_3s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_sum_3v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_sum_', 'internal error')
    end select
  end subroutine rdct_sum_3
  
#define GM_REDUCE_OP(R,A,Q) R = R + A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_3s1n
#define GM_SUBNAME_SPLIT    irdct_sum_3v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_3s2n
#define GM_SUBNAME_SPLIT    irdct_sum_3v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_3s1f
#define GM_SUBNAME_SPLIT    irdct_sum_3v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_3s2f
#define GM_SUBNAME_SPLIT    irdct_sum_3v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(dp)
#define GTYPE2  real(dp)
  !###########################################################################
  subroutine rdct_sum_4(a, rslt, dim_reduced, flush)
    ! rslt(i) = sum(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = sum(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = sum(rslt(i), sum(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = sum(rslt(i), sum(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_sum_4')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_sum_4')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_sum_4s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_sum_4v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_sum_4s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_sum_4v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_sum_4s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_sum_4v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_sum_4s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_sum_4v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_sum_', 'internal error')
    end select
  end subroutine rdct_sum_4
  
#define GM_REDUCE_OP(R,A,Q) R = R + A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_4s1n
#define GM_SUBNAME_SPLIT    irdct_sum_4v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_4s2n
#define GM_SUBNAME_SPLIT    irdct_sum_4v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_sum_4s1f
#define GM_SUBNAME_SPLIT    irdct_sum_4v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_sum_4s2f
#define GM_SUBNAME_SPLIT    irdct_sum_4v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_prod_1(a, rslt, dim_reduced, flush)
    ! rslt(i) = prod(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = prod(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = prod(rslt(i), prod(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = prod(rslt(i), prod(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_prod_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_prod_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_prod_1s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_prod_1v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_prod_1s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_prod_1v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_prod_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_prod_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_prod_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_prod_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_prod_', 'internal error')
    end select
  end subroutine rdct_prod_1
  
#define GM_REDUCE_OP(R,A,Q) R = R * A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_1s1n
#define GM_SUBNAME_SPLIT    irdct_prod_1v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_1s2n
#define GM_SUBNAME_SPLIT    irdct_prod_1v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_1s1f
#define GM_SUBNAME_SPLIT    irdct_prod_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_1s2f
#define GM_SUBNAME_SPLIT    irdct_prod_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_prod_2(a, rslt, dim_reduced, flush)
    ! rslt(i) = prod(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = prod(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = prod(rslt(i), prod(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = prod(rslt(i), prod(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_prod_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_prod_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_prod_2s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_prod_2v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_prod_2s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_prod_2v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_prod_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_prod_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_prod_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_prod_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_prod_', 'internal error')
    end select
  end subroutine rdct_prod_2
  
#define GM_REDUCE_OP(R,A,Q) R = R * A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_2s1n
#define GM_SUBNAME_SPLIT    irdct_prod_2v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_2s2n
#define GM_SUBNAME_SPLIT    irdct_prod_2v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_2s1f
#define GM_SUBNAME_SPLIT    irdct_prod_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_2s2f
#define GM_SUBNAME_SPLIT    irdct_prod_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(sp)
#define GTYPE2  real(sp)
  !###########################################################################
  subroutine rdct_prod_3(a, rslt, dim_reduced, flush)
    ! rslt(i) = prod(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = prod(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = prod(rslt(i), prod(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = prod(rslt(i), prod(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_prod_3')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_prod_3')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_prod_3s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_prod_3v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_prod_3s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_prod_3v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_prod_3s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_prod_3v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_prod_3s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_prod_3v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_prod_', 'internal error')
    end select
  end subroutine rdct_prod_3
  
#define GM_REDUCE_OP(R,A,Q) R = R * A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_3s1n
#define GM_SUBNAME_SPLIT    irdct_prod_3v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_3s2n
#define GM_SUBNAME_SPLIT    irdct_prod_3v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_3s1f
#define GM_SUBNAME_SPLIT    irdct_prod_3v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_3s2f
#define GM_SUBNAME_SPLIT    irdct_prod_3v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(dp)
#define GTYPE2  real(dp)
  !###########################################################################
  subroutine rdct_prod_4(a, rslt, dim_reduced, flush)
    ! rslt(i) = prod(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = prod(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = prod(rslt(i), prod(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = prod(rslt(i), prod(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_prod_4')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_prod_4')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_prod_4s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_prod_4v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_prod_4s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_prod_4v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_prod_4s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_prod_4v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_prod_4s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_prod_4v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_prod_', 'internal error')
    end select
  end subroutine rdct_prod_4
  
#define GM_REDUCE_OP(R,A,Q) R = R * A

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_4s1n
#define GM_SUBNAME_SPLIT    irdct_prod_4v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_4s2n
#define GM_SUBNAME_SPLIT    irdct_prod_4v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_prod_4s1f
#define GM_SUBNAME_SPLIT    irdct_prod_4v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_prod_4s2f
#define GM_SUBNAME_SPLIT    irdct_prod_4v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_iand_1(a, rslt, dim_reduced, flush)
    ! rslt(i) = iand(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = iand(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = iand(rslt(i), iand(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = iand(rslt(i), iand(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_iand_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_iand_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_iand_1s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_iand_1v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_iand_1s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_iand_1v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_iand_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_iand_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_iand_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_iand_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_iand_', 'internal error')
    end select
  end subroutine rdct_iand_1
  
#define GM_REDUCE_OP(R,A,Q) R = iand(R, A)

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_iand_1s1n
#define GM_SUBNAME_SPLIT    irdct_iand_1v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_iand_1s2n
#define GM_SUBNAME_SPLIT    irdct_iand_1v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_iand_1s1f
#define GM_SUBNAME_SPLIT    irdct_iand_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_iand_1s2f
#define GM_SUBNAME_SPLIT    irdct_iand_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_iand_2(a, rslt, dim_reduced, flush)
    ! rslt(i) = iand(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = iand(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = iand(rslt(i), iand(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = iand(rslt(i), iand(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_iand_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_iand_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_iand_2s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_iand_2v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_iand_2s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_iand_2v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_iand_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_iand_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_iand_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_iand_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_iand_', 'internal error')
    end select
  end subroutine rdct_iand_2
  
#define GM_REDUCE_OP(R,A,Q) R = iand(R, A)

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_iand_2s1n
#define GM_SUBNAME_SPLIT    irdct_iand_2v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_iand_2s2n
#define GM_SUBNAME_SPLIT    irdct_iand_2v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_iand_2s1f
#define GM_SUBNAME_SPLIT    irdct_iand_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_iand_2s2f
#define GM_SUBNAME_SPLIT    irdct_iand_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_ior_1(a, rslt, dim_reduced, flush)
    ! rslt(i) = ior(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = ior(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = ior(rslt(i), ior(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = ior(rslt(i), ior(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_ior_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_ior_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_ior_1s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_ior_1v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_ior_1s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_ior_1v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_ior_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_ior_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_ior_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_ior_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_ior_', 'internal error')
    end select
  end subroutine rdct_ior_1
  
#define GM_REDUCE_OP(R,A,Q) R = ior(R, A)

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ior_1s1n
#define GM_SUBNAME_SPLIT    irdct_ior_1v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ior_1s2n
#define GM_SUBNAME_SPLIT    irdct_ior_1v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ior_1s1f
#define GM_SUBNAME_SPLIT    irdct_ior_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ior_1s2f
#define GM_SUBNAME_SPLIT    irdct_ior_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_ior_2(a, rslt, dim_reduced, flush)
    ! rslt(i) = ior(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = ior(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = ior(rslt(i), ior(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = ior(rslt(i), ior(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_ior_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_ior_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_ior_2s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_ior_2v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_ior_2s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_ior_2v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_ior_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_ior_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_ior_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_ior_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_ior_', 'internal error')
    end select
  end subroutine rdct_ior_2
  
#define GM_REDUCE_OP(R,A,Q) R = ior(R, A)

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ior_2s1n
#define GM_SUBNAME_SPLIT    irdct_ior_2v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ior_2s2n
#define GM_SUBNAME_SPLIT    irdct_ior_2v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ior_2s1f
#define GM_SUBNAME_SPLIT    irdct_ior_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ior_2s2f
#define GM_SUBNAME_SPLIT    irdct_ior_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_ieor_1(a, rslt, dim_reduced, flush)
    ! rslt(i) = ieor(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = ieor(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = ieor(rslt(i), ieor(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = ieor(rslt(i), ieor(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_ieor_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_ieor_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_ieor_1s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_ieor_1v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_ieor_1s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_ieor_1v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_ieor_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_ieor_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_ieor_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_ieor_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_ieor_', 'internal error')
    end select
  end subroutine rdct_ieor_1
  
#define GM_REDUCE_OP(R,A,Q) R = ieor(R, A)

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ieor_1s1n
#define GM_SUBNAME_SPLIT    irdct_ieor_1v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ieor_1s2n
#define GM_SUBNAME_SPLIT    irdct_ieor_1v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ieor_1s1f
#define GM_SUBNAME_SPLIT    irdct_ieor_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ieor_1s2f
#define GM_SUBNAME_SPLIT    irdct_ieor_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_ieor_2(a, rslt, dim_reduced, flush)
    ! rslt(i) = ieor(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = ieor(a(i,:))                 if dim_reduced == 2
    ! rslt(i) = ieor(rslt(i), ieor(a(:,i))) if dim_reduced == 1 .and. flush == .false.
    ! rslt(i) = ieor(rslt(i), ieor(a(i,:))) if dim_reduced == 2 .and. flush == .false.
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    logical, intent(in), optional :: flush
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    if (present(flush)) flush_act = flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_ieor_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_ieor_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (flush_act) then
        if (dim_reduced == 1) then
          rslt(1 : iend_scan) = a(1, 1 : iend_scan)
        else
          rslt(1 : iend_scan) = a(1 : iend_scan, 1)
        end if
      end if
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_no_flush)   ! 0
      call irdct_ieor_2s1n(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_no_flush)   ! 1
      call irdct_ieor_2v1n(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_no_flush)   ! 2
      call irdct_ieor_2s2n(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_no_flush)   ! 3
      call irdct_ieor_2v2n(a, rslt)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_ieor_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_ieor_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_ieor_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_ieor_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_ieor_', 'internal error')
    end select
  end subroutine rdct_ieor_2
  
#define GM_REDUCE_OP(R,A,Q) R = ieor(R, A)

#undef GM_NEED_Q
#define GM_INIT_R           r = rslt(j)
#define GM_INIT_V           v(1 : iend_vloop) = rslt(k : k + iend_vloop - 1)
#define GM_START_I          1

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ieor_2s1n
#define GM_SUBNAME_SPLIT    irdct_ieor_2v1n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ieor_2s2n
#define GM_SUBNAME_SPLIT    irdct_ieor_2v2n
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GM_NEED_Q
#define GM_INIT_R           r = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_ieor_2s1f
#define GM_SUBNAME_SPLIT    irdct_ieor_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_ieor_2s2f
#define GM_SUBNAME_SPLIT    irdct_ieor_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_minloc_1(a, rslt, dim_reduced)
    ! rslt(i) = minloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minloc_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minloc_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minloc_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minloc_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minloc_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minloc_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minloc_', 'internal error')
    end select
  end subroutine rdct_minloc_1
  
#define GM_REDUCE_OP(R,A,Q) if (Q > A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minloc_1s1f
#define GM_SUBNAME_SPLIT    irdct_minloc_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minloc_1s2f
#define GM_SUBNAME_SPLIT    irdct_minloc_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_minloc_2(a, rslt, dim_reduced)
    ! rslt(i) = minloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minloc_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minloc_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minloc_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minloc_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minloc_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minloc_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minloc_', 'internal error')
    end select
  end subroutine rdct_minloc_2
  
#define GM_REDUCE_OP(R,A,Q) if (Q > A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minloc_2s1f
#define GM_SUBNAME_SPLIT    irdct_minloc_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minloc_2s2f
#define GM_SUBNAME_SPLIT    irdct_minloc_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(sp)
#define GTYPE2  real(sp)
  !###########################################################################
  subroutine rdct_minloc_3(a, rslt, dim_reduced)
    ! rslt(i) = minloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minloc_3')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minloc_3')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minloc_3s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minloc_3v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minloc_3s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minloc_3v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minloc_', 'internal error')
    end select
  end subroutine rdct_minloc_3
  
#define GM_REDUCE_OP(R,A,Q) if (Q > A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minloc_3s1f
#define GM_SUBNAME_SPLIT    irdct_minloc_3v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minloc_3s2f
#define GM_SUBNAME_SPLIT    irdct_minloc_3v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(dp)
#define GTYPE2  real(dp)
  !###########################################################################
  subroutine rdct_minloc_4(a, rslt, dim_reduced)
    ! rslt(i) = minloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = minloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_minloc_4')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_minloc_4')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_minloc_4s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_minloc_4v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_minloc_4s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_minloc_4v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_minloc_', 'internal error')
    end select
  end subroutine rdct_minloc_4
  
#define GM_REDUCE_OP(R,A,Q) if (Q > A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_minloc_4s1f
#define GM_SUBNAME_SPLIT    irdct_minloc_4v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_minloc_4s2f
#define GM_SUBNAME_SPLIT    irdct_minloc_4v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int32)
#define GTYPE2  integer(int32)
  !###########################################################################
  subroutine rdct_maxloc_1(a, rslt, dim_reduced)
    ! rslt(i) = maxloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxloc_1')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxloc_1')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxloc_1s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxloc_1v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxloc_1s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxloc_1v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxloc_', 'internal error')
    end select
  end subroutine rdct_maxloc_1
  
#define GM_REDUCE_OP(R,A,Q) if (Q < A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxloc_1s1f
#define GM_SUBNAME_SPLIT    irdct_maxloc_1v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxloc_1s2f
#define GM_SUBNAME_SPLIT    irdct_maxloc_1v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  integer(int64)
#define GTYPE2  integer(int64)
  !###########################################################################
  subroutine rdct_maxloc_2(a, rslt, dim_reduced)
    ! rslt(i) = maxloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxloc_2')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxloc_2')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxloc_2s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxloc_2v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxloc_2s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxloc_2v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxloc_', 'internal error')
    end select
  end subroutine rdct_maxloc_2
  
#define GM_REDUCE_OP(R,A,Q) if (Q < A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxloc_2s1f
#define GM_SUBNAME_SPLIT    irdct_maxloc_2v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxloc_2s2f
#define GM_SUBNAME_SPLIT    irdct_maxloc_2v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(sp)
#define GTYPE2  real(sp)
  !###########################################################################
  subroutine rdct_maxloc_3(a, rslt, dim_reduced)
    ! rslt(i) = maxloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxloc_3')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxloc_3')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxloc_3s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxloc_3v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxloc_3s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxloc_3v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxloc_', 'internal error')
    end select
  end subroutine rdct_maxloc_3
  
#define GM_REDUCE_OP(R,A,Q) if (Q < A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxloc_3s1f
#define GM_SUBNAME_SPLIT    irdct_maxloc_3v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxloc_3s2f
#define GM_SUBNAME_SPLIT    irdct_maxloc_3v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
#define GTYPE1  real(dp)
#define GTYPE2  real(dp)
  !###########################################################################
  subroutine rdct_maxloc_4(a, rslt, dim_reduced)
    ! rslt(i) = maxloc(a(:,i))                 if dim_reduced == 1
    ! rslt(i) = maxloc(a(i,:))                 if dim_reduced == 2
    
    GTYPE1, intent(in) :: a(1:, 1:)
    GTYPE2, intent(inout) :: rslt(1:)
    integer, intent(in) :: dim_reduced
    
    !=------------------------------------------------------------------------
    integer :: dim_scan
    integer :: iend_reduced, iend_scan
    integer :: routine_flags
    integer, parameter :: flag_sequential     = 0 ! use sequential method ?
    integer, parameter :: flag_split          = 1 ! use splitting method ?
    
    integer, parameter :: flag_dim_reduced_1  = 0 ! dim_reduced == 1 ?
    integer, parameter :: flag_dim_reduced_2  = 2 ! dim_reduced == 2 ?
    
    logical :: flush_act
    integer, parameter :: flag_no_flush      = 0 ! .not. flush_act ?
    integer, parameter :: flag_flush          = 4 ! flush_act ?
    
    flush_act = .true. ! default of flush
    
    dim_scan = 1
    if (dim_reduced == 1) dim_scan = 2
    
    !#=verify [dim_reduced] == 1 .or. dim_reduced == 2
    /*||*/  if (.not.(dim_reduced == 1 .or. dim_reduced == 2)) then
    /*||*/    print *, 'dim_reduced = ', dim_reduced
    /*||*/    print *, 'Required: dim_reduced == 1 .or. dim_reduced == 2'
    /*||*/    call gadg_abort('rdct_maxloc_4')
    /*||*/  end if
    
    !#=verify [size(a, dim_scan)] <= [size(rslt)]
    /*||*/  if (.not.(size(a, dim_scan) <= size(rslt))) then
    /*||*/    print *, 'size(a, dim_scan) = ', size(a, dim_scan)
    /*||*/    print *, 'size(rslt) = ', size(rslt)
    /*||*/    print *, 'Required: size(a, dim_scan) <= size(rslt)'
    /*||*/    call gadg_abort('rdct_maxloc_4')
    /*||*/  end if
    
    iend_reduced = ubound(a, dim_reduced)
    iend_scan = ubound(a, dim_scan)
    
    ! quick return
    if (iend_reduced == 1) then
      if (iend_reduced == 1) rslt(1 : iend_scan) = 1
      return
    end if
    
    routine_flags = merge(flag_flush, flag_no_flush, flush_act)
    routine_flags = routine_flags + &
      & merge(flag_dim_reduced_1, flag_dim_reduced_2, dim_reduced == 1)
    if (.not. GADG_CONF_VECTORIZE) then
      ! Scalar machine 
      routine_flags = routine_flags + flag_sequential
    else
      ! Vector machine: chooses a way with longer vector length.
      routine_flags = routine_flags + &
        & merge(flag_split, flag_sequential, iend_reduced < iend_scan * 4)
      ! Factor 4 ... macro operations have less performance than standard operations.
    end if
    
    select case (routine_flags)

    case (flag_sequential + flag_dim_reduced_1 + flag_flush)      ! 4
      call irdct_maxloc_4s1f(a, rslt)
    case (flag_split      + flag_dim_reduced_1 + flag_flush)      ! 5
      call irdct_maxloc_4v1f(a, rslt)
    case (flag_sequential + flag_dim_reduced_2 + flag_flush)      ! 6
      call irdct_maxloc_4s2f(a, rslt)
    case (flag_split      + flag_dim_reduced_2 + flag_flush)      ! 7
      call irdct_maxloc_4v2f(a, rslt)
    
    case default
      call gadg_abort('gadg_reduce_maxloc_', 'internal error')
    end select
  end subroutine rdct_maxloc_4
  
#define GM_REDUCE_OP(R,A,Q) if (Q < A) then; Q = A; R = i; end if

#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#define GM_NEED_Q
#define GM_INIT_R           r = 1; q = GM_A(1, j)
#define GM_INIT_V           v(1 : iend_vloop) = 1; q(1 : iend_vloop) = GM_A(1, k : k + iend_vloop - 1)
#define GM_START_I          2

#define GM_DIM_REDUCED      1
#define GM_SUBNAME_SEQ      irdct_maxloc_4s1f
#define GM_SUBNAME_SPLIT    irdct_maxloc_4v1f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#define GM_DIM_REDUCED      2
#define GM_SUBNAME_SEQ      irdct_maxloc_4s2f
#define GM_SUBNAME_SPLIT    irdct_maxloc_4v2f
#include "gadg_reduction_detail.inc"
#undef GM_DIM_REDUCED
#undef GM_SUBNAME_SEQ
#undef GM_SUBNAME_SPLIT

#undef GM_REDUCE_OP
#undef GM_NEED_Q
#undef GM_INIT_R
#undef GM_INIT_V
#undef GM_START_I

#undef GTYPE1
#undef GTYPE2
!end module m_gadg_reduction
end module gadg_reduction
