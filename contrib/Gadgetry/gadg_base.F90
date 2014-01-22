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

!module m_gadg_base
module gadg_base
  implicit none
  private
  
  ! Platform-specific constants
  
  integer, parameter :: GADG_CONF_VEC_REG_LEN = PP_CONF_VEC_REG_LEN
# ifdef PP_CONF_VECTORIZE
  logical, parameter :: GADG_CONF_VECTORIZE = .true.
# else
  logical, parameter :: GADG_CONF_VECTORIZE = .false.
# endif

# ifdef PP_CONF_USE_MPI
  logical, parameter :: GADG_CONF_USE_MPI = .true.
# else
  logical, parameter :: GADG_CONF_USE_MPI = .false.
# endif

# ifdef PP_CONF_64BIT_CPU
  logical, parameter :: GADG_CONF_64BIT_CPU = .true.
# else
  logical, parameter :: GADG_CONF_64BIT_CPU = .false.
# endif
  
#ifdef PP_COMPILER_DEPENDENT_A
  logical, parameter :: GADG_COMPILER_DEPENDENT_A = .true.
#else
  logical, parameter :: GADG_COMPILER_DEPENDENT_A = .false.
#endif
  
# ifdef NDEBUG
  logical, parameter :: GADG_NDEBUG = .true.
# else
  logical, parameter :: GADG_NDEBUG = .false.
#endif
  
  ! Kind parameters
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.D0)
  integer, parameter :: int16 = selected_int_kind(4)
  integer, parameter :: int32 = selected_int_kind(9)
  integer, parameter :: int64 = selected_int_kind(12)
# ifdef PP_CONF_64BIT_CPU
  integer, parameter :: intopt = int64
# else
  integer, parameter :: intopt = int32
# endif
  
  ! Work arrays
  integer(int32),   target, allocatable :: GADG_WORK_I32(:)
  integer(int64),   target, allocatable :: GADG_WORK_I64(:)
  integer,          target, allocatable :: GADG_WORK_I(:)
  integer(intopt),  target, allocatable :: GADG_WORK_IOPT(:)
  real(sp),         target, allocatable :: GADG_WORK_SP(:)
  real(dp),         target, allocatable :: GADG_WORK_DP(:)
  
  public :: sp, dp, int16, int32, int64, intopt
  public :: GADG_CONF_VEC_REG_LEN
  public :: GADG_CONF_VECTORIZE
  public :: GADG_CONF_USE_MPI
  public :: GADG_CONF_64BIT_CPU
  public :: GADG_COMPILER_DEPENDENT_A
  public :: GADG_NDEBUG
  
  public :: GADG_WORK_I32
  public :: GADG_WORK_I64
  public :: GADG_WORK_I
  public :: GADG_WORK_IOPT
  public :: GADG_WORK_SP
  public :: GADG_WORK_DP
  
  public :: gadg_reserve_work_i32
  public :: gadg_reserve_work_i64
  public :: gadg_reserve_work_i
  public :: gadg_reserve_work_iopt
  public :: gadg_reserve_work_sp
  public :: gadg_reserve_work_dp
  public :: gadg_dealloc_work_i32
  public :: gadg_dealloc_work_i64
  public :: gadg_dealloc_work_i
  public :: gadg_dealloc_work_iopt
  public :: gadg_dealloc_work_sp
  public :: gadg_dealloc_work_dp
  public :: gadg_dealloc_all_works
  public :: gadg_finalize
  public :: gadg_abort
  
contains

#define GADG_RESERVE_WORK_X gadg_reserve_work_i32
#define GADG_DEALLOC_WORK_X gadg_dealloc_work_i32
#define GADG_WORK_X         GADG_WORK_I32
#include "gadg_base.inc"

#define GADG_RESERVE_WORK_X gadg_reserve_work_i64
#define GADG_DEALLOC_WORK_X gadg_dealloc_work_i64
#define GADG_WORK_X         GADG_WORK_I64
#include "gadg_base.inc"

#define GADG_RESERVE_WORK_X gadg_reserve_work_i
#define GADG_DEALLOC_WORK_X gadg_dealloc_work_i
#define GADG_WORK_X         GADG_WORK_I
#include "gadg_base.inc"

#define GADG_RESERVE_WORK_X gadg_reserve_work_iopt
#define GADG_DEALLOC_WORK_X gadg_dealloc_work_iopt
#define GADG_WORK_X         GADG_WORK_IOPT
#include "gadg_base.inc"

#define GADG_RESERVE_WORK_X gadg_reserve_work_sp
#define GADG_DEALLOC_WORK_X gadg_dealloc_work_sp
#define GADG_WORK_X         GADG_WORK_SP
#include "gadg_base.inc"

#define GADG_RESERVE_WORK_X gadg_reserve_work_dp
#define GADG_DEALLOC_WORK_X gadg_dealloc_work_dp
#define GADG_WORK_X         GADG_WORK_DP
#include "gadg_base.inc"
  
  !###########################################################################
  subroutine gadg_dealloc_all_works()
    
    !=------------------------------------------------------------------------
    call gadg_dealloc_work_i32
    call gadg_dealloc_work_i64
    call gadg_dealloc_work_i
    call gadg_dealloc_work_iopt
    call gadg_dealloc_work_sp
    call gadg_dealloc_work_dp
  end subroutine gadg_dealloc_all_works
  
  !###########################################################################
  subroutine gadg_finalize()
    
    !=------------------------------------------------------------------------
    call gadg_dealloc_all_works
  end subroutine gadg_finalize
  
  !###########################################################################
  subroutine gadg_abort(sub_name, msg)
    !  outputs massage for fatal error and terminates the program.
    !  This is called from preprocessor macros.
    
#ifdef PP_CONF_USE_MPI
    use MPI
#else
    PP_CONF_USE_FOR_ABORT
#endif
    
    character(*), intent(in), optional :: sub_name ! subprogram name
    character(*), intent(in), optional :: msg ! error message
    
    !=------------------------------------------------------------------------
#ifdef PP_CONF_USE_MPI
    integer :: ierr
    logical :: flag
#endif
    
    if (present(sub_name)) print *, "In ", sub_name
    if (present(msg)) print *, msg
    print *, "Abnormal termination"
#ifdef PP_CONF_USE_MPI
    call MPI_Initialized (flag, ierr)  
    if (flag) then
      call MPI_Abort(MPI_COMM_WORLD, -1)
    end if
#endif
    PP_CONF_ABORT
    
  end subroutine gadg_abort
  
!end module m_gadg_base
end module gadg_base
