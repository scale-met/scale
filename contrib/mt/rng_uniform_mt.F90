!=COPYRIGHT
! Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! 3. The names of its contributors may not be used to endorse or promote
! products derived from this software without specific prior written
! permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!=DESCRIPTION
! Module for the MT19937, a fast and quality pseudorandom number
! generator. This is a variant of the twisted generalized feedback
! shift-register algorithm, and is known as the "Mersenne Twister".
! It has a Mersenne prime period of 2^19937 - 1 (about 10^6000)
! and is equi-distributed in 623 dimensions. It has passed the DIEHARD
! statistical tests. It uses 624 words of state per generator and is
! comparable in speed to the other generators.
!
!=NOTES
! The original code included the comment: "When you use this, send an
! email to: matumoto@math.keio.ac.jp with an appropriate reference to
! your work".
!
!=SEE ALSO
! Makoto Matsumoto has a web page with more information about the
! generator, http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/mt.html
!
! The paper below has details of the algorithm.
!
! From: Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
! 623-dimensionally equidistributerd uniform pseudorandom number
! generator". ACM Transactions on Modeling and Computer Simulation,
! Vol. 8, No. 1 (Jan. 1998), Pages 3-30
!
! You can obtain the paper directly from Makoto Matsumoto's web page.
!
!=HISTORY
! Jan-13-1999 Hiroshi Takano::  Fortran translation.
! Feb-??-2004 Akio Kawano::     Translation into Fortran90 compatible with F.
! May-16-2004 Akio Kawano::     Changed to use the second revision of the
!                               seeding procedure published in 2002.
! Mar-04-2006 Akio Kawano::     Modified for Fortran90/ES
!

! Uncommenting the following might provide a faster object code for 64-bit CPU's.
! #define CONF_USE_64BIT_ARITHMETIC

module rng_uniform_mt
  implicit none
  private

  integer, parameter :: int32 = selected_int_kind(9)

#ifdef CONF_USE_64BIT_ARITHMETIC
  integer, parameter :: intopt = selected_int_kind(15) ! optimal integer kind
#else
  integer, parameter :: intopt = int32 ! optimal integer kind
#endif

  integer, parameter :: dp = selected_real_kind(12)

  integer, parameter :: n = 624, m = 397  ! Period parameters

  type :: c_rng_uniform_mt
    integer(intopt) :: mt(0:n-1)  ! The array for the state vector
    integer :: mti = n+1          ! mti == n+1 means mt[n] is not initialized
  end type

  ! Although non-decimal constants are permitted only in DATA statements
  ! in the Fortran95, most compilers and Fortran2003 extend the use of
  ! non-decimal constants.

  integer(intopt), parameter :: mata   = z"9908b0df"  ! constant vector a
  integer(intopt), parameter :: umask  = z"80000000"  ! most significant w-r bits
  integer(intopt), parameter :: lmask  = z"7fffffff"  ! least significant r bits
  integer(intopt), parameter :: tmaskb = z"9d2c5680"  ! tempering parameter b
  integer(intopt), parameter :: tmaskc = z"efc60000"  ! tempering parameter c

  integer, parameter :: default_seed = 4357   ! default random seed

  real(dp), parameter :: factor_int2real = 1.0D0 / (2.0D0 ** 32 + 1.0D0)
#ifdef CONF_USE_64BIT_ARITHMETIC
  real(dp), parameter :: shift_const = 0.5D0 * factor_int2real
#else
  real(dp), parameter :: shift_const = 0.5D0 + 0.5D0 * factor_int2real
#endif

  interface rng_init
    module procedure rng_init_
  end interface

  interface rng_load_state
    module procedure rng_load_state_
  end interface

  interface rng_generate
    module procedure rng_generate_
  end interface

  interface rng_generate_array
    module procedure rng_generate_array_
  end interface

  interface rng_save_state
    module procedure rng_save_state_
  end interface

  interface rng_final
    module procedure rng_final_
  end interface

  public :: c_rng_uniform_mt

  public :: rng_init
  public :: rng_load_state
  public :: rng_save_state
  public :: rng_generate
  public :: rng_generate_array
  public :: rng_final

  ! Obsolete interfaces
  interface gen_rand
    module procedure rng_generate_
  end interface

  interface gen_rand_array
    module procedure rng_generate_array_
  end interface

  public :: gen_rand, gen_rand_array

  character(len=61) :: ERR_FMT = &
  &   '(2X, "Error", I4, " in ", A, "in m_rng_uniform_mt -- ", A)'

contains

  !###########################################################################
  subroutine rng_init_(self, seed)
    ! Sets initial seeds to mt(:). Before call rng_generate(), this must be
    ! called once.

    type(c_rng_uniform_mt), intent(inout) :: self
    integer, intent(in) :: seed

    !=------------------------------------------------------------------------
    integer :: i
    integer(int32) ::  kmt
    integer(intopt), parameter :: mask = z"ffffffff"

    kmt = int(seed, int32)
    self%mt(0) = kmt
    do i = 1, n - 1
      !!  See Knuth's "Art of Computer Programming" Vol. 2, 3rd Ed.
      !!  p.106 for multiplier.
      kmt = 1812433253_int32 * ieor(kmt, ishft(kmt, -30)) + i
      self%mt(i) = kmt
    end do
    self%mt(:) = iand(self%mt(:), mask)

    self%mti = n
    ! Warm up
    do i = 1, 100
      call reload(self)
    end do

    return
  end subroutine rng_init_

  !###########################################################################
  function rng_generate_(self) result(rand)
    ! Generates one pseudorandom real number
    ! uniformly distributed on (0,1)-interval.

    type(c_rng_uniform_mt), intent(inout) :: self
    real(dp) :: rand

    !=------------------------------------------------------------------------
    integer(intopt) :: y

    if (self%mti >= n) call reload(self)
    y = self%mt(self%mti)
    self%mti = self%mti + 1
    y = ieor( y, ishft(y, -11) )
    y = ieor( y, iand(ishft(y, 7), tmaskb))
    y = ieor( y, iand(ishft(y, 15), tmaskc))
    y = ieor( y, ishft(y, -18) )
    rand = y * factor_int2real + shift_const

    return
  end function rng_generate_

  !###########################################################################
  subroutine rng_generate_array_(self, arr)
    ! Generates pseudorandom real number array
    ! uniformly distributed on (0,1)-interval.

    type(c_rng_uniform_mt), intent(inout) :: self
    real(dp), intent(inout) :: arr(1:)

    !=------------------------------------------------------------------------
    integer(intopt) :: y
    integer :: iend, i, j, k, nrest_arr, nrest_mt, nloop

    i = 1
    iend = ubound(arr, 1)
    do while (i <= iend)
      if (self%mti > n - 1) call reload(self)
      nrest_arr = iend - i
      nrest_mt = (n - 1) - self%mti
      nloop = min(nrest_arr, nrest_mt)
      k = self%mti
      do j = 0, nloop
        y = self%mt(k + j)
        y = ieor(y, ishft(y, -11))
        y = ieor(y, iand(ishft(y, 7), tmaskb))
        y = ieor(y, iand(ishft(y, 15), tmaskc))
        y = ieor(y, ishft(y, -18))
        arr(i + j) =  y * factor_int2real + shift_const
      end do
      self%mti = self%mti + (nloop + 1)
      i = i + (nloop + 1)
    end do

    return
  end subroutine rng_generate_array_

  !###########################################################################
  subroutine rng_final_(self)
    ! rng_finalizes the generator

    type(c_rng_uniform_mt), intent(inout) :: self

    !=------------------------------------------------------------------------
    self%mti = n+1

    return
  end subroutine rng_final_

  !###########################################################################
  subroutine rng_save_state_(self, filename, unit)
    ! Save the state of the object into "filename"
    ! If unit is specified, filename is ignored and
    ! the object is written in the unit.

    type(c_rng_uniform_mt), intent(in) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: unit

    !=-----------------------------------------------------------------------
    integer :: io_unit

    io_unit = 2993
    if (present(unit)) then
      io_unit = unit
    else
      open (io_unit, file = filename, action = "write", &
      &     access = "sequential", status = "replace", form = "unformatted", &
      &     err = 999)
    end if
    write (io_unit, err = 999) "rng_uniform_mt"
    write (io_unit, err = 999) self%mt(:), self%mti
    if (.not. present(unit)) close (io_unit, err = 999)
    return

    999 continue
    print ERR_FMT, 1, "rng_save_state", "Write error: " // filename
    stop
  end subroutine rng_save_state_

  !###########################################################################
  subroutine rng_load_state_(self, filename, unit)
    ! Initialize the state of the object by "filename"
    ! If unit is specified, filename is ignored and
    ! the object is initialized by the unit.

    type(c_rng_uniform_mt), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: unit

    !=-----------------------------------------------------------------------
    integer :: io_unit
    character(len=14) :: s

    io_unit = 2993
    if (present(unit)) then
      io_unit = unit
    else
      open (io_unit, file = filename, action = "read", &
      &     access = "sequential", status = "old", form = "unformatted", &
      &     err = 999)
    end if
    read(io_unit, err = 999) s
    if (s /= "rng_uniform_mt") then
      print ERR_FMT, 2, "rng_load_state", "Wrong format: " // filename
      stop
    end if
    read (io_unit, err = 999) self%mt(:), self%mti
    if (.not. present(unit)) close (io_unit, err = 999)
    return

    999 continue
    print ERR_FMT, 3, "rng_load_state", "Read error: " // filename
    stop
  end subroutine rng_load_state_

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  subroutine reload(self)

    type(c_rng_uniform_mt), intent(inout) :: self

    !=------------------------------------------------------------------------
    integer(intopt) :: y, z, x0, x1, x2
    integer :: i

    ! Generate N words at one time
    if(self%mti == n + 1) then
        ! A default initial seed is used if mt_srand() has not been called
        call rng_init(self, default_seed)
    end if

    ! the original code:
    !   for (kk=0;kk<N-M;kk++) {
    !       y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    !       mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    !   }
    !   for (;kk<N-1;kk++) {
    !       y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    !       mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    !   }
    !   y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    !   mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    do i = 0, n - m - 1
      x0 = self%mt(i)
      x1 = self%mt(i + 1)
      x2 = self%mt(i + m)
      z = iand(x1, 1_intopt) * mata
      y = ishft(iand(x0, umask) + iand(x1, lmask), -1)
      self%mt(i) = ieor(ieor(x2, y), z)
    end do

    do i = n - m, 2 * (n - m) - 1
      x0 = self%mt(i)
      x1 = self%mt(i + 1)
      x2 = self%mt(i + m - n)
      z = iand(x1, 1_intopt) * mata
      y = ishft(iand(x0, umask) + iand(x1, lmask), -1)
      self%mt(i) = ieor(ieor(x2, y), z)
    end do

    do i = 2 * (n - m), n - 2
      x0 = self%mt(i)
      x1 = self%mt(i + 1)
      x2 = self%mt(i + m - n)
      z = iand(x1, 1_intopt) * mata
      y = ishft(iand(x0, umask) + iand(x1, lmask), -1)
      self%mt(i) = ieor(ieor(x2, y), z)
    end do

    x0 = self%mt(n - 1)
    x1 = self%mt(0)
    x2 = self%mt(m - 1)
    z = iand(x1, 1_intopt) * mata
    y = ishft(iand(x0, umask) + iand(x1, lmask), -1)
    self%mt(n - 1) = ieor(ieor(x2, y), z)
    self%mti = 0

    return
  end subroutine reload

end module rng_uniform_mt
