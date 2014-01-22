!fpx3_header(0.20)
!
!dependencies
!   dir: ~/work/scale3_140107/Gadgetry 
!   sources: gadg_constants.F90_fpx
!   includes: 
!   uses: 
!   provides: m_gadg_constants
!end dependencies
!
!end fpx3_header
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

!module m_gadg_constants
module gadg_constants
  ! Physical and Mathematical Constants
  implicit none
  private

  integer, parameter :: dp = selected_real_kind(12)

  ! Mathematical constants
  !
  !   MATH_E         -- e,
  !   MATH_LOG2E     -- log2 e,
  !   MATH_LOG10E    -- log10 e,
  !   MATH_PI        -- pi,
  !   MATH_LN2       -- ln 2,
  !   MATH_LN10      -- ln 10,
  !   MATH_SQRT_PI   -- sqrt(pi),
  !   MATH_2PI       -- 2 * pi,
  !   MATH_PI_ON_2   -- pi / 2,
  !   MATH_PI_ON_4   -- pi / 4,
  !   MATH_1_ON_PI   -- 1 / pi,
  !   MATH_1_ON_2PI  -- 1 / (2pi),
  !   MATH_SQRT2     -- sqrt(2),
  !   MATH_SQRT3     -- sqrt(3),
  !   MATH_SQRT5     -- sqrt(5),
  !   MATH_SQRT6     -- sqrt(6),
  !   MATH_SQRT7     -- sqrt(7),
  !   MATH_SQRT8     -- sqrt(8),
  !   MATH_SQRT10    -- sqrt(10).
  real(dp), parameter, public :: MATH_E        = 2.718281828459045235360287471352662497757d0
  real(dp), parameter, public :: MATH_LOG2E    = 1.442695040888963407359924681001892137427d0
  real(dp), parameter, public :: MATH_LOG10E   = 0.434294481903251827651128918916605082294d0
  real(dp), parameter, public :: MATH_PI       = 3.141592653589793238462643383279502884197d0
  real(dp), parameter, public :: MATH_LN2      = 0.693147180559945309417232121458176568076d0
  real(dp), parameter, public :: MATH_LN10     = 2.302585092994045684017991454684364207601d0
  real(dp), parameter, public :: MATH_SQRT_PI  = 1.772453850905516027298167483341145182798d0
  real(dp), parameter, public :: MATH_2PI      = 6.283185307179586476925286766559005768394d0
  real(dp), parameter, public :: MATH_PI_ON_2  = 1.570796326794896619231321691639751442099d0
  real(dp), parameter, public :: MATH_PI_ON_4  = 0.785398163397448309615660845819875721049d0
  real(dp), parameter, public :: MATH_1_ON_PI  = 0.318309886183790671537767526745028724069d0
  real(dp), parameter, public :: MATH_1_ON_2PI = 0.159154943091895335768883763372514362034d0
  real(dp), parameter, public :: MATH_SQRT2    = 1.414213562373095048801688724209698078570d0
  real(dp), parameter, public :: MATH_SQRT3    = 1.732050807568877293527446341505872366943d0
  real(dp), parameter, public :: MATH_SQRT5    = 2.236067977499789696409173668731276235441d0
  real(dp), parameter, public :: MATH_SQRT6    = 2.449489742783178098197284074705891391966d0
  real(dp), parameter, public :: MATH_SQRT7    = 2.645751311064590590501615753639260425710d0
  real(dp), parameter, public :: MATH_SQRT8    = 2.828427124746190097603377448419396157139d0
  real(dp), parameter, public :: MATH_SQRT10   = 3.162277660168379331998893544432718533720d0

  ! Physical constants in MKSA
  !  This contains
  !   PHYS_PLANCKS_H        -- Planck's constant h            /  Js,
  !   PHYS_SPEED_OF_LIGHT   -- The speed of light in vacuum   /  m s^-1,
  !   PHYS_AVOGADRO         -- Avogadro's number              /  mol^-1,
  !   PHYS_BOHR_RADIUS      -- The Bohr radius                /  m,
  !   PHYS_MASS_OF_ELECTRON -- The mass of the electron       /  kg,
  !   PHYS_PLANCKS_HBAR     -- Planck's constant hbar         /  Js.
  !   PHYS_BOLTZMANN        -- Boltzmann's constant           / J K^-1
  !   PHYS_GAS_CONSTANTS    -- Ideal gas constants            / J K^-1 mol^-1
  real(dp), parameter, public :: PHYS_SPEED_OF_LIGHT   = 299792458.0d0                    ! m s^-1
  real(dp), parameter, public :: PHYS_AVOGADRO         = 6.02214199d23                    ! mol^-1
  real(dp), parameter, public :: PHYS_BOHR_RADIUS      = 0.5291772083d-10                 ! m 
  real(dp), parameter, public :: PHYS_MASS_OF_ELECTRON = 9.10938188d-31                   ! kg
  real(dp), parameter, public :: PHYS_PLANCKS_H        = 6.62606876d-34                   ! J s
  real(dp), parameter, public :: PHYS_PLANCKS_HBAR     = PHYS_PLANCKS_H / MATH_2PI        ! J s
  real(dp), parameter, public :: PHYS_BOLTZMANN        = 1.3806503d-23                    ! J K^-1
  real(dp), parameter, public :: PHYS_GAS_CONSTANT     = PHYS_BOLTZMANN * PHYS_AVOGADRO   ! J K^-1 mol^-1

  !=== Conversion factors ===

  ! Dimensionless prefixes
  real(dp), parameter, public :: CONV_PF_PETA  = 1.0d15  ! P
  real(dp), parameter, public :: CONV_PF_TERA  = 1.0d12  ! T
  real(dp), parameter, public :: CONV_PF_GIGA  = 1.0d9   ! G
  real(dp), parameter, public :: CONV_PF_MEGA  = 1.0d6   ! M
  real(dp), parameter, public :: CONV_PF_KILO  = 1.0d3   ! K
  real(dp), parameter, public :: CONV_PF_HECTO = 1.0d2   ! h
  real(dp), parameter, public :: CONV_PF_DECA  = 1.0d1   ! da
  real(dp), parameter, public :: CONV_PF_UNITY = 1.0d0
  real(dp), parameter, public :: CONV_PF_DECI  = 1.0d-1  ! d
  real(dp), parameter, public :: CONV_PF_CENTI = 1.0d-2  ! c
  real(dp), parameter, public :: CONV_PF_MILLI = 1.0d-3  ! m
  real(dp), parameter, public :: CONV_PF_MICRO = 1.0d-6  ! u
  real(dp), parameter, public :: CONV_PF_NANO  = 1.0d-9  ! n
  real(dp), parameter, public :: CONV_PF_PICO  = 1.0d-12 ! p
  real(dp), parameter, public :: CONV_PF_FEMTO = 1.0d-15 ! f

  ! Conversion factors of energy units.
  ! CONV_(A)_TO_JOULE,
  ! (A) is ones of the following items:
  !   JOULE      -- joule (J),
  !   EV         -- electron volt (eV),
  !   AUENERGY   -- au, or hartree,
  !   WAVENUM    -- wave number (cm^-1),
  !   CALORIE    -- calorie (cal),
  !   KCALPERMOL -- (kcal mol^-1).
  real(dp), parameter, public :: CONV_EV_TO_JOULE        = 1.602176462d-19
  real(dp), parameter, public :: CONV_AUENERGY_TO_JOULE  = 4.35974381d-18
  real(dp), parameter, public :: CONV_WAVENUM_TO_JOULE   = PHYS_PLANCKS_H * PHYS_SPEED_OF_LIGHT * 100.d0
  real(dp), parameter, public :: CONV_CALORIE_TO_JOULE   = 4.1840d0
  real(dp), parameter, public :: CONV_KCALPERMOL_TO_JOULE= 1000.d0 * CONV_CALORIE_TO_JOULE / PHYS_AVOGADRO

  real(dp), parameter, public :: CONV_JOULE_TO_EV = 1.0d0 / CONV_EV_TO_JOULE
  real(dp), parameter, public :: CONV_JOULE_TO_AUENERGY = 1.0d0 / CONV_AUENERGY_TO_JOULE
  real(dp), parameter, public :: CONV_JOULE_TO_WAVENUM = 1.0d0 / CONV_WAVENUM_TO_JOULE
  real(dp), parameter, public :: CONV_JOULE_TO_CALORIE = 1.0d0 / CONV_CALORIE_TO_JOULE
  real(dp), parameter, public :: CONV_JOULE_TO_KCALPERMOL = 1.0d0 / CONV_KCALPERMOL_TO_JOULE
  real(dp), parameter, public :: CONV_EV_TO_AUENERGY = CONV_EV_TO_JOULE / CONV_JOULE_TO_AUENERGY
  real(dp), parameter, public :: CONV_EV_TO_WAVENUM = CONV_EV_TO_JOULE / CONV_JOULE_TO_WAVENUM
  real(dp), parameter, public :: CONV_EV_TO_CALORIE = CONV_EV_TO_JOULE / CONV_JOULE_TO_CALORIE
  real(dp), parameter, public :: CONV_EV_TO_KCALPERMOL = CONV_EV_TO_JOULE / CONV_JOULE_TO_KCALPERMOL
  real(dp), parameter, public :: CONV_AUENERGY_TO_EV = CONV_AUENERGY_TO_JOULE / CONV_JOULE_TO_EV
  real(dp), parameter, public :: CONV_AUENERGY_TO_WAVENUM = CONV_AUENERGY_TO_JOULE / CONV_JOULE_TO_WAVENUM
  real(dp), parameter, public :: CONV_AUENERGY_TO_CALORIE = CONV_AUENERGY_TO_JOULE / CONV_JOULE_TO_CALORIE
  real(dp), parameter, public :: CONV_AUENERGY_TO_KCALPERMOL = CONV_AUENERGY_TO_JOULE / CONV_JOULE_TO_KCALPERMOL
  real(dp), parameter, public :: CONV_WAVENUM_TO_EV = CONV_WAVENUM_TO_JOULE / CONV_JOULE_TO_EV
  real(dp), parameter, public :: CONV_WAVENUM_TO_AUENERGY = CONV_WAVENUM_TO_JOULE / CONV_JOULE_TO_AUENERGY
  real(dp), parameter, public :: CONV_WAVENUM_TO_CALORIE = CONV_WAVENUM_TO_JOULE / CONV_JOULE_TO_CALORIE
  real(dp), parameter, public :: CONV_WAVENUM_TO_KCALPERMOL = CONV_WAVENUM_TO_JOULE / CONV_JOULE_TO_KCALPERMOL
  real(dp), parameter, public :: CONV_CALORIE_TO_EV = CONV_CALORIE_TO_JOULE / CONV_JOULE_TO_EV
  real(dp), parameter, public :: CONV_CALORIE_TO_AUENERGY = CONV_CALORIE_TO_JOULE / CONV_JOULE_TO_AUENERGY
  real(dp), parameter, public :: CONV_CALORIE_TO_WAVENUM = CONV_CALORIE_TO_JOULE / CONV_JOULE_TO_WAVENUM
  real(dp), parameter, public :: CONV_CALORIE_TO_KCALPERMOL = CONV_CALORIE_TO_JOULE / CONV_JOULE_TO_KCALPERMOL
  real(dp), parameter, public :: CONV_KCALPERMOL_TO_EV = CONV_KCALPERMOL_TO_JOULE / CONV_JOULE_TO_EV
  real(dp), parameter, public :: CONV_KCALPERMOL_TO_AUENERGY = CONV_KCALPERMOL_TO_JOULE / CONV_JOULE_TO_AUENERGY
  real(dp), parameter, public :: CONV_KCALPERMOL_TO_WAVENUM = CONV_KCALPERMOL_TO_JOULE / CONV_JOULE_TO_WAVENUM
  real(dp), parameter, public :: CONV_KCALPERMOL_TO_CALORIE = CONV_KCALPERMOL_TO_JOULE / CONV_JOULE_TO_CALORIE

  ! Conversion factors of length units, quadruple precision.
  ! CONV_(A)_TO_METER,
  ! (A) is ones of the following items:
  !   METER     -- meter (m),
  !   ANGSTROM  -- angstrom,
  !   AULENGTH  -- au, or bohr,
  real(dp), parameter, public :: CONV_ANGSTROM_TO_METER  = 1.0d-10
  real(dp), parameter, public :: CONV_AULENGTH_TO_METER  = PHYS_BOHR_RADIUS

  real(dp), parameter, public :: CONV_METER_TO_ANGSTROM = 1.0d0 / CONV_ANGSTROM_TO_METER
  real(dp), parameter, public :: CONV_METER_TO_AULENGTH = 1.0d0 / CONV_AULENGTH_TO_METER
  real(dp), parameter, public :: CONV_ANGSTROM_TO_AULENGTH = CONV_ANGSTROM_TO_METER / CONV_METER_TO_AULENGTH
  real(dp), parameter, public :: CONV_AULENGTH_TO_ANGSTROM = CONV_AULENGTH_TO_METER / CONV_METER_TO_ANGSTROM

  ! Conversion factors of time units, quadruple precision.
  ! CONV_(A)_TO_SECOND,
  ! (A) is ones of the following items:
  !   SECOND -- second (s).
  !   AUTIME -- atomic unit of time,
  real(dp), parameter, public :: CONV_AUTIME_TO_SECOND = 2.418884326500d-17;

  real(dp), parameter, public :: CONV_SECOND_TO_AUTIME = 1.0d0 / CONV_AUTIME_TO_SECOND

!end module m_gadg_constants
end module gadg_constants
