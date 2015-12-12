!-------------------------------------------------------------------------------
!> Pre-defined grid index
!!
!! @par Description
!!          Include file of grid index
!!          If you set environment "ENABLE_FIXEDINDEX=T", this file is used
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  !++ grid parameters (common)
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: KHALO = 2               ! # of halo cells: z
  integer, public, parameter :: IHALO = 2               ! # of halo cells: x
  integer, public, parameter :: JHALO = 2               ! # of halo cells: y

  integer, public, parameter :: KA   = KMAX + KHALO * 2 ! # of z whole cells (local, with HALO)
  integer, public, parameter :: IA   = IMAX + IHALO * 2 ! # of x whole cells (local, with HALO)
  integer, public, parameter :: JA   = JMAX + JHALO * 2 ! # of y whole cells (local, with HALO)

  integer, public, parameter :: KS   = 1    + KHALO     ! start point of inner domain: z, local
  integer, public, parameter :: KE   = KMAX + KHALO     ! end   point of inner domain: z, local
  integer, public, parameter :: IS   = 1    + IHALO     ! start point of inner domain: x, local
  integer, public, parameter :: IE   = IMAX + IHALO     ! end   point of inner domain: x, local
  integer, public, parameter :: JS   = 1    + JHALO     ! start point of inner domain: y, local
  integer, public, parameter :: JE   = JMAX + JHALO     ! end   point of inner domain: y, local

  integer, public, parameter :: KIJMAX = KMAX * IMAX * JMAX ! # of computational cells: z*x*y
