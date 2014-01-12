
  !-----------------------------------------------------------------------------
  !
  !++ scale-les grid index parameters (common)
  !
  !-----------------------------------------------------------------------------

  integer, private, parameter :: KHALO = 2               ! # of halo cells: z
  integer, private, parameter :: IHALO = 2               ! # of halo cells: x
  integer, private, parameter :: JHALO = 2               ! # of halo cells: y

  integer, private, parameter :: KA   = KMAX + KHALO * 2 ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA   = IMAX + IHALO * 2 ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA   = JMAX + JHALO * 2 ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS   = 1    + KHALO     ! start point of inner domain: z, local
  integer, private, parameter :: KE   = KMAX + KHALO     ! end   point of inner domain: z, local
  integer, private, parameter :: IS   = 1    + IHALO     ! start point of inner domain: x, local
  integer, private, parameter :: IE   = IMAX + IHALO     ! end   point of inner domain: x, local
  integer, private, parameter :: JS   = 1    + JHALO     ! start point of inner domain: y, local
  integer, private, parameter :: JE   = JMAX + JHALO     ! end   point of inner domain: y, local
