  integer, private, parameter :: KHALO = 2 ! # of halo cells: z
  integer, private, parameter :: IHALO = 2 ! # of halo cells: x
  integer, private, parameter :: JHALO = 2 ! # of halo cells: y

  integer, private, parameter :: KA   = KMAX+KHALO*2 ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA   = IMAX+IHALO*2 ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA   = JMAX+JHALO*2 ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS   = KHALO+1    ! start point of inner domain: z, local
  integer, private, parameter :: KE   = KMAX+KHALO ! end   point of inner domain: z, local
  integer, private, parameter :: IS   = IHALO+1    ! start point of inner domain: x, local
  integer, private, parameter :: IE   = IMAX+IHALO ! end   point of inner domain: x, local
  integer, private, parameter :: JS   = JHALO+1    ! start point of inner domain: y, local
  integer, private, parameter :: JE   = JMAX+JHALO ! end   point of inner domain: y, local

  integer, private, parameter :: IJA  = IMAX*JMAX ! # merged inner region(tentative usage)
  integer, private, parameter :: IJS  = 1   !
  integer, private, parameter :: IJE  = IJA !
