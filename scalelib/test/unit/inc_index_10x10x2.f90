
  !-----------------------------------------------------------------------------
  !
  !++ scale-les grid parameters (very small for unit tests)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: KHALO = 2 ! # of halo cells: z
  integer, private, parameter :: IHALO = 2 ! # of halo cells: x
  integer, private, parameter :: JHALO = 2 ! # of halo cells: y

  real(RP), private, save :: DX  = 200.0_RP ! length in the main region [m]: x
  real(RP), private, save :: DY  = 200.0_RP ! length in the main region [m]: y
  real(RP), private, save :: DZ  =  10.0_RP ! length in the main region [m]: z

  real(RP), private, save :: BUFFER_DZ =  50.E0_RP ! thickness of buffer region [m]: z
  real(RP), private, save :: BUFFER_DX =  0.0E0_RP ! thickness of buffer region [m]: x
  real(RP), private, save :: BUFFER_DY =  0.0E0_RP ! thickness of buffer region [m]: y
  real(RP), private, save :: BUFFFACT  =  1.2E0_RP ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX =  10 ! # of computational cells: z
  integer, private, parameter :: IMAX =  10 ! # of computational cells: x
  integer, private, parameter :: JMAX =   2 ! # of computational cells: y

  integer, private, parameter :: KA   =  14 ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA   =  14 ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA   =   6 ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS   =   3 ! start point of inner domain: z, local
  integer, private, parameter :: KE   =  12 ! end   point of inner domain: z, local
  integer, private, parameter :: IS   =   3 ! start point of inner domain: x, local
  integer, private, parameter :: IE   =  12 ! end   point of inner domain: x, local
  integer, private, parameter :: JS   =   3 ! start point of inner domain: y, local
  integer, private, parameter :: JE   =   4 ! end   point of inner domain: y, local

  integer, private, parameter :: IJA  =  20 ! # merged inner region(tentative usage)
  integer, private, parameter :: IJS  =   1 !
  integer, private, parameter :: IJE  =  20 !

  integer, private, parameter :: IBLOCK =  5 !
  integer, private, parameter :: JBLOCK =  1 !
