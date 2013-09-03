
  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: KHALO =   2 ! # of halo cells: z
  integer, private, parameter :: IHALO =   2 ! # of halo cells: x
  integer, private, parameter :: JHALO =   2 ! # of halo cells: y

  real(RP), private, save :: DZ  = 5.0_RP ! length in the main region [m]: x,y,z
  real(RP), private, save :: DX  = 50.0_RP ! length in the main region [m]: x,y,z
  real(RP), private, save :: DY  = 50.0_RP ! length in the main region [m]: x,y,z

  real(RP), private, save :: BUFFER_DZ =6.E2_RP! thickness of buffer region [m]: z
  real(RP), private, save :: BUFFER_DX = 0.0_RP ! thickness of buffer region [m]: x
  real(RP), private, save :: BUFFER_DY = 0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, save :: BUFFFACT  = 1.1_RP ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX  = 200 ! # of computational cells: z
  integer, private, parameter :: IMAX  =  11 ! # of computational cells: x
  integer, private, parameter :: JMAX  =  11 ! # of computational cells: y

  integer, private, parameter :: KA    = 204 ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA    =  15 ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA    =  15 ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS    =   3 ! start point of inner domain: z, local
  integer, private, parameter :: KE    = 202 ! end   point of inner domain: z, local
  integer, private, parameter :: IS    =   3 ! start point of inner domain: x, local
  integer, private, parameter :: IE    =  13 ! end   point of inner domain: x, local
  integer, private, parameter :: JS    =   3 ! start point of inner domain: y, local
  integer, private, parameter :: JE    =  13 ! end   point of inner domain: y, local

  integer, private, parameter :: IJA   = 121 ! # merged inner region(tentative usage)
  integer, private, parameter :: IJS   =   1 !
  integer, private, parameter :: IJE   = 121 !

  integer, private, parameter :: IBLOCK = 11 !
  integer, private, parameter :: JBLOCK = 11 !
