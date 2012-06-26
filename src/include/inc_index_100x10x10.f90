
  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters (200m res., 20km model top)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: KHALO = 2 ! # of halo cells: z
  integer, private, parameter :: IHALO = 2 ! # of halo cells: x
  integer, private, parameter :: JHALO = 2 ! # of halo cells: y

  real(8), private, parameter :: DX  = 200 ! length in the main region [m]: x
  real(8), private, parameter :: DY  = 200 ! length in the main region [m]: y
  real(8), private, parameter :: DZ  = 200 ! length in the main region [m]: z

  real(8), private, parameter :: BUFFER_DZ =  5.D3 ! thickness of buffer region [m]: z
  real(8), private, parameter :: BUFFER_DX =  0.0D0 ! thickness of buffer region [m]: x
  real(8), private, parameter :: BUFFER_DY =  0.0D0 ! thickness of buffer region [m]: y
  real(8), private, parameter :: BUFFFACT  =  1.0D0 ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX = 100 ! # of computational cells: z
  integer, private, parameter :: IMAX =  10 ! # of computational cells: x
  integer, private, parameter :: JMAX =  10 ! # of computational cells: y

  integer, private, parameter :: KA   = 104 ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA   =  14 ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA   =  14 ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS   =   3 ! start point of inner domain: z, local
  integer, private, parameter :: KE   = 102 ! end   point of inner domain: z, local
  integer, private, parameter :: IS   =   3 ! start point of inner domain: x, local
  integer, private, parameter :: IE   =  12 ! end   point of inner domain: x, local
  integer, private, parameter :: JS   =   3 ! start point of inner domain: y, local
  integer, private, parameter :: JE   =  12 ! end   point of inner domain: y, local

  integer, private, parameter :: IJA  = 100 ! # merged inner region(tentative usage)
  integer, private, parameter :: IJS  =    1 !
  integer, private, parameter :: IJE  = 100 !

  integer, private, parameter :: IBLOCK = 10 !
  integer, private, parameter :: JBLOCK = 10 !
