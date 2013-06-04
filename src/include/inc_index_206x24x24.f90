
  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters (50m res., 9km isotropic, 15km model top)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: KHALO = 2 ! # of halo cells: z
  integer, private, parameter :: IHALO = 2 ! # of halo cells: x
  integer, private, parameter :: JHALO = 2 ! # of halo cells: y

  real(RP), private, parameter :: DX  = 50.0_RP ! length in the main region [m]: x
  real(RP), private, parameter :: DY  = 50.0_RP ! length in the main region [m]: y
  real(RP), private, parameter :: DZ  = 50.0_RP ! length in the main region [m]: z

  real(RP), private, parameter :: BUFFER_DZ =  6.0E3_RP ! thickness of buffer region [m]: z
  real(RP), private, parameter :: BUFFER_DX =  0.0E0_RP ! thickness of buffer region [m]: x
  real(RP), private, parameter :: BUFFER_DY =  0.0E0_RP ! thickness of buffer region [m]: y
  real(RP), private, parameter :: BUFFFACT  =  1.1E0_RP ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX = 206 ! # of computational cells: z
  integer, private, parameter :: IMAX =  24 ! # of computational cells: x
  integer, private, parameter :: JMAX =  24 ! # of computational cells: y

  integer, private, parameter :: KA   = 210 ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA   =  28 ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA   =  28 ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS   =   3 ! start point of inner domain: z, local
  integer, private, parameter :: KE   = 208 ! end   point of inner domain: z, local
  integer, private, parameter :: IS   =   3 ! start point of inner domain: x, local
  integer, private, parameter :: IE   =  26 ! end   point of inner domain: x, local
  integer, private, parameter :: JS   =   3 ! start point of inner domain: y, local
  integer, private, parameter :: JE   =  26 ! end   point of inner domain: y, local

  integer, private, parameter :: IJA  =  576 ! # merged inner region(tentative usage)
  integer, private, parameter :: IJS  =    1 !
  integer, private, parameter :: IJE  =  576 !

  integer, private, parameter :: IBLOCK = 8 !
  integer, private, parameter :: JBLOCK = 8 !
