
  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters (5m res., 2km isotropic, 3.5km model top)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: KHALO = 2 ! # of halo cells: z
  integer, private, parameter :: IHALO = 2 ! # of halo cells: x
  integer, private, parameter :: JHALO = 2 ! # of halo cells: y

  real(8), private, parameter :: DXYZ  = 5 ! length in the main region [m]: x,y,z

  real(8), private, parameter :: BUFFER_DZ =  1.5D3 ! thickness of buffer region [m]: z
  real(8), private, parameter :: BUFFER_DX =  0.0D0 ! thickness of buffer region [m]: x
  real(8), private, parameter :: BUFFER_DY =  0.0D0 ! thickness of buffer region [m]: y
  real(8), private, parameter :: BUFFFACT  =  1.1D0 ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX = 436 ! # of computational cells: z
  integer, private, parameter :: IMAX =  14 ! # of computational cells: x
  integer, private, parameter :: JMAX =  14 ! # of computational cells: y

  integer, private, parameter :: KA   = 440 ! # of z whole cells (local, with HALO)
  integer, private, parameter :: IA   =  18 ! # of x whole cells (local, with HALO)
  integer, private, parameter :: JA   =  18 ! # of y whole cells (local, with HALO)

  integer, private, parameter :: KS   =   3 ! start point of inner domain: z, local
  integer, private, parameter :: KE   = 438 ! end   point of inner domain: z, local
  integer, private, parameter :: IS   =   3 ! start point of inner domain: x, local
  integer, private, parameter :: IE   =  16 ! end   point of inner domain: x, local
  integer, private, parameter :: JS   =   3 ! start point of inner domain: y, local
  integer, private, parameter :: JE   =  16 ! end   point of inner domain: y, local

  integer, private, parameter :: IJA  =  196 ! # merged inner region(tentative usage)
  integer, private, parameter :: IJS  =    1 !
  integer, private, parameter :: IJE  =  196 !

  integer, private, parameter :: IBLOCK = 7 !
  integer, private, parameter :: JBLOCK = 7 !