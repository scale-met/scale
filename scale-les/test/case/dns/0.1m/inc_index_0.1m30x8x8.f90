
  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters (0.1m res., 3m isotropic)
  !
  !-----------------------------------------------------------------------------
! 30grids x 4node = 120 grids
! 3m depth / 0.1 m = 30 (8grids x 4node)
! 3m depth / 0.01 m = 300
  integer,  private, parameter :: KMAX =   30 ! # of computational cells: z
  integer,  private, parameter :: IMAX =    8 ! # of computational cells: x
  integer,  private, parameter :: JMAX =    8 ! # of computational cells: y

  integer,  private, parameter :: IBLOCK =  8 ! block size for cache blocking: x
  integer,  private, parameter :: JBLOCK =  8 ! block size for cache blocking: y

  real(RP), private, save      :: DZ        =  0.1_RP ! length in the main region [m]: z
  real(RP), private, save      :: DX        =  0.1_RP ! length in the main region [m]: x
  real(RP), private, save      :: DY        =  0.1_RP ! length in the main region [m]: y

! assume 5 layers for sponge region (0.025*5=0.125)
! real(RP), private, save      :: BUFFER_DZ =  0.125_RP ! thickness of buffer region [m]: z
  real(RP), private, save      :: BUFFER_DZ =    0.5_RP ! thickness of buffer region [m]: z
  real(RP), private, save      :: BUFFER_DX =    0.0_RP ! thickness of buffer region [m]: x
  real(RP), private, save      :: BUFFER_DY =    0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, save      :: BUFFFACT  =    1.1_RP ! strech factor for dx/dy/dz of buffer region

  include 'inc_index_all.h'
