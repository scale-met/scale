
  !-----------------------------------------------------------------------------
  !
  !++ SCALE-LES grid parameters (100m for vertical, 1000m for horizontal res)
  !
  !-----------------------------------------------------------------------------
  integer,  private, parameter :: KMAX =   80 ! # of computational cells: z
  integer,  private, parameter :: IMAX =   15 ! # of computational cells: x
  integer,  private, parameter :: JMAX =   15 ! # of computational cells: y

  integer,  private, parameter :: IBLOCK = 15 ! block size for cache blocking: x
  integer,  private, parameter :: JBLOCK = 15 ! block size for cache blocking: y

  real(RP), private, save      :: DZ        =  100.0_RP ! length in the main region [m]: z
  real(RP), private, save      :: DX        = 1000.0_RP ! length in the main region [m]: x
  real(RP), private, save      :: DY        = 1000.0_RP ! length in the main region [m]: y

  real(RP), private, save      :: BUFFER_DZ = 2000.0_RP ! thickness of buffer region [m]: z
  real(RP), private, save      :: BUFFER_DX = 2000.0_RP ! thickness of buffer region [m]: x
  real(RP), private, save      :: BUFFER_DY =    0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, save      :: BUFFFACT  =   1.05_RP ! strech factor for dx/dy/dz of buffer region

  include 'inc_index_all.h'
