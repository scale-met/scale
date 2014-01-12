
  !-----------------------------------------------------------------------------
  !
  !++ scale-les grid parameters (500m res., 17km isotropic, 22km model top)
  !
  !-----------------------------------------------------------------------------
  integer,  private, parameter :: KMAX =   40 ! # of computational cells: z
  integer,  private, parameter :: IMAX =   40 ! # of computational cells: x
  integer,  private, parameter :: JMAX =   40 ! # of computational cells: y

  integer,  private, parameter :: IBLOCK = 40 ! block size for cache blocking: x
  integer,  private, parameter :: JBLOCK = 40 ! block size for cache blocking: y

  real(RP), private, save      :: DZ        =  500.0_RP ! length in the main region [m]: z
  real(RP), private, save      :: DX        =  500.0_RP ! length in the main region [m]: x
  real(RP), private, save      :: DY        =  500.0_RP ! length in the main region [m]: y

  real(RP), private, save      :: BUFFER_DZ = 5000.0_RP ! thickness of buffer region [m]: z
  real(RP), private, save      :: BUFFER_DX =    0.0_RP ! thickness of buffer region [m]: x
  real(RP), private, save      :: BUFFER_DY =    0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, save      :: BUFFFACT  =    1.1_RP ! strech factor for dx/dy/dz of buffer region

  include 'inc_index_all.h'
