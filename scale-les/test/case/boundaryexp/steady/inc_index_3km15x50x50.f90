
  !-----------------------------------------------------------------------------
  !
  !++ SCALE-LES grid parameters (3km res., 300km isotropic, 7.5km model top)
  !
  !-----------------------------------------------------------------------------
  integer,  private, parameter :: KMAX =   15 ! # of computational cells: z
  integer,  private, parameter :: IMAX =   50 ! # of computational cells: x
  integer,  private, parameter :: JMAX =   50 ! # of computational cells: y

  integer,  private, parameter :: IBLOCK = 50 ! block size for cache blocking: x
  integer,  private, parameter :: JBLOCK = 50 ! block size for cache blocking: y

  real(RP), private, save      :: DZ        =  500.0_RP ! length in the main region [m]: z
  real(RP), private, save      :: DX        =  3000.0_RP ! length in the main region [m]: x
  real(RP), private, save      :: DY        =  3000.0_RP ! length in the main region [m]: y

  real(RP), private, save      :: BUFFER_DZ =  2000.0_RP ! thickness of buffer region [m]: z
  real(RP), private, save      :: BUFFER_DX =  9000.0_RP ! thickness of buffer region [m]: x
  real(RP), private, save      :: BUFFER_DY =  9000.0_RP ! thickness of buffer region [m]: y
  real(RP), private, save      :: BUFFFACT  =  1.0_RP    ! strech factor for dx/dy/dz of buffer region

  include 'inc_index_all.h'
