  !-----------------------------------------------------------------------------
  !
  !++ SCALE-LES grid parameters (3.125m res)
  !
  !-----------------------------------------------------------------------------
  real(RP), private, parameter :: DX  = 3.125_RP ! length in the main region [m]: x
  real(RP), private, parameter :: DY  = 3.125_RP ! length in the main region [m]: y
  real(RP), private, parameter :: DZ  = 3.125_RP ! length in the main region [m]: z

  real(RP), private, parameter :: BUFFER_DZ =  0.0_RP ! thickness of buffer region [m]: z
  real(RP), private, parameter :: BUFFER_DX =  0.0_RP ! thickness of buffer region [m]: x
  real(RP), private, parameter :: BUFFER_DY =  0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, parameter :: BUFFFACT  =  1.0_RP ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX = 2048 ! # of computational cells: z
  integer, private, parameter :: IMAX =    8 ! # of computational cells: x
  integer, private, parameter :: JMAX =    2 ! # of computational cells: y

  integer, private, parameter :: IBLOCK = IMAX !
  integer, private, parameter :: JBLOCK = JMAX !

  include 'inc_index_all.h'
