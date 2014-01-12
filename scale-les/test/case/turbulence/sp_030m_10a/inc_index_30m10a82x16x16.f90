  !-----------------------------------------------------------------------------
  !
  !++ SCALE-LES grid parameters (res. 30m, asp. 10)
  !
  !-----------------------------------------------------------------------------
  real(RP), private, parameter :: DZ  =  30 ! length in the main region [m]: z
  real(RP), private, parameter :: DX  = 300 ! length in the main region [m]: x
  real(RP), private, parameter :: DY  = 300 ! length in the main region [m]: y

  real(RP), private, parameter :: BUFFER_DZ = 1.E3_RP ! thickness of buffer region [m]: z
  real(RP), private, parameter :: BUFFER_DX = 0.0_RP ! thickness of buffer region [m]: x
  real(RP), private, parameter :: BUFFER_DY = 0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, parameter :: BUFFFACT  = 1.1_RP ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX  = 82 ! # of computational cells: z
  integer, private, parameter :: IMAX  = 16 ! # of computational cells: x
  integer, private, parameter :: JMAX  = 16 ! # of computational cells: y

  integer, private, parameter :: IBLOCK = IMAX !
  integer, private, parameter :: JBLOCK = JMAX !

  include 'inc_index_all.h'
