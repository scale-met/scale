  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters (res. 200m)
  !
  !-----------------------------------------------------------------------------
  real(RP), private, parameter :: DZ  = 200 ! length in the main region [m]: z
  real(RP), private, parameter :: DX  = 200 ! length in the main region [m]: x
  real(RP), private, parameter :: DY  = 200 ! length in the main region [m]: y

  real(RP), private, parameter :: BUFFER_DZ = 1.E3_RP ! thickness of buffer region [m]: z
  real(RP), private, parameter :: BUFFER_DX = 0.0_RP ! thickness of buffer region [m]: x
  real(RP), private, parameter :: BUFFER_DY = 0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, parameter :: BUFFFACT  = 1.1_RP ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX  =  12 ! # of computational cells: z
  integer, private, parameter :: IMAX  =  15 ! # of computational cells: x
  integer, private, parameter :: JMAX  =  15 ! # of computational cells: y

  integer, private, parameter :: IBLOCK = IMAX !
  integer, private, parameter :: JBLOCK = JMAX !

  include 'inc_index_all.h'
