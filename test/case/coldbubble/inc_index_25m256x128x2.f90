  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters (25m res)
  !
  !-----------------------------------------------------------------------------
  real(RP), private, parameter :: DX  = 25 ! length in the main region [m]: x
  real(RP), private, parameter :: DY  = 25 ! length in the main region [m]: y
  real(RP), private, parameter :: DZ  = 25 ! length in the main region [m]: z

  real(RP), private, parameter :: BUFFER_DZ =  0.0_RP ! thickness of buffer region [m]: z
  real(RP), private, parameter :: BUFFER_DX =  0.0_RP ! thickness of buffer region [m]: x
  real(RP), private, parameter :: BUFFER_DY =  0.0_RP ! thickness of buffer region [m]: y
  real(RP), private, parameter :: BUFFFACT  =  1.0_RP ! strech factor for dx/dy/dz of buffer region

  integer, private, parameter :: KMAX = 256 ! # of computational cells: z
  integer, private, parameter :: IMAX = 128 ! # of computational cells: x
  integer, private, parameter :: JMAX =   2 ! # of computational cells: y

  include 'inc_index_all.h'
