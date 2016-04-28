  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !------------------------------------------------------------
  !
  !++ scale-les grid parameter for inc_aetracer_none
  !
  !------------------------------------------------------------
  character(len=H_SHORT), public, parameter :: AETRACER_TYPE = "DUMMY"

  integer, public :: NCTG = 0
  integer, public :: NKAP = 0
  integer, public :: NSIZ = 0
  integer, public :: QA_AE = 0
  integer, public :: QAES = 0
  integer, public :: QAEE = 0

  integer, public, parameter :: IC_MIX = 0
  integer, public, parameter :: IC_SEA = 0
  integer, public, parameter :: IC_DUS = 0

  character(len=H_SHORT), public, allocatable :: AQ_AE_NAME(:)
  character(len=H_MID)  , public, allocatable :: AQ_AE_DESC(:)
  character(len=H_SHORT), public, allocatable :: AQ_AE_UNIT(:)
  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_binf+AE_none+RD_mstrnx)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_none = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! none

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! none => MSTRN_nptype=3: dust
  integer :: m, n, ierr
