  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !------------------------------------------------------------
  !
  !++ grid parameter for inc_aetracer_none
  !
  !------------------------------------------------------------
  character(len=H_SHORT), public, parameter :: AETRACER_TYPE = "NONE"

  integer, public, parameter :: AE_CTG = 0
  integer, public, parameter :: GAS_CTG = 0
  integer, public, parameter :: N_ATR = 0
  integer, public            :: NKAP(AE_CTG)
  integer, public            :: NSIZ(AE_CTG)
  integer, public, parameter :: QA_AE = 0
  integer, public, parameter :: QAES = 0
  integer, public, parameter :: QAEE = 0

  integer, public, parameter :: IC_MIX = 0
  integer, public, parameter :: IC_SEA = 0
  integer, public, parameter :: IC_DUS = 0

  integer, public, parameter :: IG_H2SO4 =  0
  integer, public, parameter :: IG_CGAS  =  0

  character(len=H_SHORT), public :: AQ_AE_NAME(QA_AE)
  character(len=H_MID)  , public :: AQ_AE_DESC(QA_AE)
  character(len=H_SHORT), public :: AQ_AE_UNIT(QA_AE)
  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_binf+AE_none+RD_mstrnx)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! none

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! none => MSTRN_nptype=3: dust
