  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !------------------------------------------------------------
  !
  !++ grid parameter for inc_aetracer_kajino13
  !
  !------------------------------------------------------------
  character(len=H_SHORT), public, parameter :: AETRACER_TYPE = "KAJINO13"

  integer, public :: AE_CTG = 1
  integer, public, parameter :: GAS_CTG = 2
  integer, public, parameter :: N_ATR = 5
  integer, public, allocatable :: NKAP(:)
  integer, public, allocatable :: NSIZ(:)
  integer, public :: QA_AE
  integer, public :: QAES
  integer, public :: QAEE

  integer, public, parameter :: IC_MIX   =  1
  integer, public, parameter :: IC_SEA   =  2
  integer, public, parameter :: IC_DUS   =  3

  integer, public, parameter :: IG_H2SO4 =  1
  integer, public, parameter :: IG_CGAS  =  2

  character(len=H_SHORT), public, allocatable :: AQ_AE_NAME(:)
  character(len=H_MID)  , public, allocatable :: AQ_AE_DESC(:)
  character(len=H_SHORT), public, allocatable :: AQ_AE_UNIT(:)
  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_binf+AE_dummy+RD_mstrnx)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust
