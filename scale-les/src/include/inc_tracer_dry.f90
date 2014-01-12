  !-----------------------------------------------------------------------------
  !
  !++ scale-les tracer parameters (vapor only)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: QA = 1

  integer, private, parameter :: I_QV =  1
  integer, private, parameter :: I_QC =  0
  integer, private, parameter :: I_QR =  0
  integer, private, parameter :: I_QI =  0
  integer, private, parameter :: I_QS =  0
  integer, private, parameter :: I_QG =  0
  integer, private, parameter :: I_NC =  0
  integer, private, parameter :: I_NR =  0
  integer, private, parameter :: I_NI =  0
  integer, private, parameter :: I_NS =  0
  integer, private, parameter :: I_NG =  0

  integer, private, parameter :: QQA =  1 ! mass tracer (water)
  integer, private, parameter :: QQS =  1 ! start index for mass tracer
  integer, private, parameter :: QQE =  1 ! end   index for mass tracer

  integer, private, parameter :: QWS =  0 ! start index for water tracer
  integer, private, parameter :: QWE =  0 ! end   index for water tracer
  integer, private, parameter :: QIS =  0 ! start index for ice tracer
  integer, private, parameter :: QIE =  0 ! end   index for ice tracer

  character(len=16), private, save :: AQ_NAME(QA)
  character(len=64), private, save :: AQ_DESC(QA)
  character(len=16), private, save :: AQ_UNIT(QA)

  data AQ_NAME / 'QV' /

  data AQ_DESC / 'Water Vapor mixing ratio' /

  data AQ_UNIT / 'kg/kg' /

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_dry+AE_dummy+RD_mstrnX)
  !
  !-----------------------------------------------------------------------------

  integer, private, parameter :: MP_QA = 1 ! number of hydrometeor tracer
  integer, private, parameter :: I_mp_dummy = 1

  integer, private, save :: I_MP2ALL(MP_QA)
  data I_MP2ALL / -999 / ! dummy

  integer, private, save :: I_MP2RD(MP_QA)
  data I_MP2RD / 1    / ! dummy

  integer, private, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, private, parameter :: I_ae_dummy = 1

  integer, private, save :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, private, save :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust
