  !-----------------------------------------------------------------------------
  !
  !++ scale-les tracer parameters (1-moment bulk 3 category)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: QA = 3

  integer, private, parameter :: I_QV =  1
  integer, private, parameter :: I_QC =  2
  integer, private, parameter :: I_QR =  3
  integer, private, parameter :: I_QI =  0
  integer, private, parameter :: I_QS =  0
  integer, private, parameter :: I_QG =  0
  integer, private, parameter :: I_NC =  0
  integer, private, parameter :: I_NR =  0
  integer, private, parameter :: I_NI =  0
  integer, private, parameter :: I_NS =  0
  integer, private, parameter :: I_NG =  0

  integer, private, parameter :: QQA =  3 ! mass tracer (water)
  integer, private, parameter :: QQS =  1 ! start index for mass tracer
  integer, private, parameter :: QQE =  3 ! end   index for mass tracer

  integer, private, parameter :: QWS =  2 ! start index for water tracer
  integer, private, parameter :: QWE =  3 ! end   index for water tracer
  integer, private, parameter :: QIS =  0 ! start index for ice tracer
  integer, private, parameter :: QIE =  0 ! end   index for ice tracer

  character(len=16), private, save :: AQ_NAME(QA)
  character(len=64), private, save :: AQ_DESC(QA)
  character(len=16), private, save :: AQ_UNIT(QA)

  data AQ_NAME / 'QV', &
                 'QC', &
                 'QR'  /

  data AQ_DESC / 'Water Vapor mixing ratio', &
                 'Cloud Water mixing ratio', &
                 'Rain Water mixing ratio'   /

  data AQ_UNIT / 'kg/kg', &
                 'kg/kg', &
                 'kg/kg'  /

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_kessler+AE_dummy+RD_mstrnX)
  !
  !-----------------------------------------------------------------------------

  integer, private, parameter :: MP_QA = 2 ! number of hydrometeor tracer
  integer, private, parameter :: I_mp_QC = 1
  integer, private, parameter :: I_mp_QR = 2

  integer, private, save :: I_MP2ALL(MP_QA)
  data I_MP2ALL / I_QC, & ! I_mp_QC => I_QC
                  I_QR  / ! I_mp_QR => I_QR

  integer, private, save :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1,    & ! I_mp_QC => MSTRN_nptype=1: water cloud
                  1     / ! I_mp_QR => MSTRN_nptype=1: water cloud

  integer, private, save :: I_MP_BIN_NUM(MP_QA) !-- bin number
  data I_MP_BIN_NUM     &
                / 1,    &
                  1     /

  integer, private, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, private, parameter :: I_ae_dummy = 1

  integer, private, save :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, private, save :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust
