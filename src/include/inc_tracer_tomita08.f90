  !-----------------------------------------------------------------------------
  !
  !++ scale3 tracer parameters (1-moment bulk 6 category)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: QA = 6

  integer, private, parameter :: I_QV =  1
  integer, private, parameter :: I_QC =  2
  integer, private, parameter :: I_QR =  3
  integer, private, parameter :: I_QI =  4
  integer, private, parameter :: I_QS =  5
  integer, private, parameter :: I_QG =  6
  integer, private, parameter :: I_NC =  0
  integer, private, parameter :: I_NR =  0
  integer, private, parameter :: I_NI =  0
  integer, private, parameter :: I_NS =  0
  integer, private, parameter :: I_NG =  0

  integer, private, parameter :: QQA =  6 ! mass tracer (water)
  integer, private, parameter :: QQS =  1 ! start index for mass tracer
  integer, private, parameter :: QQE =  6 ! end   index for mass tracer

  integer, private, parameter :: QWS =  2 ! start index for water tracer
  integer, private, parameter :: QWE =  3 ! end   index for water tracer
  integer, private, parameter :: QIS =  4 ! start index for ice tracer
  integer, private, parameter :: QIE =  6 ! end   index for ice tracer

  character(len=16), private, save :: AQ_NAME(QA)
  character(len=64), private, save :: AQ_DESC(QA)
  character(len=16), private, save :: AQ_UNIT(QA)

  data AQ_NAME / 'QV', &
                 'QC', &
                 'QR', &
                 'QI', &
                 'QS', &
                 'QG'  /

  data AQ_DESC / 'Water Vapor mixing ratio',   &
                 'Cloud Water mixing ratio',   &
                 'Rain Water mixing ratio',    &
                 'Cloud Ice mixing ratio',     &
                 'Snow mixing ratio',          &
                 'Graupel mixing ratio'        /

  data AQ_UNIT / 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg'   /

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_tomita08+AE_dummy+RD_mstrnX)
  !
  !-----------------------------------------------------------------------------

  integer, private, parameter :: MP_QA = 5 ! number of hydrometeor tracer
  integer, private, parameter :: I_mp_QC = 1
  integer, private, parameter :: I_mp_QR = 2
  integer, private, parameter :: I_mp_QI = 3
  integer, private, parameter :: I_mp_QS = 4
  integer, private, parameter :: I_mp_QG = 5

  integer, private, save :: I_MP2ALL(MP_QA)
  data I_MP2ALL / I_QC, & ! I_mp_QC => I_QC
                  I_QR, & ! I_mp_QR => I_QR
                  I_QI, & ! I_mp_QI => I_QI
                  I_QS, & ! I_mp_QS => I_QS
                  I_QG  / ! I_mp_QG => I_QG

  integer, private, save :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1,    & ! I_mp_QC => MSTRN_nptype=1: water cloud
                  1,    & ! I_mp_QR => MSTRN_nptype=1: water cloud
                  2,    & ! I_mp_QI => MSTRN_nptype=2: ice cloud
                  2,    & ! I_mp_QS => MSTRN_nptype=2: ice cloud
                  2     / ! I_mp_QG => MSTRN_nptype=2: ice cloud

  integer, private, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, private, parameter :: I_ae_dummy = 1

  integer, private, save :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, private, save :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust
