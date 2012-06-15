  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters
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

  character(len=16), private, save :: AQ_NAME(QA)
  character(len=64), private, save :: AQ_DESC(QA)
  character(len=16), private, save :: AQ_UNIT(QA)

  data AQ_NAME / 'QV' /

  data AQ_DESC / 'Water Vapor mixing ratio' /

  data AQ_UNIT / 'kg/kg' /
