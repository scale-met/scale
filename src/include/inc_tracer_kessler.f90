  !-----------------------------------------------------------------------------
  !
  !++ scale3 tracer parameters (1-moment bulk 3 category)
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
