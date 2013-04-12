  !-----------------------------------------------------------------------------
  !
  !++ scale3 tracer parameters (2-moment bulk 6 category)
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: QA = 11

  integer, private, parameter :: I_QV =  1
  integer, private, parameter :: I_QC =  2
  integer, private, parameter :: I_QR =  3
  integer, private, parameter :: I_QI =  4
  integer, private, parameter :: I_QS =  5
  integer, private, parameter :: I_QG =  6
  integer, private, parameter :: I_NC =  7
  integer, private, parameter :: I_NR =  8
  integer, private, parameter :: I_NI =  9
  integer, private, parameter :: I_NS = 10
  integer, private, parameter :: I_NG = 11

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
                 'QG', &
                 'NC', &
                 'NR', &
                 'NI', &
                 'NS', &
                 'NG'  /

  data AQ_DESC / 'Water Vapor mixing ratio',   &
                 'Cloud Water mixing ratio',   &
                 'Rain Water mixing ratio',    &
                 'Cloud Ice mixing ratio',     &
                 'Snow mixing ratio',          &
                 'Graupel mixing ratio',       &
                 'Cloud Water Number Density', &
                 'Rain Water Number Density',  &
                 'Cloud Ice Number Density',   &
                 'Snow Number Density',        &
                 'Graupel Number Density'      /

  data AQ_UNIT / 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'num/kg', &
                 'num/kg', &
                 'num/kg', &
                 'num/kg', &
                 'num/kg'  /
