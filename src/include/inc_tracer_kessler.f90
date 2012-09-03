
  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters
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
