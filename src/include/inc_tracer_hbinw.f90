
  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters for hbinw
  !
  !-----------------------------------------------------------------------------
  integer, private, parameter :: nbin = 33
  integer, private, parameter :: nccn = 20
  integer, private, parameter :: QA = 54

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

  integer, private, parameter :: QQA = 34 ! mass tracer (water)
  integer, private, parameter :: QQS = 1            ! start index for mass tracer
  integer, private, parameter :: QQE = 34  ! end   index for mass tracer

  integer, private, parameter :: QWS =  2 ! start index for water tracer
  integer, private, parameter :: QWE =  34 ! end   index for water tracer
  integer, private, parameter :: QIS =  0 ! start index for ice tracer
  integer, private, parameter :: QIE =  0 ! end   index for ice tracer

  character(len=16), private, save :: AQ_NAME(QA)
  character(len=64), private, save :: AQ_DESC(QA)
  character(len=16), private, save :: AQ_UNIT(QA)

  data AQ_NAME(1)  / 'QV' /
  data AQ_NAME(2)  / 'Qc1' /
  data AQ_NAME(3)  / 'Qc2' /
  data AQ_NAME(4)  / 'Qc3' /
  data AQ_NAME(5)  / 'Qc4' /
  data AQ_NAME(6)  / 'Qc5' /
  data AQ_NAME(7)  / 'Qc6' /
  data AQ_NAME(8)  / 'Qc7' /
  data AQ_NAME(9)  / 'Qc8' /
  data AQ_NAME(10) / 'Qc9' /
  data AQ_NAME(11) / 'Qc10' /
  data AQ_NAME(12) / 'Qc11' /
  data AQ_NAME(13) / 'Qc12' /
  data AQ_NAME(14) / 'Qc13' /
  data AQ_NAME(15) / 'Qc14' /
  data AQ_NAME(16) / 'Qc15' /
  data AQ_NAME(17) / 'Qc16' /
  data AQ_NAME(18) / 'Qc17' /
  data AQ_NAME(19) / 'Qc18' /
  data AQ_NAME(20) / 'Qc19' /
  data AQ_NAME(21) / 'Qc20' /
  data AQ_NAME(22) / 'Qc21' /
  data AQ_NAME(23) / 'Qc22' /
  data AQ_NAME(24) / 'Qc23' /
  data AQ_NAME(25) / 'Qc24' /
  data AQ_NAME(26) / 'Qc25' /
  data AQ_NAME(27) / 'Qc26' /
  data AQ_NAME(28) / 'Qc27' /
  data AQ_NAME(29) / 'Qc28' /
  data AQ_NAME(30) / 'Qc29' /
  data AQ_NAME(31) / 'Qc30' /
  data AQ_NAME(32) / 'Qc31' /
  data AQ_NAME(33) / 'Qc32' /
  data AQ_NAME(34) / 'Qc33' /
  data AQ_NAME(35) / 'Qaer1' /
  data AQ_NAME(36) / 'Qaer2' /
  data AQ_NAME(37) / 'Qaer3' /
  data AQ_NAME(38) / 'Qaer4' /
  data AQ_NAME(39) / 'Qaer5' /
  data AQ_NAME(40) / 'Qaer6' /
  data AQ_NAME(41) / 'Qaer7' /
  data AQ_NAME(42) / 'Qaer8' /
  data AQ_NAME(43) / 'Qaer9' /
  data AQ_NAME(44) / 'Qaer10' /
  data AQ_NAME(45) / 'Qaer11' /
  data AQ_NAME(46) / 'Qaer12' /
  data AQ_NAME(47) / 'Qaer13' /
  data AQ_NAME(48) / 'Qaer14' /
  data AQ_NAME(49) / 'Qaer15' /
  data AQ_NAME(50) / 'Qaer16' /
  data AQ_NAME(51) / 'Qaer17' /
  data AQ_NAME(52) / 'Qaer18' /
  data AQ_NAME(53) / 'Qaer19' /
  data AQ_NAME(54) / 'Qaer20' /

  data AQ_DESC(1) / 'Water Vapor mixing ratio' / 
  data AQ_DESC(2) / 'Water mixing ratio of bin1' / 
  data AQ_DESC(3) / 'Water mixing ratio of bin2' / 
  data AQ_DESC(4) / 'Water mixing ratio of bin3' / 
  data AQ_DESC(5) / 'Water mixing ratio of bin4' / 
  data AQ_DESC(6) / 'Water mixing ratio of bin5' / 
  data AQ_DESC(7) / 'Water mixing ratio of bin6' / 
  data AQ_DESC(8) / 'Water mixing ratio of bin7' / 
  data AQ_DESC(9) / 'Water mixing ratio of bin8' / 
  data AQ_DESC(10)/ 'Water mixing ratio of bin9' / 
  data AQ_DESC(11)/ 'Water mixing ratio of bin10' / 
  data AQ_DESC(12)/ 'Water mixing ratio of bin11' / 
  data AQ_DESC(13)/ 'Water mixing ratio of bin12' / 
  data AQ_DESC(14)/ 'Water mixing ratio of bin13' / 
  data AQ_DESC(15)/ 'Water mixing ratio of bin14' / 
  data AQ_DESC(16)/ 'Water mixing ratio of bin15' / 
  data AQ_DESC(17)/ 'Water mixing ratio of bin16' / 
  data AQ_DESC(18)/ 'Water mixing ratio of bin17' / 
  data AQ_DESC(19)/ 'Water mixing ratio of bin18' / 
  data AQ_DESC(20)/ 'Water mixing ratio of bin19' / 
  data AQ_DESC(21)/ 'Water mixing ratio of bin20' / 
  data AQ_DESC(22)/ 'Water mixing ratio of bin21' / 
  data AQ_DESC(23)/ 'Water mixing ratio of bin22' / 
  data AQ_DESC(24)/ 'Water mixing ratio of bin23' / 
  data AQ_DESC(25)/ 'Water mixing ratio of bin24' / 
  data AQ_DESC(26)/ 'Water mixing ratio of bin25' / 
  data AQ_DESC(27)/ 'Water mixing ratio of bin26' / 
  data AQ_DESC(28)/ 'Water mixing ratio of bin27' / 
  data AQ_DESC(29)/ 'Water mixing ratio of bin28' / 
  data AQ_DESC(30)/ 'Water mixing ratio of bin29' / 
  data AQ_DESC(31)/ 'Water mixing ratio of bin30' / 
  data AQ_DESC(32)/ 'Water mixing ratio of bin31' / 
  data AQ_DESC(33)/ 'Water mixing ratio of bin32' / 
  data AQ_DESC(34)/ 'Water mixing ratio of bin33' / 
  data AQ_DESC(35)/ 'aerosol mixing ratio of bin1' / 
  data AQ_DESC(36)/ 'aerosol mixing ratio of bin2' / 
  data AQ_DESC(37)/ 'aerosol mixing ratio of bin3' / 
  data AQ_DESC(38)/ 'aerosol mixing ratio of bin4' / 
  data AQ_DESC(39)/ 'aerosol mixing ratio of bin5' / 
  data AQ_DESC(40)/ 'aerosol mixing ratio of bin6' / 
  data AQ_DESC(41)/ 'aerosol mixing ratio of bin7' / 
  data AQ_DESC(42)/ 'aerosol mixing ratio of bin8' / 
  data AQ_DESC(43)/ 'aerosol mixing ratio of bin9' / 
  data AQ_DESC(44)/ 'aerosol mixing ratio of bin10' / 
  data AQ_DESC(45)/ 'aerosol mixing ratio of bin11' / 
  data AQ_DESC(46)/ 'aerosol mixing ratio of bin12' / 
  data AQ_DESC(47)/ 'aerosol mixing ratio of bin13' / 
  data AQ_DESC(48)/ 'aerosol mixing ratio of bin14' / 
  data AQ_DESC(49)/ 'aerosol mixing ratio of bin15' / 
  data AQ_DESC(50)/ 'aerosol mixing ratio of bin16' / 
  data AQ_DESC(51)/ 'aerosol mixing ratio of bin17' / 
  data AQ_DESC(52)/ 'aerosol mixing ratio of bin18' / 
  data AQ_DESC(53)/ 'aerosol mixing ratio of bin19' / 
  data AQ_DESC(54)/ 'aerosol mixing ratio of bin20' / 

  data AQ_UNIT(1)  / 'kg/kg' /
  data AQ_UNIT(2)  / 'kg/kg/unit logr' /
  data AQ_UNIT(3)  / 'kg/kg/unit logr' /
  data AQ_UNIT(4)  / 'kg/kg/unit logr' /
  data AQ_UNIT(5)  / 'kg/kg/unit logr' /
  data AQ_UNIT(6)  / 'kg/kg/unit logr' /
  data AQ_UNIT(7)  / 'kg/kg/unit logr' /
  data AQ_UNIT(8)  / 'kg/kg/unit logr' /
  data AQ_UNIT(9)  / 'kg/kg/unit logr' /
  data AQ_UNIT(10)  / 'kg/kg/unit logr' /
  data AQ_UNIT(11)  / 'kg/kg/unit logr' /
  data AQ_UNIT(12)  / 'kg/kg/unit logr' /
  data AQ_UNIT(13)  / 'kg/kg/unit logr' /
  data AQ_UNIT(14)  / 'kg/kg/unit logr' /
  data AQ_UNIT(15)  / 'kg/kg/unit logr' /
  data AQ_UNIT(16)  / 'kg/kg/unit logr' /
  data AQ_UNIT(17)  / 'kg/kg/unit logr' /
  data AQ_UNIT(18)  / 'kg/kg/unit logr' /
  data AQ_UNIT(19)  / 'kg/kg/unit logr' /
  data AQ_UNIT(20)  / 'kg/kg/unit logr' /
  data AQ_UNIT(21)  / 'kg/kg/unit logr' /
  data AQ_UNIT(22)  / 'kg/kg/unit logr' /
  data AQ_UNIT(23)  / 'kg/kg/unit logr' /
  data AQ_UNIT(24)  / 'kg/kg/unit logr' /
  data AQ_UNIT(25)  / 'kg/kg/unit logr' /
  data AQ_UNIT(26)  / 'kg/kg/unit logr' /
  data AQ_UNIT(27)  / 'kg/kg/unit logr' /
  data AQ_UNIT(28)  / 'kg/kg/unit logr' /
  data AQ_UNIT(29)  / 'kg/kg/unit logr' /
  data AQ_UNIT(30)  / 'kg/kg/unit logr' /
  data AQ_UNIT(31)  / 'kg/kg/unit logr' /
  data AQ_UNIT(32)  / 'kg/kg/unit logr' /
  data AQ_UNIT(33)  / 'kg/kg/unit logr' /
  data AQ_UNIT(34)  / 'kg/kg/unit logr' /
  data AQ_UNIT(35) / 'kg/kg/unit logra' /
  data AQ_UNIT(36) / 'kg/kg/unit logra' /
  data AQ_UNIT(37) / 'kg/kg/unit logra' /
  data AQ_UNIT(38) / 'kg/kg/unit logra' /
  data AQ_UNIT(39) / 'kg/kg/unit logra' /
  data AQ_UNIT(40) / 'kg/kg/unit logra' /
  data AQ_UNIT(41) / 'kg/kg/unit logra' /
  data AQ_UNIT(42) / 'kg/kg/unit logra' /
  data AQ_UNIT(43) / 'kg/kg/unit logra' /
  data AQ_UNIT(44) / 'kg/kg/unit logra' /
  data AQ_UNIT(45) / 'kg/kg/unit logra' /
  data AQ_UNIT(46) / 'kg/kg/unit logra' /
  data AQ_UNIT(47) / 'kg/kg/unit logra' /
  data AQ_UNIT(48) / 'kg/kg/unit logra' /
  data AQ_UNIT(49) / 'kg/kg/unit logra' /
  data AQ_UNIT(50) / 'kg/kg/unit logra' /
  data AQ_UNIT(51) / 'kg/kg/unit logra' /
  data AQ_UNIT(52) / 'kg/kg/unit logra' /
  data AQ_UNIT(53) / 'kg/kg/unit logra' /
  data AQ_UNIT(54) / 'kg/kg/unit logra' /
