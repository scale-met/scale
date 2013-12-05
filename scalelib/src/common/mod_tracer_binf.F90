!-------------------------------------------------------------------------------
!> module TRACER / binf
!!
!! @par Description
!!          Tracer binf module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_binf.f90
!!
!<
!-------------------------------------------------------------------------------
module mod_tracer_binf
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TRACER_binf_setup
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !------------------------------------------------------------
  !
  !++ scale3 grid parameter for inc_tracer_hbinf_h33_c20.f90
  !
  !------------------------------------------------------------
  integer, public, parameter :: nbin =  33
  integer, public, parameter :: nccn =  0
  integer, public, parameter :: nspc =   7
  integer, public, parameter :: QA   = 232

  integer, public, parameter :: I_QV =  1
  integer, public, parameter :: I_QC =  2
  integer, public, parameter :: I_QR =  3
  integer, public, parameter :: I_QI =  4
  integer, public, parameter :: I_QS =  5
  integer, public, parameter :: I_QG =  6
  integer, public, parameter :: I_NC =  8
  integer, public, parameter :: I_NR =  9
  integer, public, parameter :: I_NI =  10
  integer, public, parameter :: I_NS =  11
  integer, public, parameter :: I_NG =  12

  integer, public, parameter :: QQA  = 232
  integer, public, parameter :: QQS  = 1
  integer, public, parameter :: QQE  = 232

  integer, public, parameter :: QWS =  2 ! start index for water tracer
  integer, public, parameter :: QWE = 34 ! end   index for water tracer
  integer, public, parameter :: QIS = 35 ! start index for ice tracer
  integer, public, parameter :: QIE = 232 ! end   index for ice tracer

  character(len=16), public :: AQ_NAME(QA)
  character(len=64), public :: AQ_DESC(QA)
  character(len=16), public :: AQ_UNIT(QA)

  data AQ_NAME(  1)  / 'QV' /
  data AQ_NAME(  2)  / 'Qc1' /
  data AQ_NAME(  3)  / 'Qc2' /
  data AQ_NAME(  4)  / 'Qc3' /
  data AQ_NAME(  5)  / 'Qc4' /
  data AQ_NAME(  6)  / 'Qc5' /
  data AQ_NAME(  7)  / 'Qc6' /
  data AQ_NAME(  8)  / 'Qc7' /
  data AQ_NAME(  9)  / 'Qc8' /
  data AQ_NAME( 10)  / 'Qc9' /
  data AQ_NAME( 11)  / 'Qc10' /
  data AQ_NAME( 12)  / 'Qc11' /
  data AQ_NAME( 13)  / 'Qc12' /
  data AQ_NAME( 14)  / 'Qc13' /
  data AQ_NAME( 15)  / 'Qc14' /
  data AQ_NAME( 16)  / 'Qc15' /
  data AQ_NAME( 17)  / 'Qc16' /
  data AQ_NAME( 18)  / 'Qc17' /
  data AQ_NAME( 19)  / 'Qc18' /
  data AQ_NAME( 20)  / 'Qc19' /
  data AQ_NAME( 21)  / 'Qc20' /
  data AQ_NAME( 22)  / 'Qc21' /
  data AQ_NAME( 23)  / 'Qc22' /
  data AQ_NAME( 24)  / 'Qc23' /
  data AQ_NAME( 25)  / 'Qc24' /
  data AQ_NAME( 26)  / 'Qc25' /
  data AQ_NAME( 27)  / 'Qc26' /
  data AQ_NAME( 28)  / 'Qc27' /
  data AQ_NAME( 29)  / 'Qc28' /
  data AQ_NAME( 30)  / 'Qc29' /
  data AQ_NAME( 31)  / 'Qc30' /
  data AQ_NAME( 32)  / 'Qc31' /
  data AQ_NAME( 33)  / 'Qc32' /
  data AQ_NAME( 34)  / 'Qc33' /
  data AQ_NAME( 35)  / 'Qi1' /
  data AQ_NAME( 36)  / 'Qi2' /
  data AQ_NAME( 37)  / 'Qi3' /
  data AQ_NAME( 38)  / 'Qi4' /
  data AQ_NAME( 39)  / 'Qi5' /
  data AQ_NAME( 40)  / 'Qi6' /
  data AQ_NAME( 41)  / 'Qi7' /
  data AQ_NAME( 42)  / 'Qi8' /
  data AQ_NAME( 43)  / 'Qi9' /
  data AQ_NAME( 44)  / 'Qi10' /
  data AQ_NAME( 45)  / 'Qi11' /
  data AQ_NAME( 46)  / 'Qi12' /
  data AQ_NAME( 47)  / 'Qi13' /
  data AQ_NAME( 48)  / 'Qi14' /
  data AQ_NAME( 49)  / 'Qi15' /
  data AQ_NAME( 50)  / 'Qi16' /
  data AQ_NAME( 51)  / 'Qi17' /
  data AQ_NAME( 52)  / 'Qi18' /
  data AQ_NAME( 53)  / 'Qi19' /
  data AQ_NAME( 54)  / 'Qi20' /
  data AQ_NAME( 55)  / 'Qi21' /
  data AQ_NAME( 56)  / 'Qi22' /
  data AQ_NAME( 57)  / 'Qi23' /
  data AQ_NAME( 58)  / 'Qi24' /
  data AQ_NAME( 59)  / 'Qi25' /
  data AQ_NAME( 60)  / 'Qi26' /
  data AQ_NAME( 61)  / 'Qi27' /
  data AQ_NAME( 62)  / 'Qi28' /
  data AQ_NAME( 63)  / 'Qi29' /
  data AQ_NAME( 64)  / 'Qi30' /
  data AQ_NAME( 65)  / 'Qi31' /
  data AQ_NAME( 66)  / 'Qi32' /
  data AQ_NAME( 67)  / 'Qi33' /
  data AQ_NAME( 68)  / 'Qp1' /
  data AQ_NAME( 69)  / 'Qp2' /
  data AQ_NAME( 70)  / 'Qp3' /
  data AQ_NAME( 71)  / 'Qp4' /
  data AQ_NAME( 72)  / 'Qp5' /
  data AQ_NAME( 73)  / 'Qp6' /
  data AQ_NAME( 74)  / 'Qp7' /
  data AQ_NAME( 75)  / 'Qp8' /
  data AQ_NAME( 76)  / 'Qp9' /
  data AQ_NAME( 77)  / 'Qp10' /
  data AQ_NAME( 78)  / 'Qp11' /
  data AQ_NAME( 79)  / 'Qp12' /
  data AQ_NAME( 80)  / 'Qp13' /
  data AQ_NAME( 81)  / 'Qp14' /
  data AQ_NAME( 82)  / 'Qp15' /
  data AQ_NAME( 83)  / 'Qp16' /
  data AQ_NAME( 84)  / 'Qp17' /
  data AQ_NAME( 85)  / 'Qp18' /
  data AQ_NAME( 86)  / 'Qp19' /
  data AQ_NAME( 87)  / 'Qp20' /
  data AQ_NAME( 88)  / 'Qp21' /
  data AQ_NAME( 89)  / 'Qp22' /
  data AQ_NAME( 90)  / 'Qp23' /
  data AQ_NAME( 91)  / 'Qp24' /
  data AQ_NAME( 92)  / 'Qp25' /
  data AQ_NAME( 93)  / 'Qp26' /
  data AQ_NAME( 94)  / 'Qp27' /
  data AQ_NAME( 95)  / 'Qp28' /
  data AQ_NAME( 96)  / 'Qp29' /
  data AQ_NAME( 97)  / 'Qp30' /
  data AQ_NAME( 98)  / 'Qp31' /
  data AQ_NAME( 99)  / 'Qp32' /
  data AQ_NAME(100)  / 'Qp33' /
  data AQ_NAME(101)  / 'Qd1' /
  data AQ_NAME(102)  / 'Qd2' /
  data AQ_NAME(103)  / 'Qd3' /
  data AQ_NAME(104)  / 'Qd4' /
  data AQ_NAME(105)  / 'Qd5' /
  data AQ_NAME(106)  / 'Qd6' /
  data AQ_NAME(107)  / 'Qd7' /
  data AQ_NAME(108)  / 'Qd8' /
  data AQ_NAME(109)  / 'Qd9' /
  data AQ_NAME(110)  / 'Qd10' /
  data AQ_NAME(111)  / 'Qd11' /
  data AQ_NAME(112)  / 'Qd12' /
  data AQ_NAME(113)  / 'Qd13' /
  data AQ_NAME(114)  / 'Qd14' /
  data AQ_NAME(115)  / 'Qd15' /
  data AQ_NAME(116)  / 'Qd16' /
  data AQ_NAME(117)  / 'Qd17' /
  data AQ_NAME(118)  / 'Qd18' /
  data AQ_NAME(119)  / 'Qd19' /
  data AQ_NAME(120)  / 'Qd20' /
  data AQ_NAME(121)  / 'Qd21' /
  data AQ_NAME(122)  / 'Qd22' /
  data AQ_NAME(123)  / 'Qd23' /
  data AQ_NAME(124)  / 'Qd24' /
  data AQ_NAME(125)  / 'Qd25' /
  data AQ_NAME(126)  / 'Qd26' /
  data AQ_NAME(127)  / 'Qd27' /
  data AQ_NAME(128)  / 'Qd28' /
  data AQ_NAME(129)  / 'Qd29' /
  data AQ_NAME(130)  / 'Qd30' /
  data AQ_NAME(131)  / 'Qd31' /
  data AQ_NAME(132)  / 'Qd32' /
  data AQ_NAME(133)  / 'Qd33' /
  data AQ_NAME(134)  / 'Qs1' /
  data AQ_NAME(135)  / 'Qs2' /
  data AQ_NAME(136)  / 'Qs3' /
  data AQ_NAME(137)  / 'Qs4' /
  data AQ_NAME(138)  / 'Qs5' /
  data AQ_NAME(139)  / 'Qs6' /
  data AQ_NAME(140)  / 'Qs7' /
  data AQ_NAME(141)  / 'Qs8' /
  data AQ_NAME(142)  / 'Qs9' /
  data AQ_NAME(143)  / 'Qs10' /
  data AQ_NAME(144)  / 'Qs11' /
  data AQ_NAME(145)  / 'Qs12' /
  data AQ_NAME(146)  / 'Qs13' /
  data AQ_NAME(147)  / 'Qs14' /
  data AQ_NAME(148)  / 'Qs15' /
  data AQ_NAME(149)  / 'Qs16' /
  data AQ_NAME(150)  / 'Qs17' /
  data AQ_NAME(151)  / 'Qs18' /
  data AQ_NAME(152)  / 'Qs19' /
  data AQ_NAME(153)  / 'Qs20' /
  data AQ_NAME(154)  / 'Qs21' /
  data AQ_NAME(155)  / 'Qs22' /
  data AQ_NAME(156)  / 'Qs23' /
  data AQ_NAME(157)  / 'Qs24' /
  data AQ_NAME(158)  / 'Qs25' /
  data AQ_NAME(159)  / 'Qs26' /
  data AQ_NAME(160)  / 'Qs27' /
  data AQ_NAME(161)  / 'Qs28' /
  data AQ_NAME(162)  / 'Qs29' /
  data AQ_NAME(163)  / 'Qs30' /
  data AQ_NAME(164)  / 'Qs31' /
  data AQ_NAME(165)  / 'Qs32' /
  data AQ_NAME(166)  / 'Qs33' /
  data AQ_NAME(167)  / 'Qg1' /
  data AQ_NAME(168)  / 'Qg2' /
  data AQ_NAME(169)  / 'Qg3' /
  data AQ_NAME(170)  / 'Qg4' /
  data AQ_NAME(171)  / 'Qg5' /
  data AQ_NAME(172)  / 'Qg6' /
  data AQ_NAME(173)  / 'Qg7' /
  data AQ_NAME(174)  / 'Qg8' /
  data AQ_NAME(175)  / 'Qg9' /
  data AQ_NAME(176)  / 'Qg10' /
  data AQ_NAME(177)  / 'Qg11' /
  data AQ_NAME(178)  / 'Qg12' /
  data AQ_NAME(179)  / 'Qg13' /
  data AQ_NAME(180)  / 'Qg14' /
  data AQ_NAME(181)  / 'Qg15' /
  data AQ_NAME(182)  / 'Qg16' /
  data AQ_NAME(183)  / 'Qg17' /
  data AQ_NAME(184)  / 'Qg18' /
  data AQ_NAME(185)  / 'Qg19' /
  data AQ_NAME(186)  / 'Qg20' /
  data AQ_NAME(187)  / 'Qg21' /
  data AQ_NAME(188)  / 'Qg22' /
  data AQ_NAME(189)  / 'Qg23' /
  data AQ_NAME(190)  / 'Qg24' /
  data AQ_NAME(191)  / 'Qg25' /
  data AQ_NAME(192)  / 'Qg26' /
  data AQ_NAME(193)  / 'Qg27' /
  data AQ_NAME(194)  / 'Qg28' /
  data AQ_NAME(195)  / 'Qg29' /
  data AQ_NAME(196)  / 'Qg30' /
  data AQ_NAME(197)  / 'Qg31' /
  data AQ_NAME(198)  / 'Qg32' /
  data AQ_NAME(199)  / 'Qg33' /
  data AQ_NAME(200)  / 'Qh1' /
  data AQ_NAME(201)  / 'Qh2' /
  data AQ_NAME(202)  / 'Qh3' /
  data AQ_NAME(203)  / 'Qh4' /
  data AQ_NAME(204)  / 'Qh5' /
  data AQ_NAME(205)  / 'Qh6' /
  data AQ_NAME(206)  / 'Qh7' /
  data AQ_NAME(207)  / 'Qh8' /
  data AQ_NAME(208)  / 'Qh9' /
  data AQ_NAME(209)  / 'Qh10' /
  data AQ_NAME(210)  / 'Qh11' /
  data AQ_NAME(211)  / 'Qh12' /
  data AQ_NAME(212)  / 'Qh13' /
  data AQ_NAME(213)  / 'Qh14' /
  data AQ_NAME(214)  / 'Qh15' /
  data AQ_NAME(215)  / 'Qh16' /
  data AQ_NAME(216)  / 'Qh17' /
  data AQ_NAME(217)  / 'Qh18' /
  data AQ_NAME(218)  / 'Qh19' /
  data AQ_NAME(219)  / 'Qh20' /
  data AQ_NAME(220)  / 'Qh21' /
  data AQ_NAME(221)  / 'Qh22' /
  data AQ_NAME(222)  / 'Qh23' /
  data AQ_NAME(223)  / 'Qh24' /
  data AQ_NAME(224)  / 'Qh25' /
  data AQ_NAME(225)  / 'Qh26' /
  data AQ_NAME(226)  / 'Qh27' /
  data AQ_NAME(227)  / 'Qh28' /
  data AQ_NAME(228)  / 'Qh29' /
  data AQ_NAME(229)  / 'Qh30' /
  data AQ_NAME(230)  / 'Qh31' /
  data AQ_NAME(231)  / 'Qh32' /
  data AQ_NAME(232)  / 'Qh33' /

  data AQ_DESC(  1)  / 'Water Vapor mixing ratio' /
  data AQ_DESC(  2)  / 'Mixing ratio of cloud bin1' /
  data AQ_DESC(  3)  / 'Mixing ratio of cloud bin2' /
  data AQ_DESC(  4)  / 'Mixing ratio of cloud bin3' /
  data AQ_DESC(  5)  / 'Mixing ratio of cloud bin4' /
  data AQ_DESC(  6)  / 'Mixing ratio of cloud bin5' /
  data AQ_DESC(  7)  / 'Mixing ratio of cloud bin6' /
  data AQ_DESC(  8)  / 'Mixing ratio of cloud bin7' /
  data AQ_DESC(  9)  / 'Mixing ratio of cloud bin8' /
  data AQ_DESC( 10)  / 'Mixing ratio of cloud bin9' /
  data AQ_DESC( 11)  / 'Mixing ratio of cloud bin10' /
  data AQ_DESC( 12)  / 'Mixing ratio of cloud bin11' /
  data AQ_DESC( 13)  / 'Mixing ratio of cloud bin12' /
  data AQ_DESC( 14)  / 'Mixing ratio of cloud bin13' /
  data AQ_DESC( 15)  / 'Mixing ratio of cloud bin14' /
  data AQ_DESC( 16)  / 'Mixing ratio of cloud bin15' /
  data AQ_DESC( 17)  / 'Mixing ratio of cloud bin16' /
  data AQ_DESC( 18)  / 'Mixing ratio of cloud bin17' /
  data AQ_DESC( 19)  / 'Mixing ratio of cloud bin18' /
  data AQ_DESC( 20)  / 'Mixing ratio of cloud bin19' /
  data AQ_DESC( 21)  / 'Mixing ratio of cloud bin20' /
  data AQ_DESC( 22)  / 'Mixing ratio of cloud bin21' /
  data AQ_DESC( 23)  / 'Mixing ratio of cloud bin22' /
  data AQ_DESC( 24)  / 'Mixing ratio of cloud bin23' /
  data AQ_DESC( 25)  / 'Mixing ratio of cloud bin24' /
  data AQ_DESC( 26)  / 'Mixing ratio of cloud bin25' /
  data AQ_DESC( 27)  / 'Mixing ratio of cloud bin26' /
  data AQ_DESC( 28)  / 'Mixing ratio of cloud bin27' /
  data AQ_DESC( 29)  / 'Mixing ratio of cloud bin28' /
  data AQ_DESC( 30)  / 'Mixing ratio of cloud bin29' /
  data AQ_DESC( 31)  / 'Mixing ratio of cloud bin30' /
  data AQ_DESC( 32)  / 'Mixing ratio of cloud bin31' /
  data AQ_DESC( 33)  / 'Mixing ratio of cloud bin32' /
  data AQ_DESC( 34)  / 'Mixing ratio of cloud bin33' /
  data AQ_DESC( 35)  / 'Mixing ratio of colmn bin1' /
  data AQ_DESC( 36)  / 'Mixing ratio of colmn bin2' /
  data AQ_DESC( 37)  / 'Mixing ratio of colmn bin3' /
  data AQ_DESC( 38)  / 'Mixing ratio of colmn bin4' /
  data AQ_DESC( 39)  / 'Mixing ratio of colmn bin5' /
  data AQ_DESC( 40)  / 'Mixing ratio of colmn bin6' /
  data AQ_DESC( 41)  / 'Mixing ratio of colmn bin7' /
  data AQ_DESC( 42)  / 'Mixing ratio of colmn bin8' /
  data AQ_DESC( 43)  / 'Mixing ratio of colmn bin9' /
  data AQ_DESC( 44)  / 'Mixing ratio of colmn bin10' /
  data AQ_DESC( 45)  / 'Mixing ratio of colmn bin11' /
  data AQ_DESC( 46)  / 'Mixing ratio of colmn bin12' /
  data AQ_DESC( 47)  / 'Mixing ratio of colmn bin13' /
  data AQ_DESC( 48)  / 'Mixing ratio of colmn bin14' /
  data AQ_DESC( 49)  / 'Mixing ratio of colmn bin15' /
  data AQ_DESC( 50)  / 'Mixing ratio of colmn bin16' /
  data AQ_DESC( 51)  / 'Mixing ratio of colmn bin17' /
  data AQ_DESC( 52)  / 'Mixing ratio of colmn bin18' /
  data AQ_DESC( 53)  / 'Mixing ratio of colmn bin19' /
  data AQ_DESC( 54)  / 'Mixing ratio of colmn bin20' /
  data AQ_DESC( 55)  / 'Mixing ratio of colmn bin21' /
  data AQ_DESC( 56)  / 'Mixing ratio of colmn bin22' /
  data AQ_DESC( 57)  / 'Mixing ratio of colmn bin23' /
  data AQ_DESC( 58)  / 'Mixing ratio of colmn bin24' /
  data AQ_DESC( 59)  / 'Mixing ratio of colmn bin25' /
  data AQ_DESC( 60)  / 'Mixing ratio of colmn bin26' /
  data AQ_DESC( 61)  / 'Mixing ratio of colmn bin27' /
  data AQ_DESC( 62)  / 'Mixing ratio of colmn bin28' /
  data AQ_DESC( 63)  / 'Mixing ratio of colmn bin29' /
  data AQ_DESC( 64)  / 'Mixing ratio of colmn bin30' /
  data AQ_DESC( 65)  / 'Mixing ratio of colmn bin31' /
  data AQ_DESC( 66)  / 'Mixing ratio of colmn bin32' /
  data AQ_DESC( 67)  / 'Mixing ratio of colmn bin33' /
  data AQ_DESC( 68)  / 'Mixing ratio of plate bin1' /
  data AQ_DESC( 69)  / 'Mixing ratio of plate bin2' /
  data AQ_DESC( 70)  / 'Mixing ratio of plate bin3' /
  data AQ_DESC( 71)  / 'Mixing ratio of plate bin4' /
  data AQ_DESC( 72)  / 'Mixing ratio of plate bin5' /
  data AQ_DESC( 73)  / 'Mixing ratio of plate bin6' /
  data AQ_DESC( 74)  / 'Mixing ratio of plate bin7' /
  data AQ_DESC( 75)  / 'Mixing ratio of plate bin8' /
  data AQ_DESC( 76)  / 'Mixing ratio of plate bin9' /
  data AQ_DESC( 77)  / 'Mixing ratio of plate bin10' /
  data AQ_DESC( 78)  / 'Mixing ratio of plate bin11' /
  data AQ_DESC( 79)  / 'Mixing ratio of plate bin12' /
  data AQ_DESC( 80)  / 'Mixing ratio of plate bin13' /
  data AQ_DESC( 81)  / 'Mixing ratio of plate bin14' /
  data AQ_DESC( 82)  / 'Mixing ratio of plate bin15' /
  data AQ_DESC( 83)  / 'Mixing ratio of plate bin16' /
  data AQ_DESC( 84)  / 'Mixing ratio of plate bin17' /
  data AQ_DESC( 85)  / 'Mixing ratio of plate bin18' /
  data AQ_DESC( 86)  / 'Mixing ratio of plate bin19' /
  data AQ_DESC( 87)  / 'Mixing ratio of plate bin20' /
  data AQ_DESC( 88)  / 'Mixing ratio of plate bin21' /
  data AQ_DESC( 89)  / 'Mixing ratio of plate bin22' /
  data AQ_DESC( 90)  / 'Mixing ratio of plate bin23' /
  data AQ_DESC( 91)  / 'Mixing ratio of plate bin24' /
  data AQ_DESC( 92)  / 'Mixing ratio of plate bin25' /
  data AQ_DESC( 93)  / 'Mixing ratio of plate bin26' /
  data AQ_DESC( 94)  / 'Mixing ratio of plate bin27' /
  data AQ_DESC( 95)  / 'Mixing ratio of plate bin28' /
  data AQ_DESC( 96)  / 'Mixing ratio of plate bin29' /
  data AQ_DESC( 97)  / 'Mixing ratio of plate bin30' /
  data AQ_DESC( 98)  / 'Mixing ratio of plate bin31' /
  data AQ_DESC( 99)  / 'Mixing ratio of plate bin32' /
  data AQ_DESC(100)  / 'Mixing ratio of plate bin33' /
  data AQ_DESC(101)  / 'Mixing ratio of dend. bin1' /
  data AQ_DESC(102)  / 'Mixing ratio of dend. bin2' /
  data AQ_DESC(103)  / 'Mixing ratio of dend. bin3' /
  data AQ_DESC(104)  / 'Mixing ratio of dend. bin4' /
  data AQ_DESC(105)  / 'Mixing ratio of dend. bin5' /
  data AQ_DESC(106)  / 'Mixing ratio of dend. bin6' /
  data AQ_DESC(107)  / 'Mixing ratio of dend. bin7' /
  data AQ_DESC(108)  / 'Mixing ratio of dend. bin8' /
  data AQ_DESC(109)  / 'Mixing ratio of dend. bin9' /
  data AQ_DESC(110)  / 'Mixing ratio of dend. bin10' /
  data AQ_DESC(111)  / 'Mixing ratio of dend. bin11' /
  data AQ_DESC(112)  / 'Mixing ratio of dend. bin12' /
  data AQ_DESC(113)  / 'Mixing ratio of dend. bin13' /
  data AQ_DESC(114)  / 'Mixing ratio of dend. bin14' /
  data AQ_DESC(115)  / 'Mixing ratio of dend. bin15' /
  data AQ_DESC(116)  / 'Mixing ratio of dend. bin16' /
  data AQ_DESC(117)  / 'Mixing ratio of dend. bin17' /
  data AQ_DESC(118)  / 'Mixing ratio of dend. bin18' /
  data AQ_DESC(119)  / 'Mixing ratio of dend. bin19' /
  data AQ_DESC(120)  / 'Mixing ratio of dend. bin20' /
  data AQ_DESC(121)  / 'Mixing ratio of dend. bin21' /
  data AQ_DESC(122)  / 'Mixing ratio of dend. bin22' /
  data AQ_DESC(123)  / 'Mixing ratio of dend. bin23' /
  data AQ_DESC(124)  / 'Mixing ratio of dend. bin24' /
  data AQ_DESC(125)  / 'Mixing ratio of dend. bin25' /
  data AQ_DESC(126)  / 'Mixing ratio of dend. bin26' /
  data AQ_DESC(127)  / 'Mixing ratio of dend. bin27' /
  data AQ_DESC(128)  / 'Mixing ratio of dend. bin28' /
  data AQ_DESC(129)  / 'Mixing ratio of dend. bin29' /
  data AQ_DESC(130)  / 'Mixing ratio of dend. bin30' /
  data AQ_DESC(131)  / 'Mixing ratio of dend. bin31' /
  data AQ_DESC(132)  / 'Mixing ratio of dend. bin32' /
  data AQ_DESC(133)  / 'Mixing ratio of dend. bin33' /
  data AQ_DESC(134)  / 'Mixing ratio of snow  bin1' /
  data AQ_DESC(135)  / 'Mixing ratio of snow  bin2' /
  data AQ_DESC(136)  / 'Mixing ratio of snow  bin3' /
  data AQ_DESC(137)  / 'Mixing ratio of snow  bin4' /
  data AQ_DESC(138)  / 'Mixing ratio of snow  bin5' /
  data AQ_DESC(139)  / 'Mixing ratio of snow  bin6' /
  data AQ_DESC(140)  / 'Mixing ratio of snow  bin7' /
  data AQ_DESC(141)  / 'Mixing ratio of snow  bin8' /
  data AQ_DESC(142)  / 'Mixing ratio of snow  bin9' /
  data AQ_DESC(143)  / 'Mixing ratio of snow  bin10' /
  data AQ_DESC(144)  / 'Mixing ratio of snow  bin11' /
  data AQ_DESC(145)  / 'Mixing ratio of snow  bin12' /
  data AQ_DESC(146)  / 'Mixing ratio of snow  bin13' /
  data AQ_DESC(147)  / 'Mixing ratio of snow  bin14' /
  data AQ_DESC(148)  / 'Mixing ratio of snow  bin15' /
  data AQ_DESC(149)  / 'Mixing ratio of snow  bin16' /
  data AQ_DESC(150)  / 'Mixing ratio of snow  bin17' /
  data AQ_DESC(151)  / 'Mixing ratio of snow  bin18' /
  data AQ_DESC(152)  / 'Mixing ratio of snow  bin19' /
  data AQ_DESC(153)  / 'Mixing ratio of snow  bin20' /
  data AQ_DESC(154)  / 'Mixing ratio of snow  bin21' /
  data AQ_DESC(155)  / 'Mixing ratio of snow  bin22' /
  data AQ_DESC(156)  / 'Mixing ratio of snow  bin23' /
  data AQ_DESC(157)  / 'Mixing ratio of snow  bin24' /
  data AQ_DESC(158)  / 'Mixing ratio of snow  bin25' /
  data AQ_DESC(159)  / 'Mixing ratio of snow  bin26' /
  data AQ_DESC(160)  / 'Mixing ratio of snow  bin27' /
  data AQ_DESC(161)  / 'Mixing ratio of snow  bin28' /
  data AQ_DESC(162)  / 'Mixing ratio of snow  bin29' /
  data AQ_DESC(163)  / 'Mixing ratio of snow  bin30' /
  data AQ_DESC(164)  / 'Mixing ratio of snow  bin31' /
  data AQ_DESC(165)  / 'Mixing ratio of snow  bin32' /
  data AQ_DESC(166)  / 'Mixing ratio of snow  bin33' /
  data AQ_DESC(167)  / 'Mixing ratio of grpl  bin1' /
  data AQ_DESC(168)  / 'Mixing ratio of grpl  bin2' /
  data AQ_DESC(169)  / 'Mixing ratio of grpl  bin3' /
  data AQ_DESC(170)  / 'Mixing ratio of grpl  bin4' /
  data AQ_DESC(171)  / 'Mixing ratio of grpl  bin5' /
  data AQ_DESC(172)  / 'Mixing ratio of grpl  bin6' /
  data AQ_DESC(173)  / 'Mixing ratio of grpl  bin7' /
  data AQ_DESC(174)  / 'Mixing ratio of grpl  bin8' /
  data AQ_DESC(175)  / 'Mixing ratio of grpl  bin9' /
  data AQ_DESC(176)  / 'Mixing ratio of grpl  bin10' /
  data AQ_DESC(177)  / 'Mixing ratio of grpl  bin11' /
  data AQ_DESC(178)  / 'Mixing ratio of grpl  bin12' /
  data AQ_DESC(179)  / 'Mixing ratio of grpl  bin13' /
  data AQ_DESC(180)  / 'Mixing ratio of grpl  bin14' /
  data AQ_DESC(181)  / 'Mixing ratio of grpl  bin15' /
  data AQ_DESC(182)  / 'Mixing ratio of grpl  bin16' /
  data AQ_DESC(183)  / 'Mixing ratio of grpl  bin17' /
  data AQ_DESC(184)  / 'Mixing ratio of grpl  bin18' /
  data AQ_DESC(185)  / 'Mixing ratio of grpl  bin19' /
  data AQ_DESC(186)  / 'Mixing ratio of grpl  bin20' /
  data AQ_DESC(187)  / 'Mixing ratio of grpl  bin21' /
  data AQ_DESC(188)  / 'Mixing ratio of grpl  bin22' /
  data AQ_DESC(189)  / 'Mixing ratio of grpl  bin23' /
  data AQ_DESC(190)  / 'Mixing ratio of grpl  bin24' /
  data AQ_DESC(191)  / 'Mixing ratio of grpl  bin25' /
  data AQ_DESC(192)  / 'Mixing ratio of grpl  bin26' /
  data AQ_DESC(193)  / 'Mixing ratio of grpl  bin27' /
  data AQ_DESC(194)  / 'Mixing ratio of grpl  bin28' /
  data AQ_DESC(195)  / 'Mixing ratio of grpl  bin29' /
  data AQ_DESC(196)  / 'Mixing ratio of grpl  bin30' /
  data AQ_DESC(197)  / 'Mixing ratio of grpl  bin31' /
  data AQ_DESC(198)  / 'Mixing ratio of grpl  bin32' /
  data AQ_DESC(199)  / 'Mixing ratio of grpl  bin33' /
  data AQ_DESC(200)  / 'Mixing ratio of hail  bin1' /
  data AQ_DESC(201)  / 'Mixing ratio of hail  bin2' /
  data AQ_DESC(202)  / 'Mixing ratio of hail  bin3' /
  data AQ_DESC(203)  / 'Mixing ratio of hail  bin4' /
  data AQ_DESC(204)  / 'Mixing ratio of hail  bin5' /
  data AQ_DESC(205)  / 'Mixing ratio of hail  bin6' /
  data AQ_DESC(206)  / 'Mixing ratio of hail  bin7' /
  data AQ_DESC(207)  / 'Mixing ratio of hail  bin8' /
  data AQ_DESC(208)  / 'Mixing ratio of hail  bin9' /
  data AQ_DESC(209)  / 'Mixing ratio of hail  bin10' /
  data AQ_DESC(210)  / 'Mixing ratio of hail  bin11' /
  data AQ_DESC(211)  / 'Mixing ratio of hail  bin12' /
  data AQ_DESC(212)  / 'Mixing ratio of hail  bin13' /
  data AQ_DESC(213)  / 'Mixing ratio of hail  bin14' /
  data AQ_DESC(214)  / 'Mixing ratio of hail  bin15' /
  data AQ_DESC(215)  / 'Mixing ratio of hail  bin16' /
  data AQ_DESC(216)  / 'Mixing ratio of hail  bin17' /
  data AQ_DESC(217)  / 'Mixing ratio of hail  bin18' /
  data AQ_DESC(218)  / 'Mixing ratio of hail  bin19' /
  data AQ_DESC(219)  / 'Mixing ratio of hail  bin20' /
  data AQ_DESC(220)  / 'Mixing ratio of hail  bin21' /
  data AQ_DESC(221)  / 'Mixing ratio of hail  bin22' /
  data AQ_DESC(222)  / 'Mixing ratio of hail  bin23' /
  data AQ_DESC(223)  / 'Mixing ratio of hail  bin24' /
  data AQ_DESC(224)  / 'Mixing ratio of hail  bin25' /
  data AQ_DESC(225)  / 'Mixing ratio of hail  bin26' /
  data AQ_DESC(226)  / 'Mixing ratio of hail  bin27' /
  data AQ_DESC(227)  / 'Mixing ratio of hail  bin28' /
  data AQ_DESC(228)  / 'Mixing ratio of hail  bin29' /
  data AQ_DESC(229)  / 'Mixing ratio of hail  bin30' /
  data AQ_DESC(230)  / 'Mixing ratio of hail  bin31' /
  data AQ_DESC(231)  / 'Mixing ratio of hail  bin32' /
  data AQ_DESC(232)  / 'Mixing ratio of hail  bin33' /

  data AQ_UNIT(  1)  / 'kg/kg' /
  data AQ_UNIT(  2)  / 'kg/kg/unit logr' /
  data AQ_UNIT(  3)  / 'kg/kg/unit logr' /
  data AQ_UNIT(  4)  / 'kg/kg/unit logr' /
  data AQ_UNIT(  5)  / 'kg/kg/unit logr' /
  data AQ_UNIT(  6)  / 'kg/kg/unit logr' /
  data AQ_UNIT(  7)  / 'kg/kg/unit logr' /
  data AQ_UNIT(  8)  / 'kg/kg/unit logr' /
  data AQ_UNIT(  9)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 10)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 11)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 12)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 13)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 14)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 15)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 16)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 17)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 18)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 19)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 20)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 21)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 22)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 23)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 24)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 25)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 26)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 27)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 28)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 29)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 30)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 31)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 32)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 33)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 34)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 35)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 36)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 37)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 38)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 39)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 40)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 41)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 42)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 43)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 44)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 45)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 46)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 47)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 48)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 49)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 50)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 51)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 52)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 53)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 54)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 55)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 56)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 57)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 58)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 59)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 60)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 61)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 62)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 63)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 64)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 65)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 66)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 67)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 68)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 69)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 70)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 71)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 72)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 73)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 74)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 75)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 76)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 77)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 78)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 79)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 80)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 81)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 82)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 83)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 84)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 85)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 86)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 87)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 88)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 89)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 90)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 91)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 92)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 93)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 94)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 95)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 96)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 97)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 98)  / 'kg/kg/unit logr' /
  data AQ_UNIT( 99)  / 'kg/kg/unit logr' /
  data AQ_UNIT(100)  / 'kg/kg/unit logr' /
  data AQ_UNIT(101)  / 'kg/kg/unit logr' /
  data AQ_UNIT(102)  / 'kg/kg/unit logr' /
  data AQ_UNIT(103)  / 'kg/kg/unit logr' /
  data AQ_UNIT(104)  / 'kg/kg/unit logr' /
  data AQ_UNIT(105)  / 'kg/kg/unit logr' /
  data AQ_UNIT(106)  / 'kg/kg/unit logr' /
  data AQ_UNIT(107)  / 'kg/kg/unit logr' /
  data AQ_UNIT(108)  / 'kg/kg/unit logr' /
  data AQ_UNIT(109)  / 'kg/kg/unit logr' /
  data AQ_UNIT(110)  / 'kg/kg/unit logr' /
  data AQ_UNIT(111)  / 'kg/kg/unit logr' /
  data AQ_UNIT(112)  / 'kg/kg/unit logr' /
  data AQ_UNIT(113)  / 'kg/kg/unit logr' /
  data AQ_UNIT(114)  / 'kg/kg/unit logr' /
  data AQ_UNIT(115)  / 'kg/kg/unit logr' /
  data AQ_UNIT(116)  / 'kg/kg/unit logr' /
  data AQ_UNIT(117)  / 'kg/kg/unit logr' /
  data AQ_UNIT(118)  / 'kg/kg/unit logr' /
  data AQ_UNIT(119)  / 'kg/kg/unit logr' /
  data AQ_UNIT(120)  / 'kg/kg/unit logr' /
  data AQ_UNIT(121)  / 'kg/kg/unit logr' /
  data AQ_UNIT(122)  / 'kg/kg/unit logr' /
  data AQ_UNIT(123)  / 'kg/kg/unit logr' /
  data AQ_UNIT(124)  / 'kg/kg/unit logr' /
  data AQ_UNIT(125)  / 'kg/kg/unit logr' /
  data AQ_UNIT(126)  / 'kg/kg/unit logr' /
  data AQ_UNIT(127)  / 'kg/kg/unit logr' /
  data AQ_UNIT(128)  / 'kg/kg/unit logr' /
  data AQ_UNIT(129)  / 'kg/kg/unit logr' /
  data AQ_UNIT(130)  / 'kg/kg/unit logr' /
  data AQ_UNIT(131)  / 'kg/kg/unit logr' /
  data AQ_UNIT(132)  / 'kg/kg/unit logr' /
  data AQ_UNIT(133)  / 'kg/kg/unit logr' /
  data AQ_UNIT(134)  / 'kg/kg/unit logr' /
  data AQ_UNIT(135)  / 'kg/kg/unit logr' /
  data AQ_UNIT(136)  / 'kg/kg/unit logr' /
  data AQ_UNIT(137)  / 'kg/kg/unit logr' /
  data AQ_UNIT(138)  / 'kg/kg/unit logr' /
  data AQ_UNIT(139)  / 'kg/kg/unit logr' /
  data AQ_UNIT(140)  / 'kg/kg/unit logr' /
  data AQ_UNIT(141)  / 'kg/kg/unit logr' /
  data AQ_UNIT(142)  / 'kg/kg/unit logr' /
  data AQ_UNIT(143)  / 'kg/kg/unit logr' /
  data AQ_UNIT(144)  / 'kg/kg/unit logr' /
  data AQ_UNIT(145)  / 'kg/kg/unit logr' /
  data AQ_UNIT(146)  / 'kg/kg/unit logr' /
  data AQ_UNIT(147)  / 'kg/kg/unit logr' /
  data AQ_UNIT(148)  / 'kg/kg/unit logr' /
  data AQ_UNIT(149)  / 'kg/kg/unit logr' /
  data AQ_UNIT(150)  / 'kg/kg/unit logr' /
  data AQ_UNIT(151)  / 'kg/kg/unit logr' /
  data AQ_UNIT(152)  / 'kg/kg/unit logr' /
  data AQ_UNIT(153)  / 'kg/kg/unit logr' /
  data AQ_UNIT(154)  / 'kg/kg/unit logr' /
  data AQ_UNIT(155)  / 'kg/kg/unit logr' /
  data AQ_UNIT(156)  / 'kg/kg/unit logr' /
  data AQ_UNIT(157)  / 'kg/kg/unit logr' /
  data AQ_UNIT(158)  / 'kg/kg/unit logr' /
  data AQ_UNIT(159)  / 'kg/kg/unit logr' /
  data AQ_UNIT(160)  / 'kg/kg/unit logr' /
  data AQ_UNIT(161)  / 'kg/kg/unit logr' /
  data AQ_UNIT(162)  / 'kg/kg/unit logr' /
  data AQ_UNIT(163)  / 'kg/kg/unit logr' /
  data AQ_UNIT(164)  / 'kg/kg/unit logr' /
  data AQ_UNIT(165)  / 'kg/kg/unit logr' /
  data AQ_UNIT(166)  / 'kg/kg/unit logr' /
  data AQ_UNIT(167)  / 'kg/kg/unit logr' /
  data AQ_UNIT(168)  / 'kg/kg/unit logr' /
  data AQ_UNIT(169)  / 'kg/kg/unit logr' /
  data AQ_UNIT(170)  / 'kg/kg/unit logr' /
  data AQ_UNIT(171)  / 'kg/kg/unit logr' /
  data AQ_UNIT(172)  / 'kg/kg/unit logr' /
  data AQ_UNIT(173)  / 'kg/kg/unit logr' /
  data AQ_UNIT(174)  / 'kg/kg/unit logr' /
  data AQ_UNIT(175)  / 'kg/kg/unit logr' /
  data AQ_UNIT(176)  / 'kg/kg/unit logr' /
  data AQ_UNIT(177)  / 'kg/kg/unit logr' /
  data AQ_UNIT(178)  / 'kg/kg/unit logr' /
  data AQ_UNIT(179)  / 'kg/kg/unit logr' /
  data AQ_UNIT(180)  / 'kg/kg/unit logr' /
  data AQ_UNIT(181)  / 'kg/kg/unit logr' /
  data AQ_UNIT(182)  / 'kg/kg/unit logr' /
  data AQ_UNIT(183)  / 'kg/kg/unit logr' /
  data AQ_UNIT(184)  / 'kg/kg/unit logr' /
  data AQ_UNIT(185)  / 'kg/kg/unit logr' /
  data AQ_UNIT(186)  / 'kg/kg/unit logr' /
  data AQ_UNIT(187)  / 'kg/kg/unit logr' /
  data AQ_UNIT(188)  / 'kg/kg/unit logr' /
  data AQ_UNIT(189)  / 'kg/kg/unit logr' /
  data AQ_UNIT(190)  / 'kg/kg/unit logr' /
  data AQ_UNIT(191)  / 'kg/kg/unit logr' /
  data AQ_UNIT(192)  / 'kg/kg/unit logr' /
  data AQ_UNIT(193)  / 'kg/kg/unit logr' /
  data AQ_UNIT(194)  / 'kg/kg/unit logr' /
  data AQ_UNIT(195)  / 'kg/kg/unit logr' /
  data AQ_UNIT(196)  / 'kg/kg/unit logr' /
  data AQ_UNIT(197)  / 'kg/kg/unit logr' /
  data AQ_UNIT(198)  / 'kg/kg/unit logr' /
  data AQ_UNIT(199)  / 'kg/kg/unit logr' /
  data AQ_UNIT(200)  / 'kg/kg/unit logr' /
  data AQ_UNIT(201)  / 'kg/kg/unit logr' /
  data AQ_UNIT(202)  / 'kg/kg/unit logr' /
  data AQ_UNIT(203)  / 'kg/kg/unit logr' /
  data AQ_UNIT(204)  / 'kg/kg/unit logr' /
  data AQ_UNIT(205)  / 'kg/kg/unit logr' /
  data AQ_UNIT(206)  / 'kg/kg/unit logr' /
  data AQ_UNIT(207)  / 'kg/kg/unit logr' /
  data AQ_UNIT(208)  / 'kg/kg/unit logr' /
  data AQ_UNIT(209)  / 'kg/kg/unit logr' /
  data AQ_UNIT(210)  / 'kg/kg/unit logr' /
  data AQ_UNIT(211)  / 'kg/kg/unit logr' /
  data AQ_UNIT(212)  / 'kg/kg/unit logr' /
  data AQ_UNIT(213)  / 'kg/kg/unit logr' /
  data AQ_UNIT(214)  / 'kg/kg/unit logr' /
  data AQ_UNIT(215)  / 'kg/kg/unit logr' /
  data AQ_UNIT(216)  / 'kg/kg/unit logr' /
  data AQ_UNIT(217)  / 'kg/kg/unit logr' /
  data AQ_UNIT(218)  / 'kg/kg/unit logr' /
  data AQ_UNIT(219)  / 'kg/kg/unit logr' /
  data AQ_UNIT(220)  / 'kg/kg/unit logr' /
  data AQ_UNIT(221)  / 'kg/kg/unit logr' /
  data AQ_UNIT(222)  / 'kg/kg/unit logr' /
  data AQ_UNIT(223)  / 'kg/kg/unit logr' /
  data AQ_UNIT(224)  / 'kg/kg/unit logr' /
  data AQ_UNIT(225)  / 'kg/kg/unit logr' /
  data AQ_UNIT(226)  / 'kg/kg/unit logr' /
  data AQ_UNIT(227)  / 'kg/kg/unit logr' /
  data AQ_UNIT(228)  / 'kg/kg/unit logr' /
  data AQ_UNIT(229)  / 'kg/kg/unit logr' /
  data AQ_UNIT(230)  / 'kg/kg/unit logr' /
  data AQ_UNIT(231)  / 'kg/kg/unit logr' /
  data AQ_UNIT(232)  / 'kg/kg/unit logr' /

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_binf+AE_dummy+RD_mstrnX)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: MP_QA = 7 ! number of hydrometeor tracer
  integer, public, parameter :: I_mp_QC = 1
  integer, public, parameter :: I_mp_QP = 2
  integer, public, parameter :: I_mp_QCL= 3
  integer, public, parameter :: I_mp_QD = 4
  integer, public, parameter :: I_mp_QS = 5
  integer, public, parameter :: I_mp_QG = 6
  integer, public, parameter :: I_mp_QH = 7

  integer, public :: I_MP2ALL(MP_QA)  !--- dummy for Bin model
  data I_MP2ALL / I_mp_QC,   & ! start of Cloud bin
                  I_mp_QP,   & ! start of plate bin
                  I_mp_QCL,  & ! start of columner bin
                  I_mp_QD,   & ! start of dendrite bin
                  I_mp_QS,   & ! start of snow bin
                  I_mp_QG,   & ! start of graupel bin
                  I_mp_QH    / ! start of hail bin

  integer, public :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1,    & ! I_mp_QC => MSTRN_nptype=1: water cloud
                  2,    & ! I_mp_QP => MSTRN_nptype=2: ice cloud (plate)
                  2,    & ! I_mp_QCL=> MSTRN_nptype=2: ice cloud (columner)
                  2,    & ! I_mp_QD => MSTRN_nptype=2: ice cloud (dendrite)
                  2,    & ! I_mp_QS => MSTRN_nptype=2: ice cloud (snow)
                  2,    & ! I_mp_QG => MSTRN_nptype=2: ice cloud (graupel)
                  2     / ! I_mp_QH => MSTRN_nptype=2: ice cloud (hail)
  
  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust


contains

  subroutine TRACER_binf_setup
    implicit none

    return
  end subroutine TRACER_binf_setup

end module mod_tracer_binf
