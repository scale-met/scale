!-------------------------------------------------------------------------------
!> module TRACER / binw
!!
!! @par Description
!!          Tracer binw module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_binw.f90
!!
!<
!-------------------------------------------------------------------------------
module mod_tracer_binw
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
  public :: TRACER_binw_setup
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ scale3 grid parameters for hbinw
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: nbin = 33
  integer, public, parameter :: nccn = 0
  integer, public, parameter :: QA = 34

  integer, public, parameter :: I_QV =  1
  integer, public, parameter :: I_QC =  2
  integer, public, parameter :: I_QR =  3
  integer, public, parameter :: I_QI =  4
  integer, public, parameter :: I_QS =  5
  integer, public, parameter :: I_QG =  6
  integer, public, parameter :: I_NC =  7
  integer, public, parameter :: I_NR =  8
  integer, public, parameter :: I_NI =  9
  integer, public, parameter :: I_NS = 10
  integer, public, parameter :: I_NG = 11

  integer, public, parameter :: QQA = 34 ! mass tracer (water)
  integer, public, parameter :: QQS = 1            ! start index for mass tracer
  integer, public, parameter :: QQE = 34  ! end   index for mass tracer

  integer, public, parameter :: QWS =  2 ! start index for water tracer
  integer, public, parameter :: QWE =  34 ! end   index for water tracer
  integer, public, parameter :: QIS =  0 ! start index for ice tracer
  integer, public, parameter :: QIE =  0 ! end   index for ice tracer

  character(len=16), public :: AQ_NAME(QA)
  character(len=64), public :: AQ_DESC(QA)
  character(len=16), public :: AQ_UNIT(QA)

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

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_binw+AE_dummy+RD_mstrnX)
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: MP_QA = 1 ! number of hydrometeor tracer
  integer, public, parameter :: I_mp_QC = 1

  integer, public :: I_MP2ALL(MP_QA)
  data I_MP2ALL / I_mp_QC / ! dummy (meaningless) for bin model

  integer, public :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1    / ! I_mp_QC => MSTRN_nptype=1: water cloud

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust

contains

  subroutine TRACER_binw_setup
    implicit none

    return
  end subroutine TRACER_binw_setup


end module mod_tracer_binw
