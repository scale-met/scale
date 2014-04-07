!-------------------------------------------------------------------------------
!> module TRACER / sn14
!!
!! @par Description
!!          Tracer sn14 module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_sn14.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_sn14
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: TRACER_sn14_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ scale-les tracer parameters (2-moment bulk 6 category)
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: QA = 11

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

  integer, public, parameter :: QQA =  6 ! mass tracer (water)
  integer, public, parameter :: QQS =  1 ! start index for mass tracer
  integer, public, parameter :: QQE =  6 ! end   index for mass tracer

  integer, public, parameter :: QWS =  2 ! start index for water tracer
  integer, public, parameter :: QWE =  3 ! end   index for water tracer
  integer, public, parameter :: QIS =  4 ! start index for ice tracer
  integer, public, parameter :: QIE =  6 ! end   index for ice tracer

  character(len=H_SHORT), public :: AQ_NAME(QA)
  character(len=H_MID)  , public :: AQ_DESC(QA)
  character(len=H_SHORT), public :: AQ_UNIT(QA)

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

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_sn14+AE_dummy+RD_mstrnX)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: MP_QA = 5 ! number of hydrometeor tracer
  integer, public, parameter :: I_mp_QC = 1
  integer, public, parameter :: I_mp_QR = 2
  integer, public, parameter :: I_mp_QI = 3
  integer, public, parameter :: I_mp_QS = 4
  integer, public, parameter :: I_mp_QG = 5

  integer, public :: I_MP2ALL(MP_QA)
  data I_MP2ALL / I_QC, & ! I_mp_QC => I_QC
                  I_QR, & ! I_mp_QR => I_QR
                  I_QI, & ! I_mp_QI => I_QI
                  I_QS, & ! I_mp_QS => I_QS
                  I_QG  / ! I_mp_QG => I_QG

  integer, public :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1,    & ! I_mp_QC => MSTRN_nptype=1: water cloud
                  1,    & ! I_mp_QR => MSTRN_nptype=1: water cloud
                  2,    & ! I_mp_QI => MSTRN_nptype=2: ice cloud
                  2,    & ! I_mp_QS => MSTRN_nptype=2: ice cloud
                  2     / ! I_mp_QG => MSTRN_nptype=2: ice cloud

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_sn14_setup
    implicit none

    return
  end subroutine TRACER_sn14_setup

end module scale_tracer_sn14
