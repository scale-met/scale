!-------------------------------------------------------------------------------
!> module TRACER / kessler
!!
!! @par Description
!!          Tracer kessler module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_kessler.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_kessler
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
  public :: TRACER_kessler_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ scale-les tracer parameters (1-moment bulk 3 category)
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: QA = 3

  integer, public, parameter :: I_QV =  1
  integer, public, parameter :: I_QC =  2
  integer, public, parameter :: I_QR =  3
  integer, public, parameter :: I_QI =  0
  integer, public, parameter :: I_QS =  0
  integer, public, parameter :: I_QG =  0
  integer, public, parameter :: I_NC =  0
  integer, public, parameter :: I_NR =  0
  integer, public, parameter :: I_NI =  0
  integer, public, parameter :: I_NS =  0
  integer, public, parameter :: I_NG =  0

  integer, public, parameter :: QQA =  3 ! mass tracer (water)
  integer, public, parameter :: QQS =  1 ! start index for mass tracer
  integer, public, parameter :: QQE =  3 ! end   index for mass tracer

  integer, public, parameter :: QWS =  2 ! start index for water tracer
  integer, public, parameter :: QWE =  3 ! end   index for water tracer
  integer, public, parameter :: QIS =  0 ! start index for ice tracer
  integer, public, parameter :: QIE = -1 ! end   index for ice tracer

  character(len=H_SHORT), public :: AQ_NAME(QA)
  character(len=H_MID)  , public :: AQ_DESC(QA)
  character(len=H_SHORT), public :: AQ_UNIT(QA)

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
  !++ tracer index & relationship (MP_kessler+AE_dummy+RD_mstrnx)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: MP_QA = 2 ! number of hydrometeor tracer
  integer, public, parameter :: I_mp_QC = 1
  integer, public, parameter :: I_mp_QR = 2

  integer, public :: I_MP2ALL(MP_QA)
  data I_MP2ALL / I_QC, & ! I_mp_QC => I_QC
                  I_QR  / ! I_mp_QR => I_QR

  integer, public :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1,    & ! I_mp_QC => MSTRN_nptype=1: water cloud
                  1     / ! I_mp_QR => MSTRN_nptype=1: water cloud

  integer, public :: I_MP_BIN_NUM(MP_QA) !-- bin number
  data I_MP_BIN_NUM     &
                / 1,    &
                  1     /

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_kessler_setup
    implicit none

    return
  end subroutine TRACER_kessler_setup

end module scale_tracer_kessler
