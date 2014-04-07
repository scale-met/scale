!-------------------------------------------------------------------------------
!> module TRACER / dry
!!
!! @par Description
!!          Tracer dry module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)   [new] imported from inc_tracer_dry.f90
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_dry
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
  public :: TRACER_dry_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ scale-les tracer parameters (vapor only)
  !
  !-----------------------------------------------------------------------------
  integer, public, parameter :: QA = 1

  integer, public, parameter :: I_QV =  1
  integer, public, parameter :: I_QC =  0
  integer, public, parameter :: I_QR =  0
  integer, public, parameter :: I_QI =  0
  integer, public, parameter :: I_QS =  0
  integer, public, parameter :: I_QG =  0
  integer, public, parameter :: I_NC =  0
  integer, public, parameter :: I_NR =  0
  integer, public, parameter :: I_NI =  0
  integer, public, parameter :: I_NS =  0
  integer, public, parameter :: I_NG =  0

  integer, public, parameter :: QQA =  1 ! mass tracer (water)
  integer, public, parameter :: QQS =  1 ! start index for mass tracer
  integer, public, parameter :: QQE =  1 ! end   index for mass tracer

  integer, public, parameter :: QWS =  0 ! start index for water tracer
  integer, public, parameter :: QWE = -1 ! end   index for water tracer
  integer, public, parameter :: QIS =  0 ! start index for ice tracer
  integer, public, parameter :: QIE = -1 ! end   index for ice tracer

  character(len=H_SHORT), public :: AQ_NAME(QA)
  character(len=H_MID)  , public :: AQ_DESC(QA)
  character(len=H_SHORT), public :: AQ_UNIT(QA)

  data AQ_NAME / 'QV' /

  data AQ_DESC / 'Water Vapor mixing ratio' /

  data AQ_UNIT / 'kg/kg' /

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_dry+AE_dummy+RD_mstrnx)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: MP_QA = 1 ! number of hydrometeor tracer
  integer, public, parameter :: I_mp_dummy = 1

  integer, public :: I_MP2ALL(MP_QA)
  data I_MP2ALL / -999 / ! dummy

  integer, public :: I_MP2RD(MP_QA)
  data I_MP2RD / 1    / ! dummy

  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine TRACER_dry_setup
    implicit none

    return
  end subroutine TRACER_dry_setup

end module scale_tracer_dry
