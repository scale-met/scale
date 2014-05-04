!-----------------------------------------------------------------------------
!!
!! @par Description
!!          Tracer parameters for Super Droplet Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-09-18 (Y.Sato)  [new]
!! @li      2014-01-22 (Y.Sato)  [new] revised for develop 4.0.0
!!
!<
!-------------------------------------------------------------------------------
module scale_tracer_sdm
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
  public :: TRACER_sdm_setup
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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

  character(len=H_SHORT), public, save :: AQ_NAME(QA)
  character(len=H_MID), public, save :: AQ_DESC(QA)
  character(len=H_SHORT), public, save :: AQ_UNIT(QA)

  data AQ_NAME / 'QV', 'QC', 'QR' /

  data AQ_DESC / 'Water Vapor mixing ratio', &
                 'cloud water mixing ratio', &  
                 'rain water mixing ratio'  /

  data AQ_UNIT / 'kg/kg', 'kg/kg', 'kg/kg' /

  !-----------------------------------------------------------------------------
  !
  !++ tracer index & relationship (MP_sn13+AE_dummy+RD_mstrnX)
  !
  !-----------------------------------------------------------------------------

  integer, public, parameter :: MP_QA = 2 ! number of hydrometeor tracer
  integer, public, parameter :: I_mp_QC = 1
  integer, private, parameter :: I_m_QR = 2

  integer, public, save :: I_MP2ALL(MP_QA)
  data I_MP2ALL / I_QC, & ! I_mp_QC => I_QC
                  I_QR  / ! I_mp_QR => I_QR

  integer, public, save :: I_MP2RD(MP_QA)
  data I_MP2RD  / 1,    & ! I_mp_QC => MSTRN_nptype=1: water cloud
                  1     / ! I_mp_QR => MSTRN_nptype=1: water cloud
  
  integer, public, parameter :: AE_QA = 1 ! number of aerosol tracer
  integer, public, parameter :: I_ae_dummy = 1

  integer, public, save :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 / ! dummy

  integer, public, save :: I_AE2RD(AE_QA)
  data I_AE2RD  / 3    / ! dummy => MSTRN_nptype=3: dust

contains
  subroutine TRACER_sdm_setup

   return
  end subroutine TRACER_sdm_setup

end module scale_tracer_sdm
