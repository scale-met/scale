!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics Subset
!!
!! @par Description
!!          Subsets for Aerosol
!!          OFFLINE scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-06 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_aesub
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AEsub_setup
  public :: ATMOS_PHY_AEsub_AE2RD
  public :: ATMOS_PHY_AEsub_EffectiveRadius

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public, parameter :: AE_QA = 1

  integer,  public, parameter :: I_ae_seasalt = 1
  !integer,  public, parameter :: I_ae_dust    = 2
  !integer,  public, parameter :: I_ae_bc      = 3
  !integer,  public, parameter :: I_ae_oc      = 4
  !integer,  public, parameter :: I_ae_sulfate = 5

  integer,  public, save :: I_AE2ALL(AE_QA)
  data I_AE2ALL / -999 /

  real(RP), public, save :: AE_DENS(AE_QA) ! aerosol density [kg/m3]=[g/L]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine ATMOS_PHY_AEsub_setup
    implicit none
    !---------------------------------------------------------------------------

    AE_DENS(I_ae_seasalt) = 2200.0_RP
!    AE_DENS(I_ae_dust   ) = 2500.0_RP
!    AE_DENS(I_ae_bc     ) = 1250.0_RP
!    AE_DENS(I_ae_oc     ) = 1500.0_RP
!    AE_DENS(I_ae_sulfate) = 1770.0_RP

    return
  end subroutine ATMOS_PHY_AEsub_setup

  !-----------------------------------------------------------------------------
  !> Make look-up table between hydrometeor tracer and particle type in radiation scheme
  subroutine ATMOS_PHY_AEsub_AE2RD( &
       I_AE2RD )
    implicit none

    integer, intent(out) :: I_AE2RD(AE_QA)
    !---------------------------------------------------------------------------

    I_AE2RD(I_ae_seasalt) = 8
!    I_AE2RD(I_ae_dust   ) = 3
!    I_AE2RD(I_ae_bc     ) = 4
!    I_AE2RD(I_ae_oc     ) = 8
!    I_AE2RD(I_ae_sulfate) = 8

    return
  end subroutine ATMOS_PHY_AEsub_AE2RD

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_AEsub_EffectiveRadius( &
       Re,   &
       QTRC, &
       RH    )
    use mod_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    real(RP), intent(out) :: Re  (KA,IA,JA,AE_QA) ! effective radius
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: RH  (KA,IA,JA)       ! relative humidity         [0-1]
    !---------------------------------------------------------------------------

    Re(:,:,:,:) = UNDEF

    Re(:,:,:,I_ae_seasalt) = 2.E-4_RP
!    Re(:,:,:,I_ae_dust   ) = 4.E-6_RP
!    Re(:,:,:,I_ae_bc     ) = 4.D-8_RP
!    Re(:,:,:,I_ae_oc     ) = RH(:,:,:)
!    Re(:,:,:,I_ae_sulfate) = RH(:,:,:)

    return
  end subroutine ATMOS_PHY_AEsub_EffectiveRadius

end module mod_atmos_phy_aesub 
