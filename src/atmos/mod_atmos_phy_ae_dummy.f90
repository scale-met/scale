!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          dummy code
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-06 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_ae
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_setup
  public :: ATMOS_PHY_AE

  public :: ATMOS_PHY_AE_EffectiveRadius

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: AE_DENS(AE_QA) = 0.0_RP ! aerosol density [kg/m3]=[g/L]

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
  !> Setup
  subroutine ATMOS_PHY_AE_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_AE
    implicit none

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-AE]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ dummy aerosol process'

    if ( ATMOS_TYPE_PHY_AE /= 'NONE' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_AE is not NONE. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_AE_setup

  !-----------------------------------------------------------------------------
  !> Aerosol Microphysics
  subroutine ATMOS_PHY_AE
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Aerosol(dummy)'

    return
  end subroutine ATMOS_PHY_AE

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_AE_EffectiveRadius( &
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

!    Re(:,:,:,I_ae_seasalt) = 2.E-4_RP
!    Re(:,:,:,I_ae_dust   ) = 4.E-6_RP
!    Re(:,:,:,I_ae_bc     ) = 4.D-8_RP
!    Re(:,:,:,I_ae_oc     ) = RH(:,:,:)
!    Re(:,:,:,I_ae_sulfate) = RH(:,:,:)

    return
  end subroutine ATMOS_PHY_AE_EffectiveRadius

end module mod_atmos_phy_ae 
