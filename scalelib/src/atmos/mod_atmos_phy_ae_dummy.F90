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
module mod_atmos_phy_ae_dummy
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_tracer
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
  public :: ATMOS_PHY_AE_dummy_setup
  public :: ATMOS_PHY_AE_dummy

  public :: ATMOS_PHY_AE_dummy_EffectiveRadius

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: AE_DENS(:) ! aerosol density [kg/m3]=[g/L]

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
  subroutine ATMOS_PHY_AE_dummy_setup( AE_TYPE )
    use mod_process, only: &
       PRC_MPIstop
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    implicit none
    character(len=IO_SYSCHR), intent(in) :: AE_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-AE]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ dummy aerosol process'

    if ( AE_TYPE /= 'DUMMY' .and. AE_TYPE /= 'NONE' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_AE_TYPE is not DUMMY. Check!'
       call PRC_MPIstop
    endif

    allocate( AE_DENS(AE_QA) )
    AE_DENS(:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_AE_dummy_setup

  !-----------------------------------------------------------------------------
  !> Aerosol Microphysics
  subroutine ATMOS_PHY_AE_dummy
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Aerosol(dummy)'

    return
  end subroutine ATMOS_PHY_AE_dummy

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_AE_dummy_EffectiveRadius( &
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
  end subroutine ATMOS_PHY_AE_dummy_EffectiveRadius

end module mod_atmos_phy_ae_dummy
