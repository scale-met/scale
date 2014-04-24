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
module scale_atmos_phy_ae_dummy
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
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
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: AE_TYPE
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
  subroutine ATMOS_PHY_AE_dummy( &
          DENS, &
          MOMZ, &
          MOMX, &
          MOMY, &
          RHOT, &
          QTRC  )
    use scale_tracer, only: &
       QAD => QA
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Aerosol(dummy)'

    return
  end subroutine ATMOS_PHY_AE_dummy

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_AE_dummy_EffectiveRadius( &
       Re,   &
       QTRC, &
       RH    )
    use scale_const, only: &
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

end module scale_atmos_phy_ae_dummy
