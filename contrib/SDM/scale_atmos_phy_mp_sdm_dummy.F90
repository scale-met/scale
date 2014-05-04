!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Dummy Cloud Microphysics for sdm atmosphere
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-10-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_mp_sdm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_tracer_sdm
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_sdm_setup
  public :: ATMOS_PHY_MP_sdm

  public :: ATMOS_PHY_MP_sdm_CloudFraction
  public :: ATMOS_PHY_MP_sdm_EffectiveRadius
  public :: ATMOS_PHY_MP_sdm_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]

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
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** SDM not supported.'
    if( IO_L ) write(IO_FID_LOG,*) '*** Please contact SCALE developers'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm( &
          DENS, &
          MOMZ, &
          MOMX, &
          MOMY, &
          RHOT, &
          QTRC  )
    use scale_process, only: &
       PRC_MPIstop
    use scale_tracer, only: &
       QAD => QA
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** SDM not supported.'
    if( IO_L ) write(IO_FID_LOG,*) '*** Please contact SCALE developers'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_MP_sdm

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sdm_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_tracer, only: &
       QAD => QA
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QAD)
    !---------------------------------------------------------------------------

    cldfrac(:,:,:) = 0.0_RP ! dummy

    return
  end subroutine ATMOS_PHY_MP_sdm_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_sdm_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0  )
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
    !---------------------------------------------------------------------------

    Re(:,:,:,:) = 8.E-6_RP ! dummy

    return
  end subroutine ATMOS_PHY_MP_sdm_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_Mixingratio( &
       Qe,    &
       QTRC0  )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,MP_QAD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]

    integer  :: ihydro
    !---------------------------------------------------------------------------

    do ihydro = 1, MP_QA
       Qe(:,:,:,ihydro) = 8.E-6_RP ! dummy
    enddo

    return

  end subroutine ATMOS_PHY_MP_sdm_Mixingratio
  !-----------------------------------------------------------------------------

end module scale_atmos_phy_mp_sdm
