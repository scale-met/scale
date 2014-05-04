!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Dummy Cloud Microphysics for dry atmosphere
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-10-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_mp_dry
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_tracer_dry
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_dry_setup
  public :: ATMOS_PHY_MP_dry

  public :: ATMOS_PHY_MP_dry_CloudFraction
  public :: ATMOS_PHY_MP_dry_EffectiveRadius
  public :: ATMOS_PHY_MP_dry_Mixingratio

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
  subroutine ATMOS_PHY_MP_dry_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_tracer, only: &
       QAD => QA
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Dry Atmosphere'

    ATMOS_PHY_MP_DENS(:) = 0.0_RP

    if ( MP_TYPE /= 'DRY' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_MP_TYPE is not DRY. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_MP_dry_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_dry( &
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

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics(dummy)'

    return
  end subroutine ATMOS_PHY_MP_dry

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_dry_CloudFraction( &
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
  end subroutine ATMOS_PHY_MP_dry_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_dry_EffectiveRadius( &
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
  end subroutine ATMOS_PHY_MP_dry_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_dry_Mixingratio( &
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

  end subroutine ATMOS_PHY_MP_dry_Mixingratio
  !-----------------------------------------------------------------------------

end module scale_atmos_phy_mp_dry
