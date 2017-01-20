!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Super Droplet Method (SDM), dummy interface
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-10-18 (S.Nishizawa) [new]
!! @li      2015-09-08 (Y.Sato)  [mod] update for version SCALE 0.2.4
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
  use scale_atmos_hydrometeor, only: &
     N_HYD
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_sdm_config
  public :: ATMOS_PHY_MP_sdm_setup
  public :: ATMOS_PHY_MP_sdm

  public :: ATMOS_PHY_MP_sdm_CloudFraction
  public :: ATMOS_PHY_MP_sdm_EffectiveRadius
  public :: ATMOS_PHY_MP_sdm_MixingRatio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, parameter :: QA_MP  = 3

  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_sdm_NAME(QA_MP)
  character(len=H_MID),   public, target :: ATMOS_PHY_MP_sdm_DESC(QA_MP)
  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_sdm_UNIT(QA_MP)

  real(RP), public, target :: ATMOS_PHY_MP_sdm_DENS(N_HYD) ! hydrometeor density [kg/m3]=[g/L]

  data ATMOS_PHY_MP_sdm_NAME / 'QV', 'QC', 'QR' /

  data ATMOS_PHY_MP_sdm_DESC / &
       'Ratio of Water Vapor mass to total mass (Specific humidity)',   &
       'cloud water mixing ratio', &
       'rain water mixing ratio'  /

  data ATMOS_PHY_MP_sdm_UNIT / 'kg/kg', 'kg/kg', 'kg/kg' /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: I_mp_QC =  1
  integer, private, parameter :: I_mp_QR =  2
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Confif
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_config( &
       MP_TYPE, &
       QA, QS   )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in)  :: MP_TYPE
    integer,          intent(out) :: QA
    integer,          intent(out) :: QS
    !---------------------------------------------------------------------------

    write(*,*) '*** SDM not supported.'
    write(*,*) '*** Please contact SCALE developers'
    call PRC_MPIstop

    QA = 0
    QS = 0

    return
  end subroutine ATMOS_PHY_MP_sdm_config

  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none
    !---------------------------------------------------------------------------

    write(*,*) '*** SDM not supported.'
    write(*,*) '*** Please contact SCALE developers'
    call PRC_MPIstop

    ATMOS_PHY_MP_sdm_DENS(:) = UNDEF

    return
  end subroutine ATMOS_PHY_MP_sdm_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm( &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN,       &
       EVAPORATE, &
       SFLX_rain, &
       SFLX_snow  )
    use scale_grid_index
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_tracer, only: &
       QA
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)    :: CCN (KA,IA,JA)
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)   !---- evaporated cloud number concentration [/m3]
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)
    !---------------------------------------------------------------------------

    write(*,*) '*** SDM not supported.'
    write(*,*) '*** Please contact SCALE developers'
    call PRC_MPIstop

    EVAPORATE = UNDEF
    SFLX_rain = UNDEF
    SFLX_snow = UNDEF

    return
  end subroutine ATMOS_PHY_MP_sdm

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sdm_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA)
    !---------------------------------------------------------------------------

    cldfrac(:,:,:) = 0.0_RP ! dummy

    return
  end subroutine ATMOS_PHY_MP_sdm_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_sdm_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0, &
       TEMP0  )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! Density [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)       ! Temperatuer [K]
    !---------------------------------------------------------------------------

    Re(:,:,:,:) = 8.E-6_RP ! dummy

    return
  end subroutine ATMOS_PHY_MP_sdm_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_Mixingratio( &
       Qe,   &
       QTRC0 )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,N_HYD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]

    integer  :: ihydro
    !---------------------------------------------------------------------------

    Qe(:,:,:,:) = 8.E-6_RP ! dummy

    return
  end subroutine ATMOS_PHY_MP_sdm_Mixingratio

end module scale_atmos_phy_mp_sdm
