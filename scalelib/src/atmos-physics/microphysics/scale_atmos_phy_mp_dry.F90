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
!! @li      2015-09-08 (Y.Sato)      [add] Add evaporated cloud number concentration
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
  use scale_atmos_hydrometeor, only: &
     N_HYD
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_dry_config
  public :: ATMOS_PHY_MP_dry_setup
  public :: ATMOS_PHY_MP_dry

  public :: ATMOS_PHY_MP_dry_CloudFraction
  public :: ATMOS_PHY_MP_dry_EffectiveRadius
  public :: ATMOS_PHY_MP_dry_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: QA_MP  = 0

  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_dry_NAME(QA_MP)
  character(len=H_MID)  , public, target :: ATMOS_PHY_MP_dry_DESC(QA_MP)
  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_dry_UNIT(QA_MP)

  real(RP), public, target :: ATMOS_PHY_MP_dry_DENS(N_HYD)

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
  !> Configure
  subroutine ATMOS_PHY_MP_dry_config( &
       MP_TYPE, &
       QA, QS   )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry
    implicit none

    character(len=*), intent(in) :: MP_TYPE
    integer, intent(out) :: QA
    integer, intent(out) :: QS

    !---------------------------------------------------------------------------

    if ( MP_TYPE /= 'DRY' ) then
       write(*,*) 'xxx ATMOS_PHY_MP_TYPE is not DRY. Check!'
       call PRC_MPIstop
    endif

    QS = -1
    QA = QA_MP

    return
  end subroutine ATMOS_PHY_MP_dry_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_dry_setup
    implicit none

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Cloud Microphysics] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** dummy process (dry Atmosphere)'

    return
  end subroutine ATMOS_PHY_MP_dry_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  subroutine ATMOS_PHY_MP_dry( &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN ,      &
       EVAPORATE, &
       SFLX_rain, &
       SFLX_snow  )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)    :: CCN(KA,IA,JA)
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cloud microphysics(dummy)'

    ATMOS_PHY_MP_dry_DENS(:) = UNDEF

    EVAPORATE(:,:,:) = 0.0_RP
    SFLX_rain(:,:) = 0.0_RP
    SFLX_snow(:,:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_dry

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_dry_CloudFraction( &
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
  end subroutine ATMOS_PHY_MP_dry_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_dry_EffectiveRadius( &
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

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius          [cm]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)       ! temperature               [K]
    !---------------------------------------------------------------------------

    Re = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_dry_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_dry_Mixingratio( &
       Qe,   &
       QTRC0 )
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

    Qe = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_dry_Mixingratio

end module scale_atmos_phy_mp_dry
