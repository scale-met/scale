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
  !> Setup
  subroutine ATMOS_PHY_MP_dry_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: MP_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Cloud Microphysics] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** dummy process (dry Atmosphere)'

    if ( MP_TYPE /= 'DRY' ) then
       write(*,*) 'xxx ATMOS_PHY_MP_TYPE is not DRY. Check!'
       call PRC_MPIstop
    endif

    ATMOS_PHY_MP_DENS(:) = 0.0_RP

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
       SFLX_rain, &
       SFLX_snow  )
    use scale_tracer, only: &
       QAD => QA
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cloud microphysics(dummy)'

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

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius          [cm]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)        ! density                   [kg/m3]

    real(RP), parameter :: um2cm = 100.0_RP
    !---------------------------------------------------------------------------

    Re(:,:,:,:) = 8.E-6_RP * um2cm ! dummy

    return
  end subroutine ATMOS_PHY_MP_dry_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_dry_Mixingratio( &
       Qe,    &
       QTRC0  )
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
