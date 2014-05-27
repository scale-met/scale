!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Ocean Driver
!!
!! @par Description
!!          Coupler driver: atmosphere-ocean
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-26 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_ocean_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_AtmOcn_driver_setup
  public :: CPL_AtmOcn_driver

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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

  subroutine CPL_AtmOcn_driver_setup
    use mod_atmos_driver, only: &
       ATMOS_SURFACE_SET
    use scale_ocean_roughness, only: &
       OCEAN_roughness_setup
    use mod_ocean_phy_slab, only: &
       OCEAN_PHY_driver_final
    use mod_cpl_vars, only: &
       CPL_TYPE_AtmOcn
    use scale_cpl_atmos_ocean, only: &
       CPL_AtmOcn_setup
    implicit none
    !---------------------------------------------------------------------------

    call ATMOS_SURFACE_SET
    call OCEAN_PHY_driver_final

    !--- set up roughness length of sea surface
    call OCEAN_roughness_setup

    call CPL_AtmOcn_setup( CPL_TYPE_AtmOcn )
    call CPL_AtmOcn_driver( .false. )

    return
  end subroutine CPL_AtmOcn_driver_setup

  subroutine CPL_AtmOcn_driver( update_flag )
    use scale_const, only: &
       LH0  => CONST_LH0,  &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_ocean_roughness, only: &
       OCEAN_roughness
    use scale_cpl_atmos_ocean, only: &
       CPL_AtmOcn
    use mod_cpl_vars, only: &
       SST,               &
       ALBW,              &
       Z0W,               &
       RHOA => CPL_RHOA,  &
       UA   => CPL_UA,    &
       VA   => CPL_VA,    &
       WA   => CPL_WA,    &
       TMPA => CPL_TMPA,  &
       PRSA => CPL_PRSA,  &
       QVA  => CPL_QVA,   &
       PRSS => CPL_PRSS,  &
       PREC => CPL_PREC,  &
       SWD  => CPL_SWD,   &
       LWD  => CPL_LWD,   &
       TW   => CPL_TW,    &
       CPL_AtmOcn_XMFLX,  &
       CPL_AtmOcn_YMFLX,  &
       CPL_AtmOcn_ZMFLX,  &
       CPL_AtmOcn_SHFLX,  &
       CPL_AtmOcn_LHFLX,  &
       CPL_AtmOcn_QVFLX,  &
       Ocn_WHFLX,         &
       Ocn_PRECFLX,       &
       Ocn_QVFLX,         &
       CNT_Atm_Ocn,       &
       CNT_Ocn
    implicit none

    ! arguments
    logical, intent(in) :: update_flag

    ! works
    real(RP) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

    real(RP) :: Z0M(IA,JA) ! roughness length of momentum [m]
    real(RP) :: Z0H(IA,JA) ! roughness length of heat [m]
    real(RP) :: Z0E(IA,JA) ! roughness length of vapor [m]
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Ocean'

    call OCEAN_roughness( &
      Z0W(:,:), & ! (inout)
      Z0M(:,:), & ! (out)
      Z0H(:,:), & ! (out)
      Z0E(:,:), & ! (out)
      UA (:,:), & ! (in)
      VA (:,:), & ! (in)
      WA (:,:)  ) ! (in)

    call CPL_AtmOcn( &
      SST  (:,:),      & ! (inout)
      XMFLX(:,:),      & ! (out)
      YMFLX(:,:),      & ! (out)
      ZMFLX(:,:),      & ! (out)
      SHFLX(:,:),      & ! (out)
      LHFLX(:,:),      & ! (out)
      WHFLX(:,:),      & ! (out)
      update_flag,     & ! (in)
      RHOA (:,:),      & ! (in)
      UA   (:,:),      & ! (in)
      VA   (:,:),      & ! (in)
      WA   (:,:),      & ! (in)
      TMPA (:,:),      & ! (in)
      PRSA (:,:),      & ! (in)
      QVA  (:,:),      & ! (in)
      PRSS (:,:),      & ! (in)
      SWD  (:,:),      & ! (in)
      LWD  (:,:),      & ! (in)
      TW   (:,:),      & ! (in)
      ALBW (:,:,I_SW), & ! (in)
      ALBW (:,:,I_LW), & ! (in)
      Z0M  (:,:),      & ! (in)
      Z0H  (:,:),      & ! (in)
      Z0E  (:,:)       ) ! (in)

    ! temporal average flux
    CPL_AtmOcn_XMFLX(:,:) = ( CPL_AtmOcn_XMFLX(:,:) * CNT_Atm_Ocn + XMFLX(:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    CPL_AtmOcn_YMFLX(:,:) = ( CPL_AtmOcn_YMFLX(:,:) * CNT_Atm_Ocn + YMFLX(:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    CPL_AtmOcn_ZMFLX(:,:) = ( CPL_AtmOcn_ZMFLX(:,:) * CNT_Atm_Ocn + ZMFLX(:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    CPL_AtmOcn_SHFLX(:,:) = ( CPL_AtmOcn_SHFLX(:,:) * CNT_Atm_Ocn + SHFLX(:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    CPL_AtmOcn_LHFLX(:,:) = ( CPL_AtmOcn_LHFLX(:,:) * CNT_Atm_Ocn + LHFLX(:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    CPL_AtmOcn_QVFLX(:,:) = ( CPL_AtmOcn_QVFLX(:,:) * CNT_Atm_Ocn + LHFLX(:,:)/LH0 ) / ( CNT_Atm_Ocn + 1.0_RP )

    Ocn_WHFLX  (:,:) = ( Ocn_WHFLX  (:,:) * CNT_Ocn + WHFLX(:,:)     ) / ( CNT_Ocn + 1.0_RP )
    Ocn_PRECFLX(:,:) = ( Ocn_PRECFLX(:,:) * CNT_Ocn + PREC (:,:)     ) / ( CNT_Ocn + 1.0_RP )
    Ocn_QVFLX  (:,:) = ( Ocn_QVFLX  (:,:) * CNT_Ocn - LHFLX(:,:)/LH0 ) / ( CNT_Ocn + 1.0_RP )

    CNT_Atm_Ocn = CNT_Atm_Ocn + 1.0_RP
    CNT_Ocn     = CNT_Ocn     + 1.0_RP

    return
  end subroutine CPL_AtmOcn_driver

end module mod_cpl_atmos_ocean_driver
