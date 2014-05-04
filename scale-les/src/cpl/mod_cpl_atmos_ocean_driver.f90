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
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_driver_final
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

    call ATMOS_PHY_SF_driver_final
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
    use scale_grid_real, only: &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
    use scale_ocean_roughness, only: &
       OCEAN_roughness
    use scale_cpl_atmos_ocean, only: &
       CPL_AtmOcn
    use mod_cpl_vars, only: &
       SST,              &
       ALBW,             &
       Z0W,              &
       DENS => CPL_DENS, &
       MOMX => CPL_MOMX, &
       MOMY => CPL_MOMY, &
       MOMZ => CPL_MOMZ, &
       RHOS => CPL_RHOS, &
       PRES => CPL_PRES, &
       ATMP => CPL_ATMP, &
       QV   => CPL_QV  , &
       PREC => CPL_PREC, &
       SWD  => CPL_SWD , &
       LWD  => CPL_LWD , &
       TW   => CPL_TW,   &
       AtmOcn_XMFLX,     &
       AtmOcn_YMFLX,     &
       AtmOcn_ZMFLX,     &
       AtmOcn_SWUFLX,    &
       AtmOcn_LWUFLX,    &
       AtmOcn_SHFLX,     &
       AtmOcn_LHFLX,     &
       AtmOcn_QVFLX,     &
       Ocn_WHFLX,        &
       Ocn_PRECFLX,      &
       Ocn_QVFLX,        &
       CNT_Atm_Ocn,      &
       CNT_Ocn
    implicit none

    ! argument
    logical, intent(in) :: update_flag

    ! work
    real(RP) :: XMFLX (IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP) :: YMFLX (IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP) :: ZMFLX (IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP) :: SWUFLX(IA,JA) ! upward shortwave flux at the surface [W/m2]
    real(RP) :: LWUFLX(IA,JA) ! upward longwave flux at the surface [W/m2]
    real(RP) :: SHFLX (IA,JA) ! sensible heat flux at the surface [W/m2]
    real(RP) :: LHFLX (IA,JA) ! latent heat flux at the surface [W/m2]
    real(RP) :: WHFLX (IA,JA) ! water heat flux at the surface [W/m2]

    real(RP) :: DZ    (IA,JA) ! height from the surface to the lowest atmospheric layer [m]

    ! work
    integer :: i, j

    real(RP) :: tmpX(IA,JA) ! temporary XMFLX [kg/m2/s]
    real(RP) :: tmpY(IA,JA) ! temporary YMFLX [kg/m2/s]

    real(RP) :: Z0M(IA,JA) ! roughness length of momentum [m]
    real(RP) :: Z0H(IA,JA) ! roughness length of heat [m]
    real(RP) :: Z0E(IA,JA) ! roughness length of vapor [m]

    real(RP) :: Uabs       ! absolute velocity at the lowest atmospheric layer [m/s]
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Ocean'

    DZ(:,:) = CZ(KS,:,:) - FZ(KS-1,:,:)

    do j = JS-1, JE+1
    do i = IS-1, IE+1
      ! at cell center
      Uabs = sqrt( &
             ( MOMZ(i,j)               )**2 &
           + ( MOMX(i-1,j) + MOMX(i,j) )**2 &
           + ( MOMY(i,j-1) + MOMY(i,j) )**2 &
           ) / DENS(i,j) * 0.5_RP

      call OCEAN_roughness( &
        Z0W(i,j),                     & ! (inout)
        Z0M(i,j), Z0H(i,j), Z0E(i,j), & ! (out)
        Uabs, DZ(i,j)                 ) ! (in)
    end do
    end do

    call CPL_AtmOcn( &
      SST,                                 & ! (inout)
      XMFLX, YMFLX, ZMFLX,                 & ! (out)
      SWUFLX, LWUFLX, SHFLX, LHFLX, WHFLX, & ! (out)
      update_flag,                         & ! (in)
      DZ, DENS, MOMX, MOMY, MOMZ,          & ! (in)
      RHOS, PRES, ATMP, QV, SWD, LWD,      & ! (in)
      TW, ALBW(:,:,I_SW), ALBW(:,:,I_LW),  & ! (in)
      Z0M, Z0H, Z0E                        ) ! (in)

    ! interpolate momentum fluxes
    do j = JS, JE
    do i = IS, IE
      tmpX(i,j) = ( XMFLX(i,j) + XMFLX(i+1,j  ) ) * 0.5_RP ! at u/y-layer
      tmpY(i,j) = ( YMFLX(i,j) + YMFLX(i,  j+1) ) * 0.5_RP ! at x/v-layer
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
      XMFLX(i,j) = tmpX(i,j)
      YMFLX(i,j) = tmpY(i,j)
      ZMFLX(i,j) = ZMFLX(i,j) * 0.5_RP ! at w-layer
    enddo
    enddo

    ! temporal average flux
    AtmOcn_XMFLX (:,:) = ( AtmOcn_XMFLX (:,:) * CNT_Atm_Ocn + XMFLX (:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    AtmOcn_YMFLX (:,:) = ( AtmOcn_YMFLX (:,:) * CNT_Atm_Ocn + YMFLX (:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    AtmOcn_ZMFLX (:,:) = ( AtmOcn_ZMFLX (:,:) * CNT_Atm_Ocn + ZMFLX (:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    AtmOcn_SWUFLX(:,:) = ( AtmOcn_SWUFLX(:,:) * CNT_Atm_Ocn + SWUFLX(:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    AtmOcn_LWUFLX(:,:) = ( AtmOcn_LWUFLX(:,:) * CNT_Atm_Ocn + LWUFLX(:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    AtmOcn_SHFLX (:,:) = ( AtmOcn_SHFLX (:,:) * CNT_Atm_Ocn + SHFLX (:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    AtmOcn_LHFLX (:,:) = ( AtmOcn_LHFLX (:,:) * CNT_Atm_Ocn + LHFLX (:,:)     ) / ( CNT_Atm_Ocn + 1.0_RP )
    AtmOcn_QVFLX (:,:) = ( AtmOcn_QVFLX (:,:) * CNT_Atm_Ocn + LHFLX (:,:)/LH0 ) / ( CNT_Atm_Ocn + 1.0_RP )

    Ocn_WHFLX  (:,:) = ( Ocn_WHFLX  (:,:) * CNT_Ocn + WHFLX(:,:)     ) / ( CNT_Ocn + 1.0_RP )
    Ocn_PRECFLX(:,:) = ( Ocn_PRECFLX(:,:) * CNT_Ocn + PREC (:,:)     ) / ( CNT_Ocn + 1.0_RP )
    Ocn_QVFLX  (:,:) = ( Ocn_QVFLX  (:,:) * CNT_Ocn - LHFLX(:,:)/LH0 ) / ( CNT_Ocn + 1.0_RP )

    CNT_Atm_Ocn = CNT_Atm_Ocn + 1.0_RP
    CNT_Ocn     = CNT_Ocn     + 1.0_RP

    return
  end subroutine CPL_AtmOcn_driver

end module mod_cpl_atmos_ocean_driver
