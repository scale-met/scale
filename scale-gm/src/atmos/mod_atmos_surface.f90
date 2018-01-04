!-------------------------------------------------------------------------------
!> module ATMOSPHERE driver
!!
!! @par Description
!!          Atmosphere module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2017-12-15 (K.Kikuchi) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_surface
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
!  use scale_prof
!  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_SURFACE_GET
  public :: ATMOS_SURFACE_SET

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
  !-----------------------------------------------------------------------------
  !> Get surface boundary condition
  subroutine ATMOS_SURFACE_GET
    use mod_atmos_phy_sf_vars, only: &
       SFC_TEMP   => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo => ATMOS_PHY_SF_SFC_albedo, &
       SFC_Z0M    => ATMOS_PHY_SF_SFC_Z0M,    &
       SFC_Z0H    => ATMOS_PHY_SF_SFC_Z0H,    &
       SFC_Z0E    => ATMOS_PHY_SF_SFC_Z0E,    &
       SFLX_MW    => ATMOS_PHY_SF_SFLX_MW,    &
       SFLX_MU    => ATMOS_PHY_SF_SFLX_MU,    &
       SFLX_MV    => ATMOS_PHY_SF_SFLX_MV,    &
       SFLX_SH    => ATMOS_PHY_SF_SFLX_SH,    &
       SFLX_LH    => ATMOS_PHY_SF_SFLX_LH,    &
       SFLX_GH    => ATMOS_PHY_SF_SFLX_GH,    &
       SFLX_QTRC  => ATMOS_PHY_SF_SFLX_QTRC,  &
       U10        => ATMOS_PHY_SF_U10,        &
       V10        => ATMOS_PHY_SF_V10,        &
       T2         => ATMOS_PHY_SF_T2,         &
       Q2         => ATMOS_PHY_SF_Q2
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_getSFC_ATM
    implicit none
    !---------------------------------------------------------------------------

    if ( CPL_sw ) then
       call CPL_getSFC_ATM( SFC_TEMP  (:,:),   & ! [OUT]
                            SFC_albedo(:,:,:), & ! [OUT]
                            SFC_Z0M   (:,:),   & ! [OUT]
                            SFC_Z0H   (:,:),   & ! [OUT]
                            SFC_Z0E   (:,:),   & ! [OUT]
                            SFLX_MW   (:,:),   & ! [OUT]
                            SFLX_MU   (:,:),   & ! [OUT]
                            SFLX_MV   (:,:),   & ! [OUT]
                            SFLX_SH   (:,:),   & ! [OUT]
                            SFLX_LH   (:,:),   & ! [OUT]
                            SFLX_GH   (:,:),   & ! [OUT]
                            SFLX_QTRC (:,:,:), & ! [OUT]
                            U10       (:,:),   & ! [OUT]
                            V10       (:,:),   & ! [OUT]
                            T2        (:,:),   & ! [OUT]
                            Q2        (:,:)    ) ! [OUT]
    endif

    return
  end subroutine ATMOS_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Set surface boundary condition
  subroutine ATMOS_SURFACE_SET( countup )
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_Z1
    use scale_topography, only: &
       TOPO_Zsfc
    use scale_atmos_bottom, only: &
       BOTTOM_estimate => ATMOS_BOTTOM_estimate
    use mod_atmos_vars, only: &
       DENS, &
       QTRC, &
       TEMP, &
       PRES, &
       W,    &
       U,    &
       V
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_rd_vars, only: &
       SFLX_rad_dn => ATMOS_PHY_RD_SFLX_downall, &
       cosSZA      => ATMOS_PHY_RD_cosSZA
    use mod_cpl_admin, only: &
       CPL_sw
    use mod_cpl_vars, only: &
       CPL_putATM
    use scale_grid_index, only: &
       IA,    &
       JA,    &
       KS
    implicit none

    ! arguments
    logical, intent(in) :: countup

    ! works
    real(RP) :: SFC_DENS(IA,JA)
    real(RP) :: SFC_PRES(IA,JA)
    real(RP) :: ATM_PBL (IA,JA)
    !---------------------------------------------------------------------------

    if ( CPL_sw ) then
       ! planetary boundary layer
       ATM_PBL(:,:) = 100.0_RP ! tentative

       call BOTTOM_estimate( DENS     (:,:,:), & ! [IN]
                             PRES     (:,:,:), & ! [IN]
                             REAL_CZ  (:,:,:), & ! [IN]
                             TOPO_Zsfc(:,:),   & ! [IN]
                             REAL_Z1  (:,:),   & ! [IN]
                             SFC_DENS (:,:),   & ! [OUT]
                             SFC_PRES (:,:)    ) ! [OUT]

       call CPL_putATM( TEMP       (KS,:,:),   & ! [IN]
                        PRES       (KS,:,:),   & ! [IN]
                        W          (KS,:,:),   & ! [IN]
                        U          (KS,:,:),   & ! [IN]
                        V          (KS,:,:),   & ! [IN]
                        DENS       (KS,:,:),   & ! [IN]
                        QTRC       (KS,:,:,:), & ! [IN]
                        ATM_PBL    (:,:),      & ! [IN]
                        SFC_PRES   (:,:),      & ! [IN]
                        SFLX_rad_dn(:,:,:,:),  & ! [IN]
                        cosSZA     (:,:),      & ! [IN]
                        SFLX_rain  (:,:),      & ! [IN]
                        SFLX_snow  (:,:),      & ! [IN]
                        countup                ) ! [IN]
    endif

    return
  end subroutine ATMOS_SURFACE_SET

end module mod_atmos_surface
