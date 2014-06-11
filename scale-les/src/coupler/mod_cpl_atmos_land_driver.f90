!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Land Driver
!!
!! @par Description
!!          Coupler driver: atmosphere-land
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-02-25 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_land_driver
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
  public :: CPL_AtmLnd_driver_setup
  public :: CPL_AtmLnd_driver

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
  subroutine CPL_AtmLnd_driver_setup
    use mod_cpl_admin, only: &
       CPL_sw_AtmLnd,   &
       CPL_TYPE_AtmLnd
    use scale_cpl_atmos_land, only: &
       CPL_AtmLnd_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[CPL AtmLnd] / Origin[SCALE-LES]'

    if ( CPL_sw_AtmLnd ) then

       call CPL_AtmLnd_setup( CPL_TYPE_AtmLnd )

       call CPL_AtmLnd_driver( sfc_temp_update=.false. )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine CPL_AtmLnd_driver_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmLnd_driver( sfc_temp_update )
    use scale_const, only: &
       LH0  => CONST_LH0,  &
       I_SW => CONST_I_SW, &
       I_LW => CONST_I_LW
    use scale_cpl_atmos_land, only: &
       CPL_AtmLnd
    use mod_cpl_vars, only: &
       ATM_DENS   => CPL_fromAtm_ATM_DENS,   &
       ATM_U      => CPL_fromAtm_ATM_U,      &
       ATM_V      => CPL_fromAtm_ATM_V,      &
       ATM_W      => CPL_fromAtm_ATM_W,      &
       ATM_TEMP   => CPL_fromAtm_ATM_TEMP,   &
       ATM_PRES   => CPL_fromAtm_ATM_PRES,   &
       ATM_QV     => CPL_fromAtm_ATM_QV,     &
       SFC_PRES   => CPL_fromAtm_SFC_PRES,   &
       FLX_precip => CPL_fromAtm_FLX_precip, &
       FLX_LW_dn  => CPL_fromAtm_FLX_LW_dn,  &
       FLX_SW_dn  => CPL_fromAtm_FLX_SW_dn,  &
       SFC_TEMP   => CPL_fromLnd_SFC_TEMP,   &
       SFC_albedo => CPL_fromLnd_SFC_albedo, &
       LND_TCS    => CPL_fromLnd_LND_TCS,    &
       LND_DZ     => CPL_fromLnd_LND_DZ,     &
       SFC_Z0M    => CPL_fromLnd_SFC_Z0M,    &
       SFC_Z0H    => CPL_fromLnd_SFC_Z0H,    &
       SFC_Z0E    => CPL_fromLnd_SFC_Z0E,    &
       LND_TEMP   => CPL_fromLnd_LND_TEMP,   &
       LND_BETA   => CPL_fromLnd_LND_BETA,   &
       CPL_AtmLnd_ATM_FLX_MW,                &
       CPL_AtmLnd_ATM_FLX_MU,                &
       CPL_AtmLnd_ATM_FLX_MV,                &
       CPL_AtmLnd_ATM_FLX_SH,                &
       CPL_AtmLnd_ATM_FLX_LH,                &
       CPL_AtmLnd_ATM_FLX_evap,              &
       CPL_AtmLnd_ATM_U10,                   &
       CPL_AtmLnd_ATM_V10,                   &
       CPL_AtmLnd_ATM_T2,                    &
       CPL_AtmLnd_ATM_Q2,                    &
       CPL_AtmLnd_LND_FLX_heat,              &
       CPL_AtmLnd_LND_FLX_precip,            &
       CPL_AtmLnd_LND_FLX_evap,              &
       CNT_AtmLnd,                           &
       CNT_Lnd
    implicit none

    logical, intent(in) :: sfc_temp_update

    real(RP) :: ATM_FLX_MW  (IA,JA)
    real(RP) :: ATM_FLX_MU  (IA,JA)
    real(RP) :: ATM_FLX_MV  (IA,JA)
    real(RP) :: ATM_FLX_SH  (IA,JA)
    real(RP) :: ATM_FLX_LH  (IA,JA)
    real(RP) :: LND_FLX_heat(IA,JA)
    real(RP) :: ATM_U10     (IA,JA)
    real(RP) :: ATM_V10     (IA,JA)
    real(RP) :: ATM_T2      (IA,JA)
    real(RP) :: ATM_Q2      (IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Land'

    call CPL_AtmLnd( SFC_TEMP    (:,:),      & ! [INOUT]
                     ATM_FLX_MU  (:,:),      & ! [OUT]
                     ATM_FLX_MV  (:,:),      & ! [OUT]
                     ATM_FLX_MW  (:,:),      & ! [OUT]
                     ATM_FLX_SH  (:,:),      & ! [OUT]
                     ATM_FLX_LH  (:,:),      & ! [OUT]
                     LND_FLX_heat(:,:),      & ! [OUT]
                     ATM_U10     (:,:),      & ! [OUT]
                     ATM_V10     (:,:),      & ! [OUT]
                     ATM_T2      (:,:),      & ! [OUT]
                     ATM_Q2      (:,:),      & ! [OUT]
                     sfc_temp_update,            & ! [IN]
                     ATM_DENS    (:,:),      & ! [IN]
                     ATM_U       (:,:),      & ! [IN]
                     ATM_V       (:,:),      & ! [IN]
                     ATM_W       (:,:),      & ! [IN]
                     ATM_TEMP    (:,:),      & ! [IN]
                     ATM_PRES    (:,:),      & ! [IN]
                     ATM_QV      (:,:),      & ! [IN]
                     SFC_PRES    (:,:),      & ! [IN]
                     FLX_SW_dn   (:,:),      & ! [IN]
                     FLX_LW_dn   (:,:),      & ! [IN]
                     LND_TEMP    (:,:),      & ! [IN]
                     LND_BETA    (:,:),      & ! [IN]
                     SFC_albedo  (:,:,I_SW), & ! [IN]
                     SFC_albedo  (:,:,I_LW), & ! [IN]
                     LND_TCS     (:,:),      & ! [IN]
                     LND_DZ      (:,:),      & ! [IN]
                     SFC_Z0M     (:,:),      & ! [IN]
                     SFC_Z0H     (:,:),      & ! [IN]
                     SFC_Z0E     (:,:)       ) ! [IN]

    ! temporal average flux
    CPL_AtmLnd_ATM_FLX_MW  (:,:) = ( CPL_AtmLnd_ATM_FLX_MW  (:,:) * CNT_AtmLnd + ATM_FLX_MW(:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_FLX_MU  (:,:) = ( CPL_AtmLnd_ATM_FLX_MU  (:,:) * CNT_AtmLnd + ATM_FLX_MU(:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_FLX_MV  (:,:) = ( CPL_AtmLnd_ATM_FLX_MV  (:,:) * CNT_AtmLnd + ATM_FLX_MV(:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_FLX_SH  (:,:) = ( CPL_AtmLnd_ATM_FLX_SH  (:,:) * CNT_AtmLnd + ATM_FLX_SH(:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_FLX_LH  (:,:) = ( CPL_AtmLnd_ATM_FLX_LH  (:,:) * CNT_AtmLnd + ATM_FLX_LH(:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_FLX_evap(:,:) = ( CPL_AtmLnd_ATM_FLX_evap(:,:) * CNT_AtmLnd + ATM_FLX_LH(:,:)/LH0 ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_U10     (:,:) = ( CPL_AtmLnd_ATM_U10     (:,:) * CNT_AtmLnd + ATM_U10   (:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_V10     (:,:) = ( CPL_AtmLnd_ATM_V10     (:,:) * CNT_AtmLnd + ATM_V10   (:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_T2      (:,:) = ( CPL_AtmLnd_ATM_T2      (:,:) * CNT_AtmLnd + ATM_T2    (:,:)     ) / ( CNT_AtmLnd + 1.0_RP )
    CPL_AtmLnd_ATM_Q2      (:,:) = ( CPL_AtmLnd_ATM_Q2      (:,:) * CNT_AtmLnd + ATM_Q2    (:,:)     ) / ( CNT_AtmLnd + 1.0_RP )

    CPL_AtmLnd_LND_FLX_heat  (:,:) = ( CPL_AtmLnd_LND_FLX_heat  (:,:) * CNT_Lnd + LND_FLX_heat(:,:)     ) / ( CNT_Lnd + 1.0_RP )
    CPL_AtmLnd_LND_FLX_precip(:,:) = ( CPL_AtmLnd_LND_FLX_precip(:,:) * CNT_Lnd + FLX_precip  (:,:)     ) / ( CNT_Lnd + 1.0_RP )
    CPL_AtmLnd_LND_FLX_evap  (:,:) = ( CPL_AtmLnd_LND_FLX_evap  (:,:) * CNT_Lnd - ATM_FLX_LH  (:,:)/LH0 ) / ( CNT_Lnd + 1.0_RP )

    CNT_AtmLnd = CNT_AtmLnd + 1.0_RP
    CNT_Lnd    = CNT_Lnd    + 1.0_RP

    return
  end subroutine CPL_AtmLnd_driver

end module mod_cpl_atmos_land_driver
