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
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_tracer

  use scale_const, only: &
     I_SW  => CONST_I_SW, &
     I_LW  => CONST_I_LW
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
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmOcn_driver_setup
    use mod_cpl_admin, only: &
       CPL_sw_AtmOcn,  &
       CPL_TYPE_AtmOcn
    use scale_ocean_roughness, only: &
       OCEAN_roughness_setup
    use scale_cpl_atmos_ocean, only: &
       CPL_AtmOcn_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[CPL AtmOcn] / Origin[SCALE-LES]'

    if ( CPL_sw_AtmOcn ) then

       !--- set up roughness length of sea surface
       call OCEAN_roughness_setup

       call CPL_AtmOcn_setup( CPL_TYPE_AtmOcn )

       call CPL_AtmOcn_driver( sfc_temp_update=.false. )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine CPL_AtmOcn_driver_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmOcn_driver( sfc_temp_update )
    use scale_const, only: &
       LH0  => CONST_LH0
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_ocean_roughness, only: &
       OCEAN_roughness
    use scale_cpl_atmos_ocean, only: &
       CPL_AtmOcn
    use mod_cpl_vars, only: &
       ATM_DENS   => CPL_fromAtm_ATM_DENS,   &
       ATM_U      => CPL_fromAtm_ATM_U,      &
       ATM_V      => CPL_fromAtm_ATM_V,      &
       ATM_W      => CPL_fromAtm_ATM_W,      &
       ATM_TEMP   => CPL_fromAtm_ATM_TEMP,   &
       ATM_PRES   => CPL_fromAtm_ATM_PRES,   &
       ATM_QV     => CPL_fromAtm_ATM_QV,     &
       ATM_PBL    => CPL_fromAtm_ATM_PBL,    &
       SFC_PRES   => CPL_fromAtm_SFC_PRES,   &
       FLX_precip => CPL_fromAtm_FLX_precip, &
       FLX_LW_dn  => CPL_fromAtm_FLX_LW_dn,  &
       FLX_SW_dn  => CPL_fromAtm_FLX_SW_dn,  &
       SFC_TEMP   => CPL_fromOcn_SFC_TEMP,   &
       SFC_albedo => CPL_fromOcn_SFC_albedo, &
       SFC_Z0     => CPL_fromOcn_SFC_Z0M,    &
       OCN_TEMP   => CPL_fromOcn_OCN_TEMP,   &
       CPL_AtmOcn_ATM_FLX_MW,                &
       CPL_AtmOcn_ATM_FLX_MU,                &
       CPL_AtmOcn_ATM_FLX_MV,                &
       CPL_AtmOcn_ATM_FLX_SH,                &
       CPL_AtmOcn_ATM_FLX_LH,                &
       CPL_AtmOcn_ATM_FLX_evap,              &
       CPL_AtmOcn_ATM_U10,                   &
       CPL_AtmOcn_ATM_V10,                   &
       CPL_AtmOcn_ATM_T2,                    &
       CPL_AtmOcn_ATM_Q2,                    &
       CPL_AtmOcn_OCN_FLX_heat,              &
       CPL_AtmOcn_OCN_FLX_precip,            &
       CPL_AtmOcn_OCN_FLX_evap,              &
       CNT_AtmOcn,                           &
       CNT_Ocn
    implicit none

    logical, intent(in) :: sfc_temp_update

    real(RP) :: ATM_FLX_MW  (IA,JA)
    real(RP) :: ATM_FLX_MU  (IA,JA)
    real(RP) :: ATM_FLX_MV  (IA,JA)
    real(RP) :: ATM_FLX_SH  (IA,JA)
    real(RP) :: ATM_FLX_LH  (IA,JA)
    real(RP) :: OCN_FLX_heat(IA,JA)
    real(RP) :: ATM_U10     (IA,JA)
    real(RP) :: ATM_V10     (IA,JA)
    real(RP) :: ATM_T2      (IA,JA)
    real(RP) :: ATM_Q2      (IA,JA)

    real(RP) :: SFC_Z0M     (IA,JA)
    real(RP) :: SFC_Z0H     (IA,JA)
    real(RP) :: SFC_Z0E     (IA,JA)

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Ocean'

    call OCEAN_roughness(  SFC_Z0 (:,:), & ! [INOUT]
                           SFC_Z0M(:,:), & ! [OUT]
                           SFC_Z0H(:,:), & ! [OUT]
                           SFC_Z0E(:,:), & ! [OUT]
                           ATM_U  (:,:), & ! [IN]
                           ATM_V  (:,:), & ! [IN]
                           ATM_W  (:,:)  ) ! [IN]

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, SFC_Z0M(:,:), 'SFC_Z0M' )
       call STAT_total( total, SFC_Z0H(:,:), 'SFC_Z0H' )
       call STAT_total( total, SFC_Z0E(:,:), 'SFC_Z0E' )
    endif

    call CPL_AtmOcn( SFC_TEMP    (:,:),      & ! [INOUT]
                     ATM_FLX_MU  (:,:),      & ! [OUT]
                     ATM_FLX_MV  (:,:),      & ! [OUT]
                     ATM_FLX_MW  (:,:),      & ! [OUT]
                     ATM_FLX_SH  (:,:),      & ! [OUT]
                     ATM_FLX_LH  (:,:),      & ! [OUT]
                     OCN_FLX_heat(:,:),      & ! [OUT]
                     ATM_U10     (:,:),      & ! [OUT]
                     ATM_V10     (:,:),      & ! [OUT]
                     ATM_T2      (:,:),      & ! [OUT]
                     ATM_Q2      (:,:),      & ! [OUT]
                     sfc_temp_update,        & ! [IN]
                     ATM_DENS    (:,:),      & ! [IN]
                     ATM_U       (:,:),      & ! [IN]
                     ATM_V       (:,:),      & ! [IN]
                     ATM_W       (:,:),      & ! [IN]
                     ATM_TEMP    (:,:),      & ! [IN]
                     ATM_PRES    (:,:),      & ! [IN]
                     ATM_QV      (:,:),      & ! [IN]
                     ATM_PBL     (:,:),      & ! [IN]
                     SFC_PRES    (:,:),      & ! [IN]
                     FLX_SW_dn   (:,:),      & ! [IN]
                     FLX_LW_dn   (:,:),      & ! [IN]
                     OCN_TEMP    (:,:),      & ! [IN]
                     SFC_albedo  (:,:,I_SW), & ! [IN]
                     SFC_albedo  (:,:,I_LW), & ! [IN]
                     SFC_Z0M     (:,:),      & ! [IN]
                     SFC_Z0H     (:,:),      & ! [IN]
                     SFC_Z0E     (:,:)       ) ! [IN]

    CPL_AtmOcn_ATM_FLX_MW  (:,:) = ( CPL_AtmOcn_ATM_FLX_MW  (:,:) * CNT_AtmOcn + ATM_FLX_MW(:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_FLX_MU  (:,:) = ( CPL_AtmOcn_ATM_FLX_MU  (:,:) * CNT_AtmOcn + ATM_FLX_MU(:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_FLX_MV  (:,:) = ( CPL_AtmOcn_ATM_FLX_MV  (:,:) * CNT_AtmOcn + ATM_FLX_MV(:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_FLX_SH  (:,:) = ( CPL_AtmOcn_ATM_FLX_SH  (:,:) * CNT_AtmOcn + ATM_FLX_SH(:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_FLX_LH  (:,:) = ( CPL_AtmOcn_ATM_FLX_LH  (:,:) * CNT_AtmOcn + ATM_FLX_LH(:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_FLX_evap(:,:) = ( CPL_AtmOcn_ATM_FLX_evap(:,:) * CNT_AtmOcn + ATM_FLX_LH(:,:)/LH0 ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_U10     (:,:) = ( CPL_AtmOcn_ATM_U10     (:,:) * CNT_AtmOcn + ATM_U10   (:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_V10     (:,:) = ( CPL_AtmOcn_ATM_V10     (:,:) * CNT_AtmOcn + ATM_V10   (:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_T2      (:,:) = ( CPL_AtmOcn_ATM_T2      (:,:) * CNT_AtmOcn + ATM_T2    (:,:)     ) / ( CNT_AtmOcn + 1.0_RP )
    CPL_AtmOcn_ATM_Q2      (:,:) = ( CPL_AtmOcn_ATM_Q2      (:,:) * CNT_AtmOcn + ATM_Q2    (:,:)     ) / ( CNT_AtmOcn + 1.0_RP )

    CPL_AtmOcn_OCN_FLX_heat  (:,:) = ( CPL_AtmOcn_OCN_FLX_heat  (:,:) * CNT_Ocn + OCN_FLX_heat(:,:)     ) / ( CNT_Ocn + 1.0_RP )
    CPL_AtmOcn_OCN_FLX_precip(:,:) = ( CPL_AtmOcn_OCN_FLX_precip(:,:) * CNT_Ocn + FLX_precip  (:,:)     ) / ( CNT_Ocn + 1.0_RP )
    CPL_AtmOcn_OCN_FLX_evap  (:,:) = ( CPL_AtmOcn_OCN_FLX_evap  (:,:) * CNT_Ocn - ATM_FLX_LH  (:,:)/LH0 ) / ( CNT_Ocn + 1.0_RP )

    CNT_AtmOcn = CNT_AtmOcn + 1.0_RP
    CNT_Ocn    = CNT_Ocn    + 1.0_RP

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, SFC_TEMP                 (:,:), 'OCN_SFC_TEMP  ' )
       call STAT_total( total, CPL_AtmOcn_ATM_FLX_MW    (:,:), 'ATM_FLX_MW    ' )
       call STAT_total( total, CPL_AtmOcn_ATM_FLX_MU    (:,:), 'ATM_FLX_MU    ' )
       call STAT_total( total, CPL_AtmOcn_ATM_FLX_MV    (:,:), 'ATM_FLX_MV    ' )
       call STAT_total( total, CPL_AtmOcn_ATM_FLX_SH    (:,:), 'ATM_FLX_SH    ' )
       call STAT_total( total, CPL_AtmOcn_ATM_FLX_LH    (:,:), 'ATM_FLX_LH    ' )
       call STAT_total( total, CPL_AtmOcn_ATM_FLX_evap  (:,:), 'ATM_FLX_evap  ' )
       call STAT_total( total, CPL_AtmOcn_ATM_U10       (:,:), 'ATM_U10       ' )
       call STAT_total( total, CPL_AtmOcn_ATM_V10       (:,:), 'ATM_V10       ' )
       call STAT_total( total, CPL_AtmOcn_ATM_T2        (:,:), 'ATM_T2        ' )
       call STAT_total( total, CPL_AtmOcn_ATM_Q2        (:,:), 'ATM_Q2        ' )
       call STAT_total( total, CPL_AtmOcn_OCN_FLX_heat  (:,:), 'OCN_FLX_heat  ' )
       call STAT_total( total, CPL_AtmOcn_OCN_FLX_precip(:,:), 'OCN_FLX_precip' )
       call STAT_total( total, CPL_AtmOcn_OCN_FLX_evap  (:,:), 'OCN_FLX_evap  ' )
    endif

    return
  end subroutine CPL_AtmOcn_driver

end module mod_cpl_atmos_ocean_driver