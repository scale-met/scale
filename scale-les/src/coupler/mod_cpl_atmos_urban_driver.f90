!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Urban Driver
!!
!! @par Description
!!          Coupler driver: atmosphere-urban
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_urban_driver
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
  public :: CPL_AtmUrb_driver_setup
  public :: CPL_AtmUrb_driver

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
  logical, allocatable, private :: is_FLX(:,:) ! is urban coupler run?

contains
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmUrb_driver_setup
    use mod_cpl_admin, only: &
       CPL_sw_AtmUrb,  &
       CPL_TYPE_AtmUrb
    use scale_cpl_atmos_urban, only: &
       CPL_AtmUrb_setup
   use scale_landuse, only: &
       LANDUSE_fact_urban
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[CPL AtmUrb] / Origin[SCALE-LES]'

    if ( CPL_sw_AtmUrb ) then

       !--- judge to run atmos-urban coupler
       allocate( is_FLX(IA,JA) )

       do j = 1, JA
       do i = 1, IA
         if( LANDUSE_fact_urban(i,j) > 0.0_RP ) then
           is_FLX(i,j) = .true.
         else
           is_FLX(i,j) = .false.
         end if
       end do
       end do

       call CPL_AtmUrb_setup( CPL_TYPE_AtmUrb )

       call CPL_AtmUrb_driver( sfc_temp_update=.false. )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine CPL_AtmUrb_driver_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmUrb_driver( sfc_temp_update )
    use scale_const, only: &
       LHV0 => CONST_LHV0
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_cpl_atmos_urban, only: &
       CPL_AtmUrb
    use scale_cpl_atmos_urban_bulk, only: &
       CPL_AtmUrb_bulk_restart,    &
       CPL_AtmUrb_bulk_momentum
    use mod_cpl_vars, only: &
       ATM_DENS       => CPL_fromAtm_ATM_DENS,      &
       ATM_U          => CPL_fromAtm_ATM_U,         &
       ATM_V          => CPL_fromAtm_ATM_V,         &
       ATM_W          => CPL_fromAtm_ATM_W,         &
       ATM_TEMP       => CPL_fromAtm_ATM_TEMP,      &
       ATM_QV         => CPL_fromAtm_ATM_QV,        &
       ATM_PRES       => CPL_fromAtm_ATM_PRES,      &
       ATM_PBL        => CPL_fromAtm_ATM_PBL,       &
       SFC_PRES       => CPL_fromAtm_SFC_PRES,      &
       FLX_precip     => CPL_fromAtm_FLX_precip,    &
       FLX_LW_dn      => CPL_fromAtm_FLX_LW_dn,     &
       FLX_SW_dn      => CPL_fromAtm_FLX_SW_dn,     &
       SFC_TEMP       => CPL_fromUrb_SFC_TEMP,      &
       SFC_albedo     => CPL_fromUrb_SFC_albedo,    &
       SFC_Z0         => CPL_fromUrb_SFC_Z0M,       &
       ATM_FLX_MW     => CPL_AtmUrb_ATM_FLX_MW,     &
       ATM_FLX_MU     => CPL_AtmUrb_ATM_FLX_MU,     &
       ATM_FLX_MV     => CPL_AtmUrb_ATM_FLX_MV,     &
       ATM_FLX_SH     => CPL_AtmUrb_ATM_FLX_SH,     &
       ATM_FLX_LH     => CPL_AtmUrb_ATM_FLX_LH,     &
       ATM_FLX_evap   => CPL_AtmUrb_ATM_FLX_evap,   &
       ATM_U10        => CPL_AtmUrb_ATM_U10,        &
       ATM_V10        => CPL_AtmUrb_ATM_V10,        &
       ATM_T2         => CPL_AtmUrb_ATM_T2,         &
       ATM_Q2         => CPL_AtmUrb_ATM_Q2,         &
       URB_FLX_heat   => CPL_AtmUrb_URB_FLX_heat,   &
       URB_FLX_precip => CPL_AtmUrb_URB_FLX_precip, &
       URB_FLX_evap   => CPL_AtmUrb_URB_FLX_evap,   &
       CNT_getAtmUrb,                               &
       CNT_getUrb
    use scale_grid_real, only: &
       Z1  => REAL_Z1
    use scale_mapproj,   only: &
       CLON => MPRJ_basepoint_lon, &
       CLAT => MPRJ_basepoint_lat
    use mod_urban_vars, only: &
       TR_URB,    &
       TG_URB,    &
       TB_URB,    &
       TC_URB,    &
       QC_URB,    &
       UC_URB,    &
       TS_URB,    &
       TRL_URB,   &
       TGL_URB,   &
       TBL_URB,   &
       SHR_URB,   &
       SHB_URB,   &
       SHG_URB,   &
       LHR_URB,   &
       LHB_URB,   &
       LHG_URB,   &
       GHR_URB,   &
       GHB_URB,   &
       GHG_URB,   &
       RnR_URB,   &
       RnB_URB,   &
       RnG_URB,   &
       RAINR_URB, &
       RAINB_URB, &
       RAING_URB, &
       ROFF_URB,  &
       AH_URB,    &
       ALH_URB,   &
       Rngrd_URB, &
       SHFLX_URB, &
       LHFLX_URB, &
       GHFLX_URB

    implicit none

    logical, intent(in) :: sfc_temp_update

    real(RP) :: tmp_ATM_FLX_MW  (IA,JA)
    real(RP) :: tmp_ATM_FLX_MU  (IA,JA)
    real(RP) :: tmp_ATM_FLX_MV  (IA,JA)
    real(RP) :: tmp_ATM_FLX_SH  (IA,JA)
    real(RP) :: tmp_ATM_FLX_LH  (IA,JA)
    real(RP) :: tmp_URB_FLX_heat(IA,JA)
    real(RP) :: tmp_ATM_U10     (IA,JA)
    real(RP) :: tmp_ATM_V10     (IA,JA)
    real(RP) :: tmp_ATM_T2      (IA,JA)
    real(RP) :: tmp_ATM_Q2      (IA,JA)

    logical  :: LSOLAR = .false.    ! logical [true=both, false=SSG only]
    real(RP) :: Uabs                ! wind speed at the lowest atmospheric level    [m/s]
    integer  :: i, j

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Urban'

    if ( sfc_temp_update ) then  ! temperature update

      do j = 1, JA
      do i = 1, IA

      if( is_FLX(i,j) ) then

        Uabs = max( sqrt( ATM_U(i,j)**2 + ATM_V(i,j)**2 + ATM_W(i,j)**2 ), 0.1_RP)

        call CPL_AtmUrb( TR_URB          (i,j),      & ! [INOUT]
                         TB_URB          (i,j),      & ! [INOUT]
                         TG_URB          (i,j),      & ! [INOUT]
                         TC_URB          (i,j),      & ! [INOUT]
                         QC_URB          (i,j),      & ! [INOUT]
                         UC_URB          (i,j),      & ! [INOUT]
                         TRL_URB         (:,i,j),    & ! [INOUT]
                         TBL_URB         (:,i,j),    & ! [INOUT]
                         TGL_URB         (:,i,j),    & ! [INOUT]
                         RAINR_URB       (i,j),      & ! [INOUT]
                         RAINB_URB       (i,j),      & ! [INOUT]
                         RAING_URB       (i,j),      & ! [INOUT]
                         ROFF_URB        (i,j),      & ! [INOUT]
                         AH_URB          (i,j),      & ! [INOUT]
                         ALH_URB         (i,j),      & ! [INOUT]
                         SFC_albedo      (i,j,I_LW), & ! [INOUT]
                         SFC_albedo      (i,j,I_SW), & ! [INOUT]
                         SHR_URB         (i,j),      & ! [OUT]
                         SHB_URB         (i,j),      & ! [OUT]
                         SHG_URB         (i,j),      & ! [OUT]
                         LHR_URB         (i,j),      & ! [INOUT]
                         LHB_URB         (i,j),      & ! [INOUT]
                         LHG_URB         (i,j),      & ! [INOUT]
                         GHR_URB         (i,j),      & ! [OUT]
                         GHB_URB         (i,j),      & ! [OUT]
                         GHG_URB         (i,j),      & ! [OUT]
                         RnR_URB         (i,j),      & ! [OUT]
                         RnB_URB         (i,j),      & ! [OUT]
                         RnG_URB         (i,j),      & ! [OUT]
                         SFC_TEMP        (i,j),      & ! [OUT]
                         Rngrd_URB       (i,j),      & ! [OUT]
                         tmp_ATM_FLX_SH  (i,j),      & ! [OUT]
                         tmp_ATM_FLX_LH  (i,j),      & ! [OUT]
                         tmp_URB_FLX_heat(i,j),      & ! [OUT]
                         tmp_ATM_U10     (i,j),      & ! [OUT]
                         tmp_ATM_V10     (i,j),      & ! [OUT]
                         tmp_ATM_T2      (i,j),      & ! [OUT]
                         tmp_ATM_Q2      (i,j),      & ! [OUT]
                         LSOLAR,                     & ! [IN]
                         SFC_PRES        (i,j),      & ! [IN]
                         ATM_TEMP        (i,j),      & ! [IN]
                         ATM_QV          (i,j),      & ! [IN]
                         Uabs,                       & ! [IN]
                         ATM_U           (i,j),      & ! [IN]
                         ATM_V           (i,j),      & ! [IN]
                         Z1              (i,j),      & ! [IN]
                         FLX_SW_dn       (i,j),      & ! [IN]
                         FLX_LW_dn       (i,j),      & ! [IN]
                         FLX_precip      (i,j),      & ! [IN]
                         ATM_DENS        (i,j),      & ! [IN]
                         CLON,                       & ! [IN]
                         CLAT                        ) ! [IN]

                         TS_URB      (i,j) = SFC_TEMP        (i,j)
                         SHFLX_URB   (i,j) = tmp_ATM_FLX_SH  (i,j)
                         LHFLX_URB   (i,j) = tmp_ATM_FLX_LH  (i,j)
                         GHFLX_URB   (i,j) = tmp_URB_FLX_heat(i,j)

        call CPL_AtmUrb_bulk_momentum( &
                         tmp_ATM_FLX_MU  (i,j), & ! (out)
                         tmp_ATM_FLX_MV  (i,j), & ! (out)
                         tmp_ATM_FLX_MW  (i,j), & ! (out)
                         SFC_Z0          (i,j), & ! (out)
                         Z1              (i,j), & ! (in)
                         ATM_DENS        (i,j), & ! (in)
                         ATM_U           (i,j), & ! (in)
                         ATM_V           (i,j), & ! (in)
                         ATM_W           (i,j), & ! (in)
                         ATM_TEMP        (i,j), & ! (in)
                         ATM_PRES        (i,j), & ! (in)
                         ATM_QV          (i,j), & ! (in)
                         ATM_PBL         (i,j), & ! (in)
                         SFC_PRES        (i,j), & ! (in)
                         SFC_TEMP        (i,j)  ) ! (in)

      else
             ! not calculate urban scheme
                         SHR_URB         (i,j) = 0.0_RP
                         SHB_URB         (i,j) = 0.0_RP
                         SHG_URB         (i,j) = 0.0_RP
                         LHR_URB         (i,j) = 0.0_RP
                         LHB_URB         (i,j) = 0.0_RP
                         LHG_URB         (i,j) = 0.0_RP
                         GHR_URB         (i,j) = 0.0_RP
                         GHB_URB         (i,j) = 0.0_RP
                         GHG_URB         (i,j) = 0.0_RP
                         RnR_URB         (i,j) = 0.0_RP
                         RnB_URB         (i,j) = 0.0_RP
                         RnG_URB         (i,j) = 0.0_RP
                         AH_URB          (i,j) = 0.0_RP
                         ALH_URB         (i,j) = 0.0_RP
                         Rngrd_URB       (i,j) = 0.0_RP
                         SHFLX_URB       (i,j) = 0.0_RP
                         LHFLX_URB       (i,j) = 0.0_RP
                         GHFLX_URB       (i,j) = 0.0_RP
                         tmp_ATM_FLX_SH  (i,j) = 0.0_RP
                         tmp_ATM_FLX_LH  (i,j) = 0.0_RP
                         tmp_URB_FLX_heat(i,j) = 0.0_RP
                         tmp_ATM_U10     (i,j) = 0.0_RP
                         tmp_ATM_V10     (i,j) = 0.0_RP
                         tmp_ATM_T2      (i,j) = 0.0_RP
                         tmp_ATM_Q2      (i,j) = 0.0_RP
                         SFC_TEMP        (i,j) = 300.0_RP
                         TS_URB          (i,j) = SFC_TEMP(i,j)
                         tmp_ATM_FLX_MU  (i,j) = 0.0_RP
                         tmp_ATM_FLX_MV  (i,j) = 0.0_RP
                         tmp_ATM_FLX_MW  (i,j) = 0.0_RP
      endif
      enddo
      enddo

    else     ! temperature not update

      do j = 1, JA
      do i = 1, IA

      if( is_FLX(i,j) ) then

        Uabs = max( sqrt( ATM_U(i,j)**2 + ATM_V(i,j)**2 + ATM_W(i,j)**2 ), 0.1_RP)

        call CPL_AtmUrb_bulk_restart(  &
                         TR_URB          (i,j),      & ! [IN]
                         TB_URB          (i,j),      & ! [IN]
                         TG_URB          (i,j),      & ! [IN]
                         TC_URB          (i,j),      & ! [IN]
                         QC_URB          (i,j),      & ! [IN]
                         UC_URB          (i,j),      & ! [IN]
                         TRL_URB         (:,i,j),    & ! [IN]
                         TBL_URB         (:,i,j),    & ! [IN]
                         TGL_URB         (:,i,j),    & ! [IN]
                         RAINR_URB       (i,j),      & ! [IN]
                         RAINB_URB       (i,j),      & ! [IN]
                         RAING_URB       (i,j),      & ! [IN]
                         AH_URB          (i,j),      & ! [IN]
                         ALH_URB         (i,j),      & ! [IN]
                         SFC_albedo      (i,j,I_LW), & ! [INOUT]
                         SFC_albedo      (i,j,I_SW), & ! [INOUT]
                         SHR_URB         (i,j),      & ! [OUT]
                         SHB_URB         (i,j),      & ! [OUT]
                         SHG_URB         (i,j),      & ! [OUT]
                         LHR_URB         (i,j),      & ! [OUT]
                         LHB_URB         (i,j),      & ! [OUT]
                         LHG_URB         (i,j),      & ! [OUT]
                         GHR_URB         (i,j),      & ! [OUT]
                         GHB_URB         (i,j),      & ! [OUT]
                         GHG_URB         (i,j),      & ! [OUT]
                         RnR_URB         (i,j),      & ! [OUT]
                         RnB_URB         (i,j),      & ! [OUT]
                         RnG_URB         (i,j),      & ! [OUT]
                         SFC_TEMP        (i,j),      & ! [OUT]
                         Rngrd_URB       (i,j),      & ! [OUT]
                         tmp_ATM_FLX_SH  (i,j),      & ! [OUT]
                         tmp_ATM_FLX_LH  (i,j),      & ! [OUT]
                         tmp_URB_FLX_heat(i,j),      & ! [OUT]
                         tmp_ATM_U10     (i,j),      & ! [OUT]
                         tmp_ATM_V10     (i,j),      & ! [OUT]
                         tmp_ATM_T2      (i,j),      & ! [OUT]
                         tmp_ATM_Q2      (i,j),      & ! [OUT]
                         LSOLAR,                     & ! [IN]
                         SFC_PRES        (i,j),      & ! [IN]
                         ATM_TEMP        (i,j),      & ! [IN]
                         ATM_QV          (i,j),      & ! [IN]
                         Uabs,                       & ! [IN]
                         ATM_U           (i,j),      & ! [IN]
                         ATM_V           (i,j),      & ! [IN]
                         Z1              (i,j),      & ! [IN]
                         FLX_SW_dn       (i,j),      & ! [IN]
                         FLX_LW_dn       (i,j),      & ! [IN]
                         ATM_DENS        (i,j),      & ! [IN]
                         CLON,                       & ! [IN]
                         CLAT                        ) ! [IN]

                         TS_URB      (i,j) = SFC_TEMP (i,j)
                         SHFLX_URB   (i,j) = ATM_FLX_SH  (i,j)
                         LHFLX_URB   (i,j) = ATM_FLX_LH  (i,j)
                         GHFLX_URB   (i,j) = URB_FLX_heat(i,j)

        call CPL_AtmUrb_bulk_momentum( &
                         tmp_ATM_FLX_MU  (i,j), & ! (out)
                         tmp_ATM_FLX_MV  (i,j), & ! (out)
                         tmp_ATM_FLX_MW  (i,j), & ! (out)
                         SFC_Z0          (i,j), & ! (out)
                         Z1              (i,j), & ! (in)
                         ATM_DENS        (i,j), & ! (in)
                         ATM_U           (i,j), & ! (in)
                         ATM_V           (i,j), & ! (in)
                         ATM_W           (i,j), & ! (in)
                         ATM_TEMP        (i,j), & ! (in)
                         ATM_PRES        (i,j), & ! (in)
                         ATM_QV          (i,j), & ! (in)
                         ATM_PBL         (i,j), & ! (in)
                         SFC_PRES        (i,j), & ! (in)
                         SFC_TEMP        (i,j)  ) ! (in)
      else
             ! not calculate surface flux
                         SHR_URB         (i,j) = 0.0_RP
                         SHB_URB         (i,j) = 0.0_RP
                         SHG_URB         (i,j) = 0.0_RP
                         LHR_URB         (i,j) = 0.0_RP
                         LHB_URB         (i,j) = 0.0_RP
                         LHG_URB         (i,j) = 0.0_RP
                         GHR_URB         (i,j) = 0.0_RP
                         GHB_URB         (i,j) = 0.0_RP
                         GHG_URB         (i,j) = 0.0_RP
                         RnR_URB         (i,j) = 0.0_RP
                         RnB_URB         (i,j) = 0.0_RP
                         RnG_URB         (i,j) = 0.0_RP
                         AH_URB          (i,j) = 0.0_RP
                         ALH_URB         (i,j) = 0.0_RP
                         Rngrd_URB       (i,j) = 0.0_RP
                         SHFLX_URB       (i,j) = 0.0_RP
                         LHFLX_URB       (i,j) = 0.0_RP
                         GHFLX_URB       (i,j) = 0.0_RP
                         tmp_ATM_FLX_SH  (i,j) = 0.0_RP
                         tmp_ATM_FLX_LH  (i,j) = 0.0_RP
                         tmp_URB_FLX_heat(i,j) = 0.0_RP
                         tmp_ATM_U10     (i,j) = 0.0_RP
                         tmp_ATM_V10     (i,j) = 0.0_RP
                         tmp_ATM_T2      (i,j) = 0.0_RP
                         tmp_ATM_Q2      (i,j) = 0.0_RP
                         SFC_TEMP        (i,j) = 300.0_RP
                         TS_URB          (i,j) = SFC_TEMP(i,j)
                         tmp_ATM_FLX_MU  (i,j) = 0.0_RP
                         tmp_ATM_FLX_MV  (i,j) = 0.0_RP
                         tmp_ATM_FLX_MW  (i,j) = 0.0_RP
      endif
      enddo
      enddo

    endif

    ! temporal average flux
    ATM_FLX_MW  (:,:) = ( ATM_FLX_MW  (:,:) * CNT_getAtmUrb + tmp_ATM_FLX_MW(:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_FLX_MU  (:,:) = ( ATM_FLX_MU  (:,:) * CNT_getAtmUrb + tmp_ATM_FLX_MU(:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_FLX_MV  (:,:) = ( ATM_FLX_MV  (:,:) * CNT_getAtmUrb + tmp_ATM_FLX_MV(:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_FLX_SH  (:,:) = ( ATM_FLX_SH  (:,:) * CNT_getAtmUrb + tmp_ATM_FLX_SH(:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_FLX_LH  (:,:) = ( ATM_FLX_LH  (:,:) * CNT_getAtmUrb + tmp_ATM_FLX_LH(:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_FLX_evap(:,:) = ( ATM_FLX_evap(:,:) * CNT_getAtmUrb + tmp_ATM_FLX_LH(:,:)/LHV0 ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_U10     (:,:) = ( ATM_U10     (:,:) * CNT_getAtmUrb + tmp_ATM_U10   (:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_V10     (:,:) = ( ATM_V10     (:,:) * CNT_getAtmUrb + tmp_ATM_V10   (:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_T2      (:,:) = ( ATM_T2      (:,:) * CNT_getAtmUrb + tmp_ATM_T2    (:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )
    ATM_Q2      (:,:) = ( ATM_Q2      (:,:) * CNT_getAtmUrb + tmp_ATM_Q2    (:,:)      ) / ( CNT_getAtmUrb + 1.0_RP )

    URB_FLX_heat  (:,:) = ( URB_FLX_heat  (:,:) * CNT_getUrb + tmp_URB_FLX_heat(:,:)      ) / ( CNT_getUrb + 1.0_RP )
    URB_FLX_precip(:,:) = ( URB_FLX_precip(:,:) * CNT_getUrb + FLX_precip      (:,:)      ) / ( CNT_getUrb + 1.0_RP )
    URB_FLX_evap  (:,:) = ( URB_FLX_evap  (:,:) * CNT_getUrb - tmp_ATM_FLX_LH  (:,:)/LHV0 ) / ( CNT_getUrb + 1.0_RP )

    CNT_getAtmUrb = CNT_getAtmUrb + 1.0_RP
    CNT_getUrb    = CNT_getUrb    + 1.0_RP

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, SFC_TEMP      (:,:), 'URB_SFC_TEMP  ' )
       call STAT_total( total, ATM_FLX_MW    (:,:), 'ATM_FLX_MW    ' )
       call STAT_total( total, ATM_FLX_MU    (:,:), 'ATM_FLX_MU    ' )
       call STAT_total( total, ATM_FLX_MV    (:,:), 'ATM_FLX_MV    ' )
       call STAT_total( total, ATM_FLX_SH    (:,:), 'ATM_FLX_SH    ' )
       call STAT_total( total, ATM_FLX_LH    (:,:), 'ATM_FLX_LH    ' )
       call STAT_total( total, ATM_FLX_evap  (:,:), 'ATM_FLX_evap  ' )
       call STAT_total( total, ATM_U10       (:,:), 'ATM_U10       ' )
       call STAT_total( total, ATM_V10       (:,:), 'ATM_V10       ' )
       call STAT_total( total, ATM_T2        (:,:), 'ATM_T2        ' )
       call STAT_total( total, ATM_Q2        (:,:), 'ATM_Q2        ' )
       call STAT_total( total, URB_FLX_heat  (:,:), 'URB_FLX_heat  ' )
       call STAT_total( total, URB_FLX_precip(:,:), 'URB_FLX_precip' )
       call STAT_total( total, URB_FLX_evap  (:,:), 'URB_FLX_evap  ' )
    endif

    return
  end subroutine CPL_AtmUrb_driver

end module mod_cpl_atmos_urban_driver
