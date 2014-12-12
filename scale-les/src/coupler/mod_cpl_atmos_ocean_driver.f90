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
    use mod_admin_restart, only: &
       RESTART_RUN
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

       if( .NOT. RESTART_RUN ) then
          call CPL_AtmOcn_driver( sfc_temp_update=.false. )
       end if

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine CPL_AtmOcn_driver_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmOcn_driver( sfc_temp_update )
    use scale_const, only: &
       LHV0 => CONST_LHV0
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_ocean_roughness, only: &
       OCEAN_roughness
    use scale_cpl_atmos_ocean, only: &
       CPL_AtmOcn
    use mod_cpl_vars, only: &
       ATM_DENS       => CPL_fromAtm_ATM_DENS,      &
       ATM_U          => CPL_fromAtm_ATM_U,         &
       ATM_V          => CPL_fromAtm_ATM_V,         &
       ATM_W          => CPL_fromAtm_ATM_W,         &
       ATM_TEMP       => CPL_fromAtm_ATM_TEMP,      &
       ATM_PRES       => CPL_fromAtm_ATM_PRES,      &
       ATM_QV         => CPL_fromAtm_ATM_QV,        &
       ATM_PBL        => CPL_fromAtm_ATM_PBL,       &
       SFC_PRES       => CPL_fromAtm_SFC_PRES,      &
       FLX_precip     => CPL_fromAtm_FLX_precip,    &
       FLX_LW_dn      => CPL_fromAtm_FLX_LW_dn,     &
       FLX_SW_dn      => CPL_fromAtm_FLX_SW_dn,     &
       SFC_TEMP       => CPL_fromOcn_SFC_TEMP,      &
       SFC_albedo     => CPL_fromOcn_SFC_albedo,    &
       SFC_Z0M        => CPL_fromOcn_SFC_Z0M,       &
       SFC_Z0H        => CPL_fromOcn_SFC_Z0H,       &
       SFC_Z0E        => CPL_fromOcn_SFC_Z0E,       &
       OCN_TEMP       => CPL_fromOcn_OCN_TEMP,      &
       ATM_FLX_MW     => CPL_AtmOcn_ATM_FLX_MW,     &
       ATM_FLX_MU     => CPL_AtmOcn_ATM_FLX_MU,     &
       ATM_FLX_MV     => CPL_AtmOcn_ATM_FLX_MV,     &
       ATM_FLX_SH     => CPL_AtmOcn_ATM_FLX_SH,     &
       ATM_FLX_LH     => CPL_AtmOcn_ATM_FLX_LH,     &
       ATM_FLX_evap   => CPL_AtmOcn_ATM_FLX_evap,   &
       ATM_U10        => CPL_AtmOcn_ATM_U10,        &
       ATM_V10        => CPL_AtmOcn_ATM_V10,        &
       ATM_T2         => CPL_AtmOcn_ATM_T2,         &
       ATM_Q2         => CPL_AtmOcn_ATM_Q2,         &
       OCN_FLX_heat   => CPL_AtmOcn_OCN_FLX_heat,   &
       OCN_FLX_precip => CPL_AtmOcn_OCN_FLX_precip, &
       OCN_FLX_evap   => CPL_AtmOcn_OCN_FLX_evap,   &
       CNT_putAtm,                                  &
       CNT_putOcn,                                  &
       CNT_getAtmOcn,                               &
       CNT_getOcn
    implicit none

    logical, intent(in) :: sfc_temp_update

    real(RP) :: tmp_ATM_FLX_MW  (IA,JA)
    real(RP) :: tmp_ATM_FLX_MU  (IA,JA)
    real(RP) :: tmp_ATM_FLX_MV  (IA,JA)
    real(RP) :: tmp_ATM_FLX_SH  (IA,JA)
    real(RP) :: tmp_ATM_FLX_LH  (IA,JA)
    real(RP) :: tmp_OCN_FLX_heat(IA,JA)
    real(RP) :: tmp_ATM_U10     (IA,JA)
    real(RP) :: tmp_ATM_V10     (IA,JA)
    real(RP) :: tmp_ATM_T2      (IA,JA)
    real(RP) :: tmp_ATM_Q2      (IA,JA)

    real(RP) :: total ! dummy
    integer :: i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler: Atmos-Ocean'

    ! not used at initial setup
    if( sfc_temp_update ) then
       call OCEAN_roughness(  SFC_Z0M(:,:),   & ! [INOUT]
                              SFC_Z0H(:,:),   & ! [INOUT]
                              SFC_Z0E(:,:),   & ! [INOUT]
                              ATM_U  (:,:),   & ! [IN]
                              ATM_V  (:,:)    ) ! [IN]
    end if

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, SFC_Z0M(:,:), 'SFC_Z0M' )
       call STAT_total( total, SFC_Z0H(:,:), 'SFC_Z0H' )
       call STAT_total( total, SFC_Z0E(:,:), 'SFC_Z0E' )
    endif

    call CPL_AtmOcn( SFC_TEMP        (:,:),      & ! [INOUT]
                     tmp_ATM_FLX_MU  (:,:),      & ! [OUT]
                     tmp_ATM_FLX_MV  (:,:),      & ! [OUT]
                     tmp_ATM_FLX_MW  (:,:),      & ! [OUT]
                     tmp_ATM_FLX_SH  (:,:),      & ! [OUT]
                     tmp_ATM_FLX_LH  (:,:),      & ! [OUT]
                     tmp_OCN_FLX_heat(:,:),      & ! [OUT]
                     tmp_ATM_U10     (:,:),      & ! [OUT]
                     tmp_ATM_V10     (:,:),      & ! [OUT]
                     tmp_ATM_T2      (:,:),      & ! [OUT]
                     tmp_ATM_Q2      (:,:),      & ! [OUT]
                     sfc_temp_update,            & ! [IN]
                     ATM_DENS        (:,:),      & ! [IN]
                     ATM_U           (:,:),      & ! [IN]
                     ATM_V           (:,:),      & ! [IN]
                     ATM_W           (:,:),      & ! [IN]
                     ATM_TEMP        (:,:),      & ! [IN]
                     ATM_PRES        (:,:),      & ! [IN]
                     ATM_QV          (:,:),      & ! [IN]
                     ATM_PBL         (:,:),      & ! [IN]
                     SFC_PRES        (:,:),      & ! [IN]
                     FLX_SW_dn       (:,:),      & ! [IN]
                     FLX_LW_dn       (:,:),      & ! [IN]
                     OCN_TEMP        (:,:),      & ! [IN]
                     SFC_albedo      (:,:,I_SW), & ! [IN]
                     SFC_albedo      (:,:,I_LW), & ! [IN]
                     SFC_Z0M         (:,:),      & ! [IN]
                     SFC_Z0H         (:,:),      & ! [IN]
                     SFC_Z0E         (:,:)       ) ! [IN]

    do j = JS, JE
    do i = IS, IE
       ATM_FLX_MW  (i,j) = ( ATM_FLX_MW  (i,j) * CNT_getAtmOcn + tmp_ATM_FLX_MW(i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_FLX_MU  (i,j) = ( ATM_FLX_MU  (i,j) * CNT_getAtmOcn + tmp_ATM_FLX_MU(i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_FLX_MV  (i,j) = ( ATM_FLX_MV  (i,j) * CNT_getAtmOcn + tmp_ATM_FLX_MV(i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_FLX_SH  (i,j) = ( ATM_FLX_SH  (i,j) * CNT_getAtmOcn + tmp_ATM_FLX_SH(i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_FLX_LH  (i,j) = ( ATM_FLX_LH  (i,j) * CNT_getAtmOcn + tmp_ATM_FLX_LH(i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_FLX_evap(i,j) = ( ATM_FLX_evap(i,j) * CNT_getAtmOcn + tmp_ATM_FLX_LH(i,j)/LHV0 ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_U10     (i,j) = ( ATM_U10     (i,j) * CNT_getAtmOcn + tmp_ATM_U10   (i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_V10     (i,j) = ( ATM_V10     (i,j) * CNT_getAtmOcn + tmp_ATM_V10   (i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_T2      (i,j) = ( ATM_T2      (i,j) * CNT_getAtmOcn + tmp_ATM_T2    (i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       ATM_Q2      (i,j) = ( ATM_Q2      (i,j) * CNT_getAtmOcn + tmp_ATM_Q2    (i,j)      ) / ( CNT_getAtmOcn + 1.0_RP )
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       OCN_FLX_heat  (i,j) = ( OCN_FLX_heat  (i,j) * CNT_getOcn + tmp_OCN_FLX_heat(i,j)      ) / ( CNT_getOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       OCN_FLX_precip(i,j) = ( OCN_FLX_precip(i,j) * CNT_getOcn + FLX_precip      (i,j)      ) / ( CNT_getOcn + 1.0_RP )
    end do
    end do
    do j = JS, JE
    do i = IS, IE
       OCN_FLX_evap  (i,j) = ( OCN_FLX_evap  (i,j) * CNT_getOcn - tmp_ATM_FLX_LH  (i,j)/LHV0 ) / ( CNT_getOcn + 1.0_RP )
    end do
    end do

    CNT_putAtm    = 0.0_RP
    CNT_putOcn    = 0.0_RP
    CNT_getAtmOcn = CNT_getAtmOcn + 1.0_RP
    CNT_getOcn    = CNT_getOcn    + 1.0_RP

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, SFC_TEMP      (:,:), 'OCN_SFC_TEMP  ' )
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
       call STAT_total( total, OCN_FLX_heat  (:,:), 'OCN_FLX_heat  ' )
       call STAT_total( total, OCN_FLX_precip(:,:), 'OCN_FLX_precip' )
       call STAT_total( total, OCN_FLX_evap  (:,:), 'OCN_FLX_evap  ' )
    endif

    return
  end subroutine CPL_AtmOcn_driver

end module mod_cpl_atmos_ocean_driver
