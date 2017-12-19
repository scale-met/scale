!-------------------------------------------------------------------------------
!> module LAND / Physics
!!
!! @par Description
!!          land physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_land_phy_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_debug
  use scale_grid_index
  use scale_land_grid_index

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
  public :: LAND_PHY_driver_setup
  public :: LAND_PHY_driver_resume
  public :: LAND_PHY_driver

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
  !> Setup
  subroutine LAND_PHY_driver_setup
    use scale_land_phy, only: &
       LAND_PHY_setup
    use scale_land_sfc, only: &
       LAND_SFC_setup
    use mod_land_admin, only: &
       LAND_TYPE, &
       LAND_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[LAND PHY] / Origin[SCALE-RM]'

    if ( LAND_sw ) then

       ! setup library component
       call LAND_PHY_setup( LAND_TYPE )
       call LAND_SFC_setup( LAND_TYPE )

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine LAND_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine LAND_PHY_driver_resume
    use mod_admin_restart, only: &
       RESTART_RUN
    use mod_land_admin, only: &
       LAND_sw
    implicit none

    if ( LAND_sw ) then

       if ( .NOT. RESTART_RUN ) then ! tentative
          ! run once (only for the diagnostic value)
          call PROF_rapstart('LND_Physics', 1)
          call LAND_PHY_driver( update_flag = .true. )
          call PROF_rapend  ('LND_Physics', 1)
       end if
    end if

    return
  end subroutine LAND_PHY_driver_resume
  !-----------------------------------------------------------------------------
  !> Driver
  subroutine LAND_PHY_driver( update_flag )
    use scale_const, only: &
       PI => CONST_PI
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_grid_real, only: &
       REAL_Z1
    use scale_land_grid, only: &
       GRID_LCDZ
    use scale_land_phy, only: &
       LAND_PHY
    use scale_land_sfc, only: &
       LAND_SFC
    use scale_land_phy_snow_ky90, only: &
       LAND_PHY_SNOW_KY90
    use scale_land_phy_snow_diagnos, only: &
       LAND_PHY_SNOW_DIAGS
    use mod_land_admin, only: &
       LAND_TYPE
    use mod_land_vars, only: &
       LAND_PROPERTY,     &
       I_WaterLimit,      &
       I_WaterCritical,   &
       I_StomataResist,   &
       I_ThermalCond,     &
       I_HeatCapacity,    &
       I_WaterDiff,       &
       I_Z0M,             &
       I_Z0H,             &
       I_Z0E,             &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_SFC_TEMP,     &
       LAND_SFC_albedo,   &
       LAND_TEMP_t,       &
       LAND_WATER_t,      &
       LAND_SFC_TEMP_t,   &
       LAND_SFC_albedo_t, &
       LAND_SFLX_MW,      &
       LAND_SFLX_MU,      &
       LAND_SFLX_MV,      &
       LAND_SFLX_SH,      &
       LAND_SFLX_LH,      &
       LAND_SFLX_GH,      &
       LAND_SFLX_evap,    &
       LAND_U10,          &
       LAND_V10,          &
       LAND_T2,           &
       LAND_Q2,           &
       SNOW_TEMP,         &
       SNOW_SWE,          &
       SNOW_Depth,        &
       SNOW_Dzero,        &
       SNOW_nosnowsec,    &
       ATMOS_TEMP,        &
       ATMOS_PRES,        &
       ATMOS_W,           &
       ATMOS_U,           &
       ATMOS_V,           &
       ATMOS_DENS,        &
       ATMOS_QV,          &
       ATMOS_PBL,         &
       ATMOS_SFC_DENS,    &
       ATMOS_SFC_PRES,    &
       ATMOS_SFLX_LW,     &
       ATMOS_SFLX_SW,     &
       ATMOS_SFLX_prec,   & ! rain + snow
       ATMOS_SFLX_rain,   &
       ATMOS_SFLX_snow
    implicit none

    ! parameters
    real(RP), parameter :: BETA_MAX = 1.0_RP
    !real(RP) :: sw

    ! arguments
    logical, intent(in) :: update_flag

    ! works
    real(RP) :: LAND_QVEF(IA,JA)
    real(RP) :: LAND_DZ1 (IA,JA)

    real(RP) :: LHV      (IA,JA) ! latent heat of vaporization [J/kg]
    real(RP) :: total            ! dummy

    character(len=2) :: sk

    ! for snow
    real(RP) :: SNOW_TEMP_t         (IA,JA)
    real(RP) :: SNOW_albedo         (IA,JA)
    real(RP) :: SNOW_ATMO_SFLX_SH   (IA,JA)
    real(RP) :: SNOW_ATMO_SFLX_LH   (IA,JA)
    real(RP) :: SNOW_ATMO_SFLX_GH   (IA,JA)
    real(RP) :: SNOW_LAND_SFLX_GH   (IA,JA)
    real(RP) :: SNOW_LAND_SFLX_Water(IA,JA)
    real(RP) :: SNOW_frac           (IA,JA)

    real(RP) :: SNOW_ATMO_SFLX_MW   (IA,JA)
    real(RP) :: SNOW_ATMO_SFLX_MU   (IA,JA)
    real(RP) :: SNOW_ATMO_SFLX_MV   (IA,JA)
    real(RP) :: SNOW_U10            (IA,JA)
    real(RP) :: SNOW_V10            (IA,JA)
    real(RP) :: SNOW_T2             (IA,JA)
    real(RP) :: SNOW_Q2             (IA,JA)

    real(RP) :: SNOW_LAND_TEMP_t    (LKMAX,IA,JA)
    real(RP) :: SNOW_LAND_WATER_t   (LKMAX,IA,JA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

    !------------------------------------------------------------------------
    !> snow area

    if ( LAND_TYPE == 'SLAB'      .or.  &
         LAND_TYPE == 'THIN-SLAB' .or.  &
         LAND_TYPE == 'THICK-SLAB' )then

       ! accumulation and melt of snow if there is snow
       call LAND_PHY_SNOW_KY90(  SNOW_TEMP           (:,:),   & ! [INOUT]
                                 SNOW_SWE            (:,:),   & ! [INOUT]
                                 SNOW_Depth          (:,:),   & ! [INOUT]
                                 SNOW_Dzero          (:,:),   & ! [INOUT]
                                 SNOW_albedo         (:,:),   & ! [INOUT]
                                 SNOW_nosnowsec      (:,:),   & ! [INOUT]
                                 SNOW_TEMP_t         (:,:),   & ! [OUT]
                                 SNOW_ATMO_SFLX_SH   (:,:),   & ! [OUT]
                                 SNOW_ATMO_SFLX_LH   (:,:),   & ! [OUT]
                                 SNOW_ATMO_SFLX_GH   (:,:),   & ! [OUT]
                                 SNOW_LAND_SFLX_GH   (:,:),   & ! [OUT]
                                 SNOW_LAND_SFLX_Water(:,:),   & ! [OUT]
                                 SNOW_frac           (:,:),   & ! [OUT]
                                 ATMOS_SFLX_rain     (:,:),   & ! [IN]
                                 ATMOS_SFLX_snow     (:,:),   & ! [IN]
                                 ATMOS_PRES          (:,:),   & ! [IN]
                                 ATMOS_TEMP          (:,:),   & ! [IN]
                                 ATMOS_W             (:,:),   & ! [IN]
                                 ATMOS_U             (:,:),   & ! [IN]
                                 ATMOS_V             (:,:),   & ! [IN]
                                 ATMOS_QV            (:,:),   & ! [IN]
                                 ATMOS_SFLX_SW       (:,:),   & ! [IN]
                                 ATMOS_SFLX_LW       (:,:),   & ! [IN]
                                 dt                           ) ! [IN]

       ! momentum fluxes and diagnostic variables above snowpack
       call LAND_PHY_SNOW_DIAGS( SNOW_ATMO_SFLX_MW   (:,:),  & ! [OUT]
                                 SNOW_ATMO_SFLX_MU   (:,:),  & ! [OUT]
                                 SNOW_ATMO_SFLX_MV   (:,:),  & ! [OUT]
                                 SNOW_U10       (:,:),       & ! [OUT]
                                 SNOW_V10       (:,:),       & ! [OUT]
                                 SNOW_T2        (:,:),       & ! [OUT]
                                 SNOW_Q2        (:,:),       & ! [OUT]
                                 SNOW_frac      (:,:),       & ! [IN]
                                 ATMOS_TEMP     (:,:),       & ! [IN]
                                 ATMOS_PRES     (:,:),       & ! [IN]
                                 ATMOS_W        (:,:),       & ! [IN]
                                 ATMOS_U        (:,:),       & ! [IN]
                                 ATMOS_V        (:,:),       & ! [IN]
                                 ATMOS_DENS     (:,:),       & ! [IN]
                                 ATMOS_QV       (:,:),       & ! [IN]
                                 REAL_Z1        (:,:),       & ! [IN]
                                 ATMOS_PBL      (:,:),       & ! [IN]
                                 ATMOS_SFC_DENS (:,:),       & ! [IN]
                                 ATMOS_SFC_PRES (:,:),       & ! [IN]
                                 SNOW_TEMP      (:,:),       & ! [IN]
                                 LAND_QVEF      (:,:),       & ! [IN]  ! ? beta?
                                 LAND_PROPERTY  (:,:,I_Z0M), & ! [IN]
                                 LAND_PROPERTY  (:,:,I_Z0H), & ! [IN]
                                 LAND_PROPERTY  (:,:,I_Z0E)  ) ! [IN]

       ! update land temp and land water under snowpack
       do j = JS, JE
       do i = IS, IE
          ATMOS_SFLX_prec(i,j) = SNOW_LAND_SFLX_Water(i,j)
          LAND_SFLX_evap (i,j) = 0.0_RP
       enddo
       enddo

       call LAND_PHY( SNOW_LAND_TEMP_t (:,:,:),              & ! [OUT]
                      SNOW_LAND_WATER_t(:,:,:),              & ! [OUT]
                      LAND_TEMP        (:,:,:),              & ! [IN]
                      LAND_WATER       (:,:,:),              & ! [IN]
                      LAND_PROPERTY    (:,:,I_WaterLimit),   & ! [IN]
                      LAND_PROPERTY    (:,:,I_ThermalCond),  & ! [IN]
                      LAND_PROPERTY    (:,:,I_HeatCapacity), & ! [IN]
                      LAND_PROPERTY    (:,:,I_WaterDiff),    & ! [IN]
                      SNOW_LAND_SFLX_GH(:,:),                & ! [IN]
                      ATMOS_SFLX_prec  (:,:),                & ! [IN]
                      LAND_SFLX_evap   (:,:),                & ! [IN]
                      GRID_LCDZ        (:),                  & ! [IN]
                      dt                                     ) ! [IN]

    endif

    !------------------------------------------------------------------------
    !> no snow area

!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          LAND_QVEF(i,j) = min( LAND_WATER(LKS,i,j) / LAND_PROPERTY(i,j,I_WaterCritical), BETA_MAX )

          ! eq.(12) in Merlin et al.(2011) but simplified P=0.5 used
          !sw = 0.5_RP + sign(0.5_RP,LAND_WATER(LKS,i,j)-LAND_PROPERTY(i,j,I_WaterCritical)) ! if W > Wc, sw = 1
          !LAND_QVEF(i,j) = (        sw ) * 1.0_RP &
          !               + ( 1.0_RP-sw ) * sqrt( 0.5_RP - 0.5_RP * cos( PI * LAND_WATER(LKS,i,j) / LAND_PROPERTY(i,j,I_WaterCritical) ) )

          LAND_DZ1 (i,j) = GRID_LCDZ(LKS)
       end do
       end do


       call LAND_SFC( LAND_SFC_TEMP_t(:,:),                 & ! [OUT]
                      LAND_SFLX_MW   (:,:),                 & ! [OUT]
                      LAND_SFLX_MU   (:,:),                 & ! [OUT]
                      LAND_SFLX_MV   (:,:),                 & ! [OUT]
                      LAND_SFLX_SH   (:,:),                 & ! [OUT]
                      LAND_SFLX_LH   (:,:),                 & ! [OUT]
                      LAND_SFLX_GH   (:,:),                 & ! [OUT]
                      LAND_U10       (:,:),                 & ! [OUT]
                      LAND_V10       (:,:),                 & ! [OUT]
                      LAND_T2        (:,:),                 & ! [OUT]
                      LAND_Q2        (:,:),                 & ! [OUT]
                      ATMOS_TEMP     (:,:),                 & ! [IN]
                      ATMOS_PRES     (:,:),                 & ! [IN]
                      ATMOS_W        (:,:),                 & ! [IN]
                      ATMOS_U        (:,:),                 & ! [IN]
                      ATMOS_V        (:,:),                 & ! [IN]
                      ATMOS_DENS     (:,:),                 & ! [IN]
                      ATMOS_QV       (:,:),                 & ! [IN]
                      REAL_Z1        (:,:),                 & ! [IN]
                      ATMOS_PBL      (:,:),                 & ! [IN]
                      ATMOS_SFC_DENS (:,:),                 & ! [IN]
                      ATMOS_SFC_PRES (:,:),                 & ! [IN]
                      ATMOS_SFLX_LW  (:,:),                 & ! [IN]
                      ATMOS_SFLX_SW  (:,:),                 & ! [IN]
                      LAND_TEMP      (LKS,:,:),             & ! [IN]
                      LAND_SFC_TEMP  (:,:),                 & ! [IN]
                      LAND_QVEF      (:,:),                 & ! [IN]
                      LAND_SFC_albedo(:,:,I_LW),            & ! [IN]
                      LAND_SFC_albedo(:,:,I_SW),            & ! [IN]
                      LAND_DZ1       (:,:),                 & ! [IN]
                      LAND_PROPERTY  (:,:,I_StomataResist), & ! [IN]
                      LAND_PROPERTY  (:,:,I_ThermalCond),   & ! [IN]
                      LAND_PROPERTY  (:,:,I_Z0M),           & ! [IN]
                      LAND_PROPERTY  (:,:,I_Z0H),           & ! [IN]
                      LAND_PROPERTY  (:,:,I_Z0E),           & ! [IN]
                      dt                                    ) ! [IN]

       call HYDROMETEOR_LHV( LHV(:,:), ATMOS_TEMP(:,:) )

!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          LAND_SFLX_evap(i,j) = LAND_SFLX_LH(i,j) / LHV(i,j)
       end do
       end do

       call LAND_PHY( LAND_TEMP_t    (:,:,:),              & ! [OUT]
                      LAND_WATER_t   (:,:,:),              & ! [OUT]
                      LAND_TEMP      (:,:,:),              & ! [IN]
                      LAND_WATER     (:,:,:),              & ! [IN]
                      LAND_PROPERTY  (:,:,I_WaterLimit),   & ! [IN]
                      LAND_PROPERTY  (:,:,I_ThermalCond),  & ! [IN]
                      LAND_PROPERTY  (:,:,I_HeatCapacity), & ! [IN]
                      LAND_PROPERTY  (:,:,I_WaterDiff),    & ! [IN]
                      LAND_SFLX_GH   (:,:),                & ! [IN]
                      ATMOS_SFLX_prec(:,:),                & ! [IN]
                      LAND_SFLX_evap (:,:),                & ! [IN]
                      GRID_LCDZ      (:),                  & ! [IN]
                      dt                                   ) ! [IN]

       ! no albedo update (tentative)
!OCL XFILL
       LAND_SFC_albedo_t(:,:,:) = 0.0_RP


       ! marge land surface and snow surface !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( LAND_TYPE == 'SLAB'      .or.  &
         LAND_TYPE == 'THIN-SLAB' .or.  &
         LAND_TYPE == 'THICK-SLAB' )then

!OCL XFILL
       do j = JS, JE
       do i = IS, IE
           LAND_SFC_TEMP_t(i,j) = SNOW_frac(i,j)*SNOW_TEMP_t(i,j)        + (1.0_RP-SNOW_frac(i,j))*LAND_SFC_TEMP_t(i,j)
           !LAND_SFC_albedo_t(:,:,:) = 0.0_RP  ! currently not considered, but this is very important

           LAND_SFLX_MW(i,j)    = SNOW_frac(i,j)*SNOW_ATMO_SFLX_MW(i,j)  + (1.0_RP-SNOW_frac(i,j))*LAND_SFLX_MW(i,j)
           LAND_SFLX_MU(i,j)    = SNOW_frac(i,j)*SNOW_ATMO_SFLX_MU(i,j)  + (1.0_RP-SNOW_frac(i,j))*LAND_SFLX_MU(i,j)
           LAND_SFLX_MV(i,j)    = SNOW_frac(i,j)*SNOW_ATMO_SFLX_MV(i,j)  + (1.0_RP-SNOW_frac(i,j))*LAND_SFLX_MV(i,j)
           LAND_SFLX_SH(i,j)    = SNOW_frac(i,j)*SNOW_ATMO_SFLX_SH(i,j)  + (1.0_RP-SNOW_frac(i,j))*LAND_SFLX_SH(i,j)
           LAND_SFLX_LH(i,j)    = SNOW_frac(i,j)*SNOW_ATMO_SFLX_LH(i,j)  + (1.0_RP-SNOW_frac(i,j))*LAND_SFLX_LH(i,j)
           LAND_SFLX_GH(i,j)    = SNOW_frac(i,j)*SNOW_LAND_SFLX_GH(i,j)  + (1.0_RP-SNOW_frac(i,j))*LAND_SFLX_GH(i,j)
           LAND_U10(i,j)        = SNOW_frac(i,j)*SNOW_U10(i,j)           + (1.0_RP-SNOW_frac(i,j))*LAND_U10(i,j)
           LAND_V10(i,j)        = SNOW_frac(i,j)*SNOW_V10(i,j)           + (1.0_RP-SNOW_frac(i,j))*LAND_V10(i,j)
           LAND_T2 (i,j)        = SNOW_frac(i,j)*SNOW_T2 (i,j)           + (1.0_RP-SNOW_frac(i,j))*LAND_T2(i,j)
           LAND_Q2 (i,j)        = SNOW_frac(i,j)*SNOW_Q2 (i,j)           + (1.0_RP-SNOW_frac(i,j))*LAND_Q2(i,j)

         do k = LKS+1, LKE-1
           LAND_TEMP_t (k,i,j)  = SNOW_frac(i,j)*SNOW_LAND_TEMP_t (k,i,j) + (1.0_RP-SNOW_frac(i,j))*LAND_TEMP_t (k,i,j)
           LAND_WATER_t(k,i,j)  = SNOW_frac(i,j)*SNOW_LAND_WATER_t(k,i,j) + (1.0_RP-SNOW_frac(i,j))*LAND_WATER_t(k,i,j)
         enddo
       enddo
       enddo
    endif

       call FILE_HISTORY_in( LAND_TEMP_t (:,:,:), 'LAND_TEMP_t',  'tendency of LAND_TEMP',  'K',     dim_type='LXY' )
       call FILE_HISTORY_in( LAND_WATER_t(:,:,:), 'LAND_WATER_t', 'tendency of LAND_WATER', 'm3/m3', dim_type='LXY' )

       call FILE_HISTORY_in( LAND_SFC_TEMP_t  (:,:),      'LAND_SFC_TEMP_t', 'tendency of LAND_SFC_TEMP', 'K', dim_type='XY' )
       call FILE_HISTORY_in( LAND_SFC_albedo_t(:,:,I_LW), 'LAND_ALB_LW_t',   'tendency of LAND_ALB_LW',   '1', dim_type='XY' )
       call FILE_HISTORY_in( LAND_SFC_albedo_t(:,:,I_SW), 'LAND_ALB_SW_t',   'tendency of LAND_ALB_SW',   '1', dim_type='XY' )

       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_WaterLimit),    'LP_WaterLimit',    'LAND PROPERTY, WaterLimit',    'm3/m3',  dim_type='XY' )
       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_WaterCritical), 'LP_WaterCritical', 'LAND PROPERTY, WaterCritical', 'm3/m3',  dim_type='XY' )
       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_ThermalCond),   'LP_ThermalCond',   'LAND PROPERTY, ThermalCond',   'W/K/m',  dim_type='XY' )
       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_HeatCapacity),  'LP_HeatCapacity',  'LAND PROPERTY, HeatCapacity',  'J/K/m3', dim_type='XY' )
       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_WaterDiff),     'LP_WaterDiff',     'LAND PROPERTY, WaterDiff',     'm2/s',   dim_type='XY' )
       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_Z0M),           'LP_Z0M',           'LAND PROPERTY, Z0M',           'm',      dim_type='XY' )
       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_Z0H),           'LP_Z0H',           'LAND PROPERTY, Z0H',           'm',      dim_type='XY' )
       call FILE_HISTORY_in( LAND_PROPERTY(:,:,I_Z0E),           'LP_Z0E',           'LAND PROPERTY, Z0E',           'm',      dim_type='XY' )

    endif

    if ( STATISTICS_checktotal ) then
       do k = LKS, LKE
          write(sk,'(I2.2)') k

          call STAT_total( total, LAND_TEMP_t (k,:,:), 'LAND_TEMP_t'//sk  )
          call STAT_total( total, LAND_WATER_t(k,:,:), 'LAND_WATER_t'//sk )
       enddo

       call STAT_total( total, LAND_SFC_TEMP_t  (:,:),      'LAND_SFC_TEMP_t'  )
       call STAT_total( total, LAND_SFC_albedo_t(:,:,I_LW), 'LAND_ALB_LW_t'    )
       call STAT_total( total, LAND_SFC_albedo_t(:,:,I_SW), 'LAND_ALB_SW_t'    )
    endif

    return
  end subroutine LAND_PHY_driver

end module mod_land_phy_driver
