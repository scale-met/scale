!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land model driver
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_land_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_land_grid_cartesC_index

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
  public :: LAND_driver_setup
  public :: LAND_driver_calc_tendency
  public :: LAND_driver_update
  public :: LAND_SURFACE_GET
  public :: LAND_SURFACE_SET

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
  logical, private :: snow_flag
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_driver_setup
    use scale_prc, only: &
       PRC_abort
    use mod_land_admin, only: &
       LAND_do, &
       LAND_DYN_TYPE, &
       LAND_SFC_TYPE, &
       SNOW_TYPE
    use scale_land_dyn_bucket, only: &
       LAND_DYN_BUCKET_setup
    use scale_land_phy_snow_ky90, only: &
       LAND_PHY_SNOW_KY90_setup
    use scale_cpl_phy_sfc_skin, only: &
       CPL_PHY_SFC_skin_setup
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER Setup] / Categ[LAND] / Origin[SCALE-RM]'

    snow_flag = .false.

    if ( LAND_do ) then

       select case ( LAND_DYN_TYPE )
       case ( 'BUCKET' )
          call LAND_DYN_BUCKET_setup
       case ( 'CONST' )
          ! do nothing
       case default
          write(*,*) 'xxx LAND_DYN_TYPE is invalid: ', trim(LAND_DYN_TYPE)
          call PRC_abort
       end select

       select case ( LAND_SFC_TYPE )
       case ( 'SKIN' )
          call CPL_PHY_SFC_skin_setup
       case ( 'FIXED-TEMP' )
          call CPL_PHY_SFC_fixed_temp_setup
       case default
          write(*,*) 'xxx LAND_SFC_TYPE is invalid: ', trim(LAND_SFC_TYPE)
          call PRC_abort
       end select

       select case ( SNOW_TYPE )
       case ( 'NONE', 'OFF' )
       case ( 'KY90' )
          if ( IO_L ) write(IO_FID_LOG,*) '*** SNOW model is enabled'
          if ( IO_L ) write(IO_FID_LOG,*) '*** SNOW model is on experimental stage.'
          if ( IO_L ) write(IO_FID_LOG,*) '*** Use this with your own risk.'
          call LAND_PHY_SNOW_KY90_setup
          snow_flag = .true.
       case default
          write(*,*) 'xxx SNOW_TYPE is invalid: ', trim(SNOW_TYPE)
          call PRC_abort
       end select

    end if

    return
  end subroutine LAND_driver_setup

  !-----------------------------------------------------------------------------
  !> Calculate tendency
  subroutine LAND_driver_calc_tendency( force )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_land_grid_cartesC, only: &
       LCDZ => LAND_GRID_CARTESC_CDZ
    use scale_land_phy_snow_ky90, only: &
       LAND_PHY_SNOW_KY90
    use scale_land_phy_snow_diagnos, only: &
       LAND_PHY_SNOW_DIAGS
    use scale_cpl_phy_sfc_skin, only: &
       CPL_PHY_SFC_skin
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp
    use mod_land_admin, only: &
       LAND_SFC_TYPE, &
       SNOW_TYPE
    use mod_land_vars, only: &
       I_WaterLimit,      &
       I_WaterCritical,   &
       I_StomataResist,   &
       I_ThermalCond,     &
       I_HeatCapacity,    &
       I_WaterDiff,       &
       I_Z0M,             &
       I_Z0H,             &
       I_Z0E,             &
       LAND_PROPERTY,     &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_SFC_TEMP,     &
       LAND_SFC_albedo,   &
       LAND_type_albedo,  &
       LAND_TEMP_t,     &
       LAND_WATER_t,    &
       LAND_SFLX_MW,      &
       LAND_SFLX_MU,      &
       LAND_SFLX_MV,      &
       LAND_SFLX_SH,      &
       LAND_SFLX_LH,      &
       LAND_SFLX_GH,      &
       LAND_SFLX_evap,    &
       LAND_SFLX_water,   &
       LAND_SFLX_ice,     &
       LAND_U10,          &
       LAND_V10,          &
       LAND_T2,           &
       LAND_Q2,           &
       SNOW_SFC_TEMP,     &
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
       ATMOS_SFLX_rain,   &
       ATMOS_SFLX_snow
    use scale_landuse, only: &
       LANDUSE_fact_land
    implicit none
    logical, intent(in) :: force

    ! parameters
    real(RP), parameter :: BETA_MAX = 1.0_RP

    ! works
    real(RP) :: SNOW_QVEF(LIA,LJA)
    real(RP) :: LAND_QVEF(LIA,LJA)
    real(RP) :: LAND_DZ1 (LIA,LJA)
    real(RP) :: SFLX_GH  (LIA,LJA)
    real(RP) :: LHV      (LIA,LJA) ! latent heat of vaporization [J/kg]

    ! for snow
    real(RP) :: SNOW_albedo         (LIA,LJA,2)
    real(RP) :: SNOW_albedo_t       (LIA,LJA,2)
    real(RP) :: SNOW_ATMOS_SFLX_SH  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_LH  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_GH  (LIA,LJA)
    !real(RP) :: SNOW_ATMOS_SFLX_evap(LIA,LJA)
    real(RP) :: SNOW_LAND_SFLX_GH   (LIA,LJA)
    real(RP) :: SNOW_LAND_SFLX_water(LIA,LJA)
    real(RP) :: SNOW_LAND_SFLX_ice  (LIA,LJA)
    real(RP) :: SNOW_frac           (LIA,LJA)

    real(RP) :: SNOW_ATMOS_SFLX_MW  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_MU  (LIA,LJA)
    real(RP) :: SNOW_ATMOS_SFLX_MV  (LIA,LJA)
    real(RP) :: SNOW_U10            (LIA,LJA)
    real(RP) :: SNOW_V10            (LIA,LJA)
    real(RP) :: SNOW_T2             (LIA,LJA)
    real(RP) :: SNOW_Q2             (LIA,LJA)

    ! monitor
    !real(RP) :: MONIT_WCONT0        (LIA,LJA)
    !real(RP) :: MONIT_WCONT1        (LIA,LJA)
    !real(RP) :: MONIT_ENG0          (LIA,LJA)
    !real(RP) :: MONIT_ENG1          (LIA,LJA)
    !
    !real(RP) :: MONIT_SNOW_heat     (LIA,LJA)
    !real(RP) :: MONIT_SNOW_water    (LIA,LJA)
    !real(RP) :: MONIT_LAND_heat     (LIA,LJA)
    !real(RP) :: MONIT_LAND_water    (LIA,LJA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_CalcTend', 1)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER CALC TEND] / Categ[LAND] / Origin[SCALE-RM]'

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET

    !########## reset tendencies ##########
!OCL XFILL
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_TEMP_t (k,i,j) = 0.0_RP
       LAND_WATER_t(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    call HYDROMETEOR_LHV( LIA, LIS, LIE, LJA, LJS, LJE, &
                          ATMOS_TEMP(:,:), LHV(:,:) )

    if ( snow_flag ) then
       !------------------------------------------------------------------------
       !> snow area

!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          ! This is for debug---adachi start
          !if(( int(SNOW_frac(i,j)) == 1 ).and.( abs(SNOW_SFC_TEMP(i,j)-LAND_SFC_TEMP(i,j))/=0 ))then
          !   write(*,*) "xxx Error please check SNOW_SFC_TEMP routine"
          !endif
          ! This is for debug---adachi end
          SNOW_SFC_TEMP(i,j) = LAND_SFC_TEMP(i,j)
       end do
       end do

       select case ( SNOW_TYPE )
       case ( 'KY90' )
          ! accumulation and melt of snow if there is snow

          !MONIT_WCONT0 = 0.0_RP
          !call monitor_snow_water(  SNOW_Depth          (:,:),   & ! [IN]
          !                          SNOW_Dzero          (:,:),   & ! [IN]
          !                          MONIT_WCONT0        (:,:)    ) ! [OUT]

          call LAND_PHY_SNOW_KY90( LIA, LIS, LIE, LJA, LJS, LJE, &
                                   ATMOS_SFLX_rain(:,:), ATMOS_SFLX_snow(:,:),       & ! [IN]
                                   ATMOS_PRES(:,:), ATMOS_TEMP(:,:), ATMOS_QV(:,:),  & ! [IN]
                                   ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),         & ! [IN]
                                   ATMOS_SFC_DENS(:,:),                              & ! [IN]
                                   ATMOS_SFLX_SW(:,:), ATMOS_SFLX_LW(:,:),           & ! [IN]
                                   LANDUSE_fact_land(:,:), dt,                       & ! [IN]
                                   SNOW_SFC_TEMP(:,:), SNOW_SWE(:,:),                & ! [INOUT]
                                   SNOW_Depth(:,:), SNOW_Dzero(:,:),                 & ! [INOUT]
                                   SNOW_nosnowsec(:,:),                              & ! [INOUT]
                                   SNOW_albedo(:,:,:),                               & ! [OUT]
                                   SNOW_ATMOS_SFLX_SH(:,:), SNOW_ATMOS_SFLX_LH(:,:), & ! [OUT]
                                   SNOW_ATMOS_SFLX_GH(:,:), SNOW_LAND_SFLX_GH(:,:),  & ! [OUT]
                                   SNOW_LAND_SFLX_water(:,:),                        & ! [OUT]
                                   SNOW_frac           (:,:)                         ) ! [OUT]

!OCL XFILL
          !$omp parallel do
          do j = LJS, LJE
          do i = LIS, LIE
             SNOW_LAND_SFLX_ice(i,j) = 0.0_RP
             SNOW_albedo_t(i,j,I_SW) = ( SNOW_albedo(i,j,I_SW) - LAND_SFC_albedo(i,j,I_SW) ) / dt
             SNOW_albedo_t(i,j,I_LW) = ( SNOW_albedo(i,j,I_LW) - LAND_SFC_albedo(i,j,I_LW) ) / dt
          enddo
          enddo

       end select

!OCL XFILL
       !!$omp parallel do
       !do j = LJS, LJE
       !do i = LIS, LIE
       !   SNOW_ATMOS_SFLX_evap (i,j) = - SNOW_ATMOS_SFLX_LH(i,j) / LHV(i,j)
       !end do
       !end do
       !call monitor_snow_water(  SNOW_Depth          (:,:),   & ! [IN]
       !                          SNOW_Dzero          (:,:),   & ! [IN]
       !                          MONIT_WCONT1        (:,:)    ) ! [OUT]

       !call monitor_land_regidual( ATMOS_SFLX_rain     (:,:),   & ! [IN] ! downward at surface
       !                            ATMOS_SFLX_snow     (:,:),   & ! [IN] ! downward at surface
       !                            SNOW_ATMOS_SFLX_evap(:,:),   & ! [IN] ! upward   at surface
       !                            SNOW_LAND_SFLX_water(:,:),   & ! [IN] ! downward at bottom
       !                            SNOW_LAND_SFLX_ice  (:,:),   & ! [IN] ! downward at bottom
       !                            MONIT_WCONT0        (:,:),   & ! [IN]
       !                            MONIT_WCONT1        (:,:),   & ! [IN]
       !                            MONIT_SNOW_water    (:,:)    ) ! [OUT]

!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          SNOW_QVEF(i,j) = 1.0_RP ! tentative
       end do
       end do

       ! momentum fluxes and diagnostic variables above snowpack
       call LAND_PHY_SNOW_DIAGS( LIA, LIS, LIE, LJA, LJS, LJE, &
                                 SNOW_frac(:,:),                                                & ! [IN]
                                 ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                              & ! [IN]
                                 ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                      & ! [IN]
                                 ATMOS_DENS(:,:), ATMOS_QV(:,:),                                & ! [IN]
                                 REAL_Z1(:,:), ATMOS_PBL(:,:),                                  & ! [IN]
                                 ATMOS_SFC_DENS (:,:), ATMOS_SFC_PRES(:,:), SNOW_SFC_TEMP(:,:), & ! [IN]
                                 SNOW_QVEF(:,:),                                                & ! [IN]
                                 LAND_PROPERTY(:,:,I_Z0M),                                      & ! [IN]
                                 LAND_PROPERTY(:,:,I_Z0H),                                      & ! [IN]
                                 LAND_PROPERTY(:,:,I_Z0E),                                      & ! [IN]
                                 SNOW_ATMOS_SFLX_MW(:,:),                                       & ! [OUT]
                                 SNOW_ATMOS_SFLX_MU(:,:),                                       & ! [OUT]
                                 SNOW_ATMOS_SFLX_MV(:,:),                                       & ! [OUT]
                                 SNOW_U10(:,:), SNOW_V10(:,:),                                  & ! [OUT]
                                 SNOW_T2(:,:), SNOW_Q2(:,:)                                     ) ! [OUT]


       call FILE_HISTORY_in( SNOW_albedo         (:,:,I_SW), 'SNOW_ALB_SW',          'Snow surface albedo (short wave)',         '0-1',    dim_type='XY' )
       call FILE_HISTORY_in( SNOW_albedo         (:,:,I_LW), 'SNOW_ALB_LW',          'Snow surface albedo (long wave)',          '0-1',    dim_type='XY' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_SH  (:,:),      'SNOW_ATMOS_SFLX_SH',   'Snow surface sensible heat flux',          'J/m2/s', dim_type='XY' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_LH  (:,:),      'SNOW_ATMOS_SFLX_LH',   'Snow surface latent heat flux',            'J/m2/s', dim_type='XY' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_GH  (:,:),      'SNOW_ATMOS_SFLX_GH',   'Snowpack received heat flux',              'J/m2/s', dim_type='XY' )
       call FILE_HISTORY_in( SNOW_LAND_SFLX_GH   (:,:),      'SNOW_LAND_SFLX_GH',    'land surface ground heat flux under snow', 'J/m2/s', dim_type='XY' )
       call FILE_HISTORY_in( SNOW_LAND_SFLX_water(:,:),      'SNOW_LAND_SFLX_water', 'land surface liquid water flux under snow', 'kg/m2/s',dim_type='XY' )
       call FILE_HISTORY_in( SNOW_LAND_SFLX_water(:,:),      'SNOW_LAND_SFLX_ice',   'land surface ice water flux under snow',    'kg/m2/s',dim_type='XY' )
       call FILE_HISTORY_in( SNOW_frac           (:,:),      'SNOW_frac',            'Snow fraction on land subgrid',            '0-1',    dim_type='XY' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_MW  (:,:),      'SNOW_ATMOS_SFLX_MW',   'Snow surface w-momentum flux',             'J/m2/s', dim_type='XY' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_MU  (:,:),      'SNOW_ATMOS_SFLX_MU',   'Snow surface u-momentum flux',             'J/m2/s', dim_type='XY' )
       call FILE_HISTORY_in( SNOW_ATMOS_SFLX_MV  (:,:),      'SNOW_ATMOS_SFLX_MV',   'Snow surface v-momentum flux',             'J/m2/s', dim_type='XY' )
       call FILE_HISTORY_in( SNOW_U10            (:,:),      'SNOW_U10',             'Wind velocity u at 10 m on snow surface',  'm/s',    dim_type='XY' )
       call FILE_HISTORY_in( SNOW_V10            (:,:),      'SNOW_V10',             'Wind velocity v at 10 m on snow surface',  'm/s',    dim_type='XY' )
       call FILE_HISTORY_in( SNOW_T2             (:,:),      'SNOW_T2',              'Air temperature at 2m on snow surface',    'K',      dim_type='XY' )
       call FILE_HISTORY_in( SNOW_Q2             (:,:),      'SNOW_Q2',              'Specific humidity at 2m on snow surface',  'kg/kg',  dim_type='XY' )

    endif


!OCL XFILL
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       LAND_QVEF(i,j) = min( LAND_WATER(LKS,i,j) / LAND_PROPERTY(i,j,I_WaterCritical), BETA_MAX )

       ! eq.(12) in Merlin et al.(2011) but simplified P=0.5 used
       !sw = 0.5_RP + sign(0.5_RP,LAND_WATER(LKS,i,j)-LAND_PROPERTY(i,j,I_WaterCritical)) ! if W > Wc, sw = 1
       !LAND_QVEF(i,j) = (        sw ) * 1.0_RP &
       !               + ( 1.0_RP-sw ) * sqrt( 0.5_RP - 0.5_RP * cos( PI * LAND_WATER(LKS,i,j) / LAND_PROPERTY(i,j,I_WaterCritical) ) )

       LAND_DZ1(i,j) = LCDZ(LKS)
    end do
    end do


    !------------------------------------------------------------------------
    !> all land area without snow model or no snow area with snow model


    select case ( LAND_SFC_TYPE )
    case ( 'SKIN' )

       call CPL_PHY_SFC_skin( LIA, LIS, LIE, LJA, LJS, LJE, &
                              ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                        & ! [IN]
                              ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                & ! [IN]
                              ATMOS_DENS(:,:), ATMOS_QV(:,:), LHV(:,:),                & ! [IN]
                              REAL_Z1(:,:), ATMOS_PBL(:,:),                            & ! [IN]
                              ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                & ! [IN]
                              ATMOS_SFLX_LW(:,:), ATMOS_SFLX_SW(:,:),                  & ! [IN]
                              LAND_TEMP(LKS,:,:), LAND_QVEF(:,:),                      & ! [IN]
                              LAND_type_albedo(:,:,I_LW), LAND_type_albedo(:,:,I_SW),  & ! [IN]
                              LAND_DZ1(:,:),                                           & ! [IN]
                              LAND_PROPERTY(:,:,I_StomataResist),                      & ! [IN]
                              LAND_PROPERTY(:,:,I_ThermalCond),                        & ! [IN]
                              LAND_PROPERTY(:,:,I_Z0M),                                & ! [IN]
                              LAND_PROPERTY(:,:,I_Z0H),                                & ! [IN]
                              LAND_PROPERTY(:,:,I_Z0E),                                & ! [IN]
                              LANDUSE_fact_land, dt,                                   & ! [IN]
                              'LAND',                                                  & ! [IN]
                              LAND_SFC_TEMP(:,:),                                      & ! [INOUT]
                              LAND_SFLX_MW(:,:), LAND_SFLX_MU(:,:), LAND_SFLX_MV(:,:), & ! [OUT]
                              LAND_SFLX_SH(:,:), LAND_SFLX_LH(:,:), SFLX_GH(:,:),      & ! [OUT]
                              LAND_U10(:,:), LAND_V10(:,:), LAND_T2(:,:), LAND_Q2(:,:) ) ! [OUT]

!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFC_albedo(i,j,I_LW) = LAND_type_albedo(i,j,I_LW)
          LAND_SFC_albedo(i,j,I_SW) = LAND_type_albedo(i,j,I_SW)
       end do
       end do

    case ( 'FIXED-TEMP' )

       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFC_TEMP(i,j) = LAND_TEMP(LKS,i,j)
       end do
       end do

       call CPL_PHY_SFC_fixed_temp( LIA, LIS, LIE, LJA, LJS, LJE, &
                                    ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                        & ! [IN]
                                    ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                & ! [IN]
                                    ATMOS_DENS(:,:), ATMOS_QV(:,:), LHV(:,:),                & ! [IN]
                                    REAL_Z1(:,:), ATMOS_PBL(:,:),                            & ! [IN]
                                    ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                & ! [IN]
                                    ATMOS_SFLX_LW(:,:), ATMOS_SFLX_SW(:,:),                  & ! [IN]
                                    LAND_SFC_TEMP(:,:), LAND_QVEF(:,:),                      & ! [IN]
                                    LAND_type_albedo(:,:,I_LW), LAND_type_albedo(:,:,I_SW),  & ! [IN]
                                    LAND_PROPERTY(:,:,I_StomataResist),                      & ! [IN]
                                    LAND_PROPERTY(:,:,I_Z0M),                                & ! [IN]
                                    LAND_PROPERTY(:,:,I_Z0H),                                & ! [IN]
                                    LAND_PROPERTY(:,:,I_Z0E),                                & ! [IN]
                                    LANDUSE_fact_land, dt,                                   & ! [IN]
                                    LAND_SFLX_MW(:,:), LAND_SFLX_MU(:,:), LAND_SFLX_MV(:,:), & ! [OUT]
                                    LAND_SFLX_SH(:,:), LAND_SFLX_LH(:,:), SFLX_GH(:,:),      & ! [OUT]
                                    LAND_U10(:,:), LAND_V10(:,:),                            & ! [OUT]
                                    LAND_T2(:,:), LAND_Q2(:,:)                               ) ! [OUT]
    end select

    ! LAND_SFLX_* are positive for downward
!OCL XFILL
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       LAND_SFLX_GH   (i,j) = - SFLX_GH(i,j) ! inverse sign ( positive for upward to downward )
       LAND_SFLX_evap (i,j) = LAND_SFLX_LH(i,j) / LHV(i,j)
       LAND_SFLX_water(i,j) = ATMOS_SFLX_rain(i,j) - LAND_SFLX_evap(i,j)
       LAND_SFLX_ice  (i,j) = ATMOS_SFLX_snow(i,j)
    end do
    end do


    if ( snow_flag ) then

       ! marge land surface and snow surface !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OCL XFILL
       !$omp parallel do
       do j = LJS, LJE
       do i = LIS, LIE
          LAND_SFC_TEMP(i,j) =          SNOW_frac(i,j)   * SNOW_SFC_TEMP(i,j) &
                             + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_SFC_TEMP(i,j)
          LAND_SFC_albedo(i,j,I_LW) =   SNOW_frac(i,j)   * SNOW_albedo    (i,j,I_LW) &
                             + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_SFC_albedo(i,j,I_LW)
          LAND_SFC_albedo(i,j,I_SW) =   SNOW_frac(i,j)   * SNOW_albedo    (i,j,I_SW) &
                             + ( 1.0_RP-SNOW_frac(i,j) ) * LAND_SFC_albedo(i,j,I_SW)

          ! flux to the soil
          LAND_SFLX_GH   (i,j) =          SNOW_frac(i,j)   * SNOW_LAND_SFLX_GH   (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_SFLX_GH   (i,j)
          LAND_SFLX_water(i,j) =          SNOW_frac(i,j)   * SNOW_LAND_SFLX_water(i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_SFLX_water(i,j)
          LAND_SFLX_ice  (i,j) =          SNOW_frac(i,j)   * SNOW_LAND_SFLX_ice  (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_SFLX_ice  (i,j)
          ! flux to the atmosphere
          LAND_SFLX_SH   (i,j) =          SNOW_frac(i,j)   * SNOW_ATMOS_SFLX_SH  (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *       LAND_SFLX_SH  (i,j)
          LAND_SFLX_LH   (i,j) =          SNOW_frac(i,j)   * SNOW_ATMOS_SFLX_LH  (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *       LAND_SFLX_LH  (i,j)
          LAND_SFLX_evap (i,j) =  LAND_SFLX_LH(i,j) / LHV(i,j)
          ! diagnostics
          LAND_U10       (i,j) =          SNOW_frac(i,j)   *      SNOW_U10       (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_U10       (i,j)
          LAND_V10       (i,j) =          SNOW_frac(i,j)   *      SNOW_V10       (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_V10       (i,j)
          LAND_T2        (i,j) =          SNOW_frac(i,j)   *      SNOW_T2        (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_T2        (i,j)
          LAND_Q2        (i,j) =          SNOW_frac(i,j)   *      SNOW_Q2        (i,j) &
                               + ( 1.0_RP-SNOW_frac(i,j) ) *      LAND_Q2        (i,j)
       enddo
       enddo
             
    end if

    !########## Set Surface Boundary to coupler ##########
    call LAND_SURFACE_SET( countup=.true. )

    call PROF_rapend  ('LND_CalcTend', 1)

    return
  end subroutine LAND_driver_calc_tendency

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_driver_update
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use mod_land_vars, only: &
       LAND_PROPERTY,     &
       I_WaterLimit,      &
       I_ThermalCond,     &
       I_HeatCapacity,    &
       I_WaterDiff,       &
       LAND_TEMP,         &
       LAND_WATER,        &
       LAND_SFLX_GH,      &
       LAND_SFLX_water,   &
       LAND_SFLX_ice,     &
       LAND_SFC_TEMP,     &
       LAND_SFC_albedo,   &
       LAND_TEMP_t,       &
       LAND_WATER_t,      &
       LAND_vars_total,   &
       LAND_vars_history
    use scale_land_grid_cartesC, only: &
       LCDZ => LAND_GRID_CARTESC_CDZ
    use scale_land_dyn_bucket, only: &
       LAND_DYN_bucket
    use scale_landuse, only: &
       LANDUSE_fact_land
    use scale_time, only: &
       NOWDAYSEC => TIME_NOWDAYSEC
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_land_admin, only: &
       LAND_DYN_TYPE
    implicit none

    real(RP) :: RUNOFF(LIA,LJA)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_Update', 2)

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER UPDATE] / Categ[LAND] / Origin[SCALE-RM]'

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET

    !########## Dynamics / Update variables ##########
    select case ( LAND_DYN_TYPE )
    case ( 'CONST' )
       ! do nothing
    case ( 'BUCKET' )
       call LAND_DYN_bucket( LKMAX, LKS, LKE, LIA, LIS, LIE, LJA, LJS, LJE, &
                             LAND_TEMP_t(:,:,:), LAND_WATER_t(:,:,:),   & ! [IN]
                             LAND_PROPERTY(:,:,I_WaterLimit),           & ! [IN]
                             LAND_PROPERTY(:,:,I_ThermalCond),          & ! [IN]
                             LAND_PROPERTY(:,:,I_HeatCapacity),         & ! [IN]
                             LAND_PROPERTY(:,:,I_WaterDiff),            & ! [IN]
                             LAND_SFLX_GH(:,:),                         & ! [IN]
                             LAND_SFLX_water(:,:), LAND_SFLX_ice(:,:),  & ! [IN]
                             LANDUSE_fact_land(:,:), LCDZ(:),           & ! [IN]
                             dt, NOWDAYSEC,                             & ! [IN]
                             LAND_TEMP(:,:,:), LAND_WATER(:,:,:),       & ! [INOUT]
                             RUNOFF(:,:)                                ) ! [OUT]
       call FILE_HISTORY_in( RUNOFF(:,:), 'RUNOFF', 'runoff water', 'kg', dim_type='XY' )
    end select

    !########## Negative Fixer ##########
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_WATER(k,i,j) = max( LAND_WATER(k,i,j), 0.0_RP )
    enddo
    enddo
    enddo

    call LAND_vars_total

    !########## History & Monitor ##########
    call LAND_vars_history


    call PROF_rapend  ('LND_Update', 1)

    return
  end subroutine LAND_driver_update

  !-----------------------------------------------------------------------------
  !> Get surface boundary from other model
  subroutine LAND_SURFACE_GET
    use mod_land_admin, only: &
       LAND_do
    use mod_land_vars, only: &
       ATMOS_TEMP,      &
       ATMOS_PRES,      &
       ATMOS_W,         &
       ATMOS_U,         &
       ATMOS_V,         &
       ATMOS_DENS,      &
       ATMOS_QV,        &
       ATMOS_PBL,       &
       ATMOS_SFC_DENS,  &
       ATMOS_SFC_PRES,  &
       ATMOS_SFLX_LW,   &
       ATMOS_SFLX_SW,   &
       ATMOS_cosSZA,    &
       ATMOS_SFLX_prec, &
       ATMOS_SFLX_rain, &
       ATMOS_SFLX_snow
    use mod_cpl_vars, only: &
       CPL_getATM_LND
    implicit none

    real(RP) :: ATMOS_SFLX_rad_dn(LIA,LJA,2,2)

    integer  :: i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_SfcExch', 2)

    if ( LAND_do ) then
       call CPL_getATM_LND( ATMOS_TEMP       (:,:),     & ! [OUT]
                            ATMOS_PRES       (:,:),     & ! [OUT]
                            ATMOS_W          (:,:),     & ! [OUT]
                            ATMOS_U          (:,:),     & ! [OUT]
                            ATMOS_V          (:,:),     & ! [OUT]
                            ATMOS_DENS       (:,:),     & ! [OUT]
                            ATMOS_QV         (:,:),     & ! [OUT]
                            ATMOS_PBL        (:,:),     & ! [OUT]
                            ATMOS_SFC_DENS   (:,:),     & ! [OUT]
                            ATMOS_SFC_PRES   (:,:),     & ! [OUT]
                            ATMOS_SFLX_rad_dn(:,:,:,:), & ! [OUT]
                            ATMOS_cosSZA     (:,:),     & ! [OUT]
                            ATMOS_SFLX_rain  (:,:),     & ! [OUT]
                            ATMOS_SFLX_snow  (:,:)      ) ! [OUT]
    endif

!OCL XFILL
    !$omp parallel do
    do j = LJS, LJE
    do i = LIS, LIE
       ATMOS_SFLX_SW  (i,j) = ATMOS_SFLX_rad_dn(i,j,I_SW,1) + ATMOS_SFLX_rad_dn(i,j,I_SW,2) ! direct+diffuse
       ATMOS_SFLX_LW  (i,j) = ATMOS_SFLX_rad_dn(i,j,I_LW,1) + ATMOS_SFLX_rad_dn(i,j,I_LW,2) ! direct+diffuse

       ATMOS_SFLX_prec(i,j) = ATMOS_SFLX_rain(i,j) + ATMOS_SFLX_snow(i,j) ! liquid+ice
    enddo
    enddo

    call PROF_rapend  ('LND_SfcExch', 2)

    return
  end subroutine LAND_SURFACE_GET

  !-----------------------------------------------------------------------------
  !> Put surface boundary to other model
  subroutine LAND_SURFACE_SET( countup )
    use mod_land_admin, only: &
       LAND_do
    use mod_land_vars, only: &
       LAND_PROPERTY,   &
       I_Z0M,           &
       I_Z0H,           &
       I_Z0E,           &
       LAND_SFC_TEMP,   &
       LAND_SFC_albedo, &
       LAND_SFLX_MW,    &
       LAND_SFLX_MU,    &
       LAND_SFLX_MV,    &
       LAND_SFLX_SH,    &
       LAND_SFLX_LH,    &
       LAND_SFLX_GH,    &
       LAND_SFLX_evap,  &
       LAND_U10,        &
       LAND_V10,        &
       LAND_T2,         &
       LAND_Q2
    use mod_cpl_vars, only: &
       CPL_putLND
    implicit none

    ! arguments
    logical, intent(in) :: countup
    !---------------------------------------------------------------------------

    call PROF_rapstart('LND_SfcExch', 2)

    if ( LAND_do ) then
       call CPL_putLND( LAND_SFC_TEMP  (:,:),       & ! [IN]
                        LAND_SFC_albedo(:,:,:),     & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0M), & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0H), & ! [IN]
                        LAND_PROPERTY  (:,:,I_Z0E), & ! [IN]
                        LAND_SFLX_MW   (:,:),       & ! [IN]
                        LAND_SFLX_MU   (:,:),       & ! [IN]
                        LAND_SFLX_MV   (:,:),       & ! [IN]
                        LAND_SFLX_SH   (:,:),       & ! [IN]
                        LAND_SFLX_LH   (:,:),       & ! [IN]
                        LAND_SFLX_GH   (:,:),       & ! [IN]
                        LAND_SFLX_evap (:,:),       & ! [IN]
                        LAND_U10       (:,:),       & ! [IN]
                        LAND_V10       (:,:),       & ! [IN]
                        LAND_T2        (:,:),       & ! [IN]
                        LAND_Q2        (:,:),       & ! [IN]
                        countup                     ) ! [IN]
    endif

    call PROF_rapend  ('LND_SfcExch', 2)

    return
  end subroutine LAND_SURFACE_SET

end module mod_land_driver
