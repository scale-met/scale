!-------------------------------------------------------------------------------
!> module LAND driver
!!
!! @par Description
!!          Land model driver
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
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
  public :: LAND_driver_resume
  public :: LAND_driver
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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_driver_setup
    use scale_prc, only: &
       PRC_abort
    use mod_land_phy_driver, only: &
       LAND_PHY_driver_setup
    use mod_land_admin, only: &
       LAND_do, &
       LAND_DYN_TYPE
    use scale_land_dyn_bucket, only: &
       LAND_DYN_BUCKET_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER Setup] / Categ[LAND] / Origin[SCALE-RM]'

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

       call LAND_PHY_driver_setup

    end if

    return
  end subroutine LAND_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine LAND_driver_resume
    use mod_land_phy_driver, only: &
       LAND_PHY_driver_resume
    use mod_land_vars, only: &
       LAND_vars_history, &
       LAND_TEMP_t,     &
       LAND_WATER_t,    &
       LAND_SFC_TEMP_t, &
       LAND_SFC_albedo_t
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[LAND] / Origin[SCALE-RM]'

    !########## Get Surface Boundary from coupler ##########
    call LAND_SURFACE_GET

    !########## initilize tendency ##########
!OCL XFILL
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_TEMP_t (k,i,j) = 0.0_RP
       LAND_WATER_t(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = LJS, LJE
    do i = LIS, LIE
       LAND_SFC_TEMP_t  (i,j)      = 0.0_RP
       LAND_SFC_albedo_t(i,j,I_LW) = 0.0_RP
       LAND_SFC_albedo_t(i,j,I_SW) = 0.0_RP
    enddo
    enddo

    ! resume each component
    call LAND_PHY_driver_resume

    !########## Set Surface Boundary to coupler ##########
    call LAND_SURFACE_SET( countup=.true. )

    !########## History & Monitor ##########
    call PROF_rapstart('LND_History', 1)
    call LAND_vars_history
    call PROF_rapend  ('LND_History', 1)

    return
  end subroutine LAND_driver_resume

  !-----------------------------------------------------------------------------
  !> Land step
  subroutine LAND_driver
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
       LAND_SFC_TEMP_t,   &
       LAND_SFC_albedo_t, &
       LAND_vars_total,   &
       LAND_vars_history
    use mod_land_phy_driver, only: &
       LAND_PHY_driver
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

    !########## Get Surface Boundary from coupler ##########
    call PROF_rapstart('LND_SfcExch', 2)
    call LAND_SURFACE_GET
    call PROF_rapend  ('LND_SfcExch', 2)

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

    do j = LJS, LJE
    do i = LIS, LIE
       LAND_SFC_TEMP  (i,j)      = LAND_SFC_TEMP  (i,j)      + LAND_SFC_TEMP_t  (i,j)      * dt
       LAND_SFC_albedo(i,j,I_LW) = LAND_SFC_albedo(i,j,I_LW) + LAND_SFC_albedo_t(i,j,I_LW) * dt
       LAND_SFC_albedo(i,j,I_SW) = LAND_SFC_albedo(i,j,I_SW) + LAND_SFC_albedo_t(i,j,I_SW) * dt
    enddo
    enddo

    !########## Negative Fixer ##########
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_WATER(k,i,j) = max( LAND_WATER(k,i,j), 0.0_RP )
    enddo
    enddo
    enddo

    call LAND_vars_total

    !########## Set Surface Boundary to coupler ##########
    call PROF_rapstart('LND_SfcExch', 2)
    call LAND_SURFACE_SET( countup=.true. )
    call PROF_rapend  ('LND_SfcExch', 2)

    !########## reset tendencies ##########
!OCL XFILL
    do j = LJS, LJE
    do i = LIS, LIE
    do k = LKS, LKE
       LAND_TEMP_t (k,i,j) = 0.0_RP
       LAND_WATER_t(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

!OCL XFILL
    do j = LJS, LJE
    do i = LIS, LIE
       LAND_SFC_TEMP_t  (i,j)      = 0.0_RP
       LAND_SFC_albedo_t(i,j,I_LW) = 0.0_RP
       LAND_SFC_albedo_t(i,j,I_SW) = 0.0_RP
    enddo
    enddo

    !########## Physics ##########
    call PROF_rapstart('LND_Physics', 1)
    call LAND_PHY_driver( update_flag = .true. )
    call PROF_rapend  ('LND_Physics', 1)

    !########## History & Monitor ##########
    call PROF_rapstart('LND_History', 1)
    call LAND_vars_history
    call PROF_rapend  ('LND_History', 1)

    return
  end subroutine LAND_driver

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
    do j = LJS, LJE
    do i = LIS, LIE
       ATMOS_SFLX_SW  (i,j) = ATMOS_SFLX_rad_dn(i,j,I_SW,1) + ATMOS_SFLX_rad_dn(i,j,I_SW,2) ! direct+diffuse
       ATMOS_SFLX_LW  (i,j) = ATMOS_SFLX_rad_dn(i,j,I_LW,1) + ATMOS_SFLX_rad_dn(i,j,I_LW,2) ! direct+diffuse

       ATMOS_SFLX_prec(i,j) = ATMOS_SFLX_rain(i,j) + ATMOS_SFLX_snow(i,j) ! liquid+ice
    enddo
    enddo

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

    return
  end subroutine LAND_SURFACE_SET

end module mod_land_driver
