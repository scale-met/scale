!-------------------------------------------------------------------------------
!> module OCEAN / Physics
!!
!! @par Description
!!          ocean physics module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module mod_ocean_phy_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_ocean_grid_cartesC_index

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
  public :: OCEAN_PHY_driver_setup
  public :: OCEAN_PHY_driver_resume
  public :: OCEAN_PHY_driver

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
  real(RP), allocatable :: QVEF(:,:)
  real(RP), allocatable :: SR(:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_driver_setup
    use scale_process, only: &
       PRC_abort
    use mod_ocean_admin, only: &
       OCEAN_SFC_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE, &
       OCEAN_do
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp_setup
    use scale_ocean_phy_albedo_nakajima00, only: &
       OCEAN_PHY_ALBEDO_nakajima00_setup
    use scale_ocean_phy_roughness_miller92, only: &
       OCEAN_PHY_ROUGHNESS_miller92_setup
    use scale_ocean_phy_roughness_moon07, only: &
       OCEAN_PHY_ROUGHNESS_moon07_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN PHY] / Origin[SCALE-RM]'

    if ( OCEAN_do ) then

       select case ( OCEAN_SFC_TYPE )
       case ( 'FIXED-TEMP' )
          call CPL_PHY_SFC_fixed_temp_setup
       case default
          write(*,*) 'xxx OCEAN_SFC_TYPE is invalid: ', trim(OCEAN_SFC_TYPE)
          call PRC_abort
       end select

       select case ( OCEAN_ALB_TYPE )
       case ( 'CONST' )
          ! do nothing
       case ( 'NAKAJIMA00' )
          call OCEAN_PHY_ALBEDO_nakajima00_setup
       case default
          write(*,*) 'xxx OCEAN_ALB_TYPE is invalid: ', trim(OCEAN_ALB_TYPE)
          call PRC_abort
       end select

       select case ( OCEAN_RGN_TYPE )
       case ( 'CONST' )
          ! do nothing
       case ( 'MILLER92' )
          call OCEAN_PHY_ROUGHNESS_miller92_setup
       case ( 'MOON07' )
          call OCEAN_PHY_ROUGHNESS_moon07_setup
       case default
          write(*,*) 'xxx OCEAN_RGN_TYPE is invalid: ', trim(OCEAN_RGN_TYPE)
          call PRC_abort
       end select

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    allocate( QVEF(OIA,OJA) )
    allocate( SR(OIA,OJA) )
    QVEF(:,:) = 1.0_RP
    SR(:,:) = 0.0_RP

    return
  end subroutine OCEAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine OCEAN_PHY_driver_resume
    use mod_admin_restart, only: &
       RESTART_RUN
    use mod_ocean_admin, only: &
       OCEAN_do
    implicit none

    if ( OCEAN_do ) then

       if ( .NOT. RESTART_RUN ) then ! tentative
          ! run once (only for the diagnostic value)
          call PROF_rapstart('OCN_Physics', 1)
          call OCEAN_PHY_driver( update_flag = .true. )
          call PROF_rapend  ('OCN_Physics', 1)
       end if

    end if

    return
  end subroutine OCEAN_PHY_driver_resume

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine OCEAN_PHY_driver( update_flag )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_VOL, &
       OCEAN_GRID_CARTESC_REAL_TOTVOL, &
       OCEAN_GRID_CARTESC_REAL_AREA, &
       OCEAN_GRID_CARTESC_REAL_TOTAREA
    use scale_atmos_grid_cartesC_real, only: &
       REAL_Z1 => ATMOS_GRID_CARTESC_REAL_Z1
    use scale_ocean_phy_albedo_nakajima00, only: &
       OCEAN_PHY_albedo_nakajima00
    use scale_ocean_phy_roughness_miller92, only: &
       OCEAN_PHY_ROUGHNESS_miller92
    use scale_ocean_phy_roughness_moon07, only: &
       OCEAN_PHY_ROUGHNESS_moon07
    use scale_cpl_phy_sfc_fixed_temp, only: &
       CPL_PHY_SFC_fixed_temp
    use mod_ocean_admin, only: &
       OCEAN_SFC_TYPE, &
       OCEAN_ALB_TYPE, &
       OCEAN_RGN_TYPE
    use mod_ocean_vars, only: &
       OCEAN_TEMP,         &
       OCEAN_SALT,         &
       OCEAN_UVEL,         &
       OCEAN_VVEL,         &
       OCEAN_SFC_TEMP,     &
       OCEAN_SFC_albedo,   &
       OCEAN_SFC_Z0M,      &
       OCEAN_SFC_Z0H,      &
       OCEAN_SFC_Z0E,      &
!!$       OCEAN_TEMP_t,       &
!!$       OCEAN_SALT_t,       &
!!$       OCEAN_UVEL_t,       &
!!$       OCEAN_VVEL_t,       &
       OCEAN_SFC_TEMP_t,   &
       OCEAN_SFC_albedo_t, &
       OCEAN_SFC_Z0M_t,    &
       OCEAN_SFC_Z0H_t,    &
       OCEAN_SFC_Z0E_t,    &
       OCEAN_SFLX_MW,      &
       OCEAN_SFLX_MU,      &
       OCEAN_SFLX_MV,      &
       OCEAN_SFLX_SH,      &
       OCEAN_SFLX_LH,      &
       OCEAN_SFLX_evap,    &
       OCEAN_SFLX_WH,      &
       OCEAN_SFLX_water,   &
       OCEAN_SFLX_ice,     &
       OCEAN_U10,          &
       OCEAN_V10,          &
       OCEAN_T2,           &
       OCEAN_Q2,           &
       ATMOS_TEMP,         &
       ATMOS_PRES,         &
       ATMOS_W,            &
       ATMOS_U,            &
       ATMOS_V,            &
       ATMOS_DENS,         &
       ATMOS_QV,           &
       ATMOS_PBL,          &
       ATMOS_cosSZA,       &
       ATMOS_SFC_DENS,     &
       ATMOS_SFC_PRES,     &
       ATMOS_SFLX_LW,      &
       ATMOS_SFLX_SW,      &
       ATMOS_SFLX_rain,    &
       ATMOS_SFLX_snow
    use scale_landuse, only: &
       LANDUSE_fact_ocean
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: SFC_TEMP(OIA,OJA)
    real(RP) :: Z0M(OIA,OJA), Z0H(OIA,OJA), Z0E(OIA,OJA)
    real(RP) :: albedo(OIA,OJA)
    real(RP) :: LHV(OIA,OJA) ! latent heat of vaporization [J/kg]
    real(RP) :: ATMOS_Uabs(OIA,OJA)

    integer  :: i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          ATMOS_Uabs(i,j) = sqrt( ATMOS_U(i,j)**2 + ATMOS_V(i,j)**2 )
       end do
       end do

       select case ( OCEAN_RGN_TYPE )
       case ( 'CONST' )
          ! do nothing
       case ( 'MILLER92' )
          call OCEAN_PHY_ROUGHNESS_miller92( OIA, OIS, OIE, OJA, OJS, OJE, &
                                             ATMOS_Uabs(:,:),             & ! [IN]
                                             Z0M(:,:), Z0H(:,:), Z0E(:,:) ) ! [OUT]
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             OCEAN_SFC_Z0M_t(i,j) = ( Z0M(i,j) - OCEAN_SFC_Z0M(i,j) ) / dt
             OCEAN_SFC_Z0H_t(i,j) = ( Z0H(i,j) - OCEAN_SFC_Z0H(i,j) ) / dt
             OCEAN_SFC_Z0E_t(i,j) = ( Z0E(i,j) - OCEAN_SFC_Z0E(i,j) ) / dt
          end do
          end do
       case ( 'MOON07' )
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             Z0M(i,j) = OCEAN_SFC_Z0M(i,j)
          end do
          end do
          call OCEAN_PHY_ROUGHNESS_moon07( OIA, OIS, OIE, OJA, OJS, OJE, &
                                           ATMOS_Uabs(:,:), REAL_Z1(:,:), & ! [IN]
                                           Z0M(:,:),                      & ! [INOUT]
                                           Z0H(:,:), Z0E(:,:)             ) ! [OUT]
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             OCEAN_SFC_Z0M_t(i,j) = ( Z0M(i,j) - OCEAN_SFC_Z0M(i,j) ) / dt
             OCEAN_SFC_Z0H_t(i,j) = ( Z0H(i,j) - OCEAN_SFC_Z0H(i,j) ) / dt
             OCEAN_SFC_Z0E_t(i,j) = ( Z0E(i,j) - OCEAN_SFC_Z0E(i,j) ) / dt
          end do
          end do
       end select

       select case ( OCEAN_ALB_TYPE )
       case ( 'CONST' )
          ! do nothing
       case ( 'NAKAJIMA00' )
          call OCEAN_PHY_albedo_nakajima00( OIA, OIS, OIE, OJA, OJS, OJE, &
                                            ATMOS_cosSZA(:,:), & ! [IN]
                                            albedo(:,:)        ) ! [OUT]
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             OCEAN_SFC_albedo_t(i,j,I_SW) = ( albedo(i,j) - OCEAN_SFC_albedo(i,j,I_SW) ) / dt
          end do
          end do
       end select


       call HYDROMETEOR_LHV( OIA, OIS, OIE, OJA, OJS, OJE, &
                             ATMOS_TEMP(:,:), LHV(:,:) )

       select case ( OCEAN_SFC_TYPE )
       case ( 'FIXED-TEMP' )
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             SFC_TEMP(i,j) = OCEAN_TEMP(OKS,i,j)
          end do
          end do

          call CPL_PHY_SFC_fixed_temp( OIA, OIS, OIE, OJA, OJS, OJE, &
                                       ATMOS_TEMP(:,:), ATMOS_PRES(:,:),                           & ! [IN]
                                       ATMOS_W(:,:), ATMOS_U(:,:), ATMOS_V(:,:),                   & ! [IN]
                                       ATMOS_DENS(:,:), ATMOS_QV(:,:), LHV(:,:),                   & ! [IN]
                                       REAL_Z1(:,:), ATMOS_PBL(:,:),                               & ! [IN]
                                       ATMOS_SFC_DENS(:,:), ATMOS_SFC_PRES(:,:),                   & ! [IN]
                                       ATMOS_SFLX_LW(:,:), ATMOS_SFLX_SW(:,:),                     & ! [IN]
                                       SFC_TEMP(:,:), QVEF(:,:),                                   & ! [IN]
                                       OCEAN_SFC_albedo(:,:,I_LW), OCEAN_SFC_albedo(:,:,I_SW),     & ! [IN]
                                       SR(:,:),                                                    & ! [IN]
                                       OCEAN_SFC_Z0M(:,:), OCEAN_SFC_Z0H(:,:), OCEAN_SFC_Z0E(:,:), & ! [IN]
                                       LANDUSE_fact_ocean, dt,                                     & ! [IN]
                                       OCEAN_SFLX_MW(:,:), OCEAN_SFLX_MU(:,:), OCEAN_SFLX_MV(:,:), & ! [OUT]
                                       OCEAN_SFLX_SH(:,:), OCEAN_SFLX_LH(:,:), OCEAN_SFLX_WH(:,:), & ! [OUT]
                                       OCEAN_U10(:,:), OCEAN_V10(:,:),                             & ! [OUT]
                                       OCEAN_T2(:,:), OCEAN_Q2(:,:)                                ) ! [OUT]
          !OCL XFILL
          !$omp parallel do
          do j = OJS, OJE
          do i = OIS, OIE
             OCEAN_SFC_TEMP_t(i,j) = ( SFC_TEMP(i,j) - OCEAN_SFC_TEMP(i,j) ) / dt
             OCEAN_SFLX_WH   (i,j) = - OCEAN_SFLX_WH(i,j) ! upward to downward
          end do
          end do

       end select

!OCL XFILL
       !$omp parallel do
       do j = OJS, OJE
       do i = OIS, OIE
          OCEAN_SFLX_evap(i,j) = OCEAN_SFLX_LH(i,j) / LHV(i,j)
          OCEAN_SFLX_water(i,j) = ATMOS_SFLX_rain(i,j) - OCEAN_SFLX_evap(i,j)
          OCEAN_SFLX_ice  (i,j) = ATMOS_SFLX_snow(i,j)
       end do
       end do


!       call FILE_HISTORY_in( OCEAN_TEMP_t      (:,:,:),    'OCEAN_TEMP_t',     'tendency of OCEAN_TEMP',     'K', dim_type='OXY' )
       call FILE_HISTORY_in( OCEAN_SFC_TEMP_t  (:,:),      'OCEAN_SFC_TEMP_t', 'tendency of OCEAN_SFC_TEMP', 'K', dim_type='XY' )
       call FILE_HISTORY_in( OCEAN_SFC_albedo_t(:,:,I_LW), 'OCEAN_ALB_LW_t',   'tendency of OCEAN_ALB_LW',   '1', dim_type='XY' )
       call FILE_HISTORY_in( OCEAN_SFC_albedo_t(:,:,I_SW), 'OCEAN_ALB_SW_t',   'tendency of OCEAN_ALB_SW',   '1', dim_type='XY' )
       call FILE_HISTORY_in( OCEAN_SFC_Z0M_t   (:,:),      'OCEAN_SFC_Z0M_t',  'tendency of OCEAN_SFC_Z0M',  'm', dim_type='XY' )
       call FILE_HISTORY_in( OCEAN_SFC_Z0H_t   (:,:),      'OCEAN_SFC_Z0H_t',  'tendency of OCEAN_SFC_Z0H',  'm', dim_type='XY' )
       call FILE_HISTORY_in( OCEAN_SFC_Z0E_t   (:,:),      'OCEAN_SFC_Z0E_t',  'tendency of OCEAN_SFC_Z0E',  'm', dim_type='XY' )

    end if

    if ( STATISTICS_checktotal ) then
!!$       call STATISTICS_total( OKA, OKS, OKE, OIA, OIS, OIE, OJA, OJS, OJE, &
!!$                              OCEAN_TEMP_t      (:,:,:),    'OCEAN_TEMP_t', &
!!$                              OCEAN_GRID_CARTESC_REAL_VOL(:,:,:),           &
!!$                              OCEAN_GRID_CARTESC_REAL_TOTVOL                )

       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_TEMP_t  (:,:),      'OCEAN_SFC_TEMP_t', &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   )
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_albedo_t(:,:,I_LW), 'OCEAN_ALB_LW_t',   &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   )
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_albedo_t(:,:,I_SW), 'OCEAN_ALB_SW_t',   &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   )
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_Z0M_t   (:,:),      'OCEAN_SFC_Z0M_t',  &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   )
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_Z0H_t   (:,:),      'OCEAN_SFC_Z0H_t',  &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   )
       call STATISTICS_total( OIA, OIS, OIE, OJA, OJS, OJE, &
                              OCEAN_SFC_Z0E_t   (:,:),      'OCEAN_SFC_Z0E_t',  &
                              OCEAN_GRID_CARTESC_REAL_AREA(:,:),                &
                              OCEAN_GRID_CARTESC_REAL_TOTAREA                   )
    end if

    return
  end subroutine OCEAN_PHY_driver

end module mod_ocean_phy_driver
