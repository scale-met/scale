!-------------------------------------------------------------------------------
!> module URBAN / Physics Urban Canopy Model (UCM)
!!
!! @par Description
!!          Urban physics module
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_urban_phy_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_urban_grid_cartesC_index

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
  public :: URBAN_PHY_driver_setup
  public :: URBAN_PHY_driver_resume
  public :: URBAN_PHY_driver

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
  subroutine URBAN_PHY_driver_setup
    use scale_urban_phy, only: &
       URBAN_PHY_setup
    use mod_urban_admin, only: &
       URBAN_TYPE, &
       URBAN_sw
    use mod_urban_vars, only: &
       URBAN_Z0M, &
       URBAN_Z0H, &
       URBAN_Z0E
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[URBAN PHY] / Origin[SCALE-RM]'

    if ( URBAN_sw ) then

       ! setup library component
       call URBAN_PHY_setup( URBAN_TYPE,     & ! [IN]
                             URBAN_Z0M(:,:), & ! [OUT]
                             URBAN_Z0H(:,:), & ! [OUT]
                             URBAN_Z0E(:,:)  ) ! [OUT]

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine URBAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Resume
  subroutine URBAN_PHY_driver_resume
    use mod_admin_restart, only: &
       RESTART_RUN
    use mod_urban_admin, only: &
       URBAN_sw
    implicit none

    if ( URBAN_sw ) then

       if ( .NOT. RESTART_RUN ) then ! tentative
          ! run once (only for the diagnostic value)
          call PROF_rapstart('URB_Physics', 1)
          call URBAN_PHY_driver( update_flag = .true. )
          call PROF_rapend  ('URB_Physics', 1)
       end if

    end if

    return
  end subroutine URBAN_PHY_driver_resume
  !-----------------------------------------------------------------------------
  !> Driver
  subroutine URBAN_PHY_driver( update_flag )
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE,     &
       dt      => TIME_DTSEC_URBAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_urban_grid_cartesC_real, only: &
       URBAN_GRID_CARTESC_REAL_VOL,    &
       URBAN_GRID_CARTESC_REAL_TOTVOL, &
       URBAN_GRID_CARTESC_REAL_AREA,   &
       URBAN_GRID_CARTESC_REAL_TOTAREA
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_mapprojection, only: &
       BASE_LON => MAPPROJECTION_basepoint_lon, &
       BASE_LAT => MAPPROJECTION_basepoint_lat
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_Z1
    use scale_urban_phy, only: &
       URBAN_PHY
    use mod_urban_vars, only: &
       URBAN_TR,         &
       URBAN_TB,         &
       URBAN_TG,         &
       URBAN_TC,         &
       URBAN_QC,         &
       URBAN_UC,         &
       URBAN_TRL,        &
       URBAN_TBL,        &
       URBAN_TGL,        &
       URBAN_RAINR,      &
       URBAN_RAINB,      &
       URBAN_RAING,      &
       URBAN_ROFF,       &
       URBAN_TR_t,       &
       URBAN_TB_t,       &
       URBAN_TG_t,       &
       URBAN_TC_t,       &
       URBAN_QC_t,       &
       URBAN_UC_t,       &
       URBAN_TRL_t,      &
       URBAN_TBL_t,      &
       URBAN_TGL_t,      &
       URBAN_RAINR_t,    &
       URBAN_RAINB_t,    &
       URBAN_RAING_t,    &
       URBAN_ROFF_t,     &
       URBAN_SFC_TEMP,   &
       URBAN_SFC_albedo, &
       URBAN_SFLX_MW,    &
       URBAN_SFLX_MU,    &
       URBAN_SFLX_MV,    &
       URBAN_SFLX_SH,    &
       URBAN_SFLX_LH,    &
       URBAN_SFLX_GH,    &
       URBAN_SFLX_evap,  &
       URBAN_Z0M,        &
       URBAN_Z0H,        &
       URBAN_Z0E,        &
       URBAN_U10,        &
       URBAN_V10,        &
       URBAN_T2,         &
       URBAN_Q2,         &
       ATMOS_TEMP,       &
       ATMOS_PRES,       &
       ATMOS_W,          &
       ATMOS_U,          &
       ATMOS_V,          &
       ATMOS_DENS,       &
       ATMOS_QV,         &
       ATMOS_PBL,        &
       ATMOS_SFC_DENS,   &
       ATMOS_SFC_PRES,   &
       ATMOS_SFLX_LW,    &
       ATMOS_SFLX_SW,    &
       ATMOS_SFLX_prec
    implicit none

    ! arguments
    logical, intent(in) :: update_flag

    ! works
    real(RP) :: LHV(UIA,UJA) ! latent heat of vaporization [J/kg]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       call URBAN_PHY( URBAN_TR_t      (:,:),      & ! [OUT]
                       URBAN_TB_t      (:,:),      & ! [OUT]
                       URBAN_TG_t      (:,:),      & ! [OUT]
                       URBAN_TC_t      (:,:),      & ! [OUT]
                       URBAN_QC_t      (:,:),      & ! [OUT]
                       URBAN_UC_t      (:,:),      & ! [OUT]
                       URBAN_TRL_t     (:,:,:),    & ! [OUT]
                       URBAN_TBL_t     (:,:,:),    & ! [OUT]
                       URBAN_TGL_t     (:,:,:),    & ! [OUT]
                       URBAN_RAINR_t   (:,:),      & ! [OUT]
                       URBAN_RAINB_t   (:,:),      & ! [OUT]
                       URBAN_RAING_t   (:,:),      & ! [OUT]
                       URBAN_ROFF_t    (:,:),      & ! [OUT]
                       URBAN_SFC_TEMP  (:,:),      & ! [OUT]
                       URBAN_SFC_albedo(:,:,I_LW), & ! [OUT]
                       URBAN_SFC_albedo(:,:,I_SW), & ! [OUT]
                       URBAN_SFLX_MW   (:,:),      & ! [OUT]
                       URBAN_SFLX_MU   (:,:),      & ! [OUT]
                       URBAN_SFLX_MV   (:,:),      & ! [OUT]
                       URBAN_SFLX_SH   (:,:),      & ! [OUT]
                       URBAN_SFLX_LH   (:,:),      & ! [OUT]
                       URBAN_SFLX_GH   (:,:),      & ! [OUT]
                       URBAN_Z0M       (:,:),      & ! [OUT]
                       URBAN_Z0H       (:,:),      & ! [OUT]
                       URBAN_Z0E       (:,:),      & ! [OUT]
                       URBAN_U10       (:,:),      & ! [OUT]
                       URBAN_V10       (:,:),      & ! [OUT]
                       URBAN_T2        (:,:),      & ! [OUT]
                       URBAN_Q2        (:,:),      & ! [OUT]
                       ATMOS_TEMP      (:,:),      & ! [IN]
                       ATMOS_PRES      (:,:),      & ! [IN]
                       ATMOS_W         (:,:),      & ! [IN]
                       ATMOS_U         (:,:),      & ! [IN]
                       ATMOS_V         (:,:),      & ! [IN]
                       ATMOS_DENS      (:,:),      & ! [IN]
                       ATMOS_QV        (:,:),      & ! [IN]
                       ATMOS_GRID_CARTESC_REAL_Z1         (:,:),      & ! [IN]
                       ATMOS_PBL       (:,:),      & ! [IN]
                       ATMOS_SFC_DENS  (:,:),      & ! [IN]
                       ATMOS_SFC_PRES  (:,:),      & ! [IN]
                       ATMOS_SFLX_LW   (:,:,:),    & ! [IN]
                       ATMOS_SFLX_SW   (:,:,:),    & ! [IN]
                       ATMOS_SFLX_prec (:,:),      & ! [IN]
                       URBAN_TR        (:,:),      & ! [IN]
                       URBAN_TB        (:,:),      & ! [IN]
                       URBAN_TG        (:,:),      & ! [IN]
                       URBAN_TC        (:,:),      & ! [IN]
                       URBAN_QC        (:,:),      & ! [IN]
                       URBAN_UC        (:,:),      & ! [IN]
                       URBAN_TRL       (:,:,:),    & ! [IN]
                       URBAN_TBL       (:,:,:),    & ! [IN]
                       URBAN_TGL       (:,:,:),    & ! [IN]
                       URBAN_RAINR     (:,:),      & ! [IN]
                       URBAN_RAINB     (:,:),      & ! [IN]
                       URBAN_RAING     (:,:),      & ! [IN]
                       URBAN_ROFF      (:,:),      & ! [IN]
                       BASE_LON,                   & ! [IN]
                       BASE_LAT,                   & ! [IN]
                       NOWDATE         (:),        & ! [IN]
                       dt                          ) ! [IN]

       call HYDROMETEOR_LHV( LHV(:,:), ATMOS_TEMP(:,:) )

!OCL XFILL
       do j = UJS, UJE
       do i = UIS, UIE
          URBAN_SFLX_evap(i,j) = URBAN_SFLX_LH(i,j) / LHV(i,j)
       end do
       end do

       call FILE_HISTORY_in( URBAN_TR_t(:,:), 'URBAN_TR_t', 'tendency of URBAN_TR', 'K',     dim_type='XY' )
       call FILE_HISTORY_in( URBAN_TB_t(:,:), 'URBAN_TB_t', 'tendency of URBAN_TB', 'K',     dim_type='XY' )
       call FILE_HISTORY_in( URBAN_TG_t(:,:), 'URBAN_TG_t', 'tendency of URBAN_TG', 'K',     dim_type='XY' )
       call FILE_HISTORY_in( URBAN_TC_t(:,:), 'URBAN_TC_t', 'tendency of URBAN_TC', 'K',     dim_type='XY' )
       call FILE_HISTORY_in( URBAN_QC_t(:,:), 'URBAN_QC_t', 'tendency of URBAN_QC', 'kg/kg', dim_type='XY' )
       call FILE_HISTORY_in( URBAN_UC_t(:,:), 'URBAN_UC_t', 'tendency of URBAN_UC', 'm/s',   dim_type='XY' )

       call FILE_HISTORY_in( URBAN_TRL_t(:,:,:), 'URBAN_TRL_t', 'tendency of URBAN_TRL', 'K', dim_type='UXY' )
       call FILE_HISTORY_in( URBAN_TBL_t(:,:,:), 'URBAN_TBL_t', 'tendency of URBAN_TBL', 'K', dim_type='UXY' )
       call FILE_HISTORY_in( URBAN_TGL_t(:,:,:), 'URBAN_TGL_t', 'tendency of URBAN_TGL', 'K', dim_type='UXY' )

       call FILE_HISTORY_in( URBAN_RAINR_t(:,:), 'URBAN_RAINR_t', 'tendency of URBAN_RAINR', 'K', dim_type='XY' )
       call FILE_HISTORY_in( URBAN_RAINB_t(:,:), 'URBAN_RAINB_t', 'tendency of URBAN_RAINB', 'K', dim_type='XY' )
       call FILE_HISTORY_in( URBAN_RAING_t(:,:), 'URBAN_RAING_t', 'tendency of URBAN_RAING', 'K', dim_type='XY' )
       call FILE_HISTORY_in( URBAN_ROFF_t (:,:), 'URBAN_ROFF_t',  'tendency of URBAN_ROFF',  'K', dim_type='XY' )

    endif

    if ( STATISTICS_checktotal ) then

       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TRL_t (:,:,:), 'URBAN_TRL_t', &
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTVOL      )
       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TBL_t (:,:,:), 'URBAN_TBL_t', &
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTVOL      )
       call STATISTICS_total( UKA, UKS, UKE, UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TGL_t (:,:,:), 'URBAN_TGL_t', &
                              URBAN_GRID_CARTESC_REAL_VOL(:,:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTVOL      )

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TR_t(:,:), 'URBAN_TR_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TB_t(:,:), 'URBAN_TB_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TG_t(:,:), 'URBAN_TG_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_TC_t(:,:), 'URBAN_TC_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_QC_t(:,:), 'URBAN_QC_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_UC_t(:,:), 'URBAN_UC_t',     &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:), &
                              URBAN_GRID_CARTESC_REAL_TOTAREA    )

       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAINR_t(:,:), 'URBAN_RAINR_t', &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAINB_t(:,:), 'URBAN_RAINB_t', &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_RAING_t(:,:), 'URBAN_RAING_t', &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
       call STATISTICS_total( UIA, UIS, UIE, UJA, UJS, UJE, &
                              URBAN_ROFF_t (:,:), 'URBAN_ROFF_t',  &
                              URBAN_GRID_CARTESC_REAL_AREA(:,:),   &
                              URBAN_GRID_CARTESC_REAL_TOTAREA      )
    endif

    return
  end subroutine URBAN_PHY_driver

end module mod_urban_phy_driver
