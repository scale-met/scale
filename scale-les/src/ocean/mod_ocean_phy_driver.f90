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
  use scale_debug
  use scale_grid_index

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OCEAN_PHY_driver_setup
    use scale_ocean_phy, only: &
       OCEAN_PHY_setup
    use scale_ocean_sfc, only: &
       OCEAN_SFC_setup
    use mod_admin_restart, only: &
       RESTART_RUN
    use mod_ocean_admin, only: &
       OCEAN_TYPE, &
       OCEAN_sw
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[OCEAN PHY] / Origin[SCALE-LES]'

    if ( OCEAN_sw ) then

       ! setup library component
       call OCEAN_PHY_setup( OCEAN_TYPE )
       call OCEAN_SFC_setup( OCEAN_TYPE )

       if( .NOT. RESTART_RUN ) then
          ! run once (only for the diagnostic value)
          call PROF_rapstart('OCN Physics', 1)
          call OCEAN_PHY_driver( update_flag = .true. )
          call PROF_rapend  ('OCN Physics', 1)
       else
          ! no update in order to use restart value
       end if

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine OCEAN_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine OCEAN_PHY_driver( update_flag )
    use scale_const, only: &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL,    &
       LHV0  => CONST_LHV0,  &
       TEM00 => CONST_TEM00
    use scale_time, only: &
       dt => TIME_DTSEC_OCEAN
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_history, only: &
       HIST_in
    use scale_grid_real, only: &
       REAL_Z1
    use scale_roughness, only: &
       ROUGHNESS
    use scale_ocean_phy, only: &
       OCEAN_PHY
    use scale_ocean_sfc, only: &
       OCEAN_SFC
    use mod_ocean_vars, only: &
       OCEAN_TEMP,         &
       OCEAN_SFC_TEMP,     &
       OCEAN_SFC_albedo,   &
       OCEAN_SFC_Z0M,      &
       OCEAN_SFC_Z0H,      &
       OCEAN_SFC_Z0E,      &
       OCEAN_TEMP_t,       &
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
       OCEAN_SFLX_WH,      &
       OCEAN_SFLX_evap,    &
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
       ATMOS_SFC_PRES,     &
       ATMOS_SFLX_LW,      &
       ATMOS_SFLX_SW,      &
       ATMOS_SFLX_prec
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: total ! dummy
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       call ROUGHNESS( OCEAN_SFC_Z0M_t(:,:), & ! [OUT]
                       OCEAN_SFC_Z0H_t(:,:), & ! [OUT]
                       OCEAN_SFC_Z0E_t(:,:), & ! [OUT]
                       OCEAN_SFC_Z0M  (:,:), & ! [IN]
                       OCEAN_SFC_Z0H  (:,:), & ! [IN]
                       OCEAN_SFC_Z0E  (:,:), & ! [IN]
                       ATMOS_U        (:,:), & ! [IN]
                       ATMOS_V        (:,:), & ! [IN]
                       REAL_Z1        (:,:), & ! [IN]
                       dt                    ) ! [IN]

       call OCEAN_SFC( OCEAN_SFC_TEMP_t(:,:),      & ! [OUT]
                       OCEAN_SFLX_MW   (:,:),      & ! [OUT]
                       OCEAN_SFLX_MU   (:,:),      & ! [OUT]
                       OCEAN_SFLX_MV   (:,:),      & ! [OUT]
                       OCEAN_SFLX_SH   (:,:),      & ! [OUT]
                       OCEAN_SFLX_LH   (:,:),      & ! [OUT]
                       OCEAN_SFLX_WH   (:,:),      & ! [OUT]
                       OCEAN_U10       (:,:),      & ! [OUT]
                       OCEAN_V10       (:,:),      & ! [OUT]
                       OCEAN_T2        (:,:),      & ! [OUT]
                       OCEAN_Q2        (:,:),      & ! [OUT]
                       ATMOS_TEMP      (:,:),      & ! [IN]
                       ATMOS_PRES      (:,:),      & ! [IN]
                       ATMOS_W         (:,:),      & ! [IN]
                       ATMOS_U         (:,:),      & ! [IN]
                       ATMOS_V         (:,:),      & ! [IN]
                       ATMOS_DENS      (:,:),      & ! [IN]
                       ATMOS_QV        (:,:),      & ! [IN]
                       REAL_Z1         (:,:),      & ! [IN]
                       ATMOS_PBL       (:,:),      & ! [IN]
                       ATMOS_SFC_PRES  (:,:),      & ! [IN]
                       ATMOS_SFLX_LW   (:,:),      & ! [IN]
                       ATMOS_SFLX_SW   (:,:),      & ! [IN]
                       OCEAN_TEMP      (:,:),      & ! [IN]
                       OCEAN_SFC_TEMP  (:,:),      & ! [IN]
                       OCEAN_SFC_albedo(:,:,I_LW), & ! [IN]
                       OCEAN_SFC_albedo(:,:,I_SW), & ! [IN]
                       OCEAN_SFC_Z0M   (:,:),      & ! [IN]
                       OCEAN_SFC_Z0H   (:,:),      & ! [IN]
                       OCEAN_SFC_Z0E   (:,:),      & ! [IN]
                       dt                          ) ! [IN]

       OCEAN_SFLX_evap(:,:) = OCEAN_SFLX_LH(:,:) / ( LHV0 + ( CPvap-CL ) * ( ATMOS_TEMP(:,:)-TEM00 ) )

       call OCEAN_PHY( OCEAN_TEMP_t   (:,:), & ! [OUT]
                       OCEAN_TEMP     (:,:), & ! [IN]
                       OCEAN_SFLX_WH  (:,:), & ! [IN]
                       ATMOS_SFLX_prec(:,:), & ! [IN]
                       OCEAN_SFLX_evap(:,:), & ! [IN]
                       dt                    ) ! [IN]

       ! no albedo update (tentative)
       OCEAN_SFC_albedo_t(:,:,:) = 0.0_RP

       call HIST_in( OCEAN_TEMP_t      (:,:),      'OCEAN_TEMP_t',     'tendency of OCEAN_TEMP',     'K'   )
       call HIST_in( OCEAN_SFC_TEMP_t  (:,:),      'OCEAN_SFC_TEMP_t', 'tendency of OCEAN_SFC_TEMP', 'K'   )
       call HIST_in( OCEAN_SFC_albedo_t(:,:,I_LW), 'OCEAN_ALB_LW_t',   'tendency of OCEAN_ALB_LW',   '0-1' )
       call HIST_in( OCEAN_SFC_albedo_t(:,:,I_SW), 'OCEAN_ALB_SW_t',   'tendency of OCEAN_ALB_SW',   '0-1' )
       call HIST_in( OCEAN_SFC_Z0M_t   (:,:),      'OCEAN_SFC_Z0M_t',  'tendency of OCEAN_SFC_Z0M',  'm'   )
       call HIST_in( OCEAN_SFC_Z0H_t   (:,:),      'OCEAN_SFC_Z0H_t',  'tendency of OCEAN_SFC_Z0H',  'm'   )
       call HIST_in( OCEAN_SFC_Z0E_t   (:,:),      'OCEAN_SFC_Z0E_t',  'tendency of OCEAN_SFC_Z0E',  'm'   )

    end if

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, OCEAN_TEMP_t      (:,:),      'OCEAN_TEMP_t'     )
       call STAT_total( total, OCEAN_SFC_TEMP_t  (:,:),      'OCEAN_SFC_TEMP_t' )
       call STAT_total( total, OCEAN_SFC_albedo_t(:,:,I_LW), 'OCEAN_ALB_LW_t'   )
       call STAT_total( total, OCEAN_SFC_albedo_t(:,:,I_SW), 'OCEAN_ALB_SW_t'   )
       call STAT_total( total, OCEAN_SFC_Z0M_t   (:,:),      'OCEAN_SFC_Z0M_t'  )
       call STAT_total( total, OCEAN_SFC_Z0H_t   (:,:),      'OCEAN_SFC_Z0H_t'  )
       call STAT_total( total, OCEAN_SFC_Z0E_t   (:,:),      'OCEAN_SFC_Z0E_t'  )
    end if

    return
  end subroutine OCEAN_PHY_driver

end module mod_ocean_phy_driver
