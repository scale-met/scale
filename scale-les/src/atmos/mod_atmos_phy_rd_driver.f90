!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process wrapper
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-06 (S.Nishizawa)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_stdio
  use mod_prof
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_driver_setup
  public :: ATMOS_PHY_RD_driver

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
  real(RP), private, allocatable :: RHOT_t(:,:,:)

  integer, private :: ktop
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_driver_setup( RD_TYPE )
    use mod_stdio, only: &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_phy_rd, only: &
       ATMOS_PHY_RD_setup
    implicit none

    character(len=IO_SYSCHR), intent(in) :: RD_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'

    allocate( RHOT_t(KA,IA,JA) )

    call ATMOS_PHY_RD_setup( RD_TYPE )

    call ATMOS_PHY_RD_driver( .true., .false. )

    return
  end subroutine ATMOS_PHY_RD_driver_setup

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_driver( update_flag, history_flag )
    use mod_history, only: &
       HIST_in
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_RD
    use mod_grid, only: &
       CZ => GRID_CZ, &
       FZ => GRID_FZ, &
       CDZ => GRID_CDZ, &
       RCDZ => GRID_RCDZ
    use mod_time, only: &
       TIME_NOWDATE
    use mod_grid_real, only: &
       REAL_lon, &
       REAL_lat
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total, &
       DENS, &
       RHOT, &
       QTRC, &
       RHOT_tp
    use mod_atmos_vars_sf, only: &
       ATMOS_vars_sf_fillhalo, &
       SWD,                    &
       LWD
    use mod_cpl_vars, only: &
       sw_AtmLnd => CPL_sw_AtmLnd, &
       SkinT
    use mod_atmos_thermodyn, only: &
       THERMODYN_qd        => ATMOS_THERMODYN_qd,       &
       THERMODYN_cv        => ATMOS_THERMODYN_cv,       &
       THERMODYN_rhoe      => ATMOS_THERMODYN_rhoe,     &
       THERMODYN_rhot      => ATMOS_THERMODYN_rhot,     &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_atmos_phy_rd, only: &
       ATMOS_PHY_RD
    use mod_atmos_phy_rd_common, only: &
       RD_heating => ATMOS_PHY_RD_heating, &
       I_SW, &
       I_LW, &
       I_dn, &
       I_up
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in), optional :: history_flag

    real(RP) :: RHOE    (KA,IA,JA)
    real(RP) :: RHOE_t  (KA,IA,JA,2)
    real(RP) :: RHOT_tmp(KA,IA,JA)
    real(RP) :: TEMP_t  (KA,IA,JA)
    real(RP) :: QDRY    (KA,IA,JA)
    real(RP) :: CVtot   (KA,IA,JA)
    real(RP) :: temp    (KA,IA,JA)
    real(RP) :: pres    (KA,IA,JA)

    real(RP) :: flux_rad(KA,IA,JA,2,2)
    real(RP) :: flux_net(KA,IA,JA,2)
    real(RP) :: flux_up (KA,IA,JA,2)
    real(RP) :: flux_dn (KA,IA,JA,2)
    real(RP) :: flux_rad_top(IA,JA,2)
    real(RP) :: flux_rad_sfc(IA,JA,2)

    real(RP) :: solins  (IA,JA)
    real(RP) :: cosSZA  (IA,JA)
    real(RP) :: temp_sfc(IA,JA)

    real(RP) :: param_sfc(5) !< surface parameter

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       if( sw_AtmLnd ) then
          temp_sfc(:,:) = SkinT(:,:)
       else
          call THERMODYN_temp_pres( temp(:,:,:),  & ! [OUT]
                                    pres(:,:,:),  & ! [OUT]
                                    DENS(:,:,:),  & ! [IN]
                                    RHOT(:,:,:),  & ! [IN]
                                    QTRC(:,:,:,:) ) ! [IN]
          temp_sfc(:,:) = temp(KS,:,:)
       end if

       ! surface parameter
       param_sfc(1)   = 2.D0 ! land surface only, tentative
       param_sfc(2:5) = 0.D0

       call ATMOS_PHY_RD( &
            flux_rad, flux_rad_top, & ! [out]
            solins, cosSZA,         & ! [out]
            DENS, RHOT, QTRC,       & ! [in]
            temp_sfc, param_sfc,    & ! [in]
            CZ, FZ, CDZ, RCDZ,      & ! [in]
            REAL_lon, REAL_lat,     & ! [in]
            TIME_NOWDATE            ) ! [in]

       call THERMODYN_rhoe( RHOE(:,:,:),  & ! [OUT]
                            RHOT(:,:,:),  & ! [IN]
                            QTRC(:,:,:,:) ) ! [IN]

       ! apply radiative flux convergence -> heating rate
       call RD_heating( RHOE_t  (:,:,:,:),  & ! [OUT]
                        RHOE    (:,:,:),    & ! [INOUT]
                        flux_rad(:,:,:,:,:) ) ! [IN]

       ! update rhot
       call THERMODYN_rhot( RHOT_tmp(:,:,:),& ! [OUT]
                            RHOE(:,:,:),    & ! [IN]
                            QTRC(:,:,:,:)   ) ! [IN]

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          RHOT_t(k,i,j) = ( RHOT_tmp(k,i,j)-RHOT(k,i,j) ) / dt
       enddo
       enddo
       enddo

       if( sw_AtmLnd ) then
          do j = JS, JE
          do i = IS, IE
             LWD(i,j) = flux_rad(KS-1,i,j,I_LW,I_dn)
             SWD(i,j) = flux_rad(KS-1,i,j,I_SW,I_dn)
          enddo
          enddo
       end if

     if ( present(history_flag) ) then
     if ( history_flag ) then

        do j = JS, JE
        do i = IS, IE
           do k = KS, KE
              flux_net(k,i,j,I_LW) = &
                   ( ( flux_rad(k  ,i,j,I_LW,I_up) - flux_rad(k-1,i,j,I_LW,I_dn) ) &
                   + ( flux_rad(k  ,i,j,I_LW,I_up) - flux_rad(k  ,i,j,I_LW,I_dn) ) ) * 0.5_RP
              flux_net(k,i,j,I_SW) = &
                   ( ( flux_rad(k-1,i,j,I_SW,I_up) - flux_rad(k-1,i,j,I_SW,I_dn) ) &
                   + ( flux_rad(k  ,i,j,I_SW,I_up) - flux_rad(k  ,i,j,I_SW,I_dn) ) ) * 0.5_RP

              flux_up (k,i,j,I_LW) = &
                   ( flux_rad(k-1,i,j,I_LW,I_up) + flux_rad(k,i,j,I_LW,I_up) ) * 0.5_RP
              flux_up (k,i,j,I_SW) = &
                   ( flux_rad(k-1,i,j,I_SW,I_up) + flux_rad(k,i,j  ,I_SW,I_up) ) * 0.5_RP
              flux_dn (k,i,j,I_LW) = &
                   ( flux_rad(k-1,i,j,I_LW,I_dn) + flux_rad(k,i,j  ,I_LW,I_dn) ) * 0.5_RP
              flux_dn (k,i,j,I_SW) = &
                   ( flux_rad(k-1,i,j,I_SW,I_dn) + flux_rad(k,i,j  ,I_SW,I_dn) ) * 0.5_RP
           enddo
           flux_rad_sfc(i,j,I_LW) = flux_rad(KS-1,i,j,I_LW,I_up)-flux_rad(KS-1,i,j,I_LW,I_dn)
           flux_rad_sfc(i,j,I_SW) = flux_rad(KS-1,i,j,I_SW,I_up)-flux_rad(KS-1,i,j,I_SW,I_dn)

        enddo
        enddo

        call THERMODYN_qd( QDRY(:,:,:),  & ! [OUT]
                           QTRC(:,:,:,:) ) ! [IN]
        call THERMODYN_cv( CVtot(:,:,:),   & ! [OUT]
                           QTRC (:,:,:,:), & ! [IN]
                           QDRY (:,:,:)    ) ! [IN]

        TEMP_t(:,:,:) = RHOE_t(:,:,:,I_LW) / CVtot(:,:,:) * 86400.0_RP
        call HIST_in( TEMP_t(:,:,:), 'TEMP_t_rd_LW', 'tendency of temp in rd(LW)', 'K/day', dt )
        TEMP_t(:,:,:) = RHOE_t(:,:,:,I_SW) / CVtot(:,:,:) * 86400.0_RP
        call HIST_in( TEMP_t(:,:,:), 'TEMP_t_rd_SW', 'tendency of temp in rd(SW)', 'K/day', dt )
        TEMP_t(:,:,:) = ( RHOE_t(:,:,:,I_LW) + RHOE_t(:,:,:,I_SW) ) / CVtot(:,:,:) * 86400.0_RP
        call HIST_in( TEMP_t(:,:,:), 'TEMP_t_rd', 'tendency of temp in rd', 'K/day', dt )

        call HIST_in( flux_net(:,:,:,I_LW), 'RADFLUX_LW', 'net radiation flux(LW)', 'W/m2', dt )
        call HIST_in( flux_net(:,:,:,I_SW), 'RADFLUX_SW', 'net radiation flux(SW)', 'W/m2', dt )
        call HIST_in( flux_up (:,:,:,I_LW), 'RADFLUX_LWUP', 'up radiation flux(LW)', 'W/m2', dt )
        call HIST_in( flux_up (:,:,:,I_SW), 'RADFLUX_SWUP', 'up radiation flux(SW)', 'W/m2', dt )
        call HIST_in( flux_dn (:,:,:,I_LW), 'RADFLUX_LWDN', 'dn radiation flux(LW)', 'W/m2', dt )
        call HIST_in( flux_dn (:,:,:,I_SW), 'RADFLUX_SWDN', 'dn radiation flux(SW)', 'W/m2', dt )

        call HIST_in( flux_rad_top(:,:,I_LW), 'OLR', 'TOA     longwave  radiation', 'W/m2', dt )
        call HIST_in( flux_rad_top(:,:,I_SW), 'OSR', 'TOA     shortwave radiation', 'W/m2', dt )
        call HIST_in( flux_rad_sfc(:,:,I_LW), 'SLR', 'Surface longwave  radiation', 'W/m2', dt )
        call HIST_in( flux_rad_sfc(:,:,I_LW), 'SSR', 'Surface shortwave radiation', 'W/m2', dt )

        call HIST_in( solins(:,:), 'SOLINS', 'solar insolation', 'W/m2', dt )
        call HIST_in( cosSZA(:,:), 'COSZ', 'cos(solar zenith angle)', '0-1', dt )
      endif
      endif

      ! fill halo
      call ATMOS_vars_fillhalo
      call ATMOS_vars_sf_fillhalo

      ! check total (optional)
      call ATMOS_vars_total

    endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + RHOT_t(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_driver

end module mod_atmos_phy_rd_driver
