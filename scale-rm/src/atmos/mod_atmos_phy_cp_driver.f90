!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cumulus
!!
!! @par Description
!!          Cumulus parameterization driver
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-05-04 (H.Yashiro)    [new]
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module mod_atmos_phy_cp_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_driver_setup
  public :: ATMOS_PHY_CP_driver_resume
  public :: ATMOS_PHY_CP_driver

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
  subroutine ATMOS_PHY_CP_driver_setup
    use scale_atmos_phy_cp, only: &
         ATMOS_PHY_CP_setup
    use scale_atmos_phy_cp_common, only: &
         ATMOS_PHY_CP_common_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_CP_TYPE, &
       ATMOS_sw_phy_cp
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[DRIVER] / Categ[ATMOS PHY_CP] / Origin[SCALE-RM]'

    if ( ATMOS_sw_phy_cp ) then

       ! setup library component
       call ATMOS_PHY_CP_setup( ATMOS_PHY_CP_TYPE )
       call ATMOS_PHY_CP_common_setup

    else
       if( IO_L ) write(IO_FID_LOG,*) '*** this component is never called.'
    endif

    return
  end subroutine ATMOS_PHY_CP_driver_setup

  !-----------------------------------------------------------------------------
  !> Redume
  subroutine ATMOS_PHY_CP_driver_resume
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_cp
    implicit none

    if ( ATMOS_sw_phy_cp ) then

       ! run once (only for the diagnostic value)
       call PROF_rapstart('ATM_Cumulus', 1)
       call ATMOS_PHY_CP_driver( update_flag = .true. )
       call PROF_rapend  ('ATM_Cumulus', 1)

    endif

    return
  end subroutine ATMOS_PHY_CP_driver_resume

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_CP_driver( update_flag )
    use scale_rm_statistics, only: &
       STATISTICS_checktotal, &
       STAT_total
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_diagnose_number_concentration
    use scale_atmos_phy_cp, only: &
       ATMOS_PHY_CP
    use scale_atmos_phy_cp_common, only: &
       ATMOS_PHY_CP_wmean
    use scale_atmos_phy_mp, only: &
       QS_MP, &
       QE_MP
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use mod_atmos_vars, only: &
       DENS   => DENS_av, &
       MOMZ   => MOMZ_av, &
       MOMX   => MOMX_av, &
       MOMY   => MOMY_av, &
       RHOT   => RHOT_av, &
       QTRC   => QTRC_av, &
       DENS_t => DENS_tp, &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_cp_vars, only: &
       DENS_t_CP      => ATMOS_PHY_CP_DENS_t,         &
       MOMZ_t_CP      => ATMOS_PHY_CP_MOMZ_t,         &
       MOMX_t_CP      => ATMOS_PHY_CP_MOMX_t,         &
       MOMY_t_CP      => ATMOS_PHY_CP_MOMY_t,         &
       RHOT_t_CP      => ATMOS_PHY_CP_RHOT_t,         &
       RHOQ_t_CP      => ATMOS_PHY_CP_RHOQ_t,         &
       MFLX_cloudbase => ATMOS_PHY_CP_MFLX_cloudbase, &
       SFLX_rain      => ATMOS_PHY_CP_SFLX_rain,      &  ! convective rain [kg/m2/s]
       cloudtop       => ATMOS_PHY_CP_cloudtop,       &  ! cloud top height [m]
       cloudbase      => ATMOS_PHY_CP_cloudbase,      &  ! cloud base height [m]
       cldfrac_dp     => ATMOS_PHY_CP_cldfrac_dp,     &  ! cloud fraction (deep convection) (0-1)
       cldfrac_sh     => ATMOS_PHY_CP_cldfrac_sh,     &  ! cloud fraction (shallow convection) (0-1)
       w0avg          => ATMOS_PHY_CP_w0avg,          &  ! running mean vertical wind velocity [m/s]
       kf_nca         => ATMOS_PHY_CP_kf_nca             ! advection/cumulus convection timescale/dt for KF [step]

    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: total ! dummy

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    ! temporal running mean of vertical velocity
    call ATMOS_PHY_CP_wmean( w0avg(:,:,:), & ! [OUT]
                             DENS (:,:,:), & ! [IN]
                             MOMZ (:,:,:)  ) ! [IN]
    call FILE_HISTORY_in( w0avg(:,:,:), 'w0avg', 'running mean vertical wind velocity', 'kg/m2/s', fill_halo=.true. )

    if ( update_flag ) then ! update

       call ATMOS_PHY_CP( DENS,           & ! [IN]
                          MOMZ,           & ! [IN]
                          MOMX,           & ! [IN]
                          MOMY,           & ! [IN]
                          RHOT,           & ! [IN]
                          QTRC,           & ! [IN]
                          w0avg,          & ! [IN]
                          DENS_t_CP,      & ! [INOUT]
                          MOMZ_t_CP,      & ! [INOUT]
                          MOMX_t_CP,      & ! [INOUT]
                          MOMY_t_CP,      & ! [INOUT]
                          RHOT_t_CP,      & ! [INOUT]
                          RHOQ_t_CP,      & ! [INOUT]
                          MFLX_cloudbase, & ! [INOUT]
                          SFLX_rain,      & ! [OUT]
                          cloudtop,       & ! [OUT]
                          cloudbase,      & ! [OUT]
                          cldfrac_dp,     & ! [OUT]
                          cldfrac_sh,     & ! [OUT]
                          kf_nca          ) ! [OUT]

       ! tentative reset
!OCL XFILL
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          MOMZ_t_CP(k,i,j) = 0.0_RP
          MOMX_t_CP(k,i,j) = 0.0_RP
          MOMY_t_CP(k,i,j) = 0.0_RP
       enddo
       enddo
       enddo

!OCL XFILL
       do j  = JS, JE
       do i  = IS, IE
          MFLX_cloudbase(i,j) = 0.0_RP
       enddo
       enddo

       ! diagnose tendency of number concentration
       call ATMOS_HYDROMETEOR_diagnose_number_concentration( RHOQ_t_CP(:,:,:,:) ) ! [INOUT]

       call FILE_HISTORY_in( MFLX_cloudbase(:,:),   'CBMFX',     'cloud base mass flux',             'kg/m2/s', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rain     (:,:),   'RAIN_CP',   'surface rain rate by CP',          'kg/m2/s', fill_halo=.true. )
       call FILE_HISTORY_in( SFLX_rain     (:,:),   'PREC_CP',   'surface precipitation rate by CP', 'kg/m2/s', fill_halo=.true. )
       call FILE_HISTORY_in( cloudtop      (:,:),   'CUMHGT',    'CP cloud top height',              'm',       fill_halo=.true. )
       call FILE_HISTORY_in( cloudbase     (:,:),   'CUBASE',    'CP cloud base height',             'm',       fill_halo=.true. )
       call FILE_HISTORY_in( cldfrac_dp    (:,:,:), 'CUMFRC_DP', 'CP cloud fraction (deep)',         '1',       fill_halo=.true. )
       call FILE_HISTORY_in( cldfrac_sh    (:,:,:), 'CUMFRC_SH', 'CP cloud fraction (shallow)',      '1',       fill_halo=.true. )
       call FILE_HISTORY_in( kf_nca        (:,:),   'kf_nca',    'advection or cumulus convection timescale for KF', 's',       fill_halo=.true. )

       call FILE_HISTORY_in( DENS_t_CP(:,:,:), 'DENS_t_CP', 'tendency DENS in CP', 'kg/m3/s'  , fill_halo=.true. )
       call FILE_HISTORY_in( MOMZ_t_CP(:,:,:), 'MOMZ_t_CP', 'tendency MOMZ in CP', 'kg/m2/s2' , fill_halo=.true. )
       call FILE_HISTORY_in( MOMX_t_CP(:,:,:), 'MOMX_t_CP', 'tendency MOMX in CP', 'kg/m2/s2' , fill_halo=.true. )
       call FILE_HISTORY_in( MOMY_t_CP(:,:,:), 'MOMY_t_CP', 'tendency MOMY in CP', 'kg/m2/s2' , fill_halo=.true. )
       call FILE_HISTORY_in( RHOT_t_CP(:,:,:), 'RHOT_t_CP', 'tendency RHOT in CP', 'K*kg/m3/s', fill_halo=.true. )

       do iq = QS_MP, QE_MP
          call FILE_HISTORY_in( RHOQ_t_CP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_CP', &
                        'tendency rho*'//trim(TRACER_NAME(iq))//' in CP', 'kg/m3/s', fill_halo=.true. )
       enddo

    endif ! update

    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       DENS_t(k,i,j) = DENS_t(k,i,j) + DENS_t_CP(k,i,j)
       MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_CP(k,i,j)
       MOMX_t(k,i,j) = MOMX_t(k,i,j) + MOMX_t_CP(k,i,j)
       MOMY_t(k,i,j) = MOMY_t(k,i,j) + MOMY_t_CP(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_CP(k,i,j)
    enddo
    enddo
    enddo

    do iq = QS_MP,  QE_MP
    !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(3)
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_CP(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    if ( STATISTICS_checktotal ) then
       call STAT_total( total, DENS_t_CP(:,:,:), 'DENS_t_CP' )
       call STAT_total( total, MOMZ_t_CP(:,:,:), 'MOMZ_t_CP' )
       call STAT_total( total, MOMX_t_CP(:,:,:), 'MOMX_t_CP' )
       call STAT_total( total, MOMY_t_CP(:,:,:), 'MOMY_t_CP' )
       call STAT_total( total, RHOT_t_CP(:,:,:), 'RHOT_t_CP' )

       do iq = QS_MP, QE_MP
          call STAT_total( total, RHOQ_t_CP(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_CP' )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_CP_driver

end module mod_atmos_phy_cp_driver
