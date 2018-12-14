!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_tb_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
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
  public :: ATMOS_PHY_TB_driver_tracer_setup
  public :: ATMOS_PHY_TB_driver_setup
  public :: ATMOS_PHY_TB_driver_calc_tendency

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
  integer, private :: monit_west
  integer, private :: monit_east
  integer, private :: monit_south
  integer, private :: monit_north

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine ATMOS_PHY_TB_driver_tracer_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_tb_d1980, only: &
       ATMOS_PHY_TB_d1980_config
    use scale_atmos_phy_tb_dns, only: &
       ATMOS_PHY_TB_dns_config
    use mod_atmos_admin, only: &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_sw_phy_tb
    use mod_atmos_phy_tb_vars, only: &
       I_TKE
    implicit none
    !---------------------------------------------------------------------------

    if ( .not. ATMOS_sw_phy_tb ) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_driver_tracer_setup",*) 'Setup'

    select case( ATMOS_PHY_TB_TYPE )
    case( 'SMAGORINSKY' )
       I_TKE = -1
    case( 'D1980' )
       call ATMOS_PHY_TB_d1980_config( ATMOS_PHY_TB_TYPE, & ! [IN]
                                       I_TKE              ) ! [OUT]
    case( 'DNS' )
       call ATMOS_PHY_TB_dns_config( ATMOS_PHY_TB_TYPE, & ! [IN]
                                     I_TKE              ) ! [OUT]
    case('OFF')
       ! do nothing
    case default
       LOG_ERROR("ATMOS_PHY_TB_driver_tracer_setup",*) 'ATMOS_PHY_TB_TYPE is invalid: ', ATMOS_PHY_TB_TYPE
       call PRC_abort
    end select

    return
  end subroutine ATMOS_PHY_TB_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_TB_driver_setup
    use scale_atmos_grid_cartesC, only: &
       CDZ => ATMOS_GRID_CARTESC_CDZ, &
       CDX => ATMOS_GRID_CARTESC_CDX, &
       CDY => ATMOS_GRID_CARTESC_CDY
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ
    use scale_atmos_grid_cartesC_metric, only: &
       MAPF  => ATMOS_GRID_CARTESC_METRIC_MAPF
    use scale_atmos_phy_tb_smg, only: &
       ATMOS_PHY_TB_smg_setup
    use scale_atmos_phy_tb_d1980, only: &
       ATMOS_PHY_TB_d1980_setup
    use scale_atmos_phy_tb_dns, only: &
       ATMOS_PHY_TB_dns_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_TB_TYPE, &
       ATMOS_sw_phy_tb
    use mod_atmos_phy_tb_vars, only: &
       MOMZ_t_TB => ATMOS_PHY_TB_MOMZ_t
    use scale_monitor, only: &
       MONITOR_reg
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    if ( .not. ATMOS_sw_phy_tb ) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_driver_setup",*) 'Setup'

    ! initialize
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       MOMZ_t_TB(KS-1,i,j) = 0.0_RP
       MOMZ_t_TB(KE  ,i,j) = 0.0_RP
    enddo
    enddo

    ! setup library component
    select case( ATMOS_PHY_TB_TYPE )
    case( 'SMAGORINSKY' )
       call ATMOS_PHY_TB_smg_setup( REAL_FZ, REAL_CZ, CDX, CDY, MAPF(:,:,:,I_XY) ) ! [IN]
    case( 'D1980' )
       call ATMOS_PHY_TB_d1980_setup( CDZ, CDX, CDY, REAL_CZ ) ! [IN]
    case( 'DNS' )
       call ATMOS_PHY_TB_dns_setup( CDZ, CDX, CDY, REAL_CZ ) ! [IN]
    end select


    ! monitor
    call MONITOR_reg( "QTOTFLX_TB_WEST",  "water mass flux at the western boundary",  "kg/s", & ! [IN]
                      monit_west,                                                             & ! [OUT]
                      dim_type="ZY-W", is_tendency=.true.                                     ) ! [IN]
    call MONITOR_reg( "QTOTFLX_TB_EAST",  "water mass flux at the eastern boundary",  "kg/s", & ! [IN]
                      monit_east,                                                             & ! [OUT]
                      dim_type="ZY-E", is_tendency=.true.                                     ) ! [IN]
    call MONITOR_reg( "QTOTFLX_TB_SOUTH", "water mass flux at the southern boundary", "kg/s", & ! [IN]
                      monit_south,                                                            & ! [OUT]
                      dim_type="ZX-S", is_tendency=.true.                                     ) ! [IN]
    call MONITOR_reg( "QTOTFLX_TB_NORTH", "water mass flux at the northern boundary", "kg/s", & ! [IN]
                      monit_north,                                                            & ! [OUT]
                      dim_type="ZX-N", is_tendency=.true.                                     ) ! [IN]


    return
  end subroutine ATMOS_PHY_TB_driver_setup

  !-----------------------------------------------------------------------------
  !> calclate tendency
  subroutine ATMOS_PHY_TB_driver_calc_tendency( update_flag )
    use scale_prc_cartesC, only: &
       PRC_HAS_W, &
       PRC_HAS_E, &
       PRC_HAS_S, &
       PRC_HAS_N
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC, only: &
       FDZ  => ATMOS_GRID_CARTESC_FDZ, &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
       CDX  => ATMOS_GRID_CARTESC_CDX, &
       FDX  => ATMOS_GRID_CARTESC_FDX, &
       CDY  => ATMOS_GRID_CARTESC_CDY, &
       FDY  => ATMOS_GRID_CARTESC_FDY
    use scale_atmos_grid_cartesC_real, only: &
       FZ  => ATMOS_GRID_CARTESC_REAL_FZ, &
       ATMOS_GRID_CARTESC_REAL_VOL,       &
       ATMOS_GRID_CARTESC_REAL_TOTVOL,    &
       ATMOS_GRID_CARTESC_REAL_VOLWXY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLWXY, &
       ATMOS_GRID_CARTESC_REAL_VOLZUY,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZUY, &
       ATMOS_GRID_CARTESC_REAL_VOLZXV,    &
       ATMOS_GRID_CARTESC_REAL_TOTVOLZXV
    use scale_atmos_grid_cartesC_metric, only: &
       GSQRT => ATMOS_GRID_CARTESC_METRIC_GSQRT, &
       J13G  => ATMOS_GRID_CARTESC_METRIC_J13G,  &
       J23G  => ATMOS_GRID_CARTESC_METRIC_J23G,  &
       J33G  => ATMOS_GRID_CARTESC_METRIC_J33G,  &
       MAPF  => ATMOS_GRID_CARTESC_METRIC_MAPF
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_monitor, only: &
       MONITOR_put
    use scale_time, only: &
       dt_TB => TIME_DTSEC_ATMOS_PHY_TB
    use scale_atmos_phy_tb_smg, only: &
       ATMOS_PHY_TB_smg
    use scale_atmos_phy_tb_d1980, only: &
       ATMOS_PHY_TB_d1980
    use scale_atmos_phy_tb_dns, only: &
       ATMOS_PHY_TB_dns
    use scale_atmos_phy_tb_common, only: &
       calc_tend_momz => ATMOS_PHY_TB_calc_tend_momz, &
       calc_tend_momx => ATMOS_PHY_TB_calc_tend_momx, &
       calc_tend_momy => ATMOS_PHY_TB_calc_tend_momy, &
       calc_tend_phi  => ATMOS_PHY_TB_calc_tend_phi
    use mod_atmos_admin, only: &
       ATMOS_PHY_TB_TYPE
    use mod_atmos_vars, only: &
       ATMOS_vars_get_diagnostic, &
       DENS => DENS_av,   &
       MOMZ => MOMZ_av,   &
       MOMX => MOMX_av,   &
       MOMY => MOMY_av,   &
       RHOT => RHOT_av,   &
       POTT,              &
       QTRC => QTRC_av,   &
       MOMZ_t => MOMZ_tp, &
       MOMX_t => MOMX_tp, &
       MOMY_t => MOMY_tp, &
       RHOT_t => RHOT_tp, &
       RHOQ_t => RHOQ_tp
    use mod_atmos_phy_tb_vars, only: &
       I_TKE, &
       MOMZ_t_TB => ATMOS_PHY_TB_MOMZ_t, &
       MOMX_t_TB => ATMOS_PHY_TB_MOMX_t, &
       MOMY_t_TB => ATMOS_PHY_TB_MOMY_t, &
       RHOT_t_TB => ATMOS_PHY_TB_RHOT_t, &
       RHOQ_t_TB => ATMOS_PHY_TB_RHOQ_t
    use mod_atmos_phy_sf_vars, only: &
       SFLX_MW => ATMOS_PHY_SF_SFLX_MW, &
       SFLX_MU => ATMOS_PHY_SF_SFLX_MU, &
       SFLX_MV => ATMOS_PHY_SF_SFLX_MV, &
       SFLX_SH => ATMOS_PHY_SF_SFLX_SH, &
       SFLX_Q  => ATMOS_PHY_SF_SFLX_QTRC
    implicit none

    logical, intent(in) :: update_flag

    ! eddy viscosity/diffusion flux
    real(RP) :: QFLX_MOMZ(KA,IA,JA,3)
    real(RP) :: QFLX_MOMX(KA,IA,JA,3)
    real(RP) :: QFLX_MOMY(KA,IA,JA,3)
    real(RP) :: QFLX_RHOT(KA,IA,JA,3)
    real(RP) :: QFLX_RHOQ(KA,IA,JA,3,QA)

    real(RP) :: Nu(KA,IA,JA) ! eddy viscosity
    real(RP) :: Ri(KA,IA,JA) ! Richardson number
    real(RP) :: Pr(KA,IA,JA) ! Prandtl number

    real(RP) :: N2(KA,IA,JA)

    real(RP) :: tend(KA,IA,JA)

    ! for monitor
    real(RP) :: qflx_x(KA,JA)
    real(RP) :: qflx_y(KA,IA)

    integer  :: JJS, JJE
    integer  :: IIS, IIE

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       call ATMOS_vars_get_diagnostic( "N2", N2 )

       select case( ATMOS_PHY_TB_TYPE )
       case( 'SMAGORINSKY' )
          call ATMOS_PHY_TB_smg( QFLX_MOMZ, QFLX_MOMX, QFLX_MOMY,         & ! [OUT]
                                 QFLX_RHOT, QFLX_RHOQ,                    & ! [OUT]
                                 MOMZ_t_TB, MOMX_t_TB, MOMY_t_TB,         & ! [OUT]
                                 RHOT_t_TB, RHOQ_t_TB,                    & ! [OUT]
                                 Nu, Ri, Pr,                              & ! [OUT]
                                 MOMZ, MOMX, MOMY, POTT, DENS, QTRC, N2,  & ! [IN]
                                 FZ, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY, & ! [IN]
                                 GSQRT, J13G, J23G, J33G, MAPF,           & ! [IN]
                                 dt_TB                                    ) ! [IN]
       case( 'D1980' )
          MOMZ_t_TB(:,:,:)   = 0.0_RP
          MOMX_t_TB(:,:,:)   = 0.0_RP
          MOMY_t_TB(:,:,:)   = 0.0_RP
          RHOT_t_TB(:,:,:)   = 0.0_RP
          call ATMOS_PHY_TB_d1980( QFLX_MOMZ, QFLX_MOMX, QFLX_MOMY,        & ! [OUT]
                                   QFLX_RHOT, QFLX_RHOQ,                   & ! [OUT]
                                   RHOQ_t_TB,                              & ! [OUT]
                                   Nu, Ri, Pr,                             & ! [OUT]
                                   MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2, & ! [IN]
                                   SFLX_MW, SFLX_MU, SFLX_MV,              & ! [IN]
                                   SFLX_SH, SFLX_Q,                        & ! [IN]
                                   GSQRT, J13G, J23G, J33G, MAPF,          & ! [IN]
                                   dt_TB                                   ) ! [IN]
       case( 'DNS' )
          MOMZ_t_TB(:,:,:)   = 0.0_RP
          MOMX_t_TB(:,:,:)   = 0.0_RP
          MOMY_t_TB(:,:,:)   = 0.0_RP
          RHOT_t_TB(:,:,:)   = 0.0_RP
          call ATMOS_PHY_TB_dns( QFLX_MOMZ, QFLX_MOMX, QFLX_MOMY,        & ! [OUT]
                                 QFLX_RHOT, QFLX_RHOQ,                   & ! [OUT]
                                 RHOQ_t_TB,                              & ! [OUT]
                                 Nu, Ri, Pr,                             & ! [OUT]
                                 MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2, & ! [IN]
                                 SFLX_MW, SFLX_MU, SFLX_MV,              & ! [IN]
                                 SFLX_SH, SFLX_Q,                        & ! [IN]
                                 GSQRT, J13G, J23G, J33G, MAPF,          & ! [IN]
                                 dt_TB                                   ) ! [IN]
       end select


       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          call calc_tend_momz( tend(:,:,:), & ! (out)
                               QFLX_MOMZ,   & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
          !$omp parallel do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMZ_t_TB(k,i,j) = MOMZ_t_TB(k,i,j) + tend(k,i,j)
          end do
          end do
          end do

          call calc_tend_momx( tend(:,:,:), & ! (out)
                               QFLX_MOMX,   & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
          !$omp parallel do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             MOMX_t_TB(k,i,j) = MOMX_t_TB(k,i,j) + tend(k,i,j)
          end do
          end do
          end do

          call calc_tend_momy( tend(:,:,:), & ! (out)
                               QFLX_MOMY,   & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
          !$omp parallel do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             MOMY_t_TB(k,i,j) = MOMY_t_TB(k,i,j) + tend(k,i,j)
          end do
          end do
          end do

          call calc_tend_phi ( tend(:,:,:), & ! (out)
                               QFLX_RHOT,   & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)
          !$omp parallel do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             RHOT_t_TB(k,i,j) = RHOT_t_TB(k,i,j) + tend(k,i,j)
          end do
          end do
          end do

          do iq = 1, QA
             if ( iq == I_TKE .or. .not. TRACER_ADVC(iq) ) cycle

             call calc_tend_phi( tend(:,:,:),                   & ! (out)
                                 QFLX_RHOQ(:,:,:,:,iq),         & ! (in)
                                 GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                                 IIS, IIE, JJS, JJE ) ! (in)

             !$omp parallel do
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                RHOQ_t_TB(k,i,j,iq) = RHOQ_t_TB(k,i,j,iq) + tend(k,i,j)
             end do
             end do
             end do

          end do

       end do
       end do


       call FILE_HISTORY_in( NU (:,:,:), 'NU',  'eddy viscosity',           'm2/s' , fill_halo=.true. )
       call FILE_HISTORY_in( Ri (:,:,:), 'Ri',  'Richardson number',        'NIL'  , fill_halo=.true. )
       call FILE_HISTORY_in( Pr (:,:,:), 'Pr',  'Prantle number',           'NIL'  , fill_halo=.true. )

       call FILE_HISTORY_in( MOMZ_t_TB(:,:,:), 'MOMZ_t_TB', 'MOMZ tendency (TB)', 'kg/m2/s2',  fill_halo=.true. )
       call FILE_HISTORY_in( MOMX_t_TB(:,:,:), 'MOMX_t_TB', 'MOMX tendency (TB)', 'kg/m2/s2',  fill_halo=.true. )
       call FILE_HISTORY_in( MOMY_t_TB(:,:,:), 'MOMY_t_TB', 'MOMY tendency (TB)', 'kg/m2/s2',  fill_halo=.true. )
       call FILE_HISTORY_in( RHOT_t_TB(:,:,:), 'RHOT_t_TB', 'RHOT tendency (TB)', 'K.kg/m3/s', fill_halo=.true. )

       do iq = 1, QA
          call FILE_HISTORY_in( RHOQ_t_TB(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_TB',                      &
                        'tendency rho*'//trim(TRACER_NAME(iq))//' in TB', 'kg/m3/s', fill_halo=.true. )
       enddo

       call FILE_HISTORY_in( QFLX_MOMZ(:,:,:,ZDIR), 'SGS_ZFLX_MOMZ', 'SGS Z FLUX of MOMZ', 'kg/m/s2', &
                     fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_MOMZ(:,:,:,XDIR), 'SGS_XFLX_MOMZ', 'SGS X FLUX of MOMZ', 'kg/m/s2', &
                     dim_type='ZHXHY', fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_MOMZ(:,:,:,YDIR), 'SGS_YFLX_MOMZ', 'SGS Y FLUX of MOMZ', 'kg/m/s2', &
                     dim_type='ZHXYH', fill_halo=.true.)

       call FILE_HISTORY_in( QFLX_MOMX(:,:,:,ZDIR), 'SGS_ZFLX_MOMX', 'SGS Z FLUX of MOMX', 'kg/m/s2', &
                     dim_type='ZHXHY', fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_MOMX(:,:,:,XDIR), 'SGS_XFLX_MOMX', 'SGS X FLUX of MOMX', 'kg/m/s2', &
                     fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_MOMX(:,:,:,YDIR), 'SGS_YFLX_MOMX', 'SGS Y FLUX of MOMX', 'kg/m/s2', &
                     dim_type='ZXHYH', fill_halo=.true.)

       call FILE_HISTORY_in( QFLX_MOMY(:,:,:,ZDIR), 'SGS_ZFLX_MOMY', 'SGS Z FLUX of MOMY', 'kg/m/s2', &
                     dim_type='ZHXYH', fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_MOMY(:,:,:,XDIR), 'SGS_XFLX_MOMY', 'SGS X FLUX of MOMY', 'kg/m/s2', &
                     dim_type='ZXHYH', fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_MOMY(:,:,:,YDIR), 'SGS_YFLX_MOMY', 'SGS Y FLUX of MOMY', 'kg/m/s2', &
                     fill_halo=.true.)

       call FILE_HISTORY_in( QFLX_RHOT(:,:,:,ZDIR), 'SGS_ZFLX_RHOT', 'SGS Z FLUX of RHOT', 'K*kg/m2/s', &
                     dim_type='ZHXY', fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_RHOT(:,:,:,XDIR), 'SGS_XFLX_RHOT', 'SGS X FLUX of RHOT', 'K*kg/m2/s', &
                     dim_type='ZXHY', fill_halo=.true.)
       call FILE_HISTORY_in( QFLX_RHOT(:,:,:,YDIR), 'SGS_YFLX_RHOT', 'SGS Y FLUX of RHOT', 'K*kg/m2/s', &
                     dim_type='ZXYH', fill_halo=.true.)


       do iq = 1, QA
          if ( iq == I_TKE .or. .not. TRACER_ADVC(iq) ) cycle

          call FILE_HISTORY_in( QFLX_RHOQ(:,:,:,ZDIR,iq), &
               'SGS_ZFLX_'//trim(TRACER_NAME(iq)), 'SGS Z FLUX of '//trim(TRACER_NAME(iq)), 'kg/m2/s', &
               dim_type='ZHXY', fill_halo=.true.)
          call FILE_HISTORY_in( QFLX_RHOQ(:,:,:,XDIR,iq), &
               'SGS_XFLX_'//trim(TRACER_NAME(iq)), 'SGS X FLUX of '//trim(TRACER_NAME(iq)), 'kg/m2/s', &
               dim_type='ZXHY', fill_halo=.true.)
          call FILE_HISTORY_in( QFLX_RHOQ(:,:,:,YDIR,iq), &
               'SGS_YFLX_'//trim(TRACER_NAME(iq)), 'SGS Y FLUX of '//trim(TRACER_NAME(iq)), 'kg/m2/s', &
               dim_type='ZXYH', fill_halo=.true.)
       end do

       if ( STATISTICS_checktotal ) then
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 MOMZ_t_TB(:,:,:), 'MOMZ_t_TB',         &
                                 ATMOS_GRID_CARTESC_REAL_VOLWXY(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOLWXY      )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 MOMX_t_TB(:,:,:), 'MOMX_t_TB',         &
                                 ATMOS_GRID_CARTESC_REAL_VOLZUY(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOLZUY      )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 MOMY_t_TB(:,:,:), 'MOMY_t_TB',         &
                                 ATMOS_GRID_CARTESC_REAL_VOLZXV(:,:,:), &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOLZXV      )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOT_t_TB(:,:,:), 'RHOT_t_TB',       &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL       )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 Nu(:,:,:),        'Nu',              &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL       )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 Ri(:,:,:),        'Ri',              &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL       )
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 Pr(:,:,:),        'Pr',              &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL       )

          do iq = 1, QA
             call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                    RHOQ_t_TB(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_TB', &
                                    ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                    ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
          enddo
       endif

    endif

    !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,KS,KE,MOMZ_t,MOMZ_t_TB,MOMX_t,MOMX_t_TB,MOMY_t,MOMY_t_TB,RHOT_t,RHOT_t_TB)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMZ_t(k,i,j) = MOMZ_t(k,i,j) + MOMZ_t_TB(k,i,j)
       MOMX_t(k,i,j) = MOMX_t(k,i,j) + MOMX_t_TB(k,i,j)
       MOMY_t(k,i,j) = MOMY_t(k,i,j) + MOMY_t_TB(k,i,j)
       RHOT_t(k,i,j) = RHOT_t(k,i,j) + RHOT_t_TB(k,i,j)
    enddo
    enddo
    enddo

    do iq = 1, QA

       if ( .not. ( iq == I_TKE .or. TRACER_ADVC(iq) ) ) cycle

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_TB(k,i,j,iq)
       enddo
       enddo
      enddo

    enddo

    if ( monit_west > 0 ) then
       qflx_x(:,:) = 0.0_RP
       if ( .not. PRC_HAS_W ) then
          do iq = 1, QA
             if ( TRACER_ADVC(iq) .and. TRACER_MASS(iq) == 1.0_RP ) then
                do j = JS, JE
                do k = KS, KE
                   qflx_x(k,j) = qflx_x(k,j) + QFLX_RHOQ(k,IS-1,j,XDIR,iq)
                end do
                end do
             end if
          end do
       end if
       call MONITOR_put( monit_west, qflx_x(:,:) )
    end if
    if ( monit_east > 0 ) then
       qflx_x(:,:) = 0.0_RP
       if ( .not. PRC_HAS_E ) then
          do iq = 1, QA
             if ( TRACER_ADVC(iq) .and. TRACER_MASS(iq) == 1.0_RP ) then
                do j = JS, JE
                do k = KS, KE
                   qflx_x(k,j) = qflx_x(k,j) + QFLX_RHOQ(k,IE,j,XDIR,iq)
                end do
                end do
             end if
          end do
       end if
       call MONITOR_put( monit_east, qflx_x(:,:) )
    end if
    if ( monit_south > 0 ) then
       qflx_y(:,:) = 0.0_RP
       if ( .not. PRC_HAS_S ) then
          do iq = 1, QA
             if ( TRACER_ADVC(iq) .and. TRACER_MASS(iq) == 1.0_RP ) then
                do i = IS, IE
                do k = KS, KE
                   qflx_y(k,i) = qflx_y(k,i) + QFLX_RHOQ(k,i,JS-1,YDIR,iq)
                end do
                end do
             end if
          end do
       end if
       call MONITOR_put( monit_south, qflx_y(:,:) )
    end if
    if ( monit_north > 0 ) then
       qflx_y(:,:) = 0.0_RP
       if ( .not. PRC_HAS_N ) then
          do iq = 1, QA
             if ( TRACER_ADVC(iq) .and. TRACER_MASS(iq) == 1.0_RP ) then
                do i = IS, IE
                do k = KS, KE
                   qflx_y(k,i) = qflx_y(k,i) + QFLX_RHOQ(k,i,JE,YDIR,iq)
                end do
                end do
             end if
          end do
       end if
       call MONITOR_put( monit_north, qflx_y(:,:) )
    end if

    return
  end subroutine ATMOS_PHY_TB_driver_calc_tendency

end module mod_atmos_phy_tb_driver
