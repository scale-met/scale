!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          HEVE FVM scheme for large time step in Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tstep_large_fvm_heve
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer

#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tstep_large_fvm_heve_setup
  public :: ATMOS_DYN_Tstep_large_fvm_heve_finalize
  public :: ATMOS_DYN_Tstep_large_fvm_heve

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
  ! tendency
  real(RP), private, allocatable :: DENS_t(:,:,:)
  real(RP), private, allocatable :: MOMZ_t(:,:,:)
  real(RP), private, allocatable :: MOMX_t(:,:,:)
  real(RP), private, allocatable :: MOMY_t(:,:,:)
  real(RP), private, allocatable :: RHOT_t(:,:,:)
  real(RP), private, allocatable, target :: RHOQ_t(:,:,:,:)
  real(RP), private, allocatable, target :: ZERO(:,:,:)
  real(RP), private, pointer     :: RHOQ_tn(:,:,:)

  real(RP), private, allocatable :: DENS_damp(:,:,:)

  ! flux
  real(RP), private, allocatable, target :: mflx(:,:,:,:) ! rho * vel(x,y,z) * GSQRT / mapf

  ! for communication
  integer :: I_COMM_DENS
  integer :: I_COMM_MOMZ
  integer :: I_COMM_MOMX
  integer :: I_COMM_MOMY
  integer :: I_COMM_RHOT
  integer, allocatable :: I_COMM_PROG(:)

  integer :: I_COMM_DENS_t
  integer :: I_COMM_MOMZ_t
  integer :: I_COMM_MOMX_t
  integer :: I_COMM_MOMY_t
  integer :: I_COMM_RHOT_t

  integer :: I_COMM_DENS_damp

  integer, allocatable :: I_COMM_RHOQ_t(:)
  integer, allocatable :: I_COMM_QTRC(:)

  integer :: I_COMM_mflx_z
  integer :: I_COMM_mflx_x
  integer :: I_COMM_mflx_y

  ! for history
  integer :: HIST_mflx(3)
  integer :: HIST_tflx(3)
  integer :: HIST_phys(5)
  integer :: HIST_damp(5)
  integer, allocatable :: HIST_qflx(:,:)
  integer, allocatable :: HIST_phys_QTRC(:)
  integer, allocatable :: HIST_damp_QTRC(:)

  ! for monitor
  real(RP), allocatable, target :: zero_x(:,:)
  real(RP), allocatable, target :: zero_y(:,:)
  integer :: MONIT_damp_mass
  integer :: MONIT_damp_qtot
  integer :: MONIT_mflx_west
  integer :: MONIT_mflx_east
  integer :: MONIT_mflx_south
  integer :: MONIT_mflx_north
  integer :: MONIT_qflx_west
  integer :: MONIT_qflx_east
  integer :: MONIT_qflx_south
  integer :: MONIT_qflx_north

  character(len=H_SHORT) :: EVAL_TYPE_NUMFILTER = 'TENDENCY'

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_large_fvm_heve_setup( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD,  &
       PRC_HAS_E, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_const, only: &
       OHM => CONST_OHM, &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    use scale_file_history, only: &
       FILE_HISTORY_reg, &
       FILE_HISTORY_put
    use scale_monitor, only: &
       MONITOR_reg, &
       MONITOR_put
    implicit none

    ! MPI_RECV_INIT requires intent(inout)
    real(RP),               intent(inout) :: DENS(KA,IA,JA)
    real(RP),               intent(inout) :: MOMZ(KA,IA,JA)
    real(RP),               intent(inout) :: MOMX(KA,IA,JA)
    real(RP),               intent(inout) :: MOMY(KA,IA,JA)
    real(RP),               intent(inout) :: RHOT(KA,IA,JA)
    real(RP),               intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP),               intent(inout) :: PROG(KA,IA,JA,VA)

    integer :: iv, iq

    namelist /ATMOS_DYN_TSTEP_LARGE_FVM_HEVE/ &
      EVAL_TYPE_NUMFILTER
   
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=ATMOS_DYN_TSTEP_LARGE_FVM_HEVE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_DYN_Tstep_large_fvm_heve_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_DYN_Tstep_large_fvm_heve_setup",*) 'Not appropriate names in namelist ATMOS_DYN_TSTEP_LARGE_FVM_HEVE. Check!'
       call PRC_abort
    endif
    LOG_NML(ATMOS_DYN_TSTEP_LARGE_FVM_HEVE)

    select case( EVAL_TYPE_NUMFILTER ) 
    case( 'TENDENCY', 'FILTER' )
    case default
      LOG_ERROR("ATMOS_DYN_Tstep_large_fvm_heve_setup",*) 'The specfied value of EVAL_TYPE_NUMFILTER is not appropriate. Check!'
      call PRC_abort
    end select
    !--

    allocate( DENS_t(KA,IA,JA) )
    allocate( MOMZ_t(KA,IA,JA) )
    allocate( MOMX_t(KA,IA,JA) )
    allocate( MOMY_t(KA,IA,JA) )
    allocate( RHOT_t(KA,IA,JA) )
    allocate( RHOQ_t(KA,IA,JA,QA) )

    allocate( DENS_damp(KA,IA,JA) )

    allocate( mflx(KA,IA,JA,3) )

    allocate( I_COMM_PROG    (max(VA,1)) )
    allocate( I_COMM_QTRC(QA) )
    allocate( I_COMM_RHOQ_t(QA) )

    I_COMM_DENS = 1
    I_COMM_MOMZ = 2
    I_COMM_MOMX = 3
    I_COMM_MOMY = 4
    I_COMM_RHOT = 5
    call COMM_vars8_init( 'DENS', DENS, I_COMM_DENS )
    call COMM_vars8_init( 'MOMZ', MOMZ, I_COMM_MOMZ )
    call COMM_vars8_init( 'MOMX', MOMX, I_COMM_MOMX )
    call COMM_vars8_init( 'MOMY', MOMY, I_COMM_MOMY )
    call COMM_vars8_init( 'RHOT', RHOT, I_COMM_RHOT )
    do iv = 1, VA
       I_COMM_PROG(iv) = 5 + iv
       call COMM_vars8_init( 'PROG', PROG(:,:,:,iv), I_COMM_PROG(iv) )
    end do

    I_COMM_DENS_t = 1
    I_COMM_MOMZ_t = 2
    I_COMM_MOMX_t = 3
    I_COMM_MOMY_t = 4
    I_COMM_RHOT_t = 5
    I_COMM_DENS_damp = 6
    call COMM_vars8_init( 'DENS_t', DENS_t, I_COMM_DENS_t )
    call COMM_vars8_init( 'MOMZ_t', MOMZ_t, I_COMM_MOMZ_t )
    call COMM_vars8_init( 'MOMX_t', MOMX_t, I_COMM_MOMX_t )
    call COMM_vars8_init( 'MOMY_t', MOMY_t, I_COMM_MOMY_t )
    call COMM_vars8_init( 'RHOT_t', RHOT_t, I_COMM_RHOT_t )
    call COMM_vars8_init( 'DENS_t', DENS_t, I_COMM_DENS_damp )

    do iq = 1, QA
       I_COMM_RHOQ_t(iq) = 5 + VA + iq
       I_COMM_QTRC(iq) = 5 + VA + iq

       call COMM_vars8_init( 'RHOQ_t', RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq) )
       call COMM_vars8_init( 'QTRC',   QTRC  (:,:,:,iq), I_COMM_QTRC(iq) )
    end do

    I_COMM_mflx_z = 1
    I_COMM_mflx_x = 2
    I_COMM_mflx_y = 3
    call COMM_vars8_init( 'mflx_Z', mflx(:,:,:,ZDIR), I_COMM_mflx_z )
    call COMM_vars8_init( 'mflx_X', mflx(:,:,:,XDIR), I_COMM_mflx_x )
    call COMM_vars8_init( 'mflx_Y', mflx(:,:,:,YDIR), I_COMM_mflx_y )

    allocate( ZERO(KA,IA,JA) )
    ZERO(:,:,:) = 0.0_RP

    mflx(:,:,:,:) = UNDEF
    if ( PRC_TwoD ) then
       mflx(:,:,:,XDIR) = 0.0_RP
    end if


    ! history
    call FILE_HISTORY_reg( 'ZFLX_MOM', 'momentum flux of z-direction', 'kg/m2/s', & ! [IN]
                           HIST_mflx(1),                                          & ! [OUT]
                           dim_type='ZHXY'                                        ) ! [IN]
    call FILE_HISTORY_reg( 'XFLX_MOM', 'momentum flux of x-direction', 'kg/m2/s', & ! [IN]
                           HIST_mflx(2),                                          & ! [OUT]
                           dim_type='ZXHY'                                        ) ! [IN]
    call FILE_HISTORY_reg( 'YFLX_MOM', 'momentum flux of y-direction', 'kg/m2/s', & ! [IN]
                           HIST_mflx(3),                                          & ! [OUT]
                           dim_type='ZXYH'                                        ) ! [IN]

    call FILE_HISTORY_reg( 'ZFLX_RHOT', 'potential temperature flux of z-direction', 'K*kg/m2/s', & ! [IN]
                           HIST_tflx(1),                                                          & ! [OUT]
                           dim_type='ZHXY'                                                        ) ! [IN]
    call FILE_HISTORY_reg( 'XFLX_RHOT', 'potential temperature flux of x-direction', 'K*kg/m2/s', & ! [IN]
                           HIST_tflx(2),                                                          & ! [OUT]
                           dim_type='ZXHY'                                                        ) ! [IN]
    call FILE_HISTORY_reg( 'YFLX_RHOT', 'potential temperature flux of y-direction', 'K*kg/m2/s', & ! [IN]
                           HIST_tflx(3),                                                          & ! [OUT]
                           dim_type='ZXYH'                                                        ) ! [IN]

    call FILE_HISTORY_reg( 'DENS_t_phys', 'tendency of dencity due to physics', 'kg/m3/s', & ! [IN]
                           HIST_phys(1)                                                    ) ! [OUT]
    call FILE_HISTORY_reg( 'MOMZ_t_phys', 'tendency of momentum z due to physics', 'kg/m2/s2', & ! [IN]
                           HIST_phys(2),                                                       & ! [OUT]
                           dim_type='ZHXY'                                                     ) ! [IN]
    call FILE_HISTORY_reg( 'MOMX_t_phys', 'tendency of momentum x due to physics', 'kg/m2/s2', & ! [IN]
                           HIST_phys(3),                                                       & ! [OUT]
                           dim_type='ZXHY'                                                     ) ! [IN]
    call FILE_HISTORY_reg( 'MOMY_t_phys', 'tendency of momentum y due to physics', 'kg/m2/s2', & ! [IN]
                           HIST_phys(4),                                                       & ! [OUT]
                           dim_type='ZXYH'                                                     ) ! [IN]
    call FILE_HISTORY_reg( 'RHOT_t_phys', 'tendency of rho*theta temperature due to physics', 'K*kg/m3/s', & ! [IN]
                           HIST_phys(5)                                                                    ) ! [OUT]

    call FILE_HISTORY_reg( 'DENS_t_damp', 'tendency of dencity due to damping', 'kg/m3/s', & ! [IN]
                           HIST_damp(1)                                                    ) ! [OUT]
    call FILE_HISTORY_reg( 'MOMZ_t_damp', 'tendency of momentum z due to damping', 'kg/m2/s2', & ! [IN]
                           HIST_damp(2),                                                       & ! [OUT]
                           dim_type='ZHXY'                                                     ) ! [IN]
    call FILE_HISTORY_reg( 'MOMX_t_damp', 'tendency of momentum x due to damping', 'kg/m2/s2', & ! [IN]
                           HIST_damp(3),                                                       & ! [OUT]
                           dim_type='ZXHY'                                                     ) ! [IN]
    call FILE_HISTORY_reg( 'MOMY_t_damp', 'tendency of momentum y due to damping', 'kg/m2/s2', & ! [IN]
                           HIST_damp(4),                                                       & ! [OUT]
                           dim_type='ZXYH'                                                     ) ! [IN]
    call FILE_HISTORY_reg( 'RHOT_t_damp', 'tendency of rho*theta temperature due to damping', 'K kg/m3/s', & ! [IN]
                           HIST_damp(5)                                                                    ) ! [OUT]

    allocate( HIST_qflx(3,QA) )
    allocate( HIST_phys_QTRC(QA) )
    allocate( HIST_damp_QTRC(QA) )
    do iq = 1, QA
       call FILE_HISTORY_reg( 'ZFLX_'//trim(TRACER_NAME(iq)), trim(TRACER_NAME(iq))//' flux of z-direction', 'kg/m2/s', & ! [IN]
                              HIST_qflx(1,iq),                                                                          & ! [OUT]
                              dim_type='ZHXY'                                                                           ) ! [IN]
       call FILE_HISTORY_reg( 'XFLX_'//trim(TRACER_NAME(iq)), trim(TRACER_NAME(iq))//' flux of x-direction', 'kg/m2/s', & ! [IN]
                              HIST_qflx(2,iq),                                                                          & ! [OUT]
                              dim_type='ZXHY'                                                                           ) ! [IN]
       call FILE_HISTORY_reg( 'YFLX_'//trim(TRACER_NAME(iq)), trim(TRACER_NAME(iq))//' flux of y-direction', 'kg/m2/s', & ! [IN]
                              HIST_qflx(3,iq),                                                                          & ! [OUT]
                              dim_type='ZXYH'                                                                           ) ! [IN]
       call FILE_HISTORY_reg( trim(TRACER_NAME(iq))//'_t_phys', 'tendency of '//trim(TRACER_NAME(iq))//' due to physics', 'kg/m3/s', &
                              HIST_phys_QTRC(iq) )
       call FILE_HISTORY_reg( trim(TRACER_NAME(iq))//'_t_damp', 'tendency of '//trim(TRACER_NAME(iq))//' due to damping', 'kg/m3/s', &
                              HIST_damp_QTRC(iq) )
    end do

    ! for history at t=0
    do iv = 1, 3
       call FILE_HISTORY_put( HIST_mflx(iv), ZERO(:,:,:) )
       call FILE_HISTORY_put( HIST_tflx(iv), ZERO(:,:,:) )
    end do
    do iv = 1, 5
       call FILE_HISTORY_put( HIST_phys(iv), ZERO(:,:,:) )
       call FILE_HISTORY_put( HIST_damp(iv), ZERO(:,:,:) )
    end do
    do iq = 1, QA
       do iv = 1, 3
          call FILE_HISTORY_put( HIST_qflx(iv,iq), ZERO(:,:,:) )
       end do
       call FILE_HISTORY_put( HIST_phys_QTRC(iq), ZERO(:,:,:) )
       call FILE_HISTORY_put( HIST_damp_QTRC(iq), ZERO(:,:,:) )
    end do


    ! for monitor
    allocate( zero_x(KA,JA) )
    allocate( zero_y(KA,IA) )
    zero_x(:,:) = 0.0_RP
    zero_y(:,:) = 0.0_RP

    call MONITOR_reg( "MASSTND_DAMP", "mass tendency by the damping", "kg", & ! [IN]
                      MONIT_damp_mass,                                      & ! [OUT]
                      is_tendency=.true.                                    ) ! [IN]

    call MONITOR_reg( "MASSFLX_WEST",  "mass flux at the western boundary",  "kg", & ! [IN]
                      MONIT_mflx_west,                                             & ! [OUT]
                      dim_type="ZY-W", is_tendency=.true.                          ) ! [IN]
    call MONITOR_reg( "MASSFLX_EAST",  "mass flux at the eastern boundary",  "kg", & ! [IN]
                      MONIT_mflx_east,                                             & ! [OUT]
                      dim_type="ZY-E", is_tendency=.true.                          ) ! [IN]
    call MONITOR_reg( "MASSFLX_SOUTH", "mass flux at the southern boundary", "kg", & ! [IN]
                      MONIT_mflx_south,                                            & ! [OUT]
                      dim_type="ZX-S", is_tendency=.true.                          ) ! [IN]
    call MONITOR_reg( "MASSFLX_NORTH", "mass flux at the northern boundary", "kg", & ! [IN]
                      MONIT_mflx_north,                                            & ! [OUT]
                      dim_type="ZX-N", is_tendency=.true.                          ) ! [IN]

    call MONITOR_reg( "QTOTTND_DAMP", "water mass tendency by the damping", "kg", & ! [IN]
                      MONIT_damp_qtot,                                            & ! [OUT]
                      is_tendency=.true.                                          ) ! [IN]

    call MONITOR_reg( "QTOTFLX_WEST",  "water mass flux at the western boundary",  "kg", & ! [IN]
                      MONIT_qflx_west,                                                   & ! [OUT]
                      dim_type="ZY-W", is_tendency=.true.                                ) ! [IN]
    call MONITOR_reg( "QTOTFLX_EAST",  "water mass flux at the eastern boundary",  "kg", & ! [IN]
                      MONIT_qflx_east,                                                   & ! [OUT]
                      dim_type="ZY-E", is_tendency=.true.                                ) ! [IN]
    call MONITOR_reg( "QTOTFLX_SOUTH", "water mass flux at the southern boundary", "kg", & ! [IN]
                      MONIT_qflx_south,                                                  & ! [OUT]
                      dim_type="ZX-S", is_tendency=.true.                                ) ! [IN]
    call MONITOR_reg( "QTOTFLX_NORTH", "water mass flux at the northern boundary", "kg", & ! [IN]
                      MONIT_qflx_north,                                                  & ! [OUT]
                      dim_type="ZX-N", is_tendency=.true.                                ) ! [IN]

    ! at t=0
    call MONITOR_put( MONIT_damp_mass,  ZERO(:,:,:) )
    call MONITOR_put( MONIT_mflx_west,  zero_x(:,:) )
    call MONITOR_put( MONIT_mflx_east,  zero_x(:,:) )
    call MONITOR_put( MONIT_mflx_south, zero_y(:,:) )
    call MONITOR_put( MONIT_mflx_north, zero_y(:,:) )

    call MONITOR_put( MONIT_damp_qtot,  ZERO(:,:,:) )
    call MONITOR_put( MONIT_qflx_west,  zero_x(:,:) )
    call MONITOR_put( MONIT_qflx_east,  zero_x(:,:) )
    call MONITOR_put( MONIT_qflx_south, zero_y(:,:) )
    call MONITOR_put( MONIT_qflx_north, zero_y(:,:) )

    return
  end subroutine ATMOS_DYN_Tstep_large_fvm_heve_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_DYN_Tstep_large_fvm_heve_finalize

    deallocate( DENS_t )
    deallocate( MOMZ_t )
    deallocate( MOMX_t )
    deallocate( MOMY_t )
    deallocate( RHOT_t )
    deallocate( RHOQ_t )

    deallocate( DENS_damp )

    deallocate( mflx )

    deallocate( I_COMM_PROG )
    deallocate( I_COMM_QTRC )
    deallocate( I_COMM_RHOQ_t )

    deallocate( ZERO )

    deallocate( HIST_qflx )
    deallocate( HIST_phys_QTRC )
    deallocate( HIST_damp_QTRC )

    deallocate( zero_x )
    deallocate( zero_y )

    return
  end subroutine ATMOS_DYN_Tstep_large_fvm_heve_finalize

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  subroutine ATMOS_DYN_Tstep_large_fvm_heve( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG,             &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       num_diff, num_diff_q,                                 &
       QTRC0,                                                &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, &
       CORIOLI,                                              &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       PHI, GSQRT,                                           &
       J13G, J23G, J33G, MAPF,                               &
       AQ_R, AQ_CV, AQ_CP, AQ_MASS,                          &
       REF_dens, REF_pott, REF_qv, REF_pres,                 &
       BND_W, BND_E, BND_S, BND_N, TwoD,                     &
       ND_COEF, ND_COEF_Q, ND_LAPLACIAN_NUM,                 &
       ND_SFC_FACT, ND_USE_RS,                               &
       BND_QA, BND_IQ, BND_SMOOTHER_FACT,                    &
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          &
       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    &
       DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    &
       MFLUX_OFFSET_X, MFLUX_OFFSET_Y,                       &
       wdamp_coef, divdmp_coef,                              &
       FLAG_TRACER_SPLIT_TEND,                               &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       &
       FLAG_FCT_ALONG_STREAM,                                &
       USE_AVERAGE,                                          &
       I_QV,                                                 &
       DTLS, DTSS, Llast                                     )
    use scale_const, only: &
       EPS    => CONST_EPS
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_fill_halo,               &
       ATMOS_DYN_prep_pres_linearization, &
       ATMOS_DYN_divergence
    use scale_atmos_dyn_fvm_numfilter, only: &
       ATMOS_DYN_FVM_numfilter_flux,   &
       ATMOS_DYN_FVM_numfilter_flux_q, &
       ATMOS_DYN_fvm_apply_numfilter
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_put,   &
       FILE_HISTORY_set_disable
    use scale_monitor, only: &
       MONITOR_put
    use scale_atmos_dyn_tinteg_short, only: &
       ATMOS_DYN_tinteg_short
    use scale_atmos_dyn_tinteg_tracer, only: &
       ATMOS_DYN_tinteg_tracer
    use scale_spnudge, only: &
       SPNUDGE_uv,         &
       SPNUDGE_uv_divfree, &
       SPNUDGE_pt,         &
       SPNUDGE_qv,         &
       SPNUDGE_u_alpha,    &
       SPNUDGE_v_alpha,    &
       SPNUDGE_pt_alpha,   &
       SPNUDGE_qv_alpha,   &
       SPNUDGE_uv_lm,      &
       SPNUDGE_uv_mm,      &
       SPNUDGE_pt_lm,      &
       SPNUDGE_pt_mm,      &
       SPNUDGE_qv_lm,      &
       SPNUDGE_qv_mm
    use scale_dft, only: &
       DFT_g2g_divfree, &
       DFT_g2g
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

    real(RP), intent(inout) :: DENS_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_av(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_av(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_av(KA,IA,JA,QA)

    real(RP), intent(out)   :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(out)   :: num_diff_q(KA,IA,JA,3)

    real(RP), intent(in)    :: QTRC0(KA,IA,JA,QA)

    real(RP), intent(in)    :: DENS_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMX_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMY_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOT_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOQ_tp(KA,IA,JA,QA)

    real(RP), intent(in)    :: CORIOLI(IA,JA)

    real(RP), intent(in)    :: CDZ (KA)
    real(RP), intent(in)    :: CDX (IA)
    real(RP), intent(in)    :: CDY (JA)
    real(RP), intent(in)    :: FDZ (KA-1)
    real(RP), intent(in)    :: FDX (IA-1)
    real(RP), intent(in)    :: FDY (JA-1)
    real(RP), intent(in)    :: RCDZ(KA)
    real(RP), intent(in)    :: RCDX(IA)
    real(RP), intent(in)    :: RCDY(JA)
    real(RP), intent(in)    :: RFDZ(KA-1)
    real(RP), intent(in)    :: RFDX(IA-1)
    real(RP), intent(in)    :: RFDY(JA-1)

    real(RP), intent(in)    :: PHI  (KA,IA,JA)   !< geopotential
    real(RP), intent(in)    :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)    :: J13G (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)    :: J33G              !< (3,3) element of Jacobian matrix
    real(RP), intent(in)    :: MAPF (IA,JA,2,4)  !< map factor

    real(RP), intent(in)    :: AQ_R   (QA)
    real(RP), intent(in)    :: AQ_CV  (QA)
    real(RP), intent(in)    :: AQ_CP  (QA)
    real(RP), intent(in)    :: AQ_MASS(QA)

    real(RP), intent(in)    :: REF_dens(KA,IA,JA)
    real(RP), intent(in)    :: REF_pott(KA,IA,JA)
    real(RP), intent(in)    :: REF_qv  (KA,IA,JA)
    real(RP), intent(in)    :: REF_pres(KA,IA,JA)   !< reference pressure

    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N
    logical,  intent(in)    :: TwoD

    real(RP), intent(in)    :: ND_COEF
    real(RP), intent(in)    :: ND_COEF_Q
    integer,  intent(in)    :: ND_LAPLACIAN_NUM
    real(RP), intent(in)    :: ND_SFC_FACT
    logical,  intent(in)    :: ND_USE_RS

    integer,  intent(in)    :: BND_QA
    integer,  intent(in)    :: BND_IQ(QA)
    real(RP), intent(in)    :: BND_SMOOTHER_FACT

    real(RP), intent(in)    :: DAMP_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_POTT(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_QTRC(KA,IA,JA,BND_QA)

    real(RP), intent(in)    :: DAMP_alpha_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_POTT(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_QTRC(KA,IA,JA,BND_QA)
    real(RP), intent(in)    :: MFLUX_OFFSET_X(KA,JA,2)
    real(RP), intent(in)    :: MFLUX_OFFSET_Y(KA,IA,2)

    real(RP), intent(in)    :: wdamp_coef(KA)
    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_TRACER_SPLIT_TEND
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    logical,  intent(in)    :: USE_AVERAGE

    integer,  intent(in)    :: I_QV

    real(DP), intent(in)    :: DTLS
    real(DP), intent(in)    :: DTSS
    logical , intent(in)    :: Llast

    ! for time integartion
    real(RP) :: DENS00  (KA,IA,JA) ! saved density before small step loop
    real(RP) :: qflx    (KA,IA,JA,3)

    ! diagnostic variables
    real(RP) :: DDIV    (KA,IA,JA) ! 3 dimensional divergence
    real(RP) :: DPRES0  (KA,IA,JA) ! pressure deviation
    real(RP) :: RT2P    (KA,IA,JA) ! factor of RHOT to PRES
    real(RP) :: REF_rhot(KA,IA,JA) ! reference of RHOT

    real(RP) :: DENS_tq(KA,IA,JA)
    real(RP) :: diff(KA,IA,JA), diff2(KA,IA,JA), diff3(KA,IA,JA)
    real(RP) :: damp
    real(RP) :: damp_t_DENS(KA,IA,JA)
    real(RP) :: damp_t_MOMZ(KA,IA,JA)
    real(RP) :: damp_t_MOMX(KA,IA,JA)
    real(RP) :: damp_t_MOMY(KA,IA,JA)
    real(RP) :: damp_t_RHOT(KA,IA,JA)
    real(RP) :: damp_t_QTRC(KA,IA,JA)

    real(RP) :: tflx(KA,IA,JA,3)

    ! For tracer advection
    real(RP) :: mflx_av(KA,IA,JA,3)  ! rho * vel(x,y,z) @ (u,v,w)-face average

    real(RP) :: dtl
    real(RP) :: dts
    integer  :: nstep

    ! for history
    logical :: do_put

    ! for monitor
    real(RP), target  :: qflx_west (KA,JA)
    real(RP), target  :: qflx_east (KA,JA)
    real(RP), target  :: qflx_south(KA,IA)
    real(RP), target  :: qflx_north(KA,IA)
    logical :: MONIT_lateral_flag(3)

    integer  :: i, j, k, iq, iqb, step
    integer  :: iv
    integer  :: n
    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_Large_Preparation", 2)

    dts   = real(DTSS, kind=RP)            ! short time step
    dtl   = real(DTLS, kind=RP)            ! large time step
    nstep = ceiling( ( dtl - eps ) / dts )
    dts   = dtl / nstep                    ! dts is divisor of dtl and smaller or equal to dtss

    MONIT_lateral_flag(ZDIR) = .false.
    MONIT_lateral_flag(XDIR) = MONIT_mflx_west > 0 .or. MONIT_mflx_east > 0 
    MONIT_lateral_flag(YDIR) = MONIT_mflx_south > 0 .or. MONIT_mflx_north > 0 

#ifdef DEBUG
    LOG_INFO("ATMOS_DYN_Tstep_large_fvm_heve",*)                         'Dynamics large time step'
    LOG_INFO_CONT('(1x,A,F0.2,A,F0.2,A,I0)') &
    '-> DT_large, DT_small, DT_large/DT_small : ', dtl, ', ', dts, ', ', nstep

    DENS00  (:,:,:) = UNDEF
    num_diff (:,:,:,:,:) = UNDEF
#endif

!OCL XFILL
    DENS00(:,:,:) = DENS(:,:,:)

    if ( USE_AVERAGE ) then
!OCL XFILL
       !$omp parallel do
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          DENS_av(k,i,j) = 0.0_RP
          MOMZ_av(k,i,j) = 0.0_RP
          MOMX_av(k,i,j) = 0.0_RP
          MOMY_av(k,i,j) = 0.0_RP
          RHOT_av(k,i,j) = 0.0_RP
       end do
       end do
       end do
    endif

!OCL XFILL
    mflx_av(:,:,:,:) = 0.0_RP


!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do k = KS, KE
       qflx_west(k,j) = 0.0_RP
       qflx_east(k,j) = 0.0_RP
    end do
    end do
!OCL XFILL
    !$omp parallel do
    do i = IS, IE
    do k = KS, KE
       qflx_south(k,i) = 0.0_RP
       qflx_north(k,i) = 0.0_RP
    end do
    end do

    !- prepare some variables for pressure linearization 
    call ATMOS_DYN_prep_pres_linearization( &
       DPRES0, RT2P, REF_rhot,                            & ! (out)
       RHOT, QTRC, REF_pres, AQ_R, AQ_CV, AQ_CP, AQ_MASS  ) ! (in)
   
    call PROF_rapend  ("DYN_Large_Preparation", 2)

    !###########################################################################
    ! Update DENS,MONZ,MOMX,MOMY,MOMZ,RHOT
    !###########################################################################

    call PROF_rapstart("DYN_Large_Tendency", 2)

!OCL XFILL
    DENS_tq(:,:,:) = 0.0_RP

    do iq = 1, QA

       iqb = BND_IQ(iq)

       if ( iqb > 0 ) then

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = JS-1, JE+1
          do i = max(IS-1,1), min(IE+1,IA)
          do k = KS, KE
             diff(k,i,j) = QTRC(k,i,j,iq) - DAMP_QTRC(k,i,j,iqb)
          enddo
          enddo
          enddo

          call FILE_HISTORY_query( HIST_damp_QTRC(iq), do_put )
          if ( TwoD ) then
             !$omp parallel do default(none) OMP_SCHEDULE_ &
             !$omp private(j,k,damp) &
             !$omp shared(IS,JS,JE,KS,KE,iq,iqb) &
             !$omp shared(RHOQ_t,RHOQ_tp,DENS_tq,DAMP_alpha_QTRC,diff,BND_SMOOTHER_FACT,DENS00,TRACER_MASS,I_QV) &
             !$omp shared(damp_t_QTRC,do_put)
!OCL XFILL
             do j = JS, JE
             do k = KS, KE
                damp = - DAMP_alpha_QTRC(k,IS,j,iqb) &
                     * ( diff(k,IS,j) & ! rayleigh damping
                       - ( diff(k,IS,j-1) + diff(k,IS,j+1) - diff(k,IS,j)*2.0_RP ) &
                       * 0.25_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
                damp = damp * DENS00(k,IS,j)
                if ( do_put ) damp_t_QTRC(k,IS,j) = damp
                RHOQ_t(k,IS,j,iq) = RHOQ_tp(k,IS,j,iq) + damp
                DENS_tq(k,IS,j) = DENS_tq(k,IS,j) + damp * TRACER_MASS(iq) ! only for mass tracer
             enddo
             enddo
          else
             if( iq == I_QV .and. SPNUDGE_qv ) then

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff(k,i,j)
                enddo
                enddo
                enddo

                call DFT_g2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,SPNUDGE_qv_lm,SPNUDGE_qv_mm,diff2)

             else

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = 0
                enddo
                enddo
                enddo

             endif

             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,damp) &
             !$omp shared(JS,JE,IS,IE,KS,KE,iq,iqb) &
             !$omp shared(RHOQ_t,RHOQ_tp,DENS_tq,DAMP_alpha_QTRC,diff,diff2,BND_SMOOTHER_FACT,SPNUDGE_qv_alpha,DENS00,TRACER_MASS) &
             !$omp shared(damp_t_QTRC,do_put)
!OCL XFILL
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                damp = - DAMP_alpha_QTRC(k,i,j,iqb) &
                     * ( diff(k,i,j) & ! rayleigh damping
                       - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                       * 0.125_RP * BND_SMOOTHER_FACT ) & ! horizontal smoother
                       - SPNUDGE_qv_alpha(k,i,j) * diff2(k,i,j)
                damp = damp * DENS00(k,i,j)
                if ( do_put ) damp_t_QTRC(k,i,j) = damp
                RHOQ_t(k,i,j,iq) = RHOQ_tp(k,i,j,iq) + damp
                DENS_tq(k,i,j) = DENS_tq(k,i,j) + damp * TRACER_MASS(iq) ! only for mass tracer
             enddo
             enddo
             enddo
          end if

          if ( Llast ) then
             if ( do_put ) call FILE_HISTORY_put( HIST_damp_QTRC(iq), damp_t_QTRC(:,:,:) )
             call FILE_HISTORY_put( HIST_phys_QTRC(iq), RHOQ_tp(:,:,:,iq) )
          end if

          call ATMOS_DYN_fill_halo( RHOQ_t(:,:,:,iq), 0.0_RP, .false., .true. )
          call COMM_vars8( RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq) )
          call COMM_wait ( RHOQ_t(:,:,:,iq), I_COMM_RHOQ_t(iq), .false. )

       else

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = 1, JA
          do i = 1, IA
          do k = 1, KA
             RHOQ_t(k,i,j,iq) = RHOQ_tp(k,i,j,iq)
          enddo
          enddo
          enddo

       end if

    end do

    call PROF_rapend  ("DYN_Large_Tendency", 2)

    call PROF_rapstart("DYN_Large_Boundary", 2)

    if ( BND_W .and. (.not. TwoD) ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_
       do j = JS, JE
       do k = KS, KE
          mflx(k,IS-1,j,XDIR) = ( MOMX(k,IS-1,j) + MFLUX_OFFSET_X(k,j,1) ) &
                              * GSQRT(k,IS-1,j,I_UYZ)  / MAPF(IS-1,j,2,I_UY)
       enddo
       enddo
    end if
    if ( BND_E .and. (.not. TwoD) ) then
       !$omp parallel do private(j,k) OMP_SCHEDULE_
       do j = JS, JE
       do k = KS, KE
          mflx(k,IE,j,XDIR) = ( MOMX(k,IE,j) + MFLUX_OFFSET_X(k,j,2) ) &
                            * GSQRT(k,IE,j,I_UYZ) / MAPF(IE,j,2,I_UY)
       enddo
       enddo
    end if
    if ( BND_S ) then
       !$omp parallel do private(i,k) OMP_SCHEDULE_
       do i = IS, IE
       do k = KS, KE
          mflx(k,i,JS-1,YDIR) = ( MOMY(k,i,JS-1) + MFLUX_OFFSET_Y(k,i,1) ) &
                              * GSQRT(k,i,JS-1,I_XVZ) / MAPF(i,JS-1,1,I_XV)
       enddo
       enddo
    end if
    if ( BND_N ) then
       !$omp parallel do private(i,k) OMP_SCHEDULE_
       do i = IS, IE
       do k = KS, KE
          mflx(k,i,JE,YDIR) = ( MOMY(k,i,JE) + MFLUX_OFFSET_Y(k,i,2) ) &
                            * GSQRT(k,i,JE,I_XVZ) / MAPF(i,JE,1,I_XV)
       enddo
       enddo
    end if

    call PROF_rapend  ("DYN_Large_Boundary", 2)

!OCL XFILL
    !$omp parallel do collapse(2)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       damp_t_DENS(k,i,j) = 0.0_RP
       damp_t_MOMZ(k,i,j) = 0.0_RP
       damp_t_MOMX(k,i,j) = 0.0_RP
       damp_t_MOMY(k,i,j) = 0.0_RP
       damp_t_RHOT(k,i,j) = 0.0_RP
    end do
    end do
    end do

    call ATMOS_DYN_fill_halo( DENS_damp, 0.0_RP, .true., .true. )
    call ATMOS_DYN_fill_halo( DENS_t, 0.0_RP, .false., .true. )

    do step = 1, nstep

       !-----< prepare tendency >-----

       call PROF_rapstart("DYN_Large_Tendency", 2)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+1
       do i = max(IS-1,1), min(IE+1,IA)
       do k = KS, KE
          diff(k,i,j) = DENS(k,i,j) - DAMP_DENS(k,i,j)
       enddo
       enddo
       enddo

       call FILE_HISTORY_query( HIST_damp(1), do_put )
       if ( TwoD ) then
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k) &
          !$omp shared(IS,JS,JE,KS,KE) &
          !$omp shared(DAMP_alpha_DENS,diff,DENS_tq,DENS_t,DENS_tp,BND_SMOOTHER_FACT,EPS) &
          !$omp shared(damp_t_DENS,DENS_damp,do_put,nstep)
!OCL XFILL
          do j = JS, JE
          do k = KS, KE
             DENS_damp(k,IS,j) = - DAMP_alpha_DENS(k,IS,j) &
                  * ( diff(k,IS,j) & ! rayleigh damping
                    - ( diff(k,IS,j-1) + diff(k,IS,j+1) - diff(k,IS,j)*2.0_RP ) &
                    * 0.25_RP * BND_SMOOTHER_FACT ) & ! horizontal smoother
                  + DENS_tq(k,IS,j) * ( 0.5_RP - sign( 0.5_RP, DAMP_alpha_DENS(k,IS,j)-EPS ) ) ! dencity change due to rayleigh damping for tracers
             DENS_t(k,IS,j) = DENS_tp(k,IS,j) & ! tendency from physical step
                           + DENS_damp(k,IS,j)
             if ( do_put ) damp_t_DENS(k,IS,j) = damp_t_DENS(k,IS,j) + DENS_damp(k,IS,j) / nstep
          enddo
          enddo
       else
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k) &
          !$omp shared(JS,JE,IS,IE,KS,KE) &
          !$omp shared(DAMP_alpha_DENS,diff,DENS_tq,DENS_t,DENS_tp,BND_SMOOTHER_FACT,EPS) &
          !$omp shared(damp_t_DENS,DENS_damp,do_put,nstep)
!OCL XFILL
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             DENS_damp(k,i,j) = - DAMP_alpha_DENS(k,i,j) &
                  * ( diff(k,i,j) & ! rayleigh damping
                    - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                    * 0.125_RP * BND_SMOOTHER_FACT ) & ! horizontal smoother
                  + DENS_tq(k,i,j) * ( 0.5_RP - sign( 0.5_RP, DAMP_alpha_DENS(k,i,j)-EPS ) ) ! dencity change due to rayleigh damping for tracers
             DENS_t(k,i,j) = DENS_tp(k,i,j) & ! tendency from physical step
                           + DENS_damp(k,i,j)
             if ( do_put ) damp_t_DENS(k,i,j) = damp_t_DENS(k,i,j) + DENS_damp(k,i,j) / nstep
          enddo
          enddo
          enddo
       end if

       call COMM_vars8( DENS_damp(:,:,:), I_COMM_DENS_damp )
       call COMM_vars8( DENS_t(:,:,:), I_COMM_DENS_t )

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+1
       do i = max(IS-1,1), min(IE+1,IA)
       do k = KS, KE-1
          diff(k,i,j) = MOMZ(k,i,j) - DAMP_VELZ(k,i,j) * ( DENS(k,i,j)+DENS(k+1,i,j) ) * 0.5_RP
       enddo
       enddo
       enddo

       call FILE_HISTORY_query( HIST_damp(2), do_put )
       if ( TwoD ) then
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,damp) &
          !$omp shared(IS,JS,JE,KS,KE) &
          !$omp shared(DENS_damp,DENS,MOMZ) &
          !$omp shared(DAMP_alpha_VELZ,diff,BND_SMOOTHER_FACT,MOMZ_t,MOMZ_tp) &
          !$omp shared(damp_t_MOMZ,do_put,nstep)
!OCL XFILL
          do j = JS, JE
          do k = KS, KE-1
             damp = - DAMP_alpha_VELZ(k,IS,j) &
                  * ( diff(k,IS,j) & ! rayleigh damping
                    - ( diff(k,IS,j-1) + diff(k,IS,j+1) - diff(k,IS,j)*2.0_RP ) &
                    * 0.25_RP * BND_SMOOTHER_FACT ) ! horizontal smoother

             MOMZ_t(k,IS,j) = MOMZ_tp(k,IS,j) & ! tendency from physical step
                           + damp &
                           + ( DENS_damp(k,IS,j) + DENS_damp(k+1,IS,j) ) * MOMZ(k,IS,j) / ( DENS(k,IS,j) + DENS(k+1,IS,j) )
             if ( do_put ) damp_t_MOMZ(k,IS,j) = damp_t_MOMZ(k,IS,j) + damp / nstep
          enddo
          enddo
       else
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,damp) &
          !$omp shared(JS,JE,IS,IE,KS,KE) &
          !$omp shared(DENS_damp,DENS,MOMZ) &
          !$omp shared(DAMP_alpha_VELZ,diff,BND_SMOOTHER_FACT,MOMZ_t,MOMZ_tp) &
          !$omp shared(damp_t_MOMZ,do_put,nstep)
!OCL XFILL
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE-1
             damp = - DAMP_alpha_VELZ(k,i,j) &
                  * ( diff(k,i,j) & ! rayleigh damping
                    - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                    * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother

             MOMZ_t(k,i,j) = MOMZ_tp(k,i,j) & ! tendency from physical step
                           + damp &
                           + ( DENS_damp(k,i,j) + DENS_damp(k+1,i,j) ) * MOMZ(k,i,j) / ( DENS(k,i,j) + DENS(k+1,i,j) )
             if ( do_put ) damp_t_MOMZ(k,i,j) = damp_t_MOMZ(k,i,j) + damp / nstep
          enddo
          enddo
          enddo
       end if
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
          MOMZ_t( 1:KS-1,i,j) = 0.0_RP
          MOMZ_t(KE:KA  ,i,j) = 0.0_RP
       enddo
       enddo
       call COMM_vars8( MOMZ_t(:,:,:), I_COMM_MOMZ_t )

       call COMM_wait( DENS_damp(:,:,:), I_COMM_DENS_damp )

       if ( TwoD ) then
          !$omp parallel do private(j,k) OMP_SCHEDULE_
!OCL XFILL
          do j = JS-1, JE+1
          do k = KS, KE
             diff(k,IS,j) = MOMX(k,IS,j) - DAMP_VELX(k,IS,j) * DENS(k,IS,j)
          enddo
          enddo
       else
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             diff(k,i,j) = MOMX(k,i,j) - DAMP_VELX(k,i,j) * ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP
          enddo
          enddo
          enddo
       end if


       call FILE_HISTORY_query( HIST_damp(3), do_put )
       if ( TwoD ) then
!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,damp) &
          !$omp shared(IS,JS,JE,KS,KE) &
          !$omp shared(DENS_damp,DENS,MOMX) &
          !$omp shared(DAMP_alpha_VELX,diff,BND_SMOOTHER_FACT,MOMX_tp,MOMX_t) &
          !$omp shared(damp_t_MOMX,do_put,nstep)
          do j = JS, JE
          do k = KS, KE
             damp = - DAMP_alpha_VELX(k,IS,j) &
                  * ( diff(k,IS,j) & ! rayleigh damping
                    - ( diff(k,IS,j-1) + diff(k,IS,j+1) - diff(k,IS,j)*2.0_RP ) &
                    * 0.25_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
             MOMX_t(k,IS,j) = MOMX_tp(k,IS,j) & ! tendency from physical step
                           + damp &
                           + DENS_damp(k,IS,j) * MOMX(k,IS,j) / DENS(k,IS,j)
             if ( do_put ) damp_t_MOMX(k,IS,j) = damp_t_MOMX(k,IS,j) + damp / nstep
          enddo
          enddo
       else
          if( SPNUDGE_uv ) then
             if( SPNUDGE_uv_divfree ) then

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS-1, JE+1
                do i = IS-1, IE+1
                do k = KS, KE
                   diff3(k,i,j) = MOMY(k,i,j) - DAMP_VELY(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP
                enddo
                enddo
                enddo

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff(k,i,j) / ( ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP )
                enddo
                enddo
                enddo

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff3(k,i,j) = diff3(k,i,j) / ( ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP )
                enddo
                enddo
                enddo

                call DFT_g2g_divfree(KA,KS,KE,IA,IS,IE,JA,JS,JE,SPNUDGE_uv_lm,SPNUDGE_uv_mm,diff2,diff3)

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff2(k,i,j) * ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP
                enddo
                enddo
                enddo

             else

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff(k,i,j) / ( ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP )
                enddo
                enddo
                enddo

                call DFT_g2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,SPNUDGE_uv_lm,SPNUDGE_uv_mm,diff2)

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff2(k,i,j) * ( DENS(k,i,j)+DENS(k,i+1,j) ) * 0.5_RP
                enddo
                enddo
                enddo

             endif
          else

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                diff2(k,i,j) = 0
             enddo
             enddo
             enddo

          endif

!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,damp) &
          !$omp shared(JS,JE,IS,IE,KS,KE) &
          !$omp shared(DENS_damp,DENS,MOMX) &
          !$omp shared(DAMP_alpha_VELX,diff,diff2,BND_SMOOTHER_FACT,SPNUDGE_u_alpha,MOMX_tp,MOMX_t) &
          !$omp shared(damp_t_MOMX,do_put,nstep)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             damp = - DAMP_alpha_VELX(k,i,j) &
                  * ( diff(k,i,j) & ! rayleigh damping
                    - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                    * 0.125_RP * BND_SMOOTHER_FACT ) & ! horizontal smoother
                  - SPNUDGE_u_alpha(k,i,j) * diff2(k,i,j)
             MOMX_t(k,i,j) = MOMX_tp(k,i,j) & ! tendency from physical step
                           + damp &
                           + ( DENS_damp(k,i,j) + DENS_damp(k,i+1,j) ) * MOMX(k,i,j) / ( DENS(k,i,j) + DENS(k,i+1,j) )
             if ( do_put ) damp_t_MOMX(k,i,j) = damp_t_MOMX(k,i,j) + damp / nstep
          enddo
          enddo
          enddo
       end if

       call ATMOS_DYN_fill_halo( MOMX_t, 0.0_RP, .false., .true. )
       call COMM_vars8( MOMX_t(:,:,:), I_COMM_MOMX_t )

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+1
       do i = max(IS-1,1), min(IE+1,IA)
       do k = KS, KE
          diff(k,i,j) = MOMY(k,i,j) - DAMP_VELY(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP
       enddo
       enddo
       enddo

       call FILE_HISTORY_query( HIST_damp(4), do_put )
       if ( TwoD ) then
!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,damp) &
          !$omp shared(IS,JS,JE,KS,KE) &
          !$omp shared(DENS_damp,DENS,MOMY) &
          !$omp shared(DAMP_alpha_VELY,diff,BND_SMOOTHER_FACT,MOMY_tp,MOMY_t) &
          !$omp shared(damp_t_MOMY,do_put,nstep)
          do j = JS, JE
          do k = KS, KE
             damp = - DAMP_alpha_VELY(k,IS,j) &
                  * ( diff(k,IS,j) & ! rayleigh damping
                    - ( diff(k,IS,j-1) + diff(k,IS,j+1) - diff(k,IS,j)*2.0_RP ) &
                    * 0.25_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
             MOMY_t(k,IS,j) = MOMY_tp(k,IS,j) & ! tendency from physical step
                           + damp &
                           + ( DENS_damp(k,IS,j) + DENS_damp(k,IS,j+1) ) * MOMY(k,IS,j) / ( DENS(k,IS,j) + DENS(k,IS,j+1) )
             if ( do_put ) damp_t_MOMY(k,IS,j) = damp_t_MOMY(k,IS,j) + damp / nstep
          enddo
          enddo
       else
          if( SPNUDGE_uv ) then
             if( SPNUDGE_uv_divfree ) then

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff3(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP
                enddo
                enddo
                enddo

             else
                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff(k,i,j) / ( ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP )
                enddo
                enddo
                enddo

                call DFT_g2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,SPNUDGE_uv_lm,SPNUDGE_uv_mm,diff2)

                !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
                do j = JS, JE
                do i = IS, IE
                do k = KS, KE
                   diff2(k,i,j) = diff2(k,i,j) * ( DENS(k,i,j)+DENS(k,i,j+1) ) * 0.5_RP
                enddo
                enddo
                enddo
             endif
          endif

!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,damp) &
          !$omp shared(JS,JE,IS,IE,KS,KE) &
          !$omp shared(DENS_damp,DENS,MOMY) &
          !$omp shared(DAMP_alpha_VELY,diff,diff2,BND_SMOOTHER_FACT,SPNUDGE_v_alpha,MOMY_tp,MOMY_t) &
          !$omp shared(damp_t_MOMY,do_put,nstep)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             damp = - DAMP_alpha_VELY(k,i,j) &
                  * ( diff(k,i,j) & ! rayleigh damping
                    - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                    * 0.125_RP * BND_SMOOTHER_FACT ) & ! horizontal smoother
                  - SPNUDGE_v_alpha(k,i,j) * diff2(k,i,j)
             MOMY_t(k,i,j) = MOMY_tp(k,i,j) & ! tendency from physical step
                           + damp &
                           + ( DENS_damp(k,i,j) + DENS_damp(k,i,j+1) ) * MOMY(k,i,j) / ( DENS(k,i,j) + DENS(k,i,j+1) )
             if ( do_put ) damp_t_MOMY(k,i,j) = damp_t_MOMY(k,i,j) + damp / nstep
          enddo
          enddo
          enddo
       end if
!OCL XFILL

       call ATMOS_DYN_fill_halo( MOMY_t, 0.0_RP, .false., .true. )
       call COMM_vars8( MOMY_t(:,:,:), I_COMM_MOMY_t )

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
       do j = JS-1, JE+2
       do i = max(IS-1,1), min(IE+2,IA)
       do k = KS, KE
          diff(k,i,j) = RHOT(k,i,j) - DAMP_POTT(k,i,j) * DENS(k,i,j)
       enddo
       enddo
       enddo

       call FILE_HISTORY_query( HIST_damp(5), do_put )
       if ( TwoD ) then
!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ &
          !$omp private(j,k,damp) &
          !$omp shared(IS,JS,JE,KS,KE) &
          !$omp shared(DENS_damp,DENS,RHOT) &
          !$omp shared(DAMP_alpha_POTT,diff,BND_SMOOTHER_FACT,RHOT_t,RHOT_tp) &
          !$omp shared(damp_t_RHOT,do_put,nstep)
          do j = JS, JE
          do k = KS, KE
             damp = - DAMP_alpha_POTT(k,IS,j) &
                  * ( diff(k,IS,j) & ! rayleigh damping
                    - ( diff(k,IS,j-1) + diff(k,IS,j+1) - diff(k,IS,j)*2.0_RP ) &
                    * 0.25_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
             RHOT_t(k,IS,j) = RHOT_tp(k,IS,j) & ! tendency from physical step
                           + damp &
                           + DENS_damp(k,IS,j) * RHOT(k,IS,j) / DENS(k,IS,j)
             if ( do_put ) damp_t_RHOT(k,IS,j) = damp_t_RHOT(k,IS,j) + damp / nstep
          enddo
          enddo
       else
          if( SPNUDGE_pt ) then

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                diff2(k,i,j) = diff(k,i,j) / DENS(k,i,j)
             enddo
             enddo
             enddo

             call DFT_g2g(KA,KS,KE,IA,IS,IE,JA,JS,JE,SPNUDGE_pt_lm,SPNUDGE_pt_mm,diff2)

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                diff2(k,i,j) = diff2(k,i,j) * DENS(k,i,j)
             enddo
             enddo
             enddo

          else

             !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
!OCL XFILL
             do j = JS, JE
             do i = IS, IE
             do k = KS, KE
                diff2(k,i,j) = 0
             enddo
             enddo
             enddo

          endif

!OCL XFILL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,damp) &
          !$omp shared(JS,JE,IS,IE,KS,KE) &
          !$omp shared(DENS_damp,DENS,RHOT) &
          !$omp shared(DAMP_alpha_POTT,diff,diff2,BND_SMOOTHER_FACT,SPNUDGE_pt_alpha,RHOT_t,RHOT_tp) &
          !$omp shared(damp_t_RHOT,do_put,nstep)
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             damp = - DAMP_alpha_POTT(k,i,j) &
                  * ( diff(k,i,j) & ! rayleigh damping
                    - ( diff(k,i-1,j) + diff(k,i+1,j) + diff(k,i,j-1) + diff(k,i,j+1) - diff(k,i,j)*4.0_RP ) &
                    * 0.125_RP * BND_SMOOTHER_FACT ) & ! horizontal smoother
                  - SPNUDGE_pt_alpha(k,i,j) * diff2(k,i,j)
             RHOT_t(k,i,j) = RHOT_tp(k,i,j) & ! tendency from physical step
                           + damp &
                           + DENS_damp(k,i,j) * RHOT(k,i,j) / DENS(k,i,j)
             if ( do_put ) damp_t_RHOT(k,i,j) = damp_t_RHOT(k,i,j) + damp / nstep
          enddo
          enddo
          enddo
       end if
!OCL XFILL

       call ATMOS_DYN_fill_halo( RHOT_t, 0.0_RP, .false., .true. )
       call COMM_vars8( RHOT_t(:,:,:), I_COMM_RHOT_t )

       call COMM_wait ( DENS_t(:,:,:), I_COMM_DENS_t, .false. )
       call COMM_wait ( MOMZ_t(:,:,:), I_COMM_MOMZ_t, .false. )
       call COMM_wait ( MOMX_t(:,:,:), I_COMM_MOMX_t, .false. )
       call COMM_wait ( MOMY_t(:,:,:), I_COMM_MOMY_t, .false. )
       call COMM_wait ( RHOT_t(:,:,:), I_COMM_RHOT_t, .false. )

       call PROF_rapend  ("DYN_Large_Tendency", 2)

       call PROF_rapstart("DYN_Large_Numfilter", 2)

       !-----< prepare the fluxes of explicit numerical diffusion >-----
       if ( ND_COEF == 0.0_RP .or. EVAL_TYPE_NUMFILTER == 'FILTER' ) then
!OCL XFILL
          num_diff(:,:,:,:,:) = 0.0_RP
       else
          call ATMOS_DYN_FVM_numfilter_flux( num_diff(:,:,:,:,:),                          & ! [OUT]
                                         DENS, MOMZ, MOMX, MOMY, RHOT,                     & ! [IN]
                                         CDZ, CDX, CDY, FDZ, FDX, FDY, TwoD, dts,          & ! [IN]
                                         REF_dens, REF_pott,                               & ! [IN]
                                         ND_COEF, ND_LAPLACIAN_NUM, ND_SFC_FACT, ND_USE_RS ) ! [IN]
       endif

       call PROF_rapend  ("DYN_Large_Numfilter", 2)

       !-----< calculate the divegence term >-----

       if ( divdmp_coef > 0.0_RP ) then

          call ATMOS_DYN_divergence( DDIV,    & ! (out)
               MOMZ, MOMX, MOMY,              & ! (in)
               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
               TwoD,                          & ! (in)
               RCDZ, RCDX, RCDY, RFDZ, FDZ    ) ! (in)
      
       else

!XFILL
          DDIV = 0.0_RP

       end if

       !------------------------------------------------------------------------
       ! Start short time integration
       !------------------------------------------------------------------------

       call FILE_HISTORY_set_disable( .not. ( Llast .and. ( step==nstep ) ) )

       call PROF_rapstart("DYN_Short_Tinteg", 2)

       call ATMOS_DYN_tinteg_short( DENS, MOMZ, MOMX, MOMY, RHOT, PROG,      & ! (inout)
                                    mflx, tflx,                              & ! (inout, out)
                                    DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,  & ! (in)
                                    DPRES0, RT2P, CORIOLI,                   & ! (in)
                                    num_diff, wdamp_coef, divdmp_coef, DDIV, & ! (in)
                                    FLAG_FCT_MOMENTUM, FLAG_FCT_T,           & ! (in)
                                    FLAG_FCT_ALONG_STREAM,                   & ! (in)
                                    CDZ, FDZ, FDX, FDY,                      & ! (in)
                                    RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,      & ! (in)
                                    PHI, GSQRT, J13G, J23G, J33G, MAPF,      & ! (in)
                                    REF_dens, REF_rhot,                      & ! (in)
                                    BND_W, BND_E, BND_S, BND_N, TwoD,        & ! (in)
                                    dts                                      ) ! (in)

       call PROF_rapend  ("DYN_Short_Tinteg", 2)


       call FILE_HISTORY_set_disable( .false. )

       if ( ND_COEF > 0.0_RP .and. EVAL_TYPE_NUMFILTER == 'FILTER' ) then
         call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
         call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
         call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
         call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
         call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )
         call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
         call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
         call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
         call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
         call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )

         call ATMOS_DYN_fvm_apply_numfilter( &
            num_diff,                                          & ! (out)
            DENS, MOMZ, MOMX, MOMY, RHOT,                      & ! (inout)
            CDZ, CDX, CDY, FDZ, FDX, FDY,                      & ! (in)
            RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                & ! (in)
            TwoD, dts,                                         & ! (in)
            GSQRT, MAPF, REF_dens, REF_pott,                   & ! (in)
            ND_COEF, ND_LAPLACIAN_NUM, ND_SFC_FACT, ND_USE_RS  ) ! (in)         
       endif

       !$omp parallel do default(none) private(i,j,iv) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JSB,JEB,ISB,IEB,KS,KA,DENS,MOMZ,MOMX,MOMY,RHOT,VA,PROG,KE)
       do j  = JSB, JEB
       do i  = ISB, IEB
          DENS(   1:KS-1,i,j) = DENS(KS,i,j)
          MOMZ(   1:KS-1,i,j) = 0.0_RP
          MOMX(   1:KS-1,i,j) = MOMX(KS,i,j)
          MOMY(   1:KS-1,i,j) = MOMY(KS,i,j)
          RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
          do iv = 1, VA
             PROG(   1:KS-1,i,j,iv) = PROG(KS,i,j,iv)
          end do
          DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
          MOMZ(KE+1:KA,  i,j) = 0.0_RP
          MOMX(KE+1:KA,  i,j) = MOMX(KE,i,j)
          MOMY(KE+1:KA,  i,j) = MOMY(KE,i,j)
          RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
          do iv = 1, VA
             PROG(KE+1:KA,  i,j,iv) = PROG(KE,i,j,iv)
          end do
       enddo
       enddo

       call COMM_vars8( DENS(:,:,:), I_COMM_DENS )
       call COMM_vars8( MOMZ(:,:,:), I_COMM_MOMZ )
       call COMM_vars8( MOMX(:,:,:), I_COMM_MOMX )
       call COMM_vars8( MOMY(:,:,:), I_COMM_MOMY )
       call COMM_vars8( RHOT(:,:,:), I_COMM_RHOT )
       do iv = 1, VA
          call COMM_vars8( PROG(:,:,:,iv), I_COMM_PROG(iv) )
       end do
       call COMM_wait ( DENS(:,:,:), I_COMM_DENS, .false. )
       call COMM_wait ( MOMZ(:,:,:), I_COMM_MOMZ, .false. )
       call COMM_wait ( MOMX(:,:,:), I_COMM_MOMX, .false. )
       call COMM_wait ( MOMY(:,:,:), I_COMM_MOMY, .false. )
       call COMM_wait ( RHOT(:,:,:), I_COMM_RHOT, .false. )
       do iv = 1, VA
          call COMM_wait ( PROG(:,:,:,iv), I_COMM_PROG(iv), .false. )
       end do

       if ( USE_AVERAGE ) then
          !$omp parallel do
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             DENS_av(k,i,j) = DENS_av(k,i,j) + DENS(k,i,j) / nstep
             MOMZ_av(k,i,j) = MOMZ_av(k,i,j) + MOMZ(k,i,j) / nstep
             MOMX_av(k,i,j) = MOMX_av(k,i,j) + MOMX(k,i,j) / nstep
             MOMY_av(k,i,j) = MOMY_av(k,i,j) + MOMY(k,i,j) / nstep
             RHOT_av(k,i,j) = RHOT_av(k,i,j) + RHOT(k,i,j) / nstep
          end do
          end do
          end do
       endif

       !$omp parallel
       do n = 1, 3
       !$omp do
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          mflx_av(k,i,j,n) = mflx_av(k,i,j,n) + mflx(k,i,j,n)
       end do
       end do
       end do
       !$omp end do nowait
       end do
       !$omp end parallel

    enddo ! dynamical steps


    !###########################################################################
    ! Update Tracers
    !###########################################################################

    !$omp parallel
    do n = 1, 3
!OCL XFILL
       !$omp do
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          mflx(k,i,j,n) = mflx_av(k,i,j,n) / nstep
       end do
       end do
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    call COMM_vars8( mflx(:,:,:,ZDIR), I_COMM_mflx_z )
    if ( .not. TwoD ) call COMM_vars8( mflx(:,:,:,XDIR), I_COMM_mflx_x )
    call COMM_vars8( mflx(:,:,:,YDIR), I_COMM_mflx_y )
    call COMM_wait ( mflx(:,:,:,ZDIR), I_COMM_mflx_z, .false. )
    if ( .not. TwoD ) call COMM_wait ( mflx(:,:,:,XDIR), I_COMM_mflx_x, .false. )
    call COMM_wait ( mflx(:,:,:,YDIR), I_COMM_mflx_y, .false. )

#ifndef DRY

    !------------------------------------------------------------------------
    ! Update each tracer
    !------------------------------------------------------------------------

    do iq = 1, QA

       if ( TRACER_ADVC(iq) ) then

          call PROF_rapstart("DYN_Large_Numfilter", 2)

          if ( ND_COEF_Q == 0.0_RP ) then
!OCL XFILL
             num_diff_q(:,:,:,:) = 0.0_RP
          else
             call ATMOS_DYN_FVM_numfilter_flux_q( num_diff_q(:,:,:,:),                & ! [OUT]
                                              DENS00, QTRC(:,:,:,iq), iq==I_QV,       & ! [IN]
                                              CDZ, CDX, CDY, TwoD, dtl,               & ! [IN]
                                              REF_qv, iq,                             & ! [IN]
                                              ND_COEF_Q, ND_LAPLACIAN_NUM, ND_SFC_FACT, ND_USE_RS ) ! [IN]
          endif

          call PROF_rapend  ("DYN_Large_Numfilter", 2)

          call PROF_rapstart("DYN_Tracer_Tinteg", 2)

          if ( FLAG_TRACER_SPLIT_TEND ) then
             !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
             !$omp shared(KA,IA,JA,iq,QTRC,RHOQ_t,DENS,dtl)
             do j = 1, JA
             do i = 1, IA
             do k = 1, KA
                QTRC(k,i,j,iq) = QTRC(k,i,j,iq) &
                               + RHOQ_t(k,i,j,iq) * dtl / DENS(k,i,j)
             end do
             end do
             end do
             RHOQ_tn => ZERO
          else
             RHOQ_tn => RHOQ_t(:,:,:,iq)
          end if

          call ATMOS_DYN_tinteg_tracer( &
               QTRC(:,:,:,iq),              & ! (inout)
               qflx(:,:,:,:),               & ! (out)
               QTRC0(:,:,:,iq), RHOQ_tn,    & ! (in)
               DENS00, DENS,                & ! (in)
               mflx, num_diff_q,            & ! (in)
               GSQRT, MAPF(:,:,:,I_XY),     & ! (in)
               CDZ, RCDZ, RCDX, RCDY,       & ! (in)
               BND_W, BND_E, BND_S, BND_N,  & ! (in)
               TwoD,                        & ! (in)
               dtl,                         & ! (in)
               Llast .AND. FLAG_FCT_TRACER, & ! (in)
               FLAG_FCT_ALONG_STREAM        ) ! (in)

          call PROF_rapend  ("DYN_Tracer_Tinteg", 2)

          if ( Llast ) then
             do iv = 1, 3
              call FILE_HISTORY_query( HIST_qflx(iv,iq), do_put )
              if ( do_put .or. MONIT_lateral_flag(iv) ) then
               call multiply_flux_by_metric_xyz( iv, qflx, GSQRT, MAPF )               
               call FILE_HISTORY_put( HIST_qflx(iv,iq), qflx(:,:,:,iv) )
              end if 
             end do

             if ( TRACER_MASS(iq) == 1.0_RP ) then
                if ( BND_W .and. MONIT_qflx_west > 0 ) then
                   !$omp parallel do
                   do j = JS, JE
                   do k = KS, KE
                      qflx_west(k,j) = qflx_west(k,j) + qflx(k,IS-1,j,XDIR)
                   end do
                   end do
                end if
                if ( BND_E .and. MONIT_qflx_east > 0 ) then
                   !$omp parallel do
                   do j = JS, JE
                   do k = KS, KE
                      qflx_east(k,j) = qflx_east(k,j) + qflx(k,IE,j,XDIR)
                   end do
                   end do
                end if
                if ( BND_S .and. MONIT_qflx_south > 0 ) then
                   !$omp parallel do
                   do i = IS, IE
                   do k = KS, KE
                      qflx_south(k,i) = qflx_south(k,i) + qflx(k,i,JS-1,YDIR)
                   end do
                   end do
                end if
                if ( BND_N .and. MONIT_qflx_north > 0 ) then
                   !$omp parallel do
                   do i = IS, IE
                   do k = KS, KE
                      qflx_north(k,i) = qflx_north(k,i) + qflx(k,i,JE,YDIR)
                   end do
                   end do
                end if

             end if

          end if

       else

          !$omp parallel do
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq) = ( QTRC0(k,i,j,iq) * DENS00(k,i,j) &
                              + RHOQ_t(k,i,j,iq) * dtl          ) / DENS(k,i,j)
          end do
          end do
          end do

       end if

       call COMM_vars8( QTRC(:,:,:,iq), I_COMM_QTRC(iq) )

       if ( USE_AVERAGE ) then
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             QTRC_av(k,i,j,iq) = QTRC(k,i,j,iq)
          end do
          end do
          end do
       endif

    enddo ! scalar quantities loop

    do iq = 1, QA
       call COMM_wait ( QTRC(:,:,:,iq), I_COMM_QTRC(iq), .false. )
    enddo
#endif

    if ( Llast ) then
       !- history ---------------------------------------
       call FILE_HISTORY_put( HIST_phys(1), DENS_tp )
       call FILE_HISTORY_put( HIST_phys(2), MOMZ_tp )
       call FILE_HISTORY_put( HIST_phys(3), MOMX_tp )
       call FILE_HISTORY_put( HIST_phys(4), MOMY_tp )
       call FILE_HISTORY_put( HIST_phys(5), RHOT_tp )

       call FILE_HISTORY_put( HIST_damp(1), damp_t_DENS )
       call FILE_HISTORY_put( HIST_damp(2), damp_t_MOMZ )
       call FILE_HISTORY_put( HIST_damp(3), damp_t_MOMX )
       call FILE_HISTORY_put( HIST_damp(4), damp_t_MOMY )
       call FILE_HISTORY_put( HIST_damp(5), damp_t_RHOT )

       do iv = 1, 3
         call FILE_HISTORY_query( HIST_mflx(iv), do_put )
         if ( do_put .or. MONIT_lateral_flag(iv) ) then
           call multiply_flux_by_metric_xyz( iv, mflx, GSQRT, MAPF )
           call FILE_HISTORY_put( HIST_mflx(iv), mflx(:,:,:,iv) )
         end if   
       end do

       do iv = 1, 3
         call FILE_HISTORY_query( HIST_tflx(iv), do_put )
         if ( do_put ) then
           call multiply_flux_by_metric_xyz( iv, tflx, GSQRT, MAPF )
           call FILE_HISTORY_put( HIST_tflx(iv), tflx(:,:,:,iv) )
         end if   
       end do

       !- monitor mass budget ------------------------------------
       call MONITOR_put( MONIT_damp_mass, damp_t_DENS(:,:,:) )
       if ( IS>0 ) &
       call MONITOR_put_lateral_flux( MONIT_mflx_west, BND_W .and. MONIT_mflx_west > 0, mflx(:,IS-1,:,XDIR), zero_x )
       call MONITOR_put_lateral_flux( MONIT_mflx_east, BND_E .and. MONIT_mflx_east > 0, mflx(:,IE,:,XDIR), zero_x )
       if ( JS>0 ) &
       call MONITOR_put_lateral_flux( MONIT_mflx_south, BND_S .and. MONIT_mflx_south > 0, mflx(:,:,JS-1,YDIR), zero_y )
       call MONITOR_put_lateral_flux( MONIT_mflx_north, BND_N .and. MONIT_mflx_north > 0, mflx(:,:,JE,YDIR), zero_y )

       call MONITOR_put( MONIT_damp_qtot, DENS_tq(:,:,:) )
       call MONITOR_put_lateral_flux( MONIT_qflx_west, BND_W, qflx_west, zero_x )
       call MONITOR_put_lateral_flux( MONIT_qflx_east, BND_E, qflx_east, zero_x )
       call MONITOR_put_lateral_flux( MONIT_qflx_south, BND_S, qflx_south, zero_y )       
       call MONITOR_put_lateral_flux( MONIT_qflx_north, BND_N, qflx_north, zero_y )
    end if

    return
  end subroutine ATMOS_DYN_Tstep_large_fvm_heve

  !-- private subroutines --------------------------------------------
  
  subroutine MONITOR_put_lateral_flux( MONIT_ID, BND_flag, flx, flx_zero )
    use scale_monitor, only: &
      MONITOR_put   
    implicit none

    integer, intent(in) :: MONIT_ID
    logical, intent(in) :: BND_flag
    real(RP), target, intent(in) :: flx(:,:)
    real(RP), target, intent(in) :: flx_zero(:,:)

    real(RP), pointer :: flx_ptr(:,:)
    !------------------------------------------
    if ( BND_flag ) then
      flx_ptr => flx
    else
      flx_ptr => flx_zero
    end if
    call MONITOR_put( MONIT_ID, flx_ptr(:,:) )

    return
  end subroutine MONITOR_put_lateral_flux

  subroutine multiply_flux_by_metric_xyz( I_DIR, flx, GSQRT, MAPF )
    implicit none
    integer, intent(in) :: I_DIR
    real(RP), intent(inout) :: flx(KA,IA,JA,3)
    real(RP), intent(in)    :: GSQRT(KA,IA,JA,7)
    real(RP), intent(in)    :: MAPF (IA,JA,2,4)

    integer :: i, j, k
    !---------------------------------------------------

    if (I_DIR == ZDIR) then
      !$omp parallel do
      do j = JS, JE
      do i = IS, IE
      do k = KS-1, KE
         flx(k,i,j,ZDIR) = flx(k,i,j,ZDIR) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY)
      end do
      end do
      end do
    end if 

    if (I_DIR == XDIR) then
      !$omp parallel do
      do j = JS, JE
      do i = ISB, IEB
      do k = KS, KE
         flx(k,i,j,XDIR) = flx(k,i,j,XDIR) * MAPF(i,j,2,I_UY) / GSQRT(k,i,j,I_UYZ)
      end do
      end do
      end do
    end if 

    if (I_DIR == YDIR) then
      !$omp parallel do
      do j = JSB, JEB
      do i = IS, IE
      do k = KS, KE
        flx(k,i,j,YDIR) = flx(k,i,j,YDIR) * MAPF(i,j,1,I_XV) / GSQRT(k,i,j,I_XVZ)
      end do
      end do
      end do
    end if 

    return
  end subroutine multiply_flux_by_metric_xyz

end module scale_atmos_dyn_tstep_large_fvm_heve
