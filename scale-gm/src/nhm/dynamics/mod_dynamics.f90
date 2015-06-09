!-------------------------------------------------------------------------------
!>
!! Dynamical core
!!
!! @par Description
!!         This module is the core module of fluid dynamics
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita ) Imported from igdc-4.33
!! @li      2006-04-17 (H.Tomita ) Add IN_LARGE_STEP2
!! @li      2007-05-08 (H.Tomita ) Change the treatment of I_TKE.
!! @li      2008-01-24 (Y.Niwa   ) add revised MIURA2004 for tracer advection
!! @li      2008-05-24 (T.Mitsui ) fix miss-conditioning for frcvar
!! @li      2008-09-09 (Y.Niwa   ) move nudging routine here
!! @li      2009-09-08 (S.Iga    ) frhog and frhog_pl in ndg are deleted ( suggested by ES staff)
!! @li      2010-05-06 (M.Satoh  ) define QV_conv only if CP_TYPE='TDK' .or. 'KUO'
!! @li      2010-07-16 (A.T.Noda ) bug fix for TDK
!! @li      2010-08-16 (A.T.Noda ) Bug fix (Qconv not diveded by density)
!! @li      2010-08-20 (A.T.Noda ) Bug fix (Qconv should be TEND, and not be multiplied by DT)
!! @li      2010-11-29 (A.T.Noda ) Introduce the Smagorinsky model
!! @li      2011-08-16 (M.Satoh  ) Bug fix for TDK: conv => TEND; qv_dyn_tend = v grad q = ( div(rho v q) - div(rho v)*q )/rho
!! @li      2011-08-16 (M.Satoh  ) Move codes related to CP_TYPE below the tracer calculation
!! @li      2013-08-23 (H.Yashiro) Change arguments from character to index/switch
!! @li      2013-06-13 (R.Yoshida) Add tracer advection mode
!!
!<
module mod_dynamics
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: dynamics_setup
  public :: dynamics_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: I_RHOG     = 1 ! Density x G^1/2 x gamma^2
  integer, private, parameter :: I_RHOGVX   = 2 ! Density x G^1/2 x gamma^2 x Horizontal velocity (X-direction)
  integer, private, parameter :: I_RHOGVY   = 3 ! Density x G^1/2 x gamma^2 x Horizontal velocity (Y-direction)
  integer, private, parameter :: I_RHOGVZ   = 4 ! Density x G^1/2 x gamma^2 x Horizontal velocity (Z-direction)
  integer, private, parameter :: I_RHOGW    = 5 ! Density x G^1/2 x gamma^2 x Vertical   velocity
  integer, private, parameter :: I_RHOGE    = 6 ! Density x G^1/2 x gamma^2 x Internal Energy
  integer, private, parameter :: I_RHOGETOT = 7 ! Density x G^1/2 x gamma^2 x Total Energy

  integer, private :: num_of_iteration_lstep    ! number of large steps ( 0-4 )
  integer, private :: num_of_iteration_sstep(4) ! number of small steps in each of large steps

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine dynamics_setup
    use mod_adm, only: &
       ADM_proc_stop
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_SSTEP_MAX
    use mod_runconf, only: &
       TRC_ADV_TYPE
    use mod_bndcnd, only: &
       bndcnd_setup
    use mod_bsstate, only: &
       bsstate_setup
    use mod_numfilter, only: &
       numfilter_setup
    use mod_vi, only: &
       vi_setup
!    use mod_sgs, only: &
!       sgs_setup
    use mod_nudge, only: &
       NDG_setup
    implicit none
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*)
    write(ADM_LOG_FID,*) '+++ Module[dynamics]/Category[nhm]'

    write(ADM_LOG_FID,*) '+++ Time integration type: ', trim(TIME_INTEG_TYPE)
    select case(TIME_INTEG_TYPE)
    case('RK2')
       write(ADM_LOG_FID,*) '+++ 2-stage Runge-Kutta'

       num_of_iteration_lstep    = 2
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX

    case('RK3')
       write(ADM_LOG_FID,*) '+++ 3-stage Runge-Kutta'

       num_of_iteration_lstep    = 3
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX

    case('RK4')
       write(ADM_LOG_FID,*) '+++ 4-stage Runge-Kutta'

       num_of_iteration_lstep    = 4
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 4
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(4) = TIME_SSTEP_MAX

    case('TRCADV')
       write(ADM_LOG_FID,*) '+++ Offline tracer experiment'

       num_of_iteration_lstep    = 0

       if ( TRC_ADV_TYPE == 'DEFAULT' ) then
          write(ADM_LOG_FID,*) 'xxx unsupported advection scheme for TRCADV test! STOP.'
          call ADM_proc_stop
       endif

    case default
       write(ADM_LOG_FID,*) 'xxx unsupported integration type! STOP.'
       call ADM_proc_stop
    endselect

    !---< boundary condition module setup >---
    call bndcnd_setup

    !---< basic state module setup >---
    call bsstate_setup

    !---< numerical filter module setup >---
    call numfilter_setup

    !---< vertical implicit module setup >---
    call vi_setup

    !---< sub-grid scale dynamics module setup >---
!    call sgs_setup

    !---< nudging module setup >---
    call NDG_setup

    return
  end subroutine dynamics_setup

  !-----------------------------------------------------------------------------
  subroutine dynamics_step
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_gall_1d, &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_cnst, only: &
       Rdry  => CNST_RAIR,  &
       Rvap  => CNST_RVAP,  &
       CVdry => CNST_CV,    &
       KAPPA => CNST_KAPPA, &
       PRE00 => CNST_PRE00
    use mod_comm, only: &
       COMM_data_transfer
    use mod_vmtr, only: &
       VMTR_GSGAM2,       &
       VMTR_GSGAM2_pl,    &
       VMTR_PHI,          &
       VMTR_PHI_pl,       &
       VMTR_C2WfactGz,    &
       VMTR_C2WfactGz_pl, &
       VMTR_C2Wfact,      &
       VMTR_C2Wfact_pl
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_DTL,        &
       TIME_DTS,        &
       TIME_SPLIT,      &
       TIME_CTIME
    use mod_runconf, only: &
       TRC_VMAX,       &
       I_QV,           &
       I_TKE,          &
       NQW_STR,        &
       NQW_END,        &
       CVW,            &
       DYN_DIV_NUM,    &
       NDIFF_LOCATION, &
       TRC_ADV_TYPE,   &
       FLAG_NUDGING,   &
!       TB_TYPE,        &
       THUBURN_LIM
    use mod_prgvar, only: &
       prgvar_set, &
       prgvar_get
    use mod_bndcnd, only: &
       bndcnd_all
    use mod_bsstate, only: &
       rho_bs,    &
       rho_bs_pl, &
       pre_bs,    &
       pre_bs_pl
    use mod_thrmdyn, only: &
       THRMDYN_th, &
       THRMDYN_eth
    use mod_numfilter, only: &
       NUMFILTER_DOrayleigh,       &
       NUMFILTER_DOverticaldiff,   &
       numfilter_rayleigh_damping, &
       numfilter_hdiffusion,       &
       numfilter_vdiffusion
    use mod_src, only: &
       src_advection_convergence_momentum, &
       src_advection_convergence,          &
       I_SRC_default
    use mod_vi, only: &
       vi_small_step
    use mod_src_tracer, only: &
       src_tracer_advection
    use mod_forcing_driver, only: &
       forcing_update
!    use mod_sgs, only: &
!       sgs_smagorinsky
    use mod_nudge, only: &
       NDG_update_reference, &
       NDG_apply_uvtp
    !##### OpenACC (for data copy) #####
    use mod_prgvar, only: &
       PRG_var,  &
       PRG_var1
    !##### OpenACC #####
    implicit none

    real(RP) :: PROG         (ADM_gall,   ADM_kall,ADM_lall,   6)         ! prognostic variables
    real(RP) :: PROG_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: PROGq        (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  ! tracer variables
    real(RP) :: PROGq_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: g_TEND       (ADM_gall,   ADM_kall,ADM_lall,   7)         ! tendency of prognostic variables
    real(RP) :: g_TEND_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,7)
    real(RP) :: g_TENDq      (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  ! tendency of tracer variables
    real(RP) :: g_TENDq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: f_TEND       (ADM_gall,   ADM_kall,ADM_lall,   7)         ! forcing tendency of prognostic variables
    real(RP) :: f_TEND_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,7)
    real(RP) :: f_TENDq      (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  ! forcing tendency of tracer variables
    real(RP) :: f_TENDq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: PROG0        (ADM_gall,   ADM_kall,ADM_lall,   6)         ! prognostic variables (save)
    real(RP) :: PROG0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: PROGq0       (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)  ! tracer variables (save)
    real(RP) :: PROGq0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: PROG_split   (ADM_gall,   ADM_kall,ADM_lall,   6)         ! prognostic variables (split)
    real(RP) :: PROG_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP) :: PROG_mean    (ADM_gall,   ADM_kall,ADM_lall   ,5)
    real(RP) :: PROG_mean_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,5)

    !--- density ( physical )
    real(RP) :: rho   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rho_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_x  ( physical )
    real(RP) :: vx   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_y  ( physical )
    real(RP) :: vy   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- horizontal velocity_z  ( physical )
    real(RP) :: vz   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- vertical velocity ( physical )
    real(RP) :: w   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: w_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- internal energy  ( physical )
    real(RP) :: ein   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: ein_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- mass concentration of water substance ( physical )
    real(RP) :: q   (ADM_gall,   ADM_kall,ADM_lall,   TRC_VMAX)
    real(RP) :: q_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- enthalpy ( physical )
    real(RP) :: eth   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: eth_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure ( physical )
    real(RP) :: pre   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: pre_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temperature ( physical )
    real(RP) :: tem   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: tem_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- potential temperature ( physical )
    real(RP) :: th   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: th_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- density deviation from the base state ( G^1/2 X gamma2 )
    real(RP) :: rhogd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: rhogd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure deviation from the base state ( G^1/2 X gamma2 )
    real(RP) :: pregd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: pregd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temporary variables
    real(RP) :: qd   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: qd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: cv   (ADM_gall,   ADM_kall,ADM_lall   )
    real(RP) :: cv_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: pre0_kappa

    real(RP), parameter :: TKE_MIN = 0.0_RP
    real(RP)            :: TKEg_corr

    integer :: small_step_ite
    real(RP) :: large_step_dt
    real(RP) :: small_step_dt

    logical :: ndg_TEND_out
    logical :: do_tke_correction

    integer :: g, k ,l, nq, nl, ndyn, m

    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------
    !$acc wait

    call DEBUG_rapstart('__Dynamics')

    !$acc  data &
    !$acc& create(PROG,PROGq,g_TEND,g_TENDq,f_TEND,f_TENDq,PROG0,PROGq0,PROG_split,PROG_mean) &
    !$acc& create(rho,vx,vy,vz,w,ein,tem,pre,eth,th,rhogd,pregd,q,qd,cv) &
    !$acc& pcopy(PRG_var,PRG_var1)

    call DEBUG_rapstart('___Pre_Post')

    large_step_dt = TIME_DTL / real(DYN_DIV_NUM,kind=RP)

    !--- get from prg0
    call prgvar_get( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & ! [OUT]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & ! [OUT]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & ! [OUT]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & ! [OUT]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & ! [OUT]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & ! [OUT]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:),       & ! [OUT]
                     0                                              ) ! [IN]

    call DEBUG_rapend  ('___Pre_Post')

    do ndyn = 1, DYN_DIV_NUM

    call DEBUG_rapstart('___Pre_Post')

    !--- save
    !$acc kernels pcopy(PROG0) pcopyin(PROG) async(0)
    PROG0   (:,:,:,:) = PROG   (:,:,:,:)
    !$acc end kernels

    if ( ADM_have_pl ) then
       PROG0_pl(:,:,:,:) = PROG_pl(:,:,:,:)
    endif

    if ( TRC_ADV_TYPE == 'DEFAULT' ) then
       !$acc kernels pcopy(PROGq0) pcopyin(PROGq) async(0)
       PROGq0   (:,:,:,:) = PROGq   (:,:,:,:)
       !$acc end kernels

       if ( ADM_have_pl ) then
          PROGq0_pl(:,:,:,:) = PROGq_pl(:,:,:,:)
       endif
    endif

    call DEBUG_rapend  ('___Pre_Post')

    if ( TIME_INTEG_TYPE == 'TRCADV' ) then  ! TRC-ADV Test Bifurcation

       call DEBUG_rapstart('___Tracer_Advection')

       !$acc kernels pcopy(f_TEND) async(0)
       f_TEND   (:,:,:,:) = 0.0_RP
       f_TEND_pl(:,:,:,:) = 0.0_RP
       !$acc end kernels

       call src_tracer_advection( TRC_VMAX,                                          & ! [IN]
                                  PROGq (:,:,:,:),        PROGq_pl (:,:,:,:),        & ! [INOUT]
                                  PROG0 (:,:,:,I_RHOG),   PROG0_pl (:,:,:,I_RHOG),   & ! [IN]
                                  PROG  (:,:,:,I_RHOG  ), PROG_pl  (:,:,:,I_RHOG  ), & ! [IN]
                                  PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                  PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                  PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                  PROG  (:,:,:,I_RHOGW ), PROG_pl  (:,:,:,I_RHOGW ), & ! [IN]
                                  f_TEND(:,:,:,I_RHOG),   f_TEND_pl(:,:,:,I_RHOG),   & ! [IN]
                                  large_step_dt,                                     & ! [IN]
                                  THUBURN_LIM                                        ) ! [IN]

       call DEBUG_rapend  ('___Tracer_Advection')

       call forcing_update( PROG(:,:,:,:), PROG_pl(:,:,:,:) ) ! [INOUT]
    endif

    !---------------------------------------------------------------------------
    !
    !> Start large time step integration
    !
    !---------------------------------------------------------------------------
    do nl = 1, num_of_iteration_lstep

       call DEBUG_rapstart('___Pre_Post')

       !---< Generate diagnostic values and set the boudary conditions
       !$acc  kernels pcopy(rho,vx,vy,vz,ein) pcopyin(PROG,VMTR_GSGAM2) async(0)
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          rho(g,k,l) = PROG(g,k,l,I_RHOG  ) / VMTR_GSGAM2(g,k,l)
          vx (g,k,l) = PROG(g,k,l,I_RHOGVX) / PROG(g,k,l,I_RHOG)
          vy (g,k,l) = PROG(g,k,l,I_RHOGVY) / PROG(g,k,l,I_RHOG)
          vz (g,k,l) = PROG(g,k,l,I_RHOGVZ) / PROG(g,k,l,I_RHOG)
          ein(g,k,l) = PROG(g,k,l,I_RHOGE ) / PROG(g,k,l,I_RHOG)
       enddo
       enddo
       enddo
       !$acc end kernels

       !$acc kernels pcopy(q) pcopyin(PROGq,PROG) async(0)
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          !$acc loop seq
          do nq = 1, TRC_VMAX
             q(g,k,l,nq) = PROGq(g,k,l,nq) / PROG(g,k,l,I_RHOG)
          enddo
          !$acc end loop
       enddo
       enddo
       enddo
       !$acc end kernels

       !$acc kernels pcopy(cv,qd,tem,pre) pcopyin(q,ein,rho,CVW) async(0)
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          cv(g,k,l) = 0.0_RP
          qd(g,k,l) = 1.0_RP

          !$acc loop seq
          do nq = NQW_STR, NQW_END
             cv(g,k,l) = cv(g,k,l) + q(g,k,l,nq) * CVW(nq)
             qd(g,k,l) = qd(g,k,l) - q(g,k,l,nq)
          enddo
          !$acc end loop

          cv(g,k,l) = cv(g,k,l) + qd(g,k,l) * CVdry

          tem(g,k,l) = ein(g,k,l) / cv(g,k,l)
          pre(g,k,l) = rho(g,k,l) * tem(g,k,l) * ( qd(g,k,l)*Rdry + q(g,k,l,I_QV)*Rvap )
       enddo
       enddo
       enddo
       !$acc end kernels

       !$acc kernels pcopy(w) pcopyin(PROG,VMTR_C2Wfact) async(0)
       do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          w(g,k,l) = PROG(g,k,l,I_RHOGW) / ( VMTR_C2Wfact(g,k,1,l) * PROG(g,k  ,l,I_RHOG) &
                                           + VMTR_C2Wfact(g,k,2,l) * PROG(g,k-1,l,I_RHOG) )
       enddo
       enddo
       enddo
       !$acc end kernels

       do l = 1, ADM_lall
          !--- boundary conditions
          call bndcnd_all( ADM_gall,                & ! [IN]
                           rho (:,:,l),             & ! [INOUT]
                           vx  (:,:,l),             & ! [INOUT]
                           vy  (:,:,l),             & ! [INOUT]
                           vz  (:,:,l),             & ! [INOUT]
                           w   (:,:,l),             & ! [INOUT]
                           ein (:,:,l),             & ! [INOUT]
                           tem (:,:,l),             & ! [INOUT]
                           pre (:,:,l),             & ! [INOUT]
                           PROG(:,:,l,I_RHOG),      & ! [INOUT]
                           PROG(:,:,l,I_RHOGVX),    & ! [INOUT]
                           PROG(:,:,l,I_RHOGVY),    & ! [INOUT]
                           PROG(:,:,l,I_RHOGVZ),    & ! [INOUT]
                           PROG(:,:,l,I_RHOGW),     & ! [INOUT]
                           PROG(:,:,l,I_RHOGE),     & ! [INOUT]
                           VMTR_GSGAM2   (:,:,l),   & ! [IN]
                           VMTR_PHI      (:,:,l),   & ! [IN]
                           VMTR_C2Wfact  (:,:,:,l), & ! [IN]
                           VMTR_C2WfactGz(:,:,:,l)  ) ! [IN]
       enddo

       pre0_kappa = PRE00**KAPPA

       !$acc kernels pcopy(th,eth) pcopyin(tem,pre,rho,ein) async(0)
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          th (g,k,l) = tem(g,k,l) + pre(g,k,l)**KAPPA * pre0_kappa
          eth(g,k,l) = ein(g,k,l) + pre(g,k,l) / rho(g,k,l)
       enddo
       enddo
       enddo
       !$acc end kernels

       !--- perturbations ( pre, rho with metrics )
       !$acc  kernels pcopy(pregd,rhogd) pcopyin(pre,pre_bs,rho,rho_bs,VMTR_GSGAM2) async(0)
       pregd(:,:,:) = ( pre(:,:,:) - pre_bs(:,:,:) ) * VMTR_GSGAM2(:,:,:)
       rhogd(:,:,:) = ( rho(:,:,:) - rho_bs(:,:,:) ) * VMTR_GSGAM2(:,:,:)
       !$acc end kernels

       if ( ADM_have_pl ) then

          rho_pl(:,:,:) = PROG_pl(:,:,:,I_RHOG  ) / VMTR_GSGAM2_pl(:,:,:)
          vx_pl (:,:,:) = PROG_pl(:,:,:,I_RHOGVX) / PROG_pl(:,:,:,I_RHOG)
          vy_pl (:,:,:) = PROG_pl(:,:,:,I_RHOGVY) / PROG_pl(:,:,:,I_RHOG)
          vz_pl (:,:,:) = PROG_pl(:,:,:,I_RHOGVZ) / PROG_pl(:,:,:,I_RHOG)
          ein_pl(:,:,:) = PROG_pl(:,:,:,I_RHOGE ) / PROG_pl(:,:,:,I_RHOG)

          do nq = 1, TRC_VMAX
             q_pl(:,:,:,nq) = PROGq_pl(:,:,:,nq) / PROG_pl(:,:,:,I_RHOG)
          enddo

          cv_pl(:,:,:) = 0.0_RP
          qd_pl(:,:,:) = 1.0_RP
          do nq = NQW_STR, NQW_END
             cv_pl(:,:,:) = cv_pl(:,:,:) + q_pl(:,:,:,nq) * CVW(nq)
             qd_pl(:,:,:) = qd_pl(:,:,:) - q_pl(:,:,:,nq)
          enddo
          cv_pl(:,:,:) = cv_pl(:,:,:) + qd_pl(:,:,:) * CVdry

          tem_pl(:,:,:) = ein_pl(:,:,:) / cv_pl(:,:,:)
          pre_pl(:,:,:) = rho_pl(:,:,:) * tem_pl(:,:,:) * ( qd_pl(:,:,:)*Rdry + q_pl(:,:,:,I_QV)*Rvap )

          do l = 1, ADM_lall_pl
             do k = ADM_kmin+1, ADM_kmax
             do g = 1, ADM_gall_pl
                w_pl(g,k,l) = PROG_pl(g,k,l,I_RHOGW) / ( VMTR_C2Wfact_pl(g,k,1,l) * PROG_pl(g,k  ,l,I_RHOG) &
                                                       + VMTR_C2Wfact_pl(g,k,2,l) * PROG_pl(g,k-1,l,I_RHOG) )
             enddo
             enddo

             !--- boundary conditions
             call bndcnd_all( ADM_gall_pl,                & ! [IN]
                              rho_pl (:,:,l),             & ! [INOUT]
                              vx_pl  (:,:,l),             & ! [INOUT]
                              vy_pl  (:,:,l),             & ! [INOUT]
                              vz_pl  (:,:,l),             & ! [INOUT]
                              w_pl   (:,:,l),             & ! [INOUT]
                              ein_pl (:,:,l),             & ! [INOUT]
                              tem_pl (:,:,l),             & ! [INOUT]
                              pre_pl (:,:,l),             & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOG  ),    & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGVX),    & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGVY),    & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGVZ),    & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGW ),    & ! [INOUT]
                              PROG_pl(:,:,l,I_RHOGE ),    & ! [INOUT]
                              VMTR_GSGAM2_pl   (:,:,l),   & ! [IN]
                              VMTR_PHI_pl      (:,:,l),   & ! [IN]
                              VMTR_C2Wfact_pl  (:,:,:,l), & ! [IN]
                              VMTR_C2WfactGz_pl(:,:,:,l)  ) ! [IN]
          enddo

          call THRMDYN_th ( ADM_gall_pl,   & ! [IN]
                            ADM_kall,      & ! [IN]
                            ADM_lall_pl,   & ! [IN]
                            tem_pl(:,:,:), & ! [IN]
                            pre_pl(:,:,:), & ! [IN]
                            th_pl (:,:,:)  ) ! [OUT]

          call THRMDYN_eth( ADM_gall_pl,   & ! [IN]
                            ADM_kall,      & ! [IN]
                            ADM_lall_pl,   & ! [IN]
                            ein_pl(:,:,:), & ! [IN]
                            pre_pl(:,:,:), & ! [IN]
                            rho_pl(:,:,:), & ! [IN]
                            eth_pl(:,:,:)  ) ! [OUT]

          pregd_pl(:,:,:) = ( pre_pl(:,:,:) - pre_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
          rhogd_pl(:,:,:) = ( rho_pl(:,:,:) - rho_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
       else

          PROG_pl  (:,:,:,:) = 0.0_RP

          rho_pl(:,:,:) = 0.0_RP
          vx_pl (:,:,:) = 0.0_RP
          vy_pl (:,:,:) = 0.0_RP
          vz_pl (:,:,:) = 0.0_RP
          w_pl  (:,:,:) = 0.0_RP
          ein_pl(:,:,:) = 0.0_RP

          q_pl(:,:,:,:) = 0.0_RP

          tem_pl(:,:,:) = 0.0_RP
          pre_pl(:,:,:) = 0.0_RP
          th_pl (:,:,:) = 0.0_RP
          eth_pl(:,:,:) = 0.0_RP

          pregd_pl(:,:,:) = 0.0_RP
          rhogd_pl(:,:,:) = 0.0_RP

       endif

       !$acc wait
       call DEBUG_rapend  ('___Pre_Post')
       !------------------------------------------------------------------------
       !> LARGE step
       !------------------------------------------------------------------------
       call DEBUG_rapstart('___Large_step')

       !--- calculation of advection tendency including Coriolis force
       call src_advection_convergence_momentum( vx,                     vx_pl,                     & ! [IN]
                                                vy,                     vy_pl,                     & ! [IN]
                                                vz,                     vz_pl,                     & ! [IN]
                                                w,                      w_pl,                      & ! [IN]
                                                PROG  (:,:,:,I_RHOG  ), PROG_pl  (:,:,:,I_RHOG  ), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                                PROG  (:,:,:,I_RHOGW ), PROG_pl  (:,:,:,I_RHOGW ), & ! [IN]
                                                g_TEND(:,:,:,I_RHOGVX), g_TEND_pl(:,:,:,I_RHOGVX), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVY), g_TEND_pl(:,:,:,I_RHOGVY), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVZ), g_TEND_pl(:,:,:,I_RHOGVZ), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGW ), g_TEND_pl(:,:,:,I_RHOGW )  ) ! [OUT]

       !$acc kernels pcopy(g_TEND) async(0)
       g_TEND   (:,:,:,I_RHOG)     = 0.0_RP
       g_TEND   (:,:,:,I_RHOGE)    = 0.0_RP
       g_TEND   (:,:,:,I_RHOGETOT) = 0.0_RP
       !$acc end kernels
       g_TEND_pl(:,:,:,I_RHOG)     = 0.0_RP
       g_TEND_pl(:,:,:,I_RHOGE)    = 0.0_RP
       g_TEND_pl(:,:,:,I_RHOGETOT) = 0.0_RP

       !---< numerical diffusion term
       if ( NDIFF_LOCATION == 'IN_LARGE_STEP' ) then

          if ( nl == 1 ) then ! only first step
             !------ numerical diffusion
             call numfilter_hdiffusion( PROG(:,:,:,I_RHOG), PROG_pl(:,:,:,I_RHOG), & ! [IN]
                                        rho,                rho_pl,                & ! [IN]
                                        vx,                 vx_pl,                 & ! [IN]
                                        vy,                 vy_pl,                 & ! [IN]
                                        vz,                 vz_pl,                 & ! [IN]
                                        w,                  w_pl,                  & ! [IN]
                                        tem,                tem_pl,                & ! [IN]
                                        q,                  q_pl,                  & ! [IN]
                                        f_TEND,             f_TEND_pl,             & ! [OUT]
                                        f_TENDq,            f_TENDq_pl             ) ! [OUT]

             if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
                call numfilter_vdiffusion( PROG(:,:,:,I_RHOG), PROG_pl(:,:,:,I_RHOG), & ! [IN]
                                           rho,                rho_pl,                & ! [IN]
                                           vx,                 vx_pl,                 & ! [IN]
                                           vy,                 vy_pl,                 & ! [IN]
                                           vz,                 vz_pl,                 & ! [IN]
                                           w,                  w_pl,                  & ! [IN]
                                           tem,                tem_pl,                & ! [IN]
                                           q,                  q_pl,                  & ! [IN]
                                           f_TEND,             f_TEND_pl,             & ! [INOUT]
                                           f_TENDq,            f_TENDq_pl             ) ! [INOUT]
             endif

             if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
                call numfilter_rayleigh_damping( PROG  (:,:,:,I_RHOG  ), PROG_pl  (:,:,:,I_RHOG  ), & ! [IN]
                                                 vx,                     vx_pl,                     & ! [IN]
                                                 vy,                     vy_pl,                     & ! [IN]
                                                 vz,                     vz_pl,                     & ! [IN]
                                                 w,                      w_pl,                      & ! [IN]
                                                 f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGW ), f_TEND_pl(:,:,:,I_RHOGW )  ) ! [INOUT]
             endif
          endif

       elseif( NDIFF_LOCATION == 'IN_LARGE_STEP2' ) then

          !------ numerical diffusion
          call numfilter_hdiffusion( PROG(:,:,:,I_RHOG), PROG_pl(:,:,:,I_RHOG), & ! [IN]
                                     rho,                rho_pl,                & ! [IN]
                                     vx,                 vx_pl,                 & ! [IN]
                                     vy,                 vy_pl,                 & ! [IN]
                                     vz,                 vz_pl,                 & ! [IN]
                                     w,                  w_pl,                  & ! [IN]
                                     tem,                tem_pl,                & ! [IN]
                                     q,                  q_pl,                  & ! [IN]
                                     f_TEND,             f_TEND_pl,             & ! [OUT]
                                     f_TENDq,            f_TENDq_pl             ) ! [OUT]

          if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
             call numfilter_vdiffusion( PROG(:,:,:,I_RHOG), PROG_pl(:,:,:,I_RHOG), & ! [IN]
                                        rho,                rho_pl,                & ! [IN]
                                        vx,                 vx_pl,                 & ! [IN]
                                        vy,                 vy_pl,                 & ! [IN]
                                        vz,                 vz_pl,                 & ! [IN]
                                        w,                  w_pl,                  & ! [IN]
                                        tem,                tem_pl,                & ! [IN]
                                        q,                  q_pl,                  & ! [IN]
                                        f_TEND,             f_TEND_pl,             & ! [INOUT]
                                        f_TENDq,            f_TENDq_pl             ) ! [INOUT]
          endif

          if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
             call numfilter_rayleigh_damping( PROG  (:,:,:,I_RHOG  ), PROG_pl  (:,:,:,I_RHOG  ), & ! [IN]
                                              vx,                     vx_pl,                     & ! [IN]
                                              vy,                     vy_pl,                     & ! [IN]
                                              vz,                     vz_pl,                     & ! [IN]
                                              w,                      w_pl,                      & ! [IN]
                                              f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGW ), f_TEND_pl(:,:,:,I_RHOGW )  ) ! [INOUT]
          endif

       endif

!       if ( TB_TYPE == 'SMG' ) then ! Smagorinksy-type SGS model
!          call sgs_smagorinsky( nl,                                                      &
!                                rho,                       rho_pl,                       & ! [IN]
!                                PROG(:,:,:,I_RHOG  ),      PROG_pl(:,:,:,I_RHOG  ),      & ! [IN]
!                                PROGq(:,:,:,:),            PROGq_pl(:,:,:,:  ),          & ! [IN]
!                                vx,                        vx_pl,                        & ! [IN]
!                                vy,                        vy_pl,                        & ! [IN]
!                                vz,                        vz_pl,                        & ! [IN]
!                                w,                         w_pl,                         & ! [IN]
!                                tem,                       tem_pl,                       & ! [IN]
!                                q,                         q_pl,                         & ! [IN]
!                                th,                        th_pl,                        & ! [IN]
!                                f_TEND (:,:,:,I_RHOG    ), f_TEND_pl (:,:,:,I_RHOG    ), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGVX  ), f_TEND_pl (:,:,:,I_RHOGVX  ), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGVY  ), f_TEND_pl (:,:,:,I_RHOGVY  ), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGVZ  ), f_TEND_pl (:,:,:,I_RHOGVZ  ), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGW   ), f_TEND_pl (:,:,:,I_RHOGW   ), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGE   ), f_TEND_pl (:,:,:,I_RHOGE   ), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGETOT), f_TEND_pl (:,:,:,I_RHOGETOT), & ! [INOUT]
!                                f_TENDq(:,:,:,:),          f_TENDq_pl(:,:,:,:)           ) ! [INOUT]
!       endif

       if ( FLAG_NUDGING ) then

          if ( nl == 1 ) then
             call NDG_update_reference( TIME_CTIME )
          endif

          if ( nl == num_of_iteration_lstep ) then
             ndg_TEND_out = .true.
          else
             ndg_TEND_out = .false.
          endif

          call NDG_apply_uvtp( rho,                      rho_pl,                      & ! [IN]
                               vx,                       vx_pl,                       & ! [IN]
                               vy,                       vy_pl,                       & ! [IN]
                               vz,                       vz_pl,                       & ! [IN]
                               w,                        w_pl,                        & ! [IN]
                               tem,                      tem_pl,                      & ! [IN]
                               pre,                      pre_pl,                      & ! [IN]
                               f_TEND(:,:,:,I_RHOGVX  ), f_TEND_pl(:,:,:,I_RHOGVX  ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVY  ), f_TEND_pl(:,:,:,I_RHOGVY  ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVZ  ), f_TEND_pl(:,:,:,I_RHOGVZ  ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGW   ), f_TEND_pl(:,:,:,I_RHOGW   ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGE   ), f_TEND_pl(:,:,:,I_RHOGE   ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGETOT), f_TEND_pl(:,:,:,I_RHOGETOT), & ! [INOUT]
                               ndg_TEND_out                                           ) ! [IN] ( TEND out )
       endif

       !--- sum the large step TEND ( advection + coriolis + num.diff.,SGS,nudge )
       !$acc kernels pcopy(g_TEND) pcopyin(f_TEND) async(0)
       g_TEND   (:,:,:,:) = g_TEND   (:,:,:,:) + f_TEND   (:,:,:,:)
       !$acc end kernels
       g_TEND_pl(:,:,:,:) = g_TEND_pl(:,:,:,:) + f_TEND_pl(:,:,:,:)

       !$acc wait
       call DEBUG_rapend  ('___Large_step')
       !------------------------------------------------------------------------
       !> SMALL step
       !------------------------------------------------------------------------
       call DEBUG_rapstart('___Small_step')

       if ( nl /= 1 ) then ! update split values
          !$acc kernels pcopy(PROG_split) pcopyin(PROG0,PROG) async(0)
          PROG_split   (:,:,:,:) = PROG0   (:,:,:,:) - PROG   (:,:,:,:)
          !$acc end kernels
          PROG_split_pl(:,:,:,:) = PROG0_pl(:,:,:,:) - PROG_pl(:,:,:,:)
       else
          !$acc kernels pcopy(PROG_split) async(0)
          PROG_split   (:,:,:,:) = 0.0_RP
          !$acc end kernels
          PROG_split_pl(:,:,:,:) = 0.0_RP
       endif

       !------ Core routine for small step
       !------    1. By this subroutine, prognostic variables ( rho,.., rhoge ) are calculated through
       !------       the 'num_of_iteration_sstep(nl)'-th times small step.
       !------    2. grho, grhogvx, ..., and  grhoge has the large step
       !------       tendencies initially, however, they are re-used in this subroutine.
       !------
       if ( TIME_SPLIT ) then
          small_step_ite = num_of_iteration_sstep(nl)
          small_step_dt  = TIME_DTS
       else
          small_step_ite = 1
          small_step_dt  = large_step_dt / (num_of_iteration_lstep-nl+1)
       endif

       call vi_small_step( PROG(:,:,:,I_RHOG  ),       PROG_pl(:,:,:,I_RHOG  ),       & ! [INOUT] prognostic variables
                           PROG(:,:,:,I_RHOGVX),       PROG_pl(:,:,:,I_RHOGVX),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGVY),       PROG_pl(:,:,:,I_RHOGVY),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGVZ),       PROG_pl(:,:,:,I_RHOGVZ),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGW ),       PROG_pl(:,:,:,I_RHOGW ),       & ! [INOUT]
                           PROG(:,:,:,I_RHOGE ),       PROG_pl(:,:,:,I_RHOGE ),       & ! [INOUT]
                           vx,                         vx_pl,                         & ! [IN] diagnostic value
                           vy,                         vy_pl,                         & ! [IN]
                           vz,                         vz_pl,                         & ! [IN]
                           eth,                        eth_pl,                        & ! [IN]
                           rhogd,                      rhogd_pl,                      & ! [IN]
                           pregd,                      pregd_pl,                      & ! [IN]
                           g_TEND(:,:,:,I_RHOG    ),   g_TEND_pl(:,:,:,I_RHOG    ),   & ! [IN] large step TEND
                           g_TEND(:,:,:,I_RHOGVX  ),   g_TEND_pl(:,:,:,I_RHOGVX  ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGVY  ),   g_TEND_pl(:,:,:,I_RHOGVY  ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGVZ  ),   g_TEND_pl(:,:,:,I_RHOGVZ  ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGW   ),   g_TEND_pl(:,:,:,I_RHOGW   ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGE   ),   g_TEND_pl(:,:,:,I_RHOGE   ),   & ! [IN]
                           g_TEND(:,:,:,I_RHOGETOT),   g_TEND_pl(:,:,:,I_RHOGETOT),   & ! [IN]
                           PROG_split(:,:,:,I_RHOG  ), PROG_split_pl(:,:,:,I_RHOG  ), & ! [INOUT] split value
                           PROG_split(:,:,:,I_RHOGVX), PROG_split_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGVY), PROG_split_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGVZ), PROG_split_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGW ), PROG_split_pl(:,:,:,I_RHOGW ), & ! [INOUT]
                           PROG_split(:,:,:,I_RHOGE ), PROG_split_pl(:,:,:,I_RHOGE ), & ! [INOUT]
                           PROG_mean,                  PROG_mean_pl,                  & ! [OUT] mean value
                           small_step_ite,                                            & ! [IN]
                           small_step_dt                                              ) ! [IN]

       !$acc wait
       call DEBUG_rapend  ('___Small_step')
       !------------------------------------------------------------------------
       !>  Tracer advection
       !------------------------------------------------------------------------
       call DEBUG_rapstart('___Tracer_Advection')
       do_tke_correction = .false.

       if ( TRC_ADV_TYPE == 'MIURA2004' ) then

          if ( nl == num_of_iteration_lstep ) then

             call src_tracer_advection( TRC_VMAX,                                                & ! [IN]
                                        PROGq    (:,:,:,:),        PROGq_pl    (:,:,:,:),        & ! [INOUT]
                                        PROG0    (:,:,:,I_RHOG  ), PROG0_pl    (:,:,:,I_RHOG  ), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOG  ), PROG_mean_pl(:,:,:,I_RHOG  ), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVX), PROG_mean_pl(:,:,:,I_RHOGVX), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVY), PROG_mean_pl(:,:,:,I_RHOGVY), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVZ), PROG_mean_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGW ), PROG_mean_pl(:,:,:,I_RHOGW ), & ! [IN]
                                        f_TEND   (:,:,:,I_RHOG  ), f_TEND_pl   (:,:,:,I_RHOG  ), & ! [IN]
                                        large_step_dt,                                           & ! [IN]
                                        THUBURN_LIM                                              ) ! [IN]

             !$acc kernels pcopy(PROGq) pcopyin(f_TENDq) async(0)
             PROGq(:,:,:,:) = PROGq(:,:,:,:) + large_step_dt * f_TENDq(:,:,:,:) ! update rhogq by viscosity
             !$acc end kernels

             !$acc kernels pcopy(PROGq) async(0)
             PROGq(:,ADM_kmin-1,:,:) = 0.0_RP
             PROGq(:,ADM_kmax+1,:,:) = 0.0_RP
             !$acc end kernels

             if ( ADM_have_pl ) then
                PROGq_pl(:,:,:,:) = PROGq_pl(:,:,:,:) + large_step_dt * f_TENDq_pl(:,:,:,:)

                PROGq_pl(:,ADM_kmin-1,:,:) = 0.0_RP
                PROGq_pl(:,ADM_kmax+1,:,:) = 0.0_RP
             endif

             ! [comment] H.Tomita: I don't recommend adding the hyperviscosity term because of numerical instability in this case.
             if( I_TKE >= 0 ) do_tke_correction = .true.

          endif ! Last large step only

       elseif( TRC_ADV_TYPE == 'DEFAULT' ) then

          do nq = 1, TRC_VMAX

             call src_advection_convergence( PROG_mean(:,:,:,I_RHOGVX), PROG_mean_pl(:,:,:,I_RHOGVX), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGVY), PROG_mean_pl(:,:,:,I_RHOGVY), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGVZ), PROG_mean_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGW ), PROG_mean_pl(:,:,:,I_RHOGW ), & ! [IN]
                                             q        (:,:,:,nq),       q_pl        (:,:,:,nq),       & ! [IN]
                                             g_TENDq  (:,:,:,nq),       g_TENDq_pl  (:,:,:,nq),       & ! [OUT]
                                             I_SRC_default                                            ) ! [IN]

          enddo ! tracer LOOP

          !$acc kernels pcopy(PROGq) pcopyin(PROGq0,g_TENDq,f_TENDq) async(0)
          PROGq(:,:,:,:) = PROGq0(:,:,:,:)                                                                   &
                         + ( num_of_iteration_sstep(nl) * TIME_DTS ) * ( g_TENDq(:,:,:,:) + f_TENDq(:,:,:,:) )

          PROGq(:,ADM_kmin-1,:,:) = 0.0_RP
          PROGq(:,ADM_kmax+1,:,:) = 0.0_RP
          !$acc end kernels

          if ( ADM_have_pl ) then
             PROGq_pl(:,:,:,:) = PROGq0_pl(:,:,:,:)                                                                      &
                               + ( num_of_iteration_sstep(nl) * TIME_DTS ) * ( g_TENDq_pl(:,:,:,:) + f_TENDq_pl(:,:,:,:) )

             PROGq_pl(:,ADM_kmin-1,:,:) = 0.0_RP
             PROGq_pl(:,ADM_kmax+1,:,:) = 0.0_RP
          endif

          if( I_TKE >= 0 ) do_tke_correction = .true.

       endif

       call DEBUG_rapend  ('___Tracer_Advection')

       call DEBUG_rapstart('___Pre_Post')

       !--- TKE fixer [comment] 2011/08/16 M.Satoh: this fixer is needed for every small time steps
       if ( do_tke_correction ) then
          !$acc kernels pcopy(PROG,PROGq) pcopyin(VMTR_GSGAM2) async(0)
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             TKEG_corr = max( TKE_MIN * VMTR_GSGAM2(g,k,l) - PROGq(g,k,l,I_TKE), 0.0_RP )

             PROG (g,k,l,I_RHOGE) = PROG (g,k,l,I_RHOGE) - TKEG_corr
             PROGq(g,k,l,I_TKE)   = PROGq(g,k,l,I_TKE)   + TKEG_corr
          enddo
          enddo
          enddo
          !$acc end kernels

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                TKEg_corr = max( TKE_MIN * VMTR_GSGAM2_pl(g,k,l) - PROGq_pl(g,k,l,I_TKE), 0.0_RP )

                PROG_pl (g,k,l,I_RHOGE) = PROG_pl (g,k,l,I_RHOGE) - TKEG_corr
                PROGq_pl(g,k,l,I_TKE)   = PROGq_pl(g,k,l,I_TKE)   + TKEG_corr
             enddo
             enddo
             enddo
          endif
       endif

       !------ Update
       if ( nl /= num_of_iteration_lstep ) then
          call COMM_data_transfer( PROG, PROG_pl )

          !$acc kernels pcopy(PROG) async(0)
          do m = 1, 6
          do l = 1, ADM_lall
          do k = 1, ADM_kall
             PROG(suf(ADM_gmax+1,ADM_gmin-1),k,l,m) = PROG(suf(ADM_gmax+1,ADM_gmin),k,l,m)
             PROG(suf(ADM_gmin-1,ADM_gmax+1),k,l,m) = PROG(suf(ADM_gmin,ADM_gmax+1),k,l,m)
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       call DEBUG_rapend  ('___Pre_Post')

    enddo !--- large step

    enddo !--- divided step for dynamics

    call DEBUG_rapstart('___Pre_Post')

    call prgvar_set( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & ! [IN]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & ! [IN]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & ! [IN]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & ! [IN]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & ! [IN]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & ! [IN]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:),       & ! [IN]
                     0                                              ) ! [IN]

    call DEBUG_rapend  ('___Pre_Post')

    !$acc end data

    !$acc wait
    call DEBUG_rapend  ('__Dynamics')

    return
  end subroutine dynamics_step

end module mod_dynamics
