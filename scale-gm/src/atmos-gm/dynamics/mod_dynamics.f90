!-------------------------------------------------------------------------------
!> Module dynamics
!!
!! @par Description
!!          This module contains the core component of fluid dynamics on icosahedral grid system
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_dynamics
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
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
  integer, private :: num_of_iteration_lstep    ! number of large steps ( 0-4 )
  integer, private :: num_of_iteration_sstep(4) ! number of small steps in each of large steps

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine dynamics_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_SSTEP_MAX
    use mod_runconf, only: &
       TRC_ADV_TYPE
    use mod_bndcnd, only: &
       BNDCND_setup
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

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[dynamics]/Category[nhm]'

    if( IO_L ) write(IO_FID_LOG,*) '+++ Time integration type: ', trim(TIME_INTEG_TYPE)
    select case(TIME_INTEG_TYPE)
    case('RK2')
       if( IO_L ) write(IO_FID_LOG,*) '+++ 2-stage Runge-Kutta'

       num_of_iteration_lstep    = 2
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX

    case('RK3')
       if( IO_L ) write(IO_FID_LOG,*) '+++ 3-stage Runge-Kutta'

       num_of_iteration_lstep    = 3
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX

    case('RK4')
       if( IO_L ) write(IO_FID_LOG,*) '+++ 4-stage Runge-Kutta'

       num_of_iteration_lstep    = 4
       num_of_iteration_sstep(1) = TIME_SSTEP_MAX / 4
       num_of_iteration_sstep(2) = TIME_SSTEP_MAX / 3
       num_of_iteration_sstep(3) = TIME_SSTEP_MAX / 2
       num_of_iteration_sstep(4) = TIME_SSTEP_MAX

    case('TRCADV')
       if( IO_L ) write(IO_FID_LOG,*) '+++ Offline tracer experiment'

       num_of_iteration_lstep    = 0

       if ( TRC_ADV_TYPE == 'DEFAULT' ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx unsupported advection scheme for TRCADV test! STOP.'
          call PRC_MPIstop
       endif

    case default
       if( IO_L ) write(IO_FID_LOG,*) 'xxx unsupported integration type! STOP.'
       call PRC_MPIstop
    endselect

    !---< boundary condition module setup >---
    call BNDCND_setup

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
    use scale_const, only: &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       CVdry => CONST_CVdry
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kmax,    &
       ADM_kmin
    use mod_comm, only: &
       COMM_data_transfer
    use mod_vmtr, only: &
       VMTR_getIJ_GSGAM2,    &
       VMTR_getIJ_C2Wfact,   &
       VMTR_getIJ_C2WfactGz, &
       VMTR_getIJ_PHI
    use mod_time, only: &
       TIME_INTEG_TYPE, &
       TIME_DTL,        &
       TIME_DTS,        &
       TIME_SPLIT,      &
       TIME_CTIME
    use mod_runconf, only: &
       I_RHOG,         &
       I_RHOGVX,       &
       I_RHOGVY,       &
       I_RHOGVZ,       &
       I_RHOGW,        &
       I_RHOGE,        &
       I_pre,          &
       I_tem,          &
       I_vx,           &
       I_vy,           &
       I_vz,           &
       I_w,            &
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
       THUBURN_LIM
    use mod_prgvar, only: &
       prgvar_get, &
       prgvar_set
    use mod_bndcnd, only: &
       BNDCND_all
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
       PRG_var
    !##### OpenACC #####
    implicit none

    real(RP) :: PROG         (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! prognostic variables
    real(RP) :: PROG_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: PROGq        (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tracer variables
    real(RP) :: PROGq_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: g_TEND       (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! tendency of prognostic variables
    real(RP) :: g_TEND_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: g_TENDq      (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tendency of tracer variables
    real(RP) :: g_TENDq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: f_TEND       (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! forcing tendency of prognostic variables
    real(RP) :: f_TEND_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: f_TENDq      (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! forcing tendency of tracer variables
    real(RP) :: f_TENDq_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: PROG0        (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! prognostic variables (save)
    real(RP) :: PROG0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: PROGq0       (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tracer variables (save)
    real(RP) :: PROGq0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    real(RP) :: PROG_split   (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! prognostic variables (split)
    real(RP) :: PROG_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,6)

    real(RP) :: PROG_mean    (ADM_gall   ,ADM_kall,ADM_lall   ,5)
    real(RP) :: PROG_mean_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,5)

    real(RP) :: DIAG         (ADM_gall   ,ADM_kall,ADM_lall   ,6)        ! diagnostic variables
    real(RP) :: DIAG_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl,6)
    real(RP) :: q            (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_VMAX) ! tracer variables
    real(RP) :: q_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_VMAX)

    !--- density
    real(RP) :: rho   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- internal energy  ( physical )
    real(RP) :: ein   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ein_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- enthalpy ( physical )
    real(RP) :: eth   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: eth_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- potential temperature ( physical )
    real(RP) :: th   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: th_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- density deviation from the base state ( G^1/2 X gamma2 )
    real(RP) :: rhogd   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- pressure deviation from the base state ( G^1/2 X gamma2 )
    real(RP) :: pregd   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pregd_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    !--- temporary variables
    real(RP) :: qd      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qd_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: cv      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: cv_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: VMTR_GSGAM2      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_GSGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: VMTR_C2Wfact     (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP) :: VMTR_C2Wfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl)
    real(RP) :: VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   )
    real(RP) :: VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl)
    real(RP) :: VMTR_PHI         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_PHI_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: TKEg_corr

    integer  :: small_step_ite
    real(RP) :: large_step_dt
    real(RP) :: small_step_dt

    logical  :: ndg_TEND_out
    logical  :: do_tke_correction

    integer  :: g, k ,l, nq, nl, ndyn
    !---------------------------------------------------------------------------
    !$acc wait

    call PROF_rapstart('__Dynamics',1)
    !$acc  data &
    !$acc& create(PROG,PROGq,g_TEND,g_TENDq,f_TEND,f_TENDq,PROG0,PROGq0,PROG_split,PROG_mean) &
    !$acc& create(rho,vx,vy,vz,w,ein,tem,pre,eth,th,rhogd,pregd,q,qd,cv) &
    !$acc& pcopy(PRG_var)

    call PROF_rapstart('___Pre_Post',1)

    call VMTR_getIJ_GSGAM2   ( VMTR_GSGAM2,    VMTR_GSGAM2_pl    )
    call VMTR_getIJ_C2Wfact  ( VMTR_C2Wfact,   VMTR_C2Wfact_pl   )
    call VMTR_getIJ_C2WfactGz( VMTR_C2WfactGz, VMTR_C2WfactGz_pl )
    call VMTR_getIJ_PHI      ( VMTR_PHI,       VMTR_PHI_pl       )

    large_step_dt = TIME_DTL / real(DYN_DIV_NUM,kind=DP)

    !--- get from prg0
    call prgvar_get( PROG (:,:,:,I_RHOG),   PROG_pl (:,:,:,I_RHOG),   & ! [OUT]
                     PROG (:,:,:,I_RHOGVX), PROG_pl (:,:,:,I_RHOGVX), & ! [OUT]
                     PROG (:,:,:,I_RHOGVY), PROG_pl (:,:,:,I_RHOGVY), & ! [OUT]
                     PROG (:,:,:,I_RHOGVZ), PROG_pl (:,:,:,I_RHOGVZ), & ! [OUT]
                     PROG (:,:,:,I_RHOGW),  PROG_pl (:,:,:,I_RHOGW),  & ! [OUT]
                     PROG (:,:,:,I_RHOGE),  PROG_pl (:,:,:,I_RHOGE),  & ! [OUT]
                     PROGq(:,:,:,:),        PROGq_pl(:,:,:,:)         ) ! [OUT]

    call PROF_rapend  ('___Pre_Post',1)

    do ndyn = 1, DYN_DIV_NUM

    call PROF_rapstart('___Pre_Post',1)

    !--- save
    !$acc kernels pcopy(PROG0) pcopyin(PROG) async(0)
    PROG0(:,:,:,:) = PROG(:,:,:,:)
    !$acc end kernels

    if ( ADM_have_pl ) then
       PROG0_pl(:,:,:,:) = PROG_pl(:,:,:,:)
    endif

    if ( TRC_ADV_TYPE == 'DEFAULT' ) then
       !$acc kernels pcopy(PROGq0) pcopyin(PROGq) async(0)
       PROGq0(:,:,:,:) = PROGq(:,:,:,:)
       !$acc end kernels

       if ( ADM_have_pl ) then
          PROGq0_pl(:,:,:,:) = PROGq_pl(:,:,:,:)
       endif
    endif

    call PROF_rapend  ('___Pre_Post',1)

    if ( TIME_INTEG_TYPE == 'TRCADV' ) then  ! TRC-ADV Test Bifurcation

       call PROF_rapstart('___Tracer_Advection',1)

       !$acc kernels pcopy(f_TEND) async(0)
       f_TEND   (:,:,:,:) = 0.0_RP
       f_TEND_pl(:,:,:,:) = 0.0_RP
       !$acc end kernels

       call src_tracer_advection( TRC_VMAX,                                          & ! [IN]
                                  PROGq (:,:,:,:),        PROGq_pl (:,:,:,:),        & ! [INOUT]
                                  PROG0 (:,:,:,I_RHOG),   PROG0_pl (:,:,:,I_RHOG),   & ! [IN]
                                  PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                  PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                  PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                  PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                  PROG  (:,:,:,I_RHOGW),  PROG_pl  (:,:,:,I_RHOGW),  & ! [IN]
                                  f_TEND(:,:,:,I_RHOG),   f_TEND_pl(:,:,:,I_RHOG),   & ! [IN]
                                  large_step_dt,                                     & ! [IN]
                                  THUBURN_LIM                                        ) ! [IN]

       call PROF_rapend  ('___Tracer_Advection',1)

       call forcing_update( PROG(:,:,:,:), PROG_pl(:,:,:,:) ) ! [INOUT]
    endif

    !---------------------------------------------------------------------------
    !
    !> Start large time step integration
    !
    !---------------------------------------------------------------------------
    do nl = 1, num_of_iteration_lstep

       call PROF_rapstart('___Pre_Post',1)

       !---< Generate diagnostic values and set the boudary conditions
       !$acc kernels pcopy(rho,vx,vy,vz,ein) pcopyin(PROG,VMTR_GSGAM2) async(0)
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          rho (g,k,l)      = PROG(g,k,l,I_RHOG)   / VMTR_GSGAM2(g,k,l)
          DIAG(g,k,l,I_vx) = PROG(g,k,l,I_RHOGVX) / PROG(g,k,l,I_RHOG)
          DIAG(g,k,l,I_vy) = PROG(g,k,l,I_RHOGVY) / PROG(g,k,l,I_RHOG)
          DIAG(g,k,l,I_vz) = PROG(g,k,l,I_RHOGVZ) / PROG(g,k,l,I_RHOG)
          ein (g,k,l)      = PROG(g,k,l,I_RHOGE)  / PROG(g,k,l,I_RHOG)
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

          DIAG(g,k,l,I_tem) = ein(g,k,l) / cv(g,k,l)
          DIAG(g,k,l,I_pre) = rho(g,k,l) * DIAG(g,k,l,I_tem) * ( qd(g,k,l)*Rdry + q(g,k,l,I_QV)*Rvap )
       enddo
       enddo
       enddo
       !$acc end kernels

       !$acc kernels pcopy(w) pcopyin(PROG,VMTR_C2Wfact) async(0)
       do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          DIAG(g,k,l,I_w) = PROG(g,k,l,I_RHOGW) / ( VMTR_C2Wfact(g,k,1,l) * PROG(g,k  ,l,I_RHOG) &
                                                  + VMTR_C2Wfact(g,k,2,l) * PROG(g,k-1,l,I_RHOG) )
       enddo
       enddo
       enddo
       !$acc end kernels

       call BNDCND_all( ADM_gall,                       & ! [IN]
                        ADM_kall,                       & ! [IN]
                        ADM_lall,                       & ! [IN]
                        rho           (:,:,:),          & ! [INOUT]
                        DIAG          (:,:,:,I_vx),     & ! [INOUT]
                        DIAG          (:,:,:,I_vy),     & ! [INOUT]
                        DIAG          (:,:,:,I_vz),     & ! [INOUT]
                        DIAG          (:,:,:,I_w),      & ! [INOUT]
                        ein           (:,:,:),          & ! [INOUT]
                        DIAG          (:,:,:,I_tem),    & ! [INOUT]
                        DIAG          (:,:,:,I_pre),    & ! [INOUT]
                        PROG          (:,:,:,I_RHOG),   & ! [INOUT]
                        PROG          (:,:,:,I_RHOGVX), & ! [INOUT]
                        PROG          (:,:,:,I_RHOGVY), & ! [INOUT]
                        PROG          (:,:,:,I_RHOGVZ), & ! [INOUT]
                        PROG          (:,:,:,I_RHOGW),  & ! [INOUT]
                        PROG          (:,:,:,I_RHOGE),  & ! [INOUT]
                        VMTR_GSGAM2   (:,:,:),          & ! [IN]
                        VMTR_PHI      (:,:,:),          & ! [IN]
                        VMTR_C2Wfact  (:,:,:,:),        & ! [IN]
                        VMTR_C2WfactGz(:,:,:,:)         ) ! [IN]

       call THRMDYN_th ( ADM_gall,          & ! [IN]
                         ADM_kall,          & ! [IN]
                         ADM_lall,          & ! [IN]
                         DIAG(:,:,:,I_tem), & ! [IN]
                         DIAG(:,:,:,I_pre), & ! [IN]
                         th  (:,:,:)        ) ! [OUT]

       call THRMDYN_eth( ADM_gall,          & ! [IN]
                         ADM_kall,          & ! [IN]
                         ADM_lall,          & ! [IN]
                         ein (:,:,:),       & ! [IN]
                         DIAG(:,:,:,I_pre), & ! [IN]
                         rho (:,:,:),       & ! [IN]
                         eth (:,:,:)        ) ! [OUT]

       ! perturbations ( pre, rho with metrics )
       !$acc  kernels pcopy(pregd,rhogd) pcopyin(pre,pre_bs,rho,rho_bs,VMTR_GSGAM2) async(0)
       pregd(:,:,:) = ( DIAG(:,:,:,I_pre) - pre_bs(:,:,:) ) * VMTR_GSGAM2(:,:,:)
       rhogd(:,:,:) = ( rho (:,:,:)       - rho_bs(:,:,:) ) * VMTR_GSGAM2(:,:,:)
       !$acc end kernels

       if ( ADM_have_pl ) then

          rho_pl (:,:,:)      = PROG_pl(:,:,:,I_RHOG)   / VMTR_GSGAM2_pl(:,:,:)
          DIAG_pl(:,:,:,I_vx) = PROG_pl(:,:,:,I_RHOGVX) / PROG_pl(:,:,:,I_RHOG)
          DIAG_pl(:,:,:,I_vy) = PROG_pl(:,:,:,I_RHOGVY) / PROG_pl(:,:,:,I_RHOG)
          DIAG_pl(:,:,:,I_vz) = PROG_pl(:,:,:,I_RHOGVZ) / PROG_pl(:,:,:,I_RHOG)
          ein_pl (:,:,:)      = PROG_pl(:,:,:,I_RHOGE)  / PROG_pl(:,:,:,I_RHOG)

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

          DIAG_pl(:,:,:,I_tem) = ein_pl(:,:,:) / cv_pl(:,:,:)
          DIAG_pl(:,:,:,I_pre) = rho_pl(:,:,:) * DIAG_pl(:,:,:,I_tem) * ( qd_pl(:,:,:)*Rdry + q_pl(:,:,:,I_QV)*Rvap )

          do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             DIAG_pl(g,k,l,I_w) = PROG_pl(g,k,l,I_RHOGW) / ( VMTR_C2Wfact_pl(g,k,1,l) * PROG_pl(g,k  ,l,I_RHOG) &
                                                           + VMTR_C2Wfact_pl(g,k,2,l) * PROG_pl(g,k-1,l,I_RHOG) )
          enddo
          enddo
          enddo

          call BNDCND_all( ADM_gall_pl,                       & ! [IN]
                           ADM_kall,                          & ! [IN]
                           ADM_lall_pl,                       & ! [IN]
                           rho_pl           (:,:,:),          & ! [INOUT]
                           DIAG_pl          (:,:,:,I_vx),     & ! [INOUT]
                           DIAG_pl          (:,:,:,I_vy),     & ! [INOUT]
                           DIAG_pl          (:,:,:,I_vz),     & ! [INOUT]
                           DIAG_pl          (:,:,:,I_w),      & ! [INOUT]
                           ein_pl           (:,:,:),          & ! [INOUT]
                           DIAG_pl          (:,:,:,I_tem),    & ! [INOUT]
                           DIAG_pl          (:,:,:,I_pre),    & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOG),   & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGVX), & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGVY), & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGVZ), & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGW),  & ! [INOUT]
                           PROG_pl          (:,:,:,I_RHOGE),  & ! [INOUT]
                           VMTR_GSGAM2_pl   (:,:,:),          & ! [IN]
                           VMTR_PHI_pl      (:,:,:),          & ! [IN]
                           VMTR_C2Wfact_pl  (:,:,:,:),        & ! [IN]
                           VMTR_C2WfactGz_pl(:,:,:,:)         ) ! [IN]

          call THRMDYN_th ( ADM_gall_pl,          & ! [IN]
                            ADM_kall,             & ! [IN]
                            ADM_lall_pl,          & ! [IN]
                            DIAG_pl(:,:,:,I_tem), & ! [IN]
                            DIAG_pl(:,:,:,I_pre), & ! [IN]
                            th_pl  (:,:,:)        ) ! [OUT]

          call THRMDYN_eth( ADM_gall_pl,          & ! [IN]
                            ADM_kall,             & ! [IN]
                            ADM_lall_pl,          & ! [IN]
                            ein_pl (:,:,:),       & ! [IN]
                            DIAG_pl(:,:,:,I_pre), & ! [IN]
                            rho_pl (:,:,:),       & ! [IN]
                            eth_pl (:,:,:)        ) ! [OUT]

          pregd_pl(:,:,:) = ( DIAG_pl(:,:,:,I_pre) - pre_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
          rhogd_pl(:,:,:) = ( rho_pl (:,:,:)       - rho_bs_pl(:,:,:) ) * VMTR_GSGAM2_pl(:,:,:)
       else

          DIAG_pl (:,:,:,:) = 0.0_RP

          th_pl   (:,:,:) = 0.0_RP
          eth_pl  (:,:,:) = 0.0_RP
          pregd_pl(:,:,:) = 0.0_RP
          rhogd_pl(:,:,:) = 0.0_RP

       endif

       !$acc wait
       call PROF_rapend  ('___Pre_Post',1)
       !------------------------------------------------------------------------
       !> LARGE step
       !------------------------------------------------------------------------
       call PROF_rapstart('___Large_step',1)

       !--- calculation of advection tendency including Coriolis force
       call src_advection_convergence_momentum( DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                                                DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                                                DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                                                DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                                                PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                                PROG  (:,:,:,I_RHOGVX), PROG_pl  (:,:,:,I_RHOGVX), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVY), PROG_pl  (:,:,:,I_RHOGVY), & ! [IN]
                                                PROG  (:,:,:,I_RHOGVZ), PROG_pl  (:,:,:,I_RHOGVZ), & ! [IN]
                                                PROG  (:,:,:,I_RHOGW),  PROG_pl  (:,:,:,I_RHOGW),  & ! [IN]
                                                g_TEND(:,:,:,I_RHOGVX), g_TEND_pl(:,:,:,I_RHOGVX), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVY), g_TEND_pl(:,:,:,I_RHOGVY), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGVZ), g_TEND_pl(:,:,:,I_RHOGVZ), & ! [OUT]
                                                g_TEND(:,:,:,I_RHOGW),  g_TEND_pl(:,:,:,I_RHOGW)   ) ! [OUT]

       !$acc kernels pcopy(g_TEND) async(0)
       g_TEND   (:,:,:,I_RHOG)  = 0.0_RP
       g_TEND   (:,:,:,I_RHOGE) = 0.0_RP
       !$acc end kernels
       g_TEND_pl(:,:,:,I_RHOG)  = 0.0_RP
       g_TEND_pl(:,:,:,I_RHOGE) = 0.0_RP

       !---< numerical diffusion term
       if ( NDIFF_LOCATION == 'IN_LARGE_STEP' ) then

          if ( nl == 1 ) then ! only first step
             !------ numerical diffusion
             call numfilter_hdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                        rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                        DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                        DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                        DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                        DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                        DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                        q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                        f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [OUT]
                                        f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [OUT]

             if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
                call numfilter_vdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                           rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                           DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                           DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                           DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                           DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                           DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                           q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                           f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [INOUT]
                                           f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [INOUT]
             endif

             if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
                call numfilter_rayleigh_damping( PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                                 DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                                                 DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                                                 DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                                                 DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                                                 f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                                                 f_TEND(:,:,:,I_RHOGW),  f_TEND_pl(:,:,:,I_RHOGW)   ) ! [INOUT]
             endif
          endif

       elseif( NDIFF_LOCATION == 'IN_LARGE_STEP2' ) then

          !------ numerical diffusion
          call numfilter_hdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                     rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                     DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                     DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                     DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                     DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                     DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                     q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                     f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [OUT]
                                     f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [OUT]

          if ( NUMFILTER_DOverticaldiff ) then ! numerical diffusion (vertical)
             call numfilter_vdiffusion( PROG   (:,:,:,I_RHOG), PROG_pl   (:,:,:,I_RHOG), & ! [IN]
                                        rho    (:,:,:),        rho_pl    (:,:,:),        & ! [IN]
                                        DIAG   (:,:,:,I_vx),   DIAG_pl   (:,:,:,I_vx),   & ! [IN]
                                        DIAG   (:,:,:,I_vy),   DIAG_pl   (:,:,:,I_vy),   & ! [IN]
                                        DIAG   (:,:,:,I_vz),   DIAG_pl   (:,:,:,I_vz),   & ! [IN]
                                        DIAG   (:,:,:,I_w),    DIAG_pl   (:,:,:,I_w),    & ! [IN]
                                        DIAG   (:,:,:,I_tem),  DIAG_pl   (:,:,:,I_tem),  & ! [IN]
                                        q      (:,:,:,:),      q_pl      (:,:,:,:),      & ! [IN]
                                        f_TEND (:,:,:,:),      f_TEND_pl (:,:,:,:),      & ! [INOUT]
                                        f_TENDq(:,:,:,:),      f_TENDq_pl(:,:,:,:)       ) ! [INOUT]
          endif

          if ( NUMFILTER_DOrayleigh ) then ! rayleigh damping
             call numfilter_rayleigh_damping( PROG  (:,:,:,I_RHOG),   PROG_pl  (:,:,:,I_RHOG),   & ! [IN]
                                              DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                                              DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                                              DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                                              DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                                              f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                                              f_TEND(:,:,:,I_RHOGW),  f_TEND_pl(:,:,:,I_RHOGW)   ) ! [INOUT]
          endif

       endif

!       if ( TB_TYPE == 'SMG' ) then ! Smagorinksy-type SGS model
!          call sgs_smagorinsky( nl,                                                  &
!                                rho    (:,:,:),          rho_pl    (:,:,:),          & ! [IN]
!                                PROG   (:,:,:,I_RHOG),   PROG_pl   (:,:,:,I_RHOG),   & ! [IN]
!                                PROGq  (:,:,:,:),        PROGq_pl  (:,:,:,:  ),      & ! [IN]
!                                DIAG   (:,:,:,I_vx),     DIAG_pl   (:,:,:,I_vx),     & ! [IN]
!                                DIAG   (:,:,:,I_vy),     DIAG_pl   (:,:,:,I_vy),     & ! [IN]
!                                DIAG   (:,:,:,I_vz),     DIAG_pl   (:,:,:,I_vz),     & ! [IN]
!                                DIAG   (:,:,:,I_w),      DIAG_pl   (:,:,:,I_w),      & ! [IN]
!                                DIAG   (:,:,:,I_tem),    DIAG_pl   (:,:,:,I_tem),    & ! [IN]
!                                q      (:,:,:,:),        q_pl      (:,:,:,:),        & ! [IN]
!                                th     (:,:,:),          th_pl     (:,:,:),          & ! [IN]
!                                f_TEND (:,:,:,I_RHOG),   f_TEND_pl (:,:,:,I_RHOG),   & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGVX), f_TEND_pl (:,:,:,I_RHOGVX), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGVY), f_TEND_pl (:,:,:,I_RHOGVY), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGVZ), f_TEND_pl (:,:,:,I_RHOGVZ), & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGW),  f_TEND_pl (:,:,:,I_RHOGW),  & ! [INOUT]
!                                f_TEND (:,:,:,I_RHOGE),  f_TEND_pl (:,:,:,I_RHOGE),  & ! [INOUT]
!                                f_TENDq(:,:,:,:),        f_TENDq_pl(:,:,:,:)         ) ! [INOUT]
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

          call NDG_apply_uvtp( rho   (:,:,:),          rho_pl   (:,:,:),          & ! [IN]
                               DIAG  (:,:,:,I_vx),     DIAG_pl  (:,:,:,I_vx),     & ! [IN]
                               DIAG  (:,:,:,I_vy),     DIAG_pl  (:,:,:,I_vy),     & ! [IN]
                               DIAG  (:,:,:,I_vz),     DIAG_pl  (:,:,:,I_vz),     & ! [IN]
                               DIAG  (:,:,:,I_w),      DIAG_pl  (:,:,:,I_w),      & ! [IN]
                               DIAG  (:,:,:,I_tem),    DIAG_pl  (:,:,:,I_tem),    & ! [IN]
                               DIAG  (:,:,:,I_pre),    DIAG_pl  (:,:,:,I_pre),    & ! [IN]
                               f_TEND(:,:,:,I_RHOGVX), f_TEND_pl(:,:,:,I_RHOGVX), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVY), f_TEND_pl(:,:,:,I_RHOGVY), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGVZ), f_TEND_pl(:,:,:,I_RHOGVZ), & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGW),  f_TEND_pl(:,:,:,I_RHOGW),  & ! [INOUT]
                               f_TEND(:,:,:,I_RHOGE),  f_TEND_pl(:,:,:,I_RHOGE),  & ! [INOUT]
                               ndg_TEND_out                                       ) ! [IN] ( TEND out )
       endif

       !--- sum the large step TEND ( advection + coriolis + num.diff.,SGS,nudge )
       !$acc kernels pcopy(g_TEND) pcopyin(f_TEND) async(0)
       g_TEND   (:,:,:,:) = g_TEND   (:,:,:,:) + f_TEND   (:,:,:,:)
       !$acc end kernels
       g_TEND_pl(:,:,:,:) = g_TEND_pl(:,:,:,:) + f_TEND_pl(:,:,:,:)

       !$acc wait
       call PROF_rapend  ('___Large_step',1)
       !------------------------------------------------------------------------
       !> SMALL step
       !------------------------------------------------------------------------
       call PROF_rapstart('___Small_step',1)

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

       call vi_small_step( PROG      (:,:,:,:),    PROG_pl      (:,:,:,:),    & ! [INOUT] prognostic variables
                           DIAG      (:,:,:,I_vx), DIAG_pl      (:,:,:,I_vx), & ! [IN] diagnostic value
                           DIAG      (:,:,:,I_vy), DIAG_pl      (:,:,:,I_vy), & ! [IN]
                           DIAG      (:,:,:,I_vz), DIAG_pl      (:,:,:,I_vz), & ! [IN]
                           eth       (:,:,:),      eth_pl       (:,:,:),      & ! [IN]
                           rhogd     (:,:,:),      rhogd_pl     (:,:,:),      & ! [IN]
                           pregd     (:,:,:),      pregd_pl     (:,:,:),      & ! [IN]
                           g_TEND    (:,:,:,:),    g_TEND_pl    (:,:,:,:),    & ! [IN] large step TEND
                           PROG_split(:,:,:,:),    PROG_split_pl(:,:,:,:),    & ! [INOUT] split value
                           PROG_mean (:,:,:,:),    PROG_mean_pl(:,:,:,:),     & ! [OUT] mean value
                           small_step_ite,                                    & ! [IN]
                           small_step_dt                                      ) ! [IN]

       !$acc wait
       call PROF_rapend  ('___Small_step',1)
       !------------------------------------------------------------------------
       !>  Tracer advection
       !------------------------------------------------------------------------
       call PROF_rapstart('___Tracer_Advection',1)
       do_tke_correction = .false.

       if ( TRC_ADV_TYPE == 'MIURA2004' ) then

          if ( nl == num_of_iteration_lstep ) then

             call src_tracer_advection( TRC_VMAX,                                                & ! [IN]
                                        PROGq    (:,:,:,:),        PROGq_pl    (:,:,:,:),        & ! [INOUT]
                                        PROG0    (:,:,:,I_RHOG),   PROG0_pl    (:,:,:,I_RHOG),   & ! [IN]
                                        PROG_mean(:,:,:,I_RHOG),   PROG_mean_pl(:,:,:,I_RHOG),   & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVX), PROG_mean_pl(:,:,:,I_RHOGVX), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVY), PROG_mean_pl(:,:,:,I_RHOGVY), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGVZ), PROG_mean_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                        PROG_mean(:,:,:,I_RHOGW),  PROG_mean_pl(:,:,:,I_RHOGW),  & ! [IN]
                                        f_TEND   (:,:,:,I_RHOG),   f_TEND_pl   (:,:,:,I_RHOG),   & ! [IN]
                                        large_step_dt,                                           & ! [IN]
                                        THUBURN_LIM                                              ) ! [IN]

             !$acc kernels pcopy(PROGq) pcopyin(f_TENDq) async(0)
             PROGq(:,:,:,:) = PROGq(:,:,:,:) + large_step_dt * f_TENDq(:,:,:,:) ! update rhogq by viscosity
             !$acc end kernels

             if ( ADM_have_pl ) then
                PROGq_pl(:,:,:,:) = PROGq_pl(:,:,:,:) + large_step_dt * f_TENDq_pl(:,:,:,:)
             endif

             ! [comment] H.Tomita: I don't recommend adding the hyperviscosity term because of numerical instability in this case.
             if( I_TKE >= 0 ) do_tke_correction = .true.

          endif ! Last large step only

       elseif( TRC_ADV_TYPE == 'DEFAULT' ) then

          do nq = 1, TRC_VMAX

             call src_advection_convergence( PROG_mean(:,:,:,I_RHOGVX), PROG_mean_pl(:,:,:,I_RHOGVX), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGVY), PROG_mean_pl(:,:,:,I_RHOGVY), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGVZ), PROG_mean_pl(:,:,:,I_RHOGVZ), & ! [IN]
                                             PROG_mean(:,:,:,I_RHOGW),  PROG_mean_pl(:,:,:,I_RHOGW),  & ! [IN]
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

       call PROF_rapend  ('___Tracer_Advection',1)

       call PROF_rapstart('___Pre_Post',1)

       ! TKE fixer
       if ( do_tke_correction ) then
          !$acc kernels pcopy(PROG,PROGq) pcopyin(VMTR_GSGAM2) async(0)
          do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             TKEG_corr = max( -PROGq(g,k,l,I_TKE), 0.0_RP )

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
                TKEg_corr = max( -PROGq_pl(g,k,l,I_TKE), 0.0_RP )

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
       endif

       call PROF_rapend  ('___Pre_Post',1)

    enddo !--- large step

    enddo !--- divided step for dynamics

    call PROF_rapstart('___Pre_Post',1)

    call prgvar_set( PROG(:,:,:,I_RHOG),   PROG_pl(:,:,:,I_RHOG),   & ! [IN]
                     PROG(:,:,:,I_RHOGVX), PROG_pl(:,:,:,I_RHOGVX), & ! [IN]
                     PROG(:,:,:,I_RHOGVY), PROG_pl(:,:,:,I_RHOGVY), & ! [IN]
                     PROG(:,:,:,I_RHOGVZ), PROG_pl(:,:,:,I_RHOGVZ), & ! [IN]
                     PROG(:,:,:,I_RHOGW),  PROG_pl(:,:,:,I_RHOGW),  & ! [IN]
                     PROG(:,:,:,I_RHOGE),  PROG_pl(:,:,:,I_RHOGE),  & ! [IN]
                     PROGq(:,:,:,:),       PROGq_pl(:,:,:,:)        ) ! [IN]

    call PROF_rapend  ('___Pre_Post',1)

    !$acc end data

    !$acc wait
    call PROF_rapend  ('__Dynamics',1)

    return
  end subroutine dynamics_step

end module mod_dynamics
