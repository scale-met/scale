!-------------------------------------------------------------------------------
!> Module history variables
!!
!! @par Description
!!          This module prepares diagnostic/prognostic variables for history output
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_history_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: history_vars_setup
  public :: history_vars

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
  logical, private :: out_uv_cos   = .false.
  logical, private :: out_omg      = .false.
  logical, private :: out_rha      = .false.
  logical, private :: out_rh       = .false.
  logical, private :: out_rhi      = .false.
  logical, private :: out_th       = .false.
  logical, private :: out_th_prime = .false.
  logical, private :: out_mse      = .false.
  logical, private :: out_850hPa   = .false.
  logical, private :: out_500hPa   = .false.
  logical, private :: out_250hPa   = .false.
  logical, private :: out_100hPa   = .false.

  logical, private :: out_pw       = .false.
  logical, private :: out_lwp      = .false.
  logical, private :: out_iwp      = .false.
  logical, private :: out_duvw     = .false.
  logical, private :: out_dtem     = .false.
  logical, private :: out_dq       = .false.

  real(RP), private, allocatable :: u_old  (:,:,:)
  real(RP), private, allocatable :: v_old  (:,:,:)
  real(RP), private, allocatable :: wc_old (:,:,:)
  real(RP), private, allocatable :: tem_old(:,:,:)
  real(RP), private, allocatable :: qv_old (:,:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine history_vars_setup
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kmin,    &
       ADM_kmax
    use mod_vmtr, only: &
       VMTR_getIJ_GSGAM2,    &
       VMTR_getIJ_W2Cfact,   &
       VMTR_getIJ_C2Wfact,   &
       VMTR_getIJ_C2WfactGz, &
       VMTR_getIJ_PHI
    use mod_runconf, only: &
       TRC_VMAX, &
       I_QV,     &
       AF_TYPE
    use mod_prgvar, only: &
       prgvar_get_withdiag
    use mod_bndcnd, only: &
       BNDCND_all
    use mod_cnvvar, only: &
       cnvvar_vh2uv
    use mod_history, only: &
       HIST_req_nmax, &
       item_save,     &
       history_in
    implicit none

    real(RP) :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhoge    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogq    (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: pre      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tem      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: w        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q        (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(RP) :: u        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: u_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: v        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: v_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: wc       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ein      (ADM_gall   ,ADM_kall,ADM_lall   )

    real(RP) :: tmp3d(ADM_gall,ADM_kall)
    real(RP) :: tmp2d(ADM_gall,ADM_KNONE)

    real(RP) :: VMTR_GSGAM2      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_GSGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: VMTR_W2Cfact     (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP) :: VMTR_W2Cfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl)
    real(RP) :: VMTR_C2Wfact     (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP) :: VMTR_C2Wfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl)
    real(RP) :: VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   )
    real(RP) :: VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl)
    real(RP) :: VMTR_PHI         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_PHI_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    character(len=H_SHORT) :: varname

    integer  :: n, k, l, nq
    !---------------------------------------------------------------------------

    call VMTR_getIJ_GSGAM2   ( VMTR_GSGAM2,    VMTR_GSGAM2_pl    )
    call VMTR_getIJ_W2Cfact  ( VMTR_W2Cfact,   VMTR_W2Cfact_pl   )
    call VMTR_getIJ_C2Wfact  ( VMTR_C2Wfact,   VMTR_C2Wfact_pl   )
    call VMTR_getIJ_C2WfactGz( VMTR_C2WfactGz, VMTR_C2WfactGz_pl )
    call VMTR_getIJ_PHI      ( VMTR_PHI,       VMTR_PHI_pl       )

    do n = 1, HIST_req_nmax
       if(      item_save(n) == 'ml_ucos'     &
           .OR. item_save(n) == 'ml_vcos'     ) out_uv_cos   = .true.
       if(      item_save(n) == 'ml_omg'      ) out_omg      = .true.
       if(      item_save(n) == 'ml_rha'      ) out_rha      = .true.
       if(      item_save(n) == 'ml_rh'       ) out_rh       = .true.
       if(      item_save(n) == 'ml_rhi'      ) out_rhi      = .true.
       if(      item_save(n) == 'ml_th'       &
           .OR. item_save(n) == 'ml_thv'      ) out_th       = .true.
       if(      item_save(n) == 'ml_th_prime' ) out_th_prime = .true.
       if(      item_save(n) == 'ml_mse'      ) out_mse      = .true.
       if(      item_save(n) == 'sl_u850'     &
           .OR. item_save(n) == 'sl_v850'     &
           .OR. item_save(n) == 'sl_w850'     &
           .OR. item_save(n) == 'sl_t850'     ) out_850hPa   = .true.
       if(      item_save(n) == 'sl_u500'     &
           .OR. item_save(n) == 'sl_v500'     &
           .OR. item_save(n) == 'sl_w500'     &
           .OR. item_save(n) == 'sl_t500'     ) out_500hPa   = .true.
       if(      item_save(n) == 'sl_u250'     &
           .OR. item_save(n) == 'sl_v250'     &
           .OR. item_save(n) == 'sl_w250'     &
           .OR. item_save(n) == 'sl_t250'     ) out_250hPa   = .true.
       if(      item_save(n) == 'sl_u100'     &
           .OR. item_save(n) == 'sl_v100'     &
           .OR. item_save(n) == 'sl_w100'     &
           .OR. item_save(n) == 'sl_t100'     ) out_100hPa   = .true.

       if(      item_save(n) == 'sl_pw'       ) out_pw       = .true.
       if(      item_save(n) == 'sl_lwp'      ) out_lwp      = .true.
       if(      item_save(n) == 'sl_iwp'      ) out_iwp      = .true.

       if(      item_save(n) == 'ml_du'       &
           .OR. item_save(n) == 'ml_dv'       &
           .OR. item_save(n) == 'ml_dw'       ) out_duvw     = .true.
       if(      item_save(n) == 'ml_dtem'     ) out_dtem     = .true.
       if(      item_save(n) == 'ml_dq'       ) out_dq       = .true.

    enddo

    !--- get variables
    call prgvar_get_withdiag( rhog,   rhog_pl,   & ! [OUT]
                              rhogvx, rhogvx_pl, & ! [OUT]
                              rhogvy, rhogvy_pl, & ! [OUT]
                              rhogvz, rhogvz_pl, & ! [OUT]
                              rhogw,  rhogw_pl,  & ! [OUT]
                              rhoge,  rhoge_pl,  & ! [OUT]
                              rhogq,  rhogq_pl,  & ! [OUT]
                              rho,    rho_pl,    & ! [OUT]
                              pre,    pre_pl,    & ! [OUT]
                              tem,    tem_pl,    & ! [OUT]
                              vx,     vx_pl,     & ! [OUT]
                              vy,     vy_pl,     & ! [OUT]
                              vz,     vz_pl,     & ! [OUT]
                              w,      w_pl,      & ! [OUT]
                              q,      q_pl       ) ! [OUT]

    ein(:,:,:) = rhoge(:,:,:) / rhog(:,:,:)

    ! boundary condition
    call BNDCND_all( ADM_gall,                & ! [IN]
                     ADM_kall,                & ! [IN]
                     ADM_lall,                & ! [IN]
                     rho           (:,:,:),   & ! [INOUT]
                     vx            (:,:,:),   & ! [INOUT]
                     vy            (:,:,:),   & ! [INOUT]
                     vz            (:,:,:),   & ! [INOUT]
                     w             (:,:,:),   & ! [INOUT]
                     ein           (:,:,:),   & ! [INOUT]
                     tem           (:,:,:),   & ! [INOUT]
                     pre           (:,:,:),   & ! [INOUT]
                     rhog          (:,:,:),   & ! [INOUT]
                     rhogvx        (:,:,:),   & ! [INOUT]
                     rhogvy        (:,:,:),   & ! [INOUT]
                     rhogvz        (:,:,:),   & ! [INOUT]
                     rhogw         (:,:,:),   & ! [INOUT]
                     rhoge         (:,:,:),   & ! [INOUT]
                     VMTR_GSGAM2   (:,:,:),   & ! [IN]
                     VMTR_PHI      (:,:,:),   & ! [IN]
                     VMTR_C2Wfact  (:,:,:,:), & ! [IN]
                     VMTR_C2WfactGz(:,:,:,:)  ) ! [IN]

    ! zonal and meridonal wind
    call cnvvar_vh2uv( u (:,:,:), u_pl (:,:,:), & ! [OUT]
                       v (:,:,:), v_pl (:,:,:), & ! [OUT]
                       vx(:,:,:), vx_pl(:,:,:), & ! [IN]
                       vy(:,:,:), vy_pl(:,:,:), & ! [IN]
                       vz(:,:,:), vz_pl(:,:,:), & ! [IN]
                       withcos = .false.        ) ! [IN]

    ! vertical wind at cell center
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          wc(:,k,l) = ( VMTR_W2Cfact(:,k,1,l) * w(:,k+1,l) &
                      + VMTR_W2Cfact(:,k,2,l) * w(:,k  ,l) )
       enddo
       wc(:,ADM_kmin-1,l) = w(:,ADM_kmin  ,l)
       wc(:,ADM_kmax+1,l) = w(:,ADM_kmax+1,l)
    enddo

    if (out_duvw) then
       allocate( u_old (ADM_gall,ADM_kall,ADM_lall) )
       allocate( v_old (ADM_gall,ADM_kall,ADM_lall) )
       allocate( wc_old(ADM_gall,ADM_kall,ADM_lall) )

       u_old (:,:,:) = u (:,:,:)
       v_old (:,:,:) = v (:,:,:)
       wc_old(:,:,:) = wc(:,:,:)
    endif

    if (out_dtem) then
       allocate( tem_old(ADM_gall,ADM_kall,ADM_lall) )

       tem_old(:,:,:) = tem(:,:,:)
    endif

    if (out_dq) then
       allocate( qv_old(ADM_gall,ADM_kall,ADM_lall) )

       qv_old(:,:,:) = q(:,:,:,I_QV)
    endif

    tmp2d(:,:) = 0.0_RP
    tmp3d(:,:) = 0.0_RP

    select case(AF_TYPE)
    case('HELD-SUAREZ')
       do l = 1, ADM_lall
          call history_in( 'ml_af_fvx', tmp3d(:,:) )
          call history_in( 'ml_af_fvy', tmp3d(:,:) )
          call history_in( 'ml_af_fvz', tmp3d(:,:) )
          call history_in( 'ml_af_fe',  tmp3d(:,:) )
       enddo
    case('DCMIP')
       do l = 1, ADM_lall
          call history_in( 'ml_af_fvx', tmp3d(:,:) )
          call history_in( 'ml_af_fvy', tmp3d(:,:) )
          call history_in( 'ml_af_fvz', tmp3d(:,:) )
          call history_in( 'ml_af_fe',  tmp3d(:,:) )

          do nq = 1, TRC_VMAX
             write(varname,'(A,I2.2)') 'ml_af_fq', nq
             call history_in( varname, tmp3d(:,:) )
          enddo

          call history_in( 'sl_af_prcp', tmp2d(:,:) )
       enddo
    end select

    return
  end subroutine history_vars_setup

  !----------------------------------------------------------------------
  subroutine history_vars
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       GRAV  => CONST_GRAV,  &
       CPdry => CONST_CPdry, &
       Rvap  => CONST_Rvap
    use scale_atmos_hydrometer, only: &
       LHV
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_dgz,  &
       GRD_ZSFC, &
       GRD_Z,    &
       GRD_zs,   &
       GRD_vz
    use mod_vmtr, only: &
       VMTR_getIJ_GSGAM2,    &
       VMTR_getIJ_W2Cfact,   &
       VMTR_getIJ_C2Wfact,   &
       VMTR_getIJ_C2WfactGz, &
       VMTR_getIJ_PHI
    use mod_gm_statistics, only: &
       GTL_global_sum_eachlayer
    use mod_time, only: &
       TIME_DTL
    use mod_runconf, only: &
       TRC_VMAX,  &
       TRC_name,  &
       NQW_STR,   &
       NQW_END,   &
       I_QV,      &
       I_QC,      &
       I_QR,      &
       I_QI,      &
       I_QS,      &
       I_QG,      &
       AF_TYPE,   &
       NCHEM_STR, &
       NCHEM_END
    use mod_prgvar, only: &
       prgvar_get_withdiag
    use mod_thrmdyn, only: &
       THRMDYN_th
    use mod_saturation, only: &
       SATURATION_psat_all, &
       SATURATION_psat_liq, &
       SATURATION_psat_ice
    use mod_bndcnd, only: &
       BNDCND_all
    use mod_cnvvar, only: &
       cnvvar_vh2uv
    use mod_history, only: &
       history_in
    implicit none

    real(RP) :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhoge    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhoge_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogq    (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: rhogq_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)
    real(RP) :: rho      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rho_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: pre      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: pre_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: tem      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tem_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: w        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: q        (ADM_gall   ,ADM_kall,ADM_lall   ,TRC_vmax)
    real(RP) :: q_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl,TRC_vmax)

    real(RP) :: u        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: u_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: v        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: v_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ucos     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vcos     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: wc       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ein      (ADM_gall   ,ADM_kall,ADM_lall   )

    real(RP) :: omg      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: psat     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rh       (ADM_gall   ,ADM_kall,ADM_lall   )

    real(RP) :: u_slice  (ADM_gall   ,ADM_KNONE,ADM_lall  )
    real(RP) :: v_slice  (ADM_gall   ,ADM_KNONE,ADM_lall  )
    real(RP) :: w_slice  (ADM_gall   ,ADM_KNONE,ADM_lall  )
    real(RP) :: t_slice  (ADM_gall   ,ADM_KNONE,ADM_lall  )
    real(RP) :: rho_sfc  (ADM_gall   ,ADM_KNONE,ADM_lall  )
    real(RP) :: pre_sfc  (ADM_gall   ,ADM_KNONE,ADM_lall  )

    real(RP) :: th_prime (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: one      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: one_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: area_prof(ADM_kall)
    real(RP) :: th       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: th_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: th_prof  (ADM_kall)
    real(RP) :: thv      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: mse      (ADM_gall   ,ADM_kall,ADM_lall   )

    real(RP) :: q_clw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: q_cli    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: qtot     (ADM_gall   ,ADM_kall,ADM_lall   )

    real(RP) :: tmp2d    (ADM_gall   ,ADM_KNONE,ADM_lall  )
    real(RP) :: rhodz    (ADM_gall   ,ADM_KNONE,ADM_lall  )

    real(RP) :: VMTR_GSGAM2      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_GSGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: VMTR_W2Cfact     (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP) :: VMTR_W2Cfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl)
    real(RP) :: VMTR_C2Wfact     (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP) :: VMTR_C2Wfact_pl  (ADM_gall_pl,ADM_kall,2,ADM_lall_pl)
    real(RP) :: VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   )
    real(RP) :: VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl)
    real(RP) :: VMTR_PHI         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_PHI_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: mxval, mnval

    integer  :: g, k, l, nq, K0
    !---------------------------------------------------------------------------

    call VMTR_getIJ_GSGAM2   ( VMTR_GSGAM2,    VMTR_GSGAM2_pl    )
    call VMTR_getIJ_W2Cfact  ( VMTR_W2Cfact,   VMTR_W2Cfact_pl   )
    call VMTR_getIJ_C2Wfact  ( VMTR_C2Wfact,   VMTR_C2Wfact_pl   )
    call VMTR_getIJ_C2WfactGz( VMTR_C2WfactGz, VMTR_C2WfactGz_pl )
    call VMTR_getIJ_PHI      ( VMTR_PHI,       VMTR_PHI_pl       )

    K0 = ADM_KNONE

    !--- get variables
    call prgvar_get_withdiag( rhog,   rhog_pl,   & ! [OUT]
                              rhogvx, rhogvx_pl, & ! [OUT]
                              rhogvy, rhogvy_pl, & ! [OUT]
                              rhogvz, rhogvz_pl, & ! [OUT]
                              rhogw,  rhogw_pl,  & ! [OUT]
                              rhoge,  rhoge_pl,  & ! [OUT]
                              rhogq,  rhogq_pl,  & ! [OUT]
                              rho,    rho_pl,    & ! [OUT]
                              pre,    pre_pl,    & ! [OUT]
                              tem,    tem_pl,    & ! [OUT]
                              vx,     vx_pl,     & ! [OUT]
                              vy,     vy_pl,     & ! [OUT]
                              vz,     vz_pl,     & ! [OUT]
                              w,      w_pl,      & ! [OUT]
                              q,      q_pl       ) ! [OUT]

    ein(:,:,:) = rhoge(:,:,:) / rhog(:,:,:)

    ! boundary condition
    call BNDCND_all( ADM_gall,                & ! [IN]
                     ADM_kall,                & ! [IN]
                     ADM_lall,                & ! [IN]
                     rho           (:,:,:),   & ! [INOUT]
                     vx            (:,:,:),   & ! [INOUT]
                     vy            (:,:,:),   & ! [INOUT]
                     vz            (:,:,:),   & ! [INOUT]
                     w             (:,:,:),   & ! [INOUT]
                     ein           (:,:,:),   & ! [INOUT]
                     tem           (:,:,:),   & ! [INOUT]
                     pre           (:,:,:),   & ! [INOUT]
                     rhog          (:,:,:),   & ! [INOUT]
                     rhogvx        (:,:,:),   & ! [INOUT]
                     rhogvy        (:,:,:),   & ! [INOUT]
                     rhogvz        (:,:,:),   & ! [INOUT]
                     rhogw         (:,:,:),   & ! [INOUT]
                     rhoge         (:,:,:),   & ! [INOUT]
                     VMTR_GSGAM2   (:,:,:),   & ! [IN]
                     VMTR_PHI      (:,:,:),   & ! [IN]
                     VMTR_C2Wfact  (:,:,:,:), & ! [IN]
                     VMTR_C2WfactGz(:,:,:,:)  ) ! [IN]

    ! value check
    mxval = maxval( pre(:,ADM_kmin:ADM_kmax,:) )
    mnval = minval( pre(:,ADM_kmin:ADM_kmax,:) )

    if (      mxval >= 2000.E+2_RP .OR. mxval <= 0.0_RP &
         .OR. mnval >= 2000.E+2_RP .OR. mnval <= 0.0_RP ) then ! > 2000hPa or negative?

       if( IO_L ) write(IO_FID_LOG,*) 'xxx Numerical instability occurs! STOP.', mxval, mnval
       write(*,*)          'xxx Numerical instability occurs! STOP.', mxval, mnval
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          if ( pre(g,k,l) > 2000.E+2_RP .OR. pre(g,k,l) < 0.0_RP ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx Invalid pressure value=', pre(g,k,l), ' at ', g, k, l
          endif
       enddo
       enddo
       enddo
       call PRC_MPIstop
    endif

    ! zonal and meridonal wind
    call cnvvar_vh2uv( u (:,:,:), u_pl (:,:,:), & ! [OUT]
                       v (:,:,:), v_pl (:,:,:), & ! [OUT]
                       vx(:,:,:), vx_pl(:,:,:), & ! [IN]
                       vy(:,:,:), vy_pl(:,:,:), & ! [IN]
                       vz(:,:,:), vz_pl(:,:,:), & ! [IN]
                       withcos = .false.        ) ! [IN]

    ! vertical wind at cell center
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax
          wc(:,k,l) = ( VMTR_W2Cfact(:,k,1,l) * w(:,k+1,l) &
                      + VMTR_W2Cfact(:,k,2,l) * w(:,k  ,l) )
       enddo
       wc(:,ADM_kmin-1,l) = w(:,ADM_kmin  ,l)
       wc(:,ADM_kmax+1,l) = w(:,ADM_kmax+1,l)
    enddo

    do l = 1, ADM_lall
       call history_in( 'ml_rho',  rho(:,:,l) )
       call history_in( 'ml_tem',  tem(:,:,l) )
       call history_in( 'ml_pres', pre(:,:,l) )

       call history_in( 'ml_u',    u  (:,:,l) )
       call history_in( 'ml_v',    v  (:,:,l) )
       call history_in( 'ml_w',    wc (:,:,l) )

       call history_in( 'ml_hgt',  GRD_vz(:,:,l,GRD_Z) ) ! geopotential height : Hydrostatic assumption
    enddo

    if (out_duvw) then
       u_old (:,:,:) = ( u_old (:,:,:) - u (:,:,:) ) / TIME_DTL * 86400.0_RP ! [m/s/day]
       v_old (:,:,:) = ( v_old (:,:,:) - v (:,:,:) ) / TIME_DTL * 86400.0_RP ! [m/s/day]
       wc_old(:,:,:) = ( wc_old(:,:,:) - wc(:,:,:) ) / TIME_DTL * 86400.0_RP ! [m/s/day]

       do l = 1, ADM_lall
          call history_in( 'ml_du', u_old (:,:,l) )
          call history_in( 'ml_dv', v_old (:,:,l) )
          call history_in( 'ml_dw', wc_old(:,:,l) )
       enddo

       u_old (:,:,:) = u (:,:,:)
       v_old (:,:,:) = v (:,:,:)
       wc_old(:,:,:) = wc(:,:,:)
    endif

    if (out_dtem) then
       tem_old(:,:,:) = ( tem_old(:,:,:) - tem(:,:,:) ) / TIME_DTL * 86400.0_RP ! [K/day]

       do l = 1, ADM_lall
          call history_in( 'ml_dtem', tem_old(:,:,l) )
       enddo

       tem_old(:,:,:) = tem(:,:,:)
    endif

    if (out_dq) then
       qv_old(:,:,:) = ( qv_old(:,:,:) - q(:,:,:,I_QV) ) * 1.E3_RP / TIME_DTL * 86400.0_RP ! [g/kg/day]

       do l = 1, ADM_lall
          call history_in( 'ml_dq', qv_old(:,:,l) )
       enddo

       qv_old(:,:,:) = q(:,:,:,I_QV)
    endif

    ! zonal and meridonal wind with cos(phi)
    if (out_uv_cos) then
       call cnvvar_vh2uv( ucos(:,:,:), u_pl (:,:,:), & ! [OUT]
                          vcos(:,:,:), v_pl (:,:,:), & ! [OUT]
                          vx  (:,:,:), vx_pl(:,:,:), & ! [IN]
                          vy  (:,:,:), vy_pl(:,:,:), & ! [IN]
                          vz  (:,:,:), vz_pl(:,:,:), & ! [IN]
                          withcos = .true.           ) ! [IN]

       do l = 1, ADM_lall
          call history_in( 'ml_ucos', ucos(:,:,l) )
          call history_in( 'ml_vcos', vcos(:,:,l) )
       enddo
    endif

    ! omega
    if (out_omg) then
       do l = 1, ADM_lall
          omg(:,:,l) = -GRAV * rho(:,:,l) * wc(:,:,l)

          call history_in( 'ml_omg', omg(:,:,l) )
       enddo
    endif

    ! relative humidity (liq+ice, liq, ice)

    if ( out_rha ) then
       call SATURATION_psat_all( ADM_gall,    & ! [IN]
                                 ADM_kall,    & ! [IN]
                                 ADM_lall,    & ! [IN]
                                 tem (:,:,:), & ! [IN]
                                 psat(:,:,:)  ) ! [OUT]

       do l = 1, ADM_lall
          rh(:,:,l) = q(:,:,l,I_QV) * rho(:,:,l) * Rvap * tem(:,:,l) &
                    / psat(:,:,l) &
                    * 100.0_RP

          call history_in( 'ml_rha', rh(:,:,l) )
       enddo
    endif

    if ( out_rh ) then
       call SATURATION_psat_liq( ADM_gall,    & ! [IN]
                                 ADM_kall,    & ! [IN]
                                 ADM_lall,    & ! [IN]
                                 tem (:,:,:), & ! [IN]
                                 psat(:,:,:)  ) ! [OUT]

       do l = 1, ADM_lall
          rh(:,:,l) = q(:,:,l,I_QV) * rho(:,:,l) * Rvap * tem(:,:,l) &
                    / psat(:,:,l) &
                    * 100.0_RP

          call history_in( 'ml_rh', rh(:,:,l) )
       enddo
    endif

    if ( out_rhi ) then
       call SATURATION_psat_ice( ADM_gall,    & ! [IN]
                                 ADM_kall,    & ! [IN]
                                 ADM_lall,    & ! [IN]
                                 tem (:,:,:), & ! [IN]
                                 psat(:,:,:)  ) ! [OUT]

       do l = 1, ADM_lall
          rh(:,:,l) = q(:,:,l,I_QV) * rho(:,:,l) * Rvap * tem(:,:,l) &
                    / psat(:,:,l) &
                    * 100.0_RP

          call history_in( 'ml_rhi', rh(:,:,l) )
       enddo
    endif

    ! potential temperature
    if ( out_th ) then
       call THRMDYN_th( ADM_gall,   & ! [IN]
                        ADM_kall,   & ! [IN]
                        ADM_lall,   & ! [IN]
                        tem(:,:,:), & ! [IN]
                        pre(:,:,:), & ! [IN]
                        th (:,:,:)  ) ! [OUT]

       thv(:,:,:) = th(:,:,:) * ( 1.D0 + 0.61D0 * q(:,:,:,I_QV) )

       do l = 1, ADM_lall
          call history_in( 'ml_th',  th (:,:,l) )
          call history_in( 'ml_thv', thv(:,:,l) )
       enddo
    endif

    ! moist static energy
    if ( out_mse ) then
       do l = 1, ADM_lall
          do k = 1, ADM_kall
             mse(:,k,l) = CPdry * tem(:,k,l)      &
                        + GRAV  * ( GRD_vz(:,k,l,GRD_Z) - GRD_zs (:,K0,l,GRD_ZSFC) ) &
                        + LHV   * q(:,k,l,I_QV)
          enddo

          call history_in( 'ml_mse',  mse(:,:,l) )
       enddo
    endif

    if ( out_th_prime ) then
       one   (:,:,:) = 1.0_RP
       one_pl(:,:,:) = 1.0_RP

       call GTL_global_sum_eachlayer( one, one_pl, area_prof )

       call THRMDYN_th( ADM_gall,   & ! [IN]
                        ADM_kall,   & ! [IN]
                        ADM_lall,   & ! [IN]
                        tem(:,:,:), & ! [IN]
                        pre(:,:,:), & ! [IN]
                        th (:,:,:)  ) ! [OUT]

       if ( ADM_have_pl ) then
          call THRMDYN_th( ADM_gall_pl,   & ! [IN]
                           ADM_kall,      & ! [IN]
                           ADM_lall_pl,   & ! [IN]
                           tem_pl(:,:,:), & ! [IN]
                           pre_pl(:,:,:), & ! [IN]
                           th_pl (:,:,:)  ) ! [OUT]
       endif

       call GTL_global_sum_eachlayer( th, th_pl, th_prof )

       do l = 1, ADM_lall
          do k = 1, ADM_kall
             th_prime(:,k,l) = th(:,k,l) - th_prof(k) / area_prof(k)
          enddo

          call history_in( 'ml_th_prime', th_prime(:,:,l) )
       enddo
    endif

    if (out_850hPa) then   ! [add] 20130705 R.Yoshida
       do l = 1, ADM_lall
          call sv_plev_uvwt( ADM_gall,        & ! [IN]
                             pre    (:,:,l),  & ! [IN]
                             u      (:,:,l),  & ! [IN]
                             v      (:,:,l),  & ! [IN]
                             w      (:,:,l),  & ! [IN]
                             tem    (:,:,l),  & ! [IN]
                             850.E2_RP,       & ! [IN]
                             u_slice(:,K0,l), & ! [OUT]
                             v_slice(:,K0,l), & ! [OUT]
                             w_slice(:,K0,l), & ! [OUT]
                             t_slice(:,K0,l)  ) ! [OUT]

          call history_in( 'sl_u850', u_slice(:,:,l) )
          call history_in( 'sl_v850', v_slice(:,:,l) )
          call history_in( 'sl_w850', w_slice(:,:,l) )
          call history_in( 'sl_t850', t_slice(:,:,l) )
       enddo
    endif

    if (out_500hPa) then
       do l = 1, ADM_lall
          call sv_plev_uvwt( ADM_gall,        & ! [IN]
                             pre    (:,:,l),  & ! [IN]
                             u      (:,:,l),  & ! [IN]
                             v      (:,:,l),  & ! [IN]
                             w      (:,:,l),  & ! [IN]
                             tem    (:,:,l),  & ! [IN]
                             500.E2_RP,       & ! [IN]
                             u_slice(:,K0,l), & ! [OUT]
                             v_slice(:,K0,l), & ! [OUT]
                             w_slice(:,K0,l), & ! [OUT]
                             t_slice(:,K0,l)  ) ! [OUT]

          call history_in( 'sl_u500', u_slice(:,:,l) )
          call history_in( 'sl_v500', v_slice(:,:,l) )
          call history_in( 'sl_w500', w_slice(:,:,l) )
          call history_in( 'sl_t500', t_slice(:,:,l) )
       enddo
    endif

    if (out_250hPa) then
       do l = 1, ADM_lall
          call sv_plev_uvwt( ADM_gall,        & ! [IN]
                             pre    (:,:,l),  & ! [IN]
                             u      (:,:,l),  & ! [IN]
                             v      (:,:,l),  & ! [IN]
                             w      (:,:,l),  & ! [IN]
                             tem    (:,:,l),  & ! [IN]
                             250.E2_RP,       & ! [IN]
                             u_slice(:,K0,l), & ! [OUT]
                             v_slice(:,K0,l), & ! [OUT]
                             w_slice(:,K0,l), & ! [OUT]
                             t_slice(:,K0,l)  ) ! [OUT]

          call history_in( 'sl_u250', u_slice(:,:,l) )
          call history_in( 'sl_v250', v_slice(:,:,l) )
          call history_in( 'sl_w250', w_slice(:,:,l) )
          call history_in( 'sl_t250', t_slice(:,:,l) )
       enddo
    endif

    if (out_100hPa) then
       do l = 1, ADM_lall
          call sv_plev_uvwt( ADM_gall,        & ! [IN]
                             pre    (:,:,l),  & ! [IN]
                             u      (:,:,l),  & ! [IN]
                             v      (:,:,l),  & ! [IN]
                             w      (:,:,l),  & ! [IN]
                             tem    (:,:,l),  & ! [IN]
                             100.E2_RP,       & ! [IN]
                             u_slice(:,K0,l), & ! [OUT]
                             v_slice(:,K0,l), & ! [OUT]
                             w_slice(:,K0,l), & ! [OUT]
                             t_slice(:,K0,l)  ) ! [OUT]

          call history_in( 'sl_u100', u_slice(:,:,l) )
          call history_in( 'sl_v100', v_slice(:,:,l) )
          call history_in( 'sl_w100', w_slice(:,:,l) )
          call history_in( 'sl_t100', t_slice(:,:,l) )
       enddo
    endif

    do l = 1, ADM_lall
       call sv_pre_sfc( ADM_gall,                 & ! [IN]
                        rho    (:,:,l),           & ! [IN]
                        pre    (:,:,l),           & ! [IN]
                        GRD_vz (:,:,l,GRD_Z),     & ! [IN]
                        GRD_zs (:,K0,l,GRD_ZSFC), & ! [IN]
                        rho_sfc(:,K0,l),          & ! [OUT]
                        pre_sfc(:,K0,l)           ) ! [OUT]

       call history_in( 'sl_ps', pre_sfc(:,:,l) )
    enddo

    !### tracers ###

    ! tracers
    do nq = 1, TRC_vmax
    do l  = 1, ADM_lall
       call history_in( 'ml_'//TRC_name(nq), q(:,:,l,nq) )
    enddo
    enddo

    ! hydrometeors
    do l  = 1, ADM_lall
       q_clw(:,:,l) = 0.0_RP
       q_cli(:,:,l) = 0.0_RP
       do nq = NQW_STR, NQW_END
          if ( nq == I_QC ) then
             q_clw(:,:,l) = q_clw(:,:,l) + q(:,:,l,nq)
          elseif( nq == I_QR ) then
             q_clw(:,:,l) = q_clw(:,:,l) + q(:,:,l,nq)
          elseif( nq == I_QI ) then
             q_cli(:,:,l) = q_cli(:,:,l) + q(:,:,l,nq)
          elseif( nq == I_QS ) then
             q_cli(:,:,l) = q_cli(:,:,l) + q(:,:,l,nq)
          elseif( nq == I_QG ) then
             q_cli(:,:,l) = q_cli(:,:,l) + q(:,:,l,nq)
          endif
       enddo
    enddo

    do l = 1, ADM_lall
       qtot(:,:,l) = q_clw(:,:,l) + q_cli(:,:,l)

       call history_in( 'ml_qtot', qtot(:,:,l) )
    enddo

    if (out_pw) then
       do l = 1, ADM_lall
          do k = ADM_kmin, ADM_kmax
             tmp2d(:,K0,l) = tmp2d(:,K0,l) + rho(:,k,l) * q(:,k,l,I_QV) * VMTR_GSGAM2(:,k,l) * GRD_dgz(k)
          enddo
          call history_in( 'sl_pw', tmp2d(:,:,l) )
       enddo
    endif

    if (out_lwp) then
       do l = 1, ADM_lall
          tmp2d(:,K0,l) = 0.0_RP
          do k = ADM_kmin, ADM_kmax
             tmp2d(:,K0,l) = tmp2d(:,K0,l) + rho(:,k,l) * q_clw(:,k,l) * VMTR_GSGAM2(:,k,l) * GRD_dgz(k)
          enddo
          call history_in( 'sl_lwp', tmp2d(:,:,l) )
       enddo
    endif

    if (out_iwp) then
       do l  = 1, ADM_lall
          tmp2d(:,K0,l) = 0.0_RP
          do k = ADM_kmin, ADM_kmax
             tmp2d(:,K0,l) = tmp2d(:,K0,l) + rho(:,k,l) * q_cli(:,k,l) * VMTR_GSGAM2(:,k,l) * GRD_dgz(k)
          enddo
          call history_in( 'sl_iwp', tmp2d(:,:,l) )
       enddo
    endif

    if ( AF_TYPE == 'DCMIP' .AND. NCHEM_STR > 0 .AND. NCHEM_END > 0 ) then
       do l  = 1, ADM_lall

          rhodz(:,k0,l) = 0.0_RP
          do k = ADM_kmin, ADM_kmax
             rhodz(:,k0,l) = rhodz(:,k0,l) + ( rho(:,k,l) * VMTR_GSGAM2(:,k,l) * GRD_dgz(k) )
          enddo

          tmp2d(:,k0,l) = 0.0_RP
          do k = ADM_kmin, ADM_kmax
             tmp2d(:,k0,l) = tmp2d(:,k0,l) + rho(:,k,l) * q(:,k,l,NCHEM_STR) * VMTR_GSGAM2(:,k,l) * GRD_dgz(k)
          enddo
          tmp2d(:,k0,l) = tmp2d(:,k0,l) / rhodz(:,k0,l)
          call history_in( 'sl_cl', tmp2d(:,:,l) )

          tmp2d(:,k0,l) = 0.0_RP
          do k = ADM_kmin, ADM_kmax
             tmp2d(:,k0,l) = tmp2d(:,k0,l) + rho(:,k,l) * q(:,k,l,NCHEM_END) * VMTR_GSGAM2(:,k,l) * GRD_dgz(k)
          enddo
          tmp2d(:,k0,l) = tmp2d(:,k0,l) / rhodz(:,k0,l)
          call history_in( 'sl_cl2', tmp2d(:,:,l) )

          tmp2d(:,k0,l) = 0.0_RP
          do k = ADM_kmin, ADM_kmax
             tmp2d(:,k0,l) = tmp2d(:,k0,l) + rho(:,k,l) * ( q(:,k,l,NCHEM_STR) + 2.0_RP * q(:,k,l,NCHEM_END) ) &
                                                        * VMTR_GSGAM2(:,k,l) * GRD_dgz(k)
          enddo
          tmp2d(:,k0,l) = tmp2d(:,k0,l) / rhodz(:,k0,l)
          call history_in( 'sl_cly', tmp2d(:,:,l) )

       enddo
    endif

    return
  end subroutine history_vars

  !-----------------------------------------------------------------------------
  subroutine sv_pre_sfc( &
       ijdim,   &
       rho,     &
       pre,     &
       z,       &
       z_srf,   &
       rho_srf, &
       pre_srf  )
    use scale_const, only: &
       GRAV => CONST_GRAV
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: rho    (ijdim,kdim)
    real(RP), intent(in)  :: pre    (ijdim,kdim)
    real(RP), intent(in)  :: z      (ijdim,kdim)
    real(RP), intent(in)  :: z_srf  (ijdim)
    real(RP), intent(out) :: rho_srf(ijdim)
    real(RP), intent(out) :: pre_srf(ijdim)

    integer :: ij
    !---------------------------------------------------------------------------

    !--- surface density ( extrapolation )
    do ij = 1, ijdim
       rho_srf(ij) = rho(ij,kmin) &
                   - ( rho(ij,kmin+1)-rho(ij,kmin) ) * ( z(ij,kmin)-z_srf(ij) ) / ( z(ij,kmin+1)-z(ij,kmin) )
    enddo

    !--- surface pressure ( hydrostatic balance )
    do ij = 1, ijdim
       pre_srf(ij) = pre(ij,kmin) + rho(ij,kmin) * GRAV * ( z(ij,kmin)-z_srf(ij) )
    enddo

    return
  end subroutine sv_pre_sfc

  !-----------------------------------------------------------------------------
  ! [add] 20130705 R.Yoshida
  subroutine sv_plev_uvwt( &
       ijdim, &
       pre,   &
       u_z,   &
       v_z,   &
       w_z,   &
       t_z,   &
       plev,  &
       u_p,   &
       v_p,   &
       w_p,   &
       t_p    )
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       kdim => ADM_kall, &
       kmin => ADM_kmin
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: pre(ijdim,kdim)
    real(RP), intent(in)  :: u_z(ijdim,kdim)
    real(RP), intent(in)  :: v_z(ijdim,kdim)
    real(RP), intent(in)  :: w_z(ijdim,kdim)
    real(RP), intent(in)  :: t_z(ijdim,kdim)
    real(RP), intent(in)  :: plev
    real(RP), intent(out) :: u_p(ijdim)
    real(RP), intent(out) :: v_p(ijdim)
    real(RP), intent(out) :: w_p(ijdim)
    real(RP), intent(out) :: t_p(ijdim)

    integer  :: kl(ijdim)
    integer  :: ku(ijdim)
    real(RP) :: wght_l, wght_u

    integer :: ij, k
    !---------------------------------------------------------------------------

    ! search z-level
    do ij = 1, ijdim
       do k = kmin, kdim
          if( pre(ij,k) < plev ) exit
       enddo
       if ( k >= kdim ) then
          write(*         ,*) 'xxx internal error! [sv_plev_uvwt/mod_history_vars] STOP.'
          if( IO_L ) write(IO_FID_LOG,*) 'xxx internal error! [sv_plev_uvwt/mod_history_vars] STOP.',kdim,k,plev,ij,pre(ij,:)
          call PRC_MPIstop
       endif

       ku(ij) = k
       kl(ij) = k - 1
    enddo

    ! interpolate
    do ij = 1, ijdim
       wght_l = ( log(plev)           - log(pre(ij,ku(ij))) ) &
              / ( log(pre(ij,kl(ij))) - log(pre(ij,ku(ij))) )

       wght_u = ( log(pre(ij,kl(ij))) - log(plev)           ) &
              / ( log(pre(ij,kl(ij))) - log(pre(ij,ku(ij))) )

       u_p(ij) = wght_l * u_z(ij,kl(ij)) + wght_u * u_z(ij,ku(ij))
       v_p(ij) = wght_l * v_z(ij,kl(ij)) + wght_u * v_z(ij,ku(ij))
       w_p(ij) = wght_l * w_z(ij,kl(ij)) + wght_u * w_z(ij,ku(ij))
       t_p(ij) = wght_l * t_z(ij,kl(ij)) + wght_u * t_z(ij,ku(ij))
    enddo

    return
  end subroutine sv_plev_uvwt

end module mod_history_vars
