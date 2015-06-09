!-------------------------------------------------------------------------------
!
!+  Vertical Implicit module
!
!-------------------------------------------------------------------------------
module mod_vi
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module is for the vertical implicit scheme of non-hydorostatic
  !       model.
  !
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                06-08-11   Add averaged rhog for tracer advection.
  !                11-05-07   Y.Yamada: Implementation of ES tuning cord by NEC.
  !                             Modified line: (20110405 NEC)
  !                                                or (!ftr< vi_small_step.r???)
  !                11-11-28   Y.Yamada: Merge Terai-san timer code
  !                                                    into the original code.
  !                11-12-29   Y.Yamada: Delete ES tuning and merge fjtimer
  !                12-3-9    S.Iga: tuned (phase4-1)
  !
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID
  use mod_adm, only: &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_gall,    &
     ADM_gall_pl, &
     ADM_kall
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: vi_setup
  public :: vi_small_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: vi_main
  private :: vi_rhow_update_matrix
  private :: vi_rhow_solver

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
#ifdef _FIXEDINDEX_
  real(RP), public               :: Mc   (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), private              :: Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public               :: Ml   (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), private              :: Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public               :: Mu   (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), private              :: Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
#else
  real(RP), public,  allocatable :: Mc   (:,:,:)
  real(RP), private, allocatable :: Mc_pl(:,:,:)
  real(RP), public,  allocatable :: Ml   (:,:,:)
  real(RP), private, allocatable :: Ml_pl(:,:,:)
  real(RP), public,  allocatable :: Mu   (:,:,:)
  real(RP), private, allocatable :: Mu_pl(:,:,:)
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine vi_setup
    use mod_adm, only: &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    implicit none
    !---------------------------------------------------------------------------

#ifndef _FIXEDINDEX_
    allocate( Mc   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Mu   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mu_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Ml   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Ml_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
#endif

    return
  end subroutine vi_setup

  !-----------------------------------------------------------------------------
  subroutine vi_small_step( &
       rhog,         rhog_pl,         &
       rhogvx,       rhogvx_pl,       &
       rhogvy,       rhogvy_pl,       &
       rhogvz,       rhogvz_pl,       &
       rhogw,        rhogw_pl,        &
       rhoge,        rhoge_pl,        &
       vx,           vx_pl,           &
       vy,           vy_pl,           &
       vz,           vz_pl,           &
       eth,          eth_pl,          &
       rhog_prim,    rhog_prim_pl,    &
       preg_prim,    preg_prim_pl,    &
       grhog0,       grhog0_pl,       &
       grhogvx0,     grhogvx0_pl,     &
       grhogvy0,     grhogvy0_pl,     &
       grhogvz0,     grhogvz0_pl,     &
       grhogw0,      grhogw0_pl,      &
       grhoge0,      grhoge0_pl,      &
       grhogetot0,   grhogetot0_pl,   &
       rhog_split,   rhog_split_pl,   &
       rhogvx_split, rhogvx_split_pl, &
       rhogvy_split, rhogvy_split_pl, &
       rhogvz_split, rhogvz_split_pl, &
       rhogw_split,  rhogw_split_pl,  &
       rhoge_split,  rhoge_split_pl,  &
       PROG_mean,     PROG_mean_pl,   &
       num_of_itr,                    &
       dt                             )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_gall_1d, &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_kmax,    &
       ADM_kmin
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       Rdry  => CNST_RAIR,  &
       CVdry => CNST_CV
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_afac, &
       GRD_bfac
    use mod_oprt, only: &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl, &
       VMTR_W2Cfact,    &
       VMTR_W2Cfact_pl
    use mod_time, only: &
       TIME_SPLIT
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    use mod_bndcnd, only: &
       BNDCND_rhovxvyvz
    use mod_numfilter, only: &
       numfilter_divdamp,   &
       numfilter_divdamp_2d
    use mod_src, only: &
       src_advection_convergence, &
       src_flux_convergence,      &
       src_pres_gradient,         &
       src_buoyancy,              &
       I_SRC_default,             &
       I_SRC_horizontal
    implicit none

    real(RP), intent(inout) :: rhog           (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho     ( G^1/2 x gam2 )
    real(RP), intent(inout) :: rhog_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogvx         (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vx  ( G^1/2 x gam2 )
    real(RP), intent(inout) :: rhogvx_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogvy         (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vy  ( G^1/2 x gam2 )
    real(RP), intent(inout) :: rhogvy_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogvz         (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vz  ( G^1/2 x gam2 )
    real(RP), intent(inout) :: rhogvz_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogw          (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w   ( G^1/2 x gam2 )
    real(RP), intent(inout) :: rhogw_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhoge          (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Ein ( G^1/2 x gam2 )
    real(RP), intent(inout) :: rhoge_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(RP), intent(in)    :: vx             (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_x
    real(RP), intent(in)    :: vx_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vy             (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_y
    real(RP), intent(in)    :: vy_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vz             (ADM_gall   ,ADM_kall,ADM_lall   ) ! Vh_z
    real(RP), intent(in)    :: vz_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: eth            (ADM_gall   ,ADM_kall,ADM_lall   ) ! enthalpy
    real(RP), intent(in)    :: eth_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog_prim      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog_prim_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg_prim      (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg_prim_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: grhog0         (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term (large step)
    real(RP), intent(in)    :: grhog0_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: grhogvx0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: grhogvx0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: grhogvy0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: grhogvy0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: grhogvz0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: grhogvz0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: grhogw0        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: grhogw0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: grhoge0        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: grhoge0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: grhogetot0     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: grhogetot0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(inout) :: rhog_split     (ADM_gall   ,ADM_kall,ADM_lall   ) ! split value
    real(RP), intent(inout) :: rhog_split_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogvx_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: rhogvx_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogvy_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: rhogvy_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogvz_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: rhogvz_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhogw_split    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: rhogw_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: rhoge_split    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: rhoge_split_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(out)   :: PROG_mean      (ADM_gall   ,ADM_kall,ADM_lall   ,5) ! mean_flux for tracer advection
    real(RP), intent(out)   :: PROG_mean_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl,5)

    integer, intent(in)    :: num_of_itr
    real(RP), intent(in)    :: dt

    ! merged array for communication
    real(RP) :: diff_vh   (ADM_gall   ,ADM_kall,ADM_lall   ,3)
    real(RP) :: diff_vh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)
    real(RP) :: diff_we   (ADM_gall   ,ADM_kall,ADM_lall   ,3)
    real(RP) :: diff_we_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,3)

    ! tendency term (large step + small step)
    real(RP) :: grhog        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhog_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhogvx      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhogvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhogvy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhogvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhogvz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhogvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhogw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhogw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhoge       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhoge_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    ! tendency term 2
    real(RP) :: drhog        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhog_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhogvx
    real(RP) :: drhogvx_pl
    real(RP) :: drhogvy
    real(RP) :: drhogvy_pl
    real(RP) :: drhogvz
    real(RP) :: drhogvz_pl
    real(RP) :: drhogw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhogw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhoge       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    ! divergence damping
    real(RP) :: ddivdvx      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: ddivdvx_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvx_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvy_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvy_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdvz_2d   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdvz_2d_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    ! pressure gradient force
    real(RP) :: dpgrad       (ADM_gall   ,ADM_kall,ADM_lall   ,GRD_XDIR:GRD_ZDIR)
    real(RP) :: dpgrad_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl,GRD_XDIR:GRD_ZDIR)
    real(RP) :: dpgradw      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dpgradw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    ! buoyancy force
    real(RP) :: dbuoiw       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dbuoiw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    ! pressure work
    real(RP) :: drhoge_pw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhoge_pwh   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pwh_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: gz_tilde     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: gz_tilde_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhog_h       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhog_h_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: eth_h        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: eth_h_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: preg_prim_split   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: preg_prim_split_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rweight_itr

    integer :: g, k, l, ns

    integer :: i, j, suf
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____vi_path0')

    ! full level -> half level
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall
          rhog_h(g,k,l) = ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k,  l) &
                          + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )
          eth_h (g,k,l) = 0.5_RP * ( GRD_afac(k) * eth(g,k,  l) &
                                  + GRD_bfac(k) * eth(g,k-1,l) )
       enddo
       enddo
       do g = 1, ADM_gall
          rhog_h(g,ADM_kmin-1,l) = rhog_h(g,ADM_kmin,l)
          eth_h (g,ADM_kmin-1,l) = eth_h (g,ADM_kmin,l)
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k,  l) &
                                + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
             eth_h_pl (g,k,l) = 0.5_RP * ( GRD_afac(k) * eth_pl(g,k,  l) &
                                        + GRD_bfac(k) * eth_pl(g,k-1,l) )
          enddo
          enddo
          do g = 1, ADM_gall_pl
             rhog_h_pl(g,ADM_kmin-1,l) = rhog_h_pl(g,ADM_kmin,l)
             eth_h_pl (g,ADM_kmin-1,l) = eth_h_pl (g,ADM_kmin,l)
          enddo
       enddo
    endif

    !---< Calculation of source term for rhog >

    call src_flux_convergence( rhogvx, rhogvx_pl, & ! [IN]
                               rhogvy, rhogvy_pl, & ! [IN]
                               rhogvz, rhogvz_pl, & ! [IN]
                               rhogw,  rhogw_pl,  & ! [IN]
                               drhog,  drhog_pl,  & ! [OUT]
                               I_SRC_default      ) ! [IN]

    !---< Calculation of source term for Vh(vx,vy,vz) and W >

    ! divergence damping
    call numfilter_divdamp( rhogvx,  rhogvx_pl,  & ! [IN]
                            rhogvy,  rhogvy_pl,  & ! [IN]
                            rhogvz,  rhogvz_pl,  & ! [IN]
                            rhogw,   rhogw_pl,   & ! [IN]
                            ddivdvx, ddivdvx_pl, & ! [OUT]
                            ddivdvy, ddivdvy_pl, & ! [OUT]
                            ddivdvz, ddivdvz_pl, & ! [OUT]
                            ddivdw,  ddivdw_pl   ) ! [OUT]

    call numfilter_divdamp_2d( rhogvx,     rhogvx_pl,     & ! [IN]
                               rhogvy,     rhogvy_pl,     & ! [IN]
                               rhogvz,     rhogvz_pl,     & ! [IN]
                               ddivdvx_2d, ddivdvx_2d_pl, & ! [OUT]
                               ddivdvy_2d, ddivdvy_2d_pl, & ! [OUT]
                               ddivdvz_2d, ddivdvz_2d_pl  ) ! [OUT]

    ! pressure force
    call src_pres_gradient( preg_prim, preg_prim_pl, & ! [IN]
                            dpgrad,    dpgrad_pl,    & ! [OUT]
                            dpgradw,   dpgradw_pl,   & ! [OUT]
                            I_SRC_default            ) ! [IN]

    ! buoyancy force
    call src_buoyancy( rhog_prim, rhog_prim_pl, & ! [IN]
                       dbuoiw,    dbuoiw_pl     ) ! [OUT]

    !---< Calculation of source term for rhoge >

    ! advection convergence for eth
    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    eth,    eth_pl,    & ! [IN]
                                    drhoge, drhoge_pl, & ! [OUT]
                                    I_SRC_default      ) ! [IN]

    ! pressure work
    do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          gz_tilde  (g,k,l) = GRAV - ( dpgradw(g,k,l)-dbuoiw(g,k,l) ) / rhog_h(g,k,l)
          drhoge_pwh(g,k,l) = -gz_tilde(g,k,l) * rhogw(g,k,l)
       enddo
       enddo

       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          drhoge_pw(g,k,l) = ( vx(g,k,l) * dpgrad(g,k,l,GRD_XDIR)          &
                             + vy(g,k,l) * dpgrad(g,k,l,GRD_YDIR)          &
                             + vz(g,k,l) * dpgrad(g,k,l,GRD_ZDIR)          ) &
                           + ( VMTR_W2Cfact(g,k,1,l) * drhoge_pwh(g,k+1,l) &
                             + VMTR_W2Cfact(g,k,2,l) * drhoge_pwh(g,k,  l) )
       enddo
       enddo
       do g = 1, ADM_gall
          drhoge_pw(g,ADM_kmin-1,l) = 0.0_RP
          drhoge_pw(g,ADM_kmax+1,l) = 0.0_RP
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             gz_tilde_pl  (g,k,l) = GRAV - ( dpgradw_pl(g,k,l)-dbuoiw_pl(g,k,l) ) / rhog_h_pl(g,k,l)
             drhoge_pwh_pl(g,k,l) = -gz_tilde_pl(g,k,l) * rhogw_pl(g,k,l)
          enddo
          enddo

          do k  = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             drhoge_pw_pl(g,k,l) = ( vx_pl(g,k,l) * dpgrad_pl(g,k,l,GRD_XDIR)          &
                                   + vy_pl(g,k,l) * dpgrad_pl(g,k,l,GRD_YDIR)          &
                                   + vz_pl(g,k,l) * dpgrad_pl(g,k,l,GRD_ZDIR)          ) &
                                 + ( VMTR_W2Cfact_pl(g,k,1,l) * drhoge_pwh_pl(g,k+1,l) &
                                   + VMTR_W2Cfact_pl(g,k,2,l) * drhoge_pwh_pl(g,k,  l) )
          enddo
          enddo
          do g = 1, ADM_gall_pl
             drhoge_pw_pl(g,ADM_kmin-1,l) = 0.0_RP
             drhoge_pw_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    endif

    !---< sum of tendencies ( large step + pres-grad + div-damp + div-damp_2d + buoyancy ) >

!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       grhog  (g,k,l) = grhog0  (g,k,l) + drhog  (g,k,l)
       grhogvx(g,k,l) = grhogvx0(g,k,l) - dpgrad (g,k,l,GRD_XDIR) + ddivdvx(g,k,l) + ddivdvx_2d(g,k,l)
       grhogvy(g,k,l) = grhogvy0(g,k,l) - dpgrad (g,k,l,GRD_YDIR) + ddivdvy(g,k,l) + ddivdvy_2d(g,k,l)
       grhogvz(g,k,l) = grhogvz0(g,k,l) - dpgrad (g,k,l,GRD_ZDIR) + ddivdvz(g,k,l) + ddivdvz_2d(g,k,l)
       grhogw (g,k,l) = grhogw0 (g,k,l)                           + ddivdw (g,k,l) * NON_HYDRO_ALPHA &
                                        - dpgradw(g,k,l) + dbuoiw(g,k,l)
       grhoge (g,k,l) = grhoge0 (g,k,l) + drhoge (g,k,l) + drhoge_pw(g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l  = 1, ADM_lall_pl
       do k  = 1, ADM_kall
       do g = 1, ADM_gall_pl
          grhog_pl  (g,k,l) = grhog0_pl  (g,k,l) + drhog_pl  (g,k,l)
          grhogvx_pl(g,k,l) = grhogvx0_pl(g,k,l) - dpgrad_pl (g,k,l,GRD_XDIR) + ddivdvx_pl(g,k,l) + ddivdvx_2d_pl(g,k,l)
          grhogvy_pl(g,k,l) = grhogvy0_pl(g,k,l) - dpgrad_pl (g,k,l,GRD_YDIR) + ddivdvy_pl(g,k,l) + ddivdvy_2d_pl(g,k,l)
          grhogvz_pl(g,k,l) = grhogvz0_pl(g,k,l) - dpgrad_pl (g,k,l,GRD_ZDIR) + ddivdvz_pl(g,k,l) + ddivdvz_2d_pl(g,k,l)
          grhogw_pl (g,k,l) = grhogw0_pl (g,k,l)                              + ddivdw_pl (g,k,l) * NON_HYDRO_ALPHA &
                                                 - dpgradw_pl(g,k,l) + dbuoiw_pl(g,k,l)
          grhoge_pl (g,k,l) = grhoge0_pl (g,k,l) + drhoge_pl (g,k,l) + drhoge_pw_pl(g,k,l)
       enddo
       enddo
       enddo
    endif

    ! initialization of mean mass flux
    rweight_itr = 1.0_RP / real(num_of_itr,kind=RP)

    PROG_mean(:,:,:,1) = rhogvx(:,:,:)
    PROG_mean(:,:,:,2) = rhogvy(:,:,:)
    PROG_mean(:,:,:,3) = rhogvz(:,:,:)
    PROG_mean(:,:,:,4) = rhog  (:,:,:)
    PROG_mean(:,:,:,5) = rhogw (:,:,:)

    PROG_mean_pl(:,:,:,1) = rhogvx_pl(:,:,:)
    PROG_mean_pl(:,:,:,2) = rhogvy_pl(:,:,:)
    PROG_mean_pl(:,:,:,3) = rhogvz_pl(:,:,:)
    PROG_mean_pl(:,:,:,4) = rhog_pl  (:,:,:)
    PROG_mean_pl(:,:,:,5) = rhogw_pl (:,:,:)

    ! update working matrix for vertical implicit solver
    call vi_rhow_update_matrix( eth_h,    eth_h_pl,    & ! [IN]
                                gz_tilde, gz_tilde_pl, & ! [IN]
                                dt                     ) ! [IN]

    call DEBUG_rapend  ('____vi_path0')

    !---------------------------------------------------------------------------
    !
    !> Start small step iteration
    !
    !---------------------------------------------------------------------------
    do ns = 1, num_of_itr

       call DEBUG_rapstart('____vi_path1')

       !---< calculation of preg_prim(*) from rhog(*) & rhoge(*) >

       do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             preg_prim_split(g,k,l) = rhoge_split(g,k,l) * Rdry / CVdry
          enddo
          enddo
          do g = 1, ADM_gall
             preg_prim_split(g,ADM_kmin-1,l) = preg_prim_split(g,ADM_kmin,l)
             preg_prim_split(g,ADM_kmax+1,l) = preg_prim_split(g,ADM_kmax,l)
          enddo

          do g = 1, ADM_gall
             rhoge_split(g,ADM_kmin-1,l) = rhoge_split(g,ADM_kmin,l)
             rhoge_split(g,ADM_kmax+1,l) = rhoge_split(g,ADM_kmax,l)
          enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                preg_prim_split_pl(g,k,l) = rhoge_split_pl(g,k,l) * Rdry / CVdry
             enddo
             enddo
             do g = 1, ADM_gall_pl
                preg_prim_split_pl(g,ADM_kmin-1,l) = preg_prim_split_pl(g,ADM_kmin,l)
                preg_prim_split_pl(g,ADM_kmax+1,l) = preg_prim_split_pl(g,ADM_kmax,l)
             enddo

             do g = 1, ADM_gall_pl
                rhoge_split_pl(g,ADM_kmin-1,l) = rhoge_split_pl(g,ADM_kmin,l)
                rhoge_split_pl(g,ADM_kmax+1,l) = rhoge_split_pl(g,ADM_kmax,l)
             enddo
          enddo
       endif

       if ( TIME_SPLIT ) then

          !---< Calculation of source term for Vh(vx,vy,vz) and W (split) >

          ! divergence damping
          call numfilter_divdamp( rhogvx_split, rhogvx_split_pl, & ! [IN]
                                  rhogvy_split, rhogvy_split_pl, & ! [IN]
                                  rhogvz_split, rhogvz_split_pl, & ! [IN]
                                  rhogw_split,  rhogw_split_pl,  & ! [IN]
                                  ddivdvx,      ddivdvx_pl,      & ! [OUT]
                                  ddivdvy,      ddivdvy_pl,      & ! [OUT]
                                  ddivdvz,      ddivdvz_pl,      & ! [OUT]
                                  ddivdw,       ddivdw_pl        ) ! [OUT]

          ! 2d divergence damping
          call numfilter_divdamp_2d( rhogvx_split, rhogvx_split_pl, & ! [IN]
                                     rhogvy_split, rhogvy_split_pl, & ! [IN]
                                     rhogvz_split, rhogvz_split_pl, & ! [IN]
                                     ddivdvx_2d,   ddivdvx_2d_pl,   & ! [OUT]
                                     ddivdvy_2d,   ddivdvy_2d_pl,   & ! [OUT]
                                     ddivdvz_2d,   ddivdvz_2d_pl    ) ! [OUT]

          ! pressure force
          ! dpgradw=0.0_RP becaude of f_type='HORIZONTAL'.
          call src_pres_gradient( preg_prim_split, preg_prim_split_pl, & ! [IN]
                                  dpgrad,          dpgrad_pl,          & ! [OUT]
                                  dpgradw,         dpgradw_pl,         & ! [OUT]
                                  I_SRC_horizontal                     ) ! [IN]

          ! buoyancy force
          ! not calculated, because this term is implicit.

          !---< sum of tendencies ( large step + split{ pres-grad + div-damp + div-damp_2d } ) >
!OCL SERIAL
          do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             drhogvx       = grhogvx(g,k,l) - dpgrad(g,k,l,GRD_XDIR) + ddivdvx(g,k,l) + ddivdvx_2d(g,k,l)
             drhogvy       = grhogvy(g,k,l) - dpgrad(g,k,l,GRD_YDIR) + ddivdvy(g,k,l) + ddivdvy_2d(g,k,l)
             drhogvz       = grhogvz(g,k,l) - dpgrad(g,k,l,GRD_ZDIR) + ddivdvz(g,k,l) + ddivdvz_2d(g,k,l)
             drhogw(g,k,l) = grhogw (g,k,l)                          + ddivdw (g,k,l) * NON_HYDRO_ALPHA

             diff_vh(g,k,l,1) = rhogvx_split(g,k,l) + drhogvx * dt
             diff_vh(g,k,l,2) = rhogvy_split(g,k,l) + drhogvy * dt
             diff_vh(g,k,l,3) = rhogvz_split(g,k,l) + drhogvz * dt
          enddo
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                drhogvx_pl       = grhogvx_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_XDIR) + ddivdvx_pl(g,k,l) + ddivdvx_2d_pl(g,k,l)
                drhogvy_pl       = grhogvy_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_YDIR) + ddivdvy_pl(g,k,l) + ddivdvy_2d_pl(g,k,l)
                drhogvz_pl       = grhogvz_pl(g,k,l) - dpgrad_pl(g,k,l,GRD_ZDIR) + ddivdvz_pl(g,k,l) + ddivdvz_2d_pl(g,k,l)
                drhogw_pl(g,k,l) = grhogw_pl (g,k,l)                             + ddivdw_pl (g,k,l) * NON_HYDRO_ALPHA

                diff_vh_pl(g,k,l,1) = rhogvx_split_pl(g,k,l) + drhogvx_pl * dt
                diff_vh_pl(g,k,l,2) = rhogvy_split_pl(g,k,l) + drhogvy_pl * dt
                diff_vh_pl(g,k,l,3) = rhogvz_split_pl(g,k,l) + drhogvz_pl * dt
             enddo
             enddo
             enddo
          endif

       else ! NO-SPLITING

          !---< sum of tendencies ( large step ) >
!OCL SERIAL
          do l = 1, ADM_lall
!OCL PARALLEL
          do k = 1, ADM_kall
          do g = 1, ADM_gall
             drhogw(g,k,l) = grhogw(g,k,l)

             diff_vh(g,k,l,1) = rhogvx_split(g,k,l) + grhogvx(g,k,l) * dt
             diff_vh(g,k,l,2) = rhogvy_split(g,k,l) + grhogvy(g,k,l) * dt
             diff_vh(g,k,l,3) = rhogvz_split(g,k,l) + grhogvz(g,k,l) * dt
          enddo
          enddo
          enddo

          if ( ADM_have_pl ) then
             do l = 1, ADM_lall_pl
             do k = 1, ADM_kall
             do g = 1, ADM_gall_pl
                drhogw_pl(g,k,l) = grhogw_pl(g,k,l)

                diff_vh_pl(g,k,l,1) = rhogvx_split_pl(g,k,l) + grhogvx_pl(g,k,l) * dt
                diff_vh_pl(g,k,l,2) = rhogvy_split_pl(g,k,l) + grhogvy_pl(g,k,l) * dt
                diff_vh_pl(g,k,l,3) = rhogvz_split_pl(g,k,l) + grhogvz_pl(g,k,l) * dt
             enddo
             enddo
             enddo
          endif

       endif ! Split/Non-split

       ! treatment for boundary condition
       do l = 1, ADM_lall
          call BNDCND_rhovxvyvz( ADM_gall,         & ! [IN]
                                 rhog(:,:,l),      & ! [IN]
                                 diff_vh(:,:,l,1), & ! [INOUT]
                                 diff_vh(:,:,l,2), & ! [INOUT]
                                 diff_vh(:,:,l,3)  ) ! [INOUT]
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             call BNDCND_rhovxvyvz( ADM_gall_pl,         & ! [IN]
                                    rhog_pl(:,:,l),      & ! [IN]
                                    diff_vh_pl(:,:,l,1), & ! [INOUT]
                                    diff_vh_pl(:,:,l,2), & ! [INOUT]
                                    diff_vh_pl(:,:,l,3)  ) ! [INOUT]
          enddo
       endif

       call COMM_data_transfer( diff_vh, diff_vh_pl )

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
          diff_vh(suf(ADM_gmax+1,ADM_gmin-1),k,l,1) = diff_vh(suf(ADM_gmax+1,ADM_gmin),k,l,1)
          diff_vh(suf(ADM_gmin-1,ADM_gmax+1),k,l,1) = diff_vh(suf(ADM_gmin,ADM_gmax+1),k,l,1)
          diff_vh(suf(ADM_gmax+1,ADM_gmin-1),k,l,2) = diff_vh(suf(ADM_gmax+1,ADM_gmin),k,l,2)
          diff_vh(suf(ADM_gmin-1,ADM_gmax+1),k,l,2) = diff_vh(suf(ADM_gmin,ADM_gmax+1),k,l,2)
          diff_vh(suf(ADM_gmax+1,ADM_gmin-1),k,l,3) = diff_vh(suf(ADM_gmax+1,ADM_gmin),k,l,3)
          diff_vh(suf(ADM_gmin-1,ADM_gmax+1),k,l,3) = diff_vh(suf(ADM_gmin,ADM_gmax+1),k,l,3)
       enddo
       enddo

       call DEBUG_rapend  ('____vi_path1')
       call DEBUG_rapstart('____vi_path2')

       !---< vertical implicit scheme >

       call vi_main( diff_we(:,:,:,1), diff_we_pl(:,:,:,1), & ! [OUT]
                     diff_we(:,:,:,2), diff_we_pl(:,:,:,2), & ! [OUT]
                     diff_we(:,:,:,3), diff_we_pl(:,:,:,3), & ! [OUT]
                     diff_vh(:,:,:,1), diff_vh_pl(:,:,:,1), & ! [IN]
                     diff_vh(:,:,:,2), diff_vh_pl(:,:,:,2), & ! [IN]
                     diff_vh(:,:,:,3), diff_vh_pl(:,:,:,3), & ! [IN]
                     rhog_split,       rhog_split_pl,       & ! [IN]
                     rhogvx_split,     rhogvx_split_pl,     & ! [IN]
                     rhogvy_split,     rhogvy_split_pl,     & ! [IN]
                     rhogvz_split,     rhogvz_split_pl,     & ! [IN]
                     rhogw_split,      rhogw_split_pl,      & ! [IN]
                     rhoge_split,      rhoge_split_pl,      & ! [IN]
                     preg_prim_split,  preg_prim_split_pl,  & ! [IN]
                     rhog,             rhog_pl,             & ! [IN]
                     rhogvx,           rhogvx_pl,           & ! [IN]
                     rhogvy,           rhogvy_pl,           & ! [IN]
                     rhogvz,           rhogvz_pl,           & ! [IN]
                     rhogw,            rhogw_pl,            & ! [IN]
                     eth,              eth_pl,              & ! [IN]
                     grhog,            grhog_pl,            & ! [IN]
                     drhogw,           drhogw_pl,           & ! [IN]
                     grhoge,           grhoge_pl,           & ! [IN]
                     grhogetot0,       grhogetot0_pl,       & ! [IN]
                     dt                                     ) ! [IN]

       ! treatment for boundary condition
       call COMM_data_transfer( diff_we, diff_we_pl )

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
          diff_we(suf(ADM_gmax+1,ADM_gmin-1),k,l,1) = diff_we(suf(ADM_gmax+1,ADM_gmin),k,l,1)
          diff_we(suf(ADM_gmin-1,ADM_gmax+1),k,l,1) = diff_we(suf(ADM_gmin,ADM_gmax+1),k,l,1)
          diff_we(suf(ADM_gmax+1,ADM_gmin-1),k,l,2) = diff_we(suf(ADM_gmax+1,ADM_gmin),k,l,2)
          diff_we(suf(ADM_gmin-1,ADM_gmax+1),k,l,2) = diff_we(suf(ADM_gmin,ADM_gmax+1),k,l,2)
          diff_we(suf(ADM_gmax+1,ADM_gmin-1),k,l,3) = diff_we(suf(ADM_gmax+1,ADM_gmin),k,l,3)
          diff_we(suf(ADM_gmin-1,ADM_gmax+1),k,l,3) = diff_we(suf(ADM_gmin,ADM_gmax+1),k,l,3)
       enddo
       enddo

       ! update split value and mean mass flux
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          rhogvx_split(g,k,l) = diff_vh(g,k,l,1)
          rhogvy_split(g,k,l) = diff_vh(g,k,l,2)
          rhogvz_split(g,k,l) = diff_vh(g,k,l,3)
          rhog_split  (g,k,l) = diff_we(g,k,l,1)
          rhogw_split (g,k,l) = diff_we(g,k,l,2)
          rhoge_split (g,k,l) = diff_we(g,k,l,3)

          PROG_mean(g,k,l,1) = PROG_mean(g,k,l,1) + rhogvx_split(g,k,l) * rweight_itr
          PROG_mean(g,k,l,2) = PROG_mean(g,k,l,2) + rhogvy_split(g,k,l) * rweight_itr
          PROG_mean(g,k,l,3) = PROG_mean(g,k,l,3) + rhogvz_split(g,k,l) * rweight_itr
          PROG_mean(g,k,l,4) = PROG_mean(g,k,l,4) + rhog_split  (g,k,l) * rweight_itr
          PROG_mean(g,k,l,5) = PROG_mean(g,k,l,5) + rhogw_split (g,k,l) * rweight_itr
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             rhogvx_split_pl(g,k,l) = diff_vh_pl(g,k,l,1)
             rhogvy_split_pl(g,k,l) = diff_vh_pl(g,k,l,2)
             rhogvz_split_pl(g,k,l) = diff_vh_pl(g,k,l,3)
             rhog_split_pl  (g,k,l) = diff_we_pl(g,k,l,1)
             rhogw_split_pl (g,k,l) = diff_we_pl(g,k,l,2)
             rhoge_split_pl (g,k,l) = diff_we_pl(g,k,l,3)

             PROG_mean_pl(g,k,l,1) = PROG_mean_pl(g,k,l,1) + rhogvx_split_pl(g,k,l) * rweight_itr
             PROG_mean_pl(g,k,l,2) = PROG_mean_pl(g,k,l,2) + rhogvy_split_pl(g,k,l) * rweight_itr
             PROG_mean_pl(g,k,l,3) = PROG_mean_pl(g,k,l,3) + rhogvz_split_pl(g,k,l) * rweight_itr
             PROG_mean_pl(g,k,l,4) = PROG_mean_pl(g,k,l,4) + rhog_split_pl  (g,k,l) * rweight_itr
             PROG_mean_pl(g,k,l,5) = PROG_mean_pl(g,k,l,5) + rhogw_split_pl (g,k,l) * rweight_itr
          enddo
          enddo
          enddo
       endif

       call DEBUG_rapend  ('____vi_path2')

    enddo  ! small step end

    !---------------------------------------------------------------------------
    !
    !
    !
    !---------------------------------------------------------------------------

    ! update prognostic variables
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       rhog  (g,k,l) = rhog  (g,k,l) + rhog_split  (g,k,l)
       rhogvx(g,k,l) = rhogvx(g,k,l) + rhogvx_split(g,k,l)
       rhogvy(g,k,l) = rhogvy(g,k,l) + rhogvy_split(g,k,l)
       rhogvz(g,k,l) = rhogvz(g,k,l) + rhogvz_split(g,k,l)
       rhogw (g,k,l) = rhogw (g,k,l) + rhogw_split (g,k,l)
       rhoge (g,k,l) = rhoge (g,k,l) + rhoge_split (g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          rhog_pl  (g,k,l) = rhog_pl  (g,k,l) + rhog_split_pl  (g,k,l)
          rhogvx_pl(g,k,l) = rhogvx_pl(g,k,l) + rhogvx_split_pl(g,k,l)
          rhogvy_pl(g,k,l) = rhogvy_pl(g,k,l) + rhogvy_split_pl(g,k,l)
          rhogvz_pl(g,k,l) = rhogvz_pl(g,k,l) + rhogvz_split_pl(g,k,l)
          rhogw_pl (g,k,l) = rhogw_pl (g,k,l) + rhogw_split_pl (g,k,l)
          rhoge_pl (g,k,l) = rhoge_pl (g,k,l) + rhoge_split_pl (g,k,l)
       enddo
       enddo
       enddo
    endif

    call OPRT_horizontalize_vec( rhogvx, rhogvx_pl, & ! [INOUT]
                                 rhogvy, rhogvy_pl, & ! [INOUT]
                                 rhogvz, rhogvz_pl  ) ! [INOUT]

    ! communication of mean velocity
    call COMM_data_transfer( PROG_mean, PROG_mean_pl )

    return
  end subroutine vi_small_step

  !-----------------------------------------------------------------------------
  subroutine vi_main( &
       rhog_split1,      rhog_split1_pl,      &
       rhogw_split1,     rhogw_split1_pl,     &
       rhoge_split1,     rhoge_split1_pl,     &
       rhogvx_split1,    rhogvx_split1_pl,    &
       rhogvy_split1,    rhogvy_split1_pl,    &
       rhogvz_split1,    rhogvz_split1_pl,    &
       rhog_split0,      rhog_split0_pl,      &
       rhogvx_split0,    rhogvx_split0_pl,    &
       rhogvy_split0,    rhogvy_split0_pl,    &
       rhogvz_split0,    rhogvz_split0_pl,    &
       rhogw_split0,     rhogw_split0_pl,     &
       rhoge_split0,     rhoge_split0_pl,     &
       preg_prim_split0, preg_prim_split0_pl, &
       rhog0,            rhog0_pl,            &
       rhogvx0,          rhogvx0_pl,          &
       rhogvy0,          rhogvy0_pl,          &
       rhogvz0,          rhogvz0_pl,          &
       rhogw0,           rhogw0_pl,           &
       eth0,             eth0_pl,             &
       grhog,            grhog_pl,            &
       grhogw,           grhogw_pl,           &
       grhoge,           grhoge_pl,           &
       grhogetot,        grhogetot_pl,        &
       dt                                     )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_cnst, only: &
       Rdry  => CNST_RAIR, &
       CVdry => CNST_CV
    use mod_vmtr, only: &
       VMTR_PHI,          &
       VMTR_PHI_pl,       &
       VMTR_C2WfactGz,    &
       VMTR_C2WfactGz_pl
    use mod_time, only: &
       TIME_SPLIT
    use mod_bndcnd, only: &
       BNDCND_rhow
    use mod_cnvvar, only: &
       cnvvar_rhogkin
    use mod_src, only: &
       src_flux_convergence,      &
       src_advection_convergence, &
       I_SRC_horizontal,          &
       I_SRC_default
    implicit none

    real(RP), intent(out) :: rhog_split1        (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars (split, at n+1 step)
    real(RP), intent(out) :: rhog_split1_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhogw_split1       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhogw_split1_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: rhoge_split1       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: rhoge_split1_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx_split1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx_split1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy_split1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy_split1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz_split1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz_split1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: rhog_split0        (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars (split, at n step)
    real(RP), intent(in)  :: rhog_split0_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx_split0      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx_split0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy_split0      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy_split0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz_split0      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz_split0_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw_split0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogw_split0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhoge_split0       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhoge_split0_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: preg_prim_split0   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: preg_prim_split0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: rhog0              (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars ( previous )
    real(RP), intent(in)  :: rhog0_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx0            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx0_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy0            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy0_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz0            (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz0_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw0             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogw0_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: eth0               (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: eth0_pl            (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: grhog              (ADM_gall   ,ADM_kall,ADM_lall   ) ! large step tendency
    real(RP), intent(in)  :: grhog_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: grhogw             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: grhogw_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: grhoge             (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: grhoge_pl          (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: grhogetot          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: grhogetot_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: dt

    real(RP) :: drhog       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term at t=n+1
    real(RP) :: drhog_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhoge      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhoge_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: drhogetot   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: drhogetot_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: grhog1      (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term ( large step + t=n+1 )
    real(RP) :: grhog1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhoge1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhoge1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: gpre        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: gpre_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhog1       (ADM_gall   ,ADM_kall,ADM_lall   ) ! prognostic vars ( previous + t=n,t=n+1 )
    real(RP) :: rhog1_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz1     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz1_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw1      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw1_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: rhogkin0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! kinetic energy ( previous                )
    real(RP) :: rhogkin0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogkin10   (ADM_gall   ,ADM_kall,ADM_lall   ) ! kinetic energy ( previous + split(t=n)   )
    real(RP) :: rhogkin10_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogkin11   (ADM_gall   ,ADM_kall,ADM_lall   ) ! kinetic energy ( previous + split(t=n+1) )
    real(RP) :: rhogkin11_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ethtot0     (ADM_gall   ,ADM_kall,ADM_lall   ) ! total enthalpy ( h + v^{2}/2 + phi, previous )
    real(RP) :: ethtot0_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    !---< update grhog & grhoge >

    if ( TIME_SPLIT ) then
       ! horizontal flux convergence
       call src_flux_convergence( rhogvx_split1, rhogvx_split1_pl, & ! [IN]
                                  rhogvy_split1, rhogvy_split1_pl, & ! [IN]
                                  rhogvz_split1, rhogvz_split1_pl, & ! [IN]
                                  rhogw_split0,  rhogw_split0_pl,  & ! [IN]
                                  drhog,         drhog_pl,         & ! [OUT]
                                  I_SRC_horizontal                 ) ! [IN]

       ! horizontal advection convergence
       call src_advection_convergence( rhogvx_split1, rhogvx_split1_pl, & ! [IN]
                                       rhogvy_split1, rhogvy_split1_pl, & ! [IN]
                                       rhogvz_split1, rhogvz_split1_pl, & ! [IN]
                                       rhogw_split0,  rhogw_split0_pl,  & ! [IN]
                                       eth0,          eth0_pl,          & ! [IN]
                                       drhoge,        drhoge_pl,        & ! [OUT]
                                       I_SRC_horizontal                 ) ! [IN]
    else
       drhog    (:,:,:) = 0.0_RP
       drhog_pl (:,:,:) = 0.0_RP
       drhoge   (:,:,:) = 0.0_RP
       drhoge_pl(:,:,:) = 0.0_RP
    endif

    ! update grhog, grhoge and calc source term of pressure
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       grhog1 (g,k,l) = grhog  (g,k,l) + drhog (g,k,l)
       grhoge1(g,k,l) = grhoge (g,k,l) + drhoge(g,k,l)
       gpre   (g,k,l) = grhoge1(g,k,l) * Rdry / CVdry
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          grhog1_pl (g,k,l) = grhog_pl  (g,k,l) + drhog_pl (g,k,l)
          grhoge1_pl(g,k,l) = grhoge_pl (g,k,l) + drhoge_pl(g,k,l)
          gpre_pl   (g,k,l) = grhoge1_pl(g,k,l) * Rdry / CVdry
       enddo
       enddo
       enddo
    endif

    !---------------------------------------------------------------------------
    ! verical implict calculation core
    !---------------------------------------------------------------------------

    ! boundary condition for rhogw_split1
    do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          rhogw_split1(g,k,l) = 0.0_RP
       enddo
       enddo

       call BNDCND_rhow( ADM_gall,               & ! [IN]
                         rhogvx_split1 (:,:,l),  & ! [IN]
                         rhogvy_split1 (:,:,l),  & ! [IN]
                         rhogvz_split1 (:,:,l),  & ! [IN]
                         rhogw_split1  (:,:,l),  & ! [INOUT]
                         VMTR_C2WfactGz(:,:,:,l) ) ! [IN]
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             rhogw_split1_pl(g,k,l) = 0.0_RP
          enddo
          enddo

          call BNDCND_rhow( ADM_gall_pl,               & ! [IN]
                            rhogvx_split1_pl (:,:,l),  & ! [IN]
                            rhogvy_split1_pl (:,:,l),  & ! [IN]
                            rhogvz_split1_pl (:,:,l),  & ! [IN]
                            rhogw_split1_pl  (:,:,l),  & ! [INOUT]
                            VMTR_C2WfactGz_pl(:,:,:,l) ) ! [IN]
       enddo
    endif

    ! update rhogw_split1
    call vi_rhow_solver( rhogw_split1,     rhogw_split1_pl,     & ! [INOUT]
                         rhogw_split0,     rhogw_split0_pl,     & ! [IN]
                         preg_prim_split0, preg_prim_split0_pl, & ! [IN]
                         rhog_split0,      rhog_split0_pl,      & ! [IN]
                         grhog1,           grhog1_pl,           & ! [IN]
                         grhogw,           grhogw_pl,           & ! [IN]
                         gpre,             gpre_pl,             & ! [IN]
                         dt                                     ) ! [IN]

    ! update rhog_split1
    call src_flux_convergence( rhogvx_split1, rhogvx_split1_pl, & ! [IN]
                               rhogvy_split1, rhogvy_split1_pl, & ! [IN]
                               rhogvz_split1, rhogvz_split1_pl, & ! [IN]
                               rhogw_split1,  rhogw_split1_pl,  & ! [IN]
                               drhog,         drhog_pl,         & ! [OUT]
                               I_SRC_default                    ) ! [IN]

    rhog_split1(:,:,:) = rhog_split0(:,:,:) + ( grhog(:,:,:) + drhog(:,:,:) ) * dt

    if ( ADM_have_pl ) then
       rhog_split1_pl(:,:,:) = rhog_split0_pl(:,:,:) + ( grhog_pl(:,:,:) + drhog_pl(:,:,:) ) * dt
    endif

    !---------------------------------------------------------------------------
    ! energy correction by Etotal (Satoh,2002)
    !---------------------------------------------------------------------------

    ! calc rhogkin ( previous )
    call cnvvar_rhogkin( rhog0,    rhog0_pl,   & ! [IN]
                         rhogvx0,  rhogvx0_pl, & ! [IN]
                         rhogvy0,  rhogvy0_pl, & ! [IN]
                         rhogvz0,  rhogvz0_pl, & ! [IN]
                         rhogw0,   rhogw0_pl,  & ! [IN]
                         rhogkin0, rhogkin0_pl ) ! [OUT]

    ! prognostic variables ( previous + split (t=n) )
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       rhog1  (g,k,l) = rhog0  (g,k,l) + rhog_split0  (g,k,l)
       rhogvx1(g,k,l) = rhogvx0(g,k,l) + rhogvx_split0(g,k,l)
       rhogvy1(g,k,l) = rhogvy0(g,k,l) + rhogvy_split0(g,k,l)
       rhogvz1(g,k,l) = rhogvz0(g,k,l) + rhogvz_split0(g,k,l)
       rhogw1 (g,k,l) = rhogw0 (g,k,l) + rhogw_split0 (g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          rhog1_pl  (g,k,l) = rhog0_pl  (g,k,l) + rhog_split0_pl  (g,k,l)
          rhogvx1_pl(g,k,l) = rhogvx0_pl(g,k,l) + rhogvx_split0_pl(g,k,l)
          rhogvy1_pl(g,k,l) = rhogvy0_pl(g,k,l) + rhogvy_split0_pl(g,k,l)
          rhogvz1_pl(g,k,l) = rhogvz0_pl(g,k,l) + rhogvz_split0_pl(g,k,l)
          rhogw1_pl (g,k,l) = rhogw0_pl (g,k,l) + rhogw_split0_pl (g,k,l)
       enddo
       enddo
       enddo
    endif

    ! calc rhogkin ( previous + split(t=n) )
    call cnvvar_rhogkin( rhog1,     rhog1_pl,    & ! [IN]
                         rhogvx1,   rhogvx1_pl,  & ! [IN]
                         rhogvy1,   rhogvy1_pl,  & ! [IN]
                         rhogvz1,   rhogvz1_pl,  & ! [IN]
                         rhogw1,    rhogw1_pl,   & ! [IN]
                         rhogkin10, rhogkin10_pl ) ! [OUT]

    ! prognostic variables ( previous + split (t=n+1) )
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       rhog1  (g,k,l) = rhog0  (g,k,l) + rhog_split1  (g,k,l)
       rhogvx1(g,k,l) = rhogvx0(g,k,l) + rhogvx_split1(g,k,l)
       rhogvy1(g,k,l) = rhogvy0(g,k,l) + rhogvy_split1(g,k,l)
       rhogvz1(g,k,l) = rhogvz0(g,k,l) + rhogvz_split1(g,k,l)
       rhogw1 (g,k,l) = rhogw0 (g,k,l) + rhogw_split1 (g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          rhog1_pl  (g,k,l) = rhog0_pl  (g,k,l) + rhog_split1_pl  (g,k,l)
          rhogvx1_pl(g,k,l) = rhogvx0_pl(g,k,l) + rhogvx_split1_pl(g,k,l)
          rhogvy1_pl(g,k,l) = rhogvy0_pl(g,k,l) + rhogvy_split1_pl(g,k,l)
          rhogvz1_pl(g,k,l) = rhogvz0_pl(g,k,l) + rhogvz_split1_pl(g,k,l)
          rhogw1_pl (g,k,l) = rhogw0_pl (g,k,l) + rhogw_split1_pl (g,k,l)
       enddo
       enddo
       enddo
    endif

    ! calc rhogkin ( previous + split(t=n+1) )
    call cnvvar_rhogkin( rhog1,     rhog1_pl,    & ! [IN]
                         rhogvx1,   rhogvx1_pl,  & ! [IN]
                         rhogvy1,   rhogvy1_pl,  & ! [IN]
                         rhogvz1,   rhogvz1_pl,  & ! [IN]
                         rhogw1,    rhogw1_pl,   & ! [IN]
                         rhogkin11, rhogkin11_pl ) ! [OUT]

    ! calculate total enthalpy ( h + v^{2}/2 + phi, previous )
    ethtot0(:,:,:) = eth0(:,:,:)                    &
                   + rhogkin0(:,:,:) / rhog0(:,:,:) &
                   + VMTR_PHI(:,:,:)

    if ( ADM_have_pl ) then
       ethtot0_pl(:,:,:) = eth0_pl(:,:,:)                       &
                         + rhogkin0_pl(:,:,:) / rhog0_pl(:,:,:) &
                         + VMTR_PHI_pl(:,:,:)
    endif

    ! advection convergence for eth + kin + phi
    call src_advection_convergence( rhogvx1,   rhogvx1_pl,   & ! [IN]
                                    rhogvy1,   rhogvy1_pl,   & ! [IN]
                                    rhogvz1,   rhogvz1_pl,   & ! [IN]
                                    rhogw1,    rhogw1_pl,    & ! [IN]
                                    ethtot0,   ethtot0_pl,   & ! [IN]
                                    drhogetot, drhogetot_pl, & ! [OUT]
                                    I_SRC_default            ) ! [IN]

    rhoge_split1(:,:,:) = rhoge_split0 (:,:,:)                     & ! t=n
                        + ( grhogetot  (:,:,:)                     & ! tendency of total energy (num.diff+smg+nudge)
                          + drhogetot  (:,:,:) ) * dt              & ! tendency of total energy (adv.conv.)
                        + ( rhogkin10  (:,:,:)                     & ! kinetic   energy (t=n)
                          - rhogkin11  (:,:,:) )                   & ! kinetic   energy (t=n+1)
                        + ( rhog_split0(:,:,:)                     & ! potential energy (diff,t=n)
                          - rhog_split1(:,:,:) ) * VMTR_PHI(:,:,:)   ! potential energy (diff,t=n+1)

    if ( ADM_have_pl ) then
       rhoge_split1_pl(:,:,:) = rhoge_split0_pl (:,:,:)                        &
                              + ( grhogetot_pl  (:,:,:)                        &
                                + drhogetot_pl  (:,:,:) ) * dt                 &
                              + ( rhogkin10_pl  (:,:,:)                        &
                                - rhogkin11_pl  (:,:,:) )                      &
                              + ( rhog_split0_pl(:,:,:)                        &
                                - rhog_split1_pl(:,:,:) ) * VMTR_PHI_pl(:,:,:)
    endif

    return
  end subroutine vi_main

  !-----------------------------------------------------------------------------
  subroutine vi_rhow_update_matrix( &
       eth,     eth_pl,     &
       g_tilde, g_tilde_pl, &
       dt                   )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       Rdry  => CNST_RAIR,  &
       CVdry => CNST_CV
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_rdgz,  &
       GRD_cfac,  &
       GRD_dfac
    use mod_vmtr, only: &
       VMTR_RGSGAM2,    &
       VMTR_RGSGAM2_pl, &
       VMTR_GAM2H,      &
       VMTR_GAM2H_pl,   &
       VMTR_RGAMH,      &
       VMTR_RGAMH_pl,   &
       VMTR_RGSQRTH,    &
       VMTR_RGSQRTH_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    implicit none

    real(RP), intent(in) :: eth       (ADM_gall   ,ADM_kall,ADM_lall   ) ! enthalpy at the half lev
    real(RP), intent(in) :: eth_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: g_tilde   (ADM_gall   ,ADM_kall,ADM_lall   ) ! effective gravitation at the half lev
    real(RP), intent(in) :: g_tilde_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: dt

    real(RP) :: GCVovR   ! g * Cv / R
    real(RP) :: ACVovRt2 ! alfa * Cv / R / dt**2

    integer :: g, k, l
    !---------------------------------------------------------------------------
    ! Original concept
    !
    ! A_o(:,:,:) = VMTR_RGSGAM2(:,:,:)
    ! A_i(:,:,:) = VMTR_GAM2H(:,:,:) * eth(:,:,:) ! [debug] 20120727 H.Yashiro
    ! B  (:,:,:) = g_tilde(:,:,:)
    ! C_o(:,:,:) = VMTR_RGAM2H (:,:,:) * ( CNST_CV / CNST_RAIR * CNST_EGRAV )
    ! C_i(:,:,:) = 1.0_RP / VMTR_RGAM2H(:,:,:)
    ! D  (:,:,:) = CNST_CV / CNST_RAIR / ( dt*dt ) / VMTR_RGSQRTH(:,:,:)
    !
    ! do k = ADM_kmin+1, ADM_kmax
    !    Mc(:,k,:) = dble(NON_HYDRO_ALPHA) *D(:,k,:)              &
    !              + GRD_rdgzh(k)                                 &
    !              * ( GRD_rdgz (k)   * A_o(:,k  ,:) * A_i(:,k,:) &
    !                + GRD_rdgz (k-1) * A_o(:,k-1,:) * A_i(:,k,:) &
    !                - 0.5_RP * ( GRD_dfac(k) - GRD_cfac(k-1) )    &
    !                * ( B(:,k,:) + C_o(:,k,:) * C_i(:,k,:) )     &
    !                )
    !    Mu(:,k,:) = - GRD_rdgzh(k) * GRD_rdgz(k) * A_o(:,k,:) * A_i(:,k+1,:) &
    !                - GRD_rdgzh(k) * 0.5_RP * GRD_cfac(k)                     &
    !                * ( B(:,k+1,:) + C_o(:,k,:) * C_i(:,k+1,:) )
    !    Ml(:,k,:) = - GRD_rdgzh(k) * GRD_rdgz(k) * A_o(:,k,:) * A_i(:,k-1,:) &
    !                + GRD_rdgzh(k) * 0.5_RP * GRD_dfac(k-1)                   &
    !                * ( B(:,k-1,:) + C_o(:,k,:) * C_i(:,k-1,:) )
    ! enddo

    call DEBUG_rapstart('____vi_rhow_update_matrix')

    GCVovR   = GRAV * CVdry / Rdry
    ACVovRt2 = real(NON_HYDRO_ALPHA,kind=RP) * CVdry / Rdry / ( dt*dt )

    do l = 1, ADM_lall
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Mc(g,k,l) = ACVovRt2 / VMTR_RGSQRTH(g,k,l)                             &
                    + GRD_rdgzh(k) * ( ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k  )   &
                                       + VMTR_RGSGAM2(g,k-1,l) * GRD_rdgz(k-1) ) &
                                       * VMTR_GAM2H  (g,k  ,l) * eth(g,k,l)      &
                                     - 0.5_RP * ( GRD_dfac(k) - GRD_cfac(k-1) )   &
                                     * ( g_tilde(g,k,l) + GCVovR )               )
       enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Mu(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k)                     &
                                      * VMTR_GAM2H  (g,k+1,l) * eth(g,k+1,l)                    &
                                      + 0.5_RP * GRD_cfac(k)                                     &
                                      * ( g_tilde   (g,k+1,l)                                   &
                                        + VMTR_GAM2H(g,k+1,l)* VMTR_RGAMH(g,k,l)**2 * GCVovR  ) )
       enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Ml(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2(g,k  ,l) * GRD_rdgz(k)                     &
                                      * VMTR_GAM2H  (g,k-1,l) * eth(g,k-1,l)                    &
                                      - 0.5_RP * GRD_dfac(k-1)                                   &
                                      * ( g_tilde   (g,k-1,l)                                   &
                                        + VMTR_GAM2H(g,k-1,l) * VMTR_RGAMH(g,k,l)**2 * GCVovR ) )
       enddo
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Mc_pl(g,k,l) = ACVovRt2 / VMTR_RGSQRTH_pl(g,k,l)                             &
                          + GRD_rdgzh(k) * ( ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k  )   &
                                             + VMTR_RGSGAM2_pl(g,k-1,l) * GRD_rdgz(k-1) ) &
                                             * VMTR_GAM2H_pl  (g,k  ,l) * eth_pl(g,k,l)   &
                                           - 0.5_RP * ( GRD_dfac(k) - GRD_cfac(k-1) )      &
                                           * ( g_tilde_pl(g,k,l) + GCVovR )               )
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Mu_pl(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k)                        &
                                            * VMTR_GAM2H_pl  (g,k+1,l) * eth_pl(g,k+1,l)                    &
                                            + 0.5_RP * GRD_cfac(k)                                           &
                                            * ( g_tilde_pl   (g,k+1,l)                                      &
                                              + VMTR_GAM2H_pl(g,k+1,l)* VMTR_RGAMH_pl(g,k,l)**2 * GCVovR  ) )
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Ml_pl(g,k,l) = -GRD_rdgzh(k) * ( VMTR_RGSGAM2_pl(g,k  ,l) * GRD_rdgz(k)                        &
                                            * VMTR_GAM2H_pl  (g,k-1,l) * eth_pl(g,k-1,l)                    &
                                            - 0.5_RP * GRD_dfac(k-1)                                         &
                                            * ( g_tilde_pl   (g,k-1,l)                                      &
                                              + VMTR_GAM2H_pl(g,k-1,l) * VMTR_RGAMH_pl(g,k,l)**2 * GCVovR ) )
          enddo
          enddo
       enddo
    endif

    call DEBUG_rapend('____vi_rhow_update_matrix')

    return
  end subroutine vi_rhow_update_matrix

  !-----------------------------------------------------------------------------
  subroutine vi_rhow_solver( &
       rhogw,  rhogw_pl,  &
       rhogw0, rhogw0_pl, &
       preg0,  preg0_pl,  &
       rhog0,  rhog0_pl,  &
       Sr,     Sr_pl,     &
       Sw,     Sw_pl,     &
       Sp,     Sp_pl,     &
       dt                 )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_cnst, only: &
       GRAV  => CNST_EGRAV, &
       Rdry  => CNST_RAIR,  &
       CVdry => CNST_CV
    use mod_grd, only: &
       GRD_rdgzh, &
       GRD_afac,  &
       GRD_bfac
    use mod_vmtr, only: &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl,  &
       VMTR_RGSGAM2H,    &
       VMTR_RGSGAM2H_pl, &
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_GSGAM2H,     &
       VMTR_GSGAM2H_pl
    use mod_runconf, only: &
       NON_HYDRO_ALPHA
    implicit none

    real(RP), intent(inout) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 ), n+1
    real(RP), intent(inout) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho            ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sr       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rho  at the full level
    real(RP), intent(in)    :: Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sw       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rhow at the half level
    real(RP), intent(in)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sp       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for pres at the full level
    real(RP), intent(in)    :: Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt

    real(RP) :: Sall    (ADM_gall,   ADM_kall)
    real(RP) :: Sall_pl (ADM_gall_pl,ADM_kall)
    real(RP) :: beta    (ADM_gall   )
    real(RP) :: beta_pl (ADM_gall_pl)
    real(RP) :: gamma   (ADM_gall,   ADM_kall)
    real(RP) :: gamma_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: alfa
    real(RP) :: CVovRt2 ! Cv / R / dt**2

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____vi_rhow_solver')

    alfa = real(NON_HYDRO_ALPHA,kind=RP)
    CVovRt2 = CVdry / Rdry / (dt*dt)

    do l = 1, ADM_lall
       ! calc Sall
       do k  = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Sall(g,k) = (   ( rhogw0(g,k,  l)*alfa + dt * Sw(g,k,  l) ) * VMTR_RGAMH  (g,k,  l)**2            &
                      - ( ( preg0 (g,k,  l)      + dt * Sp(g,k,  l) ) * VMTR_RGSGAM2(g,k,  l)               &
                        - ( preg0 (g,k-1,l)      + dt * Sp(g,k-1,l) ) * VMTR_RGSGAM2(g,k-1,l)               &
                        ) * dt * GRD_rdgzh(k)                                                               &
                      - ( ( rhog0 (g,k,  l)      + dt * Sr(g,k,  l) ) * VMTR_RGAM(g,k,  l)**2 * GRD_afac(k) &
                        + ( rhog0 (g,k-1,l)      + dt * Sr(g,k-1,l) ) * VMTR_RGAM(g,k-1,l)**2 * GRD_bfac(k) &
                        ) * dt * 0.5_RP * GRAV                                                               &
                      ) * CVovRt2
       enddo
       enddo

       ! boundary conditions
       do g = 1, ADM_gall
          rhogw(g,ADM_kmin,  l) = rhogw(g,ADM_kmin,  l) * VMTR_RGSGAM2H(g,ADM_kmin,  l)
          rhogw(g,ADM_kmax+1,l) = rhogw(g,ADM_kmax+1,l) * VMTR_RGSGAM2H(g,ADM_kmax+1,l)
          Sall (g,ADM_kmin+1)   = Sall (g,ADM_kmin+1) - Ml(g,ADM_kmin+1,l) * rhogw(g,ADM_kmin,  l)
          Sall (g,ADM_kmax  )   = Sall (g,ADM_kmax  ) - Mu(g,ADM_kmax,  l) * rhogw(g,ADM_kmax+1,l)
       enddo

       !---< solve tri-daigonal matrix >

       ! condition at ADM_kmin+1
       k = ADM_kmin+1
       do g = 1, ADM_gall
          beta (g)     = Mc(g,k,l)
          rhogw(g,k,l) = Sall(g,k) / beta(g)
       enddo

       ! forward
       do k = ADM_kmin+2, ADM_kmax
       do g = 1, ADM_gall
          gamma(g,k)   = Mu(g,k-1,l) / beta(g)
          beta (g)     = Mc(g,k,l) - Ml(g,k,l) * gamma(g,k) ! update beta
          rhogw(g,k,l) = ( Sall(g,k) - Ml(g,k,l) * rhogw(g,k-1,l) ) / beta(g)
       enddo
       enddo

       ! backward
       do k = ADM_kmax-1, ADM_kmin+1, -1
       do g = 1, ADM_gall
          rhogw(g,k  ,l) = rhogw(g,k  ,l) - gamma(g,k+1) * rhogw(g,k+1,l)
          rhogw(g,k+1,l) = rhogw(g,k+1,l) * VMTR_GSGAM2H(g,k+1,l) ! return value ( G^1/2 x gam2 )
       enddo
       enddo

       ! boundary treatment
       do g = 1, ADM_gall
          rhogw(g,ADM_kmin  ,l) = rhogw(g,ADM_kmin  ,l) * VMTR_GSGAM2H(g,ADM_kmin  ,l)
          rhogw(g,ADM_kmin+1,l) = rhogw(g,ADM_kmin+1,l) * VMTR_GSGAM2H(g,ADM_kmin+1,l)
          rhogw(g,ADM_kmax+1,l) = rhogw(g,ADM_kmax+1,l) * VMTR_GSGAM2H(g,ADM_kmax+1,l)
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k  = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Sall_pl(g,k) = (   ( rhogw0_pl(g,k,  l)*alfa + dt * Sw_pl(g,k,  l) ) * VMTR_RGAMH_pl  (g,k,  l)**2            &
                            - ( ( preg0_pl (g,k,  l)      + dt * Sp_pl(g,k,  l) ) * VMTR_RGSGAM2_pl(g,k,  l)               &
                              - ( preg0_pl (g,k-1,l)      + dt * Sp_pl(g,k-1,l) ) * VMTR_RGSGAM2_pl(g,k-1,l)               &
                              ) * dt * GRD_rdgzh(k)                                                                        &
                            - ( ( rhog0_pl (g,k,  l)      + dt * Sr_pl(g,k,  l) ) * VMTR_RGAM_pl(g,k,  l)**2 * GRD_afac(k) &
                              + ( rhog0_pl (g,k-1,l)      + dt * Sr_pl(g,k-1,l) ) * VMTR_RGAM_pl(g,k-1,l)**2 * GRD_bfac(k) &
                              ) * dt * 0.5_RP * GRAV                                                                        &
                            ) * CVovRt2
          enddo
          enddo

          do g = 1, ADM_gall_pl
             rhogw_pl(g,ADM_kmin,  l) = rhogw_pl(g,ADM_kmin,  l) * VMTR_RGSGAM2H_pl(g,ADM_kmin,  l)
             rhogw_pl(g,ADM_kmax+1,l) = rhogw_pl(g,ADM_kmax+1,l) * VMTR_RGSGAM2H_pl(g,ADM_kmax+1,l)
             Sall_pl (g,ADM_kmin+1)   = Sall_pl (g,ADM_kmin+1) - Ml_pl(g,ADM_kmin+1,l) * rhogw_pl(g,ADM_kmin,  l)
             Sall_pl (g,ADM_kmax  )   = Sall_pl (g,ADM_kmax  ) - Mu_pl(g,ADM_kmax,  l) * rhogw_pl(g,ADM_kmax+1,l)
          enddo

          k = ADM_kmin+1
          do g = 1, ADM_gall_pl
             beta_pl (g)     = Mc_pl(g,k,l)
             rhogw_pl(g,k,l) = Sall_pl(g,k) / beta_pl(g)
          enddo

          do k = ADM_kmin+2, ADM_kmax
          do g = 1, ADM_gall_pl
             gamma_pl(g,k)   = Mu_pl(g,k-1,l) / beta_pl(g)
             beta_pl (g)     = Mc_pl(g,k,l) - Ml_pl(g,k,l) * gamma_pl(g,k) ! update beta
             rhogw_pl(g,k,l) = ( Sall_pl(g,k) - Ml_pl(g,k,l) * rhogw_pl(g,k-1,l) ) / beta_pl(g)
          enddo
          enddo

          do k = ADM_kmax-1, ADM_kmin+1, -1
          do g = 1, ADM_gall_pl
             rhogw_pl(g,k,l) = rhogw_pl(g,k,l) - gamma_pl(g,k+1) * rhogw_pl(g,k+1,l)
          enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhogw_pl(g,k,l) = rhogw_pl(g,k,l) * VMTR_GSGAM2H_pl(g,k,l)
          enddo
          enddo
       enddo
    endif

    call DEBUG_rapend('____vi_rhow_solver')

    return
  end subroutine vi_rhow_solver

end module mod_vi
!-------------------------------------------------------------------------------

