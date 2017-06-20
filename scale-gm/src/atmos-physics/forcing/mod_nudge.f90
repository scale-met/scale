!-------------------------------------------------------------------------------
!> Module nudging
!!
!! @par Description
!!          This module is for the artificial forcing by nudging
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_nudge
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
  public :: NDG_setup
  public :: NDG_update_reference
  public :: NDG_apply_uvtp
  public :: NDG_apply_q

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: calc_wgt_horizontal

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: NDG_fact   (:,:,:,:)
  real(RP), private, allocatable :: NDG_fact_pl(:,:,:,:)

  real(RP), private, allocatable :: NDG_ref   (:,:,:,:)
  real(RP), private, allocatable :: NDG_ref_pl(:,:,:,:)

  integer,  private :: NDG_VMAX = -1
  integer,  private :: I_vx  = -1
  integer,  private :: I_vy  = -1
  integer,  private :: I_vz  = -1
  integer,  private :: I_w   = -1
  integer,  private :: I_tem = -1
  integer,  private :: I_pre = -1
  integer,  private :: I_qv  = -1

  real(RP), private :: NDG_tau_vxvyvz = -999.0_RP
  real(RP), private :: NDG_tau_w      = -999.0_RP
  real(RP), private :: NDG_tau_tem    = -999.0_RP
  real(RP), private :: NDG_tau_pre    = -999.0_RP
  real(RP), private :: NDG_tau_qv     = -999.0_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine NDG_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI    => CONST_PI,    &
       UNDEF => CONST_UNDEF, &
       Rdry  => CONST_Rdry,  &
       CVdry => CONST_CVdry
    use mod_adm, only: &
       ADM_have_pl,   &
       ADM_KNONE,     &
       ADM_lall,      &
       ADM_lall_pl,   &
       ADM_gall,      &
       ADM_gall_pl,   &
       ADM_kall,      &
       ADM_kmin,      &
       ADM_vlayer
    use mod_grd, only: &
       GRD_gz
    use mod_vmtr, only : &
       VMTR_getIJ_GSGAM2
    use mod_runconf, only: &
       FLAG_NUDGING
    implicit none

    integer :: NDG_kmin0 = -1
    integer :: NDG_kmax0 = 1000
    integer :: NDG_kmin1 = -1
    integer :: NDG_kmax1 = 1000

    logical :: NDG_hwgt            = .false.     ! weighted nudging option, depending on the distance from the pole
    real(RP) :: NDG_hwgt_center_lat =  35.0_RP      ! lat. of the pole ( -90<=v<=90 )
    real(RP) :: NDG_hwgt_center_lon = 135.0_RP      ! lon. of the pole (-180<=v<=180)
    real(RP) :: NDG_hwgt_halo1_dist =   0.0_RP      ! distance from the pole to the halo1 in [m] (0<=v<=NDG_hwgt_halo2_dist)
    real(RP) :: NDG_hwgt_halo2_dist = 2.0015778E7_RP ! distance from the pole to the halo2 in [m] (wt_ngd_halo1<=v<=pi*r_e)
    real(RP) :: NDG_hwgt_halo1_coef =   0.0_RP      ! min. coefficient (0<=v<=wt_ngd_max)
    real(RP) :: NDG_hwgt_halo2_coef =   1.0_RP      ! max. coefficient (wt_ngd_min<=v<=1)

    namelist /NUDGEPARAM/ &
       NDG_tau_vxvyvz,      &
       NDG_tau_w,           &
       NDG_tau_tem,         &
       NDG_tau_pre,         &
       NDG_tau_qv,          &
       NDG_kmin0,           &
       NDG_kmax0,           &
       NDG_kmin1,           &
       NDG_kmax1,           &
       NDG_hwgt,            &
       NDG_hwgt_center_lat, &
       NDG_hwgt_center_lon, &
       NDG_hwgt_halo1_dist, &
       NDG_hwgt_halo2_dist, &
       NDG_hwgt_halo1_coef, &
       NDG_hwgt_halo2_coef

    real(RP) :: NDG_rtau_vxvyvz
    real(RP) :: NDG_rtau_w
    real(RP) :: NDG_rtau_tem
    real(RP) :: NDG_rtau_pre
    real(RP) :: NDG_rtau_qv

    real(RP) :: wgt_vertical     (ADM_kall)
    real(RP) :: wgt_horizontal   (ADM_gall   ,ADM_KNONE,ADM_lall   ,1) ! 2008/09/10 [Add] M.Hara
    real(RP) :: wgt_horizontal_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,1)

    real(RP) :: VMTR_GSGAM2      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_GSGAM2_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer  :: k0, k1
    integer  :: g, k, l

    integer  :: ierr
    !---------------------------------------------------------------------------

    if( .NOT. FLAG_NUDGING ) return

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[nudging]/Category[nhm forcing]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=NUDGEPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** NUDGEPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist NUDGEPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=NUDGEPARAM)

    NDG_VMAX = 0

    if ( NDG_tau_vxvyvz > 0.0_RP ) then
       NDG_VMAX = NDG_VMAX + 3
       I_vx     = 1
       I_vy     = 2
       I_vz     = 3
       NDG_rtau_vxvyvz = 1.0_RP / NDG_tau_vxvyvz
    else
       NDG_rtau_vxvyvz = 0.0_RP
    endif

    if ( NDG_tau_w > 0.0_RP ) then
       NDG_VMAX = NDG_VMAX + 1
       I_w      = NDG_VMAX
       NDG_rtau_w = 1.0_RP / NDG_tau_w
    else
       NDG_rtau_w = 0.0_RP
    endif

    if ( NDG_tau_tem > 0.0_RP ) then
       NDG_VMAX = NDG_VMAX + 1
       I_tem    = NDG_VMAX
       NDG_rtau_tem = 1.0_RP / NDG_tau_tem
    else
       NDG_rtau_tem = 0.0_RP
    endif

    if ( NDG_tau_pre > 0.0_RP ) then
       NDG_VMAX = NDG_VMAX + 1
       I_pre    = NDG_VMAX
       NDG_rtau_pre = 1.0_RP / NDG_tau_pre
    else
       NDG_rtau_pre = 0.0_RP
    endif

    if ( NDG_tau_qv > 0.0_RP ) then
       NDG_VMAX = NDG_VMAX + 1
       I_qv     = NDG_VMAX
       NDG_rtau_qv = 1.0_RP / NDG_tau_qv
    else
       NDG_rtau_qv = 0.0_RP
    endif

    allocate( NDG_fact   (ADM_gall   ,ADM_kall,ADM_lall   ,NDG_VMAX) )
    allocate( NDG_fact_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,NDG_VMAX) )
    NDG_fact   (:,:,:,:) = 0.0_RP
    NDG_fact_pl(:,:,:,:) = 0.0_RP

    allocate( NDG_ref   (ADM_gall   ,ADM_kall,ADM_lall   ,NDG_VMAX) )
    allocate( NDG_ref_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,NDG_VMAX) )
    NDG_ref   (:,:,:,:) = UNDEF
    NDG_ref_pl(:,:,:,:) = UNDEF

    !---< vertical weight >--- [Add] M.Hara 2008/09/10
    NDG_kmin0 = max( min( NDG_kmin0, ADM_vlayer ), 1 )
    NDG_kmin1 = max( min( NDG_kmin1, ADM_vlayer ), 1 )
    NDG_kmax0 = max( min( NDG_kmax0, ADM_vlayer ), 1 )
    NDG_kmax1 = max( min( NDG_kmax1, ADM_vlayer ), 1 )

    k0 = NDG_kmin0
    k1 = NDG_kmin1
    NDG_kmin0 = min(k0,k1)
    NDG_kmin1 = max(k0,k1)
    k0 = NDG_kmax0
    k1 = NDG_kmax1
    NDG_kmax0 = min(k0,k1)
    NDG_kmax1 = max(k0,k1)

    if ( NDG_kmin1 > NDG_kmax0 ) then
       write(*,*) 'xxx Invalid vertical layers! STOP', NDG_kmin1, NDG_kmax0
       call PRC_MPIstop
    endif

    NDG_kmin0 = NDG_kmin0 + ADM_kmin - 1
    NDG_kmin1 = NDG_kmin1 + ADM_kmin - 1
    NDG_kmax0 = NDG_kmax0 + ADM_kmin - 1
    NDG_kmax1 = NDG_kmax1 + ADM_kmin - 1

    do k = 1, ADM_kall
       if ( k < NDG_kmin0 ) then
          wgt_vertical(k) = 0.0_RP
       elseif( k >= NDG_kmin0 .AND. k < NDG_kmin1 ) then
          wgt_vertical(k) = 0.5_RP * ( 1.0_RP + cos( PI * (GRD_gz(k)        -GRD_gz(NDG_kmin1) ) &
                                                        / (GRD_gz(NDG_kmin0)-GRD_gz(NDG_kmin1) ) ) )
       elseif( k >= NDG_kmax0 .AND. k < NDG_kmax1 ) then
          wgt_vertical(k) = 0.5_RP * ( 1.0_RP + cos( PI * (GRD_gz(k)        -GRD_gz(NDG_kmax0) ) &
                                                        / (GRD_gz(NDG_kmax1)-GRD_gz(NDG_kmax0) ) ) )
       elseif( NDG_kmax1 <= k ) then
          wgt_vertical(k) = 0.0_RP
       else
          wgt_vertical(k) = 1.0_RP
       endif
    enddo

    !---< horizontal weight >--- [Add] M.Hara 2008/09/10
    if ( NDG_hwgt ) then
       call calc_wgt_horizontal( NDG_hwgt_center_lon,        & ! [IN]
                                 NDG_hwgt_center_lat,        & ! [IN]
                                 NDG_hwgt_halo1_dist,        & ! [IN]
                                 NDG_hwgt_halo2_dist,        & ! [IN]
                                 NDG_hwgt_halo1_coef,        & ! [IN]
                                 NDG_hwgt_halo2_coef,        & ! [IN]
                                 wgt_horizontal   (:,:,:,:), & ! [OUT]
                                 wgt_horizontal_pl(:,:,:,:)  ) ! [OUT]
    else
       wgt_horizontal   (:,:,:,:) = 1.0_RP
       wgt_horizontal_pl(:,:,:,:) = 1.0_RP
    endif

    call VMTR_getIJ_GSGAM2( VMTR_GSGAM2, VMTR_GSGAM2_pl )

    !---< calc factor >---
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       NDG_fact(g,k,l,I_vx ) = wgt_vertical(k) * wgt_horizontal(g,ADM_KNONE,l,1) * NDG_rtau_vxvyvz
       NDG_fact(g,k,l,I_vy ) = wgt_vertical(k) * wgt_horizontal(g,ADM_KNONE,l,1) * NDG_rtau_vxvyvz
       NDG_fact(g,k,l,I_vz ) = wgt_vertical(k) * wgt_horizontal(g,ADM_KNONE,l,1) * NDG_rtau_vxvyvz
       NDG_fact(g,k,l,I_w  ) = wgt_vertical(k) * wgt_horizontal(g,ADM_KNONE,l,1) * NDG_rtau_w
       NDG_fact(g,k,l,I_tem) = wgt_vertical(k) * wgt_horizontal(g,ADM_KNONE,l,1) * NDG_rtau_tem &
                             * CVdry
       NDG_fact(g,k,l,I_pre) = wgt_vertical(k) * wgt_horizontal(g,ADM_KNONE,l,1) * NDG_rtau_pre &
                             * CVdry / Rdry * VMTR_GSGAM2(g,k,l)
       NDG_fact(g,k,l,I_qv ) = wgt_vertical(k) * wgt_horizontal(g,ADM_KNONE,l,1) * NDG_rtau_qv
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          NDG_fact_pl(g,k,l,I_vx ) = wgt_vertical(k) * wgt_horizontal_pl(g,ADM_KNONE,l,1) * NDG_rtau_vxvyvz
          NDG_fact_pl(g,k,l,I_vy ) = wgt_vertical(k) * wgt_horizontal_pl(g,ADM_KNONE,l,1) * NDG_rtau_vxvyvz
          NDG_fact_pl(g,k,l,I_vz ) = wgt_vertical(k) * wgt_horizontal_pl(g,ADM_KNONE,l,1) * NDG_rtau_vxvyvz
          NDG_fact_pl(g,k,l,I_w  ) = wgt_vertical(k) * wgt_horizontal_pl(g,ADM_KNONE,l,1) * NDG_rtau_w
          NDG_fact_pl(g,k,l,I_tem) = wgt_vertical(k) * wgt_horizontal_pl(g,ADM_KNONE,l,1) * NDG_rtau_tem &
                                   * CVdry
          NDG_fact_pl(g,k,l,I_pre) = wgt_vertical(k) * wgt_horizontal_pl(g,ADM_KNONE,l,1) * NDG_rtau_pre &
                                   * CVdry / Rdry * VMTR_GSGAM2_pl(g,k,l)
          NDG_fact_pl(g,k,l,I_qv ) = wgt_vertical(k) * wgt_horizontal_pl(g,ADM_KNONE,l,1) * NDG_rtau_qv
       enddo
       enddo
       enddo
    endif

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Nudging factors'
    if( IO_L ) write(IO_FID_LOG,*) '*** Tau(vxvyvz) [sec]: ', NDG_tau_vxvyvz
    if( IO_L ) write(IO_FID_LOG,*) '*** Tau(w)      [sec]: ', NDG_tau_w
    if( IO_L ) write(IO_FID_LOG,*) '*** Tau(tem)    [sec]: ', NDG_tau_tem
    if( IO_L ) write(IO_FID_LOG,*) '*** Tau(pre)    [sec]: ', NDG_tau_pre

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** vertical weight'
    if( IO_L ) write(IO_FID_LOG,*) '   k       z   weight(0-1)'
    do k = 1, ADM_kall
       if( IO_L ) write(IO_FID_LOG,'(1x,I3,1x,F7.1,1x,E14.6)') k, GRD_gz(k), wgt_vertical(k)
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** horizontal weight: ', NDG_hwgt
    if ( NDG_hwgt ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** center point(lat) [deg]                      : ', NDG_hwgt_center_lat
       if( IO_L ) write(IO_FID_LOG,*) '*** center point(lon) [deg]                      : ', NDG_hwgt_center_lon
       if( IO_L ) write(IO_FID_LOG,*) '*** distance to inner boundary of sponge zone [m]: ', NDG_hwgt_halo1_dist
       if( IO_L ) write(IO_FID_LOG,*) '*** distance to outer boundary of sponge zone [m]: ', NDG_hwgt_halo2_dist
       if( IO_L ) write(IO_FID_LOG,*) '*** weight at inner boundary [1]                 : ', NDG_hwgt_halo1_coef
       if( IO_L ) write(IO_FID_LOG,*) '*** weight at outer boundary [1]                 : ', NDG_hwgt_halo2_coef
    endif

    return
  end subroutine NDG_setup

  !-----------------------------------------------------------------------------
  subroutine NDG_update_reference( &
       ctime )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall_in, &
       ADM_kall,    &
       ADM_gmin,    &
       ADM_gmax
    use mod_comm, only: &
       COMM_var
    use mod_extdata, only: &
       extdata_update
    implicit none

    real(DP), intent(in) :: ctime

    real(RP) :: temp(ADM_gall_in,ADM_kall)

    logical :: eflag
    integer :: i, j, n, k, l
    !---------------------------------------------------------------------------

    if ( NDG_tau_vxvyvz > 0.0_RP ) then
       do l = 1, ADM_lall
          call extdata_update(temp(:,:),'vx',l,ctime,eflag)
          if ( .NOT. eflag ) then
             if( IO_L ) write(IO_FID_LOG,*) 'Not found reference file! :  vx (at ', ctime
          endif

          do k = 1, ADM_kall
             n = 1
             do j = ADM_gmin, ADM_gmax+1
             do i = ADM_gmin, ADM_gmax+1
                NDG_ref(suf(i,j),k,l,I_vx) = temp(n,k)
                n = n + 1
             enddo
             enddo
          enddo

          call extdata_update(temp(:,:),'vy',l,ctime,eflag)
          if ( .NOT. eflag ) then
             if( IO_L ) write(IO_FID_LOG,*) 'Not found reference file! :  vy (at ', ctime
          endif

          do k = 1, ADM_kall
             n = 1
             do j = ADM_gmin, ADM_gmax+1
             do i = ADM_gmin, ADM_gmax+1
                NDG_ref(suf(i,j),k,l,I_vy) = temp(n,k)
                n = n + 1
             enddo
             enddo
          enddo

          call extdata_update(temp(:,:),'vz',l,ctime,eflag)
          if ( .NOT. eflag ) then
             if( IO_L ) write(IO_FID_LOG,*) 'Not found reference file! :  vz (at ', ctime
          endif

          do k = 1, ADM_kall
             n = 1
             do j = ADM_gmin, ADM_gmax+1
             do i = ADM_gmin, ADM_gmax+1
                NDG_ref(suf(i,j),k,l,I_vz) = temp(n,k)
                n = n + 1
             enddo
             enddo
          enddo
       enddo
    endif

    if ( NDG_tau_w > 0.0_RP ) then
       do l = 1, ADM_lall
          call extdata_update(temp(:,:),'w',l,ctime,eflag)
          if ( .NOT. eflag ) then
             if( IO_L ) write(IO_FID_LOG,*) 'Not found reference file! :   w (at ', ctime
          endif

          do k = 1, ADM_kall
             n = 1
             do j = ADM_gmin, ADM_gmax+1
             do i = ADM_gmin, ADM_gmax+1
                NDG_ref(suf(i,j),k,l,I_w) = temp(n,k)
                n = n + 1
             enddo
             enddo
          enddo
       enddo
    endif

    if ( NDG_tau_tem > 0.0_RP ) then
       do l = 1, ADM_lall
          call extdata_update(temp(:,:),'tem',l,ctime,eflag)
          if ( .NOT. eflag ) then
             if( IO_L ) write(IO_FID_LOG,*) 'Not found reference file! : tem (at ', ctime
          endif

          do k = 1, ADM_kall
             n = 1
             do j = ADM_gmin, ADM_gmax+1
             do i = ADM_gmin, ADM_gmax+1
                NDG_ref(suf(i,j),k,l,I_tem) = temp(n,k)
                n = n + 1
             enddo
             enddo
          enddo
       enddo
    endif

    if ( NDG_tau_pre > 0.0_RP ) then
       do l = 1, ADM_lall
          call extdata_update(temp(:,:),'pre',l,ctime,eflag)
          if ( .NOT. eflag ) then
             if( IO_L ) write(IO_FID_LOG,*) 'Not found reference file! : pre (at ', ctime
          endif

          do k = 1, ADM_kall
             n = 1
             do j = ADM_gmin, ADM_gmax+1
             do i = ADM_gmin, ADM_gmax+1
                NDG_ref(suf(i,j),k,l,I_pre) = temp(n,k)
                n = n + 1
             enddo
             enddo
          enddo
       enddo
    endif

    if ( NDG_tau_qv > 0.0_RP ) then
       do l = 1, ADM_lall
          call extdata_update(temp(:,:),'qv',l,ctime,eflag)
          if ( .NOT. eflag ) then
             if( IO_L ) write(IO_FID_LOG,*) 'Not found reference file! :  qv (at ', ctime
          endif

          do k = 1, ADM_kall
             n = 1
             do j = ADM_gmin, ADM_gmax+1
             do i = ADM_gmin, ADM_gmax+1
                NDG_ref(suf(i,j),k,l,I_qv) = temp(n,k)
                n = n + 1
             enddo
             enddo
          enddo
       enddo
    endif

    call COMM_var( NDG_ref(:,:,:,:), NDG_ref_pl(:,:,:,:), ADM_kall, NDG_VMAX )

    return
  end subroutine NDG_update_reference

  !-----------------------------------------------------------------------------
  subroutine NDG_apply_uvtp( &
       rhog,      rhog_pl,      &
       vx,        vx_pl,        &
       vy,        vy_pl,        &
       vz,        vz_pl,        &
       w,         w_pl,         &
       tem,       tem_pl,       &
       pre,       pre_pl,       &
       frhogvx,   frhogvx_pl,   &
       frhogvy,   frhogvy_pl,   &
       frhogvz,   frhogvz_pl,   &
       frhogw,    frhogw_pl,    &
       frhoge,    frhoge_pl,    &
       out_tendency             )
    use scale_const, only: &
       GRAV  => CONST_GRAV, &
       CVdry => CONST_CVdry
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_oprt, only: &
       OPRT_horizontalize_vec
    use mod_vmtr, only: &
       VMTR_getIJ_GSGAM2H, &
       VMTR_getIJ_C2Wfact
    use mod_cnvvar, only: &
       cnvvar_vh2uv
    use mod_history, only: &
       history_in
    implicit none

    real(RP), intent(in)    :: rhog        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: rhog_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vx          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vx_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vy          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vy_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: vz          (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: vz_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: w           (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: w_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: tem         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: tem_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: pre         (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)    :: pre_pl      (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogvx     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogvx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogvy     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogvy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogvz     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhogw      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhogw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: frhoge      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: frhoge_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    logical,  intent(in)    :: out_tendency

    real(RP) :: dvx    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dvx_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dvy    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dvy_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dvz    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dvz_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dw     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dw_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dein   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dein_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: du     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: du_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dv     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dv_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dtem   (ADM_gall   ,ADM_kall,ADM_lall   )

    real(RP) :: VMTR_GSGAM2H   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: VMTR_GSGAM2H_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: VMTR_C2Wfact   (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP) :: VMTR_C2Wfact_pl(ADM_gall_pl,ADM_kall,2,ADM_lall_pl)

    real(RP) :: rhog_h
    real(RP) :: NDG_ref_w

    integer  :: g, k, l
    !---------------------------------------------------------------------------

    call VMTR_getIJ_GSGAM2H( VMTR_GSGAM2H, VMTR_GSGAM2H_pl )
    call VMTR_getIJ_C2Wfact( VMTR_C2Wfact, VMTR_C2Wfact_pl )

    dvx (:,:,:) = NDG_fact(:,:,:,I_vx ) * ( NDG_ref(:,:,:,I_vx ) - vx (:,:,:) )
    dvy (:,:,:) = NDG_fact(:,:,:,I_vy ) * ( NDG_ref(:,:,:,I_vy ) - vy (:,:,:) )
    dvz (:,:,:) = NDG_fact(:,:,:,I_vz ) * ( NDG_ref(:,:,:,I_vz ) - vz (:,:,:) )
    dein(:,:,:) = NDG_fact(:,:,:,I_tem) * ( NDG_ref(:,:,:,I_tem) - tem(:,:,:) ) &
                + NDG_fact(:,:,:,I_pre) * ( NDG_ref(:,:,:,I_pre) - pre(:,:,:) )

    if ( ADM_have_pl ) then
       dvx_pl (:,:,:) = NDG_fact_pl(:,:,:,I_vx ) * ( NDG_ref_pl(:,:,:,I_vx ) - vx_pl (:,:,:) )
       dvy_pl (:,:,:) = NDG_fact_pl(:,:,:,I_vy ) * ( NDG_ref_pl(:,:,:,I_vy ) - vy_pl (:,:,:) )
       dvz_pl (:,:,:) = NDG_fact_pl(:,:,:,I_vz ) * ( NDG_ref_pl(:,:,:,I_vz ) - vz_pl (:,:,:) )
       dein_pl(:,:,:) = NDG_fact_pl(:,:,:,I_tem) * ( NDG_ref_pl(:,:,:,I_tem) - tem_pl(:,:,:) ) &
                      + NDG_fact_pl(:,:,:,I_pre) * ( NDG_ref_pl(:,:,:,I_pre) - pre_pl(:,:,:) )
    endif

    call OPRT_horizontalize_vec( dvx, dvx_pl, & ! [INOUT]
                                 dvy, dvy_pl, & ! [INOUT]
                                 dvz, dvz_pl  ) ! [INOUT]

    frhogvx  (:,:,:) = frhogvx  (:,:,:) + dvx (:,:,:) * rhog(:,:,:)
    frhogvy  (:,:,:) = frhogvy  (:,:,:) + dvy (:,:,:) * rhog(:,:,:)
    frhogvz  (:,:,:) = frhogvz  (:,:,:) + dvz (:,:,:) * rhog(:,:,:)
    frhoge   (:,:,:) = frhoge   (:,:,:) + dein(:,:,:) * rhog(:,:,:)

    if ( ADM_have_pl ) then
       frhogvx_pl  (:,:,:) = frhogvx_pl  (:,:,:) + dvx_pl (:,:,:) * rhog_pl(:,:,:)
       frhogvy_pl  (:,:,:) = frhogvy_pl  (:,:,:) + dvy_pl (:,:,:) * rhog_pl(:,:,:)
       frhogvz_pl  (:,:,:) = frhogvz_pl  (:,:,:) + dvz_pl (:,:,:) * rhog_pl(:,:,:)
       frhoge_pl   (:,:,:) = frhoge_pl   (:,:,:) + dein_pl(:,:,:) * rhog_pl(:,:,:)
    endif

    if ( NDG_tau_w > 0.0_RP ) then
       do l = 1, ADM_lall
       do k = ADM_kmin,ADM_kmax+1
       do g = 1, ADM_gall
          rhog_h = ( VMTR_C2Wfact(g,k,1,l) * rhog(g,k  ,l) &
                   + VMTR_C2Wfact(g,k,2,l) * rhog(g,k-1,l) )

          NDG_ref_w = NDG_ref(g,k,l,I_w) * VMTR_GSGAM2H(g,k,l) / ( rhog_h * GRAV )

          dw(g,k,l) = NDG_fact(g,k,l,I_w) * ( NDG_ref_w - w(g,k,l) )

          frhogw(g,k,l) = frhogw(g,k,l) + dw(g,k,l) * rhog_h
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
          do k = ADM_kmin,ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhog_h = ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
                      + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )

             NDG_ref_w = NDG_ref_pl(g,k,l,I_w) * VMTR_GSGAM2H_pl(g,k,l) / ( rhog_h * GRAV )

             dw_pl(g,k,l) = NDG_fact_pl(g,k,l,I_w) * ( NDG_ref_w - w_pl(g,k,l) )

             frhogw_pl(g,k,l) = frhogw_pl(g,k,l) + dw_pl(g,k,l) * rhog_h
          enddo
          enddo
          enddo
       endif
    endif

    !--- output tendency
    if ( out_tendency ) then
       call cnvvar_vh2uv( du,  du_pl,  & ! [OUT]
                          dv,  dv_pl,  & ! [OUT]
                          dvx, dvx_pl, & ! [IN]
                          dvy, dvy_pl, & ! [IN]
                          dvz, dvz_pl  ) ! [IN]

       dtem(:,:,:) = dein(:,:,:) / CVdry

       do l = 1, ADM_lall
          call history_in('nudge_du',   du  (:,:,l))
          call history_in('nudge_dv',   dv  (:,:,l))
          call history_in('nudge_dtem', dtem(:,:,l))
          if ( NDG_tau_w > 0.0_RP ) then
             call history_in('nudge_dw', dw(:,:,l))
          endif
       enddo
    endif

    return
  end subroutine NDG_apply_uvtp

  !-----------------------------------------------------------------------------
  !> nudge q (only for specific humidity)
  !> note: treatment of energy/mass conservation?
  subroutine NDG_apply_q( &
       rhog,    &
       rhogq,   &
       dt       )
    use mod_adm, only: &
       ADM_lall,    &
       ADM_gall_in, &
       ADM_kall
    use mod_gtl, only: &
       GTL_clip_region
    use mod_runconf, only: &
       TRC_VMAX, &
       I_QV
    use mod_history, only: &
       history_in
    implicit none

    real(RP), intent(inout) :: rhog (ADM_gall_in,ADM_kall,ADM_lall)
    real(RP), intent(inout) :: rhogq(ADM_gall_in,ADM_kall,ADM_lall,TRC_VMAX)
    real(RP), intent(in)    :: dt

    real(RP) :: NDG_ref_qv_in(ADM_gall_in,ADM_kall,ADM_lall) ! trimmed
    real(RP) :: dqv          (ADM_gall_in,ADM_kall,ADM_lall) ! tendency of qv     [kg/kg/s]

    real(RP) :: drhogqv ! tendency of rhogqv [kg/m3/s]

    integer :: g, k, l
    !---------------------------------------------------------------------------

    if ( NDG_tau_qv > 0.0_RP ) then
       call GTL_clip_region(NDG_ref(:,:,:,I_qv),NDG_ref_qv_in(:,:,:),1,ADM_kall)

       do l = 1, ADM_lall
          do k = 1, ADM_kall
          do g = 1, ADM_gall_in
             drhogqv = NDG_fact(g,k,l,I_qv) * ( NDG_ref_qv_in(g,k,l)*rhog(g,k,l) - rhogq(g,k,l,I_QV) )
             ! limiter
             drhogqv = max( drhogqv, -rhogq(g,k,l,I_QV)/dt )
             dqv    (g,k,l) = drhogqv/rhog(g,k,l)

             rhogq(g,k,l,I_QV) = rhogq(g,k,l,I_QV) + drhogqv * dt
             rhog (g,k,l)      = rhog (g,k,l)      + drhogqv * dt
          enddo
          enddo

          call history_in('nudge_dqv',dqv(:,:,l))
       enddo
    endif

    return
  end subroutine NDG_apply_q

  !-----------------------------------------------------------------------------
  subroutine calc_wgt_horizontal( &
       center_lon, &
       center_lat, &
       halo1_dist, &
       halo2_dist, &
       halo1_coef, &
       halo2_coef, &
       weight,     &
       weight_pl   )
    use scale_const, only: &
       D2R    => CONST_D2R,   &
       UNDEF  => CONST_UNDEF, &
       RADIUS => CONST_RADIUS
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_KNONE
    use scale_vector, only: &
       VECTR_distance
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_LAT,  &
       GRD_LON,  &
       GRD_s,    &
       GRD_s_pl
    use mod_history, only: &
       history_in
    implicit none

    real(RP), intent(in)  :: center_lon ! nudging center longitude [degree]
    real(RP), intent(in)  :: center_lat ! nudging center latitude  [degree]
    real(RP), intent(in)  :: halo1_dist ! distance from the pole to the halo1 in [m] (0<=v<=halo2_dist)
    real(RP), intent(in)  :: halo2_dist ! distance from the pole to the halo2 in [m] (wt_ngd_halo1<=v<=pi*r_e)
    real(RP), intent(in)  :: halo1_coef ! coefficient (0<=v<=wt_ngd_max)
    real(RP), intent(in)  :: halo2_coef ! coefficient (wt_ngd_min<=v<=1)
    real(RP), intent(out) :: weight   (ADM_gall   ,ADM_KNONE,ADM_lall   ,1)
    real(RP), intent(out) :: weight_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,1)

    real(RP) :: center_lon_rad, center_lat_rad ! [rad]
    real(RP) :: dist, fact
    integer :: g, k0, l
    !---------------------------------------------------------------------------

    k0 = ADM_KNONE

    center_lon_rad = center_lon * D2R
    center_lat_rad = center_lat * D2R

    do l = 1, ADM_lall
    do g = 1, ADM_gall

       call VECTR_distance( RADIUS,                & ! [IN]
                            center_lon_rad,        & ! [IN]
                            center_lat_rad,        & ! [IN]
                            GRD_s(g,k0,l,GRD_LON), & ! [IN]
                            GRD_s(g,k0,l,GRD_LAT), & ! [IN]
                            dist                   ) ! [OUT]

       if ( dist < halo1_dist ) then
          fact = 0.0_RP
       elseif( dist >= halo1_dist .AND. dist <= halo2_dist ) then
          fact = (dist-halo1_dist) / (halo2_dist-halo1_dist)
       elseif( dist > halo2_dist ) then
          fact = 1.0_RP
       endif

       weight(g,k0,l,1) = ( 1.0_RP-fact ) * halo1_coef &
                        + (        fact ) * halo2_coef

    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl

          call VECTR_distance( RADIUS,                   & ! [IN]
                               center_lon,               & ! [IN]
                               center_lat,               & ! [IN]
                               GRD_s_pl(g,k0,l,GRD_LON), & ! [IN]
                               GRD_s_pl(g,k0,l,GRD_LAT), & ! [IN]
                               dist                      ) ! [OUT]

          if ( dist < halo1_dist ) then
             fact = 0.0_RP
          elseif( dist >= halo1_dist .AND. dist <= halo2_dist ) then
             fact = (dist-halo1_dist) / (halo2_dist-halo1_dist)
          elseif( dist > halo2_dist ) then
             fact = 1.0_RP
          endif

          weight_pl(g,k0,l,1) = ( 1.0_RP-fact ) * halo1_coef &
                              + (        fact ) * halo2_coef
       enddo
       enddo
    else
       weight_pl(:,:,:,:) = UNDEF
    endif

    call COMM_data_transfer(weight,weight_pl)

    do l = 1, ADM_lall
       call history_in('nudge_weight',weight(:,:,l,1))
    enddo

    return
  end subroutine calc_wgt_horizontal

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
    use mod_adm, only: &
       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_nudge
