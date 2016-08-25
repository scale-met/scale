!-------------------------------------------------------------------------------
!> Module vertical metrics
!!
!! @par Description
!!          In this module, the vertical metrics is calculated for the icoshaedral model
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_vmtr
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio

  use mod_adm, only: &
     ADM_lall,    &
     ADM_lall_pl, &
     ADM_iall,    &
     ADM_jall,    &
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
  public :: VMTR_setup

  public :: VMTR_getIJ_GAM2H
  public :: VMTR_getIJ_GSGAM2
  public :: VMTR_getIJ_GSGAM2H
  public :: VMTR_getIJ_RGSQRTH
  public :: VMTR_getIJ_RGAM
  public :: VMTR_getIJ_RGAMH
  public :: VMTR_getIJ_RGSGAM2
  public :: VMTR_getIJ_RGSGAM2H
  public :: VMTR_getIJ_W2Cfact
  public :: VMTR_getIJ_C2Wfact
  public :: VMTR_getIJ_C2WfactGz
  public :: VMTR_getIJ_VOLUME
  public :: VMTR_getIJ_PHI

  public :: VMTR_getin_GSGAM2
  public :: VMTR_getin_GSGAM2H
  public :: VMTR_getin_PHI

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  integer, public, parameter :: I_a = 1      ! index for W2Cfact
  integer, public, parameter :: I_b = 2

  integer, public, parameter :: I_c = 1      ! index for C2Wfact
  integer, public, parameter :: I_d = 2

  integer, public, parameter :: I_a_GZXH = 1 ! index for C2WfactGz
  integer, public, parameter :: I_b_GZXH = 2
  integer, public, parameter :: I_a_GZYH = 3
  integer, public, parameter :: I_b_GZYH = 4
  integer, public, parameter :: I_a_GZZH = 5
  integer, public, parameter :: I_b_GZZH = 6

#ifdef _FIXEDINDEX_
  real(RP), public              :: VMTR_GAM2H       (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_GAM2H_pl    (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_GSGAM2      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_GSGAM2_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_GSGAM2H     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_GSGAM2H_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)

  real(RP), public              :: VMTR_RGSQRTH     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSQRTH_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGAM        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGAM_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGAMH       (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGAMH_pl    (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGSGAM2     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSGAM2_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGSGAM2H    (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSGAM2H_pl (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)

  real(RP), public              :: VMTR_W2Cfact     (ADM_iall,ADM_jall,ADM_kall,2,ADM_lall   )
  real(RP), public              :: VMTR_W2Cfact_pl  (ADM_gall_pl      ,ADM_kall,2,ADM_lall_pl)
  real(RP), public              :: VMTR_C2Wfact     (ADM_iall,ADM_jall,ADM_kall,2,ADM_lall   )
  real(RP), public              :: VMTR_C2Wfact_pl  (ADM_gall_pl      ,ADM_kall,2,ADM_lall_pl)
  real(RP), public              :: VMTR_C2WfactGz   (ADM_iall,ADM_jall,ADM_kall,6,ADM_lall   )
  real(RP), public              :: VMTR_C2WfactGz_pl(ADM_gall_pl      ,ADM_kall,6,ADM_lall_pl)

  real(RP), public              :: VMTR_VOLUME      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_VOLUME_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)

  real(RP), public              :: VMTR_PHI         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_PHI_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
#else
  real(RP), public, allocatable :: VMTR_GAM2H       (:,:,:,:)   ! Gamma^2 at the half level
  real(RP), public, allocatable :: VMTR_GAM2H_pl    (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2      (:,:,:,:)   ! G^1/2 X Gamma^2 at the full level
  real(RP), public, allocatable :: VMTR_GSGAM2_pl   (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2H     (:,:,:,:)   ! G^1/2 X Gamma^2 at the half level
  real(RP), public, allocatable :: VMTR_GSGAM2H_pl  (:,:,:)

  real(RP), public, allocatable :: VMTR_RGSQRTH     (:,:,:,:)   ! 1 / G^1/2 at the half level
  real(RP), public, allocatable :: VMTR_RGSQRTH_pl  (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAM        (:,:,:,:)   ! 1 / Gamma at the integer level
  real(RP), public, allocatable :: VMTR_RGAM_pl     (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAMH       (:,:,:,:)   ! 1 / Gamma at the half level
  real(RP), public, allocatable :: VMTR_RGAMH_pl    (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2     (:,:,:,:)   ! 1 / (G^1/2 X Gamma^2) at the full level
  real(RP), public, allocatable :: VMTR_RGSGAM2_pl  (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2H    (:,:,:,:)   ! 1 / (G^1/2 X Gamma^2) at the half level
  real(RP), public, allocatable :: VMTR_RGSGAM2H_pl (:,:,:)

  real(RP), public, allocatable :: VMTR_W2Cfact     (:,:,:,:,:) ! factor for half to full level
  real(RP), public, allocatable :: VMTR_W2Cfact_pl  (:,:,:,:)
  real(RP), public, allocatable :: VMTR_C2Wfact     (:,:,:,:,:) ! factor for full to half level
  real(RP), public, allocatable :: VMTR_C2Wfact_pl  (:,:,:,:)
  real(RP), public, allocatable :: VMTR_C2WfactGz   (:,:,:,:,:) ! factor for full to half level with Gz
  real(RP), public, allocatable :: VMTR_C2WfactGz_pl(:,:,:,:)

  real(RP), public, allocatable :: VMTR_VOLUME      (:,:,:,:)   ! volume at the full level
  real(RP), public, allocatable :: VMTR_VOLUME_pl   (:,:,:)

  real(RP), public, allocatable :: VMTR_PHI         (:,:,:,:)   ! geopotential at the full level
  real(RP), public, allocatable :: VMTR_PHI_pl      (:,:,:)
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: deep = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine VMTR_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_adm, only: &
       ADM_KNONE,   &
       ADM_have_pl, &
       ADM_kmin,    &
       ADM_kmax
    use scale_const, only: &
       GRAV => CONST_GRAV
    use mod_comm, only: &
       COMM_data_transfer
    use mod_grd, only: &
       GRD_Z,        &
       GRD_ZH,       &
       GRD_dgz,      &
       GRD_dgzh,     &
       ORIG_GRD_vz => GRD_vz, &
       GRD_vz_pl,    &
       GRD_afac,     &
       GRD_bfac,     &
       GRD_cfac,     &
       GRD_dfac,     &
       GRD_rscale,   &
       GRD_grid_type
    use mod_gmtr, only: &
       ORIG_GMTR_area => GMTR_area, &
       GMTR_area_pl
    use mod_oprt, only: &
       OPRT_gradient,          &
       OPRT_horizontalize_vec, &
       OPRT_coef_grad,         &
       OPRT_coef_grad_pl
    implicit none

    integer, parameter :: var_max = 6

    integer, parameter :: JXH     = 1
    integer, parameter :: JYH     = 2
    integer, parameter :: JZH     = 3
    integer, parameter :: JX      = 4
    integer, parameter :: JY      = 5
    integer, parameter :: JZ      = 6

    real(RP) :: var   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ,var_max)
    real(RP) :: var_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl,var_max)

    !--- G^1/2
    real(RP) :: GSQRT    (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GSQRT_pl (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: GSQRTH   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GSQRTH_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    !--- Gamma factor
    real(RP) :: GAM      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GAM_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: GAMH     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GAMH_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    !--- vector G^z at the full level
    real(RP) :: GZX      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GZX_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: GZY      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GZY_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: GZZ      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GZZ_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    !--- vector G^z at the half level
    real(RP) :: GZXH     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GZXH_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: GZYH     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GZYH_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: GZZH     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: GZZH_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)

    real(RP) :: var_IJ   (ADM_gall,         ADM_kall,ADM_lall,var_max)
    real(RP) :: GRD_vz   (ADM_iall,ADM_jall,ADM_kall,ADM_lall,2)
    real(RP) :: GMTR_area(ADM_iall,ADM_jall,         ADM_lall)

    namelist / VMTRPARAM / &
       deep

    integer :: ierr
    integer :: i, j, g, k, l
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[vmtr]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=VMTRPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** VMTRPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist VMTRPARAM. STOP.'
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist VMTRPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=VMTRPARAM)

#ifndef _FIXEDINDEX_
    allocate( VMTR_GAM2H       (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_GAM2H_pl    (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GSGAM2      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GSGAM2H     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2H_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_RGSQRTH     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSQRTH_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAM        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAM_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAMH       (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAMH_pl    (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2H    (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2H_pl (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_W2Cfact     (ADM_iall,ADM_jall,ADM_kall,2,ADM_lall   ) )
    allocate( VMTR_W2Cfact_pl  (ADM_gall_pl      ,ADM_kall,2,ADM_lall_pl) )
    allocate( VMTR_C2Wfact     (ADM_iall,ADM_jall,ADM_kall,2,ADM_lall   ) )
    allocate( VMTR_C2Wfact_pl  (ADM_gall_pl      ,ADM_kall,2,ADM_lall_pl) )
    allocate( VMTR_C2WfactGz   (ADM_iall,ADM_jall,ADM_kall,6,ADM_lall   ) )
    allocate( VMTR_C2WfactGz_pl(ADM_gall_pl      ,ADM_kall,6,ADM_lall_pl) )

    allocate( VMTR_VOLUME      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_VOLUME_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )

    allocate( VMTR_PHI         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
    allocate( VMTR_PHI_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
#endif

    if( IO_L ) write(IO_FID_LOG,*) '*** setup metrics for 3-D control volume'

    GRD_vz    = reshape(ORIG_GRD_vz   ,shape(GRD_vz)   )
    GMTR_area = reshape(ORIG_GMTR_area,shape(GMTR_area))

    !--- if 1 layer model( shallow water model ),
    if ( ADM_kall == ADM_KNONE ) then

       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          VMTR_VOLUME(i,j,k,l) = GMTR_area(i,j,l)
       enddo
       enddo
       enddo
       enddo

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          VMTR_VOLUME_pl(g,k,l) = GMTR_area_pl(g,l)
       enddo
       enddo
       enddo

       return
    endif

    var   (:,:,:,:,:) = 0.0_RP
    var_pl(:,:,:,:)   = 0.0_RP

    !--- calculation of Jxh, Jyh, and Jzh

    call OPRT_gradient( var           (:,:,:,:,JXH:JZH), var_pl           (:,:,:,JXH:JZH), & ! [OUT]
                        GRD_vz        (:,:,:,:,GRD_ZH),  GRD_vz_pl        (:,:,:,GRD_ZH),  & ! [IN]
                        OPRT_coef_grad(:,:,:,:,:),       OPRT_coef_grad_pl(:,:,:)          ) ! [IN]

    call OPRT_horizontalize_vec( var(:,:,:,:,JXH), var_pl(:,:,:,JXH), & ! [INOUT]
                                 var(:,:,:,:,JYH), var_pl(:,:,:,JYH), & ! [INOUT]
                                 var(:,:,:,:,JZH), var_pl(:,:,:,JZH)  ) ! [INOUT]

    !--- calculation of Jx, Jy, and Jz
    call OPRT_gradient( var           (:,:,:,:,JX:JZ), var_pl           (:,:,:,JX:JZ), & ! [OUT]
                        GRD_vz        (:,:,:,:,GRD_Z), GRD_vz_pl        (:,:,:,GRD_Z), & ! [IN]
                        OPRT_coef_grad(:,:,:,:,:),     OPRT_coef_grad_pl(:,:,:)        ) ! [IN]

    call OPRT_horizontalize_vec( var(:,:,:,:,JX), var_pl(:,:,:,JX), & ! [INOUT]
                                 var(:,:,:,:,JY), var_pl(:,:,:,JY), & ! [INOUT]
                                 var(:,:,:,:,JZ), var_pl(:,:,:,JZ)  ) ! [INOUT]

    !--- fill HALO
    var_IJ = reshape(var,shape(var_IJ))
    call COMM_data_transfer( var_IJ, var_pl )
    var = reshape(var_IJ,shape(var))

    !--- G^1/2 = dz/dgz
    do l = 1, ADM_lall
       !--- calculation of G^1/2 at full level
       do k = ADM_kmin, ADM_kmax
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          GSQRT(i,j,k,l) = ( GRD_vz(i,j,k+1,l,GRD_ZH) - GRD_vz(i,j,k,l,GRD_ZH) ) / GRD_dgz(k)
       enddo
       enddo
       enddo
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          GSQRT(i,j,ADM_kmin-1,l) = GSQRT(i,j,ADM_kmin,l)
          GSQRT(i,j,ADM_kmax+1,l) = GSQRT(i,j,ADM_kmax,l)
       enddo
       enddo

       !--- calculation of G^1/2 at half level
       do k = ADM_kmin, ADM_kmax+1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          GSQRTH(i,j,k,l) = ( GRD_vz(i,j,k,l,GRD_Z) - GRD_vz(i,j,k-1,l,GRD_Z) ) / GRD_dgzh(k)
       enddo
       enddo
       enddo
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          GSQRTH(i,j,ADM_kmin-1,l) = GSQRTH(i,j,ADM_kmin,l)
       enddo
       enddo
    enddo

    !--- Gamma = (a+z) / a
    if ( deep ) then
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          GAM (i,j,k,l) = 1.0_RP + GRD_vz(i,j,k,l,GRD_Z)  / GRD_rscale
          GAMH(i,j,k,l) = 1.0_RP + GRD_vz(i,j,k,l,GRD_ZH) / GRD_rscale
       enddo
       enddo
       enddo
       enddo
    else
       do l = 1, ADM_lall
       do k = 1, ADM_kall
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          GAM (i,j,k,l) = 1.0_RP
          GAMH(i,j,k,l) = 1.0_RP
       enddo
       enddo
       enddo
       enddo
    endif

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do j = 1, ADM_jall
    do i = 1, ADM_iall
       VMTR_GAM2H   (i,j,k,l) = GAMH(i,j,k,l) * GAMH(i,j,k,l)
       VMTR_GSGAM2  (i,j,k,l) = GAM (i,j,k,l) * GAM (i,j,k,l) * GSQRT (i,j,k,l)
       VMTR_GSGAM2H (i,j,k,l) = GAMH(i,j,k,l) * GAMH(i,j,k,l) * GSQRTH(i,j,k,l)

       VMTR_RGSQRTH (i,j,k,l) = 1.0_RP / GSQRTH(i,j,k,l)
       VMTR_RGAM    (i,j,k,l) = 1.0_RP / GAM (i,j,k,l)
       VMTR_RGAMH   (i,j,k,l) = 1.0_RP / GAMH(i,j,k,l)
       VMTR_RGSGAM2 (i,j,k,l) = 1.0_RP / VMTR_GSGAM2 (i,j,k,l)
       VMTR_RGSGAM2H(i,j,k,l) = 1.0_RP / VMTR_GSGAM2H(i,j,k,l)
    enddo
    enddo
    enddo
    enddo

    ! full level <-> half level interpolation factor
    do l = 1, ADM_lall
       do k = ADM_kmin, ADM_kmax+1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          VMTR_C2Wfact(i,j,k,I_a,l) = GRD_afac(k) * VMTR_RGSGAM2(i,j,k  ,l) * VMTR_GSGAM2H(i,j,k,l)
          VMTR_C2Wfact(i,j,k,I_b,l) = GRD_bfac(k) * VMTR_RGSGAM2(i,j,k-1,l) * VMTR_GSGAM2H(i,j,k,l)
       enddo
       enddo
       enddo
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          VMTR_C2Wfact(i,j,ADM_kmin-1,I_a,l) = 0.0_RP
          VMTR_C2Wfact(i,j,ADM_kmin-1,I_b,l) = 0.0_RP
       enddo
       enddo

       do k = ADM_kmin-1, ADM_kmax
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          VMTR_W2Cfact(i,j,k,I_c,l) = GRD_cfac(k) * VMTR_GSGAM2(i,j,k,l) * VMTR_RGSGAM2H(i,j,k+1,l)
          VMTR_W2Cfact(i,j,k,I_d,l) = GRD_dfac(k) * VMTR_GSGAM2(i,j,k,l) * VMTR_RGSGAM2H(i,j,k  ,l)
       enddo
       enddo
       enddo
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          VMTR_W2Cfact(i,j,ADM_kmax+1,I_c,l) = 0.0_RP
          VMTR_W2Cfact(i,j,ADM_kmax+1,I_d,l) = 0.0_RP
       enddo
       enddo
    enddo

    ! full level <-> half level interpolation factor with Gz

    !--- Gz(X) = - JX / G^1/2
    !--- Gz(Y) = - JY / G^1/2
    !--- Gz(Z) = - JZ / G^1/2
    do l = 1, ADM_lall
       do k = 1, ADM_kall
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          GZXH(i,j,k,l) = -var(i,j,k,l,JXH) / GSQRTH(i,j,k,l)
          GZYH(i,j,k,l) = -var(i,j,k,l,JYH) / GSQRTH(i,j,k,l)
          GZZH(i,j,k,l) = -var(i,j,k,l,JZH) / GSQRTH(i,j,k,l)
          GZX (i,j,k,l) = -var(i,j,k,l,JX)  / GSQRT (i,j,k,l)
          GZY (i,j,k,l) = -var(i,j,k,l,JY)  / GSQRT (i,j,k,l)
          GZZ (i,j,k,l) = -var(i,j,k,l,JZ)  / GSQRT (i,j,k,l)
       enddo
       enddo
       enddo

       do k = ADM_kmin, ADM_kmax+1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          VMTR_C2WfactGz(i,j,k,I_a_GZXH,l) = GRD_afac(k) * VMTR_RGSGAM2(i,j,k  ,l) * VMTR_GSGAM2H(i,j,k,l) * GZXH(i,j,k,l)
          VMTR_C2WfactGz(i,j,k,I_b_GZXH,l) = GRD_bfac(k) * VMTR_RGSGAM2(i,j,k-1,l) * VMTR_GSGAM2H(i,j,k,l) * GZXH(i,j,k,l)
          VMTR_C2WfactGz(i,j,k,I_a_GZYH,l) = GRD_afac(k) * VMTR_RGSGAM2(i,j,k  ,l) * VMTR_GSGAM2H(i,j,k,l) * GZYH(i,j,k,l)
          VMTR_C2WfactGz(i,j,k,I_b_GZYH,l) = GRD_bfac(k) * VMTR_RGSGAM2(i,j,k-1,l) * VMTR_GSGAM2H(i,j,k,l) * GZYH(i,j,k,l)
          VMTR_C2WfactGz(i,j,k,I_a_GZZH,l) = GRD_afac(k) * VMTR_RGSGAM2(i,j,k  ,l) * VMTR_GSGAM2H(i,j,k,l) * GZZH(i,j,k,l)
          VMTR_C2WfactGz(i,j,k,I_b_GZZH,l) = GRD_bfac(k) * VMTR_RGSGAM2(i,j,k-1,l) * VMTR_GSGAM2H(i,j,k,l) * GZZH(i,j,k,l)
       enddo
       enddo
       enddo

       do j = 1, ADM_jall
       do i = 1, ADM_iall
          VMTR_C2WfactGz(i,j,ADM_kmin-1,I_a_GZXH,l) = 0.0_RP
          VMTR_C2WfactGz(i,j,ADM_kmin-1,I_b_GZXH,l) = 0.0_RP
          VMTR_C2WfactGz(i,j,ADM_kmin-1,I_a_GZYH,l) = 0.0_RP
          VMTR_C2WfactGz(i,j,ADM_kmin-1,I_b_GZYH,l) = 0.0_RP
          VMTR_C2WfactGz(i,j,ADM_kmin-1,I_a_GZZH,l) = 0.0_RP
          VMTR_C2WfactGz(i,j,ADM_kmin-1,I_b_GZZH,l) = 0.0_RP
       enddo
       enddo
    enddo

    !--- calculation of volume, geopotential
    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do j = 1, ADM_jall
    do i = 1, ADM_iall
       VMTR_VOLUME(i,j,k,l) = GMTR_area(i,j,l) * VMTR_GSGAM2(i,j,k,l) * GRD_dgz(k)

       VMTR_PHI   (i,j,k,l) = GRD_vz(i,j,k,l,GRD_Z) * GRAV
    enddo
    enddo
    enddo
    enddo



    if ( ADM_have_pl ) then

       !---   G^1/2 = dz/dgz
       do l = 1, ADM_lall_pl
          !--- calculation of G^1/2 at full level
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             GSQRT_pl(g,k,l) = ( GRD_vz_pl(g,k+1,l,GRD_ZH) - GRD_vz_pl(g,k,l,GRD_ZH) ) / GRD_dgz(k)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             GSQRT_pl(g,ADM_kmin-1,l) = GSQRT_pl(g,ADM_kmin,l)
             GSQRT_pl(g,ADM_kmax+1,l) = GSQRT_pl(g,ADM_kmax,l)
          enddo
          !--- calculation of G^1/2 at half level
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             GSQRTH_pl(g,k,l) = ( GRD_vz_pl(g,k,l,GRD_Z) - GRD_vz_pl(g,k-1,l,GRD_Z) ) / GRD_dgzh(k)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             GSQRTH_pl(g,ADM_kmin-1,l) = GSQRTH_pl(g,ADM_kmin,l)
          enddo
       enddo

       !--- Gamma = (a+z) / a
       if ( deep ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             GAM_pl (g,k,l) = 1.0_RP + GRD_vz_pl(g,k,l,GRD_Z)  / GRD_rscale
             GAMH_pl(g,k,l) = 1.0_RP + GRD_vz_pl(g,k,l,GRD_ZH) / GRD_rscale
          enddo
          enddo
          enddo
       else
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             GAM_pl (g,k,l) = 1.0_RP
             GAMH_pl(g,k,l) = 1.0_RP
          enddo
          enddo
          enddo
       endif

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          VMTR_GAM2H_pl   (g,k,l) = GAMH_pl(g,k,l) * GAMH_pl(g,k,l)
          VMTR_GSGAM2_pl  (g,k,l) = GAM_pl (g,k,l) * GAM_pl (g,k,l) * GSQRT_pl (g,k,l)
          VMTR_GSGAM2H_pl (g,k,l) = GAMH_pl(g,k,l) * GAMH_pl(g,k,l) * GSQRTH_pl(g,k,l)

          VMTR_RGSQRTH_pl (g,k,l) = 1.0_RP / GSQRTH_pl(g,k,l)
          VMTR_RGAM_pl    (g,k,l) = 1.0_RP / GAM_pl (g,k,l)
          VMTR_RGAMH_pl   (g,k,l) = 1.0_RP / GAMH_pl(g,k,l)
          VMTR_RGSGAM2_pl (g,k,l) = 1.0_RP / VMTR_GSGAM2_pl (g,k,l)
          VMTR_RGSGAM2H_pl(g,k,l) = 1.0_RP / VMTR_GSGAM2H_pl(g,k,l)
       enddo
       enddo
       enddo

       ! full level <-> half level interpolation factor
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             VMTR_C2Wfact_pl(g,k,I_a,l) = GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l)
             VMTR_C2Wfact_pl(g,k,I_b,l) = GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             VMTR_C2Wfact_pl(g,ADM_kmin-1,I_a,l) = 0.0_RP
             VMTR_C2Wfact_pl(g,ADM_kmin-1,I_b,l) = 0.0_RP
          enddo

          do k = ADM_kmin-1, ADM_kmax
          do g = 1, ADM_gall_pl
             VMTR_W2Cfact_pl(g,k,I_c,l) = GRD_cfac(k) * VMTR_GSGAM2_pl(g,k,l) * VMTR_RGSGAM2H_pl(g,k+1,l)
             VMTR_W2Cfact_pl(g,k,I_d,l) = GRD_dfac(k) * VMTR_GSGAM2_pl(g,k,l) * VMTR_RGSGAM2H_pl(g,k  ,l)
          enddo
          enddo
          do g = 1, ADM_gall_pl
             VMTR_W2Cfact_pl(g,ADM_kmax+1,I_c,l) = 0.0_RP
             VMTR_W2Cfact_pl(g,ADM_kmax+1,I_d,l) = 0.0_RP
          enddo
       enddo

       ! full level <-> half level interpolation factor with Gz

       !--- Gz(X) = - JX / G^1/2
       !--- Gz(Y) = - JY / G^1/2
       !--- Gz(Z) = - JZ / G^1/2
       do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             GZXH_pl(g,k,l) = -var_pl(g,k,l,JXH) / GSQRTH_pl(g,k,l)
             GZYH_pl(g,k,l) = -var_pl(g,k,l,JYH) / GSQRTH_pl(g,k,l)
             GZZH_pl(g,k,l) = -var_pl(g,k,l,JZH) / GSQRTH_pl(g,k,l)
             GZX_pl (g,k,l) = -var_pl(g,k,l,JX)  / GSQRT_pl (g,k,l)
             GZY_pl (g,k,l) = -var_pl(g,k,l,JY)  / GSQRT_pl (g,k,l)
             GZZ_pl (g,k,l) = -var_pl(g,k,l,JZ)  / GSQRT_pl (g,k,l)
          enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             VMTR_C2WfactGz_pl(g,k,I_a_GZXH,l) = GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l) * GZXH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_b_GZXH,l) = GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l) * GZXH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_a_GZYH,l) = GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l) * GZYH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_b_GZYH,l) = GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l) * GZYH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_a_GZZH,l) = GRD_afac(k) * VMTR_RGSGAM2_pl(g,k  ,l) * VMTR_GSGAM2H_pl(g,k,l) * GZZH_pl(g,k,l)
             VMTR_C2WfactGz_pl(g,k,I_b_GZZH,l) = GRD_bfac(k) * VMTR_RGSGAM2_pl(g,k-1,l) * VMTR_GSGAM2H_pl(g,k,l) * GZZH_pl(g,k,l)
          enddo
          enddo

          do g = 1, ADM_gall_pl
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_a_GZXH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_b_GZXH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_a_GZYH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_b_GZYH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_a_GZZH,l) = 0.0_RP
             VMTR_C2WfactGz_pl(g,ADM_kmin-1,I_b_GZZH,l) = 0.0_RP
          enddo
       enddo

       !--- calculation of volume, geopotential
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          VMTR_VOLUME_pl(g,k,l) = GMTR_area_pl(g,l) * VMTR_GSGAM2_pl(g,k,l) * GRD_dgz(k)

          VMTR_PHI_pl   (g,k,l) = GRD_vz_pl(g,k,l,GRD_Z) * GRAV
       enddo
       enddo
       enddo
    else
       VMTR_GAM2H_pl    (:,:,:)   = 0.0_RP
       VMTR_GSGAM2_pl   (:,:,:)   = 0.0_RP
       VMTR_GSGAM2H_pl  (:,:,:)   = 0.0_RP
       VMTR_RGSQRTH_pl  (:,:,:)   = 0.0_RP
       VMTR_RGAM_pl     (:,:,:)   = 0.0_RP
       VMTR_RGAMH_pl    (:,:,:)   = 0.0_RP
       VMTR_RGSGAM2_pl  (:,:,:)   = 0.0_RP
       VMTR_RGSGAM2H_pl (:,:,:)   = 0.0_RP
       VMTR_W2Cfact_pl  (:,:,:,:) = 0.0_RP
       VMTR_C2Wfact_pl  (:,:,:,:) = 0.0_RP
       VMTR_C2WfactGz_pl(:,:,:,:) = 0.0_RP
       VMTR_VOLUME_pl   (:,:,:)   = 0.0_RP
       VMTR_PHI_pl      (:,:,:)   = 0.0_RP
    endif

    return
  end subroutine VMTR_setup

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_GAM2H( &
       IJ_VMTR_GAM2H,   &
       IJ_VMTR_GAM2H_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_GAM2H   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_GAM2H_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_GAM2H(g,k,l) = VMTR_GAM2H(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_GAM2H_pl(:,:,:) = VMTR_GAM2H_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_GAM2H

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_GSGAM2( &
       IJ_VMTR_GSGAM2,   &
       IJ_VMTR_GSGAM2_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_GSGAM2   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_GSGAM2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_GSGAM2(g,k,l) = VMTR_GSGAM2(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_GSGAM2_pl(:,:,:) = VMTR_GSGAM2_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_GSGAM2

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_GSGAM2H( &
       IJ_VMTR_GSGAM2H,   &
       IJ_VMTR_GSGAM2H_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_GSGAM2H   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_GSGAM2H_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_GSGAM2H(g,k,l) = VMTR_GSGAM2H(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_GSGAM2H_pl(:,:,:) = VMTR_GSGAM2H_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_GSGAM2H

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_RGSQRTH( &
       IJ_VMTR_RGSQRTH,   &
       IJ_VMTR_RGSQRTH_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_RGSQRTH   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_RGSQRTH_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_RGSQRTH(g,k,l) = VMTR_RGSQRTH(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_RGSQRTH_pl(:,:,:) = VMTR_RGSQRTH_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_RGSQRTH

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_RGAM( &
       IJ_VMTR_RGAM,   &
       IJ_VMTR_RGAM_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_RGAM   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_RGAM_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_RGAM(g,k,l) = VMTR_RGAM(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_RGAM_pl(:,:,:) = VMTR_RGAM_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_RGAM

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_RGAMH( &
       IJ_VMTR_RGAMH,   &
       IJ_VMTR_RGAMH_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_RGAMH   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_RGAMH_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_RGAMH(g,k,l) = VMTR_RGAMH(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_RGAMH_pl(:,:,:) = VMTR_RGAMH_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_RGAMH

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_RGSGAM2( &
       IJ_VMTR_RGSGAM2,   &
       IJ_VMTR_RGSGAM2_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_RGSGAM2   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_RGSGAM2_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_RGSGAM2(g,k,l) = VMTR_RGSGAM2(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_RGSGAM2_pl(:,:,:) = VMTR_RGSGAM2_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_RGSGAM2

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_RGSGAM2H( &
       IJ_VMTR_RGSGAM2H,   &
       IJ_VMTR_RGSGAM2H_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_RGSGAM2H   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_RGSGAM2H_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_RGSGAM2H(g,k,l) = VMTR_RGSGAM2H(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_RGSGAM2H_pl(:,:,:) = VMTR_RGSGAM2H_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_RGSGAM2H

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_W2Cfact( &
       IJ_VMTR_W2Cfact,   &
       IJ_VMTR_W2Cfact_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_W2Cfact   (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_W2Cfact_pl(ADM_gall_pl,ADM_kall,2,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_W2Cfact(g,k,1,l) = VMTR_W2Cfact(i,j,k,1,l)
          IJ_VMTR_W2Cfact(g,k,2,l) = VMTR_W2Cfact(i,j,k,2,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_W2Cfact_pl(:,:,:,:) = VMTR_W2Cfact_pl(:,:,:,:)

    return
  end subroutine VMTR_getIJ_W2Cfact

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_C2Wfact( &
       IJ_VMTR_C2Wfact,   &
       IJ_VMTR_C2Wfact_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_C2Wfact   (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_C2Wfact_pl(ADM_gall_pl,ADM_kall,2,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_C2Wfact(g,k,1,l) = VMTR_C2Wfact(i,j,k,1,l)
          IJ_VMTR_C2Wfact(g,k,2,l) = VMTR_C2Wfact(i,j,k,2,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_C2Wfact_pl(:,:,:,:) = VMTR_C2Wfact_pl(:,:,:,:)

    return
  end subroutine VMTR_getIJ_C2Wfact

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_C2WfactGz( &
       IJ_VMTR_C2WfactGz,   &
       IJ_VMTR_C2WfactGz_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_C2WfactGz   (ADM_gall   ,ADM_kall,6,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_C2WfactGz_pl(ADM_gall_pl,ADM_kall,6,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_C2WfactGz(g,k,1,l) = VMTR_C2WfactGz(i,j,k,1,l)
          IJ_VMTR_C2WfactGz(g,k,2,l) = VMTR_C2WfactGz(i,j,k,2,l)
          IJ_VMTR_C2WfactGz(g,k,3,l) = VMTR_C2WfactGz(i,j,k,3,l)
          IJ_VMTR_C2WfactGz(g,k,4,l) = VMTR_C2WfactGz(i,j,k,4,l)
          IJ_VMTR_C2WfactGz(g,k,5,l) = VMTR_C2WfactGz(i,j,k,5,l)
          IJ_VMTR_C2WfactGz(g,k,6,l) = VMTR_C2WfactGz(i,j,k,6,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_C2WfactGz_pl(:,:,:,:) = VMTR_C2WfactGz_pl(:,:,:,:)

    return
  end subroutine VMTR_getIJ_C2WfactGz

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_VOLUME( &
       IJ_VMTR_VOLUME,   &
       IJ_VMTR_VOLUME_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_VOLUME   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_VOLUME_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_VOLUME(g,k,l) = VMTR_VOLUME(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_VOLUME_pl(:,:,:) = VMTR_VOLUME_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_VOLUME

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getIJ_PHI( &
       IJ_VMTR_PHI,   &
       IJ_VMTR_PHI_pl )
    implicit none

    real(RP), intent(out) :: IJ_VMTR_PHI   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: IJ_VMTR_PHI_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = 1, ADM_jall
       do i = 1, ADM_iall
          IJ_VMTR_PHI(g,k,l) = VMTR_PHI(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo
    IJ_VMTR_PHI_pl(:,:,:) = VMTR_PHI_pl(:,:,:)

    return
  end subroutine VMTR_getIJ_PHI

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getin_GSGAM2( &
       in_VMTR_GSGAM2 )
    use mod_adm, only: &
       ADM_gall_in, &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_imin,    &
       ADM_imax
    implicit none

    real(RP), intent(out) :: in_VMTR_GSGAM2(ADM_gall_in,ADM_kall,ADM_lall)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = ADM_jmin, ADM_jmax+1
       do i = ADM_imin, ADM_imax+1
          in_VMTR_GSGAM2(g,k,l) = VMTR_GSGAM2(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo

    return
  end subroutine VMTR_getin_GSGAM2

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getin_GSGAM2H( &
       in_VMTR_GSGAM2H )
    use mod_adm, only: &
       ADM_gall_in, &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_imin,    &
       ADM_imax
    implicit none

    real(RP), intent(out) :: in_VMTR_GSGAM2H(ADM_gall_in,ADM_kall,ADM_lall)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = ADM_jmin, ADM_jmax+1
       do i = ADM_imin, ADM_imax+1
          in_VMTR_GSGAM2H(g,k,l) = VMTR_GSGAM2H(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo

    return
  end subroutine VMTR_getin_GSGAM2H

  !-----------------------------------------------------------------------------
  !> Tentative Converter
  subroutine VMTR_getin_PHI( &
       in_VMTR_PHI )
    use mod_adm, only: &
       ADM_gall_in, &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_imin,    &
       ADM_imax
    implicit none

    real(RP), intent(out) :: in_VMTR_PHI(ADM_gall_in,ADM_kall,ADM_lall)

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do k = 1, ADM_kall
       g = 1
       do j = ADM_jmin, ADM_jmax+1
       do i = ADM_imin, ADM_imax+1
          in_VMTR_PHI(g,k,l) = VMTR_PHI(i,j,k,l)
          g = g + 1
       enddo
       enddo
    enddo
    enddo

    return
  end subroutine VMTR_getin_PHI

end module mod_vmtr
