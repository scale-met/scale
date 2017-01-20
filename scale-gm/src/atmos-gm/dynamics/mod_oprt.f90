!-------------------------------------------------------------------------------
!> Module operator
!!
!! @par Description
!!          This module contains the subroutines for differential operators
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_oprt
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof

  use mod_adm, only: &
     ADM_nxyz,           &
     TI    => ADM_TI,    &
     TJ    => ADM_TJ,    &
     AI    => ADM_AI,    &
     AIJ   => ADM_AIJ,   &
     AJ    => ADM_AJ,    &
     K0    => ADM_KNONE, &
     vlink => ADM_vlink, &
     ADM_lall,           &
     ADM_lall_pl,        &
     ADM_kall,           &
     ADM_jall,           &
     ADM_iall,           &
     ADM_gall_pl
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OPRT_setup
  public :: OPRT_divergence
  public :: OPRT_gradient
  public :: OPRT_laplacian
  public :: OPRT_diffusion
  public :: OPRT_horizontalize_vec
  public :: OPRT_rotation
  public :: OPRT_divdamp

  interface OPRT_horizontalize_vec
     module procedure OPRT_horizontalize_vec_ij
     module procedure OPRT_horizontalize_vec_ixj
  end interface OPRT_horizontalize_vec

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
#ifdef _FIXEDINDEX_
!  real(RP), public              :: OPRT_coef_div    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
!  real(RP), public              :: OPRT_coef_div_pl (ADM_nxyz,         0:vlink,ADM_lall_pl)
!  real(RP), public              :: OPRT_coef_rot    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
!  real(RP), public              :: OPRT_coef_rot_pl (ADM_nxyz,         0:vlink,ADM_lall_pl)
!  real(RP), public              :: OPRT_coef_grad   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
!  real(RP), public              :: OPRT_coef_grad_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)
!  real(RP), public              :: OPRT_coef_lap    (         ADM_gall,0:6    ,ADM_lall   )
!  real(RP), public              :: OPRT_coef_lap_pl (                  0:vlink,ADM_lall_pl)
!  real(RP), public              :: OPRT_coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
!  real(RP), public              :: OPRT_coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
!  real(RP), public              :: OPRT_coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   )
!  real(RP), public              :: OPRT_coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_div    (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
  real(RP), public              :: OPRT_coef_div_pl (                  0:vlink,ADM_nxyz,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_rot    (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
  real(RP), public              :: OPRT_coef_rot_pl (                  0:vlink,ADM_nxyz,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_grad   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
  real(RP), public              :: OPRT_coef_grad_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl)
  real(RP), public              :: OPRT_coef_lap    (ADM_iall,ADM_jall,0:6    ,         ADM_lall   )
  real(RP), public              :: OPRT_coef_lap_pl (                  0:vlink,         ADM_lall_pl)
  real(RP), public              :: OPRT_coef_intp   (ADM_iall,ADM_jall,1:3    ,ADM_nxyz,TI:TJ,ADM_lall   )
  real(RP), public              :: OPRT_coef_intp_pl(ADM_gall_pl      ,1:3    ,ADM_nxyz,      ADM_lall_pl)
  real(RP), public              :: OPRT_coef_diff   (ADM_iall,ADM_jall,1:6    ,ADM_nxyz,ADM_lall   )
  real(RP), public              :: OPRT_coef_diff_pl(                  1:vlink,ADM_nxyz,ADM_lall_pl)
#else
  real(RP), public, allocatable :: OPRT_coef_div    (:,:,:,:,:)   ! coefficient for divergence operator
  real(RP), public, allocatable :: OPRT_coef_div_pl (:,:,:)
  real(RP), public, allocatable :: OPRT_coef_rot    (:,:,:,:,:)   ! coefficient for rotation operator
  real(RP), public, allocatable :: OPRT_coef_rot_pl (:,:,:)
  real(RP), public, allocatable :: OPRT_coef_grad   (:,:,:,:,:)   ! coefficient for gradient operator
  real(RP), public, allocatable :: OPRT_coef_grad_pl(:,:,:)
  real(RP), public, allocatable :: OPRT_coef_lap    (:,:,:,:)     ! coefficient for laplacian operator
  real(RP), public, allocatable :: OPRT_coef_lap_pl (:,:)
  real(RP), public, allocatable :: OPRT_coef_intp   (:,:,:,:,:,:) ! coefficient for interpolation operator
  real(RP), public, allocatable :: OPRT_coef_intp_pl(:,:,:,:)
  real(RP), public, allocatable :: OPRT_coef_diff   (:,:,:,:,:)
  real(RP), public, allocatable :: OPRT_coef_diff_pl(:,:,:)
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  character(len=H_LONG),  private :: OPRT_fname   = ''
  character(len=H_SHORT), private :: OPRT_io_mode = 'ADVANCED'

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OPRT_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_gmtr, only: &
       GMTR_p,    &
       GMTR_p_pl, &
       GMTR_t,    &
       GMTR_t_pl, &
       GMTR_a,    &
       GMTR_a_pl, &
       GMTR_p_nmax, &
       GMTR_t_nmax, &
       GMTR_a_nmax
    implicit none

    namelist / OPRTPARAM / &
       OPRT_io_mode, &
       OPRT_fname

    real(RP) :: IxJ_GMTR_p(ADM_iall,ADM_jall,K0,ADM_lall,      GMTR_p_nmax)
    real(RP) :: IxJ_GMTR_t(ADM_iall,ADM_jall,K0,ADM_lall,TI:TJ,GMTR_t_nmax)
    real(RP) :: IxJ_GMTR_a(ADM_iall,ADM_jall,K0,ADM_lall,AI:AJ,GMTR_a_nmax)

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[oprt]/Category[common share]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=OPRTPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** OPRTPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist OPRTPARAM. STOP.'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=OPRTPARAM)

#ifndef _FIXEDINDEX_
!    allocate( OPRT_coef_div    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   ) )
!    allocate( OPRT_coef_div_pl (ADM_nxyz,         0:vlink,ADM_lall_pl) )
!    allocate( OPRT_coef_rot    (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   ) )
!    allocate( OPRT_coef_rot_pl (ADM_nxyz,         0:vlink,ADM_lall_pl) )
!    allocate( OPRT_coef_grad   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   ) )
!    allocate( OPRT_coef_grad_pl(ADM_nxyz,         0:vlink,ADM_lall_pl) )
!    allocate( OPRT_coef_lap    (         ADM_gall,0:6    ,ADM_lall   ) )
!    allocate( OPRT_coef_lap_pl (                  0:vlink,ADM_lall_pl) )
!    allocate( OPRT_coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   ) )
!    allocate( OPRT_coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl) )
!    allocate( OPRT_coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   ) )
!    allocate( OPRT_coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl) )
    allocate( OPRT_coef_div    (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   ) )
    allocate( OPRT_coef_div_pl (                  0:vlink,ADM_nxyz,ADM_lall_pl) )
    allocate( OPRT_coef_rot    (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   ) )
    allocate( OPRT_coef_rot_pl (                  0:vlink,ADM_nxyz,ADM_lall_pl) )
    allocate( OPRT_coef_grad   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   ) )
    allocate( OPRT_coef_grad_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl) )
    allocate( OPRT_coef_lap    (ADM_iall,ADM_jall,0:6    ,         ADM_lall   ) )
    allocate( OPRT_coef_lap_pl (                  0:vlink,         ADM_lall_pl) )
    allocate( OPRT_coef_intp   (ADM_iall,ADM_jall,1:3    ,ADM_nxyz,TI:TJ,ADM_lall   ) )
    allocate( OPRT_coef_intp_pl(ADM_gall_pl      ,1:3    ,ADM_nxyz,      ADM_lall_pl) )
    allocate( OPRT_coef_diff   (ADM_iall,ADM_jall,1:6    ,ADM_nxyz,ADM_lall   ) )
    allocate( OPRT_coef_diff_pl(                  1:vlink,ADM_nxyz,ADM_lall_pl) )
#endif

    IxJ_GMTR_p = reshape(GMTR_p,shape(IxJ_GMTR_p))
    IxJ_GMTR_t = reshape(GMTR_t,shape(IxJ_GMTR_t))
    IxJ_GMTR_a = reshape(GMTR_a,shape(IxJ_GMTR_a))

    call OPRT_divergence_setup( IxJ_GMTR_p    (:,:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_t    (:,:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_a    (:,:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_div (:,:,:,:,:),   OPRT_coef_div_pl (:,:,:)    ) ! [OUT]

    call OPRT_rotation_setup  ( IxJ_GMTR_p    (:,:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_t    (:,:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_a    (:,:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_rot (:,:,:,:,:),   OPRT_coef_rot_pl (:,:,:)    ) ! [OUT]

    call OPRT_gradient_setup  ( IxJ_GMTR_p    (:,:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_t    (:,:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_a    (:,:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_grad(:,:,:,:,:),   OPRT_coef_grad_pl(:,:,:)    ) ! [OUT]

    call OPRT_laplacian_setup ( IxJ_GMTR_p    (:,:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_t    (:,:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_a    (:,:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_lap (:,:,:,:),     OPRT_coef_lap_pl (:,:)      ) ! [OUT]

    call OPRT_diffusion_setup ( IxJ_GMTR_p    (:,:,:,:,:),   GMTR_p_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_t    (:,:,:,:,:,:), GMTR_t_pl        (:,:,:,:), & ! [IN]
                                IxJ_GMTR_a    (:,:,:,:,:,:), GMTR_a_pl        (:,:,:,:), & ! [IN]
                                OPRT_coef_intp(:,:,:,:,:,:), OPRT_coef_intp_pl(:,:,:,:), & ! [OUT]
                                OPRT_coef_diff(:,:,:,:,:),   OPRT_coef_diff_pl(:,:,:)    ) ! [OUT]

    if ( OPRT_fname /= "" ) then
       call OPRT_output_coef( OPRT_fname )
    endif

    return
  end subroutine OPRT_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_div, coef_div_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_jmin,     &
       ADM_jmax,     &
       ADM_imin,     &
       ADM_imax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       W1      => GMTR_t_W1,    &
       W2      => GMTR_t_W2,    &
       W3      => GMTR_t_W3,    &
       HNX     => GMTR_a_HNX,   &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_div   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
!    real(RP), intent(out) :: coef_div_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_div   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_div_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, hn
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of divergence operator'

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       hn = d + HNX - 1

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          ! ij
!          coef_div(d,ij,0,l) = &
          coef_div(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_div(d,ij,1,l) = &
          coef_div(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_div(d,ij,2,l) = &
          coef_div(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_div(d,ij,3,l) = &
          coef_div(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP*GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_div(d,ij,4,l) = &
          coef_div(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_div(d,ij,5,l) = &
          coef_div(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_div(d,ij,6,l) = &
          coef_div(i,j,6,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          ! ij
!          coef_div(d,ij,0,l) = &
          coef_div(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_div(d,ij,1,l) = &
          coef_div(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_div(d,ij,2,l) = &
          coef_div(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_div(d,ij,3,l) = &
          coef_div(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_div(d,ij,4,l) = &
          coef_div(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_div(d,ij,5,l) = &
          coef_div(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_div(d,ij,6,l) = &
          coef_div(i,j,6,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       endif

    enddo ! loop d
    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,hn) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,hn) )
          enddo
!          coef_div_pl(d,0,l) &
          coef_div_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

!          coef_div_pl(d,v-1,l) &
             coef_div_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_divergence_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_rotation_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_rot, coef_rot_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_jmin,     &
       ADM_jmax,     &
       ADM_imin,     &
       ADM_imax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       W1      => GMTR_t_W1,    &
       W2      => GMTR_t_W2,    &
       W3      => GMTR_t_W3,    &
       HTX     => GMTR_a_HTX,   &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_rot   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
!    real(RP), intent(out) :: coef_rot_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_rot   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_rot_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, ht
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of rotation operator'

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       ht = d + HTX - 1

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          ! ij
!          coef_rot(d,ij,0,l) &
          coef_rot(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_rot(d,ij,1,l) &
          coef_rot(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_rot(d,ij,2,l) &
          coef_rot(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_rot(d,ij,3,l) &
          coef_rot(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                ) * 0.5_RP*GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_rot(d,ij,4,l) &
          coef_rot(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_rot(d,ij,5,l) &
          coef_rot(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_rot(d,ij,6,l) &
          coef_rot(i,j,6,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          ! ij
!          coef_rot(d,ij,0,l) &
          coef_rot(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_rot(d,ij,1,l) &
          coef_rot(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_rot(d,ij,2,l) &
          coef_rot(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_rot(d,ij,3,l) &
          coef_rot(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_rot(d,ij,4,l) &
          coef_rot(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_rot(d,ij,5,l) &
          coef_rot(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_rot(d,ij,6,l) &
          coef_rot(i,j,6,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       endif

    enddo ! loop d
    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          ht = d + HTX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,ht) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,ht) )
          enddo
!          coef_rot_pl(d,0,l) &
          coef_rot_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

!             coef_rot_pl(d,v-1,l) &
             coef_rot_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,ht) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,ht) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_rotation_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient_setup( &
       GMTR_p,    GMTR_p_pl,   &
       GMTR_t,    GMTR_t_pl,   &
       GMTR_a,    GMTR_a_pl,   &
       coef_grad, coef_grad_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_jmin,     &
       ADM_jmax,     &
       ADM_imin,     &
       ADM_imax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       W1      => GMTR_t_W1,    &
       W2      => GMTR_t_W2,    &
       W3      => GMTR_t_W3,    &
       HNX     => GMTR_a_HNX,   &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_grad   (ADM_nxyz,ADM_gall,0:6    ,ADM_lall   )
!    real(RP), intent(out) :: coef_grad_pl(ADM_nxyz,         0:vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_grad   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_grad_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, hn
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of gradient operator'

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       hn = d + HNX - 1

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          ! ij
!          coef_grad(d,ij,0,l) &
          coef_grad(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                   - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                   - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,hn)                     & ! P0 * b1
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,hn)                     & ! P0 * b2
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,hn)                     & ! P0 * b3
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,hn)                     & ! P0 * b4
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,hn)                     & ! P0 * b5
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,hn)                     & ! P0 * b6
                                 ) * 0.5_RP * GMTR_p(i  ,j  ,k0,l,P_RAREA)
          ! ip1j
!          coef_grad(d,ij,1,l) &
          coef_grad(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_grad(d,ij,2,l) &
          coef_grad(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_grad(d,ij,3,l) &
          coef_grad(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_grad(d,ij,4,l) &
          coef_grad(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_grad(d,ij,5,l) &
          coef_grad(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_grad(d,ij,6,l) &
          coef_grad(i,j,6,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          ! ij
!          coef_grad(d,ij,0,l) &
          coef_grad(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                   - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,hn)                     & ! P0 * b1
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,hn)                     & ! P0 * b2
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,hn)                     & ! P0 * b3
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,hn)                     & ! P0 * b4
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,hn)                     & ! P0 * b6
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_grad(d,ij,1,l) &
          coef_grad(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_grad(d,ij,2,l) &
          coef_grad(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_grad(d,ij,3,l) &
          coef_grad(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_grad(d,ij,4,l) &
          coef_grad(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_grad(d,ij,5,l) &
          coef_grad(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_grad(d,ij,6,l) &
          coef_grad(i,j,6,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       endif
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl

             coef = coef + 2.0_RP * ( GMTR_t_pl(ij,k0,l,W1) - 1.0_RP ) * GMTR_a_pl(ijp1,k0,l,hn)
          enddo
!          coef_grad_pl(d,0,l) &
          coef_grad_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

!             coef_grad_pl(d,v-1,l) &
             coef_grad_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                     ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_gradient_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_lap, coef_lap_pl )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_jmin,     &
       ADM_jmax,     &
       ADM_imin,     &
       ADM_imax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       T_RAREA => GMTR_t_RAREA, &
       HNX     => GMTR_a_HNX,   &
       TNX     => GMTR_a_TNX,   &
       TN2X    => GMTR_a_TN2X,  &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_lap   (ADM_iall,ADM_jall,0:6    ,ADM_lall   )
    real(RP), intent(out) :: coef_lap_pl(                  0:vlink,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    integer  :: i, j, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of laplacian operator'

    do l = 1, ADM_lall

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          coef_lap(i,j,:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn = d + HNX - 1
             tn = d + TNX - 1

             ! ij
             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j-1,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(i,j,5,l) = coef_lap(i,j,5,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             coef_lap(i,j,5,l) = coef_lap(i,j,5,l) &
                               + GMTR_t(i-1,j-1,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) )

             ! ijm1
             coef_lap(i,j,6,l) = coef_lap(i,j,6,l) &
                               + GMTR_t(i-1,j-1,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) )

             coef_lap(i,j,6,l) = coef_lap(i,j,6,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )
          enddo
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          coef_lap(i,j,:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn = d + HNX - 1
             tn = d + TNX - 1

             ! ij
             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(i,j,5,l) = coef_lap(i,j,5,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             ! ijm1
             coef_lap(i,j,6,l) = coef_lap(i,j,6,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )
          enddo
       endif

       coef_lap(:,:,0,l) = coef_lap(:,:,0,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,1,l) = coef_lap(:,:,1,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,2,l) = coef_lap(:,:,2,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,3,l) = coef_lap(:,:,3,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,4,l) = coef_lap(:,:,4,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,5,l) = coef_lap(:,:,5,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,6,l) = coef_lap(:,:,6,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
    enddo

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1,ADM_lall_pl

          coef_lap_pl(:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                   * ( - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) )

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                   * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) )
             enddo
          enddo ! d loop

          do v = ADM_gslf_pl, ADM_gmax_pl
             coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) * GMTR_p_pl(n,k0,l,P_RAREA) / 12.0_RP
          enddo

       enddo ! l loop
    endif

    return
  end subroutine OPRT_laplacian_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion_setup( &
       GMTR_p,    GMTR_p_pl,    &
       GMTR_t,    GMTR_t_pl,    &
       GMTR_a,    GMTR_a_pl,    &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_jmin,     &
       ADM_jmax,     &
       ADM_imin,     &
       ADM_imax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_gmtr, only: &
       P_RAREA => GMTR_p_RAREA, &
       T_RAREA => GMTR_t_RAREA, &
       HNX     => GMTR_a_HNX,   &
       TNX     => GMTR_a_TNX,   &
       TN2X    => GMTR_a_TN2X,  &
       GMTR_p_nmax,             &
       GMTR_t_nmax,             &
       GMTR_a_nmax,             &
       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
!    real(RP), intent(out) :: coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
!    real(RP), intent(out) :: coef_diff   (ADM_nxyz,ADM_gall,1:6    ,ADM_lall   )
!    real(RP), intent(out) :: coef_diff_pl(ADM_nxyz,         1:vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_intp   (ADM_iall,ADM_jall,1:3    ,ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(out) :: coef_intp_pl(ADM_gall_pl      ,1:3    ,ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(out) :: coef_diff   (ADM_iall,ADM_jall,1:6    ,ADM_nxyz,      ADM_lall   )
    real(RP), intent(out) :: coef_diff_pl(                  1:vlink,ADM_nxyz,      ADM_lall_pl)

    integer  :: ij, ijp1

    integer :: i, j, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** setup coefficient of diffusion operator'

    do l = 1, ADM_lall
       do d = 1, ADM_nxyz
          hn = d + HNX - 1
          tn = d + TNX - 1

          do j = ADM_jmin-1, ADM_jmax
          do i = ADM_imin-1, ADM_imax

!             coef_intp(d,ij,1,TI,l) = + GMTR_a(ij  ,k0,l,AIJ,tn) - GMTR_a(ij  ,k0,l,AI ,tn)
!             coef_intp(d,ij,2,TI,l) = - GMTR_a(ij  ,k0,l,AI ,tn) - GMTR_a(ip1j,k0,l,AJ ,tn)
!             coef_intp(d,ij,3,TI,l) = - GMTR_a(ip1j,k0,l,AJ ,tn) + GMTR_a(ij  ,k0,l,AIJ,tn)
             coef_intp(i,j,1,d,TI,l) = + GMTR_a(i  ,j  ,k0,l,AIJ,tn) - GMTR_a(i  ,j  ,k0,l,AI ,tn)
             coef_intp(i,j,2,d,TI,l) = - GMTR_a(i  ,j  ,k0,l,AI ,tn) - GMTR_a(i+1,j  ,k0,l,AJ ,tn)
             coef_intp(i,j,3,d,TI,l) = - GMTR_a(i+1,j  ,k0,l,AJ ,tn) + GMTR_a(i  ,j  ,k0,l,AIJ,tn)

!             coef_intp(d,ij,1,TJ,l) = + GMTR_a(ij  ,k0,l,AJ ,tn) - GMTR_a(ij  ,k0,l,AIJ,tn)
!             coef_intp(d,ij,2,TJ,l) = - GMTR_a(ij  ,k0,l,AIJ,tn) + GMTR_a(ijp1,k0,l,AI ,tn)
!             coef_intp(d,ij,3,TJ,l) = + GMTR_a(ijp1,k0,l,AI ,tn) + GMTR_a(ij  ,k0,l,AJ ,tn)
             coef_intp(i,j,1,d,TJ,l) = + GMTR_a(i  ,j  ,k0,l,AJ ,tn) - GMTR_a(i  ,j  ,k0,l,AIJ,tn)
             coef_intp(i,j,2,d,TJ,l) = - GMTR_a(i  ,j  ,k0,l,AIJ,tn) + GMTR_a(i  ,j+1,k0,l,AI ,tn)
             coef_intp(i,j,3,d,TJ,l) = + GMTR_a(i  ,j+1,k0,l,AI ,tn) + GMTR_a(i  ,j  ,k0,l,AJ ,tn)

!             coef_intp(d,ij,:,TI,l) = coef_intp(d,ij,:,TI,l) * 0.5_RP * GMTR_t(ij,k0,l,TI,T_RAREA)
!             coef_intp(d,ij,:,TJ,l) = coef_intp(d,ij,:,TJ,l) * 0.5_RP * GMTR_t(ij,k0,l,TJ,T_RAREA)
             coef_intp(i,j,:,d,TI,l) = coef_intp(i,j,:,d,TI,l) * 0.5_RP * GMTR_t(i,j,k0,l,TI,T_RAREA)
             coef_intp(i,j,:,d,TJ,l) = coef_intp(i,j,:,d,TJ,l) * 0.5_RP * GMTR_t(i,j,k0,l,TJ,T_RAREA)
          enddo
          enddo

          do j = ADM_jmin, ADM_jmax
          do i = ADM_imin, ADM_imax

!             coef_diff(d,ij,1,l) = + GMTR_a(i  ,j  ,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,2,l) = + GMTR_a(i  ,j  ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,3,l) = - GMTR_a(i-1,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,4,l) = - GMTR_a(i-1,j-1,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,5,l) = - GMTR_a(i  ,j-1,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,6,l) = + GMTR_a(i  ,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
             coef_diff(i,j,1,d,l) = + GMTR_a(i  ,j  ,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,2,d,l) = + GMTR_a(i  ,j  ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,3,d,l) = - GMTR_a(i-1,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,4,d,l) = - GMTR_a(i-1,j-1,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,5,d,l) = - GMTR_a(i  ,j-1,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,6,d,l) = + GMTR_a(i  ,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          enddo
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          coef_diff(ADM_imin,ADM_jmin,5,:,l) = 0.0_RP
       endif

    enddo ! l loop

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1, ADM_lall_pl

          coef_intp_pl(:,:,:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

!                coef_intp_pl(d,v,1,l) = - GMTR_a_pl(ijp1,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn )
!                coef_intp_pl(d,v,2,l) = + GMTR_a_pl(ij  ,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn2)
!                coef_intp_pl(d,v,3,l) = + GMTR_a_pl(ij  ,k0,l,tn2) - GMTR_a_pl(ijp1,k0,l,tn )
                coef_intp_pl(v,1,d,l) = - GMTR_a_pl(ijp1,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn )
                coef_intp_pl(v,2,d,l) = + GMTR_a_pl(ij  ,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn2)
                coef_intp_pl(v,3,d,l) = + GMTR_a_pl(ij  ,k0,l,tn2) - GMTR_a_pl(ijp1,k0,l,tn )

!                coef_intp_pl(d,v,:,l) = coef_intp_pl(d,v,:,l) * 0.5_RP * GMTR_t_pl(v,k0,l,T_RAREA)
                coef_intp_pl(v,:,d,l) = coef_intp_pl(v,:,d,l) * 0.5_RP * GMTR_t_pl(v,k0,l,T_RAREA)

!                coef_diff_pl(d,v-1,l) = GMTR_a_pl(v,k0,l,hn) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
                coef_diff_pl(v-1,d,l) = GMTR_a_pl(v,k0,l,hn) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
             enddo
          enddo
       enddo ! l loop
    endif

    return
  end subroutine OPRT_diffusion_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence( &
       scl,      scl_pl,     &
       vx,       vx_pl,      &
       vy,       vy_pl,      &
       vz,       vz_pl,      &
       coef_div, coef_div_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_imin,    &
       ADM_imax,    &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: scl        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: scl_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_div   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
    real(RP), intent(in)  :: coef_div_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: imin, imax, jmin, jmax, kall, lall

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_divergence',2)

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l), &
    !$omp shared(imin,imax,jmin,jmax,kall,lall,scl,vx,vy,vz,coef_div)
    do l = 1, lall
    do k = 1, kall
       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          scl(i,j,k,l) = coef_div(i,j,0,XDIR,l) * vx(i  ,j  ,k,l) &
                       + coef_div(i,j,1,XDIR,l) * vx(i+1,j  ,k,l) &
                       + coef_div(i,j,2,XDIR,l) * vx(i+1,j+1,k,l) &
                       + coef_div(i,j,3,XDIR,l) * vx(i  ,j+1,k,l) &
                       + coef_div(i,j,4,XDIR,l) * vx(i-1,j  ,k,l) &
                       + coef_div(i,j,5,XDIR,l) * vx(i-1,j-1,k,l) &
                       + coef_div(i,j,6,XDIR,l) * vx(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          scl(i,j,k,l) = scl(i,j,k,l) &
                       + coef_div(i,j,0,YDIR,l) * vy(i  ,j  ,k,l) &
                       + coef_div(i,j,1,YDIR,l) * vy(i+1,j  ,k,l) &
                       + coef_div(i,j,2,YDIR,l) * vy(i+1,j+1,k,l) &
                       + coef_div(i,j,3,YDIR,l) * vy(i  ,j+1,k,l) &
                       + coef_div(i,j,4,YDIR,l) * vy(i-1,j  ,k,l) &
                       + coef_div(i,j,5,YDIR,l) * vy(i-1,j-1,k,l) &
                       + coef_div(i,j,6,YDIR,l) * vy(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          scl(i,j,k,l) = scl(i,j,k,l) &
                       + coef_div(i,j,0,ZDIR,l) * vz(i  ,j  ,k,l) &
                       + coef_div(i,j,1,ZDIR,l) * vz(i+1,j  ,k,l) &
                       + coef_div(i,j,2,ZDIR,l) * vz(i+1,j+1,k,l) &
                       + coef_div(i,j,3,ZDIR,l) * vz(i  ,j+1,k,l) &
                       + coef_div(i,j,4,ZDIR,l) * vz(i-1,j  ,k,l) &
                       + coef_div(i,j,5,ZDIR,l) * vz(i-1,j-1,k,l) &
                       + coef_div(i,j,6,ZDIR,l) * vz(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do

       !$omp workshare
       scl(:,jmin-1,k,l) = 0.0_RP
       scl(:,jmax+1,k,l) = 0.0_RP
       scl(imin-1,:,k,l) = 0.0_RP
       scl(imax+1,:,k,l) = 0.0_RP
       !$omp end workshare
    enddo ! k loop
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          scl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( coef_div_pl(v-1,XDIR,l) * vx_pl(v,k,l) &
                                             + coef_div_pl(v-1,YDIR,l) * vy_pl(v,k,l) &
                                             + coef_div_pl(v-1,ZDIR,l) * vz_pl(v,k,l) )
          enddo
       enddo
       enddo
    else
       scl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_divergence',2)

    return
  end subroutine OPRT_divergence

  !-----------------------------------------------------------------------------
  subroutine OPRT_rotation( &
       scl,      scl_pl,     &
       vx,       vx_pl,      &
       vy,       vy_pl,      &
       vz,       vz_pl,      &
       coef_rot, coef_rot_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_imin,    &
       ADM_imax,    &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: scl        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: scl_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_rot   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
    real(RP), intent(in)  :: coef_rot_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: imin, imax, jmin, jmax, kall, lall

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_rotation',2)

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l), &
    !$omp shared(imin,imax,jmin,jmax,kall,lall,scl,vx,vy,vz,coef_rot)
    do l = 1, lall
    do k = 1, kall
       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          scl(i,j,k,l) = coef_rot(i,j,0,XDIR,l) * vx(i  ,j  ,k,l) &
                       + coef_rot(i,j,1,XDIR,l) * vx(i+1,j  ,k,l) &
                       + coef_rot(i,j,2,XDIR,l) * vx(i+1,j+1,k,l) &
                       + coef_rot(i,j,3,XDIR,l) * vx(i  ,j+1,k,l) &
                       + coef_rot(i,j,4,XDIR,l) * vx(i-1,j  ,k,l) &
                       + coef_rot(i,j,5,XDIR,l) * vx(i-1,j-1,k,l) &
                       + coef_rot(i,j,6,XDIR,l) * vx(i  ,j-1,k,l) &
                       + coef_rot(i,j,0,YDIR,l) * vy(i  ,j  ,k,l) &
                       + coef_rot(i,j,1,YDIR,l) * vy(i+1,j  ,k,l) &
                       + coef_rot(i,j,2,YDIR,l) * vy(i+1,j+1,k,l) &
                       + coef_rot(i,j,3,YDIR,l) * vy(i  ,j+1,k,l) &
                       + coef_rot(i,j,4,YDIR,l) * vy(i-1,j  ,k,l) &
                       + coef_rot(i,j,5,YDIR,l) * vy(i-1,j-1,k,l) &
                       + coef_rot(i,j,6,YDIR,l) * vy(i  ,j-1,k,l) &
                       + coef_rot(i,j,0,ZDIR,l) * vz(i  ,j  ,k,l) &
                       + coef_rot(i,j,1,ZDIR,l) * vz(i+1,j  ,k,l) &
                       + coef_rot(i,j,2,ZDIR,l) * vz(i+1,j+1,k,l) &
                       + coef_rot(i,j,3,ZDIR,l) * vz(i  ,j+1,k,l) &
                       + coef_rot(i,j,4,ZDIR,l) * vz(i-1,j  ,k,l) &
                       + coef_rot(i,j,5,ZDIR,l) * vz(i-1,j-1,k,l) &
                       + coef_rot(i,j,6,ZDIR,l) * vz(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do

       !$omp workshare
       scl(:,jmin-1,k,l) = 0.0_RP
       scl(:,jmax+1,k,l) = 0.0_RP
       scl(imin-1,:,k,l) = 0.0_RP
       scl(imax+1,:,k,l) = 0.0_RP
       !$omp end workshare
    enddo ! k loop
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          scl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             scl_pl(n,k,l) = scl_pl(n,k,l) + ( coef_rot_pl(v-1,XDIR,l) * vx_pl(v,k,l) &
                                             + coef_rot_pl(v-1,YDIR,l) * vy_pl(v,k,l) &
                                             + coef_rot_pl(v-1,ZDIR,l) * vz_pl(v,k,l) )
          enddo
       enddo
       enddo
    else
       scl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_rotation',2)

    return
  end subroutine OPRT_rotation

  !-----------------------------------------------------------------------------
  !> horizontal gradient operator
  subroutine OPRT_gradient( &
       grad,      grad_pl,     &
       scl,       scl_pl,      &
       coef_grad, coef_grad_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_imin,    &
       ADM_imax,    &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: grad        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ,ADM_nxyz)
    real(RP), intent(out) :: grad_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl,ADM_nxyz)
    real(RP), intent(in)  :: scl         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_grad   (ADM_iall,ADM_jall,0:6    ,ADM_nxyz,ADM_lall   )
    real(RP), intent(in)  :: coef_grad_pl(                  0:vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: imin, imax, jmin, jmax, kall, lall

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_gradient',2)

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l), &
    !$omp shared(imin,imax,jmin,jmax,kall,lall,grad,scl,coef_grad)
    do l = 1, lall
    do k = 1, kall
       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          grad(i,j,k,l,XDIR) = coef_grad(i,j,0,XDIR,l) * scl(i  ,j  ,k,l) &
                             + coef_grad(i,j,1,XDIR,l) * scl(i+1,j  ,k,l) &
                             + coef_grad(i,j,2,XDIR,l) * scl(i+1,j+1,k,l) &
                             + coef_grad(i,j,3,XDIR,l) * scl(i  ,j+1,k,l) &
                             + coef_grad(i,j,4,XDIR,l) * scl(i-1,j  ,k,l) &
                             + coef_grad(i,j,5,XDIR,l) * scl(i-1,j-1,k,l) &
                             + coef_grad(i,j,6,XDIR,l) * scl(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do nowait

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          grad(i,j,k,l,YDIR) = coef_grad(i,j,0,YDIR,l) * scl(i  ,j  ,k,l) &
                             + coef_grad(i,j,1,YDIR,l) * scl(i+1,j  ,k,l) &
                             + coef_grad(i,j,2,YDIR,l) * scl(i+1,j+1,k,l) &
                             + coef_grad(i,j,3,YDIR,l) * scl(i  ,j+1,k,l) &
                             + coef_grad(i,j,4,YDIR,l) * scl(i-1,j  ,k,l) &
                             + coef_grad(i,j,5,YDIR,l) * scl(i-1,j-1,k,l) &
                             + coef_grad(i,j,6,YDIR,l) * scl(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do nowait

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          grad(i,j,k,l,ZDIR) = coef_grad(i,j,0,ZDIR,l) * scl(i  ,j  ,k,l) &
                             + coef_grad(i,j,1,ZDIR,l) * scl(i+1,j  ,k,l) &
                             + coef_grad(i,j,2,ZDIR,l) * scl(i+1,j+1,k,l) &
                             + coef_grad(i,j,3,ZDIR,l) * scl(i  ,j+1,k,l) &
                             + coef_grad(i,j,4,ZDIR,l) * scl(i-1,j  ,k,l) &
                             + coef_grad(i,j,5,ZDIR,l) * scl(i-1,j-1,k,l) &
                             + coef_grad(i,j,6,ZDIR,l) * scl(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do

       !$omp workshare
       grad(:,jmin-1,k,l,XDIR) = 0.0_RP
       grad(:,jmin-1,k,l,YDIR) = 0.0_RP
       grad(:,jmin-1,k,l,ZDIR) = 0.0_RP
       grad(:,jmax+1,k,l,XDIR) = 0.0_RP
       grad(:,jmax+1,k,l,YDIR) = 0.0_RP
       grad(:,jmax+1,k,l,ZDIR) = 0.0_RP
       grad(imin-1,:,k,l,XDIR) = 0.0_RP
       grad(imin-1,:,k,l,YDIR) = 0.0_RP
       grad(imin-1,:,k,l,ZDIR) = 0.0_RP
       grad(imax+1,:,k,l,XDIR) = 0.0_RP
       grad(imax+1,:,k,l,YDIR) = 0.0_RP
       grad(imax+1,:,k,l,ZDIR) = 0.0_RP
       !$omp end workshare
    enddo ! k loop
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          grad_pl(:,k,l,XDIR) = 0.0_RP
          grad_pl(:,k,l,YDIR) = 0.0_RP
          grad_pl(:,k,l,ZDIR) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             grad_pl(n,k,l,XDIR) = grad_pl(n,k,l,XDIR) + coef_grad_pl(v-1,XDIR,l) * scl_pl(v,k,l)
             grad_pl(n,k,l,YDIR) = grad_pl(n,k,l,YDIR) + coef_grad_pl(v-1,YDIR,l) * scl_pl(v,k,l)
             grad_pl(n,k,l,ZDIR) = grad_pl(n,k,l,ZDIR) + coef_grad_pl(v-1,ZDIR,l) * scl_pl(v,k,l)
          enddo
       enddo
       enddo
    else
       grad_pl(:,:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_gradient',2)

    return
  end subroutine OPRT_gradient

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian( &
       dscl,     dscl_pl,    &
       scl,      scl_pl,     &
       coef_lap, coef_lap_pl )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_imin,    &
       ADM_imax,    &
       ADM_jmin,    &
       ADM_jmax,    &
       ADM_gslf_pl, &
       ADM_gmax_pl
    implicit none

    real(RP), intent(out) :: dscl       (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl    (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_lap   (ADM_iall,ADM_jall,0:6    ,ADM_lall   )
    real(RP), intent(in)  :: coef_lap_pl(                  0:vlink,ADM_lall_pl)

    integer  :: imin, imax, jmin, jmax, kall, lall

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_laplacian',2)

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l), &
    !$omp shared(imin,imax,jmin,jmax,kall,lall,dscl,scl,coef_lap)
    do l = 1, lall
    do k = 1, kall
       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = coef_lap(i,j,0,l) * scl(i  ,j  ,k,l) &
                        + coef_lap(i,j,1,l) * scl(i+1,j  ,k,l) &
                        + coef_lap(i,j,2,l) * scl(i+1,j+1,k,l) &
                        + coef_lap(i,j,3,l) * scl(i  ,j+1,k,l) &
                        + coef_lap(i,j,4,l) * scl(i-1,j  ,k,l) &
                        + coef_lap(i,j,5,l) * scl(i-1,j-1,k,l) &
                        + coef_lap(i,j,6,l) * scl(i  ,j-1,k,l)
       enddo
       enddo
       !$omp end do

       !$omp workshare
       dscl(:,jmin-1,k,l) = 0.0_RP
       dscl(:,jmax+1,k,l) = 0.0_RP
       dscl(imin-1,:,k,l) = 0.0_RP
       dscl(imax+1,:,k,l) = 0.0_RP
       !$omp end workshare
    enddo ! k loop
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
          dscl_pl(:,k,l) = 0.0_RP
          do v = ADM_gslf_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + coef_lap_pl(v-1,l) * scl_pl(v,k,l)
          enddo
       enddo
       enddo
    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_laplacian',2)

    return
  end subroutine OPRT_laplacian

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion( &
       dscl,      dscl_pl,      &
       scl,       scl_pl,       &
       kh,        kh_pl,        &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_imin,     &
       ADM_imax,     &
       ADM_jmin,     &
       ADM_jmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: dscl        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh          (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl       (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_iall,ADM_jall,1:3,    ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_gall_pl      ,1:3,    ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_iall,ADM_jall,1:6,    ADM_nxyz,      ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(                  1:vlink,ADM_nxyz,      ADM_lall_pl)

    real(RP) :: vt   (ADM_iall,ADM_jall,ADM_nxyz,TI:TJ)
    real(RP) :: vt_pl(ADM_gall_pl      ,ADM_nxyz)

    integer  :: imin, imax, jmin, jmax, kall, lall
    integer  :: ij, ijp1, ijm1

    integer  :: i, j, k, l, n, v, d
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_diffusion',2)

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l,d), &
    !$omp shared(imin,imax,jmin,jmax,kall,lall,ADM_have_sgp,dscl,scl,kh,vt,coef_intp,coef_diff)
    do l = 1, lall
    do k = 1, kall
       !$omp do schedule(static) collapse(2)
       do d = XDIR, ZDIR
       do j = jmin-1, jmax
       do i = imin-1, imax
          vt(i,j,d,TI) = ( ( + 2.0_RP * coef_intp(i,j,1,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TI,l) ) * scl(i  ,j  ,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TI,l) &
                             + 2.0_RP * coef_intp(i,j,2,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TI,l) ) * scl(i+1,j  ,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TI,l) &
                             + 2.0_RP * coef_intp(i,j,3,d,TI,l) ) * scl(i+1,j+1,k,l) &
                         ) / 3.0_RP
       enddo
       enddo
       enddo
       !$omp end do nowait

       !$omp do schedule(static) collapse(2)
       do d = XDIR, ZDIR
       do j = jmin-1, jmax
       do i = imin-1, imax
          vt(i,j,d,TJ) = ( ( + 2.0_RP * coef_intp(i,j,1,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TJ,l) ) * scl(i  ,j  ,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TJ,l) &
                             + 2.0_RP * coef_intp(i,j,2,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TJ,l) ) * scl(i+1,j+1,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TJ,l) &
                             + 2.0_RP * coef_intp(i,j,3,d,TJ,l) ) * scl(i  ,j+1,k,l) &
                         ) / 3.0_RP
       enddo
       enddo
       enddo
       !$omp end do

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          vt(imin-1,jmin-1,:,TI) = vt(imin,jmin-1,:,TJ)
          !$omp end master
       endif

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = ( coef_diff(i,j,1,XDIR,l) * ( vt(i  ,j  ,XDIR,TI) + vt(i  ,j  ,XDIR,TJ) ) &
                          + coef_diff(i,j,1,YDIR,l) * ( vt(i  ,j  ,YDIR,TI) + vt(i  ,j  ,YDIR,TJ) ) &
                          + coef_diff(i,j,1,ZDIR,l) * ( vt(i  ,j  ,ZDIR,TI) + vt(i  ,j  ,ZDIR,TJ) ) &
                          ) * 0.5_RP * ( kh(i  ,j  ,k,l) + kh(i+1,j+1,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,2,XDIR,l) * ( vt(i  ,j  ,XDIR,TJ) + vt(i-1,j  ,XDIR,TI) ) &
                          + coef_diff(i,j,2,YDIR,l) * ( vt(i  ,j  ,YDIR,TJ) + vt(i-1,j  ,YDIR,TI) ) &
                          + coef_diff(i,j,2,ZDIR,l) * ( vt(i  ,j  ,ZDIR,TJ) + vt(i-1,j  ,ZDIR,TI) ) &
                          ) * 0.5_RP * ( kh(i  ,j  ,k,l) + kh(i  ,j+1,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,3,XDIR,l) * ( vt(i-1,j  ,XDIR,TI) + vt(i-1,j-1,XDIR,TJ) ) &
                          + coef_diff(i,j,3,YDIR,l) * ( vt(i-1,j  ,YDIR,TI) + vt(i-1,j-1,YDIR,TJ) ) &
                          + coef_diff(i,j,3,ZDIR,l) * ( vt(i-1,j  ,ZDIR,TI) + vt(i-1,j-1,ZDIR,TJ) ) &
                          ) * 0.5_RP * ( kh(i-1,j  ,k,l) + kh(i  ,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,4,XDIR,l) * ( vt(i-1,j-1,XDIR,TJ) + vt(i-1,j-1,XDIR,TI) ) &
                          + coef_diff(i,j,4,YDIR,l) * ( vt(i-1,j-1,YDIR,TJ) + vt(i-1,j-1,YDIR,TI) ) &
                          + coef_diff(i,j,4,ZDIR,l) * ( vt(i-1,j-1,ZDIR,TJ) + vt(i-1,j-1,ZDIR,TI) ) &
                          ) * 0.5_RP * ( kh(i-1,j-1,k,l) + kh(i  ,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,5,XDIR,l) * ( vt(i-1,j-1,XDIR,TI) + vt(i  ,j-1,XDIR,TJ) ) &
                          + coef_diff(i,j,5,YDIR,l) * ( vt(i-1,j-1,YDIR,TI) + vt(i  ,j-1,YDIR,TJ) ) &
                          + coef_diff(i,j,5,ZDIR,l) * ( vt(i-1,j-1,ZDIR,TI) + vt(i  ,j-1,ZDIR,TJ) ) &
                          ) * 0.5_RP * ( kh(i  ,j-1,k,l) + kh(i  ,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,6,XDIR,l) * ( vt(i  ,j-1,XDIR,TJ) + vt(i  ,j  ,XDIR,TI) ) &
                          + coef_diff(i,j,6,YDIR,l) * ( vt(i  ,j-1,YDIR,TJ) + vt(i  ,j  ,YDIR,TI) ) &
                          + coef_diff(i,j,6,ZDIR,l) * ( vt(i  ,j-1,ZDIR,TJ) + vt(i  ,j  ,ZDIR,TI) ) &
                          ) * 0.5_RP * ( kh(i  ,j  ,k,l) + kh(i+1,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp workshare
       dscl(:,jmin-1,k,l) = 0.0_RP
       dscl(:,jmax+1,k,l) = 0.0_RP
       dscl(imin-1,:,k,l) = 0.0_RP
       dscl(imax+1,:,k,l) = 0.0_RP
       !$omp end workshare
    enddo ! k loop
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             do d = 1, ADM_nxyz
                vt_pl(ij,d) = ( ( + 2.0_RP * coef_intp_pl(v,1,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,2,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(n   ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(v,1,d,l) &
                                  + 2.0_RP * coef_intp_pl(v,2,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(ij  ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(v,1,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,2,d,l) &
                                  + 2.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(ijp1,k,l) &
                              ) / 3.0_RP
             enddo
          enddo

          dscl_pl(:,k,l) = 0.0_RP

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

             dscl_pl(n,k,l) = dscl_pl(n,k,l) &
                            + ( coef_diff_pl(v-1,XDIR,l) * ( vt_pl(ijm1,XDIR) + vt_pl(ij,XDIR) ) &
                              + coef_diff_pl(v-1,YDIR,l) * ( vt_pl(ijm1,YDIR) + vt_pl(ij,YDIR) ) &
                              + coef_diff_pl(v-1,ZDIR,l) * ( vt_pl(ijm1,ZDIR) + vt_pl(ij,ZDIR) ) &
                              ) * 0.5_RP * ( kh_pl(n,k,l) + kh_pl(ij,k,l) )
          enddo

       enddo
       enddo
    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_diffusion',2)

    return
  end subroutine OPRT_diffusion

  !-----------------------------------------------------------------------------
  subroutine OPRT_divdamp( &
       ddivdx,    ddivdx_pl,    &
       ddivdy,    ddivdy_pl,    &
       ddivdz,    ddivdz_pl,    &
       vx,        vx_pl,        &
       vy,        vy_pl,        &
       vz,        vz_pl,        &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
    use mod_adm, only: &
       ADM_have_pl,  &
       ADM_have_sgp, &
       ADM_imin,     &
       ADM_imax,     &
       ADM_jmin,     &
       ADM_jmax,     &
       ADM_gslf_pl,  &
       ADM_gmin_pl,  &
       ADM_gmax_pl
    use mod_grd, only: &
       XDIR => GRD_XDIR, &
       YDIR => GRD_YDIR, &
       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: ddivdx      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdx_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdy      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdy_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: ddivdz      (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz_pl   (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vx          (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl       (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy          (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl       (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz          (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl       (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_iall,ADM_jall,1:3    ,ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_gall_pl      ,1:3    ,ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_iall,ADM_jall,1:6    ,ADM_nxyz,      ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(                  1:vlink,ADM_nxyz,      ADM_lall_pl)

    real(RP) :: sclt   (ADM_iall,ADM_jall,TI:TJ)
    real(RP) :: sclt_pl(ADM_gall_pl      )

    integer  :: imin, imax, jmin, jmax, kall, lall
    integer  :: ij, ijp1, ijm1

    integer  :: i, j, k, l, n, v
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_divdamp',2)

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l), &
    !$omp shared(imin,imax,jmin,jmax,kall,lall,ADM_have_sgp, &
    !$omp ddivdx,ddivdy,ddivdz,vx,vy,vz,sclt,coef_intp,coef_diff)
    do l = 1, lall
    do k = 1, kall

       !$omp do schedule(static)
       do j = jmin-1, jmax
       do i = imin-1, imax
          sclt(i,j,TI) = coef_intp(i,j,1,XDIR,TI,l) * vx(i  ,j  ,k,l) &
                       + coef_intp(i,j,2,XDIR,TI,l) * vx(i+1,j  ,k,l) &
                       + coef_intp(i,j,3,XDIR,TI,l) * vx(i+1,j+1,k,l) &
                       + coef_intp(i,j,1,YDIR,TI,l) * vy(i  ,j  ,k,l) &
                       + coef_intp(i,j,2,YDIR,TI,l) * vy(i+1,j  ,k,l) &
                       + coef_intp(i,j,3,YDIR,TI,l) * vy(i+1,j+1,k,l) &
                       + coef_intp(i,j,1,ZDIR,TI,l) * vz(i  ,j  ,k,l) &
                       + coef_intp(i,j,2,ZDIR,TI,l) * vz(i+1,j  ,k,l) &
                       + coef_intp(i,j,3,ZDIR,TI,l) * vz(i+1,j+1,k,l)
       enddo
       enddo
       !$omp end do nowait

       !$omp do schedule(static)
       do j = jmin-1, jmax
       do i = imin-1, imax
          sclt(i,j,TJ) = coef_intp(i,j,1,XDIR,TJ,l) * vx(i  ,j  ,k,l) &
                       + coef_intp(i,j,2,XDIR,TJ,l) * vx(i+1,j+1,k,l) &
                       + coef_intp(i,j,3,XDIR,TJ,l) * vx(i  ,j+1,k,l) &
                       + coef_intp(i,j,1,YDIR,TJ,l) * vy(i  ,j  ,k,l) &
                       + coef_intp(i,j,2,YDIR,TJ,l) * vy(i+1,j+1,k,l) &
                       + coef_intp(i,j,3,YDIR,TJ,l) * vy(i  ,j+1,k,l) &
                       + coef_intp(i,j,1,ZDIR,TJ,l) * vz(i  ,j  ,k,l) &
                       + coef_intp(i,j,2,ZDIR,TJ,l) * vz(i+1,j+1,k,l) &
                       + coef_intp(i,j,3,ZDIR,TJ,l) * vz(i  ,j+1,k,l)
       enddo
       enddo
       !$omp end do

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          sclt(imin-1,jmin-1,TI) = sclt(imin,jmin-1,TJ)
          !$omp end master
       endif

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          ddivdx(i,j,k,l) = coef_diff(i,j,1,XDIR,l) * ( sclt(i  ,j  ,TI) + sclt(i  ,j  ,TJ) ) &
                          + coef_diff(i,j,2,XDIR,l) * ( sclt(i  ,j  ,TJ) + sclt(i-1,j  ,TI) ) &
                          + coef_diff(i,j,3,XDIR,l) * ( sclt(i-1,j  ,TI) + sclt(i-1,j-1,TJ) ) &
                          + coef_diff(i,j,4,XDIR,l) * ( sclt(i-1,j-1,TJ) + sclt(i-1,j-1,TI) ) &
                          + coef_diff(i,j,5,XDIR,l) * ( sclt(i-1,j-1,TI) + sclt(i  ,j-1,TJ) ) &
                          + coef_diff(i,j,6,XDIR,l) * ( sclt(i  ,j-1,TJ) + sclt(i  ,j  ,TI) )
       enddo
       enddo
       !$omp end do nowait

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          ddivdy(i,j,k,l) = coef_diff(i,j,1,YDIR,l) * ( sclt(i  ,j  ,TI) + sclt(i  ,j  ,TJ) ) &
                          + coef_diff(i,j,2,YDIR,l) * ( sclt(i  ,j  ,TJ) + sclt(i-1,j  ,TI) ) &
                          + coef_diff(i,j,3,YDIR,l) * ( sclt(i-1,j  ,TI) + sclt(i-1,j-1,TJ) ) &
                          + coef_diff(i,j,4,YDIR,l) * ( sclt(i-1,j-1,TJ) + sclt(i-1,j-1,TI) ) &
                          + coef_diff(i,j,5,YDIR,l) * ( sclt(i-1,j-1,TI) + sclt(i  ,j-1,TJ) ) &
                          + coef_diff(i,j,6,YDIR,l) * ( sclt(i  ,j-1,TJ) + sclt(i  ,j  ,TI) )
       enddo
       enddo
       !$omp end do nowait

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          ddivdz(i,j,k,l) = coef_diff(i,j,1,ZDIR,l) * ( sclt(i  ,j  ,TI) + sclt(i  ,j  ,TJ) ) &
                          + coef_diff(i,j,2,ZDIR,l) * ( sclt(i  ,j  ,TJ) + sclt(i-1,j  ,TI) ) &
                          + coef_diff(i,j,3,ZDIR,l) * ( sclt(i-1,j  ,TI) + sclt(i-1,j-1,TJ) ) &
                          + coef_diff(i,j,4,ZDIR,l) * ( sclt(i-1,j-1,TJ) + sclt(i-1,j-1,TI) ) &
                          + coef_diff(i,j,5,ZDIR,l) * ( sclt(i-1,j-1,TI) + sclt(i  ,j-1,TJ) ) &
                          + coef_diff(i,j,6,ZDIR,l) * ( sclt(i  ,j-1,TJ) + sclt(i  ,j  ,TI) )
       enddo
       enddo
       !$omp end do

       !$omp workshare
       ddivdx(:,jmin-1,k,l) = 0.0_RP
       ddivdy(:,jmin-1,k,l) = 0.0_RP
       ddivdz(:,jmin-1,k,l) = 0.0_RP
       ddivdx(:,jmax+1,k,l) = 0.0_RP
       ddivdy(:,jmax+1,k,l) = 0.0_RP
       ddivdz(:,jmax+1,k,l) = 0.0_RP
       ddivdx(imin-1,:,k,l) = 0.0_RP
       ddivdy(imin-1,:,k,l) = 0.0_RP
       ddivdz(imin-1,:,k,l) = 0.0_RP
       ddivdx(imax+1,:,k,l) = 0.0_RP
       ddivdy(imax+1,:,k,l) = 0.0_RP
       ddivdz(imax+1,:,k,l) = 0.0_RP
       !$omp end workshare
    enddo ! k loop
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             sclt_pl(ij) = coef_intp_pl(v,1,XDIR,l) * vx_pl(n   ,k,l) &
                         + coef_intp_pl(v,2,XDIR,l) * vx_pl(ij  ,k,l) &
                         + coef_intp_pl(v,3,XDIR,l) * vx_pl(ijp1,k,l) &
                         + coef_intp_pl(v,1,YDIR,l) * vy_pl(n   ,k,l) &
                         + coef_intp_pl(v,2,YDIR,l) * vy_pl(ij  ,k,l) &
                         + coef_intp_pl(v,3,YDIR,l) * vy_pl(ijp1,k,l) &
                         + coef_intp_pl(v,1,ZDIR,l) * vz_pl(n   ,k,l) &
                         + coef_intp_pl(v,2,ZDIR,l) * vz_pl(ij  ,k,l) &
                         + coef_intp_pl(v,3,ZDIR,l) * vz_pl(ijp1,k,l)
          enddo

          ddivdx_pl(:,k,l) = 0.0_RP
          ddivdy_pl(:,k,l) = 0.0_RP
          ddivdz_pl(:,k,l) = 0.0_RP

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

             ddivdx_pl(n,k,l) = ddivdx_pl(n,k,l) + coef_diff_pl(v-1,XDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
             ddivdy_pl(n,k,l) = ddivdy_pl(n,k,l) + coef_diff_pl(v-1,YDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
             ddivdz_pl(n,k,l) = ddivdz_pl(n,k,l) + coef_diff_pl(v-1,ZDIR,l) * ( sclt_pl(ijm1) + sclt_pl(ij) )
          enddo
       enddo
       enddo
    else
       ddivdx_pl(:,:,:) = 0.0_RP
       ddivdy_pl(:,:,:) = 0.0_RP
       ddivdz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_divdamp',2)

    return
  end subroutine OPRT_divdamp

  !-----------------------------------------------------------------------------
  subroutine OPRT_horizontalize_vec_ij( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl  )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall
    use mod_grd, only: &
       GRD_XDIR,   &
       GRD_YDIR,   &
       GRD_ZDIR,   &
       GRD_rscale, &
       GRD_x,      &
       GRD_x_pl
    implicit none

    real(RP), intent(inout) :: vx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: prd
    integer  :: g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_horizontalize_vec',2)

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       prd = vx(g,k,l) * GRD_x(g,k0,l,GRD_XDIR) / GRD_rscale &
           + vy(g,k,l) * GRD_x(g,k0,l,GRD_YDIR) / GRD_rscale &
           + vz(g,k,l) * GRD_x(g,k0,l,GRD_ZDIR) / GRD_rscale

       vx(g,k,l) = vx(g,k,l) - prd * GRD_x(g,k0,l,GRD_XDIR) / GRD_rscale
       vy(g,k,l) = vy(g,k,l) - prd * GRD_x(g,k0,l,GRD_YDIR) / GRD_rscale
       vz(g,k,l) = vz(g,k,l) - prd * GRD_x(g,k0,l,GRD_ZDIR) / GRD_rscale
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          prd = vx_pl(g,k,l) * GRD_x_pl(g,k0,l,GRD_XDIR) / GRD_rscale &
              + vy_pl(g,k,l) * GRD_x_pl(g,k0,l,GRD_YDIR) / GRD_rscale &
              + vz_pl(g,k,l) * GRD_x_pl(g,k0,l,GRD_ZDIR) / GRD_rscale

          vx_pl(g,k,l) = vx_pl(g,k,l) - prd * GRD_x_pl(g,k0,l,GRD_XDIR) / GRD_rscale
          vy_pl(g,k,l) = vy_pl(g,k,l) - prd * GRD_x_pl(g,k0,l,GRD_YDIR) / GRD_rscale
          vz_pl(g,k,l) = vz_pl(g,k,l) - prd * GRD_x_pl(g,k0,l,GRD_ZDIR) / GRD_rscale
       enddo
       enddo
       enddo
    else
       vx_pl(:,:,:) = 0.0_RP
       vy_pl(:,:,:) = 0.0_RP
       vz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_horizontalize_vec',2)

    return
  end subroutine OPRT_horizontalize_vec_ij

  !-----------------------------------------------------------------------------
  subroutine OPRT_horizontalize_vec_ixj( &
       vx, vx_pl, &
       vy, vy_pl, &
       vz, vz_pl  )
    use mod_adm, only: &
       ADM_have_pl
    use mod_grd, only: &
       GRD_XDIR,   &
       GRD_YDIR,   &
       GRD_ZDIR,   &
       GRD_rscale, &
       ORIG_GRD_x => GRD_x, &
       GRD_x_pl
    implicit none

    real(RP), intent(inout) :: vx   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vx_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vy   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vy_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(inout) :: vz   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: vz_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)

    real(RP) :: GRD_x(ADM_iall,ADM_jall,k0,ADM_lall,ADM_nxyz)

    real(RP) :: prd
    integer  :: i, j, g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('OPRT_horizontalize_vec',2)

    GRD_x = reshape(ORIG_GRD_x,shape(GRD_x))

    do l = 1, ADM_lall
    do k = 1, ADM_kall
    do j = 1, ADM_jall
    do i = 1, ADM_iall
       prd = vx(i,j,k,l) * GRD_x(i,j,k0,l,GRD_XDIR) / GRD_rscale &
           + vy(i,j,k,l) * GRD_x(i,j,k0,l,GRD_YDIR) / GRD_rscale &
           + vz(i,j,k,l) * GRD_x(i,j,k0,l,GRD_ZDIR) / GRD_rscale

       vx(i,j,k,l) = vx(i,j,k,l) - prd * GRD_x(i,j,k0,l,GRD_XDIR) / GRD_rscale
       vy(i,j,k,l) = vy(i,j,k,l) - prd * GRD_x(i,j,k0,l,GRD_YDIR) / GRD_rscale
       vz(i,j,k,l) = vz(i,j,k,l) - prd * GRD_x(i,j,k0,l,GRD_ZDIR) / GRD_rscale
    enddo
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          prd = vx_pl(g,k,l) * GRD_x_pl(g,k0,l,GRD_XDIR) / GRD_rscale &
              + vy_pl(g,k,l) * GRD_x_pl(g,k0,l,GRD_YDIR) / GRD_rscale &
              + vz_pl(g,k,l) * GRD_x_pl(g,k0,l,GRD_ZDIR) / GRD_rscale

          vx_pl(g,k,l) = vx_pl(g,k,l) - prd * GRD_x_pl(g,k0,l,GRD_XDIR) / GRD_rscale
          vy_pl(g,k,l) = vy_pl(g,k,l) - prd * GRD_x_pl(g,k0,l,GRD_YDIR) / GRD_rscale
          vz_pl(g,k,l) = vz_pl(g,k,l) - prd * GRD_x_pl(g,k0,l,GRD_ZDIR) / GRD_rscale
       enddo
       enddo
       enddo
    else
       vx_pl(:,:,:) = 0.0_RP
       vy_pl(:,:,:) = 0.0_RP
       vz_pl(:,:,:) = 0.0_RP
    endif

    call PROF_rapend('OPRT_horizontalize_vec',2)

    return
  end subroutine OPRT_horizontalize_vec_ixj

  !-----------------------------------------------------------------------------
  subroutine OPRT_output_coef( &
       basename )
    use scale_process, only: &
       PRC_MPIstop
    use mod_io_param, only: &
       IO_REAL8, &
       IO_REAL4
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_gall,    &
       ADM_gall_1d
    use mod_comm, only: &
       COMM_data_transfer
    use mod_fio, only: &
       FIO_output
    implicit none

    character(len=*), intent(in) :: basename

    character(len=H_MID) :: desc = 'Coefficients info'

    real(RP) :: tmp   (ADM_gall   ,106,ADM_lall   ,1)
    real(RP) :: tmp_pl(ADM_gall_pl,106,ADM_lall_pl,1)

    integer :: dtype
    integer :: i, j, g, l
    !---------------------------------------------------------------------------

    if    ( RP == SP ) then
       dtype = IO_REAL4
    elseif( RP == DP ) then
       dtype = IO_REAL8
    endif

    do l = 1, ADM_lall
    do j = 1, ADM_jall
    do i = 1, ADM_iall
       g = (j-1) * ADM_gall_1d + i

       tmp(g,  1,l,1) = abs(OPRT_coef_div (i,j,0,1,l))
       tmp(g,  2,l,1) = abs(OPRT_coef_div (i,j,0,2,l))
       tmp(g,  3,l,1) = abs(OPRT_coef_div (i,j,0,3,l))
       tmp(g,  4,l,1) = abs(OPRT_coef_div (i,j,1,1,l))
       tmp(g,  5,l,1) = abs(OPRT_coef_div (i,j,1,2,l))
       tmp(g,  6,l,1) = abs(OPRT_coef_div (i,j,1,3,l))
       tmp(g,  7,l,1) = abs(OPRT_coef_div (i,j,2,1,l))
       tmp(g,  8,l,1) = abs(OPRT_coef_div (i,j,2,2,l))
       tmp(g,  9,l,1) = abs(OPRT_coef_div (i,j,2,3,l))
       tmp(g, 10,l,1) = abs(OPRT_coef_div (i,j,3,1,l))
       tmp(g, 11,l,1) = abs(OPRT_coef_div (i,j,3,2,l))
       tmp(g, 12,l,1) = abs(OPRT_coef_div (i,j,3,3,l))
       tmp(g, 13,l,1) = abs(OPRT_coef_div (i,j,4,1,l))
       tmp(g, 14,l,1) = abs(OPRT_coef_div (i,j,4,2,l))
       tmp(g, 15,l,1) = abs(OPRT_coef_div (i,j,4,3,l))
       tmp(g, 16,l,1) = abs(OPRT_coef_div (i,j,5,1,l))
       tmp(g, 17,l,1) = abs(OPRT_coef_div (i,j,5,2,l))
       tmp(g, 18,l,1) = abs(OPRT_coef_div (i,j,5,3,l))
       tmp(g, 19,l,1) = abs(OPRT_coef_div (i,j,6,1,l))
       tmp(g, 20,l,1) = abs(OPRT_coef_div (i,j,6,2,l))
       tmp(g, 21,l,1) = abs(OPRT_coef_div (i,j,6,3,l))
       tmp(g, 22,l,1) = abs(OPRT_coef_rot (i,j,0,1,l))
       tmp(g, 23,l,1) = abs(OPRT_coef_rot (i,j,0,2,l))
       tmp(g, 24,l,1) = abs(OPRT_coef_rot (i,j,0,3,l))
       tmp(g, 25,l,1) = abs(OPRT_coef_rot (i,j,1,1,l))
       tmp(g, 26,l,1) = abs(OPRT_coef_rot (i,j,1,2,l))
       tmp(g, 27,l,1) = abs(OPRT_coef_rot (i,j,1,3,l))
       tmp(g, 28,l,1) = abs(OPRT_coef_rot (i,j,2,1,l))
       tmp(g, 29,l,1) = abs(OPRT_coef_rot (i,j,2,2,l))
       tmp(g, 30,l,1) = abs(OPRT_coef_rot (i,j,2,3,l))
       tmp(g, 31,l,1) = abs(OPRT_coef_rot (i,j,3,1,l))
       tmp(g, 32,l,1) = abs(OPRT_coef_rot (i,j,3,2,l))
       tmp(g, 33,l,1) = abs(OPRT_coef_rot (i,j,3,3,l))
       tmp(g, 34,l,1) = abs(OPRT_coef_rot (i,j,4,1,l))
       tmp(g, 35,l,1) = abs(OPRT_coef_rot (i,j,4,2,l))
       tmp(g, 36,l,1) = abs(OPRT_coef_rot (i,j,4,3,l))
       tmp(g, 37,l,1) = abs(OPRT_coef_rot (i,j,5,1,l))
       tmp(g, 38,l,1) = abs(OPRT_coef_rot (i,j,5,2,l))
       tmp(g, 39,l,1) = abs(OPRT_coef_rot (i,j,5,3,l))
       tmp(g, 40,l,1) = abs(OPRT_coef_rot (i,j,6,1,l))
       tmp(g, 41,l,1) = abs(OPRT_coef_rot (i,j,6,2,l))
       tmp(g, 42,l,1) = abs(OPRT_coef_rot (i,j,6,3,l))
       tmp(g, 43,l,1) = abs(OPRT_coef_grad(i,j,0,1,l))
       tmp(g, 44,l,1) = abs(OPRT_coef_grad(i,j,0,2,l))
       tmp(g, 45,l,1) = abs(OPRT_coef_grad(i,j,0,3,l))
       tmp(g, 46,l,1) = abs(OPRT_coef_grad(i,j,1,1,l))
       tmp(g, 47,l,1) = abs(OPRT_coef_grad(i,j,1,2,l))
       tmp(g, 48,l,1) = abs(OPRT_coef_grad(i,j,1,3,l))
       tmp(g, 49,l,1) = abs(OPRT_coef_grad(i,j,2,1,l))
       tmp(g, 50,l,1) = abs(OPRT_coef_grad(i,j,2,2,l))
       tmp(g, 51,l,1) = abs(OPRT_coef_grad(i,j,2,3,l))
       tmp(g, 52,l,1) = abs(OPRT_coef_grad(i,j,3,1,l))
       tmp(g, 53,l,1) = abs(OPRT_coef_grad(i,j,3,2,l))
       tmp(g, 54,l,1) = abs(OPRT_coef_grad(i,j,3,3,l))
       tmp(g, 55,l,1) = abs(OPRT_coef_grad(i,j,4,1,l))
       tmp(g, 56,l,1) = abs(OPRT_coef_grad(i,j,4,2,l))
       tmp(g, 57,l,1) = abs(OPRT_coef_grad(i,j,4,3,l))
       tmp(g, 58,l,1) = abs(OPRT_coef_grad(i,j,5,1,l))
       tmp(g, 59,l,1) = abs(OPRT_coef_grad(i,j,5,2,l))
       tmp(g, 60,l,1) = abs(OPRT_coef_grad(i,j,5,3,l))
       tmp(g, 61,l,1) = abs(OPRT_coef_grad(i,j,6,1,l))
       tmp(g, 62,l,1) = abs(OPRT_coef_grad(i,j,6,2,l))
       tmp(g, 63,l,1) = abs(OPRT_coef_grad(i,j,6,3,l))
       tmp(g, 64,l,1) = abs(OPRT_coef_lap (i,j,0,  l))
       tmp(g, 65,l,1) = abs(OPRT_coef_lap (i,j,1,  l))
       tmp(g, 66,l,1) = abs(OPRT_coef_lap (i,j,2,  l))
       tmp(g, 67,l,1) = abs(OPRT_coef_lap (i,j,3,  l))
       tmp(g, 68,l,1) = abs(OPRT_coef_lap (i,j,4,  l))
       tmp(g, 69,l,1) = abs(OPRT_coef_lap (i,j,5,  l))
       tmp(g, 70,l,1) = abs(OPRT_coef_lap (i,j,6,  l))
       tmp(g, 71,l,1) = abs(OPRT_coef_intp(i,j,1,1,TI,l))
       tmp(g, 72,l,1) = abs(OPRT_coef_intp(i,j,1,2,TI,l))
       tmp(g, 73,l,1) = abs(OPRT_coef_intp(i,j,1,3,TI,l))
       tmp(g, 74,l,1) = abs(OPRT_coef_intp(i,j,2,1,TI,l))
       tmp(g, 75,l,1) = abs(OPRT_coef_intp(i,j,2,2,TI,l))
       tmp(g, 76,l,1) = abs(OPRT_coef_intp(i,j,2,3,TI,l))
       tmp(g, 77,l,1) = abs(OPRT_coef_intp(i,j,3,1,TI,l))
       tmp(g, 78,l,1) = abs(OPRT_coef_intp(i,j,3,2,TI,l))
       tmp(g, 79,l,1) = abs(OPRT_coef_intp(i,j,3,3,TI,l))
       tmp(g, 80,l,1) = abs(OPRT_coef_intp(i,j,1,1,TJ,l))
       tmp(g, 81,l,1) = abs(OPRT_coef_intp(i,j,1,2,TJ,l))
       tmp(g, 82,l,1) = abs(OPRT_coef_intp(i,j,1,3,TJ,l))
       tmp(g, 83,l,1) = abs(OPRT_coef_intp(i,j,2,1,TJ,l))
       tmp(g, 84,l,1) = abs(OPRT_coef_intp(i,j,2,2,TJ,l))
       tmp(g, 85,l,1) = abs(OPRT_coef_intp(i,j,2,3,TJ,l))
       tmp(g, 86,l,1) = abs(OPRT_coef_intp(i,j,3,1,TJ,l))
       tmp(g, 87,l,1) = abs(OPRT_coef_intp(i,j,3,2,TJ,l))
       tmp(g, 88,l,1) = abs(OPRT_coef_intp(i,j,3,3,TJ,l))
       tmp(g, 89,l,1) = abs(OPRT_coef_diff(i,j,1,1,l))
       tmp(g, 90,l,1) = abs(OPRT_coef_diff(i,j,1,2,l))
       tmp(g, 91,l,1) = abs(OPRT_coef_diff(i,j,1,3,l))
       tmp(g, 92,l,1) = abs(OPRT_coef_diff(i,j,2,1,l))
       tmp(g, 93,l,1) = abs(OPRT_coef_diff(i,j,2,2,l))
       tmp(g, 94,l,1) = abs(OPRT_coef_diff(i,j,2,3,l))
       tmp(g, 95,l,1) = abs(OPRT_coef_diff(i,j,3,1,l))
       tmp(g, 96,l,1) = abs(OPRT_coef_diff(i,j,3,2,l))
       tmp(g, 97,l,1) = abs(OPRT_coef_diff(i,j,3,3,l))
       tmp(g, 98,l,1) = abs(OPRT_coef_diff(i,j,4,1,l))
       tmp(g, 99,l,1) = abs(OPRT_coef_diff(i,j,4,2,l))
       tmp(g,100,l,1) = abs(OPRT_coef_diff(i,j,4,3,l))
       tmp(g,101,l,1) = abs(OPRT_coef_diff(i,j,5,1,l))
       tmp(g,102,l,1) = abs(OPRT_coef_diff(i,j,5,2,l))
       tmp(g,103,l,1) = abs(OPRT_coef_diff(i,j,5,3,l))
       tmp(g,104,l,1) = abs(OPRT_coef_diff(i,j,6,1,l))
       tmp(g,105,l,1) = abs(OPRT_coef_diff(i,j,6,2,l))
       tmp(g,106,l,1) = abs(OPRT_coef_diff(i,j,6,3,l))
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) Then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          tmp_pl(g,  1,l,1) = abs(OPRT_coef_div_pl (0,1,l))
          tmp_pl(g,  2,l,1) = abs(OPRT_coef_div_pl (0,2,l))
          tmp_pl(g,  3,l,1) = abs(OPRT_coef_div_pl (0,3,l))
          tmp_pl(g,  4,l,1) = abs(OPRT_coef_div_pl (1,1,l))
          tmp_pl(g,  5,l,1) = abs(OPRT_coef_div_pl (1,2,l))
          tmp_pl(g,  6,l,1) = abs(OPRT_coef_div_pl (1,3,l))
          tmp_pl(g,  7,l,1) = abs(OPRT_coef_div_pl (2,1,l))
          tmp_pl(g,  8,l,1) = abs(OPRT_coef_div_pl (2,2,l))
          tmp_pl(g,  9,l,1) = abs(OPRT_coef_div_pl (2,3,l))
          tmp_pl(g, 10,l,1) = abs(OPRT_coef_div_pl (3,1,l))
          tmp_pl(g, 11,l,1) = abs(OPRT_coef_div_pl (3,2,l))
          tmp_pl(g, 12,l,1) = abs(OPRT_coef_div_pl (3,3,l))
          tmp_pl(g, 13,l,1) = abs(OPRT_coef_div_pl (4,1,l))
          tmp_pl(g, 14,l,1) = abs(OPRT_coef_div_pl (4,2,l))
          tmp_pl(g, 15,l,1) = abs(OPRT_coef_div_pl (4,3,l))
          tmp_pl(g, 16,l,1) = abs(OPRT_coef_div_pl (5,1,l))
          tmp_pl(g, 17,l,1) = abs(OPRT_coef_div_pl (5,2,l))
          tmp_pl(g, 18,l,1) = abs(OPRT_coef_div_pl (5,3,l))
          tmp_pl(g, 19,l,1) = abs(OPRT_coef_div_pl (1,1,l))
          tmp_pl(g, 20,l,1) = abs(OPRT_coef_div_pl (1,2,l))
          tmp_pl(g, 21,l,1) = abs(OPRT_coef_div_pl (1,3,l))
          tmp_pl(g, 22,l,1) = abs(OPRT_coef_rot_pl (0,1,l))
          tmp_pl(g, 23,l,1) = abs(OPRT_coef_rot_pl (0,2,l))
          tmp_pl(g, 24,l,1) = abs(OPRT_coef_rot_pl (0,3,l))
          tmp_pl(g, 25,l,1) = abs(OPRT_coef_rot_pl (1,1,l))
          tmp_pl(g, 26,l,1) = abs(OPRT_coef_rot_pl (1,2,l))
          tmp_pl(g, 27,l,1) = abs(OPRT_coef_rot_pl (1,3,l))
          tmp_pl(g, 28,l,1) = abs(OPRT_coef_rot_pl (2,1,l))
          tmp_pl(g, 29,l,1) = abs(OPRT_coef_rot_pl (2,2,l))
          tmp_pl(g, 30,l,1) = abs(OPRT_coef_rot_pl (2,3,l))
          tmp_pl(g, 31,l,1) = abs(OPRT_coef_rot_pl (3,1,l))
          tmp_pl(g, 32,l,1) = abs(OPRT_coef_rot_pl (3,2,l))
          tmp_pl(g, 33,l,1) = abs(OPRT_coef_rot_pl (3,3,l))
          tmp_pl(g, 34,l,1) = abs(OPRT_coef_rot_pl (4,1,l))
          tmp_pl(g, 35,l,1) = abs(OPRT_coef_rot_pl (4,2,l))
          tmp_pl(g, 36,l,1) = abs(OPRT_coef_rot_pl (4,3,l))
          tmp_pl(g, 37,l,1) = abs(OPRT_coef_rot_pl (5,1,l))
          tmp_pl(g, 38,l,1) = abs(OPRT_coef_rot_pl (5,2,l))
          tmp_pl(g, 39,l,1) = abs(OPRT_coef_rot_pl (5,3,l))
          tmp_pl(g, 40,l,1) = abs(OPRT_coef_rot_pl (1,1,l))
          tmp_pl(g, 41,l,1) = abs(OPRT_coef_rot_pl (1,2,l))
          tmp_pl(g, 42,l,1) = abs(OPRT_coef_rot_pl (1,3,l))
          tmp_pl(g, 43,l,1) = abs(OPRT_coef_grad_pl(0,1,l))
          tmp_pl(g, 44,l,1) = abs(OPRT_coef_grad_pl(0,2,l))
          tmp_pl(g, 45,l,1) = abs(OPRT_coef_grad_pl(0,3,l))
          tmp_pl(g, 46,l,1) = abs(OPRT_coef_grad_pl(1,1,l))
          tmp_pl(g, 47,l,1) = abs(OPRT_coef_grad_pl(1,2,l))
          tmp_pl(g, 48,l,1) = abs(OPRT_coef_grad_pl(1,3,l))
          tmp_pl(g, 49,l,1) = abs(OPRT_coef_grad_pl(2,1,l))
          tmp_pl(g, 50,l,1) = abs(OPRT_coef_grad_pl(2,2,l))
          tmp_pl(g, 51,l,1) = abs(OPRT_coef_grad_pl(2,3,l))
          tmp_pl(g, 52,l,1) = abs(OPRT_coef_grad_pl(3,1,l))
          tmp_pl(g, 53,l,1) = abs(OPRT_coef_grad_pl(3,2,l))
          tmp_pl(g, 54,l,1) = abs(OPRT_coef_grad_pl(3,3,l))
          tmp_pl(g, 55,l,1) = abs(OPRT_coef_grad_pl(4,1,l))
          tmp_pl(g, 56,l,1) = abs(OPRT_coef_grad_pl(4,2,l))
          tmp_pl(g, 57,l,1) = abs(OPRT_coef_grad_pl(4,3,l))
          tmp_pl(g, 58,l,1) = abs(OPRT_coef_grad_pl(5,1,l))
          tmp_pl(g, 59,l,1) = abs(OPRT_coef_grad_pl(5,2,l))
          tmp_pl(g, 60,l,1) = abs(OPRT_coef_grad_pl(5,3,l))
          tmp_pl(g, 61,l,1) = abs(OPRT_coef_grad_pl(1,1,l))
          tmp_pl(g, 62,l,1) = abs(OPRT_coef_grad_pl(1,2,l))
          tmp_pl(g, 63,l,1) = abs(OPRT_coef_grad_pl(1,3,l))
          tmp_pl(g, 64,l,1) = abs(OPRT_coef_lap_pl (0,  l))
          tmp_pl(g, 65,l,1) = abs(OPRT_coef_lap_pl (1,  l))
          tmp_pl(g, 66,l,1) = abs(OPRT_coef_lap_pl (2,  l))
          tmp_pl(g, 67,l,1) = abs(OPRT_coef_lap_pl (3,  l))
          tmp_pl(g, 68,l,1) = abs(OPRT_coef_lap_pl (4,  l))
          tmp_pl(g, 69,l,1) = abs(OPRT_coef_lap_pl (5,  l))
          tmp_pl(g, 70,l,1) = abs(OPRT_coef_lap_pl (1,  l))
          tmp_pl(g, 71,l,1) = abs(OPRT_coef_intp_pl(g,1,1,l))
          tmp_pl(g, 72,l,1) = abs(OPRT_coef_intp_pl(g,1,2,l))
          tmp_pl(g, 73,l,1) = abs(OPRT_coef_intp_pl(g,1,3,l))
          tmp_pl(g, 74,l,1) = abs(OPRT_coef_intp_pl(g,2,1,l))
          tmp_pl(g, 75,l,1) = abs(OPRT_coef_intp_pl(g,2,2,l))
          tmp_pl(g, 76,l,1) = abs(OPRT_coef_intp_pl(g,2,3,l))
          tmp_pl(g, 77,l,1) = abs(OPRT_coef_intp_pl(g,3,1,l))
          tmp_pl(g, 78,l,1) = abs(OPRT_coef_intp_pl(g,3,2,l))
          tmp_pl(g, 79,l,1) = abs(OPRT_coef_intp_pl(g,3,3,l))
          tmp_pl(g, 80,l,1) = abs(OPRT_coef_intp_pl(g,1,1,l))
          tmp_pl(g, 81,l,1) = abs(OPRT_coef_intp_pl(g,1,2,l))
          tmp_pl(g, 82,l,1) = abs(OPRT_coef_intp_pl(g,1,3,l))
          tmp_pl(g, 83,l,1) = abs(OPRT_coef_intp_pl(g,2,1,l))
          tmp_pl(g, 84,l,1) = abs(OPRT_coef_intp_pl(g,2,2,l))
          tmp_pl(g, 85,l,1) = abs(OPRT_coef_intp_pl(g,2,3,l))
          tmp_pl(g, 86,l,1) = abs(OPRT_coef_intp_pl(g,3,1,l))
          tmp_pl(g, 87,l,1) = abs(OPRT_coef_intp_pl(g,3,2,l))
          tmp_pl(g, 88,l,1) = abs(OPRT_coef_intp_pl(g,3,3,l))
          tmp_pl(g, 89,l,1) = abs(OPRT_coef_diff_pl(1,1,l))
          tmp_pl(g, 90,l,1) = abs(OPRT_coef_diff_pl(1,2,l))
          tmp_pl(g, 91,l,1) = abs(OPRT_coef_diff_pl(1,3,l))
          tmp_pl(g, 92,l,1) = abs(OPRT_coef_diff_pl(2,1,l))
          tmp_pl(g, 93,l,1) = abs(OPRT_coef_diff_pl(2,2,l))
          tmp_pl(g, 94,l,1) = abs(OPRT_coef_diff_pl(2,3,l))
          tmp_pl(g, 95,l,1) = abs(OPRT_coef_diff_pl(3,1,l))
          tmp_pl(g, 96,l,1) = abs(OPRT_coef_diff_pl(3,2,l))
          tmp_pl(g, 97,l,1) = abs(OPRT_coef_diff_pl(3,3,l))
          tmp_pl(g, 98,l,1) = abs(OPRT_coef_diff_pl(4,1,l))
          tmp_pl(g, 99,l,1) = abs(OPRT_coef_diff_pl(4,2,l))
          tmp_pl(g,100,l,1) = abs(OPRT_coef_diff_pl(4,3,l))
          tmp_pl(g,101,l,1) = abs(OPRT_coef_diff_pl(5,1,l))
          tmp_pl(g,102,l,1) = abs(OPRT_coef_diff_pl(5,2,l))
          tmp_pl(g,103,l,1) = abs(OPRT_coef_diff_pl(5,3,l))
          tmp_pl(g,104,l,1) = abs(OPRT_coef_diff_pl(1,1,l))
          tmp_pl(g,105,l,1) = abs(OPRT_coef_diff_pl(1,2,l))
          tmp_pl(g,106,l,1) = abs(OPRT_coef_diff_pl(1,3,l))
       enddo
       enddo
    endif

    call COMM_data_transfer( tmp, tmp_pl )

    if ( OPRT_io_mode == 'ADVANCED' ) then

       call FIO_output( tmp(:,:,:,1), basename, desc, "",               & ! [IN]
                        "oprtcoef", "oprt coef", "",                    & ! [IN]
                        "", dtype, "LAYERNM", 1, 106, 1, 0.0_DP, 0.0_DP ) ! [IN]

    else
       if( IO_L ) write(IO_FID_LOG,*) 'Invalid io_mode!'
       call PRC_MPIstop
    endif

    return
  end subroutine OPRT_output_coef

end module mod_oprt
