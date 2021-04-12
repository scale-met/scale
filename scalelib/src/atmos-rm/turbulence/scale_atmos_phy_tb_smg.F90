!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!          Smagolinsky-type
!!
!! @author Team SCALE
!!
!! @par Reference
!!  - Brown et al., 1994:
!!    Large-eddy simulaition of stable atmospheric boundary layers with a revised stochastic subgrid model.
!!    Roy. Meteor. Soc., 120, 1485-1512
!!  - Scotti et al., 1993:
!!    Generalized Smagorinsky model for anisotropic grids.
!!    Phys. Fluids A, 5, 2306-2308
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_tb_smg
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer

#if defined DEBUG || defined QUICKDEBUG
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
  public :: ATMOS_PHY_TB_smg_setup
  public :: ATMOS_PHY_TB_smg_finalize
  public :: ATMOS_PHY_TB_smg

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
  real(RP), private, parameter   :: OneOverThree  = 1.0_RP / 3.0_RP
  real(RP), private, parameter   :: twoOverThree  = 2.0_RP / 3.0_RP
  real(RP), private, parameter   :: FourOverThree = 4.0_RP / 3.0_RP

  real(RP), private              :: Cs            = 0.13_RP ! Smagorinsky constant (Scotti et al. 1993)
  real(RP), private, parameter   :: PrN           = 0.7_RP  ! Prandtl number in neutral conditions
  real(RP), private, parameter   :: RiC           = 0.25_RP ! critical Richardson number
  real(RP), private, parameter   :: FmC           = 16.0_RP ! fum = sqrt(1 - c*Ri)
  real(RP), private, parameter   :: FhB           = 40.0_RP ! fuh = sqrt(1 - b*Ri)/PrN
  real(RP), private              :: RPrN                    ! 1 / PrN
  real(RP), private              :: RRiC                    ! 1 / RiC
  real(RP), private              :: PrNovRiC                ! PrN / RiC

  ! for backscatter
  real(RP), private, parameter   :: CB   = 1.4_RP
  real(RP), private, parameter   :: CBt  = 0.45_RP
  real(RP), private, parameter   :: aN   = 0.47958315233127197_RP ! a_N = sqrt(0.23)
  real(RP), private, parameter   :: atN4 = 0.09_RP                 ! a_{\theta N}^4 = 0.3**2
  real(RP), private, parameter   :: C1o  = aN**3
  real(RP), private, parameter   :: D1o  = PrN * atN4 / aN

  real(RP), private, allocatable :: lambda0(:,:,:)
  real(RP), private, allocatable :: lambda (:,:,:)

  real(RP), private              :: ATMOS_PHY_TB_SMG_NU_MAX      = 10000.0_RP
  logical,  private              :: ATMOS_PHY_TB_SMG_backscatter = .false.
  logical,  private              :: ATMOS_PHY_TB_SMG_bottom      = .true.
  logical,  private              :: ATMOS_PHY_TB_SMG_implicit    = .false.
  logical,  private              :: ATMOS_PHY_TB_SMG_horizontal  = .false.

  real(RP), private              :: tke_fact

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_TB_smg_setup( &
       FZ, CZ, CDX, CDY, MAPF, &
       horizontal )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       KARMAN  => CONST_KARMAN
    implicit none

    real(RP), intent(in) :: FZ  (0:KA,IA,JA)
    real(RP), intent(in) :: CZ  (  KA,IA,JA)
    real(RP), intent(in) :: CDX (IA)
    real(RP), intent(in) :: CDY (JA)
    real(RP), intent(in) :: MAPF(IA,JA,2)

    logical,  intent(in), optional :: horizontal

    real(RP) :: ATMOS_PHY_TB_SMG_Cs
    real(RP) :: ATMOS_PHY_TB_SMG_filter_fact
    logical  :: ATMOS_PHY_TB_SMG_consistent_tke

    namelist / PARAM_ATMOS_PHY_TB_SMG / &
       ATMOS_PHY_TB_SMG_Cs,             &
       ATMOS_PHY_TB_SMG_NU_MAX,         &
       ATMOS_PHY_TB_SMG_filter_fact,    &
       ATMOS_PHY_TB_SMG_bottom,         &
       ATMOS_PHY_TB_SMG_backscatter,    &
       ATMOS_PHY_TB_SMG_implicit,       &
       ATMOS_PHY_TB_SMG_consistent_tke, &
       ATMOS_PHY_TB_SMG_horizontal

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_smg_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_smg_setup",*) 'Smagorinsky-type Eddy Viscocity Model'

    ATMOS_PHY_TB_SMG_Cs             = Cs
    ATMOS_PHY_TB_SMG_filter_fact    = 2.0_RP
    ATMOS_PHY_TB_SMG_consistent_tke = .true.

    if ( present(horizontal) ) ATMOS_PHY_TB_SMG_horizontal = horizontal

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_SMG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_smg_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_smg_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_SMG. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_SMG)

    Cs = ATMOS_PHY_TB_SMG_Cs

    RPrN     = 1.0_RP / PrN
    RRiC     = 1.0_RP / RiC
    PrNovRiC = ( 1.0_RP - PrN ) * RRiC

    allocate( lambda0(KA,IA,JA) )
    allocate( lambda (KA,IA,JA) )

#ifdef DEBUG
    lambda0(:,:,:) = UNDEF
    lambda (:,:,:) = UNDEF
#endif
    if ( ATMOS_PHY_TB_SMG_horizontal ) then
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          lambda0(k,i,j) = Cs * sqrt( CDX(i) * CDY(j) / ( MAPF(i,j,1) * MAPF(i,j,2) ) )
          lambda (k,i,j) = lambda0(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ATMOS_PHY_TB_SMG_consistent_tke = .false.
       ATMOS_PHY_TB_SMG_implicit       = .false. ! flux in the z-direction is not necessary
       ATMOS_PHY_TB_SMG_backscatter    = .false.
    else
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          lambda0(k,i,j) = Cs * mixlen( FZ(k,i,j) - FZ(k-1,i,j),     &
                                        CDX(i) / MAPF(i,j,1),        &
                                        CDY(j) / MAPF(i,j,2),        &
                                        ATMOS_PHY_TB_SMG_filter_fact )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       if ( ATMOS_PHY_TB_SMG_bottom ) then
          !$omp parallel do collapse(2)
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             lambda(k,i,j) = sqrt( 1.0_RP / ( 1.0_RP / lambda0(k,i,j)**2 + 1.0_RP / ( KARMAN*(CZ(k,i,j)-FZ(KS-1,i,j)) )**2 ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       else
          !$omp parallel do collapse(2)
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             lambda(k,i,j) = lambda0(k,i,j)
          enddo
          enddo
          enddo
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       end if
    end if

    ! flag for isotropic stress tensor
    if ( ATMOS_PHY_TB_SMG_consistent_tke ) then
       tke_fact = 1.0_RP
    else
       tke_fact = 0.0_RP ! neglect
    end if

    return
  end subroutine ATMOS_PHY_TB_smg_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_TB_smg_finalize

    deallocate( lambda0 )
    deallocate( lambda  )

    return
  end subroutine ATMOS_PHY_TB_smg_finalize

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_smg( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, &
       qflx_sgs_rhot, qflx_sgs_rhoq,                &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, RHOQ_t,      &
       Nu, Ri, Pr,                                  &
       MOMZ, MOMX, MOMY, POTT, DENS, QTRC, N2,      &
       FZ, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY,     &
       GSQRT, J13G, J23G, J33G, MAPF, dt            )
    use scale_precision
    use scale_atmos_grid_cartesC_index
    use scale_tracer
    use scale_const, only: &
       EPS  => CONST_EPS, &
       GRAV => CONST_GRAV
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_phy_tb_common, only: &
       calc_strain_tensor => ATMOS_PHY_TB_calc_strain_tensor, &
       diffusion_solver   => ATMOS_PHY_TB_diffusion_solver,   &
       calc_tend_momz     => ATMOS_PHY_TB_calc_tend_momz,     &
       calc_tend_momx     => ATMOS_PHY_TB_calc_tend_momx,     &
       calc_tend_momy     => ATMOS_PHY_TB_calc_tend_momy,     &
       calc_flux_phi      => ATMOS_PHY_TB_calc_flux_phi
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_random, only: &
       RANDOM_normal
    implicit none

    ! SGS flux
    real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

    real(RP), intent(out) :: MOMZ_t       (KA,IA,JA)
    real(RP), intent(out) :: MOMX_t       (KA,IA,JA)
    real(RP), intent(out) :: MOMY_t       (KA,IA,JA)
    real(RP), intent(out) :: RHOT_t       (KA,IA,JA)
    real(RP), intent(out) :: RHOQ_t       (KA,IA,JA,QA) ! tendency of rho * QTRC

    real(RP), intent(out) :: nu           (KA,IA,JA)    ! eddy viscosity (center)
    real(RP), intent(out) :: Ri           (KA,IA,JA)    ! Richardson number
    real(RP), intent(out) :: Pr           (KA,IA,JA)    ! Prantle number

    real(RP), intent(in)  :: MOMZ         (KA,IA,JA)
    real(RP), intent(in)  :: MOMX         (KA,IA,JA)
    real(RP), intent(in)  :: MOMY         (KA,IA,JA)
    real(RP), intent(in)  :: POTT         (KA,IA,JA)
    real(RP), intent(in)  :: DENS         (KA,IA,JA)
    real(RP), intent(in)  :: QTRC         (KA,IA,JA,QA)
    real(RP), intent(in)  :: N2           (KA,IA,JA)

    real(RP), intent(in)  :: FZ           (0:KA,IA,JA)
    real(RP), intent(in)  :: FDZ          (KA-1)
    real(RP), intent(in)  :: RCDZ         (KA)
    real(RP), intent(in)  :: RFDZ         (KA-1)
    real(RP), intent(in)  :: CDX          (IA)
    real(RP), intent(in)  :: FDX          (IA-1)
    real(RP), intent(in)  :: CDY          (JA)
    real(RP), intent(in)  :: FDY          (JA-1)

    real(RP), intent(in)  :: GSQRT         (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G          (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G          (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                       !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF(IA,JA,2,4)            !< map factor
    real(DP), intent(in)  :: dt

    ! tke
    real(RP) :: TKE(KA,IA,JA)

    ! deformation rate tensor
    real(RP) :: S33_C(KA,IA,JA) ! (cell center)
    real(RP) :: S11_C(KA,IA,JA)
    real(RP) :: S22_C(KA,IA,JA)
    real(RP) :: S31_C(KA,IA,JA)
    real(RP) :: S12_C(KA,IA,JA)
    real(RP) :: S23_C(KA,IA,JA)
    real(RP) :: S12_Z(KA,IA,JA) ! (z edge or x-y plane)
    real(RP) :: S23_X(KA,IA,JA) ! (x edge or y-z plane)
    real(RP) :: S31_Y(KA,IA,JA) ! (y edge or z-x plane)
    real(RP) :: S2   (KA,IA,JA) ! |S|^2

    real(RP) :: Kh(KA,IA,JA) ! eddy diffusion
    real(RP) :: fm(KA)
    real(RP) :: e(KA)
    real(RP) :: et
    real(RP) :: lambda_r(KA)
    real(RP) :: Rf
    real(RP) :: C1(KA)
    real(RP) :: C2
    real(RP) :: D2

    ! implicit scheme
    real(RP) :: TEND(KA,IA,JA)
    real(RP) :: a   (KA,IA,JA)
    real(RP) :: b   (KA,IA,JA)
    real(RP) :: c   (KA,IA,JA)
    real(RP) :: d   (KA)
    real(RP) :: ap

    ! backscatter
    real(RP) :: random   (KA,IA,JA)
    real(RP) :: random_mz(KA,IA,JA)
    real(RP) :: random_mx(KA,IA,JA)
    real(RP) :: random_my(KA,IA,JA)
    real(RP) :: random_qz(KA,IA,JA)
    real(RP) :: random_qx(KA,IA,JA)
    real(RP) :: random_qy(KA,IA,JA)
    real(RP) :: dd       (KA,IA,JA)
    real(RP) :: leOvleo5
    real(RP) :: dz, dx, dy
    real(RP) :: fact
    real(RP) :: flxz(KA)

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / turbulence / Smagorinsky'

#ifdef DEBUG
    qflx_sgs_momz(:,:,:,:)   = UNDEF
    qflx_sgs_momx(:,:,:,:)   = UNDEF
    qflx_sgs_momy(:,:,:,:)   = UNDEF
    qflx_sgs_rhot(:,:,:,:)   = UNDEF
    qflx_sgs_rhoq(:,:,:,:,:) = UNDEF

    nu           (:,:,:)     = UNDEF
    tke          (:,:,:)     = UNDEF
    Pr           (:,:,:)     = UNDEF
    Ri           (:,:,:)     = UNDEF
    Kh           (:,:,:)     = UNDEF
#endif

#ifdef QUICKDEBUG
    qflx_sgs_momz(KS:KE,   1:IS-1,    :    ,:) = UNDEF
    qflx_sgs_momz(KS:KE,IE+1:IA  ,    :    ,:) = UNDEF
    qflx_sgs_momz(KS:KE,    :    ,   1:JS-1,:) = UNDEF
    qflx_sgs_momz(KS:KE,    :    ,JE+1:JA  ,:) = UNDEF
    qflx_sgs_momx(KS:KE,   1:IS-1,    :    ,:) = UNDEF
    qflx_sgs_momx(KS:KE,IE+1:IA  ,    :    ,:) = UNDEF
    qflx_sgs_momx(KS:KE,    :    ,   1:JS-1,:) = UNDEF
    qflx_sgs_momx(KS:KE,    :    ,JE+1:JA  ,:) = UNDEF
    qflx_sgs_momy(KS:KE,   1:IS-1,    :    ,:) = UNDEF
    qflx_sgs_momy(KS:KE,IE+1:IA  ,    :    ,:) = UNDEF
    qflx_sgs_momy(KS:KE,    :    ,   1:JS-1,:) = UNDEF
    qflx_sgs_momy(KS:KE,    :    ,JE+1:JA  ,:) = UNDEF
#endif


    !##### Start Upadate #####

    call calc_strain_tensor( &
         S33_C, S11_C, S22_C,          & ! (out)
         S31_C, S12_C, S23_C,          & ! (out)
         S12_Z, S23_X, S31_Y,          & ! (out)
         S2                 ,          & ! (out)
         DENS, MOMZ, MOMX, MOMY,       & ! (in)
         GSQRT, J13G, J23G, J33G, MAPF ) ! (in)

    if ( ATMOS_PHY_TB_SMG_backscatter ) then

       call RANDOM_normal( random(:,:,:) )
       ! 1:2:1 filter
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          do k = KS+1, KE-1
             random_mz(k,i,j) = ( random(k,i,j) * 2.0_RP &
                                + random(k+1,i,j) + random(k-1,i,j) &
                                + random(k,i+1,j) + random(k,i-1,j) &
                                + random(k,i,j+1) + random(k,i,j-1) ) / 8.0_RP
          end do
          random_mz(KS,i,j) = ( random(KS,i,j) * 2.0_RP &
                              + random(KS+1,i,j) &
                              + random(KS,i+1,j) + random(KS,i-1,j) &
                              + random(KS,i,j+1) + random(KS,i,j-1) ) / 7.0_RP
          random_mz(KE,i,j) = ( random(KE,i,j) * 2.0_RP &
                              + random(KE-1,i,j) &
                              + random(KE,i+1,j) + random(KE,i-1,j) &
                              + random(KE,i,j+1) + random(KE,i,j-1) ) / 7.0_RP
       end do
       end do

       call RANDOM_normal( random(:,:,:) )
       ! 1:2:1 filter
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          do k = KS+1, KE-1
             random_mx(k,i,j) = ( random(k,i,j) * 2.0_RP &
                                + random(k+1,i,j) + random(k-1,i,j) &
                                + random(k,i+1,j) + random(k,i-1,j) &
                                + random(k,i,j+1) + random(k,i,j-1) ) / 8.0_RP
          end do
          random_mx(KS,i,j) = ( random(KS,i,j) * 2.0_RP &
                              + random(KS+1,i,j) &
                              + random(KS,i+1,j) + random(KS,i-1,j) &
                              + random(KS,i,j+1) + random(KS,i,j-1) ) / 7.0_RP
          random_mx(KE,i,j) = ( random(KE,i,j) * 2.0_RP &
                              + random(KE-1,i,j) &
                              + random(KE,i+1,j) + random(KE,i-1,j) &
                              + random(KE,i,j+1) + random(KE,i,j-1) ) / 7.0_RP
       end do
       end do

       call RANDOM_normal( random(:,:,:) )
       ! 1:2:1 filter
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          do k = KS+1, KE-1
             random_my(k,i,j) = ( random(k,i,j) * 2.0_RP &
                                + random(k+1,i,j) + random(k-1,i,j) &
                                + random(k,i+1,j) + random(k,i-1,j) &
                                + random(k,i,j+1) + random(k,i,j-1) ) / 8.0_RP
          end do
          random_my(KS,i,j) = ( random(KS,i,j) * 2.0_RP &
                              + random(KS+1,i,j) &
                              + random(KS,i+1,j) + random(KS,i-1,j) &
                              + random(KS,i,j+1) + random(KS,i,j-1) ) / 7.0_RP
          random_my(KE,i,j) = ( random(KE,i,j) * 2.0_RP &
                              + random(KE-1,i,j) &
                              + random(KE,i+1,j) + random(KE,i-1,j) &
                              + random(KE,i,j+1) + random(KE,i,j-1) ) / 7.0_RP
       end do
       end do

       call RANDOM_normal( random(:,:,:) )
       ! 1:2:1 filter
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          do k = KS+1, KE-1
             random_qz(k,i,j) = ( random(k,i,j) * 2.0_RP &
                                + random(k+1,i,j) + random(k-1,i,j) &
                                + random(k,i+1,j) + random(k,i-1,j) &
                                + random(k,i,j+1) + random(k,i,j-1) ) / 8.0_RP
          end do
          random_qz(KS,i,j) = ( random(KS,i,j) * 2.0_RP &
                              + random(KS+1,i,j) &
                              + random(KS,i+1,j) + random(KS,i-1,j) &
                              + random(KS,i,j+1) + random(KS,i,j-1) ) / 7.0_RP
          random_qz(KE,i,j) = ( random(KE,i,j) * 2.0_RP &
                              + random(KE-1,i,j) &
                              + random(KE,i+1,j) + random(KE,i-1,j) &
                              + random(KE,i,j+1) + random(KE,i,j-1) ) / 7.0_RP
       end do
       end do

       call RANDOM_normal( random(:,:,:) )
       ! 1:2:1 filter
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          do k = KS+1, KE-1
             random_qx(k,i,j) = ( random(k,i,j) * 2.0_RP &
                                + random(k+1,i,j) + random(k-1,i,j) &
                                + random(k,i+1,j) + random(k,i-1,j) &
                                + random(k,i,j+1) + random(k,i,j-1) ) / 8.0_RP
          end do
          random_qx(KS,i,j) = ( random(KS,i,j) * 2.0_RP &
                              + random(KS+1,i,j) &
                              + random(KS,i+1,j) + random(KS,i-1,j) &
                              + random(KS,i,j+1) + random(KS,i,j-1) ) / 7.0_RP
          random_qx(KE,i,j) = ( random(KE,i,j) * 2.0_RP &
                              + random(KE-1,i,j) &
                              + random(KE,i+1,j) + random(KE,i-1,j) &
                              + random(KE,i,j+1) + random(KE,i,j-1) ) / 7.0_RP
       end do
       end do

       call RANDOM_normal( random(:,:,:) )
       ! 1:2:1 filter
       !$omp parallel do collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
          do k = KS+1, KE-1
             random_qy(k,i,j) = ( random(k,i,j) * 2.0_RP &
                                + random(k+1,i,j) + random(k-1,i,j) &
                                + random(k,i+1,j) + random(k,i-1,j) &
                                + random(k,i,j+1) + random(k,i,j-1) ) / 8.0_RP
          end do
          random_qy(KS,i,j) = ( random(KS,i,j) * 2.0_RP &
                              + random(KS+1,i,j) &
                              + random(KS,i+1,j) + random(KS,i-1,j) &
                              + random(KS,i,j+1) + random(KS,i,j-1) ) / 7.0_RP
          random_qy(KE,i,j) = ( random(KE,i,j) * 2.0_RP &
                              + random(KE-1,i,j) &
                              + random(KE,i+1,j) + random(KE,i-1,j) &
                              + random(KE,i,j+1) + random(KE,i,j-1) ) / 7.0_RP
       end do
       end do


    end if


    !$omp parallel do collapse(2) &
    !$omp private(fm,Rf,lambda_r,leOvleo5,C1,C2,D2,e,et,fact,dz,dx,dy)
    do j = JS-1, JE+1
    do i = IS-1, IE+1

       ! Ri = N^2 / |S|^2, N^2 = g / theta * dtheta/dz
       do k = KS, KE
          Ri(k,i,j) = N2(k,i,j) / max(S2(k,i,j), EPS)
       enddo

       ! Nu
       ! Pr = Nu / Kh = fm / fh
       do k = KS, KE
          if ( Ri(k,i,j) < 0.0_RP ) then ! stable
             fm(k) = sqrt( 1.0_RP - FmC * Ri(k,i,j) )
             nu(k,i,j) = lambda(k,i,j)**2 * sqrt( S2(k,i,j) ) * fm(k)
             Pr(k,i,j) = fm(k) / sqrt( 1.0_RP - FhB*Ri(k,i,j) ) * PrN
          else if ( Ri(k,i,j) < RiC ) then ! weakly stable
             fm(k) = ( 1.0_RP - Ri(k,i,j)*RRiC )**4
             nu(k,i,j) = lambda(k,i,j)**2 * sqrt( S2(k,i,j) ) * fm(k)
             Pr(k,i,j) = PrN / ( 1.0_RP - PrNovRiC * Ri(k,i,j) )
          else ! strongly stable
             fm(k) = 0.0_RP
             nu(k,i,j) = 0.0_RP
             Kh(k,i,j) = 0.0_RP
             Pr(k,i,j) = 1.0_RP
          endif

          if ( Ri(k,i,j) < RiC ) then
             Kh(k,i,j) = max( min( nu(k,i,j) / Pr(k,i,j), ATMOS_PHY_TB_SMG_NU_MAX ), EPS )
             nu(k,i,j) = max( min( nu(k,i,j), ATMOS_PHY_TB_SMG_NU_MAX ), EPS )
             Pr(k,i,j) = nu(k,i,j) / Kh(k,i,j)
             Rf = Ri(k,i,j) / Pr(k,i,j)
             lambda_r(k) = lambda(k,i,j) * sqrt( fm(k) / sqrt( 1.0_RP - Rf ) )
          else
             lambda_r(k) = 0.0_RP
          end if

       enddo

       if ( ATMOS_PHY_TB_SMG_backscatter ) then

          do k = KS, KE
             lambda_r(k) = min( 1.8_RP * lambda0(k,i,j), lambda_r(k) )
             leOvleo5 = ( lambda_r(k) / lambda0(k,i,j) )**5
             C2 = CB  * leOvleo5 / ( 1.0_RP + CB * leOvleo5 )
             D2 = CBt  * leOvleo5 / ( 1.0_RP + CBt * leOvleo5 )
             C1(k) = C1o / sqrt( 1.0_RP - C2 )
             ! D1 = D1o / sqrt( 1.0_RP - C2 )
             e(k) = nu(k,i,j)**3 * ( 1.0_RP - C2 ) / ( lambda_r(k)**4 + EPS )
             et = Kh(k,i,j) * ( 1.0_RP - D2 ) ! epsilon_t / D^2

             dz = FZ(k,i,j) - FZ(k-1,i,j)
             dx = CDX(i) / MAPF(i,j,1,I_XY)
             dy = CDY(j) / MAPF(i,j,2,I_XY)

             fact = sqrt( CB * leovleo5 * e(k) / ( 3.0_RP * dt ) ) * DENS(k,i,j)
             random_mz(k,i,j) = random_mz(k,i,j) * fact * dx * dy / sqrt( dx**2 + dy**2 )
             random_mx(k,i,j) = random_mx(k,i,j) * fact * dy * dz / sqrt( dy**2 + dz**2 )
             random_my(k,i,j) = random_my(k,i,j) * fact * dz * dx / sqrt( dz**2 + dx**2 )

             fact = sqrt( CBt * leovleo5 * et / ( 3.0_RP * dt ) ) * DENS(k,i,j)
             random_qz(k,i,j) = random_qz(k,i,j) * fact * dz
             random_qx(k,i,j) = random_qx(k,i,j) * fact * dx
             random_qy(k,i,j) = random_qy(k,i,j) * fact * dy
          end do

       else

          do k = KS, KE
             e(k) = nu(k,i,j)**3 / ( lambda_r(k)**4 + EPS )
             C1(k) = C1o
          end do

       end if

       ! TKE
       do k = KS, KE
          TKE(k,i,j) = ( e(k) * lambda_r(k) / C1(k) )**TwoOverThree
       enddo

    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif



    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (z) #####

       !$omp parallel private(i,j,k)

       ! (cell center)
       if ( ATMOS_PHY_TB_SMG_horizontal ) then
          !$omp workshare
          qflx_sgs_momz(:,:,:,ZDIR) = 0.0_RP
          !$omp end workshare nowait
       else
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS+1, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, DENS(k,i,j) )
             call CHECK( __LINE__, nu(k,i,j) )
             call CHECK( __LINE__, S33_C(k,i,j) )
             call CHECK( __LINE__, S11_C(k,i,j) )
             call CHECK( __LINE__, S22_C(k,i,j) )
             call CHECK( __LINE__, tke(k,i,j) )
#endif
             qflx_sgs_momz(k,i,j,ZDIR) = DENS(k,i,j) * ( &
                  - 2.0_RP * nu(k,i,j) &
                  * ( S33_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) &
                  + twoOverThree * tke(k,i,j) * tke_fact )
          enddo
          enddo
          enddo
          !$omp end do nowait
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
          !$omp do
          do j = JJS, JJE
          do i = IIS, IIE
             ! momentum will not be conserved
             qflx_sgs_momz(KS,i,j,ZDIR) = DENS(KS,i,j) * twoOverThree * tke(KS,i,j) * tke_fact
             qflx_sgs_momz(KE,i,j,ZDIR) = DENS(KE,i,j) * twoOverThree * tke(KE,i,j) * tke_fact
             ! anti-isotropic stress is calculated by the surface scheme
          enddo
          enddo
          !$omp end do nowait
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       end if
       ! (y edge)
       !$omp do collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i+1,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, nu(k+1,i+1,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,XDIR) = - 0.125_RP & ! 2.0 / 4 / 4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
               * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i+1,j)+nu  (k+1,i+1,j)) &
               * S31_Y(k,i,j)
       enddo
       enddo
       enddo
       !$omp end do nowait
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       !$omp do collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j+1) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, nu(k+1,i,j+1) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,YDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
               * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i,j+1)+nu  (k+1,i,j+1) ) &
               * S23_X(k,i,j)
       enddo
       enddo
       enddo
       !$omp end do nowait
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !$omp end parallel

       if ( ATMOS_PHY_TB_SMG_implicit ) then

          call calc_tend_MOMZ( TEND, & ! (out)
                               qflx_sgs_momz, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          !$omp parallel do collapse(2) &
          !$omp private (ap,d)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - FourOverThree * dt &
                  * DENS(KS+1,i,j)*Nu(KS+1,i,j) &
                  * RCDZ(KS+1) / GSQRT(KS+1,i,j,I_XYZ)
             a(KS,i,j) = ap * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS+1,i,j) )
             do k = KS+1, KE-2
                c(k,i,j) = ap * RFDZ(k+1) / GSQRT(k+1,i,j,I_XYW)
                ap = - FourOverThree * dt &
                     * DENS(k+1,i,j)*Nu(k+1,i,j) &
                     * RCDZ(k+1) / GSQRT(k+1,i,j,I_XYZ)
                a(k,i,j) = ap * RFDZ(k) / GSQRT(k,i,j,I_XYW)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k+1,i,j) )
             end do
             a(KE-1,i,j) = 0.0_RP
             c(KE-1,i,j) = ap * RFDZ(KE) / GSQRT(KE,i,j,I_XYW)
             b(KE-1,i,j) = - c(KE-1,i,j) + 0.5_RP * ( DENS(KE-1,i,j)+DENS(KE,i,j) )

             do k = KS, KE-1
                d(k) = TEND(k,i,j)
             end do

             call diffusion_solver( &
                  TEND(:,i,j),                     & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d, & ! (in)
                  KE-1                             ) ! (in)

             do k = KS+1, KE-1
                qflx_sgs_momz(k,i,j,ZDIR) = qflx_sgs_momz(k,i,j,ZDIR) &
                     - FourOverThree * DENS(k,i,j) * Nu(k,i,j) * dt &
                     * ( TEND(k,i,j) - TEND(k-1,i,j) ) * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
             end do

          end do
          end do

       end if

       !##### momentum equation (x) #####

       !$omp parallel private(i,j,k)

       ! (y edge)
       if ( ATMOS_PHY_TB_SMG_horizontal ) then
          !$omp workshare
          qflx_sgs_momx(:,:,:,ZDIR) = 0.0_RP
          !$omp end workshare nowait
       else
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, DENS(k,i,j) )
             call CHECK( __LINE__, DENS(k,i+1,j) )
             call CHECK( __LINE__, DENS(k+1,i,j) )
             call CHECK( __LINE__, DENS(k+1,i+1,j) )
             call CHECK( __LINE__, nu(k,i,j) )
             call CHECK( __LINE__, nu(k,i+1,j) )
             call CHECK( __LINE__, nu(k+1,i,j) )
             call CHECK( __LINE__, nu(k+1,i+1,j) )
             call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
             qflx_sgs_momx(k,i,j,ZDIR) = - 0.125_RP & ! 2/4/4
                  * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
                  * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i+1,j)+nu  (k+1,i+1,j) ) &
                  * S31_Y(k,i,j)
          enddo
          enddo
          enddo
          !$omp end do nowait
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
          !$omp do
          do j = JJS, JJE
          do i = IIS, IIE
             qflx_sgs_momx(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
             qflx_sgs_momx(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
          enddo
          enddo
          !$omp end do nowait
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       end if

       ! (cell center)
       !$omp do collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, TKE(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,XDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu(k,i,j) &
               * ( S11_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) &
             + twoOverThree * TKE(k,i,j) * tke_fact )
       enddo
       enddo
       enddo
       !$omp end do nowait
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z edge)
       !$omp do collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i+1,j+1) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, nu(k,i+1,j+1) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,YDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ( nu  (k,i,j)+nu  (k,i+1,j)+nu  (k,i,j+1)+nu  (k,i+1,j+1) ) &
               * S12_Z(k,i,j)
       enddo
       enddo
       enddo
       !$omp end do nowait
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !$omp end parallel

       if ( ATMOS_PHY_TB_SMG_implicit ) then
          call calc_tend_MOMX( TEND, & ! (out)
                               qflx_sgs_momx, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          !$omp parallel do collapse(2) &
          !$omp private(ap,d)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS  ,i  ,j)*Nu(KS  ,i  ,j) &
                                   + DENS(KS+1,i  ,j)*Nu(KS+1,i  ,j) &
                                   + DENS(KS  ,i+1,j)*Nu(KS  ,i+1,j) &
                                   + DENS(KS+1,i+1,j)*Nu(KS+1,i+1,j) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_UYW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_UYZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS,i+1,j) )
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_UYZ)
                ap = - dt * 0.25_RP * ( DENS(k  ,i  ,j)*Nu(k  ,i  ,j) &
                                      + DENS(k+1,i  ,j)*Nu(k+1,i  ,j) &
                                      + DENS(k  ,i+1,j)*Nu(k  ,i+1,j) &
                                      + DENS(k+1,i+1,j)*Nu(k+1,i+1,j) ) &
                                    * RFDZ(k) / GSQRT(k,i,j,I_UYW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_UYZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) )
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_UYZ)
             b(KE,i,j) = - c(KE,i,j) + 0.5_RP * ( DENS(KE,i,j)+DENS(KE,i+1,j) )

             do k = KS, KE
                d(k) = TEND(k,i,j)
             end do

             call diffusion_solver( &
                  TEND(:,i,j),                     & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d, & ! (in)
                  KE                               ) ! (in)

             do k = KS, KE-1
                qflx_sgs_momx(k,i,j,ZDIR) = qflx_sgs_momx(k,i,j,ZDIR) &
                     - 0.25_RP * ( DENS(k  ,i  ,j)*Nu(k  ,i  ,j) &
                                 + DENS(k+1,i  ,j)*Nu(k+1,i  ,j) &
                                 + DENS(k  ,i+1,j)*Nu(k  ,i+1,j) &
                                 + DENS(k+1,i+1,j)*Nu(k+1,i+1,j) ) &
                     * dt * ( TEND(k+1,i,j) - TEND(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_UYW)
             end do

          end do
          end do

       end if

       !##### momentum equation (y) #####
       ! (x edge)

       !$omp parallel private(i,j,k)

       if ( ATMOS_PHY_TB_SMG_horizontal ) then
          !$omp workshare
          qflx_sgs_momy(:,:,:,ZDIR) = 0.0_RP
          !$omp end workshare nowait
       else
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, DENS(k,i,j) )
             call CHECK( __LINE__, DENS(k,i,j+1) )
             call CHECK( __LINE__, DENS(k+1,i,j) )
             call CHECK( __LINE__, DENS(k+1,i,j+1) )
             call CHECK( __LINE__, nu(k,i,j) )
             call CHECK( __LINE__, nu(k,i,j+1) )
             call CHECK( __LINE__, nu(k+1,i,j) )
             call CHECK( __LINE__, nu(k+1,i,j+1) )
             call CHECK( __LINE__, S23_X(k,i,j) )
#endif
             qflx_sgs_momy(k,i,j,ZDIR) = - 0.125_RP & ! 2/4/4
                  * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
                  * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i,j+1)+nu  (k+1,i,j+1) ) &
                  * S23_X(k,i,j)
          enddo
          enddo
          enddo
          !$omp end do nowait
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
          !$omp do
          do j = JJS, JJE
          do i = IIS, IIE
             qflx_sgs_momy(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
             qflx_sgs_momy(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
          enddo
          enddo
          !$omp end do nowait
#ifdef DEBUG
          i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       end if

       ! (z edge)
       !$omp do collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i+1,j+1) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, nu(k,i+1,j+1) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,XDIR) = - 0.125_RP & !
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ( nu  (k,i,j)+nu  (k,i+1,j)+nu  (k,i,j+1)+nu  (k,i+1,j+1) ) &
               * S12_Z(k,i,j)
       enddo
       enddo
       enddo
       !$omp end do nowait
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z-x plane)
       !$omp do collapse(2)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, TKE(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,YDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu(k,i,j) &
               * ( S22_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) &
             + twoOverThree * TKE(k,i,j) * tke_fact)
       enddo
       enddo
       enddo
       !$omp end do nowait
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !$omp end parallel

       if ( ATMOS_PHY_TB_SMG_implicit ) then
          call calc_tend_MOMY( TEND, & ! (out)
                               qflx_sgs_momy, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          !$omp parallel do collapse(2) &
          !$omp private(ap,d)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS  ,i,j  )*Nu(KS  ,i,j  ) &
                                   + DENS(KS+1,i,j  )*Nu(KS+1,i,j  ) &
                                   + DENS(KS  ,i,j+1)*Nu(KS  ,i,j+1) &
                                   + DENS(KS+1,i,j+1)*Nu(KS+1,i,j+1) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_XVW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_XVZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS,i,j+1) )
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XVZ)
                ap = - dt * 0.25_RP * ( DENS(k  ,i,j  )*Nu(k  ,i,j  ) &
                                      + DENS(k+1,i,j  )*Nu(k+1,i,j  ) &
                                      + DENS(k  ,i,j+1)*Nu(k  ,i,j+1) &
                                      + DENS(k+1,i,j+1)*Nu(k+1,i,j+1) ) &
                                    * RFDZ(k) / GSQRT(k,i,j,I_XVW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XVZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) )
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_XVZ)
             b(KE,i,j) = - c(KE,i,j) + 0.5_RP * ( DENS(KE,i,j)+DENS(KE,i,j+1) )

             do k = KS, KE
                d(k) = TEND(k,i,j)
             end do

             call diffusion_solver( &
                  TEND(:,i,j),                     & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d, & ! (in)
                  KE                               ) ! (in)

             do k = KS, KE-1
                qflx_sgs_momy(k,i,j,ZDIR) = qflx_sgs_momy(k,i,j,ZDIR) &
                     - 0.25_RP * ( DENS(k  ,i,j  )*Nu(k  ,i,j  ) &
                                 + DENS(k+1,i,j  )*Nu(k+1,i,j  ) &
                                 + DENS(k  ,i,j+1)*Nu(k  ,i,j+1) &
                                 + DENS(k+1,i,j+1)*Nu(k+1,i,j+1) ) &
                     * dt * ( TEND(k+1,i,j) - TEND(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_XVW)
             end do

          end do
          end do

       end if

       if ( ATMOS_PHY_TB_SMG_backscatter ) then

#define f2h(k,i,j,p) ( ( FZ(k+p-1,i,j) - FZ(k+p-2,i,j) ) / ( FZ(k+1,i,j) - FZ(k-1,i,j) ) )

          !$omp parallel private(flxz)

          ! MOMZ : dfy/dx - dfx/dy
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS+1, KE-1
                flxz(k) = J13G(k,i,j,I_XYZ) * random_my(k,i,j) &
                        - J23G(k,i,j,I_XYZ) * random_mx(k,i,j)
             end do
             flxz(KS) = 0.0_RP
             flxz(KE) = 0.0_RP
             do k = KS, KE-1
                MOMZ_t(k,i,j) = ( ( GSQRT(k,i+1,j,I_XYW) * ( f2h(k,i+1,j,1) * random_my(k+1,i+1,j) &
                                                           + f2h(k,i+1,j,2) * random_my(k,i+1,j) ) &
                                  - GSQRT(k,i-1,j,I_XYW) * ( f2h(k,i-1,j,1) * random_my(k+1,i-1,j) &
                                                           + f2h(k,i-1,j,2) * random_my(k,i-1,j) ) &
                                  ) / ( FDX(i) + FDX(i-1) ) * MAPF(i,j,1,I_XY) &
                                - ( GSQRT(k,i,j+1,I_XYW) * ( f2h(k,i,j+1,1) * random_mx(k+1,i,j+1) &
                                                           + f2h(k,i,j+1,2) * random_mx(k,i,j+1) ) &
                                  - GSQRT(k,i,j-1,I_XYW) * ( f2h(k,i,j-1,1) * random_mx(k+1,i,j-1) &
                                                           + f2h(k,i,j+1,2) * random_mx(k,i,j-1) ) &
                                  ) / ( FDY(j) + FDY(j-1) ) * MAPF(i,j,2,I_XY) &
                             + ( flxz(k+1) - flxz(k) ) * RFDZ(k) &
                                ) / GSQRT(k,i,j,I_XYW)
             end do
             MOMZ_t(KE,i,j) = 0.0_RP
          end do
          end do
          !$omp end do nowait

          ! MOMX : dfz/dy - dfy/dz
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS, KE-1
                flxz(k) = J23G(k,i,j,I_UYW) * ( f2h(k,i+1,j,1) * random_mz(k+1,i+1,j) + f2h(k,i,j+1,2) * random_mz(k,i+1,j) &
                                              + f2h(k,i  ,j,1) * random_mz(k+1,i  ,j) + f2h(k,i  ,j,2) * random_mz(k,i  ,j) ) * 0.5_RP &
                        - J33G * ( f2h(k,i+1,j,1) * random_my(k+1,i+1,j) + f2h(k,i+1,j,2) * random_my(k,i+1,j) &
                                 + f2h(k,i  ,j,1) * random_my(k+1,i  ,j) + f2h(k,i  ,j,2) * random_my(k,i  ,j) ) * 0.5_RP
             end do
             flxz(KS-1) = 0.0_RP
             flxz(KE  ) = 0.0_RP
             do k = KS, KE
                MOMX_t(k,i,j) = ( ( GSQRT(k,i,j+1,I_UYZ) * ( random_mz(k,i+1,j+1) + random_mz(k,i,j+1) ) * 0.5_RP &
                                  - GSQRT(k,i,j-1,I_UYZ) * ( random_mz(k,i+1,j-1) + random_mz(k,i,j-1) ) * 0.5_RP ) / ( FDY(j) + FDY(j-1) ) * MAPF(i,j,2,I_UY) &
                                  + ( flxz(k) - flxz(k-1) ) * RCDZ(k) &
                                ) / GSQRT(k,i,j,I_UYZ)
             end do
          end do
          end do
          !$omp end do nowait

          ! MOMY : dfx/dz - dfz/dx
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS, KE-1
                flxz(k) = J33G * ( f2h(k,i,j+1,1) * random_mx(k+1,i,j+1) + f2h(k,i,j+1,2) * random_mx(k,i,j+1) &
                                 + f2h(k,i,j  ,1) * random_mx(k+1,i,j  ) + f2h(k,i,j  ,2) * random_mx(k,i,j  ) ) * 0.5_RP &
                        - J13G(k,i,j,I_XVW) * ( f2h(k,i,j+1,1) * random_mz(k+1,i,j+1) + f2h(k,i,j+1,2) * random_mz(k,i,j+1) &
                                              + f2h(k,i,j  ,1) * random_mz(k+1,i,j  ) + f2h(k,i,j  ,2) * random_mz(k,i,j  ) ) * 0.5_RP
             end do
             flxz(KS-1) = 0.0_RP
             flxz(KE  ) = 0.0_RP
             do k = KS, KE
                MOMY_t(k,i,j) = ( ( flxz(k) - flxz(k-1) ) * RCDZ(k) &
                                  - ( GSQRT(k,i+1,j,I_XVZ) * ( random_mz(k,i+1,j+1) + random_mz(k,i+1,j) ) * 0.5_RP &
                                    - GSQRT(k,i-1,j,I_XYZ) * ( random_mz(k,i-1,j+1) + random_mz(k,i-1,j) ) * 0.5_RP ) / ( FDX(i) + FDX(i-1) ) * MAPF(i,j,1,I_XV) &
                                ) / GSQRT(k,i,j,I_XVZ)
             end do
          end do
          end do
          !$omp end do nowait

          !$omp end parallel

       else

          !$omp parallel

!OCL XFILL
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             MOMZ_t(k,i,j) = 0.0_RP
          end do
          end do
          end do
          !$omp end do nowait

!OCL XFILL
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             MOMX_t(k,i,j) = 0.0_RP
          end do
          end do
          end do
          !$omp end do nowait

!OCL XFILL
          !$omp do collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             MOMY_t(k,i,j) = 0.0_RP
          end do
          end do
          end do
          !$omp end do nowait

          !$omp end parallel

       end if

       !##### Thermodynamic Equation #####

       if ( ATMOS_PHY_TB_SMG_implicit ) then

          !$omp parallel do collapse(2) &
          !$omp private (ap,d)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS,i,j)+DENS(KS+1,i,j) ) &
                                 * ( Kh(KS,i,j)+Kh(KS+1,i,j) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + DENS(KS,i,j)
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
                ap = - dt * 0.25_RP * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                                    * ( Kh(k,i,j)+Kh(k+1,i,j) ) &
                                   * RFDZ(k) / GSQRT(k,i,j,I_XYW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + DENS(k,i,j)
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ)
             b(KE,i,j) = - c(KE,i,j) + DENS(KE,i,j)

          end do
          end do

       end if

       call calc_flux_phi( &
            qflx_sgs_rhot, &
            DENS, POTT, Kh, 1.0_RP, &
            GSQRT, J13G, J23G, J33G, MAPF, &
            ATMOS_PHY_TB_SMG_horizontal, &
            ATMOS_PHY_TB_SMG_implicit, &
            a, b, c, dt, &
            IIS, IIE, JJS, JJE )

       if ( ATMOS_PHY_TB_SMG_backscatter ) then

          !$omp parallel do collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1

             do k = KS+1, KE-1
                dd(k,i,j) = sqrt( ( ( POTT(k+1,i,j) - POTT(k-1,i,j) ) * J33G / ( FDZ(k) + FDZ(k-1) ) )**2 &
                                + ( ( ( GSQRT(k,i+1,j,I_XYZ)*POTT(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ)*POTT(k,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                                    + ( J13G(k+1,i,j,I_XYZ)*POTT(k+1,i,j) - J13G(k-1,i,j,I_XYZ)*POTT(k-1,i,j) ) / ( FDZ(k) + FDZ(k-1) ) ) * MAPF(i,j,1,I_XY) )**2 &
                                + ( ( ( GSQRT(k,i,j+1,I_XYZ)*POTT(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ)*POTT(k,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
                                    + ( J23G(k+1,i,j,I_XYZ)*POTT(k+1,i,j) - J23G(k-1,i,j,I_XYZ)*POTT(k-1,i,j) ) / ( FDZ(k) + FDZ(k-1) ) ) * MAPF(i,j,2,I_XY) )**2 &
                            ) / GSQRT(k,i,j,I_XYZ)
             end do
             dd(KS,i,j) = sqrt( ( ( POTT(KS+1,i,j) - POTT(KS,i,j) ) * J33G * RFDZ(KS) )**2 &
                              + ( ( ( GSQRT(KS,i+1,j,I_XYZ)*POTT(KS,i+1,j) - GSQRT(KS,i-1,j,I_XYZ)*POTT(KS,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                                  + ( J13G(KS+1,i,j,I_XYZ)*POTT(KS+1,i,j) - J13G(KS,i,j,I_XYZ)*POTT(KS,i,j) ) * RFDZ(KS) ) * MAPF(i,j,1,I_XY) )**2 &
                              + ( ( ( GSQRT(KS,i,j+1,I_XYZ)*POTT(KS,i,j+1) - GSQRT(KS,i,j-1,I_XYZ)*POTT(KS,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
                                  + ( J23G(KS+1,i,j,I_XYZ)*POTT(KS+1,i,j) - J23G(KS,i,j,I_XYZ)*POTT(KS,i,j) ) * RFDZ(KS) ) * MAPF(i,j,2,I_XY) )**2 &
                              ) / GSQRT(KS,i,j,I_XYZ)
             dd(KE,i,j) = sqrt( ( ( POTT(KE,i,j) - POTT(KE-1,i,j) ) * J33G * RFDZ(KE-1) )**2 &
                              + ( ( ( GSQRT(KE,i+1,j,I_XYZ)*POTT(KE,i+1,j) - GSQRT(KE,i-1,j,I_XYZ)*POTT(KE,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                                  + ( J13G(KE,i,j,I_XYZ)*POTT(KE,i,j) - J13G(KE-1,i,j,I_XYZ)*POTT(KE-1,i,j) ) * RFDZ(KE-1) ) * MAPF(i,j,1,I_XY) )**2 &
                              + ( ( ( GSQRT(KE,i,j+1,I_XYZ)*POTT(KE,i,j+1) - GSQRT(KE,i,j-1,I_XYZ)*POTT(KE,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
                                  + ( J23G(KE,i,j,I_XYZ)*POTT(KE,i,j) - J23G(KE-1,i,j,I_XYZ)*POTT(KE-1,i,j) ) * RFDZ(KE-1) ) * MAPF(i,j,2,I_XY) )**2 &
                              ) / GSQRT(KE,i,j,I_XYZ)

          end do
          end do

          !$omp parallel do collapse(2) &
          !$omp private(flxz)
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS, KE-1
                flxz(k) = J33G * ( f2h(k,i,j,1) * random_qz(k+1,i,j) * dd(k+1,i,j) + f2h(k,i,j,2) * random_qz(k,i,j) * dd(k,i,j) ) &
                        + J13G(k,i,j,I_XYW) * ( f2h(k,i,j,1) * random_qx(k+1,i,j) * dd(k+1,i,j) + f2h(k,i,j,2) * random_qx(k,i,j) * dd(k,i,j) ) &
                        + J23G(k,i,j,I_XYW) * ( f2h(k,i,j,1) * random_qy(k+1,i,j) * dd(k+1,i,j) + f2h(k,i,j,2) * random_qy(k,i,j) * dd(k,i,j) )
             end do
             flxz(KS-1) = 0.0_RP
             flxz(KE  ) = 0.0_RP
             do k = KS, KE
                RHOT_t(k,i,j) = ( ( flxz(k) - flxz(k-1) ) * RCDZ(k) &
                                + ( GSQRT(k,i+1,j,I_XYZ) * random_qx(k,i+1,j) * dd(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ) * random_qx(k,i-1,j) * dd(k,i-1,j) ) / ( FDX(i) + FDX(i-1) ) * MAPF(i,j,1,I_XY) &
                                + ( GSQRT(k,i,j+1,I_XYZ) * random_qy(k,i,j+1) * dd(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ) * random_qy(k,i,j-1) * dd(k,i,j-1) ) / ( FDY(j) + FDY(j-1) ) * MAPF(i,j,2,I_XY) &
                                ) / GSQRT(k,i,j,I_XYZ)
             end do
          end do
          end do

       else

          !$omp parallel do collapse(2)
!OCL XFILL
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             RHOT_t(k,i,j) = 0.0_RP
          end do
          end do
          end do

       end if

    enddo
    enddo

    !##### Tracers #####
    do iq = 1, QA

       if ( .not. TRACER_ADVC(iq) ) cycle

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          call calc_flux_phi( &
               qflx_sgs_rhoq(:,:,:,:,iq), &
               DENS, QTRC(:,:,:,iq), Kh, 1.0_RP, &
               GSQRT, J13G, J23G, J33G, MAPF, &
               ATMOS_PHY_TB_SMG_horizontal, &
               ATMOS_PHY_TB_SMG_implicit, &
               a, b, c, dt, &
               IIS, IIE, JJS, JJE )


          if ( ATMOS_PHY_TB_SMG_backscatter .and. iq == I_QV ) then

             !$omp parallel do collapse(2)
             do j = JJS-1, JJE+1
             do i = IIS-1, IIE+1

                do k = KS+1, KE-1
                   dd(k,i,j) = sqrt( ( ( QTRC(k+1,i,j,iq) - QTRC(k-1,i,j,iq) ) * J33G / ( FDZ(k) + FDZ(k-1) ) )**2 &
                                   + ( ( ( GSQRT(k,i+1,j,I_XYZ)*QTRC(k,i+1,j,iq) - GSQRT(k,i-1,j,I_XYZ)*QTRC(k,i-1,j,iq) ) / ( FDX(i) + FDX(i-1) ) &
                                       + ( J13G(k+1,i,j,I_XYZ)*QTRC(k+1,i,j,iq) - J13G(k-1,i,j,I_XYZ)*QTRC(k-1,i,j,iq) ) / ( FDZ(k) + FDZ(k-1) ) ) * MAPF(i,j,1,I_XY) )**2 &
                                   + ( ( ( GSQRT(k,i,j+1,I_XYZ)*QTRC(k,i,j+1,iq) - GSQRT(k,i,j-1,I_XYZ)*QTRC(k,i,j-1,iq) ) / ( FDY(j) + FDY(j-1) ) &
                                       + ( J23G(k+1,i,j,I_XYZ)*QTRC(k+1,i,j,iq) - J23G(k-1,i,j,I_XYZ)*QTRC(k-1,i,j,iq) ) / ( FDZ(k) + FDZ(k-1) ) ) * MAPF(i,j,2,I_XY) )**2 &
                               ) / GSQRT(k,i,j,I_XYZ)
                end do
                dd(KS,i,j) = sqrt( ( ( QTRC(KS+1,i,j,iq) - QTRC(KS,i,j,iq) ) * J33G * RFDZ(KS) )**2 &
                                 + ( ( ( GSQRT(KS,i+1,j,I_XYZ)*QTRC(KS,i+1,j,iq) - GSQRT(KS,i-1,j,I_XYZ)*QTRC(KS,i-1,j,iq) ) / ( FDX(i) + FDX(i-1) ) &
                                     + ( J13G(KS+1,i,j,I_XYZ)*QTRC(KS+1,i,j,iq) - J13G(KS,i,j,I_XYZ)*QTRC(KS,i,j,iq) ) * RFDZ(KS) ) * MAPF(i,j,1,I_XY) )**2 &
                                 + ( ( ( GSQRT(KS,i,j+1,I_XYZ)*QTRC(KS,i,j+1,iq) - GSQRT(KS,i,j-1,I_XYZ)*QTRC(KS,i,j-1,iq) ) / ( FDY(j) + FDY(j-1) ) &
                                     + ( J23G(KS+1,i,j,I_XYZ)*QTRC(KS+1,i,j,iq) - J23G(KS,i,j,I_XYZ)*QTRC(KS,i,j,iq) ) * RFDZ(KS) ) * MAPF(i,j,2,I_XY) )**2 &
                                 ) / GSQRT(KS,i,j,I_XYZ)
                dd(KE,i,j) = sqrt( ( ( QTRC(KE,i,j,iq) - QTRC(KE-1,i,j,iq) ) * J33G * RFDZ(KE-1) )**2 &
                                 + ( ( ( GSQRT(KE,i+1,j,I_XYZ)*QTRC(KE,i+1,j,iq) - GSQRT(KE,i-1,j,I_XYZ)*QTRC(KE,i-1,j,iq) ) / ( FDX(i) + FDX(i-1) ) &
                                     + ( J13G(KE,i,j,I_XYZ)*QTRC(KE,i,j,iq) - J13G(KE-1,i,j,I_XYZ)*QTRC(KE-1,i,j,iq) ) * RFDZ(KE-1) ) * MAPF(i,j,1,I_XY) )**2 &
                                 + ( ( ( GSQRT(KE,i,j+1,I_XYZ)*QTRC(KE,i,j+1,iq) - GSQRT(KE,i,j-1,I_XYZ)*QTRC(KE,i,j-1,iq) ) / ( FDY(j) + FDY(j-1) ) &
                                     + ( J23G(KE,i,j,I_XYZ)*QTRC(KE,i,j,iq) - J23G(KE-1,i,j,I_XYZ)*QTRC(KE-1,i,j,iq) ) * RFDZ(KE-1) ) * MAPF(i,j,2,I_XY) )**2 &
                                 ) / GSQRT(KE,i,j,I_XYZ)

             end do
             end do

             !$omp parallel do collapse(2) &
             !$omp private (flxz)
             do j = JJS, JJE
             do i = IIS, IIE
                do k = KS, KE-1
                   flxz(k) = J33G * ( f2h(k,i,j,1) * random_qz(k+1,i,j) * dd(k+1,i,j) + f2h(k,i,j,2) * random_qz(k,i,j) * dd(k,i,j) ) &
                           + J13G(k,i,j,I_XYW) * ( f2h(k,i,j,1) * random_qx(k+1,i,j) * dd(k+1,i,j) + f2h(k,i,j,2) * random_qx(k,i,j) * dd(k,i,j) ) &
                           + J23G(k,i,j,I_XYW) * ( f2h(k,i,j,1) * random_qy(k+1,i,j) * dd(k+1,i,j) + f2h(k,i,j,2) * random_qy(k,i,j) * dd(k,i,j) )
                end do
                flxz(KS-1) = 0.0_RP
                flxz(KE  ) = 0.0_RP
                do k = KS, KE
                   RHOQ_t(k,i,j,iq) = ( ( flxz(k) - flxz(k-1) ) * RCDZ(k) &
                                      + ( GSQRT(k,i+1,j,I_XYZ) * random_qx(k,i+1,j) * dd(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ) * random_qx(k,i-1,j) * dd(k,i-1,j) ) / ( FDX(i) + FDX(i-1) ) * MAPF(i,j,1,I_XY) &
                                      + ( GSQRT(k,i,j+1,I_XYZ) * random_qy(k,i,j+1) * dd(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ) * random_qy(k,i,j-1) * dd(k,i,j-1) ) / ( FDY(j) + FDY(j-1) ) * MAPF(i,j,2,I_XY) &
                                      ) / GSQRT(k,i,j,I_XYZ)
                end do
             end do
             end do

          else

             !$omp parallel do collapse(2)
!OCL XFILL
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                RHOQ_t(k,i,j,iq) = 0.0_RP
             end do
             end do
             end do

          end if

       enddo
       enddo

    enddo ! scalar quantities loop
#ifdef DEBUG
    iq = IUNDEF
#endif

    call FILE_HISTORY_in( TKE(:,:,:), 'TKE_SMG', 'turbulent kinetic energy (Smagorinsky)', 'm2/s2', fill_halo=.true. )


    return
  end subroutine ATMOS_PHY_TB_smg


  function mixlen(dz, dx, dy, filter_fact)
    implicit none
    real(RP), intent(in) :: dz
    real(RP), intent(in) :: dx
    real(RP), intent(in) :: dy
    real(RP), intent(in) :: filter_fact
    real(RP) :: mixlen ! (out)

    mixlen = fact(dz, dx, dy) * filter_fact * ( dz * dx * dy )**OneOverThree ! Scotti et al. (1993)

    return
  end function mixlen

  function fact(dz, dx, dy)
    real(RP), intent(in) :: dz
    real(RP), intent(in) :: dx
    real(RP), intent(in) :: dy
    real(RP) :: fact ! (out)

    real(RP), parameter :: oot = -1.0_RP/3.0_RP
    real(RP), parameter :: fot =  5.0_RP/3.0_RP
    real(RP), parameter :: eot = 11.0_RP/3.0_RP
    real(RP), parameter :: tof = -3.0_RP/4.0_RP
    real(RP) :: a1, a2, b1, b2, dmax


    dmax = max(dz, dx, dy)
    if ( dz == dmax ) then
       a1 = dx / dmax
       a2 = dy / dmax
    else if ( dx == dmax ) then
       a1 = dz / dmax
       a2 = dy / dmax
    else ! dy == dmax
       a1 = dz / dmax
       a2 = dx / dmax
    end if
    b1 = atan( a1/a2 )
    b2 = atan( a2/a1 )

   fact = 1.736_RP * (a1*a2)**oot &
         * ( 4.0_RP*p1(b1)*a1**oot + 0.222_RP*p2(b1)*a1**fot + 0.077*p3(b1)*a1**eot - 3.0_RP*b1 &
           + 4.0_RP*p1(b2)*a2**oot + 0.222_RP*p2(b2)*a2**fot + 0.077*p3(b2)*a2**eot - 3.0_RP*b2 &
           )**tof
   return
  end function fact
  function p1(z)
    real(RP), intent(in) :: z
    real(RP) :: p1 ! (out)

    p1 = 2.5_RP * p2(z) - 1.5_RP * sin(z) * cos(z)**TwoOverThree
    return
  end function p1
  function p2(z)
    real(RP), intent(in) :: z
    real(RP) :: p2 ! (out)

    p2 = 0.986_RP * z + 0.073_RP * z**2 - 0.418_RP * z**3 + 0.120_RP * z**4
    return
  end function p2
  function p3(z)
    real(RP), intent(in) :: z
    real(RP) :: p3 ! (out)

    p3 = 0.976_RP * z + 0.188_RP * z**2 - 1.169_RP * z**3 + 0.755_RP * z**4 - 0.151_RP * z**5
    return
  end function p3

end module scale_atmos_phy_tb_smg
