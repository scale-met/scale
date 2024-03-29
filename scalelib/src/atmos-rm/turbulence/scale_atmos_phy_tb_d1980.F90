!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!          1.5th TKE model Deardorff (1980)
!!
!! @author Team SCALE
!!
!! @par Reference
!!  - Deardorff, 1980:
!!    Stratocumulus-capped mixed layers derived from a three-dimensional model.
!!    Bound.-Layer Meteor., 18, 495-527
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_tb_d1980
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
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
  public :: ATMOS_PHY_TB_d1980_config
  public :: ATMOS_PHY_TB_d1980_setup
  public :: ATMOS_PHY_TB_d1980_finalize
  public :: ATMOS_PHY_TB_d1980

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
  real(RP), private, parameter   :: C1 = 0.19_RP
  real(RP), private, parameter   :: C2 = 0.51_RP

  real(RP), private, parameter   :: OneOverThree  = 1.0_RP / 3.0_RP
  real(RP), private, parameter   :: TwoOverThree  = 2.0_RP / 3.0_RP
  real(RP), private, parameter   :: FourOverThree = 4.0_RP / 3.0_RP

  real(RP), private, allocatable :: delta(:,:,:) ! (dx*dy*dz)^(1/3)

  integer,  private :: I_TKE

  real(RP), private :: ATMOS_PHY_TB_D1980_Km_MAX = 1000.0_RP
  logical,  private :: ATMOS_PHY_TB_D1980_implicit = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_TB_d1980_config( &
       TYPE_TB,  &
       I_TKE_out )
    use scale_prc, only: &
       PRC_abort
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    character(len=*), intent(in)  :: TYPE_TB
    integer,          intent(out) :: I_TKE_out
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_d1980_config",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_d1980_config",*) 'Tracers for Deardorff (1980) 1.5th TKE Model'

    if ( TYPE_TB /= 'D1980' ) then
       LOG_ERROR("ATMOS_PHY_TB_d1980_config",*) 'ATMOS_PHY_TB_TYPE is not D1980. Check!'
       call PRC_abort
    endif

    call TRACER_regist( I_TKE,                                    & ! [OUT]
                        1,                                        & ! [IN]
                        (/ 'TKE_D1980' /),                        & ! [IN]
                        (/ 'turbulent kinetic energy (D1980)' /), & ! [IN]
                        (/ 'm2/s2' /)                             ) ! [IN]

    I_TKE_out = I_TKE

    return
  end subroutine ATMOS_PHY_TB_d1980_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_TB_d1980_setup( &
       CDZ, CDX, CDY, CZ )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA,IA,JA)

    namelist / PARAM_ATMOS_PHY_TB_D1980 / &
       ATMOS_PHY_TB_D1980_Km_MAX,   &
       ATMOS_PHY_TB_D1980_implicit

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_d1980_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_d1980_setup",*) 'Deardorff (1980) 1.5th TKE Model'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_D1980,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_d1980_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_d1980_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_D1980. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_D1980)

    allocate( delta(KA,IA,JA) )

#ifdef DEBUG
    delta(:,:,:) = UNDEF
#endif
    !$omp parallel do
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
       delta(k,i,j) = ( CDZ(k) * CDX(i) * CDY(j) )**(1.0_RP/3.0_RP)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB_d1980_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_TB_d1980_finalize

    deallocate( delta )

    return
  end subroutine ATMOS_PHY_TB_d1980_finalize

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_d1980( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, &
       qflx_sgs_rhot, qflx_sgs_rhoq,                &
       RHOQ_t,                                      &
       Km, Ri, Pr,                                  &
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2,      &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_Q,  &
       GSQRT, J13G, J23G, J33G, MAPF, dt            )
    use scale_precision
    use scale_atmos_grid_cartesC_index
    use scale_tracer
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_atmos_grid_cartesC, only: &
       FDZ  => ATMOS_GRID_CARTESC_FDZ,  &
       RCDZ => ATMOS_GRID_CARTESC_RCDZ, &
       RCDX => ATMOS_GRID_CARTESC_RCDX, &
       RCDY => ATMOS_GRID_CARTESC_RCDY, &
       RFDZ => ATMOS_GRID_CARTESC_RFDZ, &
       RFDX => ATMOS_GRID_CARTESC_RFDX, &
       RFDY => ATMOS_GRID_CARTESC_RFDY
    use scale_matrix, only: &
       MATRIX_SOLVER_tridiagonal
    use scale_atmos_phy_tb_common, only: &
       calc_strain_tensor => ATMOS_PHY_TB_calc_strain_tensor, &
       calc_tend_momz     => ATMOS_PHY_TB_calc_tend_momz,     &
       calc_tend_momx     => ATMOS_PHY_TB_calc_tend_momx,     &
       calc_tend_momy     => ATMOS_PHY_TB_calc_tend_momy,     &
       calc_tend_phi      => ATMOS_PHY_TB_calc_tend_phi,      &
       calc_flux_phi      => ATMOS_PHY_TB_calc_flux_phi
    implicit none

    ! SGS flux
    real(RP), intent(out)   :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

    real(RP), intent(inout) :: RHOQ_t       (KA,IA,JA,QA) ! tendency of rho * QTRC

    real(RP), intent(out)   :: Km           (KA,IA,JA)    ! eddy viscosity (center)
    real(RP), intent(out)   :: Ri           (KA,IA,JA)    ! Richardson number
    real(RP), intent(out)   :: Pr           (KA,IA,JA)    ! Prantle number

    real(RP), intent(in)    :: MOMZ         (KA,IA,JA)
    real(RP), intent(in)    :: MOMX         (KA,IA,JA)
    real(RP), intent(in)    :: MOMY         (KA,IA,JA)
    real(RP), intent(in)    :: RHOT         (KA,IA,JA)
    real(RP), intent(in)    :: DENS         (KA,IA,JA)
    real(RP), intent(in)    :: QTRC         (KA,IA,JA,QA)
    real(RP), intent(in)    :: N2           (KA,IA,JA)

    real(RP), intent(in)    :: SFLX_MW      (IA,JA)
    real(RP), intent(in)    :: SFLX_MU      (IA,JA)
    real(RP), intent(in)    :: SFLX_MV      (IA,JA)
    real(RP), intent(in)    :: SFLX_SH      (IA,JA)
    real(RP), intent(in)    :: SFLX_Q       (IA,JA,QA)

    real(RP), intent(in)    :: GSQRT        (KA,IA,JA,7)  !< vertical metrics {G}^1/2
    real(RP), intent(in)    :: J13G         (KA,IA,JA,7)  !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G         (KA,IA,JA,7)  !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J33G                       !< (3,3) element of Jacobian matrix
    real(RP), intent(in)    :: MAPF         (IA,JA,2,4)   !< map factor
    real(DP), intent(in)    :: dt

    ! diagnostic variables
    real(RP) :: POTT   (KA,IA,JA)

    ! deformation rate tensor
    real(RP) :: S33_C (KA,IA,JA) ! (cell center)
    real(RP) :: S11_C (KA,IA,JA)
    real(RP) :: S22_C (KA,IA,JA)
    real(RP) :: S31_C (KA,IA,JA)
    real(RP) :: S12_C (KA,IA,JA)
    real(RP) :: S23_C (KA,IA,JA)
    real(RP) :: S12_Z (KA,IA,JA) ! (z edge or x-y plane)
    real(RP) :: S23_X (KA,IA,JA) ! (x edge or y-z plane)
    real(RP) :: S31_Y (KA,IA,JA) ! (y edge or z-x plane)
    real(RP) :: S2    (KA,IA,JA) ! |S|^2

    real(RP) :: Kh    (KA,IA,JA) ! eddy diffusion
    real(RP) :: l     (KA,IA,JA) ! mixing length

    real(RP) :: qflx_tke(KA,IA,JA,3)

    real(RP) :: TEND(KA,IA,JA)
    real(RP) :: TEND2(KA,IA,JA)
    real(RP) :: a(KA,IA,JA)
    real(RP) :: b(KA,IA,JA)
    real(RP) :: c(KA,IA,JA)
    real(RP) :: ap

    integer  :: IIS, IIE
    integer  :: JJS, JJE

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / turbulence / D1980'

#ifdef DEBUG
    qflx_sgs_momz(:,:,:,:)   = UNDEF
    qflx_sgs_momx(:,:,:,:)   = UNDEF
    qflx_sgs_momy(:,:,:,:)   = UNDEF
    qflx_sgs_rhot(:,:,:,:)   = UNDEF
    qflx_sgs_rhoq(:,:,:,:,:) = UNDEF

    Pr (:,:,:) = UNDEF
    Ri (:,:,:) = UNDEF
    Kh (:,:,:) = UNDEF
    Km (:,:,:) = UNDEF

    POTT   (:,:,:) = UNDEF
#endif

    ! potential temperature
    !$omp parallel do
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, RHOT(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
    enddo
    enddo
    enddo
#ifdef DEBUG
    i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    !##### Start Upadate #####

    call calc_strain_tensor( &
         S33_C, S11_C, S22_C,          & ! (out)
         S31_C, S12_C, S23_C,          & ! (out)
         S12_Z, S23_X, S31_Y,          & ! (out)
         S2                 ,          & ! (out)
         DENS, MOMZ, MOMX, MOMY,       & ! (in)
         GSQRT, J13G, J23G, J33G, MAPF ) ! (in)

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! Ri = N^2 / |S|^2, N^2 = g / theta * dtheta/dz
       !$omp parallel do
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          Ri(k,i,j) = N2(k,i,j) / S2(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! mixing length
       !$omp parallel do
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          if ( N2(k,i,j) > 1e-10_RP ) then
             l(k,i,j) = max( min( 0.76_RP * sqrt(QTRC(k,i,j,I_TKE)/N2(k,i,j)), delta(k,i,j) ), 1e-10_RP )
          else
             l(k,i,j) = delta(k,i,j)
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !$omp parallel do
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          Km(k,i,j) = min( 0.10_RP * l(k,i,j) * sqrt(QTRC(k,i,j,I_TKE)), ATMOS_PHY_TB_D1980_Km_MAX )
          Pr(k,i,j) = 1.0_RP / ( 1.0_RP + 2.0_RP * l(k,i,j) / delta(k,i,j) )
          Kh(k,i,j) = Km(k,i,j) / Pr(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !##### momentum equation (z) #####
       ! (cell center)
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,ZDIR) = DENS(k,i,j) * ( - 2.0_RP * Km(k,i,j) * S33_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momz(KS,i,j,ZDIR) = 0.0_RP ! just above bottom boundary
          qflx_sgs_momz(KE,i,j,ZDIR) = 0.0_RP ! just below top boundary
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       !$omp parallel do
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i+1,j) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, Km(k,i+1,j) )
       call CHECK( __LINE__, Km(k+1,i,j) )
       call CHECK( __LINE__, Km(k+1,i+1,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,XDIR) = - 0.125_RP & ! 2.0 / 4 / 4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
               * ( Km  (k,i,j)+Km  (k+1,i,j)+Km  (k,i+1,j)+Km  (k+1,i+1,j)) &
               * S31_Y(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       !$omp parallel do
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j+1) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, Km(k,i,j+1) )
       call CHECK( __LINE__, Km(k+1,i,j) )
       call CHECK( __LINE__, Km(k+1,i,j+1) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,YDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
               * ( Km  (k,i,j)+Km  (k+1,i,j)+Km  (k,i,j+1)+Km  (k+1,i,j+1) ) &
               * S23_X(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       if ( ATMOS_PHY_TB_D1980_implicit ) then
          call calc_tend_MOMZ( TEND, & ! (out)
                               qflx_sgs_momz, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          !$omp parallel do &
          !$omp private(ap)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - FourOverThree * dt &
                  * DENS(KS+1,i,j)*Km(KS+1,i,j) &
                  * RCDZ(KS+1) / GSQRT(KS+1,i,j,I_XYZ)
             a(KS,i,j) = ap * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS+1,i,j) )
             do k = KS+1, KE-2
                c(k,i,j) = ap * RFDZ(k+1) / GSQRT(k+1,i,j,I_XYW)
                ap = - FourOverThree * dt &
                     * DENS(k+1,i,j)*Km(k+1,i,j) &
                     * RCDZ(k+1) / GSQRT(k+1,i,j,I_XYZ)
                a(k,i,j) = ap * RFDZ(k) / GSQRT(k,i,j,I_XYW)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k+1,i,j) )
             end do
             a(KE-1,i,j) = 0.0_RP
             c(KE-1,i,j) = ap * RFDZ(KE) / GSQRT(KE,i,j,I_XYW)
             b(KE-1,i,j) = - c(KE-1,i,j) + 0.5_RP * ( DENS(KE-1,i,j)+DENS(KE,i,j) )
          end do
          end do

          call MATRIX_SOLVER_tridiagonal( KA, KS, KE-1, IA, IIS, IIE, JA, JJS, JJE, &
                                          a(:,:,:), b(:,:,:), c(:,:,:), TEND(:,:,:), & ! (in)
                                          TEND2(:,:,:)                               ) ! (out)

          !$omp parallel do
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS+1, KE-1
                qflx_sgs_momz(k,i,j,ZDIR) = qflx_sgs_momz(k,i,j,ZDIR) &
                     - FourOverThree * DENS(k,i,j) * Km(k,i,j) * dt &
                     * ( TEND2(k,i,j) - TEND2(k-1,i,j) ) * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
             end do

          end do
          end do

       end if

       !##### momentum equation (x) #####
       ! (y edge)
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i+1,j) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, Km(k,i+1,j) )
       call CHECK( __LINE__, Km(k+1,i,j) )
       call CHECK( __LINE__, Km(k+1,i+1,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,ZDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
               * ( Km  (k,i,j)+Km  (k+1,i,j)+Km  (k,i+1,j)+Km  (k+1,i+1,j) ) &
               * S31_Y(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momx(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_momx(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (cell center)
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,XDIR) = DENS(k,i,j) * ( - 2.0_RP * Km(k,i,j) * S11_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       !$omp parallel do
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i+1,j+1) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, Km(k,i+1,j) )
       call CHECK( __LINE__, Km(k,i,j+1) )
       call CHECK( __LINE__, Km(k,i+1,j+1) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,YDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ( Km  (k,i,j)+Km  (k,i+1,j)+Km  (k,i,j+1)+Km  (k,i+1,j+1) ) &
               * S12_Z(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       if ( ATMOS_PHY_TB_D1980_implicit ) then
          call calc_tend_MOMX( TEND, & ! (out)
                               qflx_sgs_momx, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          !$omp parallel do &
          !$omp private(ap)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS  ,i  ,j)*Km(KS  ,i  ,j) &
                                   + DENS(KS+1,i  ,j)*Km(KS+1,i  ,j) &
                                   + DENS(KS  ,i+1,j)*Km(KS  ,i+1,j) &
                                   + DENS(KS+1,i+1,j)*Km(KS+1,i+1,j) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_UYW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_UYZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS,i+1,j) )
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_UYZ)
                ap = - dt * 0.25_RP * ( DENS(k  ,i  ,j)*Km(k  ,i  ,j) &
                                      + DENS(k+1,i  ,j)*Km(k+1,i  ,j) &
                                      + DENS(k  ,i+1,j)*Km(k  ,i+1,j) &
                                      + DENS(k+1,i+1,j)*Km(k+1,i+1,j) ) &
                                    * RFDZ(k) / GSQRT(k,i,j,I_UYW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_UYZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) )
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_UYZ)
             b(KE,i,j) = - c(KE,i,j) + 0.5_RP * ( DENS(KE,i,j)+DENS(KE,i+1,j) )
          end do
          end do

          call MATRIX_SOLVER_tridiagonal( KA, KS, KE, IA, IIS, IIE, JA, JJS, JJE, &
                                          a(:,:,:), b(:,:,:), c(:,:,:), TEND(:,:,:), & ! (in)
                                          TEND2(:,:,:)                               ) ! (out)

          !$omp parallel do
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS, KE-1
                qflx_sgs_momx(k,i,j,ZDIR) = qflx_sgs_momx(k,i,j,ZDIR) &
                     - 0.25_RP * ( DENS(k  ,i  ,j)*Km(k  ,i  ,j) &
                                 + DENS(k+1,i  ,j)*Km(k+1,i  ,j) &
                                 + DENS(k  ,i+1,j)*Km(k  ,i+1,j) &
                                 + DENS(k+1,i+1,j)*Km(k+1,i+1,j) ) &
                     * dt * ( TEND2(k+1,i,j) - TEND2(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_UYW)
             end do

          end do
          end do

       end if

       !##### momentum equation (y) #####
       ! (x edge)
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j+1) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, Km(k,i,j+1) )
       call CHECK( __LINE__, Km(k+1,i,j) )
       call CHECK( __LINE__, Km(k+1,i,j+1) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,ZDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
               * ( Km  (k,i,j)+Km  (k+1,i,j)+Km  (k,i,j+1)+Km  (k+1,i,j+1) ) &
               * S23_X(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momy(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_momy(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z edge)
       !$omp parallel do
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i+1,j+1) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, Km(k,i+1,j) )
       call CHECK( __LINE__, Km(k,i,j+1) )
       call CHECK( __LINE__, Km(k,i+1,j+1) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,XDIR) = - 0.125_RP & !
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ( Km  (k,i,j)+Km  (k,i+1,j)+Km  (k,i,j+1)+Km  (k,i+1,j+1) ) &
               * S12_Z(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z-x plane)
       !$omp parallel do
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, Km(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,YDIR) = DENS(k,i,j) * ( - 2.0_RP * Km(k,i,j) * S22_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       if ( ATMOS_PHY_TB_D1980_implicit ) then
          call calc_tend_MOMY( TEND, & ! (out)
                               qflx_sgs_momy, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          !$omp parallel do &
          !$omp private(ap)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS  ,i,j  )*Km(KS  ,i,j  ) &
                                   + DENS(KS+1,i,j  )*Km(KS+1,i,j  ) &
                                   + DENS(KS  ,i,j+1)*Km(KS  ,i,j+1) &
                                   + DENS(KS+1,i,j+1)*Km(KS+1,i,j+1) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_XVW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_XVZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS,i,j+1) )
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XVZ)
                ap = - dt * 0.25_RP * ( DENS(k  ,i,j  )*Km(k  ,i,j  ) &
                                      + DENS(k+1,i,j  )*Km(k+1,i,j  ) &
                                      + DENS(k  ,i,j+1)*Km(k  ,i,j+1) &
                                      + DENS(k+1,i,j+1)*Km(k+1,i,j+1) ) &
                                    * RFDZ(k) / GSQRT(k,i,j,I_XVW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XVZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) )
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_XVZ)
             b(KE,i,j) = - c(KE,i,j) + 0.5_RP * ( DENS(KE,i,j)+DENS(KE,i,j+1) )
          end do
          end do

          call MATRIX_SOLVER_tridiagonal( KA, KS, KE, IA, IIS, IIE, JA, JJS, JJE, &
                                          a(:,:,:), b(:,:,:), c(:,:,:), TEND(:,:,:), & ! (in)
                                          TEND2(:,:,:)                               ) ! (out)

          !$omp parallel do
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS, KE-1
                qflx_sgs_momy(k,i,j,ZDIR) = qflx_sgs_momy(k,i,j,ZDIR) &
                     - 0.25_RP * ( DENS(k  ,i,j  )*Km(k  ,i,j  ) &
                                 + DENS(k+1,i,j  )*Km(k+1,i,j  ) &
                                 + DENS(k  ,i,j+1)*Km(k  ,i,j+1) &
                                 + DENS(k+1,i,j+1)*Km(k+1,i,j+1) ) &
                     * dt * ( TEND2(k+1,i,j) - TEND2(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_XVW)
             end do

          end do
          end do

       end if

       !##### Thermodynamic Equation #####

       if ( ATMOS_PHY_TB_D1980_implicit ) then

          !$omp parallel do &
          !$omp private(ap)
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
            .false., &
            ATMOS_PHY_TB_D1980_implicit, &
            a, b, c, dt, &
            IIS, IIE, JJS, JJE )

    enddo
    enddo


    !##### Tracers #####
    do iq = 1, QA

       if ( iq == I_TKE ) then
          qflx_sgs_rhoq(:,:,:,:,iq) = 0.0_RP
          cycle
       end if
       if ( .not. TRACER_ADVC(iq) ) cycle

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          call calc_flux_phi( &
               qflx_sgs_rhoq(:,:,:,:,iq), &
               DENS, QTRC(:,:,:,iq), Kh, 1.0_RP, &
               GSQRT, J13G, J23G, J33G, MAPF, &
               .false., &
               ATMOS_PHY_TB_D1980_implicit, &
               a, b, c, dt, &
               IIS, IIE, JJS, JJE )

       enddo
       enddo
#ifdef DEBUG
       IIS = IUNDEF; IIE = IUNDEF; JJS = IUNDEF; JJE = IUNDEF
#endif

    enddo ! scalar quantities loop
#ifdef DEBUG
    iq = IUNDEF
#endif

    ! TKE
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       if ( ATMOS_PHY_TB_D1980_implicit ) then

          !$omp parallel do &
          !$omp private(ap)
          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS,i,j)+DENS(KS+1,i,j) ) &
                                 * 2.0_RP * ( Km(KS,i,j)+Km(KS+1,i,j) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + DENS(KS,i,j)
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
                ap = - dt * 0.25_RP * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                                    * 2.0_RP * ( Km(k,i,j)+Km(k+1,i,j) ) &
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
            qflx_tke,                            & ! (out)
            DENS, QTRC(:,:,:,I_TKE), Km, 2.0_RP, & ! (in)
            GSQRT, J13G, J23G, J33G, MAPF,       & ! (in)
            .false.,                             & ! (in)
            ATMOS_PHY_TB_D1980_implicit,         & ! (in)
            a, b, c, dt,                         & ! (in)
            IIS, IIE, JJS, JJE                   ) ! (in)

       call calc_tend_phi( RHOQ_t(:,:,:,I_TKE),           & ! (out)
                           qflx_tke,                      & ! (in)
                           GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                           IIS, IIE, JJS, JJE             ) ! (in)
       !$omp parallel do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOQ_t(k,i,j,I_TKE) = max( &
                 RHOQ_t(k,i,j,I_TKE) &
               + ( Km(k,i,j) * S2(k,i,j) &
                 - Kh(k,i,j) * N2(k,i,j) &
                 - ( C1 + C2*l(k,i,j)/delta(k,i,j) ) * QTRC(k,i,j,I_TKE)**(1.5_RP) / l(k,i,j) ) * DENS(k,i,j), &
                 - QTRC(k,i,j,I_TKE) * DENS(k,i,j) / real(dt,kind=RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    end do
    end do


    return
  end subroutine ATMOS_PHY_TB_D1980

end module scale_atmos_phy_tb_d1980
