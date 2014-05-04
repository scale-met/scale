!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          Runge-Kutta for Atmospheric dynamical process
!!          HEVI,
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-04-15 (S.Nishizawa) [new] newly impremented
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
#define HIVI_BICGSTAB 1
module scale_atmos_dyn_rk_hivi
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

#ifdef DEBUG
  use scale_debug, only: &
       CHECK
  use scale_const, only: &
       UNDEF  => CONST_UNDEF,  &
       IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_rk_hivi_setup
  public :: ATMOS_DYN_rk_hivi

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
  integer,  private :: ITMAX
  real(RP), private :: epsilon

#ifdef DRY
  real(RP), private :: kappa
#endif
  integer, private :: mtype ! MPI DATATYPE
  !-----------------------------------------------------------------------------


contains

  subroutine ATMOS_DYN_rk_hivi_setup( ATMOS_TYPE_DYN )
    use scale_process, only: &
       PRC_MPIstop
#ifdef DRY
    use scale_const, only: &
       CVdry  => CONST_CVdry,  &
       CPdry  => CONST_CPdry
#endif
    implicit none

    character(len=H_SHORT), intent(in) :: ATMOS_TYPE_DYN

    integer :: ierr

    namelist / PARAM_ATMOS_DYN_RK_HIVI / &
         ITMAX, &
         EPSILON
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HIVI'

    if ( ATMOS_TYPE_DYN .ne. 'HIVI' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_DYN is not HIVI. Check!'
       call PRC_MPIstop
    end if

#ifdef HIVI_BICGSTAB
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING Bi-CGSTAB'
#else
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING Multi-Grid'
    if ( IO_L ) write(IO_FID_LOG,*) 'xxx Not Implemented yet'
    call PRC_MPIstop
#endif

    ITMAX = 100
    epsilon = 0.1_RP ** (RP*2)
    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_RK_HIVI,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_DYN_RK_HIVI. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_DYN_RK_HIVI)


#ifdef DRY
    kappa = CPdry / CVdry
#endif

    if ( RP == DP ) then
       mtype = MPI_DOUBLE_PRECISION
    else if ( RP == SP ) then
       mtype = MPI_REAL
    else
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx Unsupported precision'
       call PRC_MPIstop
    end if

    return
  end subroutine ATMOS_DYN_rk_hivi_setup


  subroutine ATMOS_DYN_rk_hivi( &
    DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
    mflx_hi,                                     &
    DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
    DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
    DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
    Rtot, CVtot, CORIOLI,                        &
    num_diff, divdmp_coef,                       &
    FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, &
    CDZ, FDZ, FDX, FDY,                          &
    RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
    PHI, GSQRT, J13G, J23G, J33G,                &
    REF_pres, REF_dens,                          &
    dtrk                                         )
    use scale_const, only: &
#ifdef DRY
       Rdry   => CONST_Rdry,   &
       CVdry  => CONST_CVdry,  &
       CPdry  => CONST_CPdry,  &
#endif
       GRAV   => CONST_GRAV,   &
       P00    => CONST_PRE00
    use scale_comm, only: &
#ifdef _USE_RDMA
       COMM_rdma_vars8, &
#endif
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_common, only: &
       FACT_N, &
       FACT_F, &
       ATMOS_DYN_fct
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ
    implicit none

    real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(out) :: mflx_hi(KA,IA,JA,3) ! rho * vel(x,y,z)

    real(RP), intent(in),target :: DENS0(KA,IA,JA) ! prognostic variables
    real(RP), intent(in),target :: MOMZ0(KA,IA,JA) ! at previous dynamical time step
    real(RP), intent(in),target :: MOMX0(KA,IA,JA) !
    real(RP), intent(in),target :: MOMY0(KA,IA,JA) !
    real(RP), intent(in),target :: RHOT0(KA,IA,JA) !

    real(RP), intent(in) :: DENS(KA,IA,JA)   ! prognostic variables
    real(RP), intent(in) :: MOMZ(KA,IA,JA)   ! at previous RK step
    real(RP), intent(in) :: MOMX(KA,IA,JA)   !
    real(RP), intent(in) :: MOMY(KA,IA,JA)   !
    real(RP), intent(in) :: RHOT(KA,IA,JA)   !

    real(RP), intent(in) :: DENS_t(KA,IA,JA)
    real(RP), intent(in) :: MOMZ_t(KA,IA,JA)
    real(RP), intent(in) :: MOMX_t(KA,IA,JA)
    real(RP), intent(in) :: MOMY_t(KA,IA,JA)
    real(RP), intent(in) :: RHOT_t(KA,IA,JA)

    real(RP), intent(in) :: Rtot(KA,IA,JA) ! R for dry air + vapor
    real(RP), intent(in) :: CVtot(KA,IA,JA) ! CV
    real(RP), intent(in) :: CORIOLI(1,IA,JA)
    real(RP), intent(in) :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in) :: divdmp_coef

    logical,  intent(in) :: FLAG_FCT_RHO
    logical,  intent(in) :: FLAG_FCT_MOMENTUM
    logical,  intent(in) :: FLAG_FCT_T

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: FDX(IA-1)
    real(RP), intent(in) :: FDY(JA-1)
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: RFDZ(KA-1)
    real(RP), intent(in) :: RFDX(IA-1)
    real(RP), intent(in) :: RFDY(JA-1)

    real(RP), intent(in)  :: PHI     (KA,IA,JA)   !< geopotential
    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: REF_pres(KA,IA,JA)   !< reference pressure
    real(RP), intent(in)  :: REF_dens(KA,IA,JA)   !< reference density

    real(RP), intent(in)  :: dtrk


    ! diagnostic variables (work space)
    real(RP) :: PRES(KA,IA,JA) ! pressure [Pa]
    real(RP) :: VELZ(KA,IA,JA) ! velocity w [m/s]
    real(RP) :: VELX(KA,IA,JA) ! velocity u [m/s]
    real(RP) :: VELY(KA,IA,JA) ! velocity v [m/s]
    real(RP) :: POTT(KA,IA,JA) ! potential temperature [K]
    real(RP) :: DDIV(KA,IA,JA) ! divergence
    real(RP) :: DPRES(KA,IA,JA) ! pressure deviation from reference
    real(RP) :: DDENS(KA,IA,JA) ! density deviation from reference

    real(RP) :: Sr(KA,IA,JA)
    real(RP) :: Sw(KA,IA,JA)
    real(RP) :: Su(KA,IA,JA)
    real(RP) :: Sv(KA,IA,JA)
    real(RP) :: St(KA,IA,JA)
    real(RP) :: RCs2(KA,IA,JA)
    real(RP) :: B(KA,IA,JA)

#ifndef DRY
    real(RP) :: kappa(KA,IA,JA)
#endif

    real(RP) :: duvw

    real(RP) :: qflx_hi(KA,IA,JA,3)
    real(RP) :: qflx_J (KA,IA,JA)
    real(RP) :: mflx_hi2(KA,IA,JA,3)

    ! for implicit solver
    real(RP) :: M(7,KA,IA,JA)
    real(RP) :: r(KA,IA,JA)
    real(RP) :: p(KA,IA,JA)

    real(RP) :: rdt

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j

#ifdef DEBUG
    PRES(:,:,:) = UNDEF
    VELZ(:,:,:) = UNDEF
    VELX(:,:,:) = UNDEF
    VELY(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF
    DDIV(:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    mflx_hi(KS-1,:,:,ZDIR) = 0.0_RP

    qflx_hi(:,:,:,:) = UNDEF
    qflx_J (:,:,:)   = UNDEF
    mflx_hi2(:,:,:,:)   = UNDEF

    Sr(:,:,:) = UNDEF
    Sw(:,:,:) = UNDEF
    Su(:,:,:) = UNDEF
    Sv(:,:,:) = UNDEF
    St(:,:,:) = UNDEF
    RCs2(:,:,:) = UNDEF

    B(:,:,:) = UNDEF

    r(:,:,:) = UNDEF
    p(:,:,:) = UNDEF

    kappa(:,:,:) = UNDEF
#endif

    rdt = 1.0_RP / dtrk

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+2
       do i = IIS-1, IIE+2
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k,i,j) )
          call CHECK( __LINE__, DENS(k  ,i,j) )
          call CHECK( __LINE__, DENS(k+1,i,j) )
#endif
          VELZ(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+2
       do i = IIS-2, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMX(k,i,j) )
          call CHECK( __LINE__, DENS(k,i  ,j) )
          call CHECK( __LINE__, DENS(k,i+1,j) )
#endif
          VELX(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+1
       do i = IIS-1, IIE+2
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMY(k,i,j) )
          call CHECK( __LINE__, DENS(k,i,j  ) )
          call CHECK( __LINE__, DENS(k,i,j+1) )
#endif
          VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


#ifndef DRY
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          kappa(k,i,j) = ( CVtot(k,i,j) + Rtot(k,i,j) ) / CVtot(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, RHOT(k,i,j) )
          call CHECK( __LINE__, Rtot(k,i,j) )
#ifndef DRY
          call CHECK( __LINE__, kappa(k,i,j) )
#endif
#endif
#ifdef DRY
          PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rdry / P00 )**kappa
#else
          PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rtot(k,i,j) / P00 )**kappa(k,i,j)
#endif
       enddo
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
          call CHECK( __LINE__, PRES(KS+1,i,j) )
          call CHECK( __LINE__, DENS(KS,i,j) )
          call CHECK( __LINE__, PHI(KS-1,i,j) )
          call CHECK( __LINE__, PHI(KS+1,i,j) )
          call CHECK( __LINE__, PRES(KE-1,i,j) )
          call CHECK( __LINE__, DENS(KE,i,j) )
          call CHECK( __LINE__, PHI(KE+1,i,j) )
          call CHECK( __LINE__, PHI(KE-1,i,j) )
#endif
          PRES(KS-1,i,j) = PRES(KS+1,i,j) - DENS(KS,i,j) * ( PHI(KS-1,i,j) - PHI(KS+1,i,j) )
          PRES(KE+1,i,j) = PRES(KE-1,i,j) - DENS(KE,i,j) * ( PHI(KE+1,i,j) - PHI(KE-1,i,j) )
       end do
       end do
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS-1, KE+1
#ifdef DEBUG
          call CHECK( __LINE__, PRES(k,i,j) )
          call CHECK( __LINE__, REF_PRES(k,i,j) )
#endif
          DPRES(k,i,j) = PRES(k,i,j) - REF_pres(k,i,j)
       end do
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS(k,i,j) )
          call CHECK( __LINE__, REF_dens(k,i,j) )
#endif
          DDENS(k,i,j) = DENS(k,i,j) - REF_dens(k,i,j)
       end do
       enddo
       enddo
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
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
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! 3D divergence for damping
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k  ,i  ,j  ) )
          call CHECK( __LINE__, MOMZ(k-1,i  ,j  ) )
          call CHECK( __LINE__, MOMX(k  ,i  ,j  ) )
          call CHECK( __LINE__, MOMX(k  ,i-1,j  ) )
          call CHECK( __LINE__, MOMY(k  ,i  ,j  ) )
          call CHECK( __LINE__, MOMY(k  ,i  ,j-1) )
#endif
            DDIV(k,i,j) = ( MOMZ(k,i,j) - MOMZ(k-1,i  ,j  ) ) * RCDZ(k) &
                        + ( MOMX(k,i,j) - MOMX(k  ,i-1,j  ) ) * RCDX(i) &
                        + ( MOMY(k,i,j) - MOMY(k  ,i  ,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KS  ,i  ,j  ) )
          call CHECK( __LINE__, MOMX(KS  ,i  ,j  ) )
          call CHECK( __LINE__, MOMX(KS  ,i-1,j  ) )
          call CHECK( __LINE__, MOMY(KS  ,i  ,j  ) )
          call CHECK( __LINE__, MOMY(KS  ,i  ,j-1) )
          call CHECK( __LINE__, MOMZ(KE-1,i  ,j  ) )
          call CHECK( __LINE__, MOMX(KE  ,i  ,j  ) )
          call CHECK( __LINE__, MOMX(KE  ,i-1,j  ) )
          call CHECK( __LINE__, MOMY(KE  ,i  ,j  ) )
          call CHECK( __LINE__, MOMY(KE  ,i  ,j-1) )
#endif
            DDIV(KS,i,j) = ( MOMZ(KS,i,j)                      ) * RCDZ(KS) &
                         + ( MOMX(KS,i,j) - MOMX(KS  ,i-1,j  ) ) * RCDX(i) &
                         + ( MOMY(KS,i,j) - MOMY(KS  ,i  ,j-1) ) * RCDY(j)
            DDIV(KE,i,j) = (              - MOMZ(KE-1,i  ,j  ) ) * RCDZ(KE) &
                         + ( MOMX(KE,i,j) - MOMX(KE  ,i-1,j  ) ) * RCDX(i) &
                         + ( MOMY(KE,i,j) - MOMY(KE  ,i  ,j-1) ) * RCDY(j)
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### continuity equation #####

       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMX(k+1,i  ,j) )
          call CHECK( __LINE__, MOMX(k+1,i-1,j) )
          call CHECK( __LINE__, MOMX(k  ,i  ,j) )
          call CHECK( __LINE__, MOMX(k  ,i+1,j) )
          call CHECK( __LINE__, MOMY(k+1,i,j) )
          call CHECK( __LINE__, MOMY(k+1,i,j-1) )
          call CHECK( __LINE__, MOMY(k  ,i,j) )
          call CHECK( __LINE__, MOMY(k  ,i,j-1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,ZDIR) )
          call CHECK( __LINE__, GSQRT(k,i,j,I_XYW) )
#endif
          mflx_hi2(k,i,j,ZDIR) = J13G(k,i,j,I_XYW) * 0.25_RP * ( MOMX(k+1,i,j)+MOMX(k+1,i-1,j) &
                                                               + MOMX(k  ,i,j)+MOMX(k  ,i-1,j) ) &
                               + J23G(k,i,j,I_XYW) * 0.25_RP * ( MOMY(k+1,i,j)+MOMY(k+1,i,j-1) &
                                                               + MOMY(k  ,i,j)+MOMY(k  ,i,j-1) ) &
                               + GSQRT(k,i,j,I_XYW) * num_diff(k,i,j,I_DENS,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          mflx_hi2(KS-1,i,j,ZDIR) = 0.0_RP
          mflx_hi2(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo

       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,XDIR) )
#endif
          mflx_hi2(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) * num_diff(k,i,j,I_DENS,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
          mflx_hi2(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) * num_diff(k,i,j,I_DENS,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (z) #####

       ! at (x, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k  ,i,j) )
          call CHECK( __LINE__, VELZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k-2,i,j) )
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMZ(k+1,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,ZDIR) )
          call CHECK( __LINE__, GSQRT(k,i,j,I_XYZ) )
#endif
          qflx_hi(k,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(k,i,j)+VELZ(k-1,i,j) ) &
                              * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                                + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) ) &
                              + GSQRT(k,i,j,I_XYZ) * num_diff(k,i,j,I_MOMZ,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, MOMZ(KS+1,i,j) )
          call CHECK( __LINE__, num_diff(KS+1,i,j,I_MOMZ,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, MOMZ(KE-2,i,j) )
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMZ,ZDIR) )
#endif
          ! k = KS
          qflx_hi(KS,i,j,ZDIR) = 0.0_RP
          ! k = KS+1
          qflx_hi(KS+1,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(KS+1,i,j)+VELZ(KS,i,j) ) &
                                        * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS,i,j) ) &
                                          + FACT_F * ( MOMZ(KS+2,i,j)            ) ) &
                                 + GSQRT(KS+1,i,j,I_XYZ) * num_diff(KS+1,i,j,I_MOMZ,ZDIR)
          ! k = KE-1
          qflx_hi(KE-1,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(KE-1,i,j)+VELZ(KE-2,i,j) ) &
                                        * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                                          + FACT_F * (                MOMZ(KE-3,i,j) ) ) &
                                 + GSQRT(KE-1,i,j,I_XYZ) * num_diff(KE-1,i,j,I_MOMZ,ZDIR)
          ! k = KE
          qflx_hi(KE,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k  ,i,j) )
          call CHECK( __LINE__, VELZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k-2,i,j) )
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMZ(k+1,i,j) )
#endif
          qflx_J(k,i,j) = 0.5_RP &
                        * ( J13G(k,i,j,I_XYZ) * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                          + J23G(k,i,j,I_XYZ) * ( VELY(k,i,j)+VELY(k,i,j-1) ) ) &
                        * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                          + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j) )
          call CHECK( __LINE__, VELZ(KS+1,i,j) )
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, MOMZ(KS+1,i,j) )
          call CHECK( __LINE__, VELZ(KE-2,i,j) )
          call CHECK( __LINE__, VELZ(KE-1,i,j) )
          call CHECK( __LINE__, MOMZ(KE-2,i,j) )
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
#endif
          ! k = KS
          qflx_J(KS,i,j) = 0.0_RP
          ! k = KS+1
          qflx_J(KS+1,i,j) = 0.5_RP &
                           * ( J13G(KS+1,i,j,I_XYZ) * ( VELX(KS+1,i,j)+VELX(KS+1,i-1,j) ) &
                             + J23G(KS+1,i,j,I_XYZ) * ( VELY(KS+1,i,j)+VELY(KS+1,i,j-1) ) ) &
                           * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS  ,i,j) ) &
                             + FACT_F * ( MOMZ(KS+2,i,j)                ) )
          ! k = KE-1
          qflx_J(KE-1,i,j) = 0.5_RP &
                           * ( J13G(KE-1,i,j,I_XYZ) * ( VELX(KE-1,i,j)+VELX(KE-1,i-1,j) ) &
                             + J23G(KE-1,i,j,I_XYZ) * ( VELY(KE-1,i,j)+VELY(KE-1,i,j-1) ) ) &
                           * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                             + FACT_F * (                MOMZ(KE-3,i,j) ) )
          ! k = KE
          qflx_J(KE,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (u, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k  ,i,j) )
          call CHECK( __LINE__, VELX(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k,i-1,j) )
          call CHECK( __LINE__, MOMZ(k,i  ,j) )
          call CHECK( __LINE__, MOMZ(k,i+1,j) )
          call CHECK( __LINE__, MOMZ(k,i+2,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYW) &
                              * ( 0.5_RP * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                                * ( FACT_N * ( MOMZ(k,i+1,j)+MOMZ(k,i  ,j) ) &
                                  + FACT_F * ( MOMZ(k,i+2,j)+MOMZ(k,i-1,j) ) ) &
                                + num_diff(k,i,j,I_MOMZ,XDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (x, v, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k  ,i,j) )
          call CHECK( __LINE__, VELY(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k,i,j-1) )
          call CHECK( __LINE__, MOMZ(k,i,j  ) )
          call CHECK( __LINE__, MOMZ(k,i,j+1) )
          call CHECK( __LINE__, MOMZ(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVW) &
                              * ( 0.5_RP * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                                * ( FACT_N * ( MOMZ(k,i,j+1)+MOMZ(k,i,j  ) ) &
                                  + FACT_F * ( MOMZ(k,i,j+2)+MOMZ(k,i,j-1) ) ) &
                                + num_diff(k,i,j,I_MOMZ,YDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k+1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_J (k+1,i  ,j)        )
          call CHECK( __LINE__, qflx_J (k  ,i  ,j)        )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DDIV(k  ,i,j) )
          call CHECK( __LINE__, DDIV(k+1,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, MOMZ_t(k,i,j) )
#endif
          Sw(k,i,j) = &
               - ( ( qflx_hi(k+1,i,j,ZDIR) - qflx_hi(k,i  ,j  ,ZDIR) ) * RFDZ(k) &
                 + ( qflx_J (k+1,i,j)      - qflx_J (k,i  ,j  )      ) * RFDZ(k) &
                 + ( qflx_hi(k  ,i,j,XDIR) - qflx_hi(k,i-1,j  ,XDIR) ) * RCDX(i) &
                 + ( qflx_hi(k  ,i,j,YDIR) - qflx_hi(k,i  ,j-1,YDIR) ) * RCDY(j) &
                 ) / GSQRT(k,i,j,I_XYW) &
               + divdmp_coef * dtrk  * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) & ! divergence damping
               + MOMZ_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          Sw(KS-1,i,j) = 0.0_RP
          Sw(KE  ,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


       !##### momentum equation (x) #####

       ! at (u, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k,i  ,j) )
          call CHECK( __LINE__, VELZ(k,i+1,j) )
          call CHECK( __LINE__, DENS(k+1,i+1,j) )
          call CHECK( __LINE__, DENS(k+1,i  ,j) )
          call CHECK( __LINE__, DENS(k  ,i+1,j) )
          call CHECK( __LINE__, DENS(k  ,i  ,j) )
          call CHECK( __LINE__, MOMX(k-1,i,j) )
          call CHECK( __LINE__, MOMX(k  ,i,j) )
          call CHECK( __LINE__, MOMX(k+1,i,j) )
          call CHECK( __LINE__, MOMX(k+2,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,ZDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(k,i+1,j) + VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                                + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) ) &
                              + GSQRT(k,i,j,I_UYW) * num_diff(k,i,j,I_MOMX,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS,i  ,j) )
          call CHECK( __LINE__, VELZ(KS,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i  ,j) )
          call CHECK( __LINE__, DENS(KS  ,i+1,j) )
          call CHECK( __LINE__, DENS(KS  ,i  ,j) )
          call CHECK( __LINE__, MOMX(KS+1,i,j) )
          call CHECK( __LINE__, MOMX(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_MOMX,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-1,i  ,j) )
          call CHECK( __LINE__, VELZ(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i  ,j) )
          call CHECK( __LINE__, DENS(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE-1,i  ,j) )
          call CHECK( __LINE__, MOMX(KE-1,i,j) )
          call CHECK( __LINE__, MOMX(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMX,ZDIR) )
#endif
          ! at the bottom boundary
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          ! just above the bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = J33G * 0.25_RP * ( VELZ(KS,i+1,j) + VELZ(KS,i,j) ) &
                                 * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) ) &
                                 + GSQRT(KS,i,j,I_UYW) * num_diff(KS  ,i,j,I_MOMX,ZDIR)
          ! just below the top boundary
          qflx_hi(KE-1,i,j,ZDIR) = J33G * 0.25_RP * ( VELZ(KE-1,i+1,j) + VELZ(KE-1,i,j) ) &
                                 * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) ) &
                                 + GSQRT(KE-1,i,j,I_UYW) * num_diff(KE-1,i,j,I_MOMX,ZDIR)
          ! at the top boundary
          qflx_hi(KE,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
          qflx_J(k,i,j) = &
               ( J13G(k,i,j,I_UYW) * 0.5_RP   * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
               + J23G(k,i,j,I_UYW) * 0.125_RP * ( VELY(k+1,i+1,j  )+VELY(k,i+1,j  ) &
                                                + VELY(k+1,i  ,j  )+VELY(k,i  ,j  ) &
                                                + VELY(k+1,i+1,j-1)+VELY(k,i+1,j-1) &
                                                + VELY(k+1,i  ,j-1)+VELY(k,i  ,j-1) ) ) &
               * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                 + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) )

       enddo
       enddo
       enddo
       do j = JJS  , JJE
       do i = IIS, IIE
          qflx_J(KS-1,i,j) = 0.0_RP
          qflx_J(KS  ,i,j) = &
               ( J13G(KS,i,j,I_UYW) * 0.5_RP   * ( VELX(KS+1,i,j)+VELX(KS,i,j) ) &
               + J23G(KS,i,j,I_UYW) * 0.125_RP * ( VELY(KS+1,i+1,j  )+VELY(KS,i+1,j  ) &
                                                 + VELY(KS+1,i  ,j  )+VELY(KS,i  ,j  ) &
                                                 + VELY(KS+1,i+1,j-1)+VELY(KS,i+1,j-1) &
                                                 + VELY(KS+1,i  ,j-1)+VELY(KS,i  ,j-1) ) ) &
               * 0.5_RP * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) )
          qflx_J(KE-1,i,j) = &
               ( J13G(KE-1,i,j,I_UYW) * 0.5_RP    * ( VELX(KE,i,j)+VELX(KE-1,i,j) ) &
               + J23G(KE-1,i,j,I_UYW) * 0.125_RP * ( VELY(KE,i+1,j  )+VELY(KE-1,i+1,j  ) &
                                                   + VELY(KE,i  ,j  )+VELY(KE-1,i  ,j  ) &
                                                   + VELY(KE,i+1,j-1)+VELY(KE-1,i+1,j-1) &
                                                   + VELY(KE,i  ,j-1)+VELY(KE-1,i  ,j-1) ) ) &
               * 0.5_RP * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) )
          qflx_J(KE  ,i,j) = 0.0_RP
       enddo
       enddo

       ! at (x, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k,i  ,j) )
          call CHECK( __LINE__, VELX(k,i-1,j) )
          call CHECK( __LINE__, MOMX(k,i-2,j) )
          call CHECK( __LINE__, MOMX(k,i-1,j) )
          call CHECK( __LINE__, MOMX(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_XYZ) &
                              * ( 0.5_RP * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                                * ( FACT_N * ( MOMX(k,i  ,j)+MOMX(k,i-1,j) ) &
                                  + FACT_F * ( MOMX(k,i+1,j)+MOMX(k,i-2,j) ) ) &
                                + num_diff(k,i,j,I_MOMX,XDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k,i+1,j) )
          call CHECK( __LINE__, VELY(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i,j-1) )
          call CHECK( __LINE__, MOMX(k,i,j  ) )
          call CHECK( __LINE__, MOMX(k,i,j+1) )
          call CHECK( __LINE__, MOMX(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_UVZ) &
                              * ( 0.5_RP * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                                * ( FACT_N * ( MOMX(k,i,j+1)+MOMX(k,i,j  ) ) &
                                  + FACT_F * ( MOMX(k,i,j+2)+MOMX(k,i,j-1) ) ) &
                                + num_diff(k,i,j,I_MOMX,YDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(x)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i+1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DPRES(k,i+1,j) )
          call CHECK( __LINE__, DPRES(k,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i+1,j) )
          call CHECK( __LINE__, VELY(k,i  ,j  ) )
          call CHECK( __LINE__, VELY(k,i+1,j  ) )
          call CHECK( __LINE__, VELY(k,i  ,j-1) )
          call CHECK( __LINE__, VELY(k,i+1,j-1) )
          call CHECK( __LINE__, DDIV(k,i+1,j) )
          call CHECK( __LINE__, DDIV(k,i  ,j) )
          call CHECK( __LINE__, MOMX0(k,i,j) )
#endif
          Su(k,i,j) = ( &
               - ( ( qflx_hi(k,i  ,j,ZDIR) - qflx_hi(k-1,i,j  ,ZDIR) ) * RCDZ(k) &
                 + ( qflx_J (k,i  ,j)      - qflx_J (k-1,i,j)        ) * RCDZ(k) &
                 + ( qflx_hi(k,i+1,j,XDIR) - qflx_hi(k  ,i,j  ,XDIR) ) * RFDX(i) &
                 + ( qflx_hi(k,i  ,j,YDIR) - qflx_hi(k  ,i,j-1,YDIR) ) * RCDY(j) ) &
               - ( J13G(k+1,i,j,I_UYZ) * ( DPRES(k+1,i+1,j)+DPRES(k+1,i,j) ) &
                 - J13G(k-1,i,j,I_UYZ) * ( DPRES(k-1,i+1,j)+DPRES(k-1,i,j) ) ) &
                 * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) & 
               ) / GSQRT(k,i,j,I_UYZ) &
               + 0.0625_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i+1,j) ) &
                           * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                           * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) & ! coriolis force
               + divdmp_coef * dtrk * ( DDIV(k,i+1,j)-DDIV(k,i,j) ) * FDX(i) & ! divergence damping
               + MOMX_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (y) #####

       ! at (x, v, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS  , IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k,i,j+1) )
          call CHECK( __LINE__, VELZ(k,i,j  ) )
          call CHECK( __LINE__, DENS(k+1,i+1,j) )
          call CHECK( __LINE__, DENS(k+1,i  ,j) )
          call CHECK( __LINE__, DENS(k  ,i+1,j) )
          call CHECK( __LINE__, DENS(k  ,i  ,j) )
          call CHECK( __LINE__, MOMY(k-1,i,j) )
          call CHECK( __LINE__, MOMY(k  ,i,j) )
          call CHECK( __LINE__, MOMY(k+1,i,j) )
          call CHECK( __LINE__, MOMY(k+2,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,ZDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = J33G * 0.5_RP * ( VELZ(k,i,j+1) + VELZ(k,i,j) ) &
                              * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                                + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) ) &
                              + GSQRT(k,i,j,I_XVW) * num_diff(k,i,j,I_MOMY,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k,vel) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS  , IIE
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(KS  ,i,j+1) )
          call CHECK( __LINE__, VELZ(KS  ,i,j  ) )
          call CHECK( __LINE__, DENS(KS+1,i+1,j) )
          call CHECK( __LINE__, DENS(KS+1,i  ,j) )
          call CHECK( __LINE__, DENS(KS  ,i+1,j) )
          call CHECK( __LINE__, DENS(KS  ,i  ,j) )
          call CHECK( __LINE__, MOMY(KS+1,i,j) )
          call CHECK( __LINE__, MOMY(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_MOMY,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-1,i,j+1) )
          call CHECK( __LINE__, VELZ(KE-1,i,j  ) )
          call CHECK( __LINE__, DENS(KE  ,i+1,j) )
          call CHECK( __LINE__, DENS(KE  ,i  ,j) )
          call CHECK( __LINE__, DENS(KE-1,i+1,j) )
          call CHECK( __LINE__, DENS(KE-1,i  ,j) )
          call CHECK( __LINE__, MOMY(KE-1,i,j) )
          call CHECK( __LINE__, MOMY(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMY,ZDIR) )
#endif
          ! at the bottom boundary
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          ! just above the bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = J33G * 0.25_RP &
                                 * ( VELZ(KS,i,j+1) + VELZ(KS,i,j) ) &
                                 * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) ) &
                                 + GSQRT(KS,i,j,I_XVW) * num_diff(KS  ,i,j,I_MOMY,ZDIR)
          ! just below the top boundary
          qflx_hi(KE-1,i,j,ZDIR) = J33G * 0.25_RP &
                                 * ( VELZ(KE-1,i,j+1) + VELZ(KE-1,i,j) ) &
                                 * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) ) &
                                 + GSQRT(KE-1,i,j,I_XVW) * num_diff(KE-1,i,j,I_MOMY,ZDIR)
          ! at the top boundary
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS  , IIE
       do k = KS+1, KE-2
          qflx_J(k,i,j) = &
               ( J13G(k,i,j,I_XVW) * 0.125_RP * ( VELX(k+1,i  ,j+1)+VELX(k,i  ,j+1) &
                                                + VELX(k+1,i-1,j+1)+VELX(k,i-1,j+1) &
                                                + VELX(k+1,i  ,j  )+VELX(k,i  ,j  ) &
                                                + VELX(k+1,i-1,j  )+VELX(k,i-1,j  ) ) &
               + J23G(k,i,j,I_XVW) * 0.5_RP   * ( VELY(k+1,i,j)+VELY(k,i,j) ) ) &
               * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                 + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) )
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J(KS-1,i,j) = 0.0_RP
          qflx_J(KS  ,i,j) = &
               ( J13G(KS,i,j,I_XVW) * 0.125_RP * ( VELX(KS+1,i  ,j+1)+VELX(KS,i  ,j+1) &
                                                 + VELX(KS+1,i-1,j+1)+VELX(KS,i-1,j+1) &
                                                 + VELX(KS+1,i  ,j  )+VELX(KS,i  ,j  ) &
                                                 + VELX(KS+1,i-1,j  )+VELX(KS,i-1,j  ) ) &
               + J23G(KS,i,j,I_XVW) * 0.5_RP   * ( VELY(KS+1,i,j)+VELY(KS,i,j) ) ) &
               * 0.5_RP * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) )
          qflx_J(KE-1,i,j) = &
               ( J23G(KE-1,i,j,I_XVW) * 0.5_RP   * ( VELY(KE,i,j)+VELY(KE-1,i,j) ) &
               + J13G(KE-1,i,j,I_XVW) * 0.125_RP * ( VELX(KE,i  ,j+1)+VELX(KE-1,i  ,j+1) &
                                                   + VELX(KE,i-1,j+1)+VELX(KE-1,i-1,j+1) &
                                                   + VELX(KE,i  ,j  )+VELX(KE-1,i  ,j  ) &
                                                   + VELX(KE,i-1,j  )+VELX(KE-1,i-1,j  ) ) ) &
               * 0.5_RP * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) )
          qflx_J(KE  ,i,j) = 0.0_RP
       enddo
       enddo

       ! at (u, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELX(k,i,j+1) )
          call CHECK( __LINE__, VELX(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i-1,j) )
          call CHECK( __LINE__, MOMY(k,i  ,j) )
          call CHECK( __LINE__, MOMY(k,i+1,j) )
          call CHECK( __LINE__, MOMY(k,i+2,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UVZ) &
                              * ( 0.5_RP * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                                * ( FACT_N * ( MOMY(k,i+1,j)+MOMY(k,i  ,j) ) &
                                  + FACT_F * ( MOMY(k,i+2,j)+MOMY(k,i-1,j) ) ) &
                                + num_diff(k,i,j,I_MOMY,XDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE+1
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, VELY(k,i,j  ) )
          call CHECK( __LINE__, VELY(k,i,j-1) )
          call CHECK( __LINE__, MOMY(k,i,j-2) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i,j+1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMY,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XYZ) &
                              * ( 0.5_RP * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                                * ( FACT_N * ( MOMY(k,i,j  )+MOMY(k,i,j-1) ) &
                                  + FACT_F * ( MOMY(k,i,j+1)+MOMY(k,i,j-2) ) ) &
                                + num_diff(k,i,j,I_MOMY,YDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(y)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j+1,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j,YDIR) )
          call CHECK( __LINE__, DPRES(k,i,j  ) )
          call CHECK( __LINE__, DPRES(k,i,j+1) )
          call CHECK( __LINE__, CORIOLI(1,i,j  ) )
          call CHECK( __LINE__, CORIOLI(1,i,j+1) )
          call CHECK( __LINE__, VELX(k,i  ,j  ) )
          call CHECK( __LINE__, VELX(k,i  ,j+1) )
          call CHECK( __LINE__, VELX(k,i-1,j  ) )
          call CHECK( __LINE__, VELX(k,i-1,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j  ) )
          call CHECK( __LINE__, MOMY_t(k,i,j) )
#endif
          Sv(k,i,j) = &
               ( - ( ( qflx_hi(k,i,j  ,ZDIR) - qflx_hi(k-1,i  ,j,ZDIR) ) * RCDZ(k) &
                   + ( qflx_J (k,i,j  )      - qflx_J (k-1,i  ,j)      ) * RCDZ(k) &
                   + ( qflx_hi(k,i,j  ,XDIR) - qflx_hi(k  ,i-1,j,XDIR) ) * RCDX(i) &
                   + ( qflx_hi(k,i,j+1,YDIR) - qflx_hi(k  ,i  ,j,YDIR) ) * RFDY(j) ) &
                 - ( J23G(k+1,i,j,I_XVZ) * ( DPRES(k+1,i,j+1)+DPRES(k+1,i,j) ) &
                   - J23G(k-1,i,j,I_XVZ) * ( DPRES(k-1,i,j+1)+DPRES(k-1,i,j) ) ) &
                   * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) & ! pressure gradient force
               ) / GSQRT(k,i,j,I_XVZ) &
               - 0.0625_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) &
                           * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
                           * ( VELX(k,i,j+1)+VELX(k,i,j)+VELX(k,i-1,j+1)+VELX(k,i-1,j) ) & ! coriolis force
               + divdmp_coef * dtrk * ( DDIV(k,i,j+1)-DDIV(k,i,j) ) * FDY(j) & ! divergence damping
               + MOMY_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### Thermodynamic Equation #####


       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi2(k,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = mflx_hi2(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j) + POTT(k  ,i,j) ) &
                                + FACT_F * ( POTT(k+2,i,j) + POTT(k-1,i,j) ) ) &
                              + GSQRT(k,i,j,I_XYW) * num_diff(k,i,j,I_RHOT,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi2(KS  ,i,j,ZDIR) )
          call CHECK( __LINE__, mflx_hi2(KE-1,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, POTT(KS+1,i,j) )
          call CHECK( __LINE__, POTT(KE-1,i,j) )
          call CHECK( __LINE__, POTT(KE  ,i,j) )
          call CHECK( __LINE__, GSQRT(KS  ,i,j,I_XYW) )
          call CHECK( __LINE__, GSQRT(KE-1,i,j,I_XYW) )
          call CHECK( __LINE__, num_diff(KS  ,i,j,I_RHOT,ZDIR) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_RHOT,ZDIR) )
#endif

          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = mflx_hi2(KS  ,i,j,ZDIR) * 0.5_RP * ( POTT(KS+1,i,j) + POTT(KS  ,i,j) ) &
                                 + GSQRT(KS,i,j,I_XYW) * num_diff(KS  ,i,j,I_RHOT,ZDIR)
          qflx_hi(KE-1,i,j,ZDIR) = mflx_hi2(KE-1,i,j,ZDIR) * 0.5_RP * ( POTT(KE  ,i,j) + POTT(KE-1,i,j) ) &
                                 + GSQRT(KE-1,i,j,I_XYW) * num_diff(KE-1,i,j,I_RHOT,ZDIR)
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi2(k,i,j,XDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
#endif
          qflx_hi(k,i,j,XDIR) = mflx_hi2(k,i,j,XDIR) &
                                * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                  + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi2(k,i,j,YDIR) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j+2) )
#endif
          qflx_hi(k,i,j,YDIR) = mflx_hi2(k,i,j,YDIR) &
                                * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) ) &
                                  + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, RHOT_t(k,i,j) )
#endif
          St(k,i,j) = &
               - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                 + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                 + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) &
                 ) / GSQRT(k,i,j,I_XYZ) &
               + RHOT_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum flux #####

       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_XYW) )
          call CHECK( __LINE__, MOMZ(k,i,j) )
#endif
          mflx_hi(k,i,j,ZDIR) = J33G * MOMZ(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          mflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
          call CHECK( __LINE__, MOMX(k,i,j) )
#endif
          mflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) * MOMX(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
          call CHECK( __LINE__, MOMY(k,i,j) )
#endif
          mflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) * MOMY(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### Thermodynamic Equation #####

       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j) + POTT(k  ,i,j) ) &
                                + FACT_F * ( POTT(k+2,i,j) + POTT(k-1,i,j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = mflx_hi(KS  ,i,j,ZDIR) * 0.5_RP * ( POTT(KS+1,i,j) + POTT(KS  ,i,j) )
          qflx_hi(KE-1,i,j,ZDIR) = mflx_hi(KE-1,i,j,ZDIR) * 0.5_RP * ( POTT(KE  ,i,j) + POTT(KE-1,i,j) )
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) &
                                * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                  + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) &
                                * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) ) &
                                  + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT_RK(k,i,j) = RHOT0(k,i,j) &
               + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                            + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                            / GSQRT(k,i,j,I_XYZ) &
                          + St(k,i,j) )
       end do
       end do
       end do

       !##### continuous equation #####

       !--- update density
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS0(k,i,j) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, mflx_hi2(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi2(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi2(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, mflx_hi2(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DENS_t(k,i,j) )
#endif
          DENS_RK(k,i,j) = DENS0(k,i,j) &
               + dtrk * ( - ( ( mflx_hi (k,i,j,ZDIR)-mflx_hi (k-1,i  ,j  ,ZDIR) &
                              + mflx_hi2(k,i,j,ZDIR)-mflx_hi2(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                            + ( mflx_hi (k,i,j,XDIR)-mflx_hi (k  ,i-1,j  ,XDIR) &
                              + mflx_hi2(k,i,j,XDIR)-mflx_hi2(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( mflx_hi (k,i,j,YDIR)-mflx_hi (k  ,i  ,j-1,YDIR) &
                              + mflx_hi2(k,i,j,YDIR)-mflx_hi2(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                            / GSQRT(k,i,j,I_XYZ) & ! divergence
                          + DENS_t(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    ! implicit solver
    !$omp parallel do private(i,j,k,PT,A,B,C) OMP_SCHEDULE_ collapse(2)
#ifdef HIVI_BICGSTAB
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, RHOT(k,i,j) )
          call CHECK( __LINE__, PRES(k,i,j) )
#ifndef DRY
          call CHECK( __LINE__, kappa(k,i,j) )
#endif
#endif
          RCs2(k,i,j) = RHOT(k,i,j) / ( PRES(k,i,j) &
#ifdef DRY
               * kappa &
#else
               * kappa(k,i,j) &
#endif
               )
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    end do
    end do

    call COMM_vars8( Su, 1 )
    call COMM_vars8( Sv, 2 )
    call COMM_wait ( Su, 1 )
    call COMM_wait ( Sv, 2 )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
       ! r = b - M x0
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMX(k,i-1,j) )
          call CHECK( __LINE__, MOMX(k,i  ,j) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, POTT(k-2,i,j) )
          call CHECK( __LINE__, POTT(k-1,i,j) )
          call CHECK( __LINE__, POTT(k  ,i,j) )
          call CHECK( __LINE__, POTT(k+1,i,j) )
          call CHECK( __LINE__, POTT(k+2,i,j) )
          call CHECK( __LINE__, POTT(k,i-2,j) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+2,j) )
          call CHECK( __LINE__, POTT(k,i,j-2) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j+2) )
          call CHECK( __LINE__, Sw(k-1,i,j) )
          call CHECK( __LINE__, Sw(k  ,i,j) )
          call CHECK( __LINE__, Su(k,i-1,j) )
          call CHECK( __LINE__, Su(k,i  ,j) )
          call CHECK( __LINE__, Sv(k,i,j-1) )
          call CHECK( __LINE__, Sv(k,i,j  ) )
          call CHECK( __LINE__, GSQRT(k,i ,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(k,i-1,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(k,i,j-1,I_XVZ) )
          call CHECK( __LINE__, GSQRT(k,i,j  ,I_XVZ) )
          call CHECK( __LINE__, GSQRT(k,i,j,I_XYZ) )
          call CHECK( __LINE__, St(k,i,j) )
          call CHECK( __LINE__, DPRES(k-1,i,j) )
          call CHECK( __LINE__, DPRES(k  ,i,j) )
          call CHECK( __LINE__, DPRES(k+1,i,j) )
          call CHECK( __LINE__, RCs2(k-1,i,j) )
          call CHECK( __LINE__, RCs2(k  ,i,j) )
          call CHECK( __LINE__, RCs2(k+1,i,j) )
          call CHECK( __LINE__, RHOT(k-1,i,j) )
          call CHECK( __LINE__, RHOT(k+1,i,j) )
          call CHECK( __LINE__, DDENS(k-1,i,j) )
          call CHECK( __LINE__, DDENS(k+1,i,j) )
#endif
          B(k,i,j) = ( &
               ( J33G * ( MOMZ(k  ,i,j) + dtrk*Sw(k  ,i,j) ) &
                 * ( FACT_N*(POTT(k+1,i,j)+POTT(k  ,i,j)) &
                   + FACT_F*(POTT(k+2,i,j)+POTT(k-1,i,j)) ) &
               - J33G * ( MOMZ(k-1,i,j) + dtrk*Sw(k-1,i,j) ) &
                 * ( FACT_N*(POTT(k  ,i,j)+POTT(k-1,i,j)) &
                   + FACT_F*(POTT(k+1,i,j)+POTT(k-2,i,j)) ) ) * RCDZ(k) &
             + ( GSQRT(k,i  ,j,I_UYZ) * ( MOMX(k,i  ,j) + dtrk*Su(k,i  ,j) ) &
                 * ( FACT_N*(POTT(k,i+1,j)+POTT(k,i  ,j)) &
                   + FACT_F*(POTT(k,i+2,j)+POTT(k,i-1,j)) ) &
               - GSQRT(k,i-1,j,I_UYZ) * ( MOMX(k,i-1,j) + dtrk*Su(k,i-1,j) ) &
                 * ( FACT_N*(POTT(k,i  ,j)+POTT(k,i-1,j)) &
                   + FACT_F*(POTT(k,i+1,j)+POTT(k,i-2,j)) ) ) * RCDX(i) &
             + ( GSQRT(k,i,j  ,I_XVZ) * ( MOMY(k,i,j  ) + dtrk*Sv(k,i,j  ) ) &
                 * ( FACT_N*(POTT(k,i,j+1)+POTT(k,i,j  )) &
                   + FACT_F*(POTT(k,i,j+2)+POTT(k,i,j-1)) ) &
               - GSQRT(k,i,j-1,I_XVZ) * ( MOMY(k,i,j-1) + dtrk*Sv(k,i,j-1) ) &
                 * ( FACT_N*(POTT(k,i,j  )+POTT(k,i,j-1)) &
                   + FACT_F*(POTT(k,i,j+1)+POTT(k,i,j-2)) ) ) * RCDY(j) &
             + GSQRT(k,i,j,I_XYZ) * ( St(k,i,j) - DPRES(k,i,j) * RCs2(k,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( ( DPRES(k+1,i,j)*RCs2(k+1,i,j) &
                               - DPRES(k-1,i,j)*RCs2(k-1,i,j) ) &
                             - ( RHOT(k+1,i,j)*DDENS(k+1,i,j)/DENS(k+1,i,j) &
                               - RHOT(k-1,i,j)*DDENS(k-1,i,j)/DENS(k-1,i,j) ) &
                             - ( DENS(k+1,i,j)*(RHOT_RK(k+1,i,j)/DENS_RK(k+1,i,j) - POTT(k+1,i,j) ) &
                               - DENS(k-1,i,j)*(RHOT_RK(k-1,i,j)/DENS_RK(k-1,i,j) - POTT(k-1,i,j) ) ) &
                             ) / ( FDZ(k) + FDZ(k-1) )
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KS,i,j) )
          call CHECK( __LINE__, MOMX(KS,i,j) )
          call CHECK( __LINE__, MOMY(KS,i,j) )
          call CHECK( __LINE__, Sw(KS,i,j) )
          call CHECK( __LINE__, Su(KS,i-1,j) )
          call CHECK( __LINE__, Su(KS,i  ,j) )
          call CHECK( __LINE__, Sv(KS,i,j-1) )
          call CHECK( __LINE__, Sv(KS,i,j  ) )
          call CHECK( __LINE__, St(KS,i,j) )
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, POTT(KS+1,i,j) )
          call CHECK( __LINE__, POTT(KS,i-2,j) )
          call CHECK( __LINE__, POTT(KS,i-1,j) )
          call CHECK( __LINE__, POTT(KS,i  ,j) )
          call CHECK( __LINE__, POTT(KS,i+1,j) )
          call CHECK( __LINE__, POTT(KS,i+2,j) )
          call CHECK( __LINE__, POTT(KS,i,j-2) )
          call CHECK( __LINE__, POTT(KS,i,j-1) )
          call CHECK( __LINE__, POTT(KS,i,j  ) )
          call CHECK( __LINE__, POTT(KS,i,j+1) )
          call CHECK( __LINE__, POTT(KS,i,j+2) )
          call CHECK( __LINE__, GSQRT(KS,i,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(KS,i,j,I_XVZ) )
          call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
          call CHECK( __LINE__, DPRES(KS  ,i,j) )
          call CHECK( __LINE__, DPRES(KS+1,i,j) )
          call CHECK( __LINE__, RCs2(KS  ,i,j) )
          call CHECK( __LINE__, RCs2(KS+1,i,j) )
          call CHECK( __LINE__, RHOT(KS+1,i,j) )
          call CHECK( __LINE__, DDENS(KS+1,i,j) )
#endif
          B(KS,i,j) = ( &
               ( J33G * ( MOMZ(KS  ,i,j) + dtrk*Sw(KS  ,i,j) ) &
                 * 0.5_RP*(POTT(KS+1,i,j)+POTT(KS  ,i,j)) ) * RCDZ(KS) &
             + ( GSQRT(KS,i  ,j,I_UYZ) * ( MOMX(KS,i  ,j) + dtrk*Su(KS,i  ,j) ) &
                 * ( FACT_N*(POTT(KS,i+1,j)+POTT(KS,i  ,j)) &
                   + FACT_F*(POTT(KS,i+2,j)+POTT(KS,i-1,j)) ) &
               - GSQRT(KS,i-1,j,I_UYZ) * ( MOMX(KS,i-1,j) + dtrk*Su(KS,i-1,j) ) &
                 * ( FACT_N*(POTT(KS,i  ,j)+POTT(KS,i-1,j)) &
                   + FACT_F*(POTT(KS,i+1,j)+POTT(KS,i-2,j)) ) ) * RCDX(i) &
             + ( GSQRT(KS,i,j  ,I_XVZ) * ( MOMY(KS,i,j  ) + dtrk*Sv(KS,i,j  ) ) &
                 * ( FACT_N*(POTT(KS,i,j+1)+POTT(KS,i,j  )) &
                   + FACT_F*(POTT(KS,i,j+2)+POTT(KS,i,j-1)) ) &
               - GSQRT(KS,i,j-1,I_XVZ) * ( MOMY(KS,i,j-1) + dtrk*Sv(KS,i,j-1) ) &
                 * ( FACT_N*(POTT(KS,i,j  )+POTT(KS,i,j-1)) &
                   + FACT_F*(POTT(KS,i,j+1)+POTT(KS,i,j-2)) ) ) * RCDY(j) &
             + GSQRT(KS,i,j,I_XYZ) * ( St(KS,i,j) - DPRES(KS,i,j) * RCs2(KS,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( DPRES(KS+1,i,j)*RCs2(KS+1,i,j) &
                             - RHOT(KS+1,i,j)*DDENS(KS+1,i,j)/DENS(KS+1,i,j) &
                             - DENS(KS+1,i,j)*(RHOT_RK(KS+1,i,j)/DENS_RK(KS+1,i,j) - POTT(KS+1,i,j) ) ) &
                             / ( FDZ(KS) + FDZ(KS-1) )
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, MOMZ(KS+1  ,i,j) )
          call CHECK( __LINE__, MOMX(KS+1,i-1,j) )
          call CHECK( __LINE__, MOMX(KS+1,i  ,j) )
          call CHECK( __LINE__, MOMY(KS+1,i,j-1) )
          call CHECK( __LINE__, MOMY(KS+1,i,j  ) )
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, POTT(KS+1  ,i,j) )
          call CHECK( __LINE__, POTT(KS+1+1,i,j) )
          call CHECK( __LINE__, POTT(KS+1+2,i,j) )
          call CHECK( __LINE__, POTT(KS+1,i-2,j) )
          call CHECK( __LINE__, POTT(KS+1,i-1,j) )
          call CHECK( __LINE__, POTT(KS+1,i  ,j) )
          call CHECK( __LINE__, POTT(KS+1,i+1,j) )
          call CHECK( __LINE__, POTT(KS+1,i+2,j) )
          call CHECK( __LINE__, POTT(KS+1,i,j-2) )
          call CHECK( __LINE__, POTT(KS+1,i,j-1) )
          call CHECK( __LINE__, POTT(KS+1,i,j  ) )
          call CHECK( __LINE__, POTT(KS+1,i,j+1) )
          call CHECK( __LINE__, POTT(KS+1,i,j+2) )
          call CHECK( __LINE__, Sw(KS  ,i,j) )
          call CHECK( __LINE__, Sw(KS+1,i,j) )
          call CHECK( __LINE__, Su(KS+1,i-1,j) )
          call CHECK( __LINE__, Su(KS+1,i  ,j) )
          call CHECK( __LINE__, Sv(KS+1,i,j-1) )
          call CHECK( __LINE__, Sv(KS+1,i,j  ) )
          call CHECK( __LINE__, GSQRT(KS+1,i ,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(KS+1,i-1,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(KS+1,i,j-1,I_XVZ) )
          call CHECK( __LINE__, GSQRT(KS+1,i,j  ,I_XVZ) )
          call CHECK( __LINE__, GSQRT(KS+1,i,j,I_XYZ) )
          call CHECK( __LINE__, St(KS+1,i,j) )
          call CHECK( __LINE__, DPRES(KS  ,i,j) )
          call CHECK( __LINE__, DPRES(KS+1,i,j) )
          call CHECK( __LINE__, DPRES(KS+2,i,j) )
          call CHECK( __LINE__, RCs2(KS  ,i,j) )
          call CHECK( __LINE__, RCs2(KS+1,i,j) )
          call CHECK( __LINE__, RCs2(KS+2,i,j) )
          call CHECK( __LINE__, RHOT(KS  ,i,j) )
          call CHECK( __LINE__, RHOT(KS+2,i,j) )
          call CHECK( __LINE__, DDENS(KS  ,i,j) )
          call CHECK( __LINE__, DDENS(KS+2,i,j) )
#endif
          B(KS+1,i,j) = ( &
               ( J33G * ( MOMZ(KS+1,i,j) + dtrk*Sw(KS+1,i,j) ) &
                 * ( FACT_N*(POTT(KS+2,i,j)+POTT(KS+1,i,j)) &
                   + FACT_F*(POTT(KS+3,i,j)+POTT(KS  ,i,j)) ) &
               - J33G * ( MOMZ(KS+1-1,i,j) + dtrk*Sw(KS+1-1,i,j) ) &
                 * ( 0.5_RP*(POTT(KS+1,i,j)+POTT(KS,i,j)) ) ) * RCDZ(KS+1) &
             + ( GSQRT(KS+1,i  ,j,I_UYZ) * ( MOMX(KS+1,i  ,j) + dtrk*Su(KS+1,i  ,j) ) &
                 * ( FACT_N*(POTT(KS+1,i+1,j)+POTT(KS+1,i  ,j)) &
                   + FACT_F*(POTT(KS+1,i+2,j)+POTT(KS+1,i-1,j)) ) &
               - GSQRT(KS+1,i-1,j,I_UYZ) * ( MOMX(KS+1,i-1,j) + dtrk*Su(KS+1,i-1,j) ) &
                 * ( FACT_N*(POTT(KS+1,i  ,j)+POTT(KS+1,i-1,j)) &
                   + FACT_F*(POTT(KS+1,i+1,j)+POTT(KS+1,i-2,j)) ) ) * RCDX(i) &
             + ( GSQRT(KS+1,i,j  ,I_XVZ) * ( MOMY(KS+1,i,j  ) + dtrk*Sv(KS+1,i,j  ) ) &
                 * ( FACT_N*(POTT(KS+1,i,j+1)+POTT(KS+1,i,j  )) &
                   + FACT_F*(POTT(KS+1,i,j+2)+POTT(KS+1,i,j-1)) ) &
               - GSQRT(KS+1,i,j-1,I_XVZ) * ( MOMY(KS+1,i,j-1) + dtrk*Sv(KS+1,i,j-1) ) &
                 * ( FACT_N*(POTT(KS+1,i,j  )+POTT(KS+1,i,j-1)) &
                   + FACT_F*(POTT(KS+1,i,j+1)+POTT(KS+1,i,j-2)) ) ) * RCDY(j) &
             + GSQRT(KS+1,i,j,I_XYZ) * ( St(KS+1,i,j) - DPRES(KS+1,i,j) * RCs2(KS+1,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( ( DPRES(KS+2,i,j)*RCs2(KS+2,i,j) &
                               - DPRES(KS  ,i,j)*RCs2(KS  ,i,j) ) &
                             - ( RHOT(KS+2,i,j)*DDENS(KS+2,i,j)/DENS(KS+2,i,j) &
                               - RHOT(KS  ,i,j)*DDENS(KS  ,i,j)/DENS(KS  ,i,j) ) &
                             - ( DENS(KS+2,i,j)*(RHOT_RK(KS+2,i,j)/DENS_RK(KS+2,i,j) - POTT(KS+2,i,j) ) &
                               - DENS(KS  ,i,j)*(RHOT_RK(KS  ,i,j)/DENS_RK(KS  ,i,j) - POTT(KS  ,i,j) ) ) &
                             ) / ( FDZ(KS+1) + FDZ(KS) )
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KE-2,i,j) )
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
          call CHECK( __LINE__, MOMX(KE-1,i-1,j) )
          call CHECK( __LINE__, MOMX(KE-1,i  ,j) )
          call CHECK( __LINE__, MOMY(KE-1,i,j-1) )
          call CHECK( __LINE__, MOMY(KE-1,i,j  ) )
          call CHECK( __LINE__, POTT(KE-3,i,j) )
          call CHECK( __LINE__, POTT(KE-2,i,j) )
          call CHECK( __LINE__, POTT(KE-1,i,j) )
          call CHECK( __LINE__, POTT(KE  ,i,j) )
          call CHECK( __LINE__, POTT(KE-1,i-2,j) )
          call CHECK( __LINE__, POTT(KE-1,i-1,j) )
          call CHECK( __LINE__, POTT(KE-1,i  ,j) )
          call CHECK( __LINE__, POTT(KE-1,i+1,j) )
          call CHECK( __LINE__, POTT(KE-1,i+2,j) )
          call CHECK( __LINE__, POTT(KE-1,i,j-2) )
          call CHECK( __LINE__, POTT(KE-1,i,j-1) )
          call CHECK( __LINE__, POTT(KE-1,i,j  ) )
          call CHECK( __LINE__, POTT(KE-1,i,j+1) )
          call CHECK( __LINE__, POTT(KE-1,i,j+2) )
          call CHECK( __LINE__, Sw(KE-2,i,j) )
          call CHECK( __LINE__, Sw(KE-1,i,j) )
          call CHECK( __LINE__, Su(KE-1,i-1,j) )
          call CHECK( __LINE__, Su(KE-1,i  ,j) )
          call CHECK( __LINE__, Sv(KE-1,i,j-1) )
          call CHECK( __LINE__, Sv(KE-1,i,j  ) )
          call CHECK( __LINE__, GSQRT(KE-1,i ,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(KE-1,i-1,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(KE-1,i,j-1,I_XVZ) )
          call CHECK( __LINE__, GSQRT(KE-1,i,j  ,I_XVZ) )
          call CHECK( __LINE__, GSQRT(KE-1,i,j,I_XYZ) )
          call CHECK( __LINE__, St(KE-1,i,j) )
          call CHECK( __LINE__, DPRES(KE-2,i,j) )
          call CHECK( __LINE__, DPRES(KE-1,i,j) )
          call CHECK( __LINE__, DPRES(KE  ,i,j) )
          call CHECK( __LINE__, RCs2(KE-2,i,j) )
          call CHECK( __LINE__, RCs2(KE-1,i,j) )
          call CHECK( __LINE__, RCs2(KE  ,i,j) )
          call CHECK( __LINE__, RHOT(KE-2,i,j) )
          call CHECK( __LINE__, RHOT(KE,i,j) )
          call CHECK( __LINE__, DDENS(KE-2,i,j) )
          call CHECK( __LINE__, DDENS(KE,i,j) )
#endif
          B(KE-1,i,j) = ( &
               ( J33G * ( MOMZ(KE-1,i,j) + dtrk*Sw(KE-1,i,j) ) &
                 * ( 0.5_RP*(POTT(KE  ,i,j)+POTT(KE-1,i,j)) ) &
               - J33G * ( MOMZ(KE-2,i,j) + dtrk*Sw(KE-2,i,j) ) &
                 * ( FACT_N*(POTT(KE-1,i,j)+POTT(KE-2,i,j)) &
                   + FACT_F*(POTT(KE  ,i,j)+POTT(KE-3,i,j)) ) ) * RCDZ(KE-1) &
             + ( GSQRT(KE-1,i  ,j,I_UYZ) * ( MOMX(KE-1,i  ,j) + dtrk*Su(KE-1,i  ,j) ) &
                 * ( FACT_N*(POTT(KE-1,i+1,j)+POTT(KE-1,i  ,j)) &
                   + FACT_F*(POTT(KE-1,i+2,j)+POTT(KE-1,i-1,j)) ) &
               - GSQRT(KE-1,i-1,j,I_UYZ) * ( MOMX(KE-1,i-1,j) + dtrk*Su(KE-1,i-1,j) ) &
                 * ( FACT_N*(POTT(KE-1,i  ,j)+POTT(KE-1,i-1,j)) &
                   + FACT_F*(POTT(KE-1,i+1,j)+POTT(KE-1,i-2,j)) ) ) * RCDX(i) &
             + ( GSQRT(KE-1,i,j  ,I_XVZ) * ( MOMY(KE-1,i,j  ) + dtrk*Sv(KE-1,i,j  ) ) &
                 * ( FACT_N*(POTT(KE-1,i,j+1)+POTT(KE-1,i,j  )) &
                   + FACT_F*(POTT(KE-1,i,j+2)+POTT(KE-1,i,j-1)) ) &
               - GSQRT(KE-1,i,j-1,I_XVZ) * ( MOMY(KE-1,i,j-1) + dtrk*Sv(KE-1,i,j-1) ) &
                 * ( FACT_N*(POTT(KE-1,i,j  )+POTT(KE-1,i,j-1)) &
                   + FACT_F*(POTT(KE-1,i,j+1)+POTT(KE-1,i,j-2)) ) ) * RCDY(j) &
             + GSQRT(KE-1,i,j,I_XYZ) * ( St(KE-1,i,j) - DPRES(KE-1,i,j) * RCs2(KE-1,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( ( DPRES(KE  ,i,j)*RCs2(KE  ,i,j) &
                               - DPRES(KE-2,i,j)*RCs2(KE-2,i,j) ) &
                             - ( RHOT(KE  ,i,j)*DDENS(KE  ,i,j)/DENS(KE  ,i,j) &
                               - RHOT(KE-2,i,j)*DDENS(KE-2,i,j)/DENS(KE-2,i,j) )&
                             - ( DENS(KE  ,i,j)*(RHOT_RK(KE  ,i,j)/DENS_RK(KE  ,i,j) - POTT(KE  ,i,j) ) &
                               - DENS(KE-2,i,j)*(RHOT_RK(KE-2,i,j)/DENS_RK(KE-2,i,j) - POTT(KE-2,i,j) ) ) &
                             ) / ( FDZ(KE-1) + FDZ(KE-1-1) )
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
          call CHECK( __LINE__, MOMX(KE,i-1,j) )
          call CHECK( __LINE__, MOMX(KE,i  ,j) )
          call CHECK( __LINE__, MOMY(KE,i,j-1) )
          call CHECK( __LINE__, MOMY(KE,i,j  ) )
          call CHECK( __LINE__, Sw(KE-1,i,j) )
          call CHECK( __LINE__, Su(KE,i-1,j) )
          call CHECK( __LINE__, Su(KE,i  ,j) )
          call CHECK( __LINE__, Sv(KE,i,j-1) )
          call CHECK( __LINE__, Sv(KE,i,j  ) )
          call CHECK( __LINE__, GSQRT(KE,i-1,j,I_UYZ) )
          call CHECK( __LINE__, POTT(KE-1,i,j) )
          call CHECK( __LINE__, POTT(KE  ,i,j) )
          call CHECK( __LINE__, POTT(KE,i-2,j) )
          call CHECK( __LINE__, POTT(KE,i-1,j) )
          call CHECK( __LINE__, POTT(KE,i  ,j) )
          call CHECK( __LINE__, POTT(KE,i+1,j) )
          call CHECK( __LINE__, POTT(KE,i+2,j) )
          call CHECK( __LINE__, POTT(KE,i,j-2) )
          call CHECK( __LINE__, POTT(KE,i,j-1) )
          call CHECK( __LINE__, POTT(KE,i,j  ) )
          call CHECK( __LINE__, POTT(KE,i,j+1) )
          call CHECK( __LINE__, POTT(KE,i,j+2) )
          call CHECK( __LINE__, GSQRT(KE,i-1,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(KE,i  ,j,I_UYZ) )
          call CHECK( __LINE__, GSQRT(KE,i,j-1,I_XVZ) )
          call CHECK( __LINE__, GSQRT(KE,i,j  ,I_XVZ) )
          call CHECK( __LINE__, GSQRT(KE,i,j,I_XYZ) )
          call CHECK( __LINE__, St(KE,i,j) )
          call CHECK( __LINE__, DPRES(KE-1,i,j) )
          call CHECK( __LINE__, DPRES(KE  ,i,j) )
          call CHECK( __LINE__, RCs2(KE-1,i,j) )
          call CHECK( __LINE__, RCs2(KE  ,i,j) )
          call CHECK( __LINE__, RHOT(KE-1,i,j) )
          call CHECK( __LINE__, DDENS(KE-1,i,j) )
#endif
          B(KE,i,j) = ( &
               ( &
               - J33G * ( MOMZ(KE-1,i,j) + dtrk*Sw(KE-1,i,j) ) &
                 * 0.5_RP*(POTT(KE  ,i,j)+POTT(KE-1,i,j)) ) * RCDZ(KE) &
             + ( GSQRT(KE,i  ,j,I_UYZ) * ( MOMX(KE,i  ,j) + dtrk*Su(KE,i  ,j) ) &
                 * ( FACT_N*(POTT(KE,i+1,j)+POTT(KE,i  ,j)) &
                   + FACT_F*(POTT(KE,i+2,j)+POTT(KE,i-1,j)) ) &
               - GSQRT(KE,i-1,j,I_UYZ) * ( MOMX(KE,i-1,j) + dtrk*Su(KE,i-1,j) ) &
                 * ( FACT_N*(POTT(KE,i  ,j)+POTT(KE,i-1,j)) &
                   + FACT_F*(POTT(KE,i+1,j)+POTT(KE,i-2,j)) ) ) * RCDX(i) &
             + ( GSQRT(KE,i,j  ,I_XVZ) * ( MOMY(KE,i,j  ) + dtrk*Sv(KE,i,j  ) ) &
                 * ( FACT_N*(POTT(KE,i,j+1)+POTT(KE,i,j  )) &
                   + FACT_F*(POTT(KE,i,j+2)+POTT(KE,i,j-1)) ) &
               - GSQRT(KE,i,j-1,I_XVZ) * ( MOMY(KE,i,j-1) + dtrk*Sv(KE,i,j-1) ) &
                 * ( FACT_N*(POTT(KE,i,j  )+POTT(KE,i,j-1)) &
                   + FACT_F*(POTT(KE,i,j+1)+POTT(KE,i,j-2)) ) ) * RCDY(j) &
             + GSQRT(KE,i,j,I_XYZ) * ( St(KE,i,j) - DPRES(KE,i,j) * RCs2(KE,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( - DPRES(KE-1,i,j)*RCs2(KE-1,i,j) &
                               + RHOT(KE-1,i,j)*DDENS(KE-1,i,j)/DENS(KE-1,i,j) &
                               - DENS(KE-1,i,j)*(RHOT_RK(KE-1,i,j)/DENS_RK(KE-1,i,j) - POTT(KE-1,i,j) ) ) &
                             / ( FDZ(KE) + FDZ(KE-1) )
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       call make_matrix( &
            M, & ! (out)
            POTT, RCs2, GRAV, &
            GSQRT, J33G, &
            RCDZ, RFDZ, RCDX, RFDX, RCDY,RFDY, FDZ, &
            rdt, &
            FACT_N, FACT_F, &
            I_XYZ, I_XYW, &
            IIS, IIE, JJS, JJE )

    end do
    end do

    call solve_bicgstab( &
       DPRES, & ! (inout)
       M, B   ) ! (in)

#ifdef DEBUG
    call check_solver( DPRES, M, B )
#endif

    call COMM_vars8( DPRES, 1 )
    call COMM_wait ( DPRES, 1 )


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#else
    end do
    end do

    call solve_multigrid

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#endif

       !##### momentum equation (z) #####

       !--- update momentum(z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, DPRES(k+1,i,j) )
          call CHECK( __LINE__, DPRES(k  ,i,j) )
          call CHECK( __LINE__, PRES(k+1,i,j) )
          call CHECK( __LINE__, PRES(k  ,i,j) )
          call CHECK( __LINE__, DDENS(k+1,i,j) )
          call CHECK( __LINE__, DDENS(k  ,i,j) )
          call CHECK( __LINE__, REF_DENS(k+1,i,j) )
          call CHECK( __LINE__, REF_DENS(k,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
#endif
          duvw =  dtrk * ( ( &
                   - J33G * ( DPRES(k+1,i,j) - DPRES(k,i,j) ) * RFDZ(k) &
                   ) / GSQRT(k,i,j,I_UYZ) &
                 - 0.5_RP * GRAV &
                 * ( DDENS(k+1,i,j) &
                   + ( DPRES(k+1,i,j) - (PRES(k+1,i,j)-REF_pres(k+1,i,j)) ) &
                   * DENS(k+1,i,j) / ( kappa(k+1,i,j)*PRES(k+1,i,j) ) &
                   + DDENS(k  ,i,j) &
                   + ( DPRES(k  ,i,j) - (PRES(k  ,i,j)-REF_pres(k  ,i,j)) ) &
                   * DENS(k  ,i,j) / ( kappa(k  ,i,j)*PRES(k  ,i,j) ) ) &
                 + Sw(k,i,j) )
          MOMZ_RK(k,i,j) = MOMZ0(k,i,j) + duvw
          mflx_hi(k,i,j,ZDIR) = J33G * ( MOMZ(k,i,j) + duvw )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          mflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


       !##### momentum equation (x) #####

       !--- update momentum(x)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DPRES(k,i+1,j) )
          call CHECK( __LINE__, DPRES(k,i  ,j) )
          call CHECK( __LINE__, GSQRT(k,i+1,j,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i  ,j,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
          call CHECK( __LINE__, Su(k,i,j) )
          call CHECK( __LINE__, MOMX0(k,i,j) )
#endif
          duvw = dtrk * ( ( &
                 - ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) &
                   - GSQRT(k,i  ,j,I_XYZ) * DPRES(k,i  ,j) ) * RFDX(i) &  ! pressure gradient force
                   ) / GSQRT(k,i,j,I_UYZ) &
                 + Su(k,i,j) )

          MOMX_RK(k,i,j) = MOMX0(k,i,j) + duvw
          mflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) * ( MOMX(k,i,j) + duvw )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (y) #####

       !--- update momentum(y)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DPRES(k,i,j  ) )
          call CHECK( __LINE__, DPRES(k,i,j+1) )
          call CHECK( __LINE__, GSQRT(k,i,j  ,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i,j+1,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
          call CHECK( __LINE__, Sv(k,i,j) )
          call CHECK( __LINE__, MOMY0(k,i,j) )
#endif
          duvw = dtrk * ( ( &
                            - ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) &
                              - GSQRT(k,i,j  ,I_XYZ) * DPRES(k,i,j  ) ) * RFDY(j) &
                          ) / GSQRT(k,i,j,I_XVZ) &
                          + Sv(k,i,j) )
          MOMY_RK(k,i,j) = MOMY0(k,i,j) + duvw
          mflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) * ( MOMY(k,i,j) + duvw )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


       !##### Thermodynamic Equation #####

       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j) + POTT(k  ,i,j) ) &
                                + FACT_F * ( POTT(k+2,i,j) + POTT(k-1,i,j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = mflx_hi(KS  ,i,j,ZDIR) * 0.5_RP * ( POTT(KS+1,i,j) + POTT(KS  ,i,j) )
          qflx_hi(KE-1,i,j,ZDIR) = mflx_hi(KE-1,i,j,ZDIR) * 0.5_RP * ( POTT(KE  ,i,j) + POTT(KE-1,i,j) )
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) &
                                * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                  + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j+2) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,YDIR) )
#endif
          qflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) &
                                * ( FACT_N * ( POTT(k,i,j+1)+POTT(k,i,j  ) ) &
                                  + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT_RK(k,i,j) = RHOT0(k,i,j) &
               + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                            + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                            / GSQRT(k,i,j,I_XYZ) &
                          + St(k,i,j) )
       end do
       end do
       end do

       !##### continuous equation #####

       ! total momentum flux
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + mflx_hi2(k,i,j,ZDIR)
       end do
       end do
       end do
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
          mflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) + mflx_hi2(k,i,j,XDIR)
       end do
       end do
       end do
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
          mflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) + mflx_hi2(k,i,j,YDIR)
       end do
       end do
       end do

       !--- update density
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS0(k,i,j) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DENS_t(k,i,j) )
#endif
          DENS_RK(k,i,j) = DENS0(k,i,j) &
               + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i  ,j,  ZDIR) ) * RCDZ(k) &
                            + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                            + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) &
                            / GSQRT(k,i,j,I_XYZ) & ! divergence
                          + DENS_t(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif



#ifdef DEBUG
       call check_pres( &
            DPRES, REF_pres, PRES, RHOT_RK, RHOT &
#ifdef DRY
            , kappa
#endif
    )
#endif

    enddo
    enddo

  end subroutine ATMOS_DYN_rk_hivi

#ifdef HIVI_BICGSTAB
  subroutine solve_bicgstab( &
       DPRES, &
       M, B )
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none
    real(RP), intent(inout) :: DPRES(KA,IA,JA)
    real(RP), intent(in) :: M(7,KA,IA,JA)
    real(RP), intent(in) :: B(KA,IA,JA)

    real(RP) :: r0(KA,IA,JA)

    real(RP) :: p(KA,IA,JA)
    real(RP) :: Mp(KA,IA,JA)
    real(RP) :: s(KA,IA,JA)
    real(RP) :: Ms(KA,IA,JA)
    real(RP) :: al, be, w

    real(RP), pointer :: r(:,:,:)
    real(RP), pointer :: rn(:,:,:)
    real(RP), pointer :: swap(:,:,:)
    real(RP), target :: v0(KA,IA,JA)
    real(RP), target :: v1(KA,IA,JA)
    real(RP) :: r0r
    real(RP) :: norm, error, error2

    real(RP) :: iprod(2)
    real(RP) :: buf(2)

    integer :: k, i, j
    integer :: iis, iie, jjs, jje
    integer :: iter
    integer :: ierror

#ifdef DEBUG
    r0(:,:,:) = UNDEF
    p (:,:,:) = UNDEF
    Mp(:,:,:) = UNDEF
    s (:,:,:) = UNDEF
    Ms(:,:,:) = UNDEF

    v0(:,:,:) = UNDEF
    v1(:,:,:) = UNDEF
#endif

    r  => v0
    rn => v1

    call mul_matrix( v1, M, DPRES ) ! v1 = M x0

    norm = 0.0_RP
    r0r  = 0.0_RP

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          norm = norm + B(k,i,j)**2
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! r = b - M x0
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, B(k,i,j) )
          call CHECK( __LINE__, v1(k,i,j) )
#endif
          r(k,i,j) = B(k,i,j) - v1(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          r0(k,i,j) = r(k,i,j)
          p(k,i,j) = r(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, r(k,i,j) )
          call CHECK( __LINE__, r0(k,i,j) )
#endif
          r0r = r0r + r0(k,i,j) * r(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    end do
    end do

    iprod(1) = r0r
    iprod(2) = norm
    call MPI_AllReduce(iprod, buf, 2, mtype, MPI_SUM, MPI_COMM_WORLD, ierror)
    r0r = buf(1)
    norm = buf(2)

    error2 = norm

    do iter = 1, ITMAX

       error = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, r(k,i,j) )
#endif
          error = error + r(k,i,j)**2
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       call MPI_AllReduce(error, buf, 1, mtype, MPI_SUM, MPI_COMM_WORLD, ierror)
       error = buf(1)

#ifdef DEBUG
       if (IO_L) write(*,*) iter, error/norm
#endif
       if ( sqrt(error/norm) < epsilon .or. error > error2 ) then
#ifdef DEBUG
         IF ( IO_L ) write(*,*) "Bi-CGSTAB converged:", iter
#endif
          exit
       end if
       error2 = error

       call COMM_vars8( p, 1 )
       call COMM_wait ( p, 1 )
       call mul_matrix( Mp, M, p )

       iprod(1) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, r0(k,i,j) )
          call CHECK( __LINE__, Mp(k,i,j) )
#endif
          iprod(1) = iprod(1) + r0(k,i,j) * Mp(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       call MPI_AllReduce(iprod, buf, 1, mtype, MPI_SUM, MPI_COMM_WORLD, ierror)
       al = r0r / buf(1) ! (r0,r) / (r0,Mp)

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, r (k,i,j) )
          call CHECK( __LINE__, Mp(k,i,j) )
#endif
          s(k,i,j) = r(k,i,j) - al*Mp(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       call COMM_vars8( s, 1 )
       call COMM_wait ( s, 1 )
       call mul_matrix( Ms, M, s )
       iprod(1) = 0.0_RP
       iprod(2) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, Ms(k,i,j) )
          call CHECK( __LINE__,  s(k,i,j) )
#endif
          iprod(1) = iprod(1) + Ms(k,i,j) *  s(k,i,j)
          iprod(2) = iprod(2) + Ms(k,i,j) * Ms(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       call MPI_AllReduce(iprod, buf, 2, mtype, MPI_SUM, MPI_COMM_WORLD, ierror)
       w = buf(1) / buf(2) ! (Ms,s) / (Ms,Ms)

       iprod(1) = 0.0_RP
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DPRES(k,i,j) )
          call CHECK( __LINE__, p(k,i,j) )
          call CHECK( __LINE__, s(k,i,j) )
#endif
          DPRES(k,i,j) = DPRES(k,i,j) + al*p(k,i,j) + w*s(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, s(k,i,j) )
          call CHECK( __LINE__, Ms(k,i,j) )
#endif
          rn(k,i,j) = s(k,i,j) - w*Ms(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       be = al/w / r0r
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, r0(k,i,j) )
          call CHECK( __LINE__, rn(k,i,j) )
#endif
          iprod(1) = iprod(1) + r0(k,i,j) * rn(k,i,j)
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    end do
    end do

       call MPI_AllReduce(iprod, r0r, 1, mtype, MPI_SUM, MPI_COMM_WORLD, ierror)

       be = be * r0r ! al/w * (r0,rn)/(r0,r)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, rn(k,i,j) )
          call CHECK( __LINE__, p(k,i,j) )
          call CHECK( __LINE__, Mp(k,i,j) )
#endif
          p(k,i,j) = rn(k,i,j) + be * ( p(k,i,j) - w*Mp(k,i,j) )
       end do
       end do
       end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       swap => rn
       rn => r
       r => swap
#ifdef DEBUG
       rn(:,:,:) = UNDEF
#endif
    end do

    if ( iter >= ITMAX ) then
       write(*,*) 'xxx [atmos_dyn_hivi] Bi-CGSTAB'
       write(*,*) 'xxx not converged', error, norm
       call PRC_MPIstop
    end if

    return
  end subroutine solve_bicgstab

  subroutine make_matrix(&
       M, &
       POTT, RCs2, GRAV, &
       G, J33G, &
       RCDZ, RFDZ, RCDX, RFDX, RCDY, RFDY, FDZ, &
       rdt, &
       FACT_N, FACT_F, &
       I_XYZ, I_XYW, &
       IIS, IIE, JJS, JJE )
    implicit none
    real(RP), intent(out) :: M(7,KA,IA,JA)
    real(RP), intent(in) :: POTT(KA,IA,JA)
    real(RP), intent(in) :: RCs2(KA,IA,JA)
    real(RP), intent(in) :: GRAV
    real(RP), intent(in) :: G(KA,IA,JA,7)
    real(RP), intent(in) :: J33G
    real(RP), intent(in) :: RCDZ(KA)
    real(RP), intent(in) :: RFDZ(KA-1)
    real(RP), intent(in) :: RFDX(IA-1)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: RFDY(JA-1)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: rdt
    real(RP), intent(in) :: FACT_N
    real(RP), intent(in) :: FACT_F
    integer, intent(in) :: I_XYZ
    integer, intent(in) :: I_XYW
    integer, intent(in) :: IIS
    integer, intent(in) :: IIE
    integer, intent(in) :: JJS
    integer, intent(in) :: JJE

    integer :: k, i, j

#ifdef DEBUG
    M(:,:,:,:) = UNDEF
#endif

    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, POTT(k-2,i,j) )
          call CHECK( __LINE__, POTT(k-1,i,j) )
          call CHECK( __LINE__, POTT(k  ,i,j) )
          call CHECK( __LINE__, POTT(k+1,i,j) )
          call CHECK( __LINE__, POTT(k+2,i,j) )
          call CHECK( __LINE__, POTT(k,i-2,j) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+2,j) )
          call CHECK( __LINE__, POTT(k,i,j-2) )
          call CHECK( __LINE__, POTT(k,i,j-1) )
          call CHECK( __LINE__, POTT(k,i,j  ) )
          call CHECK( __LINE__, POTT(k,i,j+1) )
          call CHECK( __LINE__, POTT(k,i,j+2) )
          call CHECK( __LINE__, G(k-1,i,j,I_XYW) )
          call CHECK( __LINE__, G(k+1,i,j,I_XYW) )
          call CHECK( __LINE__, G(k,i,j,I_XYZ) )
          call CHECK( __LINE__, RCs2(k-1,i,j) )
          call CHECK( __LINE__, RCs2(k  ,i,j) )
          call CHECK( __LINE__, RCs2(k+1,i,j) )
#endif
       ! k,i,j
       M(1,k,i,j) = &
           - ( ( FACT_N * (POTT(k+1,i,j)+POTT(k  ,i,j)) &
               + FACT_F * (POTT(k+2,i,j)+POTT(k-1,i,j)) ) * RFDZ(k  ) / G(k  ,i,j,I_XYW) &
             + ( FACT_N * (POTT(k  ,i,j)+POTT(k-1,i,j)) &
               + FACT_F * (POTT(k+1,i,j)+POTT(k-2,i,j)) ) * RFDZ(k-1) / G(k-1,i,j,I_XYW) &
              ) * J33G * J33G * RFDZ(k) &
           - ( ( FACT_N * (POTT(k,i+1,j)+POTT(k,i  ,j)) &
               + FACT_F * (POTT(k,i+2,j)+POTT(k,i-1,j)) ) * RFDX(i  ) &
             + ( FACT_N * (POTT(k,i  ,j)+POTT(k,i-1,j)) &
               + FACT_F * (POTT(k,i+1,j)+POTT(k,i-2,j)) ) * RFDX(i-1) &
             ) * G(k,i,j,I_XYZ) * RFDX(i) &
           - ( ( FACT_N * (POTT(k,i,j+1)+POTT(k,i,j  )) &
               + FACT_F * (POTT(k,i,j+2)+POTT(k,i,j-1)) ) * RFDY(j  ) &
             + ( FACT_N * (POTT(k,i,j  )+POTT(k,i,j-1)) &
               + FACT_F * (POTT(k,i,j-1)+POTT(k,i,j-2)) ) * RFDY(j-1) &
             ) * G(k,i,j,I_XYZ) * RFDY(j) &
           - G(k,i,j,I_XYZ) * RCs2(k,i,j) * rdt * rdt
       ! k-1
       M(2,k,i,j) = J33G * J33G / G(k-1,i,j,I_XYW) &
                  * ( FACT_N * (POTT(k  ,i,j)+POTT(k-1,i,j)) &
                    + FACT_F * (POTT(k+1,i,j)+POTT(k-2,i,j)) ) &
                  * RFDZ(k-1) * RCDZ(k) &
                  - GRAV * J33G * RCs2(k-1,i,j) / ( FDZ(k)+FDZ(k-1) )
       ! k+1
       M(3,k,i,j) = J33G * J33G / G(k+1,i,j,I_XYW) &
                  * ( FACT_N * (POTT(k+1,i,j)+POTT(k  ,i,j)) &
                    + FACT_F * (POTT(k+2,i,j)+POTT(k-1,i,j)) ) &
                  * RFDZ(k  ) * RCDZ(k) &
                  + GRAV * J33G * RCs2(k+1,i,j) / ( FDZ(k)+FDZ(k-1) )
    end do
    end do
    end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
    do j = JJS, JJE
    do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, POTT(KS+1,i,j) )
          call CHECK( __LINE__, POTT(KS,i-2,j) )
          call CHECK( __LINE__, POTT(KS,i-1,j) )
          call CHECK( __LINE__, POTT(KS,i  ,j) )
          call CHECK( __LINE__, POTT(KS,i+1,j) )
          call CHECK( __LINE__, POTT(KS,i+2,j) )
          call CHECK( __LINE__, POTT(KS,i,j-2) )
          call CHECK( __LINE__, POTT(KS,i,j-1) )
          call CHECK( __LINE__, POTT(KS,i,j  ) )
          call CHECK( __LINE__, POTT(KS,i,j+1) )
          call CHECK( __LINE__, POTT(KS,i,j+2) )
          call CHECK( __LINE__, G(KS,i,j,I_XYZ) )
          call CHECK( __LINE__, RCs2(KS,i,j) )
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, POTT(KS+1,i,j) )
          call CHECK( __LINE__, G(KS+1,i,j,I_XYW) )
          call CHECK( __LINE__, RCs2(KS+1,i,j) )
#endif
       ! k,i,j
       M(1,KS,i,j) = &
            - ( 0.5_RP * (POTT(KS+1,i,j)+POTT(KS  ,i,j)) * RFDZ(KS  ) &
              ) * J33G * J33G / G(KS  ,i,j,I_XYW) * RFDZ(KS) &
            - ( ( FACT_N * (POTT(KS,i+1,j)+POTT(KS,i  ,j)) &
                + FACT_F * (POTT(KS,i+2,j)+POTT(KS,i-1,j)) ) * RFDX(i  ) &
              + ( FACT_N * (POTT(KS,i  ,j)+POTT(KS,i-1,j)) &
                + FACT_F * (POTT(KS,i+1,j)+POTT(KS,i-2,j)) ) * RFDX(i-1) &
              ) * G(KS,i,j,I_XYZ) * RFDX(i) &
            - ( ( FACT_N * (POTT(KS,i,j+1)+POTT(KS,i,j  )) &
                + FACT_F * (POTT(KS,i,j+2)+POTT(KS,i,j-1)) ) * RFDY(j  ) &
              + ( FACT_N * (POTT(KS,i,j  )+POTT(KS,i,j-1)) &
                + FACT_F * (POTT(KS,i,j-1)+POTT(KS,i,j-2)) ) * RFDY(j-1) &
              ) * G(KS,i,j,I_XYZ) * RFDY(j) &
            - G(KS,i,j,I_XYZ) * RCs2(KS,i,j) * rdt * rdt
       ! k+1
       M(3,KS,i,j) = J33G * J33G / G(KS+1,i,j,I_XYW) &
                   * 0.5_RP * (POTT(KS+1,i,j)+POTT(KS  ,i,j)) &
                   * RFDZ(KS  ) * RCDZ(KS) &
                   + GRAV * J33G * RCs2(KS+1,i,j) / ( FDZ(KS)+FDZ(KS-1) )
#ifdef DEBUG
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, POTT(KS+1,i,j) )
          call CHECK( __LINE__, POTT(KS+2,i,j) )
          call CHECK( __LINE__, POTT(KS+3,i,j) )
          call CHECK( __LINE__, POTT(KS+1,i-2,j) )
          call CHECK( __LINE__, POTT(KS+1,i-1,j) )
          call CHECK( __LINE__, POTT(KS+1,i  ,j) )
          call CHECK( __LINE__, POTT(KS+1,i+1,j) )
          call CHECK( __LINE__, POTT(KS+1,i+2,j) )
          call CHECK( __LINE__, POTT(KS+1,i,j-2) )
          call CHECK( __LINE__, POTT(KS+1,i,j-1) )
          call CHECK( __LINE__, POTT(KS+1,i,j  ) )
          call CHECK( __LINE__, POTT(KS+1,i,j+1) )
          call CHECK( __LINE__, POTT(KS+1,i,j+2) )
          call CHECK( __LINE__, G(KS  ,i,j,I_XYW) )
          call CHECK( __LINE__, G(KS+2,i,j,I_XYW) )
          call CHECK( __LINE__, G(KS+1,i,j,I_XYZ) )
          call CHECK( __LINE__, RCs2(KS  ,i,j) )
          call CHECK( __LINE__, RCs2(KS+1  ,i,j) )
          call CHECK( __LINE__, RCs2(KS+2,i,j) )
#endif
       ! k,i,j
       M(1,KS+1,i,j) = &
           - ( ( FACT_N * (POTT(KS+2,i,j)+POTT(KS+1,i,j)) &
               + FACT_F * (POTT(KS+3,i,j)+POTT(KS  ,i,j)) ) * RFDZ(KS+1) / G(KS+1,i,j,I_XYW) &
             +   0.5_RP * (POTT(KS+1,i,j)+POTT(KS  ,i,j))   * RFDZ(KS  ) / G(KS  ,i,j,I_XYW) &
              ) * J33G * J33G * RFDZ(KS+1) &
           - ( ( FACT_N * (POTT(KS+1,i+1,j)+POTT(KS+1,i  ,j)) &
               + FACT_F * (POTT(KS+1,i+2,j)+POTT(KS+1,i-1,j)) ) * RFDX(i  ) &
             + ( FACT_N * (POTT(KS+1,i  ,j)+POTT(KS+1,i-1,j)) &
               + FACT_F * (POTT(KS+1,i+1,j)+POTT(KS+1,i-2,j)) ) * RFDX(i-1) &
             ) * G(KS+1,i,j,I_XYZ) * RFDX(i) &
           - ( ( FACT_N * (POTT(KS+1,i,j+1)+POTT(KS+1,i,j  )) &
               + FACT_F * (POTT(KS+1,i,j+2)+POTT(KS+1,i,j-1)) ) * RFDY(j  ) &
             + ( FACT_N * (POTT(KS+1,i,j  )+POTT(KS+1,i,j-1)) &
               + FACT_F * (POTT(KS+1,i,j-1)+POTT(KS+1,i,j-2)) ) * RFDY(j-1) &
             ) * G(KS+1,i,j,I_XYZ) * RFDY(j) &
           - G(KS+1,i,j,I_XYZ) * RCs2(KS+1,i,j) * rdt * rdt
       ! k-1
       M(2,KS+1,i,j) = J33G * J33G / G(KS,i,j,I_XYW) &
                  * 0.5_RP * (POTT(KS+1,i,j)+POTT(KS  ,i,j)) &
                  * RFDZ(KS  ) * RCDZ(KS+1) &
                  - GRAV * J33G * RCs2(KS  ,i,j) / ( FDZ(KS+1)+FDZ(KS) )
       ! k+1
       M(3,KS+1,i,j) = J33G * J33G / G(KS+2,i,j,I_XYW) &
                  * ( FACT_N * (POTT(KS+2,i,j)+POTT(KS+1,i,j)) &
                    + FACT_F * (POTT(KS+3,i,j)+POTT(KS  ,i,j)) ) &
                  * RFDZ(KS+1) * RCDZ(KS+1) &
                  + GRAV * J33G * RCs2(KS+2,i,j) / ( FDZ(KS+1)+FDZ(KS) )
#ifdef DEBUG
          call CHECK( __LINE__, POTT(KE-3,i,j) )
          call CHECK( __LINE__, POTT(KE-2,i,j) )
          call CHECK( __LINE__, POTT(KE-1,i,j) )
          call CHECK( __LINE__, POTT(KE  ,i,j) )
          call CHECK( __LINE__, POTT(KE-1,i-2,j) )
          call CHECK( __LINE__, POTT(KE-1,i-1,j) )
          call CHECK( __LINE__, POTT(KE-1,i  ,j) )
          call CHECK( __LINE__, POTT(KE-1,i+1,j) )
          call CHECK( __LINE__, POTT(KE-1,i+2,j) )
          call CHECK( __LINE__, POTT(KE-1,i,j-2) )
          call CHECK( __LINE__, POTT(KE-1,i,j-1) )
          call CHECK( __LINE__, POTT(KE-1,i,j  ) )
          call CHECK( __LINE__, POTT(KE-1,i,j+1) )
          call CHECK( __LINE__, POTT(KE-1,i,j+2) )
          call CHECK( __LINE__, G(KE-2,i,j,I_XYW) )
          call CHECK( __LINE__, G(KE  ,i,j,I_XYW) )
          call CHECK( __LINE__, G(KE-1,i,j,I_XYZ) )
          call CHECK( __LINE__, RCs2(KE-2,i,j) )
          call CHECK( __LINE__, RCs2(KE-1  ,i,j) )
          call CHECK( __LINE__, RCs2(KE  ,i,j) )
#endif
       ! k,i,j
       M(1,KE-1,i,j) = &
           - (   0.5_RP * (POTT(KE  ,i,j)+POTT(KE-1,i,j))   * RFDZ(KE-1) / G(KE-1,i,j,I_XYW) &
             + ( FACT_N * (POTT(KE-1,i,j)+POTT(KE-2,i,j)) &
               + FACT_F * (POTT(KE  ,i,j)+POTT(KE-3,i,j)) ) * RFDZ(KE-2) / G(KE-2,i,j,I_XYW) &
              ) * J33G * J33G * RFDZ(KE-1) &
           - ( ( FACT_N * (POTT(KE-1,i+1,j)+POTT(KE-1,i  ,j)) &
               + FACT_F * (POTT(KE-1,i+2,j)+POTT(KE-1,i-1,j)) ) * RFDX(i  ) &
             + ( FACT_N * (POTT(KE-1,i  ,j)+POTT(KE-1,i-1,j)) &
               + FACT_F * (POTT(KE-1,i+1,j)+POTT(KE-1,i-2,j)) ) * RFDX(i-1) &
             ) * G(KE-1,i,j,I_XYZ) * RFDX(i) &
           - ( ( FACT_N * (POTT(KE-1,i,j+1)+POTT(KE-1,i,j  )) &
               + FACT_F * (POTT(KE-1,i,j+2)+POTT(KE-1,i,j-1)) ) * RFDY(j  ) &
             + ( FACT_N * (POTT(KE-1,i,j  )+POTT(KE-1,i,j-1)) &
               + FACT_F * (POTT(KE-1,i,j-1)+POTT(KE-1,i,j-2)) ) * RFDY(j-1) &
             ) * G(KE-1,i,j,I_XYZ) * RFDY(j) &
           - G(KE-1,i,j,I_XYZ) * RCs2(KE-1,i,j) * rdt * rdt
       ! k-1
       M(2,KE-1,i,j) = J33G * J33G / G(KE-2,i,j,I_XYW) &
                  * ( FACT_N * (POTT(KE-1,i,j)+POTT(KE-2,i,j)) &
                    + FACT_F * (POTT(KE  ,i,j)+POTT(KE-3,i,j)) ) &
                  * RFDZ(KE-2) * RCDZ(KE-1) &
                  - GRAV * J33G * RCs2(KE-2,i,j) / ( FDZ(KE-1)+FDZ(KE-2) )
       ! k+1
       M(3,KE-1,i,j) = J33G * J33G / G(KE  ,i,j,I_XYW) &
                  * 0.5_RP * (POTT(KE  ,i,j)+POTT(KE-1,i,j)) &
                  * RFDZ(KE-1) * RCDZ(KE-1) &
                  + GRAV * J33G * RCs2(KE  ,i,j) / ( FDZ(KE-1)+FDZ(KE-2) )
#ifdef DEBUG
          call CHECK( __LINE__, POTT(KE-1,i,j) )
          call CHECK( __LINE__, POTT(KE  ,i,j) )
          call CHECK( __LINE__, POTT(KE,i-2,j) )
          call CHECK( __LINE__, POTT(KE,i-1,j) )
          call CHECK( __LINE__, POTT(KE,i  ,j) )
          call CHECK( __LINE__, POTT(KE,i+1,j) )
          call CHECK( __LINE__, POTT(KE,i+2,j) )
          call CHECK( __LINE__, POTT(KE,i,j-2) )
          call CHECK( __LINE__, POTT(KE,i,j-1) )
          call CHECK( __LINE__, POTT(KE,i,j  ) )
          call CHECK( __LINE__, POTT(KE,i,j+1) )
          call CHECK( __LINE__, POTT(KE,i,j+2) )
          call CHECK( __LINE__, G(KE-1,i,j,I_XYW) )
          call CHECK( __LINE__, G(KE,i,j,I_XYZ) )
          call CHECK( __LINE__, RCs2(KE-1,i,j) )
          call CHECK( __LINE__, RCs2(KE,i,j) )
#endif
       ! k,i,j
       M(1,KE,i,j) = &
           - ( &
             + 0.5_RP * (POTT(KE  ,i,j)+POTT(KE-1,i,j)) * RFDZ(KE-1) / G(KE-1,i,j,I_XYW) &
              ) * J33G * J33G * RFDZ(KE) &
           - ( ( FACT_N * (POTT(KE,i+1,j)+POTT(KE,i  ,j)) &
               + FACT_F * (POTT(KE,i+2,j)+POTT(KE,i-1,j)) ) * RFDX(i  ) &
             + ( FACT_N * (POTT(KE,i  ,j)+POTT(KE,i-1,j)) &
               + FACT_F * (POTT(KE,i+1,j)+POTT(KE,i-2,j)) ) * RFDX(i-1) &
             ) * G(KE,i,j,I_XYZ) * RFDX(i) &
           - ( ( FACT_N * (POTT(KE,i,j+1)+POTT(KE,i,j  )) &
               + FACT_F * (POTT(KE,i,j+2)+POTT(KE,i,j-1)) ) * RFDY(j  ) &
             + ( FACT_N * (POTT(KE,i,j  )+POTT(KE,i,j-1)) &
               + FACT_F * (POTT(KE,i,j-1)+POTT(KE,i,j-2)) ) * RFDY(j-1) &
             ) * G(KE,i,j,I_XYZ) * RFDY(j) &
           - G(KE,i,j,I_XYZ) * RCs2(KE,i,j) * rdt * rdt
       ! k-1
       M(2,KE,i,j) = J33G * J33G / G(KE-1,i,j,I_XYW) &
                  * 0.5_RP * (POTT(KE  ,i,j)+POTT(KE-1,i,j)) &
                  * RFDZ(KE-1) * RCDZ(KE) &
                  - GRAV * J33G * RCs2(KE-1,i,j) / ( FDZ(KE)+FDZ(KE-1) )
    end do
    end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    do j = JJS, JJE
    do i = IIS, IIE
    do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, G(k,i-1,j,I_XYZ) )
          call CHECK( __LINE__, G(k,i+1,j,I_XYZ) )
          call CHECK( __LINE__, G(k,i,j-1,I_XYZ) )
          call CHECK( __LINE__, G(k,i,j+1,I_XYZ) )
          call CHECK( __LINE__, POTT(k,i-2,j  ) )
          call CHECK( __LINE__, POTT(k,i-1,j  ) )
          call CHECK( __LINE__, POTT(k,i  ,j  ) )
          call CHECK( __LINE__, POTT(k,i+1,j  ) )
          call CHECK( __LINE__, POTT(k,i+2,j  ) )
          call CHECK( __LINE__, POTT(k,i  ,j-2) )
          call CHECK( __LINE__, POTT(k,i  ,j-1) )
          call CHECK( __LINE__, POTT(k,i  ,j  ) )
          call CHECK( __LINE__, POTT(k,i  ,j+1) )
          call CHECK( __LINE__, POTT(k,i  ,j+2) )
#endif
       ! i-1
       M(4,k,i,j) = G(k,i-1,j,I_XYZ) &
                  * ( FACT_N * (POTT(k,i  ,j)+POTT(k,i-1,j)) &
                    + FACT_F * (POTT(k,i+1,j)+POTT(k,i-2,j)) ) &
                  * RFDX(i-1) * RCDX(i)
       ! i+1
       M(5,k,i,j) = G(k,i+1,j,I_XYZ) &
                  * ( FACT_N * (POTT(k,i+1,j)+POTT(k,i  ,j)) &
                    + FACT_F * (POTT(k,i+2,j)+POTT(k,i-1,j)) ) &
                  * RFDX(i  ) * RCDX(i)
       ! j-1
       M(6,k,i,j) = G(k,i,j-1,I_XYZ) &
                  * ( FACT_N * (POTT(k,i,j  )+POTT(k,i,j-1)) &
                    + FACT_F * (POTT(k,i,j+1)+POTT(k,i,j-2)) ) &
                  * RFDY(j-1) * RCDY(j)
       ! j+1
       M(7,k,i,j) = G(k,i,j+1,I_XYZ) &
                  * ( FACT_N * (POTT(k,i,j+1)+POTT(k,i,j  )) &
                    + FACT_F * (POTT(k,i,j+2)+POTT(k,i,j-1)) ) &
                  * RFDY(j  ) * RCDY(j)
    end do
    end do
    end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    return
  end subroutine make_matrix

  subroutine mul_matrix(V, M, C)
    implicit none
    real(RP), intent(out) :: V(KA,IA,JA)
    real(RP), intent(in)  :: M(7,KA,IA,JA)
    real(RP), intent(in)  :: C(KA,IA,JA)

    integer :: k, i, j

#ifdef DEBUG
    V(:,:,:) = UNDEF
#endif

    do j = JS, JE
    do i = IS, IE
    do k = KS+1, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, M(1,k,i,j) )
          call CHECK( __LINE__, M(2,k,i,j) )
          call CHECK( __LINE__, M(3,k,i,j) )
          call CHECK( __LINE__, M(4,k,i,j) )
          call CHECK( __LINE__, M(5,k,i,j) )
          call CHECK( __LINE__, M(6,k,i,j) )
          call CHECK( __LINE__, M(7,k,i,j) )
          call CHECK( __LINE__, C(k  ,i  ,j  ) )
          call CHECK( __LINE__, C(k-1,i  ,j  ) )
          call CHECK( __LINE__, C(k+1,i  ,j  ) )
          call CHECK( __LINE__, C(k  ,i-1,j  ) )
          call CHECK( __LINE__, C(k  ,i+1,j  ) )
          call CHECK( __LINE__, C(k  ,i  ,j-1) )
          call CHECK( __LINE__, C(k  ,i  ,j+1) )
#endif
       V(k,i,j) = M(1,k,i,j) * C(k  ,i  ,j  ) &
                + M(2,k,i,j) * C(k-1,i  ,j  ) &
                + M(3,k,i,j) * C(k+1,i  ,j  ) &
                + M(4,k,i,j) * C(k  ,i-1,j  ) &
                + M(5,k,i,j) * C(k  ,i+1,j  ) &
                + M(6,k,i,j) * C(k  ,i  ,j-1) &
                + M(7,k,i,j) * C(k  ,i  ,j+1)
    end do
    end do
    end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
    do j = JS, JE
    do i = IS, IE
#ifdef DEBUG
          call CHECK( __LINE__, M(1,KS,i,j) )
          call CHECK( __LINE__, M(3,KS,i,j) )
          call CHECK( __LINE__, M(4,KS,i,j) )
          call CHECK( __LINE__, M(5,KS,i,j) )
          call CHECK( __LINE__, M(6,KS,i,j) )
          call CHECK( __LINE__, M(7,KS,i,j) )
          call CHECK( __LINE__, C(KS  ,i  ,j  ) )
          call CHECK( __LINE__, C(KS+1,i  ,j  ) )
          call CHECK( __LINE__, C(KS  ,i-1,j  ) )
          call CHECK( __LINE__, C(KS  ,i+1,j  ) )
          call CHECK( __LINE__, C(KS  ,i  ,j-1) )
          call CHECK( __LINE__, C(KS  ,i  ,j+1) )
          call CHECK( __LINE__, M(1,KE,i,j) )
          call CHECK( __LINE__, M(2,KE,i,j) )
          call CHECK( __LINE__, M(4,KE,i,j) )
          call CHECK( __LINE__, M(5,KE,i,j) )
          call CHECK( __LINE__, M(6,KE,i,j) )
          call CHECK( __LINE__, M(7,KE,i,j) )
          call CHECK( __LINE__, C(KE  ,i  ,j  ) )
          call CHECK( __LINE__, C(KE-1,i  ,j  ) )
          call CHECK( __LINE__, C(KE  ,i-1,j  ) )
          call CHECK( __LINE__, C(KE  ,i+1,j  ) )
          call CHECK( __LINE__, C(KE  ,i  ,j-1) )
          call CHECK( __LINE__, C(KE  ,i  ,j+1) )
#endif
       V(KS,i,j) = M(1,KS,i,j) * C(KS  ,i  ,j  ) &
                 + M(3,KS,i,j) * C(KS+1,i  ,j  ) &
                 + M(4,KS,i,j) * C(KS  ,i-1,j  ) &
                 + M(5,KS,i,j) * C(KS  ,i+1,j  ) &
                 + M(6,KS,i,j) * C(KS  ,i  ,j-1) &
                 + M(7,KS,i,j) * C(KS  ,i  ,j+1)
       V(KE,i,j) = M(1,KE,i,j) * C(KE  ,i  ,j  ) &
                 + M(2,KE,i,j) * C(KE-1,i  ,j  ) &
                 + M(4,KE,i,j) * C(KE  ,i-1,j  ) &
                 + M(5,KE,i,j) * C(KE  ,i+1,j  ) &
                 + M(6,KE,i,j) * C(KE  ,i  ,j-1) &
                 + M(7,KE,i,j) * C(KE  ,i  ,j+1)
    end do
    end do
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    return
  end subroutine mul_matrix

#else
  subroutine solve_multigrid
  end subroutine solve_multigrid
#endif

#ifdef DEBUG
  subroutine check_solver( &
       DPRES, M, B )
    real(RP), intent(inout) :: DPRES(KA,IA,JA)
    real(RP), intent(in) :: M(7,KA,IA,JA)
    real(RP), intent(in) :: B(KA,IA,JA)

    real(RP) :: B2(KA,IA,JA)
    integer :: k, i, j

    call mul_matrix(B2, M, DPRES)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if ( abs(B2(k,i,j) - B(k,i,j) ) / B(k,i,j) > 1.e-5_RP ) then
          write(*,*) "solver error is too large: ", k,i,j, B(k,i,j) , B2(k,i,j)
          call abort
       end if
    end do
    end do
    end do

  end subroutine check_solver

  subroutine check_pres( &
       DPRES, REF_pres, PRES, &
       RHOT_RK, RHOT &
#ifdef DRY
       , kappa
#endif
       )
    use scale_const, only: &
       Rdry  => CONST_Rdry, &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       P00   => CONST_PRE00
    real(RP), intent(in) :: DPRES(KA,IA,JA)
    real(RP), intent(in) :: REF_pres(KA,IA,JA)
    real(RP), intent(in) :: PRES(KA,IA,JA)
    real(RP), intent(in) :: RHOT_RK(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
#ifdef DRY
    real(RP), intent(in) :: kappa(KA,IA,JA)
#endif

    real(RP) :: lhs, rhs
    real(RP) :: r, kapp
    integer :: k,i,j

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
#ifdef DRY
       R = Rtot(k,i,j)
       kapp = kappa(k,i,j)
#else
       R = Rdry
       kapp  = CPdry / CVdry
#endif
       lhs = DPRES(k,i,j) - ( PRES(k,i,j) - REF_pres(k,i,j) )
       rhs = kapp * PRES(k,i,j) * ( RHOT_RK(k,i,j) - RHOT(k,i,j) ) / RHOT(k,i,j)
       if ( abs( lhs - rhs ) / lhs > 1e-15 ) then
          write(*,*) "error is too large: ", k,i,j, lhs, rhs
          call abort
       end if
    end do
    end do
    end do
  end subroutine check_pres
#endif


end module scale_atmos_dyn_rk_hivi
