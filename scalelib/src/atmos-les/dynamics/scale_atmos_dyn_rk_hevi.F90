!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          Runge-Kutta for Atmospheric dynamical process
!!          HEVI, no terrain + tracer FCT limiter
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-18 (S.Nishizawa) [new] newly impremented
!! @li      2014-04-04 (S.Nishizawa) [mod] support terrain-following coordinate
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_rk_hevi
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

#ifdef DEBUG
  use scale_debug, only: &
       CHECK
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_rk_hevi_setup
  public :: ATMOS_DYN_rk_hevi

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
  integer, private, parameter :: NB = 1
  real(RP), private :: kappa
  !-----------------------------------------------------------------------------


contains

  subroutine ATMOS_DYN_rk_hevi_setup( ATMOS_TYPE_DYN )
    use scale_process, only: &
       PRC_MPIstop
#ifdef DRY
    use scale_const, only: &
       CVdry  => CONST_CVdry,  &
       CPdry  => CONST_CPdry
#endif
    implicit none

    character(len=H_SHORT), intent(in) :: ATMOS_TYPE_DYN
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HEVI'

    if ( ATMOS_TYPE_DYN .ne. 'HEVI' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_DYN is not HEVI. Check!'
       call PRC_MPIstop
    end if

#ifdef HEVI_BICGSTAB
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING Bi-CGSTAB'
#else
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING LAPACK'
#endif

#ifdef DRY
    kappa = CPdry / CVdry
#endif

    return
  end subroutine ATMOS_DYN_rk_hevi_setup


  subroutine ATMOS_DYN_rk_hevi( &
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
       UNDEF  => CONST_UNDEF,  &
       IUNDEF => CONST_UNDEF2, &
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
    real(RP) :: DPRES(KA,IA,JA) ! pressure deviation from reference pressure

    real(RP) :: qflx_hi(KA,IA,JA,3)
    real(RP) :: qflx_J (KA,IA,JA)

    ! for implicit solver
    real(RP) :: A(KA)
    real(RP) :: B(KA)
    real(RP) :: Sr(KA,IA,JA)
    real(RP) :: Sw(KA,IA,JA)
    real(RP) :: St(KA,IA,JA)
    real(RP) :: PT(KA)
    real(RP) :: CPRES(KA) ! kappa * PRES / RHOT
    real(RP) :: CPtot(KA,IA,JA)
    real(RP) :: C(KMAX-1)

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j

#ifdef DEBUG
    PRES(:,:,:) = UNDEF
    VELZ(:,:,:) = UNDEF
    VELX(:,:,:) = UNDEF
    VELY(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF
    DPRES(:,:,:) = UNDEF

    PT(:) = UNDEF
    CPRES(:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    mflx_hi(KS-1,:,:,ZDIR) = 0.0_RP

    qflx_hi(:,:,:,:) = UNDEF
    qflx_J (:,:,:)   = UNDEF
#endif

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
       ! pressure, pott. temp.
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
          CPtot(k,i,j) = CVtot(k,i,j) + Rtot(k,i,j)
       enddo
       enddo
       enddo
#endif
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, RHOT(k,i,j) )
             call CHECK( __LINE__, Rtot(k,i,j) )
             call CHECK( __LINE__, CVtot(k,i,j) )
#endif
#ifdef DRY
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rdry / P00 )**kappa
#else
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rtot(k,i,j) / P00 )**(CPtot(k,i,j)/CVtot(k,i,j))
#endif
          enddo
          PRES(KS-1,i,j) = PRES(KS+1,i,j) - DENS(KS,i,j) * ( PHI(KS-1,i,j) - PHI(KS+1,i,j) )
          PRES(KE+1,i,j) = PRES(KE-1,i,j) - DENS(KE,i,j) * ( PHI(KE+1,i,j) - PHI(KE-1,i,j) )

          do k = KS-1, KE+1
             DPRES(k,i,j) = PRES(k,i,j) - REF_pres(k,i,j)
          enddo
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
       do j = JJS, JJE+1
       do i = IIS, IIE+1
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
       do j = JJS, JJE+1
       do i = IIS, IIE+1
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
#endif
          mflx_hi(k,i,j,ZDIR) = J13G(k,i,j,I_XYW) * 0.25_RP * ( MOMX(k+1,i,j)+MOMX(k+1,i-1,j) &
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
          mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          mflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo


       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
          call CHECK( __LINE__, MOMX(k,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,XDIR) )
#endif
          mflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) &
               * ( MOMX(k,i,j) + num_diff(k,i,j,I_DENS,XDIR) )
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
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) &
               * ( MOMY(k,i,j) + num_diff(k,i,j,I_DENS,YDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

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
          Sr(k,i,j) =  &
                 - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i  ,j,  ZDIR) ) * RCDZ(k) &
                   + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                   + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) / GSQRT(k,i,j,I_XYZ) & ! divergence
                 + DENS_t(k,i,j)
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
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = mflx_hi(KS  ,i,j,ZDIR) * 0.5_RP * ( POTT(KS+1,i,j) + POTT(KS  ,i,j) ) &
                                 + GSQRT(KS,i,j,I_XYW) * num_diff(KS  ,i,j,I_RHOT,ZDIR)
          qflx_hi(KE-1,i,j,ZDIR) = mflx_hi(KE-1,i,j,ZDIR) * 0.5_RP * ( POTT(KE  ,i,j) + POTT(KE-1,i,j) ) &
                                 + GSQRT(KE-1,i,j,ZDIR) * num_diff(KE-1,i,j,I_RHOT,ZDIR)
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
                                  + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) ) &
                                + GSQRT(k,i,j,I_UYZ) * num_diff(k,i,j,I_RHOT,XDIR)
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
                                  + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) ) &
                                + GSQRT(k,i,j,I_XVZ) * num_diff(k,i,j,I_RHOT,YDIR)
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
                 + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) / GSQRT(k,i,j,I_XYZ) &
               + RHOT_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! implicit solver
       !$omp parallel do private(i,j,k,PT,A,B,C) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE

          do k = KS+1, KE-2
             PT(k) = FACT_N * ( POTT(k+1,i,j) + POTT(k  ,i,j) ) &
                   + FACT_F * ( POTT(k+2,i,j) + POTT(k-1,i,j) )
          end do
          PT(KS  ) = ( POTT(KS+1,i,j) + POTT(KS  ,i,j) ) * 0.5_RP
          PT(KE-1) = ( POTT(KE  ,i,j) + POTT(KE-1,i,j) ) * 0.5_RP

          do k = KS, KE
             CPRES(k) = &
#ifdef DRY
             kappa        * PRES(k,i,j) /                  RHOT(k,i,j)
#else
             CPtot(k,i,j) * PRES(k,i,j) / ( CVtot(k,i,j) * RHOT(k,i,j) )
#endif
          end do

          do k = KS, KE
             A(k) = dtrk**2 * J33G * RCDZ(k) * CPRES(k) * J33G / GSQRT(k,i,j,I_XYZ)
          end do
          do k = KS, KE
             B(k) = GRAV * dtrk**2 * J33G / ( CDZ(k+1) + CDZ(k) )
          end do

          ! vector
          do k = KS, KE-1
             C(k-KS+1) = MOMZ(k,i,j) &
                  + dtrk * ( &
                      - ( DPRES(k+1,i,j) + CPRES(k+1)*dtrk*St(k+1,i,j) &
                        - DPRES(k  ,i,j) - CPRES(k  )*dtrk*St(k  ,i,j) ) &
                      * RFDZ(k) * J33G / GSQRT(k,i,j,I_XYW) &
                      - 0.5_RP * GRAV &
                      * ( ( DENS(k+1,i,j) - REF_dens(k+1,i,j) ) &
                        + ( DENS(k  ,i,j) - REF_dens(k  ,i,j) ) &
                        + dtrk * ( Sr(k+1,i,j) + Sr(k,i,j) ) ) &
                      + Sw(k,i,j) )
          end do

#ifdef HEVI_BICGSTAB
          call solve_bicgstab( &
               C,                      & ! (inout)
               A, B, PT,               & ! (in)
               GSQRT(:,i,j,I_XYW), RFDZ) ! (in)
#else
          call solve_lapack( &
               C,                       & ! (inout)
               A, B, PT,                & ! (in)
#ifdef DEBUG
               i, j,                    & ! (in)
#endif
               GSQRT(:,i,j,I_XYW), RFDZ ) ! (in)
#endif

          ! z-momentum flux
          do k = KS, KE-1
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + J33G * C(k-KS+1)
          end do

          ! z-momentum
          do k = KS, KE-1
             MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                  + ( C(k-KS+1) - MOMZ(k,i,j) )
          end do

#ifdef DEBUG_HEVI2HEVE
          ! for debug (change to explicit integration)
          do k = KS, KE-1
             C(k-KS+1) = MOMZ(k,i,j)
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + J33G * MOMZ(k,i,j)
             MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                  + dtrk*( &
                  - J33G * ( DPRES(k+1,i,j)-DPRES(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,i_XYW) &
                  - 0.5_RP * GRAV * ( (DENS(k,i,j)-REF_dens(k,i,j))+(DENS(k+1,i,j)-REF_dens(k+1,i,j)) ) &
                  + Sw(k,i,j) )
          end do
#endif

          ! density
          do k = KS+1, KE-1
             DENS_RK(k,i,j) = DENS0(k,i,j) &
                  + dtrk * ( - J33G * ( C(k-KS+1) - C(k-KS) ) * RCDZ(k) / GSQRT(k,i,j,I_XYZ) &
                             + Sr(k,i,j) )
          end do
          DENS_RK(KS,i,j) = DENS0(KS,i,j) &
                  + dtrk * ( - J33G * C(1) * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) & ! C(0) = 0
                             + Sr(KS,i,j) )
          DENS_RK(KE,i,j) = DENS0(KE,i,j) &
                  + dtrk * (   J33G * C(KE-KS) * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) & ! C(KE-KS+1) = 0
                             + Sr(KE,i,j) )

          ! rho*theta
          do k = KS+1, KE-1
             RHOT_RK(k,i,j) = RHOT0(k,i,j) &
                  + dtrk * ( - ( C(k-KS+1) * PT(k) - C(k-KS) * PT(k-1) ) &
                               * J33G * RCDZ(k) / GSQRT(k,i,j,I_XYZ) &
                             + St(k,i,j) )
          end do
          RHOT_RK(KS,i,j) = RHOT0(KS,i,j) &
                  + dtrk * ( - C(1) * PT(KS) & ! C(0) = 0
                               * J33G * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) &
                             + St(KS,i,j) )
          RHOT_RK(KE,i,j) = RHOT0(KE,i,j) &
                  + dtrk * (   C(KE-KS) * PT(KE-1) & ! C(KE-KS+1) = 0
                               * J33G * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) &
                             + St(KE,i,j) )


#ifdef DEBUG
       call check_equation( &
            C, &
            DENS(:,i,j), MOMZ(:,i,j), RHOT(:,i,j), PRES(:,i,j), &
            Sr(:,i,j), Sw(:,i,j), St(:,i,j), &
            J33G, GSQRT(:,i,j,:), &
#ifdef DRY
            kappa, &
#else
            CPtot(:,i,j), CVtot(:,i,j), &
#endif
            dtrk, i, j )
#endif
       end do
       end do



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
       do j = JJS, JJE
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
       do j = JJS, JJE
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
       do j = JJS, JJE
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
       do i = IIS,   IIE
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
          call CHECK( __LINE__, PRES(k,i+1,j) )
          call CHECK( __LINE__, PRES(k,i  ,j) )
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
          MOMX_RK(k,i,j) = MOMX0(k,i,j) &
               + dtrk * ( ( - ( ( qflx_hi(k,i  ,j,ZDIR) - qflx_hi(k-1,i,j  ,ZDIR) ) * RCDZ(k) &
                              + ( qflx_J (k,i  ,j)      - qflx_J (k-1,i,j)        ) * RCDZ(k) &
                              + ( qflx_hi(k,i+1,j,XDIR) - qflx_hi(k  ,i,j  ,XDIR) ) * RFDX(i) &
                              + ( qflx_hi(k,i  ,j,YDIR) - qflx_hi(k  ,i,j-1,YDIR) ) * RCDY(j) ) &
                            - ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) &
                              - GSQRT(k,i  ,j,I_XYZ) * DPRES(k,i  ,j) ) * RFDX(i) &
                            - ( J13G(k+1,i,j,I_UYZ) * ( DPRES(k+1,i+1,j)+DPRES(k+1,i,j) ) &
                              - J13G(k-1,i,j,I_UYZ) * ( DPRES(k-1,i+1,j)+DPRES(k-1,i,j) ) ) &
                            * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) &  ! pressure gradient force
                          ) / GSQRT(k,i,j,I_UYZ) &
                          + 0.0625_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i+1,j) ) &
                                      * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                                      * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) & ! coriolis force
                          + divdmp_coef * dtrk * ( DDIV(k,i+1,j)-DDIV(k,i,j) ) * FDX(i) & ! divergence damping
                          + MOMX_t(k,i,j) )
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
       do i = IIS, IIE
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
       do i = IIS, IIE
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
       do i = IIS, IIE
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
       do j = JJS,   JJE
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
       do i = IIS, IIE
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
          call CHECK( __LINE__, PRES(k,i,j  ) )
          call CHECK( __LINE__, PRES(k,i,j+1) )
          call CHECK( __LINE__, CORIOLI(1,i,j  ) )
          call CHECK( __LINE__, CORIOLI(1,i,j+1) )
          call CHECK( __LINE__, VELX(k,i  ,j  ) )
          call CHECK( __LINE__, VELX(k,i  ,j+1) )
          call CHECK( __LINE__, VELX(k,i-1,j  ) )
          call CHECK( __LINE__, VELX(k,i-1,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j  ) )
          call CHECK( __LINE__, MOMY_t(k,i,j) )
          call CHECK( __LINE__, MOMY0(k,i,j) )
#endif
          MOMY_RK(k,i,j) = MOMY0(k,i,j) &
               + dtrk * ( ( - ( ( qflx_hi(k,i,j  ,ZDIR) - qflx_hi(k-1,i  ,j,ZDIR) ) * RCDZ(k) &
                              + ( qflx_J (k,i,j  )      - qflx_J (k-1,i  ,j)      ) * RCDZ(k) &
                              + ( qflx_hi(k,i,j  ,XDIR) - qflx_hi(k  ,i-1,j,XDIR) ) * RCDX(i) &
                              + ( qflx_hi(k,i,j+1,YDIR) - qflx_hi(k  ,i  ,j,YDIR) ) * RFDY(j) ) &
                            - ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) &
                              - GSQRT(k,i,j  ,I_XYZ) * DPRES(k,i,j  ) ) * RFDY(j) &
                            - ( J23G(k+1,i,j,I_XVZ) * ( DPRES(k+1,i,j+1)+DPRES(k+1,i,j) ) &
                              - J23G(k-1,i,j,I_XVZ) * ( DPRES(k-1,i,j+1)+DPRES(k-1,i,j) ) ) &
                            * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) & ! pressure gradient force
                          ) / GSQRT(k,i,j,I_XVZ) &
                          - 0.0625_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) &
                                      * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
                                      * ( VELX(k,i,j+1)+VELX(k,i,j)+VELX(k,i-1,j+1)+VELX(k,i-1,j) ) & ! coriolis force
                          + divdmp_coef * dtrk * ( DDIV(k,i,j+1)-DDIV(k,i,j) ) * FDY(j) & ! divergence damping
                          + MOMY_t(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

  end subroutine ATMOS_DYN_rk_hevi

#ifdef HEVI_BICGSTAB
  subroutine solve_bicgstab( &
       C,        & ! (inout)
       A, B, PT, & ! (in)
       G, RFDZ   ) ! (in)

    use scale_process, only: &
       PRC_MPIstop
    implicit none
    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: A(KA)
    real(RP), intent(in)    :: B(KA)
    real(RP), intent(in)    :: PT(KA)
    real(RP), intent(in)    :: G(KA)
    real(RP), intent(in)    :: RFDZ(KA-1)

    real(RP) :: r0(KMAX-1)

    real(RP) :: M(3,KMAX-1)
    real(RP) :: p(KMAX-1)
    real(RP) :: ap(KMAX-1)
    real(RP) :: s(KMAX-1)
    real(RP) :: as(KMAX-1)
    real(RP) :: al, be, w

    real(RP), pointer :: r(:)
    real(RP), pointer :: rn(:)
    real(RP), pointer :: swap(:)
    real(RP), target :: v0(KMAX-1)
    real(RP), target :: v1(KMAX-1)
    real(RP) :: r0r
    real(RP) :: norm, error, epsilon

    integer :: k, iter

    epsilon = 0.1_RP**(RP-1)

    call make_matrix(M,   & ! (out)
         A, B, PT, G, RFDZ) ! (in)

    norm = 0.0_RP
    do k = 1, KMAX-1
       norm = norm + C(k)**2
    end do

    r  => v0
    rn => v1

    call mul_matrix( v1, M, C )

    do k = 1, KMAX-1
       r(k) = C(k) - v1(k)
       r0(k) = r(k)
       p(k) = r(k)
    end do

    r0r = r0(1) * r(1)
    do k = 2, KMAX-1
       r0r = r0r + r0(k)*r(k)
    enddo
    do iter = 1, KMAX-1
       error = 0.0_RP
       do k = 1, KMAX-1
          error = error + r(k)**2
       end do

!       if ( error < epsilon .or. error / norm < epsilon ) then
       if ( error/norm < epsilon ) then
#ifdef DEBUG
!          write(*,*) "Bi-CGSTAB converged:", iter
#endif
          exit
       end if

       call mul_matrix( ap, M, p )
       al = r0(1) * ap(1)
       do k = 2, KMAX-1
          al = al + r0(k)*ap(k)
       end do
       al = r0r / al ! (r0,r) / (r0,Mp)
       s(:) = r(:) - al*ap(:)
       call mul_matrix( as, M, s )
       be = as(1) * s(1)  ! be is used as just work variable here
       w =  as(1) * as(1)
       do k = 2, KMAX-1
          be = be + as(k)*s(k)
          w  = w  + as(k)*as(k)
       end do
       w = be / w ! (as,s) / (as,as)

       c(:) = c(:) + al*p(:) + w*s(:)
       rn(:) = s(:) - w*as(:)
       be = al/w / r0r
       r0r = r0(1) * rn(1)
       do k = 2, KMAX-1
          r0r = r0r + r0(k)*rn(k)
       end do
       be = be * r0r ! al/w * (r0,rn)/(r0,r)
       p(:) = rn(:) + be * ( p(:) - w*ap(:) )

       swap => rn
       rn => r
       r => swap
    end do

    if ( iter >= KMAX-1 ) then
       write(*,*) 'xxx [atmos_dyn_hevi] Bi-CGSTAB'
       write(*,*) 'xxx not converged', error, norm
       call PRC_MPIstop
    end if

    return
  end subroutine solve_bicgstab

  subroutine make_matrix(M, &
       A, B, PT, &
       G, RFDZ)
    implicit none
    real(RP), intent(out) :: M(3, KMAX-1)
    real(RP), intent(in)  :: A(KA)
    real(RP), intent(in)  :: B(KA)
    real(RP), intent(in)  :: PT(KA)
    real(RP), intent(in)  :: G(KA)
    real(RP), intent(in)  :: RFDZ(KA-1)

    integer :: k

    ! k = KS
    M(3,1) =        - ( PT(KS+1)*RFDZ(KS)* A(KS+1) + B(KS) ) / G(KS)
    M(2,1) = 1.0_RP +   PT(KS  )*RFDZ(KS)*(A(KS+1) + A(KS) ) / G(KS)
    do k = KS+1, KE-2
       M(3,k-KS+1) =        - ( PT(k+1)*RFDZ(k)* A(k+1) + B(k) ) / G(k)
       M(2,k-KS+1) = 1.0_RP +   PT(k  )*RFDZ(k)*(A(k+1) + A(k) ) / G(k)
       M(1,k-KS+1) =        - ( PT(k-1)*RFDZ(k)* A(k  ) - B(k) ) / G(k)
    enddo
    ! k = KE-1
    M(2,KE-KS) = 1.0_RP +   PT(KE-1)*RFDZ(KE-1)*(A(KE  ) + A(KE-1) ) / G(KE-1)
    M(1,KE-KS) =        - ( PT(KE-2)*RFDZ(KE-1)* A(KE-1) - B(KE-1) ) / G(KE-1)

    return
  end subroutine make_matrix

  subroutine mul_matrix(V, M, C)
    implicit none
    real(RP), intent(out) :: V(KMAX-1)
    real(RP), intent(in)  :: M(3, KMAX-1)
    real(RP), intent(in)  :: C(KMAX-1)

    integer :: k

    ! k = KS
    V(1) = M(3,1)*C(2) + M(2,1)*C(1)
    do k = 2, KMAX-2
       V(k) = M(3,k)*C(k+1) + M(2,k)*C(k) + M(1,k)*C(k-1)
    enddo
    ! k = KE-1
    V(KMAX-1) = M(2,KMAX-1)*C(KMAX-1) + M(1,KMAX-1)*C(KMAX-2)

    return
  end subroutine mul_matrix

#else

  subroutine solve_lapack( &
       C,   & ! (inout)
       A, B, PT, &
#ifdef DEBUG
       i, j, &
#endif
       G, RFDZ ) ! (in)
    use scale_process, only: &
         PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: A(KA)
    real(RP), intent(in)    :: B(KA)
    real(RP), intent(in)    :: PT(KA)
#ifdef DEBUG
    integer , intent(in)    :: i
    integer , intent(in)    :: j
#endif
    real(RP), intent(in)    :: G(KA)
    real(RP), intent(in)    :: RFDZ(KA-1)

    real(RP) :: M(NB*3+1,KMAX-1)
    integer  :: IPIV(KMAX-1)
    integer  :: INFO

    integer :: k

#ifdef DEBUG
    real(RP) :: M2(KMAX-1,KMAX-1)
    real(RP) :: C2(KMAX-1)
    real(RP) :: sum
#endif


    ! band matrix
    do k = KS+1, KE-2
       ! k-1, +1
       M(NB+1,k-KS+1) =        - ( PT(k)*RFDZ(k-1)*  A(k  ) + B(k-1) ) / G(k-1)
       ! k, 0
       M(NB+2,k-KS+1) = 1.0_RP +   PT(k)*RFDZ(k  )*( A(k+1) + A(k  ) ) / G(k)
       ! k+1, -1
       M(NB+3,k-KS+1) =        - ( PT(k)*RFDZ(k+1)*  A(k+1) - B(k+1) ) / G(k+1)
    end do
    ! KS, 0
    M(NB+2,1    ) = 1.0_RP +   PT(KS  )*RFDZ(KS  )*( A(KS+1) + A(KS  ) ) / G(KS)
    ! KS+1, -1
    M(NB+3,1    ) =        - ( PT(KS  )*RFDZ(KS+1)*  A(KS+1) - B(KS+1) ) / G(KS+1)
    ! KE-2, +1
    M(NB+1,KE-KS) =        - ( PT(KE-1)*RFDZ(KE-2)*  A(KE-1) + B(KE-2) ) / G(KE-2)
    ! KE-1, 0
    M(NB+2,KE-KS) = 1.0_RP +   PT(KE-1)*RFDZ(KE-1)*( A(KE  ) + A(KE-1) ) / G(KE-1)

#ifdef DEBUG
          M2(:,:) = 0.0_RP
          do k = 1, KMAX-1
             if (k>1) M2(k-1,k) = M(NB+3,k-1)
             M2(k,k) = M(NB+2,k)
             if (k<kmax-1) M2(k+1,k) = M(NB+1,k+1)
             c2(k) = c(k)
          end do
#endif
          if ( RP == DP ) then
             call DGBSV( KMAX-1, NB, NB, 1, M, NB*3+1, IPIV, C, KMAX-1, INFO)
          else
             call SGBSV( KMAX-1, NB, NB, 1, M, NB*3+1, IPIV, C, KMAX-1, INFO)
          end if
          ! C is (\rho w)^{n+1}
#ifdef DEBUG
          if ( INFO .ne. 0 ) then
             write(*,*) "DGBSV was failed", info
             call PRC_MPIstop
          end if

          do k = 1, KMAX-1
             sum = 0.0_RP
             if (k>1) sum = sum + M2(k-1,k)*c(k-1)
             sum = sum + M2(k,k)*c(k)
             if (k<kmax-1) sum = sum + M2(k+1,k)*c(k+1)
             if ( abs(sum-c2(k)) > 1E-10_RP ) then
                write(*,*) "sum is different"
                write(*,*) k+2, i, j, sum, c2(k)
!                write(*,*) M2(k-1:k+1,k), c(k-1:k+1)
                call PRC_MPIstop
             end if
          end do
#endif
  end subroutine solve_lapack
#endif

#ifdef DEBUG
  subroutine check_equation( &
       VECT, &
       DENS, MOMZ, RHOT, PRES, &
       Sr, Sw, St, &
       J33G, G, &
#ifdef DRY
       kappa, &
#else
       CPtot, CVtot, &
#endif
       dt, i, j )
    use scale_const, only: &
         EPS => CONST_EPS, &
         GRAV => CONST_GRAV
    use scale_process, only: &
         PRC_MPIstop
    use scale_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW
    implicit none
    real(RP), intent(in) :: VECT(KMAX-1)
    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: MOMZ(KA)
    real(RP), intent(in) :: RHOT(KA)
    real(RP), intent(in) :: PRES(KA)
    real(RP), intent(in) :: Sr(KA)
    real(RP), intent(in) :: Sw(KA)
    real(RP), intent(in) :: St(KA)
    real(RP), intent(in) :: J33G
    real(RP), intent(in) :: G(KA,8)
#ifdef DRY
    real(RP), intent(in) :: kappa
#else
    real(RP), intent(in) :: CPtot(KA)
    real(RP), intent(in) :: CVtot(KA)
#endif
    real(RP), intent(in) :: dt
    integer , intent(in) :: i
    integer , intent(in) :: j

    real(RP), parameter :: small = 1e-6_RP

    real(RP) :: MOMZ_N(KA)
    real(RP) :: DENS_N(KA)
    real(RP) :: RHOT_N(KA)
    real(RP) :: PRES_N(KA)

    real(RP) :: POTT(KA)
    real(RP) :: PT(KA)

#ifndef DRY
    real(RP) :: kappa
#endif

    real(RP) :: error, lhs, rhs
    real(RP) :: a0, a1, b
    integer :: k


    do k = KS, KE-1
       MOMZ_N(k) = VECT(k-KS+1)
    end do
    MOMZ_N(:KS-1) = 0.0_RP
    MOMZ_N(KE:) = 0.0_RP

    ! density
    do k = KS+1, KE-1
       DENS_N(k) = DENS(k) &
            + dt * ( - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k) )
    end do
    DENS_N(KS) = DENS(KS) &
         + dt * ( - J33G * MOMZ_N(KS) * RCDZ(KS) / G(KS,I_XYZ) + Sr(KS) )
    DENS_N(KE) = DENS(KE) &
         + dt * ( J33G * MOMZ_N(KE-1) * RCDZ(KE) / G(KE,I_XYZ) + Sr(KE) )

    ! rho*theta
    do k = KS, KE
       POTT(k) = RHOT(k) / DENS(k)
    end do
    do k = KS+1, KE-2
       PT(k) = ( 7.0_RP * ( POTT(k+1) + POTT(k  ) ) &
                 -        ( POTT(k+2) + POTT(k-1) ) ) / 12.0_RP
    end do
    PT(KS-1) = 0.0_RP
    PT(KS  ) = ( POTT(KS+1) + POTT(KS  ) ) * 0.5_RP
    PT(KE-1) = ( POTT(KE  ) + POTT(KE-1) ) * 0.5_RP
    PT(KE  ) = 0.0_RP
    do k = KS+1, KE-1
       RHOT_N(k) = RHOT(k) &
            + dt * ( - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) &
                     + St(k) )
    end do
    RHOT_N(KS) = RHOT(KS) &
         + dt * ( - J33G * MOMZ_N(KS)*PT(KS) * RCDZ(KS) / G(KS,I_XYZ) + St(KS) )
    RHOT_N(KE) = RHOT(KE) &
         + dt * ( J33G * MOMZ_N(KE-1)*PT(KE-1) * RCDZ(KE) / G(KE-1,I_XYZ) + St(KE) )


    do k = KS, KE
#ifndef DRY
       kappa = CPtot(k) / CVtot(k)
#endif
       PRES_N(k) = PRES(k) * ( 1.0_RP + kappa * ( RHOT_N(k) - RHOT(k) ) / RHOT(k) )
    end do

    do k = KS, KE
       lhs = ( DENS_N(k) - DENS(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       end if
       if ( abs(error) > small ) then
          write(*,*)"HEVI: DENS error", k, i, j, error, lhs, rhs
          write(*,*)eps
          call PRC_MPIstop
       end if
    end do

    do k = KS, KE-1
       lhs = ( MOMZ_N(k) - MOMZ(k) ) / dt
       rhs = - J33G * ( PRES_N(k+1) - PRES_N(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( DENS_N(k+1) + DENS_N(k) ) * 0.5_RP &
             + Sw(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       end if
       if ( abs(error) > small ) then
          write(*,*)"HEVI: MOMZ error", k, i, j, error, lhs, rhs
          write(*,*) MOMZ_N(k), MOMZ(k), dt
          write(*,*) - J33G * ( PRES(k+1) - PRES(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( DENS(k+1) + DENS(k) ) * 0.5_RP &
             + Sw(k)
          call PRC_MPIstop
       end if
    end do

    do k = KS, KE
       lhs = ( RHOT_N(k) - RHOT(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) + St(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       end if
       if ( abs(error) > small ) then
          write(*,*)"HEVI: RHOT error", k, i, j, error, lhs, rhs
          call PRC_MPIstop
       end if
    end do

    return
  end subroutine check_equation
#endif

end module scale_atmos_dyn_rk_hevi
