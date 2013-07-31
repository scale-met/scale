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
!!
!<
!-------------------------------------------------------------------------------
#define OMP_SCHEDULE_ schedule(static)
module mod_atmos_dyn_rk
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
#ifdef DEBUG
  use mod_debug, only: &
       CHECK
#endif
  use dc_types, only: &
       DP
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_rk_setup
  public :: ATMOS_DYN_rk

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  integer, private, parameter :: NB = 2
  real(RP), private :: kappa
  !-----------------------------------------------------------------------------


contains

  subroutine ATMOS_DYN_rk_setup
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_TYPE_DYN
#ifdef DRY
    use mod_const, only : &
       CVdry  => CONST_CVdry,  &
       CPdry  => CONST_CPdry
#endif

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
  end subroutine ATMOS_DYN_rk_setup


#define __BOUND 1

  subroutine ATMOS_DYN_rk( &
    DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
    DDIV,                                        &
    mflx_hi,                                     &
    DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
    DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
    DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
    Rtot, CVtot, CORIOLI,                        &
    num_diff, divdmp_coef,                       &
    FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, &
    CDZ, FDZ, FDX, FDY,                          &
    RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
    dtrk, rk, rko,                               &
    VELZ, VELX, VELY, PRES, POTT                 )
    use mod_const, only : &
       UNDEF  => CONST_UNDEF,  &
       IUNDEF => CONST_UNDEF2, &
       GRAV   => CONST_GRAV,   &
       P00    => CONST_PRE00,  &
       Rdry   => CONST_Rdry,   &
       CVdry  => CONST_CVdry,  &
       CPdry  => CONST_CPdry
    use mod_comm, only: &
#ifdef _USE_RDMA
       COMM_rdma_vars8, &
#endif
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
       ZDIR, &
       XDIR, &
       YDIR, &
       I_DENS, &
       I_MOMZ, &
       I_MOMX, &
       I_MOMY, &
       I_RHOT
    use mod_atmos_dyn_common, only: &
       FACT_N, &
       FACT_F, &
       ATMOS_DYN_fct
    implicit none

    real(RP), intent(inout) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(inout) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(inout) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(inout) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(inout) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(out)   :: DDIV(KA,IA,JA)      ! divergence

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! rho * vel(x,y,z)

    real(RP), intent(inout),target :: DENS0(KA,IA,JA)   ! prognostic variables
    real(RP), intent(inout),target :: MOMZ0(KA,IA,JA)   ! at previous dynamical time step
    real(RP), intent(inout),target :: MOMX0(KA,IA,JA)   !
    real(RP), intent(inout),target :: MOMY0(KA,IA,JA)   !
    real(RP), intent(inout),target :: RHOT0(KA,IA,JA)   !

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
    real(RP), intent(in) :: dtrk
    integer , intent(in) :: rk
    integer , intent(in) :: rko

    ! diagnostic variables (work space)
    real(RP), intent(out) :: PRES(KA,IA,JA) ! pressure [Pa]
    real(RP), intent(out) :: VELZ(KA,IA,JA) ! velocity w [m/s]
    real(RP), intent(out) :: VELX(KA,IA,JA) ! velocity u [m/s]
    real(RP), intent(out) :: VELY(KA,IA,JA) ! velocity v [m/s]
    real(RP), intent(out) :: POTT(KA,IA,JA) ! potential temperature [K]

    ! for implicit solver
    real(RP) :: A(KA)
    real(RP) :: B
    real(RP) :: PT(KA)
    real(RP) :: Sr(KA,IA,JA)
    real(RP) :: Sw(KA,IA,JA)
    real(RP) :: St(KA,IA,JA)
    real(RP) :: CPtot(KA,IA,JA)
    real(RP) :: C(KMAX-1)

    real(RP) :: qflx_hi(KA,IA,JA,3)
    real(RP) :: vel
    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j

#ifdef DEBUG
    PRES(:,:,:) = UNDEF
    VELZ(:,:,:) = UNDEF
    VELX(:,:,:) = UNDEF
    VELY(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    mflx_hi(KS-1,:,:,ZDIR) = 0.0_RP

    qflx_hi(:,:,:,:) = UNDEF
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
       do k = KS, KE
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

       !##### continuity equation #####

       ! at (x, y, interface)
       ! at (u, y, layer)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMX(k,i+1,j) )
          call CHECK( __LINE__, MOMX(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i-1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,XDIR) )
#endif
!          mflx_hi(k,i,j,XDIR) = ( -MOMX(k,i+1,j) + 8.0_RP*MOMX(k,i,j) -MOMX(k,i-1,j) ) / 6.0_RP &
          mflx_hi(k,i,j,XDIR) = MOMX(k,i,j) &
                              + num_diff(k,i,j,I_DENS,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, layer)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMY(k,i,j+1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
!          mflx_hi(k,i,j,YDIR) = ( -MOMY(k,i,j+1) + 8.0_RP*MOMY(k,i,j) -MOMY(k,i,j-1) ) / 6.0_RP &
          mflx_hi(k,i,j,YDIR) = MOMY(k,i,j) &
                              + num_diff(k,i,j,I_DENS,YDIR)
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
                 - ( ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                   + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) & ! divergence
                 + DENS_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (z) #####

       ! high order flux
       ! at (x, y, layer)
       ! note than z-index is added by -1
       !$omp parallel do private(i,j,k,vel) OMP_SCHEDULE_ collapse(2)
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
          vel = 0.5_RP * ( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) / DENS(k,i,j)
          qflx_hi(k-1,i,j,ZDIR) = vel  &
                              * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                                + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMZ,ZDIR)
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
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          ! k = KS+1
          vel = 0.5_RP * ( MOMZ(KS+1,i,j)+MOMZ(KS,i,j) ) / DENS(KS+1,i,j)
          qflx_hi(KS,i,j,ZDIR) = vel  &
                              * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS,i,j) ) &
                                + FACT_F * ( MOMZ(KS+2,i,j)            ) ) &
                              + num_diff(KS+1,i,j,I_MOMZ,ZDIR)
          ! k = KE-1
          vel = 0.5_RP * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) / DENS(KE-1,i,j)
          qflx_hi(KE-2,i,j,ZDIR) = vel  &
                              * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                                + FACT_F * (                MOMZ(KE-3,i,j) ) ) &
                              + num_diff(KE-1,i,j,I_MOMZ,ZDIR)
          ! k = KE
          qflx_hi(KE-1,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, interface)
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
          qflx_hi(k,i,j,XDIR) = 0.5_RP * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMZ(k,i+1,j)+MOMZ(k,i  ,j) ) &
                                + FACT_F * ( MOMZ(k,i+2,j)+MOMZ(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_MOMZ,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
          qflx_hi(KE,i,j,XDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, interface)
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
          qflx_hi(k,i,j,YDIR) = 0.5_RP * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMZ(k,i,j+1)+MOMZ(k,i,j  ) ) &
                                + FACT_F * ( MOMZ(k,i,j+2)+MOMZ(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_MOMZ,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
          qflx_hi(KE,i,j,YDIR) = 0.0_RP
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
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, PRES(k  ,i,j) )
          call CHECK( __LINE__, PRES(k+1,i,j) )
          call CHECK( __LINE__, DENS(k  ,i,j) )
          call CHECK( __LINE__, DENS(k+1,i,j) )
          call CHECK( __LINE__, DDIV(k  ,i,j) )
          call CHECK( __LINE__, DDIV(k+1,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, MOMZ_t(k,i,j) )
#endif
          Sw(k,i,j) = &
                        ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
                            + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                          + divdmp_coef * dtrk  * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) & ! divergence damping
                          + MOMZ_t(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif



       !##### Thermodynamic Equation #####

       ! at (u, y, layer)
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
                              + num_diff(k,i,j,I_RHOT,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, layer)
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
                              + num_diff(k,i,j,I_RHOT,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update rho*theta
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, RHOT_t(k,i,j) )
          call CHECK( __LINE__, RHOT0(k,i,j) )
#endif
          St(k,i,j) = &
               - ( ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                 + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
               + RHOT_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! implicit solver
       B = GRAV * dtrk**2 / 12.0_RP
       !$omp parallel do private(i,j,k,PT,A,C) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE

          do k = KS+1, KE-2
             PT(k) = ( - ( POTT(k+2,i,j) + POTT(k-1,i,j) ) &
                       + ( POTT(k+1,i,j) + POTT(k  ,i,j) ) * 7.0_RP ) / 12.0_RP
          end do
          PT(KS  ) = ( POTT(KS+1,i,j) + POTT(KS  ,i,j) ) * 0.5_RP
          PT(KE-1) = ( POTT(KE  ,i,j) + POTT(KE-1,i,j) ) * 0.5_RP

          do k = KS, KE
#ifdef DRY
             A(k) = CPdry * dtrk**2 * PRES(k,i,j) &
                  / ( 6.0_RP * CVdry * CDZ(k) * RHOT(k,i,j) )
#else
             A(k) = CPtot(k,i,j) * dtrk**2 * PRES(k,i,j) &
                  / ( 6.0_RP * CVtot(k,i,j) * CDZ(k) * RHOT(k,i,j) )
#endif
          end do

          ! vector
          do k = KS, KE-1
             C(k-KS+1) = MOMZ(k,i,j) &
                  + dtrk * ( &
                      - RFDZ(k) * ( &
#ifdef DRY
                          PRES(k+1,i,j)*( 1.0_RP + dtrk*kappa*St(k+1,i,j)/RHOT(k+1,i,j) ) &
                        - PRES(k  ,i,j)*( 1.0_RP + dtrk*kappa*St(k  ,i,j)/RHOT(k  ,i,j) ) )  &
#else
                          PRES(k+1,i,j)*( 1.0_RP + dtrk*CPtot(k+1,i,j)*St(k+1,i,j)/(CVtot(k+1,i,j)*RHOT(k+1,i,j)) ) &
                        - PRES(k  ,i,j)*( 1.0_RP + dtrk*CPtot(k  ,i,j)*St(k  ,i,j)/(CVtot(k  ,i,j)*RHOT(k  ,i,j)) ) )  &
#endif
                      - 0.5_RP * GRAV * ( DENS(k+1,i,j) + DENS(k,i,j) + dtrk * ( Sr(k+1,i,j) + Sr(k,i,j) ) ) &
                      + Sw(k,i,j) )
          end do

#ifdef HEVI_BICGSTAB
          call solve_bicgstab( &
               C,         & ! (out)
               PT, A, B ) ! (in)
#else
          call solve_lapack( &
               C,         & ! (out)
               PT, A, B ) ! (in)
#endif

          ! z momentum flux
          do k = KS+1, KE-2
             mflx_hi(k,i,j,ZDIR) = ( - C(k-KS+2) + 8.0_RP * C(k-KS+1) - C(k-KS) ) / 6.0_RP
          end do
#if __BOUND
          mflx_hi(KS  ,i,j,ZDIR) = ( - C(2) + 8.0_RP * C(1) ) / 6.0_RP
          mflx_hi(KE-1,i,j,ZDIR) = ( 8.0_RP * C(KE-KS) - C(KE-KS-1) ) / 6.0_RP
#else
          mflx_hi(KS,i,j,ZDIR) = C(1)
          mflx_hi(KE-1,i,j,ZDIR) = C(KE-KS)
#endif

          ! z-momentum
          do k = KS, KE-1
             MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                  + ( C(k-KS+1) - MOMZ(k,i,j) )
          end do

          ! density
          do k = KS+1, KE-1
             DENS_RK(k,i,j) = DENS0(k,i,j) &
                  + dtrk * ( - ( mflx_hi(k,i,j,ZDIR) - mflx_hi(k-1,i,j,ZDIR) ) * RCDZ(k) &
                             + Sr(k,i,j) )
          end do
          DENS_RK(KS,i,j) = DENS0(KS,i,j) &
                  + dtrk * ( - mflx_hi(KS,i,j,ZDIR) * RCDZ(KS) + Sr(KS,i,j) )
          DENS_RK(KE,i,j) = DENS0(KE,i,j) &
                  + dtrk * ( mflx_hi(KE-1,i,j,ZDIR) * RCDZ(KE) + Sr(KE,i,j) )

          ! rho*theta
          do k = KS+1, KE-1
             RHOT_RK(k,i,j) = RHOT0(k,i,j) &
                  + dtrk * ( - ( mflx_hi(k,i,j,ZDIR) * PT(k) - mflx_hi(k-1,i,j,ZDIR) * PT(k-1) ) * RCDZ(k) &
                             + St(k,i,j) )
          end do
          RHOT_RK(KS,i,j) = RHOT0(KS,i,j) &
                  + dtrk * ( - mflx_hi(KS,i,j,ZDIR) * PT(KS) * RCDZ(KS) &
                             + St(KS,i,j) )
          RHOT_RK(KE,i,j) = RHOT0(KE,i,j) &
                  + dtrk * ( mflx_hi(KE-1,i,j,ZDIR) * PT(KE-1) * RCDZ(KE) &
                             + St(KE,i,j) )


#ifdef DEBUG
       call check_equation( &
            C, &
            DENS(:,i,j), MOMZ(:,i,j), RHOT(:,i,j), PRES(:,i,j), &
            Sr(:,i,j), Sw(:,i,j), St(:,i,j), &
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

       ! high order flux
       ! at (u, y, interface)
       !$omp parallel do private(i,j,k,vel) OMP_SCHEDULE_ collapse(2)
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
          vel = 0.5_RP * ( VELZ(k,i+1,j) + VELZ(k,i,j) )
          qflx_hi(k,i,j,ZDIR) = vel &
                              * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                                + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMX,ZDIR)
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
          ! just above the bottom boundary
          vel = 0.5_RP * ( VELZ(KS,i+1,j) + VELZ(KS,i,j) )
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * vel &
                                 * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR)
          ! just below the top boundary
          vel = 0.5_RP * ( VELZ(KE-1,i+1,j) + VELZ(KE-1,i,j) )
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * vel &
                                 * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_MOMX,ZDIR)
          qflx_hi(KE,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, y, layer)
       ! note that x-index is added by -1
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
          qflx_hi(k,i-1,j,XDIR) = 0.5_RP * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                              * ( FACT_N * ( MOMX(k,i  ,j)+MOMX(k,i-1,j) ) &
                                + FACT_F * ( MOMX(k,i+1,j)+MOMX(k,i-2,j) ) ) &
                              + num_diff(k,i,j,I_MOMX,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, v, layer)
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
          qflx_hi(k,i,j,YDIR) = 0.5_RP * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                              * ( FACT_N * ( MOMX(k,i,j+1)+MOMX(k,i,j  ) ) &
                                + FACT_F * ( MOMX(k,i,j+2)+MOMX(k,i,j-1) ) ) &
                              + num_diff(k,i,j,I_MOMX,YDIR)
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
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
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
                          + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                                       + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                                     - ( PRES(k,i+1,j)-PRES(k,i,j) ) * RFDX(i)                     & ! pressure gradient force
                                     + 0.125_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i+1,j) )              &
                                     * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) & ! coriolis force
                                     + divdmp_coef * dtrk * ( DDIV(k,i+1,j)-DDIV(k,i,j) ) * FDX(i) & ! divergence damping
                                     + MOMX_t(k,i,j)                                               )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (y) #####
       ! high order flux
       ! at (x, v, interface)
       !$omp parallel do private(i,j,k,vel) OMP_SCHEDULE_ collapse(2)
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
          vel = 0.5_RP * ( VELZ(k,i,j+1) + VELZ(k,i,j) )
          qflx_hi(k,i,j,ZDIR) = vel &
                              * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                                + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_MOMY,ZDIR)
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
          ! just above the bottom boundary
          vel = 0.5_RP * ( VELZ(KS,i,j+1) + VELZ(KS,i,j) )
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * vel * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR)
          ! just below the top boundary
          vel = 0.5_RP * ( VELZ(KE-1,i,j+1) + VELZ(KE-1,i,j) )
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * vel * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR)
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, v, layer)
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
          qflx_hi(k,i,j,XDIR) = 0.5_RP * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                              * ( FACT_N * ( MOMY(k,i+1,j)+MOMY(k,i  ,j) ) &
                                + FACT_F * ( MOMY(k,i+2,j)+MOMY(k,i-1,j) ) ) &
                              + num_diff(k,i,j,I_MOMY,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, y, layer)
       ! note that y-index is added by -1
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
          qflx_hi(k,i,j-1,YDIR) = 0.5_RP * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                              * ( FACT_N * ( MOMY(k,i,j  )+MOMY(k,i,j-1) ) &
                                + FACT_F * ( MOMY(k,i,j+1)+MOMY(k,i,j-2) ) ) &
                              + num_diff(k,i,j,I_MOMY,YDIR)
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
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
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
                          + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                       + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RFDY(j) ) &
                                     - ( PRES(k,i,j+1)-PRES(k,i,j) ) * RFDY(j)                     & ! pressure gradient force
                                     - 0.125_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i,j+1) )              &
                                     * ( VELX(k,i,j)+VELX(k,i,j+1)+VELX(k,i-1,j)+VELX(k,i-1,j+1) ) & ! coriolis force
                                     + divdmp_coef * dtrk * ( DDIV(k,i,j+1)-DDIV(k,i,j) ) * FDY(j) & ! divergence damping
                                     + MOMY_t(k,i,j)                                               )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

  end subroutine ATMOS_DYN_rk

#ifdef HEVI_BICGSTAB
  subroutine solve_bicgstab( &
       C,       & ! (inout)
       PT, A, B ) ! (in)
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    implicit none
    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: PT(KA)
    real(RP), intent(in)    :: A(KA)
    real(RP), intent(in)    :: B

    real(RP) :: r0(KMAX-1)

    real(RP) :: M(5,KMAX-1)
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

    call make_matrix(M, PT, A, B)

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
          write(*,*) "Bi-CGSTAB converged:", iter
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

  subroutine make_matrix(M, PT, A, B)
    use mod_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    implicit none
    real(RP), intent(out) :: M(5, KMAX-1)
    real(RP), intent(in)  :: PT(KA)
    real(RP), intent(in)  :: A(KA)
    real(RP), intent(in)  :: B

    integer :: k

#ifdef __BOUND
    ! k = KS
    M(5,KS-2) = A(KS+1)*pt(KS+1)*rfdz(KS) + B *rcdz(KS+1)
    M(4,KS-2) = -(A(KS+1)*(8.0_RP*pt(KS+1)+pt(KS))+A(KS)*pt(KS))*rfdz(KS)-B*(9.0_RP*rcdz(KS+1)-rcdz(KS))
    M(3,KS-2) = (A(KS+1)*(pt(KS+1)+8.0_RP*pt(KS))+A(KS)*(8.0_RP*pt(KS)))*rfdz(KS) + B*(9.0_RP*rcdz(KS+1)-8.0_RP*rcdz(KS)) + 1.0_RP
    ! k = KS+1
    M(5,KS-1) = A(KS+2)*pt(KS+2)*rfdz(KS+1) + B *rcdz(KS+2)
    M(4,KS-1) = -(A(KS+2)*(8.0_RP*pt(KS+2)+pt(KS+1))+A(KS+1)*pt(KS+1))*rfdz(KS+1)-B*(9.0_RP*rcdz(KS+2)-rcdz(KS+1))
    M(3,KS-1) = (A(KS+2)*(pt(KS+2)+8.0_RP*pt(KS+1))+A(KS+1)*(8.0_RP*pt(KS+1)+pt(KS)))*rfdz(KS+1) &
         + 9.0_RP*B*(rcdz(KS+2)-rcdz(KS+1)) + 1.0_RP
    M(2,KS-1) = -(A(KS+2)*pt(KS+1)+a(KS+1)*(pt(KS+1)+8.0_RP*pt(KS)))*rfdz(KS+1) -B*(rcdz(KS+2)-9.0_RP*rcdz(KS+1))
#else
    ! k = KS
    M(5,KS-2) =  A(KS+1)* pt(KS+1)                                    *rfdz(KS) &
         + B        *rcdz(KS+1)
    M(4,KS-2) = -A(KS+1)* pt(KS+1)                                    *rfdz(KS) &
         - B* 8.0_RP*rcdz(KS+1)
    M(3,KS-2) = (A(KS+1)*(pt(KS+1)+6.0_RP*pt(KS))+6.0_RP*A(KS)*pt(KS))*rfdz(KS) &
         + B*(7.0_RP*rcdz(KS+1)-6.0_RP*rcdz(KS)) + 1.0_RP
    ! k = KS+1
    M(5,KS-1) =   A(KS+2)*        pt(KS+2)                                   *rfdz(KS+1) &
         + B*rcdz(KS+2)
    M(4,KS-1) = -(A(KS+2)*(8.0_RP*pt(KS+2)+pt(KS+1))+       A(KS+1)*pt(KS+1)                )*rfdz(KS+1) &
         - B*(9.0_RP*rcdz(KS+2)-       rcdz(KS+1))
    M(3,KS-1) =  (A(KS+2)*(pt(KS+2)+8.0_RP*pt(KS+1))+8.0_RP*A(KS+1)*pt(KS+1)                )*rfdz(KS+1) &
         + B*(9.0_RP*rcdz(KS+2)-8.0_RP*rcdz(KS+1)) + 1.0_RP
    M(2,KS-1) = -(A(KS+2)* pt(KS+1)                 +       A(KS+1)*(pt(KS+1)+6.0_RP*pt(KS)))*rfdz(KS+1) &
         - B*(       rcdz(KS+2)-7.0_RP*rcdz(KS+1))
    end if
#endif
    do k = KS+2, KE-3
       M(5,k-2) =   A(k+1)*        pt(k+1)                                           *rfdz(k) &
            + B*        rcdz(k+1)
       M(4,k-2) = -(A(k+1)*(8.0_RP*pt(k+1)+       pt(k))+A(k)        *pt(k)         )*rfdz(k) &
            - B*(9.0_RP*rcdz(k+1)-rcdz(k))
       M(3,k-2) =  (A(k+1)*(       pt(k+1)+8.0_RP*pt(k))+A(k)*(8.0_RP*pt(k)+pt(k-1)))*rfdz(k) &
            + B*9.0_RP*(rcdz(k+1)-rcdz(k)) + 1.0_RP
       M(2,k-2) = -(A(k+1)*        pt(k)                +A(k)*(pt(k)+8.0_RP*pt(k-1)))*rfdz(k) &
            - B*       (rcdz(k+1)-9.0_RP*rcdz(k))
       M(1,k-2) =                                        A(k)*pt(k-1)                *rfdz(k) &
            - B*                         rcdz(k)
    enddo
#ifdef __BOUND
    ! k = KE-2
    M(4,KE-4) = -(A(KE-1)*(8.0_RP*pt(KE-1)+pt(KE-2))+A(KE-2)*        pt(KE-2)          )*rfdz(KE-2) &
         - B*(9.0_RP*rcdz(KE-1)-rcdz(KE-2))
    M(3,KE-4) =  (A(KE-1)*(pt(KE-1)+8.0_RP*pt(KE-2))+A(KE-2)*(8.0_RP*pt(KE-2)+pt(KE-3)))*rfdz(KE-2) &
         + B*9.0_RP*(rcdz(KE-1)-rcdz(KE-2)) + 1.0_RP
    M(2,KE-4) = -(A(KE-1)*pt(KE-2)                  +A(KE-2)*(pt(KE-2)+8.0_RP*pt(KE-3)))*rfdz(KE-2) &
         - B*       (rcdz(KE-1)-9.0_RP*rcdz(KE-2))
    M(1,KE-4) =                                      A(KE-2)*                 pt(KE-3  )*rfdz(KE-2) &
         - B*                          rcdz(KE-2)
    ! k = KE-1
    M(3,KE-3) =  (A(KE)*(8.0_RP*pt(KE-1))+A(KE-1)*(8.0_RP*pt(KE-1)+pt(KE-2)))*rfdz(KE-1) &
           + B*(8.0_RP*rcdz(KE)-9.0_RP*rcdz(KE-1)) + 1.0_RP
    M(2,KE-3) = -(A(KE)*pt(KE-1)         +A(KE-1)*(pt(KE-1)+8.0_RP*pt(KE-2)))*rfdz(KE-1) &
           - B*(       rcdz(KE)-9.0_RP*rcdz(KE-1))
    M(1,KE-3) =                           A(KE-1)*                 pt(KE-2)  *rfdz(KE-1) &
           - B*                        rcdz(KE-1)
#else
    ! k = KE-2
    M(4,KE-2) = -(A(KE-1)*(6.0_RP*pt(KE-1)+pt(KE-2))+A(KE-2)*        pt(KE-2)          )*rfdz(KE-2) &
         - B*(7.0_RP*rcdz(KE-1)       -rcdz(KE-2))
    M(3,KE-2) =  (A(KE-1)* 8.0_RP*pt(KE-2)          +A(KE-2)*(8.0_RP*pt(KE-2)+pt(KE-3)))*rfdz(KE-2) &
         + B*(8.0_RP*rcdz(KE-1)-9.0_RP*rcdz(KE-2)) + 1.0_RP
    M(2,KE-2) = -(A(KE-1)*        pt(KE-2)          +A(KE-2)*(pt(KE-2)+8.0_RP*pt(KE-3)))*rfdz(KE-2) &
         - B*(       rcdz(KE-1)-9.0_RP*rcdz(KE-2))
    M(1,KE-2) =                                      A(KE-2)*                 pt(KE-3)  *rfdz(KE-2) &
         - B*                          rcdz(KE-2)
    ! k = KE-1
    M(3,KE-3) = (6.0_RP*A(KE)*pt(KE-1)+A(KE-1)*(6.0_RP*pt(KE-1)+pt(KE-2)))*rfdz(KE-1) &
         + B*(6.0_RP*rcdz(KE)-7.0_RP*rcdz(KE-1)) + 1.0_RP
    M(2,KE-3) = -                      A(KE-1)*                 pt(KE-2)  *rfdz(KE-1) &
         + B*                 8.0_RP*rcdz(KE-1)
    M(1,KE-3) =                        A(KE-1)*                 pt(KE-2)  *rfdz(KE-1) &
         - B*                        rcdz(KE-1)
#endif

    return
  end subroutine make_matrix

  subroutine mul_matrix(V, M, C)
    implicit none
    real(RP), intent(out) :: V(KMAX-1)
    real(RP), intent(in)  :: M(5, KMAX-1)
    real(RP), intent(in)  :: C(KMAX-1)

    integer :: k

    V(1) = M(5,1)*C(3) + M(4,1)*C(2) + M(3,1)*C(1)
    V(2) = M(5,2)*C(4) + M(4,2)*C(3) + M(3,2)*C(2) + M(2,2)*C(1)
    do k = 3, KMAX-3
       V(k) = M(5,k)*C(k+2) + M(4,k)*C(k+1) + M(3,k)*C(k) + M(2,k)*C(k-1) + M(1,k)*C(k-2)
    enddo
    ! k = KE-2
    V(KMAX-2) = M(4,KMAX-2)*C(KMAX-1) + M(3,KMAX-2)*C(KMAX-2) + M(2,KMAX-2)*C(KMAX-3) + M(1,KMAX-2)*C(KMAX-4)
    ! k = KE-1
    V(KMAX-1) = M(3,KMAX-1)*C(KMAX-1) + M(2,KMAX-1)*C(KMAX-2) + M(1,KMAX-1)*C(KMAX-3)

    return
  end subroutine mul_matrix

#else

  subroutine solve_lapack( &
       C,       & ! (inout)
       PT, A, B ) ! (in)
    use mod_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    implicit none

    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: PT(KA)
    real(RP), intent(in)    :: A(KA)
    real(RP), intent(in)    :: B

    real(RP) :: M(NB*3+1,KMAX-1)
    integer  :: IPIV(KMAX-1)
    integer  :: INFO

    integer :: k

#ifdef DEBUG
    real(RP) :: M2(KMAX-1,KMAX-1)
    real(RP) :: C2(KMAX-1)
    real(RP) :: sum
#endif


    M(NB+1,1) = 0.0_RP
    M(NB+2,1) = 0.0_RP
    M(NB+1,2) = 0.0_RP
    M(NB*2  ,KMAX-2) = 0.0_RP
    M(NB*2+1,KMAX-1) = 0.0_RP
    M(NB*2+1,KMAX-1) = 0.0_RP


          ! band matrix
#ifdef __BOUND
          do k = KS+2, KE-3
#else
          do k = KS+3, KE-4
#endif
             M(NB+1,k-KS+1) =   A(k-1) * PT(k-1) * RFDZ(k-2) &
                              + B * RCDZ(k-1)
             M(NB+2,k-KS+1) = - ( A(k  ) * ( 8.0_RP*PT(k) + PT(k-1) ) &
                                + A(k-1) *   PT(k-1)                  ) * RFDZ(k-1) &
                              - B * ( 9.0_RP*RCDZ(k) - RCDZ(k-1) )
             M(NB+3,k-KS+1) =   ( A(k+1) * ( PT(k+1) + 8.0_RP*PT(k) ) &
                                + A(k  ) * ( 8.0_RP*PT(k) + PT(k-1) ) ) * RFDZ(k  ) &
                              + 9.0_RP * B * ( RCDZ(k+1) - RCDZ(k) ) &
                              + 1.0_RP
             M(NB+4,k-KS+1) = - ( A(k+2) *   PT(k+1) &
                                + A(k+1) * ( PT(k+1) + 8.0_RP*PT(k) ) ) * RFDZ(k+1) &
                              - B * ( RCDZ(k+2) - 9.0_RP * RCDZ(k+1) )
             M(NB+5,k-KS+1) =   A(k+2) * PT(k+1) * RFDZ(k+2) &
                              - B * RCDZ(k+2)
          end do
#ifdef __BOUND
          M(NB+3,1) =   ( A(KS+1) * ( PT(KS+1) + 8.0_RP*PT(KS) ) &
                        + A(KS  ) *   8.0_RP*PT(KS)              ) * RFDZ(KS  ) &
                      + B * ( 9.0_RP*RCDZ(KS+1) - 8.0_RP*RCDZ(KS) ) &
                      + 1.0_RP
          M(NB+4,1) = - ( A(KS+2) *   PT(KS+1) &
                        + A(KS+1) * ( PT(KS+1) + 8.0_RP*PT(KS) ) ) * RFDZ(KS+1) &
                      - B * ( RCDZ(KS+2) - 9.0_RP * RCDZ(KS+1) )
          M(NB+5,1) =   A(KS+2) * PT(KS+1) * RFDZ(KS+2) &
                      - B * RCDZ(KS+2)
          M(NB+2,2) = - ( A(KS+1) * ( 8.0_RP*PT(KS+1) + PT(KS) ) &
                        + A(KS  ) *   PT(KS)                     ) * RFDZ(KS) &
                      - B * ( 9.0_RP*RCDZ(KS+1) - RCDZ(KS) )
          M(NB+3,2) =   ( A(KS+2) * ( PT(KS+2) + 8.0_RP*PT(KS+1) ) &
                        + A(KS+1) * ( 8.0_RP*PT(KS+1) + PT(KS  ) ) ) * RFDZ(KS+1) &
                      + 9.0_RP * B * ( RCDZ(KS+2) - RCDZ(KS+1) ) &
                      + 1.0_RP
          M(NB+4,2) = - ( A(KS+3) *   PT(KS+2) &
                        + A(KS+2) * ( PT(KS+2) + 8.0_RP*PT(KS+1) ) ) * RFDZ(KS+2) &
                      - B * ( RCDZ(KS+3) - 9.0_RP * RCDZ(KS+3) )
          M(NB+5,2) =   A(KS+3) * PT(KS+2) * RFDZ(KS+3) &
                      - B * RCDZ(KS+3)
          M(NB+1,KE-KS-1) =   A(KE-3) * PT(KE-3) * RFDZ(KE-4) &
                            + B * RCDZ(KE-3)
          M(NB+2,KE-KS-1) = - ( A(KE-2) * ( 8.0_RP*PT(KE-2) + PT(KE-3) ) &
                              + A(KE-3) *   PT(KE-3)                   ) * RFDZ(KE-3) &
                            - B * ( 9.0_RP*RCDZ(KE-2) - RCDZ(KE-3) )
          M(NB+3,KE-KS-1) =   ( A(KE-1) * ( PT(KE-1) + 8.0_RP*PT(KE-2) ) &
                              + A(KE-2) * ( 8.0_RP*PT(KE-2) + PT(KE-3) ) ) * RFDZ(KE-2) &
                            + 9.0_RP * B * ( RCDZ(KE-1) - RCDZ(KE-2) ) &
                            + 1.0_RP
          M(NB+4,KE-KS-1) = - ( A(KE  ) *   PT(KE-1) &
                              + A(KE-1) * ( PT(KE-1) + 8.0_RP*PT(KE-2) ) ) * RFDZ(KE-1) &
                            - B * ( RCDZ(KE) - 9.0_RP * RCDZ(KE-1) )
          M(NB+1,KE-KS) =   A(KE-2) * PT(KE-2) * RFDZ(KE-3) &
                          + B * RCDZ(KE-2)
          M(NB+2,KE-KS) = - ( A(KE-1) * ( 8.0_RP*PT(KE-1) + PT(KE-2) ) &
                            + A(KE-2) *   PT(KE-2)                   ) * RFDZ(KE-2) &
                          - B * ( 9.0_RP*RCDZ(KE-1) - RCDZ(KE-2) )
          M(NB+3,KE-KS) =   ( A(KE  ) * (              8.0_RP*PT(KE-1) ) &
                            + A(KE-1) * ( 8.0_RP*PT(KE-1) + PT(KE-2) ) ) * RFDZ(KE-1) &
                          + B * ( 8.0_RP*RCDZ(KE) - 9.0_RP*RCDZ(KE-1) ) &
                          + 1.0_RP
#else
          M(NB+3,1) =   ( A(KS+1) * ( PT(KS+1) + 6.0_RP*PT(KS) ) &
                        + 6.0_RP * A(KS) * PT(KS)              ) * RFDZ(KS) &
                      + B * ( 7.0_RP*RCDZ(KS+1) - 6.0_RP*RCDZ(KS) ) &
                      + 1.0_RP
          M(NB+4,1) = - ( A(KS+2) *   PT(KS+1) &
                        + A(KS+1) * ( PT(KS+1) + 6.0_RP*PT(KS) ) ) * RFDZ(KS+1) &
                      - B * ( RCDZ(KS+2) - 7.0_RP*RCDZ(KS+1) )
          M(NB+5,1) =   A(KS+2) * PT(KS+1) * RFDZ(KS+2) &
                      - B * RCDZ(KS+2)
          M(NB+2,2) = - A(KS+1) * PT(KS+1) * RFDZ(KS) &
                      - 8.0_RP * B * RCDZ(KS+1)
          M(NB+3,2) =   ( A(KS+2) * ( PT(KS+2) + 8.0_RP*PT(KS+1) ) &
                        + 8.0_RP*A(KS+1) * PT(KS+1)              ) * RFDZ(KS+1) &
                      + B * ( 9.0_RP*RCDZ(KS+2) - 8.0_RP*RCDZ(KS+1) ) &
                      + 1.0_RP
          M(NB+4,2) = - ( A(KS+3) *   PT(KS+2) &
                        + A(KS+2) * ( PT(KS+2) + 8.0_RP*PT(KS+1) ) ) * RFDZ(KS+2) &
                      - B * ( RCDZ(KS+3) - 9.0_RP * RCDZ(KS+2) )
          M(NB+5,2) =   A(KS+3) * PT(KS+2) * RFDZ(KS+3) &
                      - B * RCDZ(KS+3)
          M(NB+1,3) =   A(KS+1) * PT(KS+1) *RFDZ(KS) &
                      + B * RCDZ(KS+1)
          M(NB+2,3) = - ( A(KS+2) * ( 8.0_RP*PT(KS+2) + PT(KS+1) ) &
                        + A(KS+1) *   PT(KS+1)                     ) * RFDZ(KS+1) &
                      - B * ( 9.0_RP*RCDZ(KS+2) - RCDZ(KS+1) )
          M(NB+3,3) =   ( A(KS+3) * ( PT(KS+3) + 8.0_RP*PT(KS+2) ) &
                        + A(KS+2) * ( 8.0_RP*PT(KS+2) + PT(KS+1) ) ) * RFDZ(KS+2) &
                      + 9.0_RP * B * ( RCDZ(KS+3) - RCDZ(KS+2) ) &
                      + 1.0_RP
          M(NB+4,3) = - ( A(KS+4) *   PT(KS+3) &
                        + A(KS+3) * ( PT(KS+3) + 8.0_RP*PT(KS+2) ) ) * RFDZ(KS+3) &
                      - B * ( RCDZ(KS+4) - 9.0_RP * RCDZ(KS+3) )
          M(NB+5,3) =   A(KS+4) * PT(KS+3) * RFDZ(KS+4) &
                      - B * RCDZ(KS+4)
          M(NB+1,KE-KS-2) =   A(KE-4) * PT(KE-4) * RFDZ(KE-5) &
                            + B * RCDZ(KE-4)
          M(NB+2,KE-KS-2) = - ( A(KE-3) * ( 8.0_RP*PT(KE-3) + PT(KE-4) ) &
                              + A(KE-4) *   PT(KE-4)                     ) * RFDZ(KE-4) &
                            - B * ( 9.0_RP*RCDZ(KE-3) - RCDZ(KE-4) )
          M(NB+3,KE-KS-2) =   ( A(KE-2) * ( PT(KE-2) + 8.0_RP*PT(KE-3) ) &
                              + A(KE-3) * ( 8.0_RP*PT(KE-3) + PT(KE-4) ) ) * RFDZ(KE-3) &
                            + 9.0_RP * B * ( RCDZ(KE-2) - RCDZ(KE-3) ) &
                            + 1.0_RP
          M(NB+4,KE-KS-2) = - ( A(KE-1) *   PT(KE-2) &
                              + A(KE-2) * ( PT(KE-2) + 8.0_RP*PT(KE-3) ) ) * RFDZ(KE-2) &
                            - B * ( RCDZ(KE-1) - 9.0_RP * RCDZ(KE-2) )
          M(NB+5,KE-KS-2) =   A(KE-1) * PT(KE-2) * RFDZ(KE-1) &
                            - B * RCDZ(KE-1)
          M(NB+1,KE-KS-1) =   A(KE-3) * PT(KE-3) * RFDZ(KE-4) &
                            + B * RCDZ(KE-3)
          M(NB+2,KE-KS-1) = - ( A(KE-2) * ( 8.0_RP*PT(KE-2) + PT(KE-3) ) &
                              + A(KE-3) *   PT(KE-3)                     ) * RFDZ(KE-3) &
                            - B * ( 9.0_RP*RCDZ(KE-2) - RCDZ(KE-3) )
          M(NB+3,KE-KS-1) =   ( 8.0_RP * A(KE-1) * PT(KE-2) &
                              + A(KE-2) * ( 8.0_RP*PT(KE-2) + PT(KE-3) ) ) * RFDZ(KE-2) &
                            + B * ( 8.0_RP*RCDZ(KE-1) - 9.0_RP*RCDZ(KE-2) ) &
                            + 1.0_RP
          M(NB+4,KE-KS-1) = - A(KE-1) * PT(KE-2) * RFDZ(KE-1) &
                            + 8.0_RP * B * RCDZ(KE-1)
          M(NB+1,KE-KS  ) =   A(KE-2) * PT(KE-2) * RFDZ(KE-3) &
                            + B * RCDZ(KE-2)
          M(NB+2,KE-KS  ) = - ( A(KE-1) * ( 6.0_RP*PT(KE-1) + PT(KE-2) ) &
                              + A(KE-2) *   PT(KE-2)                     ) * RFDZ(KE-2) &
                            - B * ( 7.0_RP*RCDZ(KE-1) - RCDZ(KE-2) )
          M(NB+3,KE-KS  ) =   ( 6.0_RP * A(KE) * PT(KE-1) &
                              + A(KE-1) * ( 6.0_RP*PT(KE-1) + PT(KE-2) ) ) * RFDZ(KE-1) &
                            + B * ( 6.0_RP*RCDZ(KE) - 7.0_RP*RCDZ(KE-1) ) &
                            + 1.0_RP
#endif


#ifdef DEBUG
          call check_matrix(M, A, B, PT)

          M2(:,:) = 0.0_RP
          do k = 1, KMAX-1
             if (k>2) M2(k-2,k) = M(NB+5,k-2)
             if (k>1) M2(k-1,k) = M(NB+4,k-1)
             M2(k,k) = M(NB+3,k)
             if (k<kmax-1) M2(k+1,k) = M(NB+2,k+1)
             if (k<kmax-2) M2(k+2,k) = M(NB+1,k+2)
             C2(k) = C(k)
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
             call abort
          end if

          do k = 1, KMAX-1
             sum = 0.0_RP
             if (k>2) sum = sum + M2(k-2,k)*c(k-2)
             if (k>1) sum = sum + M2(k-1,k)*c(k-1)
             sum = sum + M2(k,k)*c(k)
             if (k<kmax-1) sum = sum + M2(k+1,k)*c(k+1)
             if (k<kmax-2) sum = sum + M2(k+2,k)*c(k+2)
             if ( abs(sum-c2(k)) > 1E-10_RP ) then
                write(*,*) k+2, i, j, sum, c2(k)
             end if
          end do
#endif
  end subroutine solve_lapack

  subroutine check_matrix( M, A, B, PT )
    use mod_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    implicit none
    real(RP), intent(in) :: M(NB*3+1,KMAX-1)
    real(RP), intent(in) :: A(KA)
    real(RP), intent(in) :: B
    real(RP), intent(in) :: PT(KA)

    integer :: k

#ifdef __BOUND
    k = KS
    if ( M(NB+3,k-2) .ne. (A(k+1)*(pt(k+1)+8.0_RP*pt(k))+A(k)*(8.0_RP*pt(k)))*rfdz(k) + B*(9.0_RP*rcdz(k+1)-8.0_RP*rcdz(k))+1.0_RP ) then
       write(*,*)k,0
    end if
#else
    if ( M(NB+1,3) .ne. A(k+1)*pt(k+1)*rfdz(k) + B*rcdz(k+1) ) then
       write(*,*)k, -2
    end if
    if ( M(NB+2,2) .ne. -A(k+1)*pt(k+1)*rfdz(k) -8.0_RP*B*rcdz(k+1) ) then
       write(*,*)k, -1
    end if
    if ( M(NB+3,1) .ne. (A(k+1)*(pt(k+1)+6.0_RP*pt(k))+6.0_RP*A(k)*pt(k))*rfdz(k) + B*(7.0_RP*rcdz(k+1)-6.0_RP*rcdz(k))+1.0_RP ) then
       write(*,*)k, 0
    end if
    k = KS+1
    if ( M(NB+1,4) .ne. A(k+1)*pt(k+1)*rfdz(k) + B*rcdz(k+1) ) then
       write(*,*)k, -2
    end if
    if ( M(NB+2,3) .ne. -(A(k+1)*(8.0_RP*pt(k+1)+pt(k))+A(k)*pt(k))*rfdz(k) -B*(9.0_RP*rcdz(k+1)-rcdz(k)) ) then
       write(*,*)k, -1
    end if
    if ( M(NB+3,2) .ne. (A(k+1)*(pt(k+1)+8.0_RP*pt(k))+8.0_RP*A(k)*pt(k))*rfdz(k) + B*(9.0_RP*rcdz(k+1)-8.0_RP*rcdz(k))+1.0_RP ) then
       write(*,*)k, 0
    end if
    if ( M(NB+4,1) .ne. -(A(k+1)*pt(k)+A(k)*(pt(k)+6.0_RP*pt(k-1)))*rfdz(k) - B*(rcdz(k+1)-7.0_RP*rcdz(k)) ) then
       write(*,*)k, 1
    end if
#endif

#if __BOUND
    do k = KS, KE-1
#else
    do k = KS+2, KE-3
#endif
       if ( k <= KE-3 ) then
          if( M(NB+1,k+2-2) .ne. A(k+1)*pt(k+1)*rfdz(k) + B *rcdz(k+1) ) then
             write(*,*)k,-2
          end if
       end if
       if ( k <= KE-2 ) then
          if( M(NB+2,k+1-2) .ne. -(A(k+1)*(8.0_RP*pt(k+1)+pt(k))+A(k)*pt(k))*rfdz(k)-B*(9.0_RP*rcdz(k+1)-rcdz(k)) ) then
             write(*,*)k,-1
          end if
       end if
       if ( k >= KS+1 .and. k <= KE-2 ) then
          if ( M(NB+3,k-2) .ne. (A(k+1)*(pt(k+1)+8.0_RP*pt(k))+A(k)*(8.0_RP*pt(k)+pt(k-1)))*rfdz(k) + 9.0_RP*B*(rcdz(k+1)-rcdz(k))+1.0_RP ) then
             write(*,*)k,0
          end if
       end if
       if ( k >= KS+1 ) then
          if ( M(NB+4,k-1-2) .ne. -(A(k+1)*pt(k)+a(k)*(pt(k)+8.0_RP*pt(k-1)))*rfdz(k) -B*(rcdz(k+1)-9.0_RP*rcdz(k)) ) then
             write(*,*)k,1
          end if
       end if
       if ( k >= KS+2 ) then
          if (M(NB+5,k-2-2) .ne. A(k)*pt(k-1)*rfdz(k)-B*rcdz(k) ) then
             write(*,*)k,2
          end if
       end if
    enddo

#ifdef __BOUND
    k = KE-1
    if ( M(NB+3,k-2) .ne. (A(k+1)*(8.0_RP*pt(k))+A(k)*(8.0_RP*pt(k)+pt(k-1)))*rfdz(k) + B*(8.0_RP*rcdz(k+1)-9.0_RP*rcdz(k))+1.0_RP ) then
       write(*,*)k,0
    end if
#else
    k = KE-2
    if ( M(NB+2,k+1-2) .ne. -(A(k+1)*(6.0_RP*pt(k+1)+pt(k))+A(k)*pt(k))*rfdz(k) -B*(7.0_RP*rcdz(k+1)-rcdz(k)) ) then
       write(*,*)k, -1
    end if
    if ( M(NB+3,k-2) .ne. (8.0_RP*A(k+1)*pt(k)+A(k)*(8.0_RP*pt(k)+pt(k-1)))*rfdz(k) + B*(8.0_RP*rcdz(k+1)-9.0_RP*rcdz(k))+1.0_RP ) then
       write(*,*)k, 0
    end if
    if ( M(NB+4,k-1-2) .ne. -(A(k+1)*pt(k)+A(k)*(pt(k)+8.0_RP*pt(k-1)))*rfdz(k) - B*(rcdz(k+1)-9.0_RP*rcdz(k)) ) then
       write(*,*)k,1
    end if
    if ( M(NB+5,k-2-2) .ne. A(k)*pt(k-1)*rfdz(k) - B*rcdz(k) ) then
       write(*,*)k,2
    end if

    k = KE-1
    if ( M(NB+3,k-2) .ne. (6.0_RP*A(k+1)*pt(k)+A(k)*(6.0_RP*pt(k)+pt(k-1)))*rfdz(k) + B*(6.0_RP*rcdz(k+1)-7.0_RP*rcdz(k))+1.0_RP ) then
       write(*,*)k, 0
    end if
    if ( M(NB+4,k-1-2) .ne. -A(k)*pt(k-1)*rfdz(k) + 8.0_RP*B*rcdz(k) ) then
       write(*,*)k, 1
    end if
    if ( M(NB+5,k-2-2) .ne. A(k)*pt(k-1)*rfdz(k) - B*rcdz(k) ) then
       write(*,*)k, 2
    end if
#endif
  end subroutine check_matrix
#endif

  subroutine check_equation( &
       VECT, &
       DENS, MOMZ, RHOT, PRES, &
       Sr, Sw, St, &
#ifdef DRY
       kappa, &
#else
       CPtot, CVtot, &
#endif
       dt, i, j )
    use mod_const, only: &
         EPS => CONST_EPS, &
         GRAV => CONST_GRAV
    use mod_process, only: &
         PRC_MPIstop
    use mod_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    implicit none
    real(RP), intent(in) :: VECT(KMAX-1)
    real(RP), intent(in) :: DENS(KA)
    real(RP), intent(in) :: MOMZ(KA)
    real(RP), intent(in) :: RHOT(KA)
    real(RP), intent(in) :: PRES(KA)
    real(RP), intent(in) :: Sr(KA)
    real(RP), intent(in) :: Sw(KA)
    real(RP), intent(in) :: St(KA)
#ifdef DRY
    real(RP), intent(in) :: kappa
#else
    real(RP), intent(in) :: CPtot(KA)
    real(RP), intent(in) :: CVtot(KA)
#endif
    real(RP), intent(in) :: dt
    integer , intent(in) :: i
    integer , intent(in) :: j

    real(RP), parameter :: small = 1e-5_RP

    real(RP) :: MOMZ_N(KA)
    real(RP) :: DENS_N(KA)
    real(RP) :: RHOT_N(KA)
    real(RP) :: PRES_N(KA)

    real(RP) :: mflx(KA)
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
    MOMZ_N(KE-2:) = 0.0_RP

    ! z momentum flux
    do k = KS+1, KE-2
       mflx(k) = ( - MOMZ_N(k+1) + 8.0_RP * MOMZ_N(k) - MOMZ_N(k-1) ) / 6.0_RP
    end do
    mflx(KS-1) = 0.0_RP
#if __BOUND
    mflx(KS  ) = ( - MOMZ_N(KS+1) + 8.0_RP * MOMZ_N(KS) ) / 6.0_RP
    mflx(KE-1) = ( 8.0_RP * MOMZ_N(KE-1) - MOMZ_N(KE-2) ) / 6.0_RP
#else
    mflx(KS  ) = MOMZ_N(KS)
    mflx(KE-1) = MOMZ_N(KE-1)
#endif
    mflx(KE  ) = 0.0_RP

    ! density
    do k = KS+1, KE-1
       DENS_N(k) = DENS(k) &
            + dt * ( - ( mflx(k) - mflx(k-1) ) * RCDZ(k) + Sr(k) )
    end do
    DENS_N(KS) = DENS(KS) &
         + dt * ( - mflx(KS) * RCDZ(KS) + Sr(KS) )
    DENS_N(KE) = DENS(KE) &
         + dt * ( mflx(KE-1) * RCDZ(KE) + Sr(KE) )

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
            + dt * ( - ( mflx(k)*PT(k) - mflx(k-1)*PT(k-1) ) * RCDZ(k) &
                     + St(k) )
    end do
    RHOT_N(KS) = RHOT(KS) &
         + dt * ( - mflx(KS)*PT(KS) * RCDZ(KS) + St(KS) )
    RHOT_N(KE) = RHOT(KE) &
         + dt * ( mflx(KE-1)*PT(KE-1) * RCDZ(KE) + St(KE) )


    do k = KS, KE
#ifndef DRY
       kappa = CPtot(k) / CVtot(k)
#endif
       PRES_N(k) = PRES(k) * ( 1.0_RP + kappa * ( RHOT_N(k) - RHOT(k) ) / RHOT(k) )
    end do

    do k = KS, KE
       lhs = ( DENS_N(k) - DENS(k) ) / dt
       rhs = - ( mflx(k) - mflx(k-1) ) * RCDZ(k) + Sr(k)
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
       rhs = - ( PRES_N(k+1) - PRES_N(k) ) * RFDZ(k) &
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
          write(*,*) (PRES_N(k+1)-PRES_N(k))*RFDZ(k), GRAV*(DENS_N(k+1)+DENS_N(k))*0.5_RP, Sw(k)
          lhs = MOMZ(k) - dt*RFDZ(k)*( PRES(k+1)*(1.0_RP+kappa*dt*St(k+1)/RHOT(k+1)) &
                                      -PRES(k  )*(1.0_RP+kappa*dt*St(k  )/RHOT(k  ))) &
                        - dt*GRAV*0.5_RP*( DENS(k+1)+DENS(k) + dt*(Sr(k+1)+Sr(k)) ) &
                        + dt*Sw(k)
          a1 = kappa * dt**2 * PRES(k+1) * RCDZ(k+1) / ( 6.0_RP * RHOT(k+1) )
          a0 = kappa * dt**2 * PRES(k  ) * RCDZ(k  ) / ( 6.0_RP * RHOT(k  ) )
          b = GRAV * dt**2 / 12.0_RP
          rhs = ( a1*PT(k+1)*RFDZ(k) + B*RFDZ(k+1) ) * MOMZ_N(k+2) &
              - ( (A1*(8.0_RP*PT(k+1)+PT(k))+A0*PT(k))*RFDZ(k) + B*(9.0_RP*RCDZ(k+1)-RCDZ(k)) ) * MOMZ_N(k+1) &
              + ( (A1*(PT(k+1)+8.0_RP*PT(k))+A0*(8.0_RP*PT(k)+PT(k-1)))*RFDZ(k) + 9.0_RP*B*(RCDZ(k+1)-RCDZ(k)) &
                 + 1.0_RP ) * MOMZ_N(k) &
              - ( (A1*PT(k)+A0*(PT(k)+8.0_RP*PT(k-1)))*RFDZ(k) + B*(RFDZ(k+1)-9.0_RP*RFDZ(k)) ) * MOMZ_N(k-1) &
              + ( A0*PT(k-1)*RFDZ(k) - B*RFDZ(k) ) * MOMZ_N(k-2)
          write(*,*) lhs, rhs
          call PRC_MPIstop
       end if
    end do

    do k = KS, KE
       lhs = ( RHOT_N(k) - RHOT(k) ) / dt
       rhs = - ( mflx(k)*PT(k) - mflx(k-1)*PT(k-1) ) * RCDZ(k) + St(k)
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

end module mod_atmos_dyn_rk
