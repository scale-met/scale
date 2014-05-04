!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          Runge-Kutta for Atmospheric dynamical process
!!          HEVE, no terrain + tracer FCT limiter
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-18 (S.Nishizawa) [new] splited from dynamical core
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_rk_heve
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
  public :: ATMOS_DYN_rk_heve_setup
  public :: ATMOS_DYN_rk_heve

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_rk_heve_setup( ATMOS_TYPE_DYN )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: ATMOS_TYPE_DYN
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HEVE'

    if ( ATMOS_TYPE_DYN .ne. 'HEVE' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_DYN is not HEVE. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_DYN_rk_heve_setup

  !-----------------------------------------------------------------------------
  !> Runge-Kutta loop
  subroutine ATMOS_DYN_rk_heve( &
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
       GRAV   => CONST_GRAV,   &
       P00    => CONST_PRE00,  &
#ifdef DRY
       CPovCV => CONST_CPovCV
#endif
       Rdry   => CONST_Rdry
    use scale_comm, only: &
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

    real(RP), intent(out) :: mflx_hi(KA,IA,JA,3) ! mass flux

    real(RP), intent(in),target :: DENS0(KA,IA,JA) ! prognostic variables at previous dynamical time step
    real(RP), intent(in),target :: MOMZ0(KA,IA,JA) !
    real(RP), intent(in),target :: MOMX0(KA,IA,JA) !
    real(RP), intent(in),target :: MOMY0(KA,IA,JA) !
    real(RP), intent(in),target :: RHOT0(KA,IA,JA) !

    real(RP), intent(in)  :: DENS(KA,IA,JA)      ! prognostic variables at previous RK step
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)      !
    real(RP), intent(in)  :: MOMX(KA,IA,JA)      !
    real(RP), intent(in)  :: MOMY(KA,IA,JA)      !
    real(RP), intent(in)  :: RHOT(KA,IA,JA)      !

    real(RP), intent(in)  :: DENS_t(KA,IA,JA)    ! tendency
    real(RP), intent(in)  :: MOMZ_t(KA,IA,JA)    !
    real(RP), intent(in)  :: MOMX_t(KA,IA,JA)    !
    real(RP), intent(in)  :: MOMY_t(KA,IA,JA)    !
    real(RP), intent(in)  :: RHOT_t(KA,IA,JA)    !

    real(RP), intent(in)  :: Rtot    (KA,IA,JA)  ! total R
    real(RP), intent(in)  :: CVtot   (KA,IA,JA)  ! total CV
    real(RP), intent(in)  :: CORIOLI (1, IA,JA)
    real(RP), intent(in)  :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in)  :: divdmp_coef

    logical,  intent(in)  :: FLAG_FCT_RHO
    logical,  intent(in)  :: FLAG_FCT_MOMENTUM
    logical,  intent(in)  :: FLAG_FCT_T

    real(RP), intent(in)  :: CDZ (KA)
    real(RP), intent(in)  :: FDZ (KA-1)
    real(RP), intent(in)  :: FDX (IA-1)
    real(RP), intent(in)  :: FDY (JA-1)
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: RCDX(IA)
    real(RP), intent(in)  :: RCDY(JA)
    real(RP), intent(in)  :: RFDZ(KA-1)
    real(RP), intent(in)  :: RFDX(IA-1)
    real(RP), intent(in)  :: RFDY(JA-1)

    real(RP), intent(in)  :: PHI     (KA,IA,JA)   !< geopotential
    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: REF_pres(KA,IA,JA)   !< reference pressure
    real(RP), intent(in)  :: REF_dens(KA,IA,JA)   !< reference density

    real(RP), intent(in)  :: dtrk

    ! diagnostic variables
    real(RP) :: PRES (KA,IA,JA) ! pressure [Pa]
    real(RP) :: VELZ (KA,IA,JA) ! velocity w [m/s]
    real(RP) :: VELX (KA,IA,JA) ! velocity u [m/s]
    real(RP) :: VELY (KA,IA,JA) ! velocity v [m/s]
    real(RP) :: POTT (KA,IA,JA) ! potential temperature [K]
    real(RP) :: DDIV (KA,IA,JA) ! divergence
    real(RP) :: DPRES(KA,IA,JA) ! pressure - reference pressure

    real(RP) :: qflx_J13(KA,IA,JA)
    real(RP) :: qflx_J23(KA,IA,JA)
    real(RP) :: pgf     (KA,IA,JA)  ! pressure gradient force
    real(RP) :: buoy    (KA,IA,JA)  ! buoyancy force
    real(RP) :: cor     (KA,IA,JA)  ! Coriolis force

    ! flux
    real(RP) :: qflx_hi  (KA,IA,JA,3)
    real(RP) :: qflx_lo  (KA,IA,JA,3)
    real(RP) :: qflx_anti(KA,IA,JA,3)

    real(RP) :: sw

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: k, i, j
    !---------------------------------------------------------------------------

#ifdef DEBUG
    PRES(:,:,:) = UNDEF
    VELZ(:,:,:) = UNDEF
    VELX(:,:,:) = UNDEF
    VELY(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF

    DPRES(:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    qflx_hi(:,:,:,:) = UNDEF

#ifndef NO_FCT_DYN
    qflx_lo  (:,:,:,:) = UNDEF
    qflx_anti(:,:,:,:) = UNDEF
#endif
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

       ! pressure, pott. temp.

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, RHOT(k,i,j) )
             call CHECK( __LINE__, Rtot(k,i,j) )
             call CHECK( __LINE__, CVtot(k,i,j) )
#endif
#ifdef DRY
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rdry / P00 )**CPovCV &
#else
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rtot(k,i,j) / P00 )**((CVtot(k,i,j)+Rtot(k,i,j))/CVtot(k,i,j))
#endif
          enddo

          PRES(KS-1,i,j) = PRES(KS+1,i,j) - DENS(KS,i,j) * ( PHI(KS-1,i,j) - PHI(KS+1,i,j) )
          PRES(KE+1,i,j) = PRES(KE-1,i,j) - DENS(KE,i,j) * ( PHI(KE+1,i,j) - PHI(KE-1,i,j) )

          do k = KS-1, KE+1
             DPRES(k,i,j) = PRES(k,i,j) - REF_pres(k,i,j)
          enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

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

#ifndef NO_FCT_DYN
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
       call COMM_vars8( VELZ(:,:,:), 1 )
       call COMM_vars8( VELX(:,:,:), 2 )
       call COMM_vars8( VELY(:,:,:), 3 )
       call COMM_wait ( VELZ(:,:,:), 1 )
       call COMM_wait ( VELX(:,:,:), 2 )
       call COMM_wait ( VELY(:,:,:), 3 )
    endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#endif

       !########################################################################
       ! continuity equation (total rho)
       !########################################################################

       !-----< high order flux >-----

       ! at (x, y, w)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,ZDIR) )
#endif
          mflx_hi(k,i,j,ZDIR) = J33G * MOMZ(k,i,j) &
                              + J13G(k,i,j,I_XYW) * 0.25_RP * ( MOMX(k+1,i,j)+MOMX(k+1,i-1,j) &
                                                              + MOMX(k  ,i,j)+MOMX(k  ,i-1,j) ) & ! [{u,y,z->x,y,w}]
                              + J23G(k,i,j,I_XYW) * 0.25_RP * ( MOMY(k+1,i,j)+MOMY(k+1,i,j-1) &
                                                              + MOMY(k  ,i,j)+MOMY(k  ,i,j-1) ) & ! [{x,v,z->x,y,w}]
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
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (u, y, z)

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
          mflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) * ( MOMX(k,i,j) + num_diff(k,i,j,I_DENS,XDIR) )
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
          call CHECK( __LINE__, MOMY(k,i,j+1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) * ( MOMY(k,i,j) + num_diff(k,i,j,I_DENS,YDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !-----< update density >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS0(k,i,j) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, mflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, mflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DENS_t(k,i,j) )
#endif
          DENS_RK(k,i,j) = DENS0(k,i,j) &
                         + dtrk * ( - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                                      + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                                      + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) & ! divergence
                                  ) / GSQRT(k,i,j,I_XYZ) &
                         + dtrk * DENS_t(k,i,j) ! physics tendency
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifndef NO_FCT_DYN
       if ( FLAG_FCT_RHO ) then

          ! at (x, y, interface)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k,i,j) )
             call CHECK( __LINE__, DENS(k  ,i,j) )
             call CHECK( __LINE__, DENS(k+1,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * (     VELZ(k,i,j)  * ( DENS(k+1,i,j)+DENS(k,i,j) ) &
                                            - abs(VELZ(k,i,j)) * ( DENS(k+1,i,j)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (u, y, layer)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELX(k,i,j) )
             call CHECK( __LINE__, DENS(k,i  ,j) )
             call CHECK( __LINE__, DENS(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.5_RP * (     VELX(k,i,j)  * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                                            - abs(VELX(k,i,j)) * ( DENS(k,i+1,j)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (x, v, layer)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELY(k,i,j) )
             call CHECK( __LINE__, DENS(k,i,j  ) )
             call CHECK( __LINE__, DENS(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.5_RP * (     VELY(k,i,j)  * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
                                            - abs(VELY(k,i,j)) * ( DENS(k,i,j+1)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       endif
    enddo
    enddo

    if ( FLAG_FCT_RHO ) then
       call ATMOS_DYN_fct( qflx_anti,               & ! (out)
                           DENS0, mflx_hi, qflx_lo, & ! (in)
                           RCDZ, RCDX, RCDY, dtrk   ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update rho
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, DENS_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             DENS_RK(k,i,j) = DENS_RK(k,i,j) &
                            + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                         + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                         + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       enddo
       enddo

       !--- update mflx_hi
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k,i,j,ZDIR) )
#endif
          mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + qflx_anti(k,i,j,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS  , JE
       do i = IS-1, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, qflx_anti(k,i,j,XDIR) )
#endif
          mflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) + qflx_anti(k,i,j,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE
       do i = IS  , IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, qflx_anti(k,i,j,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) + qflx_anti(k,i,j,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifdef DEBUG
       qflx_lo  (:,:,:,:) = UNDEF
       qflx_anti(:,:,:,:) = UNDEF
#endif

    endif ! FLAG_FCT_RHO?

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#endif

       !########################################################################
       ! momentum equation (z)
       !########################################################################

       !-----< high order flux >-----

       ! at (x, y, z)
       ! note than z-index is added by -1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+2, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k-2,i,j) )
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMZ(k+1,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMZ,ZDIR) )
#endif
          qflx_hi(k-1,i,j,ZDIR) = J33G &
                                * 0.5_RP * ( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) / DENS(k,i,j) & ! [x,y,w->x,y,z]
                                * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                                  + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) ) &         ! {x,y,w->x,y,z}
                                + GSQRT(k,i,j,I_XYZ) * num_diff(k,i,j,I_MOMZ,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k,sw) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, MOMZ(KS+1,i,j) )
          call CHECK( __LINE__, num_diff(KS+1,i,j,I_MOMZ,ZDIR) )
          call CHECK( __LINE__, MOMZ(KE-2,i,j) )
          call CHECK( __LINE__, MOMZ(KE-1,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMZ,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP ! k = KS

          qflx_hi(KS  ,i,j,ZDIR) = J33G &
                                 * 0.5_RP * ( MOMZ(KS+1,i,j)+MOMZ(KS  ,i,j) ) / DENS(KS+1,i,j) & ! [x,y,w->x,y,z]
                                 * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS  ,i,j) ) &
                                   + FACT_F * ( MOMZ(KS+2,i,j)                ) ) &              ! {x,y,w->x,y,z}
                                 + GSQRT(KS+1,i,j,I_XYZ) * num_diff(KS+1,i,j,I_MOMZ,ZDIR) ! k = KS+1

          ! if w>0; min(f,w*dz/dt)
          ! else  ; max(f,w*dz/dt) = -min(-f,-w*dz/dt)
          sw = sign( 1.0_RP, MOMZ(KS,i,j) )
          qflx_hi(KS  ,i,j,ZDIR) = sw * min( sw*qflx_hi(KS,i,j,ZDIR), sw*MOMZ(KS,i,j)*GSQRT(KS,i,j,I_XYZ)*FDZ(KS)/dtrk )

          qflx_hi(KE-2,i,j,ZDIR) = J33G &
                                 * 0.5_RP * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) / DENS(KE-1,i,j) & ! [x,y,w->x,y,z]
                                 * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                                   + FACT_F * (                MOMZ(KE-3,i,j) ) ) &              ! {x,y,w->x,y,z}
                                 + GSQRT(KE-1,i,j,I_XYZ) * num_diff(KE-1,i,j,I_MOMZ,ZDIR) ! k = KE-1

          qflx_hi(KE-1,i,j,ZDIR) = 0.0_RP ! k = KE
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+2, KE-2
          qflx_J13(k-1,i,j) = J13G(k,i,j,I_XYZ) &
                            * 0.5_RP * ( MOMX(k,i,j)+MOMX(k,i-1,j) ) / DENS(k,i,j) & ! [u,y,z->x,y,z]
                            * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                              + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) )           ! {x,y,w->x,y,z}
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J13(KS-1,i,j) = 0.0_RP

          qflx_J13(KS  ,i,j) = J13G(KS+1,i,j,I_XYZ) &
                             * 0.5_RP * ( MOMX(KS+1,i,j)+MOMX(KS+1,i-1,j) ) / DENS(KS+1,i,j) & ! [u,y,z->x,y,z]
                             * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS  ,i,j) ) &
                               + FACT_F * ( MOMZ(KS+2,i,j)                ) )                  ! {x,y,w->x,y,z}

          qflx_J13(KE-2,i,j) = J13G(KE-1,i,j,I_XYZ) &
                             * 0.5_RP * ( MOMX(KE-1,i,j)+MOMX(KE-1,i-1,j) ) / DENS(KE-1,i,j) & ! [u,y,z->x,y,z]
                             * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                               + FACT_F * (                MOMZ(KE-3,i,j) ) )                  ! {x,y,w->x,y,z}

          qflx_J13(KE-1,i,j) = 0.0_RP
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+2, KE-2
          qflx_J23(k-1,i,j) = J23G(k,i,j,I_XYZ) &
                            * 0.5_RP * ( MOMY(k,i,j)+MOMY(k,i,j-1) ) / DENS(k,i,j) & ! [x,v,z->x,y,z]
                            * ( FACT_N * ( MOMZ(k  ,i,j)+MOMZ(k-1,i,j) ) &
                              + FACT_F * ( MOMZ(k+1,i,j)+MOMZ(k-2,i,j) ) )           ! {x,y,w->x,y,z}
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J23(KS-1,i,j) = 0.0_RP

          qflx_J23(KS  ,i,j) = J23G(KS+1,i,j,I_XYZ) &
                             * 0.5_RP * ( MOMY(KS+1,i,j)+MOMY(KS+1,i,j-1) ) / DENS(KS+1,i,j) & ! [x,v,z->x,y,z]
                             * ( FACT_N * ( MOMZ(KS+1,i,j)+MOMZ(KS  ,i,j) ) &
                               + FACT_F * ( MOMZ(KS+2,i,j)                ) )            ! {x,y,w->x,y,z}

          qflx_J23(KE-2,i,j) = J23G(KE-1,i,j,I_XYZ) &
                             * 0.5_RP * ( MOMY(KE-1,i,j)+MOMY(KE-1,i,j-1) ) / DENS(KE-1,i,j) & ! [x,v,z->x,y,z]
                             * ( FACT_N * ( MOMZ(KE-1,i,j)+MOMZ(KE-2,i,j) ) &
                               + FACT_F * (                MOMZ(KE-3,i,j) ) )                  ! {x,y,w->x,y,z}

          qflx_J23(KE-1,i,j) = 0.0_RP
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
                              * 0.5_RP * ( VELX(k+1,i,j)+VELX(k,i,j) ) &       ! [u,y,z->u,y,w]
                              * ( FACT_N * ( MOMZ(k,i+1,j)+MOMZ(k,i  ,j) ) &
                                + FACT_F * ( MOMZ(k,i+2,j)+MOMZ(k,i-1,j) ) ) & ! {x,y,w->u,y,w}
                              + GSQRT(k,i,j,I_UYW) * num_diff(k,i,j,I_MOMZ,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
          qflx_hi(KE,i,j,XDIR) = 0.0_RP
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
                              * 0.5_RP * ( VELY(k+1,i,j)+VELY(k,i,j) ) &       ! [x,v,z->x,v,w]
                              * ( FACT_N * ( MOMZ(k,i,j+1)+MOMZ(k,i,j  ) ) &
                                + FACT_F * ( MOMZ(k,i,j+2)+MOMZ(k,i,j-1) ) ) & ! {x,y,w->x,v,w}
                              + GSQRT(k,i,j,I_XVW) * num_diff(k,i,j,I_MOMZ,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
          qflx_hi(KE,i,j,YDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! pressure gradient force at (x, y, w)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          pgf(k,i,j) = J33G * ( DPRES(k+1,i,j)-DPRES(k,i,j) ) * RFDZ(k) ! [x,y,z]
       enddo
       enddo
       enddo

       ! buoyancy force at (x, y, w)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          buoy(k,i,j) = GRAV * 0.5_RP * ( GSQRT(k+1,i,j,I_XYZ) * ( DENS(k+1,i,j)-REF_dens(k+1,i,j) ) &
                                        + GSQRT(k  ,i,j,I_XYZ) * ( DENS(k  ,i,j)-REF_dens(k  ,i,j) ) ) ! [x,y,z]
       enddo
       enddo
       enddo

       !-----< update momentum (z) -----

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
          call CHECK( __LINE__, DDIV(k  ,i,j) )
          call CHECK( __LINE__, DDIV(k+1,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, MOMZ_t(k,i,j) )
#endif
          MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                         + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
                                      + ( qflx_J13(k,i,j) - qflx_J13(k-1,i,j)             ) * RFDZ(k) &
                                      + ( qflx_J23(k,i,j) - qflx_J23(k-1,i,j)             ) * RFDZ(k) &
                                      + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                      + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) & ! advection
                                    - pgf (k,i,j)                                                       & ! pressure gradient force
                                    - buoy(k,i,j)                                                       & ! buoyancy force
                                  ) / GSQRT(k,i,j,I_XYW) &
                         + dtrk * ( divdmp_coef * dtrk * FDZ(k) * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) & ! divergence damping
                                  + MOMZ_t(k,i,j)                                               ) ! physics tendency
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          MOMZ_RK(KS-1,i,j) = 0.0_RP
          MOMZ_RK(KE  ,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifndef NO_FCT_DYN
       if ( FLAG_FCT_MOMENTUM ) then

          ! monotonic flux
          ! at (x, y, layer)
          ! note than z-index is added by -1

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS+1, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k-1,i,j) )
             call CHECK( __LINE__, MOMZ(k  ,i,j) )
             call CHECK( __LINE__, VELZ(k-1,i,j) )
             call CHECK( __LINE__, VELZ(k  ,i,j) )
#endif
             qflx_lo(k-1,i,j,ZDIR) = 0.25_RP * (    ( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) * ( VELZ(k,i,j)+VELZ(k-1,i,j) ) &
                                               - abs( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) * ( VELZ(k,i,j)-VELZ(k-1,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP ! k = KS
             qflx_lo(KE-1,i,j,ZDIR) = 0.0_RP ! k = KE
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP ! k = KE+1
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (u, y, interface)
          !$omp parallel do private(i,j,k,vel) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMX(k  ,i,j) )
             call CHECK( __LINE__, MOMX(k+1,i,j) )
             call CHECK( __LINE__, VELZ(k,i,j) )
             call CHECK( __LINE__, VELZ(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.25_RP * (    ( MOMX(k+1,i,j)+MOMX(k,i,j) ) * ( VELZ(k,i+1,j)+VELZ(k,i,j) ) &
                                             - abs( MOMX(k+1,i,j)+MOMX(k,i,j) ) * ( VELZ(k,i+1,j)-VELZ(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
             qflx_lo(KE,i,j,XDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (x, v, interface)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i,j  ) )
             call CHECK( __LINE__, MOMY(k+1,i,j) )
             call CHECK( __LINE__, VELZ(k,i,j) )
             call CHECK( __LINE__, VELZ(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.25_RP * (    (MOMY(k+1,i,j)+MOMY(k,i,j)) * ( VELZ(k,i,j+1)+VELZ(k,i,j) ) &
                                             - abs(MOMY(k+1,i,j)+MOMY(k,i,j)) * ( VELZ(k,i,j+1)-VELZ(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE,i,j,YDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       endif
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then

       call ATMOS_DYN_fct( qflx_anti,               & ! (out)
                           MOMZ0, qflx_hi, qflx_lo, & ! (in)
                           RFDZ, RCDX, RCDY, dtrk   ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(z)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMZ_RK(k,i,j) = MOMZ_RK(k,i,j) &
                            + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
                                         + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                         + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo

       enddo
       enddo

    endif ! FLAG_FCT_MOMENTUM

#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#endif

       !########################################################################
       ! momentum equation (x)
       !########################################################################

       !-----< high order flux >-----

       ! at (u, y, w)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, VELZ(k,i  ,j) )
          call CHECK( __LINE__, VELZ(k,i+1,j) )
          call CHECK( __LINE__, MOMX(k-1,i,j) )
          call CHECK( __LINE__, MOMX(k  ,i,j) )
          call CHECK( __LINE__, MOMX(k+1,i,j) )
          call CHECK( __LINE__, MOMX(k+2,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_MOMX,ZDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = J33G &
                              * 0.5_RP * ( VELZ(k,i+1,j)+VELZ(k,i,j) )     &   ! [x,y,w->u,y,w]
                              * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                                + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) ) & ! {u,y,z->u,y,w}
                              + GSQRT(k,i,j,I_UYW) * num_diff(k,i,j,I_MOMX,ZDIR)
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
          call CHECK( __LINE__, MOMX(KS+1,i,j) )
          call CHECK( __LINE__, MOMX(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_MOMX,ZDIR) )
          call CHECK( __LINE__, VELZ(KE-1,i  ,j) )
          call CHECK( __LINE__, VELZ(KE-1,i+1,j) )
          call CHECK( __LINE__, MOMX(KE-1,i,j) )
          call CHECK( __LINE__, MOMX(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_MOMX,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP

          qflx_hi(KS  ,i,j,ZDIR) = J33G &
                                 * 0.25_RP * ( VELZ(KS,i+1,j)+VELZ(KS,i,j) ) & ! [x,y,w->u,y,w]
                                           * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) ) & ! [u,y,z->u,y,w]
                                 + GSQRT(KS  ,i,j,I_UYW) * num_diff(KS  ,i,j,I_MOMX,ZDIR)

          qflx_hi(KE-1,i,j,ZDIR) = J33G &
                                 * 0.25_RP * ( VELZ(KE-1,i+1,j)+VELZ(KE-1,i,j) ) & ! [x,y,w->u,y,w]
                                           * ( MOMX(KE  ,i  ,j)+MOMX(KE-1,i,j) ) & ! [u,y,z->u,y,w]
                                 + GSQRT(KE-1,i,j,I_UYW) * num_diff(KE-1,i,j,I_MOMX,ZDIR)

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
          qflx_J13(k,i,j) = J13G(k,i,j,I_UYW) &
                          * 0.5_RP * ( VELX(k+1,i,j)+VELX(k,i,j) ) &     ! [u,y,z->u,y,w]
                          * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                            + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) ) ! {u,y,z->u,y,w}
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J13(KS-1,i,j) = 0.0_RP

          qflx_J13(KS  ,i,j) = J13G(KS,i,j,I_UYW) &
                             * 0.25_RP * ( VELX(KS+1,i,j)+VELX(KS,i,j) ) & ! [u,y,z->u,y,w]
                                       * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) )   ! [u,y,z->u,y,w]

          qflx_J13(KE-1,i,j) = J13G(KE-1,i,j,I_UYW) &
                             * 0.25_RP * ( VELX(KE,i,j)+VELX(KE-1,i,j) ) & ! [u,y,z->u,y,w]
                                       * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) )   ! [u,y,z->u,y,w]

          qflx_J13(KE  ,i,j) = 0.0_RP
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
          qflx_J23(k,i,j) = J23G(k,i,j,I_UYW) &
                          * 0.125_RP * ( VELY(k+1,i+1,j  )+VELY(k,i+1,j  ) &
                                       + VELY(k+1,i  ,j  )+VELY(k,i  ,j  ) &
                                       + VELY(k+1,i+1,j-1)+VELY(k,i+1,j-1) &
                                       + VELY(k+1,i  ,j-1)+VELY(k,i  ,j-1) ) & ! [x,u,z->u,y,w]
                          * ( FACT_N * ( MOMX(k+1,i,j)+MOMX(k  ,i,j) ) &
                            + FACT_F * ( MOMX(k+2,i,j)+MOMX(k-1,i,j) ) )       ! {u,y,z->u,y,w}
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J23(KS-1,i,j) = 0.0_RP

          qflx_J23(KS  ,i,j) = J23G(KS,i,j,I_UYW) &
                             * 0.0625_RP * ( VELY(KS+1,i+1,j  )+VELY(KS,i+1,j  ) &
                                           + VELY(KS+1,i  ,j  )+VELY(KS,i  ,j  ) &
                                           + VELY(KS+1,i+1,j-1)+VELY(KS,i+1,j-1) &
                                           + VELY(KS+1,i  ,j-1)+VELY(KS,i  ,j-1) ) & ! [u,y,z->u,y,w]
                                         * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) )           ! [u,y,z->u,y,w]

          qflx_J23(KE-1,i,j) = J23G(KE-1,i,j,I_UYW) &
                             * 0.0625_RP * ( VELY(KE,i+1,j  )+VELY(KE-1,i+1,j  ) &
                                           + VELY(KE,i  ,j  )+VELY(KE-1,i  ,j  ) &
                                           + VELY(KE,i+1,j-1)+VELY(KE-1,i+1,j-1) &
                                           + VELY(KE,i  ,j-1)+VELY(KE-1,i  ,j-1) ) & ! [u,y,z->u,y,w]
                                         * ( MOMX(KE,i,j)+MOMX(KE-1,i,j) )           ! [u,y,z->u,y,w]

          qflx_J23(KE  ,i,j) = 0.0_RP
       enddo
       enddo

       ! at (x, y, z)
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
          qflx_hi(k,i-1,j,XDIR) = GSQRT(k,i,j,I_XYZ) &
                                * 0.5_RP * ( VELX(k,i,j)+VELX(k,i-1,j) ) &       ! [u,y,z->x,y,z]
                                * ( FACT_N * ( MOMX(k,i  ,j)+MOMX(k,i-1,j) ) &
                                  + FACT_F * ( MOMX(k,i+1,j)+MOMX(k,i-2,j) ) ) & ! {u,y,z->x,y,z}
                                + GSQRT(k,i,j,I_XYZ) * num_diff(k,i,j,I_MOMX,XDIR)
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
                              * 0.5_RP * ( VELY(k,i+1,j)+VELY(k,i,j) ) &       ! [x,v,z->u,v,z]
                              * ( FACT_N * ( MOMX(k,i,j+1)+MOMX(k,i,j  ) ) &
                                + FACT_F * ( MOMX(k,i,j+2)+MOMX(k,i,j-1) ) ) & ! {u,y,z->u,v,z}
                              + GSQRT(k,i,j,I_UVZ) * num_diff(k,i,j,I_MOMX,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! pressure gradient force at (u, y, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          pgf(k,i,j) = ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) & ! [x,y,z]
                       - GSQRT(k,i  ,j,I_XYZ) * DPRES(k,i  ,j) & ! [x,y,z]
                       ) * RFDX(i) &
                     + ( J13G(k  ,i,j,I_UYW) &
                       * 0.25_RP * ( DPRES(k+1,i+1,j)+DPRES(k,i+1,j) &
                                   + DPRES(k+1,i  ,j)+DPRES(k,i  ,j) ) & ! [x,y,z->u,y,w]
                       - J13G(k-1,i,j,I_UYW) &
                       * 0.25_RP * ( DPRES(k,i+1,j)+DPRES(k-1,i+1,j) &
                                   + DPRES(k,i  ,j)+DPRES(k-1,i  ,j) ) & ! [x,y,z->u,y,w]
                       ) * RCDZ(k)
       enddo
       enddo
       enddo

       ! coriolis force at (u, y, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          cor(k,i,j) = 0.0625_RP * ( CORIOLI(1,i+1,j  )+CORIOLI(1,i,j  ) ) & ! [x,y,z->u,y,z]
                                 * ( DENS   (k,i+1,j  )+DENS   (k,i,j  ) ) & ! [x,y,z->u,y,z]
                                 * ( VELY   (k,i+1,j  )+VELY   (k,i,j  ) &
                                   + VELY   (k,i+1,j-1)+VELY   (k,i,j-1) )   ! [x,v,z->u,y,z]
       enddo
       enddo
       enddo

       !-----< update momentum (x) >-----

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
                                      + ( qflx_J13(k,i,j) - qflx_J13(k-1,i,j)             ) * RCDZ(k) &
                                      + ( qflx_J23(k,i,j) - qflx_J23(k-1,i,j)             ) * RCDZ(k) &
                                      + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                                      + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) & ! advection
                                    - pgf(k,i,j)                                                        & ! pressure gradient force
                                  ) / GSQRT(k,i,j,I_UYZ) &
                         + dtrk * ( + cor(k,i,j)                                                        & ! coriolis force
                                    + divdmp_coef * dtrk * FDX(i) * ( DDIV(k,i+1,j)-DDIV(k,i,j) )       & ! divergence damping
                                    + MOMX_t(k,i,j)                                                     ) ! physics tendency
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifndef NO_FCT_DYN
       if ( FLAG_FCT_MOMENTUM ) then

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k,i  ,j) )
             call CHECK( __LINE__, MOMZ(k,i+1,j) )
             call CHECK( __LINE__, VELX(k  ,i,j) )
             call CHECK( __LINE__, VELX(k+1,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.25_RP * (    ( MOMZ(k,i+1,j)+MOMZ(k,i,j) ) * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
                                             - abs( MOMZ(k,i+1,j)+MOMZ(k,i,j) ) * ( VELX(k+1,i,j)-VELX(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! note that x-index is added by -1
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+2
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMX(k,i-1,j) )
             call CHECK( __LINE__, MOMX(k,i  ,j) )
             call CHECK( __LINE__, VELX(k,i-1,j) )
             call CHECK( __LINE__, VELX(k,i  ,j) )
#endif
             qflx_lo(k,i-1,j,XDIR) = 0.25_RP * (    ( MOMX(k,i,j)+MOMX(k,i-1,j) ) * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
                                               - abs( MOMX(k,i,j)+MOMX(k,i-1,j) ) * ( VELX(k,i,j)-VELX(k,i-1,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i+1,j) )
             call CHECK( __LINE__, MOMY(k,i  ,j) )
             call CHECK( __LINE__, VELX(k,i,j  ) )
             call CHECK( __LINE__, VELX(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.25_RP * (    ( MOMY(k,i+1,j)+MOMY(k,i,j) ) * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
                                             - abs( MOMY(k,i+1,j)+MOMY(k,i,j) ) * ( VELX(k,i,j+1)-VELX(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       endif
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then

       call ATMOS_DYN_fct( qflx_anti,               & ! (out)
                           MOMX0, qflx_hi, qflx_lo, & ! (in)
                           RCDZ, RFDX, RCDY, dtrk   ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(x)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMX_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             MOMX_RK(k,i,j) = MOMX_RK(k,i,j) &
                            + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                         + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                                         + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       enddo
       enddo

#ifdef DEBUG
       qflx_lo  (:,:,:,:) = UNDEF
       qflx_anti(:,:,:,:) = UNDEF
#endif

    endif ! FLAG_FCT_MOMENTUM

#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#endif

       !########################################################################
       ! momentum equation (y)
       !########################################################################

       !-----< high order flux >-----

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
          qflx_hi(k,i,j,ZDIR) = J33G &
                              * 0.5_RP * ( VELZ(k,i,j+1)+VELZ(k,i,j) ) &       ! [x,y,w->x,v,w]
                              * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                                + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) ) & ! {x,v,z->x,v,w}
                              + GSQRT(k,i,j,I_XVW) * num_diff(k,i,j,I_MOMY,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

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
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP

          qflx_hi(KS  ,i,j,ZDIR) = J33G &
                                 * 0.25_RP * ( VELZ(KS  ,i,j+1)+VELZ(KS,i,j) ) & ! [x,y,w->x,v,w]
                                           * ( MOMY(KS+1,i,j  )+MOMY(KS,i,j) ) & ! [x,v,z->x,v,w]
                                 + GSQRT(KS  ,i,j,I_XVW) * num_diff(KS  ,i,j,I_MOMY,ZDIR)

          qflx_hi(KE-1,i,j,ZDIR) = J33G &
                                 * 0.25_RP * ( VELZ(KE-1,i,j+1)+VELZ(KE-1,i,j) ) & ! [x,y,w->x,v,w]
                                           * ( MOMY(KE  ,i,j  )+MOMY(KE-1,i,j) ) & ! [x,v,z->x,v,w]
                                 + GSQRT(KE-1,i,j,I_XVW) * num_diff(KE-1,i,j,I_MOMY,ZDIR)

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
          qflx_J13(k,i,j) = J13G(k,i,j,I_XVW) &
                          * 0.125_RP * ( VELX(k+1,i  ,j+1)+VELX(k,i  ,j+1) &
                                       + VELX(k+1,i-1,j+1)+VELX(k,i-1,j+1) &
                                       + VELX(k+1,i  ,j  )+VELX(k,i  ,j  ) &
                                       + VELX(k+1,i-1,j  )+VELX(k,i-1,j  ) ) & ! [u,y,z->x,v,w]
                          * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                            + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) )       ! {x,v,z->x,v,w}
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J13(KS-1,i,j) = 0.0_RP

          qflx_J13(KS  ,i,j) = J13G(KS,i,j,I_XVW) &
                             * 0.0625_RP * ( VELX(KS+1,i  ,j+1)+VELX(KS,i  ,j+1) &
                                           + VELX(KS+1,i-1,j+1)+VELX(KS,i-1,j+1) &
                                           + VELX(KS+1,i  ,j  )+VELX(KS,i  ,j  ) &
                                           + VELX(KS+1,i-1,j  )+VELX(KS,i-1,j  ) ) & ! [u,y,z->x,v,w]
                                         * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) )           ! [x,v,z->x,v,w]

          qflx_J13(KE-1,i,j) = J13G(KE-1,i,j,I_XVW) &
                             * 0.0625_RP * ( VELX(KE,i  ,j+1)+VELX(KE-1,i  ,j+1) &
                                           + VELX(KE,i-1,j+1)+VELX(KE-1,i-1,j+1) &
                                           + VELX(KE,i  ,j  )+VELX(KE-1,i  ,j  ) &
                                           + VELX(KE,i-1,j  )+VELX(KE-1,i-1,j  ) ) & ! [u,y,z->x,v,w]
                                         * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) )           ! [x,v,z->x,v,w]

          qflx_J13(KE  ,i,j) = 0.0_RP
       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-2
          qflx_J23(k,i,j) = J23G(k,i,j,I_XVW) &
                          * 0.5_RP * ( VELY(k+1,i,j)+VELY(k,i,j) ) &     ! [x,v,z->x,v,w]
                          * ( FACT_N * ( MOMY(k+1,i,j)+MOMY(k  ,i,j) ) &
                            + FACT_F * ( MOMY(k+2,i,j)+MOMY(k-1,i,j) ) ) ! {x,v,z->x,v,w}
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_J23(KS-1,i,j) = 0.0_RP

          qflx_J23(KS  ,i,j) = J23G(KS,i,j,I_XVW) &
                             * 0.25_RP * ( VELY(KS+1,i,j)+VELY(KS,i,j) ) & ! [x,v,z->x,v,w]
                                       * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) )   ! [x,v,z->x,v,w]

          qflx_J23(KE-1,i,j) = J23G(KE-1,i,j,I_XVW) &
                             * 0.25_RP * ( VELY(KE,i,j)+VELY(KE-1,i,j) ) & ! [x,v,z->x,v,w]
                                       * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) )   ! [x,v,z->x,v,w]

          qflx_J23(KE  ,i,j) = 0.0_RP
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
                              * 0.5_RP * ( VELX(k,i,j+1)+VELX(k,i,j) ) &       ! [u,y,z->u,v,z]
                              * ( FACT_N * ( MOMY(k,i+1,j)+MOMY(k,i  ,j) ) &
                                + FACT_F * ( MOMY(k,i+2,j)+MOMY(k,i-1,j) ) ) & ! {x,v,z->u,v,z}
                              + GSQRT(k,i,j,I_UVZ) * num_diff(k,i,j,I_MOMY,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (x, y, z)
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
          qflx_hi(k,i,j-1,YDIR) = GSQRT(k,i,j,I_XYZ) &
                                * 0.5_RP * ( VELY(k,i,j)+VELY(k,i,j-1) ) &       ! [x,v,z->x,y,z]
                                * ( FACT_N * ( MOMY(k,i,j  )+MOMY(k,i,j-1) ) &
                                  + FACT_F * ( MOMY(k,i,j+1)+MOMY(k,i,j-2) ) ) & ! {x,v,z->x,y,z}
                                + GSQRT(k,i,j,I_XYZ) * num_diff(k,i,j,I_MOMY,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! pressure gradient force at (x, v, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          pgf(k,i,j) = ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) & ! [x,y,z]
                       - GSQRT(k,i,j  ,I_XYZ) * DPRES(k,i,j  ) & ! [x,y,z]
                       ) * RFDY(j) &
                     + ( J23G(k  ,i,j,I_XVW) &
                       * 0.25_RP * ( DPRES(k+1,i,j+1)+DPRES(k,i,j+1) &
                                   + DPRES(k+1,i,j  )+DPRES(k,i,j  ) ) & ! [x,y,z->x,v,w]
                       - J23G(k-1,i,j,I_XVW) &
                       * 0.25_RP * ( DPRES(k,i,j+1)+DPRES(k-1,i,j+1) &
                                   + DPRES(k,i,j  )+DPRES(k-1,i,j  ) ) & ! [x,y,z->x,v,w]
                       ) * RCDZ(k)
       enddo
       enddo
       enddo

       ! coriolis force at (x, v, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          cor(k,i,j) = 0.0625_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) & ! [x,y,z->x,v,z]
                                 * ( DENS   (k,i  ,j+1)+DENS   (k,i  ,j) ) & ! [x,y,z->x,v,z]
                                 * ( VELX   (k,i  ,j+1)+VELX   (k,i  ,j) &
                                   + VELX   (k,i-1,j+1)+VELX   (k,i-1,j) )   ! [u,y,z->x,v,z]
       enddo
       enddo
       enddo

       !-----< update momentum (y) >-----

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
                                      + ( qflx_J13(k,i,j) - qflx_J13(k-1,i,j)             ) * RCDZ(k) &
                                      + ( qflx_J23(k,i,j) - qflx_J23(k-1,i,j)             ) * RCDZ(k) &
                                      + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                      + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RFDY(j) ) & ! advection
                                    - pgf(k,i,j)                                                        & ! pressure gradient force
                                  ) / GSQRT(k,i,j,I_XVZ) &
                         + dtrk * ( + cor(k,i,j)                                                        & ! coriolis force
                                    + divdmp_coef * dtrk * FDY(j) * ( DDIV(k,i,j+1)-DDIV(k,i,j) )       & ! divergence damping
                                    + MOMY_t(k,i,j)                                                     ) ! physics tendency
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifndef NO_FCT_DYN
       if ( FLAG_FCT_MOMENTUM ) then

          ! monotonic flux
          ! at (x, v, interface)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k,i,j+1) )
             call CHECK( __LINE__, MOMZ(k,i,j  ) )
             call CHECK( __LINE__, VELY(k+1,i,j) )
             call CHECK( __LINE__, VELY(k  ,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.25_RP * (    ( MOMZ(k,i,j+1)+MOMZ(k,i,j) ) * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
                                             - abs( MOMZ(k,i,j+1)+MOMZ(k,i,j) ) * ( VELY(k+1,i,j)-VELY(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (u, v, layer)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELX(k,i,j) )
             call CHECK( __LINE__, VELX(k,i+1,j) )
             call CHECK( __LINE__, MOMY(k,i,j) )
             call CHECK( __LINE__, MOMY(k,i,j+1) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.25_RP * (    ( MOMX(k,i,j+1)+MOMX(k,i,j) ) * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
                                             - abs( MOMX(k,i,j+1)+MOMX(k,i,j) ) * ( VELY(k,i+1,j)-VELY(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (x, y, layer)
          ! note that y-index is added by -1
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+2
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i,j-1) )
             call CHECK( __LINE__, MOMY(k,i,j  ) )
             call CHECK( __LINE__, VELY(k,i,j-1) )
             call CHECK( __LINE__, VELY(k,i,j  ) )
#endif
             qflx_lo(k,i,j-1,YDIR) = 0.25_RP * (    ( MOMY(k,i,j)+MOMY(k,i,j-1) ) * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
                                               - abs( MOMY(k,i,j)+MOMY(k,i,j-1) ) * ( VELY(k,i,j)-VELY(k,i,j-1) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       endif
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then

       call ATMOS_DYN_fct( qflx_anti,               & ! (out)
                           MOMY0, qflx_hi, qflx_lo, & ! (in)
                           RCDZ, RCDX, RFDY, dtrk   ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(y)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             MOMY_RK(k,i,j) = MOMY_RK(k,i,j) &
                            + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                         + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                         + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RFDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       enddo
       enddo

#ifdef DEBUG
       qflx_lo  (:,:,:,:) = UNDEF
       qflx_anti(:,:,:,:) = UNDEF
#endif
    endif ! FLAG_FCT_MOMENTUM

#ifdef DEBUG
    qflx_hi(KS:,:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#endif

       !########################################################################
       ! Thermodynamic equation
       !########################################################################

       !-----< high order flux >-----

       ! at (x, y, w)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS,   JJE
       do i = IIS,   IIE
       do k = KS+1, KE-2
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(k-1,i,j) )
          call CHECK( __LINE__, POTT(k  ,i,j) )
          call CHECK( __LINE__, POTT(k+1,i,j) )
          call CHECK( __LINE__, POTT(k+2,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,ZDIR) )
#endif
          qflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                              * ( FACT_N * ( POTT(k+1,i,j)+POTT(k  ,i,j) ) &
                                + FACT_F * ( POTT(k+2,i,j)+POTT(k-1,i,j) ) ) & ! [{x,y,z->x,y,w}]
                              + GSQRT(k,i,j,I_XYW) * num_diff(k,i,j,I_RHOT,ZDIR)
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
          call CHECK( __LINE__, mflx_hi(KS,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(KS+1,i,j) )
          call CHECK( __LINE__, POTT(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_RHOT,ZDIR) )
          call CHECK( __LINE__, mflx_hi(KE-1,i,j,ZDIR) )
          call CHECK( __LINE__, POTT(KE-1,i,j) )
          call CHECK( __LINE__, POTT(KE  ,i,j) )
          call CHECK( __LINE__, num_diff(KE-1,i,j,I_RHOT,ZDIR) )
#endif
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_hi(KS  ,i,j,ZDIR) = mflx_hi(KS  ,i,j,ZDIR) * 0.5_RP * ( POTT(KS+1,i,j)+POTT(KS,i,j) ) &
                                 + GSQRT(KS  ,i,j,I_XYW) * num_diff(KS  ,i,j,I_RHOT,ZDIR)
          qflx_hi(KE-1,i,j,ZDIR) = mflx_hi(KE-1,i,j,ZDIR) * 0.5_RP * ( POTT(KE,i,j)+POTT(KE-1,i,j) ) &
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
          call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, POTT(k,i-1,j) )
          call CHECK( __LINE__, POTT(k,i  ,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, POTT(k,i+1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_RHOT,XDIR) )
#endif
          qflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) &
                              * ( FACT_N * ( POTT(k,i+1,j)+POTT(k,i  ,j) ) &
                                + FACT_F * ( POTT(k,i+2,j)+POTT(k,i-1,j) ) ) & ! [{x,y,z->u,y,z}]
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
                                + FACT_F * ( POTT(k,i,j+2)+POTT(k,i,j-1) ) ) & ! [{x,y,z->x,v,z}]
                              + GSQRT(k,i,j,I_XVZ) * num_diff(k,i,j,I_RHOT,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !-----< update rho*theta >-----

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
          call CHECK( __LINE__, RHOT0(k,i,j) )
#endif
          RHOT_RK(k,i,j) = RHOT0(k,i,j) &
                         + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                      + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                      + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                                  ) / GSQRT(k,i,j,I_XYZ) &
                         + dtrk * RHOT_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifndef NO_FCT_DYN
       if ( FLAG_FCT_T ) then

          call COMM_vars8( mflx_hi(:,:,:,ZDIR), 1 )
          call COMM_vars8( mflx_hi(:,:,:,XDIR), 2 )
          call COMM_vars8( mflx_hi(:,:,:,YDIR), 3 )
          call COMM_wait ( mflx_hi(:,:,:,ZDIR), 1 )
          call COMM_wait ( mflx_hi(:,:,:,XDIR), 2 )
          call COMM_wait ( mflx_hi(:,:,:,YDIR), 3 )

          ! monotonic flux
          ! at (x, y, interface)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
             call CHECK( __LINE__, POTT(k  ,i,j) )
             call CHECK( __LINE__, POTT(k+1,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * (     mflx_hi(k,i,j,ZDIR)  * ( POTT(k+1,i,j)+POTT(k,i,j) ) &
                                            - abs(mflx_hi(k,i,j,ZDIR)) * ( POTT(k+1,i,j)-POTT(k,i,j) ) )
          enddo
          enddo
          enddo

          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo

#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (u, y, layer)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
             call CHECK( __LINE__, POTT(k,i  ,j) )
             call CHECK( __LINE__, POTT(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.5_RP * (     mflx_hi(k,i,j,XDIR)  * ( POTT(k,i+1,j)+POTT(k,i,j) ) &
                                            - abs(mflx_hi(k,i,j,XDIR)) * ( POTT(k,i+1,j)-POTT(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

          ! at (x, v, layer)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
             call CHECK( __LINE__, POTT(k,i,j  ) )
             call CHECK( __LINE__, POTT(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.5_RP * (     mflx_hi(k,i,j,YDIR)  * ( POTT(k,i,j+1)+POTT(k,i,j) ) &
                                            - abs(mflx_hi(k,i,j,YDIR)) * ( POTT(k,i,j+1)-POTT(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       endif
#endif

    enddo
    enddo

#ifndef NO_FCT_DYN
    if ( FLAG_FCT_T ) then
       call ATMOS_DYN_fct( qflx_anti,               & ! (out)
                           RHOT0, qflx_hi, qflx_lo, & ! (in)
                           RCDZ, RCDX, RCDY, dtrk   ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update rho*theta
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, RHOT_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             RHOT_RK(k,i,j) = RHOT_RK(k,i,j) &
                            + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                         + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                         + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       enddo
       enddo
    endif
#endif

    return
  end subroutine ATMOS_DYN_rk_heve

end module scale_atmos_dyn_rk_heve
