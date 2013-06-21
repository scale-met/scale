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
  real(RP) :: kappa
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

    if( IO_L ) write(IO_FID_LOG,*) '*** HEVE'

    if ( ATMOS_TYPE_DYN .ne. 'HEVE' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_DYN is not HEVE. Check!'
       call PRC_MPIstop
    end if

#ifdef DRY
    kappa = CPdry / CVdry
#endif

    return
  end subroutine ATMOS_DYN_rk_setup


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
    FDZ, FDX, FDY,                               &
    RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
    dtrk, rk, rko,                               &
    VELZ, VELX, VELY, PRES, POTT                 )
    use mod_const, only : &
       UNDEF  => CONST_UNDEF,  &
       IUNDEF => CONST_UNDEF2, &
       GRAV   => CONST_GRAV,   &
       P00    => CONST_PRE00,  &
       Rdry   => CONST_Rdry
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

    ! flux (work space)
    real(RP) :: qflx_hi  (KA,IA,JA,3)
    real(RP) :: qflx_lo  (KA,IA,JA,3)
    real(RP) :: qflx_anti(KA,IA,JA,3)

#ifndef NO_FCT_DYN
    real(RP), pointer :: org(:,:,:)
    real(RP), target  :: work(KA,IA,JA)
#endif

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

    qflx_lo(KS:,:,:,:) = UNDEF
    qflx_hi(KS:,:,:,:) = UNDEF

#ifndef NO_FCT_DYN
    work(:,:,:) = UNDEF
#endif
#endif

#ifndef NO_FCT_DYN
    if ( FLAG_FCT_RHO ) then
       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = DENS0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => DENS0
       end if
    end if
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
       enddo
       enddo

#ifndef NO_FCT_DYN

    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
#ifdef _USE_RDMA
       call COMM_rdma_vars8( 5+QA+16, 3 )
#else
       call COMM_vars8( VELZ(:,:,:), 1 )
       call COMM_vars8( VELX(:,:,:), 2 )
       call COMM_vars8( VELY(:,:,:), 3 )
       call COMM_wait ( VELZ(:,:,:), 1 )
       call COMM_wait ( VELX(:,:,:), 2 )
       call COMM_wait ( VELY(:,:,:), 3 )
#endif
    end if


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1
#endif
       ! pressure, pott. temp.
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
          PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rtot(k,i,j) / P00 )**((CVtot(k,i,j)+Rtot(k,i,j))/CVtot(k,i,j))
#endif
       enddo
       enddo
       enddo
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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

#ifndef NO_FCT_DYN
       ! monotonic flux
       if ( FLAG_FCT_RHO ) then
          ! at (x, y, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, VELZ(k,i,j) )
             call CHECK( __LINE__, DENS(k  ,i,j) )
             call CHECK( __LINE__, DENS(k+1,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  VELZ(k,i,j)  * ( DENS(k+1,i,j)+DENS(k,i,j) ) &
             -abs(VELZ(k,i,j)) * ( DENS(k+1,i,j)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
          ! at (u, y, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELX(k,i,j) )
             call CHECK( __LINE__, DENS(k,i  ,j) )
             call CHECK( __LINE__, DENS(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.5_RP * ( &
                  VELX(k,i,j)  * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
             -abs(VELX(k,i,j)) * ( DENS(k,i+1,j)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELY(k,i,j) )
             call CHECK( __LINE__, DENS(k,i,j  ) )
             call CHECK( __LINE__, DENS(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.5_RP * ( &
                  VELY(k,i,j)  * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
             -abs(VELY(k,i,j)) * ( DENS(k,i,j+1)-DENS(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       endif
#endif

       ! at (x, y, interface)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(k+1,i,j) )
          call CHECK( __LINE__, MOMZ(k  ,i,j) )
          call CHECK( __LINE__, MOMZ(k-1,i,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,ZDIR) )
#endif
!          mflx_hi(k,i,j,ZDIR) = ( -MOMZ(k+1,i,j) + 8.0_RP*MOMZ(k,i,j) - MOMZ(k-1,i,j) ) / 6.0_RP &
          mflx_hi(k,i,j,ZDIR) = MOMZ(k,i,j) &
                              + num_diff(k,i,j,I_DENS,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ(KS+1,i,j) )
          call CHECK( __LINE__, MOMZ(KS  ,i,j) )
          call CHECK( __LINE__, num_diff(KS,i,j,I_DENS,ZDIR) )
#endif
          mflx_hi(KE,i,j,ZDIR) = 0.0_RP
!          mflx_hi(KS,i,j,ZDIR) = ( -MOMZ(KS+1,i,j) + 8.0_RP*MOMZ(KS,i,j) ) / 6.0_RP & ! MOMZ(KS-1,i,j) == 0
          mflx_hi(KS,i,j,ZDIR) = MOMZ(KS,i,j) &
                               + num_diff(KS,i,j,I_DENS,ZDIR)
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
                          + DENS_t(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
#ifndef NO_FCT_DYN
    enddo
    enddo

    if ( FLAG_FCT_RHO ) then
       call ATMOS_DYN_fct( &
            qflx_anti,          & ! (out)
            org, mflx_hi, qflx_lo,  & ! (in)
            RCDZ, RCDX, RCDY, dtrk  ) ! (in)
#ifdef DEBUG
       qflx_lo(KS:,:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update rho
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
          !--- update mflx_hi
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
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
       enddo
       enddo

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#ifdef DEBUG
    qflx_hi(KS:,:,:,:) = UNDEF
#endif


#ifdef DEBUG
    work(:,:,:) = UNDEF
#endif
    if ( FLAG_FCT_MOMENTUM ) then
       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = MOMZ0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => MOMZ0
       end if
    end if


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (z) #####

       if ( FLAG_FCT_MOMENTUM ) then

          ! monotonic flux
          ! at (x, y, layer)
          ! note than z-index is added by -1
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS+1, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k-1,i,j) )
             call CHECK( __LINE__, MOMZ(k  ,i,j) )
             call CHECK( __LINE__, VELZ(k-1,i,j) )
             call CHECK( __LINE__, VELZ(k  ,i,j) )
#endif
             vel = ( MOMZ(k,i,j)+MOMZ(k-1,i,j) ) * 0.5_RP
             qflx_lo(k-1,i,j,ZDIR) = 0.5_RP * ( &
                       ( vel ) * ( VELZ(k,i,j)+VELZ(k-1,i,j) ) &
                  - abs( vel ) * ( VELZ(k,i,j)-VELZ(k-1,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             ! k = KE
             qflx_lo(KE-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, y, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMX(k  ,i,j) )
             call CHECK( __LINE__, MOMX(k+1,i,j) )
             call CHECK( __LINE__, VELZ(k,i,j) )
             call CHECK( __LINE__, VELZ(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.25_RP * ( &
                       ( MOMX(k+1,i,j)+MOMX(k,i,j) ) * ( VELZ(k,i+1,j)+VELZ(k,i,j) ) &
                  - abs( MOMX(k+1,i,j)+MOMX(k,i,j) ) * ( VELZ(k,i+1,j)-VELZ(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
             qflx_lo(KE,i,j,XDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, v, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i,j  ) )
             call CHECK( __LINE__, MOMY(k+1,i,j) )
             call CHECK( __LINE__, VELZ(k,i,j) )
             call CHECK( __LINE__, VELZ(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.25_RP * ( &
                      ( MOMY(k+1,i,j)+MOMY(k,i,j) ) * ( VELZ(k,i,j+1)+VELZ(k,i,j) ) &
                  -abs( MOMY(k+1,i,j)+MOMY(k,i,j) ) * ( VELZ(k,i,j+1)-VELZ(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE,i,j,YDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if
#endif

       ! high order flux
       ! at (x, y, layer)
       ! note than z-index is added by -1
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
          qflx_hi(KE,i,j,XDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (x, v, interface)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
          qflx_hi(KE,i,j,YDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update momentum(z)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
          MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
               + dtrk * ( - ( ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
                            + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                          - ( PRES(k+1,i,j)-PRES(k,i,j) ) * RFDZ(k)                      & ! pressure gradient force
                          - ( DENS(k+1,i,j)+DENS(k,i,j) ) * 0.5_RP * GRAV                & ! gravity force
                          + divdmp_coef * dtrk  * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) & ! divergence damping
                          + MOMZ_t(k,i,j) )
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
       call ATMOS_DYN_fct( &
            qflx_anti,              & ! (out)
            org, qflx_hi, qflx_lo,  & ! (in)
            RFDZ, RCDX, RCDY, dtrk  ) ! (in)
#ifdef DEBUG
       qflx_lo(KS:,:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(z)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ_RK(k,i,j) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             MOMZ_RK(k,i,j) = MOMZ_RK(k,i,j) &
                  + dtrk * ( - ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
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
#ifdef DEBUG
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if

#ifdef DEBUG
    qflx_hi(KS:,:,:,:) = UNDEF
    work(:,:,:) = UNDEF
#endif
    if ( FLAG_FCT_MOMENTUM ) then

       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = MOMX0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => MOMX0
       end if
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (x) #####

       ! monotonic flux
       if ( FLAG_FCT_MOMENTUM ) then
          ! at (u, y, interface)
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k,i  ,j) )
             call CHECK( __LINE__, MOMZ(k,i+1,j) )
             call CHECK( __LINE__, VELX(k  ,i,j) )
             call CHECK( __LINE__, VELX(k+1,i,j) )
#endif
             vel = 0.5_RP * ( MOMZ(k,i+1,j) + MOMZ(k,i,j) )
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  ( vel ) * ( VELX(k+1,i,j)+VELX(k,i,j) ) &
              -abs( vel ) * ( VELX(k+1,i,j)-VELX(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KS-1,i,j,ZDIR) = 0.0_RP
             qflx_lo(KE,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, y, layer)
          ! note that x-index is added by -1
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+2
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMX(k,i-1,j) )
             call CHECK( __LINE__, MOMX(k,i  ,j) )
             call CHECK( __LINE__, VELX(k,i-1,j) )
             call CHECK( __LINE__, VELX(k,i  ,j) )
#endif
             vel = MOMX(k,i,j)+MOMX(k,i-1,j)
             qflx_lo(k,i-1,j,XDIR) = 0.25_RP * ( &
                  ( vel ) * ( VELX(k,i,j)+VELX(k,i-1,j) ) &
              -abs( vel ) * ( VELX(k,i,j)-VELX(k,i-1,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i+1,j) )
             call CHECK( __LINE__, MOMY(k,i  ,j) )
             call CHECK( __LINE__, VELX(k,i,j  ) )
             call CHECK( __LINE__, VELX(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.25_RP * ( &
                  ( MOMY(k,i+1,j)+MOMY(k,i,j) ) * ( VELX(k,i,j+1)+VELX(k,i,j) ) &
              -abs( MOMY(k,i+1,j)+MOMY(k,i,j) ) * ( VELX(k,i,j+1)-VELX(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if
#endif

       ! high order flux
       ! at (u, y, interface)
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
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
          vel = 0.5_RP * ( VELZ(KS,i+1,j) + VELZ(KS,i,j) )
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * vel & ! just above the bottom boundary
                                 * ( MOMX(KS+1,i,j)+MOMX(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_MOMX,ZDIR)
          vel = 0.5_RP * ( VELZ(KE-1,i+1,j) + VELZ(KE-1,i,j) )
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * vel & ! just below the top boundary
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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


#ifndef NO_FCT_DYN
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
       call ATMOS_DYN_fct( &
            qflx_anti,              & ! (out)
            org, qflx_hi, qflx_lo,  & ! (in)
            RCDZ, RFDX, RCDY, dtrk  ) ! (in)
#ifdef DEBUG
       qflx_lo(KS:,:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(x)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#ifdef DEBUG
    qflx_hi(KS:,:,:,:) = UNDEF
    work(:,:,:) = UNDEF
#endif
    if ( FLAG_FCT_MOMENTUM ) then

       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = MOMY0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => MOMY0
       end if
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (y) #####
       ! monotonic flux
       if ( FLAG_FCT_MOMENTUM ) then
          ! at (x, v, interface)
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, MOMZ(k,i,j+1) )
             call CHECK( __LINE__, MOMZ(k,i,j  ) )
             call CHECK( __LINE__, VELY(k+1,i,j) )
             call CHECK( __LINE__, VELY(k  ,i,j) )
#endif
             vel = 0.5_RP * ( MOMZ(k,i,j+1) + MOMZ(k,i,j) )
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  ( vel ) * ( VELY(k+1,i,j)+VELY(k,i,j) ) &
              -abs( vel ) * ( VELY(k+1,i,j)-VELY(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, VELX(k,i,j) )
             call CHECK( __LINE__, VELX(k,i+1,j) )
             call CHECK( __LINE__, MOMY(k,i,j) )
             call CHECK( __LINE__, MOMY(k,i,j+1) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.25_RP * ( &
                  ( MOMX(k,i,j+1)+MOMX(k,i,j) ) * ( VELY(k,i+1,j)+VELY(k,i,j) ) &
              -abs( MOMX(k,i,j+1)+MOMX(k,i,j) ) * ( VELY(k,i+1,j)-VELY(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, y, layer)
          ! note that y-index is added by -1
          !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+2
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, MOMY(k,i,j-1) )
             call CHECK( __LINE__, MOMY(k,i,j  ) )
             call CHECK( __LINE__, VELY(k,i,j-1) )
             call CHECK( __LINE__, VELY(k,i,j  ) )
#endif
             vel = MOMY(k,i,j)+MOMY(k,i,j-1)
             qflx_lo(k,i,j-1,YDIR) = 0.25_RP * ( &
                  ( vel ) * ( VELY(k,i,j)+VELY(k,i,j-1) ) &
              -abs( vel ) * ( VELY(k,i,j)-VELY(k,i,j-1) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if
#endif

       ! high order flux
       ! at (x, v, interface)
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k,vel) schedule(static,1) collapse(2)
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
          vel = 0.5_RP * ( VELZ(KS,i,j+1) + VELZ(KS,i,j) )
          ! just above the bottom boundary
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * vel * ( MOMY(KS+1,i,j)+MOMY(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_MOMY,ZDIR)
          vel = 0.5_RP * ( VELZ(KE-1,i,j+1) + VELZ(KE-1,i,j) )
          ! just below the top boundary
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * vel * ( MOMY(KE,i,j)+MOMY(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_MOMY,ZDIR)
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, v, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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

#ifndef NO_FCT_DYN
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
       call ATMOS_DYN_fct( &
            qflx_anti,              & ! (out)
            org, qflx_hi, qflx_lo,  & ! (in)
            RCDZ, RCDX, RFDY, dtrk  ) ! (in)
#ifdef DEBUG
       qflx_lo(KS:,:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(y)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#ifdef DEBUG
    qflx_hi(KS:,:,:,:) = UNDEF
    work(:,:,:) = UNDEF
#endif

    if ( FLAG_FCT_T ) then

#ifdef _USE_RDMA
       call COMM_rdma_vars8( 5+QA+11, 3 )
#else
       call COMM_vars8( mflx_hi(:,:,:,ZDIR), 1 )
       call COMM_vars8( mflx_hi(:,:,:,XDIR), 2 )
       call COMM_vars8( mflx_hi(:,:,:,YDIR), 3 )
       call COMM_wait ( mflx_hi(:,:,:,ZDIR), 1 )
       call COMM_wait ( mflx_hi(:,:,:,XDIR), 2 )
       call COMM_wait ( mflx_hi(:,:,:,YDIR), 3 )
#endif
       if ( rko == RK ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          !OCL XFILL
          do j = JS-1, JE+1
          do i = IS-1, IE+1
          do k = KS, KE
             work(k,i,j) = RHOT0(k,i,j)
          enddo
          enddo
          enddo
          org => work
       else
          org => RHOT0
       end if
    end if

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### Thermodynamic Equation #####

       ! monotonic flux
       if ( FLAG_FCT_T ) then
          ! at (x, y, interface)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE-1
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,ZDIR) )
             call CHECK( __LINE__, POTT(k  ,i,j) )
             call CHECK( __LINE__, POTT(k+1,i,j) )
#endif
             qflx_lo(k,i,j,ZDIR) = 0.5_RP * ( &
                  mflx_hi(k,i,j,ZDIR)  * ( POTT(k+1,i,j)+POTT(k,i,j) ) &
             -abs(mflx_hi(k,i,j,ZDIR)) * ( POTT(k+1,i,j)-POTT(k,i,j) ) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-1, IIE+1
             qflx_lo(KE  ,i,j,ZDIR) = 0.0_RP
          enddo
          enddo

#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (u, y, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE+1
          do i = IIS-2, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
             call CHECK( __LINE__, POTT(k,i  ,j) )
             call CHECK( __LINE__, POTT(k,i+1,j) )
#endif
             qflx_lo(k,i,j,XDIR) = 0.5_RP * ( &
                  mflx_hi(k,i,j,XDIR)  * ( POTT(k,i+1,j)+POTT(k,i,j) ) &
             -abs(mflx_hi(k,i,j,XDIR)) * ( POTT(k,i+1,j)-POTT(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
          ! at (x, v, layer)
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-2, JJE+1
          do i = IIS-1, IIE+1
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
             call CHECK( __LINE__, POTT(k,i,j  ) )
             call CHECK( __LINE__, POTT(k,i,j+1) )
#endif
             qflx_lo(k,i,j,YDIR) = 0.5_RP * ( &
                  mflx_hi(k,i,j,YDIR)  * ( POTT(k,i,j+1)+POTT(k,i,j) ) &
             -abs(mflx_hi(k,i,j,YDIR)) * ( POTT(k,i,j+1)-POTT(k,i,j) ) )
          enddo
          enddo
          enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       end if
#endif

       ! at (x, y, interface)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
                                + FACT_F * ( POTT(k+2,i,j)+POTT(k-1,i,j) ) ) &
                              + num_diff(k,i,j,I_RHOT,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
          qflx_hi(KS  ,i,j,ZDIR) = 0.5_RP * mflx_hi(KS  ,i,j,ZDIR)  &      ! just above the bottom boundary
                                 * ( POTT(KS+1,i,j)+POTT(KS,i,j) ) &
                                 + num_diff(KS  ,i,j,I_RHOT,ZDIR)
          qflx_hi(KE-1,i,j,ZDIR) = 0.5_RP * mflx_hi(KE-1,i,j,ZDIR)  &      ! just below the top boundary
                                 * ( POTT(KE,i,j)+POTT(KE-1,i,j) ) &
                                 + num_diff(KE-1,i,j,I_RHOT,ZDIR)
          qflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       ! at (u, y, layer)
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
                            + RHOT_t(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifndef NO_FCT_DYN
    enddo
    enddo

    if ( FLAG_FCT_T ) then
       call ATMOS_DYN_fct( &
            qflx_anti,              & ! (out)
            org, qflx_hi, qflx_lo,  & ! (in)
            RCDZ, RCDX, RCDY, dtrk  ) ! (in)
#ifdef DEBUG
       qflx_lo(KS:,:,:,:) = UNDEF
#endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update rho*theta
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
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
#ifdef DEBUG
       qflx_anti(:,:,:,:) = UNDEF
#endif
    end if
#endif
#ifdef DEBUG
    qflx_hi(KS:,:,:,:) = UNDEF
#endif


  end subroutine ATMOS_DYN_rk


end module mod_atmos_dyn_rk
