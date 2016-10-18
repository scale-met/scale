!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          HEVE FVM scheme for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-18 (S.Nishizawa) [new] splited from dynamical core
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_dyn_tstep_short_fvm_heve
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_index
  use scale_tracer
  use scale_process
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
  public :: ATMOS_DYN_Tstep_short_fvm_heve_regist
  public :: ATMOS_DYN_Tstep_short_fvm_heve_setup
  public :: ATMOS_DYN_Tstep_short_fvm_heve

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
#if 1
#define F2H(k,p,idx) (CDZ(k+p-1)*GSQRT(k+p-1,i,j,idx)/(CDZ(k)*GSQRT(k,i,j,idx)+CDZ(k+1)*GSQRT(k+1,i,j,idx)))
#else
#define F2H(k,p,idx) 0.5_RP
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: VA_FVM_HEVE = 0

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tstep_short_fvm_heve_regist( &
       ATMOS_DYN_TYPE, &
       VA_out,         &
       VAR_NAME,       &
       VAR_DESC,       &
       VAR_UNIT        )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*),       intent(in)  :: ATMOS_DYN_TYPE
    integer,                intent(out) :: VA_out         !< number of prognostic variables
    character(len=H_SHORT), intent(out) :: VAR_NAME(:)    !< name   of the variables
    character(len=H_MID),   intent(out) :: VAR_DESC(:)    !< desc.  of the variables
    character(len=H_SHORT), intent(out) :: VAR_UNIT(:)    !< unit   of the variables
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HEVE Register'

    if ( ATMOS_DYN_TYPE /= 'FVM-HEVE' .AND. ATMOS_DYN_TYPE /= 'HEVE' ) then
       write(*,*) 'xxx ATMOS_DYN_TYPE is not FVM-HEVE. Check!'
       call PRC_MPIstop
    endif

    VA_out      = VA_FVM_HEVE
    VAR_NAME(:) = ""
    VAR_DESC(:) = ""
    VAR_UNIT(:) = ""

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_heve_regist

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_short_fvm_heve_setup
    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_heve_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_Tstep_short_fvm_heve( &
       DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
       PROG_RK,                                     &
       mflx_hi, tflx_hi,                            &
       DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
       DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
       PROG0, PROG,                                 &
       DPRES0, RT2P, CORIOLI,                       &
       num_diff, divdmp_coef, DDIV,                 &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T,               &
       FLAG_FCT_ALONG_STREAM,                       &
       CDZ, FDZ, FDX, FDY,                          &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
       PHI, GSQRT, J13G, J23G, J33G, MAPF,          &
       REF_dens, REF_rhot,                          &
       BND_W, BND_E, BND_S, BND_N,                  &
       dtrk, dt                                     )
    use scale_grid_index
    use scale_const, only: &
       GRAV   => CONST_GRAV,  &
       P00    => CONST_PRE00
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_fct
    use scale_atmos_dyn_fvm_flux_ud1, only: &
       ATMOS_DYN_FVM_fluxZ_XYW_ud1, &
       ATMOS_DYN_FVM_fluxX_XYW_ud1, &
       ATMOS_DYN_FVM_fluxY_XYW_ud1, &
       ATMOS_DYN_FVM_fluxZ_UYZ_ud1, &
       ATMOS_DYN_FVM_fluxX_UYZ_ud1, &
       ATMOS_DYN_FVM_fluxY_UYZ_ud1, &
       ATMOS_DYN_FVM_fluxZ_XVZ_ud1, &
       ATMOS_DYN_FVM_fluxX_XVZ_ud1, &
       ATMOS_DYN_FVM_fluxY_XVZ_ud1, &
       ATMOS_DYN_FVM_fluxZ_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxX_XYZ_ud1, &
       ATMOS_DYN_FVM_fluxY_XYZ_ud1
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_fluxZ_XYZ,   &
       ATMOS_DYN_FVM_fluxX_XYZ,   &
       ATMOS_DYN_FVM_fluxY_XYZ,   &
       ATMOS_DYN_FVM_fluxZ_XYW,   &
       ATMOS_DYN_FVM_fluxJ13_XYW, &
       ATMOS_DYN_FVM_fluxJ23_XYW, &
       ATMOS_DYN_FVM_fluxX_XYW,   &
       ATMOS_DYN_FVM_fluxY_XYW,   &
       ATMOS_DYN_FVM_fluxZ_UYZ,   &
       ATMOS_DYN_FVM_fluxJ13_UYZ, &
       ATMOS_DYN_FVM_fluxJ23_UYZ, &
       ATMOS_DYN_FVM_fluxX_UYZ,   &
       ATMOS_DYN_FVM_fluxY_UYZ,   &
       ATMOS_DYN_FVM_fluxZ_XVZ,   &
       ATMOS_DYN_FVM_fluxJ13_XVZ, &
       ATMOS_DYN_FVM_fluxJ23_XVZ, &
       ATMOS_DYN_FVM_fluxX_XVZ,   &
       ATMOS_DYN_FVM_fluxY_XVZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY,  &
       I_UY,  &
       I_XV,  &
       I_UV
#ifdef HIST_TEND
    use scale_history, only: &
       HIST_in
#endif
    implicit none

    real(RP), intent(out)        :: DENS_RK (KA,IA,JA)    ! prognostic variables
    real(RP), intent(out)        :: MOMZ_RK (KA,IA,JA)    !
    real(RP), intent(out)        :: MOMX_RK (KA,IA,JA)    !
    real(RP), intent(out)        :: MOMY_RK (KA,IA,JA)    !
    real(RP), intent(out)        :: RHOT_RK (KA,IA,JA)    !
    real(RP), intent(out)        :: PROG_RK (KA,IA,JA,VA) !

    real(RP), intent(inout)      :: mflx_hi (KA,IA,JA,3)  ! mass flux
    real(RP), intent(out)        :: tflx_hi (KA,IA,JA,3)  ! internal energy flux

    real(RP), intent(in), target :: DENS0   (KA,IA,JA)    ! prognostic variables at previous dynamical time step
    real(RP), intent(in), target :: MOMZ0   (KA,IA,JA)    !
    real(RP), intent(in), target :: MOMX0   (KA,IA,JA)    !
    real(RP), intent(in), target :: MOMY0   (KA,IA,JA)    !
    real(RP), intent(in), target :: RHOT0   (KA,IA,JA)    !
    real(RP), intent(in)         :: PROG0   (KA,IA,JA,VA)

    real(RP), intent(in)         :: DENS    (KA,IA,JA)    ! prognostic variables at previous RK step
    real(RP), intent(in)         :: MOMZ    (KA,IA,JA)    !
    real(RP), intent(in)         :: MOMX    (KA,IA,JA)    !
    real(RP), intent(in)         :: MOMY    (KA,IA,JA)    !
    real(RP), intent(in)         :: RHOT    (KA,IA,JA)    !
    real(RP), intent(in)         :: PROG    (KA,IA,JA,VA)

    real(RP), intent(in)         :: DENS_t  (KA,IA,JA)    ! tendency
    real(RP), intent(in)         :: MOMZ_t  (KA,IA,JA)    !
    real(RP), intent(in)         :: MOMX_t  (KA,IA,JA)    !
    real(RP), intent(in)         :: MOMY_t  (KA,IA,JA)    !
    real(RP), intent(in)         :: RHOT_t  (KA,IA,JA)    !

    real(RP), intent(in)         :: DPRES0  (KA,IA,JA)
    real(RP), intent(in)         :: RT2P    (KA,IA,JA)
    real(RP), intent(in)         :: CORIOLI (1, IA,JA)
    real(RP), intent(in)         :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in)         :: divdmp_coef
    real(RP), intent(in)         :: DDIV    (KA,IA,JA)

    logical,  intent(in)         :: FLAG_FCT_MOMENTUM
    logical,  intent(in)         :: FLAG_FCT_T
    logical,  intent(in)         :: FLAG_FCT_ALONG_STREAM

    real(RP), intent(in)         :: CDZ (KA)
    real(RP), intent(in)         :: FDZ (KA-1)
    real(RP), intent(in)         :: FDX (IA-1)
    real(RP), intent(in)         :: FDY (JA-1)
    real(RP), intent(in)         :: RCDZ(KA)
    real(RP), intent(in)         :: RCDX(IA)
    real(RP), intent(in)         :: RCDY(JA)
    real(RP), intent(in)         :: RFDZ(KA-1)
    real(RP), intent(in)         :: RFDX(IA-1)
    real(RP), intent(in)         :: RFDY(JA-1)

    real(RP), intent(in)         :: PHI     (KA,IA,JA)   !< geopotential
    real(RP), intent(in)         :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)         :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)         :: J23G    (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)         :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)         :: MAPF    (IA,JA,2,4)  !< map factor
    real(RP), intent(in)         :: REF_dens(KA,IA,JA)   !< reference density
    real(RP), intent(in)         :: REF_rhot(KA,IA,JA)

    logical,  intent(in)         :: BND_W
    logical,  intent(in)         :: BND_E
    logical,  intent(in)         :: BND_S
    logical,  intent(in)         :: BND_N

    real(RP), intent(in)         :: dtrk
    real(RP), intent(in)         :: dt

    ! diagnostic variables
    real(RP) :: VELZ (KA,IA,JA) ! velocity w [m/s]
    real(RP) :: VELX (KA,IA,JA) ! velocity u [m/s]
    real(RP) :: VELY (KA,IA,JA) ! velocity v [m/s]
    real(RP) :: POTT (KA,IA,JA) ! potential temperature [K]
    real(RP) :: DPRES(KA,IA,JA) ! pressure - reference pressure

    real(RP) :: qflx_J13(KA,IA,JA)
    real(RP) :: qflx_J23(KA,IA,JA)
    real(RP) :: pgf     (KA,IA,JA)  ! pressure gradient force
    real(RP) :: buoy    (KA,IA,JA)  ! buoyancy force
    real(RP) :: cor     (KA,IA,JA)  ! Coriolis force

    ! flux
    real(RP) :: qflx_hi  (KA,IA,JA,3)
#ifndef NO_FCT_DYN
    real(RP) :: qflx_lo  (KA,IA,JA,3)
    real(RP) :: qflx_anti(KA,IA,JA,3)
    real(RP) :: tflx_lo  (KA,IA,JA,3)
    real(RP) :: tflx_anti(KA,IA,JA,3)
    real(RP) :: DENS0_uvw(KA,IA,JA)
    real(RP) :: DENS_uvw (KA,IA,JA)
#endif
    real(RP) :: advch ! horizontal advection
    real(RP) :: advcv ! vertical advection
    real(RP) :: div  ! divergence damping
#ifdef HIST_TEND
    real(RP) :: advch_t(KA,IA,JA,5)
    real(RP) :: advcv_t(KA,IA,JA,5)
    real(RP) :: ddiv_t(KA,IA,JA,3)
    real(RP) :: pg_t(KA,IA,JA,3)
    real(RP) :: cf_t(KA,IA,JA,2)
    logical  :: lhist
#endif

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: IFS_OFF, JFS_OFF
    integer  :: k, i, j
    !---------------------------------------------------------------------------

#ifdef DEBUG
    VELZ(:,:,:) = UNDEF
    VELX(:,:,:) = UNDEF
    VELY(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF

    DPRES(:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    tflx_hi(:,:,:,:) = UNDEF
    qflx_hi(:,:,:,:) = UNDEF

#ifndef NO_FCT_DYN
    qflx_lo  (:,:,:,:) = UNDEF
    qflx_anti(:,:,:,:) = UNDEF
    tflx_lo  (:,:,:,:) = UNDEF
    tflx_anti(:,:,:,:) = UNDEF
#endif
#endif

#ifdef HIST_TEND
    advch_t = 0.0_RP
    advcv_t = 0.0_RP
    ddiv_t = 0.0_RP
    pg_t = 0.0_RP
    cf_t = 0.0_RP

    lhist = dt .eq. dtrk
#endif

    IFS_OFF = 1
    JFS_OFF = 1
    if ( BND_W ) IFS_OFF = 0
    if ( BND_S ) JFS_OFF = 0


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! pressure, pott. temp.

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-2, JJE+2
       do i = IIS-2, IIE+2
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, DPRES0(k,i,j) )
             call CHECK( __LINE__, RT2P(k,i,j) )
             call CHECK( __LINE__, RHOT(k,i,j) )
             call CHECK( __LINE__, REF_rhot(k,i,j) )
#endif
             DPRES(k,i,j) = DPRES0(k,i,j) + RT2P(k,i,j) * ( RHOT(k,i,j) - REF_rhot(k,i,j) )
          enddo
          DPRES(KS-1,i,j) = DPRES0(KS-1,i,j) - DENS(KS,i,j) * ( PHI(KS-1,i,j) - PHI(KS+1,i,j) )
          DPRES(KE+1,i,j) = DPRES0(KE+1,i,j) - DENS(KE,i,j) * ( PHI(KE+1,i,j) - PHI(KE-1,i,j) )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       do j = JJS-JHALO, JJE+JHALO
       do i = IIS-IHALO, IIE+IHALO
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

#ifndef NO_FCT_DYN
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then
       ! momentum -> velocity
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-1, IE+2
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, DENS0(k  ,i,j) )
          call CHECK( __LINE__, DENS0(k+1,i,j) )
#endif
          VELZ(k,i,j) = 2.0_RP * MOMZ0(k,i,j) / ( DENS0(k+1,i,j)+DENS0(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-1, IE+2
          VELZ(KE,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+2
       do i = IS-2, IE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMX0(k,i,j) )
          call CHECK( __LINE__, DENS0(k,i  ,j) )
          call CHECK( __LINE__, DENS0(k,i+1,j) )
#endif
          VELX(k,i,j) = 2.0_RP * MOMX0(k,i,j) / ( DENS0(k,i+1,j)+DENS0(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-2, JE+1
       do i = IS-1, IE+2
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMY0(k,i,j) )
          call CHECK( __LINE__, DENS0(k,i,j  ) )
          call CHECK( __LINE__, DENS0(k,i,j+1) )
#endif
          VELY(k,i,j) = 2.0_RP * MOMY0(k,i,j) / ( DENS0(k,i,j+1)+DENS0(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       call COMM_vars8( VELZ(:,:,:), 4 )
       call COMM_vars8( VELX(:,:,:), 5 )
       call COMM_vars8( VELY(:,:,:), 6 )
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
          mflx_hi(k,i,j,ZDIR) = J33G * MOMZ(k,i,j) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) ) &
                              + J13G(k,i,j,I_XYW) * 0.25_RP * ( MOMX(k+1,i,j)+MOMX(k+1,i-1,j) &
                                                              + MOMX(k  ,i,j)+MOMX(k  ,i-1,j) ) &
                              / MAPF(i,j,2,I_XY) & ! [{u,y,z->x,y,w}]
                              + J23G(k,i,j,I_XYW) * 0.25_RP * ( MOMY(k+1,i,j)+MOMY(k+1,i,j-1) &
                                                              + MOMY(k  ,i,j)+MOMY(k  ,i,j-1) ) &
                              / MAPF(i,j,1,I_XY) & ! [{x,v,z->x,y,w}]
                              + GSQRT(k,i,j,I_XYW) * num_diff(k,i,j,I_DENS,ZDIR) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) )
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
       do i = IIS-IFS_OFF, min(IIE,IEH)
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMX(k,i+1,j) )
          call CHECK( __LINE__, MOMX(k,i  ,j) )
          call CHECK( __LINE__, MOMX(k,i-1,j) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,XDIR) )
#endif
          mflx_hi(k,i,j,XDIR) = GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY) &
                              * ( MOMX(k,i,j) + num_diff(k,i,j,I_DENS,XDIR) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       ! at (x, v, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-JFS_OFF, min(JJE,JEH)
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMY(k,i,j+1) )
          call CHECK( __LINE__, MOMY(k,i,j  ) )
          call CHECK( __LINE__, MOMY(k,i,j-1) )
          call CHECK( __LINE__, num_diff(k,i,j,I_DENS,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV) &
                              * ( MOMY(k,i,j) + num_diff(k,i,j,I_DENS,YDIR) )
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
          advcv = - ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i,  j,  ZDIR) ) * RCDZ(k)
          advch = - ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                  - ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j)
          DENS_RK(k,i,j) = DENS0(k,i,j) &
                         + dtrk * ( ( advcv + advch ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) &
                                  + DENS_t(k,i,j) )
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_DENS) = advcv * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
             advch_t(k,i,j,I_DENS) = advch * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
          endif
#endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


       !########################################################################
       ! momentum equation (z)
       !########################################################################

       !-----< high order flux >-----

       ! at (x, y, z)
       ! note than z-index is added by -1
       call ATMOS_DYN_FVM_fluxZ_XYW( qflx_hi(:,:,:,ZDIR), & ! (out)
            MOMZ, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), J33G, & ! (in)
            num_diff(:,:,:,I_MOMZ,ZDIR), & ! (in)
            CDZ, FDZ, dtrk, &
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ13_XYW( qflx_J13, & ! (out)
            MOMX, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), J13G(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            CDZ, &
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ23_XYW( qflx_J23, & ! (out)
            MOMY, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), J23G(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            CDZ, &
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, y, w)
       call ATMOS_DYN_FVM_fluxX_XYW( qflx_hi(:,:,:,XDIR), & ! (out)
            MOMX, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), & ! (in)
            num_diff(:,:,:,I_MOMZ,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, v, w)
       call ATMOS_DYN_FVM_fluxY_XYW( qflx_hi(:,:,:,YDIR), & ! (out)
            MOMY, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), & ! (in)
            num_diff(:,:,:,I_MOMZ,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

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
          buoy(k,i,j) = GRAV * GSQRT(k,i,j,I_XYW) &
               * ( F2H(k,1,I_XYZ) * ( DENS(k+1,i,j)-REF_dens(k+1,i,j) ) &
                 + F2H(k,2,I_XYZ) * ( DENS(k  ,i,j)-REF_dens(k  ,i,j) ) )
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
          advcv = - ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RFDZ(k)
          advch = - ( ( qflx_J13(k,i,j) - qflx_J13(k-1,i,j) &
                      + qflx_J23(k,i,j) - qflx_J23(k-1,i,j) ) * RFDZ(k) &
                    + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k,i-1,j,XDIR) ) * RCDX(i) &
                    + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k,i,j-1,YDIR) ) * RCDY(j) ) &
                  * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY)
          div = divdmp_coef / dtrk * FDZ(k) * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) ! divergence damping
          MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                         + dtrk * ( ( advcv + advch        &
                                    - pgf (k,i,j)          & ! pressure gradient force
                                    - buoy(k,i,j)          & ! buoyancy force
                                    ) / GSQRT(k,i,j,I_XYW) &
                                  + div                    &
                                  + MOMZ_t(k,i,j) )        ! physics tendency
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMZ) = advcv / GSQRT(k,i,j,I_XYW)
             advch_t(k,i,j,I_MOMZ) = advch / GSQRT(k,i,j,I_XYW)
             pg_t(k,i,j,1) = ( - pgf(k,i,j) - buoy(k,i,j) ) / GSQRT(k,i,j,I_XYW)
             ddiv_t(k,i,j,1) = div
          endif
#endif
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
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(KE,i,j,I_MOMZ) = 0.0_RP
             advch_t(KE,i,j,I_MOMZ) = 0.0_RP
             pg_t(KE,i,j,1) = 0.0_RP
             ddiv_t(KE,i,j,1) = 0.0_RP
          endif
#endif
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
          call ATMOS_DYN_FVM_fluxZ_XYW_ud1( qflx_lo(:,:,:,ZDIR), & ! (out)
               MOMZ, MOMZ, DENS, & ! (in)
               GSQRT(:,:,:,I_XYZ), J33G, & ! (in)
               num_diff(:,:,:,I_MOMZ,ZDIR), & ! (in)
               CDZ, FDZ, dtrk, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          call ATMOS_DYN_FVM_fluxX_XYW_ud1( qflx_lo(:,:,:,XDIR), & ! (out)
               MOMX, MOMZ, DENS, & ! (in)
               GSQRT(:,:,:,I_UYZ), MAPF(:,:,:,I_UY), & ! (in)
               num_diff(:,:,:,I_MOMZ,XDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          call ATMOS_DYN_FVM_fluxY_XYW_ud1( qflx_lo(:,:,:,YDIR), & ! (out)
               MOMY, MOMZ, DENS, & ! (in)
               GSQRT(:,:,:,I_XVZ), MAPF(:,:,:,I_XV), & ! (in)
               num_diff(:,:,:,I_MOMZ,YDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

       endif
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then

       call COMM_vars8( DENS_RK, 1 )
       call COMM_wait ( DENS_RK, 1, .false. )

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qflx_hi(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) ) &
                              + qflx_J13(k,i,j) + qflx_J23(k,i,j)
       enddo
       enddo
       enddo

       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE-1
          DENS0_uvw(k,i,j) = 0.5_RP * ( DENS0(k,i,j) + DENS0(k+1,i,j) )
          DENS_uvw(k,i,j) = 0.5_RP * ( DENS_RK(k,i,j) + DENS_RK(k+1,i,j) )
       enddo
       enddo
       enddo

       do j = JS-1, JE+1
       do i = IS-1, IE+1
          DENS_uvw(KE,i,j) = DENS_uvw(KE-1,i,j)
          DENS0_uvw(KE,i,j) = DENS0_uvw(KE-1,i,j)
       enddo
       enddo

       call COMM_wait ( VELZ(:,:,:), 4 )

       call ATMOS_DYN_fct( qflx_anti,                 & ! (out)
                           VELZ, DENS0_uvw, DENS_uvw, & ! (in)
                           qflx_hi, qflx_lo,          & ! (in)
                           mflx_hi,                   & ! (in)
                           RFDZ, RCDX, RCDY,          & ! (in)
                           GSQRT(:,:,:,I_XYW),        & ! (in)
                           MAPF(:,:,:,I_XY), dtrk,    & ! (in)
                           FLAG_FCT_ALONG_STREAM      ) ! (in)

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
                            + dtrk * (   ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RFDZ(k) &
                                     + ( ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                       + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                                     * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) ) &
                                     / GSQRT(k,i,j,I_XYW)
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
       call ATMOS_DYN_FVM_fluxZ_UYZ( qflx_hi(:,:,:,ZDIR), & ! (out)
            MOMZ, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UYW), J33G, & ! (in)
            num_diff(:,:,:,I_MOMX,ZDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ13_UYZ( qflx_J13, & ! (out)
            MOMX, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UYZ), J13G(:,:,:,I_UYW), MAPF(:,:,:,I_UY), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ23_UYZ( qflx_J23, & ! (out)
            MOMY, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UYZ), J23G(:,:,:,I_UYW), MAPF(:,:,:,I_UY), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, y, z)
       ! note that x-index is added by -1
       call ATMOS_DYN_FVM_fluxX_UYZ( qflx_hi(:,:,:,XDIR), & ! (out)
            MOMX, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            num_diff(:,:,:,I_MOMX,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, v, z)
       call ATMOS_DYN_FVM_fluxY_UYZ( qflx_hi(:,:,:,YDIR), & ! (out)
            MOMY, MOMX, DENS, & ! (in)
            GSQRT(:,:,:,I_UVZ), MAPF(:,:,1,I_UV), & ! (in)
            num_diff(:,:,:,I_MOMX,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! pressure gradient force at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          pgf(k,i,j) = ( ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) & ! [x,y,z]
                         - GSQRT(k,i  ,j,I_XYZ) * DPRES(k,i  ,j) & ! [x,y,z]
                         ) * RFDX(i) &
                       + ( J13G(k  ,i,j,I_UYW) &
                         * 0.5_RP * ( F2H(k,1,I_UYZ) * ( DPRES(k+1,i+1,j)+DPRES(k+1,i,j) ) &
                                    + F2H(k,2,I_UYZ) * ( DPRES(k  ,i+1,j)+DPRES(k  ,i,j) ) ) & ! [x,y,z->u,y,w]
                         - J13G(k-1,i,j,I_UYW) &
                         * 0.5_RP * ( F2H(k,1,I_UYZ) * ( DPRES(k  ,i+1,j)+DPRES(k  ,i,j) ) &
                                    + F2H(k,2,I_UYZ) * ( DPRES(k-1,i+1,j)+DPRES(k-1,i,j) ) ) & ! [x,y,z->u,y,w]
                         ) * RCDZ(k) ) &
                     * MAPF(i,j,1,I_UY)
       enddo
       enddo
       enddo

       ! coriolis force at (u, y, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMY(k,i  ,j  ) )
          call CHECK( __LINE__, MOMY(k,i+1,j  ) )
          call CHECK( __LINE__, MOMY(k,i  ,j-1) )
          call CHECK( __LINE__, MOMY(k,i+1,j-1) )
#endif
          cor(k,i,j) = 0.125_RP * ( CORIOLI(1,i+1,j  )+CORIOLI(1,i,j  ) ) & ! [x,y,z->u,y,z]
                                * ( MOMY   (k,i+1,j  )+MOMY   (k,i,j  ) &
                                  + MOMY   (k,i+1,j-1)+MOMY   (k,i,j-1) ) &  ! [x,v,z->u,y,z]
                      - 0.25_RP * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) &
                      * ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) &
                      * ( ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) * 0.25_RP &
                        * ( 1.0_RP/MAPF(i+1,j,2,I_XY) - 1.0_RP/MAPF(i,j,2,I_XY) ) * RCDX(i) &
                        - MOMX(k,i,j) &
                        * ( 1.0_RP/MAPF(i,j,1,I_UV) - 1.0_RP/MAPF(i,j-1,1,I_UV) ) * RFDY(j) ) &
                      * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) ) ! metric term
       enddo
       enddo
       enddo

       !-----< update momentum (x) >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, min(IIE, IEH)
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DDIV(k,i+1,j) )
          call CHECK( __LINE__, DDIV(k,i  ,j) )
          call CHECK( __LINE__, MOMX0(k,i,j) )
#endif
          ! advection
          advcv = - ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k)
          advch = - ( ( qflx_J13(k,i,j) - qflx_J13(k-1,i,j)             ) * RCDZ(k) &
                    + ( qflx_J23(k,i,j) - qflx_J23(k-1,i,j)             ) * RCDZ(k) &
                    + ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                    + ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                  * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY)
          div = divdmp_coef / dtrk * FDX(i) * ( DDIV(k,i+1,j)-DDIV(k,i,j) ) ! divergence damping
          MOMX_RK(k,i,j) = MOMX0(k,i,j) &
                         + dtrk * ( ( advcv + advch        & ! advection
                                    - pgf(k,i,j)           & ! pressure gradient force
                                    ) / GSQRT(k,i,j,I_UYZ) &
                                    + cor(k,i,j)           & ! coriolis force
                                    + div                  & ! divergence damping
                                    + MOMX_t(k,i,j)        ) ! physics tendency
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMX) = advcv / GSQRT(k,i,j,I_UYZ)
             advch_t(k,i,j,I_MOMX) = advch / GSQRT(k,i,j,I_UYZ)
             pg_t(k,i,j,2) = - pgf(k,i,j) / GSQRT(k,i,j,I_UYZ)
             cf_t(k,i,j,1) = cor(k,i,j)
             ddiv_t(k,i,j,2) = div
          endif
#endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

#ifndef NO_FCT_DYN
       if ( FLAG_FCT_MOMENTUM ) then

          call ATMOS_DYN_FVM_fluxZ_UYZ_ud1( qflx_lo(:,:,:,ZDIR), & ! (out)
               MOMZ, MOMX, DENS, & ! (in)
               GSQRT(:,:,:,I_UYW), J33G, & ! (in)
               num_diff(:,:,:,I_MOMX,ZDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          ! note that x-index is added by -1
          call ATMOS_DYN_FVM_fluxX_UYZ_ud1( qflx_lo(:,:,:,XDIR), & ! (out)
               MOMX, MOMX, DENS, & ! (in)
               GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_UY), & ! (in)
               num_diff(:,:,:,I_MOMX,XDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          call ATMOS_DYN_FVM_fluxY_UYZ_ud1( qflx_lo(:,:,:,YDIR), & ! (out)
               MOMY, MOMX, DENS, & ! (in)
               GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_XV), & ! (in)
               num_diff(:,:,:,I_MOMX,YDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

       endif
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qflx_hi(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) / ( MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) )&
                              + qflx_J13(k,i,j) + qflx_J23(k,i,j)
       enddo
       enddo
       enddo

       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          DENS0_uvw(k,i,j) = 0.5_RP * ( DENS0(k,i,j) + DENS0(k,i+1,j) )
          DENS_uvw(k,i,j)  = 0.5_RP * ( DENS_RK(k,i,j) + DENS_RK(k,i+1,j) )
       enddo
       enddo
       enddo

       call COMM_wait ( VELX(:,:,:), 5 )

       call ATMOS_DYN_fct( qflx_anti,                 & ! (out)
                           VELX, DENS0_uvw, DENS_uvw, & ! (in)
                           qflx_hi, qflx_lo,          & ! (in)
                           mflx_hi,                   & ! (in)
                           RCDZ, RFDX, RCDY,          & ! (in)
                           GSQRT(:,:,:,I_UYZ),        & ! (in)
                           MAPF(:,:,:,I_UY), dtrk,    & ! (in)
                           FLAG_FCT_ALONG_STREAM      ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(x)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, JJE
          do i = IIS, min(IIE,IEH)
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
                            + dtrk * ( ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                                       + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) ) &
                            * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) &
                            / GSQRT(k,i,j,I_UYZ)
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
       call ATMOS_DYN_FVM_fluxZ_XVZ( qflx_hi(:,:,:,ZDIR), & ! (out)
            MOMZ, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XVW), J33G, & ! (in)
            num_diff(:,:,:,I_MOMY,ZDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ13_XVZ( qflx_J13, & ! (out)
            MOMX, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XVZ), J13G(:,:,:,I_XVW), MAPF(:,:,:,I_XV), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       call ATMOS_DYN_FVM_fluxJ23_XVZ( qflx_J23, & ! (out)
            MOMY, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XVZ), J23G(:,:,:,I_XVW), MAPF(:,:,:,I_XV), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, v, z)
       call ATMOS_DYN_FVM_fluxX_XVZ( qflx_hi(:,:,:,XDIR), & ! (out)
            MOMX, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_UV), & ! (in)
            num_diff(:,:,:,I_MOMY,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, y, z)
       ! note that y-index is added by -1
       call ATMOS_DYN_FVM_fluxY_XVZ( qflx_hi(:,:,:,YDIR), & ! (out)
            MOMY, MOMY, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
            num_diff(:,:,:,I_MOMY,YDIR), & ! (in
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! pressure gradient force at (x, v, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          pgf(k,i,j) = ( ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) & ! [x,y,z]
                         - GSQRT(k,i,j  ,I_XYZ) * DPRES(k,i,j  ) & ! [x,y,z]
                         ) * RFDY(j) &
                       + ( J23G(k  ,i,j,I_XVW) &
                         * 0.5_RP * ( F2H(k  ,1,I_XVZ) * ( DPRES(k+1,i,j+1)+DPRES(k+1,i,j) ) &
                                    + F2H(k  ,2,I_XVZ) * ( DPRES(k  ,i,j+1)+DPRES(k  ,i,j) ) ) & ! [x,y,z->x,v,w]
                         - J23G(k-1,i,j,I_XVW) &
                         * 0.5_RP * ( F2H(k-1,1,I_XVZ) * ( DPRES(k  ,i,j+1)+DPRES(k  ,i,j) ) &
                                    + F2H(k-1,2,I_XVZ) * ( DPRES(k-1,i,j+1)+DPRES(k-1,i,j) ) ) & ! [x,y,z->x,v,w]
                         ) * RCDZ(k) ) &
                      * MAPF(i,j,2,I_XV)
       enddo
       enddo
       enddo

       ! coriolis force at (x, v, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, MOMX(k,i  ,j  ) )
          call CHECK( __LINE__, MOMX(k,i  ,j+1) )
          call CHECK( __LINE__, MOMX(k,i-1,j  ) )
          call CHECK( __LINE__, MOMX(k,i-1,j+1) )
#endif
          cor(k,i,j) = - 0.125_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) & ! [x,y,z->x,v,z]
                                  * ( MOMX   (k,i  ,j+1)+MOMX   (k,i  ,j) &
                                    + MOMX   (k,i-1,j+1)+MOMX   (k,i-1,j) ) & ! [u,y,z->x,v,z]
                     + 0.25_RP * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) &
                     * ( MOMX(k,i,j) + MOMX(k,i-1,j) + MOMX(k,i,j+1) + MOMX(k,i-1,j+1) )&
                     * ( MOMY(k,i,j) &
                       * ( 1.0_RP/MAPF(i,j,2,I_UV) - 1.0_RP/MAPF(i-1,j,2,I_UV) ) * RCDX(i) &
                       - 0.25_RP * ( MOMX(k,i,j)+MOMX(k,i-1,j)+MOMX(k,i,j+1)+MOMX(k,i-1,j+1) ) &
                       * ( 1.0_RP/MAPF(i,j+1,1,I_XY) - 1.0_RP/MAPF(i,j,1,I_XY) ) * RFDY(j) ) &
                     * 2.0_RP / ( DENS(k,i,j) + DENS(k,i,j+1) ) ! metoric term
       enddo
       enddo
       enddo

       !-----< update momentum (y) >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, min(JJE, JEH)
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DDIV(k,i,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j  ) )
          call CHECK( __LINE__, MOMY_t(k,i,j) )
          call CHECK( __LINE__, MOMY0(k,i,j) )
#endif

          advcv = - ( qflx_hi(k,i,j,ZDIR) - qflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k)
          advch = - ( qflx_J13(k,i,j) - qflx_J13(k-1,i,j)             ) * RCDZ(k) &
                  - ( qflx_J23(k,i,j) - qflx_J23(k-1,i,j)             ) * RCDZ(k) &
                  - ( qflx_hi(k,i,j,XDIR) - qflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                  - ( qflx_hi(k,i,j,YDIR) - qflx_hi(k  ,i  ,j-1,YDIR) ) * RFDY(j) &
                * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV)
          div = divdmp_coef / dtrk * FDY(j) * ( DDIV(k,i,j+1)-DDIV(k,i,j) )
          MOMY_RK(k,i,j) = MOMY0(k,i,j) &
                         + dtrk * ( ( advcv + advch        & ! advection
                                    - pgf(k,i,j)           & ! pressure gradient force
                                    ) / GSQRT(k,i,j,I_XVZ) &
                                  + cor(k,i,j)             & ! coriolis force
                                  + div                    & ! divergence damping
                                  + MOMY_t(k,i,j)          ) ! physics tendency
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMY) = advcv / GSQRT(k,i,j,I_UYZ)
             advch_t(k,i,j,I_MOMY) = advch / GSQRT(k,i,j,I_UYZ)
             pg_t(k,i,j,3) = - pgf(k,i,j) / GSQRT(k,i,j,I_UYZ)
             cf_t(k,i,j,2) = cor(k,i,j)
             ddiv_t(k,i,j,3) = div
          endif
#endif
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
          call ATMOS_DYN_FVM_fluxZ_XVZ_ud1( qflx_lo(:,:,:,ZDIR), & ! (out)
               MOMZ, MOMY, DENS, & ! (in)
               GSQRT(:,:,:,I_XVZ), J33G, & ! (in)
               num_diff(:,:,:,I_MOMY,ZDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          ! at (u, v, layer)
          call ATMOS_DYN_FVM_fluxX_XVZ_ud1( qflx_lo(:,:,:,XDIR), & ! (out)
               MOMX, MOMY, DENS, & ! (in)
               GSQRT(:,:,:,I_UVZ), MAPF(:,:,:,I_XY), & ! (in)
               num_diff(:,:,:,I_MOMY,XDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          ! at (x, y, layer)
          ! note that y-index is added by -1
          call ATMOS_DYN_FVM_fluxY_XVZ_ud1( qflx_lo(:,:,:,YDIR), & ! (out)
               MOMY, MOMY, DENS, & ! (in)
               GSQRT(:,:,:,I_XYZ), MAPF(:,:,:,I_XY), & ! (in)
               num_diff(:,:,:,I_MOMY,YDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

       endif
    enddo
    enddo

    if ( FLAG_FCT_MOMENTUM ) then

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          qflx_hi(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) / ( MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) ) &
                              + qflx_J13(k,i,j) + qflx_J23(k,i,j)
       enddo
       enddo
       enddo

       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          DENS0_uvw(k,i,j) = 0.5_RP * ( DENS0(k,i,j) + DENS0(k,i,j+1) )
          DENS_uvw(k,i,j) = 0.5_RP * ( DENS_RK(k,i,j) + DENS_RK(k,i,j+1) )
       enddo
       enddo
       enddo

       call COMM_wait ( VELY(:,:,:), 6 )

       call ATMOS_DYN_fct( qflx_anti,                 & ! (out)
                           VELY, DENS0_uvw, DENS_uvw, & ! (in)
                           qflx_hi, qflx_lo,          & ! (in)
                           mflx_hi,                   & ! (in)
                           RCDZ, RCDX, RFDY,          & ! (in)
                           GSQRT(:,:,:,I_XVZ),        & ! (in)
                           MAPF(:,:,:,I_XV), dtrk,    & ! (in)
                           FLAG_FCT_ALONG_STREAM      ) ! (in)

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          !--- update momentum(y)
          !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
          do j = JJS, min(JJE, JEH)
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
                            + dtrk * ( ( ( qflx_anti(k,i,j,ZDIR) - qflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                       + ( qflx_anti(k,i,j,XDIR) - qflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                       + ( qflx_anti(k,i,j,YDIR) - qflx_anti(k  ,i  ,j-1,YDIR) ) * RFDY(j) ) ) &
                            * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) &
                            / GSQRT(k,i,j,I_XVZ)
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
       call ATMOS_DYN_FVM_fluxZ_XYZ( tflx_hi(:,:,:,ZDIR), & ! (out)
            mflx_hi(:,:,:,ZDIR), POTT, GSQRT(:,:,:,I_XYW), & ! (in)
            num_diff(:,:,:,I_RHOT,ZDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, y, z)
       call ATMOS_DYN_FVM_fluxX_XYZ( tflx_hi(:,:,:,XDIR), & ! (out)
            mflx_hi(:,:,:,XDIR), POTT, GSQRT(:,:,:,I_UYZ), & ! (in)
            num_diff(:,:,:,I_RHOT,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, v, z)
       call ATMOS_DYN_FVM_fluxY_XYZ( tflx_hi(:,:,:,YDIR), & ! (out)
            mflx_hi(:,:,:,YDIR), POTT, GSQRT(:,:,:,I_XVZ), & ! (in)
            num_diff(:,:,:,I_RHOT,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       !-----< update rho*theta >-----

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, tflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, tflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, RHOT_t(k,i,j) )
          call CHECK( __LINE__, RHOT0(k,i,j) )
#endif
          advcv = - ( tflx_hi(k,i,j,ZDIR) - tflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k)
          advch = - ( tflx_hi(k,i,j,XDIR) - tflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                  - ( tflx_hi(k,i,j,YDIR) - tflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j)
          RHOT_RK(k,i,j) = RHOT0(k,i,j) &
                         + dtrk * ( ( advcv + advch ) * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ) &
                                  + RHOT_t(k,i,j) )
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_RHOT) = advcv * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY)/ GSQRT(k,i,j,I_XYZ)
             advch_t(k,i,j,I_RHOT) = advch * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY)/ GSQRT(k,i,j,I_XYZ)
          endif
#endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

#ifndef NO_FCT_DYN

    if ( FLAG_FCT_T ) then

       call COMM_vars8( mflx_hi(:,:,:,ZDIR), 1 )
       call COMM_vars8( mflx_hi(:,:,:,XDIR), 2 )
       call COMM_vars8( mflx_hi(:,:,:,YDIR), 3 )
       call COMM_wait ( mflx_hi(:,:,:,ZDIR), 1, .false. )
       call COMM_wait ( mflx_hi(:,:,:,XDIR), 2, .false. )
       call COMM_wait ( mflx_hi(:,:,:,YDIR), 3, .false. )

       if ( .NOT. FLAG_FCT_MOMENTUM ) then
          call COMM_vars8( DENS_RK, 1 )
          call COMM_wait ( DENS_RK, 1, .false. )
       endif

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          ! monotonic flux
          ! at (x, y, interface)
          call ATMOS_DYN_FVM_fluxZ_XYZ_ud1( tflx_lo(:,:,:,ZDIR), & ! (out)
               mflx_hi(:,:,:,ZDIR), POTT, GSQRT(:,:,:,I_XYZ), & ! (in)
               num_diff(:,:,:,I_RHOT,ZDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          ! at (u, y, layer)
          call ATMOS_DYN_FVM_fluxX_XYZ_ud1( tflx_lo(:,:,:,XDIR), & ! (out)
               mflx_hi(:,:,:,XDIR), POTT, GSQRT(:,:,:,I_UYZ), & ! (in)
               num_diff(:,:,:,I_RHOT,XDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

          ! at (x, v, layer)
          call ATMOS_DYN_FVM_fluxY_XYZ_ud1( tflx_lo(:,:,:,YDIR), & ! (out)
               mflx_hi(:,:,:,YDIR), POTT, GSQRT(:,:,:,I_XVZ), & ! (in)
               num_diff(:,:,:,I_RHOT,YDIR), & ! (in)
               CDZ, & ! (in)
               IIS-1, IIE+1, JJS-1, JJE+1 ) ! (in)

       enddo
       enddo

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          POTT(k,i,j) = RHOT0(k,i,j) / DENS0(k,i,j)
       enddo
       enddo
       enddo

       call ATMOS_DYN_fct( tflx_anti,               & ! (out)
                           POTT, DENS0, DENS_RK,    & ! (out)
                           tflx_hi, tflx_lo,        & ! (in)
                           mflx_hi,                 & ! (in)
                           RCDZ, RCDX, RCDY,        & ! (in)
                           GSQRT(:,:,:,I_XYZ),      & ! (in)
                           MAPF(:,:,:,I_XY), dtrk,  & ! (in)
                           FLAG_FCT_ALONG_STREAM    ) ! (in)

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
             call CHECK( __LINE__, tflx_anti(k  ,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, tflx_anti(k-1,i  ,j  ,ZDIR) )
             call CHECK( __LINE__, tflx_anti(k  ,i  ,j  ,XDIR) )
             call CHECK( __LINE__, tflx_anti(k  ,i-1,j  ,XDIR) )
             call CHECK( __LINE__, tflx_anti(k  ,i  ,j  ,YDIR) )
             call CHECK( __LINE__, tflx_anti(k  ,i  ,j-1,YDIR) )
#endif
             RHOT_RK(k,i,j) = RHOT_RK(k,i,j) &
                            + dtrk * ( ( tflx_anti(k,i,j,ZDIR) - tflx_anti(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                                     + ( tflx_anti(k,i,j,XDIR) - tflx_anti(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                                     + ( tflx_anti(k,i,j,YDIR) - tflx_anti(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                           * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                           / GSQRT(k,i,j,I_XYZ)
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

#ifdef HIST_TEND
    if ( lhist ) then
       call HIST_in(advcv_t(:,:,:,I_DENS), 'DENS_t_advcv', 'tendency of density    (vert. advection)',    'kg/m3/s'   )
       call HIST_in(advcv_t(:,:,:,I_MOMZ), 'MOMZ_t_advcv', 'tendency of momentum z (vert. advection)',    'kg/m2/s2', zdim='half')
       call HIST_in(advcv_t(:,:,:,I_MOMX), 'MOMX_t_advcv', 'tendency of momentum x (vert. advection)',    'kg/m2/s2', xdim='half')
       call HIST_in(advcv_t(:,:,:,I_MOMY), 'MOMY_t_advcv', 'tendency of momentum y (vert. advection)',    'kg/m2/s2', ydim='half')
       call HIST_in(advcv_t(:,:,:,I_RHOT), 'RHOT_t_advcv', 'tendency of rho*theta  (vert. advection)',    'K kg/m3/s' )

       call HIST_in(advch_t(:,:,:,I_DENS), 'DENS_t_advch', 'tendency of density    (horiz. advection)',   'kg/m3/s'   )
       call HIST_in(advch_t(:,:,:,I_MOMZ), 'MOMZ_t_advch', 'tendency of momentum z (horiz. advection)',   'kg/m2/s2', zdim='half')
       call HIST_in(advch_t(:,:,:,I_MOMX), 'MOMX_t_advch', 'tendency of momentum x (horiz. advection)',   'kg/m2/s2', xdim='half')
       call HIST_in(advch_t(:,:,:,I_MOMY), 'MOMY_t_advch', 'tendency of momentum y (horiz. advection)',   'kg/m2/s2', ydim='half')
       call HIST_in(advch_t(:,:,:,I_RHOT), 'RHOT_t_advch', 'tendency of rho*theta  (horiz. advection)',   'K kg/m3/s' )

       call HIST_in(pg_t   (:,:,:,1),      'MOMZ_t_pg',    'tendency of momentum z (pressure gradient)',  'kg/m2/s2', zdim='half')
       call HIST_in(pg_t   (:,:,:,2),      'MOMX_t_pg',    'tendency of momentum x (pressure gradient)',  'kg/m2/s2', xdim='half')
       call HIST_in(pg_t   (:,:,:,3),      'MOMY_t_pg',    'tendency of momentum y (pressure gradient)',  'kg/m2/s2', ydim='half')

       call HIST_in(ddiv_t (:,:,:,1),      'MOMZ_t_ddiv',  'tendency of momentum z (divergence damping)', 'kg/m2/s2', zdim='half')
       call HIST_in(ddiv_t (:,:,:,2),      'MOMX_t_ddiv',  'tendency of momentum x (divergence damping)', 'kg/m2/s2', xdim='half')
       call HIST_in(ddiv_t (:,:,:,3),      'MOMY_t_ddiv',  'tendency of momentum y (divergence damping)', 'kg/m2/s2', ydim='half')

       call HIST_in(cf_t   (:,:,:,1),      'MOMX_t_cf',    'tendency of momentum x (coliolis force)',     'kg/m2/s2', xdim='half')
       call HIST_in(cf_t   (:,:,:,2),      'MOMY_t_cf',    'tendency of momentum y (coliolis force)',     'kg/m2/s2', ydim='half')
    endif
#endif

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_heve

end module scale_atmos_dyn_tstep_short_fvm_heve
