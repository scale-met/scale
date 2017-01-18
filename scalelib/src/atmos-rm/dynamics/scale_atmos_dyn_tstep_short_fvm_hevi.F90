!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          HEVI FVM scheme for Atmospheric dynamical process
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

#ifdef PROFILE_FAPP
#define PROFILE_START(name) call fapp_start(name, 1, 1)
#define PROFILE_STOP(name)  call fapp_stop (name, 1, 1)
#elif defined(PROFILE_FINEPA)
#define PROFILE_START(name) call start_collection(name)
#define PROFILE_STOP(name)  call stop_collection (name)
#else
#define PROFILE_START(name)
#define PROFILE_STOP(name)
#endif

module scale_atmos_dyn_tstep_short_fvm_hevi
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
  public :: ATMOS_DYN_Tstep_short_fvm_hevi_regist
  public :: ATMOS_DYN_Tstep_short_fvm_hevi_setup
  public :: ATMOS_DYN_Tstep_short_fvm_hevi

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
# define F2H(k,p,idx) 0.5_RP
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: NB = 1
  integer,  private, parameter :: VA_FVM_HEVI = 0
  integer                      :: IFS_OFF
  integer                      :: JFS_OFF

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tstep_short_fvm_hevi_regist( &
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

    if ( ATMOS_DYN_TYPE /= 'FVM-HEVI' .AND. ATMOS_DYN_TYPE /= 'HEVI' ) then
       write(*,*) 'xxx ATMOS_DYN_TYPE is not FVM-HEVI. Check!'
       call PRC_MPIstop
    endif

    VA_out      = VA_FVM_HEVI
    VAR_NAME(:) = ""
    VAR_DESC(:) = ""
    VAR_UNIT(:) = ""

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Register additional prognostic variables (HEVI)'
    if ( VA_out < 1 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** => nothing.'
    endif

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hevi_regist

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_short_fvm_hevi_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HEVI Setup'
#ifdef HEVI_BICGSTAB
    if( IO_L ) write(IO_FID_LOG,*) '*** USING Bi-CGSTAB'
#elif defined(HEVI_LAPACK)
    if( IO_L ) write(IO_FID_LOG,*) '*** USING LAPACK'
#else
    if( IO_L ) write(IO_FID_LOG,*) '*** USING DIRECT'
#endif

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hevi_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_Tstep_short_fvm_hevi( &
       DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
       PROG_RK,                                     &
       mflx_hi, tflx_hi,                            &
       DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
       DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
       PROG0, PROG,                                 &
       DPRES0, RT2P, CORIOLI,                       &
       num_diff, wdamp_coef, divdmp_coef, DDIV,     &
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
#ifdef DRY
       Rdry   => CONST_Rdry,  &
       CVdry  => CONST_CVdry, &
       CPdry  => CONST_CPdry, &
#endif
       GRAV   => CONST_GRAV,  &
       P00    => CONST_PRE00
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_fct
    use scale_atmos_dyn_fvm_flux, only: &
       ATMOS_DYN_FVM_flux_valueW_Z, &
       ATMOS_DYN_FVM_fluxZ_XYZ,     &
       ATMOS_DYN_FVM_fluxX_XYZ,     &
       ATMOS_DYN_FVM_fluxY_XYZ,     &
       ATMOS_DYN_FVM_fluxZ_XYW,     &
       ATMOS_DYN_FVM_fluxJ13_XYW,   &
       ATMOS_DYN_FVM_fluxJ23_XYW,   &
       ATMOS_DYN_FVM_fluxX_XYW,     &
       ATMOS_DYN_FVM_fluxY_XYW,     &
       ATMOS_DYN_FVM_fluxZ_UYZ,     &
       ATMOS_DYN_FVM_fluxJ13_UYZ,   &
       ATMOS_DYN_FVM_fluxJ23_UYZ,   &
       ATMOS_DYN_FVM_fluxX_UYZ,     &
       ATMOS_DYN_FVM_fluxY_UYZ,     &
       ATMOS_DYN_FVM_fluxZ_XVZ,     &
       ATMOS_DYN_FVM_fluxJ13_XVZ,   &
       ATMOS_DYN_FVM_fluxJ23_XVZ,   &
       ATMOS_DYN_FVM_fluxX_XVZ,     &
       ATMOS_DYN_FVM_fluxY_XVZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       I_XY , &
       I_UY , &
       I_XV , &
       I_UV
#ifdef HIST_TEND
    use scale_history, only: &
       HIST_in
#endif
    implicit none

    real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(out) :: PROG_RK(KA,IA,JA,VA)  !

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! rho * vel(x,y,z)
    real(RP), intent(out)   :: tflx_hi(KA,IA,JA,3) ! rho * theta * vel(x,y,z)

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

    real(RP), intent(in)  :: PROG0(KA,IA,JA,VA)
    real(RP), intent(in)  :: PROG (KA,IA,JA,VA)

    real(RP), intent(in) :: DPRES0(KA,IA,JA)
    real(RP), intent(in) :: RT2P(KA,IA,JA)
    real(RP), intent(in) :: CORIOLI(1,IA,JA)
    real(RP), intent(in) :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in) :: wdamp_coef(KA)
    real(RP), intent(in) :: divdmp_coef
    real(RP), intent(in) :: DDIV(KA,IA,JA)

    logical,  intent(in) :: FLAG_FCT_MOMENTUM
    logical,  intent(in) :: FLAG_FCT_T
    logical,  intent(in) :: FLAG_FCT_ALONG_STREAM

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
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
    real(RP), intent(in)  :: REF_dens(KA,IA,JA)   !< reference density
    real(RP), intent(in)  :: REF_rhot(KA,IA,JA)

    logical,  intent(in)  :: BND_W
    logical,  intent(in)  :: BND_E
    logical,  intent(in)  :: BND_S
    logical,  intent(in)  :: BND_N

    real(RP), intent(in)  :: dtrk
    real(RP), intent(in)  :: dt


    ! diagnostic variables (work space)
    real(RP) :: POTT(KA,IA,JA) ! potential temperature [K]
    real(RP) :: DPRES(KA,IA,JA) ! pressure deviation from reference pressure

    real(RP) :: qflx_hi (KA,IA,JA,3)
    real(RP) :: qflx_J13(KA,IA,JA)
    real(RP) :: qflx_J23(KA,IA,JA)

    real(RP) :: advch ! horizontal advection
    real(RP) :: advcv ! vertical advection
    real(RP) :: wdmp  ! Raileight damping
    real(RP) :: div   ! divergence damping
    real(RP) :: pg    ! pressure gradient force
    real(RP) :: cf    ! colioris force
#ifdef HIST_TEND
    real(RP) :: advch_t(KA,IA,JA,5)
    real(RP) :: advcv_t(KA,IA,JA,5)
    real(RP) :: wdmp_t(KA,IA,JA)
    real(RP) :: ddiv_t(KA,IA,JA,3)
    real(RP) :: pg_t(KA,IA,JA,3)
    real(RP) :: cf_t(KA,IA,JA,2)
    logical  :: lhist
#endif

    ! for implicit solver
    real(RP) :: A(KA,IA,JA)
    real(RP) :: B
    real(RP) :: Sr(KA,IA,JA)
    real(RP) :: Sw(KA,IA,JA)
    real(RP) :: St(KA,IA,JA)
    real(RP) :: PT(KA,IA,JA)
    real(RP) :: C(KMAX-1,IA,JA)

    real(RP) :: F1(KA,IA,JA)
    real(RP) :: F2(KA,IA,JA)
    real(RP) :: F3(KA,IA,JA)

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j
    integer :: iss, iee

#ifdef DEBUG
    POTT(:,:,:) = UNDEF
    DPRES(:,:,:) = UNDEF

    PT(:,:,:) = UNDEF

    qflx_hi (:,:,:,:) = UNDEF
    qflx_J13(:,:,:)   = UNDEF
    qflx_J23(:,:,:)   = UNDEF
#endif

#ifdef HIST_TEND
    lhist = dt .eq. dtrk
#endif

#ifdef PROFILE_FIPP
    call fipp_start()
#endif

    IFS_OFF = 1
    JFS_OFF = 1
    if ( BND_W ) IFS_OFF = 0
    if ( BND_S ) JFS_OFF = 0


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       PROFILE_START("hevi_pres")
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
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
       PROFILE_STOP("hevi_pres")

       !##### continuity equation #####

       PROFILE_START("hevi_mflx_z")
       ! at (x, y, w)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS-1, IIE
          mflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
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
             mflx_hi(k,i,j,ZDIR) = J13G(k,i,j,I_XYW) * 0.25_RP / MAPF(i,j,2,I_XY) &
                                 * ( MOMX(k+1,i,j)+MOMX(k+1,i-1,j) &
                                   + MOMX(k  ,i,j)+MOMX(k  ,i-1,j) ) &
                                 + J23G(k,i,j,I_XYW) * 0.25_RP / MAPF(i,j,1,I_XY) &
                                 * ( MOMY(k+1,i,j)+MOMY(k+1,i,j-1) &
                                   + MOMY(k  ,i,j)+MOMY(k  ,i,j-1) ) &
                                 + GSQRT(k,i,j,I_XYW) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) ) * num_diff(k,i,j,I_DENS,ZDIR)
          enddo
          mflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_mflx_z")

       PROFILE_START("hevi_mflx_x")
       iss = max(IIS-1,IS-IFS_OFF)
       iee = min(IIE,IEH)
       ! at (u, y, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = iss, iee
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
          call CHECK( __LINE__, MOMX(k,i,j) )
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
       PROFILE_STOP("hevi_mflx_x")

       ! at (x, v, z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = max(JJS-1,JS-JFS_OFF), min(JJE,JEH)
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
          call CHECK( __LINE__, MOMY(k,i,j) )
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

       !--- update density
       PROFILE_START("hevi_sr")
       !$omp parallel do private(i,j,k,advch) OMP_SCHEDULE_ collapse(2)
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
          advch = - ( ( mflx_hi(k,i,j,ZDIR)-mflx_hi(k-1,i  ,j,  ZDIR) ) * RCDZ(k) &
                    + ( mflx_hi(k,i,j,XDIR)-mflx_hi(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                    + ( mflx_hi(k,i,j,YDIR)-mflx_hi(k  ,i,  j-1,YDIR) ) * RCDY(j) ) &
                * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
          Sr(k,i,j) =  advch + DENS_t(k,i,j)
#ifdef HIST_TEND
          if ( lhist ) advch_t(k,i,j,I_DENS) = advch
#endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_sr")

       !##### momentum equation (z) #####

       ! at (x, y, z)
       ! not that z-index is added by -1
       PROFILE_START("hevi_momz_qflxhi_z")
       call ATMOS_DYN_FVM_fluxZ_XYW( qflx_hi(:,:,:,ZDIR), & ! (out)
            MOMZ, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XYZ), J33G, & ! (in)
            num_diff(:,:,:,I_MOMZ,ZDIR), & ! (in)
            CDZ, FDZ, dtrk, &
            IIS, IIE, JJS, JJE ) ! (in)
       PROFILE_STOP("hevi_momz_qflxhi_z")

       PROFILE_START("hevi_momz_qflxj")
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
       PROFILE_STOP("hevi_momz_qflxj")

       ! at (u, y, w)
       PROFILE_START("hevi_momz_qflxhi_x")
       call ATMOS_DYN_FVM_fluxX_XYW( qflx_hi(:,:,:,XDIR), & ! (out)
            MOMX, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_UYW), MAPF(:,:,:,I_UY), & ! (in)
            num_diff(:,:,:,I_MOMZ,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       PROFILE_STOP("hevi_momz_qflxhi_x")

       ! at (x, v, w)
       PROFILE_START("hevi_momz_qflxhi_y")
       call ATMOS_DYN_FVM_fluxY_XYW( qflx_hi(:,:,:,YDIR), & ! (out)
            MOMY, MOMZ, DENS, & ! (in)
            GSQRT(:,:,:,I_XVW), MAPF(:,:,:,I_XV), & ! (in)
            num_diff(:,:,:,I_MOMZ,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)
       PROFILE_STOP("hevi_momz_qflxhi_y")

       !--- update momentum(z)
       PROFILE_START("hevi_sw")
       !$omp parallel do private(i,j,k,advcv,advch,cf,wdmp,div) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi (k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_J13(k  ,i  ,j)        )
          call CHECK( __LINE__, qflx_J13(k-1,i  ,j)        )
          call CHECK( __LINE__, qflx_J23(k  ,i  ,j)        )
          call CHECK( __LINE__, qflx_J23(k-1,i  ,j)        )
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi (k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DDIV(k  ,i,j) )
          call CHECK( __LINE__, DDIV(k+1,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, MOMZ_t(k,i,j) )
#endif
          advcv = -   ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) ) * RFDZ(k)
          advch = - ( ( qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j  ) &
                      + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j  )    ) * RFDZ(k) &
                    + ( qflx_hi (k,i,j,XDIR) - qflx_hi (k,i-1,j  ,XDIR) ) * RCDX(i) &
                    + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k,i  ,j-1,YDIR) ) * RCDY(j) ) &
                * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY)
          wdmp = - wdamp_coef(k) * MOMZ0(k,i,j)
          div = divdmp_coef / dtrk * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) ! divergence damping
          Sw(k,i,j) = ( advcv + advch ) / GSQRT(k,i,j,I_XYW) &
                    + wdmp + div + MOMZ_t(k,i,j)
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMZ) = advcv / GSQRT(k,i,j,I_XYW)
             advch_t(k,i,j,I_MOMZ) = advch / GSQRT(k,i,j,I_XYW)
             wdmp_t(k,i,j) = wdmp
             ddiv_t(k,i,j,1) = div
          endif
#endif
       enddo
       enddo
       enddo
       PROFILE_STOP("hevi_sw")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif


       !##### Thermodynamic Equation #####

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
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


       PROFILE_START("hevi_st")
       !$omp parallel do private(i,j,k,advch) OMP_SCHEDULE_ collapse(2)
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
#endif
          advch = - ( ( tflx_hi(k,i,j,ZDIR) - tflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                    + ( tflx_hi(k,i,j,XDIR) - tflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                    + ( tflx_hi(k,i,j,YDIR) - tflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                  * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / GSQRT(k,i,j,I_XYZ)
          St(k,i,j) = advch + RHOT_t(k,i,j)
#ifdef HIST_TEND
          if ( lhist ) then
             advch_t(k,i,j,I_RHOT) = advch
          endif
#endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       PROFILE_STOP("hevi_st")

       ! implicit solver

       PROFILE_START("hevi_solver")

!OCL INDEPENDENT
!OCL PREFETCH_SEQUENTIAL(SOFT)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE

          call ATMOS_DYN_FVM_flux_valueW_Z( PT(:,i,j), & ! (out)
               MOMZ(:,i,j), POTT(:,i,j), GSQRT(:,i,j,I_XYZ), & ! (in)
               CDZ )

          do k = KS, KE
             A(k,i,j) = dtrk**2 * J33G * RCDZ(k) * RT2P(k,i,j) * J33G / GSQRT(k,i,j,I_XYZ)
          enddo
          B = GRAV * dtrk**2 * J33G / ( CDZ(KS+1) + CDZ(KS) )
          F1(KS,i,j) =        - ( PT(KS+1,i,j) * RFDZ(KS) *   A(KS+1,i,j)             + B ) / GSQRT(KS,i,j,I_XYW)
          F2(KS,i,j) = 1.0_RP + ( PT(KS  ,i,j) * RFDZ(KS) * ( A(KS+1,i,j)+A(KS,i,j) )     ) / GSQRT(KS,i,j,I_XYW)
          do k = KS+1, KE-2
             B = GRAV * dtrk**2 * J33G / ( CDZ(k+1) + CDZ(k) )
             F1(k,i,j) =        - ( PT(k+1,i,j) * RFDZ(k) *   A(k+1,i,j)            + B ) / GSQRT(k,i,j,I_XYW)
             F2(k,i,j) = 1.0_RP + ( PT(k  ,i,j) * RFDZ(k) * ( A(k+1,i,j)+A(k,i,j) )     ) / GSQRT(k,i,j,I_XYW)
             F3(k,i,j) =        - ( PT(k-1,i,j) * RFDZ(k) *              A(k,i,j)   - B ) / GSQRT(k,i,j,I_XYW)
          enddo
          B = GRAV * dtrk**2 * J33G / ( CDZ(KE) + CDZ(KE-1) )
          F2(KE-1,i,j) = 1.0_RP + ( PT(KE-1,i,j) * RFDZ(KE-1) * ( A(KE,i,j)+A(KE-1,i,j) )    ) / GSQRT(KE-1,i,j,I_XYW)
          F3(KE-1,i,j) =        - ( PT(KE-2,i,j) * RFDZ(KE-1) *             A(KE-1,i,j)  - B ) / GSQRT(KE-1,i,j,I_XYW)
          do k = KS, KE-1
             pg = - ( DPRES(k+1,i,j) + RT2P(k+1,i,j)*dtrk*St(k+1,i,j) &
                    - DPRES(k  ,i,j) - RT2P(k  ,i,j)*dtrk*St(k  ,i,j) ) &
                    * RFDZ(k) * J33G / GSQRT(k,i,j,I_XYW) &
                  - GRAV &
                    * ( F2H(k,1,I_XYZ) * ( DENS(k+1,i,j) - REF_dens(k+1,i,j) + Sr(k+1,i,j) * dtrk ) &
                      + F2H(k,2,I_XYZ) * ( DENS(k  ,i,j) - REF_dens(k  ,i,j) + Sr(k  ,i,j) * dtrk ) )
             C(k-KS+1,i,j) = MOMZ(k,i,j) + dtrk * ( pg + Sw(k,i,j) )
#ifdef HIST_TEND
             if ( lhist ) pg_t(k,i,j,1) = pg
#endif
          enddo

#ifdef HEVI_BICGSTAB
          call solve_bicgstab( &
               C(:,i,j),         & ! (inout)
               F1(:,i,j), F2(:,i,j), F3(:,i,j) ) ! (in)
#elif defined(HEVI_LAPACK)
          call solve_lapack( &
               C(:,i,j),        & ! (inout)
#ifdef DEBUG
               i, j,      & ! (in)
#endif
               F1(:,i,j), F2(:,i,j), F3(:,i,j) ) ! (in)
#else
          call solve_direct( &
               C(:,i,j),        & ! (inout)
               F1(:,i,j), F2(:,i,j), F3(:,i,j) ) ! (in)
#endif

          do k = KS, KE-1
#ifdef DEBUG_HEVI2HEVE
          ! for debug (change to explicit integration)
             C(k-KS+1,i,j) = MOMZ(k,i,j)
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                                 + J33G * MOMZ(k,i,j) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
             MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                  + dtrk*( &
                  - J33G * ( DPRES(k+1,i,j)-DPRES(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,i_XYW) &
                  - GRAV * ( F2H(k,2,I_XYZ)*(DENS(k,i,j)-REF_dens(k,i,j))+F2H(k,1,I_XYZ)*(DENS(k+1,i,j)-REF_dens(k+1,i,j)) ) &
                  + Sw(k,i,j) )
#else
             ! z-momentum flux
             mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
                                 + J33G * C(k-KS+1,i,j) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) )
             ! z-momentum
             MOMZ_RK(k,i,j) = MOMZ0(k,i,j) &
                            + ( C(k-KS+1,i,j) - MOMZ(k,i,j) )
#endif
          enddo
          MOMZ_RK(KS-1,i,j) = 0.0_RP
          MOMZ_RK(KE  ,i,j) = 0.0_RP

          ! density and rho*theta
          advcv = - C(1,i,j)            * J33G * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) ! C(0) = 0
          DENS_RK(KS,i,j) = DENS0(KS,i,j) + dtrk * ( advcv + Sr(KS,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KS,i,j,I_DENS) = advcv
#endif
          advcv = - C(1,i,j)*PT(KS,i,j) * J33G * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ) ! C(0) = 0
          RHOT_RK(KS,i,j) = RHOT0(KS,i,j) + dtrk * ( advcv + St(KS,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KS,i,j,I_RHOT) = advcv
#endif
          do k = KS+1, KE-1
             advcv = - ( C(k-KS+1,i,j)           - C(k-KS,i,j) ) &
                   * J33G * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
             DENS_RK(k,i,j) = DENS0(k,i,j) + dtrk * ( advcv + Sr(k,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(k,i,j,I_DENS) = advcv
#endif
             advcv = - ( C(k-KS+1,i,j)*PT(k,i,j) - C(k-KS,i,j)*PT(k-1,i,j) ) &
                   * J33G * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
             RHOT_RK(k,i,j) = RHOT0(k,i,j) + dtrk * ( advcv + St(k,i,j) )
#ifdef HIST_TEND
             if ( lhist ) advcv_t(k,i,j,I_RHOT) = advcv
#endif
          enddo
          advcv = C(KE-KS,i,j)                * J33G * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) ! C(KE-KS+1) = 0
          DENS_RK(KE,i,j) = DENS0(KE,i,j) + dtrk * ( advcv + Sr(KE,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KE,i,j,I_DENS) = advcv
#endif
          advcv = C(KE-KS,i,j) * PT(KE-1,i,j) * J33G * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ) ! C(KE-KS+1) = 0
          RHOT_RK(KE,i,j) = RHOT0(KE,i,j) + dtrk * ( advcv + St(KE,i,j) )
#ifdef HIST_TEND
          if ( lhist ) advcv_t(KE,i,j,I_RHOT) = advcv
#endif

#ifdef DEBUG
          call check_equation( &
               C(:,i,j), &
               DENS(:,i,j), MOMZ(:,i,j), RHOT(:,i,j), DPRES(:,i,j), &
               REF_dens(:,i,j), &
               Sr(:,i,j), Sw(:,i,j), St(:,i,j), &
               J33G, GSQRT(:,i,j,:), &
               RT2P(:,i,j), &
               dtrk, i, j )
#endif

       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       PROFILE_STOP("hevi_solver")


       !##### momentum equation (x) #####

       PROFILE_START("hevi_momx")
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

       !--- update momentum(x)
       iee = min(IIE,IEH)
       !$omp parallel do private(i,j,k,advch,advcv,pg,cf,div) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, iee
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DPRES(k,i+1,j) )
          call CHECK( __LINE__, DPRES(k,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i  ,j) )
          call CHECK( __LINE__, CORIOLI(1,i+1,j) )
          call CHECK( __LINE__, MOMY(k,i  ,j  ) )
          call CHECK( __LINE__, MOMY(k,i+1,j  ) )
          call CHECK( __LINE__, MOMY(k,i  ,j-1) )
          call CHECK( __LINE__, MOMY(k,i+1,j-1) )
          call CHECK( __LINE__, DDIV(k,i+1,j) )
          call CHECK( __LINE__, DDIV(k,i  ,j) )
          call CHECK( __LINE__, MOMX0(k,i,j) )
#endif
          advcv = -   ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) ) * RCDZ(k)
          advch = - ( ( qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j  ) &
                      + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j  )      ) * RFDZ(k) &
                    + ( qflx_hi (k,i,j,XDIR) - qflx_hi (k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                    + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY)
          pg = ( ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) & ! [x,y,z]
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
          cf = 0.125_RP * ( CORIOLI(1,i+1,j  )+CORIOLI(1,i,j  ) ) & ! [x,y,z->u,y,z]
             * ( MOMY   (k,i+1,j  )+MOMY   (k,i,j  ) &
               + MOMY   (k,i+1,j-1)+MOMY   (k,i,j-1) ) &  ! [x,v,z->u,y,z]
             + 0.25_RP * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) &
             * ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) &
             * ( ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) * 0.25_RP &
                 * ( 1.0_RP/MAPF(i+1,j,2,I_XY) - 1.0_RP/MAPF(i,j,2,I_XY) ) * RCDX(i) &
                   - MOMX(k,i,j) &
                   * ( 1.0_RP/MAPF(i,j,1,I_UV) - 1.0_RP/MAPF(i,j-1,1,I_UV) ) * RFDY(j) ) &
             * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) ) ! metric term
          div = divdmp_coef / dtrk * ( DDIV(k,i+1,j)/MAPF(i+1,j,2,I_XY) - DDIV(k,i,j)/MAPF(i,j,1,I_XY) ) &
              * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) * FDX(i) ! divergence damping
          MOMX_RK(k,i,j) = MOMX0(k,i,j) &
               + dtrk * ( ( advcv + advch - pg ) / GSQRT(k,i,j,I_UYZ) + cf + div + MOMX_t(k,i,j) )
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMX) = advcv / GSQRT(k,i,j,I_UYZ)
             advch_t(k,i,j,I_MOMX) = advch / GSQRT(k,i,j,I_UYZ)
             pg_t(k,i,j,2) = - pg / GSQRT(k,i,j,I_UYZ)
             cf_t(k,i,j,1) = cf
             ddiv_t(k,i,j,2) = div
          endif
#endif
       enddo
       enddo
       enddo
       PROFILE_STOP("hevi_momx")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (y) #####
       PROFILE_START("hevi_momy")
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

       !--- update momentum(y)
       !$omp parallel do private(i,j,k,advch,advcv,pg,cf,div) OMP_SCHEDULE_ collapse(2)
       do j = JJS, min(JJE,JEH)
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi(k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DPRES(k,i,j  ) )
          call CHECK( __LINE__, DPRES(k,i,j+1) )
          call CHECK( __LINE__, CORIOLI(1,i,j  ) )
          call CHECK( __LINE__, CORIOLI(1,i,j+1) )
          call CHECK( __LINE__, MOMX(k,i  ,j  ) )
          call CHECK( __LINE__, MOMX(k,i  ,j+1) )
          call CHECK( __LINE__, MOMX(k,i-1,j  ) )
          call CHECK( __LINE__, MOMX(k,i-1,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j+1) )
          call CHECK( __LINE__, DDIV(k,i,j  ) )
          call CHECK( __LINE__, MOMY_t(k,i,j) )
          call CHECK( __LINE__, MOMY0(k,i,j) )
#endif
          advcv = -   ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) ) * RCDZ(k)
          advch = - ( ( qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j  ) &
                      + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j  )      ) * RCDZ(k) &
                    + ( qflx_hi (k,i,j,XDIR) - qflx_hi (k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                    + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k  ,i  ,j-1,YDIR) ) * RFDY(j) ) &
                  * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV)
          pg = ( ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) & ! [x,y,z]
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
          cf = - 0.125_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) & ! [x,y,z->x,v,z]
                          * ( MOMX   (k,i  ,j+1)+MOMX   (k,i  ,j) &
                            + MOMX   (k,i-1,j+1)+MOMX   (k,i-1,j) ) & ! [u,y,z->x,v,z]
               - 0.25_RP * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) &
               * ( MOMX(k,i,j) + MOMX(k,i-1,j) + MOMX(k,i,j+1) + MOMX(k,i-1,j+1) ) &
               * ( MOMY(k,i,j) &
                 * ( 1.0_RP/MAPF(i,j,2,I_UV) - 1.0_RP/MAPF(i-1,j,2,I_UV) ) * RCDX(i) &
                 - 0.25_RP * ( MOMX(k,i,j)+MOMX(k,i-1,j)+MOMX(k,i,j+1)+MOMX(k,i-1,j+1) ) &
                 * ( 1.0_RP/MAPF(i,j+1,1,I_XY) - 1.0_RP/MAPF(i,j,1,I_XY) ) * RFDY(j) ) &
               * 2.0_RP / ( DENS(k,i,j+1) + DENS(k,i,j) ) ! metoric term
          div = divdmp_coef / dtrk * ( DDIV(k,i,j+1)/MAPF(i,j+1,1,I_XY) - DDIV(k,i,j)/MAPF(i,j,1,I_XY) ) &
              * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) * FDY(j) ! divergence damping
          MOMY_RK(k,i,j) = MOMY0(k,i,j) &
                         + dtrk * ( ( advcv + advch - pg ) / GSQRT(k,i,j,I_XVZ) + cf + div + MOMY_t(k,i,j) )
#ifdef HIST_TEND
          if ( lhist ) then
             advcv_t(k,i,j,I_MOMY) = advcv / GSQRT(k,i,j,I_XVZ)
             advch_t(k,i,j,I_MOMY) = advch / GSQRT(k,i,j,I_XVZ)
             pg_t(k,i,j,3) = - pg / GSQRT(k,i,j,I_XVZ)
             cf_t(k,i,j,2) = cf
             ddiv_t(k,i,j,3) = div
          endif
#endif
       enddo
       enddo
       enddo
       PROFILE_STOP("hevi_momy")
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

#ifdef PROFILE_FIPP
       call fipp_stop()
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

       call HIST_in(wdmp_t (:,:,:),        'MOMZ_t_dwamp', 'tendency of momentum z (Rayleigh damping)',   'kg/m2/s2', zdim='half')

       call HIST_in(ddiv_t (:,:,:,1),      'MOMZ_t_ddiv',  'tendency of momentum z (divergence damping)', 'kg/m2/s2', zdim='half')
       call HIST_in(ddiv_t (:,:,:,2),      'MOMX_t_ddiv',  'tendency of momentum x (divergence damping)', 'kg/m2/s2', xdim='half')
       call HIST_in(ddiv_t (:,:,:,3),      'MOMY_t_ddiv',  'tendency of momentum y (divergence damping)', 'kg/m2/s2', ydim='half')

       call HIST_in(cf_t   (:,:,:,1),      'MOMX_t_cf',    'tendency of momentum x (coliolis force)',     'kg/m2/s2', xdim='half')
       call HIST_in(cf_t   (:,:,:,2),      'MOMY_t_cf',    'tendency of momentum y (coliolis force)',     'kg/m2/s2', ydim='half')
    endif
#endif

  end subroutine ATMOS_DYN_Tstep_short_fvm_hevi

#ifdef HEVI_BICGSTAB
!OCL SERIAL
  subroutine solve_bicgstab( &
       C,        & ! (inout)
       F1, F2, F3 ) ! (in)

    use scale_process, only: &
       PRC_MPIstop
    implicit none
    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: F1(KA)
    real(RP), intent(in)    :: F2(KA)
    real(RP), intent(in)    :: F3(KA)

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

    M(3,1) = F1(KS)
    M(2,1) = F2(KS)
    do k = KS+1, KE-2
       M(3,k-KS+1) = F1(k)
       M(2,k-KS+1) = F2(k)
       M(1,k-KS+1) = F3(k)
    enddo
    M(2,KE-KS) = F2(KE-1)
    M(1,KE-KS) = F3(KE-1)

    norm = 0.0_RP
    do k = 1, KMAX-1
       norm = norm + C(k)**2
    enddo

    r  => v0
    rn => v1

    call mul_matrix( v1, M, C )

    do k = 1, KMAX-1
       r(k) = C(k) - v1(k)
       r0(k) = r(k)
       p(k) = r(k)
    enddo

    r0r = r0(1) * r(1)
    do k = 2, KMAX-1
       r0r = r0r + r0(k)*r(k)
    enddo
    do iter = 1, KMAX-1
       error = 0.0_RP
       do k = 1, KMAX-1
          error = error + r(k)**2
       enddo

!       if ( error < epsilon .OR. error / norm < epsilon ) then
       if ( error/norm < epsilon ) then
#ifdef DEBUG
!          write(*,*) "Bi-CGSTAB converged:", iter
#endif
          exit
       endif

       call mul_matrix( ap, M, p )
       al = r0(1) * ap(1)
       do k = 2, KMAX-1
          al = al + r0(k)*ap(k)
       enddo
       al = r0r / al ! (r0,r) / (r0,Mp)
       s(:) = r(:) - al*ap(:)
       call mul_matrix( as, M, s )
       be = as(1) * s(1)  ! be is used as just work variable here
       w =  as(1) * as(1)
       do k = 2, KMAX-1
          be = be + as(k)*s(k)
          w  = w  + as(k)*as(k)
       enddo
       w = be / w ! (as,s) / (as,as)

       c(:) = c(:) + al*p(:) + w*s(:)
       rn(:) = s(:) - w*as(:)
       be = al/w / r0r
       r0r = r0(1) * rn(1)
       do k = 2, KMAX-1
          r0r = r0r + r0(k)*rn(k)
       enddo
       be = be * r0r ! al/w * (r0,rn)/(r0,r)
       p(:) = rn(:) + be * ( p(:) - w*ap(:) )

       swap => rn
       rn => r
       r => swap
    enddo

    if ( iter >= KMAX-1 ) then
       write(*,*) 'xxx [atmos_dyn_hevi] Bi-CGSTAB'
       write(*,*) 'xxx not converged', error, norm
       call PRC_MPIstop
    endif

    return
  end subroutine solve_bicgstab

  subroutine mul_matrix(V, M, C)
    implicit none
    real(RP), intent(out) :: V(KMAX-1)
    real(RP), intent(in)  :: M(3, KMAX-1)
    real(RP), intent(in)  :: C(KMAX-1)

    integer :: k

    ! k = 1
    V(1) = M(3,1)*C(2) + M(2,1)*C(1)
    do k = 2, KMAX-2
       V(k) = M(3,k)*C(k+1) + M(2,k)*C(k) + M(1,k)*C(k-1)
    enddo
    ! k = KE-1
    V(KMAX-1) = M(2,KMAX-1)*C(KMAX-1) + M(1,KMAX-1)*C(KMAX-2)

    return
  end subroutine mul_matrix

#elif defined(HEVI_LAPACK)
!OCL SERIAL
  subroutine solve_lapack( &
       C,   & ! (inout)
#ifdef DEBUG
       i, j, &
#endif
       F1, F2, F3 ) ! (in)
    use scale_process, only: &
         PRC_MPIstop
    implicit none

    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: F1(KA)
    real(RP), intent(in)    :: F2(KA)
    real(RP), intent(in)    :: F3(KA)
#ifdef DEBUG
    integer , intent(in)    :: i
    integer , intent(in)    :: j
#endif

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
       M(NB+1,k-KS+1) = F1(k-1)
       ! k, 0
       M(NB+2,k-KS+1) = F2(k)
       ! k+1, -1
       M(NB+3,k-KS+1) = F3(k+1)
    enddo
    ! KS, 0
    M(NB+2,1    ) = F2(KS)
    ! KS+1, -1
    M(NB+3,1    ) = F3(KS+1)
    ! KE-2, +1
    M(NB+1,KE-KS) = F1(KE-2)
    ! KE-1, 0
    M(NB+2,KE-KS) = F2(KE-1)

#ifdef DEBUG
          M2(:,:) = 0.0_RP
          do k = 1, KMAX-1
             if (k>1) M2(k-1,k) = M(NB+3,k-1)
             M2(k,k) = M(NB+2,k)
             if (k<kmax-1) M2(k+1,k) = M(NB+1,k+1)
             c2(k) = c(k)
          enddo
#endif
          if ( RP == DP ) then
             call DGBSV( KMAX-1, NB, NB, 1, M, NB*3+1, IPIV, C, KMAX-1, INFO)
          else
             call SGBSV( KMAX-1, NB, NB, 1, M, NB*3+1, IPIV, C, KMAX-1, INFO)
          endif
          ! C is (\rho w)^{n+1}
#ifdef DEBUG
          if ( INFO /= 0 ) then
             write(*,*) "DGBSV was failed", info
             call PRC_MPIstop
          endif

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
             endif
          enddo
#endif
  end subroutine solve_lapack
#else
!OCL SERIAL
  subroutine solve_direct( &
       C,         & ! (inout)
       F1, F2, F3 ) ! (in)

    use scale_process, only: &
       PRC_MPIstop
    implicit none
    real(RP), intent(inout) :: C(KMAX-1)
    real(RP), intent(in)    :: F1(KA)
    real(RP), intent(in)    :: F2(KA)
    real(RP), intent(in)    :: F3(KA)

    real(RP) :: e(KMAX-2)
    real(RP) :: f(KMAX-2)

    real(RP) :: rdenom

    integer :: k

    rdenom = 1.0_RP / F2(KS)
    e(1) = - F1(KS) * rdenom
    f(1) = C(1) * rdenom
    do k = 2, KMAX-2
       rdenom = 1.0_RP / ( F2(k+KS-1) + F3(k+KS-1) * e(k-1) )
       e(k) = - F1(k+KS-1) * rdenom
       f(k) = ( C(k) - F3(k+KS-1) * f(k-1) ) * rdenom
    enddo

    ! C = \rho w
    C(KMAX-1) = ( C(KMAX-1) - F3(KE-1) * f(KMAX-2) ) &
              / ( F2(KE-1) + F3(KE-1) * e(KMAX-2) ) ! C(KMAX-1) = f(KMAX-1)
    do k = KMAX-2, 1, -1
       C(k) = e(k) * C(k+1) + f(k)
    enddo

    return
  end subroutine solve_direct

#endif

#ifdef DEBUG
  subroutine check_equation( &
       VECT, &
       DENS, MOMZ, RHOT, DPRES, &
       REF_dens, &
       Sr, Sw, St, &
       J33G, G, &
       RT2P, &
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
    real(RP), intent(in) :: DPRES(KA)
    real(RP), intent(in) :: REF_dens(KA)
    real(RP), intent(in) :: Sr(KA)
    real(RP), intent(in) :: Sw(KA)
    real(RP), intent(in) :: St(KA)
    real(RP), intent(in) :: J33G
    real(RP), intent(in) :: G(KA,8)
    real(RP), intent(in) :: RT2P(KA)
    real(RP), intent(in) :: dt
    integer , intent(in) :: i
    integer , intent(in) :: j

    real(RP), parameter :: small = 1e-6_RP

    real(RP) :: MOMZ_N(KA)
    real(RP) :: DENS_N(KA)
    real(RP) :: RHOT_N(KA)
    real(RP) :: DPRES_N(KA)

    real(RP) :: POTT(KA)
    real(RP) :: PT(KA)

    real(RP) :: error, lhs, rhs
    integer :: k


    do k = KS, KE-1
       MOMZ_N(k) = VECT(k-KS+1)
    enddo
    MOMZ_N(:KS-1) = 0.0_RP
    MOMZ_N(KE:) = 0.0_RP

    ! density
    do k = KS+1, KE-1
       DENS_N(k) = DENS(k) &
            + dt * ( - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k) )
    enddo
    DENS_N(KS) = DENS(KS) &
         + dt * ( - J33G * MOMZ_N(KS) * RCDZ(KS) / G(KS,I_XYZ) + Sr(KS) )
    DENS_N(KE) = DENS(KE) &
         + dt * ( J33G * MOMZ_N(KE-1) * RCDZ(KE) / G(KE,I_XYZ) + Sr(KE) )

    ! rho*theta
    do k = KS, KE
       POTT(k) = RHOT(k) / DENS(k)
    enddo
    do k = KS+1, KE-2
       PT(k) = ( 7.0_RP * ( POTT(k+1) + POTT(k  ) ) &
                 -        ( POTT(k+2) + POTT(k-1) ) ) / 12.0_RP
    enddo
    PT(KS-1) = 0.0_RP
    PT(KS  ) = ( POTT(KS+1) + POTT(KS  ) ) * 0.5_RP
    PT(KE-1) = ( POTT(KE  ) + POTT(KE-1) ) * 0.5_RP
    PT(KE  ) = 0.0_RP
    do k = KS+1, KE-1
       RHOT_N(k) = RHOT(k) &
            + dt * ( - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) &
                     + St(k) )
    enddo
    RHOT_N(KS) = RHOT(KS) &
         + dt * ( - J33G * MOMZ_N(KS)*PT(KS) * RCDZ(KS) / G(KS,I_XYZ) + St(KS) )
    RHOT_N(KE) = RHOT(KE) &
         + dt * ( J33G * MOMZ_N(KE-1)*PT(KE-1) * RCDZ(KE) / G(KE-1,I_XYZ) + St(KE) )


    do k = KS, KE
       DPRES_N(k) = DPRES(k) + RT2P(k) * ( RHOT_N(k) - RHOT(k) )
    enddo

    do k = KS, KE
       lhs = ( DENS_N(k) - DENS(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k) - MOMZ_N(k-1) ) * RCDZ(k) / G(k,I_XYZ) + Sr(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       endif
       if ( abs(error) > small ) then
          write(*,*)"HEVI: DENS error", k, i, j, error, lhs, rhs
          write(*,*)eps
          call PRC_MPIstop
       endif
    enddo

    do k = KS, KE-1
       lhs = ( MOMZ_N(k) - MOMZ(k) ) / dt
       rhs = - J33G * ( DPRES_N(k+1) - DPRES_N(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( DENS_N(k+1) - REF_dens(k+1) + DENS_N(k) - REF_dens(k) ) * 0.5_RP &
             + Sw(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       endif
       if ( abs(error) > small ) then
          write(*,*)"HEVI: MOMZ error", k, i, j, error, lhs, rhs
          write(*,*) MOMZ_N(k), MOMZ(k), dt
          write(*,*) - J33G * ( DPRES(k+1) - DPRES(k) ) * RFDZ(k) / G(k,I_XYW) &
             - GRAV * ( DENS(k+1) -REF_dens(k+1) + DENS(k) -REF_dens(k) ) * 0.5_RP &
             + Sw(k)
          call PRC_MPIstop
       endif
    enddo

    do k = KS, KE
       lhs = ( RHOT_N(k) - RHOT(k) ) / dt
       rhs = - J33G * ( MOMZ_N(k)*PT(k) - MOMZ_N(k-1)*PT(k-1) ) * RCDZ(k) / G(k,I_XYZ) + St(k)
       if ( abs(lhs) < small ) then
          error = rhs
       else
          error = ( lhs - rhs ) / lhs
       endif
       if ( abs(error) > small ) then
          write(*,*)"HEVI: RHOT error", k, i, j, error, lhs, rhs
          call PRC_MPIstop
       endif
    enddo

    return
  end subroutine check_equation
#endif

end module scale_atmos_dyn_tstep_short_fvm_hevi
