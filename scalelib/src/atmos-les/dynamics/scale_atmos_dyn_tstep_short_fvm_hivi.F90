!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          HIVI FVM scheme for Atmospheric dynamical process
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
module scale_atmos_dyn_tstep_short_fvm_hivi
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
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
  public :: ATMOS_DYN_Tstep_short_fvm_hivi_regist
  public :: ATMOS_DYN_Tstep_short_fvm_hivi_setup
  public :: ATMOS_DYN_Tstep_short_fvm_hivi

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
  integer,  private, parameter :: VA_FVM_HIVI = 0

  integer,  private            :: ITMAX
  real(RP), private            :: epsilon

  integer,  private            :: mtype ! MPI DATATYPE

  ! tentative
  real(RP), private, parameter :: FACT_N =  7.0_RP / 12.0_RP
  real(RP), private, parameter :: FACT_F = -1.0_RP / 12.0_RP

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tstep_short_fvm_hivi_regist( &
       ATMOS_DYN_TYPE, &
       VA_out,         &
       VAR_NAME,       &
       VAR_DESC,       &
       VAR_UNIT        )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*),       intent(in)  :: ATMOS_DYN_TYPE
    integer,                intent(out) :: VA_out      !< number of prognostic variables
    character(len=H_SHORT), intent(out) :: VAR_NAME(:) !< name  of the variables
    character(len=H_MID),   intent(out) :: VAR_DESC(:) !< desc. of the variables
    character(len=H_SHORT), intent(out) :: VAR_UNIT(:) !< unit  of the variables
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HIVI Register'

    if ( ATMOS_DYN_TYPE .ne. 'FVM-HIVI' .or. ATMOS_DYN_TYPE .ne. 'HIVI' ) then
       write(*,*) 'xxx ATMOS_DYN_TYPE is not FVM-HIVI. Check!'
       call PRC_MPIstop
    endif

    VA_out = VA_FVM_HIVI

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hivi_regist

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_short_fvm_hivi_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_ATMOS_DYN_TSTEP_FVM_HIVI / &
         ITMAX, &
         EPSILON

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** HIVI Setup'
#ifdef HIVI_BICGSTAB
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING Bi-CGSTAB'
#else
    if ( IO_L ) write(IO_FID_LOG,*) '*** USING Multi-Grid'
    write(*,*) 'xxx Not Implemented yet'
    call PRC_MPIstop
#endif

    ! currently, vertical difference scheme for potential temperature is the CD4
    ITMAX = 100
    epsilon = 0.1_RP ** (RP*2)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_TSTEP_FVM_HIVI,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_DYN_TSTEP_FVM_HIVI. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_DYN_TSTEP_FVM_HIVI)

    if ( RP == DP ) then
       mtype = MPI_DOUBLE_PRECISION
    elseif( RP == SP ) then
       mtype = MPI_REAL
    else
       write(*,*) 'xxx Unsupported precision'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hivi_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_Tstep_short_fvm_hivi( &
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
       GRAV   => CONST_GRAV,   &
       P00    => CONST_PRE00
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
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
    real(RP) :: POTT(KA,IA,JA)  ! potential temperature [K]
    real(RP) :: DDENS(KA,IA,JA) ! density deviation from reference
    real(RP) :: DPRES(KA,IA,JA) ! pressure deviation from reference
    real(RP) :: DPRES_N(KA,IA,JA) ! pressure deviation at t=n+1 [Pa]

    real(RP) :: Sr(KA,IA,JA)
    real(RP) :: Sw(KA,IA,JA)
    real(RP) :: Su(KA,IA,JA)
    real(RP) :: Sv(KA,IA,JA)
    real(RP) :: St(KA,IA,JA)
    real(RP) :: RCs2(KA,IA,JA)
    real(RP) :: B(KA,IA,JA)

    real(RP) :: duvw

    real(RP) :: qflx_hi (KA,IA,JA,3)
    real(RP) :: qflx_J13(KA,IA,JA)
    real(RP) :: qflx_J23(KA,IA,JA)
    real(RP) :: mflx_hi2(KA,IA,JA,3)

    ! for implicit solver
    real(RP) :: M(7,KA,IA,JA)
    real(RP) :: r(KA,IA,JA)
    real(RP) :: p(KA,IA,JA)

    real(RP) :: zero(KA,IA,JA)

    real(RP) :: rdt

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j

    zero = 0.0_RP

#ifdef DEBUG
    DPRES(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF

    mflx_hi(:,:,:,:) = UNDEF
    mflx_hi(KS-1,:,:,ZDIR) = 0.0_RP

    qflx_hi(:,:,:,:) = UNDEF
    tflx_hi(:,:,:,:) = UNDEF
    qflx_J13(:,:,:)  = UNDEF
    qflx_J23(:,:,:)  = UNDEF
    mflx_hi2(:,:,:,:) = UNDEF

    Sr(:,:,:) = UNDEF
    Sw(:,:,:) = UNDEF
    Su(:,:,:) = UNDEF
    Sv(:,:,:) = UNDEF
    St(:,:,:) = UNDEF
    RCs2(:,:,:) = UNDEF

    B(:,:,:) = UNDEF
    r(:,:,:) = UNDEF
    p(:,:,:) = UNDEF
#endif

    rdt = 1.0_RP / dtrk

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
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

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, DENS(k,i,j) )
          call CHECK( __LINE__, REF_dens(k,i,j) )
#endif
          DDENS(k,i,j) = DENS(k,i,j) - REF_dens(k,i,j)
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

       !--- update momentum(z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_hi (k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_J13(k  ,i  ,j  )      )
          call CHECK( __LINE__, qflx_J13(k-1,i  ,j  )      )
          call CHECK( __LINE__, qflx_J23(k  ,i  ,j  )      )
          call CHECK( __LINE__, qflx_J23(k-1,i  ,j  )      )
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi (k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_hi (k  ,i  ,j-1,YDIR) )
          call CHECK( __LINE__, DDIV(k  ,i,j) )
          call CHECK( __LINE__, DDIV(k+1,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
          call CHECK( __LINE__, MOMZ_t(k,i,j) )
#endif
          Sw(k,i,j) = &
               - ( ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) &
                   + qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j  ) &
                   + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j  )      ) * RFDZ(k) &
                 + ( qflx_hi (k,i,j,XDIR) - qflx_hi (k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                 + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k  ,i  ,j-1,YDIR) ) * RCDY(j) &
                 ) / GSQRT(k,i,j,I_XYW) &
               + divdmp_coef * rdt  * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) * FDZ(k) & ! divergence damping
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
          Su(k,i,j) = ( &
               - ( ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) &
                   + qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j) &
                   + qflx_J13(k,i,j)      - qflx_J23(k-1,i  ,j)        ) * RCDZ(k) &
                 + ( qflx_hi (k,i,j,XDIR) - qflx_hi (k  ,i-1,j  ,XDIR) ) * RFDX(i) &
                 + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
               - ( J13G(k+1,i,j,I_UYZ) * ( DPRES(k+1,i+1,j)+DPRES(k+1,i,j) ) &
                 - J13G(k-1,i,j,I_UYZ) * ( DPRES(k-1,i+1,j)+DPRES(k-1,i,j) ) ) &
                 * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) &
               ) / GSQRT(k,i,j,I_UYZ) &
               + 0.125_RP * ( CORIOLI(1,i,j)+CORIOLI(1,i+1,j) ) &
                          * ( MOMY(k,i,j)+MOMY(k,i+1,j)+MOMY(k,i,j-1)+MOMY(k,i+1,j-1) ) & ! coriolis force
               + divdmp_coef * rdt * ( DDIV(k,i+1,j)-DDIV(k,i,j) ) * FDX(i) & ! divergence damping
               + MOMX_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### momentum equation (y) #####

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
#endif
          Sv(k,i,j) = &
               ( - ( ( qflx_hi (k,i,j,ZDIR) - qflx_hi (k-1,i  ,j  ,ZDIR) &
                     + qflx_J13(k,i,j)      - qflx_J13(k-1,i  ,j  ) &
                     + qflx_J23(k,i,j)      - qflx_J23(k-1,i  ,j  )      ) * RCDZ(k) &
                   + ( qflx_hi (k,i,j,XDIR) - qflx_hi (k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                   + ( qflx_hi (k,i,j,YDIR) - qflx_hi (k  ,i  ,j-1,YDIR) ) * RFDY(j) ) &
                 - ( J23G(k+1,i,j,I_XVZ) * ( DPRES(k+1,i,j+1)+DPRES(k+1,i,j) ) &
                   - J23G(k-1,i,j,I_XVZ) * ( DPRES(k-1,i,j+1)+DPRES(k-1,i,j) ) ) &
                   * 0.5_RP / ( FDZ(k+1)+FDZ(k) ) & ! pressure gradient force
               ) / GSQRT(k,i,j,I_XVZ) &
               - 0.125_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) &
                          * ( MOMX(k,i,j+1)+MOMX(k,i,j)+MOMX(k,i-1,j+1)+MOMX(k,i-1,j) ) & ! coriolis force
               + divdmp_coef * rdt * ( DDIV(k,i,j+1)-DDIV(k,i,j) ) * FDY(j) & ! divergence damping
               + MOMY_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !##### Thermodynamic Equation #####


       ! at (x, y, w)
       call ATMOS_DYN_FVM_fluxZ_XYZ( tflx_hi(:,:,:,ZDIR), & ! (out)
            mflx_hi2(:,:,:,ZDIR), POTT, GSQRT(:,:,:,I_XYW), & ! (in)
            num_diff(:,:,:,I_RHOT,ZDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, y, z)
       call ATMOS_DYN_FVM_fluxX_XYZ( tflx_hi(:,:,:,XDIR), & ! (out)
            mflx_hi2(:,:,:,XDIR), POTT, GSQRT(:,:,:,I_UYZ), & ! (in)
            num_diff(:,:,:,I_RHOT,XDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, v, z)
       call ATMOS_DYN_FVM_fluxY_XYZ( tflx_hi(:,:,:,YDIR), & ! (out)
            mflx_hi2(:,:,:,YDIR), POTT, GSQRT(:,:,:,I_XVZ), & ! (in)
            num_diff(:,:,:,I_RHOT,YDIR), & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

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
#endif
          St(k,i,j) = &
               - ( ( tflx_hi(k,i,j,ZDIR) - tflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                 + ( tflx_hi(k,i,j,XDIR) - tflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                 + ( tflx_hi(k,i,j,YDIR) - tflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) &
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
       ! at (x, y, w)
       call ATMOS_DYN_FVM_fluxZ_XYZ( tflx_hi(:,:,:,ZDIR), & ! (out)
            mflx_hi(:,:,:,ZDIR), POTT, GSQRT(:,:,:,I_XYW), & ! (in)
            zero, & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (u, y, z)
       call ATMOS_DYN_FVM_fluxX_XYZ( tflx_hi(:,:,:,XDIR), & ! (out)
            mflx_hi(:,:,:,XDIR), POTT, GSQRT(:,:,:,I_UYZ), & ! (in)
            zero, & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! at (x, v, z)
       call ATMOS_DYN_FVM_fluxY_XYZ( tflx_hi(:,:,:,YDIR), & ! (out)
            mflx_hi(:,:,:,YDIR), POTT, GSQRT(:,:,:,I_XVZ), & ! (in)
            zero, & ! (in)
            CDZ, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       ! rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT_RK(k,i,j) = RHOT0(k,i,j) &
               + dtrk * ( - ( ( tflx_hi(k,i,j,ZDIR) - tflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                            + ( tflx_hi(k,i,j,XDIR) - tflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( tflx_hi(k,i,j,YDIR) - tflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                            / GSQRT(k,i,j,I_XYZ) &
                          + St(k,i,j) )
       enddo
       enddo
       enddo

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

    enddo
    enddo

    ! implicit solver
#ifdef HIVI_BICGSTAB

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
          call CHECK( __LINE__, RT2P(k-1,i,j) )
          call CHECK( __LINE__, RT2P(k  ,i,j) )
          call CHECK( __LINE__, RT2P(k+1,i,j) )
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
             + GSQRT(k,i,j,I_XYZ) * ( St(k,i,j) - DPRES(k,i,j) / RT2P(k,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( ( DPRES(k+1,i,j)*RCs2(k+1,i,j) &
                               - DPRES(k-1,i,j)*RCs2(k-1,i,j) ) &
                             - ( RHOT(k+1,i,j)*DDENS(k+1,i,j)/DENS(k+1,i,j) &
                               - RHOT(k-1,i,j)*DDENS(k-1,i,j)/DENS(k-1,i,j) ) &
                             - ( DENS(k+1,i,j)*(RHOT_RK(k+1,i,j)/DENS_RK(k+1,i,j) - POTT(k+1,i,j) ) &
                               - DENS(k-1,i,j)*(RHOT_RK(k-1,i,j)/DENS_RK(k-1,i,j) - POTT(k-1,i,j) ) ) &
                             ) / ( FDZ(k) + FDZ(k-1) )
       enddo
       enddo
       enddo
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
          call CHECK( __LINE__, RT2P(KS  ,i,j) )
          call CHECK( __LINE__, RT2P(KS+1,i,j) )
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
             + GSQRT(KS,i,j,I_XYZ) * ( St(KS,i,j) - DPRES(KS,i,j) / RT2P(KS,i,j) * rdt ) &
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
          call CHECK( __LINE__, RT2P(KS  ,i,j) )
          call CHECK( __LINE__, RT2P(KS+1,i,j) )
          call CHECK( __LINE__, RT2P(KS+2,i,j) )
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
             + GSQRT(KS+1,i,j,I_XYZ) * ( St(KS+1,i,j) - DPRES(KS+1,i,j) / RT2P(KS+1,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( ( DPRES(KS+2,i,j)/RT2P(KS+2,i,j) &
                               - DPRES(KS  ,i,j)/RT2P(KS  ,i,j) ) &
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
          call CHECK( __LINE__, RT2P(KE-2,i,j) )
          call CHECK( __LINE__, RT2P(KE-1,i,j) )
          call CHECK( __LINE__, RT2P(KE  ,i,j) )
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
             + GSQRT(KE-1,i,j,I_XYZ) * ( St(KE-1,i,j) - DPRES(KE-1,i,j) / RT2P(KE-1,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( ( DPRES(KE  ,i,j)/RT2P(KE  ,i,j) &
                               - DPRES(KE-2,i,j)/RT2P(KE-2,i,j) ) &
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
          call CHECK( __LINE__, RT2P(KE-1,i,j) )
          call CHECK( __LINE__, RT2P(KE  ,i,j) )
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
             + GSQRT(KE,i,j,I_XYZ) * ( St(KE,i,j) - DPRES(KE,i,j) / RT2P(KE,i,j) * rdt ) &
             ) * rdt &
             + GRAV * J33G * ( - DPRES(KE-1,i,j)/RT2P(KE-1,i,j) &
                               + RHOT(KE-1,i,j)*DDENS(KE-1,i,j)/DENS(KE-1,i,j) &
                               - DENS(KE-1,i,j)*(RHOT_RK(KE-1,i,j)/DENS_RK(KE-1,i,j) - POTT(KE-1,i,j) ) ) &
                             / ( FDZ(KE) + FDZ(KE-1) )
       enddo
       enddo
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

    enddo
    enddo

    call solve_bicgstab( &
       DPRES_N, & ! (out)
       DPRES, & ! (in)
       M, B   ) ! (in)

#ifdef DEBUG
    call check_solver( DPRES_N, M, B )
#endif

    call COMM_vars8( DPRES_N, 1 )
    call COMM_wait ( DPRES_N, 1 )

#else

    call solve_multigrid

#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (z) #####

       !--- update momentum(z)
       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, DPRES_N(k+1,i,j) )
          call CHECK( __LINE__, DPRES_N(k  ,i,j) )
          call CHECK( __LINE__, DPRES(k+1,i,j) )
          call CHECK( __LINE__, DPRES(k  ,i,j) )
          call CHECK( __LINE__, DDENS(k+1,i,j) )
          call CHECK( __LINE__, DDENS(k  ,i,j) )
          call CHECK( __LINE__, REF_DENS(k+1,i,j) )
          call CHECK( __LINE__, REF_DENS(k,i,j) )
          call CHECK( __LINE__, MOMZ0(k,i,j) )
#endif
          duvw =  dtrk * ( ( &
                   - J33G * ( DPRES_N(k+1,i,j) - DPRES_N(k,i,j) ) * RFDZ(k) &
                   ) / GSQRT(k,i,j,I_UYZ) &
                 - 0.5_RP * GRAV &
                 * ( DDENS(k+1,i,j) &
                   + ( DPRES_N(k+1,i,j) - DPRES(k+1,i,j) ) &
                   * POTT(k+1,i,j) / RT2P(k+1,i,j) &
                   + DDENS(k  ,i,j) &
                   + ( DPRES_N(k  ,i,j) - DPRES(k  ,i,j) ) &
                   * POTT(k  ,i,j) / RT2P(k  ,i,j) ) &
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
          call CHECK( __LINE__, DPRES_N(k,i+1,j) )
          call CHECK( __LINE__, DPRES_N(k,i  ,j) )
          call CHECK( __LINE__, GSQRT(k,i+1,j,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i  ,j,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
          call CHECK( __LINE__, Su(k,i,j) )
          call CHECK( __LINE__, MOMX0(k,i,j) )
#endif
          duvw = dtrk * ( ( &
                 - ( GSQRT(k,i+1,j,I_XYZ) * DPRES_N(k,i+1,j) &
                   - GSQRT(k,i  ,j,I_XYZ) * DPRES_N(k,i  ,j) ) * RFDX(i) &  ! pressure gradient force
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
          call CHECK( __LINE__, DPRES_N(k,i,j  ) )
          call CHECK( __LINE__, DPRES_N(k,i,j+1) )
          call CHECK( __LINE__, GSQRT(k,i,j  ,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i,j+1,I_XYZ) )
          call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
          call CHECK( __LINE__, Sv(k,i,j) )
          call CHECK( __LINE__, MOMY0(k,i,j) )
#endif
          duvw = dtrk * ( ( &
                            - ( GSQRT(k,i,j+1,I_XYZ) * DPRES_N(k,i,j+1) &
                              - GSQRT(k,i,j  ,I_XYZ) * DPRES_N(k,i,j  ) ) * RFDY(j) &
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
          tflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) &
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
          tflx_hi(KS-1,i,j,ZDIR) = 0.0_RP
          tflx_hi(KS  ,i,j,ZDIR) = mflx_hi(KS  ,i,j,ZDIR) * 0.5_RP * ( POTT(KS+1,i,j) + POTT(KS  ,i,j) )
          tflx_hi(KE-1,i,j,ZDIR) = mflx_hi(KE-1,i,j,ZDIR) * 0.5_RP * ( POTT(KE  ,i,j) + POTT(KE-1,i,j) )
          tflx_hi(KE  ,i,j,ZDIR) = 0.0_RP
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
          tflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) &
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
          tflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) &
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
               + dtrk * ( - ( ( tflx_hi(k,i,j,ZDIR) - tflx_hi(k-1,i  ,j  ,ZDIR) ) * RCDZ(k) &
                            + ( tflx_hi(k,i,j,XDIR) - tflx_hi(k  ,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( tflx_hi(k,i,j,YDIR) - tflx_hi(k  ,i  ,j-1,YDIR) ) * RCDY(j) ) &
                            / GSQRT(k,i,j,I_XYZ) &
                          + St(k,i,j) )
       enddo
       enddo
       enddo

       !##### continuous equation #####

       ! total momentum flux
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) + mflx_hi2(k,i,j,ZDIR)
       enddo
       enddo
       enddo
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
          mflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) + mflx_hi2(k,i,j,XDIR)
       enddo
       enddo
       enddo
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
          mflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) + mflx_hi2(k,i,j,YDIR)
       enddo
       enddo
       enddo

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
            DPRES_N, DPRES, RHOT_RK, RHOT, &
            RT2P )
#endif

    enddo
    enddo

    return
  end subroutine ATMOS_DYN_Tstep_short_fvm_hivi

#ifdef HIVI_BICGSTAB
  subroutine solve_bicgstab( &
       DPRES_N, &
       DPRES, &
       M, B )
    use scale_process, only: &
       PRC_MPIstop
    use scale_comm, only: &
       COMM_world, &
       COMM_vars8, &
       COMM_wait
    implicit none
    real(RP), intent(out) :: DPRES_N(KA,IA,JA)
    real(RP), intent(in)  :: DPRES(KA,IA,JA)
    real(RP), intent(in)  :: M(7,KA,IA,JA)
    real(RP), intent(in)  :: B(KA,IA,JA)

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
       enddo
       enddo
       enddo
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
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          r0(k,i,j) = r(k,i,j)
          p(k,i,j) = r(k,i,j)
       enddo
       enddo
       enddo
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
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    iprod(1) = r0r
    iprod(2) = norm
    call MPI_AllReduce(iprod, buf, 2, mtype, MPI_SUM, COMM_world, ierror)
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
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       call MPI_AllReduce(error, buf, 1, mtype, MPI_SUM, COMM_world, ierror)
       error = buf(1)

#ifdef DEBUG
       if (IO_L) write(*,*) iter, error/norm
#endif
       if ( sqrt(error/norm) < epsilon .or. error > error2 ) then
#ifdef DEBUG
         IF ( IO_L ) write(*,*) "Bi-CGSTAB converged:", iter
#endif
          exit
       endif
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
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       call MPI_AllReduce(iprod, buf, 1, mtype, MPI_SUM, COMM_world, ierror)
       al = r0r / buf(1) ! (r0,r) / (r0,Mp)

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, r (k,i,j) )
          call CHECK( __LINE__, Mp(k,i,j) )
#endif
          s(k,i,j) = r(k,i,j) - al*Mp(k,i,j)
       enddo
       enddo
       enddo
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
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       call MPI_AllReduce(iprod, buf, 2, mtype, MPI_SUM, COMM_world, ierror)
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
          DPRES_N(k,i,j) = DPRES(k,i,j) + al*p(k,i,j) + w*s(k,i,j)
       enddo
       enddo
       enddo
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
       enddo
       enddo
       enddo
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
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

       call MPI_AllReduce(iprod, r0r, 1, mtype, MPI_SUM, COMM_world, ierror)

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
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       swap => rn
       rn => r
       r => swap
#ifdef DEBUG
       rn(:,:,:) = UNDEF
#endif
    enddo

    if ( iter >= ITMAX ) then
       write(*,*) 'xxx [atmos_dyn_hivi] Bi-CGSTAB'
       write(*,*) 'xxx not converged', error, norm
       call PRC_MPIstop
    endif

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
    enddo
    enddo
    enddo
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
    enddo
    enddo
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
    enddo
    enddo
    enddo
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
    enddo
    enddo
    enddo
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
    enddo
    enddo
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
       endif
    enddo
    enddo
    enddo

  end subroutine check_solver

  subroutine check_pres( &
       DPRES_N, DPRES, &
       RHOT_RK, RHOT, &
       RT2P )
    real(RP), intent(in) :: DPRES_N(KA,IA,JA)
    real(RP), intent(in) :: DPRES(KA,IA,JA)
    real(RP), intent(in) :: RHOT_RK(KA,IA,JA)
    real(RP), intent(in) :: RHOT(KA,IA,JA)
    real(RP), intent(in) :: RT2P(KA,IA,JA)

    real(RP) :: lhs, rhs
    integer :: k,i,j

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       lhs = DPRES_N(k,i,j) - DPRES(k,i,j)
       rhs = RT2P(k,i,j) * ( RHOT_RK(k,i,j) - RHOT(k,i,j) )
       if ( abs( lhs - rhs ) / lhs > 1e-15 ) then
          write(*,*) "error is too large: ", k,i,j, lhs, rhs
          call abort
       endif
    enddo
    enddo
    enddo
  end subroutine check_pres
#endif


end module scale_atmos_dyn_tstep_short_fvm_hivi
