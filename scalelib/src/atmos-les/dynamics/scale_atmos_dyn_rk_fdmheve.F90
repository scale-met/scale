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
module scale_atmos_dyn_rk_fdmheve
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

  use scale_atmos_numeric_fdm_def, only: &
       & FLXEVALTYPE_CD4, FLXEVALTYPE_UD1,  &
       & VL_ZXY, VL_ZUY, VL_ZXV, VL_WXY
  
  use scale_atmos_numeric_fdm, only: &
       & ATMOS_NUMERIC_FDM_setup,           &
       & ATMOS_NUMERIC_FDM_SetCoordMapInfo, &
       & ATMOS_NUMERIC_FDM_RhoVar2Var,      &
       & ATMOS_NUMERIC_FDM_EvalFlux,        &
       & ATMOS_NUMERIC_FDM_EvolveVar

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_rk_fdmheve_regist
  public :: ATMOS_DYN_rk_fdmheve_setup
  public :: ATMOS_DYN_rk_fdmheve

  public :: ATMOS_DYN_rk_fdmheve_SetFluxEvalType

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
  integer, parameter :: VA_HEVE = 0

  integer :: IFS_OFF
  integer :: JFS_OFF
#ifdef HIST_TEND
  real(RP), allocatable :: advch_t(:,:,:,:)
  real(RP), allocatable :: advcv_t(:,:,:,:)
  real(RP), allocatable :: ddiv_t(:,:,:,:)
  real(RP), allocatable :: pg_t(:,:,:,:)
  real(RP), allocatable :: cf_t(:,:,:,:)
#endif

  integer :: FlxEvalTypeID
  
  !-----------------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_rk_fdmheve_regist( &
       ATMOS_DYN_TYPE, &
       VA_out, &
       VAR_NAME, VAR_DESC, VAR_UNIT )
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    character(len=*),       intent(in)  :: ATMOS_DYN_TYPE
    integer,                intent(out) :: VA_out !< number of prognostic variables
    character(len=H_SHORT), intent(out) :: VAR_NAME(:) !< name  of the variables
    character(len=H_MID),   intent(out) :: VAR_DESC(:) !< desc. of the variables
    character(len=H_SHORT), intent(out) :: VAR_UNIT(:) !< unit  of the variables

    if( IO_L ) write(IO_FID_LOG,*) '*** FDM-HEVE Register'

    if ( ATMOS_DYN_TYPE .ne. 'FDM-HEVE' ) then
       write(*,*) 'xxx ATMOS_DYN_TYPE is not FDM-HEVE. Check!'
       call PRC_MPIstop
    endif

    VA_out = VA_HEVE
!!$    VAR_NAME(1) = 'NC_rk'
!!$    VAR_DESC(1) = 'NC advanced with rk schme'
!!$    VAR_UNIT(1) = 'num/kg'
    return
    
  end subroutine ATMOS_DYN_rk_fdmheve_regist

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_rk_fdmheve_setup( &
       BND_W, &
       BND_E, &
       BND_S, &
       BND_N  )
    implicit none

    logical, intent(in) :: BND_W
    logical, intent(in) :: BND_E
    logical, intent(in) :: BND_S
    logical, intent(in) :: BND_N
    !---------------------------------------------------------------------------

    IFS_OFF = 1
    JFS_OFF = 1
    if ( BND_W ) IFS_OFF = 0
    if ( BND_S ) JFS_OFF = 0


#ifdef HIST_TEND
    allocate( advch_t(KA,IA,JA,5) )
    allocate( advcv_t(KA,IA,JA,5) )
    allocate( ddiv_t(KA,IA,JA,3) )
    allocate( pg_t(KA,IA,JA,3) )
    allocate( cf_t(KA,IA,JA,2) )

    advch_t = 0.0_RP
    advcv_t = 0.0_RP
    ddiv_t = 0.0_RP
    pg_t = 0.0_RP
    cf_t = 0.0_RP
#endif

    !
    FlxEvalTypeID = FLXEVALTYPE_CD4
    call ATMOS_NUMERIC_FDM_setup(IFS_OFF, JFS_OFF)
    
    return
  end subroutine ATMOS_DYN_rk_fdmheve_setup

  !-----------------------------------------------------------------------------
  !> Specify type of finite difference scheme
  subroutine ATMOS_DYN_rk_fdmheve_SetFluxEvalType(FlxEvalTypeName)

    use scale_atmos_numeric_fdm, only: &
         & FlxEvalTypeName2ID
    
    character(len=H_SHORT), intent(in) :: FlxEvalTypeName

    FlxEvalTypeID = FlxEvalTypeName2ID(FlxEvalTypeName)
    
  end subroutine ATMOS_DYN_rk_fdmheve_SetFluxEvalType
  
  !-----------------------------------------------------------------------------
  !> Runge-Kutta loop
  subroutine ATMOS_DYN_rk_fdmheve( &
       DENS_RK, MOMZ_RK, MOMX_RK, MOMY_RK, RHOT_RK, &
       PROG_RK,                                     &
       mflx_hi, tflx_hi,                            &
       DENS0,   MOMZ0,   MOMX0,   MOMY0,   RHOT0,   &
       DENS,    MOMZ,    MOMX,    MOMY,    RHOT,    &
       DENS_t,  MOMZ_t,  MOMX_t,  MOMY_t,  RHOT_t,  &
       PROG0, PROG,                                 &
       Rtot, CVtot, CORIOLI,                        &
       num_diff, divdmp_coef,                       &
       FLAG_FCT_RHO, FLAG_FCT_MOMENTUM, FLAG_FCT_T, &
       FLAG_FCT_ALONG_STREAM,                       &
       CDZ, FDZ, FDX, FDY,                          &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,          &
       PHI, GSQRT, J13G, J23G, J33G, MAPF,          &
       REF_pres, REF_dens,                          &
       BND_W, BND_E, BND_S, BND_N,                  &
       dtrk, dt                                     )
    use scale_grid_index
    use scale_const, only: &
       GRAV   => CONST_GRAV,  &
       P00    => CONST_PRE00, &
       Rdry   => CONST_Rdry,  &
       CPdry  => CONST_CPdry, &
       CVdry  => CONST_CVdry
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

    real(RP), intent(out) :: DENS_RK(KA,IA,JA)   ! prognostic variables
    real(RP), intent(out) :: MOMZ_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMX_RK(KA,IA,JA)   !
    real(RP), intent(out) :: MOMY_RK(KA,IA,JA)   !
    real(RP), intent(out) :: RHOT_RK(KA,IA,JA)   !

    real(RP), intent(out) :: PROG_RK(KA,IA,JA,VA)  !

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! mass flux
    real(RP), intent(out)   :: tflx_hi(KA,IA,JA,3) ! internal energy flux

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

    real(RP), intent(in)  :: PROG0(KA,IA,JA,VA)
    real(RP), intent(in)  :: PROG (KA,IA,JA,VA)

    real(RP), intent(in)  :: Rtot    (KA,IA,JA)  ! total R
    real(RP), intent(in)  :: CVtot   (KA,IA,JA)  ! total CV
    real(RP), intent(in)  :: CORIOLI (1, IA,JA)
    real(RP), intent(in)  :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in)  :: divdmp_coef

    logical,  intent(in)  :: FLAG_FCT_RHO
    logical,  intent(in)  :: FLAG_FCT_MOMENTUM
    logical,  intent(in)  :: FLAG_FCT_T
    logical,  intent(in)  :: FLAG_FCT_ALONG_STREAM

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
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
    real(RP), intent(in)  :: REF_pres(KA,IA,JA)   !< reference pressure
    real(RP), intent(in)  :: REF_dens(KA,IA,JA)   !< reference density

    logical,  intent(in)  :: BND_W
    logical,  intent(in)  :: BND_E
    logical,  intent(in)  :: BND_S
    logical,  intent(in)  :: BND_N

    real(RP), intent(in)  :: dtrk
    real(RP), intent(in)  :: dt

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
#ifndef NO_FCT_DYN
    real(RP) :: qflx_lo  (KA,IA,JA,3)
    real(RP) :: qflx_anti(KA,IA,JA,3)
    real(RP) :: tflx_lo  (KA,IA,JA,3)
    real(RP) :: tflx_anti(KA,IA,JA,3)
    real(RP) :: DENS0_uvw(KA,IA,JA)
    real(RP) :: DENS_uvw (KA,IA,JA)
#endif

    real(RP) :: sw

    real(RP) :: one(KA,IA,JA)    
    real(RP) :: VARTMP(KA,IA,JA)
    real(RP) :: MOMZ_t_nadv(KA,IA,JA) 
    real(RP) :: MOMX_t_nadv(KA,IA,JA) 
    real(RP) :: MOMY_t_nadv(KA,IA,JA)    

    real(RP) :: advch ! horizontal advection
    real(RP) :: advcv ! vertical advection
    real(RP) :: div  ! divergence damping
    real(RP) :: pg   ! pressure gradient force
    real(RP) :: cf   ! colioris force
#ifdef HIST_TEND
    logical  :: lhist
#endif

#ifdef DRY
    real(RP) :: CPovCV
#endif

    integer  :: IIS, IIE
    integer  :: JJS, JJE
    integer  :: k, i, j

    real(RP), parameter ::  extdmp_coef = 0.01D0
    real(RP) :: extdiv
    
    !---------------------------------------------------------------------------

#ifdef DEBUG
    PRES(:,:,:) = UNDEF
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
    lhist = dt .eq. dtrk
#endif

#ifdef DRY
    CPovCV = CPdry / CVdry
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

      ! Momentum -> Velocity
      call ATMOS_NUMERIC_FDM_RhoVar2Var( VELX,                        &  ! (out)
         & MOMX, VL_ZUY, DENS,                                        &  ! (in)
         & IIS-IHALO, IIE+IHALO-1, JJS-JHALO+1, JJE+JHALO, KS, KE )      ! (in)

      call ATMOS_NUMERIC_FDM_RhoVar2Var( VELY,                        &  ! (out)
         & MOMY, VL_ZXV, DENS,                                        &  ! (in)
         & IIS-IHALO+1, IIE+IHALO, JJS-JHALO, JJE+JHALO-1, KS, KE )      ! (in)

      call ATMOS_NUMERIC_FDM_RhoVar2Var( VELZ,                        &  ! (out)
         & MOMZ, VL_WXY, DENS,                                        &  ! (in)
         & IIS-IHALO+1, IIE+IHALO, JJS-JHALO+1, JJE+JHALO, KS, KE-1 )    ! (in)

       ! pressure, pott. temp.

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS-JHALO, JJE+JHALO
       do i = IIS-IHALO, IIE+IHALO
          do k = KS, KE
#ifdef DEBUG
             call CHECK( __LINE__, RHOT(k,i,j) )
             call CHECK( __LINE__, Rtot(k,i,j) )
             call CHECK( __LINE__, CVtot(k,i,j) )
#endif
#ifdef DRY
             PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rdry / P00 )**CPovCV
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

      ! RHOT --> POTT 
      call ATMOS_NUMERIC_FDM_RhoVar2Var( POTT,                    &  ! (out)
         & RHOT, VL_ZXY, DENS,                                    &  ! (in)
         & IIS-IHALO, IIE+IHALO, JJS-JHALO, JJE+JHALO, KS, KE )      ! (in)
       

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

    enddo
    enddo

#ifndef NO_FCT_DYN    
    if ( FLAG_FCT_MOMENTUM ) then
       call COMM_vars8( VELZ(:,:,:), 1 )
       call COMM_vars8( VELX(:,:,:), 2 )
       call COMM_vars8( VELY(:,:,:), 3 )
       call COMM_wait ( VELZ(:,:,:), 1 )
       call COMM_wait ( VELX(:,:,:), 2 )
       call COMM_wait ( VELY(:,:,:), 3 )
    endif
#endif

    !
    one = 1.0_RP

    call ATMOS_NUMERIC_FDM_SetCoordMapInfo(      &
       & GSQRT, J13G, J23G, J33G, MAPF           & ! (in)
       & )

    !#####################################################################################
    ! continuity equation (total rho)
    !#####################################################################################
    
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

      !-----< high order flux >-----
      ! * Input (Phi * Vec_i / Fact)
      !   Phi = 1, Vec_i = MOM{X,Y,Z}, Fact = 1
      ! * Ouput
      !   [ GSQRT*MOMX/n, GSQRT*MOMY/m, GSQRT*(MOMZ*J33/(mn) + MOMX*J13/n + MOMY*J23/m) ]
      call ATMOS_NUMERIC_FDM_EvalFlux( mflx_hi,                           &  ! (inout)
         & FlxEvalTypeID, VL_ZXY, one, one, MOMX, MOMY, MOMZ, .true.,     &  ! (in)
         & IIS, IIE, JJS, JJE, KS, KE, num_diff(:,:,:,I_DENS,:))             ! (in)
    
      !-----< update density >-----
      call ATMOS_NUMERIC_FDM_EvolveVar( DENS_RK,                         & ! (out)
         & DENS0, VL_ZXY, mflx_hi, dtrk, RCDX, RCDY, RCDZ,               & ! (in)
         & IIS, IIE, JJS, JJE, KS, KE, DENS_t                            & ! (in)
#ifdef HIST_TEND
         & ,lhist, advch_t, advcv_t                                      &
#endif       
         & )
      
#ifndef NO_FCT_DYN
      if ( FLAG_FCT_RHO ) then
         !-----< low order flux >-----         
         ! * Input (Phi * Vec_i / Fact)
         !   Phi = DENS, Vec_i = MOM{X,Y,Z}, Fact = DENS
         ! * Ouput 
         !   [ GSQRT*MOMX/n, GSQRT*MOMY/m, GSQRT*(MOMZ*J33/(mn) + MOMX*J13/n + MOMY*J23/m) ]
         call ATMOS_NUMERIC_FDM_EvalFlux( qflx_lo,                                &  ! (inout)
              & FLXEVALTYPE_UD1, VL_ZXY, DENS, DENS, MOMX, MOMY, MOMZ, .true.,    &  ! (in)
              & IIS, IIE, JJS, JJE, KS, KE )                                         ! (in)
      endif

    enddo
    enddo

    if ( FLAG_FCT_RHO ) then
       if ( BND_W ) then
          do j = JS-1, JE+1
          do k = KS, KE
             qflx_lo(k,IS-1,j,XDIR) = mflx_hi(k,IS-1,j,XDIR)
          end do
          end do
       end if
       if ( BND_E ) then
          do j = JS-1, JE+1
          do k = KS, KE
             qflx_lo(k,IE,j,XDIR) = mflx_hi(k,IE,j,XDIR)
          end do
          end do
       end if
       if ( BND_S ) then
          do i = IS-1, IE+1
          do k = KS, KE
             qflx_lo(k,i,JS-1,YDIR) = mflx_hi(k,i,JS-1,YDIR)
          end do
          end do
       end if
       if ( BND_N ) then
          do i = IS-1, IE+1
          do k = KS, KE
             qflx_lo(k,i,JE,YDIR) = mflx_hi(k,i,JE,YDIR)
          end do
          end do
       end if
       call ATMOS_DYN_fct( qflx_anti,               & ! (out)
                           DENS0, one, one,         & ! (in)
                           mflx_hi, qflx_lo,        & ! (in)
                           mflx_hi,                 & ! (in)
                           RCDZ, RCDX, RCDY,        & ! (in)
                           GSQRT(:,:,:,I_XYZ),      & ! (in)
                           MAPF(:,:,:,I_XY), dtrk,  & ! (in)
                           FLAG_FCT_ALONG_STREAM    ) ! (in)

       !--- update rho       
       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

         VARTMP(:,IIS:IIE,JJS:JJE) = DENS_RK(:,IIS:IIE,JJS:JJE)
         call ATMOS_NUMERIC_FDM_EvolveVar( DENS_RK,                    & ! (out)
           & VARTMP, VL_ZXY, qflx_anti, dtrk, RCDX, RCDY, RCDZ,        & ! (in)
           & IIS, IIE, JJS, JJE, KS, KE )                              
         
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
          mflx_hi(k,i,j,ZDIR) = mflx_hi(k,i,j,ZDIR) - qflx_anti(k,i,j,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
          k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS  , JE
       do i = IS-1, IEH
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, qflx_anti(k,i,j,XDIR) )
#endif
          mflx_hi(k,i,j,XDIR) = mflx_hi(k,i,j,XDIR) - qflx_anti(k,i,j,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JS-1, JEH
       do i = IS  , IE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, mflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, qflx_anti(k,i,j,YDIR) )
#endif
          mflx_hi(k,i,j,YDIR) = mflx_hi(k,i,j,YDIR) - qflx_anti(k,i,j,YDIR)
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
#endif    ! [-- end ifndef NO_FCT_DYN ------------------------------------------]

    
    !###################################################################################
    ! momentum equation (z)
    !###################################################################################
    
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !-----< high order flux >-----
       ! * Input (Phi * Vec_i / Fact)
       !   Phi = MOMZ, Vec_i = MOM{X,Y,Z}, Fact = DENS
       ! * Ouput
       !   [ GSQRT*MOMZ*u/n, GSQRT*MOMZ*v/m, GSQRT*MOMZ*(u*J33/(mn) + u*J13/n + v*J23/m) ]
       call ATMOS_NUMERIC_FDM_EvalFlux( qflx_hi,                          &  ! (inout)
         & FlxEvalTypeID, VL_WXY, MOMZ, DENS, MOMX, MOMY, MOMZ, .true.,   &  ! (in)
         & IIS, IIE, JJS, JJE, KS, KE, num_diff(:,:,:,I_MOMZ,:) )            ! (in)

       !$omp parallel do private(i,j,k, div) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
         do k = KS, KE-1

           ! pressure gradient force at (x, y, w)             
           pgf(k,i,j) = J33G * ( DPRES(k+1,i,j)-DPRES(k,i,j) ) * RFDZ(k) ! [x,y,z]          

           ! buoyancy force at (x, y, w)
           buoy(k,i,j) = GRAV * 0.5_RP * (   GSQRT(k+1,i,j,I_XYZ) * ( DENS(k+1,i,j)-REF_dens(k+1,i,j) ) &
                                           + GSQRT(k  ,i,j,I_XYZ) * ( DENS(k  ,i,j)-REF_dens(k  ,i,j) ) ) ! [x,y,z]
           ! divergence damping
           div = divdmp_coef * FDZ(k) / dtrk * ( DDIV(k+1,i,j)-DDIV(k,i,j) ) ! divergence damping

           ! Total tendency except advection term
           MOMZ_t_nadv(k,i,j) = &
                &   ( - pgf(k,i,j) - buoy(k,i,j) ) / GSQRT(k,i,j,I_XYW) &
                & + div                                                 &
                & + MOMZ_t(k,i,j) 
#ifdef HIST_TEND
           if ( lhist ) then
              pg_t(k,i,j,1) = ( - pgf(k,i,j) - buoy(k,i,j) ) / GSQRT(k,i,j,I_XYW)
              ddiv_t(k,i,j,1) = div
           end if
#endif          
         enddo
       enddo
       enddo
    
      !-----< update momentum (z) -----
       call ATMOS_NUMERIC_FDM_EvolveVar( MOMZ_RK,                     & ! (out)
         & MOMZ0, VL_WXY, qflx_hi, dtrk, RCDX, RCDY, RFDZ,           & ! (in)
         & IIS, IIE, JJS, JJE, KS, KE-1, MOMZ_t_nadv                 & ! (in)
#ifdef HIST_TEND
         & ,lhist, advch_t, advcv_t                                  &
#endif       
         & )

       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
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
          end if
#endif
       enddo
       enddo

#ifndef NO_FCT_DYN
       if ( FLAG_FCT_MOMENTUM ) then

         !-----< high order flux >-----
         ! * Input (Phi * Vec_i / Fact)
         !   Phi = MOMZ, Vec_i = MOM{X,Y,Z}, Fact = DENS
         ! * Ouput
         !   [ GSQRT*MOMZ*u/n, GSQRT*MOMZ*v/m, GSQRT*MOMZ*(u*J33/(mn) + u*J13/n + v*J23/m) ]          
         call ATMOS_NUMERIC_FDM_EvalFlux( qflx_lo,                             &  ! (inout)
           & FLXEVALTYPE_UD1, VL_WXY, MOMZ, DENS, MOMX, MOMY, MOMZ, .true.,    &  ! (in)
           & IIS, IIE, JJS, JJE, KS, KE)                                          ! (in)
       endif
#endif       
    enddo
    enddo

#ifndef NO_FCT_DYN    
    if ( FLAG_FCT_MOMENTUM ) then

       call COMM_vars8( DENS_RK, 1 )
       call COMM_wait ( DENS_RK, 1, .false. )

       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE-1
          DENS0_uvw(k,i,j) = 0.5_RP * ( DENS0(k,i,j) + DENS0(k+1,i,j) )
          DENS_uvw(k,i,j) = 0.5_RP * ( DENS_RK(k,i,j) + DENS_RK(k+1,i,j) )
       end do
       end do
       end do

       do j = JS-1, JE+1
       do i = IS-1, IE+1
          DENS_uvw(KE,i,j) = DENS_uvw(KE-1,i,j)
          DENS0_uvw(KE,i,j) = DENS0_uvw(KE-1,i,j)
          VELZ(KE,i,j) = 0.0_RP
       end do
       end do

       call ATMOS_DYN_fct( qflx_anti,                 & ! (out)
                           VELZ, DENS0_uvw, DENS_uvw, & ! (in)
                           qflx_hi, qflx_lo,          & ! (in)
                           mflx_hi,                   & ! (in)
                           RFDZ, RCDX, RCDY,          & ! (in)
                           GSQRT(:,:,:,I_XYW),        & ! (in)
                           MAPF(:,:,:,I_XY), dtrk,    & ! (in)
                           FLAG_FCT_ALONG_STREAM      ) ! (in)

       !--- update momentum(z)       
       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

         VARTMP(KS:KE-1,IIS:IIE,JJS:JJE) = MOMZ_RK(KS:KE-1,IIS:IIE,JJS:JJE)
         call ATMOS_NUMERIC_FDM_EvolveVar( MOMZ_RK,                  & ! (out)
           & VARTMP, VL_WXY, qflx_anti, dtrk, RCDX, RCDY, RFDZ,      & ! (in)
           & IIS, IIE, JJS, JJE, KS, KE-1                            & ! (in)
           & )

       enddo
       enddo
    endif ! FLAG_FCT_MOMENTUM
#endif ! [------- End ifndef NO_FCT_DYN -------------------------------------------------]


    !#####################################################################################
    ! momentum equation (x)
    !#####################################################################################
    
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !-----< high order flux >-----
       ! * Input (Phi * Vec_i / Fact)
       !   Phi = MOMX, Vec_i = MOM{X,Y,Z}, Fact = DENS
       ! * Ouput
       !   [ GSQRT*MOMX*u/n, GSQRT*MOMX*v/m, GSQRT*MOMX*(w*J33/(mn) + u*J13/n + v*J23/m) ]          
      call ATMOS_NUMERIC_FDM_EvalFlux( qflx_hi,                             &  ! (inout)
        & FlxEvalTypeID, VL_ZUY, MOMX, DENS, MOMX, MOMY, MOMZ, .true.,      &  ! (in)
        & IIS, IIE, JJS, JJE, KS, KE, num_diff(:,:,:,I_MOMX,:) )               ! (in)
    
       ! pressure gradient force at (u, y, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          pgf(k,i,j) = ( ( GSQRT(k,i+1,j,I_XYZ) * DPRES(k,i+1,j) & ! [x,y,z]
                         - GSQRT(k,i  ,j,I_XYZ) * DPRES(k,i  ,j) & ! [x,y,z]
                         ) * RFDX(i) &
                       + ( J13G(k  ,i,j,I_UYW) &
                         * 0.25_RP * ( DPRES(k+1,i+1,j)+DPRES(k,i+1,j) &
                                     + DPRES(k+1,i  ,j)+DPRES(k,i  ,j) ) & ! [x,y,z->u,y,w]
                         - J13G(k-1,i,j,I_UYW) &
                         * 0.25_RP * ( DPRES(k,i+1,j)+DPRES(k-1,i+1,j) &
                                     + DPRES(k,i  ,j)+DPRES(k-1,i  ,j) ) & ! [x,y,z->u,y,w]
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
          cor(k,i,j) = 0.0625_RP * ( CORIOLI(1,i+1,j  )+CORIOLI(1,i,j  ) ) & ! [x,y,z->u,y,z]
                                 * ( DENS   (k,i+1,j  )+DENS   (k,i,j  ) ) & ! [x,y,z->u,y,z]
                                 * ( VELY   (k,i+1,j  )+VELY   (k,i,j  ) &
                                   + VELY   (k,i+1,j-1)+VELY   (k,i,j-1) ) &  ! [x,v,z->u,y,z]
                      + 0.25_RP * MAPF(i,j,1,I_UY) * MAPF(i,j,2,I_UY) &
                      * ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) &
                      * ( ( MOMY(k,i,j) + MOMY(k,i,j-1) + MOMY(k,i+1,j) + MOMY(k,i+1,j-1) ) * 0.25_RP &
                        * ( 1.0_RP/MAPF(i+1,j,2,I_XY) - 1.0_RP/MAPF(i,j,2,I_XY) ) * RCDX(i) &
                        - MOMX(k,i,j) &
                        * ( 1.0_RP/MAPF(i,j,1,I_UV) - 1.0_RP/MAPF(i,j-1,1,I_UV) ) * RFDY(j) ) &
                      * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) ) ! metric term
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k, div, extdiv) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, min(IIE, IEH)
         extdiv = 0.0_RP
         do k = KS, KE
           div = divdmp_coef * FDX(i) / dtrk * ( DDIV(k,i+1,j)-DDIV(k,i,j) ) ! divergence damping
!           extdiv = extdiv + CDZ(k) * extdmp_coef * div/divdmp_coef
           
           ! Total tendency except advection term
           MOMX_t_nadv(k,i,j) = &
               & - pgf(k,i,j) / GSQRT(k,i,j,I_UYZ) &
               & + cor(k,i,j)                      &
               & + div                             &
               & + MOMX_t(k,i,j)
         end do
        
!         MOMX_t_nadv(KS:KE,i,j) = MOMX_t_nadv(KS:KE,i,j) + extdiv / sum(CDZ(KS:KE))
       end do
       end do
       
       !-----< update momentum (x) >-----
       call ATMOS_NUMERIC_FDM_EvolveVar( MOMX_RK,                       & ! (out)
         & MOMX0, VL_ZUY, qflx_hi, dtrk, RFDX, RCDY, RCDZ,              & ! (in)
         & IIS, min(IIE,IEH), JJS, JJE, KS, KE, MOMX_t_nadv             & ! (in)
#ifdef HIST_TEND
         & ,lhist, advch_t, advcv_t                                     &
#endif       
         & )

       
#ifndef NO_FCT_DYN
       if ( FLAG_FCT_MOMENTUM ) then
         !-----< low order flux >-----
         ! * Input (Phi * Vec_i / Fact)
         !   Phi = MOMX, Vec_i = MOM{X,Y,Z}, Fact = DENS
         ! * Ouput
         !   [ GSQRT*MOMX*u/n, GSQRT*MOMX*v/m, GSQRT*MOMX*(w*J33/(mn) + u*J13/n + v*J23/m) ]                    
         call ATMOS_NUMERIC_FDM_EvalFlux( qflx_lo,                            &  ! (inout)
           & FLXEVALTYPE_UD1, VL_ZUY, MOMX, DENS, MOMX, MOMY, MOMZ, .true.,   &  ! (in)
           & IIS, IIE, JJS, JJE, KS, KE)                                         ! (in)
      endif
#endif      
    enddo
    enddo

#ifndef NO_FCT_DYN    
    if ( FLAG_FCT_MOMENTUM ) then

       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          DENS0_uvw(k,i,j) = 0.5_RP * ( DENS0(k,i,j) + DENS0(k,i+1,j) )
          DENS_uvw(k,i,j)  = 0.5_RP * ( DENS_RK(k,i,j) + DENS_RK(k,i+1,j) )
       end do
       end do
       end do

       call ATMOS_DYN_fct( qflx_anti,                 & ! (out)
                           VELX, DENS0_uvw, DENS_uvw, & ! (in)
                           qflx_hi, qflx_lo,          & ! (in)
                           mflx_hi,                   & ! (in)
                           RCDZ, RFDX, RCDY,          & ! (in)
                           GSQRT(:,:,:,I_UYZ),        & ! (in)
                           MAPF(:,:,:,I_UY), dtrk,    & ! (in)
                           FLAG_FCT_ALONG_STREAM      ) ! (in)

       !--- update momentum(x)       
       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

         VARTMP(:,IIS:min(IIE,IEH),JJS:JJE) = MOMX_RK(:,IIS:min(IIE,IEH),JJS:JJE)
         call ATMOS_NUMERIC_FDM_EvolveVar( MOMX_RK,                               & ! (out)
           & VARTMP, VL_ZUY, qflx_anti, dtrk, RFDX, RCDY, RCDZ,                   & ! (in)
           & IIS, min(IIE,IEH), JJS, JJE, KS, KE                                  & ! (in)
           & )
         
       enddo
       enddo

#ifdef DEBUG
       qflx_lo  (:,:,:,:) = UNDEF
       qflx_anti(:,:,:,:) = UNDEF
#endif

    endif ! FLAG_FCT_MOMENTUM
#endif   ! [-- End ifndef NO_FCT_DYN --------------------------------------------]
    
#ifdef DEBUG
    qflx_hi(:,:,:,:) = UNDEF
#endif

    !########################################################################
    ! momentum equation (y)
    !########################################################################
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

      !-----< high order flux >-----
      ! * Input (Phi * Vec_i / Fact)
      !   Phi = MOMY, Vec_i = MOM{X,Y,Z}, Fact = DENS
      ! * Ouput
      !   [ GSQRT*MOMY*u/n, GSQRT*MOMY*v/m, GSQRT*MOMY*(w*J33/(mn) + u*J13/n + v*J23/m) ]                    
      call ATMOS_NUMERIC_FDM_EvalFlux( qflx_hi,                            &  ! (inout)
        & FlxEvalTypeID, VL_ZXV, MOMY, DENS, MOMX, MOMY, MOMZ, .true.,     &  ! (in)
        & IIS, IIE, JJS, JJE, KS, KE, num_diff(:,:,:,I_MOMY,:) )              ! (in)
    
       ! pressure gradient force at (x, v, z)

       !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          pgf(k,i,j) = ( ( GSQRT(k,i,j+1,I_XYZ) * DPRES(k,i,j+1) & ! [x,y,z]
                         - GSQRT(k,i,j  ,I_XYZ) * DPRES(k,i,j  ) & ! [x,y,z]
                         ) * RFDY(j) &
                       + ( J23G(k  ,i,j,I_XVW) &
                         * 0.25_RP * ( DPRES(k+1,i,j+1)+DPRES(k,i,j+1) &
                                     + DPRES(k+1,i,j  )+DPRES(k,i,j  ) ) & ! [x,y,z->x,v,w]
                         - J23G(k-1,i,j,I_XVW) &
                         * 0.25_RP * ( DPRES(k,i,j+1)+DPRES(k-1,i,j+1) &
                                     + DPRES(k,i,j  )+DPRES(k-1,i,j  ) ) & ! [x,y,z->x,v,w]
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
          cor(k,i,j) = - 0.0625_RP * ( CORIOLI(1,i  ,j+1)+CORIOLI(1,i  ,j) ) & ! [x,y,z->x,v,z]
                                   * ( DENS   (k,i  ,j+1)+DENS   (k,i  ,j) ) & ! [x,y,z->x,v,z]
                                   * ( VELX   (k,i  ,j+1)+VELX   (k,i  ,j) &
                                     + VELX   (k,i-1,j+1)+VELX   (k,i-1,j) ) & ! [u,y,z->x,v,z]
                     - 0.25_RP * MAPF(i,j,1,I_XV) * MAPF(i,j,2,I_XV) * GSQRT(k,i,j,I_XVZ) &
                     * ( MOMX(k,i,j) + MOMX(k,i-1,j) + MOMX(k,i,j+1) + MOMX(k,i-1,j+1) )&
                     * ( MOMY(k,i,j) &
                       * ( 1.0_RP/MAPF(i,j,2,I_UV) - 1.0_RP/MAPF(i-1,j,2,I_UV) ) * RCDX(i) &
                       - 0.25_RP * ( MOMX(k,i,j)+MOMX(k,i-1,j)+MOMX(k,i,j+1)+MOMX(k,i-1,j+1) ) &
                       * ( 1.0_RP/MAPF(i,j+1,1,I_XY) - 1.0_RP/MAPF(i,j,1,I_XY) ) * RFDY(j) ) &
                     * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) ) ! metoric term
       enddo
       enddo
       enddo

       !$omp parallel do private(i,j,k, div, extdiv) OMP_SCHEDULE_ collapse(2)
       do j = JJS, min(JJE,JEH)
       do i = IIS, IIE
         extdiv = 0.0_RP 
         do k = KS, KE            
           div = divdmp_coef / dtrk * FDY(j) * ( DDIV(k,i,j+1)-DDIV(k,i,j) ) ! divergence damping
!           extdiv = extdiv + CDZ(k) * extdmp_coef * div/divdmp_coef

           ! Total tendency except advection term          
           MOMY_t_nadv(k,i,j) = &
               & - pgf(k,i,j) / GSQRT(k,i,j,I_XVZ) &
               & + cor(k,i,j)                      &
               & + div                             &
               & + MOMY_t(k,i,j)

#ifdef HIST_TEND
            if ( lhist ) then
              pg_t(k,i,j,3) = - pgf(k,i,j) / GSQRT(k,i,j,I_XVZ)
              cf_t(k,i,j,2) = cor(k,i,j)
              ddiv_t(k,i,j,3) = div
            end if
#endif
!           MOMY_t_nadv(KS:KE,i,j) = MOMY_t_nadv(KS:KE,i,j) + extdiv / sum(CDZ(KS:KE))
          
         end do
       end do
       end do
       
       !-----< update momentum (y) >-----
       call ATMOS_NUMERIC_FDM_EvolveVar( MOMY_RK,                     & ! (out)
         & MOMY0, VL_ZXV, qflx_hi, dtrk, RCDX, RFDY, RCDZ,            & ! (in)
         & IIS, IIE, JJS, min(JJE,JEH), KS, KE, MOMY_t_nadv           & ! (in)
#ifdef HIST_TEND
         & ,lhist, advch_t, advcv_t                                   &
#endif       
         & )
       
#ifndef NO_FCT_DYN
       if ( FLAG_FCT_MOMENTUM ) then
         !-----< low order flux >-----
         ! * Input (Phi * Vec_i / Fact)
         !   Phi = MOMY, Vec_i = MOM{X,Y,Z}, Fact = DENS
         ! * Ouput
         !   [ GSQRT*MOMY*u/n, GSQRT*MOMY*v/m, GSQRT*MOMY*(w*J33/(mn) + u*J13/n + v*J23/m) ]                              
         call ATMOS_NUMERIC_FDM_EvalFlux( qflx_lo,                           &  ! (inout)
           & FLXEVALTYPE_UD1, VL_ZXV, MOMY, DENS, MOMX, MOMY, MOMZ, .true.,  &  ! (in)
           & IIS, IIE, JJS, JJE, KS, KE )                                       ! (in)
       endif
#endif
       
    enddo
    enddo

#ifndef NO_FCT_DYN     
    if ( FLAG_FCT_MOMENTUM ) then

       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          DENS0_uvw(k,i,j) = 0.5_RP * ( DENS0(k,i,j) + DENS0(k,i,j+1) )
          DENS_uvw(k,i,j) = 0.5_RP * ( DENS_RK(k,i,j) + DENS_RK(k,i,j+1) )
       end do
       end do
       end do

       call ATMOS_DYN_fct( qflx_anti,                 & ! (out)
                           VELY, DENS0_uvw, DENS_uvw, & ! (in)
                           qflx_hi, qflx_lo,          & ! (in)
                           mflx_hi,                   & ! (in)
                           RCDZ, RCDX, RFDY,          & ! (in)
                           GSQRT(:,:,:,I_XVZ),        & ! (in)
                           MAPF(:,:,:,I_XV), dtrk,    & ! (in)
                           FLAG_FCT_ALONG_STREAM      ) ! (in)

       !--- update momentum (y)
       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

         VARTMP(:,IIS:IIE,JJS:min(JJE,JEH)) = MOMY_RK(:,IIS:IIE,JJS:min(JJE,JEH))
         call ATMOS_NUMERIC_FDM_EvolveVar( MOMY_RK,                      & ! (out)
           & VARTMP, VL_ZXV, qflx_anti, dtrk, RCDX, RFDY, RCDZ,          & ! (in)
           & IIS, IIE, JJS, min(JJE,JEH), KS, KE-1                       & ! (in)
           & )

       enddo
       enddo

#ifdef DEBUG
       qflx_lo  (:,:,:,:) = UNDEF
       qflx_anti(:,:,:,:) = UNDEF
#endif
    endif ! FLAG_FCT_MOMENTUM
#endif   ! [-- End ifndef NO_FCT_DYN --------------------------------------------]

#ifdef DEBUG
    qflx_hi(KS:,:,:,:) = UNDEF
#endif
    
    !####################################################################################
    ! Thermodynamic equation
    !####################################################################################
    
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

      !-----< high order flux >-----
      ! * Input (Phi * Vec_i / Fact)
      !   Phi = POTT, Vec_i = MomFlx{X,Y,Z}, Fact = 1.0
      ! * Ouput
      !   [ POTT*MOMFlxX, POTT*MOMFlxY, POTT*MOMFlxZ ]
      call ATMOS_NUMERIC_FDM_EvalFlux( tflx_hi,                                              &  ! (inout)
        & FlxEvalTypeID, VL_ZXY,                                                             &  ! (in)
        & POTT, one, mflx_hi(:,:,:,XDIR), mflx_hi(:,:,:,YDIR), mflx_hi(:,:,:,ZDIR), .false., &  ! (in)
        & IIS, IIE, JJS, JJE, KS, KE, num_diff(:,:,:,I_RHOT,:) )                                ! (in)

      !-----< update rho*theta >-----
      call ATMOS_NUMERIC_FDM_EvolveVar( RHOT_RK,                      & ! (out)
         & RHOT0, VL_ZXY, tflx_hi, dtrk, RCDX, RCDY, RCDZ,            & ! (in)
         & IIS, IIE, JJS, JJE, KS, KE, RHOT_t                         & ! (in)
#ifdef HIST_TEND
         & ,lhist, advch_t, advcv_t                                   &
#endif       
         & )

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

       if ( .not. FLAG_FCT_MOMENTUM ) then
          call COMM_vars8( DENS_RK, 1 )
          call COMM_wait ( DENS_RK, 1, .false. )
       end if

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

         !-----< low order flux >-----
         ! * Input (Phi * Vec_i / Fact)
         !   Phi = POTT, Vec_i = MomFlx{X,Y,Z}, Fact = 1.0
         ! * Ouput
         !   [ POTT*MOMFlxX, POTT*MOMFlxY, POTT*MOMFlxZ ]
         VARTMP(:,IIS:IIE,JJS:JJE) = 1.0_RP
         call ATMOS_NUMERIC_FDM_EvalFlux( tflx_hi,                                                 &  ! (inout)
           & FLXEVALTYPE_UD1, VL_ZXY,                                                              &  ! (in)
           & POTT, VARTMP, mflx_hi(:,:,:,XDIR), mflx_hi(:,:,:,YDIR), mflx_hi(:,:,:,ZDIR), .false., &  ! (in)
           & IIS, IIE, JJS, JJE, KS, KE )                                                      ! (in)
         
         !$omp parallel do private(i,j,k) OMP_SCHEDULE_ collapse(2)
         do j = JS-1, JE+1
         do i = IS-1, IE+1
         do k = KS, KE
            POTT(k,i,j) = RHOT0(k,i,j) / DENS0(k,i,j)
         enddo
         enddo
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

       !--- update rho*theta       
       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

         VARTMP(:,IIS:IIE,JJS:JJE) = RHOT_RK(:,IIS:IIE,JJS:JJE)
         call ATMOS_NUMERIC_FDM_EvolveVar( RHOT_RK,                   & ! (out)
           & VARTMP, VL_ZXY, tflx_anti, dtrk, RCDX, RCDY, RCDZ,       & ! (in)
           & IIS, IIE, JJS, JJE, KS, KE                               & ! (in)
           & )
          
       enddo
       enddo

    endif ! FLAG_FCT_T
#endif    ! [--- End ifndef NO_FCT_DYN ---------------------------------------------------------- ]

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
  end subroutine ATMOS_DYN_rk_fdmheve

end module scale_atmos_dyn_rk_fdmheve
