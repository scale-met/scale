!-------------------------------------------------------------------------------
!> module coupler / physics / surface skin
!!
!! @par Description
!!          Skin surface model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_cpl_phy_sfc_skin
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_PHY_SFC_skin_setup
  public :: CPL_PHY_SFC_skin_finalize
  public :: CPL_PHY_SFC_skin

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
  integer,  private :: CPL_PHY_SFC_SKIN_itr_max ! maximum iteration number

  real(RP), private :: CPL_PHY_SFC_SKIN_dTS_max ! maximum delta surface temperature [K/s]
  real(RP), private :: CPL_PHY_SFC_SKIN_res_min ! minimum value of residual
  real(RP), private :: CPL_PHY_SFC_SKIN_err_min ! minimum value of error

  logical,  private :: initialized = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_PHY_SFC_SKIN_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_CPL_PHY_SFC_SKIN / &
       CPL_PHY_SFC_SKIN_itr_max, &
       CPL_PHY_SFC_SKIN_dTS_max, &
       CPL_PHY_SFC_SKIN_res_min, &
       CPL_PHY_SFC_SKIN_err_min

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( initialized ) return

    LOG_NEWLINE
    LOG_INFO("CPL_PHY_SFC_SKIN_setup",*) 'Setup'

    CPL_PHY_SFC_SKIN_itr_max = 100

    CPL_PHY_SFC_SKIN_dTS_max = 5.0E-2_RP
    CPL_PHY_SFC_SKIN_res_min = 1.0E+0_RP
    CPL_PHY_SFC_SKIN_err_min = 1.0E-2_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_PHY_SFC_SKIN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CPL_PHY_SFC_SKIN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CPL_PHY_SFC_SKIN_setup",*) 'Not appropriate names in namelist PARAM_CPL_PHY_SFC_SKIN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CPL_PHY_SFC_SKIN)

    initialized = .true.

    return
  end subroutine CPL_PHY_SFC_SKIN_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine CPL_PHY_SFC_SKIN_finalize

    initialized = .false.

    return
  end subroutine CPL_PHY_SFC_SKIN_finalize

  !-----------------------------------------------------------------------------
  subroutine CPL_PHY_SFC_skin( &
       IA, IS, IE,          &
       JA, JS, JE,          &
       TMPA, PRSA,          &
       WA, UA, VA,          &
       RHOA, QVA,           &
       LH, Z1, PBL,         &
       RHOS, PRSS,          &
       RFLXD,               &
       TG, WSTR, QVEF,      &
       ALBEDO,              &
       Rb, TC_dZ,           &
       Z0M, Z0H, Z0E,       &
       calc_flag, dt,       &
       model_name,          &
       TMPS,                &
       ZMFLX, XMFLX, YMFLX, &
       SHFLX, LHFLX, QVFLX, &
       GFLX,                &
       Ustar, Tstar, Qstar, &
       Wstar,               &
       RLmo,                &
       U10, V10, T2, Q2     )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       EPS   => CONST_EPS, &
       UNDEF => CONST_UNDEF, &
       PRE00 => CONST_PRE00, &
       TEM00 => CONST_TEM00, &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       Rvap  => CONST_Rvap,  &
       STB   => CONST_STB
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_dens2qsat_all
!       qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_LHV, &
       ATMOS_HYDROMETEOR_LHS, &
       CV_WATER, &
       CV_ICE,   &
       LHF
    use scale_bulkflux, only: &
       BULKFLUX, &
       BULKFLUX_diagnose_surface
    implicit none

    integer,          intent(in)    :: IA, IS, IE
    integer,          intent(in)    :: JA, JS, JE
    real(RP),         intent(in)    :: TMPA     (IA,JA)                     ! temperature at the lowest atmospheric layer [K]
    real(RP),         intent(in)    :: PRSA     (IA,JA)                     ! pressure    at the lowest atmospheric layer [Pa]
    real(RP),         intent(in)    :: WA       (IA,JA)                     ! velocity w  at the lowest atmospheric layer [m/s]
    real(RP),         intent(in)    :: UA       (IA,JA)                     ! velocity u  at the lowest atmospheric layer [m/s]
    real(RP),         intent(in)    :: VA       (IA,JA)                     ! velocity v  at the lowest atmospheric layer [m/s]
    real(RP),         intent(in)    :: RHOA     (IA,JA)                     ! density     at the lowest atmospheric layer [kg/m3]
    real(RP),         intent(in)    :: QVA      (IA,JA)                     ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP),         intent(in)    :: LH       (IA,JA)                     ! latent heat [J/kg]
    real(RP),         intent(in)    :: Z1       (IA,JA)                     ! cell center height at the lowest atmospheric layer [m]
    real(RP),         intent(in)    :: PBL      (IA,JA)                     ! the top of atmospheric mixing layer [m]
    real(RP),         intent(in)    :: RHOS     (IA,JA)                     ! density  at the surface [kg/m3]
    real(RP),         intent(in)    :: PRSS     (IA,JA)                     ! pressure at the surface [Pa]
    real(RP),         intent(in)    :: RFLXD    (IA,JA,N_RAD_DIR,N_RAD_RGN) ! downward radiation flux at the surface (direct/diffuse,IR/near-IR/VIS) [J/m2/s]
    real(RP),         intent(in)    :: TG       (IA,JA)                     ! subsurface temperature [K]
    real(RP),         intent(in)    :: WSTR     (IA,JA)                     ! amount of water storage [kg/m2]
    real(RP),         intent(in)    :: QVEF     (IA,JA)                     ! efficiency of evaporation (0-1)
    real(RP),         intent(in)    :: ALBEDO   (IA,JA,N_RAD_DIR,N_RAD_RGN) ! surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
    real(RP),         intent(in)    :: Rb       (IA,JA)                     ! stomata resistance [1/s]
    real(RP),         intent(in)    :: TC_dZ    (IA,JA)                     ! thermal conductivity / depth between surface and subsurface [J/m2/s/K]
    real(RP),         intent(in)    :: Z0M      (IA,JA)                     ! roughness length for momemtum [m]
    real(RP),         intent(in)    :: Z0H      (IA,JA)                     ! roughness length for heat     [m]
    real(RP),         intent(in)    :: Z0E      (IA,JA)                     ! roughness length for vapor    [m]
    logical,          intent(in)    :: calc_flag(IA,JA)                     ! to decide calculate or not
    real(DP),         intent(in)    :: dt                                   ! delta time
    character(len=*), intent(in)    :: model_name

    real(RP),         intent(inout) :: TMPS     (IA,JA)                     ! surface temperature [K]

    real(RP),         intent(out)   :: ZMFLX    (IA,JA)                     ! z-momentum      flux at the surface [kg/m/s2]
    real(RP),         intent(out)   :: XMFLX    (IA,JA)                     ! x-momentum      flux at the surface [kg/m/s2]
    real(RP),         intent(out)   :: YMFLX    (IA,JA)                     ! y-momentum      flux at the surface [kg/m/s2]
    real(RP),         intent(out)   :: SHFLX    (IA,JA)                     ! sensible heat   flux at the surface [J/m2/s]
    real(RP),         intent(out)   :: LHFLX    (IA,JA)                     ! latent heat     flux at the surface [J/m2/s]
    real(RP),         intent(out)   :: QVFLX    (IA,JA)                     ! water vapor     flux at the surface [kg/m2/s]
    real(RP),         intent(out)   :: GFLX     (IA,JA)                     ! subsurface heat flux at the surface [J/m2/s]
    real(RP),         intent(out)   :: Ustar    (IA,JA)                     ! friction velocity         [m/s]
    real(RP),         intent(out)   :: Tstar    (IA,JA)                     ! temperature scale         [K]
    real(RP),         intent(out)   :: Qstar    (IA,JA)                     ! moisture scale            [kg/kg]
    real(RP),         intent(out)   :: Wstar    (IA,JA)                     ! convective velocity scale [m/s]
    real(RP),         intent(out)   :: RLmo     (IA,JA)                     ! inversed Obukhov length   [1/m]
    real(RP),         intent(out)   :: U10      (IA,JA)                     ! velocity u  at 10m [m/s]
    real(RP),         intent(out)   :: V10      (IA,JA)                     ! velocity v  at 10m [m/s]
    real(RP),         intent(out)   :: T2       (IA,JA)                     ! temperature at 2m  [K]
    real(RP),         intent(out)   :: Q2       (IA,JA)                     ! water vapor at 2m  [kg/kg]

    real(RP), parameter :: dTS0     = 1.0E-4_RP ! delta surface temp.
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0E+0_RP ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5E+0_RP ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1E+0_RP ! factor b in Tomita (2009)

    real(RP) :: TMPS1(IA,JA)

    real(RP) :: emis   ! surface longwave emission                 [J/m2/s]
    real(RP) :: LWD    ! surface downward longwave  radiation flux [J/m2/s]
    real(RP) :: LWU    ! surface upward   longwave  radiation flux [J/m2/s]
    real(RP) :: SWD    ! surface downward shortwave radiation flux [J/m2/s]
    real(RP) :: SWU    ! surface upward   shortwave radiation flux [J/m2/s]
    real(RP) :: flx_qv ! surface upward qv flux                    [kg/m2/s]
    real(RP) :: res    ! residual

    real(RP) :: dres   ! d(residual)/dTMPS
    real(RP) :: oldres ! residual in previous step
    real(RP) :: redf   ! reduced factor
    real(RP) :: dts    ! temperature change
    real(RP) :: olddts ! temperature change in previous step

    real(RP) :: dUstar      ! friction velocity difference               [m/s]
    real(RP) :: dTstar      ! friction potential temperature difference  [K]
    real(RP) :: dQstar      ! friction water vapor mass ratio difference [kg/kg]
    real(RP) :: dWstar      ! free convection velocity scale difference  [m/s]
    real(RP) :: dRLmo       ! inversed Obukhov length         [1/m]
    real(RP) :: Uabs, dUabs ! modified absolute velocity      [m/s]
    real(RP) :: Ra,   dRa   ! Aerodynamic resistance (=1/Ce)  [1/s]

    real(RP) :: QVsat, dQVsat    ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS(IA,JA), dQVS ! water vapor mixing ratio at surface            [kg/kg]
    real(RP) :: Rtot             ! total gas constant
    real(RP) :: qdry             ! dry air mass ratio [kg/kg]

    real(RP) :: FracU10(IA,JA), dFracU10 ! calculation parameter for U10 [1]
    real(RP) :: FracT2 (IA,JA), dFracT2  ! calculation parameter for T2  [1]
    real(RP) :: FracQ2 (IA,JA), dFracQ2  ! calculation parameter for Q2  [1]

    real(RP) :: MFLUX

#ifdef _OPENACC
    logical :: err_flag
#endif

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'coupler / physics / surface / SKIN'

    !$acc data copyin(TMPA,PRSA,WA,UA,VA,RHOA,QVA,LH,Z1,PBL,RHOS,PRSS,RFLXD,TG,WSTR,QVEF,ALBEDO,Rb,TC_dZ,Z0M,Z0H,Z0E,calc_flag) &
    !$acc      copy(TMPS) &
    !$acc      copyout(ZMFLX,XMFLX,YMFLX,SHFLX,LHFLX,QVFLX,GFLX,Ustar,Tstar,Qstar,Wstar,RLmo,U10,V10,T2,Q2) &
    !$acc      create(TMPS1,QVS,FracU10,FracT2,FracQ2)

    ! copy surfce temperature for iteration
    !$omp parallel do
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       TMPS1(i,j) = TMPS(i,j)
    enddo
    enddo
    !$acc end kernels

#ifdef _OPENACC
    err_flag = .false.
#endif

    ! update surface temperature
    !$omp parallel do schedule(dynamic) collapse(2) &
#ifndef __GFORTRAN__
    !$omp default(none) &
    !$omp shared(IO_UNIVERSALRANK,IO_LOCALRANK,IO_JOBID,IO_DOMAINID) &
    !$omp shared(IS,IE,JS,JE,EPS,UNDEF,Rdry,CPdry,PRC_myrank,IO_FID_LOG,IO_L,model_name,bulkflux, &
    !$omp        CPL_PHY_SFC_SKIN_itr_max,CPL_PHY_SFC_SKIN_dTS_max,CPL_PHY_SFC_SKIN_err_min,CPL_PHY_SFC_SKIN_res_min, &
    !$omp        calc_flag,dt,QVA,QVS,TMPA,TMPS,PRSA,RHOA,WA,UA,VA,LH,Z1,PBL, &
    !$omp        TG,PRSS,RHOS,TMPS1,WSTR,QVEF,Z0M,Z0H,Z0E,Rb,TC_dZ,ALBEDO,RFLXD, &
    !$omp        FracU10,FracT2,FracQ2, &
    !$omp        ZMFLX,XMFLX,YMFLX,SHFLX,LHFLX,QVFLX,GFLX,Ustar,Tstar,Qstar,Wstar,RLmo,U10,V10,T2,Q2) &
#else
    !$omp default(shared) &
#endif
    !$omp private(qdry,Rtot,flx_qv,redf,res,dts,olddts,emis,LWD,LWU,SWD,SWU,dres,oldres,dQVS, &
    !$omp         QVsat,dQVsat,dUstar,dTstar,dQstar,dWstar,dFracU10,dFracT2,dFracQ2, &
    !$omp         Uabs,dUabs,dRLmo,Ra,dRa,MFLUX)
    !$acc parallel
    !$acc loop collapse(2) reduction(.or.:err_flag) independent &
    !$acc private(qdry,Rtot,flx_qv,redf,res,dts,olddts,emis,LWD,LWU,SWD,SWU,dres,oldres,dQVS, &
    !$acc         QVsat,dQVsat,dUstar,dTstar,dQstar,dWstar,dFracU10,dFracT2,dFracQ2, &
    !$acc         Uabs,dUabs,dRLmo,Ra,dRa,MFLUX)
    do j = JS, JE
    do i = IS, IE
       if ( calc_flag(i,j) ) then

!          qdry = 1.0_RP - QVA(i,j)
!          Rtot = qdry * Rdry + QVA(i,j) * Rvap

          redf   = 1.0_RP
          oldres = huge(0.0_RP)
          olddts = CPL_PHY_SFC_SKIN_dTS_max * dt

          ! modified Newton-Raphson method (Tomita 2009)
          !$acc loop seq
          do n = 1, CPL_PHY_SFC_SKIN_itr_max

             call qsat( TMPS1(i,j),      RHOS(i,j), QVsat  )
             call qsat( TMPS1(i,j)+dTS0, RHOS(i,j), dQVsat )
!             call qsat( TMPS1(i,j),      PRSS(i,j), qdry, QVsat  )
!             call qsat( TMPS1(i,j)+dTS0, PRSS(i,j), qdry, dQVsat )

             QVS(i,j) = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
                      + (        QVEF(i,j) ) * QVsat
             dQVS     = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
                      + (        QVEF(i,j) ) * dQVsat

             Uabs = sqrt( WA(i,j)**2 + UA(i,j)**2 + VA(i,j)**2 )

             call BULKFLUX( TMPA(i,j), TMPS1(i,j),                 & ! [IN]
                            PRSA(i,j), PRSS (i,j),                 & ! [IN]
                            QVA (i,j), QVS  (i,j),                 & ! [IN]
                            Uabs, Z1(i,j), PBL(i,j),               & ! [IN]
                            Z0M(i,j), Z0H(i,j), Z0E(i,j),          & ! [IN]
                            Ustar(i,j), Tstar(i,j), Qstar(i,j),    & ! [OUT]
                            Wstar(i,j), RLmo(i,j), Ra,             & ! [OUT]
                            FracU10(i,j), FracT2(i,j), FracQ2(i,j) ) ! [OUT]

             call BULKFLUX( TMPA (i,j), TMPS1(i,j)+dTS0,  & ! [IN]
                            PRSA (i,j), PRSS (i,j),       & ! [IN]
                            QVA  (i,j), dQVS,             & ! [IN]
                            Uabs, Z1(i,j), PBL(i,j),      & ! [IN]
                            Z0M(i,j), Z0H(i,j), Z0E(i,j), & ! [IN]
                            dUstar, dTstar, dQstar,       & ! [OUT]
                            dWstar, dRLmo, dRa,           & ! [OUT] ! not used
                            dFracU10, dFracT2, dFracQ2    ) ! [OUT] ! not used

             emis = ( 1.0_RP - ALBEDO(i,j,I_R_diffuse,I_R_IR) ) * STB * TMPS1(i,j)**4

             LWD  = RFLXD(i,j,I_R_diffuse,I_R_IR)
             LWU  = RFLXD(i,j,I_R_diffuse,I_R_IR)  * ALBEDO(i,j,I_R_diffuse,I_R_IR) + emis
             SWD  = RFLXD(i,j,I_R_direct ,I_R_NIR) &
                  + RFLXD(i,j,I_R_diffuse,I_R_NIR) &
                  + RFLXD(i,j,I_R_direct ,I_R_VIS) &
                  + RFLXD(i,j,I_R_diffuse,I_R_VIS)
             SWU  = RFLXD(i,j,I_R_direct ,I_R_NIR) * ALBEDO(i,j,I_R_direct ,I_R_NIR) &
                  + RFLXD(i,j,I_R_diffuse,I_R_NIR) * ALBEDO(i,j,I_R_diffuse,I_R_NIR) &
                  + RFLXD(i,j,I_R_direct ,I_R_VIS) * ALBEDO(i,j,I_R_direct ,I_R_VIS) &
                  + RFLXD(i,j,I_R_diffuse,I_R_VIS) * ALBEDO(i,j,I_R_diffuse,I_R_VIS)

             ! calculation for residual
             flx_qv = min( - RHOS(i,j) * Ustar(i,j) * Qstar(i,j) * Ra / ( Ra+Rb(i,j) ), WSTR(i,j)/real(dt,RP) )
             res = SWD - SWU + LWD - LWU                         &
                 + CPdry   * RHOS(i,j) * Ustar(i,j) * Tstar(i,j) &
                 - LH(i,j) * flx_qv                              &
                 - TC_dZ(i,j) * ( TMPS1(i,j) - TG(i,j) )

             ! calculation for d(residual)/dTMPS
             dres = -4.0_RP * emis / TMPS1(i,j) &
                  + CPdry   * RHOS(i,j) * ( Ustar(i,j)*(dTstar-Tstar(i,j))/dTS0 + Tstar(i,j)*(dUstar-Ustar(i,j))/dTS0 )                       &
                  + LH(i,j) * RHOS(i,j) * ( Ustar(i,j)*(dQstar-Qstar(i,j))/dTS0 + Qstar(i,j)*(dUstar-Ustar(i,j))/dTS0 ) * Ra / ( Ra+Rb(i,j) ) &
                  - TC_dZ(i,j)

             ! convergence test with residual and error levels
             if (      abs(res     ) < CPL_PHY_SFC_SKIN_res_min &
                  .OR. abs(res/dres) < CPL_PHY_SFC_SKIN_err_min ) then
                exit
             endif

             ! calculate reduced factor
             if ( dres < 0.0_RP ) then
                if ( abs(res) > abs(oldres) ) then
                   redf = max( TFa*abs(redf), redf_min )
                else
                   redf = min( TFb*abs(redf), redf_max )
                endif
             else
                redf = -1.0_RP
             endif

             ! estimate next surface temperature
             dts = - redf * res / dres
             dts = sign( min( abs(dts), abs(olddts) ), dts )
             TMPS1(i,j) = TMPS1(i,j) + dts

             ! save residual in this step
             oldres = res
             olddts = dts
          enddo

          ! update surface temperature with limitation
          TMPS1(i,j) = min( max( TMPS1(i,j),                                                 &
                                 TMPS (i,j) - CPL_PHY_SFC_SKIN_dTS_max * real(dt,kind=RP) ), &
                                 TMPS (i,j) + CPL_PHY_SFC_SKIN_dTS_max * real(dt,kind=RP) )

          if ( n > CPL_PHY_SFC_SKIN_itr_max ) then
             ! surface temperature was not converged
#ifdef _OPENACC
             LOG_WARN("CPL_PHY_SFC_skin",*) 'surface tempearture was not converged. '
!             LOG_WARN("CPL_PHY_SFC_skin",*) 'surface tempearture was not converged. ', model_name
#else
             LOG_WARN("CPL_PHY_SFC_skin",*) 'surface tempearture was not converged. ', trim(model_name)
             LOG_NEWLINE
             LOG_WARN_CONT('(A,I32)'   ) 'number of i                        [no unit]  :', i
             LOG_WARN_CONT('(A,I32)'   ) 'number of j                        [no unit]  :', j
             LOG_NEWLINE
             LOG_WARN_CONT('(A,I32)'   ) 'loop number                        [no unit]  :', n
             LOG_WARN_CONT('(A,F32.16)') 'Residual                           [J/m2/s]   :', res
             LOG_WARN_CONT('(A,F32.16)') 'delta Residual                     [J/m2/s]   :', dres
             LOG_NEWLINE
             LOG_WARN_CONT('(A,F32.16)') 'temperature                        [K]        :', TMPA  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'pressure                           [Pa]       :', PRSA  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'velocity w                         [m/s]      :', WA    (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'velocity u                         [m/s]      :', UA    (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'velocity v                         [m/s]      :', VA    (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'absolute velocity                  [m/s]      :', Uabs
             LOG_WARN_CONT('(A,F32.16)') 'density                            [kg/m3]    :', RHOA  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'water vapor mass ratio             [kg/kg]    :', QVA   (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'cell center height                 [m]        :', Z1    (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'atmospheric mixing layer height    [m]        :', PBL   (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'pressure at the surface            [Pa]       :', PRSS  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'downward radiation (IR, direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_IR )
             LOG_WARN_CONT('(A,F32.16)') 'downward radiation (IR, diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_IR )
             LOG_WARN_CONT('(A,F32.16)') 'downward radiation (NIR,direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_NIR)
             LOG_WARN_CONT('(A,F32.16)') 'downward radiation (NIR,diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_NIR)
             LOG_WARN_CONT('(A,F32.16)') 'downward radiation (VIS,direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_VIS)
             LOG_WARN_CONT('(A,F32.16)') 'downward radiation (VIS,diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_VIS)
             LOG_NEWLINE
             LOG_WARN_CONT('(A,F32.16)') 'soil temperature                   [K]        :', TG    (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'soil water                         [kg/m2]    :', WSTR  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'surface temperature                [K]        :', TMPS  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'surface density                    [kg/m3]    :', RHOS  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'efficiency of evaporation          [1]        :', QVEF  (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'surface albedo (IR, direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_IR )
             LOG_WARN_CONT('(A,F32.16)') 'surface albedo (IR, diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_IR )
             LOG_WARN_CONT('(A,F32.16)') 'surface albedo (NIR,direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_NIR)
             LOG_WARN_CONT('(A,F32.16)') 'surface albedo (NIR,diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_NIR)
             LOG_WARN_CONT('(A,F32.16)') 'surface albedo (VIS,direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_VIS)
             LOG_WARN_CONT('(A,F32.16)') 'surface albedo (VIS,diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_VIS)
             LOG_WARN_CONT('(A,F32.16)') 'latent heat                        [J/kg]     :', LH    (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'stomata registance                 [1/s]      :', Rb    (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'thermal conductivity / depth       [J/m2/s/K] :', TC_dZ (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'roughness length for momemtum      [m]        :', Z0M   (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'roughness length for heat          [m]        :', Z0H   (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'roughness length for vapor         [m]        :', Z0E   (i,j)
             LOG_WARN_CONT('(A,F32.16)') 'time step                          [s]        :', dt
             LOG_NEWLINE
             LOG_WARN_CONT('(A,F32.16)') 'friction velocity                  [m/s]      :', Ustar(i,j)
             LOG_WARN_CONT('(A,F32.16)') 'friction potential temperature     [K]        :', Tstar(i,j)
             LOG_WARN_CONT('(A,F32.16)') 'friction water vapor mass ratio    [kg/kg]    :', Qstar(i,j)
             LOG_WARN_CONT('(A,F32.16)') 'free convection velocity scale     [m/s]      :', Wstar(i,j)
             LOG_WARN_CONT('(A,F32.16)') 'd(friction velocity)               [m/s]      :', dUstar
             LOG_WARN_CONT('(A,F32.16)') 'd(friction potential temperature)  [K]        :', dTstar
             LOG_WARN_CONT('(A,F32.16)') 'd(friction water vapor mass ratio) [kg/kg]    :', dQstar
             LOG_WARN_CONT('(A,F32.16)') 'd(free convection velocity scale)  [m/s]      :', dWstar
             LOG_WARN_CONT('(A,F32.16)') 'next surface temperature           [K]        :', TMPS1(i,j)
#endif

             ! check NaN
             if ( .NOT. ( res > -1.0_RP .OR. res < 1.0_RP ) ) then ! must be NaN
#ifdef _OPENACC
                err_flag = .true.
#else
                LOG_ERROR("CPL_PHY_SFC_skin",*) 'NaN is detected for surface temperature. ', trim(model_name)
                LOG_ERROR_CONT('(A,I32)'   ) 'number of i                        [no unit]  :', i
                LOG_ERROR_CONT('(A,I32)'   ) 'number of j                        [no unit]  :', j
                LOG_ERROR_CONT('(A,I32)'   ) 'loop number                        [no unit]  :', n
                LOG_ERROR_CONT('(A,F32.16)') 'temperature                        [K]        :', TMPA  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'pressure                           [Pa]       :', PRSA  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'velocity w                         [m/s]      :', WA    (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'velocity u                         [m/s]      :', UA    (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'velocity v                         [m/s]      :', VA    (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'absolute velocity                  [m/s]      :', Uabs
                LOG_ERROR_CONT('(A,F32.16)') 'density                            [kg/m3]    :', RHOA  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'water vapor mass ratio             [kg/kg]    :', QVA   (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'cell center height                 [m]        :', Z1    (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'atmospheric mixing layer height    [m]        :', PBL   (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'pressure at the surface            [Pa]       :', PRSS  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'downward radiation (IR, direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_IR )
                LOG_ERROR_CONT('(A,F32.16)') 'downward radiation (IR, diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_IR )
                LOG_ERROR_CONT('(A,F32.16)') 'downward radiation (NIR,direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_NIR)
                LOG_ERROR_CONT('(A,F32.16)') 'downward radiation (NIR,diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_NIR)
                LOG_ERROR_CONT('(A,F32.16)') 'downward radiation (VIS,direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_VIS)
                LOG_ERROR_CONT('(A,F32.16)') 'downward radiation (VIS,diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_VIS)
                LOG_ERROR_CONT('(A,F32.16)') 'soil temperature                   [K]        :', TG    (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'soil water                         [kg/m2]    :', WSTR  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'surface temperature                [K]        :', TMPS  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'surface density                    [kg/m3]    :', RHOS  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'efficiency of evaporation          [1]        :', QVEF  (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'surface albedo (IR, direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_IR )
                LOG_ERROR_CONT('(A,F32.16)') 'surface albedo (IR, diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_IR )
                LOG_ERROR_CONT('(A,F32.16)') 'surface albedo (NIR,direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_NIR)
                LOG_ERROR_CONT('(A,F32.16)') 'surface albedo (NIR,diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_NIR)
                LOG_ERROR_CONT('(A,F32.16)') 'surface albedo (VIS,direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_VIS)
                LOG_ERROR_CONT('(A,F32.16)') 'surface albedo (VIS,diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_VIS)
                LOG_ERROR_CONT('(A,F32.16)') 'latent heat                        [J/kg]     :', LH    (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'stomata registance                 [1/s]      :', Rb    (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'thermal conductivity / depth       [J/m2/s/K] :', TC_dZ (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'roughness length for momemtum      [m]        :', Z0M   (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'roughness length for heat          [m]        :', Z0H   (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'roughness length for vapor         [m]        :', Z0E   (i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'time step                          [s]        :', dt
                LOG_ERROR_CONT('(A,F32.16)') 'friction velocity                  [m/s]      :', Ustar(i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'friction potential temperature     [K]        :', Tstar(i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'friction water vapor mass ratio    [kg/kg]    :', Qstar(i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'free convection velocity scale     [m/s]      :', Wstar(i,j)
                LOG_ERROR_CONT('(A,F32.16)') 'd(friction velocity)               [m/s]      :', dUstar
                LOG_ERROR_CONT('(A,F32.16)') 'd(friction potential temperature)  [K]        :', dTstar
                LOG_ERROR_CONT('(A,F32.16)') 'd(friction water vapor mass ratio) [kg/kg]    :', dQstar
                LOG_ERROR_CONT('(A,F32.16)') 'd(free convection velocity scale)  [m/s]      :', dWstar
                LOG_ERROR_CONT('(A,F32.16)') 'next surface temperature           [K]        :', TMPS1(i,j)
                call PRC_abort
#endif
             endif
          endif

          ! calculate surface flux
          TMPS(i,j) = TMPS1(i,j)

!          qdry = 1.0_RP - QVA(i,j)
 !         Rtot = qdry * Rdry + QVA(i,j) * Rvap

          call qsat( TMPS(i,j), RHOS(i,j), QVsat )
!          call qsat( TMPS(i,j), PRSS(i,j), qdry, QVsat )

          QVS(i,j) = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
                   + (        QVEF(i,j) ) * QVsat

          call BULKFLUX( TMPA(i,j), TMPS(i,j),                  & ! [IN]
                         PRSA(i,j), PRSS(i,j),                  & ! [IN]
                         QVA (i,j), QVS (i,j),                  & ! [IN]
                         Uabs, Z1(i,j), PBL(i,j),               & ! [IN]
                         Z0M(i,j), Z0H(i,j), Z0E(i,j),          & ! [IN]
                         Ustar(i,j), Tstar(i,j), Qstar(i,j),    & ! [OUT]
                         Wstar(i,j), RLmo(i,j), Ra,             & ! [OUT]
                         FracU10(i,j), FracT2(i,j), FracQ2(i,j) ) ! [OUT]

          if ( Uabs < EPS ) then
             ZMFLX(i,j) = 0.0_RP
             XMFLX(i,j) = 0.0_RP
             YMFLX(i,j) = 0.0_RP
          else
             MFLUX = - RHOS(i,j) * Ustar(i,j)**2
             ZMFLX(i,j) = MFLUX * WA(i,j) / Uabs
             XMFLX(i,j) = MFLUX * UA(i,j) / Uabs
             YMFLX(i,j) = MFLUX * VA(i,j) / Uabs
          end if
          SHFLX(i,j) = -RHOS(i,j) * Ustar(i,j) * Tstar(i,j) * CPdry
          QVFLX(i,j) = min( - RHOS(i,j) * Ustar(i,j) * Qstar(i,j) * Ra / ( Ra+Rb(i,j) ), WSTR(i,j)/real(dt,RP) )

          emis = ( 1.0_RP-ALBEDO(i,j,I_R_diffuse,I_R_IR) ) * STB * TMPS(i,j)**4

          LWD  = RFLXD(i,j,I_R_diffuse,I_R_IR)
          LWU  = RFLXD(i,j,I_R_diffuse,I_R_IR)  * ALBEDO(i,j,I_R_diffuse,I_R_IR) + emis
          SWD  = RFLXD(i,j,I_R_direct ,I_R_NIR) &
               + RFLXD(i,j,I_R_diffuse,I_R_NIR) &
               + RFLXD(i,j,I_R_direct ,I_R_VIS) &
               + RFLXD(i,j,I_R_diffuse,I_R_VIS)
          SWU  = RFLXD(i,j,I_R_direct ,I_R_NIR) * ALBEDO(i,j,I_R_direct ,I_R_NIR) &
               + RFLXD(i,j,I_R_diffuse,I_R_NIR) * ALBEDO(i,j,I_R_diffuse,I_R_NIR) &
               + RFLXD(i,j,I_R_direct ,I_R_VIS) * ALBEDO(i,j,I_R_direct ,I_R_VIS) &
               + RFLXD(i,j,I_R_diffuse,I_R_VIS) * ALBEDO(i,j,I_R_diffuse,I_R_VIS)

          GFLX(i,j) = TC_dZ(i,j) * ( TMPS(i,j) - TG(i,j) )

          LHFLX(i,j) = QVFLX(i,j) * LH(i,j)


          ! calculation for residual
          res = SWD - SWU + LWD - LWU - SHFLX(i,j) - LHFLX(i,j) - GFLX(i,j)

          ! put residual in ground heat flux
          GFLX(i,j) = GFLX(i,j) + res

       else ! not calculate surface flux
          ZMFLX(i,j) = UNDEF
          XMFLX(i,j) = UNDEF
          YMFLX(i,j) = UNDEF
          SHFLX(i,j) = UNDEF
          LHFLX(i,j) = UNDEF
          QVFLX(i,j) = UNDEF
          GFLX (i,j) = UNDEF
          Ustar(i,j) = UNDEF
          Tstar(i,j) = UNDEF
          Qstar(i,j) = UNDEF
          Wstar(i,j) = UNDEF
          RLmo (i,j) = UNDEF
          U10  (i,j) = UNDEF
          V10  (i,j) = UNDEF
          T2   (i,j) = UNDEF
          Q2   (i,j) = UNDEF
       endif
    enddo
    enddo
    !$acc end parallel

#ifdef _OPENACC
    if ( err_flag ) then
       LOG_ERROR("CPL_PHY_SFC_skin",*) 'NaN is detected for surface temperature. ', trim(model_name)
       call PRC_abort
    end if
#endif

    call BULKFLUX_diagnose_surface( IA, IS, IE, JA, JS, JE, &
                                    UA(:,:), VA(:,:),                      & ! (in)
                                    TMPA(:,:), QVA(:,:),                   & ! (in)
                                    TMPS(:,:), QVS(:,:),                   & ! (in)
                                    Z1(:,:), Z0M(:,:), Z0H(:,:), Z0E(:,:), & ! (in)
                                    U10(:,:), V10(:,:), T2(:,:), Q2(:,:),  & ! (out)
                                    mask = calc_flag(:,:),                 & ! (in)
                                    FracU10 = FracU10(:,:),                & ! (in)
                                    FracT2 = FracT2(:,:),                  & ! (in)
                                    FracQ2 = FracQ2(:,:)                   ) ! (in)

    !$acc end data

    return
  end subroutine CPL_PHY_SFC_skin

end module scale_cpl_phy_sfc_skin
