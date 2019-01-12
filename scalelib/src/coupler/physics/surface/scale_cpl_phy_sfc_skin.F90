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
  integer,  private :: CPL_PHY_SFC_SKIN_itr_max = 100 ! maximum iteration number

  real(RP), private :: CPL_PHY_SFC_SKIN_dTS_max = 5.0E-2_RP ! maximum delta surface temperature [K/s]
  real(RP), private :: CPL_PHY_SFC_SKIN_res_min = 1.0E+0_RP ! minimum value of residual
  real(RP), private :: CPL_PHY_SFC_SKIN_err_min = 1.0E-2_RP ! minimum value of error
  real(RP), private :: CPL_PHY_SFC_SKIN_dreslim = 1.0E+2_RP ! limiter of d(residual)

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
       CPL_PHY_SFC_SKIN_err_min, &
       CPL_PHY_SFC_SKIN_dreslim

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( initialized ) return

    LOG_NEWLINE
    LOG_INFO("CPL_PHY_SFC_SKIN_setup",*) 'Setup'

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
  subroutine CPL_PHY_SFC_skin( &
       IA, IS, IE,          &
       JA, JS, JE,          &
       TMPA, PRSA,          &
       WA, UA, VA,          &
       RHOA, QVA, LH,       &
       Z1, PBL,             &
       RHOS, PRSS,          &
       RFLXD,               &
       TG, QVEF,            &
       ALBEDO,              &
       Rb, TC_dZ,           &
       Z0M, Z0H, Z0E,       &
       calc_flag, dt,       &
       model_name,          &
       TMPS,                &
       ZMFLX, XMFLX, YMFLX, &
       SHFLX, QVFLX, GFLX,  &
       U10, V10, T2, Q2     )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       EPS   => CONST_EPS, &
       PRE00 => CONST_PRE00, &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       Rvap  => CONST_Rvap,  &
       STB   => CONST_STB
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_dens2qsat_all
!       qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
       BULKFLUX
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
    real(RP),         intent(in)    :: LH       (IA,JA)                     ! latent heat at the lowest atmospheric layer [J/kg]
    real(RP),         intent(in)    :: Z1       (IA,JA)                     ! cell center height at the lowest atmospheric layer [m]
    real(RP),         intent(in)    :: PBL      (IA,JA)                     ! the top of atmospheric mixing layer [m]
    real(RP),         intent(in)    :: RHOS     (IA,JA)                     ! density  at the surface [kg/m3]
    real(RP),         intent(in)    :: PRSS     (IA,JA)                     ! pressure at the surface [Pa]
    real(RP),         intent(in)    :: RFLXD    (IA,JA,N_RAD_DIR,N_RAD_RGN) ! downward radiation flux at the surface (direct/diffuse,IR/near-IR/VIS) [J/m2/s]
    real(RP),         intent(in)    :: TG       (IA,JA)                     ! subsurface temperature [K]
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
    real(RP),         intent(out)   :: QVFLX    (IA,JA)                     ! water vapor     flux at the surface [kg/m2/s]
    real(RP),         intent(out)   :: GFLX     (IA,JA)                     ! subsurface heat flux at the surface [J/m2/s]
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

    real(RP) :: emis          ! surface longwave emission                 [J/m2/s]
    real(RP) :: LWD           ! surface downward longwave  radiation flux [J/m2/s]
    real(RP) :: LWU           ! surface upward   longwave  radiation flux [J/m2/s]
    real(RP) :: SWD           ! surface downward shortwave radiation flux [J/m2/s]
    real(RP) :: SWU           ! surface upward   shortwave radiation flux [J/m2/s]
    real(RP) :: res           ! residual

    real(RP) :: dres          ! d(residual)/dTMPS
    real(RP) :: oldres        ! residual in previous step
    real(RP) :: redf          ! reduced factor

    real(RP) :: Ustar, dUstar ! friction velocity               [m/s]
    real(RP) :: Tstar, dTstar ! friction potential temperature  [K]
    real(RP) :: Qstar, dQstar ! friction water vapor mass ratio [kg/kg]
    real(RP) :: Wstar, dWstar ! free convection velocity scale  [m/s]
    real(RP) :: Uabs,  dUabs  ! modified absolute velocity      [m/s]
    real(RP) :: Ra,    dRa    ! Aerodynamic resistance (=1/Ce)  [1/s]

    real(RP) :: QVsat, dQVsat ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS, dQVS     ! water vapor mixing ratio at surface            [kg/kg]
    real(RP) :: Rtot          ! total gas constant
    real(RP) :: qdry          ! dry air mass ratio [kg/kg]

    real(RP) :: FracU10       ! calculation parameter for U10 [1]
    real(RP) :: FracT2        ! calculation parameter for T2  [1]
    real(RP) :: FracQ2        ! calculation parameter for Q2  [1]

    real(RP) :: MFLUX

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'coupler / physics / surface / SKIN'

    ! copy surfce temperature for iteration
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       TMPS1(i,j) = TMPS(i,j)
    enddo
    enddo

    ! update surface temperature
    !$omp parallel do &
#ifndef __GFORTRAN__
    !$omp default(none) &
    !$omp shared(IS,IE,JS,JE,EPS,Rdry,CPdry,PRC_myrank,IO_FID_LOG,IO_L,model_name,bulkflux, &
    !$omp        CPL_PHY_SFC_SKIN_itr_max,CPL_PHY_SFC_SKIN_dTS_max,CPL_PHY_SFC_SKIN_dreslim,CPL_PHY_SFC_SKIN_err_min, CPL_PHY_SFC_SKIN_res_min, &
    !$omp        calc_flag,dt,QVA,TMPA,PRSA,RHOA,WA,UA,VA,LH,Z1,PBL, &
    !$omp        TG,PRSS,RHOS,TMPS1,QVEF,Z0M,Z0H,Z0E,Rb,TC_dZ,ALBEDO,RFLXD, &
    !$omp        TMPS,ZMFLX,XMFLX,YMFLX,SHFLX,QVFLX,GFLX,U10,V10,T2,Q2) &
#else
    !$omp default(shared) &
#endif
    !$omp private(qdry,Rtot,redf,res,emis,LWD,LWU,SWD,SWU,dres,oldres,QVS,dQVS, &
    !$omp         QVsat,dQVsat,Ustar,dUstar,Tstar,dTstar,Qstar,dQstar,Wstar,dWstar,&
    !$omp         Uabs,dUabs,Ra,dRa,FracU10,FracT2,FracQ2,MFLUX)
    do j = JS, JE
    do i = IS, IE
       if ( calc_flag(i,j) ) then

!          qdry = 1.0_RP - QVA(i,j)
!          Rtot = qdry * Rdry + QVA(i,j) * Rvap

          redf   = 1.0_RP
          oldres = huge(0.0_RP)

          ! modified Newton-Raphson method (Tomita 2009)
          do n = 1, CPL_PHY_SFC_SKIN_itr_max

             call qsat( TMPS1(i,j),      RHOS(i,j), QVsat  )
             call qsat( TMPS1(i,j)+dTS0, RHOS(i,j), dQVsat )
!             call qsat( TMPS1(i,j),      PRSS(i,j), qdry, QVsat  )
!             call qsat( TMPS1(i,j)+dTS0, PRSS(i,j), qdry, dQVsat )

             QVS  = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
                  + (        QVEF(i,j) ) * QVsat
             dQVS = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
                  + (        QVEF(i,j) ) * dQVsat

             Uabs = sqrt( WA(i,j)**2 + UA(i,j)**2 + VA(i,j)**2 )

             call BULKFLUX( Ustar,           & ! [OUT]
                            Tstar,           & ! [OUT]
                            Qstar,           & ! [OUT]
                            Wstar,           & ! [OUT]
                            Ra,              & ! [OUT]
                            FracU10,         & ! [OUT] ! not used
                            FracT2,          & ! [OUT] ! not used
                            FracQ2,          & ! [OUT] ! not used
                            TMPA (i,j),      & ! [IN]
                            TMPS1(i,j),      & ! [IN]
                            PRSA (i,j),      & ! [IN]
                            PRSS (i,j),      & ! [IN]
                            QVA  (i,j),      & ! [IN]
                            QVS,             & ! [IN]
                            Uabs,            & ! [IN]
                            Z1   (i,j),      & ! [IN]
                            PBL  (i,j),      & ! [IN]
                            Z0M  (i,j),      & ! [IN]
                            Z0H  (i,j),      & ! [IN]
                            Z0E  (i,j)       ) ! [IN]

             call BULKFLUX( dUstar,          & ! [OUT]
                            dTstar,          & ! [OUT]
                            dQstar,          & ! [OUT]
                            dWstar,          & ! [OUT]
                            dRa,             & ! [OUT] ! not used
                            FracU10,         & ! [OUT] ! not used
                            FracT2,          & ! [OUT] ! not used
                            FracQ2,          & ! [OUT] ! not used
                            TMPA (i,j),      & ! [IN]
                            TMPS1(i,j)+dTS0, & ! [IN]
                            PRSA (i,j),      & ! [IN]
                            PRSS (i,j),      & ! [IN]
                            QVA  (i,j),      & ! [IN]
                            dQVS,            & ! [IN]
                            Uabs,            & ! [IN]
                            Z1   (i,j),      & ! [IN]
                            PBL  (i,j),      & ! [IN]
                            Z0M  (i,j),      & ! [IN]
                            Z0H  (i,j),      & ! [IN]
                            Z0E  (i,j)       ) ! [IN]

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
             res = SWD - SWU + LWD - LWU                                     &
                 + CPdry   * RHOS(i,j) * Ustar * Tstar                       &
                 + LH(i,j) * RHOS(i,j) * Ustar * Qstar * Ra / ( Ra+Rb(i,j) ) &
                 - TC_dZ(i,j) * ( TMPS1(i,j) - TG(i,j) )

             ! calculation for d(residual)/dTMPS
             dres = -4.0_RP * emis / TMPS1(i,j)                                                                           &
                  + CPdry   * RHOS(i,j) * ( Ustar*(dTstar-Tstar)/dTS0 + Tstar*(dUstar-Ustar)/dTS0 )                       &
                  + LH(i,j) * RHOS(i,j) * ( Ustar*(dQstar-Qstar)/dTS0 + Qstar*(dUstar-Ustar)/dTS0 ) * Ra / ( Ra+Rb(i,j) ) &
                  - TC_dZ(i,j)

             ! convergence test with residual and error levels
             if (      abs(res     ) < CPL_PHY_SFC_SKIN_res_min &
                  .OR. abs(res/dres) < CPL_PHY_SFC_SKIN_err_min ) then
                exit
             endif

             ! stop iteration to prevent numerical error
             if ( abs(dres) * CPL_PHY_SFC_SKIN_dreslim < abs(res) ) then
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
             TMPS1(i,j) = TMPS1(i,j) - redf * res / dres

             ! save residual in this step
             oldres = res
          enddo

          ! update surface temperature with limitation
          TMPS1(i,j) = min( max( TMPS1(i,j),                                                 &
                                 TMPS (i,j) - CPL_PHY_SFC_SKIN_dTS_max * real(dt,kind=RP) ), &
                                 TMPS (i,j) + CPL_PHY_SFC_SKIN_dTS_max * real(dt,kind=RP) )

          if ( n > CPL_PHY_SFC_SKIN_itr_max ) then
             ! surface temperature was not converged
             LOG_WARN("CPL_PHY_SFC_skin",*) 'surface tempearture was not converged. ', trim(model_name)
             LOG_NEWLINE
             LOG_INFO_CONT('(A,I32)'   ) 'PRC_myrank                         [no unit]  :', PRC_myrank
             LOG_INFO_CONT('(A,I32)'   ) 'number of i                        [no unit]  :', i
             LOG_INFO_CONT('(A,I32)'   ) 'number of j                        [no unit]  :', j
             LOG_NEWLINE
             LOG_INFO_CONT('(A,I32)'   ) 'loop number                        [no unit]  :', n
             LOG_INFO_CONT('(A,F32.16)') 'Residual                           [J/m2/s]   :', res
             LOG_INFO_CONT('(A,F32.16)') 'delta Residual                     [J/m2/s]   :', dres
             LOG_NEWLINE
             LOG_INFO_CONT('(A,F32.16)') 'temperature                        [K]        :', TMPA  (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'pressure                           [Pa]       :', PRSA  (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'velocity w                         [m/s]      :', WA    (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'velocity u                         [m/s]      :', UA    (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'velocity v                         [m/s]      :', VA    (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'absolute velocity                  [m/s]      :', Uabs
             LOG_INFO_CONT('(A,F32.16)') 'density                            [kg/m3]    :', RHOA  (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'water vapor mass ratio             [kg/kg]    :', QVA   (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'cell center height                 [m]        :', Z1    (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'atmospheric mixing layer height    [m]        :', PBL   (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'pressure at the surface            [Pa]       :', PRSS  (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'downward radiation (IR, direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_IR )
             LOG_INFO_CONT('(A,F32.16)') 'downward radiation (IR, diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_IR )
             LOG_INFO_CONT('(A,F32.16)') 'downward radiation (NIR,direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_NIR)
             LOG_INFO_CONT('(A,F32.16)') 'downward radiation (NIR,diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_NIR)
             LOG_INFO_CONT('(A,F32.16)') 'downward radiation (VIS,direct )   [J/m2/s]   :', RFLXD (i,j,I_R_direct ,I_R_VIS)
             LOG_INFO_CONT('(A,F32.16)') 'downward radiation (VIS,diffuse)   [J/m2/s]   :', RFLXD (i,j,I_R_diffuse,I_R_VIS)
             LOG_NEWLINE
             LOG_INFO_CONT('(A,F32.16)') 'soil temperature                   [K]        :', TG    (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'surface temperature                [K]        :', TMPS  (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'efficiency of evaporation          [1]        :', QVEF  (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'surface albedo (IR, direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_IR )
             LOG_INFO_CONT('(A,F32.16)') 'surface albedo (IR, diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_IR )
             LOG_INFO_CONT('(A,F32.16)') 'surface albedo (NIR,direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_NIR)
             LOG_INFO_CONT('(A,F32.16)') 'surface albedo (NIR,diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_NIR)
             LOG_INFO_CONT('(A,F32.16)') 'surface albedo (VIS,direct )       [1]        :', ALBEDO(i,j,I_R_direct ,I_R_VIS)
             LOG_INFO_CONT('(A,F32.16)') 'surface albedo (VIS,diffuse)       [1]        :', ALBEDO(i,j,I_R_diffuse,I_R_VIS)
             LOG_INFO_CONT('(A,F32.16)') 'thermal conductivity / depth       [J/m2/s/K] :', TC_dZ (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'roughness length for momemtum      [m]        :', Z0M   (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'roughness length for heat          [m]        :', Z0H   (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'roughness length for vapor         [m]        :', Z0E   (i,j)
             LOG_NEWLINE
             LOG_INFO_CONT('(A,F32.16)') 'latent heat                        [J/kg]     :', LH   (i,j)
             LOG_INFO_CONT('(A,F32.16)') 'friction velocity                  [m/s]      :', Ustar
             LOG_INFO_CONT('(A,F32.16)') 'friction potential temperature     [K]        :', Tstar
             LOG_INFO_CONT('(A,F32.16)') 'friction water vapor mass ratio    [kg/kg]    :', Qstar
             LOG_INFO_CONT('(A,F32.16)') 'free convection velocity scale     [m/s]      :', Wstar
             LOG_INFO_CONT('(A,F32.16)') 'd(friction velocity)               [m/s]      :', dUstar
             LOG_INFO_CONT('(A,F32.16)') 'd(friction potential temperature)  [K]        :', dTstar
             LOG_INFO_CONT('(A,F32.16)') 'd(friction water vapor mass ratio) [kg/kg]    :', dQstar
             LOG_INFO_CONT('(A,F32.16)') 'd(free convection velocity scale)  [m/s]      :', dWstar
             LOG_INFO_CONT('(A,F32.16)') 'next surface temperature           [K]        :', TMPS1(i,j)

             ! check NaN
             if ( .NOT. ( res > -1.0_RP .OR. res < 1.0_RP ) ) then ! must be NaN
                LOG_ERROR("CPL_PHY_SFC_skin",*) 'NaN is detected for surface temperature. ', trim(model_name)
                call PRC_abort
             endif
          endif

          ! calculate surface flux
          TMPS(i,j) = TMPS1(i,j)

!          qdry = 1.0_RP - QVA(i,j)
 !         Rtot = qdry * Rdry + QVA(i,j) * Rvap

          call qsat( TMPS(i,j), RHOS(i,j), QVsat )
!          call qsat( TMPS(i,j), PRSS(i,j), qdry, QVsat )

          QVS = ( 1.0_RP-QVEF(i,j) ) * QVA(i,j) &
              + (        QVEF(i,j) ) * QVsat

          call BULKFLUX( Ustar,     & ! [OUT]
                         Tstar,     & ! [OUT]
                         Qstar,     & ! [OUT]
                         Wstar,     & ! [OUT]
                         Ra,        & ! [OUT]
                         FracU10,   & ! [OUT]
                         FracT2,    & ! [OUT]
                         FracQ2,    & ! [OUT]
                         TMPA(i,j), & ! [IN]
                         TMPS(i,j), & ! [IN]
                         PRSA(i,j), & ! [IN]
                         PRSS(i,j), & ! [IN]
                         QVA (i,j), & ! [IN]
                         QVS,       & ! [IN]
                         Uabs,      & ! [IN]
                         Z1  (i,j), & ! [IN]
                         PBL (i,j), & ! [IN]
                         Z0M (i,j), & ! [IN]
                         Z0H (i,j), & ! [IN]
                         Z0E (i,j)  ) ! [IN]

          if ( Uabs < EPS ) then
             ZMFLX(i,j) = 0.0_RP
             XMFLX(i,j) = 0.0_RP
             YMFLX(i,j) = 0.0_RP
          else
             MFLUX = - RHOS(i,j) * Ustar**2
             ZMFLX(i,j) = MFLUX * WA(i,j) / Uabs
             XMFLX(i,j) = MFLUX * UA(i,j) / Uabs
             YMFLX(i,j) = MFLUX * VA(i,j) / Uabs
          end if
          SHFLX(i,j) = -RHOS(i,j) * Ustar * Tstar * CPdry
          QVFLX(i,j) = -RHOS(i,j) * Ustar * Qstar * Ra / ( Ra+Rb(i,j) )

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

          GFLX(i,j) = -TC_dZ(i,j) * ( TMPS(i,j) - TG(i,j) )

          ! calculation for residual
          res = SWD - SWU + LWD - LWU - SHFLX(i,j) - QVFLX(i,j) * LH(i,j) + GFLX(i,j)

          ! put residual in ground heat flux
          GFLX(i,j) = GFLX(i,j) - res

          ! diagnostic variables considering unstable/stable state
          !U10(i,j) = FracU10 * UA(i,j)
          !V10(i,j) = FracU10 * VA(i,j)
          !T2 (i,j) = ( 1.0_RP - FracT2 ) * TMPS(i,j) + FracT2 * TMPA(i,j)
          !Q2 (i,j) = ( 1.0_RP - FracQ2 ) * QVS       + FracQ2 * QVA (i,j)

          ! diagnostic variables for neutral state
          U10(i,j) = UA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
          V10(i,j) = VA  (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
          T2 (i,j) = TMPS(i,j) + ( TMPA(i,j) - TMPS(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                           / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
          Q2 (i,j) = QVS       + (  QVA(i,j) - QVS       ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                           / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )

       else ! not calculate surface flux
          ZMFLX(i,j) = 0.0_RP
          XMFLX(i,j) = 0.0_RP
          YMFLX(i,j) = 0.0_RP
          SHFLX(i,j) = 0.0_RP
          QVFLX(i,j) = 0.0_RP
          GFLX (i,j) = 0.0_RP
          U10  (i,j) = 0.0_RP
          V10  (i,j) = 0.0_RP
          T2   (i,j) = 0.0_RP
          Q2   (i,j) = 0.0_RP
       endif
    enddo
    enddo

    return
  end subroutine CPL_PHY_SFC_skin

end module scale_cpl_phy_sfc_skin
