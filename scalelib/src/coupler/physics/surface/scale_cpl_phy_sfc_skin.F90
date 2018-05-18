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
  !-----------------------------------------------------------------------------
  integer,  private :: CPL_PHY_SFC_SKIN_itr_max = 100 ! maximum iteration number

  real(RP), private :: CPL_PHY_SFC_SKIN_dTS_max = 5.0E-2_RP ! maximum delta surface temperature [K/s]
  real(RP), private :: CPL_PHY_SFC_SKIN_res_min = 1.0E+0_RP ! minimum value of residual
  real(RP), private :: CPL_PHY_SFC_SKIN_err_min = 1.0E-2_RP ! minimum value of error
  real(RP), private :: CPL_PHY_SFC_SKIN_dreslim = 1.0E+2_RP ! limiter of d(residual)

  logical, private :: initialized = .false.

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
       IA, IS, IE, JA, JS, JE, &
       TMPA, PRSA,          &
       WA, UA, VA,          &
       RHOA, QVA, LHV,      &
       Z1, PBL,             &
       RHOS, PRSS,          &
       RFLXD,               &
       TG, QVEF,            &
       ALBEDO,              &
       DZG,                 &
       Rb, TCS,             &
       Z0M, Z0H, Z0E,       &
       fact_area, dt,       &
       model_name,          &
       TMPS,                &
       ZMFLX, XMFLX, YMFLX, &
       SHFLX, LHFLX, GHFLX, &
       U10, V10, T2, Q2     )
    use scale_prc, only: &
      PRC_myrank,  &
      PRC_abort
    use scale_const, only: &
      PRE00 => CONST_PRE00, &
      Rdry  => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      Rvap  => CONST_Rvap,  &
      STB   => CONST_STB
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
      BULKFLUX
    implicit none
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: TMPA(IA,JA) ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: PRSA(IA,JA) ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: WA  (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: UA  (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: VA  (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: RHOA(IA,JA) ! density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: QVA (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: LHV (IA,JA) ! latent heat of vaporization [J/kg]
    real(RP), intent(in) :: Z1  (IA,JA) ! cell center height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL (IA,JA) ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: RHOS(IA,JA) ! density  at the surface [kg/m3]
    real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: RFLXD (IA,JA,N_RAD_DIR,N_RAD_RGN) ! downward radiation flux at the surface (direct/diffuse,IR/near-IR/VIS) [J/m2/s]

    real(RP), intent(in) :: TG    (IA,JA) ! soil temperature [K]
    real(RP), intent(in) :: QVEF  (IA,JA) ! efficiency of evaporation (0-1)
    real(RP), intent(in) :: ALBEDO(IA,JA,N_RAD_DIR,N_RAD_RGN) ! surface albedo (direct/diffuse,IR/near-IR/VIS) (0-1)
    real(RP), intent(in) :: DZG   (IA,JA) ! soil depth [m]
    real(RP), intent(in) :: Rb    (IA,JA) ! stomata resistance [1/s]
    real(RP), intent(in) :: TCS   (IA,JA) ! thermal conductivity for soil [J/m/K/s]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
    real(RP), intent(in) :: fact_area(IA,JA)
    real(DP), intent(in) :: dt            ! delta time

    character(len=*), intent(in) :: model_name

    real(RP), intent(inout) :: TMPS(IA,JA) ! surface temperature [K]

    real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m/s2]
    real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: GHFLX(IA,JA) ! ground heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]


    ! parameters
    real(RP), parameter :: dTS0     = 1.0E-4_RP ! delta surface temp.

    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0E+0_RP ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5E+0_RP ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1E+0_RP ! factor b in Tomita (2009)

    ! works
    real(RP) :: TMPS1(IA,JA)

    real(RP) :: emis   ! surface longwave emission                 [J/m2/s]
    real(RP) :: LWD    ! surface downward longwave  radiation flux [J/m2/s]
    real(RP) :: LWU    ! surface upward   longwave  radiation flux [J/m2/s]
    real(RP) :: SWD    ! surface downward shortwave radiation flux [J/m2/s]
    real(RP) :: SWU    ! surface upward   shortwave radiation flux [J/m2/s]
    real(RP) :: res    ! residual

    real(RP) :: dres   ! d(residual)/dTMPS
    real(RP) :: oldres ! residual in previous step
    real(RP) :: redf   ! reduced factor

    real(RP) :: Ustar, dUstar ! friction velocity [m]
    real(RP) :: Tstar, dTstar ! friction potential temperature [K]
    real(RP) :: Qstar, dQstar ! friction water vapor mass ratio [kg/kg]
    real(RP) :: Uabs,  dUabs  ! modified absolute velocity [m/s]
    real(RP) :: Ra,    dRa    ! Aerodynamic resistance (=1/Ce) [1/s]

    real(RP) :: QVsat, dQVsat ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS, dQVS     ! water vapor mixing ratio at surface [kg/kg]
    real(RP) :: Rtot          ! total gas constant
    real(RP) :: qdry          ! dry air mass ratio [kg/kg]

    real(RP) :: FracU10 ! calculation parameter for U10 [-]
    real(RP) :: FracT2  ! calculation parameter for T2 [-]
    real(RP) :: FracQ2  ! calculation parameter for Q2 [-]

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'coupler / physics / surface / SKIN'

    ! copy surfce temperature for iteration
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
      TMPS1(i,j) = TMPS(i,j)
    end do
    end do

    ! update surface temperature
    !$omp parallel do default(none) &
    !$omp private(qdry,Rtot,redf,res,emis,LWD,LWU,SWD,SWU,dres,oldres,QVS,dQVS, &
    !$omp         QVsat,dQVsat,Ustar,dUstar,Tstar,dTstar,Qstar,dQstar,Uabs,dUabs,Ra,dRa, &
    !$omp         FracU10,FracT2,FracQ2) &
    !$omp shared(IS,IE,JS,JE,Rdry,CPdry,PRC_myrank,IO_FID_LOG,IO_L,model_name, &
    !$omp        bulkflux, &
    !$omp        CPL_PHY_SFC_SKIN_itr_max,CPL_PHY_SFC_SKIN_dTS_max,CPL_PHY_SFC_SKIN_dreslim,CPL_PHY_SFC_SKIN_err_min, CPL_PHY_SFC_SKIN_res_min, &
    !$omp        fact_area,DZG,dt,QVA,TMPA,PRSA,RHOA,WA,UA,VA,LHV,Z1,PBL, &
    !$omp        TG,PRSS,RHOS,TMPS1,QVEF,Z0M,Z0H,Z0E,Rb,TCS,ALBEDO,RFLXD, &
    !$omp        TMPS,ZMFLX,XMFLX,YMFLX,SHFLX,LHFLX,GHFLX,U10,V10,T2,Q2)
    do j = JS, JE
    do i = IS, IE

      if( fact_area(i,j) > 0.0_RP ) then

        qdry = 1.0_RP - QVA(i,j)
        Rtot = qdry * Rdry + QVA(i,j) * Rvap

        redf   = 1.0_RP
        oldres = huge(0.0_RP)

        ! modified Newton-Raphson method (Tomita 2009)
        do n = 1, CPL_PHY_SFC_SKIN_itr_max

          call qsat( TMPS1(i,j),      PRSS(i,j), qdry, & ! [IN]
                     QVsat                             ) ! [OUT]
          call qsat( TMPS1(i,j)+dTS0, PRSS(i,j), qdry, & ! [IN]
                     dQVsat                            ) ! [OUT]

          QVS  = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * QVsat
          dQVS = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * dQVsat

          call BULKFLUX( &
              Ustar,     & ! [OUT]
              Tstar,     & ! [OUT]
              Qstar,     & ! [OUT]
              Uabs,      & ! [OUT]
              Ra,        & ! [OUT]
              FracU10,   & ! [OUT] ! not used
              FracT2,    & ! [OUT] ! not used
              FracQ2,    & ! [OUT] ! not used
              TMPA(i,j), & ! [IN]
              TMPS1(i,j), & ! [IN]
              PRSA(i,j), & ! [IN]
              PRSS(i,j), & ! [IN]
              QVA (i,j), & ! [IN]
              QVS,       & ! [IN]
              UA  (i,j), & ! [IN]
              VA  (i,j), & ! [IN]
              Z1  (i,j), & ! [IN]
              PBL (i,j), & ! [IN]
              Z0M (i,j), & ! [IN]
              Z0H (i,j), & ! [IN]
              Z0E (i,j)  ) ! [IN]

          call BULKFLUX( &
              dUstar,         & ! [OUT]
              dTstar,         & ! [OUT]
              dQstar,         & ! [OUT]
              dUabs,          & ! [OUT]
              dRa,            & ! [OUT] ! not used
              FracU10,        & ! [OUT] ! not used
              FracT2,         & ! [OUT] ! not used
              FracQ2,         & ! [OUT] ! not used
              TMPA(i,j),      & ! [IN]
              TMPS1(i,j)+dTS0, & ! [IN]
              PRSA(i,j),      & ! [IN]
              PRSS(i,j),      & ! [IN]
              QVA (i,j),      & ! [IN]
              dQVS,           & ! [IN]
              UA  (i,j),      & ! [IN]
              VA  (i,j),      & ! [IN]
              Z1  (i,j),      & ! [IN]
              PBL (i,j),      & ! [IN]
              Z0M (i,j),      & ! [IN]
              Z0H (i,j),      & ! [IN]
              Z0E (i,j)       ) ! [IN]

           emis = ( 1.0_RP - ALBEDO(i,j,I_R_diffuse,I_R_IR) ) * STB * TMPS1(i,j)**4

           LWD  = RFLXD(i,j,I_R_diffuse,I_R_IR)
           LWU  = RFLXD(i,j,I_R_diffuse,I_R_IR) * ALBEDO(i,j,I_R_diffuse,I_R_IR) + emis
           SWD  = RFLXD(i,j,I_R_direct ,I_R_NIR) &
                + RFLXD(i,j,I_R_diffuse,I_R_NIR) &
                + RFLXD(i,j,I_R_direct ,I_R_VIS) &
                + RFLXD(i,j,I_R_diffuse,I_R_VIS)
           SWU  = RFLXD(i,j,I_R_direct ,I_R_NIR) * ALBEDO(i,j,I_R_direct ,I_R_NIR) &
                + RFLXD(i,j,I_R_diffuse,I_R_NIR) * ALBEDO(i,j,I_R_diffuse,I_R_NIR) &
                + RFLXD(i,j,I_R_direct ,I_R_VIS) * ALBEDO(i,j,I_R_direct ,I_R_VIS) &
                + RFLXD(i,j,I_R_diffuse,I_R_VIS) * ALBEDO(i,j,I_R_diffuse,I_R_VIS)

          ! calculation for residual
          res = SWD - SWU + LWD - LWU &
              + CPdry    * RHOS(i,j) * Ustar * Tstar &
              + LHV(i,j) * RHOS(i,j) * Ustar * Qstar * Ra / ( Ra + Rb(i,j) ) &
              - 2.0_RP * TCS(i,j) * ( TMPS1(i,j) - TG(i,j) ) / DZG(i,j)

          ! calculation for d(residual)/dTMPS
          dres = -4.0_RP * emis / TMPS1(i,j) &
               + CPdry    * RHOS(i,j) * ( (dUstar-Ustar)/dTS0 * Tstar + Ustar * (dTstar-Tstar)/dTS0 ) &
               + LHV(i,j) * RHOS(i,j) * ( (dUstar-Ustar)/dTS0 * Qstar + Ustar * (dQstar-Qstar)/dTS0 ) &
               * Ra / ( Ra + Rb(i,j) ) &
               - 2.0_RP * TCS(i,j) / DZG(i,j)

          ! convergence test with residual and error levels
          if( abs( res      ) < CPL_PHY_SFC_SKIN_res_min .or. &
              abs( res/dres ) < CPL_PHY_SFC_SKIN_err_min      ) then
            exit
          end if

          ! stop iteration to prevent numerical error
          if( abs(dres) * CPL_PHY_SFC_SKIN_dreslim < abs(res) ) then
            exit
          end if

          ! calculate reduced factor
          if( dres < 0.0_RP ) then
            if( abs(res) > abs(oldres) ) then
              redf = max( TFa*abs(redf), redf_min )
            else
              redf = min( TFb*abs(redf), redf_max )
            end if
          else
            redf = -1.0_RP
          end if

          ! estimate next surface temperature
          TMPS1(i,j) = TMPS1(i,j) - redf * res / dres

          ! save residual in this step
          oldres = res

        end do

        ! update surface temperature with limitation
        TMPS1(i,j) = min( max( TMPS1(i,j), &
                               TMPS(i,j) - CPL_PHY_SFC_SKIN_dTS_max * real( dt, kind=RP ) ), &
                               TMPS(i,j) + CPL_PHY_SFC_SKIN_dTS_max * real( dt, kind=RP ) )

        if( n > CPL_PHY_SFC_SKIN_itr_max ) then
          ! surface temperature was not converged
          LOG_WARN("CPL_PHY_SFC_skin",'(A)'       ) 'surface tempearture was not converged. ', trim(model_name)
          LOG_INFO_CONT('(A)'       ) ''
          LOG_INFO_CONT('(A,I32)'   ) 'DEBUG --- PRC_myrank                         [no unit] :', PRC_myrank
          LOG_INFO_CONT('(A,I32)'   ) 'DEBUG --- number of i                        [no unit] :', i
          LOG_INFO_CONT('(A,I32)'   ) 'DEBUG --- number of j                        [no unit] :', j
          LOG_INFO_CONT('(A)'       ) ''
          LOG_INFO_CONT('(A,I32)'   ) 'DEBUG --- loop number                        [no unit] :', n
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- Residual                           [J/m2/s]  :', res
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- delta Residual                     [J/m2/s]  :', dres
          LOG_INFO_CONT('(A,F32.16)') ''
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- temperature                        [K]       :', TMPA  (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- pressure                           [Pa]      :', PRSA  (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- velocity w                         [m/s]     :', WA    (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- velocity u                         [m/s]     :', UA    (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- velocity v                         [m/s]     :', VA    (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- density                            [kg/m3]   :', RHOA  (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- water vapor mass ratio             [kg/kg]   :', QVA   (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- cell center height                 [m]       :', Z1    (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- atmospheric mixing layer height    [m]       :', PBL   (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- pressure at the surface            [Pa]      :', PRSS  (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- downward radiation (IR, direct )   [J/m2/s]  :', RFLXD (i,j,I_R_direct ,I_R_IR )
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- downward radiation (IR, diffuse)   [J/m2/s]  :', RFLXD (i,j,I_R_diffuse,I_R_IR )
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- downward radiation (NIR,direct )   [J/m2/s]  :', RFLXD (i,j,I_R_direct ,I_R_NIR)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- downward radiation (NIR,diffuse)   [J/m2/s]  :', RFLXD (i,j,I_R_diffuse,I_R_NIR)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- downward radiation (VIS,direct )   [J/m2/s]  :', RFLXD (i,j,I_R_direct ,I_R_VIS)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- downward radiation (VIS,diffuse)   [J/m2/s]  :', RFLXD (i,j,I_R_diffuse,I_R_VIS)
          LOG_INFO_CONT('(A)'       ) ''
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- soil temperature                   [K]       :', TG    (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- surface temperature                [K]       :', TMPS  (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- efficiency of evaporation          [1]       :', QVEF  (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- surface albedo (IR, direct )       [1]       :', ALBEDO(i,j,I_R_direct ,I_R_IR )
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- surface albedo (IR, diffuse)       [1]       :', ALBEDO(i,j,I_R_diffuse,I_R_IR )
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- surface albedo (NIR,direct )       [1]       :', ALBEDO(i,j,I_R_direct ,I_R_NIR)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- surface albedo (NIR,diffuse)       [1]       :', ALBEDO(i,j,I_R_diffuse,I_R_NIR)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- surface albedo (VIS,direct )       [1]       :', ALBEDO(i,j,I_R_direct ,I_R_VIS)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- surface albedo (VIS,diffuse)       [1]       :', ALBEDO(i,j,I_R_diffuse,I_R_VIS)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- soil depth                         [m]       :', DZG   (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- thermal conductivity for soil      [J/m/K/s] :', TCS   (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- roughness length for momemtum      [m]       :', Z0M   (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- roughness length for heat          [m]       :', Z0H   (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- roughness length for vapor         [m]       :', Z0E   (i,j)
          LOG_INFO_CONT('(A)'       ) ''
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- latent heat                        [J/kg]    :', LHV   (i,j)
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- friction velocity                  [m]       :', Ustar
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- friction potential temperature     [K]       :', Tstar
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- friction water vapor mass ratio    [kg/kg]   :', Qstar
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- d(friction velocity)               [m]       :', dUstar
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- d(friction potential temperature)  [K]       :', dTstar
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- d(friction water vapor mass ratio) [kg/kg]   :', dQstar
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- modified absolute velocity         [m/s]     :', Uabs
          LOG_INFO_CONT('(A,F32.16)') 'DEBUG --- next surface temperature           [K]       :', TMPS1(i,j)

          ! check NaN
          if ( .NOT. ( res > -1.0_RP .OR. res < 1.0_RP ) ) then ! must be NaN
             LOG_ERROR("CPL_PHY_SFC_skin",*) 'NaN is detected for surface temperature. ', trim(model_name)
             LOG_ERROR_CONT(*) ''
             LOG_ERROR_CONT(*) 'DEBUG --- PRC_myrank                                   :', PRC_myrank
             LOG_ERROR_CONT(*) 'DEBUG --- number of i                                  :', i
             LOG_ERROR_CONT(*) 'DEBUG --- number of j                                  :', j
             LOG_ERROR_CONT(*) ''
             LOG_ERROR_CONT(*) 'DEBUG --- Residual                           [J/m2/s]  :', res
             LOG_ERROR_CONT(*) 'DEBUG --- delta Residual                     [J/m2/s]  :', dres
             LOG_ERROR_CONT(*) ''
             LOG_ERROR_CONT(*) 'DEBUG --- temperature                        [K]       :', TMPA  (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- pressure                           [Pa]      :', PRSA  (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- velocity w                         [m/s]     :', WA    (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- velocity u                         [m/s]     :', UA    (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- velocity v                         [m/s]     :', VA    (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- density                            [kg/m3]   :', RHOA  (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- water vapor mass ratio             [kg/kg]   :', QVA   (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- cell center height                 [m]       :', Z1    (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- atmospheric mixing layer height    [m]       :', PBL   (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- pressure at the surface            [Pa]      :', PRSS  (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- downward radiation (IR, direct )   [J/m2/s]  :', RFLXD (i,j,I_R_direct ,I_R_IR )
             LOG_ERROR_CONT(*) 'DEBUG --- downward radiation (IR, diffuse)   [J/m2/s]  :', RFLXD (i,j,I_R_diffuse,I_R_IR )
             LOG_ERROR_CONT(*) 'DEBUG --- downward radiation (NIR,direct )   [J/m2/s]  :', RFLXD (i,j,I_R_direct ,I_R_NIR)
             LOG_ERROR_CONT(*) 'DEBUG --- downward radiation (NIR,diffuse)   [J/m2/s]  :', RFLXD (i,j,I_R_diffuse,I_R_NIR)
             LOG_ERROR_CONT(*) 'DEBUG --- downward radiation (VIS,direct )   [J/m2/s]  :', RFLXD (i,j,I_R_direct ,I_R_VIS)
             LOG_ERROR_CONT(*) 'DEBUG --- downward radiation (VIS,diffuse)   [J/m2/s]  :', RFLXD (i,j,I_R_diffuse,I_R_VIS)
             LOG_ERROR_CONT(*) ''
             LOG_ERROR_CONT(*) 'DEBUG --- soil temperature                   [K]       :', TG    (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- surface temperature                [K]       :', TMPS  (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- efficiency of evaporation          [1]       :', QVEF  (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- surface albedo for (IR, direct )   [1]       :', ALBEDO(i,j,I_R_direct ,I_R_IR )
             LOG_ERROR_CONT(*) 'DEBUG --- surface albedo for (IR, diffuse)   [1]       :', ALBEDO(i,j,I_R_diffuse,I_R_IR )
             LOG_ERROR_CONT(*) 'DEBUG --- surface albedo for (NIR,direct )   [1]       :', ALBEDO(i,j,I_R_direct ,I_R_NIR)
             LOG_ERROR_CONT(*) 'DEBUG --- surface albedo for (NIR,diffuse)   [1]       :', ALBEDO(i,j,I_R_diffuse,I_R_NIR)
             LOG_ERROR_CONT(*) 'DEBUG --- surface albedo for (VIS,direct )   [1]       :', ALBEDO(i,j,I_R_direct ,I_R_VIS)
             LOG_ERROR_CONT(*) 'DEBUG --- surface albedo for (VIS,diffuse)   [1]       :', ALBEDO(i,j,I_R_diffuse,I_R_VIS)
             LOG_ERROR_CONT(*) 'DEBUG --- soil depth                         [m]       :', DZG   (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- thermal conductivity for soil      [J/m/K/s] :', TCS   (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- roughness length for momemtum      [m]       :', Z0M   (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- roughness length for heat          [m]       :', Z0H   (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- roughness length for vapor         [m]       :', Z0E   (i,j)
             LOG_ERROR_CONT(*) ''
             LOG_ERROR_CONT(*) 'DEBUG --- latent heat                        [J/kg]    :', LHV   (i,j)
             LOG_ERROR_CONT(*) 'DEBUG --- friction velocity                  [m]       :', Ustar
             LOG_ERROR_CONT(*) 'DEBUG --- friction potential temperature     [K]       :', Tstar
             LOG_ERROR_CONT(*) 'DEBUG --- friction water vapor mass ratio    [kg/kg]   :', Qstar
             LOG_ERROR_CONT(*) 'DEBUG --- d(friction velocity)               [m]       :', dUstar
             LOG_ERROR_CONT(*) 'DEBUG --- d(friction potential temperature)  [K]       :', dTstar
             LOG_ERROR_CONT(*) 'DEBUG --- d(friction water vapor mass ratio) [kg/kg]   :', dQstar
             LOG_ERROR_CONT(*) 'DEBUG --- modified absolute velocity         [m/s]     :', Uabs
             LOG_ERROR_CONT(*) 'DEBUG --- next surface temperature           [K]       :', TMPS1 (i,j)

             call PRC_abort
          endif

        end if


        ! calculate surface flux

        qdry = 1.0_RP - QVA(i,j)
        Rtot = qdry * Rdry + QVA(i,j) * Rvap

        call qsat( TMPS1(i,j), PRSS(i,j), qdry, & ! [IN]
                   QVsat                        ) ! [OUT]

        QVS  = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * QVsat

        call BULKFLUX( &
            Ustar,     & ! [OUT]
            Tstar,     & ! [OUT]
            Qstar,     & ! [OUT]
            Uabs,      & ! [OUT]
            Ra,        & ! [OUT]
            FracU10,   & ! [OUT]
            FracT2,    & ! [OUT]
            FracQ2,    & ! [OUT]
            TMPA(i,j), & ! [IN]
            TMPS1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVS,       & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            Z1  (i,j), & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ZMFLX(i,j) = -RHOS(i,j) * Ustar * Ustar / Uabs * WA(i,j)
        XMFLX(i,j) = -RHOS(i,j) * Ustar * Ustar / Uabs * UA(i,j)
        YMFLX(i,j) = -RHOS(i,j) * Ustar * Ustar / Uabs * VA(i,j)
        SHFLX(i,j) = -RHOS(i,j) * Ustar * Tstar * CPdry
        LHFLX(i,j) = -RHOS(i,j) * Ustar * Qstar * LHV(i,j) * Ra / ( Ra + Rb(i,j) )

        GHFLX(i,j) = -2.0_RP * TCS(i,j) * ( TMPS1(i,j) - TG(i,j) ) / DZG(i,j)

        emis = ( 1.0_RP - ALBEDO(i,j,I_R_diffuse,I_R_IR) ) * STB * TMPS(i,j)**4

        LWD  = RFLXD(i,j,I_R_diffuse,I_R_IR)
        LWU  = RFLXD(i,j,I_R_diffuse,I_R_IR) * ALBEDO(i,j,I_R_diffuse,I_R_IR) + emis
        SWD  = RFLXD(i,j,I_R_direct ,I_R_NIR) &
             + RFLXD(i,j,I_R_diffuse,I_R_NIR) &
             + RFLXD(i,j,I_R_direct ,I_R_VIS) &
             + RFLXD(i,j,I_R_diffuse,I_R_VIS)
        SWU  = RFLXD(i,j,I_R_direct ,I_R_NIR) * ALBEDO(i,j,I_R_direct ,I_R_NIR) &
             + RFLXD(i,j,I_R_diffuse,I_R_NIR) * ALBEDO(i,j,I_R_diffuse,I_R_NIR) &
             + RFLXD(i,j,I_R_direct ,I_R_VIS) * ALBEDO(i,j,I_R_direct ,I_R_VIS) &
             + RFLXD(i,j,I_R_diffuse,I_R_VIS) * ALBEDO(i,j,I_R_diffuse,I_R_VIS)

        ! calculation for residual
        res = SWD - SWU + LWD - LWU - SHFLX(i,j) - LHFLX(i,j) + GHFLX(i,j)

        ! put residual in ground heat flux
        GHFLX(i,j) = GHFLX(i,j) - res

        ! diagnostic variables considering unstable/stable state
        !U10(i,j) = FracU10 * UA(i,j)
        !V10(i,j) = FracU10 * VA(i,j)
        !T2 (i,j) = ( 1.0_RP - FracT2 ) * TMPS1(i,j) + FracT2 * TMPA(i,j)
        !Q2 (i,j) = ( 1.0_RP - FracQ2 ) * QVS        + FracQ2 * QVA (i,j)

        ! diagnostic variables for neutral state
        U10(i,j) = UA   (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        V10(i,j) = VA   (i,j) * log( 10.0_RP / Z0M(i,j) ) / log( Z1(i,j) / Z0M(i,j) )
        T2 (i,j) = TMPS1(i,j) + ( TMPA(i,j) - TMPS1(i,j) ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0H(i,j) ) ) &
                                                           / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0H(i,j) ) )
        Q2 (i,j) = QVS       + (  QVA(i,j) - QVS       ) * ( log(  2.0_RP / Z0M(i,j) ) * log(  2.0_RP / Z0E(i,j) ) ) &
                                                         / ( log( Z1(i,j) / Z0M(i,j) ) * log( Z1(i,j) / Z0E(i,j) ) )

        TMPS(i,j) = TMPS1(i,j)

      else

        ! not calculate surface flux
        ZMFLX(i,j) = 0.0_RP
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
        SHFLX(i,j) = 0.0_RP
        LHFLX(i,j) = 0.0_RP
        GHFLX(i,j) = 0.0_RP
        U10  (i,j) = 0.0_RP
        V10  (i,j) = 0.0_RP
        T2   (i,j) = 0.0_RP
        Q2   (i,j) = 0.0_RP

      end if

    end do
    end do

    return
  end subroutine CPL_PHY_SFC_skin

end module scale_cpl_phy_sfc_skin
