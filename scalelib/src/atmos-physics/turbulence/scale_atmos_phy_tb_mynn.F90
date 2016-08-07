!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Boundary layer turbulence model
!!          Mellor-Yamada Nakanishi-Niino model
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-08-27 (S.Nishizawa) [new]
!!
!! - Reference
!!  - Nakanishi and Niino, 2009:
!!    Development of an improved turbulence closure model for the atmospheric boundary layer.
!!    J. Meteorol. Soc. Japan, 87, 895-912
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_tb_mynn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

#include "macro_thermodyn.h"

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
  public :: ATMOS_PHY_TB_mynn_setup
  public :: ATMOS_PHY_TB_mynn

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
  real(RP), parameter :: OneOverThree = 1.0_RP / 3.0_RP
  real(RP), parameter :: TKE_min = 1.0E-10_RP
  real(RP), parameter :: LT_min = 1.0E-6_RP
  real(RP), parameter :: Us_min = 1.0E-6_RP

  real(RP)              :: A1
  real(RP)              :: A2
  real(RP), parameter   :: B1 = 24.0_RP
  real(RP), parameter   :: B2 = 15.0_RP
  real(RP)              :: C1
  real(RP), parameter   :: C2 = 0.75_RP
  real(RP), parameter   :: C3 = 0.352_RP
  real(RP), parameter   :: C5 = 0.2_RP
  real(RP), parameter   :: G1 = 0.235_RP
  real(RP)              :: G2
  real(RP)              :: F1
  real(RP)              :: F2
  real(RP)              :: Rf1
  real(RP)              :: Rf2
  real(RP)              :: Rfc
  real(RP)              :: AF12 !< A1 F1 / A2 F2
  real(RP), parameter   :: PrN = 0.74_RP

  real(RP)              :: SQRT_2PI
  real(RP)              :: RSQRT_2PI
  real(RP)              :: RSQRT_2

  integer :: KE_PBL
  logical :: ATMOS_PHY_TB_MYNN_TKE_INIT = .false. !< set tke with that of level 2 at the first time if .true.
  real(RP) :: ATMOS_PHY_TB_MYNN_N2_MAX = 1.E3_RP
  real(RP) :: ATMOS_PHY_TB_MYNN_NU_MIN = 1.E-6_RP
  real(RP) :: ATMOS_PHY_TB_MYNN_NU_MAX = 10000._RP
  real(RP) :: ATMOS_PHY_TB_MYNN_KH_MIN = 1.E-6_RP
  real(RP) :: ATMOS_PHY_TB_MYNN_KH_MAX = 10000._RP
  real(RP) :: ATMOS_PHY_TB_MYNN_Lt_MAX = 700._RP ! ~ 0.23 * 3 km

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_mynn_setup( &
       TYPE_TB,       &
       CDZ, CDX, CDY, &
       CZ             )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI => CONST_PI
    implicit none

    character(len=*), intent(in) :: TYPE_TB

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA,IA,JA)

    real(RP) :: ATMOS_PHY_TB_MYNN_PBL_MAX = 1e10_RP !< maximum height of the PBL

    NAMELIST / PARAM_ATMOS_PHY_TB_MYNN / &
         ATMOS_PHY_TB_MYNN_TKE_INIT, &
         ATMOS_PHY_TB_MYNN_PBL_MAX, &
         ATMOS_PHY_TB_MYNN_N2_MAX, &
         ATMOS_PHY_TB_MYNN_NU_MIN, &
         ATMOS_PHY_TB_MYNN_NU_MAX, &
         ATMOS_PHY_TB_MYNN_KH_MIN, &
         ATMOS_PHY_TB_MYNN_KH_MAX, &
         ATMOS_PHY_TB_MYNN_Lt_MAX

    integer :: k, i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[TURBULENCE] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Mellor-Yamada Nakanishi-Niino Model'

    if ( TYPE_TB /= 'MYNN' ) then
       write(*,*) 'xxx ATMOS_PHY_TB_TYPE is not MYNN. Check!'
       call PRC_MPIstop
    endif


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_MYNN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_TB_MYNN. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_TB_MYNN)

    do k = KS, KE-1
       do j = JS, JE
       do i = IS, IE
          if ( ATMOS_PHY_TB_MYNN_PBL_MAX >= CZ(k,i,j) ) then
             KE_PBL = k
          end if
       end do
       end do
    end do

    A1 = B1 * (1.0_RP - 3.0_RP * G1) / 6.0_RP
    A2 = 1.0_RP / (3.0_RP * G1 * B1**(1.0_RP/3.0_RP) * PrN )
    C1 = G1 - 1.0_RP / ( 3.0_RP * A1 * B1**(1.0_RP/3.0_RP) )
    G2 = ( 2.0_RP * A1 * (3.0_RP - 2.0_RP * C2) + B2 * (1.0_RP - C3) ) / B1
    F1 = B1 * (G1 - C1) + 2.0_RP * A1 * (3.0_RP - 2.0_RP * C2) + 3.0_RP * A2 * (1.0_RP - C2) * (1.0_RP - C5)
    F2 = B1 * (G1 + G2) - 3.0_RP * A1 * (1.0_RP - C2)

    Rf1 = B1 * (G1 - C1) / F1
    Rf2 = B1 * G1 / F2
    Rfc = G1 / (G1 + G2)

    AF12 = A1 * F1 / ( A2 * F2 )

    SQRT_2PI = sqrt( 2.0_RP * PI )
    RSQRT_2PI = 1.0_RP / SQRT_2PI
    RSQRT_2 = 1.0_RP / sqrt( 2.0_RP )

    return
  end subroutine ATMOS_PHY_TB_mynn_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_mynn( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       tke,                                         & ! (inout)
       tke_t, Nu, Ri, Pr, N2,                       & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          & ! (in)
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)
    use scale_precision
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       R      => CONST_Rdry, &
       Rvap   => CONST_Rvap, &
       CP     => CONST_CPdry, &
       EPSTvap => CONST_EPSTvap
    use scale_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use scale_grid_real, only: &
       CZ => REAL_CZ
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_XVW, &
       I_UYW, &
       I_UYZ, &
       I_XVZ, &
       I_XY, &
       I_UY, &
       I_XV
    use scale_comm, only: &
       COMM_vars8, &
       COMM_vars, &
       COMM_wait
    use scale_atmos_phy_tb_common, only: &
       diffusion_solver => ATMOS_PHY_TB_diffusion_solver
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_temp_pres
    use scale_atmos_hydrometer, only: &
       I_QV, &
       I_QC, &
       I_QI, &
       ATMOS_HYDROMETER_templhv, &
       ATMOS_HYDROMETER_templhs
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_alpha, &
       ATMOS_SATURATION_pres2qsat => ATMOS_SATURATION_pres2qsat_all
#ifdef MORE_HIST
    use scale_history, only: &
       HIST_in
#endif
    implicit none

    real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

    real(RP), intent(inout) :: tke (KA,IA,JA) ! TKE
    real(RP), intent(out) :: tke_t(KA,IA,JA) ! tendency of TKE
    real(RP), intent(out) :: Nu(KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Ri(KA,IA,JA) ! Richardson number
    real(RP), intent(out) :: Pr(KA,IA,JA) ! Plandtle number
    real(RP), intent(out) :: N2(KA,IA,JA) ! squared Brunt-Vaisala frequency

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)  :: SFLX_MW(IA,JA)
    real(RP), intent(in)  :: SFLX_MU(IA,JA)
    real(RP), intent(in)  :: SFLX_MV(IA,JA)
    real(RP), intent(in)  :: SFLX_SH(IA,JA)
    real(RP), intent(in)  :: SFLX_QV(IA,JA)

    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
    real(DP), intent(in)  :: dt

    real(RP) :: U(KA,IA,JA)    !< velocity in x-direction (full level)
    real(RP) :: V(KA,IA,JA)    !< velocity in y-direction (full level)
    real(RP) :: phiN(KA,IA,JA)


    real(RP) :: sm(KA,IA,JA) !< stability function for velocity
    real(RP) :: sh(KA,IA,JA) !< stability function for scalars
    real(RP) :: l(KA,IA,JA) !< length scale L
    real(RP) :: q(KA,IA,JA) !< q
    real(RP) :: dudz2(KA,IA,JA) !< (dudz)^2 + (dvdz)^2
    real(RP) :: q2_2(KA,IA,JA)  !< q^2 for level 2
    real(RP) :: Kh              !< eddy diffusion coefficient
    real(RP) :: RHOKh (KA,IA,JA)!< mass-weighted eddy diffusion coefficient

    real(RP) :: SFLX_PT(IA,JA) ! surface potential temperature flux

    real(RP) :: a(KA,IA,JA)
    real(RP) :: b(KA,IA,JA)
    real(RP) :: c(KA,IA,JA)
    real(RP) :: d(KA,IA,JA)
    real(RP) :: ap
    real(RP) :: tke_N(KA,IA,JA)

    real(RP) :: POTT(KA,IA,JA) !< potential temperature
    real(RP) :: POTV(KA,IA,JA) !< virtual potential temperature
    real(RP) :: POTL(KA,IA,JA) !< liquid water potential temperature
    real(RP) :: TEML(KA,IA,JA) !< liquid water temperature

    real(RP) :: Qw(KA,IA,JA)   !< total water
    real(RP) :: ql             !< liquid water
    real(RP) :: qs             !< solid water
    real(RP) :: qdry           !< dry air

    real(RP) :: temp(KA,IA,JA) !< temperature
    real(RP) :: pres(KA,IA,JA) !< pressure

    real(RP) :: lh(KA,IA,JA)
    real(RP) :: lhv(KA,IA,JA)
    real(RP) :: lhs(KA,IA,JA)
    real(RP) :: alpha(KA,IA,JA)

    real(RP) :: ac           !< \alpha_c

    real(RP) :: Q1
    real(RP) :: Qsl(KA,IA,JA)
    real(RP) :: dQsl
    real(RP) :: sigma_s
    real(RP) :: RR
    real(RP) :: Rt
    real(RP) :: betat
    real(RP) :: betaq
    real(RP) :: aa
    real(RP) :: bb
    real(RP) :: cc

    real(RP) :: flux_z(KA,IA,JA)
    real(RP) :: flux_x(KA,IA,JA)
    real(RP) :: flux_y(KA,IA,JA)
    real(RP) :: mflx
    real(RP) :: advc

#ifdef MORE_HIST
    real(RP) :: adv(KA,IA,JA)
    real(RP) :: gen(KA,IA,JA)
#endif

    integer :: k, i, j, iq
    integer :: IIS, IIE, JJS, JJE

    if ( IO_L ) write(IO_FID_LOG, *) "*** Physics step: Turbulence (MYNN)"

#ifdef DEBUG
    POTT(:,:,:) = UNDEF
    POTL(:,:,:) = UNDEF
    temp(:,:,:) = UNDEF
#endif

!OCL XFILL
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE
       qflx_sgs_momz(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS  , JE
    do i = IS  , IE+1
    do k = KS-1, KE+1
       qflx_sgs_momx(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do

!OCL XFILL
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE+1
       qflx_sgs_momy(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE+1
       qflx_sgs_rhot(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do iq = 1, QA
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE+1
       qflx_sgs_rhoq(k,i,j,XDIR,iq) = 0.0_RP
    end do
    end do
    end do
    end do

!OCL XFILL
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE
       qflx_sgs_momz(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_momx(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS  , JE+1
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_momy(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_rhot(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do iq = 1, QA
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_rhoq(k,i,j,YDIR,iq) = 0.0_RP
    end do
    end do
    end do
    end do
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE+1
       qflx_sgs_momz(k,i,j,ZDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do iq = I_QV+1, QA
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_rhoq(k,i,j,ZDIR,iq) = 0.0_RP
    end do
    end do
    end do
    end do



    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

!OCL XFILL
       do j = JJS  , JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE_PBL+1
          U(k,i,j) = 0.5_RP * ( MOMX(k,i,j) + MOMX(k,i-1,j) ) / DENS(k,i,j)
       end do
       end do
       end do

!OCL XFILL
       do j = JJS-1, JJE+1
       do i = IIS  , IIE+1
       do k = KS, KE_PBL+1
          V(k,i,j) = 0.5_RP * ( MOMY(k,i,j) + MOMY(k,i,j-1) ) / DENS(k,i,j)
       end do
       end do
       end do

    end do
    end do

    call ATMOS_THERMODYN_temp_pres( temp, pres, & ! (out)
                                    DENS, RHOT, QTRC, & ! (in)
                                    TRACER_CV, TRACER_R, TRACER_MASS ) ! (in)


    call ATMOS_HYDROMETER_templhv( lhv, temp )
    call ATMOS_HYDROMETER_templhs( lhs, temp )
    call ATMOS_SATURATION_alpha( alpha, temp )

!OCL LOOP_NOFUSION,PREFETCH_SEQUENTIAL(SOFT),SWP
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE_PBL+1
          ql = 0.0_RP
          if ( I_QC > 0 ) ql = QTRC(k,i,j,I_QC)
!          do iq = QWS, QWE
!             ql = ql + QTRC(k,i,j,iq)
!          end do
          qs = 0.0_RP
          if ( I_QI > 0 ) qs = QTRC(k,i,j,I_QI)
!          do iq = QIS, QIE
!             qs = qs + QTRC(k,i,j,iq)
!          end do
          CALC_QDRY(qdry, QTRC, TRACER_MASS, k, i, j, iq)

          Qw(k,i,j) = QTRC(k,i,j,I_QV) + ql + qs

          lh(k,i,j) = lhv(k,i,j) * alpha(k,i,j) + lhs(k,i,j) * ( 1.0_RP-alpha(k,i,j) )

          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
          ! liquid water potential temperature
          POTL(k,i,j) = POTT(k,i,j) * (1.0_RP - 1.0_RP * (lhv(k,i,j) * ql + lhs(k,i,j) * qs) / ( temp(k,i,j) * CP ) )
          TEML(k,i,j) = POTL(k,i,j) * temp(k,i,j) / POTT(k,i,j)

          ! virtual potential temperature for derivertive
!          POTV(k,i,j) = ( 1.0_RP + EPSTvap * Qw(k,i,j) ) * POTL(k,i,j)
          POTV(k,i,j) = ( qdry + (EPSTvap+1.0_RP) * QTRC(k,i,j,I_QV) ) * POTL(k,i,j)

       end do

       SFLX_PT(i,j) = SFLX_SH(i,j) / ( CP * DENS(KS,i,j) ) &
                    * POTT(KS,i,j) / TEMP(KS,i,j)

       n2(KS,i,j) = min(ATMOS_PHY_TB_MYNN_N2_MAX, &
                        GRAV * ( POTV(KS+1,i,j) - POTV(KS,i,j) ) &
                             / ( ( CZ(KS+1,i,j)-CZ(KS,i,j) ) * POTV(KS,i,j) ) )
       do k = KS+1, KE_PBL
          n2(k,i,j) = min(ATMOS_PHY_TB_MYNN_N2_MAX, &
                          GRAV * ( POTV(k+1,i,j) - POTV(k-1,i,j) ) &
                               / ( ( CZ(k+1,i,j)-CZ(k-1,i,j) ) * POTV(k,i,j) ) )
       end do

       dudz2(KS,i,j) = ( ( U(KS+1,i,j) - U(KS,i,j) )**2 + ( V(KS+1,i,j) - V(KS,i,j) )**2 ) &
                     / ( CZ(KS+1,i,j) - CZ(KS,i,j) )**2
       do k = KS+1, KE_PBL
          dudz2(k,i,j) = ( ( U(k+1,i,j) - U(k-1,i,j) )**2 + ( V(k+1,i,j) - V(k-1,i,j) )**2 ) &
                       / ( CZ(k+1,i,j) - CZ(k-1,i,j) )**2
       end do

       do k = KS, KE_PBL
          Ri(k,i,j) = n2(k,i,j) / max(dudz2(k,i,j), 1E-10_RP)
       end do

    end do
    end do

    if ( ATMOS_PHY_TB_MYNN_TKE_INIT ) then
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE_PBL
          q(k,i,j) = sqrt( TKE_MIN*2.0_RP )
       end do
       end do
       end do
    else
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE_PBL
          q(k,i,j) = sqrt( max(tke(k,i,j), TKE_MIN)*2.0_RP )
       end do
       end do
       end do
    end if

    ! length
    call get_length( &
         l, & ! (out)
         DENS, & ! (in)
         q, n2, & ! (in)
         SFLX_MU, SFLX_MV, SFLX_PT, & ! (in)
         POTT ) ! (in)

    call get_q2_level2( &
         q2_2, & ! (out)
         dudz2, Ri, l ) ! (in)

    if ( ATMOS_PHY_TB_MYNN_TKE_INIT ) then
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE_PBL
          tke(k,i,j) = max(q2_2(k,i,j) * 0.5_RP, TKE_MIN)
          q(k,i,j) = sqrt( tke(k,i,j) * 2.0_RP )
       end do
       end do
       end do
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KE_PBL+1, KE
          tke(k,i,j) = 0.0_RP
       end do
       end do
       end do
       call COMM_vars8(tke, 1)
       call COMM_wait (tke, 1)

       ATMOS_PHY_TB_MYNN_TKE_INIT = .false.
    end if

    call get_smsh( sm, sh, & ! (out)
                   q, q2_2, & ! (in)
                   l, n2, dudz2 ) ! (in)

    call ATMOS_SATURATION_pres2qsat( Qsl(KS:KE_PBL,:,:), & ! (out)
                                     TEML(KS:KE_PBL,:,:), pres(KS:KE_PBL,:,:), & ! (in)
                                     KE_PBL-KS+1 ) ! (in)

!OCL LOOP_NOFUSION,PREFETCH_SEQUENTIAL(SOFT),SWP
    do j = JS, JE
    do i = IS, IE
       do k = KS+1, KE_PBL

          dQsl = Qsl(k,i,j) * lh(k,i,j) / ( Rvap * POTL(k,i,j)**2 )
          aa = 1.0_RP / ( 1.0_RP + lh(k,i,j)/CP * dQsl )
          bb = TEMP(k,i,j) / POTT(k,i,j) * dQsl
          ac = min( q(k,i,j)/sqrt(q2_2(k,i,j)), 1.0_RP )
          sigma_s = max( sqrt( 0.25_RP * aa**2 * l(k,i,j)**2 * ac * B2 * sh(k,i,j) ) &
                       * abs( Qw(k+1,i,j) - Qw(k-1,i,j) - bb * ( POTL(k+1,i,j)-POTL(k-1,i,j) ) ) &
                       / ( CZ(k+1,i,j) - CZ(k-1,i,j) ), &
                       1e-10_RP )
          Q1 = aa * ( Qw(k,i,j) - Qsl(k,i,j) ) * 0.5_RP / sigma_s
          RR = min( max( 0.5_RP * ( 1.0_RP + erf(Q1*rsqrt_2) ), 0.0_RP ), 1.0_RP )
          Ql = min( max( 2.0_RP * sigma_s * ( RR * Q1 + rsqrt_2pi * exp(-0.5_RP*Q1**2) ), &
                    0.0_RP ), &
                    Qw(k,i,j) * 0.5_RP )
          cc = ( 1.0_RP + EPSTvap * Qw(k,i,j) - (1.0_RP+EPSTvap) * Ql ) * POTT(k,i,j)/TEMP(k,i,j) * lh(k,i,j) / CP &
               - (1.0_RP+EPSTvap) * POTT(k,i,j)
          Rt = min( max( RR - Ql / (2.0_RP*sigma_s*sqrt_2pi) * exp(-Q1**2 * 0.5_RP), 0.0_RP ), 1.0_RP )
          betat = 1.0_RP + EPSTvap * Qw(k,i,j) - (1.0_RP+EPSTvap) * Ql - Rt * aa * bb * cc
          betaq = EPSTvap * POTT(k,i,j) + Rt * aa * cc
          n2(k,i,j) = min(ATMOS_PHY_TB_MYNN_N2_MAX, &
                          GRAV * ( ( POTL(k+1,i,j) - POTL(k-1,i,j) ) * betat &
                                 + ( Qw  (k+1,i,j) - Qw  (k-1,i,j) ) * betaq ) &
                               / ( CZ(k+1,i,j) - CZ(k-1,i,j) ) / POTV(k,i,j) )
       end do
       n2(KS,i,j) = n2(KS+1,i,j)

       do k = KS, KE_PBL
          Ri(k,i,j) = n2(k,i,j) / max(dudz2(k,i,j), 1E-10_RP)
       end do
       do k = KE_PBL+1, KE
          Ri(k,i,j) = 0.0_RP
          n2(k,i,j) = 0.0_RP
       end do

    end do
    end do

    ! length
    call get_length( &
         l, & ! (out)
         DENS, & ! (in)
         q, n2, & ! (in)
         SFLX_MU, SFLX_MV, SFLX_PT, & ! (in)
         POTT ) ! (in)

    call get_q2_level2( &
         q2_2, & ! (out)
         dudz2, Ri, l ) ! (in)

    call get_smsh( &
         sm, sh, & ! (out)
         q, q2_2, & ! (in)
         l, n2, dudz2 ) ! (in)


    do j = JS, JE
    do i = IS, IE
       do k = 1, KS-1
          Nu(k,i,j) = 0.0_RP
       end do
       do k = KS, KE_PBL
          Nu(k,i,j) = max( min( l(k,i,j) * q(k,i,j) * sm(k,i,j), &
                           ATMOS_PHY_TB_MYNN_NU_MAX ), &
                           ATMOS_PHY_TB_MYNN_NU_MIN )
       end do
       do k = KE_PBL+1, KA
          Nu(k,i,j) = 0.0_RP
       end do

       do k = KS, KE_PBL
          Kh = max( min( l(k,i,j) * q(k,i,j) * sh(k,i,j), &
                    ATMOS_PHY_TB_MYNN_KH_MAX ), &
                    ATMOS_PHY_TB_MYNN_KH_MIN )
          RHOKh(k,i,j) = DENS(k,i,j) * Kh
          Pr(k,i,j) = Nu(k,i,j) / Kh
       end do
       do k = KE_PBL+1, KE
          Pr(k,i,j) = 1.0_RP
       end do
    end do
    end do

    call COMM_vars( Nu,   1 )

    ! time integration

    !  for velocities
!OCL INDEPENDENT
    do j = JS, JE
    do i = IS, IE

       ap = - dt * 0.5_RP * ( DENS(KS  ,i,j)*Nu(KS  ,i,j) &
                            + DENS(KS+1,i,j)*Nu(KS+1,i,j) ) &
                          * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
       a(KS,i,j) = ap * RCDZ(KS) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
       c(KS,i,j) = 0.0_RP
       b(KS,i,j) = - a(KS,i,j) + 1.0_RP
       do k = KS+1, KE_PBL-1
          c(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
          ap = - dt * 0.5_RP * ( DENS(k  ,i,j)*Nu(k  ,i,j) &
                               + DENS(k+1,i,j)*Nu(k+1,i,j) ) &
                             * RFDZ(k) / GSQRT(k,i,j,I_XYW)
          a(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
          b(k,i,j) = - a(k,i,j) - c(k,i,j) + 1.0_RP
       end do
       a(KE_PBL,i,j) = 0.0_RP
       c(KE_PBL,i,j) = ap * RCDZ(KE_PBL) / ( DENS(KE_PBL,i,j) * GSQRT(KE_PBL,i,j,I_XYZ) )
       b(KE_PBL,i,j) = - c(KE_PBL,i,j) + 1.0_RP

       d(KS,i,j) = U(KS,i,j) + dt * SFLX_MU(i,j) * RCDZ(KS) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
       do k = KS+1, KE_PBL
          d(k,i,j) = U(k,i,j)
       end do
       call diffusion_solver( &
            phiN(:,i,j),                     & ! (out)
            a(:,i,j), b(:,i,j), c(:,i,j), d(:,i,j), & ! (in)
            KE_PBL                           ) ! (in)
    end do
    end do

    call COMM_vars( phiN, 2 )
    call COMM_wait( Nu,   1 )
    call COMM_wait( phiN, 2 )

!OCL INDEPENDENT
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE_PBL-1
          qflx_sgs_momx(k,i,j,ZDIR) = - 0.03125_RP & ! 1/4/4/2
               * ( DENS(k,i,j) + DENS(k+1,i,j) + DENS(k,i+1,j) + DENS(k+1,i+1,j) ) &
               * ( Nu(k,i,j) + Nu(k+1,i,j) + Nu(k,i+1,j) + Nu(k+1,i+1,j) ) &
               * ( (phiN(k+1,i,j)+phiN(k+1,i+1,j)) - (phiN(k,i,j)+phiN(k,i+1,j)) ) &
               * J33G * RFDZ(k) / GSQRT(k,i,j,I_UYW)
       end do
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momx(KS-1,i,j,ZDIR) = 0.0_RP
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KE_PBL, KE
          qflx_sgs_momx(k,i,j,ZDIR) = 0.0_RP
       end do
       end do
       end do


       ! integration V
!OCL INDEPENDENT
       do j = JJS, JJE
       do i = IIS, IIE
          d(KS,i,j) = V(KS,i,j) + dt * SFLX_MV(i,j) * RCDZ(KS) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
          do k = KS+1, KE_PBL
             d(k,i,j) = V(k,i,j)
          end do
          call diffusion_solver( &
               phiN(:,i,j),                     & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d(:,i,j), & ! (in)
               KE_PBL                           ) ! (in)
       end do
       end do

    end do
    end do

    call COMM_vars( phiN, 1 )
    call COMM_wait( phiN, 1 )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momy(KS-1,i,j,ZDIR) = 0.0_RP
          do k = KS, KE_PBL-1
             qflx_sgs_momy(k,i,j,ZDIR) = - 0.03125_RP & ! 1/4/4/2
               * ( DENS(k,i,j) + DENS(k+1,i,j) + DENS(k,i,j+1) + DENS(k+1,i,j+1) ) &
               * ( Nu(k,i,j) + Nu(k+1,i,j) + Nu(k,i,j+1) + Nu(k+1,i,j+1) ) &
               * ( (phiN(k+1,i,j)+phiN(k+1,i,j+1)) - (phiN(k,i,j)+phiN(k,i,j+1)) ) &
               * J33G * RFDZ(k) / GSQRT(k,i,j,I_XVW)
          end do
          do k = KE_PBL, KE
             qflx_sgs_momy(k,i,j,ZDIR) = 0.0_RP
          end do
       end do
       end do


       ! for scalars
       ! integration POTT
!OCL INDEPENDENT
       do j = JJS, JJE
       do i = IIS, IIE
          ap = - dt * 0.5_RP * ( RHOKh(KS,i,j) + RHOKh(KS+1,i,j) ) &
             * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
          a(KS,i,j) = ap * RCDZ(KS) / (DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
          c(KS,i,j) = 0.0_RP
          b(KS,i,j) = - a(KS,i,j) + 1.0_RP
          do k = KS+1, KE_PBL-1
             c(k,i,j) = ap * RCDZ(k) / (DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             ap = - dt * 0.5_RP * ( RHOKh(k,i,j) + RHOKh(k+1,i,j) ) &
                * RFDZ(k) / GSQRT(k,i,j,I_XYW)
             a(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             b(k,i,j) = - a(k,i,j) - c(k,i,j) + 1.0_RP
          end do
          a(KE_PBL,i,j) = 0.0_RP
          c(KE_PBL,i,j) = ap * RCDZ(KE_PBL) / (DENS(KE_PBL,i,j) * GSQRT(KE_PBL,i,j,I_XYZ) )
          b(KE_PBL,i,j) = - c(KE_PBL,i,j) + 1.0_RP
          d(KS,i,j) = POTL(KS,i,j) + dt * SFLX_PT(i,j) * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ)
          do k = KS+1, KE_PBL
             d(k,i,j) = POTL(k,i,j)
          end do
          call diffusion_solver( &
               phiN(:,i,j),                     & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d(:,i,j), & ! (in)
               KE_PBL                           ) ! (in)
       end do
       end do

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_rhot(KS-1,i,j,ZDIR) = 0.0_RP
          do k = KS, KE_PBL-1
             qflx_sgs_rhot(k,i,j,ZDIR) = - 0.5_RP & ! 1/2
               * ( RHOKh(k,i,j) + RHOKh(k+1,i,j) ) &
               * J33G * ( phiN(k+1,i,j) - PhiN(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_XYW)
          end do
          do k = KE_PBL, KE
             qflx_sgs_rhot(k,i,j,ZDIR) = 0.0_RP
          end do
       end do
       end do


       ! integration QV
       iq = I_QV
!OCL INDEPENDENT
       do j = JJS, JJE
       do i = IIS, IIE
          d(KS,i,j) = Qw(KS,i,j) + dt * SFLX_QV(i,j) * RCDZ(KS) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
          do k = KS+1, KE_PBL
             d(k,i,j) = Qw(k,i,j)
          end do
          call diffusion_solver( &
               phiN(:,i,j),                     & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d(:,i,j), & ! (in)
               KE_PBL                           ) ! (in)
       end do
       end do

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_rhoq(KS-1,i,j,ZDIR,iq) = 0.0_RP
          do k = KS, KE_PBL-1
             qflx_sgs_rhoq(k,i,j,ZDIR,iq) = - 0.5_RP & ! 1/2
                  * ( RHOKh(k,i,j) + RHOKh(k+1,i,j) ) &
                  * J33G * ( phiN(k+1,i,j) - phiN(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_XYW)
          end do
          do k = KE_PBL, KE
             qflx_sgs_rhoq(k,i,j,ZDIR,iq) = 0.0_RP
          end do
       end do
       end do


       ! time integration tke

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS , KE_PBL-1
          mflx = J33G * MOMZ(k,i,j) / ( MAPF(i,j,1,I_XY)*MAPF(i,j,2,I_XY) ) &
               + J13G(k,i,j,I_XYW) * 0.25_RP * ( MOMX(k,i-1,j)+MOMX(k,i,j)+MOMX(k+1,i-1,j)+MOMX(k+1,i,j) ) / MAPF(i,j,2,I_XY) &
               + J23G(k,i,j,I_XYW) * 0.25_RP * ( MOMY(k,i,j-1)+MOMY(k,i,j)+MOMY(k+1,i,j-1)+MOMY(k+1,i,j) ) / MAPF(i,j,1,I_XY)
          flux_z(k,i,j) = 0.5_RP * (    mflx  * ( tke(k+1,i,j)+tke(k,i,j) ) &
                                   -abs(mflx) * ( tke(k+1,i,j)-tke(k,i,j) ) )
       end do
       end do
       end do
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS   , KE_PBL
          mflx = MOMX(k,i,j) * GSQRT(k,i,j,I_UYZ) / MAPF(i,j,2,I_UY)
          flux_x(k,i,j) = 0.5_RP * (    mflx  * ( tke(k,i+1,j)+tke(k,i,j) ) &
                                   -abs(mflx) * ( tke(k,i+1,j)-tke(k,i,j) ) )
       end do
       end do
       end do
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS   , KE_PBL
          mflx = MOMY(k,i,j) * GSQRT(k,i,j,I_XVZ) / MAPF(i,j,1,I_XV)
          flux_y(k,i,j) = 0.5_RP * (    mflx  * ( tke(k,i,j+1)+tke(k,i,j) ) &
                                   -abs(mflx) * ( tke(k,i,j+1)-tke(k,i,j) ) )
       end do
       end do
       end do

       do j = JJS, JJE
       do i = IIS, IIE
          advc = ( ( flux_z(KS,i,j)                    ) * RCDZ(KS) &
                 + ( flux_x(KS,i,j) - flux_x(KS,i-1,j) ) * RCDX(i) &
                 + ( flux_y(KS,i,j) - flux_y(KS,i,j-1) ) * RCDY(j) ) &
                 * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / ( GSQRT(KS,i,j,I_XYZ) * DENS(KS,i,j) )
          d(KS,i,j) = tke(KS,i,j) + dt * ( Nu(KS,i,j) * (dudz2(KS,i,j) - n2(KS,i,j)/Pr(KS,i,j)) &
                                     - advc )
#ifdef MORE_HIST
          adv(KS,i,j) = advc
          gen(KS,i,j) = Nu(KS,i,j) * (dudz2(KS,i,j) - n2(KS,i,j)/Pr(KS,i,j))
#endif
          do k = KS+1, KE_PBL-1
             advc = ( ( flux_z(k,i,j) - flux_z(k-1,i,j) ) * RCDZ(k) &
                    + ( flux_x(k,i,j) - flux_x(k,i-1,j) ) * RCDX(i) &
                    + ( flux_y(k,i,j) - flux_y(k,i,j-1) ) * RCDY(j) ) &
                    * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / ( GSQRT(k,i,j,I_XYZ) * DENS(k,i,j) )
             d(k,i,j) = tke(k,i,j) &
                  + dt * ( Nu(k,i,j) * (dudz2(k,i,j) - n2(k,i,j)/Pr(k,i,j)) &
                         - advc )
#ifdef MORE_HIST
             adv(k,i,j) = advc
             gen(k,i,j) = Nu(k,i,j) * (dudz2(k,i,j) - n2(k,i,j)/Pr(k,i,j))
#endif
          end do
          advc = ( (                    - flux_z(KE_PBL-1,i,j) ) * RCDZ(KE_PBL) &
                 + ( flux_x(KE_PBL,i,j) - flux_x(KE_PBL,i-1,j) ) * RCDX(i) &
                 + ( flux_y(KE_PBL,i,j) - flux_y(KE_PBL,i,j-1) ) * RCDY(j) ) &
                 * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) / ( GSQRT(KE_PBL,i,j,I_XYZ) * DENS(KE_PBL,i,j) )
          d(KE_PBL,i,j) = tke(KE_PBL,i,j) + dt * ( Nu(KE_PBL,i,j) * (dudz2(KE_PBL,i,j) - n2(KE_PBL,i,j)/Pr(KE_PBL,i,j)) &
                                             - advc )
#ifdef MORE_HIST
          adv(KE_PBL,i,j) = advc
          gen(KE_PBL,i,j) = Nu(KE_PBL,i,j) * (dudz2(KE_PBL,i,j) - n2(KE_PBL,i,j)/Pr(KE_PBL,i,j))
#endif
       end do
       end do

!OCL INDEPENDENT
       do j = JJS, JJE
       do i = IIS, IIE
          ap = - dt * 1.5_RP * ( DENS(KS  ,i,j)*Nu(KS  ,i,j) &
                               + DENS(KS+1,i,j)*Nu(KS+1,i,j) ) &
             * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
          a(KS,i,j) = ap * RCDZ(KS) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
          c(KS,i,j) = 0.0_RP
          b(KS,i,j) = - a(KS,i,j) + 1.0_RP + 2.0_RP * dt * q(KS,i,j) / ( B1 * l(KS,i,j) )
          do k = KS+1, KE_PBL-1
             c(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             ap = - dt * 1.5_RP * ( DENS(k  ,i,j)*Nu(k  ,i,j) &
                                  + DENS(k+1,i,j)*Nu(k+1,i,j))  &
                * RFDZ(k) / GSQRT(k,i,j,I_XYW)
             a(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             b(k,i,j) = - a(k,i,j) - c(k,i,j) + 1.0_RP + 2.0_RP * dt * q(k,i,j) / ( B1 * l(k,i,j) )
          end do

          a(KE_PBL,i,j) = 0.0_RP
          c(KE_PBL,i,j) = ap * RCDZ(KE_PBL) / ( DENS(KE_PBL,i,j) * GSQRT(KE_PBL,i,j,I_XYZ) )
          b(KE_PBL,i,j) = - c(KE_PBL,i,j) + 1.0_RP + 2.0_RP * dt * q(KE_PBL,i,j) / ( B1 * l(KE_PBL,i,j) )
          call diffusion_solver( &
               tke_N(:,i,j),     & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d(:,i,j), & ! (in)
               KE_PBL                           ) ! (in)
#ifdef DEBUG
          do k = KS+1, KE_PBL-1
             tke(1,i,j) = d(k,i,j) &
               + ( ( (tke_N(k+1,i,j)-tke_N(k,i,j)) * RFDZ(k) / GSQRT(k,i,j,I_XYW) &
                   * (NU(k+1,i,j)*DENS(k+1,i,j)+NU(k,i,j)*DENS(k,i,j))*1.5_RP &
                   - (tke_N(k,i,j)-tke_N(k-1,i,j)) * RFDZ(k-1) / GSQRT(k-1,i,j,I_XYW) &
                   * (NU(k,i,j)*DENS(k,i,j)+NU(k-1,i,j)*DENS(k-1,i,j))*1.5_RP ) &
                 * RCDZ(k) / GSQRT(k,i,j,I_XYZ) / DENS(k,i,j) &
                 - 2.0_RP*tke_N(k,i,j) * q(k,i,j) / (B1 * l(k,i,j)) ) &
                 * dt
             if ( tke_N(k,i,j) > 0.1_RP .AND. abs(tke(1,i,j) - tke_N(k,i,j))/tke_N(k,i,j) > 1e-10_RP ) then
                advc = ( tke(k,i,j) - d(k,i,j) ) / dt + Nu(k,i,j) * (dudz2(k,i,j) - n2(k,i,j)/Pr(k,i,j))
                write(*,*)k,i,j,tke(1,i,j),tke_N(k,i,j), tke(k,i,j),nu(k,i,j),dudz2(k,i,j),n2(k,i,j),pr(k,i,j),advc
                open(90, file="mynn.dat")
                write(90,*)KE_PBL-KS+1
                write(90,*)dt, B1
                write(90,*)tke(KS:KE_PBL,i,j)
                write(90,*)q(KS:KE_PBL,i,j)
                write(90,*)l(KS:KE_PBL,i,j)
                write(90,*)rcdz(KS:KE_PBL)
                write(90,*)rfdz(KS:KE_PBL)
                write(90,*)nu(KS:KE_PBL,i,j)
                write(90,*)dens(KS:KE_PBL,i,j)
                write(90,*)dudz2(KS:KE_PBL,i,j)
                write(90,*)n2(KS:KE_PBL,i,j)
                write(90,*)pr(KS:KE_PBL,i,j)
                close(90)
                call abort
             end if
          end do
#endif
          do k = KS, KE_PBL
             tke_t(k,i,j) = ( max(tke_N(k,i,j), TKE_min) - tke(k,i,j) ) / dt
          end do
          do k = KE_PBL+1, KE
             tke_t(k,i,j) = 0.0_RP
          end do

       end do
       end do
    end do
    end do

#ifdef MORE_HIST
    adv(KE_PBL+1:KE,:,:) = 0.0_RP
    gen(KE_PBL+1:KE,:,:) = 0.0_RP
    dudz2(KE_PBL+1:KE,:,:) = 0.0_RP
    l(KE_PBL+1:KE,:,:) = 0.0_RP
    POTV(KE_PBL+1:KE,:,:) = 0.0_RP
    POTL(KE_PBL+1:KE,:,:) = 0.0_RP
    call HIST_in(adv, 'TKE_advc', 'advection of TKE', 'm2/s3', nohalo=.true.)
    call HIST_in(gen, 'TKE_gen', 'generation of TKE', 'm2/s3', nohalo=.true.)
    call HIST_in(dudz2, 'dUdZ2', 'dudz2', 'm2/s2', nohalo=.true.)
    call HIST_in(l, 'L_mix', 'minxing length', 'm', nohalo=.true.)
    call HIST_in(POTV, 'POTV', 'virtual potential temperature', 'K', nohalo=.true.)
    call HIST_in(POTL, 'POTL', 'liquid potential temperature', 'K', nohalo=.true.)
#endif

    return
  end subroutine ATMOS_PHY_TB_mynn

  subroutine get_length( &
       l, &
       DENS, &
       q, n2, &
       SFLX_MU, SFLX_MV, SFLX_PT, &
       PT0 )
    use scale_grid_real, only: &
       FZ => REAL_FZ
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       KARMAN => CONST_KARMAN, &
       CP     => CONST_CPdry, &
       EPS    => CONST_EPS
    implicit none

    real(RP), intent(out) :: l(KA,IA,JA)
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: q(KA,IA,JA)
    real(RP), intent(in) :: n2(KA,IA,JA)
    real(RP), intent(in) :: SFLX_MU(IA,JA)
    real(RP), intent(in) :: SFLX_MV(IA,JA)
    real(RP), intent(in) :: SFLX_PT(IA,JA)
    real(RP), intent(in) :: PT0(KA,IA,JA)

    real(RP) :: ls           !< L_S
    real(RP) :: lt           !< L_T
    real(RP) :: lb           !< L_B
    real(RP) :: rlm          !< 1/L_M
    real(RP) :: rlt          !< 1/L_T

    real(RP) :: qc           !< q_c
    real(RP) :: int_q        !< \int q dz
    real(RP) :: int_qz       !< \int qz dz
    real(RP) :: rn2sr         !< 1/N
    real(RP) :: us           !< friction velocity
    real(RP) :: zeta         !< height normalized by the Obukhov length

    real(RP) :: z
    real(RP) :: qdz

    real(RP) :: sw
    integer :: k, i, j

!OCL LOOP_NOFUSION,PREFETCH_SEQUENTIAL(SOFT),SWP
    do j = JS, JE
    do i = IS, IE
       int_qz = 0.0_RP
       int_q = 0.0_RP
       do k = KS, KE_PBL
          qdz = q(k,i,j) * ( FZ(k,i,j) - FZ(k-1,i,j) )
          int_qz = int_qz + ((FZ(k,i,j)+FZ(k-1,i,j))*0.5_RP-FZ(KS-1,i,j)) * qdz
          int_q  = int_q + qdz
       end do
       ! LT
       lt = min( max(0.23_RP * int_qz / (int_q + EPS), &
                    LT_min), &
                    ATMOS_PHY_TB_MYNN_Lt_MAX )
       rlt = 1.0_RP / lt

       us = ( SFLX_MU(i,j)**2 + SFLX_MV(i,j)**2 )**0.25_RP / DENS(KS,i,j) ! friction velocity
       us = max(us, Us_min)
       rlm = - KARMAN * GRAV * SFLX_PT(i,j) / (PT0(KS,i,j) * us**3 )

       qc = (GRAV/PT0(KS,i,j)*max(SFLX_PT(i,j),0.0_RP)*lt)**OneOverThree

       do k = KS, KE_PBL
          z = ( FZ(k,i,j)+FZ(k-1,i,j) )*0.5_RP - FZ(KS-1,i,j)
          zeta = z * rlm

          ! LS
          sw = sign(0.5_RP, zeta) + 0.5_RP ! 1 for zeta >= 0, 0 for zeta < 0
          ls = KARMAN * z &
             * ( sw / (1.0_RP + 2.7_RP*min(zeta,1.0_RP)*sw ) &
               + ( (1.0_RP - 100.0_RP*zeta)*(1.0_RP-sw) )**0.2_RP )

          ! LB
          sw  = sign(0.5_RP, n2(k,i,j)) + 0.5_RP ! 1 for dptdz >0, 0 for dptdz < 0
          rn2sr = 1.0_RP / ( sqrt(n2(k,i,j)*sw) + 1.0_RP-sw)
          lb = (1.0_RP + 5.0_RP * sqrt(qc*rn2sr/lt)) * q(k,i,j) * rn2sr * sw & ! qc=0 when SFLX_PT < 0
             +  999.E10_RP * (1.0_RP-sw)

          ! L
          l(k,i,j) = 1.0_RP / ( 1.0_RP/ls + rlt + 1.0_RP/lb )
       end do
    end do
    end do

    return
  end subroutine get_length

  subroutine get_q2_level2( &
       q2_2,        &
       dudz2, Ri, l )
    implicit none

    real(RP), intent(out) :: q2_2(KA,IA,JA)
    real(RP), intent(in)  :: dudz2(KA,IA,JA)
    real(RP), intent(in)  :: Ri(KA,IA,JA)
    real(RP), intent(in)  :: l(KA,IA,JA)

    real(RP) :: rf           !< Rf
    real(RP) :: sm_2         !< sm for level 2
    real(RP) :: sh_2         !< sh for level 2

    real(RP) :: q2
    integer :: k, i, j

!OCL LOOP_NOFUSION,PREFETCH_SEQUENTIAL(SOFT),SWP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE_PBL
       rf = min(0.5_RP / AF12 * ( Ri(k,i,j) &
                              + AF12*Rf1 &
                              - sqrt(Ri(k,i,j)**2 + 2.0_RP*AF12*(Rf1-2.0_RP*Rf2)*Ri(k,i,j) + (AF12*Rf1)**2) ), &
                Rfc)
       sh_2 = 3.0_RP * A2 * (G1+G2) * (Rfc-rf) / (1.0_RP-rf)
       sm_2 = sh_2 * AF12 * (Rf1-rf) / (Rf2-rf)
       q2 = B1 * l(k,i,j)**2 * sm_2 * (1.0_RP-rf) * dudz2(k,i,j)
       q2_2(k,i,j) = max( q2, 1.E-10_RP )
    end do
    end do
    end do

    return
  end subroutine get_q2_level2

  subroutine get_smsh( &
         sm, sh, & ! (out)
         q, q2_2, & ! (in)
         l, n2, dudz2 ) ! (in)
    implicit none
    real(RP), intent(out) :: sm(KA,IA,JA)
    real(RP), intent(out) :: sh(KA,IA,JA)
    real(RP), intent(in)  :: q(KA,IA,JA)
    real(RP), intent(in)  :: q2_2(KA,IA,JA)
    real(RP), intent(in)  :: l(KA,IA,JA)
    real(RP), intent(in)  :: n2(KA,IA,JA)
    real(RP), intent(in)  :: dudz2(KA,IA,JA)

    real(RP) :: l2q2         !< L^2/q^2
    real(RP) :: ac          !< \alpha_c
    real(RP) :: ac2          !< \alpha_c^2
    real(RP) :: p1           !< \Phi_1
    real(RP) :: p2           !< \Phi_2
    real(RP) :: p3           !< \Phi_3
    real(RP) :: p4           !< \Phi_4
    real(RP) :: p5           !< \Phi_5
    real(RP) :: rd25         !< 1/D_2.5
    real(RP) :: gh           !< G_H

    integer :: k, i, j

!OCL LOOP_NOFUSION,PREFETCH_SEQUENTIAL(SOFT),SWP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE_PBL

       ! level 2.5
       ac = min(q(k,i,j)/sqrt(q2_2(k,i,j)), 1.0_RP)
       ac2 = ac**2
       l2q2 = ( l(k,i,j) / q(k,i,j) )**2
       gh = - n2(k,i,j) * l2q2

       p1 = 1.0_RP - 3.0_RP * ac2 * A2 * B2 * (1.0_RP-C3) * gh
       p2 = 1.0_RP - 9.0_RP * ac2 * A1 * A2 * (1.0_RP-C2) * gh
       p3 = p1 + 9.0_RP * ac2 * A2**2 * (1.0_RP-C2) * (1.0_RP-C5) * gh
       p4 = p1 - 12.0_RP * ac2 * A1 * A2 * (1.0_RP-C2) * gh
       p5 = 6.0_RP * ac2 * A1**2 * dudz2(k,i,j) * l2q2

       rd25 = 1.0_RP / max(p2 * p4 + p5 * p3, 1.E-20_RP)
       sm(k,i,j) = max( ac * A1 * (p3 - 3.0_RP * C1 * p4) * rd25, 0.0_RP )
       sh(k,i,j) = max( ac * A2 * (p2 + 3.0_RP * C1 * p5) * rd25, 0.0_RP )

    end do
    end do
    end do

    return
  end subroutine get_smsh

  function erf(x)
    real(RP), parameter :: a = 0.1400122886866665_RP     ! -8(pi-3)/(3pi(pi-4))
    real(RP), parameter :: fourpi = 1.2732395447351628_RP ! 4/pi

    real(RP), intent(in) :: x
    real(RP) :: erf

    real(RP) :: x2

    x2 = x**2
    erf = sign( sqrt( 1.0_RP - exp(-x2 * (fourpi+a*x2)/(1.0_RP+a*x2) ) ), x )

    return
  end function erf

end module scale_atmos_phy_tb_mynn
