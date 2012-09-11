!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!          Smagolinsky-type
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-29 (S.Iga)       [new]
!! @li      2011-12-11 (H.Yashiro)   [mod] integrate to SCALE3
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-03-27 (H.Yashiro)   [mod] reconstruction
!! @li      2012-07-02 (S.Nishizawa) [mod] reconstruction with Brown et al. (1994)
!!
!! - Reference
!!  - Brown et al., 1994:
!!    Large-eddy simulaition of stable atmospheric boundary layers with a revised stochastic subgrid model.
!!    Roy. Meteor. Soc., 120, 1485-1512
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_tb
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_TB_setup
  public :: ATMOS_PHY_TB
  public :: ATMOS_PHY_TB_main

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'
  include 'inc_precision.h'

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
#ifdef DEBUG
  real(RP), private, parameter :: UNDEF = -9.999E30_RP
  integer,  private, parameter :: IUNDEF = -99999
#endif

  real(RP), private,      save :: nu_factC (KA,IA,JA) ! (Cs*Delta)^2 (cell center)
  real(RP), private,      save :: nu_factXY(KA,IA,JA) !              (x-y plane)
  real(RP), private,      save :: nu_factYZ(KA,IA,JA) !              (y-z plane)
  real(RP), private,      save :: nu_factZX(KA,IA,JA) !              (z-x plane)
  real(RP), private,      save :: nu_factZ (KA,IA,JA) !              (z edge)
  real(RP), private,      save :: nu_factX (KA,IA,JA) !              (x edge)
  real(RP), private,      save :: nu_factY (KA,IA,JA) !              (y edge)

  real(RP), private, parameter :: Cs  = 0.18_RP ! (Sullivan et al.1994, Nakanishi and Niino)
  real(RP), private, parameter :: PrN = 0.7_RP  ! Prandtl number in neutral conditions
  real(RP), private, parameter :: RiC = 0.25_RP ! critical Richardson number
  real(RP), private, parameter :: FmC = 16.0_RP ! fum = sqrt(1 - c*Ri)
  real(RP), private, parameter :: FhB = 40.0_RP ! fuh = sqrt(1 - b*Ri)/PrN
  real(RP), private            :: RPrN          ! 1 / PrN
  real(RP), private            :: RRiC          ! 1 / RiC
  real(RP), private            :: PrNovRiC      ! PrN / RiC

  real(RP), private, parameter :: OneOverThree = 1.0_RP / 3.0_RP
  real(RP), private, parameter :: twoOverThree = 2.0_RP / 3.0_RP

  integer, private, parameter :: ZDIR = 1
  integer, private, parameter :: XDIR = 2
  integer, private, parameter :: YDIR = 3

  !-----------------------------------------------------------------------------
contains

#ifdef DEBUG
  subroutine CHECK( line, v )
    integer,  intent(in) :: line
    real(RP), intent(in) :: v
    if ( v == UNDEF ) then
       write(*,*) "use uninitialized value at line ", line
       stop
    end if
  end subroutine CHECK
#endif

  subroutine ATMOS_PHY_TB_setup
    use mod_grid, only : &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY, &
       FDZ => GRID_FDZ, &
       FDX => GRID_FDX, &
       FDY => GRID_FDY
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-TB]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Smagorinsky-type Eddy Viscocity Model'

    RPrN     = 1.0_RP / PrN
    RRiC     = 1.0_RP / RiC
    PrNovRiC = (1- PrN) * RRiC

#ifdef DEBUG
    nu_factC (:,:,:) = UNDEF
    nu_factXY(:,:,:) = UNDEF
    nu_factYZ(:,:,:) = UNDEF
    nu_factZX(:,:,:) = UNDEF
    nu_factZ (:,:,:) = UNDEF
    nu_factX (:,:,:) = UNDEF
    nu_factY (:,:,:) = UNDEF
#endif
    do j = JS, JE+1
    do i = IS, IE+1
    do k = KS, KE
       nu_factC (k,i,j) = ( Cs * ( CDZ(k) * CDX(i) * CDY(j) )**OneOverThree )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       nu_factXY(k,i,j) = ( Cs * ( FDZ(k) * CDX(i) * CDY(j) )**OneOverThree )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS, KE
       nu_factYZ(k,i,j) = ( Cs * ( CDZ(k) * FDX(i) * CDY(j) )**OneOverThree )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS  , KE
       nu_factZX(k,i,j) = ( Cs * ( CDZ(k) * CDX(i) * FDY(j) )**OneOverThree )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE
    do i = IS-1, IE
    do k = KS  , KE
       nu_factZ(k,i,j) = ( Cs * ( CDZ(k) * FDX(i) * FDY(j) )**OneOverThree )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS  , KE
       nu_factX(k,i,j) = ( Cs * ( FDZ(k) * CDX(i) * FDY(j) )**OneOverThree )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS  , KE
       nu_factY(k,i,j) = ( Cs * ( FDZ(k) * FDX(i) * CDY(j) )**OneOverThree )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    return
  end subroutine ATMOS_PHY_TB_setup

  !-----------------------------------------------------------------------------
  !> Smagorinsky-type turblence
  !>
  !> comment:
  !>  1, Pr is given linearly (iga)
  !>  4, heat flux is not accurate yet. (i.e. energy is not conserved, see *1)
  !>  5, stratification effect is not considered.
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB
    use mod_time, only: &
       dttb => TIME_DTSEC_ATMOS_PHY_TB
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total,    &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_POTT, &
       SFLX_QV
    implicit none

    ! tendency
    real(RP) :: MOMZ_t(KA,IA,JA)
    real(RP) :: MOMX_t(KA,IA,JA)
    real(RP) :: MOMY_t(KA,IA,JA)
    real(RP) :: RHOT_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)

    ! diagnostic variables
    real(RP) :: tke(KA,IA,JA) ! TKE
    real(RP) :: nu (KA,IA,JA) ! eddy diffusion
    real(RP) :: Ri (KA,IA,JA) ! Richardoson number
    real(RP) :: Pr (KA,IA,JA) ! Prandtle number

    integer :: k, i, j, iq
    integer :: IIS, IIE, JJS, JJE


    call ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out) tendency
       tke, nu, Ri, Pr,                        & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)  diagnostic variables
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in) surface flux
       )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          MOMZ(k,i,j) = MOMZ(k,i,j) + dttb * MOMZ_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMX(k,i,j) = MOMX(k,i,j) + dttb * MOMX_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMY(k,i,j) = MOMY(k,i,j) + dttb * MOMY_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT(k,i,j) = RHOT(k,i,j) + dttb * RHOT_t(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    enddo
    enddo
#ifdef DEBUG
       IIS = IUNDEF; IIE = IUNDEF; JJS = IUNDEF; JJE = IUNDEF
#endif


    do iq = 1, QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + dttb * QTRC_t(k,i,j,iq)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          if ( QTRC(k,i,j,iq) < 1.0E-10_RP ) then
             QTRC(k,i,j,iq) = 0.0_RP
          endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    enddo
    enddo
#ifdef DEBUG
       IIS = IUNDEF; IIE = IUNDEF; JJS = IUNDEF; JJE = IUNDEF
#endif

    enddo ! do iq = 1, QA

    call ATMOS_vars_fillhalo

    call ATMOS_vars_total

    call HIST_in( MOMZ_t(:,:,:), 'MOMZ_t_tb', 'tendency of MOMZ in tb', 'kg/m2/s2',  '3D', dttb )
    call HIST_in( MOMX_t(:,:,:), 'MOMX_t_tb', 'tendency of MOMX in tb', 'kg/m2/s2',  '3D', dttb )
    call HIST_in( MOMY_t(:,:,:), 'MOMY_t_tb', 'tendency of MOMY in tb', 'kg/m2/s2',  '3D', dttb )
    call HIST_in( RHOT_t(:,:,:), 'RHOT_t_tb', 'tendency of RHOT in tb', 'K*kg/m3/s', '3D', dttb )
    do iq = 1, QA
       call HIST_in( QTRC_t(:,:,:,iq), AQ_NAME(iq)//'_t_tb', AQ_DESC(iq), AQ_UNIT(iq)//'/s', '3D', dttb )
    enddo

    call HIST_in( tke(:,:,:), 'TKE',  'turburent kinetic energy', 'm2/s2', '3D', dttb )
    call HIST_in( nu (:,:,:), 'NU',   'eddy viscosity',           'm2/s',  '3D', dttb )
    call HIST_in( Pr (:,:,:), 'Pr',   'Prantle number',           'NIL',   '3D', dttb )
    call HIST_in( Ri (:,:,:), 'Ri',   'Richardson number',        'NIL',   '3D', dttb )

  end subroutine ATMOS_PHY_TB


  subroutine ATMOS_PHY_TB_main( &
       MOMZ_t, MOMX_t, MOMY_t, RHOT_t, QTRC_t, & ! (out) tendency
       tke, nu_C, Ri, Pr,                      & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,     & ! (in)  diagnostic variables
       SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV & ! (in) surface flux
       )
    use mod_const, only : &
       GRAV => CONST_GRAV
    use mod_grid, only : &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    implicit none

    ! tendency
    real(RP), intent(out) :: MOMZ_t(KA,IA,JA)
    real(RP), intent(out) :: MOMX_t(KA,IA,JA)
    real(RP), intent(out) :: MOMY_t(KA,IA,JA)
    real(RP), intent(out) :: RHOT_t(KA,IA,JA)
    real(RP), intent(out) :: QTRC_t(KA,IA,JA,QA)

    real(RP), intent(out) :: tke (KA,IA,JA) ! TKE
    real(RP), intent(out) :: nu_C(KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Pr  (KA,IA,JA) ! Prantle number
    real(RP), intent(out) :: Ri  (KA,IA,JA) ! Richardson number

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)  :: SFLX_MOMZ(IA,JA)
    real(RP), intent(in)  :: SFLX_MOMX(IA,JA)
    real(RP), intent(in)  :: SFLX_MOMY(IA,JA)
    real(RP), intent(in)  :: SFLX_POTT(IA,JA)
    real(RP), intent(in)  :: SFLX_QV  (IA,JA)


    ! diagnostic variables
    real(RP) :: VELZ_C (KA,IA,JA)
    real(RP) :: VELZ_XY(KA,IA,JA)
    real(RP) :: VELX_C (KA,IA,JA)
    real(RP) :: VELX_YZ(KA,IA,JA)
    real(RP) :: VELY_C (KA,IA,JA)
    real(RP) :: VELY_ZX(KA,IA,JA)
    real(RP) :: POTT(KA,IA,JA)

    ! deformation rate tensor
    ! (cell center)
    real(RP) :: S33_C(KA,IA,JA)
    real(RP) :: S11_C(KA,IA,JA)
    real(RP) :: S22_C(KA,IA,JA)
    real(RP) :: S31_C(KA,IA,JA)
    real(RP) :: S12_C(KA,IA,JA)
    real(RP) :: S23_C(KA,IA,JA)
    ! (z edge or x-y plane)
    real(RP) :: S33_Z(KA,IA,JA)
    real(RP) :: S11_Z(KA,IA,JA)
    real(RP) :: S22_Z(KA,IA,JA)
    real(RP) :: S31_Z(KA,IA,JA)
    real(RP) :: S12_Z(KA,IA,JA)
    real(RP) :: S23_Z(KA,IA,JA)
    ! (x edge or y-z plane)
    real(RP) :: S33_X(KA,IA,JA)
    real(RP) :: S11_X(KA,IA,JA)
    real(RP) :: S22_X(KA,IA,JA)
    real(RP) :: S31_X(KA,IA,JA)
    real(RP) :: S12_X(KA,IA,JA)
    real(RP) :: S23_X(KA,IA,JA)
    ! (y edge or z-x plane)
    real(RP) :: S33_Y(KA,IA,JA)
    real(RP) :: S11_Y(KA,IA,JA)
    real(RP) :: S22_Y(KA,IA,JA)
    real(RP) :: S31_Y(KA,IA,JA)
    real(RP) :: S12_Y(KA,IA,JA)
    real(RP) :: S23_Y(KA,IA,JA)

    real(RP) :: nu_Z(KA,IA,JA)  ! eddy viscosity (z edge or x-y plane)
    real(RP) :: nu_X(KA,IA,JA)  !                (x edge or y-z plane)
    real(RP) :: nu_Y(KA,IA,JA)  !                (y edge or z-x plane)

    real(RP) :: S2(KA,IA,JA)     ! |S|^2
    real(RP) :: WORK_V(KA,IA,JA) ! work space (vertex)
    real(RP) :: WORK_Z(KA,IA,JA) !            (z edge or x-y plane)
    real(RP) :: WORK_X(KA,IA,JA) !            (x edge or y-z plane)
    real(RP) :: WORK_Y(KA,IA,JA) !            (y edge or z-x plane)

    real(RP) :: qflx_sgs(KA,IA,JA,3)

    real(RP) :: TMP1, TMP2, TMP3

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

#ifdef DEBUG
    VELZ_C (:,:,:) = UNDEF
    VELZ_XY(:,:,:) = UNDEF
    VELX_C (:,:,:) = UNDEF
    VELX_YZ(:,:,:) = UNDEF
    VELY_C (:,:,:) = UNDEF
    VELY_ZX(:,:,:) = UNDEF
    POTT(:,:,:) = UNDEF

    S33_C(:,:,:) = UNDEF
    S11_C(:,:,:) = UNDEF
    S22_C(:,:,:) = UNDEF
    S31_C(:,:,:) = UNDEF
    S12_C(:,:,:) = UNDEF
    S23_C(:,:,:) = UNDEF
    S33_Z(:,:,:) = UNDEF
    S11_Z(:,:,:) = UNDEF
    S22_Z(:,:,:) = UNDEF
    S31_Z(:,:,:) = UNDEF
    S12_Z(:,:,:) = UNDEF
    S23_Z(:,:,:) = UNDEF
    S33_X(:,:,:) = UNDEF
    S11_X(:,:,:) = UNDEF
    S22_X(:,:,:) = UNDEF
    S31_X(:,:,:) = UNDEF
    S12_X(:,:,:) = UNDEF
    S23_X(:,:,:) = UNDEF
    S33_Y(:,:,:) = UNDEF
    S11_Y(:,:,:) = UNDEF
    S22_Y(:,:,:) = UNDEF
    S31_Y(:,:,:) = UNDEF
    S12_Y(:,:,:) = UNDEF
    S23_Y(:,:,:) = UNDEF

    WORK_V(:,:,:) = UNDEF
    WORK_Z(:,:,:) = UNDEF
    WORK_X(:,:,:) = UNDEF
    WORK_Y(:,:,:) = UNDEF

    S2(:,:,:) = UNDEF

    nu_C(:,:,:) = UNDEF
    nu_Z(:,:,:) = UNDEF
    nu_X(:,:,:) = UNDEF
    nu_Y(:,:,:) = UNDEF

    tke (:,:,:) = UNDEF
    Pr  (:,:,:) = UNDEF
    Ri  (:,:,:) = UNDEF

    qflx_sgs(:,:,:,:) = UNDEF
#endif


    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: SGS Parameterization'

   ! momentum -> velocity
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMZ(k,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELZ_XY(k,i,j) = 2.0_RP * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+2
    do i = IS-1, IE+2
    do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMZ(k,i,j) )
       call CHECK( __LINE__, MOMZ(k-1,i,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELZ_C(k,i,j) = 0.5_RP * ( MOMZ(k,i,j) + MOMZ(k-1,i,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+2
    do i = IS-1, IE+2
#ifdef DEBUG
       call CHECK( __LINE__, MOMZ(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i,j) )
#endif
       VELZ_C(KS,i,j) = 0.5_RP * MOMZ(KS,i,j) / DENS(KS,i,j) ! MOMZ(KS-1,i,j) = 0
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMX(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELX_YZ(k,i,j) = 2.0_RP * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       VELX_YZ(KE+1,i,j) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+2
    do i = IS-1, IE+2
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMX(k,i,j) )
       call CHECK( __LINE__, MOMX(k,i-1,j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELX_C(k,i,j) = 0.5_RP * ( MOMX(k,i,j) + MOMX(k,i-1,j) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMY(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELY_ZX(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+1
    do i = IS-1, IE+1
       VELY_ZX(KE+1,i,j) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+2
    do i = IS-1, IE+2
    do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, MOMY(k,i,j) )
       call CHECK( __LINE__, MOMY(k,i,j-1) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
       VELY_C(k,i,j) = 0.5_RP * ( MOMY(k,i,j) + MOMY(k,i,j-1) ) / DENS(k,i,j)
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    ! potential temperature
    do j = JS-1, JE+1
    do i = IS-1, IE+1
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
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    !##### Start Upadate #####

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF; WORK_V(:,:,:) = UNDEF
#endif
       ! w
       ! (x-y plane)
       ! WORK_Z = VELZ_XY
       ! (y-z plane)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i+1,j) )
       call CHECK( __LINE__, VELZ_C(k,i,j) )
#endif
          WORK_X(k,i,j) = 0.5_RP * ( VELZ_C(k,i+1,j) + VELZ_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          WORK_X(KE+1,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i,j+1) )
       call CHECK( __LINE__, VELZ_C(k,i,j) )
#endif
          WORK_Y(k,i,j) = 0.5_RP * ( VELZ_C(k,i,j+1) + VELZ_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          WORK_Y(KE+1,i,j) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (vertex)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k,i+1,j) )
       call CHECK( __LINE__, VELZ_XY(k,i,j+1) )
       call CHECK( __LINE__, VELZ_XY(k,i+1,j+1) )
#endif
          WORK_V(k,i,j) = 0.25_RP * ( VELZ_XY(k,i,j) + VELZ_XY(k,i+1,j) + VELZ_XY(k,i,j+1) + VELZ_XY(k,i+1,j+1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! dw/dz
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S33_C(k,i,j) = ( VELZ_XY(k,i,j) - VELZ_XY(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(KS,i,j) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S33_C(KS,i,j) = VELZ_XY(KS,i,j) * RCDZ(KS) ! VELZ_XY(KS-1,i,j) == 0
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S33_Z(k,i,j) = ( WORK_V(k,i,j) - WORK_V(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_V(KS,i,j) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S33_Z(KS,i,j) = WORK_V(KS,i,j) * RCDZ(KS) ! WORK_V(KS-1,i,j) == 0
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(k+1,i,j) )
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S33_X(k,i,j) = ( WORK_Y(k+1,i,j) - WORK_Y(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k+1,i,j) )
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S33_Y(k,i,j) = ( WORK_X(k+1,i,j) - WORK_X(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dx
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i+1,j) )
       call CHECK( __LINE__, VELZ_C(k,i-1,j) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S31_C(k,i,j) = 0.5_RP * ( VELZ_C(k,i+1,j) - VELZ_C(k,i-1,j) ) / ( FDX(i) + FDX(i-1) )
       enddo
       enddo
       enddo
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(k,i+1,j) )
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S31_Z(k,i,j) = 0.5_RP * ( WORK_Y(k,i+1,j) - WORK_Y(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S31_X(k,i,j) = 0.5_RP * ( WORK_V(k,i,j) - WORK_V(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i+1,j) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S31_Y(k,i,j) = 0.5_RP * ( VELZ_XY(k,i+1,j) - VELZ_XY(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dy
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i,j+1) )
       call CHECK( __LINE__, VELZ_C(k,i,j-1) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S23_C(k,i,j) = 0.5_RP * ( VELZ_C(k,i,j+1) - VELZ_C(k,i,j-1) ) / ( FDY(j) + FDY(j-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j+1) )
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S23_Z(k,i,j) = 0.5_RP * ( WORK_X(k,i,j+1) - WORK_X(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j+1) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S23_X(k,i,j) = 0.5_RP * ( VELZ_XY(k,i,j+1) - VELZ_XY(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S23_Y(k,i,j) = 0.5_RP * ( WORK_V(k,i,j) - WORK_V(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF; WORK_V(:,:,:) = UNDEF
#endif
       ! u
       ! (x-y plane)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
#endif
          WORK_Z(k,i,j) = 0.5_RP * ( VELX_C(k+1,i,j) + VELX_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       ! WORK_X = VELX_YZ
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
#endif
          WORK_Y(k,i,j) = 0.5_RP * ( VELX_C(k,i,j+1) + VELX_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (vertex)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j+1) )
#endif
          WORK_V(k,i,j) = 0.25_RP * ( VELX_YZ(k,i,j) + VELX_YZ(k,i,j+1) + VELX_YZ(k+1,i,j) + VELX_YZ(k+1,i,j+1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! du/dx
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_C(k,i,j) = ( VELX_YZ(k,i,j) - VELX_YZ(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(k,i+1,j) )
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S11_Z(k,i,j) = ( WORK_Y(k,i+1,j) - WORK_Y(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_X(k,i,j) = ( WORK_V(k,i,j) - WORK_V(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i+1,j) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S11_Y(k,i,j) = ( WORK_Z(k,i+1,j) - WORK_Z(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dz
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_C(k,i,j) )
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S31_C(k,i,j) = S31_C(k,i,j) + &
               0.5_RP * ( VELX_C(k+1,i,j) - VELX_C(k-1,i,j) ) / ( FDZ(k) + FDZ(k-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S31_C(KS,i,j) )
       call CHECK( __LINE__, VELX_C(KS+1,i,j) )
       call CHECK( __LINE__, VELX_C(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
#endif
          S31_C(KS,i,j) = S31_C(KS,i,j) + &
               0.5_RP * ( VELX_C(KS+1,i,j) - VELX_C(KS,i,j) ) * RFDZ(KS)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S31_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
#endif
          S31_C(KE,i,j) = S31_C(KE,i,j) + &
               0.5_RP * ( VELX_C(KE,i,j) - VELX_C(KE-1,i,j) ) * RFDZ(KE-1)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_Z(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S31_Z(k,i,j) = S31_Z(k,i,j) + &
               0.5_RP * ( WORK_V(k,i,j) - WORK_V(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S31_Z(KS,i,j) )
       call CHECK( __LINE__, WORK_V(KS,i,j) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S31_Z(KS,i,j) = S31_Z(KS,i,j) + &
               0.5_RP * WORK_V(KS,i,j) * RCDZ(KS) ! WORK_V(KS-1,i,j) == 0
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S31_Z(KE,i,j) )
       call CHECK( __LINE__, WORK_V(KE-1,i,j) )
       call CHECK( __LINE__, RCDZ(KE) )
#endif
          S31_Z(KE,i,j) = S31_Z(KE,i,j) - &
               0.5_RP * WORK_V(KE-1,i,j) * RCDZ(KE) ! WORK_V(KE,i,j) == 0
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_X(k,i,j) )
       call CHECK( __LINE__, VELX_C(k+1,i,j+1) )
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S31_X(k,i,j) = S31_X(k,i,j) + &
               0.25_RP * ( VELX_C(k+1,i,j+1) + VELX_C(k+1,i,j) - VELX_C(k,i,j+1) - VELX_C(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_Y(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S31_Y(k,i,j) = S31_Y(k,i,j) + &
               0.5_RP * ( VELX_YZ(k+1,i,j) - VELX_YZ(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dy
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k,i,j-1) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S12_C(k,i,j) = 0.5_RP * ( VELX_C(k,i,j+1) - VELX_C(k,i,j-1) ) / ( FDY(j) + FDY(j-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S12_Z(k,i,j) = 0.5_RP * ( VELX_YZ(k,i,j+1) - VELX_YZ(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j+1) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S12_X(k,i,j) = 0.5_RP * ( WORK_Z(k,i,j+1) - WORK_Z(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j-1) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S12_Y(k,i,j) = 0.5_RP * ( WORK_V(k,i,j) - WORK_V(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF; WORK_V(:,:,:) = UNDEF
#endif
       ! v
       ! (x-y plane)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_C(k+1,i,j) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
#endif
          WORK_Z(k,i,j) = 0.5_RP * ( VELY_C(k+1,i,j) + VELY_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELY_C(k,i+1,j) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
#endif
          WORK_X(k,i,j) = 0.5_RP * ( VELY_C(k,i+1,j) + VELY_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       ! WORK_Y = VELY_ZX
       ! (vertex)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i+1,j) )
#endif
          WORK_V(k,i,j) = 0.25_RP * ( VELY_ZX(k,i,j) + VELY_ZX(k+1,i,j) + VELY_ZX(k,i+1,j) + VELY_ZX(k+1,i+1,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! dv/dy
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S22_C(k,i,j) = ( VELY_ZX(k,i,j) - VELY_ZX(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j+1) )
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S22_Z(k,i,j) = ( WORK_X(k,i,j+1) - WORK_X(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j+1) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S22_X(k,i,j) = ( WORK_Z(k,i,j+1) - WORK_Z(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S22_Y(k,i,j) = ( WORK_V(k,i,j) - WORK_V(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! 1/2 * dv/dx
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S12_C(k,i,j) )
       call CHECK( __LINE__, VELY_C(k,i-1,j) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S12_C(k,i,j) = S12_C(k,i,j) + &
               0.5_RP * ( VELY_C(k,i+1,j) - VELY_C(k,i-1,j) ) / ( FDX(i) + FDX(i-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S12_Z(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S12_Z(k,i,j) = S12_Z(k,i,j) + &
               0.5_RP * ( VELY_ZX(k,i+1,j) - VELY_ZX(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S12_X(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S12_X(k,i,j) = S12_X(k,i,j) + &
               0.5_RP * ( WORK_V(k,i,j) - WORK_V(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S12_Y(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k,i+1,j) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S12_Y(k,i,j) = S12_Y(k,i,j) + &
               0.5_RP * ( WORK_Z(k,i+1,j) - WORK_Z(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dv/dz
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_C(k,i,j) )
       call CHECK( __LINE__, VELY_C(k+1,i,j) )
       call CHECK( __LINE__, VELY_C(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S23_C(k,i,j) = S23_C(k,i,j) + &
               0.5_RP * ( VELY_C(k+1,i,j) - VELY_C(k-1,i,j) ) / ( FDZ(k) + FDZ(k-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S23_C(KS,i,j) )
       call CHECK( __LINE__, VELY_C(KS+1,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
#endif
          S23_C(KS,i,j) = S23_C(KS,i,j) + &
               0.5_RP * ( VELY_C(KS+1,i,j) - VELY_C(KS,i,j) ) * RFDZ(KS)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S23_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
#endif
          S23_C(KE,i,j) = S23_C(KE,i,j) + &
               0.5_RP * ( VELY_C(KE,i,j) - VELY_C(KE-1,i,j) ) * RFDZ(KE-1)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_Z(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S23_Z(k,i,j) = S23_Z(k,i,j) + &
               0.5_RP * ( WORK_V(k,i,j) - WORK_V(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S23_Z(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i+1,j) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S23_Z(KS,i,j) = S23_Z(KS,i,j) + &
               0.25_RP * ( VELY_ZX(KS+1,i,j) + VELY_ZX(KS+1,i+1,j) &
                         - VELY_ZX(KS  ,i,j) - VELY_ZX(KS  ,i+1,j) ) * RCDZ(KS)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S23_Z(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i+1,j) )
       call CHECK( __LINE__, RCDZ(KE-1) )
#endif
          S23_Z(KE,i,j) = S23_Z(KE,i,j) + &
               0.25_RP * ( VELY_ZX(KE  ,i,j) + VELY_ZX(KE  ,i+1,j) &
                         - VELY_ZX(KE-1,i,j) - VELY_ZX(KE-1,i+1,j) ) * RCDZ(KE-1)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_X(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S23_X(k,i,j) = S23_X(k,i,j) + &
               0.5_RP * ( VELY_ZX(k+1,i,j) - VELY_ZX(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_Y(k,i,j) )
       call CHECK( __LINE__, WORK_X(k+1,i,j) )
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S23_Y(k,i,j) = S23_Y(k,i,j) + &
               0.5_RP * ( WORK_X(k+1,i,j) - WORK_X(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


       ! nu_SGS = (Cs * Delta)^2 * |S|, |S|^2 = 2*Sij*Sij
       ! tke = (Cs * Delta)^2 * |S|^2
#ifdef DEBUG
       S2(:,:,:) = UNDEF
       nu_C(:,:,:) = UNDEF; nu_Z(:,:,:) = UNDEF; nu_X(:,:,:) = UNDEF; nu_Y(:,:,:) = UNDEF
       tke(:,:,:) = UNDEF
       Ri(:,:,:) = UNDEF
       Pr(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (cell center)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, S31_C(k,i,j) )
       call CHECK( __LINE__, S12_C(k,i,j) )
       call CHECK( __LINE__, S23_C(k,i,j) )
#endif
          S2(k,i,j) = &
                 2.0_RP * ( S11_C(k,i,j)**2 + S22_C(k,i,j)**2 + S33_C(k,i,j)**2 ) &
               + 4.0_RP * ( S31_C(k,i,j)**2 + S12_C(k,i,j)**2 + S23_C(k,i,j)**2 )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Ri = N^2 / |S|^2, N^2 = g / theta * dtheta/dz
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          Ri(k,i,j) = GRAV * ( POTT(k+1,i,j) - POTT(k-1,i,j) ) &
               / ( ( FDZ(k) + FDZ(k-1) ) * POTT(k,i,j) * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KS+1,i,j) )
       call CHECK( __LINE__, POTT(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S2(KS,i,j) )
#endif
          Ri(KS,i,j) = GRAV * ( POTT(KS+1,i,j) - POTT(KS,i,j) ) &
               * RFDZ(KS) / (POTT(KS,i,j) * max(S2(KS,i,j),1.0E-20_RP) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KE,i,j) )
       call CHECK( __LINE__, POTT(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
       call CHECK( __LINE__, S2(KE,i,j) )
#endif
          Ri(KE,i,j) = GRAV * ( POTT(KE,i,j) - POTT(KE-1,i,j) ) &
               * RFDZ(KE-1) / (POTT(KE,i,j) * max(S2(KE,i,j),1.0E-20_RP) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, Ri(k,i,j) )
       call CHECK( __LINE__, nu_factC(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( Ri(k,i,j) < 0.0_RP ) then ! stable
             nu_C(k,i,j) = nu_factC(k,i,j) &
                  * sqrt( S2(k,i,j) * (1.0_RP - FmC*Ri(k,i,j)) )
          else if ( Ri(k,i,j) < RiC ) then ! weakly stable
             nu_C(k,i,j) = nu_factC(k,i,j) &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - Ri(k,i,j)*RRiC )**4
          else ! strongly stable
             nu_C(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! tke = 1/2 (tau_11 + tau_22 + tau_33) = (Cs * Delta)^2 * |S|^2
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, nu_factC(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          tke(k,i,j) = nu_factC(k,i,j) * S2(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Pr = nu_m / nu_h = fm / fh
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, Ri(k,i,j) )
#endif
          if ( Ri(k,i,j) < 0.0_RP ) then ! stable
             Pr(k,i,j) = sqrt( ( 1.0_RP - FmC*Ri(k,i,j) )  &
                             / ( 1.0_RP - FhB*Ri(k,i,j) ) ) * PrN
          else if ( Ri(k,i,j) < RiC ) then ! weakly stable
             Pr(k,i,j) = PrN / ( 1.0_RP - PrNovRiC * Ri(k,i,j) )
          else ! strongly stable
             Pr(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

#ifdef DEBUG
       S2(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S11_Z(k,i,j) )
       call CHECK( __LINE__, S22_Z(k,i,j) )
       call CHECK( __LINE__, S33_Z(k,i,j) )
       call CHECK( __LINE__, S31_Z(k,i,j) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
       call CHECK( __LINE__, S23_Z(k,i,j) )
#endif
          S2(k,i,j) = &
                 2.0_RP * ( S11_Z(k,i,j)**2 + S22_Z(k,i,j)**2 + S33_Z(k,i,j)**2 ) &
               + 4.0_RP * ( S31_Z(k,i,j)**2 + S12_Z(k,i,j)**2 + S23_Z(k,i,j)**2 )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Ri
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k+1,i+1,j) )
       call CHECK( __LINE__, POTT(k+1,i,j+1) )
       call CHECK( __LINE__, POTT(k+1,i+1,j+1) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k,i+1,j) )
       call CHECK( __LINE__, POTT(k,i,j+1) )
       call CHECK( __LINE__, POTT(k,i+1,j+1) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, POTT(k-1,i+1,j) )
       call CHECK( __LINE__, POTT(k-1,i,j+1) )
       call CHECK( __LINE__, POTT(k-1,i+1,j+1) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          TMP1 = ( POTT(k+1,i,j) + POTT(k+1,i+1,j) + POTT(k+1,i,j+1) + POTT(k+1,i+1,j+1) ) * 0.25_RP
          TMP2 = ( POTT(k  ,i,j) + POTT(k  ,i+1,j) + POTT(k  ,i,j+1) + POTT(k  ,i+1,j+1) ) * 0.25_RP
          TMP3 = ( POTT(k-1,i,j) + POTT(k-1,i+1,j) + POTT(k-1,i,j+1) + POTT(k-1,i+1,j+1) ) * 0.25_RP
          WORK_Z(k,i,j) = GRAV * ( TMP1 - TMP3 ) &
               / ( ( FDZ(k) + FDZ(k-1) ) * TMP2 * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KE,i,j) )
       call CHECK( __LINE__, POTT(KE,i+1,j) )
       call CHECK( __LINE__, POTT(KE,i,j+1) )
       call CHECK( __LINE__, POTT(KE,i+1,j+1) )
       call CHECK( __LINE__, POTT(KE-1,i,j) )
       call CHECK( __LINE__, POTT(KE-1,i+1,j) )
       call CHECK( __LINE__, POTT(KE-1,i,j+1) )
       call CHECK( __LINE__, POTT(KE-1,i+1,j+1) )
       call CHECK( __LINE__, RFDZ(KE-1) )
       call CHECK( __LINE__, S2(KE,i,j) )
#endif
          TMP2 = ( POTT(KE  ,i,j) + POTT(KE  ,i+1,j) + POTT(KE  ,i,j+1) + POTT(KE  ,i+1,j+1) ) * 0.25_RP
          TMP3 = ( POTT(KE-1,i,j) + POTT(KE-1,i+1,j) + POTT(KE-1,i,j+1) + POTT(KE-1,i+1,j+1) ) * 0.25_RP
          WORK_Z(KE,i,j) = 0.5_RP * GRAV * ( TMP2 - TMP3 ) &
               * RFDZ(KE-1) / ( TMP2 * max(S2(KE,i,j),1.0E-20_RP) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KS+1,i,j) )
       call CHECK( __LINE__, POTT(KS+1,i+1,j) )
       call CHECK( __LINE__, POTT(KS+1,i,j+1) )
       call CHECK( __LINE__, POTT(KS+1,i+1,j+1) )
       call CHECK( __LINE__, POTT(KS,i,j) )
       call CHECK( __LINE__, POTT(KS,i+1,j) )
       call CHECK( __LINE__, POTT(KS,i,j+1) )
       call CHECK( __LINE__, POTT(KS,i+1,j+1) )
       call CHECK( __LINE__, RFDZ(KS-1) )
       call CHECK( __LINE__, S2(KS,i,j) )
#endif
          TMP1 = ( POTT(KS+1,i,j) + POTT(KS+1,i+1,j) + POTT(KS+1,i,j+1) + POTT(KS+1,i+1,j+1) ) * 0.25_RP
          TMP2 = ( POTT(KS  ,i,j) + POTT(KS  ,i+1,j) + POTT(KS  ,i,j+1) + POTT(KS  ,i+1,j+1) ) * 0.25_RP
          WORK_Z(KS,i,j) = 0.5_RP * GRAV * ( TMP1 - TMP2 ) &
               * RFDZ(KS) / ( TMP2 * max(S2(KS,i,j),1.0E-20_RP) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, nu_factZ(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( WORK_Z(k,i,j) < 0.0_RP ) then
             nu_Z(k,i,j) = nu_factZ(k,i,j) &
                  * sqrt( S2(k,i,j) * (1.0_RP - FmC*WORK_Z(k,i,j)) )
          else if ( WORK_Z(k,i,j) < RiC ) then
             nu_Z(k,i,j) = nu_factZ(k,i,j) &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - WORK_Z(k,i,j)*RRiC )**4
          else
             nu_Z(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
#ifdef DEBUG
       S2(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S11_X(k,i,j) )
       call CHECK( __LINE__, S22_X(k,i,j) )
       call CHECK( __LINE__, S33_X(k,i,j) )
       call CHECK( __LINE__, S31_X(k,i,j) )
       call CHECK( __LINE__, S12_X(k,i,j) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          S2(k,i,j) = &
                 2.0_RP * ( S11_X(k,i,j)**2 + S22_X(k,i,j)**2 + S33_X(k,i,j)**2 ) &
               + 4.0_RP * ( S31_X(k,i,j)**2 + S12_X(k,i,j)**2 + S23_X(k,i,j)**2 )
       enddo
       enddo
       enddo

#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Ri
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k+1,i,j+1) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k,i,j+1) )
       call CHECK( __LINE__, RFDZ(k) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          TMP1 = ( POTT(k+1,i,j) + POTT(k+1,i,j+1) ) * 0.5_RP
          TMP2 = ( POTT(k  ,i,j) + POTT(k  ,i,j+1) ) * 0.5_RP
          WORK_X(k,i,j) = 0.5_RP * GRAV * ( TMP1 - TMP2 ) &
               * RFDZ(k) / ( ( TMP1 + TMP2 ) * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, nu_factX(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( WORK_X(k,i,j) < 0.0_RP ) then
             nu_X(k,i,j) = nu_factX(k,i,j) &
                  * sqrt( S2(k,i,j) * (1.0_RP - FmC*WORK_X(k,i,j)) )
          else if ( WORK_X(k,i,j) < RiC ) then
             nu_X(k,i,j) = nu_factX(k,i,j) &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - WORK_X(k,i,j)*RRiC )**4
          else
             nu_X(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
#ifdef DEBUG
       S2(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S11_Y(k,i,j) )
       call CHECK( __LINE__, S22_Y(k,i,j) )
       call CHECK( __LINE__, S33_Y(k,i,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
       call CHECK( __LINE__, S12_Y(k,i,j) )
       call CHECK( __LINE__, S23_Y(k,i,j) )
#endif
          S2(k,i,j) = &
                 2.0_RP * ( S11_Y(k,i,j)**2 + S22_Y(k,i,j)**2 + S33_Y(k,i,j)**2 ) &
               + 4.0_RP * ( S31_Y(k,i,j)**2 + S12_Y(k,i,j)**2 + S23_Y(k,i,j)**2 )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Ri
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k+1,i+1,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k,i+1,j) )
       call CHECK( __LINE__, RFDZ(k) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          TMP1 = ( POTT(k+1,i,j) + POTT(k+1,i+1,j) ) * 0.5_RP
          TMP2 = ( POTT(k  ,i,j) + POTT(k  ,i+1,j) ) * 0.5_RP
          WORK_Y(k,i,j) = 0.5_RP * GRAV * ( TMP1 - TMP2 ) &
               * RFDZ(k) / ( ( TMP1 + TMP2 ) * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, nu_factY(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( WORK_Y(k,i,j) < 0.0_RP ) then
             nu_Y(k,i,j) = nu_factY(k,i,j) &
                  * sqrt( S2(k,i,j) * (1.0_RP - FmC*WORK_Y(k,i,j)) )
          else if ( WORK_Y(k,i,j) < RiC ) then
             nu_Y(k,i,j) = nu_factY(k,i,j) &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - WORK_Y(k,i,j)*RRiC )**4
          else
             nu_Y(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


       !##### momentum equation (z) #####
       ! (cell center)
#ifdef DEBUG
       qflx_sgs(:,:,:,:) = UNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, nu_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, tke(k,i,j) )
#endif
          qflx_sgs(k,i,j,ZDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu_C(k,i,j) &
               * ( S33_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) &
             + twoOverThree * tke(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
       call CHECK( __LINE__, SFLX_MOMZ(i,j) )
#endif
          qflx_sgs(KS,i,j,ZDIR) = SFLX_MOMZ(i,j) ! bottom boundary
          qflx_sgs(KE,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i+1,j) )
       call CHECK( __LINE__, nu_Y(k,i,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
          qflx_sgs(k,i,j,XDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ) &
                               * nu_Y(k,i,j) * S31_Y(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j+1) )
       call CHECK( __LINE__, nu_X(k,i,j) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          qflx_sgs(k,i,j,YDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) &
                               * nu_X(k,i,j) * S23_X(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !--- tendency momentum(z)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, qflx_sgs(k+1,i,j,ZDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,ZDIR) )
       call CHECK( __LINE__, RFDZ(k) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,XDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i-1,j,XDIR) )
       call CHECK( __LINE__, RCDX(i) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,YDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i,j-1,YDIR) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          MOMZ_t(k,i,j) = - ( ( qflx_sgs(k+1,i,j,ZDIR)-qflx_sgs(k,i  ,j  ,ZDIR) ) * RFDZ(k) &
                            + ( qflx_sgs(k  ,i,j,XDIR)-qflx_sgs(k,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_sgs(k  ,i,j,YDIR)-qflx_sgs(k,i  ,j-1,YDIR) ) * RCDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


#ifdef DEBUG
       qflx_sgs(:,:,:,:) = UNDEF
#endif
       !##### momentum equation (x) #####
       ! (y edge)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i+1,j) )
       call CHECK( __LINE__, nu_Y(k,i,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
          qflx_sgs(k,i,j,ZDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ) &
                               * nu_Y(k,i,j) * S31_Y(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
       call CHECK( __LINE__, SFLX_MOMX(i,j) )
#endif
          qflx_sgs(KS-1,i,j,ZDIR) = SFLX_MOMX(i,j) ! bottom boundary
          qflx_sgs(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (cell center)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, nu_C(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, tke(k,i,j) )
#endif
          qflx_sgs(k,i,j,XDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu_C(k,i,j) &
               * ( S11_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) &
             + twoOverThree * tke(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i+1,j+1) )
       call CHECK( __LINE__, nu_Z(k,i,j) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs(k,i,j,YDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
                               * nu_Z(k,i,j) * S12_Z(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !--- tendency momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, qflx_sgs(k,i,j,ZDIR) )
       call CHECK( __LINE__, qflx_sgs(k-1,i,j,ZDIR) )
       call CHECK( __LINE__, RCDZ(k) )
       call CHECK( __LINE__, qflx_sgs(k,i+1,j,XDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,XDIR) )
       call CHECK( __LINE__, RFDX(i) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,YDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i,j-1,YDIR) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          MOMX_t(k,i,j) = - ( ( qflx_sgs(k,i  ,j,ZDIR)-qflx_sgs(k-1,i,j,  ZDIR) ) * RCDZ(k) &
                            + ( qflx_sgs(k,i+1,j,XDIR)-qflx_sgs(k  ,i,j,  XDIR) ) * RFDX(i) &
                            + ( qflx_sgs(k,i  ,j,YDIR)-qflx_sgs(k  ,i,j-1,YDIR) ) * RCDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


#ifdef DEBUG
       qflx_sgs(:,:,:,:) = UNDEF
#endif
       !##### momentum equation (y) #####
       ! (x edge)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j+1) )
       call CHECK( __LINE__, nu_X(k,i,j) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          qflx_sgs(k,i,j,ZDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) &
                               * nu_X(k,i,j) * S23_X(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
       call CHECK( __LINE__, SFLX_MOMY(i,j) )
#endif
          qflx_sgs(KS-1,i,j,ZDIR) = SFLX_MOMY(i,j) ! bottom boundary
          qflx_sgs(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z edge)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, DENS(k,i+1,j+1) )
       call CHECK( __LINE__, nu_Z(k,i,j) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs(k,i,j,XDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
                               * nu_Z(k,i,j) * S12_Z(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z-x plane)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, nu_C(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, tke(k,i,j) )
#endif
          qflx_sgs(k,i,j,YDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu_C(k,i,j) &
               * ( S22_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) &
             + twoOverThree * tke(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !--- tendency momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, qflx_sgs(k,i,j,ZDIR) )
       call CHECK( __LINE__, qflx_sgs(k-1,i,j,ZDIR) )
       call CHECK( __LINE__, RCDZ(k) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,XDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i-1,j,XDIR) )
       call CHECK( __LINE__, RCDX(i) )
       call CHECK( __LINE__, qflx_sgs(k,i,j+1,YDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,YDIR) )
       call CHECK( __LINE__, RFDZ(j) )
#endif
          MOMY_t(k,i,j) = - ( ( qflx_sgs(k,i,j  ,ZDIR)-qflx_sgs(k-1,i  ,j,ZDIR) ) * RCDZ(k) &
                            + ( qflx_sgs(k,i,j  ,XDIR)-qflx_sgs(k  ,i-1,j,XDIR) ) * RCDX(i) &
                            + ( qflx_sgs(k,i,j+1,YDIR)-qflx_sgs(k  ,i  ,j,YDIR) ) * RFDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !##### Thermodynamic Equation #####

#ifdef DEBUG
       S33_Z(:,:,:) = UNDEF; S33_X(:,:,:) = UNDEF; S33_Y(:,:,:) = UNDEF
       S22_Z(:,:,:) = UNDEF; S22_X(:,:,:) = UNDEF; S22_Y(:,:,:) = UNDEF
       S11_Z(:,:,:) = UNDEF; S11_X(:,:,:) = UNDEF; S11_Y(:,:,:) = UNDEF
       S31_Z(:,:,:) = UNDEF; S31_X(:,:,:) = UNDEF; S31_Y(:,:,:) = UNDEF
       S12_Z(:,:,:) = UNDEF; S12_X(:,:,:) = UNDEF; S12_Y(:,:,:) = UNDEF
       S23_Z(:,:,:) = UNDEF; S23_X(:,:,:) = UNDEF; S23_Y(:,:,:) = UNDEF
#endif


#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! w
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i,j) )
       call CHECK( __LINE__, VELZ_C(k,i+1,j) )
       call CHECK( __LINE__, VELZ_C(k,i,j+1) )
       call CHECK( __LINE__, VELZ_C(k,i+1,j+1) )
#endif
          WORK_Z(k,i,j) = 0.25_RP * ( VELZ_C(k,i,j) + VELZ_C(k,i+1,j) + VELZ_C(k,i,j+1) + VELZ_C(k,i+1,j+1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k,i,j+1) )
#endif
          WORK_X(k,i,j) = 0.5_RP * ( VELZ_XY(k,i,j) + VELZ_XY(k,i,j+1) )
       enddo
       enddo
       enddo
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k,i+1,j) )
#endif
          WORK_Y(k,i,j) = 0.5_RP * ( VELZ_XY(k,i,j) + VELZ_XY(k,i+1,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! dw/dz
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k+1,i,j) )
       call CHECK( __LINE__, VELZ_C(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S33_Z(k,i,j) = ( VELZ_C(k+1,i,j) - VELZ_C(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, WORK_Y(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S33_X(k,i,j) = ( WORK_Y(k,i,j) - WORK_Y(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(KS,i,j) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S33_X(KS,i,j) = WORK_Y(KS,i,j) * RCDZ(KS) ! WORK_Y(k-1,i,j) = 0
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, WORK_X(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S33_Y(k,i,j) = ( WORK_X(k,i,j) - WORK_X(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS  , IIE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(KS,i,j) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S33_Y(KS,i,j) = WORK_X(KS,i,j) * RCDZ(KS) ! WORK_Z(KS-1,i,j) = 0
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dx
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i+1,j) )
       call CHECK( __LINE__, VELZ_XY(k,i-1,j) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S31_Z(k,i,j) = 0.5_RP * ( VELZ_XY(k,i+1,j) - VELZ_XY(k,i-1,j) ) / ( FDX(i) + FDX(i-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i+1,j) )
       call CHECK( __LINE__, VELZ_C(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S31_X(k,i,j) = 0.5_RP * ( VELZ_C(k,i+1,j) - VELZ_C(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S31_Y(k,i,j) = 0.5_RP * ( WORK_Z(k,i,j) - WORK_Z(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dy
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, WORK_X(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S23_Z(k,i,j) = 0.5_RP * ( WORK_X(k,i,j) - WORK_X(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S23_X(k,i,j) = 0.5_RP * ( WORK_Z(k,i,j) - WORK_Z(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S23_Y(k,i,j) = 0.5_RP * ( WORK_Z(k,i,j) - WORK_Z(k,i-1,j) ) *RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! u
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
#endif
          WORK_Z(k,i,j) = 0.5_RP * ( VELX_YZ(k,i,j+1) + VELX_YZ(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i,j) )
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k+1,i,j+1) )
#endif
          WORK_X(k,i,j) = 0.25_RP * ( VELX_C(k,i,j) + VELX_C(k+1,i,j) + VELX_C(k,i,j+1) + VELX_C(k+1,i,j+1) )
       enddo
       enddo
       enddo
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
#endif
          WORK_Y(k,i,j) = 0.5_RP * ( VELX_YZ(k+1,i,j) + VELX_YZ(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! du/dx
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, WORK_Y(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_Z(k,i,j) = ( WORK_Y(k,i,j) - WORK_Y(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i+1,j) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S11_X(k,i,j) = ( VELX_C(k,i+1,j) - VELX_C(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_Y(k,i,j) = ( WORK_Z(k,i,j) - WORK_Z(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dy
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, WORK_X(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S12_Z(k,i,j) = 0.5_RP * ( WORK_X(k,i,j) - WORK_X(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k,i,j-1) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S12_X(k,i,j) = 0.5_RP * ( VELX_YZ(k,i,j+1) - VELX_YZ(k,i,j-1) ) / ( FDY(j) + FDY(j-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S12_Y(k,i,j) = 0.5_RP * ( VELX_C(k,i,j+1) - VELX_C(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dz
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_Z(k,i,j) )
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S31_Z(k,i,j) = S31_Z(k,i,j) + &
               0.5_RP * ( VELX_C(k+1,i,j) - VELX_C(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, S31_X(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S31_X(k,i,j) = S31_X(k,i,j) + &
               0.5_RP * ( VELX_YZ(k+1,i,j) - VELX_YZ(k-1,i,j) ) / ( FDZ(k) + FDZ(k-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S31_X(KS,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
#endif
          S31_X(KS,i,j) = S31_X(KS,i,j) + &
               0.5_RP * ( VELX_YZ(KS+1,i,j) - VELX_YZ(KS,i,j) ) * RFDZ(KS)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_Y(k,i,j) )
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, WORK_X(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S31_Y(k,i,j) = S31_Y(k,i,j) + &
               0.5_RP * ( WORK_X(k,i,j) - WORK_X(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
       do j = JJS-1, JJE
       do i = IIS  , IIE
#ifdef DEBUG
       call CHECK( __LINE__, S31_Y(KS,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i+1,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i+1,j+1) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S31_Y(KS,i,j) = S31_Y(KS,i,j) + &
               0.125_RP * ( VELX_YZ(KS+1,i,j) + VELX_YZ(KS+1,i+1,j) + VELX_YZ(KS+1,i,j+1) + VELX_YZ(KS+1,i+1,j+1) &
                          - VELX_YZ(KS  ,i,j) + VELX_YZ(KS  ,i+1,j) + VELX_YZ(KS  ,i,j+1) + VELX_YZ(KS  ,i+1,j+1) ) * RCDZ(KS)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS  , IIE
#ifdef DEBUG
       call CHECK( __LINE__, S31_Y(KE,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE,i+1,j) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(KE,i+1,j+1) )
       call CHECK( __LINE__, RCDZ(KE-1) )
#endif
          S31_Y(KE,i,j) = S31_Y(KE,i,j) + &
               0.125_RP * ( VELX_YZ(KE  ,i,j) + VELX_YZ(KE  ,i+1,j) + VELX_YZ(KE  ,i,j+1) + VELX_YZ(KE  ,i+1,j+1) &
                          - VELX_YZ(KE-1,i,j) + VELX_YZ(KE-1,i+1,j) + VELX_YZ(KE-1,i,j+1) + VELX_YZ(KE-1,i+1,j+1) ) * RCDZ(KE-1)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! v
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
#endif
          WORK_Z(k,i,j) = 0.5_RP * ( VELY_ZX(k,i+1,j) + VELY_ZX(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
#endif
          WORK_X(k,i,j) = 0.5_RP * ( VELY_ZX(k+1,i,j) + VELY_ZX(k,i,j) )
       enddo
       enddo
       enddo
       ! (y edge)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_C(k,i,j) )
       call CHECK( __LINE__, VELY_C(k+1,i,j) )
       call CHECK( __LINE__, VELY_C(k,i+1,j) )
       call CHECK( __LINE__, VELY_C(k+1,i+1,j) )
#endif
          WORK_Y(k,i,j) = 0.25_RP * ( VELY_C(k,i,j) + VELY_C(k+1,i,j) + VELY_C(k,i+1,j) + VELY_C(k+1,i+1,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! dv/dy
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, WORK_X(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S22_Z(k,i,j) = ( WORK_X(k,i,j) - WORK_X(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k,i,j-1) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S22_X(k,i,j) = ( WORK_Z(k,i,j) - WORK_Z(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELY_C(k,i,j+1) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S22_Y(k,i,j) = ( VELY_C(k,i,j+1) - VELY_C(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dv/dz
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_Z(k,i,j) )
       call CHECK( __LINE__, VELY_C(k+1,i,j) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S23_Z(k,i,j) = S23_Z(k,i,j) + &
               0.5_RP * ( VELY_C(k+1,i,j) - VELY_C(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_X(k,i,j) )
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, WORK_Y(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S23_X(k,i,j) = S23_X(k,i,j) + &
               0.5_RP * ( WORK_Y(k,i,j) - WORK_Y(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S23_X(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i,j+1) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i+1,j+1) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j+1) )
       call CHECK( __LINE__, VELY_ZX(KS,i+1,j+1) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S23_X(KS,i,j) = S23_X(KS,i,j) + &
               0.125_RP * ( VELY_ZX(KS+1,i,j) + VELY_ZX(KS+1,i+1,j) + VELY_ZX(KS+1,i,j+1) + VELY_ZX(KS+1,i+1,j+1) &
                           -VELY_ZX(KS  ,i,j) - VELY_ZX(KS  ,i+1,j) - VELY_ZX(KS  ,i,j+1) - VELY_ZX(KS  ,i+1,j+1) ) * RCDZ(KS)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S23_X(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j+1) )
       call CHECK( __LINE__, VELY_ZX(KE,i+1,j+1) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i,j+1) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i+1,j+1) )
       call CHECK( __LINE__, RCDZ(KE-1) )
#endif
          S23_X(KE,i,j) = S23_X(KE,i,j) + &
               0.125_RP * ( VELY_ZX(KE  ,i,j) + VELY_ZX(KE  ,i+1,j) + VELY_ZX(KE  ,i,j+1) + VELY_ZX(KE  ,i+1,j+1) &
                           -VELY_ZX(KE-1,i,j) - VELY_ZX(KE-1,i+1,j) - VELY_ZX(KE-1,i,j+1) - VELY_ZX(KE-1,i+1,j+1) ) * RCDZ(KE-1)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, S23_Y(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S23_Y(k,i,j) = S23_Y(k,i,j) + &
               0.5_RP * ( VELY_ZX(k+1,i,j) - VELY_ZX(k-1,i,j) ) / ( FDZ(k) + FDZ(k-1) )
       enddo
       enddo
       enddo
       do j = JJS-1, JJE
       do i = IIS  , IIE
#ifdef DEBUG
       call CHECK( __LINE__, S23_Y(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
#endif
          S23_Y(KS,i,j) = S23_Y(KS,i,j) + &
               0.5_RP * ( VELY_ZX(KS+1,i,j) - VELY_ZX(KS,i,j) ) * RFDZ(KS)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dv/dx
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, WORK_Y(k,i-1,j) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S12_Z(k,i,j) = S12_Z(k,i,j) + &
               0.5_RP * ( WORK_Y(k,i,j) - WORK_Y(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S12_X(k,i,j) )
       call CHECK( __LINE__, VELY_C(k,i+1,j) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S12_X(k,i,j) = S12_X(k,i,j) + &
               0.5_RP * ( VELY_C(k,i+1,j) - VELY_C(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S12_Y(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k,i-1,j) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S12_Y(k,i,j) = S12_Y(k,i,j) + &
               0.5_RP * ( VELY_ZX(k,i+1,j) - VELY_ZX(k,i-1,j) ) / ( FDX(i) + FDX(i-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


#ifdef DEBUG
       S2(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
       nu_Z(:,:,:) = UNDEF; nu_X(:,:,:) = UNDEF; nu_Y(:,:,:) = UNDEF
#endif
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S11_Z(k,i,j) )
       call CHECK( __LINE__, S22_Z(k,i,j) )
       call CHECK( __LINE__, S33_Z(k,i,j) )
       call CHECK( __LINE__, S31_Z(k,i,j) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
       call CHECK( __LINE__, S23_Z(k,i,j) )
#endif
          S2(k,i,j) = &
                 2.0_RP * ( S11_Z(k,i,j)**2 + S22_Z(k,i,j)**2 + S33_Z(k,i,j)**2 ) &
               + 4.0_RP * ( S31_Z(k,i,j)**2 + S12_Z(k,i,j)**2 + S23_Z(k,i,j)**2 )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Ri
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          WORK_Z(k,i,j) = 0.5_RP * GRAV * ( POTT(k+1,i,j) - POTT(k,i,j) ) &
               * RFDZ(k) / ( ( POTT(k+1,i,j) + POTT(k,i,j) ) * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, nu_factXY(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( WORK_Z(k,i,j) < 0.0_RP ) then
             nu_Z(k,i,j) = nu_factXY(k,i,j) * RPrN &
                  * sqrt( S2(k,i,j) * (1.0_RP - FhB*WORK_Z(k,i,j)) )
          else if ( WORK_Z(k,i,j) < RiC ) then
             nu_Z(k,i,j) = nu_factXY(k,i,j) * RPrN &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - WORK_Z(k,i,j)*RRiC )**4 &
                  * ( 1.0_RP - PrNovRiC*WORK_Z(k,i,j) )
          else
             nu_Z(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
#ifdef DEBUG
       S2(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (y-z plane)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S11_X(k,i,j) )
       call CHECK( __LINE__, S22_X(k,i,j) )
       call CHECK( __LINE__, S33_X(k,i,j) )
       call CHECK( __LINE__, S31_X(k,i,j) )
       call CHECK( __LINE__, S12_X(k,i,j) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          S2(k,i,j) = &
                 2.0_RP * ( S11_X(k,i,j)**2 + S22_X(k,i,j)**2 + S33_X(k,i,j)**2 ) &
               + 4.0_RP * ( S31_X(k,i,j)**2 + S12_X(k,i,j)**2 + S23_X(k,i,j)**2 )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Ri
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i+1,j) )
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i+1,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k-1,i+1,j) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          TMP1 = ( POTT(k+1,i+1,j) + POTT(k+1,i,j) ) * 0.5_RP
          TMP2 = ( POTT(k  ,i+1,j) + POTT(k  ,i,j) ) * 0.5_RP
          TMP3 = ( POTT(k-1,i+1,j) + POTT(k-1,i,j) ) * 0.5_RP
          WORK_X(k,i,j) = 0.5_RP * GRAV * ( TMP1 - TMP3 ) &
               / ( ( FDZ(k) + FDZ(k-1) ) * TMP2 * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KE,i+1,j) )
       call CHECK( __LINE__, POTT(KE,i,j) )
       call CHECK( __LINE__, POTT(KE-1,i+1,j) )
       call CHECK( __LINE__, POTT(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
       call CHECK( __LINE__, S2(KE,i,j) )
#endif
          TMP2 = ( POTT(KE  ,i+1,j) + POTT(KE  ,i,j) ) * 0.5_RP
          TMP3 = ( POTT(KE-1,i+1,j) + POTT(KE-1,i,j) ) * 0.5_RP
          WORK_X(KE,i,j) = 0.5_RP * GRAV * ( TMP2 - TMP3 ) &
               * RFDZ(KE-1) / ( TMP2 * max(S2(KE,i,j),1.0E-20) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KS+1,i+1,j) )
       call CHECK( __LINE__, POTT(KS+1,i,j) )
       call CHECK( __LINE__, POTT(KS,i+1,j) )
       call CHECK( __LINE__, POTT(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S2(KS,i,j) )
#endif
          TMP1 = ( POTT(KS+1,i+1,j) + POTT(KS+1,i,j) ) * 0.5_RP
          TMP2 = ( POTT(KS  ,i+1,j) + POTT(KS  ,i,j) ) * 0.5_RP
          WORK_X(KS,i,j) = 0.5_RP * GRAV * ( TMP1 - TMP2 ) &
               * RFDZ(KS) / ( TMP2 * max(S2(KS,i,j),1.0E-20) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_X(k,i,j) )
       call CHECK( __LINE__, nu_factYZ(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( WORK_X(k,i,j) < 0.0_RP ) then
             nu_X(k,i,j) = nu_factYZ(k,i,j) * RPrN &
                  * sqrt( S2(k,i,j) * (1.0_RP - FhB*WORK_X(k,i,j)) )
          else if ( WORK_X(k,i,j) < RiC ) then
             nu_X(k,i,j) = nu_factYZ(k,i,j) * RPrN &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - WORK_X(k,i,j)*RRiC )**4 &
                  * ( 1.0_RP - PrNovRiC*WORK_X(k,i,j) )
          else
             nu_X(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

#ifdef DEBUG
       S2(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, S11_Y(k,i,j) )
       call CHECK( __LINE__, S22_Y(k,i,j) )
       call CHECK( __LINE__, S33_Y(k,i,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
       call CHECK( __LINE__, S12_Y(k,i,j) )
       call CHECK( __LINE__, S23_Y(k,i,j) )
#endif
          S2(k,i,j) = &
                 2.0_RP * ( S11_Y(k,i,j)**2 + S22_Y(k,i,j)**2 + S33_Y(k,i,j)**2 ) &
               + 4.0_RP * ( S31_Y(k,i,j)**2 + S12_Y(k,i,j)**2 + S23_Y(k,i,j)**2 )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Ri
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j+1) )
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j+1) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k-1,i,j+1) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          TMP1 = ( POTT(k+1,i,j+1) + POTT(k+1,i,j) ) * 0.5_RP
          TMP2 = ( POTT(k  ,i,j+1) + POTT(k  ,i,j) ) * 0.5_RP
          TMP3 = ( POTT(k-1,i,j+1) + POTT(k-1,i,j) ) * 0.5_RP
          WORK_Y(k,i,j) = 0.5_RP * GRAV * ( TMP1 - TMP3 ) &
               / ( ( FDZ(k) + FDZ(k-1) ) * TMP2 * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS  , IIE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KE,i,j+1) )
       call CHECK( __LINE__, POTT(KE,i,j) )
       call CHECK( __LINE__, POTT(KE-1,i,j+1) )
       call CHECK( __LINE__, POTT(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
       call CHECK( __LINE__, S2(KE,i,j) )
#endif
          TMP2 = ( POTT(KE  ,i,j+1) + POTT(KE  ,i,j) ) * 0.5_RP
          TMP3 = ( POTT(KE-1,i,j+1) + POTT(KE-1,i,j) ) * 0.5_RP
          WORK_Y(KE,i,j) = 0.5_RP * GRAV * ( TMP2 - TMP3 ) &
               * RFDZ(KE-1) / ( TMP2 * max(S2(KE,i,j),1.0E-20) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS  , IIE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KS+1,i,j+1) )
       call CHECK( __LINE__, POTT(KS+1,i,j) )
       call CHECK( __LINE__, POTT(KS,i,j+1) )
       call CHECK( __LINE__, POTT(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S2(KS,i,j) )
#endif
          TMP1 = ( POTT(KS+1,i,j+1) + POTT(KS+1,i,j) ) * 0.5_RP
          TMP2 = ( POTT(KS  ,i,j+1) + POTT(KS  ,i,j) ) * 0.5_RP
          WORK_Y(KS,i,j) = 0.5_RP * GRAV * ( TMP1 - TMP2 ) &
               * RFDZ(KS) / ( TMP2 * max(S2(KS,i,j),1.0E-20) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, WORK_Y(k,i,j) )
       call CHECK( __LINE__, nu_factZX(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( WORK_Y(k,i,j) < 0.0_RP ) then
             nu_Y(k,i,j) = nu_factZX(k,i,j) * RPrN &
                  * sqrt( S2(k,i,j) * (1.0_RP - FhB*WORK_Y(k,i,j)) )
          else if ( WORK_Y(k,i,j) < RiC ) then
             nu_Y(k,i,j) = nu_factZX(k,i,j) * RPrN &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - WORK_Y(k,i,j)*RRiC )**4 &
                  * ( 1.0_RP - PrNovRiC*WORK_Y(k,i,j) )
          else
             nu_Y(k,i,j) = 0.0_RP
          end if
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif


#ifdef DEBUG
       qflx_sgs(:,:,:,:) = UNDEF
#endif
       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, nu_Z(k,i,j) )
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          qflx_sgs(k,i,j,ZDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                               * nu_Z(k,i,j) &
                               * ( POTT(k+1,i,j)-POTT(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
       call CHECK( __LINE__, SFLX_POTT(i,j) )
#endif
          qflx_sgs(KS-1,i,j,ZDIR) = SFLX_POTT(i,j)
          qflx_sgs(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (y-z plane)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, nu_X(k,i,j) )
       call CHECK( __LINE__, POTT(k,i+1,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          qflx_sgs(k,i,j,XDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                               * nu_X(k,i,j) &
                               * ( POTT(k,i+1,j)-POTT(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, nu_Y(k,i,j) )
       call CHECK( __LINE__, POTT(k,i,j+1) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          qflx_sgs(k,i,j,YDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
                               * nu_Y(k,i,j) &
                               * ( POTT(k,i,j+1)-POTT(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !--- tendency rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, qflx_sgs(k,i,j,ZDIR) )
       call CHECK( __LINE__, qflx_sgs(k-1,i,j,ZDIR) )
       call CHECK( __LINE__, RCDZ(k) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,XDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i-1,j,XDIR) )
       call CHECK( __LINE__, RCDX(i) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,YDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i,j-1,YDIR) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          RHOT_t(k,i,j) = - ( ( qflx_sgs(k,i,j,ZDIR)-qflx_sgs(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                            + ( qflx_sgs(k,i,j,XDIR)-qflx_sgs(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                            + ( qflx_sgs(k,i,j,YDIR)-qflx_sgs(k  ,i,  j-1,YDIR) ) * RCDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    enddo
    enddo
#ifdef DEBUG
       IIS = IUNDEF; IIE = IUNDEF; JJS = IUNDEF; JJE = IUNDEF
#endif

    !##### Tracers #####
    do iq = 1, QA

#ifdef DEBUG
       qflx_sgs(:,:,:,:) = UNDEF
#endif
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, nu_Z(k,i,j) )
       call CHECK( __LINE__, QTRC(k+1,i,j,iq) )
       call CHECK( __LINE__, QTRC(k,i,j,iq) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          qflx_sgs(k,i,j,ZDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                               * nu_Z(k,i,j) &
                               * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) * RFDZ(k)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_sgs(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! Surface QV Flux
       if ( iq == I_QV ) then
         do j = JJS, JJE
         do i = IIS, IIE
#ifdef DEBUG
       call CHECK( __LINE__, SFLX_QV(i,j) )
#endif
             qflx_sgs(KS-1,i,j,ZDIR) = SFLX_QV(i,j)
          enddo
          enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       endif

       ! (y-z plane)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, nu_X(k,i,j) )
       call CHECK( __LINE__, QTRC(k,i+1,j,iq) )
       call CHECK( __LINE__, QTRC(k,i,j,iq) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          qflx_sgs(k,i,j,XDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                               * nu_X(k,i,j) &
                               * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) * RFDX(i)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, nu_Y(k,i,j) )
       call CHECK( __LINE__, QTRC(k,i,j+1,iq) )
       call CHECK( __LINE__, QTRC(k,i,j,iq) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          qflx_sgs(k,i,j,YDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
                               * nu_Y(k,i,j) &
                               * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) * RFDY(j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !--- tendency tracers
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, qflx_sgs(k,i,j,ZDIR) )
       call CHECK( __LINE__, qflx_sgs(k-1,i,j,ZDIR) )
       call CHECK( __LINE__, RCDZ(k) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,XDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i-1,j,XDIR) )
       call CHECK( __LINE__, RCDX(i) )
       call CHECK( __LINE__, qflx_sgs(k,i,j,YDIR) )
       call CHECK( __LINE__, qflx_sgs(k,i,j-1,YDIR) )
       call CHECK( __LINE__, RCDY(j) )
       call CHECK( __LINE__, DENS(k,i,j) )
#endif
          QTRC_t(k,i,j,iq) = - ( ( qflx_sgs(k,i,j,ZDIR)-qflx_sgs(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                               + ( qflx_sgs(k,i,j,XDIR)-qflx_sgs(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                               + ( qflx_sgs(k,i,j,YDIR)-qflx_sgs(k  ,i,  j-1,YDIR) ) * RCDY(j) ) &
                             / DENS(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    enddo
    enddo
#ifdef DEBUG
       IIS = IUNDEF; IIE = IUNDEF; JJS = IUNDEF; JJE = IUNDEF
#endif

    enddo ! scalar quantities loop
#ifdef DEBUG
       iq = IUNDEF
#endif

    return
  end subroutine ATMOS_PHY_TB_main

end module mod_atmos_phy_tb
