!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!          Smagolinsky-type
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-29 (S.Iga)       [new]
!! @li      2011-12-11 (H.Yashiro)   [mod] integrate to SCALE-LES ver.3
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-03-27 (H.Yashiro)   [mod] reconstruction
!! @li      2012-07-02 (S.Nishizawa) [mod] reconstruction with Brown et al. (1994)
!! @li      2012-10-26 (S.Nishizawa) [mod] remove surface flux
!! @li      2013-06-13 (S.Nishizawa) [mod] change mixing length by Brown et al. (1994) and Scotti et al. (1993)
!! @li      2014-04-02 (S.Nishizawa) [mod] modified for terrrain-following
!!
!! - Reference
!!  - Brown et al., 1994:
!!    Large-eddy simulaition of stable atmospheric boundary layers with a revised stochastic subgrid model.
!!    Roy. Meteor. Soc., 120, 1485-1512
!!  - Scotti et al., 1993:
!!    Generalized Smagorinsky model for anisotropic grids.
!!    Phys. Fluids A, 5, 2306-2308
!!  - Sullivan et al., 1994:
!!    A subgrid-scale model for large-eddy simulation of planetary boundary-layer flows.
!!    Bound.-Layer Meteor., 71, 247-276
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_tb_smg
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
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
  public :: ATMOS_PHY_TB_smg_setup
  public :: ATMOS_PHY_TB_smg

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
  real(RP), private, allocatable :: nu_fact (:,:,:) ! (Cs*Delta)^2

  real(RP), private            :: Cs  = 0.13_RP ! Smagorinsky constant (Scotti et al. 1993)
  real(RP), private, parameter :: Ck  = 0.1_RP  ! SGS constant (Moeng and Wyngaard 1988)
  real(RP), private, parameter :: PrN = 0.7_RP  ! Prandtl number in neutral conditions
  real(RP), private, parameter :: RiC = 0.25_RP ! critical Richardson number
  real(RP), private, parameter :: FmC = 16.0_RP ! fum = sqrt(1 - c*Ri)
  real(RP), private, parameter :: FhB = 40.0_RP ! fuh = sqrt(1 - b*Ri)/PrN
  real(RP), private            :: RPrN          ! 1 / PrN
  real(RP), private            :: RRiC          ! 1 / RiC
  real(RP), private            :: PrNovRiC      ! PrN / RiC

  real(RP), private, parameter :: OneOverThree = 1.0_RP / 3.0_RP
  real(RP), private, parameter :: twoOverThree = 2.0_RP / 3.0_RP

  logical, private  :: ATMOS_PHY_TB_SMG_bottom = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_smg_setup( &
       TYPE_TB,       &
       CDZ, CDX, CDY, &
       CZ             )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: TYPE_TB

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA)

    real(RP) :: ATMOS_PHY_TB_SMG_Cs
    real(RP) :: ATMOS_PHY_TB_SMG_filter_fact = 2.0_RP

    NAMELIST / PARAM_ATMOS_PHY_TB_SMG / &
         ATMOS_PHY_TB_SMG_Cs, &
         ATMOS_PHY_TB_SMG_filter_fact, &
         ATMOS_PHY_TB_SMG_bottom

    integer :: k, i, j

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[TURBULENCE] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Smagorinsky-type Eddy Viscocity Model'

    if ( TYPE_TB /= 'SMAGORINSKY' ) then
       write(*,*) 'xxx ATMOS_PHY_TB_TYPE is not SMAGORINSKY. Check!'
       call PRC_MPIstop
    endif

    ATMOS_PHY_TB_SMG_Cs = Cs

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_SMG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_TB_SMG. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_TB_SMG)

    Cs = ATMOS_PHY_TB_SMG_Cs

    RPrN     = 1.0_RP / PrN
    RRiC     = 1.0_RP / RiC
    PrNovRiC = (1- PrN) * RRiC

    allocate( nu_fact(KA,IA,JA) )

#ifdef DEBUG
    nu_fact (:,:,:) = UNDEF
#endif
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
       nu_fact(k,i,j) = ( Cs * mixlen(CDZ(k),CDX(i),CDY(j),CZ(k),ATMOS_PHY_TB_SMG_filter_fact) )**2
    enddo
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    return
  end subroutine ATMOS_PHY_TB_smg_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_smg( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, &
!       qflx_sgs_rhot, qflx_sgs_rhoq,                &
       qflx_sgs_rhot, qflx_sgs_qtrc,                &
       tke, nu, Ri, Pr,                             &
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          &
       GSQRT, J13G, J23G, J33G                      )
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_grid, only: &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ
    implicit none

    ! SGS flux
    real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
!    real(RP), intent(out) :: qflx_sgs_rhoq(KA,IA,JA,QA,3)
    real(RP), intent(out) :: qflx_sgs_qtrc(KA,IA,JA,QA,3)

    real(RP), intent(out) :: tke(KA,IA,JA) ! TKE
    real(RP), intent(out) :: nu (KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Pr (KA,IA,JA) ! Prantle number
    real(RP), intent(out) :: Ri (KA,IA,JA) ! Richardson number

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix

    ! diagnostic variables
    real(RP) :: VELZ_C (KA,IA,JA)
    real(RP) :: VELZ_XY(KA,IA,JA)
    real(RP) :: VELX_C (KA,IA,JA)
    real(RP) :: VELX_YZ(KA,IA,JA)
    real(RP) :: VELY_C (KA,IA,JA)
    real(RP) :: VELY_ZX(KA,IA,JA)
    real(RP) :: POTT   (KA,IA,JA)

    ! deformation rate tensor
    real(RP) :: S33_C (KA,IA,JA) ! (cell center)
    real(RP) :: S11_C (KA,IA,JA)
    real(RP) :: S22_C (KA,IA,JA)
    real(RP) :: S31_C (KA,IA,JA)
    real(RP) :: S12_C (KA,IA,JA)
    real(RP) :: S23_C (KA,IA,JA)
    real(RP) :: S12_Z (KA,IA,JA) ! (z edge or x-y plane)
    real(RP) :: S23_X (KA,IA,JA) ! (x edge or y-z plane)
    real(RP) :: S31_Y (KA,IA,JA) ! (y edge or z-x plane)

    real(RP) :: S2    (KA,IA,JA) ! |S|^2
    real(RP) :: WORK_V(KA,IA,JA) ! work space (vertex)
    real(RP) :: WORK_Z(KA,IA,JA) !            (z edge or x-y plane)
    real(RP) :: WORK_X(KA,IA,JA) !            (x edge or y-z plane)
    real(RP) :: WORK_Y(KA,IA,JA) !            (y edge or z-x plane)

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Turbulence(smagorinsky)'

#ifdef DEBUG
    qflx_sgs_momz(:,:,:,:)   = UNDEF
    qflx_sgs_momx(:,:,:,:)   = UNDEF
    qflx_sgs_momy(:,:,:,:)   = UNDEF
    qflx_sgs_rhot(:,:,:,:)   = UNDEF
!    qflx_sgs_rhoq(:,:,:,:,:) = UNDEF
    qflx_sgs_qtrc(:,:,:,:,:) = UNDEF

    nu (:,:,:) = UNDEF
    tke(:,:,:) = UNDEF
    Pr (:,:,:) = UNDEF
    Ri (:,:,:) = UNDEF

    VELZ_C (:,:,:) = UNDEF
    VELZ_XY(:,:,:) = UNDEF
    VELX_C (:,:,:) = UNDEF
    VELX_YZ(:,:,:) = UNDEF
    VELY_C (:,:,:) = UNDEF
    VELY_ZX(:,:,:) = UNDEF
    POTT   (:,:,:) = UNDEF

    S33_C (:,:,:) = UNDEF
    S11_C (:,:,:) = UNDEF
    S22_C (:,:,:) = UNDEF
    S31_C (:,:,:) = UNDEF
    S12_C (:,:,:) = UNDEF
    S23_C (:,:,:) = UNDEF
    S12_Z (:,:,:) = UNDEF
    S23_X (:,:,:) = UNDEF
    S31_Y (:,:,:) = UNDEF

    S2    (:,:,:) = UNDEF
    WORK_V(:,:,:) = UNDEF
    WORK_Z(:,:,:) = UNDEF
    WORK_X(:,:,:) = UNDEF
    WORK_Y(:,:,:) = UNDEF
#endif

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
    do j = JS-2, JE+2
    do i = IS-2, IE+2
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
    do j = JS-2, JE+2
    do i = IS-2, IE+2
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
    do i = IS-2, IE+1
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
    do i = IS-2, IE+1
       VELX_YZ(KE+1,i,j) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-2, JE+2
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

    do j = JS-2, JE+1
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
    do j = JS-2, JE+1
    do i = IS-1, IE+1
       VELY_ZX(KE+1,i,j) = 0.0_RP
    enddo
    enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    do j = JS-1, JE+2
    do i = IS-2, IE+2
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
       ! (x-y plane; x,y,w)
       ! WORK_Z = VELZ_XY
       ! (y-z plane; u,y,z)
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
       ! (z-x plane; x,v,z)
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

       ! dw/dz
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k-1,i,j) )
       call CHECK( __LINE__, RCDZ(k) )
#endif
          S33_C(k,i,j) = ( VELZ_XY(k,i,j) - VELZ_XY(k-1,i,j) ) * RCDZ(k) &
                       * J33G / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(KS,i,j) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDZ(KS) )
#endif
          S33_C(KS,i,j) = VELZ_XY(KS,i,j) * RCDZ(KS) & ! VELZ_XY(KS-1,i,j) == 0
                        * J33G / GSQRT(KS,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dx
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i+1,j) )
       call CHECK( __LINE__, VELZ_C(k,i-1,j) )
       call CHECK( __LINE__, GSQRT(k,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k-1,i,j) )
       call CHECK( __LINE__, J13G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J13G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S31_C(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i+1,j,I_XYZ)*VELZ_C(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ)*VELZ_C(k,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
               + ( J13G(k,i,j,I_XYW)*VELZ_XY(k,i,j) - J13G(k-1,i,j,I_XYW)*VELZ_XY(k-1,i,j) ) * RCDZ(k) &
               )

       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(KS,i+1,j) )
       call CHECK( __LINE__, VELZ_C(KS,i-1,j) )
       call CHECK( __LINE__, GSQRT(KS,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KS,i,j) )
       call CHECK( __LINE__, J13G(KS,i,j,I_XYW) )
       call CHECK( __LINE__, VELZ_C(KE,i+1,j) )
       call CHECK( __LINE__, VELZ_C(KE,i-1,j) )
       call CHECK( __LINE__, GSQRT(KE,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KE,i,j) )
       call CHECK( __LINE__, J13G(KE,i,j,I_XYW) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S31_C(KS,i,j) = 0.5_RP * ( &
                 ( GSQRT(KS,i+1,j,I_XYZ)*VELZ_C(KS,i+1,j) - GSQRT(KS,i-1,j,I_XYZ)*VELZ_C(KS,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
               + ( J13G(KS,i,j,I_XYW)*VELZ_XY(KS,i,j) ) * RCDZ(KS) &
               )
          S31_C(KE,i,j) = 0.5_RP * ( &
                 ( GSQRT(KE,i+1,j,I_XYZ)*VELZ_C(KE,i+1,j) - GSQRT(KE,i-1,j,I_XYZ)*VELZ_C(KE,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
               - ( J13G(KE-1,i,j,I_XYW)*VELZ_XY(KE-1,i,j) ) * RCDZ(KE) &
               )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (y edge, u,y,w)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i+1,j) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S31_Y(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i+1,j,I_XYW)*VELZ_XY(k,i+1,j) - GSQRT(k,i,j,I_XYW)*VELZ_XY(k,i,j) ) * RFDX(i) &
               + ( J13G(k+1,i,j,I_UYZ)*WORK_X(k+1,i,j) - J13G(k,i,j,I_UYZ)*WORK_X (k,i,j)) * RFDZ(k) &
               )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dw/dy
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(k,i,j+1) )
       call CHECK( __LINE__, VELZ_C(k,i,j-1) )
       call CHECK( __LINE__, GSQRT(k,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, VELZ_XY(k-1,i,j) )
       call CHECK( __LINE__, J23G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J23G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S23_C(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_XYZ)*VELZ_C(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ)*VELZ_C(k,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(k,i,j,I_XYW)*VELZ_XY(k,i,j) - J23G(k-1,i,j,I_XYW)*VELZ_XY(k-1,i,j) ) * RCDZ(k) &
               )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_C(KS,i,j+1) )
       call CHECK( __LINE__, VELZ_C(KS,i,j-1) )
       call CHECK( __LINE__, GSQRT(KS,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KS,i,j) )
       call CHECK( __LINE__, J23G(KS,i,j,I_XYW) )
       call CHECK( __LINE__, VELZ_C(KE,i,j+1) )
       call CHECK( __LINE__, VELZ_C(KE,i,j-1) )
       call CHECK( __LINE__, GSQRT(KE,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELZ_XY(KE,i,j) )
       call CHECK( __LINE__, J23G(KE,i,j,I_XYW) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S23_C(KS,i,j) = 0.5_RP * ( &
                 ( GSQRT(KS,i,j+1,I_XYZ)*VELZ_C(KS,i,j+1) - GSQRT(KS,i,j-1,I_XYZ)*VELZ_C(KS,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(KS,i,j,I_XYW)*VELZ_XY(KS,i,j) ) * RCDZ(KS) &
               )
          S23_C(KE,i,j) = 0.5_RP * ( &
                 ( GSQRT(KE,i,j+1,I_XYZ)*VELZ_C(KE,i,j+1) - GSQRT(KE,i,j-1,I_XYZ)*VELZ_C(KE,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               - ( J23G(KE-1,i,j,I_XYW)*VELZ_XY(KE-1,i,j) ) * RCDZ(KE) &
               )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (x edge; x,v,w)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELZ_XY(k,i,j+1) )
       call CHECK( __LINE__, VELZ_XY(k,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S23_X(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_XYW)*VELZ_XY(k,i,j+1) - GSQRT(k,i,j,I_XYW)*VELZ_XY(k,i,j) ) * RFDY(j) &
               + ( J23G(k+1,i,j,I_XVZ)*WORK_Y(k+1,i,j) - J23G(k,i,j,I_XVZ)*WORK_Y (k,i,j) ) * RFDZ(k) &
               )
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
       ! (x-y plane; x,y,w)
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
       ! (y-z plane; u,y,z)
       ! WORK_X = VELX_YZ
       ! (z-x plane; x,v,z)
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
       ! (vertex; u,v,w)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j+1) )
       call CHECK( __LINE__, J23G(k  ,i,j  ,I_UVZ) )
       call CHECK( __LINE__, J23G(k+1,i,j  ,I_UVZ) )
       call CHECK( __LINE__, J23G(k  ,i,j+1,I_UVZ) )
       call CHECK( __LINE__, J23G(k+1,i,j+1,I_UVZ) )
#endif
          WORK_V(k,i,j) = 0.25_RP &
               * ( J23G(k  ,i,j  ,I_UYZ)*VELX_YZ(k  ,i,j  ) &
                 + J23G(k+1,i,j  ,I_UYZ)*VELX_YZ(k+1,i,j  ) &
                 + J23G(k  ,i,j+1,I_UYZ)*VELX_YZ(k  ,i,j+1) &
                 + J23G(k+1,i,j+1,I_UYZ)*VELX_YZ(k+1,i,j+1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! du/dx
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i-1,j) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_UYZ) )
       call CHECK( __LINE__, GSQRT(k,i-1,j,I_UYZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J13G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J13G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_C(k,i,j) = ( &
                 ( GSQRT(k,i,j,I_UYZ)*VELX_YZ(k,i,j) - GSQRT(k,i-1,j,I_UYZ)*VELX_YZ(k,i-1,j) ) * RCDX(i) &
               + ( J13G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J13G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(KS,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS,i-1,j) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_UYZ) )
       call CHECK( __LINE__, GSQRT(KS,i-1,j,I_UYZ) )
       call CHECK( __LINE__, VELX_C(KS+1,i,j) )
       call CHECK( __LINE__, VELX_C(KS,i,j) )
       call CHECK( __LINE__, J13G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE,i-1,j) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_UYZ) )
       call CHECK( __LINE__, GSQRT(KE,i-1,j,I_UYZ) )
       call CHECK( __LINE__, VELX_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE-1,i,j) )
       call CHECK( __LINE__, J13G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KE-1,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDX(i) )
#endif
          S11_C(KS,i,j) = ( &
                 ( GSQRT(KS,i,j,I_UYZ)*VELX_YZ(KS,i,j) - GSQRT(KS,i-1,j,I_UYZ)*VELX_YZ(KS,i-1,j) ) * RCDX(i) &
               + ( J13G(KS+1,i,j,I_XYZ)*VELX_C(KS+1,i,j) - J13G(KS,i,j,I_XYZ)*VELX_C(KS,i,j) ) * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S11_C(KE,i,j) = ( &
                 ( GSQRT(KE,i,j,I_UYZ)*VELX_YZ(KE,i,j) - GSQRT(KE,i-1,j,I_UYZ)*VELX_YZ(KE,i-1,j) ) * RCDX(i) &
               + ( J13G(KE,i,j,I_XYZ)*VELX_C(KE,i,j) - J13G(KE-1,i,j,I_XYZ)*VELX_C(KE-1,i,j) ) * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dz
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_C(k,i,j) )
       call CHECK( __LINE__, VELX_C(k+1,i,j) )
       call CHECK( __LINE__, VELX_C(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S31_C(k,i,j) = ( S31_C(k,i,j) & ! dw/dx
               + 0.5_RP * ( VELX_C(k+1,i,j) - VELX_C(k-1,i,j) ) * J33G / ( FDZ(k) + FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S31_C(KS,i,j) )
       call CHECK( __LINE__, VELX_C(KS+1,i,j) )
       call CHECK( __LINE__, VELX_C(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S31_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
#endif
          S31_C(KS,i,j) = ( S31_C(KS,i,j) &
               + 0.5_RP * ( VELX_C(KS+1,i,j) - VELX_C(KS,i,j) ) * J33G * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S31_C(KE,i,j) = ( S31_C(KE,i,j) &
               + 0.5_RP * ( VELX_C(KE,i,j) - VELX_C(KE-1,i,j) ) * J33G * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (y edge; u,y,w)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S31_Y(k,i,j) )
       call CHECK( __LINE__, VELX_YZ(k+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S31_Y(k,i,j) = ( S31_Y(k,i,j) & ! dw/dx
               + 0.5_RP * ( VELX_YZ(k+1,i,j) - VELX_YZ(k,i,j) ) * J33G * RFDZ(k) &
               ) / GSQRT(k,i,j,I_UYW)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * du/dy
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(k,i,j+1) )
       call CHECK( __LINE__, VELX_C(k,i,j-1) )
       call CHECK( __LINE__, GSQRT(k,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i,j-1,I_XYZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J23G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J23G(k-1,i,j,I_XYW) )
#endif
          S12_C(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_XYZ)*VELX_C(k,i,j+1) - GSQRT(k,i,j-1,I_XYZ)*VELX_C(k,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J23G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_C(KS,i,j+1) )
       call CHECK( __LINE__, VELX_C(KS,i,j-1) )
       call CHECK( __LINE__, GSQRT(KS,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELX_C(KS+1,i,j) )
       call CHECK( __LINE__, VELX_C(KS,i,j) )
       call CHECK( __LINE__, J23G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, VELX_C(KE,i,j+1) )
       call CHECK( __LINE__, VELX_C(KE,i,j-1) )
       call CHECK( __LINE__, GSQRT(KE,i,j+1,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j-1,I_XYZ) )
       call CHECK( __LINE__, VELX_C(KE,i,j) )
       call CHECK( __LINE__, VELX_C(KE-1,i,j) )
       call CHECK( __LINE__, J23G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KE-1,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, FDY(j) )
       call CHECK( __LINE__, FDY(j-1) )
#endif
          S12_C(KS,i,j) = 0.5_RP * ( &
                 ( GSQRT(KS,i,j+1,I_XYZ)*VELX_C(KS,i,j+1) - GSQRT(KS,i,j-1,I_XYZ)*VELX_C(KS,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(KS+1,i,j,I_XYZ)*VELX_C(KS+1,i,j) - J23G(KS,i,j,I_XYZ)*VELX_C(KS,i,j) ) * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S12_C(KE,i,j) = 0.5_RP * ( &
                 ( GSQRT(KE,i,j+1,I_XYZ)*VELX_C(KE,i,j+1) - GSQRT(KE,i,j-1,I_XYZ)*VELX_C(KE,i,j-1) ) / ( FDY(j) + FDY(j-1) ) &
               + ( J23G(KE,i,j,I_XYZ)*VELX_C(KE,i,j) - J23G(KE-1,i,j,I_XYZ)*VELX_C(KE-1,i,j) ) * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (z edge; u,v,z)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(k,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k-1,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          S12_Z(k,i,j) = 0.5_RP * ( &
                 ( GSQRT(k,i,j+1,I_UYZ)*VELX_YZ(k,i,j+1) - GSQRT(k,i,j,I_UYZ)*VELX_YZ(k,i,j) ) * RFDY(j) &
               + ( WORK_V(k,i,j) - WORK_V(k-1,i,j) ) * RCDZ(k) &
               )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, VELX_YZ(KS,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(KS,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i,j) )
       call CHECK( __LINE__, VELX_YZ(KS+1,i,j+1) )
       call CHECK( __LINE__, J23G(KS+1,i,j,I_UVZ) )
       call CHECK( __LINE__, J23G(KS  ,i,j,I_UVZ) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j+1) )
       call CHECK( __LINE__, VELX_YZ(KE,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE-1,i,j) )
       call CHECK( __LINE__, VELX_YZ(KE-1,i,j+1) )
       call CHECK( __LINE__, J23G(KE  ,i,j,I_UVZ) )
       call CHECK( __LINE__, J23G(KE-1,i,j,I_UVZ) )
#endif
          S12_Z(KS,i,j) = 0.25_RP * ( &
                 ( GSQRT(KS,i,j+1,I_UYZ)*VELX_YZ(KS,i,j+1) - GSQRT(KS,i,j,I_UYZ)*VELX_YZ(KS,i,j) ) * RFDY(j) &
               + ( J23G(KS+1,i,j,I_UVZ) * ( VELX_YZ(KS+1,i,j) + VELX_YZ(KS+1,i,j+1) ) &
                 - J23G(KS  ,i,j,I_UVZ) * ( VELX_YZ(KS  ,i,j) + VELX_YZ(KS  ,i,j+1) ) ) * RFDZ(KS) &
               )
          S12_Z(KE,i,j) = 0.25_RP * ( &
                 ( GSQRT(KE,i,j+1,I_UYZ)*VELX_YZ(KE,i,j+1) - GSQRT(KE,i,j,I_UYZ)*VELX_YZ(KE,i,j) ) * RFDY(j) &
               + ( J23G(KE  ,i,j,I_UVZ) * ( VELX_YZ(KE  ,i,j) + VELX_YZ(KE  ,i,j+1) ) &
                 - J23G(KE-1,i,j,I_UVZ) * ( VELX_YZ(KE-1,i,j) + VELX_YZ(KE-1,i,j+1) ) ) * RFDZ(KE-1) &
               )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

#ifdef DEBUG
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF; WORK_V(:,:,:) = UNDEF
#endif
       ! v
       ! (x-y plane; x,y,w)
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
       ! (y-z plane; u,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, VELY_C(k,i+1,j) )
       call CHECK( __LINE__, VELY_C(k,i,j) )
#endif
          WORK_X(k,i,j) = 0.5_RP * ( VELZ_C(k,i+1,j) + VELZ_C(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane; x,v,z)
       ! WORK_Y = VELY_ZX
       ! (vertex; u,v,w)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i+1,j) )
#endif
          WORK_V(k,i,j) = 0.25_RP &
               * ( J13G(k  ,i  ,j,I_XVZ)*VELY_ZX(k  ,i  ,j) &
                 + J13G(k+1,i  ,j,I_XVZ)*VELY_ZX(k+1,i  ,j) &
                 + J13G(k  ,i+1,j,I_XVZ)*VELY_ZX(k  ,i+1,j) &
                 + J13G(k+1,i+1,j,I_XVZ)*VELY_ZX(k+1,i+1,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! dv/dy
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j-1) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_XVZ) )
       call CHECK( __LINE__, GSQRT(k,i,j-1,I_XVZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J23G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J23G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, RCDY(j) )
#endif
          S22_C(k,i,j) = ( &
                 ( GSQRT(k,i,j,I_XVZ)*VELY_ZX(k,i,j) - GSQRT(k,i,j-1,I_XVZ)*VELY_ZX(k,i,j-1) ) * RCDY(j) &
               + ( J23G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J23G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, VELY_ZX(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j-1) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XVZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j-1,I_XVZ) )
       call CHECK( __LINE__, VELY_C(KS+1,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i,j) )
       call CHECK( __LINE__, J23G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, RCDY(j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j-1) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XVZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j-1,I_XVZ) )
       call CHECK( __LINE__, VELY_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE-1,i,j) )
       call CHECK( __LINE__, J23G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J23G(KE-1,i,j,I_XYZ) )
#endif
          S22_C(KS,i,j) = ( &
                 ( GSQRT(KS,i,j,I_XVZ)*VELY_ZX(KS,i,j) - GSQRT(KS,i,j-1,I_XVZ)*VELY_ZX(KS,i,j-1) ) * RCDY(j) &
               + ( J23G(KS+1,i,j,I_XYZ)*VELY_C(KS+1,i,j) - J23G(KS,i,j,I_XYZ)*VELY_C(KS,i,j) ) * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S22_C(KE,i,j) = ( &
                 ( GSQRT(KE,i,j,I_XVZ)*VELY_ZX(KE,i,j) - GSQRT(KE,i,j-1,I_XVZ)*VELY_ZX(KE,i,j-1) ) * RCDY(j) &
               + ( J23G(KE,i,j,I_XYZ)*VELY_C(KE,i,j) - J23G(KE-1,i,j,I_XYZ)*VELY_C(KE-1,i,j) ) * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dv/dx
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S12_C(k,i,j) )
       call CHECK( __LINE__, VELY_C(k,i+1,j) )
       call CHECK( __LINE__, VELY_C(k,i-1,j) )
       call CHECK( __LINE__, GSQRT(k,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(k,i-1,j,I_XYZ) )
       call CHECK( __LINE__, WORK_Z(k,i,j) )
       call CHECK( __LINE__, WORK_Z(k-1,i,j) )
       call CHECK( __LINE__, J13G(k,i,j,I_XYW) )
       call CHECK( __LINE__, J13G(k-1,i,j,I_XYW) )
       call CHECK( __LINE__, GSQRT(k,i,j,I_XYZ) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S12_C(k,i,j) = ( S12_C(k,i,j) & ! du/dy
               + 0.5_RP * ( &
                     ( GSQRT(k,i+1,j,I_XYZ)*VELY_C(k,i+1,j) - GSQRT(k,i-1,j,I_XYZ)*VELY_C(k,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                   + ( J13G(k,i,j,I_XYW)*WORK_Z(k,i,j) - J13G(k-1,i,j,I_XYW)*WORK_Z(k-1,i,j) ) * RCDZ(k) ) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S12_C(KS,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i+1,j) )
       call CHECK( __LINE__, VELY_C(KS,i-1,j) )
       call CHECK( __LINE__, GSQRT(KS,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELY_C(KS+1,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i,j) )
       call CHECK( __LINE__, J13G(KS+1,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KS,i,j,I_XYZ) )
       call CHECK( __LINE__, S12_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE,i+1,j) )
       call CHECK( __LINE__, VELY_C(KE,i-1,j) )
       call CHECK( __LINE__, GSQRT(KE,i+1,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i-1,j,I_XYZ) )
       call CHECK( __LINE__, VELY_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE-1,i,j) )
       call CHECK( __LINE__, J13G(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, J13G(KE-1,i,j,I_XYZ) )
       call CHECK( __LINE__, GSQRT(KE,i,j,I_XYZ) )
       call CHECK( __LINE__, FDX(i) )
       call CHECK( __LINE__, FDX(i-1) )
#endif
          S12_C(KS,i,j) = ( S12_C(KS,i,j) & ! du/dy
               + 0.5_RP * ( &
                     ( GSQRT(KS,i+1,j,I_XYZ)*VELY_C(KS,i+1,j) - GSQRT(KS,i-1,j,I_XYZ)*VELY_C(KS,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                   + ( J13G(KS+1,i,j,I_XYZ)*VELY_C(KS+1,i,j) - J13G(KS,i,j,I_XYZ)*VELY_C(KS,i,j) ) * RFDZ(KS) ) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S12_C(KE,i,j) = ( S12_C(KE,i,j) & ! du/dy
               + 0.5_RP * ( &
                     ( GSQRT(KE,i+1,j,I_XYZ)*VELY_C(KE,i+1,j) - GSQRT(KE,i-1,j,I_XYZ)*VELY_C(KE,i-1,j) ) / ( FDX(i) + FDX(i-1) ) &
                   + ( J13G(KE,i,j,I_XYZ)*VELY_C(KE,i,j) - J13G(KE-1,i,j,I_XYZ)*VELY_C(KE-1,i,j) ) * RFDZ(KE-1) ) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z edge; u,v,z)
       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S12_Z(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, WORK_V(k,i,j) )
       call CHECK( __LINE__, WORK_V(k-1,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S12_Z(k,i,j) = ( S12_Z(k,i,j) &
               + 0.5_RP * ( &
                     ( GSQRT(k,i+1,j,I_XVZ)*VELY_ZX(k,i+1,j) - GSQRT(k,i,j,I_XVZ)*VELY_ZX(k,i,j) ) * RFDX(i) &
                   + ( WORK_V(k,i,j) - WORK_V(k-1,i,j) ) * RCDZ(k) ) &
               ) / GSQRT(k,i,j,I_UVZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, S12_Z(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KS,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KS+1,i+1,j) )
       call CHECK( __LINE__, S12_Z(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i+1,j) )
       call CHECK( __LINE__, VELY_ZX(KE,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i,j) )
       call CHECK( __LINE__, VELY_ZX(KE-1,i+1,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          S12_Z(KS,i,j) = ( S12_Z(KS,i,j) &
               + 0.5_RP * ( &
                     ( GSQRT(KS,i+1,j,I_XVZ)*VELY_ZX(KS,i+1,j) - GSQRT(KS,i,j,I_XVZ)*VELY_ZX(KS,i,j) ) * RFDX(i) &
                   + ( J13G(KS+1,i,j,I_UVZ) * ( VELY_ZX(KS+1,i,j) + VELY_ZX(KS+1,i+1,j) ) &
                     - J13G(KS  ,i,j,I_UVZ) * ( VELY_ZX(KS  ,i,j) + VELY_ZX(KS  ,i+1,j) ) ) * RFDZ(KS) ) &
               ) / GSQRT(KS,i,j,I_UVZ)
          S12_Z(KE,i,j) = ( S12_Z(KE,i,j) &
               + 0.5_RP * ( &
                     ( GSQRT(KE,i+1,j,I_XVZ)*VELY_ZX(KE,i+1,j) - GSQRT(KE,i,j,I_XVZ)*VELY_ZX(KE,i,j) ) * RFDX(i) &
                   + ( J13G(KE  ,i,j,I_UVZ) * ( VELY_ZX(KE  ,i,j) + VELY_ZX(KE  ,i+1,j) ) &
                     - J13G(KE-1,i,j,I_UVZ) * ( VELY_ZX(KE-1,i,j) + VELY_ZX(KE-1,i+1,j) ) ) * RFDZ(KE-1) ) &
               ) / GSQRT(KE,i,j,I_UVZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! 1/2 * dv/dz
       ! (cell center; x,y,z)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_C(k,i,j) )
       call CHECK( __LINE__, VELY_C(k+1,i,j) )
       call CHECK( __LINE__, VELY_C(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
#endif
          S23_C(k,i,j) = ( S23_C(k,i,j) & ! dw/dy
               + 0.5_RP * ( VELY_C(k+1,i,j) - VELY_C(k-1,i,j) ) * J33G / ( FDZ(k) + FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_XYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, S23_C(KS,i,j) )
       call CHECK( __LINE__, VELY_C(KS+1,i,j) )
       call CHECK( __LINE__, VELY_C(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S23_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE,i,j) )
       call CHECK( __LINE__, VELY_C(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
#endif
          S23_C(KS,i,j) = ( S23_C(KS,i,j) &
               + 0.5_RP * ( VELY_C(KS+1,i,j) - VELY_C(KS,i,j) ) * J33G * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XYZ)
          S23_C(KE,i,j) = ( S23_C(KE,i,j) &
               + 0.5_RP * ( VELY_C(KE,i,j) - VELY_C(KE-1,i,j) ) * J33G * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (x edge; x,v,w)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, S23_X(k,i,j) )
       call CHECK( __LINE__, VELY_ZX(k+1,i,j) )
       call CHECK( __LINE__, VELY_ZX(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          S23_X(k,i,j) = ( S23_X(k,i,j) &
               + 0.5_RP * ( VELY_ZX(k+1,i,j) - VELY_ZX(k,i,j) ) * J33G * RFDZ(k) &
               ) / GSQRT(k,i,j,I_XVW)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif


       ! nu_SGS = (Cs * Delta)^2 * |S|, |S|^2 = 2*Sij*Sij
#ifdef DEBUG
       S2(:,:,:) = UNDEF
       WORK_Z(:,:,:) = UNDEF; WORK_X(:,:,:) = UNDEF; WORK_Y(:,:,:) = UNDEF
#endif
       ! (cell center)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
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
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          Ri(k,i,j) = GRAV * ( POTT(k+1,i,j) - POTT(k-1,i,j) ) * J33G &
               / ( ( FDZ(k) + FDZ(k-1) ) * GSQRT(k,i,j,I_XYZ) * POTT(k,i,j) * max(S2(k,i,j),1.0E-20_RP) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KS+1,i,j) )
       call CHECK( __LINE__, POTT(KS,i,j) )
       call CHECK( __LINE__, RFDZ(KS) )
       call CHECK( __LINE__, S2(KS,i,j) )
#endif
          Ri(KS,i,j) = GRAV * ( POTT(KS+1,i,j) - POTT(KS,i,j) ) * J33G &
               * RFDZ(KS) / (GSQRT(KS,i,j,I_XYZ) * POTT(KS,i,j) * max(S2(KS,i,j),1.0E-20_RP) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
#ifdef DEBUG
       call CHECK( __LINE__, POTT(KE,i,j) )
       call CHECK( __LINE__, POTT(KE-1,i,j) )
       call CHECK( __LINE__, RFDZ(KE-1) )
       call CHECK( __LINE__, S2(KE,i,j) )
#endif
          Ri(KE,i,j) = GRAV * ( POTT(KE,i,j) - POTT(KE-1,i,j) ) * J33G &
               * RFDZ(KE-1) / (GSQRT(KE,i,j,I_XYZ) * POTT(KE,i,j) * max(S2(KE,i,j),1.0E-20_RP) )
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, Ri(k,i,j) )
       call CHECK( __LINE__, nu_fact(k,i,j) )
       call CHECK( __LINE__, S2(k,i,j) )
#endif
          if ( Ri(k,i,j) < 0.0_RP ) then ! stable
             nu(k,i,j) = nu_fact(k,i,j) &
                  * sqrt( S2(k,i,j) * (1.0_RP - FmC*Ri(k,i,j)) )
          else if ( Ri(k,i,j) < RiC ) then ! weakly stable
             nu(k,i,j) = nu_fact(k,i,j) &
                  * sqrt( S2(k,i,j) ) * ( 1.0_RP - Ri(k,i,j)*RRiC )**4
          else ! strongly stable
             nu(k,i,j) = 0.0_RP
          endif
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! tke = (nu/(Ck * Delta))^2 = ( nu * Cs / Ck )^2 / ( Cs * Delta )^2
       ! Sullivan et al. (1994)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu_fact(k,i,j) )
#endif
          tke(k,i,j) = ( nu(k,i,j) * Cs / Ck )**2 / nu_fact(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! Pr = nu_m / nu_h = fm / fh
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
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
             Pr(k,i,j) = 1.0_RP
          endif
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

       !##### momentum equation (z) #####
       ! (cell center)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, tke(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,ZDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu(k,i,j) &
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
          qflx_sgs_momz(KS,i,j,ZDIR) = 0.0_RP ! just above bottom boundary
          qflx_sgs_momz(KE,i,j,ZDIR) = 0.0_RP ! just below top boundary
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, nu(k+1,i+1,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,XDIR) = - 0.125_RP & ! 2.0 / 4 / 4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
               * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i+1,j)+nu  (k+1,i+1,j)) &
               * S31_Y(k,i,j)
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, nu(k+1,i,j+1) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          qflx_sgs_momz(k,i,j,YDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
               * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i,j+1)+nu  (k+1,i,j+1) ) &
               * S23_X(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, nu(k+1,i+1,j) )
       call CHECK( __LINE__, S31_Y(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,ZDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
               * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i+1,j)+nu  (k+1,i+1,j) ) &
               * S31_Y(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momx(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_momx(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, tke(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,XDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu(k,i,j) &
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, nu(k,i+1,j+1) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,YDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ( nu  (k,i,j)+nu  (k,i+1,j)+nu  (k,i,j+1)+nu  (k,i+1,j+1) ) &
               * S12_Z(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, nu(k+1,i,j+1) )
       call CHECK( __LINE__, S23_X(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,ZDIR) = - 0.125_RP & ! 2/4/4
               * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
               * ( nu  (k,i,j)+nu  (k+1,i,j)+nu  (k,i,j+1)+nu  (k+1,i,j+1) ) &
               * S23_X(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momy(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_momy(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, nu(k,i+1,j+1) )
       call CHECK( __LINE__, S12_Z(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,XDIR) = - 0.125_RP & !
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ( nu  (k,i,j)+nu  (k,i+1,j)+nu  (k,i,j+1)+nu  (k,i+1,j+1) ) &
               * S12_Z(k,i,j)
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
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, S11_C(k,i,j) )
       call CHECK( __LINE__, S22_C(k,i,j) )
       call CHECK( __LINE__, S33_C(k,i,j) )
       call CHECK( __LINE__, tke(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,YDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu(k,i,j) &
               * ( S22_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) &
             + twoOverThree * tke(k,i,j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !##### Thermodynamic Equation #####

       ! (x-y plane; x,y,w)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
          qflx_sgs_rhot(k,i,j,ZDIR) = - 0.25_RP & ! 2/2/2/2
               * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
               * ( nu(k,i,j)/Pr(k,i,j) + nu(k+1,i,j)/Pr(k+1,i,j) ) &
               * ( POTT(k+1,i,j)-POTT(k,i,j) ) * RFDZ(k) * J33G / GSQRT(k,i,j,I_XYW)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_rhot(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_sgs_rhot(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (y-z plane; u,y,z)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, Pr(k,i,j) )
       call CHECK( __LINE__, Pr(k,i+1,j) )
       call CHECK( __LINE__, POTT(k,i+1,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k+1,i+1,j) )
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k-1,i+1,j) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          qflx_sgs_rhot(k,i,j,XDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(k,i,j) + DENS(k,i+1,j) ) &
               * ( nu(k,i,j)/Pr(k,i,j) + nu(k,i+1,j)/Pr(k,i+1,j) ) &
               * ( &
                     ( GSQRT(k,i+1,j,I_XYZ) * POTT(k,i+1,j) &
                     - GSQRT(k,i  ,j,I_XYZ) * POTT(k,i  ,j) ) * RFDX(i) &
                   + ( J13G(k+1,i,j,I_UYZ) * ( POTT(k+1,i+1,j)+POTT(k+1,i,j) ) &
                     - J13G(k-1,i,j,I_UYZ) * ( POTT(k-1,i+1,j)+POTT(k-1,i,j) ) &
                     ) * 0.5_RP / ( FDZ(k) + FDZ(k-1) ) &
                 ) / GSQRT(k,i,j,I_UYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS,   JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i+1,j) )
       call CHECK( __LINE__, nu(KS,i,j) )
       call CHECK( __LINE__, nu(KS,i+1,j) )
       call CHECK( __LINE__, POTT(KS,i+1,j) )
       call CHECK( __LINE__, POTT(KS,i,j) )
       call CHECK( __LINE__, RFDX(i) )
#endif
          qflx_sgs_rhot(KS,i,j,XDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(KS,i,j)+DENS(KS,i+1,j) ) &
               * ( nu(KS,i,j)/Pr(KS,i,j) + nu(KS,i+1,j)/Pr(KS,i+1,j) ) &
               * ( &
                     ( GSQRT(KS,i+1,j,I_XYZ) * POTT(KS,i+1,j) &
                     - GSQRT(KS,i  ,j,I_XYZ) * POTT(KS,i  ,j) ) * RFDX(i) &
                   + ( J13G(KS+1,i,j,I_UYZ) * ( POTT(KS+1,i+1,j)+POTT(KS+1,i,j) ) &
                     - J13G(KS  ,i,j,I_UYZ) * ( POTT(KS  ,i+1,j)+POTT(KS  ,i,j) ) &
                     ) * 0.5_RP * RFDZ(KS) &
                 ) / GSQRT(KS,i,j,I_UYZ)
          qflx_sgs_rhot(KE,i,j,XDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(KE,i,j)+DENS(KE,i+1,j) ) &
               * ( nu(KE,i,j)/Pr(KE,i,j) + nu(KE,i+1,j)/Pr(KE,i+1,j) ) &
               * ( &
                     ( GSQRT(KE,i+1,j,I_XYZ) * POTT(KE,i+1,j) &
                     - GSQRT(KE,i  ,j,I_XYZ) * POTT(KE,i  ,j) ) * RFDX(i) &
                   + ( J13G(KE  ,i,j,I_UYZ) * ( POTT(KE  ,i+1,j)+POTT(KE  ,i,j) ) &
                     - J13G(KE-1,i,j,I_UYZ) * ( POTT(KE-1,i+1,j)+POTT(KE-1,i,j) ) &
                     ) * 0.5_RP * RFDZ(KE-1) &
                 ) / GSQRT(KE,i,j,I_UYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane; x,v,z)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS+1, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, POTT(k,i,j+1) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k+1,i,j+1) )
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k-1,i,j+1) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          qflx_sgs_rhot(k,i,j,YDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
               * ( nu(k,i,j)/Pr(k,i,j) + nu(k,i,j+1)/Pr(k,i,j+1) ) &
               * ( &
                   ( GSQRT(k,i,j+1,I_XYZ) * POTT(k,i,j+1) &
                   - GSQRT(k,i,j  ,I_XYZ) * POTT(k,i,j  ) ) * RFDY(j) &
                 + ( J23G(k+1,i,j,I_XVZ) * ( POTT(k+1,i,j+1)+POTT(k+1,i,j) ) &
                   - J23G(k-1,i,j,I_XVZ) * ( POTT(k-1,i,j+1)+POTT(k-1,i,j) ) &
                   ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_XVZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS,   IIE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i,j+1) )
       call CHECK( __LINE__, nu(KS,i,j) )
       call CHECK( __LINE__, nu(KS,i,j+1) )
       call CHECK( __LINE__, POTT(KS,i,j+1) )
       call CHECK( __LINE__, POTT(KS,i,j) )
       call CHECK( __LINE__, RFDY(j) )
#endif
          qflx_sgs_rhot(KS,i,j,YDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(KS,i,j)+DENS(KS,i,j+1) ) &
               * ( nu(KS,i,j)/Pr(KS,i,j) + nu(KS,i,j+1)/Pr(KS,i,j+1) ) &
               * ( &
                   ( GSQRT(KS,i,j+1,I_XYZ) * POTT(KS,i,j+1) &
                   - GSQRT(KS,i,j  ,I_XYZ) * POTT(KS,i,j  ) ) * RFDY(j) &
                 + ( J23G(KS+1,i,j,I_XVZ) * ( POTT(KS+1,i,j+1)+POTT(KS+1,i,j) ) &
                   - J23G(KS  ,i,j,I_XVZ) * ( POTT(KS  ,i,j+1)+POTT(KS  ,i,j) ) &
                   ) * 0.5_RP * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XVZ)
          qflx_sgs_rhot(KE,i,j,YDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(KE,i,j)+DENS(KE,i,j+1) ) &
               * ( nu(KE,i,j)/Pr(KE,i,j) + nu(KE,i,j+1)/Pr(KE,i,j+1) ) &
               * ( &
                   ( GSQRT(KE,i,j+1,I_XYZ) * POTT(KE,i,j+1) &
                   - GSQRT(KE,i,j  ,I_XYZ) * POTT(KE,i,j  ) ) * RFDY(j) &
                 + ( J23G(KE  ,i,j,I_XVZ) * ( POTT(KE  ,i,j+1)+POTT(KE  ,i,j) ) &
                   - J23G(KE-1,i,j,I_XVZ) * ( POTT(KE-1,i,j+1)+POTT(KE-1,i,j) ) &
                   ) * 0.5_RP * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XVZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

    enddo
    enddo


    !##### Tracers #####
    do iq = 1, QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! (x-y plane; x,y,w)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k+1,i,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k+1,i,j) )
       call CHECK( __LINE__, QTRC(k+1,i,j,iq) )
       call CHECK( __LINE__, QTRC(k,i,j,iq) )
       call CHECK( __LINE__, RFDZ(k) )
#endif
!          qflx_sgs_rhoq(k,i,j,iq,ZDIR) = - 0.25_RP & ! 1/2/2
!               * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
          qflx_sgs_qtrc(k,i,j,iq,ZDIR) = - 0.5_RP & ! 1/2
               * ( nu(k,i,j)/Pr(k,i,j) + nu(k+1,i,j)/Pr(k+1,i,j) ) &
               * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) * RFDZ(k) * J33G / GSQRT(k,i,j,I_XYW)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS, JJE
       do i = IIS, IIE
!          qflx_sgs_rhoq(KS-1,i,j,iq,ZDIR) = 0.0_RP
!          qflx_sgs_rhoq(KE  ,i,j,iq,ZDIR) = 0.0_RP
          qflx_sgs_qtrc(KS-1,i,j,iq,ZDIR) = 0.0_RP
          qflx_sgs_qtrc(KE  ,i,j,iq,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       ! (y-z plane; u,y,z)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS+1,  KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i+1,j) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i+1,j) )
       call CHECK( __LINE__, QTRC(k,i+1,j,iq) )
       call CHECK( __LINE__, QTRC(k,i,j,iq) )
       call CHECK( __LINE__, RFDX(i) )
#endif
!          qflx_sgs_rhoq(k,i,j,iq,XDIR) = - 0.25_RP & ! 1/2/2
!               * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
          qflx_sgs_qtrc(k,i,j,iq,XDIR) = - 0.5_RP & ! 1/2
               * ( nu(k,i,j)/Pr(k,i,j) + nu(k,i+1,j)/Pr(k,i+1,j) ) &
               * ( &
                   ( GSQRT(k,i+1,j,I_XYZ) * QTRC(k,i+1,j,iq) &
                   - GSQRT(k,i  ,j,I_XYZ) * QTRC(k,i  ,j,iq) ) * RFDX(i) &
                 + ( J13G(k+1,i,j,I_UYZ) * ( QTRC(k+1,i+1,j,iq)+QTRC(k+1,i,j,iq) ) &
                   - J13G(k-1,i,j,I_UYZ) * ( QTRC(k-1,i+1,j,iq)+QTRC(k-1,i,j,iq) ) &
                   ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_UYZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS,   JJE
       do i = IIS-1, IIE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i+1,j) )
       call CHECK( __LINE__, nu(KS,i,j) )
       call CHECK( __LINE__, nu(KS,i+1,j) )
       call CHECK( __LINE__, QTRC(KS,i+1,j,iq) )
       call CHECK( __LINE__, QTRC(KS,i,j,iq) )
       call CHECK( __LINE__, RFDX(i) )
#endif
!          qflx_sgs_rhoq(KS,i,j,iq,XDIR) = - 0.25_RP & ! 1/2/2
!               * ( DENS(KS,i,j)+DENS(KS,i+1,j) ) &
          qflx_sgs_qtrc(KS,i,j,iq,XDIR) = - 0.5_RP & ! 1/2
               * ( nu(KS,i,j)/Pr(KS,i,j) + nu(KS,i+1,j)/Pr(KS,i+1,j) ) &
               * ( &
                   ( GSQRT(KS,i+1,j,I_XYZ) * QTRC(KS,i+1,j,iq) &
                   - GSQRT(KS,i  ,j,I_XYZ) * QTRC(KS,i  ,j,iq) ) * RFDX(i) &
                 + ( J13G(KS+1,i,j,I_UYZ) * ( QTRC(KS+1,i+1,j,iq)+QTRC(KS+1,i,j,iq) ) &
                   - J13G(KS  ,i,j,I_UYZ) * ( QTRC(KS  ,i+1,j,iq)+QTRC(KS  ,i,j,iq) ) &
                   ) * 0.5_RP * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_UYZ)
!          qflx_sgs_rhoq(KE,i,j,iq,XDIR) = - 0.25_RP & ! 1/2/2
!               * ( DENS(KE,i,j)+DENS(KE,i+1,j) ) &
          qflx_sgs_qtrc(KE,i,j,iq,XDIR) = - 0.5_RP & ! 1/2
               * ( nu(KE,i,j)/Pr(KE,i,j) + nu(KE,i+1,j)/Pr(KE,i+1,j) ) &
               * ( &
                   ( GSQRT(KE,i+1,j,I_XYZ) * QTRC(KE,i+1,j,iq) &
                   - GSQRT(KE,i  ,j,I_XYZ) * QTRC(KE,i  ,j,iq) ) * RFDX(i) &
                 + ( J13G(KE  ,i,j,I_UYZ) * ( QTRC(KE  ,i+1,j,iq)+QTRC(KE  ,i,j,iq) ) &
                   - J13G(KE-1,i,j,I_UYZ) * ( QTRC(KE-1,i+1,j,iq)+QTRC(KE-1,i,j,iq) ) &
                   ) * 0.5_RP * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_UYZ)
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ! (z-x plane; x,v,z)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS+1,  KE-1
#ifdef DEBUG
       call CHECK( __LINE__, DENS(k,i,j) )
       call CHECK( __LINE__, DENS(k,i,j+1) )
       call CHECK( __LINE__, nu(k,i,j) )
       call CHECK( __LINE__, nu(k,i,j+1) )
       call CHECK( __LINE__, QTRC(k,i,j+1,iq) )
       call CHECK( __LINE__, QTRC(k,i,j,iq) )
       call CHECK( __LINE__, RFDY(j) )
#endif
!          qflx_sgs_rhoq(k,i,j,iq,YDIR) = - 0.25_RP &
!               * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
          qflx_sgs_qtrc(k,i,j,iq,YDIR) = - 0.5_RP &
               * ( nu(k,i,j)/Pr(k,i,j) + nu(k,i,j+1)/Pr(k,i,j+1) ) &
               * ( &
                     ( GSQRT(k,i,j+1,I_XYZ) * QTRC(k,i,j+1,iq) &
                     - GSQRT(k,i,j  ,I_XYZ) * QTRC(k,i,j  ,iq) ) * RFDY(j) &
                   + ( J23G(k+1,i,j,I_XVZ) * ( QTRC(k+1,i,j+1,iq)+QTRC(k+1,i,j,iq) ) &
                     - J23G(k-1,i,j,I_XVZ) * ( QTRC(k-1,i,j+1,iq)+QTRC(k-1,i,j,iq) ) &
                     ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
               ) / GSQRT(k,i,j,I_XVZ)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       do j = JJS-1, JJE
       do i = IIS,   IIE
#ifdef DEBUG
       call CHECK( __LINE__, DENS(KS,i,j) )
       call CHECK( __LINE__, DENS(KS,i,j+1) )
       call CHECK( __LINE__, nu(KS,i,j) )
       call CHECK( __LINE__, nu(KS,i,j+1) )
       call CHECK( __LINE__, QTRC(KS,i,j+1,iq) )
       call CHECK( __LINE__, QTRC(KS,i,j,iq) )
       call CHECK( __LINE__, RFDY(j) )
#endif
!          qflx_sgs_rhoq(KS,i,j,iq,YDIR) = - 0.25_RP &
!               * ( DENS(KS,i,j)+DENS(KS,i,j+1) ) &
          qflx_sgs_qtrc(KS,i,j,iq,YDIR) = - 0.5_RP &
               * ( nu(KS,i,j)/Pr(KS,i,j) + nu(KS,i,j+1)/Pr(KS,i,j+1) ) &
               * ( &
                     ( GSQRT(KS,i,j+1,I_XYZ) * QTRC(KS,i,j+1,iq) &
                     - GSQRT(KS,i,j  ,I_XYZ) * QTRC(KS,i,j  ,iq) ) * RFDY(j) &
                   + ( J23G(KS+1,i,j,I_XVZ) * ( QTRC(KS+1,i,j+1,iq)+QTRC(KS+1,i,j,iq) ) &
                     - J23G(KS  ,i,j,I_XVZ) * ( QTRC(KS  ,i,j+1,iq)+QTRC(KS  ,i,j,iq) ) &
                     ) * 0.5_RP * RFDZ(KS) &
               ) / GSQRT(KS,i,j,I_XVZ)
!          qflx_sgs_rhoq(KE,i,j,iq,YDIR) = - 0.25_RP &
!               * ( DENS(KE,i,j)+DENS(KE,i,j+1) ) &
          qflx_sgs_qtrc(KE,i,j,iq,YDIR) = - 0.5_RP &
               * ( nu(KE,i,j)/Pr(KE,i,j) + nu(KE,i,j+1)/Pr(KE,i,j+1) ) &
               * ( &
                     ( GSQRT(KE,i,j+1,I_XYZ) * QTRC(KE,i,j+1,iq) &
                     - GSQRT(KE,i,j  ,I_XYZ) * QTRC(KE,i,j  ,iq) ) * RFDY(j) &
                   + ( J23G(KE  ,i,j,I_XVZ) * ( QTRC(KE  ,i,j+1,iq)+QTRC(KE  ,i,j,iq) ) &
                     - J23G(KE-1,i,j,I_XVZ) * ( QTRC(KE-1,i,j+1,iq)+QTRC(KE-1,i,j,iq) ) &
                     ) * 0.5_RP * RFDZ(KE-1) &
               ) / GSQRT(KE,i,j,I_XVZ)
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
  end subroutine ATMOS_PHY_TB_smg


  function mixlen(dz, dx, dy, z, filter_fact)
  use scale_const, only: &
     KARMAN  => CONST_KARMAN
    implicit none
    real(RP), intent(in) :: dz
    real(RP), intent(in) :: dx
    real(RP), intent(in) :: dy
    real(RP), intent(in) :: z
    real(RP), intent(in) :: filter_fact
    real(RP) :: mixlen ! (out)

    real(RP) :: d0

    d0 = fact(dz, dx, dy) * filter_fact * ( dz * dx * dy )**OneOverThree ! Scotti et al. (1993)
    if ( ATMOS_PHY_TB_SMG_bottom ) then
       mixlen = sqrt( 1.0_RP / ( 1.0_RP/d0**2 + 1.0_RP/(KARMAN*z)**2 ) ) ! Brown et al. (1994)
    else
       mixlen = d0
    end if

    return
  end function mixlen

  function fact(dz, dx, dy)
    real(RP), intent(in) :: dz
    real(RP), intent(in) :: dx
    real(RP), intent(in) :: dy
    real(RP) :: fact ! (out)

    real(RP), parameter :: oot = -1.0_RP/3.0_RP
    real(RP), parameter :: fot =  5.0_RP/3.0_RP
    real(RP), parameter :: eot = 11.0_RP/3.0_RP
    real(RP), parameter :: tof = -3.0_RP/4.0_RP
    real(RP) :: a1, a2, b1, b2, dmax


    dmax = max(dz, dx, dy)
    if ( dz .eq. dmax ) then
       a1 = dx / dmax
       a2 = dy / dmax
    else if ( dx .eq. dmax ) then
       a1 = dz / dmax
       a2 = dy / dmax
    else ! dy .eq. dmax
       a1 = dz / dmax
       a2 = dx / dmax
    end if
    b1 = atan( a1/a2 )
    b2 = atan( a2/a1 )

   fact = 1.736_RP * (a1*a2)**oot &
         * ( 4.0_RP*p1(b1)*a1**oot + 0.222_RP*p2(b1)*a1**fot + 0.077*p3(b1)*a1**eot - 3.0_RP*b1 &
           + 4.0_RP*p1(b2)*a2**oot + 0.222_RP*p2(b2)*a2**fot + 0.077*p3(b2)*a2**eot - 3.0_RP*b2 &
           )**tof
   return
  end function fact
  function p1(z)
    real(RP), intent(in) :: z
    real(RP) :: p1 ! (out)

    p1 = 2.5_RP * p2(z) - 1.5_RP * sin(z) * cos(z)**TwoOverThree
    return
  end function p1
  function p2(z)
    real(RP), intent(in) :: z
    real(RP) :: p2 ! (out)

    p2 = 0.986_RP * z + 0.073_RP * z**2 - 0.418_RP * z**3 + 0.120_RP * z**4
    return
  end function p2
  function p3(z)
    real(RP), intent(in) :: z
    real(RP) :: p3 ! (out)

    p3 = 0.976_RP * z + 0.188_RP * z**2 - 1.169_RP * z**3 + 0.755_RP * z**4 - 0.151_RP * z**5
    return
  end function p3

end module scale_atmos_phy_tb_smg
