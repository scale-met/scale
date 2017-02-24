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
  public :: ATMOS_PHY_TB_smg_config
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
  real(RP), private, parameter   :: OneOverThree  = 1.0_RP / 3.0_RP
  real(RP), private, parameter   :: twoOverThree  = 2.0_RP / 3.0_RP
  real(RP), private, parameter   :: FourOverThree = 4.0_RP / 3.0_RP

  real(RP), private              :: Cs            = 0.13_RP ! Smagorinsky constant (Scotti et al. 1993)
  real(RP), private, parameter   :: Ck            = 0.1_RP  ! SGS constant (Moeng and Wyngaard 1988)
  real(RP), private, parameter   :: PrN           = 0.7_RP  ! Prandtl number in neutral conditions
  real(RP), private, parameter   :: RiC           = 0.25_RP ! critical Richardson number
  real(RP), private, parameter   :: FmC           = 16.0_RP ! fum = sqrt(1 - c*Ri)
  real(RP), private, parameter   :: FhB           = 40.0_RP ! fuh = sqrt(1 - b*Ri)/PrN
  real(RP), private              :: RPrN                    ! 1 / PrN
  real(RP), private              :: RRiC                    ! 1 / RiC
  real(RP), private              :: PrNovRiC                ! PrN / RiC

  real(RP), private, allocatable :: nu_fact (:,:,:)         ! (Cs*Delta)^2

  integer,  private              :: I_TKE = -1

  real(RP), private              :: ATMOS_PHY_TB_SMG_NU_MAX     = 10000.0_RP
  logical,  private              :: ATMOS_PHY_TB_SMG_implicit   = .false.
  logical,  private              :: ATMOS_PHY_TB_SMG_bottom     = .false.
  logical,  private              :: ATMOS_PHY_TB_SMG_horizontal = .false.

  real(RP), private              :: tke_fact

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_TB_smg_config( &
       TYPE_TB,  &
       I_TKE_out )
    use scale_process, only: &
       PRC_MPIstop
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    character(len=*), intent(in)  :: TYPE_TB
    integer,          intent(out) :: I_TKE_out
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Turbulence Tracer] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Tracers for Smagorinsky-type Eddy Viscocity Model'

    if ( TYPE_TB /= 'SMAGORINSKY' ) then
       write(*,*) 'xxx ATMOS_PHY_TB_TYPE is not SMAGORINSKY. Check!'
       call PRC_MPIstop
    endif

    call TRACER_regist( I_TKE,                                          & ! [OUT]
                        1,                                              & ! [IN]
                        (/ 'TKE_SMG' /),                                & ! [IN]
                        (/ 'turbulent kinetic energy (Smagorinsky)' /), & ! [IN]
                        (/ 'm2/s2' /),                                  & ! [IN]
                        advc = (/ .false. /)                            ) ! [IN], optional

    I_TKE_out = I_TKE

    return
  end subroutine ATMOS_PHY_TB_smg_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_TB_smg_setup( &
       CDZ, CDX, CDY, CZ )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA,IA,JA)

    real(RP) :: ATMOS_PHY_TB_SMG_Cs
    real(RP) :: ATMOS_PHY_TB_SMG_filter_fact    = 2.0_RP
    logical  :: ATMOS_PHY_TB_SMG_consistent_tke = .true.

    NAMELIST / PARAM_ATMOS_PHY_TB_SMG / &
       ATMOS_PHY_TB_SMG_Cs,             &
       ATMOS_PHY_TB_SMG_NU_MAX,         &
       ATMOS_PHY_TB_SMG_filter_fact,    &
       ATMOS_PHY_TB_SMG_implicit,       &
       ATMOS_PHY_TB_SMG_consistent_tke, &
       ATMOS_PHY_TB_SMG_horizontal,     &
       ATMOS_PHY_TB_SMG_bottom

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Turbulence] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Smagorinsky-type Eddy Viscocity Model'

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
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_TB_SMG)

    Cs = ATMOS_PHY_TB_SMG_Cs

    RPrN     = 1.0_RP / PrN
    RRiC     = 1.0_RP / RiC
    PrNovRiC = ( 1.0_RP - PrN ) * RRiC

    allocate( nu_fact(KA,IA,JA) )

#ifdef DEBUG
    nu_fact (:,:,:) = UNDEF
#endif
    if ( ATMOS_PHY_TB_SMG_horizontal ) then
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          nu_fact(k,i,j) = Cs**2 * ( CDX(i) * CDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
       ATMOS_PHY_TB_SMG_consistent_tke = .false.
       ATMOS_PHY_TB_SMG_implicit       = .false. ! flux in the z-direction is not necessary
    else
       do j = JS-1, JE+1
       do i = IS-1, IE+1
       do k = KS, KE
          nu_fact(k,i,j) = ( Cs * mixlen(CDZ(k),CDX(i),CDY(j),CZ(k,i,j),ATMOS_PHY_TB_SMG_filter_fact) )**2
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif
    end if

    if ( ATMOS_PHY_TB_SMG_consistent_tke ) then
       tke_fact = 1.0_RP
    else
       tke_fact = 0.0_RP
    end if

    return
  end subroutine ATMOS_PHY_TB_smg_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_smg( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, &
       qflx_sgs_rhot, qflx_sgs_rhoq,                &
       RHOQ_t,                                      &
       Nu, Ri, Pr,                                  &
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC, N2,      &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, &
       GSQRT, J13G, J23G, J33G, MAPF, dt            )
    use scale_precision
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       EPS  => CONST_EPS, &
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
       I_UVZ, &
       I_XY,  &
       I_UY,  &
       I_XV,  &
       I_UV
    use scale_atmos_phy_tb_common, only: &
       calc_strain_tensor => ATMOS_PHY_TB_calc_strain_tensor, &
       diffusion_solver   => ATMOS_PHY_TB_diffusion_solver,   &
       calc_tend_momz     => ATMOS_PHY_TB_calc_tend_momz,     &
       calc_tend_momx     => ATMOS_PHY_TB_calc_tend_momx,     &
       calc_tend_momy     => ATMOS_PHY_TB_calc_tend_momy,     &
       calc_tend_phi      => ATMOS_PHY_TB_calc_tend_phi,      &
       calc_flux_phi      => ATMOS_PHY_TB_calc_flux_phi
    implicit none

    ! SGS flux
    real(RP), intent(out)   :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out)   :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

    real(RP), intent(inout) :: RHOQ_t       (KA,IA,JA,QA) ! tendency of rho * QTRC

    real(RP), intent(out)   :: nu           (KA,IA,JA)    ! eddy viscosity (center)
    real(RP), intent(out)   :: Ri           (KA,IA,JA)    ! Richardson number
    real(RP), intent(out)   :: Pr           (KA,IA,JA)    ! Prantle number

    real(RP), intent(in)    :: MOMZ         (KA,IA,JA)
    real(RP), intent(in)    :: MOMX         (KA,IA,JA)
    real(RP), intent(in)    :: MOMY         (KA,IA,JA)
    real(RP), intent(in)    :: RHOT         (KA,IA,JA)
    real(RP), intent(in)    :: DENS         (KA,IA,JA)
    real(RP), intent(in)    :: QTRC         (KA,IA,JA,QA)
    real(RP), intent(in)    :: N2           (KA,IA,JA)

    real(RP), intent(in)    :: SFLX_MW      (IA,JA)
    real(RP), intent(in)    :: SFLX_MU      (IA,JA)
    real(RP), intent(in)    :: SFLX_MV      (IA,JA)
    real(RP), intent(in)    :: SFLX_SH      (IA,JA)
    real(RP), intent(in)    :: SFLX_QV      (IA,JA)

    real(RP), intent(in)    :: GSQRT         (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)    :: J13G          (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G          (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J33G                       !< (3,3) element of Jacobian matrix
    real(RP), intent(in)    :: MAPF(IA,JA,2,4)            !< map factor
    real(DP), intent(in)    :: dt

    ! diagnostic variables
    real(RP) :: TKE  (KA,IA,JA)
    real(RP) :: POTT (KA,IA,JA)

    ! deformation rate tensor
    real(RP) :: S33_C(KA,IA,JA) ! (cell center)
    real(RP) :: S11_C(KA,IA,JA)
    real(RP) :: S22_C(KA,IA,JA)
    real(RP) :: S31_C(KA,IA,JA)
    real(RP) :: S12_C(KA,IA,JA)
    real(RP) :: S23_C(KA,IA,JA)
    real(RP) :: S12_Z(KA,IA,JA) ! (z edge or x-y plane)
    real(RP) :: S23_X(KA,IA,JA) ! (x edge or y-z plane)
    real(RP) :: S31_Y(KA,IA,JA) ! (y edge or z-x plane)
    real(RP) :: S2   (KA,IA,JA) ! |S|^2

    real(RP) :: Kh   (KA,IA,JA) ! eddy diffusion

    real(RP) :: TEND (KA,IA,JA)
    real(RP) :: a    (KA,IA,JA)
    real(RP) :: b    (KA,IA,JA)
    real(RP) :: c    (KA,IA,JA)
    real(RP) :: d    (KA)
    real(RP) :: ap

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Turbulence(smagorinsky)'

#ifdef DEBUG
    qflx_sgs_momz(:,:,:,:)   = UNDEF
    qflx_sgs_momx(:,:,:,:)   = UNDEF
    qflx_sgs_momy(:,:,:,:)   = UNDEF
    qflx_sgs_rhot(:,:,:,:)   = UNDEF
    qflx_sgs_rhoq(:,:,:,:,:) = UNDEF

    nu (:,:,:) = UNDEF
    tke(:,:,:) = UNDEF
    Pr (:,:,:) = UNDEF
    Ri (:,:,:) = UNDEF
    Kh (:,:,:) = UNDEF

    POTT   (:,:,:) = UNDEF
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

    call calc_strain_tensor( &
         S33_C, S11_C, S22_C,          & ! (out)
         S31_C, S12_C, S23_C,          & ! (out)
         S12_Z, S23_X, S31_Y,          & ! (out)
         S2                 ,          & ! (out)
         DENS, MOMZ, MOMX, MOMY,       & ! (in)
         GSQRT, J13G, J23G, J33G, MAPF ) ! (in)

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! Ri = N^2 / |S|^2, N^2 = g / theta * dtheta/dz
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
       call CHECK( __LINE__, POTT(k+1,i,j) )
       call CHECK( __LINE__, POTT(k,i,j) )
       call CHECK( __LINE__, POTT(k-1,i,j) )
       call CHECK( __LINE__, FDZ(k) )
       call CHECK( __LINE__, FDZ(k-1) )
       call CHECK( __LINE__, S2(k,i,j) )
       call CHECK( __LINE__, N2(k,i,j) )
#endif
          Ri(k,i,j) = N2(k,i,j) / S2(k,i,j)
       enddo
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
          Kh(k,i,j) = min( nu(k,i,j)/Pr(k,i,j), ATMOS_PHY_TB_SMG_NU_MAX )
          nu(k,i,j) = min( nu(k,i,j), ATMOS_PHY_TB_SMG_NU_MAX )
          Pr(k,i,j) = nu(k,i,j) / ( Kh(k,i,j) + ( 0.5_RP - sign(0.5_RP, Kh(k,i,j)-EPS) ) )
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
          TKE(k,i,j) = ( nu(k,i,j) * Cs / Ck )**2 / nu_fact(k,i,j)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       !##### momentum equation (z) #####
       ! (cell center)
       if ( ATMOS_PHY_TB_SMG_horizontal ) then
          qflx_sgs_momz(:,:,:,ZDIR) = 0.0_RP
       else
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
                  * ( S33_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree * tke_fact ) &
                  + twoOverThree * tke(k,i,j) * tke_fact )
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
       end if
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

       if ( ATMOS_PHY_TB_SMG_implicit ) then

          call calc_tend_MOMZ( TEND, & ! (out)
                               qflx_sgs_momz, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          do j = JJS, JJE
          do i = IIS, IIE

             ap = - FourOverThree * dt &
                  * DENS(KS+1,i,j)*Nu(KS+1,i,j) &
                  * RCDZ(KS+1) / GSQRT(KS+1,i,j,I_XYZ)
             a(KS,i,j) = ap * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS+1,i,j) )
             do k = KS+1, KE-2
                c(k,i,j) = ap * RFDZ(k+1) / GSQRT(k+1,i,j,I_XYW)
                ap = - FourOverThree * dt &
                     * DENS(k+1,i,j)*Nu(k+1,i,j) &
                     * RCDZ(k+1) / GSQRT(k+1,i,j,I_XYZ)
                a(k,i,j) = ap * RFDZ(k) / GSQRT(k,i,j,I_XYW)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k+1,i,j) )
             end do
             a(KE-1,i,j) = 0.0_RP
             c(KE-1,i,j) = ap * RFDZ(KE) / GSQRT(KE,i,j,I_XYW)
             b(KE-1,i,j) = - c(KE-1,i,j) + 0.5_RP * ( DENS(KE-1,i,j)+DENS(KE,i,j) )

             do k = KS, KE-1
                d(k) = TEND(k,i,j)
             end do

             call diffusion_solver( &
                  TEND(:,i,j),                     & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d, & ! (in)
                  KE-1                             ) ! (in)

             do k = KS+1, KE-1
                qflx_sgs_momz(k,i,j,ZDIR) = qflx_sgs_momz(k,i,j,ZDIR) &
                     - FourOverThree * DENS(k,i,j) * Nu(k,i,j) * dt &
                     * ( TEND(k,i,j) - TEND(k-1,i,j) ) * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
             end do

          end do
          end do

       end if

       !##### momentum equation (x) #####
       ! (y edge)
       if ( ATMOS_PHY_TB_SMG_horizontal ) then
          qflx_sgs_momx(:,:,:,ZDIR) = 0.0_RP
       else
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
       end if
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
       call CHECK( __LINE__, TKE(k,i,j) )
#endif
          qflx_sgs_momx(k,i,j,XDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu(k,i,j) &
               * ( S11_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree * tke_fact ) &
             + twoOverThree * TKE(k,i,j) * tke_fact )
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

       if ( ATMOS_PHY_TB_SMG_implicit ) then
          call calc_tend_MOMX( TEND, & ! (out)
                               qflx_sgs_momx, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS  ,i  ,j)*Nu(KS  ,i  ,j) &
                                   + DENS(KS+1,i  ,j)*Nu(KS+1,i  ,j) &
                                   + DENS(KS  ,i+1,j)*Nu(KS  ,i+1,j) &
                                   + DENS(KS+1,i+1,j)*Nu(KS+1,i+1,j) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_UYW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_UYZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS,i+1,j) )
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_UYZ)
                ap = - dt * 0.25_RP * ( DENS(k  ,i  ,j)*Nu(k  ,i  ,j) &
                                      + DENS(k+1,i  ,j)*Nu(k+1,i  ,j) &
                                      + DENS(k  ,i+1,j)*Nu(k  ,i+1,j) &
                                      + DENS(k+1,i+1,j)*Nu(k+1,i+1,j) ) &
                                    * RFDZ(k) / GSQRT(k,i,j,I_UYW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_UYZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) )
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_UYZ)
             b(KE,i,j) = - c(KE,i,j) + 0.5_RP * ( DENS(KE,i,j)+DENS(KE,i+1,j) )

             do k = KS, KE
                d(k) = TEND(k,i,j)
             end do

             call diffusion_solver( &
                  TEND(:,i,j),                     & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d, & ! (in)
                  KE                               ) ! (in)

             do k = KS, KE-1
                qflx_sgs_momx(k,i,j,ZDIR) = qflx_sgs_momx(k,i,j,ZDIR) &
                     - 0.25_RP * ( DENS(k  ,i  ,j)*Nu(k  ,i  ,j) &
                                 + DENS(k+1,i  ,j)*Nu(k+1,i  ,j) &
                                 + DENS(k  ,i+1,j)*Nu(k  ,i+1,j) &
                                 + DENS(k+1,i+1,j)*Nu(k+1,i+1,j) ) &
                     * dt * ( TEND(k+1,i,j) - TEND(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_UYW)
             end do

          end do
          end do

       end if

       !##### momentum equation (y) #####
       ! (x edge)
       if ( ATMOS_PHY_TB_SMG_horizontal ) then
          qflx_sgs_momy(:,:,:,ZDIR) = 0.0_RP
       else
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
       end if

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
       call CHECK( __LINE__, TKE(k,i,j) )
#endif
          qflx_sgs_momy(k,i,j,YDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * nu(k,i,j) &
               * ( S22_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree * tke_fact ) &
             + twoOverThree * TKE(k,i,j) * tke_fact)
       enddo
       enddo
       enddo
#ifdef DEBUG
       i = IUNDEF; j = IUNDEF; k = IUNDEF
#endif

       if ( ATMOS_PHY_TB_SMG_implicit ) then
          call calc_tend_MOMY( TEND, & ! (out)
                               qflx_sgs_momy, & ! (in)
                               GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
                               IIS, IIE, JJS, JJE ) ! (in)

          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS  ,i,j  )*Nu(KS  ,i,j  ) &
                                   + DENS(KS+1,i,j  )*Nu(KS+1,i,j  ) &
                                   + DENS(KS  ,i,j+1)*Nu(KS  ,i,j+1) &
                                   + DENS(KS+1,i,j+1)*Nu(KS+1,i,j+1) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_XVW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_XVZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + 0.5_RP * ( DENS(KS,i,j)+DENS(KS,i,j+1) )
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XVZ)
                ap = - dt * 0.25_RP * ( DENS(k  ,i,j  )*Nu(k  ,i,j  ) &
                                      + DENS(k+1,i,j  )*Nu(k+1,i,j  ) &
                                      + DENS(k  ,i,j+1)*Nu(k  ,i,j+1) &
                                      + DENS(k+1,i,j+1)*Nu(k+1,i,j+1) ) &
                                    * RFDZ(k) / GSQRT(k,i,j,I_XVW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XVZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) )
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_XVZ)
             b(KE,i,j) = - c(KE,i,j) + 0.5_RP * ( DENS(KE,i,j)+DENS(KE,i,j+1) )

             do k = KS, KE
                d(k) = TEND(k,i,j)
             end do

             call diffusion_solver( &
                  TEND(:,i,j),                     & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d, & ! (in)
                  KE                               ) ! (in)

             do k = KS, KE-1
                qflx_sgs_momy(k,i,j,ZDIR) = qflx_sgs_momy(k,i,j,ZDIR) &
                     - 0.25_RP * ( DENS(k  ,i,j  )*Nu(k  ,i,j  ) &
                                 + DENS(k+1,i,j  )*Nu(k+1,i,j  ) &
                                 + DENS(k  ,i,j+1)*Nu(k  ,i,j+1) &
                                 + DENS(k+1,i,j+1)*Nu(k+1,i,j+1) ) &
                     * dt * ( TEND(k+1,i,j) - TEND(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_XVW)
             end do

          end do
          end do

       end if

       !##### Thermodynamic Equation #####

       if ( ATMOS_PHY_TB_SMG_implicit ) then

          do j = JJS, JJE
          do i = IIS, IIE

             ap = - dt * 0.25_RP * ( DENS(KS,i,j)+DENS(KS+1,i,j) ) &
                                 * ( Kh(KS,i,j)+Kh(KS+1,i,j) ) &
                                 * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
             a(KS,i,j) = ap * RCDZ(KS) / GSQRT(KS,i,j,I_XYZ)
             c(KS,i,j) = 0.0_RP
             b(KS,i,j) = - a(KS,i,j) + DENS(KS,i,j)
             do k = KS+1, KE-1
                c(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
                ap = - dt * 0.25_RP * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                                    * ( Kh(k,i,j)+Kh(k+1,i,j) ) &
                                   * RFDZ(k) / GSQRT(k,i,j,I_XYW)
                a(k,i,j) = ap * RCDZ(k) / GSQRT(k,i,j,I_XYZ)
                b(k,i,j) = - a(k,i,j) - c(k,i,j) + DENS(k,i,j)
             end do
             a(KE,i,j) = 0.0_RP
             c(KE,i,j) = ap * RCDZ(KE) / GSQRT(KE,i,j,I_XYZ)
             b(KE,i,j) = - c(KE,i,j) + DENS(KE,i,j)

          end do
          end do

       end if

       call calc_flux_phi( &
            qflx_sgs_rhot, &
            DENS, POTT, Kh, 1.0_RP, &
            GSQRT, J13G, J23G, J33G, MAPF, &
            a, b, c, dt, &
            ATMOS_PHY_TB_SMG_implicit, &
            IIS, IIE, JJS, JJE )

    enddo
    enddo


    !##### Tracers #####
    do iq = 1, QA

       if ( iq == I_TKE .or. .not. TRACER_ADVC(iq) ) cycle

       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1

          call calc_flux_phi( &
               qflx_sgs_rhoq(:,:,:,:,iq), &
               DENS, QTRC(:,:,:,iq), Kh, 1.0_RP, &
               GSQRT, J13G, J23G, J33G, MAPF, &
               a, b, c, dt, &
               ATMOS_PHY_TB_SMG_implicit, &
               IIS, IIE, JJS, JJE )


       enddo
       enddo
#ifdef DEBUG
       IIS = IUNDEF; IIE = IUNDEF; JJS = IUNDEF; JJE = IUNDEF
#endif

    enddo ! scalar quantities loop
#ifdef DEBUG
    iq = IUNDEF
#endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOQ_t(k,i,j,I_TKE) = ( TKE(k,i,j) - QTRC(k,i,j,I_TKE) ) * DENS(k,i,j) / dt
    end do
    end do
    end do

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
    if ( dz == dmax ) then
       a1 = dx / dmax
       a2 = dy / dmax
    else if ( dx == dmax ) then
       a1 = dz / dmax
       a2 = dy / dmax
    else ! dy == dmax
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
