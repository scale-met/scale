!-------------------------------------------------------------------------------
!> module URBAN / Surface fluxes with Single-layer Canpoy Model
!!
!! @par Description
!!          Surface fluxes between atmosphere and urban
!!          based on Single-layer Urban Canopy Model (Kusaka et al. 2001, BLM)
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
module scale_urban_phy_slc
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: URBAN_PHY_SLC_setup
  public :: URBAN_PHY_SLC

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: SLC_main
  private :: canopy_wind
  private :: cal_beta
  private :: cal_psi
  private :: mos
  private :: multi_layer
  private :: urban_param_setup

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  logical, allocatable, private :: is_URB(:,:) ! urban tile or not.
  ! from namelist
  real(RP), private :: DTS_MAX    =    0.1_RP ! maximum dT during one minute [K/sec]
                                              ! 0.1 [K/sec] = 6.0 [K/min]
  real(RP), private :: ZR         =   10.0_RP ! roof level ( building height) [m]
  real(RP), private :: roof_width =    9.0_RP ! roof level ( building height) [m]
  real(RP), private :: road_width =   11.0_RP ! roof level ( building height) [m]
  real(RP), private :: SIGMA_ZED  =    1.0_RP ! Standard deviation of roof height [m]
  real(RP), private :: AH         =   17.5_RP ! Sensible Anthropogenic heat [W/m^2]
  real(RP), private :: ALH        =    0.0_RP ! Latent Anthropogenic heat   [W/m^2]
  real(RP), private :: BETR       =    0.0_RP ! Evaporation efficiency of roof     [-]
  real(RP), private :: BETB       =    0.0_RP !                        of building [-]
  real(RP), private :: BETG       =    0.0_RP !                        of ground   [-]
  real(RP), private :: STRGR      =    0.0_RP ! rain strage on roof     [-]
  real(RP), private :: STRGB      =    0.0_RP !             on wall     [-]
  real(RP), private :: STRGG      =    0.0_RP !             on ground   [-]
  real(RP), private :: CAPR       =  1.2E6_RP ! heat capacity of roof   [J m-3 K]
  real(RP), private :: CAPB       =  1.2E6_RP !               of wall   [J m-3 K]
  real(RP), private :: CAPG       =  1.2E6_RP !               of ground [J m-3 K]
  real(RP), private :: AKSR       =   2.28_RP ! thermal conductivity of roof   [W m-1 K]
  real(RP), private :: AKSB       =   2.28_RP !                      of wall   [W m-1 K]
  real(RP), private :: AKSG       =   2.28_RP !                      of ground [W m-1 K]
  real(RP), private :: ALBR       =    0.2_RP ! surface albedo of roof
  real(RP), private :: ALBB       =    0.2_RP ! surface albedo of wall
  real(RP), private :: ALBG       =    0.2_RP ! surface albedo of ground
  real(RP), private :: EPSR       =   0.90_RP ! Surface emissivity of roof
  real(RP), private :: EPSB       =   0.90_RP ! Surface emissivity of wall
  real(RP), private :: EPSG       =   0.90_RP ! Surface emissivity of ground
  real(RP), private :: Z0R        =   0.01_RP ! roughness length for momentum of building roof
  real(RP), private :: Z0B        = 0.0001_RP ! roughness length for momentum of building wall
  real(RP), private :: Z0G        =   0.01_RP ! roughness length for momentum of ground
  real(RP), private :: TRLEND     = 293.00_RP ! lower boundary condition of roof temperature [K]
  real(RP), private :: TBLEND     = 293.00_RP ! lower boundary condition of wall temperature [K]
  real(RP), private :: TGLEND     = 293.00_RP ! lower boundary condition of ground temperature [K]
  integer , private :: BOUND      = 1         ! Boundary Condition for Roof, Wall, Ground Layer Temp
                                              !       [1: Zero-Flux, 2: T = Constant]
  real(RP), private :: ahdiurnal(1:24)        ! AH diurnal profile

  ! calculated in subroutine urban_param_set
  real(RP), private :: R                       ! Normalized roof wight (eq. building coverage ratio)
  real(RP), private :: RW                      ! (= 1 - R)
  real(RP), private :: HGT                     ! Normalized building height
  real(RP), private :: Z0HR                    ! roughness length for heat of roof
  real(RP), private :: Z0HB                    ! roughness length for heat of building wall
  real(RP), private :: Z0HG                    ! roughness length for heat of ground
  real(RP), private :: Z0C                     ! Roughness length above canyon for momentum [m]
  real(RP), private :: Z0HC                    ! Roughness length above canyon for heat [m]
  real(RP), private :: ZDC                     ! Displacement height [m]
  real(RP), private :: SVF                     ! Sky view factor [-]

  real(RP), private, allocatable :: DZR(:)     ! thickness of each roof layer [m]
  real(RP), private, allocatable :: DZB(:)     ! thickness of each building layer [m]
  real(RP), private, allocatable :: DZG(:)     ! thickness of each road layer [m]

  real(RP), private :: XXXR    = 0.0_RP        ! Monin-Obkhov length for roof [-]
  real(RP), private :: XXXC    = 0.0_RP        ! Monin-Obkhov length for canopy [-]
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine URBAN_PHY_SLC_setup( &
       URBAN_TYPE, &
       Z0M, &
       Z0H, &
       Z0E  )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_landuse, only: &
       LANDUSE_fact_urban
    implicit none

    character(len=*), intent(in)  :: URBAN_TYPE
    real(RP)        , intent(out) :: Z0M(IA,JA)
    real(RP)        , intent(out) :: Z0H(IA,JA)
    real(RP)        , intent(out) :: Z0E(IA,JA)

    NAMELIST / PARAM_URBAN_PHY_SLC / &
       DTS_MAX,    &
       ZR,         &
       roof_width, &
       road_width, &
       SIGMA_ZED,  &
       AH,         &
       ALH,        &
       BETR,       &
       BETB,       &
       BETG,       &
       STRGR,      &
       STRGB,      &
       STRGG,      &
       CAPR,       &
       CAPB,       &
       CAPG,       &
       AKSR,       &
       AKSB,       &
       AKSG,       &
       ALBR,       &
       ALBB,       &
       ALBG,       &
       EPSR,       &
       EPSB,       &
       EPSG,       &
       Z0R,        &
       Z0B,        &
       Z0G,        &
       TRLEND,     &
       TBLEND,     &
       TGLEND,     &
       BOUND

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[SLC] / Categ[URBAN PHY] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN_PHY_SLC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_URBAN_PHY_SLC. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_URBAN_PHY_SLC)

    allocate( DZR(UKS:UKE) )
    allocate( DZB(UKS:UKE) )
    allocate( DZG(UKS:UKE) )
    DZR(UKS:UKE) = (/0.01_RP,0.01_RP,0.03_RP,0.05_RP,0.10_RP/)
    DZB(UKS:UKE) = (/0.01_RP,0.01_RP,0.03_RP,0.05_RP,0.10_RP/)
    DZG(UKS:UKE) = (/0.01_RP,0.01_RP,0.03_RP,0.05_RP,0.10_RP/)

    ahdiurnal(:) = (/ 0.356, 0.274, 0.232, 0.251, 0.375, 0.647, 0.919, 1.135, 1.249, 1.328, &
                      1.365, 1.363, 1.375, 1.404, 1.457, 1.526, 1.557, 1.521, 1.372, 1.206, &
                      1.017, 0.876, 0.684, 0.512                                            /)

    ! set other urban parameters
    call urban_param_setup

    ! judge to run slab land model
    allocate( is_URB(IA,JA) )

    do j = JS, JE
    do i = IS, IE
       if ( LANDUSE_fact_urban(i,j) > 0.0_RP ) then
          is_URB(i,j) = .true.
       else
          is_URB(i,j) = .false.
       endif
    enddo
    enddo

    Z0M(:,:) = UNDEF
    Z0H(:,:) = UNDEF
    Z0E(:,:) = UNDEF
    do j = JS, JE
    do i = IS, IE
       if ( is_URB(i,j) ) then
          Z0M(i,j) = Z0C
          Z0H(i,j) = Z0HC
          Z0E(i,j) = Z0HC
       endif
    enddo
    enddo

    return
  end subroutine URBAN_PHY_SLC_setup

  subroutine URBAN_PHY_SLC( &
        TR_URB_t,    &
        TB_URB_t,    &
        TG_URB_t,    &
        TC_URB_t,    &
        QC_URB_t,    &
        UC_URB_t,    &
        TRL_URB_t,   &
        TBL_URB_t,   &
        TGL_URB_t,   &
        RAINR_URB_t, &
        RAINB_URB_t, &
        RAING_URB_t, &
        ROFF_URB_t,  &
        SFC_TEMP,    &
        ALBD_LW,     &
        ALBD_SW,     &
        MWFLX,       &
        MUFLX,       &
        MVFLX,       &
        SHFLX,       &
        LHFLX,       &
        GHFLX,       &
        Z0M,         &
        Z0H,         &
        Z0E,         &
        U10,         &
        V10,         &
        T2,          &
        Q2,          &
        TMPA,        &
        PRSA,        &
        W1,          &
        U1,          &
        V1,          &
        DENS,        &
        QA,          &
        Z1,          &
        PBL,         &
        PRSS,        &
        LWD,         &
        SWD,         &
        PREC,        &
        TR_URB,      &
        TB_URB,      &
        TG_URB,      &
        TC_URB,      &
        QC_URB,      &
        UC_URB,      &
        TRL_URB,     &
        TBL_URB,     &
        TGL_URB,     &
        RAINR_URB,   &
        RAINB_URB,   &
        RAING_URB,   &
        ROFF_URB,    &
        LON,         &
        LAT,         &
        NOWDATE,     &
        dt           )
    use scale_grid_index
    use scale_urban_grid_index
    use scale_history, only: &
       HIST_in
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
       BULKFLUX
    implicit none

    ! parameter
    logical,  parameter :: LSOLAR = .false. ! [true=both, false=SSG only]

    real(RP), parameter :: Uabs_min = 0.1_RP

    ! arguments
    real(RP), intent(out) :: TR_URB_t   (IA,JA)
    real(RP), intent(out) :: TB_URB_t   (IA,JA)
    real(RP), intent(out) :: TG_URB_t   (IA,JA)
    real(RP), intent(out) :: TC_URB_t   (IA,JA)
    real(RP), intent(out) :: QC_URB_t   (IA,JA)
    real(RP), intent(out) :: UC_URB_t   (IA,JA)
    real(RP), intent(out) :: TRL_URB_t  (UKS:UKE,IA,JA)
    real(RP), intent(out) :: TBL_URB_t  (UKS:UKE,IA,JA)
    real(RP), intent(out) :: TGL_URB_t  (UKS:UKE,IA,JA)
    real(RP), intent(out) :: RAINR_URB_t(IA,JA)
    real(RP), intent(out) :: RAINB_URB_t(IA,JA)
    real(RP), intent(out) :: RAING_URB_t(IA,JA)
    real(RP), intent(out) :: ROFF_URB_t (IA,JA)

    real(RP), intent(out) :: SFC_TEMP(IA,JA)
    real(RP), intent(out) :: ALBD_LW (IA,JA)
    real(RP), intent(out) :: ALBD_SW (IA,JA)
    real(RP), intent(out) :: MWFLX   (IA,JA)
    real(RP), intent(out) :: MUFLX   (IA,JA)
    real(RP), intent(out) :: MVFLX   (IA,JA)
    real(RP), intent(out) :: SHFLX   (IA,JA)
    real(RP), intent(out) :: LHFLX   (IA,JA)
    real(RP), intent(out) :: GHFLX   (IA,JA)
    real(RP), intent(out) :: Z0M     (IA,JA)
    real(RP), intent(out) :: Z0H     (IA,JA)
    real(RP), intent(out) :: Z0E     (IA,JA)
    real(RP), intent(out) :: U10     (IA,JA)
    real(RP), intent(out) :: V10     (IA,JA)
    real(RP), intent(out) :: T2      (IA,JA)
    real(RP), intent(out) :: Q2      (IA,JA)

    real(RP), intent(in) :: TMPA  (IA,JA)
    real(RP), intent(in) :: PRSA  (IA,JA)
    real(RP), intent(in) :: W1    (IA,JA)
    real(RP), intent(in) :: U1    (IA,JA)
    real(RP), intent(in) :: V1    (IA,JA)
    real(RP), intent(in) :: DENS  (IA,JA)
    real(RP), intent(in) :: QA    (IA,JA)
    real(RP), intent(in) :: Z1    (IA,JA)
    real(RP), intent(in) :: PBL   (IA,JA)
    real(RP), intent(in) :: PRSS  (IA,JA)
    real(RP), intent(in) :: LWD   (IA,JA,2)
    real(RP), intent(in) :: SWD   (IA,JA,2)
    real(RP), intent(in) :: PREC  (IA,JA)

    real(RP), intent(in) :: TR_URB   (IA,JA)
    real(RP), intent(in) :: TB_URB   (IA,JA)
    real(RP), intent(in) :: TG_URB   (IA,JA)
    real(RP), intent(in) :: TC_URB   (IA,JA)
    real(RP), intent(in) :: QC_URB   (IA,JA)
    real(RP), intent(in) :: UC_URB   (IA,JA)
    real(RP), intent(in) :: TRL_URB  (UKS:UKE,IA,JA)
    real(RP), intent(in) :: TBL_URB  (UKS:UKE,IA,JA)
    real(RP), intent(in) :: TGL_URB  (UKS:UKE,IA,JA)
    real(RP), intent(in) :: RAINR_URB(IA,JA)
    real(RP), intent(in) :: RAINB_URB(IA,JA)
    real(RP), intent(in) :: RAING_URB(IA,JA)
    real(RP), intent(in) :: ROFF_URB (IA,JA)

    real(RP), intent(in) :: LON
    real(RP), intent(in) :: LAT
    integer,  intent(in) :: NOWDATE(6)
    real(DP), intent(in) :: dt

    ! work
    real(RP) :: TR
    real(RP) :: TB
    real(RP) :: TG
    real(RP) :: TC
    real(RP) :: QC
    real(RP) :: UC
    real(RP) :: TRL(UKS:UKE)
    real(RP) :: TBL(UKS:UKE)
    real(RP) :: TGL(UKS:UKE)
    real(RP) :: RAINR
    real(RP) :: RAINB
    real(RP) :: RAING
    real(RP) :: ROFF

    real(RP) :: SHR  (IA,JA)
    real(RP) :: SHB  (IA,JA)
    real(RP) :: SHG  (IA,JA)
    real(RP) :: LHR  (IA,JA)
    real(RP) :: LHB  (IA,JA)
    real(RP) :: LHG  (IA,JA)
    real(RP) :: GHR  (IA,JA)
    real(RP) :: GHB  (IA,JA)
    real(RP) :: GHG  (IA,JA)
    real(RP) :: RNR  (IA,JA)
    real(RP) :: RNB  (IA,JA)
    real(RP) :: RNG  (IA,JA)
    real(RP) :: RNgrd(IA,JA)

    real(RP) :: Ustar ! friction velocity [m]
    real(RP) :: Tstar ! friction temperature [K]
    real(RP) :: Qstar ! friction mixing rate [kg/kg]
    real(RP) :: Uabs  ! modified absolute velocity [m/s]

    real(RP) :: QVsat ! saturation water vapor mixing ratio at surface [kg/kg]

    real(RP) :: FracU10 ! calculation parameter for U10 [-]
    real(RP) :: FracT2  ! calculation parameter for T2 [-]
    real(RP) :: FracQ2  ! calculation parameter for Q2 [-]

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Urban surface physics step: Single Layer Canopy'

    do j = JS, JE
    do i = IS, IE

    if( is_URB(i,j) ) then

       Uabs = max( sqrt( U1(i,j)**2 + V1(i,j)**2 + W1(i,j)**2 ), Uabs_min )

       ! save
       TR = TR_URB(i,j)
       TB = TB_URB(i,j)
       TG = TG_URB(i,j)
       TC = TC_URB(i,j)
       QC = QC_URB(i,j)
       UC = UC_URB(i,j)

       do k = UKS, UKE
          TRL(k) = TRL_URB(k,i,j)
          TBL(k) = TBL_URB(k,i,j)
          TGL(k) = TGL_URB(k,i,j)
       end do

       RAINR = RAINR_URB(i,j)
       RAINB = RAINB_URB(i,j)
       RAING = RAING_URB(i,j)
       ROFF  = ROFF_URB (i,j)

       call SLC_main( TR,                 & ! [INOUT]
                      TB,                 & ! [INOUT]
                      TG,                 & ! [INOUT]
                      TC,                 & ! [INOUT]
                      QC,                 & ! [INOUT]
                      UC,                 & ! [INOUT]
                      TRL     (:),        & ! [INOUT]
                      TBL     (:),        & ! [INOUT]
                      TGL     (:),        & ! [INOUT]
                      RAINR,              & ! [INOUT]
                      RAINB,              & ! [INOUT]
                      RAING,              & ! [INOUT]
                      ROFF,               & ! [INOUT]
                      ALBD_LW (i,j),      & ! [OUT]
                      ALBD_SW (i,j),      & ! [OUT]
                      SHR     (i,j),      & ! [OUT]
                      SHB     (i,j),      & ! [OUT]
                      SHG     (i,j),      & ! [OUT]
                      LHR     (i,j),      & ! [OUT]
                      LHB     (i,j),      & ! [OUT]
                      LHG     (i,j),      & ! [OUT]
                      GHR     (i,j),      & ! [OUT]
                      GHB     (i,j),      & ! [OUT]
                      GHG     (i,j),      & ! [OUT]
                      RNR     (i,j),      & ! [OUT]
                      RNB     (i,j),      & ! [OUT]
                      RNG     (i,j),      & ! [OUT]
                      SFC_TEMP(i,j),      & ! [OUT]
                      RNgrd   (i,j),      & ! [OUT]
                      SHFLX   (i,j),      & ! [OUT]
                      LHFLX   (i,j),      & ! [OUT]
                      GHFLX   (i,j),      & ! [OUT]
                      U10     (i,j),      & ! [OUT]
                      V10     (i,j),      & ! [OUT]
                      T2      (i,j),      & ! [OUT]
                      Q2      (i,j),      & ! [OUT]
                      LSOLAR,             & ! [IN]
                      PRSA    (i,j),      & ! [IN]
                      PRSS    (i,j),      & ! [IN]
                      TMPA    (i,j),      & ! [IN]
                      QA      (i,j),      & ! [IN]
                      Uabs,               & ! [IN]
                      U1      (i,j),      & ! [IN]
                      V1      (i,j),      & ! [IN]
                      Z1      (i,j),      & ! [IN]
                      SWD     (i,j,:),    & ! [IN]
                      LWD     (i,j,:),    & ! [IN]
                      PREC    (i,j),      & ! [IN]
                      DENS    (i,j),      & ! [IN]
                      LON,                & ! [IN]
                      LAT,                & ! [IN]
                      NOWDATE (:),        & ! [IN]
                      dt, i, j            ) ! [IN]

       ! calculate tendency
       TR_URB_t(i,j) = ( TR - TR_URB(i,j) ) / dt
       TB_URB_t(i,j) = ( TB - TB_URB(i,j) ) / dt
       TG_URB_t(i,j) = ( TG - TG_URB(i,j) ) / dt
       TC_URB_t(i,j) = ( TC - TC_URB(i,j) ) / dt
       QC_URB_t(i,j) = ( QC - QC_URB(i,j) ) / dt
       UC_URB_t(i,j) = ( UC - UC_URB(i,j) ) / dt

       do k = UKS, UKE
          TRL_URB_t(k,i,j) = ( TRL(k) - TRL_URB(k,i,j) ) / dt
          TBL_URB_t(k,i,j) = ( TBL(k) - TBL_URB(k,i,j) ) / dt
          TGL_URB_t(k,i,j) = ( TGL(k) - TGL_URB(k,i,j) ) / dt
       end do

       RAINR_URB_t(i,j) = ( RAINR - RAINR_URB(i,j) ) / dt
       RAINB_URB_t(i,j) = ( RAINB - RAINB_URB(i,j) ) / dt
       RAING_URB_t(i,j) = ( RAING - RAING_URB(i,j) ) / dt
       ROFF_URB_t (i,j) = ( ROFF  - ROFF_URB (i,j) ) / dt

       ! saturation at the surface
       call qsat( QVsat, SFC_TEMP(i,j), PRSS(i,j) )

       call BULKFLUX( Ustar,         & ! [OUT]
                      Tstar,         & ! [OUT]
                      Qstar,         & ! [OUT]
                      Uabs,          & ! [OUT]
                      FracU10,       & ! [OUT]
                      FracT2,        & ! [OUT]
                      FracQ2,        & ! [OUT]
                      TMPA    (i,j), & ! [IN]
                      SFC_TEMP(i,j), & ! [IN]
                      PRSA    (i,j), & ! [IN]
                      PRSS    (i,j), & ! [IN]
                      QA      (i,j), & ! [IN]
                      QVsat,         & ! [IN]
                      U1      (i,j), & ! [IN]
                      V1      (i,j), & ! [IN]
                      Z1      (i,j), & ! [IN]
                      PBL     (i,j), & ! [IN]
                      Z0C,           & ! [IN]
                      Z0HC,          & ! [IN]
                      Z0HC           ) ! [IN]

       MWFLX(i,j) = -DENS(i,j) * Ustar**2 / Uabs * W1(i,j)
       MUFLX(i,j) = -DENS(i,j) * Ustar**2 / Uabs * U1(i,j)
       MVFLX(i,j) = -DENS(i,j) * Ustar**2 / Uabs * V1(i,j)

       Z0M(i,j) = Z0C
       Z0H(i,j) = Z0HC
       Z0E(i,j) = Z0HC

    else
       ! not calculate urban module
       TR_URB_t   (i,j)   = 0.0_RP
       TB_URB_t   (i,j)   = 0.0_RP
       TG_URB_t   (i,j)   = 0.0_RP
       TC_URB_t   (i,j)   = 0.0_RP
       QC_URB_t   (i,j)   = 0.0_RP
       UC_URB_t   (i,j)   = 0.0_RP
       TRL_URB_t  (:,i,j) = 0.0_RP
       TBL_URB_t  (:,i,j) = 0.0_RP
       TGL_URB_t  (:,i,j) = 0.0_RP
       RAINR_URB_t(i,j)   = 0.0_RP
       RAINB_URB_t(i,j)   = 0.0_RP
       RAING_URB_t(i,j)   = 0.0_RP
       ROFF_URB_t (i,j)   = 0.0_RP
       SFC_TEMP   (i,j)   = 300.0_RP ! constant value
       ALBD_LW    (i,j)   = 0.0_RP
       ALBD_SW    (i,j)   = 0.0_RP
       MWFLX      (i,j)   = 0.0_RP
       MUFLX      (i,j)   = 0.0_RP
       MVFLX      (i,j)   = 0.0_RP
       SHFLX      (i,j)   = 0.0_RP
       LHFLX      (i,j)   = 0.0_RP
       GHFLX      (i,j)   = 0.0_RP
       Z0M        (i,j)   = 0.0_RP
       Z0H        (i,j)   = 0.0_RP
       Z0E        (i,j)   = 0.0_RP
       U10        (i,j)   = 0.0_RP
       V10        (i,j)   = 0.0_RP
       T2         (i,j)   = 0.0_RP
       Q2         (i,j)   = 0.0_RP
       !
       SHR        (i,j)   = 0.0_RP
       SHB        (i,j)   = 0.0_RP
       SHG        (i,j)   = 0.0_RP
       LHR        (i,j)   = 0.0_RP
       LHB        (i,j)   = 0.0_RP
       LHG        (i,j)   = 0.0_RP
       GHR        (i,j)   = 0.0_RP
       GHB        (i,j)   = 0.0_RP
       GHG        (i,j)   = 0.0_RP
       RNR        (i,j)   = 0.0_RP
       RNB        (i,j)   = 0.0_RP
       RNG        (i,j)   = 0.0_RP
       RNgrd      (i,j)   = 0.0_RP
    endif

    end do
    end do

    call HIST_in( SHR  (:,:), 'URBAN_SHR',   'urban sensible heat flux on roof',    'W/m2' )
    call HIST_in( SHB  (:,:), 'URBAN_SHB',   'urban sensible heat flux on wall',    'W/m2' )
    call HIST_in( SHG  (:,:), 'URBAN_SHG',   'urban sensible heat flux on road',    'W/m2' )
    call HIST_in( LHR  (:,:), 'URBAN_LHR',   'urban latent heat flux on roof',      'W/m2' )
    call HIST_in( LHB  (:,:), 'URBAN_LHB',   'urban latent heat flux on wall',      'W/m2' )
    call HIST_in( LHG  (:,:), 'URBAN_LHG',   'urban latent heat flux on road',      'W/m2' )
    call HIST_in( GHR  (:,:), 'URBAN_GHR',   'urban ground heat flux on roof',      'W/m2' )
    call HIST_in( GHB  (:,:), 'URBAN_GHB',   'urban ground heat flux on wall',      'W/m2' )
    call HIST_in( GHG  (:,:), 'URBAN_GHG',   'urban ground heat flux on road',      'W/m2' )
    call HIST_in( RNR  (:,:), 'URBAN_RNR',   'urban net radiation on roof',         'W/m2' )
    call HIST_in( RNB  (:,:), 'URBAN_RNB',   'urban net radiation on wall',         'W/m2' )
    call HIST_in( RNG  (:,:), 'URBAN_RNG',   'urban net radiation on road',         'W/m2' )
    call HIST_in( RNgrd(:,:), 'URBAN_RNgrd', 'urban grid average of net radiation', 'W/m2' )

    return
  end subroutine URBAN_PHY_SLC

  !-----------------------------------------------------------------------------
  subroutine SLC_main( &
        TR,           & ! (inout)
        TB,           & ! (inout)
        TG,           & ! (inout)
        TC,           & ! (inout)
        QC,           & ! (inout)
        UC,           & ! (inout)
        TRL,          & ! (inout)
        TBL,          & ! (inout)
        TGL,          & ! (inout)
        RAINR,        & ! (inout)
        RAINB,        & ! (inout)
        RAING,        & ! (inout)
        ROFF,         & ! (inout)
        ALBD_LW_grid, & ! (out)
        ALBD_SW_grid, & ! (out)
        SHR,          & ! (out)
        SHB,          & ! (out)
        SHG,          & ! (out)
        LHR,          & ! (out)
        LHB,          & ! (out)
        LHG,          & ! (out)
        GHR,          & ! (out)
        GHB,          & ! (out)
        GHG,          & ! (out)
        RNR,          & ! (out)
        RNB,          & ! (out)
        RNG,          & ! (out)
        RTS,          & ! (out)
        RN,           & ! (out)
        SH,           & ! (out)
        LH,           & ! (out)
        GHFLX,        & ! (out)
        U10,          & ! (out)
        V10,          & ! (out)
        T2,           & ! (out)
        Q2,           & ! (out)
        LSOLAR,       & ! (in)
        PRSA,         & ! (in)
        PRSS,         & ! (in)
        TA,           & ! (in)
        QA,           & ! (in)
        UA,           & ! (in)
        U1,           & ! (in)
        V1,           & ! (in)
        ZA,           & ! (in)
        SSG,          & ! (in)
        LLG,          & ! (in)
        RAIN,         & ! (in)
        RHOO,         & ! (in)
        XLON,         & ! (in)
        XLAT,         & ! (in)
        NOWDATE,      & ! (in)
        dt,           & ! (in)
        i, j          ) ! (in)
    use scale_grid_index
    use scale_process, only: &
       PRC_myrank, &
       PRC_MPIstop
    use scale_const, only: &
       EPS    => CONST_EPS,     &    ! small number (machine epsilon)
       PI     => CONST_PI,      &    ! pi               [-]
       D2R    => CONST_D2R,     &    ! degree to radian
       KARMAN => CONST_KARMAN,  &    ! Kalman constant  [-]
       CPdry  => CONST_CPdry,   &    ! Heat capacity of dry air [J/K/kg]
       GRAV   => CONST_GRAV,    &    ! gravitational constant [m/s2]
       Rdry   => CONST_Rdry,    &    ! specific gas constant (dry) [J/kg/K]
       Rvap   => CONST_Rvap,    &    ! gas constant (water vapor) [J/kg/K]
       STB    => CONST_STB,     &    ! stefan-Boltzman constant [MKS unit]
       TEM00  => CONST_TEM00,   &    ! temperature reference (0 degree C) [K]
       PRE00  => CONST_PRE00         ! pressure reference [Pa]
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_pres2qsat_all
    implicit none

    !-- configuration variables
    logical , intent(in)    :: LSOLAR ! logical   [true=both, false=SSG only]

    !-- Input variables from Coupler to Urban
    real(RP), intent(in)    :: PRSA ! Pressure at 1st atmospheric layer      [Pa]
    real(RP), intent(in)    :: PRSS ! Surface Pressure                       [Pa]
    real(RP), intent(in)    :: TA   ! temp at 1st atmospheric level          [K]
    real(RP), intent(in)    :: QA   ! specific humidity at 1st atmospheric level  [kg/kg]
    real(RP), intent(in)    :: UA   ! wind speed at 1st atmospheric level    [m/s]
    real(RP), intent(in)    :: U1   ! u at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: V1   ! v at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: ZA   ! height of 1st atmospheric level        [m]
    real(RP), intent(in)    :: SSG(2) ! downward total short wave radiation    [W/m/m]
    real(RP), intent(in)    :: LLG(2) ! downward long wave radiation           [W/m/m]
    real(RP), intent(in)    :: RAIN ! precipitation flux                     [kg/m2/s]
    real(RP), intent(in)    :: RHOO ! air density                            [kg/m^3]
    real(RP), intent(in)    :: XLAT ! latitude                               [rad,-pi,pi]
    real(RP), intent(in)    :: XLON ! longitude                              [rad,0-2pi]
    integer,  intent(in)    :: NOWDATE(6)
    real(DP), intent(in)    :: dt

    !-- In/Out variables from/to Coupler to/from Urban
    real(RP), intent(inout) :: TR   ! roof temperature              [K]
    real(RP), intent(inout) :: TB   ! building wall temperature     [K]
    real(RP), intent(inout) :: TG   ! road temperature              [K]
    real(RP), intent(inout) :: TC   ! urban-canopy air temperature  [K]
    real(RP), intent(inout) :: QC   ! urban-canopy air specific humidity [kg/kg]
    real(RP), intent(inout) :: UC   ! diagnostic canopy wind        [m/s]
    real(RP), intent(inout) :: TRL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: TBL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: TGL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: RAINR ! rain amount in storage on roof     [kg/m2]
    real(RP), intent(inout) :: RAINB ! rain amount in storage on building [kg/m2]
    real(RP), intent(inout) :: RAING ! rain amount in storage on road     [kg/m2]
    real(RP), intent(inout) :: ROFF  ! runoff from urban           [kg/m2]

    !-- Output variables from Urban to Coupler
    real(RP), intent(out)   :: ALBD_SW_grid  ! grid mean of surface albedo for SW
    real(RP), intent(out)   :: ALBD_LW_grid  ! grid mean of surface albedo for LW ( 1-emiss )
    real(RP), intent(out)   :: RTS    ! radiative surface temperature    [K]
    real(RP), intent(out)   :: SH     ! sensible heat flux               [W/m/m]
    real(RP), intent(out)   :: LH     ! latent heat flux                 [W/m/m]
    real(RP), intent(out)   :: GHFLX  ! heat flux into the ground        [W/m/m]
    real(RP), intent(out)   :: RN     ! net radition                     [W/m/m]
    real(RP), intent(out)   :: U10    ! U wind at 10m                    [m/s]
    real(RP), intent(out)   :: V10    ! V wind at 10m                    [m/s]
    real(RP), intent(out)   :: T2     ! air temperature at 2m            [K]
    real(RP), intent(out)   :: Q2     ! specific humidity at 2m          [kg/kg]
    real(RP), intent(out)   :: RNR, RNB, RNG
    real(RP), intent(out)   :: SHR, SHB, SHG
    real(RP), intent(out)   :: LHR, LHB, LHG
    real(RP), intent(out)   :: GHR, GHB, GHG
    integer , intent(in)    :: i, j

    !-- parameters
!    real(RP), parameter     :: SRATIO    = 0.75_RP     ! ratio between direct/total solar [-]
!    real(RP), parameter     :: TFa       = 0.5_RP      ! factor a in Tomita (2009)
!    real(RP), parameter     :: TFb       = 1.1_RP      ! factor b in Tomita (2009)
    real(RP), parameter     :: redf_min  = 1.0E-2_RP   ! minimum reduced factor
    real(RP), parameter     :: redf_max  = 1.0_RP      ! maximum reduced factor
!    real(RP), parameter     :: CAP_water = 4.185E6_RP ! Heat capacity of water (15 deg) [J m-3 K]
!    real(RP), parameter     :: AKS_water = 0.59_RP    ! Thermal conductivity of water   [W m-1 K]

    !-- Local variables
!    logical  :: SHADOW = .false.
           !  true  = consider svf and shadow effects,
           !  false = consider svf effect only

    real(RP) :: LON,LAT  ! longitude [deg], latitude [deg]
    integer  :: tloc     ! local time (1-24h)
    real(RP) :: dsec     ! second [s]
    real(RP) :: TIME     ! absorute part of current time
    real(RP) :: tahdiurnal ! temporal AH diurnal profile
    real(RP) :: AH_t     ! Sensible Anthropogenic heat [W/m^2]
    real(RP) :: ALH_t    ! Latent Anthropogenic heat   [W/m^2]

    real(RP) :: SSGD     ! downward direct short wave radiation   [W/m/m]
    real(RP) :: SSGQ     ! downward diffuse short wave radiation  [W/m/m]

    real(RP) :: W, VFGS, VFGW, VFWG, VFWS, VFWW
    real(RP) :: SX, RX

    real(RP) :: TRP              ! TRP: at previous time step [K]
    real(RP) :: TBP              ! TBP: at previous time step [K]
    real(RP) :: TGP              ! TGP: at previous time step [K]
    real(RP) :: TCP              ! TCP: at previous time step [K]
    real(RP) :: QCP              ! QCP: at previous time step [kg/kg]
    real(RP) :: TRLP(UKS:UKE)    ! Layer temperature at previous step  [K]
    real(RP) :: TBLP(UKS:UKE)    ! Layer temperature at previous step  [K]
    real(RP) :: TGLP(UKS:UKE)    ! Layer temperature at previous step  [K]
    !
    real(RP) :: RAINRP ! at previous step, rain amount in storage on roof     [kg/m2]
    real(RP) :: RAINBP ! at previous step, rain amount in storage on building [kg/m2]
    real(RP) :: RAINGP ! at previous step, rain amount in storage on road     [kg/m2]

    !real(RP) :: UST, TST, QST

    real(RP) :: RAINT
    real(RP) :: ROFFR, ROFFB, ROFFG ! runoff [kg/m2]

    real(RP) :: LUP, LDN, RUP
    real(RP) :: SUP, SDN

    real(RP) :: LNET, SNET, FLXUV
    real(RP) :: SW                  ! shortwave radition          [W/m/m]
    real(RP) :: LW                  ! longwave radition           [W/m/m]
    real(RP) :: psim,psim2,psim10   ! similality stability shear function for momentum
    real(RP) :: psih,psih2,psih10   ! similality stability shear function for heat

    ! for shadow effect model
    ! real(RP) :: HOUI1, HOUI2, HOUI3, HOUI4, HOUI5, HOUI6, HOUI7, HOUI8
    ! real(RP) :: SLX, SLX1, SLX2, SLX3, SLX4, SLX5, SLX6, SLX7, SLX8
    ! real(RP) :: THEATAZ    ! Solar Zenith Angle [rad]
    ! real(RP) :: THEATAS    ! = PI/2. - THETAZ
    ! real(RP) :: FAI        ! Latitude [rad]

    real(RP) :: SR, SB, SG, RR, RB, RG
    real(RP) :: SR1, SB1, SB2, SG1, SG2
    real(RP) :: RB1, RB2, RG1, RG2
    real(RP) :: HR, ELER, G0R
    real(RP) :: HB, ELEB, G0B
    real(RP) :: HG, ELEG, G0G

    real(RP) :: Z
    real(RP) :: QS0R, QS0B, QS0G

    real(RP) :: RIBR, BHR, CDR
    real(RP) :: RIBC, BHC, CDC
    real(RP) :: ALPHAR, ALPHAB, ALPHAG, ALPHAC
    real(RP) :: CHR, CHB, CHG, CHC
    real(RP) :: TC1, TC2, QC1, QC2
!    real(RP) :: CAPL1, AKSL1

    real(RP) :: DTS_MAX_onestep = 0.0_RP   ! DTS_MAX * dt
    real(RP) :: resi1,resi2     ! residual
    real(RP) :: resi1p,resi2p     ! residual
    real(RP) :: G0RP,G0BP,G0GP

    real(RP) :: XXX, X, CD, CH, CHU, XXX2, XXX10
    real(RP) :: LHV                              ! latent heat of vaporization [J/kg]
    real(RP) :: THA,THC,THS,THS1,THS2
    real(RP) :: RovCP

    integer  :: iteration

    !-----------------------------------------------------------
    ! Set parameters
    !-----------------------------------------------------------

    call HYDROMETEOR_LHV( LHV, TA )

    RovCP = Rdry / CPdry
    THA   = TA * ( PRE00 / PRSA )**RovCP

    !--- Renew surface and layer temperatures

    TRP = TR
    TBP = TB
    TGP = TG
    TCP = TC
    QCP = QC
    !
    TRLP = TRL
    TBLP = TBL
    TGLP = TGL
    !
    RAINRP = RAINR
    RAINBP = RAINB
    RAINGP = RAING

    !--- local time
    LAT = XLAT / D2R
    LON = XLON / D2R

    TIME = real( NOWDATE(4)*3600.0_RP + NOWDATE(5)*60.0_RP + NOWDATE(6), kind=RP )
    tloc = mod((NOWDATE(4) + int(LON/15.0_RP)),24 )
    dsec = real( NOWDATE(5)*60.0_RP + NOWDATE(6), kind=RP ) / 3600.0_RP
    if( tloc == 0 ) tloc = 24

    !--- Calculate AH data at LST
    if ( tloc == 24 ) then
      tahdiurnal = ( 1.0_RP-dsec ) * ahdiurnal(tloc  ) &
                 + (        dsec ) * ahdiurnal(1     )
    else
      tahdiurnal = ( 1.0_RP-dsec ) * ahdiurnal(tloc  ) &
                 + (        dsec ) * ahdiurnal(tloc+1)
    endif
    AH_t  = AH  * tahdiurnal
    ALH_t = ALH * tahdiurnal

    !--- limiter for surface temp change
    DTS_MAX_onestep = DTS_MAX * dt

    if ( ZDC + Z0C + 2.0_RP >= ZA ) then
       write(*,*) 'xxx [URBAN_PHY_SLC] ZDC + Z0C + 2m is larger than the 1st level! STOP.'
       call PRC_MPIstop
    endif

    W    = 2.0_RP * 1.0_RP * HGT
    VFGS = SVF
    VFGW = 1.0_RP - SVF
    VFWG = ( 1.0_RP - SVF ) * ( 1.0_RP - R ) / W
    VFWS = VFWG
    VFWW = 1.0_RP - 2.0_RP * VFWG

    SX   = SSG(1) + SSG(2) ! downward shortwave radiation [W/m2]
    RX   = LLG(1) + LLG(2) ! downward longwave  radiation [W/m2]

    SSGD = SSG(1)          ! downward direct  shortwave radiation [W/m2]
    SSGQ = SSG(2)          ! downward diffuse shortwave radiation [W/m2]

    !--- calculate canopy wind

    call canopy_wind(ZA, UA, UC)

    !-----------------------------------------------------------
    ! Set evaporation efficiency on roof/wall/road
    !-----------------------------------------------------------

    !!--- calculate evaporation efficiency
    RAINT = 1.0_RP * ( RAIN * dt )            ! [kg/m2/s -> kg/m2]
    call cal_beta(BETR, RAINT, RAINR, STRGR, ROFFR)

    RAINT = 0.1_RP * ( RAIN * dt )
    call cal_beta(BETB, RAINT, RAINB, STRGB, ROFFB)

    RAINT = 0.9_RP * ( RAIN * dt )
    call cal_beta(BETG, RAINT, RAING, STRGG, ROFFG)

    ROFF = ROFF +  R * ROFFR  + RW * ( ROFFB + ROFFG )

    !-----------------------------------------------------------
    ! Radiation : Net Short Wave Radiation at roof/wall/road
    !-----------------------------------------------------------

    if( SX > 0.0_RP ) then !  SSG is downward short

      ! currently we use no shadow effect model
      !!     IF(.NOT.SHADOW) THEN              ! no shadow effects model

      SR1  = SX * ( 1.0_RP - ALBR )
      SG1  = SX * VFGS * ( 1.0_RP - ALBG )
      SB1  = SX * VFWS * ( 1.0_RP - ALBB )
      SG2  = SB1 * ALBB / ( 1.0_RP - ALBB ) * VFGW * ( 1.0_RP - ALBG )
      SB2  = SG1 * ALBG / ( 1.0_RP - ALBG ) * VFWG * ( 1.0_RP - ALBB )

      SR   = SR1
      SG   = SG1 + SG2
      SB   = SB1 + SB2
      SNET = R * SR + W * SB + RW * SG

    else

      SR   = 0.0_RP
      SG   = 0.0_RP
      SB   = 0.0_RP
      SNET = 0.0_RP

    end if

    !-----------------------------------------------------------
    ! Energy balance on roof/wall/road surface
    !-----------------------------------------------------------

    !--------------------------------------------------
    !   Roof
    !--------------------------------------------------

    ! new scheme

     G0RP = 0.0_RP
     XXXR = 0.0_RP
     resi1p = 0.0_RP
     do iteration = 1, 100

      THS   = TR * ( PRE00 / PRSS )**RovCP   ! potential temp

      Z    = ZA - ZDC
      BHR  = LOG(Z0R/Z0HR) / 0.4_RP
      RIBR = ( GRAV * 2.0_RP / (THA+THS) ) * (THA-THS) * (Z+Z0R) / (UA*UA)
      call mos(XXXR,CHR,CDR,BHR,RIBR,Z,Z0R,UA,THA,THS,RHOO)

      call qsat( QS0R, TR, PRSS )

      RR    = EPSR * ( RX - STB * (TR**4)  )
      !HR    = RHOO * CPdry * CHR * UA * (TR-TA)
      HR    = RHOO * CPdry * CHR * UA * (THS-THA)
      ELER  = RHOO * LHV   * CHR * UA * BETR * (QS0R-QA)
      G0R   = SR + RR - HR - ELER

    !--- calculate temperature in roof
    !  if ( STRGR /= 0.0_RP ) then
    !    CAPL1 = CAP_water * (RAINR / (DZR(1) + RAINR)) + CAPR * (DZR(1) / (DZR(1) + RAINR))
    !    AKSL1 = AKS_water * (RAINR / (DZR(1) + RAINR)) + AKSR * (DZR(1) / (DZR(1) + RAINR))
    !  else
    !    CAPL1 = CAPR
    !    AKSL1 = AKSR
    !  endif
    !! 1st layer's cap, aks are replaced.
    !! call multi_layer2(UKE,BOUND,G0R,CAPR,AKSR,TRL,DZR,dt,TRLEND,CAPL1,AKSL1)

      TRL = TRLP
      call multi_layer(UKE,BOUND,G0R,CAPR,AKSR,TRL,DZR,dt,TRLEND)
      resi1 = TRL(1) - TR

     ! write(*,'(a3,i5,f8.3,6f15.5)') "TR,",iteration,TR,G0R,SR,RR,HR,ELER,resi1

      if( abs(resi1) < sqrt(EPS) ) then
        TR = TRL(1)
        TR = max( TRP - DTS_MAX_onestep, min( TRP + DTS_MAX_onestep, TR ) )
        exit
      endif

      if ( resi1*resi1p < 0.0_RP ) then
        TR = (TR + TRL(1)) * 0.5_RP
      else
        TR = TRL(1)
      endif
      TR = max( TRP - DTS_MAX_onestep, min( TRP + DTS_MAX_onestep, TR ) )

      resi1p = resi1

     enddo

!    if( .NOT. (resi1 < sqrt(EPS)) ) then
!       write(*,*) 'xxx Warning not converged for TR in URBAN SLC', &
!            PRC_myrank, i,j, &
!            resi1, G0R
!    end if

     ! output for debug
     if ( iteration > 100 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Warning: [URBAN_PHY_SLC/SLC_main] iteration for TR was not converged',PRC_myrank,i,j
       if( IO_L ) write(IO_FID_LOG,*) '---------------------------------------------------------------------------------'
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Residual                                          [K] :', resi1
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TRP : Initial TR                                  [K] :', TRP
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TRLP: Initial TRL                                 [K] :', TRLP
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- SX  : Shortwave radiation                      [W/m2] :', SX
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- RX  : Longwave radiation                       [W/m2] :', RX
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- PRSS: Surface pressure                           [Pa] :', PRSS
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- PRSA: Pressure at 1st atmos layer                 [m] :', PRSA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- RHOO: Air density                             [kg/m3] :', RHOO
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- ZA  : Height at 1st atmos layer                   [m] :', ZA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TA  : Temperature at 1st atmos layer              [K] :', TA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- UA  : Wind speed at 1st atmos layer             [m/s] :', UA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- QA  : Specific humidity at 1st atmos layer    [kg/kg] :', QA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- DZR : Depth of surface layer                      [m] :', DZR
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- R, W, RW : Normalized height and road width       [-] :', R, W,RW
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- SVF : Sky View Factors                            [-] :', SVF
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- BETR: Evaporation efficiency                      [-] :', BETR
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- EPSR: Surface emissivity of roof                  [-] :', EPSR
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- CAPR: Heat capacity of roof                 [J m-3 K] :', CAPR
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- AKSR: Thermal conductivity of roof          [W m-1 K] :', AKSR
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- QS0R: Surface specific humidity               [kg/kg] :', QS0R
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- ZDC : Desplacement height of canopy               [m] :', ZDC
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Z0R : Momentum roughness length of roof           [m] :', Z0R
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Z0HR: Thermal roughness length of roof            [m] :', Z0HR
       if( IO_L ) write(IO_FID_LOG,*) '---------------------------------------------------------------------------------'
     endif

    !--- update only fluxes ----
     THS   = TR * ( PRE00 / PRSS )**RovCP
     RIBR = ( GRAV * 2.0_RP / (THA+THS) ) * (THA-THS) * (Z+Z0R) / (UA*UA)
     call mos(XXXR,CHR,CDR,BHR,RIBR,Z,Z0R,UA,THA,THS,RHOO)

     call qsat( QS0R, TR, PRSS )

     RR      = EPSR * ( RX - STB * (TR**4) )
     HR      = RHOO * CPdry * CHR * UA * (THS-THA)
     ELER    = RHOO * LHV   * CHR * UA * BETR * (QS0R-QA)
     G0R     = SR + RR - HR - ELER

     TRL   = TRLP
     call multi_layer(UKE,BOUND,G0R,CAPR,AKSR,TRL,DZR,dt,TRLEND)
     resi1 = TRL(1) - TR
     TR    = TRL(1)

     if ( abs(resi1) > DTS_MAX_onestep ) then
       if ( abs(resi1) > DTS_MAX_onestep*10.0_RP ) then
         write(*,*) 'xxx [URBAN_PHY_SLC/SLC_main] tendency of TR exceeded a limit! STOP.'
         write(*,*) 'xxx previous TR and updated TR(TRL(1)) is ',TR-resi1, TR
         call PRC_MPIstop
       endif
       if( IO_L ) write(IO_FID_LOG,*) '*** [URBAN_PHY_SLC/SLC_main] tendency of TR exceeded a limit'
       if( IO_L ) write(IO_FID_LOG,*) '*** previous TR and updated TR(TRL(1)) is ', TR-resi1, TR
     endif

    !--------------------------------------------------
    !   Wall and Road
    !--------------------------------------------------

    ! new scheme

    ! empirical form
      ALPHAB = 6.15_RP + 4.18_RP * UC
      if( UC > 5.0_RP ) ALPHAB = 7.51_RP * (UC**0.78_RP )
      ALPHAG = 6.15_RP + 4.18_RP * UC
      if( UC > 5.0_RP ) ALPHAG = 7.51_RP * (UC**0.78_RP )
      CHB = ALPHAB / RHOO / CPdry / UC
      CHG = ALPHAG / RHOO / CPdry / UC

     G0BP = 0.0_RP
     G0GP = 0.0_RP
     XXXC = 0.0_RP
     resi1p = 0.0_RP
     resi2p = 0.0_RP

     do iteration = 1, 200

      THS1   = TB * ( PRE00 / PRSS )**RovCP
      THS2   = TG * ( PRE00 / PRSS )**RovCP
      THC    = TC * ( PRE00 / PRSS )**RovCP

      Z    = ZA - ZDC
      BHC  = LOG(Z0C/Z0HC) / 0.4_RP
      RIBC = ( GRAV * 2.0_RP / (THA+THC) ) * (THA-THC) * (Z+Z0C) / (UA*UA)
      call mos(XXXC,CHC,CDC,BHC,RIBC,Z,Z0C,UA,THA,THC,RHOO)
      ALPHAC = CHC * RHOO * CPdry * UA

      call qsat( QS0B, TB, PRSS )
      call qsat( QS0G, TG, PRSS )

      TC1   = RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
      !TC2   = RW*ALPHAC*TA + RW*ALPHAG*TG + W*ALPHAB*TB
      TC2   = RW*ALPHAC*THA + W*ALPHAB*THS1 + RW*ALPHAG*THS2
      THC   = TC2 / TC1
      QC1   = RW*(CHC*UA)    + RW*(CHG*BETG*UC)      + W*(CHB*BETB*UC)
      QC2   = RW*(CHC*UA)*QA + RW*(CHG*BETG*UC)*QS0G + W*(CHB*BETB*UC)*QS0B
      QC    = QC2 / QC1

      RG1   = EPSG * ( RX * VFGS                  &
                     + EPSB * VFGW * STB * TB**4  &
                     - STB * TG**4                )

      RB1   = EPSB * ( RX * VFWS                  &
                     + EPSG * VFWG * STB * TG**4  &
                     + EPSB * VFWW * STB * TB**4  &
                     - STB * TB**4                )

      RG2   = EPSG * ( (1.0_RP-EPSB) * VFGW * VFWS * RX                   &
                     + (1.0_RP-EPSB) * VFGW * VFWG * EPSG * STB * TG**4  &
                     + EPSB * (1.0_RP-EPSB) * VFGW * VFWW * STB * TB**4  )

      RB2   = EPSB * ( (1.0_RP-EPSG) * VFWG * VFGS * RX                                  &
                     + (1.0_RP-EPSG) * EPSB * VFGW * VFWG * STB * TB**4                  &
                     + (1.0_RP-EPSB) * VFWS * VFWW * RX                                  &
                     + (1.0_RP-EPSB) * VFWG * VFWW * STB * EPSG * TG**4                  &
                     + EPSB * (1.0_RP-EPSB) * VFWW * (1.0_RP-2.0_RP*VFWS) * STB * TB**4  )

      RG    = RG1 + RG2
      RB    = RB1 + RB2

      HB    = RHOO * CPdry * CHB * UC * (THS1-THC)
      ELEB  = RHOO * LHV   * CHB * UC * BETB * (QS0B-QC)
      G0B   = SB + RB - HB - ELEB

      HG    = RHOO * CPdry * CHG * UC * (THS2-THC)
      ELEG  = RHOO * LHV   * CHG * UC * BETG * (QS0G-QC)
      G0G   = SG + RG - HG - ELEG

      TBL = TBLP
      call multi_layer(UKE,BOUND,G0B,CAPB,AKSB,TBL,DZB,dt,TBLEND)
      resi1 = TBL(1) - TB

      TGL = TGLP
      call multi_layer(UKE,BOUND,G0G,CAPG,AKSG,TGL,DZG,dt,TGLEND)
      resi2 = TGL(1) - TG

      !-----------
      !print *,HB, RHOO , CPdry , CHB , UC , THS1,THC
      !print *,HG, RHOO , CPdry , CHG , UC , THS2,THC
      !print *,ELEB ,RHOO , LHV , CHB , UC , BETB , QS0B , QC
      !print *,ELEG ,RHOO , LHV , CHG , UC , BETG , QS0G , QC

      !write(*,'(a3,i5,f8.3,6f15.5)') "TB,",iteration,TB,G0B,SB,RB,HB,ELEB,resi1
      !write(*,'(a3,i5,f8.3,6f15.5)') "TG,",iteration,TG,G0G,SG,RG,HG,ELEG,resi2
      !write(*,'(a3,i5,f8.3,3f15.5)') "TC,",iteration,TC,QC,QS0B,QS0G
      !--------
      !resi1  =  abs(G0B - G0BP)
      !resi2  =  abs(G0G - G0GP)
      !G0BP   = G0B
      !G0GP   = G0G

      if ( abs(resi1) < sqrt(EPS) .AND. abs(resi2) < sqrt(EPS) ) then
         TB = TBL(1)
         TG = TGL(1)
         TB = max( TBP - DTS_MAX_onestep, min( TBP + DTS_MAX_onestep, TB ) )
         TG = max( TGP - DTS_MAX_onestep, min( TGP + DTS_MAX_onestep, TG ) )
         exit
      endif

      if ( resi1*resi1p < 0.0_RP ) then
         TB = (TB + TBL(1)) * 0.5_RP
      else
         TB = TBL(1)
      endif
      if ( resi2*resi2p < 0.0_RP ) then
         TG = (TG + TGL(1)) * 0.5_RP
      else
         TG = TGL(1)
      endif
      TB = max( TBP - DTS_MAX_onestep, min( TBP + DTS_MAX_onestep, TB ) )
      TG = max( TGP - DTS_MAX_onestep, min( TGP + DTS_MAX_onestep, TG ) )

      resi1p = resi1
      resi2p = resi2

      ! this is for TC, QC
      call qsat( QS0B, TB, PRSS )
      call qsat( QS0G, TG, PRSS )

      THS1   = TB * ( PRE00 / PRSS )**RovCP
      THS2   = TG * ( PRE00 / PRSS )**RovCP

      TC1    =  RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
      !TC2    =  RW*ALPHAC*THA + RW*ALPHAG*TG + W*ALPHAB*TB
      TC2    =  RW*ALPHAC*THA + W*ALPHAB*THS1 + RW*ALPHAG*THS2
      THC    =  TC2 / TC1
      TC = THC * ( PRSS / PRE00 )**RovCP

      QC1    =  RW*(CHC*UA)    + RW*(CHG*BETG*UC)      + W*(CHB*BETB*UC)
      QC2    =  RW*(CHC*UA)*QA + RW*(CHG*BETG*UC)*QS0G + W*(CHB*BETB*UC)*QS0B
      QC     =  QC2 / QC1

    enddo

!    if( .NOT. (resi1 < sqrt(EPS) .AND. resi2 < sqrt(EPS) ) ) then
!       write(*,*) 'xxx Warning not converged for TG, TB in URBAN SLC', &
!            PRC_myrank, i,j, &
!            resi1, resi2, TB, TG, TC, G0BP, G0GP, RB, HB, RG, HG, QC, &
!            CHC, CDC, BHC, ALPHAC
!       write(*,*) TBP, TGP, TCP, &
!                  PRSS, THA, UA, QA, RHOO, UC, QCP, &
!                  ZA, ZDC, Z0C, Z0HC, &
!                  CHG, CHB, &
!                  W, RW, ALPHAG, ALPHAB, BETG, BETB, &
!                  RX, VFGS, VFGW, VFWG, VFWW, VFWS, STB, &
!                  SB, SG, LHV, TBLP, TGLP
!       write(*,*) "6",VFGS, VFGW, VFWG, VFWW, VFWS
!    end if

     ! output for debug
     if ( iteration > 200 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Warning: [URBAN_PHY_SLC/SLC_main] iteration for TB/TG was not converged',PRC_myrank,i,j
       if( IO_L ) write(IO_FID_LOG,*) '---------------------------------------------------------------------------------'
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Residual                                       [K] :', resi1,resi2
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TBP : Initial TB                               [K] :', TBP
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TBLP: Initial TBL                              [K] :', TBLP
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TGP : Initial TG                               [K] :', TGP
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TGLP: Initial TGL                              [K] :', TGLP
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TCP : Initial TC                               [K] :', TCP
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- QCP : Initial QC                               [K] :', QCP
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- UC  : Canopy wind                            [m/s] :', UC
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- SX  : Shortwave radiation                   [W/m2] :', SX
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- RX  : Longwave radiation                    [W/m2] :', RX
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- PRSS: Surface pressure                        [Pa] :', PRSS
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- PRSA: Pressure at 1st atmos layer              [m] :', PRSA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- RHOO: Air density                          [kg/m3] :', RHOO
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- ZA  : Height at 1st atmos layer                [m] :', ZA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- TA  : Temperature at 1st atmos layer           [K] :', TA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- UA  : Wind speed at 1st atmos layer          [m/s] :', UA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- QA  : Specific humidity at 1st atmos layer [kg/kg] :', QA
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- DZB : Depth of surface layer                   [m] :', DZB
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- DZG : Depth of surface layer                   [m] :', DZG
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- R, W, RW  : Normalized height and road width    [-] :', R, W,RW
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- SVF       : Sky View Factors                    [-] :', SVF
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- BETB,BETG : Evaporation efficiency              [-] :', BETB,BETG
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- EPSB,EPSG : Surface emissivity                  [-] :', EPSB,EPSG
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- CAPB,CAPG : Heat capacity                 [J m-3 K] :', CAPB,CAPG
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- AKSB,AKSG : Thermal conductivity          [W m-1 K] :', AKSB,AKSB
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- QS0B,QS0G : Surface specific humidity       [kg/kg] :', QS0B,QS0G
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- ZDC       : Desplacement height of canopy       [m] :', ZDC
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Z0C       : Momentum roughness length of canopy [m] :', Z0C
       if( IO_L ) write(IO_FID_LOG,*) 'DEBUG Message --- Z0HC      : Thermal roughness length of canopy  [m] :', Z0HC
       if( IO_L ) write(IO_FID_LOG,*) '---------------------------------------------------------------------------------'
     endif


    !--- update only fluxes ----
     RG1      = EPSG * ( RX * VFGS                  &
                       + EPSB * VFGW * STB * TB**4  &
                       - STB * TG**4                )
     RB1      = EPSB * ( RX * VFWS                  &
                       + EPSG * VFWG * STB * TG**4  &
                       + EPSB * VFWW * STB * TB**4  &
                       - STB * TB**4                )

     RG2      = EPSG * ( (1.0_RP-EPSB) * VFGW * VFWS * RX                  &
                       + (1.0_RP-EPSB) * VFGW * VFWG * EPSG * STB * TG**4  &
                       + EPSB * (1.0_RP-EPSB) * VFGW * VFWW * STB * TB**4  )
     RB2      = EPSB * ( (1.0_RP-EPSG) * VFWG * VFGS * RX                                 &
                       + (1.0_RP-EPSG) * EPSB * VFGW * VFWG * STB * TB**4                 &
                       + (1.0_RP-EPSB) * VFWS * VFWW * RX                                 &
                       + (1.0_RP-EPSB) * VFWG * VFWW * STB * EPSG * TG**4                 &
                       + EPSB * (1.0_RP-EPSB) * VFWW * (1.0_RP-2.0_RP*VFWS) * STB * TB**4 )

     RG       = RG1 + RG2
     RB       = RB1 + RB2

     THS1   = TB * ( PRE00 / PRSS )**RovCP
     THS2   = TG * ( PRE00 / PRSS )**RovCP
     THC    = TC * ( PRE00 / PRSS )**RovCP

     HB   = RHOO * CPdry * CHB * UC * (THS1-THC)
     ELEB = RHOO * LHV   * CHB * UC * BETB * (QS0B-QC)
     G0B  = SB + RB - HB - ELEB

     HG   = RHOO * CPdry * CHG * UC * (THS2-THC)
     ELEG = RHOO * LHV   * CHG * UC * BETG * (QS0G-QC)
     G0G  = SG + RG - HG - ELEG

     TBL   = TBLP
     call multi_layer(UKE,BOUND,G0B,CAPB,AKSB,TBL,DZB,dt,TBLEND)
     resi1 = TBL(1) - TB
     TB    = TBL(1)

     TGL   = TGLP
     call multi_layer(UKE,BOUND,G0G,CAPG,AKSG,TGL,DZG,dt,TGLEND)
     resi2 = TGL(1) - TG
     TG    = TGL(1)

     if ( abs(resi1) > DTS_MAX_onestep ) then
        if ( abs(resi1) > DTS_MAX_onestep*10.0_RP ) then
           write(*,*) 'xxx [URBAN_PHY_SLC/SLC_main] tendency of TB exceeded a limit! STOP.'
           write(*,*) 'xxx previous TB and updated TB(TBL(1)) is ', TB-resi1,TB
           call PRC_MPIstop
        endif
        if( IO_L ) write(IO_FID_LOG,*) '*** [URBAN_PHY_SLC/SLC_main] tendency of TB exceeded a limit'
        if( IO_L ) write(IO_FID_LOG,*) '*** previous TB and updated TB(TBL(1)) is ', TB-resi1, TB
     endif

     if ( abs(resi2) > DTS_MAX_onestep ) then
        if ( abs(resi2) > DTS_MAX_onestep*10.0_RP ) then
           write(*,*) 'xxx [URBAN_PHY_SLC/SLC_main] tendency of TG exceeded a limit! STOP.'
           write(*,*) 'xxx previous TG and updated TG(TGL(1)) is ', TG-resi2, TG, resi2
           call PRC_MPIstop
        endif
        if( IO_L ) write(IO_FID_LOG,*) '*** [URBAN_PHY_SLC/SLC_main] tendency of TG exceeded a limit'
        if( IO_L ) write(IO_FID_LOG,*) '*** previous TG and updated TG(TGL(1)) is ', TG-resi2, TG
     endif

    !-----------------------------------------------------------
    ! Total Fluxes from Urban Canopy
    !-----------------------------------------------------------

    FLXUV = ( R*CDR + RW*CDC ) * UA * UA
    SH    = ( R*HR   + W*HB   + RW*HG )              ! Sensible heat flux   [W/m/m]
    LH    = ( R*ELER + W*ELEB + RW*ELEG )            ! Latent heat flux     [W/m/m]
    GHFLX = -1.0_RP * ( R*G0R + W*G0B + RW*G0G )
    LNET  = R*RR + W*RB + RW*RG

    !-----------------------------------------------------------
    ! Grid average
    !-----------------------------------------------------------

    LW = RX - LNET     ! Upward longwave radiation   [W/m/m]
    SW = SX - SNET     ! Upward shortwave radiation  [W/m/m]
    RN = (SNET+LNET)    ! Net radiation [W/m/m]

    !--- shortwave radiation
    SDN =  R + W * (VFWS + VFGS * ALBG * VFWG)  + RW * (VFGS + VFWS * ALBB * VFGW)
    SUP =  R * ALBR                                        &
         + W *  ( VFWS * ALBB + VFGS * ALBG * VFWG *ALBB ) &
         + RW * ( VFGS * ALBG + VFWS * ALBB * VFGW *ALBG )

    ALBD_SW_grid = SUP / SDN

    !--- longwave radiation
    LDN =  R + W*VFWS + RW*VFGS
    LUP =  R * (1.0_RP-EPSR)                                                           &
         + W*( (1.0_RP-EPSB*VFWW)*(1.0_RP-EPSB)*VFWS - EPSB*VFWG*(1.0_RP-EPSG)*VFGS )  &
         + RW*( (1.0_RP-EPSG)*VFGS - EPSG*(1.0_RP-VFGS)*(1.0_RP-EPSB)*VFWS )

    RUP = (LDN - LUP) * RX - LNET
    ALBD_LW_grid = LUP / LDN


    ! RUP =  R * (EPSR * STB * TR**4 ) &
    !      + W * (EPSB*STB*(TB**4) - EPSB*EPSG*VFWG*STB*(TG**4) - EPSB*EPSB*VFWW*STB*(TB**4)  &
    !           - EPSB *(1.0_RP-EPSG) * EPSB * VFGW * VFWG * STB * (TB**4)         &
    !           - EPSB *(1.0_RP-EPSB) * VFWG * VFWW * STB * EPSG * (TG**4)         &
    !           - EPSB * EPSB * (1.0_RP-EPSB) * VFWW * VFWW * STB * (TB**4) )      &
    !      + RW * (EPSG*STB*(TG**4) - EPSG * EPSB * VFGW * STB * (TB**4)             &
    !           - EPSG * EPSB * (1.0_RP-EPSB) * (1.0_RP-SVF) * VFWW * STB * TB**4  )


    SHR = HR             ! Sensible heat flux on roof [W/m/m]
    SHB = HB             ! Sensible heat flux on wall [W/m/m]
    SHG = HG             ! Sensible heat flux on road [W/m/m]
    LHR = ELER           ! Latent heat flux on road [W/m/m]
    LHB = ELEB           ! Latent heat flux on wall [W/m/m]
    LHG = ELEG           ! Latent heat flux on road [W/m/m]
    GHR = -1.0_RP * G0R  ! Ground heat flux on roof [W/m/m]
    GHB = -1.0_RP * G0B  ! Ground heat flux on wall [W/m/m]
    GHG = -1.0_RP * G0G  ! Ground heat flux on road [W/m/m]
    RNR = SR + RR        ! Net radiation on roof [W/m/m]
    RNB = SB + RB        ! Net radiation on building [W/m/m]
    RNG = SG + RG        ! Net radiation on ground [W/m/m]

    ! calculate rain amount remaining on the surface
    RAINR = max(0.0_RP, RAINR-(LHR/LHV)*real(dt,kind=RP)) ! [kg/m/m = mm]
    RAINB = max(0.0_RP, RAINB-(LHB/LHV)*real(dt,kind=RP)) ! [kg/m/m = mm]
    RAING = max(0.0_RP, RAING-(LHG/LHV)*real(dt,kind=RP)) ! [kg/m/m = mm]

    !-----------------------------------------------------------
    !  diagnostic GRID AVERAGED TS from upward logwave
    !-----------------------------------------------------------

    RTS = ( RUP / STB / ( 1.0_RP-ALBD_LW_grid) )**0.25

    !-----------------------------------------------------------
    !  diagnostic grid average U10, V10, T2, Q2 from urban
    !  Below method would be better to be improved. This is tentative method.
    !-----------------------------------------------------------
    !UST = sqrt( FLXUV )             ! u* [m/s]
    !TST = -SH / RHOO / CPdry / UST  ! T* [K]
    !QST = -LH / RHOO / LHV   / UST    ! q* [-]
    !Z = ZA - ZDC
    !XXX = 0.4*9.81*Z*TST/TA/UST/UST

    XXX = XXXC
    call cal_psi(XXX,psim,psih)

    XXX2 = (2.0_RP/Z) * XXX
    call cal_psi(XXX2,psim2,psih2)

    XXX10 = (10.0_RP/Z) * XXX
    call cal_psi(XXX10,psim10,psih10)

    !U10 = U1 * ((log(10.0_RP/Z0C)-psim10)/(log(Z/Z0C)-psim))  ! u at 10 m [m/s]
    !V10 = V1 * ((log(10.0_RP/Z0C)-psim10)/(log(Z/Z0C)-psim))  ! v at 10 m [m/s]
    U10 = U1 * log(10.0_RP/Z0C) / log(Z/Z0C)
    V10 = V1 * log(10.0_RP/Z0C) / log(Z/Z0C)

    T2  = RTS + (TA-RTS)*((log(2.0_RP/Z0HC)-psih2)/(log(Z/Z0HC)-psih))
    Q2 = QC

    !-----------------------------------------------------------
    ! add anthropogenic heat fluxes
    !-----------------------------------------------------------

    SH     = SH + AH_t           ! Sensible heat flux          [W/m/m]
    LH     = LH + ALH_t          ! Latent heat flux            [W/m/m]

    return
  end subroutine SLC_main

  !-----------------------------------------------------------------------------
  subroutine canopy_wind(ZA, UA, UC)
    implicit none

    real(RP), intent(in)  :: ZA   ! height at 1st atmospheric level [m]
    real(RP), intent(in)  :: UA   ! wind speed at 1st atmospheric level [m/s]
    real(RP), intent(out) :: UC   ! wind speed at 1st atmospheric level [m/s]

    real(RP) :: UR,ZC,XLB,BB

    if( ZR + 2.0_RP < ZA ) then
      UR  = UA * log((ZR-ZDC)/Z0C) / log((ZA-ZDC)/Z0C)
      ZC  = 0.7_RP * ZR
      XLB = 0.4_RP * (ZR-ZDC)
      ! BB formulation from Inoue (1963)
      BB  = 0.4_RP * ZR / ( XLB * log((ZR-ZDC)/Z0C) )
      UC  = UR * exp( -BB * (1.0_RP-ZC/ZR) )
    else
      ! PRINT *, 'Warning ZR + 2m  is larger than the 1st WRF level'
      ZC  = ZA / 2.0_RP
      UC  = UA / 2.0_RP
    endif

    UC = max(UC,0.01_RP)

    return
  end subroutine canopy_wind

  !-----------------------------------------------------------------------------
  subroutine cal_beta(BET, RAIN, WATER, STRG, ROFF)
    implicit none

    real(RP), intent(out)   :: BET    ! evapolation efficiency [-]
    real(RP), intent(in)    :: RAIN   ! precipitation [mm*dt]
    real(RP), intent(inout) :: WATER  ! rain amount in strage [kg/m2]
    real(RP), intent(inout) :: STRG   ! rain strage [kg/m2]
    real(RP), intent(inout) :: ROFF   ! runoff [kg/m2]

    if ( STRG == 0.0_RP ) then ! not consider evapolation from urban
       BET   = 0.0_RP
       ROFF  = RAIN
    else
       WATER = WATER + RAIN
       ROFF  = max(0.0_RP, WATER-STRG)
       WATER = WATER - max(0.0_RP, WATER-STRG)
       BET   = min ( WATER / STRG, 1.0_RP)
    endif

    return
  end subroutine cal_beta

  !-----------------------------------------------------------------------------
  subroutine cal_psi(zeta,psim,psih)
    use scale_const, only: &
       PI     => CONST_PI
    implicit none

    real(RP), intent(inout) :: zeta  ! z/L
    real(RP), intent(out)   :: psim
    real(RP), intent(out)   :: psih
    real(RP)                :: X

    if( zeta >=  1.0_RP ) zeta =  1.0_RP
    if( zeta <= -5.0_RP ) zeta = -5.0_RP

    if( zeta > 0.0_RP ) then
      psim = -5.0_RP * zeta
      psih = -5.0_RP * zeta
    else
      X    = ( 1.0_RP - 16.0_RP * zeta )**0.25_RP
      psim = 2.0_RP * log((1.0_RP+X)/2.0_RP) + log((1.0_RP+X*X)/2.0_RP) - 2.0_RP*atan(X) + PI/2.0_RP
      psih = 2.0_RP * log((1.0_RP+X*X)/2.0_RP)
    end if

    return
  end subroutine cal_psi

  !-----------------------------------------------------------------------------
  !  XXX:   z/L (requires iteration by Newton-Rapson method)
  !  B1:    Stanton number
  !  PSIM:  = PSIX of LSM
  !  PSIH:  = PSIT of LSM
  subroutine mos(XXX,CH,CD,B1,RIB,Z,Z0,UA,TA,TSF,RHO)
    use scale_const, only: &
       CPdry => CONST_CPdry ! CPP : heat capacity of dry air [J/K/kg]
    implicit none

    real(RP), intent(in)    :: B1, Z, Z0, UA, TA, TSF, RHO
    real(RP), intent(out)   :: CD, CH
    real(RP), intent(inout) :: XXX, RIB
    real(RP)                :: XXX0, X, X0, FAIH, DPSIM, DPSIH
    real(RP)                :: F, DF, XXXP, US, TS, AL, XKB, DD, PSIM, PSIH
    integer                 :: NEWT
    integer, parameter      :: NEWT_END = 10

    real(RP)                :: lnZ

    lnZ = log( (Z+Z0)/Z0 )

    if( RIB <= -15.0_RP ) RIB = -15.0_RP

    if( RIB < 0.0_RP ) then

      do NEWT = 1, NEWT_END

        if( XXX >= 0.0_RP ) XXX = -1.0e-3_RP

        XXX0  = XXX * Z0/(Z+Z0)

        X     = (1.0_RP-16.0_RP*XXX)**0.25
        X0    = (1.0_RP-16.0_RP*XXX0)**0.25

        PSIM  = lnZ &
              - log( (X+1.0_RP)**2 * (X**2+1.0_RP) ) &
              + 2.0_RP * atan(X) &
              + log( (X+1.0_RP)**2 * (X0**2+1.0_RP) ) &
              - 2.0_RP * atan(X0)
        FAIH  = 1.0_RP / sqrt( 1.0_RP - 16.0_RP*XXX )
        PSIH  = lnZ + 0.4_RP*B1 &
              - 2.0_RP * log( sqrt( 1.0_RP - 16.0_RP*XXX ) + 1.0_RP ) &
              + 2.0_RP * log( sqrt( 1.0_RP - 16.0_RP*XXX0 ) + 1.0_RP )

        DPSIM = ( 1.0_RP - 16.0_RP*XXX )**(-0.25) / XXX &
              - ( 1.0_RP - 16.0_RP*XXX0 )**(-0.25) / XXX
        DPSIH = 1.0_RP / sqrt( 1.0_RP - 16.0_RP*XXX ) / XXX &
              - 1.0_RP / sqrt( 1.0_RP - 16.0_RP*XXX0 ) / XXX

        F     = RIB * PSIM**2 / PSIH - XXX

        DF    = RIB * ( 2.0_RP*DPSIM*PSIM*PSIH - DPSIH*PSIM**2 ) &
              / PSIH**2 - 1.0_RP

        XXXP  = XXX
        XXX   = XXXP - F / DF

        if( XXX <= -10.0_RP ) XXX = -10.0_RP

      end do

    else if( RIB >= 0.142857_RP ) then

      XXX  = 0.714_RP
      PSIM = lnZ + 7.0_RP * XXX
      PSIH = PSIM + 0.4_RP * B1

    else

      AL   = lnZ
      XKB  = 0.4_RP * B1
      DD   = -4.0_RP * RIB * 7.0_RP * XKB * AL + (AL+XKB)**2
      if( DD <= 0.0_RP ) DD = 0.0_RP

      XXX  = ( AL + XKB - 2.0_RP*RIB*7.0_RP*AL - sqrt(DD) ) / ( 2.0_RP * ( RIB*7.0_RP**2 - 7.0_RP ) )
      PSIM = lnZ + 7.0_RP * min( XXX, 0.714_RP )
      PSIH = PSIM + 0.4_RP * B1

    endif

    US = 0.4_RP * UA / PSIM             ! u*
    if( US <= 0.01_RP ) US = 0.01_RP
    TS = 0.4_RP * (TA-TSF) / PSIH       ! T*

    CD    = US * US / UA**2             ! CD
    CH    = 0.4_RP * US / PSIH / UA     ! CH
    !ALPHA = RHO * CPdry * 0.4_RP * US / PSIH  ! RHO*CP*CH*U

    return
  end subroutine mos

  !-------------------------------------------------------------------
  subroutine multi_layer(KM,BOUND,G0,CAP,AKS,TSL,DZ,DELT,TSLEND)
  !
  !  calculate temperature in roof/building/road
  !  multi-layer heat equation model
  !  Solving Heat Equation by Tri Diagonal Matrix Algorithm
  !-------------------------------------------------------------------

    implicit none

    real(RP), intent(in)    :: G0
    real(RP), intent(in)    :: CAP
    real(RP), intent(in)    :: AKS
    real(DP), intent(in)    :: DELT      ! Time step [ s ]
    real(RP), intent(in)    :: TSLEND
    integer,  intent(in)    :: KM
    integer,  intent(in)    :: BOUND
    real(RP), intent(in)    :: DZ(KM)
    real(RP), intent(inout) :: TSL(KM)
    real(RP)                :: A(KM), B(KM), C(KM), D(KM), P(KM), Q(KM)
    real(RP)                :: DZEND
    integer                 :: K

    DZEND = DZ(KM)

    A(1)  = 0.0_RP

    B(1)  = CAP * DZ(1) / DELT &
          + 2.0_RP * AKS / (DZ(1)+DZ(2))
    C(1)  = -2.0_RP * AKS / (DZ(1)+DZ(2))
    D(1)  = CAP * DZ(1) / DELT * TSL(1) + G0

    do K = 2, KM-1
      A(K) = -2.0_RP * AKS / (DZ(K-1)+DZ(K))
      B(K) = CAP * DZ(K) / DELT + 2.0_RP * AKS / (DZ(K-1)+DZ(K)) + 2.0_RP * AKS / (DZ(K)+DZ(K+1))
      C(K) = -2.0_RP * AKS / (DZ(K)+DZ(K+1))
      D(K) = CAP * DZ(K) / DELT * TSL(K)
    end do

    if( BOUND == 1 ) then ! Flux=0
      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      C(KM) = 0.0_RP
      D(KM) = CAP * DZ(KM) / DELT * TSL(KM)
    else ! T=constant
      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM)) + 2.0_RP * AKS / (DZ(KM)+DZEND)
      C(KM) = 0.0_RP
      D(KM) = CAP * DZ(KM) / DELT * TSL(KM) + 2.0_RP * AKS * TSLEND / (DZ(KM)+DZEND)
    end if

    P(1) = -C(1) / B(1)
    Q(1) =  D(1) / B(1)

    do K = 2, KM
      P(K) = -C(K) / ( A(K) * P(K-1) + B(K) )
      Q(K) = ( -A(K) * Q(K-1) + D(K) ) / ( A(K) * P(K-1) + B(K) )
    end do

    TSL(KM) = Q(KM)

    do K = KM-1, 1, -1
      TSL(K) = P(K) * TSL(K+1) + Q(K)
    end do

    return
  end subroutine multi_layer

  !-------------------------------------------------------------------
  subroutine multi_layer2(KM,BOUND,G0,CAP,AKS,TSL,DZ,DELT,TSLEND,CAP1,AKS1)
  !
  !  calculate temperature in roof/building/road
  !  multi-layer heat equation model
  !  Solving Heat Equation by Tri Diagonal Matrix Algorithm
  !-------------------------------------------------------------------

    implicit none

    real(RP), intent(in)    :: G0
    real(RP), intent(in)    :: CAP
    real(RP), intent(in)    :: AKS
    real(RP), intent(in)    :: CAP1      ! for 1st layer
    real(RP), intent(in)    :: AKS1      ! for 1st layer
    real(DP), intent(in)    :: DELT      ! Time step [ s ]
    real(RP), intent(in)    :: TSLEND
    integer,  intent(in)    :: KM
    integer,  intent(in)    :: BOUND
    real(RP), intent(in)    :: DZ(KM)
    real(RP), intent(inout) :: TSL(KM)
    real(RP)                :: A(KM), B(KM), C(KM), D(KM), X(KM), P(KM), Q(KM)
    real(RP)                :: DZEND
    integer                 :: K

    DZEND = DZ(KM)

    A(1)  = 0.0_RP

    B(1)  = CAP1 * DZ(1) / DELT &
          + 2.0_RP * AKS1 / (DZ(1)+DZ(2))
    C(1)  = -2.0_RP * AKS1 / (DZ(1)+DZ(2))
    D(1)  = CAP1 * DZ(1) / DELT * TSL(1) + G0

    do K = 2, KM-1
      A(K) = -2.0_RP * AKS / (DZ(K-1)+DZ(K))
      B(K) = CAP * DZ(K) / DELT + 2.0_RP * AKS / (DZ(K-1)+DZ(K)) + 2.0_RP * AKS / (DZ(K)+DZ(K+1))
      C(K) = -2.0_RP * AKS / (DZ(K)+DZ(K+1))
      D(K) = CAP * DZ(K) / DELT * TSL(K)
    end do

    if( BOUND == 1 ) then ! Flux=0
      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      C(KM) = 0.0_RP
      D(KM) = CAP * DZ(KM) / DELT * TSL(KM)
    else ! T=constant
      A(KM) = -2.0_RP * AKS / (DZ(KM-1)+DZ(KM))
      B(KM) = CAP * DZ(KM) / DELT + 2.0_RP * AKS / (DZ(KM-1)+DZ(KM)) + 2.0_RP * AKS / (DZ(KM)+DZEND)
      C(KM) = 0.0_RP
      D(KM) = CAP * DZ(KM) / DELT * TSL(KM) + 2.0_RP * AKS * TSLEND / (DZ(KM)+DZEND)
    end if

    P(1) = -C(1) / B(1)
    Q(1) =  D(1) / B(1)

    do K = 2, KM
      P(K) = -C(K) / ( A(K) * P(K-1) + B(K) )
      Q(K) = ( -A(K) * Q(K-1) + D(K) ) / ( A(K) * P(K-1) + B(K) )
    end do

    X(KM) = Q(KM)

    do K = KM-1, 1, -1
      X(K) = P(K) * X(K+1) + Q(K)
    end do

    do K = 1, KM
      TSL(K) = X(K)
    enddo

    return
  end subroutine multi_layer2

  !-----------------------------------------------------------------------------
  subroutine urban_param_setup
    implicit none

    real(RP) :: DHGT,THGT,VFWS,VFGS
    integer  :: k

    ! initialize
    R    = 0.0_RP
    RW   = 0.0_RP
    HGT  = 0.0_RP
    Z0HR = 0.0_RP
    Z0HB = 0.0_RP
    Z0HG = 0.0_RP
    Z0C  = 0.0_RP
    Z0HC = 0.0_RP
    ZDC  = 0.0_RP
    SVF  = 0.0_RP

    ! set up other urban parameters
    Z0HR = 0.1_RP * Z0R
    Z0HB = 0.1_RP * Z0B
    Z0HG = 0.1_RP * Z0G
    ZDC  = ZR * 0.3_RP
    Z0C  = ZR * 0.15_RP
    Z0HC = 0.1_RP * Z0C

    ! HGT:  Normalized height
    HGT  = ZR / ( ROAD_WIDTH + ROOF_WIDTH )

    ! R:  Normalized Roof Width (a.k.a. "building coverage ratio")
    R    = ROOF_WIDTH / ( ROAD_WIDTH + ROOF_WIDTH )
    RW   = 1.0_RP - R

    ! Calculate Sky View Factor:
    DHGT = HGT / 100.0_RP
    THGT = 0.0_RP
    VFWS = 0.0_RP
    THGT = HGT - DHGT / 2.0_RP
    do k = 1, 99
      THGT  = THGT - DHGT
      VFWS = VFWS + 0.25_RP * ( 1.0_RP - THGT / sqrt( THGT**2 + RW**2 ) )
    end do

    VFWS = VFWS / 99.0_RP
    VFWS = VFWS * 2.0_RP
    VFGS = 1.0_RP - 2.0_RP * VFWS * HGT / RW
    SVF  = VFGS

    return
  end subroutine urban_param_setup

end module scale_urban_phy_slc
