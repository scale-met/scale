!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Urban Surface fluxes
!!
!! @par Description
!!          Surface fluxes between atmosphere and urban
!!          based on Single-layer Urban Canopy Model (Kusaka et al. 2000, BLM)
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
! !OCL SERIAL
module scale_cpl_atmos_urban_bulk
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
  public :: CPL_AtmUrb_bulk_setup
  public :: CPL_AtmUrb_bulk_restart
  public :: CPL_AtmUrb_bulk
  public :: CPL_AtmUrb_bulk_momentum

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
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
  ! from namelist
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
  subroutine CPL_AtmUrb_bulk_setup( CPL_TYPE_AtmUrb )
    use scale_cpl_bulkflux, only: &
       CPL_bulkflux_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: CPL_TYPE_AtmUrb

    real(RP) :: URBAN_UCM_ZR
    real(RP) :: URBAN_UCM_roof_width
    real(RP) :: URBAN_UCM_road_width
    real(RP) :: URBAN_UCM_SIGMA_ZED
    real(RP) :: URBAN_UCM_AH
    real(RP) :: URBAN_UCM_ALH
    real(RP) :: URBAN_UCM_BETR
    real(RP) :: URBAN_UCM_BETB
    real(RP) :: URBAN_UCM_BETG
    real(RP) :: URBAN_UCM_STRGR
    real(RP) :: URBAN_UCM_STRGB
    real(RP) :: URBAN_UCM_STRGG
    real(RP) :: URBAN_UCM_CAPR
    real(RP) :: URBAN_UCM_CAPB
    real(RP) :: URBAN_UCM_CAPG
    real(RP) :: URBAN_UCM_AKSR
    real(RP) :: URBAN_UCM_AKSB
    real(RP) :: URBAN_UCM_AKSG
    real(RP) :: URBAN_UCM_ALBR
    real(RP) :: URBAN_UCM_ALBB
    real(RP) :: URBAN_UCM_ALBG
    real(RP) :: URBAN_UCM_EPSR
    real(RP) :: URBAN_UCM_EPSB
    real(RP) :: URBAN_UCM_EPSG
    real(RP) :: URBAN_UCM_Z0R
    real(RP) :: URBAN_UCM_Z0B
    real(RP) :: URBAN_UCM_Z0G
    real(RP) :: URBAN_UCM_TRLEND
    real(RP) :: URBAN_UCM_TBLEND
    real(RP) :: URBAN_UCM_TGLEND
    real(RP),allocatable :: URBAN_UCM_DZR(:)
    real(RP),allocatable :: URBAN_UCM_DZB(:)
    real(RP),allocatable :: URBAN_UCM_DZG(:)
    integer  :: URBAN_UCM_BOUND

    NAMELIST / PARAM_CPL_ATMURB_BULK / &
       URBAN_UCM_ZR,         &
       URBAN_UCM_roof_width, &
       URBAN_UCM_road_width, &
       URBAN_UCM_SIGMA_ZED,  &
       URBAN_UCM_AH,         &
       URBAN_UCM_ALH,        &
       URBAN_UCM_BETR,       &
       URBAN_UCM_BETB,       &
       URBAN_UCM_BETG,       &
       URBAN_UCM_STRGR,      &
       URBAN_UCM_STRGB,      &
       URBAN_UCM_STRGG,      &
       URBAN_UCM_CAPR,       &
       URBAN_UCM_CAPB,       &
       URBAN_UCM_CAPG,       &
       URBAN_UCM_AKSR,       &
       URBAN_UCM_AKSB,       &
       URBAN_UCM_AKSG,       &
       URBAN_UCM_ALBR,       &
       URBAN_UCM_ALBB,       &
       URBAN_UCM_ALBG,       &
       URBAN_UCM_EPSR,       &
       URBAN_UCM_EPSB,       &
       URBAN_UCM_EPSG,       &
       URBAN_UCM_Z0R,        &
       URBAN_UCM_Z0B,        &
       URBAN_UCM_Z0G,        &
       URBAN_UCM_TRLEND,     &
       URBAN_UCM_TBLEND,     &
       URBAN_UCM_TGLEND,     &
       !URBAN_UCM_DZR,        &
       !URBAN_UCM_DZB,        &
       !URBAN_UCM_DZG,        &
       URBAN_UCM_BOUND

    integer :: ierr
    !---------------------------------------------------------------------------

    allocate( URBAN_UCM_DZR(UKS:UKE) )
    allocate( URBAN_UCM_DZB(UKS:UKE) )
    allocate( URBAN_UCM_DZG(UKS:UKE) )

    if( IO_L ) write(IO_FID_LOG,*) ''
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[AtmUrb bulk] / Categ[COUPLER] / Origin[SCALElib]'

    if ( CPL_TYPE_AtmUrb /= 'BULK' ) then
       write(*,*) 'xxx CPL_TYPE_AtmUrb is not BULK. Check!'
       call PRC_MPIstop
    endif

    URBAN_UCM_ZR           = ZR
    URBAN_UCM_roof_width   = roof_width
    URBAN_UCM_road_width   = road_width
    URBAN_UCM_SIGMA_ZED    = SIGMA_ZED
    URBAN_UCM_AH           = AH
    URBAN_UCM_ALH          = ALH
    URBAN_UCM_BETR         = BETR
    URBAN_UCM_BETB         = BETB
    URBAN_UCM_BETG         = BETG
    URBAN_UCM_STRGR        = STRGR
    URBAN_UCM_STRGB        = STRGB
    URBAN_UCM_STRGG        = STRGG
    URBAN_UCM_CAPR         = CAPR
    URBAN_UCM_CAPB         = CAPB
    URBAN_UCM_CAPG         = CAPG
    URBAN_UCM_AKSR         = AKSR
    URBAN_UCM_AKSB         = AKSB
    URBAN_UCM_AKSG         = AKSG
    URBAN_UCM_ALBR         = ALBR
    URBAN_UCM_ALBB         = ALBB
    URBAN_UCM_ALBG         = ALBG
    URBAN_UCM_EPSR         = EPSR
    URBAN_UCM_EPSB         = EPSB
    URBAN_UCM_EPSG         = EPSG
    URBAN_UCM_Z0R          = Z0R
    URBAN_UCM_Z0B          = Z0B
    URBAN_UCM_Z0G          = Z0G
    URBAN_UCM_TRLEND       = TRLEND
    URBAN_UCM_TBLEND       = TBLEND
    URBAN_UCM_TGLEND       = TGLEND
    URBAN_UCM_BOUND        = BOUND
    URBAN_UCM_DZR(UKS:UKE) = (/0.01_RP,0.01_RP,0.03_RP,0.05_RP,0.10_RP/)
    URBAN_UCM_DZB(UKS:UKE) = (/0.01_RP,0.01_RP,0.03_RP,0.05_RP,0.10_RP/)
    URBAN_UCM_DZG(UKS:UKE) = (/0.01_RP,0.01_RP,0.03_RP,0.05_RP,0.10_RP/)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_ATMURB_BULK,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_ATMURB_BULK. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_CPL_ATMURB_BULK)

    ZR           = URBAN_UCM_ZR
    roof_width   = URBAN_UCM_roof_width
    road_width   = URBAN_UCM_road_width
    SIGMA_ZED    = URBAN_UCM_SIGMA_ZED
    AH           = URBAN_UCM_AH
    ALH          = URBAN_UCM_ALH
    BETR         = URBAN_UCM_BETR
    BETB         = URBAN_UCM_BETB
    BETG         = URBAN_UCM_BETG
    STRGR        = URBAN_UCM_STRGR
    STRGB        = URBAN_UCM_STRGB
    STRGG        = URBAN_UCM_STRGG
    CAPR         = URBAN_UCM_CAPR
    CAPB         = URBAN_UCM_CAPB
    CAPG         = URBAN_UCM_CAPG
    AKSR         = URBAN_UCM_AKSR
    AKSB         = URBAN_UCM_AKSB
    AKSG         = URBAN_UCM_AKSG
    ALBR         = URBAN_UCM_ALBR
    ALBB         = URBAN_UCM_ALBB
    ALBG         = URBAN_UCM_ALBG
    EPSR         = URBAN_UCM_EPSR
    EPSB         = URBAN_UCM_EPSB
    EPSG         = URBAN_UCM_EPSG
    Z0R          = URBAN_UCM_Z0R
    Z0B          = URBAN_UCM_Z0B
    Z0G          = URBAN_UCM_Z0G
    TRLEND       = URBAN_UCM_TRLEND
    TBLEND       = URBAN_UCM_TBLEND
    TGLEND       = URBAN_UCM_TGLEND
    BOUND        = URBAN_UCM_BOUND

    ahdiurnal(:) = (/ 0.356, 0.274, 0.232, 0.251, 0.375, 0.647, 0.919, 1.135, 1.249, 1.328, &
                      1.365, 1.363, 1.375, 1.404, 1.457, 1.526, 1.557, 1.521, 1.372, 1.206, &
                      1.017, 0.876, 0.684, 0.512                                            /)

    allocate( DZR(UKS:UKE) )
    allocate( DZB(UKS:UKE) )
    allocate( DZG(UKS:UKE) )

    DZR(UKS:UKE) = URBAN_UCM_DZR(UKS:UKE)
    DZB(UKS:UKE) = URBAN_UCM_DZB(UKS:UKE)
    DZG(UKS:UKE) = URBAN_UCM_DZG(UKS:UKE)

    ! set other urban parameters
    call urban_param_setup

    ! set up bulk coefficient function
    call CPL_bulkflux_setup

    return
  end subroutine CPL_AtmUrb_bulk_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmUrb_bulk_restart( &
        TR,           & ! (in)
        TB,           & ! (in)
        TG,           & ! (in)
        TC,           & ! (in)
        QC,           & ! (in)
        UC,           & ! (in)
        TRL,          & ! (in)
        TBL,          & ! (in)
        TGL,          & ! (in)
        RAINR,        & ! (in)
        RAINB,        & ! (in)
        RAING,        & ! (in)
        AH_t,         & ! (in)
        ALH_t,        & ! (in)
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
        PRES,         & ! (in)
        TA,           & ! (in)
        QA,           & ! (in)
        UA,           & ! (in)
        U1,           & ! (in)
        V1,           & ! (in)
        ZA,           & ! (in)
        SSG,          & ! (in)
        LLG,          & ! (in)
        RHOO,         & ! (in)
        XLON,         & ! (in)
        XLAT          ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       EPS    => CONST_EPS,     &    ! small number (machine epsilon)
       D2R    => CONST_D2R,     &    ! degree to radian
       KARMAN => CONST_KARMAN,  &    ! kalman constant  [-]
       PI     => CONST_PI,      &    ! pi               [-]
       CPdry  => CONST_CPdry,   &    ! heat capacity of dry air [J/K/kg]
       LHV0   => CONST_LHV0,    &    ! latent heat of vaporization [J/kg]
       GRAV   => CONST_GRAV,    &    !< gravitational constant [m/s2]
       Rdry   => CONST_Rdry,    &    !< specific gas constant (dry) [J/kg/K]
       Rvap   => CONST_Rvap,    &    !< gas constant (water vapor) [J/kg/K]
       STB    => CONST_STB,     &    !< stefan-Boltzman constant [MKS unit]
       TEM00  => CONST_TEM00         !< temperature reference (0 degree C) [K]
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_time, only:       &
       DELT => TIME_DTSEC_URBAN      !< time interval of urban step [sec]
    implicit none

    !-- configuration variables
    logical, intent(in) :: LSOLAR  ! logical   [true=both, false=SSG only]

    !-- Input variables from Coupler to Urban
    real(RP), intent(in)    :: PRES    ! Surface Pressure                       [Pa]
    real(RP), intent(in)    :: TA      ! temp at 1st atmospheric level          [K]
    real(RP), intent(in)    :: QA      ! specific humidity at 1st atmospheric level  [kg/kg]
    real(RP), intent(in)    :: UA      ! wind speed at 1st atmospheric level    [m/s]
    real(RP), intent(in)    :: U1      ! u at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: V1      ! v at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: ZA      ! height of 1st atmospheric level        [m]
    real(RP), intent(in)    :: SSG     ! downward total short wave radiation    [W/m/m]
    real(RP), intent(in)    :: LLG     ! downward long wave radiation           [W/m/m]
    real(RP), intent(in)    :: RHOO    ! air density                            [kg/m^3]
    real(RP), intent(in)    :: XLAT    !< latitude                              [rad,-pi,pi]
    real(RP), intent(in)    :: XLON    !< longitude                             [rad,0-2pi]

    !-- In/Out variables from/to Coupler to/from Urban
    real(RP), intent(in)    :: TR      ! roof temperature              [K]
    real(RP), intent(in)    :: TB      ! building wall temperature     [K]
    real(RP), intent(in)    :: TG      ! road temperature              [K]
    real(RP), intent(in)    :: TC      ! urban-canopy air temperature  [K]
    real(RP), intent(in)    :: QC      ! urban-canopy air specific humidity [kg/kg]
    real(RP), intent(in)    :: UC      ! diagnostic canopy wind        [m/s]
    real(RP), intent(in)    :: TRL(UKS:UKE)  ! layer temperature       [K]
    real(RP), intent(in)    :: TBL(UKS:UKE)  ! layer temperature       [K]
    real(RP), intent(in)    :: TGL(UKS:UKE)  ! layer temperature       [K]
    real(RP), intent(in)    :: RAINR   ! rain amount in storage on roof     [kg/m2]
    real(RP), intent(in)    :: RAINB   ! rain amount in storage on building [kg/m2]
    real(RP), intent(in)    :: RAING   ! rain amount in storage on road     [kg/m2]
    real(RP), intent(in)    :: AH_t    ! Sensible Anthropogenic heat        [W/m^2]
    real(RP), intent(in)    :: ALH_t   ! Latent Anthropogenic heat          [W/m^2]

    !-- Output variables from Urban to Coupler
    real(RP), intent(out)   :: ALBD_SW_grid  ! grid mean of surface albedo for SW
    real(RP), intent(out)   :: ALBD_LW_grid  ! grid mean of surface albedo for LW ( 1-emiss )
    real(RP), intent(out)   :: RTS           ! radiative surface temperature    [K]
    real(RP), intent(out)   :: SH            ! sensible heat flux               [W/m/m]
    real(RP), intent(out)   :: LH            ! latent heat flux                 [W/m/m]
    real(RP), intent(out)   :: GHFLX         ! heat flux into the ground        [W/m/m]
    real(RP), intent(out)   :: RN            ! net radition                     [W/m/m]
    real(RP), intent(out)   :: U10           ! U wind at 10m                    [m/s]
    real(RP), intent(out)   :: V10           ! V wind at 10m                    [m/s]
    real(RP), intent(out)   :: T2            ! air temperature at 2m            [K]
    real(RP), intent(out)   :: Q2            ! specific humidity at 2m          [kg/kg]
    real(RP), intent(out)   :: RNR, RNB, RNG
    real(RP), intent(out)   :: SHR, SHB, SHG
    real(RP), intent(out)   :: LHR, LHB, LHG
    real(RP), intent(out)   :: GHR, GHB, GHG

    !-- parameters
    real(RP), parameter     :: SRATIO = 0.75_RP  ! ratio between direct/total solar [-]

    !-- Local variables
    logical  :: SHADOW = .false.
           !  true  = consider svf and shadow effects,
           !  false = consider svf effect only

    real(RP) :: LON, LAT  ! longitude [deg], latitude [deg]

    real(RP) :: SSGD      ! downward direct short wave radiation   [W/m/m]
    real(RP) :: SSGQ      ! downward diffuse short wave radiation  [W/m/m]

    real(RP) :: W, VFGS, VFGW, VFWG, VFWS, VFWW
    real(RP) :: SX, RX

    !real(RP) :: UST, TST, QST

    real(RP) :: LUP, LDN, RUP
    real(RP) :: SUP, SDN

    real(RP) :: LNET, SNET, FLXUV
    real(RP) :: psim,psim2,psim10   ! similality stability shear function for momentum
    real(RP) :: psih,psih2,psih10   ! similality stability shear function for heat

    real(RP) :: SR, SB, SG, RR, RB, RG
    real(RP) :: SB1, SB2, SG1, SG2
    real(RP) :: RB1, RB2, RG1, RG2

    real(RP) :: Z
    real(RP) :: QS0R,QS0B,QS0G

    real(RP) :: RIBC, BHC, CDC
    real(RP) :: RIBR, BHR, CDR
    real(RP) :: ALPHAB, ALPHAG
    real(RP) :: CHR, CHB, CHG, CHC

    real(RP) :: XXX, X, CD, CH, CHU, XXX2, XXX10

    !-----------------------------------------------------------
    ! Set parameters
    !-----------------------------------------------------------

     ! local time at center point
     LAT = XLAT / D2R
     LON = XLON / D2R

     !if(.NOT.LSOLAR) then   ! Radiation scheme does not have SSGD and SSGQ.
       SSGD = SRATIO * SSG   ! downward direct short wave radiation
       SSGQ = SSG - SSGD     ! downward diffuse short wave radiation
     !endif

     W    = 2.0_RP * 1.0_RP * HGT
     VFGS = SVF
     VFGW = 1.0_RP - SVF
     VFWG = ( 1.0_RP - SVF ) * ( 1.0_RP - R ) / W
     VFWS = VFWG
     VFWW = 1.0_RP - 2.0_RP * VFWG

     SX  = (SSGD+SSGQ)   ! downward short wave radition [W/m2]
     RX  = LLG           ! downward long wave radiation

    !-----------------------------------------------------------
    ! Radiation : Net Short Wave Radiation at roof/wall/road
    !-----------------------------------------------------------

     if( SSG > 0.0_RP ) then !  SSG is downward short

      ! currently we use no shadow effect model
      !!     IF(.NOT.SHADOW) THEN              ! no shadow effects model

      SG1  = SX * VFGS * ( 1.0_RP - ALBG )
      SB1  = SX * VFWS * ( 1.0_RP - ALBB )
      SG2  = SB1 * ALBB / ( 1.0_RP - ALBB ) * VFGW * ( 1.0_RP - ALBG )
      SB2  = SG1 * ALBG / ( 1.0_RP - ALBG ) * VFWG * ( 1.0_RP - ALBB )

      SR   = SX * ( 1.0_RP - ALBR )
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
    ! Set evaporation efficiency on roof/wall/road
    !-----------------------------------------------------------

    if ( STRGR == 0.0_RP ) then
       BETR = 0.0_RP
    else
       BETR = min ( RAINR / STRGR, 1.0_RP)
    endif
    if ( STRGB == 0.0_RP ) then
       BETB = 0.0_RP
    else
       BETB = min ( RAINB / STRGB, 1.0_RP)
    endif
    if ( STRGG == 0.0_RP ) then
       BETG = 0.0_RP
    else
       BETG = min ( RAING / STRGG, 1.0_RP)
    endif

    !-----------------------------------------------------------
    ! Energy balance on roof/wall/road surface
    !-----------------------------------------------------------

    !--- Roof
    Z    = ZA - ZDC
    BHR  = LOG(Z0R/Z0HR) / 0.4_RP
    RIBR = ( GRAV * 2.0_RP / (TA+TR) ) * (TA-TR) * (Z+Z0R) / (UA*UA)
xxxr=0.0_RP ! temtative
    call mos(XXXR,CHR,CDR,BHR,RIBR,Z,Z0R,UA,TA,TR,RHOO)

    call qsat( QS0R, TR, PRES )

    RR   = EPSR * ( RX - STB * (TR**4) )
    SHR  = RHOO * CPdry * CHR * UA * (TR-TA)
    LHR  = RHOO * LHV0 * CHR * UA * BETR * (QS0R-QA)
    GHR  = SR + RR - SHR - LHR

    !--- Wall and Road

    ! empirical form
    !ALPHAB = RHOO * CPdry * ( 6.15_RP + 4.18_RP * UC ) / 1200.0_RP
    !if( UC > 5.0_RP ) ALPHAB = RHOO * CPdry * ( 7.51_RP * UC**0.78_RP ) / 1200.0_RP
    !ALPHAG = RHOO * CPdry * ( 6.15_RP + 4.18_RP * UC ) / 1200.0_RP
    !if( UC > 5.0_RP ) ALPHAG = RHOO * CPdry * ( 7.51_RP * UC**0.78_RP ) / 1200.0_RP

     if ( (UC-0.0) < sqrt(EPS) ) then
       write(*,*) 'UC value is too small. Check!'
       call PRC_MPIstop
     endif

     ALPHAB = 6.15_RP + 4.18_RP * UC
     if( UC > 5.0_RP ) ALPHAB = 7.51_RP * (UC**0.78_RP )
     ALPHAG = 6.15_RP + 4.18_RP * UC
     if( UC > 5.0_RP ) ALPHAG = 7.51_RP * (UC**0.78_RP )

     CHB = ALPHAB / RHOO / CPdry / UC
     CHG = ALPHAG / RHOO / CPdry / UC

     RG1      = EPSG * ( RX * VFGS                   &
                       + EPSB * VFGW * STB * TB**4   &
                       - STB * TG**4                )
     RB1      = EPSB * ( RX * VFWS                   &
                       + EPSG * VFWG * STB * TG**4   &
                       + EPSB * VFWW * STB * TB**4   &
                       - STB * TB**4                )

     RG2      = EPSG * ( (1.0_RP-EPSB) * VFGW * VFWS * RX                  &
                       + (1.0_RP-EPSB) * VFGW * VFWG * EPSG * STB * TG**4  &
                       + EPSB * (1.0_RP-EPSB) * VFGW * VFWW * STB * TB**4  )
     RB2      = EPSB * ( (1.0_RP-EPSG) * VFWG * VFGS * RX                                  &
                       + (1.0_RP-EPSG) * EPSB * VFGW * VFWG * STB * TB**4                  &
                       + (1.0_RP-EPSB) * VFWS * VFWW * RX                                  &
                       + (1.0_RP-EPSB) * VFWG * VFWW * STB * EPSG * TG**4                  &
                       + EPSB * (1.0_RP-EPSB) * VFWW * (1.0_RP-2.0_RP*VFWS) * STB * TB**4  )

     RG       = RG1 + RG2
     RB       = RB1 + RB2

     call qsat( QS0B, TB, PRES )

     SHB  = RHOO * CPdry * CHB * UC * (TB-TC)
     LHB  = RHOO * LHV0  * CHB * UC * BETB * (QS0B-QC)
     GHB  =  SB + RB - SHB - LHB

     call qsat( QS0G, TG, PRES )

     SHG  = RHOO * CPdry * CHG * UC * (TG-TC)
     LHG  = RHOO * LHV0  * CHG * UC * BETG * (QS0G-QC)
     GHG  =  SG + RG - SHG - LHG

    !-----------------------------------------------------------
    ! Total Fluxes from Urban Canopy
    !-----------------------------------------------------------

     ! comment out below line to avoid using unintialized value
     !FLXUV  = ( R*CDR + RW*CDC ) * UA * UA
     SH     =  R*SHR  + W*SHB  + RW*SHG               ! Sensible heat flux [W/m/m]
     LH     =  R*LHR  + W*LHB  + RW*LHG               ! Latent heat flux [W/m/m]
     GHFLX  =  -1.0_RP * ( R*GHR + W*GHB + RW*GHG )
     LNET   =  R*RR   + W*RB   + RW*RG

    !-----------------------------------------------------------
    ! Net Radiation and surface albedo
    !-----------------------------------------------------------
     RNR = SR + RR        ! Net radiation on roof [W/m/m]
     RNB = SB + RB        ! Net radiation on building [W/m/m]
     RNG = SG + RG        ! Net radiation on ground [W/m/m]
     RN = (SNET+LNET)     ! Net radiation [W/m/m]

     !--- shortwave radiation
     SDN =  R + W * (VFWS + VFGS * ALBG * VFWG)  + RW * (VFGS + VFWS * ALBB * VFGW)
     SUP =  R  * ALBR  &
          + W  * ( VFWS* ALBB + VFGS * ALBG * VFWG *ALBB ) &
          + RW * ( VFGS * ALBG + VFWS * ALBB * VFGW *ALBG )

     ALBD_SW_grid = SUP / SDN

     !--- longwave radiation
     LDN =  R + W*VFWS + RW*VFGS
     LUP =  R  * (1.0_RP-EPSR) &
          + W  * ( (1.0_RP-EPSB*VFWW)*(1.0_RP-EPSB)*VFWS - EPSB*VFWG*(1.0_RP-EPSG)*VFGS )  &
          + RW * ( (1.0_RP-EPSG)*VFGS - EPSG*(1.0_RP-VFGS)*(1.0_RP-EPSB)*VFWS )

     RUP = (LDN - LUP) * RX - LNET
     ALBD_LW_grid = LUP / LDN

    !-----------------------------------------------------------
    !  diagnostic GRID AVERAGED TS from upward logwave
    !-----------------------------------------------------------

     RTS = ( RUP / STB / ( 1.0_RP-ALBD_LW_grid) )**0.25

    !-----------------------------------------------------------
    !  diagnostic grid average U10, V10, T2, Q2 from urban
    !  Below method would be better to be improved. This is tentative method.
    !-----------------------------------------------------------
    ! comment out below three lines to avoid using unintialized value
    !UST = sqrt( FLXUV )              ! u* [m/s]
    !TST = -SH / RHOO / CPdry / UST   ! T* [K]
    !QST = -LH / RHOO / LHV0 / UST    ! q* [-]
     Z    = ZA - ZDC
     BHC  = LOG(Z0C/Z0HC) / 0.4_RP
     RIBC = ( GRAV * 2.0_RP / (TA+TC) ) * (TA-TC) * (Z+Z0C) / (UA*UA)
     call mos(XXXC,CHC,CDC,BHC,RIBC,Z,Z0C,UA,TA,TC,RHOO)

     XXX  = XXXC
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
    ! Add anthropogenic heat fluxes
    !-----------------------------------------------------------

     SH     = SH + AH_t     ! Sensible heat flux  [W/m/m]
     LH     = LH + ALH_t    ! Latent heat flux    [W/m/m]

    return
  end subroutine CPL_AtmUrb_bulk_restart

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmUrb_bulk( &
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
        AH_t,         & ! (inout)
        ALH_t,        & ! (inout)
        ALBD_LW_grid, & ! (out)
        ALBD_SW_grid, & ! (out)
        SHR,          & ! (out)
        SHB,          & ! (out)
        SHG,          & ! (out)
        LHR,          & ! (inout)
        LHB,          & ! (inout)
        LHG,          & ! (inout)
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
        PRES,         & ! (in)
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
        XLAT          ) ! (in)
    use scale_grid_index
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       EPS    => CONST_EPS,     &    ! small number (machine epsilon)
       PI     => CONST_PI,      &    ! pi               [-]
       D2R    => CONST_D2R,     &    ! degree to radian
       KARMAN => CONST_KARMAN,  &    ! Kalman constant  [-]
       CPdry  => CONST_CPdry,   &    ! Heat capacity of dry air [J/K/kg]
       LHV0   => CONST_LHV0,    &    ! Latent heat of vaporization [J/kg]
       GRAV   => CONST_GRAV,    &    ! gravitational constant [m/s2]
       Rdry   => CONST_Rdry,    &    ! specific gas constant (dry) [J/kg/K]
       Rvap   => CONST_Rvap,    &    ! gas constant (water vapor) [J/kg/K]
       STB    => CONST_STB,     &    ! stefan-Boltzman constant [MKS unit]
       TEM00  => CONST_TEM00         ! temperature reference (0 degree C) [K]
    use scale_atmos_saturation, only: &
       qsat   => ATMOS_SATURATION_pres2qsat_all
    use scale_time, only:       &
       NOWDATE => TIME_NOWDATE, &    ! current time [YYYY MM DD HH MM SS]
       DELT    => TIME_DTSEC_URBAN   ! time interval of urban step [sec]
    implicit none

    !-- configuration variables
    logical , intent(in)    :: LSOLAR ! logical   [true=both, false=SSG only]

    !-- Input variables from Coupler to Urban
    real(RP), intent(in)    :: PRES ! Surface Pressure                       [Pa]
    real(RP), intent(in)    :: TA   ! temp at 1st atmospheric level          [K]
    real(RP), intent(in)    :: QA   ! specific humidity at 1st atmospheric level  [kg/kg]
    real(RP), intent(in)    :: UA   ! wind speed at 1st atmospheric level    [m/s]
    real(RP), intent(in)    :: U1   ! u at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: V1   ! v at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: ZA   ! height of 1st atmospheric level        [m]
    real(RP), intent(in)    :: SSG  ! downward total short wave radiation    [W/m/m]
    real(RP), intent(in)    :: LLG  ! downward long wave radiation           [W/m/m]
    real(RP), intent(in)    :: RAIN ! precipitation flux                     [kg/m2/s]
    real(RP), intent(in)    :: RHOO ! air density                            [kg/m^3]
    real(RP), intent(in)    :: XLAT ! latitude                               [rad,-pi,pi]
    real(RP), intent(in)    :: XLON ! longitude                              [rad,0-2pi]

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
    real(RP), intent(inout) :: AH_t  ! Sensible Anthropogenic heat [W/m^2]
    real(RP), intent(inout) :: ALH_t ! Latent Anthropogenic heat   [W/m^2]

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
    real(RP), intent(inout) :: LHR, LHB, LHG
    real(RP), intent(out)   :: GHR, GHB, GHG

    !-- parameters
    real(RP), parameter     :: SRATIO    = 0.75_RP     ! ratio between direct/total solar [-]
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

    real(RP) :: SSGD     ! downward direct short wave radiation   [W/m/m]
    real(RP) :: SSGQ     ! downward diffuse short wave radiation  [W/m/m]

    real(RP) :: W, VFGS, VFGW, VFWG, VFWS, VFWW
    real(RP) :: SX, RX

    real(RP) :: TRP = 350.0_RP   ! TRP: at previous time step [K]
    real(RP) :: TBP = 350.0_RP   ! TBP: at previous time step [K]
    real(RP) :: TGP = 350.0_RP   ! TGP: at previous time step [K]
    real(RP) :: TCP = 350.0_RP   ! TCP: at previous time step [K]
    real(RP) :: QCP = 0.01_RP    ! QCP: at previous time step [kg/kg]
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

    real(RP) :: resi1,resi2     ! residual
    real(RP) :: G0RP,G0BP,G0GP

    real(RP) :: XXX, X, CD, CH, CHU, XXX2, XXX10

    integer  :: iteration

    !-----------------------------------------------------------
    ! Set parameters
    !-----------------------------------------------------------

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

    if ( ZDC + Z0C + 2.0_RP >= ZA ) then
       if( IO_L ) write(IO_FID_LOG,*) 'ZDC + Z0C + 2m is larger than the 1st WRF level' // &
                                      'Stop in subroutine urban - change ZDC and Z0C'
       call PRC_MPIstop
    endif

    !if(.NOT.LSOLAR) then   ! Radiation scheme does not have SSGD and SSGQ.
      SSGD = SRATIO * SSG   ! downward direct short wave radiation
      SSGQ = SSG - SSGD     ! downward diffuse short wave radiation
    !endif

    W    = 2.0_RP * 1.0_RP * HGT
    VFGS = SVF
    VFGW = 1.0_RP - SVF
    VFWG = ( 1.0_RP - SVF ) * ( 1.0_RP - R ) / W
    VFWS = VFWG
    VFWW = 1.0_RP - 2.0_RP * VFWG

    SX  = (SSGD+SSGQ)   ! downward short wave radition [W/m/m]
    RX  = LLG           ! downward long wave radiation

    !--- calculate canopy wind

    call canopy_wind(ZA, UA, UC)

    !-----------------------------------------------------------
    ! Set evaporation efficiency on roof/wall/road
    !-----------------------------------------------------------

    !!--- calculate rain amount remaining on the surface
    RAINR = max(0.0_RP, RAINR-(LHR/LHV0)*DELT)   ! [kg/m/m = mm]
    RAINB = max(0.0_RP, RAINB-(LHB/LHV0)*DELT)   ! [kg/m/m = mm]
    RAING = max(0.0_RP, RAING-(LHG/LHV0)*DELT)   ! [kg/m/m = mm]

    !!--- calculate evaporation efficiency
    RAINT = 1.0_RP * ( RAIN * DELT )            ! [kg/m2/s -> kg/m2]
    call cal_beta(BETR, RAINT, RAINR, STRGR, ROFFR)

    RAINT = 0.1_RP * ( RAIN * DELT )
    call cal_beta(BETB, RAINT, RAINB, STRGB, ROFFB)

    RAINT = 0.9_RP * ( RAIN * DELT )
    call cal_beta(BETG, RAINT, RAING, STRGG, ROFFG)

    ROFF = ROFF +  R * ROFFR  + RW * ( ROFFB + ROFFG )

    !-----------------------------------------------------------
    ! Radiation : Net Short Wave Radiation at roof/wall/road
    !-----------------------------------------------------------

    if( SSG > 0.0_RP ) then !  SSG is downward short

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
     do iteration = 1, 20

      Z    = ZA - ZDC
      BHR  = LOG(Z0R/Z0HR) / 0.4_RP
      RIBR = ( GRAV * 2.0_RP / (TA+TR) ) * (TA-TR) * (Z+Z0R) / (UA*UA)
      call mos(XXXR,CHR,CDR,BHR,RIBR,Z,Z0R,UA,TA,TR,RHOO)

      call qsat( QS0R, TR, PRES )

      RR    = EPSR * ( RX - STB * (TR**4)  )
      HR    = RHOO * CPdry * CHR * UA * (TR-TA)
      ELER  = RHOO * LHV0  * CHR * UA * BETR * (QS0R-QA)
      G0R   = SR + RR - HR - ELER

    !--- calculate temperature in roof
    !  if ( STRGR /= 0.0_RP ) then
    !    CAPL1 = CAP_water * (RAINR / (DZR(1) + RAINR)) + CAPR * (DZR(1) / (DZR(1) + RAINR))
    !    AKSL1 = AKS_water * (RAINR / (DZR(1) + RAINR)) + AKSR * (DZR(1) / (DZR(1) + RAINR))
    !  else
    !    CAPL1 = CAPR
    !    AKSL1 = AKSR
    !  endif

      TRL = TRLP
      call multi_layer(UKE,BOUND,G0R,CAPR,AKSR,TRL,DZR,DELT,TRLEND)
      !! 1st layer's cap, aks are replaced.
      !! call multi_layer2(UKE,BOUND,G0R,CAPR,AKSR,TRL,DZR,DELT,TRLEND,CAPL1,AKSL1)
      TR  = TRL(1)

      resi1  =  abs(G0R - G0RP)
      G0RP   =  G0R

      if( resi1 < sqrt(EPS) ) exit

     enddo

    !--- update only fluxes ----
     RIBR = ( GRAV * 2.0_RP / (TA+TR) ) * (TA-TR) * (Z+Z0R) / (UA*UA)
     call mos(XXXR,CHR,CDR,BHR,RIBR,Z,Z0R,UA,TA,TR,RHOO)

     call qsat( QS0R, TR, PRES )

     RR      = EPSR * ( RX - STB * (TR**4) )
     HR      = RHOO * CPdry * CHR * UA * (TR-TA)
     ELER    = RHOO * LHV0  * CHR * UA * BETR * (QS0R-QA)
     G0R     = SR + RR - HR - ELER

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
     do iteration = 1, 50

      Z    = ZA - ZDC
      BHC  = LOG(Z0C/Z0HC) / 0.4_RP
      RIBC = ( GRAV * 2.0_RP / (TA+TC) ) * (TA-TC) * (Z+Z0C) / (UA*UA)
      call mos(XXXC,CHC,CDC,BHC,RIBC,Z,Z0C,UA,TA,TC,RHOO)
      ALPHAC = CHC * RHOO * CPdry * UA

      call qsat( QS0B, TB, PRES )
      call qsat( QS0G, TG, PRES )

      TC1   = RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
      TC2   = RW*ALPHAC*TA + RW*ALPHAG*TG + W*ALPHAB*TB
      TC    = TC2 / TC1
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

      HB    = RHOO * CPdry * CHB * UC * (TB-TC)
      ELEB  = RHOO * LHV0  * CHB * UC * BETB * (QS0B-QC)
      G0B   = SB + RB - HB - ELEB

      HG    = RHOO * CPdry * CHG * UC * (TG-TC)
      ELEG  = RHOO * LHV0  * CHG * UC * BETG * (QS0G-QC)
      G0G   = SG + RG - HG - ELEG

    !--- Heat capacity and Thermal conductivity of 1st layer
    !                 and calculate temperature in building/road
     ! if ( STRGB /= 0.0_RP ) then
     !   CAPL1 = CAP_water * (RAINB / (DZB(1) + RAINB)) + CAPB * (DZB(1) / (DZB(1) + RAINB))
     !   AKSL1 = AKS_water * (RAINB / (DZB(1) + RAINB)) + AKSB * (DZB(1) / (DZB(1) + RAINB))
     ! else
     !   CAPL1 = CAPB
     !   AKSL1 = AKSB
     ! endif
      TBL = TBLP
      call multi_layer(UKE,BOUND,G0B,CAPB,AKSB,TBL,DZB,DELT,TBLEND)
      !call multi_layer2(UKE,BOUND,G0B,CAPB,AKSB,TBL,DZB,DELT,TBLEND,CAPL1,AKSL1)
      TB = TBL(1)

      !if ( STRGG /= 0.0_RP ) then
      !  CAPL1 = CAP_water * (RAING / (DZG(1) + RAING)) + CAPG * (DZG(1) / (DZG(1) + RAING))
      !  AKSL1 = AKS_water * (RAING / (DZG(1) + RAING)) + AKSG * (DZG(1) / (DZG(1) + RAING))
      !else
      !  CAPL1 = CAPG
      !  AKSL1 = AKSG
      !endif
      TGL = TGLP
      call multi_layer(UKE,BOUND,G0G,CAPG,AKSG,TGL,DZG,DELT,TGLEND)
      !call multi_layer2(UKE,BOUND,G0G,CAPG,AKSG,TGL,DZG,DELT,TGLEND,CAPL1,AKSL1)
      TG = TGL(1)

      call qsat( QS0B, TB, PRES )
      call qsat( QS0G, TG, PRES )

      TC1    =  RW*ALPHAC    + RW*ALPHAG    + W*ALPHAB
      TC2    =  RW*ALPHAC*TA + RW*ALPHAG*TG + W*ALPHAB*TB
      TC     =  TC2 / TC1
      QC1    =  RW*(CHC*UA)    + RW*(CHG*BETG*UC)      + W*(CHB*BETB*UC)
      QC2    =  RW*(CHC*UA)*QA + RW*(CHG*BETG*UC)*QS0G + W*(CHB*BETB*UC)*QS0B
      QC     =  QC2 / QC1

      resi1  =  abs(G0B - G0BP)
      resi2  =  abs(G0G - G0GP)
      G0BP   =  G0B
      G0GP   =  G0G

      if( resi1 < sqrt(EPS) .and. resi2 < sqrt(EPS) ) exit

    enddo

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

     HB   = RHOO * CPdry * CHB * UC * (TB-TC)
     ELEB = RHOO * LHV0  * CHB * UC * BETB * (QS0B-QC)
     G0B  = SB + RB - HB - ELEB

     HG   = RHOO * CPdry * CHG * UC * (TG-TC)
     ELEG = RHOO * LHV0  * CHG * UC * BETG * (QS0G-QC)
     G0G  = SG + RG - HG - ELEG

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

    LW = LLG - LNET     ! Upward longwave radiation   [W/m/m]
    SW = SSG - SNET     ! Upward shortwave radiation  [W/m/m]
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
  end subroutine CPL_AtmUrb_bulk

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmUrb_bulk_momentum( &
        XMFLX,      & ! (out)
        YMFLX,      & ! (out)
        ZMFLX,      & ! (out)
        Z0M,        & ! (out)
        ZA,         & ! (in)
        RHOA,       & ! (in)
        UA,         & ! (in)
        VA,         & ! (in)
        WA,         & ! (in)
        TMPA,       & ! (in)
        PRSA,       & ! (in)
        QVA,        & ! (in)
        PBL,        & ! (in)
        PRSS,       & ! (in)
        RTS         ) ! (in)
    use scale_atmos_saturation, only: &
       qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_cpl_bulkflux, only: &
       CPL_bulkflux
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    ! arguments
    real(RP), intent(out) :: XMFLX  ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX  ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: ZMFLX  ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: Z0M    ! roughness length for momentum [m]

    real(RP), intent(in)  :: RTS    ! radiative surface temperature [K]
    real(RP), intent(in)  :: RHOA   ! density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in)  :: ZA     ! height at the lowest atmospheric layer [m/s]
    real(RP), intent(in)  :: UA     ! velocity u at the lowest atmospheric layer [m/s]
    real(RP), intent(in)  :: VA     ! velocity v at the lowest atmospheric layer [m/s]
    real(RP), intent(in)  :: WA     ! velocity w at the lowest atmospheric layer [m/s]
    real(RP), intent(in)  :: TMPA   ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in)  :: PRSA   ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in)  :: QVA    ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in)  :: PBL    ! the top of atmospheric mixing layer [m]
    real(RP), intent(in)  :: PRSS   ! pressure at the surface [Pa]

    ! works
    real(RP) :: Ustar ! friction velocity [m]
    real(RP) :: Tstar ! friction temperature [K]
    real(RP) :: Qstar ! friction mixing rate [kg/kg]
    real(RP) :: Uabs  ! modified absolute velocity [m/s]
    real(RP) :: QVS   ! saturation water vapor mixing ratio at surface [kg/kg]
    !---------------------------------------------------------------------------

    ! saturation at the surface
    call qsat( QVS, RTS, PRSS )

    call CPL_bulkflux( Ustar, & ! (out)
                       Tstar, & ! (out)
                       Qstar, & ! (out)
                       Uabs,  & ! (out)
                       TMPA,  & ! (in)
                       RTS,   & ! (in)
                       PRSA,  & ! (in)
                       PRSS,  & ! (in)
                       QVA,   & ! (in)
                       QVS,   & ! (in)
                       UA,    & ! (in)
                       VA,    & ! (in)
                       ZA,    & ! (in)
                       PBL,   & ! (in)
                       Z0C,   & ! (in)
                       Z0HC,  & ! (in)
                       Z0HC   ) ! (in)

    XMFLX  = -RHOA * Ustar**2 / Uabs * UA
    YMFLX  = -RHOA * Ustar**2 / Uabs * VA
    ZMFLX  = -RHOA * Ustar**2 / Uabs * WA

    Z0M = Z0C

    return
  end subroutine CPL_AtmUrb_bulk_momentum

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
    !   BET   = min ( WATER / STRG, 0.35_RP)
    !   BET   = max ( WATER / STRG, 0.1_RP)
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
    real(RP), intent(in)    :: DELT      ! Time step [ s ]
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
    real(RP), intent(in)    :: DELT      ! Time step [ s ]
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

end module scale_cpl_atmos_urban_bulk
