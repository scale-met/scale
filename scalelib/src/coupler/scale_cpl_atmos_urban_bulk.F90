!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Urban Surface fluxes
!!
!! @par Description
!!          Surface flux between atmosphere and urban
!!          based on Urban Canopy Model (Kusaka et al. 2000)
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
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
  public :: CPL_AtmUrb_bulk

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: canopy_wind
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
  real(RP), private :: ZR         = 10.0       ! roof level ( building height) [m]
  real(RP), private :: roof_width =  9.0       ! roof level ( building height) [m]
  real(RP), private :: road_width = 11.0       ! roof level ( building height) [m]
  real(RP), private :: SIGMA_ZED  =  1.0       ! Standard deviation of roof height [m]
  real(RP), private :: AH         = 17.5       ! Sensible Anthropogenic heat [W/m^2]
  real(RP), private :: ALH        = 0.0        ! Latent Anthropogenic heat [W/m^2]
  real(RP), private :: BETR       = 0.0        ! Evaporation efficiency of roof [-]
  real(RP), private :: BETB       = 0.0        !                        of building [-]
  real(RP), private :: BETG       = 0.0        !                        of ground [-]
  real(RP), private :: STRGR      = 0.0        ! rain strage on roof [-]
  real(RP), private :: STRGB      = 0.0        !             on building [-]
  real(RP), private :: STRGG      = 0.0        !             on ground [-]
  real(RP), private :: CAPR       = 1.2E6      ! heat capacity of roof [J m-3 K]
  real(RP), private :: CAPB       = 1.2E6      !  ( units converted in code
  real(RP), private :: CAPG       = 1.2E6      !            to [ cal cm{-3} deg{-1} ] )
  real(RP), private :: AKSR       = 2.28       ! thermal conductivity of roof, wall, and ground [W m-1 K]
  real(RP), private :: AKSB       = 2.28       !  ( units converted in code
  real(RP), private :: AKSG       = 2.28       !            to [ cal cm{-1} s{-1} deg{-1} ] )
  real(RP), private :: ALBR       = 0.2        ! surface albedo of roof
  real(RP), private :: ALBB       = 0.2        ! surface albedo of wall
  real(RP), private :: ALBG       = 0.2        ! surface albedo of ground
  real(RP), private :: EPSR       = 0.90       ! Surface emissivity of roof
  real(RP), private :: EPSB       = 0.90       ! Surface emissivity of wall
  real(RP), private :: EPSG       = 0.90       ! Surface emissivity of ground
  real(RP), private :: Z0R        = 0.01       ! roughness length for momentum of building roof
  real(RP), private :: Z0B        = 0.0001     ! roughness length for momentum of building wall
  real(RP), private :: Z0G        = 0.01       ! roughness length for momentum of ground
  real(RP), private :: TRLEND     = 293.00     ! lower boundary condition of roof temperature [K]
  real(RP), private :: TBLEND     = 293.00     ! lower boundary condition of wall temperature [K]
  real(RP), private :: TGLEND     = 293.00     ! lower boundary condition of ground temperature [K]
  integer , private :: BOUND      = 1          ! Boundary Condition for Roof, Wall, Ground Layer Temp
                                                     !       [1: Zero-Flux, 2: T = Constant]
  ! calculate in subroutine urban_param_set
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

  real(RP), private :: ahdiurnal(1:24)         ! AH diurnal profile

  real(RP), private, allocatable :: DZR(:)     ! thickness of each roof layer [m]
  real(RP), private, allocatable :: DZB(:)     ! thickness of each building layer [m]
  real(RP), private, allocatable :: DZG(:)     ! thickness of each road layer [m]
                                                     ! ( units converted in code to [cm] )
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_AtmUrb_bulk_setup( CPL_TYPE_AtmUrb )
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
    real(RP) :: URBAN_UCM_DZR
    real(RP) :: URBAN_UCM_DZB
    real(RP) :: URBAN_UCM_DZG
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
       URBAN_UCM_DZR,        &
       URBAN_UCM_DZB,        &
       URBAN_UCM_DZG,        &
       URBAN_UCM_BOUND

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
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
    URBAN_UCM_DZR          = 0.05
    URBAN_UCM_DZB          = 0.05
    URBAN_UCM_DZG          = 0.50

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

    DZR(UKS:UKE) = URBAN_UCM_DZR
    DZB(UKS:UKE) = URBAN_UCM_DZB
    DZG(UKS:UKE) = URBAN_UCM_DZG

    ! set other urban parameters
    call urban_param_setup

    return
  end subroutine CPL_AtmUrb_bulk_setup

  !-----------------------------------------------------------------------------
  subroutine CPL_AtmUrb_bulk( &
        TR,      & ! (inout)
        TB,      & ! (inout)
        TG,      & ! (inout)
        TC,      & ! (inout)
        QC,      & ! (inout)
        UC,      & ! (inout)
        TRL,     & ! (inout)
        TBL,     & ! (inout)
        TGL,     & ! (inout)
        RAINR,   & ! (inout)
        RAINB,   & ! (inout)
        RAING,   & ! (inout)
        ROFF,    & ! (inout)
        TS,      & ! (out)
        SHR,     & ! (out)
        SHB,     & ! (out)
        SHG,     & ! (out)
        LHR,     & ! (out)
        LHB,     & ! (out)
        LHG,     & ! (out)
        GHR,     & ! (out)
        GHB,     & ! (out)
        GHG,     & ! (out)
        RNR,     & ! (out)
        RNB,     & ! (out)
        RNG,     & ! (out)
        RTS,     & ! (out)
        RN,      & ! (out)
        SH,      & ! (out)
        LH,      & ! (out)
        GHFLX,   & ! (out)
        LSOLAR,  & ! (in)
        TA,      & ! (in)
        QA,      & ! (in)
        UA,      & ! (in)
        U1,      & ! (in)
        V1,      & ! (in)
        ZA,      & ! (in)
        SSG,     & ! (in)
        LLG,     & ! (in)
        RAIN,    & ! (in)
        RHOO,    & ! (in)
        XLON,    & ! (in)
        XLAT     ) ! (in)
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       KARMAN => CONST_KARMAN,  &    ! AK : kalman constant  [-]
       PI     => CONST_PI,      &    ! PI : pi               [-]
       CPdry  => CONST_CPdry,   &    ! CPP : heat capacity of dry air [J/K/kg]
       LH0    => CONST_LH0,     &    ! ELL : latent heat of vaporization [J/kg]
       GRAV   => CONST_GRAV,    &    !< gravitational constant [m/s2]
       Rdry   => CONST_Rdry,    &    !< specific gas constant (dry) [J/kg/K]
       Rvap   => CONST_Rvap,    &    !< gas constant (water vapor) [J/kg/K]
       STB    => CONST_STB,     &    !< stefan-Boltzman constant [MKS unit]
       TEM00  => CONST_TEM00         !< temperature reference (0 degree C) [K]
    use scale_time, only:       &
       TIME => TIME_NOWSEC,      &   !< absolute sec
       DELT => TIME_DTSEC_URBAN      !< time interval of urban step [sec]
    implicit none

    !-- configuration variables
    logical, intent(in) :: LSOLAR  ! logical   [true=both, false=SSG only]

    !-- Input variables from Coupler to Urban
    real(RP), intent(in)    :: TA   ! temp at 1st atmospheric level          [K]
    real(RP), intent(in)    :: QA   ! specific humidity at 1st atmospheric level  [kg/kg]
    real(RP), intent(in)    :: UA   ! wind speed at 1st atmospheric level    [m/s]
    real(RP), intent(in)    :: U1   ! u at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: V1   ! v at 1st atmospheric level             [m/s]
    real(RP), intent(in)    :: ZA   ! height of 1st atmospheric level        [m]
    real(RP), intent(in)    :: SSG  ! downward total short wave radiation    [W/m/m]
    real(RP), intent(in)    :: LLG  ! downward long wave radiation           [W/m/m]
    real(RP), intent(in)    :: RAIN ! precipitation                          [kg/m2/s]
                                    !   (if you use RAIN, check unit!)
    real(RP), intent(in)    :: RHOO ! air density                            [kg/m^3]
    real(RP), intent(in)    :: XLAT !< latitude  [rad,-pi,pi]
    real(RP), intent(in)    :: XLON !< longitude [rad,0-2pi]

    !! for shadow effect model
    !!   real(RP), intent(in)    :: DECLIN ! solar declination
    ![rad]
    !!   real(RP), intent(in)    :: COSZ !
    !sin(fai)*sin(del)+cos(fai)*cos(del)*cos(omg)
    !!   real(RP), intent(in)    :: OMG  ! solar hour angle
    ![rad]

    !-- In/Out variables from/to Coupler to/from Urban
    real(RP), intent(inout) :: TR   ! roof temperature              [K]
    real(RP), intent(inout) :: TB   ! building wall temperature     [K]
    real(RP), intent(inout) :: TG   ! road temperature              [K]
    real(RP), intent(inout) :: TC   ! urban-canopy air temperature  [K]
    real(RP), intent(inout) :: QC   ! urban-canopy air mixing ratio [kg/kg]
    real(RP), intent(inout) :: UC   ! diagnostic canopy wind        [m/s]
    real(RP), intent(inout) :: TRL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: TBL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: TGL(UKS:UKE)  ! layer temperature [K]
    real(RP), intent(inout) :: RAINR ! rain amount in storage on roof     [kg/m2]
    real(RP), intent(inout) :: RAINB ! rain amount in storage on building [kg/m2]
    real(RP), intent(inout) :: RAING ! rain amount in storage on road     [kg/m2]
    real(RP), intent(inout) :: ROFF  ! runoff from urban        [kg/m2]

    !-- Output variables from Urban to Coupler
    real(RP), intent(out)   :: TS     ! Diagnostic surface temperature   [K]
    real(RP), intent(out)   :: RTS    ! radiative surface temperature    [K]
    real(RP), intent(out)   :: SH     ! sensible heat flux               [W/m/m]
    real(RP), intent(out)   :: LH     ! latent heat flux                 [W/m/m]
    real(RP), intent(out)   :: GHFLX  ! heat flux into the ground        [W/m/m]
    real(RP), intent(out)   :: RN     ! net radition                [W/m/m]
    real(RP), intent(out)   :: RNR, RNB, RNG
    real(RP), intent(out)   :: SHR, SHB, SHG
    real(RP), intent(out)   :: LHR, LHB, LHG
    real(RP), intent(out)   :: GHR, GHB, GHG

    !-- parameters
    real(RP), parameter    :: CP  = 0.24_RP     ! heat capacity of dry air  [cgs unit]
    real(RP), parameter    :: EL  = 583.0_RP    ! latent heat of vaporation [cgs unit]
    real(RP), parameter    :: SIG = 8.17E-11_RP ! stefun bolzman constant   [cgs unit]
    real(RP), parameter    :: SRATIO = 0.75_RP  ! ratio between direct/total solar [-]
    real(RP), parameter    :: TFa      = 0.5_RP      ! factor a in Tomita (2009) 
    real(RP), parameter    :: TFb      = 1.1_RP      ! factor b in Tomita (2009) 
    real(RP), parameter    :: redf_min = 1.0E-2_RP   ! minimum reduced factor
    real(RP), parameter    :: redf_max = 1.0_RP      ! maximum reduced factor

    !-- Local variables
    logical  :: SHADOW = .false.
           ! [true=consider svf and shadow effects, false=consider svf effect
           ! only]

    real(RP) :: LON,LAT  ! longitude [deg], latitude [deg]
    integer  :: tloc     ! local time (1-24h)
    real(RP) :: dsec     ! second [s]
    real(RP) :: tahdiurnal ! temporal AH diurnal profile
    real(RP) :: AH_t     ! Sensible Anthropogenic heat [W/m^2]
    real(RP) :: ALH_t    ! Latent Anthropogenic heat [W/m^2]

    real(RP) :: SSGD     ! downward direct short wave radiation   [W/m/m]
    real(RP) :: SSGQ     ! downward diffuse short wave radiation  [W/m/m]

    real(RP) :: W, VFGS, VFGW, VFWG, VFWS, VFWW
    real(RP) :: SX, SD, SQ, RX
    real(RP) :: RHO

    real(RP) :: TRP = 350.0_RP   ! TRP: at previous time step [K]
    real(RP) :: TBP = 350.0_RP   ! TBP: at previous time step [K]
    real(RP) :: TGP = 350.0_RP   ! TGP: at previous time step [K]
    real(RP) :: TCP = 350.0_RP   ! TCP: at previous time step [K]
    real(RP) :: QCP = 0.01_RP    ! QCP: at previous time step [kg/kg]
    real(RP) :: TSP = 300.0_RP   ! TSP: at previous time step [K]
    real(RP) :: UST, TST, QST

    real(RP) :: RAINT
    real(RP) :: ROFFR, ROFFB, ROFFG ! runoff [kg/m2]

    real(RP) :: PS         ! Surface Pressure [hPa]
    real(RP) :: TAV        ! Vertial Temperature [K]
    real(RP) :: qmix       ! mixing ratio [kg/kg]
    real(RP) :: ES

    real(RP) :: LUP, LDN, RUP, ALBD_LW_grid
    real(RP) :: SUP, SDN, ALBD_SW_grid

    real(RP) :: LNET, SNET, FLXUV, FLXTH, FLXHUM, FLXG
    real(RP) :: SW     ! shortwave radition           [W/m/m]
    real(RP) :: LW     ! longwave radition           [W/m/m]
    real(RP) :: PSIM   ! similality stability shear function for momentum
    real(RP) :: PSIH   ! similality stability shear function for heat

    ! for shadow effect model
    ! real(RP) :: HOUI1, HOUI2, HOUI3, HOUI4, HOUI5, HOUI6, HOUI7, HOUI8
    ! real(RP) :: SLX, SLX1, SLX2, SLX3, SLX4, SLX5, SLX6, SLX7, SLX8
    ! real(RP) :: THEATAZ    ! Solar Zenith Angle [rad]
    ! real(RP) :: THEATAS    ! = PI/2. - THETAZ
    ! real(RP) :: FAI        ! Latitude [rad]

    real(RP) :: SR, SB, SG, RR, RB, RG
    real(RP) :: SR1, SB1, SB2, SG1, SG2
    real(RP) :: RB1, RB2, RG1, RG2
    real(RP) :: HR, ELER, G0R, FLXTHR, FLXHUMR
    real(RP) :: HB, ELEB, G0B, FLXTHB, FLXHUMB
    real(RP) :: HG, ELEG, G0G, FLXTHG, FLXHUMG

    real(RP) :: Z
    real(RP) :: QS0R, DQS0RDTR, DESDT
    real(RP) :: QS0B, DQS0BDTB
    real(RP) :: QS0G, DQS0GDTG

    real(RP) :: RIBR, XXXR, BHR, ALPHAR, CHR, CDR
    real(RP) :: RIBC, XXXC, BHC
    real(RP) :: ALPHAC, ALPHAB, ALPHAG
    real(RP) :: CHC, CHB, CHG, CDC
    real(RP) :: TC1, TC2, QC1, QC2

    real(RP) :: oldF,oldGF      ! residual in previous step
    real(RP) :: redf,redfg      ! reduced factor

    real(RP) :: F
    real(RP) :: DRRDTR, DHRDTR, DELERDTR, DG0RDTR
    real(RP) :: DTR, DFDT

    real(RP) :: FX, FY, GF, GX, GY
    real(RP) :: DTCDTB, DTCDTG
    real(RP) :: DQCDTB, DQCDTG
    real(RP) :: DRBDTB1,  DRBDTG1,  DRBDTB2,  DRBDTG2
    real(RP) :: DRGDTB1,  DRGDTG1,  DRGDTB2,  DRGDTG2
    real(RP) :: DRBDTB,   DRBDTG,   DRGDTB,   DRGDTG
    real(RP) :: DHBDTB,   DHBDTG,   DHGDTB,   DHGDTG
    real(RP) :: DELEBDTB, DELEBDTG, DELEGDTG, DELEGDTB
    real(RP) :: DG0BDTB,  DG0BDTG,  DG0GDTG,  DG0GDTB
    real(RP) :: DTB, DTG, DTC

    real(RP) :: XXX, X, Z0, Z0H, CD, CH, CHU
    real(RP) :: PSIX, PSIT, PSIX2, PSIT2, PSIX10, PSIT10

    integer  :: iteration, K

    !-----------------------------------------------------------
    ! Set parameters
    !-----------------------------------------------------------

    ! local time
    ! tloc=mod(int(OMG/PI*180./15.+12.+0.5 ),24)
    LAT=XLAT/PI*180.0_RP
    LON=XLON/PI*180.0_RP

    tloc = mod( (int(TIME/(60.0_RP*60.0_RP)) + int(LON/15.0_RP)),24 )
    dsec = mod(TIME,3600.0_RP)
    if(tloc==0) tloc = 24
    ! print *,'tloc',TIME,tloc,XLON,XLAT,LON,LAT

    ! Calculate AH data at LST
     if(tloc==24)then
       tahdiurnal = (ahdiurnal(tloc)*(3600.0_RP-dsec)+ahdiurnal(1)*dsec)/3600.0_RP
     else
       tahdiurnal = (ahdiurnal(tloc)*(3600.0_RP-dsec)+ahdiurnal(tloc+1)*dsec)/3600.0_RP
     endif
     AH_t  = AH  * tahdiurnal
     ALH_t = ALH * tahdiurnal
    ! AH_t  = AH  * ahdiurnal(tloc)
    ! ALH_t = ALH * ahdiurnal(tloc)


    !if(dsec==0.) print *,'tloc',tloc,SSG,RHOO,QA,TA

    if( ZDC+Z0C+2.0_RP >= ZA) then
      if( IO_L ) write(IO_FID_LOG,*) 'ZDC + Z0C + 2m is larger than the 1st WRF level' // &
                                     'Stop in subroutine urban - change ZDC and Z0C'
      call PRC_MPIstop
    endif

    !if(.NOT.LSOLAR) then   ! Radiation scheme does not have SSGD and SSGQ.
      SSGD = SRATIO * SSG
      SSGQ = SSG - SSGD
    !endif

    W    = 2.0_RP * 1.0_RP * HGT
    VFGS = SVF
    VFGW = 1.0_RP - SVF
    VFWG = ( 1.0_RP - SVF ) * ( 1.0_RP - R ) / W
    VFWS = VFWG
    VFWW = 1.0_RP - 2.0_RP * VFWG

    !--- Convert unit from MKS to cgs

    SX  = (SSGD+SSGQ) / 697.7_RP / 60.0_RP  ! downward short wave radition [ly/min]
    SD  = SSGD        / 697.7_RP / 60.0_RP  ! downward direct short wave radiation
    SQ  = SSGQ        / 697.7_RP / 60.0_RP  ! downward diffuse short wave radiation
    RX  = LLG         / 697.7_RP / 60.0_RP  ! downward long wave radiation
    RHO = RHOO * 0.001_RP                   ! air density at first atmospheric level

    !--- Renew surface and layer temperatures

    TRP = TR
    TBP = TB
    TGP = TG
    TCP = TC
    QCP = QC
    TSP = TS

    !--- calculate canopy wind

    call canopy_wind(ZA, UA, UC)

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
    ! Set evaporation efficiency on roof/wall/road
    !-----------------------------------------------------------

    RAINT = RAIN * DELT  !!! check [kg/m2/s -> kg/m2 ?]

    call cal_beta(BETR, RAINT, RAINR, STRGR, ROFFR)
    call cal_beta(BETB, RAINT, RAINB, STRGB, ROFFB)
    call cal_beta(BETG, RAINT, RAING, STRGG, ROFFG)

    ROFF = ROFF + (R-0.05_RP)*ROFFR     &
                + 0.1_RP*ROFFB          &
                + (RW-0.05_RP)*ROFFG

    !print *, TIME, BETR, BETB, BETG, ROFF

    !-----------------------------------------------------------
    ! Energy balance on roof/wall/road surface
    !-----------------------------------------------------------

    !--- Roof

    Z    = ZA - ZDC
    BHR  = LOG(Z0R/Z0HR) / 0.4_RP
    RIBR = ( GRAV * 2.0_RP / (TA+TRP) ) * (TA-TRP) * (Z+Z0R) / (UA*UA)

    call mos(XXXR,ALPHAR,CDR,BHR,RIBR,Z,Z0R,UA,TA,TRP,RHO)

    CHR = ALPHAR / RHO / CP / UA

    qmix = QA / ( 1.0_RP - QA )
    TAV = TA * ( 1.0_RP + 0.61_RP * qmix )
    PS  = RHOO * Rdry * TAV / 100.0_RP ! [hPa]

    ! TR  Solving Non-Linear Equation by Newton-Rapson

    ! TRP = 350.0_RP

    redf   = 1.0_RP
    oldF   = 1.0E+5_RP

    do iteration = 1, 20

      ES       = 6.11_RP * exp( ( LH0/Rvap) * (TRP-TEM00) / (TEM00*TRP) )
      DESDT    = (LH0/Rvap) * ES / (TRP**2)
      QS0R     = 0.622_RP * ES / ( PS - 0.378_RP * ES )
      DQS0RDTR = DESDT * 0.622_RP * PS / ( ( PS - 0.378_RP * ES )**2 )

      RR       = EPSR * ( RX - SIG * (TRP**4) / 60.0_RP )
      HR       = RHO * CP * CHR * UA * (TRP-TA) * 100.0_RP
      ELER     = RHO * EL * CHR * UA * BETR * (QS0R-QA) * 100.0_RP
      G0R      = AKSR * (TRP-TRL(1)) / (DZR(1)/2.0_RP)

      F        = SR + RR - HR - ELER - G0R

      DRRDTR   = ( -4.0_RP * EPSR * SIG * TRP**3 ) / 60.0_RP
      DHRDTR   = RHO * CP * CHR * UA * 100.0_RP
      DELERDTR = RHO * EL * CHR * UA * BETR * DQS0RDTR * 100.0_RP
      DG0RDTR  = 2.0_RP * AKSR / DZR(1)

      DFDT     = DRRDTR - DHRDTR - DELERDTR - DG0RDTR

      DTR      = -1.0_RP * F / DFDT

      !DTR      = F / DFDT
      !TR       = TRP - DTR
      !TRP      = TR

      !!! Tomita(2009)
      if( redf < 0.0_RP ) then
      	  redf = 1.0_RP
      end if
      if( abs(F) > abs(oldF) ) then
          redf = max( TFa*redf, redf_min )
      else
	  redf = min( TFb*redf, redf_max )
      end if
      if( DFDT > 0.0_RP ) then
      	  redf = -1.0_RP
      end if
      TR       = TRP + redf * DTR
      TRP      = TR
      oldF     = F      

      if( abs(F) < 0.000001_RP .AND. abs(DTR) < 0.000001_RP ) exit

    end do

    !print *,'roof',tloc,SR,RR,HR,ELER,G0R
    !print *,'HR  ',TIME,HR,RHO,CHR,UA,(TRP-TA)
    !print *,'ELER',TIME,ELER,RHO,CHR,UA,BETR,(QS0R-QA)
    !print *,'G0R ',TIME,G0R,AKSR,(TRP-TRL(1)),DZR(1)
    !print *,'F   ',iteration-1,abs(F),abs(DTR)

    FLXTHR  = HR / RHO / CP / 100.0_RP
    FLXHUMR = ELER / RHO / EL / 100.0_RP

    !!--- calculate the rain amount remaining on the surface
    RAINR = max(0.0_RP, RAINR-(ELER/EL)*DELT*10.0_RP)   ! from cgs to MKS [kg/m/m = mm]

    !--- Wall and Road

    Z    = ZA - ZDC
    BHC  = LOG(Z0C/Z0HC) / 0.4_RP
    RIBC = ( GRAV * 2.0_RP / (TA+TCP) ) * (TA-TCP) * (Z+Z0C) / (UA*UA)

    call mos(XXXC,ALPHAC,CDC,BHC,RIBC,Z,Z0C,UA,TA,TCP,RHO)

    ! empirical form
    ALPHAB = RHO * CP * ( 6.15_RP + 4.18_RP * UC ) / 1200.0_RP
    if( UC > 5.0_RP ) ALPHAB = RHO * CP * ( 7.51_RP * UC**0.78_RP ) / 1200.0_RP

    ALPHAG = RHO * CP * ( 6.15_RP + 4.18_RP * UC ) / 1200.0_RP
    if( UC > 5.0_RP ) ALPHAG = RHO * CP * ( 7.51_RP * UC**0.78_RP ) / 1200.0_RP

    CHC = ALPHAC / RHO / CP / UA
    CHB = ALPHAB / RHO / CP / UC
    CHG = ALPHAG / RHO / CP / UC


    ! TB,TG  Solving Non-Linear Simultaneous Equation by Newton-Rapson
    ! TBP = 350.0_RP
    ! TGP = 350.0_RP

    do iteration = 1, 20

      ES       = 6.11_RP * exp( (LH0/Rvap) * (TBP-TEM00) / (TEM00*TBP) )
      DESDT    = (LH0/Rvap) * ES / (TBP**2)
      QS0B     = 0.622_RP * ES / ( PS - 0.378_RP * ES )
      DQS0BDTB = DESDT * 0.622_RP * PS / ( ( PS - 0.378_RP * ES )**2 )

      ES       = 6.11_RP * exp( (LH0/Rvap) * (TGP-TEM00) / (TEM00*TGP) )
      DESDT    = (LH0/Rvap) * ES / (TGP**2)
      QS0G     = 0.622_RP * ES / ( PS - 0.378_RP * ES )
      DQS0GDTG = DESDT * 0.622_RP * PS / ( ( PS - 0.378_RP * ES )**2 )

      RG1      = EPSG * ( RX * VFGS                             &
                        + EPSB * VFGW * SIG * TBP**4 / 60.0_RP  &
                        - SIG * TGP**4 / 60.0_RP                )

      RB1      = EPSB * ( RX * VFWS         &
                        + EPSG * VFWG * SIG * TGP**4 / 60.0_RP &
                        + EPSB * VFWW * SIG * TBP**4 / 60.0_RP &
                        - SIG * TBP**4 / 60.0_RP               )

      RG2      = EPSG * ( (1.0_RP-EPSB) * VFGW * VFWS * RX                                            &
                        + (1.0_RP-EPSB) * VFGW * VFWG * EPSG * SIG * TGP**4 /60.0_RP                  &
                        + EPSB * (1.0_RP-EPSB) * VFGW * VFWW * SIG * TBP**4 / 60.0_RP )

!      RB2      = EPSB * ( (1.0_RP-EPSG) * VFWG * VFGS * RX                                                            &
!                        + (1.0_RP-EPSG) * EPSB * VFGW * VFWG * SIG * (TBP**4) / 60.0_RP                               &
!                        + (1.0_RP-EPSB) * VFWS * (1.0_RP-2.0_RP*VFWS) * RX                                            &
!!                       + (1.0_RP-EPSB) * VFWG * (1.0_RP-2.0_RP*VFWS) * EPSG * SIG * EPSG * TGP**4 / 60.0_RP          &
!                        + EPSB * (1.0_RP-EPSB) * (1.0_RP-2.0_RP*VFWS) * (1.0_RP-2.0_RP*VFWS) * SIG * TBP**4 / 60.0_RP )

      RB2      = EPSB * ( (1.0_RP-EPSG) * VFWG * VFGS * RX                                                            &
                        + (1.0_RP-EPSG) * EPSB * VFGW * VFWG * SIG * (TBP**4) / 60.0_RP                               &
                        + (1.0_RP-EPSB) * VFWS * VFWW * RX                                                            &
                        + (1.0_RP-EPSB) * VFWG * VFWW * SIG * EPSG * TGP**4 / 60.0_RP                                 &
                        + EPSB * (1.0_RP-EPSB) * VFWW * (1.0_RP-2.0_RP*VFWS) * SIG * TBP**4 / 60.0_RP )

      RG       = RG1 + RG2
      RB       = RB1 + RB2


      DRBDTB1  = EPSB * ( 4.0_RP * EPSB * SIG * TBP**3 * VFWW - 4.0_RP * SIG * TBP**3 ) / 60.0_RP
      DRBDTG1  = EPSB * ( 4.0_RP * EPSG * SIG * TGP**3 * VFWG ) / 60.0_RP
      DRBDTB2  = EPSB * ( 4.0_RP * (1.0_RP-EPSG) * EPSB * SIG * TBP**3 * VFGW * VFWG &
                        + 4.0_RP * EPSB * (1.0_RP-EPSB) * SIG * TBP**3 * VFWW * VFWW ) / 60.0_RP
      DRBDTG2  = EPSB * ( 4.0_RP * (1.0_RP-EPSB) * EPSG * SIG * TGP**3 * VFWG * VFWW ) / 60.0_RP

      DRGDTB1  = EPSG * ( 4.0_RP * EPSB * SIG * TBP**3 * VFGW ) / 60.0_RP
      DRGDTG1  = EPSG * ( -4.0_RP * SIG * TGP**3 ) / 60.0_RP
      DRGDTB2  = EPSG * ( 4.0_RP * EPSB * (1.0_RP-EPSB) * SIG * TBP**3 * VFWW * VFGW ) / 60.0_RP
      DRGDTG2  = EPSG * ( 4.0_RP * (1.0_RP-EPSB) * EPSG * SIG * TGP**3 * VFWG * VFGW ) / 60.0_RP

      DRBDTB   = DRBDTB1 + DRBDTB2
      DRBDTG   = DRBDTG1 + DRBDTG2
      DRGDTB   = DRGDTB1 + DRGDTB2
      DRGDTG   = DRGDTG1 + DRGDTG2

      HB       = RHO * CP * CHB * UC * (TBP-TCP) * 100.0_RP
      HG       = RHO * CP * CHG * UC * (TGP-TCP) * 100.0_RP

      DTCDTB   = W  * ALPHAB / ( RW*ALPHAC + RW*ALPHAG + W*ALPHAB )
      DTCDTG   = RW * ALPHAG / ( RW*ALPHAC + RW*ALPHAG + W*ALPHAB )

      DHBDTB   = RHO * CP * CHB * UC * (1.0_RP-DTCDTB) * 100.0_RP
      DHBDTG   = RHO * CP * CHB * UC * (0.0_RP-DTCDTG) * 100.0_RP
      DHGDTG   = RHO * CP * CHG * UC * (1.0_RP-DTCDTG) * 100.0_RP
      DHGDTB   = RHO * CP * CHG * UC * (0.0_RP-DTCDTB) * 100.0_RP

      ELEB     = RHO * EL * CHB * UC * BETB * (QS0B-QCP) * 100.0_RP
      ELEG     = RHO * EL * CHG * UC * BETG * (QS0G-QCP) * 100.0_RP

      DQCDTB   = W  * ALPHAB * BETB * DQS0BDTB / ( RW*ALPHAC + RW*ALPHAG*BETG + W*ALPHAB*BETB )
      DQCDTG   = RW * ALPHAG * BETG * DQS0GDTG / ( RW*ALPHAC + RW*ALPHAG*BETG + W*ALPHAB*BETB )

      DELEBDTB = RHO * EL * CHB * UC * BETB * (DQS0BDTB-DQCDTB) * 100.0_RP
      DELEBDTG = RHO * EL * CHB * UC * BETB * (0.0_RP-DQCDTG) * 100.0_RP
      DELEGDTG = RHO * EL * CHG * UC * BETG * (DQS0GDTG-DQCDTG) * 100.0_RP
      DELEGDTB = RHO * EL * CHG * UC * BETG * (0.0_RP-DQCDTB) * 100.0_RP

      G0B      = AKSB * (TBP-TBL(1)) / (DZB(1)/2.0_RP)
      G0G      = AKSG * (TGP-TGL(1)) / (DZG(1)/2.0_RP)

      DG0BDTB  = 2.0_RP * AKSB / DZB(1)
      DG0BDTG  = 0.0_RP
      DG0GDTG  = 2.0_RP * AKSG / DZG(1)
      DG0GDTB  = 0.0_RP

      F        = SB + RB - HB - ELEB - G0B
      FX       = DRBDTB - DHBDTB - DELEBDTB - DG0BDTB
      FY       = DRBDTG - DHBDTG - DELEBDTG - DG0BDTG

      GF       = SG + RG - HG - ELEG - G0G
      GX       = DRGDTB - DHGDTB - DELEGDTB - DG0GDTB
      GY       = DRGDTG - DHGDTG - DELEGDTG - DG0GDTG

      DTB      = ( GF*FY - F*GY ) / ( FX*GY - GX*FY )
      DTG      = -( GF + GX*DTB ) / GY

      TB       = TBP + DTB
      TG       = TGP + DTG

      TBP      = TB
      TGP      = TG

      TC1      = RW*ALPHAC + RW*ALPHAG + W*ALPHAB
      TC2      = RW*ALPHAC*TA + RW*ALPHAG*TGP + W*ALPHAB*TBP
      TC       = TC2 / TC1

      QC1      = RW*ALPHAC + RW*ALPHAG*BETG + W*ALPHAB*BETB
      QC2      = RW*ALPHAC*QA + RW*ALPHAG*BETG*QS0G + W*ALPHAB*BETB*QS0B
      QC       = QC2 / QC1

      DTC      = TCP - TC
      TCP      = TC
      QCP      = QC

      if( abs(F)   < 0.000001_RP .and. &
          abs(DTB) < 0.000001_RP .and. &
          abs(GF)  < 0.000001_RP .and. &
          abs(DTG) < 0.000001_RP .and. &
          abs(DTC) < 0.000001_RP       ) exit

    end do

    FLXTHB  = HB / RHO / CP / 100.0_RP
    FLXHUMB = ELEB / RHO / EL / 100.0_RP
    FLXTHG  = HG / RHO / CP / 100.0_RP
    FLXHUMG = ELEG / RHO / EL / 100.0_RP

    !!--- calculate the rain amount remaining on the surface
    RAINB = max(0.0_RP, RAINB-(ELEB/EL)*DELT*10.0_RP)   ! from cgs to MKS [kg/m/m = mm]

    !!--- calculate the rain amount remaining on the surface
    RAING = max(0.0_RP, RAING-(ELEG/EL)*DELT*10.0_RP)   ! from cgs to MKS [kg/m/m = mm]


    !-----------------------------------------------------------
    ! Total Fluxes from Urban Canopy
    !-----------------------------------------------------------

    FLXUV  = ( R*CDR + RW*CDC ) * UA * UA
!    FLXTH  = ( R*FLXTHR  + W*FLXTHB  + RW*FLXTHG ) + AH_t/RHOO/CPdry
!    FLXHUM = ( R*FLXHUMR + W*FLXHUMB + RW*FLXHUMG ) + ALH_t/RHOO/LH0
    FLXTH  = ( R*FLXTHR  + W*FLXTHB  + RW*FLXTHG )
    FLXHUM = ( R*FLXHUMR + W*FLXHUMB + RW*FLXHUMG )
    FLXG   = ( R*G0R + W*G0B + RW*G0G )
    LNET   = R*RR + W*RB + RW*RG

!    print *,'SNET  ',SNET
!    print *,'LNET  ',LNET
!    print *,'FLXTH ',FLXTH*RHO*CP*100.0_RP
!    print *,'FLXHUM',FLXHUM*RHO*EL*100.0
!    print *,'FLXG  ',FLXG
!    print *, SNET+LNET-FLXTH*RHO*CP*100.0_RP-FLXHUM*RHO*EL*100.0_RP-FLXG

    !-----------------------------------------------------------
    ! Grid average
    !-----------------------------------------------------------
    !--- shortwave radiation
    SDN =  R + W * (VFWS + VFGS * ALBG * VFWG)  + RW * (VFGS + VFWS * ALBB * VFGW)
    SUP =  R * ALBR  &
         + W * ( VFWS* ALBB + VFGS * ALBG * VFWG *ALBB ) &
         + RW * ( VFGS * ALBG + VFWS * ALBB * VFGW *ALBG )

    ALBD_SW_grid = SUP / SDN

    !--- longwave radiation
    LDN = R + W*VFWS + RW*VFGS
    LUP =  R * (1.0_RP-EPSR) &
         + W*( (1.0_RP-EPSB*VFWW)*(1.0_RP-EPSB)*VFWS - EPSB*VFWG*(1.0_RP-EPSG)*VFGS )  &
         + RW*( (1.0_RP-EPSG)*VFGS - EPSG*(1.0_RP-VFGS)*(1.0_RP-EPSB)*VFWS )

    RUP = (LDN - LUP) * RX - LNET
    ALBD_LW_grid = LUP / LDN

!
!    RUP =  R * (EPSR * SIG * (TR**4) / 60.0_RP) &
!          + W * (EPSB*SIG*(TB**4)/60.0_RP - EPSB*EPSG*VFWG*SIG*(TG**4)/60.0_RP - EPSB*EPSB*VFWW*SIG*(TB**4)/60.0_RP   &
!                 - EPSB *(1.0_RP-EPSG) * EPSB * VFGW * VFWG * SIG * (TB**4) / 60.0_RP                                 &
!                 - EPSB *(1.0_RP-EPSB) * VFWG * VFWW * SIG * EPSG * (TG**4) / 60.0_RP                                 &
!                 - EPSB * EPSB * (1.0_RP-EPSB) * VFWW * VFWW * SIG * (TB**4) / 60.0_RP )                              &
!          + RW * (EPSG*SIG*(TG**4)/60.0_RP - EPSG * EPSB * VFGW * SIG * (TB**4) / 60.0_RP                             &
!                 - EPSG * (1.0_RP-EPSB) * (1.0_RP-SVF) * VFWG * EPSG * SIG * TG**4 /60.0_RP                           &
!                 - EPSG * EPSB * (1.0_RP-EPSB) * (1.0_RP-SVF) * VFWW * SIG * TB**4 / 60.0_RP )
!
    !-----------------------------------------------------------
    ! Convert Unit: FLUXES and u* T* q*  [cgs] --> [MKS]
    !-----------------------------------------------------------

    SH = FLXTH  * RHOO * CPdry               ! Sensible heat flux          [W/m/m]
    LH = FLXHUM * RHOO * LH0                 ! Latent heat flux            [W/m/m]
    LW = LLG - ( LNET * 697.7_RP * 60.0_RP ) ! Upward longwave radiation   [W/m/m]
    SW = SSG - ( SNET * 697.7_RP * 60.0_RP ) ! Upward shortwave radiation  [W/m/m]

    GHFLX  = -FLXG * 697.7_RP * 60.0_RP          ! [W/m/m]

    RN = (SNET+LNET) * 697.7_RP * 60.0_RP    ! Net radiation [W/m/m]

    RUP = RUP * 697.7_RP * 60.0_RP        ! Upward longwave = sig*T**4 [W/m/m]

    SHR = FLXTHR * RHOO * CPdry           ! Sensible heat flux on roof [W/m/m]
    SHB = FLXTHB * RHOO * CPdry           ! Sensible heat flux on wall [W/m/m]
    SHG = FLXTHG * RHOO * CPdry           ! Sensible heat flux on road [W/m/m]
    LHR = FLXHUMR * RHOO * LH0           ! Latent heat flux on road [W/m/m]
    LHB = FLXHUMB * RHOO * LH0           ! Latent heat flux on wall [W/m/m]
    LHG = FLXHUMG * RHOO * LH0           ! Latent heat flux on road [W/m/m]
    GHR = -1.0_RP * G0R * 697.7_RP * 60.0_RP          ! Ground heat flux on roof [W/m/m]
    GHB = -1.0_RP * G0B * 697.7_RP * 60.0_RP          ! Ground heat flux on wall [W/m/m]
    GHG = -1.0_RP * G0G * 697.7_RP * 60.0_RP          ! Ground heat flux on road [W/m/m]
    RNR = (SR + RR) * 697.7_RP * 60.0_RP     ! Net radiation on roof [W/m/m]
    RNB = (SB + RB) * 697.7_RP * 60.0_RP     ! Net radiation on building [W/m/m]
    RNG = (SG + RG) * 697.7_RP * 60.0_RP     ! Net radiation on ground [W/m/m]

    UST = sqrt( FLXUV )                      ! u* [m/s]
    TST = -FLXTH / UST                       ! T* [K]
    QST = -FLXHUM / UST                      ! q* [-]

    !-----------------------------------------------------------
    !  diagnostic GRID AVERAGED TS from heat flux
    !-----------------------------------------------------------

    Z    = ZA - ZDC
    BHC  = LOG(Z0C/Z0HC) / 0.4_RP
    RIBC = ( GRAV * 2.0_RP / (TA+TSP) ) * (TA-TSP) * (Z+Z0C) / (UA*UA)

    call mos(XXX,ALPHAC,CD,BHC,RIBC,Z,Z0C,UA,TA,TS,RHO)
    CHU = (ALPHAC / CP / RHO) * (CPdry * RHOO)  ! [cgs -> MKS]
    TS = TA + SH / CHU    ! surface potential temp (flux temp)

!    CH = ALPHAC / RHO / CP / UA
!    print *,'total    ',TS, CH,SH/CPdry/RHOO*CP*RHO*100.,(SH/CPdry/RHOO*CP*RHO*100.)/RHO/CP/CH/UA/100.0_RP,TA
!    print *,'roof     ',TR,CHR,HR,HR/RHO/CP/CHR/UA/100.0_RP,TA
!    print *,'building ',TB,CHB,HB,HB/RHO/CP/CHB/UC/100.0_RP,TC
!    print *,'ground   ',TG,CHG,HG,HG/RHO/CP/CHG/UC/100.0_RP,TC

    !-----------------------------------------------------------
    !  diagnostic GRID AVERAGED TS from upward logwave
    !-----------------------------------------------------------

    RTS = ( RUP / STB / ( 1.0_RP-ALBD_LW_grid) )**0.25

    !-----------------------------------------------------------
    ! add anthropogenic heat fluxes
    !-----------------------------------------------------------

    FLXTH  = FLXTH + AH_t/RHOO/CPdry
    FLXHUM = FLXHUM + ALH_t/RHOO/LH0
    SH = FLXTH  * RHOO * CPdry               ! Sensible heat flux          [W/m/m]
    LH = FLXHUM * RHOO * LH0                 ! Latent heat flux            [W/m/m]

    !-----------------------------------------------------------
    !  calculate temperature in building/road
    !  multi-layer heat equation model
    !  Solving Heat Equation by Tri Diagonal Matrix Algorithm
    !-----------------------------------------------------------

    call multi_layer(UKE,BOUND,G0R,CAPR,AKSR,TRL,DZR,DELT,TRLEND)
    call multi_layer(UKE,BOUND,G0B,CAPB,AKSB,TBL,DZB,DELT,TBLEND)
    call multi_layer(UKE,BOUND,G0G,CAPG,AKSG,TGL,DZG,DELT,TGLEND)

!-----
!    Z0  = Z0C
!    Z0H = Z0HC
!    Z   = ZA - ZDC
!
!    XXX = KARMAN * GRAV * Z * TST / TA / UST / UST ! M-O theory (z/L)
!
!    if( XXX >=  1.0_RP ) XXX =  1.0_RP
!    if( XXX <= -5.0_RP ) XXX = -5.0_RP
!
!    if( XXX > 0.0_RP ) then
!      PSIM = -5.0_RP * XXX
!      PSIH = -5.0_RP * XXX
!    else
!      X    = ( 1.0_RP - 16.0_RP * XXX )**0.25_RP
!      PSIM = 2.0_RP * log((1.0_RP+X)/2.0_RP) + log((1.0_RP+X*X)/2.0_RP) - 2.0_RP*atan(X) + PI/2.0_RP
!      PSIH = 2.0_RP * log((1.0_RP+X*X)/2.0_RP)
!    end if
!
!    CD = KARMAN**2.0_RP / ( LOG(Z/Z0) - PSIM )**2
!
!    CH = 0.4_RP**2 / ( LOG(Z/Z0) - PSIM ) / ( LOG(Z/Z0H) - PSIH )
!    TS = TA + FLXTH / CH / UA   ! surface potential temp (flux temp)
!
    !   TS = TA + FLXTH/CHS    ! surface potential temp (flux temp)
    !   QS = QA + FLXHUM/CHS   ! surface humidity
    !   TS = (LW/STB/0.88)**0.25       ! Radiative temperature [K]

    return
  end subroutine CPL_AtmUrb_bulk

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

    if ( STRG == 0.0_RP ) then
       BET=0.0_RP
       ROFF  = RAIN
    else
       WATER = WATER + RAIN
       ROFF  = max(0.0_RP, WATER-STRG)
       WATER = WATER - max(0.0_RP, WATER-STRG)
       BET   = WATER / STRG
    endif

!    print *,BET, RAIN, WATER, STRG, ROFF

    return
  end subroutine cal_beta

  !-----------------------------------------------------------------------------
  subroutine mos(XXX,ALPHA,CD,B1,RIB,Z,Z0,UA,TA,TSF,RHO)

    !  XXX:   z/L (requires iteration by Newton-Rapson method)
    !  B1:    Stanton number
    !  PSIM:  = PSIX of LSM
    !  PSIH:  = PSIT of LSM

    implicit none

    real(RP), parameter     :: CP = 0.24_RP
    real(RP), intent(in)    :: B1, Z, Z0, UA, TA, TSF, RHO
    real(RP), intent(out)   :: ALPHA, CD
    real(RP), intent(inout) :: XXX, RIB
    real(RP)                :: XXX0, X, X0, FAIH, DPSIM, DPSIH
    real(RP)                :: F, DF, XXXP, US, TS, AL, XKB, DD, PSIM, PSIH
    integer                 :: NEWT
    integer, parameter      :: NEWT_END = 10

    if( RIB <= -15.0_RP ) RIB = -15.0_RP

    if( RIB < 0.0_RP ) then

      do NEWT = 1, NEWT_END

        if( XXX >= 0.0_RP ) XXX = -1.0e-3_RP

        XXX0  = XXX * Z0/(Z+Z0)

        X     = (1.0_RP-16.0_RP*XXX)**0.25
        X0    = (1.0_RP-16.0_RP*XXX0)**0.25

        PSIM  = log( (Z+Z0)/Z0 ) &
              - log( (X+1.0_RP)**2 * (X**2+1.0_RP) ) &
              + 2.0_RP * atan(X) &
              + log( (X+1.0_RP)**2 * (X0**2+1.0_RP) ) &
              - 2.0_RP * atan(X0)
        FAIH  = 1.0_RP / sqrt( 1.0_RP - 16.0_RP*XXX )
        PSIH  = log( (Z+Z0)/Z0 ) + 0.4_RP*B1 &
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
      PSIM = log( (Z+Z0)/Z0 ) + 7.0_RP * XXX
      PSIH = PSIM + 0.4_RP * B1

    else

      AL   = LOG( (Z+Z0)/Z0 )
      XKB  = 0.4_RP * B1
      DD   = -4.0_RP * RIB * 7.0_RP * XKB * AL + (AL+XKB)**2
      if( DD <= 0.0_RP ) DD = 0.0_RP

      XXX  = ( AL + XKB - 2.0_RP*RIB*7.0_RP*AL - sqrt(DD) ) / ( 2.0_RP * ( RIB*7.0_RP**2 - 7.0_RP ) )
      PSIM = log( (Z+Z0)/Z0 ) + 7.0_RP * min( XXX, 0.714_RP )
      PSIH = PSIM + 0.4_RP * B1

    endif

    US = 0.4_RP * UA / PSIM             ! u*
    if( US <= 0.01_RP ) US = 0.01_RP
    TS = 0.4_RP * (TA-TSF) / PSIH       ! T*

    CD    = US * US / UA**2                 ! CD
    ALPHA = RHO * CP * 0.4_RP * US / PSIH   ! RHO*CP*CH*U

    return
  end subroutine mos

  !-----------------------------------------------------------------------------
  subroutine multi_layer(KM,BOUND,G0,CAP,AKS,TSL,DZ,DELT,TSLEND)
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
    real(RP)                :: A(KM), B(KM), C(KM), D(KM), X(KM), P(KM), Q(KM)
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

    X(KM) = Q(KM)

    do K = KM-1, 1, -1
      X(K) = P(K) * X(K+1) + Q(K)
    end do

    do K = 1, KM
      TSL(K) = X(K)
    enddo

    return
  end subroutine multi_layer

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

    ! Thickness of roof, building wall, ground layers
    ! DZR(UKS:UKE) = 0.05 ! [m]
    ! DZB(UKS:UKE) = 0.05 ! [m]
    ! DZG(1)       = 0.05 ! [m]
    ! DZG(2)       = 0.25 ! [m]
    ! DZG(3)       = 0.50 ! [m]
    ! DZG(4)       = 0.75 ! [m]

    ! convert unit
    DZR(UKS:UKE) = DZR(UKS:UKE) * 100.0_RP ! [m]-->[cm]
    DZB(UKS:UKE) = DZB(UKS:UKE) * 100.0_RP ! [m]-->[cm]
    DZG(UKS:UKE) = DZG(UKS:UKE) * 100.0_RP ! [m]-->[cm]

    CAPR = CAPR * ( 1.0_RP / 4.1868_RP ) * 1.E-6_RP   ! [J/m^3/K] --> [cal/cm^3/deg]
    CAPB = CAPB * ( 1.0_RP / 4.1868_RP ) * 1.E-6_RP   ! [J/m^3/K] --> [cal/cm^3/deg]
    CAPG = CAPG * ( 1.0_RP / 4.1868_RP ) * 1.E-6_RP   ! [J/m^3/K] --> [cal/cm^3/deg]
    AKSR = AKSR * ( 1.0_RP / 4.1868_RP ) * 1.E-2_RP   ! [J/m/s/K] --> [cal/cm/s/deg]
    AKSB = AKSB * ( 1.0_RP / 4.1868_RP ) * 1.E-2_RP   ! [J/m/s/K] --> [cal/cm/s/deg]
    AKSG = AKSG * ( 1.0_RP / 4.1868_RP ) * 1.E-2_RP   ! [J/m/s/K] --> [cal/cm/s/deg]

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
