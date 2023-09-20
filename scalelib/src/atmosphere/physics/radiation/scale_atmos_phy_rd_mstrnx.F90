!-------------------------------------------------------------------------------
!> module atmosphere / physics / radiation / mstrnX
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          mstrnX
!!          Ref: Nakajima and Tanaka(1986)
!!               Nakajima et al.(2000)
!!               Sekiguchi and Nakajima(2008)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_rd_mstrnx
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_cpl_sfc_index

  use scale_atmos_phy_rd_common, only: &
     I_SW, &
     I_LW, &
     I_dn, &
     I_up

  use scale_atmos_hydrometeor, only: &
     N_HYD
  use scale_atmos_aerosol, only: &
     N_AE
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_mstrnx_setup
  public :: ATMOS_PHY_RD_mstrnx_finalize
  public :: ATMOS_PHY_RD_mstrnx_flux

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: RD_MSTRN_setup
  private :: RD_MSTRN_DTRN3
  private :: RD_MSTRN_two_stream

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: RD_cosSZA_min = 0.017_RP ! minimum SZA (>89.0)
  real(RP), private, parameter :: RD_EPS        = 1.E-4_RP

  real(RP), private :: RD_TOA  = 100.0_RP !< top of atmosphere [km]
  integer,  private :: RD_KADD = 10       !< RD_KMAX = KMAX + RD_KADD
  integer,  private, parameter :: RD_naero      = N_HYD + N_AE ! # of cloud/aerosol species
  integer,  private, parameter :: RD_hydro_str  = 1            ! start index for cloud
  integer,  private, parameter :: RD_hydro_end  = N_HYD        ! end   index for cloud
  integer,  private, parameter :: RD_aero_str   = N_HYD + 1    ! start index for aerosol
  integer,  private, parameter :: RD_aero_end   = N_HYD + N_AE ! end   index for aerosol

  integer,  private :: RD_KMAX      ! # of computational cells: z for radiation scheme

  real(RP), private, allocatable :: RD_zh          (:)   ! altitude    at the interface [km]
  real(RP), private, allocatable :: RD_z           (:)   ! altitude    at the center    [km]
  real(RP), private, allocatable :: RD_rhodz       (:)   ! density * delta z            [kg/m2]
  real(RP), private, allocatable :: RD_pres        (:)   ! pressure    at the center    [hPa]
  real(RP), private, allocatable :: RD_presh       (:)   ! pressure    at the interface [hPa]
  real(RP), private, allocatable :: RD_temp        (:)   ! temperature at the center    [K]
  real(RP), private, allocatable :: RD_temph       (:)   ! temperature at the interface [K]
  real(RP), private, allocatable :: RD_gas         (:,:) ! gas species   volume mixing ratio [ppmv]
  real(RP), private, allocatable :: RD_cfc         (:,:) ! CFCs          volume mixing ratio [ppmv]
  real(RP), private, allocatable :: RD_aerosol_conc(:,:) ! cloud/aerosol volume mixing ratio [ppmv]
  real(RP), private, allocatable :: RD_aerosol_radi(:,:) ! cloud/aerosol effective radius    [cm]
  real(RP), private, allocatable :: RD_cldfrac     (:)   ! cloud fraction    (0-1)

  integer,  private :: I_MPAE2RD(RD_naero) ! look-up table between input aerosol category and MSTRN particle type

  character(len=H_LONG), private :: MSTRN_GASPARA_INPUTFILE   = 'PARAG.29'     !< input file (gas parameter)
  character(len=H_LONG), private :: MSTRN_AEROPARA_INPUTFILE  = 'PARAPC.29'    !< input file (particle parameter)
  character(len=H_LONG), private :: MSTRN_HYGROPARA_INPUTFILE = 'VARDATA.RM29' !< input file (hygroscopic parameter)

  integer,  private            :: MSTRN_nband    = 29 !< # of wave bands

  integer,  private, parameter :: MSTRN_nstream  =  1 !< # of streams
  integer,  private, parameter :: MSTRN_ch_limit = 10 !< max # of subintervals
  integer,  private, parameter :: MSTRN_nflag    =  7 ! # of optical properties flag
  integer,  private, parameter :: MSTRN_nfitP    = 26 ! # of fitting point for log(pressure)
  integer,  private, parameter :: MSTRN_nfitT    =  3 ! # of fitting point for temperature

  integer,  private, parameter :: MSTRN_ngas     =  7 !< # of gas species
                                                      !   1: H2O
                                                      !   2: CO2
                                                      !   3: O3
                                                      !   4: N2O
                                                      !   5: CO
                                                      !   6: CH4
                                                      !   7: O2
  integer,  private, parameter :: MSTRN_ncfc     = 28 !< # of CFC species
                                                      !   1: CFC-11
                                                      !   2: CFC-12
                                                      !   3: CFC-13
                                                      !   4: CFC-14
                                                      !   5: CFC-113
                                                      !   6: CFC-114
                                                      !   7: CFC-115
                                                      !   8: HCFC-21
                                                      !   9: HCFC-22
                                                      !  10: HCFC-123
                                                      !  11: HCFC-124
                                                      !  12: HCFC-141b
                                                      !  13: HCFC-142b
                                                      !  14: HCFC-225ca
                                                      !  15: HCFC-225cb
                                                      !  16: HFC-32
                                                      !  17: HFC-125
                                                      !  18: HFC-134
                                                      !  19: HFC-134a
                                                      !  20: HFC-143a
                                                      !  21: HFC-152a
                                                      !  22: SF6
                                                      !  23: ClONO2
                                                      !  24: CCl4
                                                      !  25: N2O5
                                                      !  26: C2F6
                                                      !  27: HNO4
  integer,  private            :: MSTRN_nptype   =  9 !< # of particle species
                                                      !  1: water cloud
                                                      !  2: ice cloud
                                                      !  3: Soil dust
                                                      !  4: Carbonacerous (BC/OC=0.3)
                                                      !  5: Carbonacerous (BC/OC=0.15)
                                                      !  6: Carbonacerous (BC/OC=0.)
                                                      !  7: Black carbon
                                                      !  8: Sulfate
                                                      !  9: Sea salt

  integer,  private, parameter :: MSTRN_nsfc     =  7 !< # of surface types
                                                      !  1: ocean
                                                      !  2: wet land
                                                      !  3: dry land
                                                      !  4: low plants
                                                      !  5: forest
                                                      !  6: snow
                                                      !  7: ice
  integer,  private, parameter :: MSTRN_nfitPLK  =  5 !< # of fitting point for planck function
  integer,  private, parameter :: MSTRN_nplkord  =  3 !< # of orders for planck function
  integer,  private, parameter :: MSTRN_nmoment  =  6 !< absorption + # of moments for scattering phase function
  integer,  private            :: MSTRN_nradius  =  8 !< # of radius mode for hygroscopic parameter
  integer,  private, parameter :: MSTRN_ncloud   =  2 !< # of cloud types [ClearSky/Cloud]

  logical,  private            :: ATMOS_PHY_RD_MSTRN_ONLY_QCI        = .false.
  logical,  private            :: ATMOS_PHY_RD_MSTRN_ONLY_TROPOCLOUD = .false.
  logical,  private            :: ATMOS_PHY_RD_MSTRN_USE_AERO        = .false.



  real(RP), private, allocatable :: waveh   (:)         ! wavenumbers at band boundary [1/cm]

  real(RP), private, allocatable :: logfitP (:)         ! fitting point for log10(pressure)
  real(RP), private, allocatable :: fitT    (:)         ! fitting point for temperature
  real(RP), private, allocatable :: logfitT (:)         ! fitting point for log10(temperature)
  integer,  private, allocatable :: iflgb   (:,:)       ! optical properties flag   in each band
  integer,  private, allocatable :: nch     (:)         ! number  of subintervals   in each band
  real(RP), private, allocatable :: wgtch   (:,:)       ! weights of subintervals   in each band
  integer,  private, allocatable :: ngasabs (:)         ! number  of absorbers(gas) in each band
  integer,  private, allocatable :: igasabs (:,:)       ! index   of absorbers(gas) in each band

  real(RP), private, allocatable :: akd     (:,:,:,:,:) ! absorption coefficient table
  real(RP), private, allocatable :: skd     (:,:,:,:)   ! absorption coefficient table for H2O self broadening
  real(RP), private, allocatable :: acfc_pow(:,:)       ! 10 to the power of absorption coefficient table for CFC

  real(RP), private, allocatable :: fitPLK  (:,:)       ! fitting point for planck function
  real(RP), private, allocatable :: fsol    (:)         ! solar insolation    in each band
  real(RP), private              :: fsol_tot            ! total solar insolation
  real(RP), private, allocatable :: sfc     (:,:)       ! surface condition   in each band
  real(RP), private, allocatable :: rayleigh(:)         ! rayleigh scattering in each band
  real(RP), private, allocatable :: qmol    (:,:)       ! moments for rayleigh scattering phase function
  real(RP), private, allocatable :: q       (:,:,:,:)   ! moments for aerosol  scattering phase function

  integer,  private, allocatable :: hygro_flag(:)       ! flag for hygroscopic enlargement
  real(RP), private, allocatable :: radmode   (:,:)     ! radius mode for hygroscopic parameter

  integer,  private, allocatable :: ptype_nradius(:)    ! number limit of radius mode for cloud / aerosol particle type


  ! index for optical flag iflgb
  integer,  private, parameter :: I_SWLW          = 4
  integer,  private, parameter :: I_H2O_continuum = 5
  integer,  private, parameter :: I_CFC_continuum = 7
  ! for cloud type
  integer,  private, parameter :: I_ClearSky = 1
  integer,  private, parameter :: I_Cloud    = 2

  ! pre-calc
  real(RP), private :: RHO_std          ! rho(0C,1atm) [kg/m3]

  real(RP), private :: M(2)             ! discrete quadrature mu for two-stream approximation
  real(RP), private :: W(2)             ! discrete quadrature w  for two-stream approximation
  real(RP), private :: Wmns(2), Wpls(2) ! W-, W+
  real(RP), private :: Wbar(2), Wscale(2)
  !$acc declare create(M, W, Wmns, Wpls, Wscale)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_mstrnx_setup( &
       KA, KS, KE, &
       CZ, FZ )
    use scale_prc, only: &
       PRC_abort
    use scale_time, only: &
       TIME_NOWDATE
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT
    use scale_atmos_phy_rd_profile, only: &
       RD_PROFILE_setup       => ATMOS_PHY_RD_PROFILE_setup,       &
       RD_PROFILE_setup_zgrid => ATMOS_PHY_RD_PROFILE_setup_zgrid, &
       RD_PROFILE_read        => ATMOS_PHY_RD_PROFILE_read
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use scale_atmos_aerosol, only: &
       N_AE
    implicit none
    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in) :: CZ(KA)
    real(RP), intent(in) :: FZ(0:KA)

    real(RP)              :: ATMOS_PHY_RD_MSTRN_TOA
    integer               :: ATMOS_PHY_RD_MSTRN_KADD
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME
    integer               :: ATMOS_PHY_RD_MSTRN_nband
    integer               :: ATMOS_PHY_RD_MSTRN_nptype
    integer               :: ATMOS_PHY_RD_MSTRN_nradius
    integer               :: ATMOS_PHY_RD_MSTRN_nradius_cloud
    integer               :: ATMOS_PHY_RD_MSTRN_nradius_aero

    namelist / PARAM_ATMOS_PHY_RD_MSTRN / &
       ATMOS_PHY_RD_MSTRN_TOA,                   &
       ATMOS_PHY_RD_MSTRN_KADD,                  &
       ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME,   &
       ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME,  &
       ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME, &
       ATMOS_PHY_RD_MSTRN_nband,                 &
       ATMOS_PHY_RD_MSTRN_nptype,                &
       ATMOS_PHY_RD_MSTRN_nradius,               &
       ATMOS_PHY_RD_MSTRN_nradius_cloud,         &
       ATMOS_PHY_RD_MSTRN_nradius_aero,          &
       ATMOS_PHY_RD_MSTRN_ONLY_QCI,              &
       ATMOS_PHY_RD_MSTRN_ONLY_TROPOCLOUD,       &
       ATMOS_PHY_RD_MSTRN_USE_AERO

    integer :: KMAX
    integer :: ngas, ncfc
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_mstrnx_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_RD_mstrnx_setup",*) 'Sekiguchi and Nakajima (2008) mstrnX radiation process'

    !--- read namelist
    ATMOS_PHY_RD_MSTRN_TOA                   = RD_TOA
    ATMOS_PHY_RD_MSTRN_KADD                  = RD_KADD
    ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME   = MSTRN_GASPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME  = MSTRN_AEROPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = MSTRN_HYGROPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_nband                 = MSTRN_nband
    ATMOS_PHY_RD_MSTRN_nptype                = MSTRN_nptype
    ATMOS_PHY_RD_MSTRN_nradius               = MSTRN_nradius
    ATMOS_PHY_RD_MSTRN_nradius_cloud         = -1
    ATMOS_PHY_RD_MSTRN_nradius_aero          = -1

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_MSTRN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_RD_mstrnx_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_RD_mstrnx_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD_MSTRN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_RD_MSTRN)

    RD_TOA                    = ATMOS_PHY_RD_MSTRN_TOA
    RD_KADD                   = ATMOS_PHY_RD_MSTRN_KADD
    MSTRN_GASPARA_INPUTFILE   = ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME
    MSTRN_AEROPARA_INPUTFILE  = ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME
    MSTRN_HYGROPARA_INPUTFILE = ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME
    MSTRN_nband               = ATMOS_PHY_RD_MSTRN_nband
    MSTRN_nptype              = ATMOS_PHY_RD_MSTRN_nptype

    if ( ATMOS_PHY_RD_MSTRN_nradius_cloud < 0 ) then ! set default
       ATMOS_PHY_RD_MSTRN_nradius_cloud = ATMOS_PHY_RD_MSTRN_nradius
    endif
    if ( ATMOS_PHY_RD_MSTRN_nradius_aero  < 0 ) then ! set default
       ATMOS_PHY_RD_MSTRN_nradius_aero  = ATMOS_PHY_RD_MSTRN_nradius
    endif

    MSTRN_nradius = max( ATMOS_PHY_RD_MSTRN_nradius,       &
                         ATMOS_PHY_RD_MSTRN_nradius_cloud, &
                         ATMOS_PHY_RD_MSTRN_nradius_aero   )

    !--- setup MP/AE - RD(particle) coupling parameter
    allocate( ptype_nradius(MSTRN_nptype) )

    if ( MSTRN_nptype == 9 ) then

       I_MPAE2RD( 1) =  1 ! cloud water Qc -> Water
       I_MPAE2RD( 2) =  1 ! rain        Qr -> Water
       I_MPAE2RD( 3) =  2 ! cloud ice   Qi -> Ice
       I_MPAE2RD( 4) =  2 ! snow        Qs -> Ice
       I_MPAE2RD( 5) =  2 ! graupel     Qg -> Ice
       I_MPAE2RD( 6) =  2 ! hail        Qh -> Ice
       I_MPAE2RD( 7) =  3 ! aerosol type1 -> Soil dust
       I_MPAE2RD( 8) =  4 ! aerosol type2 -> Carbonacerous (BC/OC=0.3)
       I_MPAE2RD( 9) =  5 ! aerosol type3 -> Carbonacerous (BC/OC=0.15)
       I_MPAE2RD(10) =  6 ! aerosol type4 -> Carbonacerous (BC/OC=0.)
       I_MPAE2RD(11) =  7 ! aerosol type5 -> Black carbon
       I_MPAE2RD(12) =  8 ! aerosol type6 -> Sulfate
       I_MPAE2RD(13) =  9 ! aerosol type7 -> Sea salt

       ptype_nradius( 1) =  ATMOS_PHY_RD_MSTRN_nradius_cloud
       ptype_nradius( 2) =  ATMOS_PHY_RD_MSTRN_nradius_cloud
       ptype_nradius( 3) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 4) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 5) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 6) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 7) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 8) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 9) =  ATMOS_PHY_RD_MSTRN_nradius_aero

    elseif( MSTRN_nptype == 12 ) then

       I_MPAE2RD( 1) =  1 ! cloud water Qc -> CLOUD
       I_MPAE2RD( 2) =  2 ! rain        Qr -> RAIN
       I_MPAE2RD( 3) =  3 ! cloud ice   Qi -> ICE
       I_MPAE2RD( 4) =  4 ! snow        Qs -> SNOW
       I_MPAE2RD( 5) =  5 ! graupel     Qg -> GRAUPEL
       I_MPAE2RD( 6) =  5 ! hail        Qh -> GRAUPEL
       I_MPAE2RD( 7) =  6 ! aerosol type1 -> Soil dust
       I_MPAE2RD( 8) =  7 ! aerosol type2 -> Carbonacerous (BC/OC=0.3)
       I_MPAE2RD( 9) =  8 ! aerosol type3 -> Carbonacerous (BC/OC=0.15)
       I_MPAE2RD(10) =  9 ! aerosol type4 -> Carbonacerous (BC/OC=0.)
       I_MPAE2RD(11) = 10 ! aerosol type5 -> Black carbon
       I_MPAE2RD(12) = 11 ! aerosol type6 -> Sulfate
       I_MPAE2RD(13) = 12 ! aerosol type7 -> Sea salt

       ptype_nradius( 1) =  ATMOS_PHY_RD_MSTRN_nradius_cloud
       ptype_nradius( 2) =  ATMOS_PHY_RD_MSTRN_nradius_cloud
       ptype_nradius( 3) =  ATMOS_PHY_RD_MSTRN_nradius_cloud
       ptype_nradius( 4) =  ATMOS_PHY_RD_MSTRN_nradius_cloud
       ptype_nradius( 5) =  ATMOS_PHY_RD_MSTRN_nradius_cloud
       ptype_nradius( 6) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 7) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 8) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius( 9) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius(10) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius(11) =  ATMOS_PHY_RD_MSTRN_nradius_aero
       ptype_nradius(12) =  ATMOS_PHY_RD_MSTRN_nradius_aero

    endif
    !$acc enter data copyin(I_MPAE2RD, ptype_nradius)

    !--- setup MSTRN parameter
    call RD_MSTRN_setup( ngas, & ! [OUT]
                         ncfc  ) ! [OUT]

    !--- setup climatological profile
    call RD_PROFILE_setup

    KMAX = KE - KS + 1
    RD_KMAX = KMAX + RD_KADD

    !--- allocate arrays
    ! input
    allocate( RD_zh   (RD_KMAX+1) )
    allocate( RD_z    (RD_KMAX  ) )

    allocate( RD_rhodz(RD_KMAX  ) )
    allocate( RD_pres (RD_KMAX  ) )
    allocate( RD_presh(RD_KMAX+1) )
    allocate( RD_temp (RD_KMAX  ) )
    allocate( RD_temph(RD_KMAX+1) )

    allocate( RD_gas         (RD_KMAX,ngas    ) )
    allocate( RD_cfc         (RD_KMAX,ncfc    ) )
    allocate( RD_aerosol_conc(RD_KMAX,RD_naero) )
    allocate( RD_aerosol_radi(RD_KMAX,RD_naero) )
    allocate( RD_cldfrac     (RD_KMAX         ) )

    !--- setup vartical grid for radiation (larger TOA than Model domain)
    call RD_PROFILE_setup_zgrid( &
         KA, KS, KE, &
         RD_KMAX, RD_KADD,     & ! [IN]
         RD_TOA, CZ(:), FZ(:), & ! [IN]
         RD_zh(:), RD_z(:)     ) ! [OUT]

    !--- read climatological profile
    call RD_PROFILE_read( RD_KMAX,                & ! [IN]
                          ngas,                   & ! [IN]
                          ncfc,                   & ! [IN]
                          RD_naero,               & ! [IN]
                          ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT,     & ! [IN]
                          TIME_NOWDATE   (:),     & ! [IN]
                          RD_zh          (:),     & ! [IN]
                          RD_z           (:),     & ! [IN]
                          RD_rhodz       (:),     & ! [OUT]
                          RD_pres        (:),     & ! [OUT]
                          RD_presh       (:),     & ! [OUT]
                          RD_temp        (:),     & ! [OUT]
                          RD_temph       (:),     & ! [OUT]
                          RD_gas         (:,:),   & ! [OUT]
                          RD_cfc         (:,:),   & ! [OUT]
                          RD_aerosol_conc(:,:),   & ! [OUT]
                          RD_aerosol_radi(:,:),   & ! [OUT]
                          RD_cldfrac     (:)      ) ! [OUT]

    !$acc enter data copyin(RD_rhodz,RD_pres,RD_presh,RD_temp,RD_temph,RD_gas,RD_cfc,RD_aerosol_conc,RD_aerosol_radi,RD_cldfrac)

    return
  end subroutine ATMOS_PHY_RD_mstrnx_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_RD_mstrnx_finalize

    !$acc exit data delete(I_MPAE2RD, ptype_nradius)

    !$acc exit data &
    !$acc delete(waveh, &
    !$acc        wgtch, fitPLK, logfitP, logfitT, fitT, &
    !$acc        radmode, ngasabs, igasabs, &
    !$acc        fsol, q, qmol, rayleigh, acfc_pow, nch, AKD, SKD )

    !$acc exit data delete(RD_rhodz,RD_pres,RD_presh,RD_temp,RD_temph,RD_gas,RD_cfc,RD_aerosol_conc,RD_aerosol_radi,RD_cldfrac)

    deallocate( ptype_nradius )
    deallocate( RD_zh    )
    deallocate( RD_z     )

    deallocate( RD_rhodz )
    deallocate( RD_pres  )
    deallocate( RD_presh )
    deallocate( RD_temp  )
    deallocate( RD_temph )

    deallocate( RD_gas          )
    deallocate( RD_cfc          )
    deallocate( RD_aerosol_conc )
    deallocate( RD_aerosol_radi )
    deallocate( RD_cldfrac      )

    deallocate( waveh   )
    deallocate( logfitP )
    deallocate( fitT    )
    deallocate( logfitT )

    deallocate( iflgb   )
    deallocate( nch     )
    deallocate( wgtch   )
    deallocate( ngasabs )
    deallocate( igasabs )

    deallocate( akd      )
    deallocate( skd      )
    deallocate( acfc_pow )

    deallocate( fitPLK   )
    deallocate( fsol     )
    deallocate( sfc      )
    deallocate( rayleigh )

    deallocate( qmol )
    deallocate( q    )

    deallocate( hygro_flag )
    deallocate( radmode    )

    return
  end subroutine ATMOS_PHY_RD_mstrnx_finalize

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_mstrnx_flux( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, TEMP, PRES, QV,  &
       CZ, FZ,                &
       fact_ocean,            &
       fact_land,             &
       fact_urban,            &
       temp_sfc, albedo_sfc,  &
       solins, cosSZA,        &
       CLDFRAC, MP_Re, MP_Qe, &
       AE_Re, AE_Qe,          &
       flux_rad,              &
       flux_rad_top,          &
       flux_rad_sfc_dn,       &
       dtau_s, dem_s          )
       !Jval                   )
#ifdef _OPENACC
    use openacc
#endif
    use scale_const, only: &
       EPS  => CONST_EPS, &
       Mdry => CONST_Mdry, &
       Mvap => CONST_Mvap, &
       PPM  => CONST_PPM
    use scale_time, only: &
       TIME_NOWDATE
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HI, &
       HYD_DENS
    use scale_atmos_aerosol, only: &
       N_AE, &
       AE_DENS
    use scale_atmos_phy_rd_profile, only: &
       RD_PROFILE_read            => ATMOS_PHY_RD_PROFILE_read, &
       RD_PROFILE_use_climatology => ATMOS_PHY_RD_PROFILE_use_climatology
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: DENS           (KA,IA,JA)
    real(RP), intent(in)  :: TEMP           (KA,IA,JA)
    real(RP), intent(in)  :: PRES           (KA,IA,JA)
    real(RP), intent(in)  :: QV             (KA,IA,JA)
    real(RP), intent(in)  :: CZ             (  KA,IA,JA)
    real(RP), intent(in)  :: FZ             (0:KA,IA,JA)
    real(RP), intent(in)  :: fact_ocean     (IA,JA)
    real(RP), intent(in)  :: fact_land      (IA,JA)
    real(RP), intent(in)  :: fact_urban     (IA,JA)
    real(RP), intent(in)  :: temp_sfc       (IA,JA)
    real(RP), intent(in)  :: albedo_sfc     (IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(in)  :: solins         (IA,JA)
    real(RP), intent(in)  :: cosSZA         (IA,JA)
    real(RP), intent(in)  :: cldfrac        (KA,IA,JA)
    real(RP), intent(in)  :: MP_Re          (KA,IA,JA,N_HYD)
    real(RP), intent(in)  :: MP_Qe          (KA,IA,JA,N_HYD)
    real(RP), intent(in)  :: AE_Re          (KA,IA,JA,N_AE)
    real(RP), intent(in)  :: AE_Qe          (KA,IA,JA,N_AE)
    real(RP), intent(out) :: flux_rad       (KA,IA,JA,2,2,2)
    real(RP), intent(out) :: flux_rad_top   (IA,JA,2,2,2)
    real(RP), intent(out) :: flux_rad_sfc_dn(IA,JA,N_RAD_DIR,N_RAD_RGN)

    real(RP), intent(out), optional :: dtau_s(KA,IA,JA) ! 0.67 micron cloud optical depth
    real(RP), intent(out), optional :: dem_s (KA,IA,JA) ! 10.5 micron cloud emissivity
!    real(RP), intent(out), optional :: Jval  (KA,IA,JA,CH_QA_photo)

    integer  :: tropopause(IA,JA)
    real(RP) :: gamma

    real(RP), parameter :: min_cldfrac = 1.E-8_RP

    real(RP) :: rhodz_merge       (RD_KMAX,IA,JA)
    real(RP) :: pres_merge        (RD_KMAX,IA,JA)
    real(RP) :: temp_merge        (RD_KMAX,IA,JA)
    real(RP) :: temph_merge       (RD_KMAX+1,IA,JA)

    real(RP) :: gas_merge         (RD_KMAX,IA,JA,MSTRN_ngas)
    real(RP) :: cfc_merge         (RD_KMAX,IA,JA,MSTRN_ncfc)
    real(RP) :: aerosol_conc_merge(RD_KMAX,IA,JA,RD_naero  )
    real(RP) :: aerosol_radi_merge(RD_KMAX,IA,JA,RD_naero  )
    real(RP) :: cldfrac_merge     (RD_KMAX,IA,JA)

    ! output
    real(RP) :: flux_rad_merge(RD_KMAX+1,IA,JA,2,2,MSTRN_ncloud)
    real(RP) :: tauCLD_067u   (RD_KMAX,IA,JA) ! 0.67 micron cloud optical depth
    real(RP) :: emisCLD_105u  (RD_KMAX,IA,JA) ! 10.5 micron cloud emissivity

    real(RP) :: zerosw

    integer :: ihydro, iaero
    integer :: RD_k, k, i, j, v, ic
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / radiation / mstrnX'

    call PROF_rapstart('RD_Profile', 3)

    !$acc data copyin(DENS,TEMP,PRES,QV,CZ,FZ,fact_ocean,fact_land,fact_urban,temp_sfc,albedo_sfc,solins,cosSZA,cldfrac,MP_Re,AE_Re,AE_Qe) &
    !$acc      copyout(flux_rad,flux_rad_top,flux_rad_sfc_dn) &
    !$acc      create(tropopause,rhodz_merge,pres_merge,temp_merge,temph_merge,gas_merge,cfc_merge,aerosol_conc_merge,aerosol_radi_merge,cldfrac_merge,flux_rad_merge,tauCLD_067u,emisCLD_105u)

    !$acc data copyout(dtau_s) if(acc_is_present(dtau_s))
    !$acc data copyout(dem_s) if(acc_is_present(dem_s))

    if ( ATMOS_PHY_RD_MSTRN_ONLY_TROPOCLOUD ) then
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j, &
       !$omp         gamma) &
       !$omp shared (tropopause,pres,temp,CZ, &
       !$omp         KS,KE,IS,IE,JS,JE)
       !$acc kernels
       do j  = JS, JE
       do i  = IS, IE
          tropopause(i,j) = KE+1
          do k  = KE, KS, -1
             if ( pres(k,i,j) >= 300.E+2_RP ) then
                exit
             elseif( pres(k,i,j) <  50.E+2_RP ) then
                tropopause(i,j) = k
             else
                gamma = ( temp(k+1,i,j) - temp(k,i,j) ) / ( CZ(k+1,i,j) - CZ(k,i,j) )
                if( gamma > 1.E-4_RP ) tropopause(i,j) = k
             endif
          enddo
       enddo
       enddo
       !$acc end kernels
    else
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,j) &
       !$omp shared (tropopause, &
       !$omp         KE,IS,IE,JS,JE)
!OCL XFILL
       !$acc kernels
       do j  = JS, JE
       do i  = IS, IE
          tropopause(i,j) = KE+1
       end do
       end do
       !$acc end kernels
    endif

    ! marge basic profile and value in model domain

    if ( RD_PROFILE_use_climatology ) then
       call RD_PROFILE_read( RD_KMAX,                               & ! [IN]
                             MSTRN_ngas,                            & ! [IN]
                             MSTRN_ncfc,                            & ! [IN]
                             RD_naero,                              & ! [IN]
                             ATMOS_GRID_CARTESC_REAL_BASEPOINT_LAT, & ! [IN]
                             TIME_NOWDATE   (:),                    & ! [IN]
                             RD_zh          (:),                    & ! [IN]
                             RD_z           (:),                    & ! [IN]
                             RD_rhodz       (:),                    & ! [OUT]
                             RD_pres        (:),                    & ! [OUT]
                             RD_presh       (:),                    & ! [OUT]
                             RD_temp        (:),                    & ! [OUT]
                             RD_temph       (:),                    & ! [OUT]
                             RD_gas         (:,:),                  & ! [OUT]
                             RD_cfc         (:,:),                  & ! [OUT]
                             RD_aerosol_conc(:,:),                  & ! [OUT]
                             RD_aerosol_radi(:,:),                  & ! [OUT]
                             RD_cldfrac     (:)                     ) ! [OUT]
    endif

!OCL XFILL
    !$omp parallel do default(none)                                           &
    !$omp shared(JS,JE,IS,IE,RD_KADD,temph_merge,RD_temph,KE,RD_KMAX,KS,temp,CZ,FZ) &
    !$omp private(i,j,k,RD_k) OMP_SCHEDULE_ collapse(2)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KADD
          temph_merge(RD_k,i,j) = RD_temph(RD_k)
       enddo

       temph_merge(RD_KADD+1,i,j) = temp(KE,i,j)
       do RD_k = RD_KADD+2, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis

          temph_merge(RD_k,i,j) = ( temp(k  ,i,j) * (CZ(k+1,i,j)-FZ(k,i,j)) / (CZ(k+1,i,j)-CZ(k,i,j)) &
                                  + temp(k+1,i,j) * (FZ(k  ,i,j)-CZ(k,i,j)) / (CZ(k+1,i,j)-CZ(k,i,j)) )
       enddo
       temph_merge(RD_KMAX+1,i,j) = temp(KS,i,j)
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(k,i,j,RD_k) &
    !$omp shared (RD_temp,dens,FZ,pres,temp, &
    !$omp         rhodz_merge,RD_rhodz,pres_merge,RD_pres,temp_merge, &
    !$omp         KS,JS,JE,IS,IE,RD_KMAX,RD_KADD)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KADD
          rhodz_merge(RD_k,i,j) = RD_rhodz(RD_k)
          pres_merge (RD_k,i,j) = RD_pres (RD_k)
          temp_merge (RD_k,i,j) = RD_temp (RD_k)
       enddo

       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis

          rhodz_merge(RD_k,i,j) = dens(k,i,j) * ( FZ(k,i,j)-FZ(k-1,i,j) ) ! [kg/m2]
          pres_merge (RD_k,i,j) = pres(k,i,j) * 1.E-2_RP ! [hPa]
          temp_merge (RD_k,i,j) = temp(k,i,j)
       enddo
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$omp parallel default(none) &
    !$omp private(v,i,j,RD_k) &
    !$omp shared (gas_merge,RD_gas, &
    !$omp         IS,IE,JS,JE,RD_KMAX)
    !$acc kernels
    do v = 1,  MSTRN_ngas
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do RD_k = 1, RD_KMAX
          gas_merge(RD_k,i,j,v) = RD_gas(RD_k,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
    enddo
    !$acc end kernels
    !$omp end parallel

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(k,i,j,RD_k, &
    !$omp         zerosw ) &
    !$omp shared (gas_merge,QV,EPS,Mvap,Mdry, &
    !$omp         KS,IS,IE,JS,JE,RD_KADD,RD_KMAX)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do RD_k = RD_KADD+1, RD_KMAX
       k = KS + RD_KMAX - RD_k ! reverse axis
       zerosw = sign(0.5_RP, QV(k,i,j)-EPS) + 0.5_RP
       gas_merge(RD_k,i,j,1) = QV(k,i,j) / Mvap * Mdry / PPM * zerosw ! [PPM]
    enddo
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$omp parallel default(none) &
    !$omp private(v,i,j,RD_k) &
    !$omp shared (cfc_merge,RD_cfc, &
    !$omp         IS,IE,JS,JE,RD_KMAX)
    !$acc kernels
    do v = 1,  MSTRN_ncfc
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do RD_k = 1, RD_KMAX
          cfc_merge(RD_k,i,j,v) = RD_cfc(RD_k,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
    enddo
    !$acc end kernels
    !$omp end parallel

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
    !$omp private(k,i,j,RD_k) &
    !$omp shared (cldfrac_merge,RD_cldfrac,cldfrac, &
    !$omp         KS,IS,IE,JS,JE,RD_KMAX,RD_KADD)
    !$acc kernels
    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KADD
          cldfrac_merge(RD_k,i,j) = RD_cldfrac(RD_k)
       enddo

       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis

          cldfrac_merge(RD_k,i,j) = min( max( cldfrac(k,i,j), 0.0_RP ), 1.0_RP )
!          cldfrac_merge(RD_k,i,j) = 0.5_RP + sign( 0.5_RP, cldfrac(k,i,j)-min_cldfrac )
       enddo
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$omp parallel default(none) &
    !$omp private(v,i,j,RD_k) &
    !$omp shared (aerosol_conc_merge,RD_aerosol_conc,aerosol_radi_merge,RD_aerosol_radi, &
    !$omp         IS,IE,JS,JE,RD_KADD)
    !$acc kernels
    do v = 1,  RD_naero
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
       do RD_k = 1, RD_KADD
          aerosol_conc_merge(RD_k,i,j,v) = RD_aerosol_conc(RD_k,v)
          aerosol_radi_merge(RD_k,i,j,v) = RD_aerosol_radi(RD_k,v)
       enddo
       enddo
       enddo
       !$omp end do nowait
    enddo
    !$acc end kernels
    !$omp end parallel

    !$omp parallel default(none) &
    !$omp private(ihydro,k,i,j,RD_k) &
    !$omp shared (aerosol_conc_merge,tropopause,MP_Qe,HYD_DENS,RHO_std, &
    !$omp         ATMOS_PHY_RD_MSTRN_ONLY_QCI, &
    !$omp         KS,KE,IS,IE,JS,JE,RD_KMAX,RD_KADD)
    do ihydro = 1, N_HYD
       if ( ATMOS_PHY_RD_MSTRN_ONLY_QCI .and. &
            ( ihydro /= I_HC .and. ihydro /= I_HI ) ) then
!OCL XFILL
          !$omp do OMP_SCHEDULE_ collapse(2)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do RD_k = RD_KADD+1, RD_KMAX
             aerosol_conc_merge(RD_k,i,j,ihydro) = 0.0_RP
          end do
          end do
          end do
          !$acc end kernels
          !$omp end do nowait
       else
          !$omp do OMP_SCHEDULE_ collapse(2)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
             do RD_k = RD_KADD+1, RD_KADD+1 + KE - tropopause(i,j)
                aerosol_conc_merge(RD_k,i,j,ihydro) = 0.0_RP
             end do
             !$acc loop independent
             do RD_k = RD_KADD+1 + KE - tropopause(i,j) + 1, RD_KMAX
                k = KS + RD_KMAX - RD_k ! reverse axis
                aerosol_conc_merge(RD_k,i,j,ihydro) = max( MP_Qe(k,i,j,ihydro), 0.0_RP ) &
                                                    / HYD_DENS(ihydro) * RHO_std / PPM ! [PPM to standard air]
             enddo
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       end if
    enddo
    !$omp end parallel

    !$omp parallel do default(none) OMP_SCHEDULE_ collapse(3) &
    !$omp private(ihydro,k,i,j,RD_k) &
    !$omp shared (aerosol_radi_merge,MP_Re, &
    !$omp         KS,IS,IE,JS,JE,RD_KMAX,RD_KADD)
    !$acc kernels
    do ihydro = 1, N_HYD
    do j = JS, JE
    do i = IS, IE
    do RD_k = RD_KADD+1, RD_KMAX
       k = KS + RD_KMAX - RD_k ! reverse axis
       aerosol_radi_merge(RD_k,i,j,ihydro) = MP_Re(k,i,j,ihydro)
    enddo
    enddo
    enddo
    enddo
    !$acc end kernels

    !$omp parallel default(none) &
    !$omp private(iaero,k,i,j,RD_k) &
    !$omp shared (aerosol_conc_merge,RD_aerosol_conc,aerosol_radi_merge,RD_aerosol_radi, &
    !$omp         ATMOS_PHY_RD_MSTRN_USE_AERO, &
    !$omp         AE_Qe,RHO_std,AE_Re, &
    !$omp         KS,IS,IE,JS,JE,RD_KMAX,RD_KADD)
    do iaero = 1, N_AE

       if ( ATMOS_PHY_RD_MSTRN_USE_AERO ) then
          !$omp do OMP_SCHEDULE_ collapse(2)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do RD_k = RD_KADD+1, RD_KMAX
             k = KS + RD_KMAX - RD_k ! reverse axis
             aerosol_conc_merge(RD_k,i,j,N_HYD+iaero) = max( AE_Qe(k,i,j,iaero), 0.0_RP ) &
                                                      / AE_DENS(iaero) * RHO_std / PPM ! [PPM to standard air]
             aerosol_radi_merge(RD_k,i,j,N_HYD+iaero) = AE_Re(k,i,j,iaero)
          enddo
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       else
          !$omp do OMP_SCHEDULE_ collapse(2)
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do RD_k = RD_KADD+1, RD_KMAX
             aerosol_conc_merge(RD_k,i,j,N_HYD+iaero) = RD_aerosol_conc(RD_k,N_HYD+iaero)
             aerosol_radi_merge(RD_k,i,j,N_HYD+iaero) = RD_aerosol_radi(RD_k,N_HYD+iaero)
          enddo
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       endif

    enddo
    !$omp end parallel

    call PROF_rapend  ('RD_Profile', 3)
    call PROF_rapstart('RD_MSTRN_DTRN3', 3)

    ! calc radiative transfer
#ifndef _OPENACC
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
#endif
    call RD_MSTRN_DTRN3( RD_KMAX, &
#ifdef _OPENACC
                         IA, IS, IE, JA, JS, JE, &
#else
                         IA, i, i, JA, j, j, &
                         i, j, &
#endif
                         MSTRN_ngas, MSTRN_ncfc, RD_naero,                         & ! [IN]
                         RD_hydro_str, RD_hydro_end, RD_aero_str, RD_aero_end,     & ! [IN]
                         solins(:,:), cosSZA(:,:),                                 & ! [IN]
                         rhodz_merge(:,:,:), pres_merge(:,:,:),                    & ! [IN]
                         temp_merge(:,:,:), temph_merge(:,:,:), temp_sfc(:,:),     & ! [IN]
                         gas_merge(:,:,:,:), cfc_merge(:,:,:,:),                   & ! [IN]
                         aerosol_conc_merge(:,:,:,:), aerosol_radi_merge(:,:,:,:), & ! [IN]
                         I_MPAE2RD(:),                                             & ! [IN]
                         cldfrac_merge(:,:,:),                                     & ! [IN]
                         albedo_sfc(:,:,:,:),                                      & ! [IN]
                         fact_ocean(:,:), fact_land(:,:), fact_urban(:,:),         & ! [IN]
                         flux_rad_merge(:,:,:,:,:,:), flux_rad_sfc_dn(:,:,:,:),    & ! [OUT]
                         tauCLD_067u(:,:,:), emisCLD_105u(:,:,:)                   ) ! [OUT]
#ifndef _OPENACC
    end do
    end do
#endif


    call PROF_rapend  ('RD_MSTRN_DTRN3', 3)

    ! return to grid coordinate of model domain
    !$omp parallel default(none) &
    !$omp private(ic,k,i,j,RD_k) &
    !$omp shared (flux_rad,flux_rad_merge, &
    !$omp         KS,IS,IE,JS,JE,RD_KMAX,RD_KADD)
    !$acc kernels
    !$acc loop independent collapse(4)
    do ic = 1, 2
       !$omp do OMP_SCHEDULE_ collapse(2)
       do j  = JS, JE
       do i  = IS, IE
       do RD_k = RD_KADD+1, RD_KMAX+1
          k = KS + RD_KMAX - RD_k ! reverse axis

          flux_rad(k,i,j,I_LW,I_up,ic) = flux_rad_merge(RD_k,i,j,I_LW,I_up,ic)
          flux_rad(k,i,j,I_LW,I_dn,ic) = flux_rad_merge(RD_k,i,j,I_LW,I_dn,ic)
          flux_rad(k,i,j,I_SW,I_up,ic) = flux_rad_merge(RD_k,i,j,I_SW,I_up,ic)
          flux_rad(k,i,j,I_SW,I_dn,ic) = flux_rad_merge(RD_k,i,j,I_SW,I_dn,ic)
       enddo
       enddo
       enddo
       !$omp end do nowait
    enddo
    !$acc end kernels
    !$omp end parallel

!OCL XFILL
    !$omp parallel default(none) &
    !$omp private(ic,i,j) &
    !$omp shared (flux_rad_top,flux_rad_merge, &
    !$omp         IS,IE,JS,JE)
    !$acc kernels
    do ic = 1, 2
       !$omp do OMP_SCHEDULE_
       do j  = JS, JE
       do i  = IS, IE
          flux_rad_top(i,j,I_LW,I_up,ic) = flux_rad_merge(1,i,j,I_LW,I_up,ic)
          flux_rad_top(i,j,I_LW,I_dn,ic) = flux_rad_merge(1,i,j,I_LW,I_dn,ic)
          flux_rad_top(i,j,I_SW,I_up,ic) = flux_rad_merge(1,i,j,I_SW,I_up,ic)
          flux_rad_top(i,j,I_SW,I_dn,ic) = flux_rad_merge(1,i,j,I_SW,I_dn,ic)
       enddo
       enddo
       !$omp end do nowait
    enddo
    !$acc end kernels
    !$omp end parallel

    if ( present( dtau_s ) ) then
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j,RD_k) &
       !$omp shared (dtau_s,tauCLD_067u, &
       !$omp         KS,IS,IE,JS,JE,RD_KMAX,RD_KADD)
       !$acc kernels
       !$acc loop independent collapse(3)
       do j = JS, JE
       do i = IS, IE
       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis
          dtau_s(k,i,j) = tauCLD_067u(RD_k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
    end if

    if ( present( dem_s ) ) then
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(k,i,j,RD_k) &
       !$omp shared (dem_s,emisCLD_105u, &
       !$omp         KS,IS,IE,JS,JE,RD_KMAX,RD_KADD)
       !$acc kernels
       !$acc loop independent collapse(3)
       do j = JS, JE
       do i = IS, IE
       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis
          dem_s(k,i,j) = emisCLD_105u(RD_k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
    end if

    !$acc end data
    !$acc end data
    !$acc end data

    return
  end subroutine ATMOS_PHY_RD_mstrnx_flux

  !-----------------------------------------------------------------------------
  !> Setup MSTRN parameter table
  subroutine RD_MSTRN_setup( &
       ngas, &
       ncfc  )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       GRAV  => CONST_GRAV, &
       Rdry  => CONST_Rdry, &
       Pstd  => CONST_Pstd, &
       TEM00 => CONST_TEM00
    implicit none

    integer, intent(out) :: ngas
    integer, intent(out) :: ncfc

    real(RP) :: acfc(MSTRN_ncfc)                   !< absorption coefficient table for CFC

    integer :: nband, nstream, nfitP, nfitT, nflag !< gas             parameters for check
    integer :: nsfc, nptype, nplkord, nfitPLK      !< aerosol/surface parameters for check
    integer :: nradius

    character(len=H_LONG) :: dummy, fname

    integer :: fid, ierr
    integer :: iw, ich, ip, it, igas, icfc, iptype, im
    !---------------------------------------------------------------------------

    !---< gas absorption parameter input >---
    ngas = MSTRN_ngas

    ! allocate arrays
    allocate( waveh  (MSTRN_nband+1) )
    allocate( logfitP(MSTRN_nfitP)   )
    allocate( fitT   (MSTRN_nfitT)   )
    allocate( logfitT(MSTRN_nfitT)   )

    allocate( iflgb  (MSTRN_nflag,   MSTRN_nband) )
    allocate( nch    (               MSTRN_nband) )
    allocate( wgtch  (MSTRN_ch_limit,MSTRN_nband) )
    allocate( ngasabs(               MSTRN_nband) )
    allocate( igasabs(MSTRN_ngas,    MSTRN_nband) )

    allocate( akd (MSTRN_nfitP,MSTRN_ch_limit,MSTRN_nfitT,MSTRN_ngas,MSTRN_nband) )
    allocate( skd (MSTRN_nfitP,MSTRN_ch_limit,MSTRN_nfitT,           MSTRN_nband) )
    allocate( acfc_pow(MSTRN_ncfc,MSTRN_nband) )

    fid = IO_get_available_fid()
    call IO_get_fname(fname, MSTRN_GASPARA_INPUTFILE)
    open( fid,                  &
          file   = fname,       &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Input data file does not found! ', trim(fname)
          call PRC_abort
       endif

       ! read gas parameters for check
       read(fid,*) nband, nstream, nfitP, nfitT, nflag, ncfc

       if ( nband /= MSTRN_nband ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nband(given,file)=', MSTRN_nband, nband
          call PRC_abort
       endif
       if ( nstream /= MSTRN_nstream ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nstream(given,file)=', MSTRN_nstream, nstream
          call PRC_abort
       endif
       if ( nfitP /= MSTRN_nfitP ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nfitP(given,file)=', MSTRN_nfitP, nfitP
          call PRC_abort
       endif
       if ( nfitT /= MSTRN_nfitT ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nfitT(given,file)=', MSTRN_nfitT, nfitT
          call PRC_abort
       endif
       if ( nflag /= MSTRN_nflag ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nflag(given,file)=', MSTRN_nflag, nflag
          call PRC_abort
       endif
       if ( ncfc /= MSTRN_ncfc ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! ncfc(given,file)=', MSTRN_ncfc, ncfc
          call PRC_abort
       endif

       ! wave band boundaries
       read(fid,*) dummy
       read(fid,*) waveh(:)
       ! fitting point for log(pressure)
       read(fid,*) dummy
       read(fid,*) logfitP(:)
       ! fitting point for temperature
       read(fid,*) dummy
       read(fid,*) fitT(:)

       logfitT(:) = log10( fitT(:) )

       ! for each band
       do iw = 1, MSTRN_nband

          ! optical properties flag
          read(fid,*) dummy
          read(fid,*) iflgb(:,iw)
          ! number of subintervals
          read(fid,*) dummy
          read(fid,*) nch(iw)
          ! weights for subintervals
          read(fid,*) dummy
          read(fid,*) (wgtch(ich,iw),ich=1,nch(iw))
          ! number of considering gases
          read(fid,*) dummy
          read(fid,*) ngasabs(iw)

          ! major gas absorption
          if ( ngasabs(iw) > 0 ) then
             do igas = 1, ngasabs(iw)
                read(fid,*) igasabs(igas,iw)
                do it  = 1, MSTRN_nfitT
                do ip  = 1, MSTRN_nfitP
                   read(fid,*) (akd(ip,ich,it,igasabs(igas,iw),iw),ich=1,nch(iw))
                enddo
                enddo
             enddo
          endif

          ! H2O continuum
          if ( iflgb(I_H2O_continuum,iw) > 0 ) then
             read(fid,*) dummy
             do it = 1, MSTRN_nfitT
             do ip = 1, MSTRN_nfitP
                read(fid,*) (skd(ip,ich,it,iw),ich=1,nch(iw))
             enddo
             enddo
          endif

          ! CFC absorption
          if ( iflgb(I_CFC_continuum,iw) > 0 ) then
             read(fid,*) dummy
             read(fid,*) (acfc(icfc),icfc=1,MSTRN_ncfc)

             do icfc = 1, MSTRN_ncfc
                acfc_pow(icfc,iw) = 10.0_RP**acfc(icfc)
             end do
          endif

       enddo ! band loop

    close(fid)

    !---< aerosol(particle) parameter input >---

    ! allocate arrays
    allocate( fitPLK  (MSTRN_nfitPLK,MSTRN_nband) )
    allocate( fsol    (              MSTRN_nband) )
    allocate( sfc     (MSTRN_nsfc,   MSTRN_nband) )
    allocate( rayleigh(              MSTRN_nband) )

    allocate( qmol    (                                MSTRN_nmoment,MSTRN_nband) )
    allocate( q       (-1:MSTRN_nradius+1,MSTRN_nptype,MSTRN_nmoment,MSTRN_nband) )
    do iptype = 1, MSTRN_nptype
       nradius = ptype_nradius(iptype)

       q(       -1:0              ,iptype,:,:) = 0.D0 ! dummy for NaN
       q(nradius+1:MSTRN_nradius+1,iptype,:,:) = 0.D0 ! dummy for extrapolation
    enddo


    call IO_get_fname(fname, MSTRN_AEROPARA_INPUTFILE)
    open( fid,                  &
          file   = fname,       &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Input data file does not found! ', trim(fname)
          call PRC_abort
       endif

       ! read aerosol/surface parameters for check
       read(fid,*) nband, nsfc, nptype, nstream, nplkord, nfitPLK

       if ( nband /= MSTRN_nband ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nband(given,file)=', MSTRN_nband, nband
          call PRC_abort
       endif
       if ( nsfc /= MSTRN_nsfc ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nsfc(given,file)=', MSTRN_nsfc, nsfc
          call PRC_abort
       endif
       if ( nptype /= MSTRN_nptype ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nptype(given,file)=', MSTRN_nptype, nptype
          call PRC_abort
       endif
       if ( nstream /= MSTRN_nstream ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nstream(given,file)=', MSTRN_nstream, nstream
          call PRC_abort
       endif
       if ( nplkord /= MSTRN_nplkord ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nplkord(given,file)=', MSTRN_nplkord, nplkord
          call PRC_abort
       endif
       if ( nfitPLK /= MSTRN_nfitPLK ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nfitPLK(given,file)=', MSTRN_nfitPLK, nfitPLK
          call PRC_abort
       endif

       ! wave band boundaries
       read(fid,*) dummy
       read(fid,*) waveh(:)

       ! for each band
       do iw = 1, MSTRN_nband

          ! fitting point for planck functions
          read(fid,*) dummy
          read(fid,*) fitPLK(:,iw)
          ! solar insolation
          read(fid,*) dummy
          read(fid,*) fsol(iw)
          ! surface properties
          read(fid,*) dummy
          read(fid,*) sfc(:,iw)
          ! rayleigh scattering
          read(fid,*) dummy
          read(fid,*) rayleigh(iw)

          ! moments
          read(fid,*) dummy
          do im = 1, MSTRN_nmoment
             ! for rayleigh scattering phase function
             read(fid,*) qmol(im,iw)
             ! for aerosol scattering phase function
             do iptype = 1, MSTRN_nptype
                nradius = ptype_nradius(iptype)
                read(fid,*) q(1:nradius,iptype,im,iw)
             enddo
          enddo

       enddo

    close(fid)

    fsol_tot = 0.0_RP
    do iw = 1, MSTRN_nband
       fsol_tot = fsol_tot + fsol(iw)
    enddo

    !---< radius mode & hygroscopic parameter input >---

    ! allocate arrays
    allocate( hygro_flag(MSTRN_nptype)               )
    allocate( radmode   (MSTRN_nptype,MSTRN_nradius) )

    call IO_get_fname(fname, MSTRN_HYGROPARA_INPUTFILE)
    open( fid,                  &
          file   = fname,       &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Input data file does not found! ', trim(fname)
          call PRC_abort
       endif

       read(fid,*) nptype

       if ( nptype /= MSTRN_nptype ) then
          LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nptype(given,file)=', MSTRN_nptype, nptype
          call PRC_abort
       endif

       do iptype = 1, MSTRN_nptype
          read(fid,*) dummy
          read(fid,*) hygro_flag(iptype), nradius

          if ( nradius /= ptype_nradius(iptype) ) then
             LOG_ERROR("RD_MSTRN_setup",*) 'Inconsistent parameter value! nradius(given,file)=', ptype_nradius(iptype), nradius
             call PRC_abort
          endif

          read(fid,*) radmode(iptype,1:nradius)
       enddo

    close(fid)

    !----- report data -----
    LOG_NEWLINE
    LOG_INFO("RD_MSTRN_setup",'(1x,A,F12.7)') 'Baseline of total solar insolation : ', fsol_tot

    !---< constant parameter for main scheme >---
    RHO_std = Pstd / ( Rdry * TEM00 ) ! [kg/m3]

    M   (I_SW) = 1.0_RP / sqrt(3.0_RP)
    W   (I_SW) = 1.0_RP
    Wbar(I_SW) = 1.0_RP

    M   (I_LW) = 1.0_RP / 1.66_RP
    W   (I_LW) = 1.0_RP
    Wbar(I_LW) = 2.0_RP * M(I_LW)

    do im = 1, 2
       Wmns  (im) = sqrt( W(im) / M(im) )
       Wpls  (im) = sqrt( W(im) * M(im) )
       Wscale(im) = Wpls(im) / Wbar(im)
    end do
    !$acc update device(M, W, Wmns, Wpls, Wscale)

    !$acc enter data &
    !$acc copyin(waveh, &
    !$acc        wgtch, fitPLK, logfitP, logfitT, fitT, &
    !$acc        radmode, ngasabs, igasabs, ptype_nradius, &
    !$acc        fsol, q, qmol, rayleigh, acfc_pow, nch, AKD, SKD )

    return
  end subroutine RD_MSTRN_setup

  !-----------------------------------------------------------------------------
  !> DTRN v3.2
!OCL SERIAL
  subroutine RD_MSTRN_DTRN3( &
       rd_kmax, IA, IS, IE, JA, JS, JE, &
#ifndef _OPENACC
       i, j, &
#endif
       ngas,         &
       ncfc,         &
       naero,        &
       hydro_str,    &
       hydro_end,    &
       aero_str,     &
       aero_end,     &
       solins,       &
       cosSZA0,      &
       rhodz,        &
       pres,         &
       temp,         &
       temph,        &
       temp_sfc,     &
       gas,          &
       cfc,          &
       aerosol_conc, &
       aerosol_radi, &
       aero2ptype,   &
       cldfrac,      &
       albedo_sfc,   &
       fact_ocean,   &
       fact_land,    &
       fact_urban,   &
       rflux,        &
       rflux_sfc_dn, &
       tauCLD_067u,  &
       emisCLD_105u  )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       GRAV => CONST_GRAV, &
       Pstd => CONST_Pstd, &
       PPM  => CONST_PPM
    implicit none
    integer, intent(in) :: rd_kmax
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
#ifndef _OPENACC
    integer, intent(in) :: i, j
#endif

    integer,  intent(in)  :: ngas
    integer,  intent(in)  :: ncfc
    integer,  intent(in)  :: naero
    integer,  intent(in)  :: hydro_str
    integer,  intent(in)  :: hydro_end
    integer,  intent(in)  :: aero_str
    integer,  intent(in)  :: aero_end
    real(RP), intent(in)  :: solins      (IA,JA)
    real(RP), intent(in)  :: cosSZA0     (IA,JA)
    real(RP), intent(in)  :: rhodz       (rd_kmax  ,IA,JA)
    real(RP), intent(in)  :: pres        (rd_kmax  ,IA,JA)
    real(RP), intent(in)  :: temp        (rd_kmax  ,IA,JA)
    real(RP), intent(in)  :: temph       (rd_kmax+1,IA,JA)
    real(RP), intent(in)  :: temp_sfc    (IA,JA)
    real(RP), intent(in)  :: gas         (rd_kmax,IA,JA,ngas )
    real(RP), intent(in)  :: cfc         (rd_kmax,IA,JA,ncfc )
    real(RP), intent(in)  :: aerosol_conc(rd_kmax,IA,JA,naero)
    real(RP), intent(in)  :: aerosol_radi(rd_kmax,IA,JA,naero)
    integer,  intent(in)  :: aero2ptype  (naero)
    real(RP), intent(in)  :: cldfrac     (rd_kmax,IA,JA)
    real(RP), intent(in)  :: albedo_sfc  (IA,JA,N_RAD_DIR,N_RAD_RGN)
    real(RP), intent(in)  :: fact_ocean  (IA,JA)
    real(RP), intent(in)  :: fact_land   (IA,JA)
    real(RP), intent(in)  :: fact_urban  (IA,JA)
    real(RP), intent(inout) :: rflux       (rd_kmax+1,IA,JA,2,2,MSTRN_ncloud)
    real(RP), intent(inout) :: rflux_sfc_dn(IA,JA,N_RAD_DIR,N_RAD_RGN)        ! surface downward radiation flux (direct/diffuse,IR/NIR/VIS)
    real(RP), intent(inout) :: tauCLD_067u (rd_kmax,IA,JA)                    ! 0.67 micron cloud optical depth
    real(RP), intent(inout) :: emisCLD_105u(rd_kmax,IA,JA)                    ! 10.5 micron cloud emissivity

#ifdef _OPENACC
#define LOOP_INNER do j=JS,JE; do i=IS,IE
#define LOOP_END_INNER end do; end do
#else
#define LOOP_INNER
#define LOOP_END_INNER
#endif

    ! for P-T fitting
    real(RP) :: dz_std (rd_kmax,IS:IE,JS:JE)       ! layer thickness at 0C, 1atm [cm]
    real(RP) :: logP                                     ! log10(pres)
    real(RP) :: logT                                     ! log10(temp)
    integer  :: indexP (rd_kmax,IS:IE,JS:JE)       ! index for interpolation in P-fitting
    real(RP) :: factP  (rd_kmax,IS:IE,JS:JE)       ! interpolation factor    in P-fitting
    real(RP) :: factT32(rd_kmax,IS:IE,JS:JE)       ! interpolation factor    in T-fitting
    real(RP) :: factT21(rd_kmax,IS:IE,JS:JE)       ! interpolation factor    in T-fitting
    integer  :: indexR (rd_kmax,naero,IS:IE,JS:JE) ! index for interpolation in R-fitting
    real(RP) :: factR  (rd_kmax,naero,IS:IE,JS:JE) ! interpolation factor    in R-fitting

    ! for optical thickness by gas
    real(RP) :: tauGAS(rd_kmax,MSTRN_ch_limit,IS:IE,JS:JE) ! optical thickness by gas absorption (total)
    real(RP) :: A1, A2, A3, factPT
    real(RP) :: qv, length
    integer  :: gasno

    ! for optical thickness by particles
    real(RP) :: tauPR   (rd_kmax,MSTRN_ncloud)               ! optical thickness        by Rayleigh/cloud/aerosol
    real(RP) :: omgPR   (rd_kmax,MSTRN_ncloud)               ! single scattering albedo by Rayleigh/cloud/aerosol
    real(DP) :: optparam(rd_kmax,MSTRN_nmoment,MSTRN_ncloud,IS:IE,JS:JE) ! optical parameters
    real(RP) :: q_fit, dp_P

    ! for planck functions
    real(RP) :: bbar (rd_kmax  ) ! planck functions for thermal source at the interface
    real(RP) :: bbarh(rd_kmax+1) ! planck functions for thermal source at the center
    real(RP) :: b_sfc            ! planck functions for thermal source at the surface
    real(RP) :: wl, beta

    ! for two-stream
    real(RP) :: tau(rd_kmax,    MSTRN_ncloud,MSTRN_ch_limit) ! total optical thickness
    real(RP) :: omg(rd_kmax,    MSTRN_ncloud,MSTRN_ch_limit) ! single scattering albedo
    real(RP) :: g  (rd_kmax,0:2,MSTRN_ncloud               ) ! two-stream approximation factors
                                                             ! 0: always 1
                                                             ! 1: asymmetry factor
                                                             ! 2: truncation factor
    real(RP) :: b  (rd_kmax,0:2,MSTRN_ncloud,MSTRN_ch_limit) ! planck expansion coefficients (zero if SW)
    real(RP) :: fsol_rgn                                     ! solar insolation              (zero if LW)

    real(RP) :: flux       (rd_kmax+1,2,MSTRN_ncloud,MSTRN_ch_limit) ! upward/downward flux
    real(RP) :: flux_direct(rd_kmax+1,  MSTRN_ncloud,MSTRN_ch_limit) ! downward flux (direct solar)

    ! for satellite simulator
    real(RP) :: emisCLD(rd_kmax,MSTRN_ch_limit) ! cloud emissivity

    ! work
    real(RP) :: cosSZA(IA,JA)

    ! works for two_stream
    ! main factors
    real(RP) :: Tdir0(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! transmission factor for solar direct (clear-sky/cloud)
    real(RP) :: R0   (rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! reflection   factor                  (clear-sky/cloud)
    real(RP) :: T0   (rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! transmission factor                  (clear-sky/cloud)
    real(RP) :: Em_LW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! thermal source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Em_SW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! solar   source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Ep_LW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! thermal source (TOA->sfc)            (clear-sky/cloud)
    real(RP) :: Ep_SW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! solar   source (TOA->sfc)            (clear-sky/cloud)
    ! Averaged factors, considering cloud overwrap
    real(RP) :: cf         (rd_kmax  ,MSTRN_ncloud*MSTRN_ch_limit)       ! cloud fraction
    real(RP) :: tau_bar_sol(rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit)       ! solar insolation through accumulated optical thickness at each layer
    real(RP) :: R          (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! reflection   factor
    real(RP) :: T          (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! transmission factor
    real(RP) :: Em         (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! source (sfc->TOA)
    real(RP) :: Ep         (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! source (TOA->sfc)
    real(RP) :: Em0        (MSTRN_ncloud)
    ! Doubling-Adding
    real(RP) :: R12mns(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! reflection factor in doubling method
    real(RP) :: R12pls(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! reflection factor in doubling method
    real(RP) :: E12mns(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! source function   in doubling method
    real(RP) :: E12pls(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! source function   in doubling method


    real(RP) :: zerosw
    real(RP) :: valsum(rd_kmax)
    integer  :: chmax
    integer  :: ip, ir, irgn, irgn_alb
    integer  :: igas, icfc, iaero, iptype, nradius
    integer  :: iw, ich, iplk, icloud, im, idir
    integer  :: k
#ifdef _OPENACC
    integer :: i, j
#endif
    !---------------------------------------------------------------------------

    !$acc data &
    !$acc copyin(solins, cosSZA0, rhodz, pres, temp, temph, temp_sfc, gas, cfc, &
    !$acc        aerosol_conc, aerosol_radi, aero2ptype, cldfrac, albedo_sfc, &
    !$acc        fact_ocean, fact_land, fact_urban) &
    !$acc copyout(rflux, rflux_sfc_dn, tauCLD_067u, emisCLD_105u) &
    !$acc create(dz_std, indexP, factP, factT32, factT21, indexR, factR, cosSZA, optparam, tauGAS)

    !$acc kernels
    LOOP_INNER
       cosSZA(i,j) = max( cosSZA0(i,j), RD_cosSZA_min )
    LOOP_END_INNER
    !$acc end kernels

    !$acc parallel loop collapse(2)
    LOOP_INNER
    do k = 1, rd_kmax
       dz_std(k,i,j) = rhodz(k,i,j) / RHO_std * 100.0_RP ! [cm]

       logP = log10( pres(k,i,j) )
       logT = log10( temp(k,i,j) )

       indexP(k,i,j) = MSTRN_nfitP
       !$acc loop seq
       do ip = MSTRN_nfitP, 2, -1
          if( logP >= logfitP(ip) ) indexP(k,i,j) = ip
       enddo

       ip = indexP(k,i,j)

       factP(k,i,j) = ( logP - logfitP(ip-1) ) / ( logfitP(ip) - logfitP(ip-1) )

       factT32(k,i,j) = ( logT        - logfitT(2)  ) / ( logfitT(3) - logfitT(2) ) &
                      * ( temp(k,i,j) - fitT(1)     ) / ( fitT(3)    - fitT(1)    )
       factT21(k,i,j) = ( logT        - logfitT(2)  ) / ( logfitT(2) - logfitT(1) ) &
                      * ( fitT(3)     - temp(k,i,j) ) / ( fitT(3)    - fitT(1)    )
    enddo
    LOOP_END_INNER
    !$acc end parallel

    !---< interpolation of mode radius & hygroscopic parameter (R-fitting) >---
    !$acc parallel loop collapse(3)
    LOOP_INNER
    do iaero = 1, naero
    do k = 1, rd_kmax
       iptype  = aero2ptype(iaero)
       nradius = ptype_nradius(iptype)

       ! [Note] Extrapolation sometimes makes unexpected value
       if ( aerosol_radi(k,i,j,iaero) <= radmode(iptype,1) ) then
          indexR(k,iaero,i,j) = 1
          factR (k,iaero,i,j) = 0.0_RP
       elseif( aerosol_radi(k,i,j,iaero) > radmode(iptype,nradius) ) then
          indexR(k,iaero,i,j) = nradius
          factR (k,iaero,i,j) = 1.0_RP
       else
          indexR(k,iaero,i,j) = -1
          !$acc loop seq
          do ir = 1, nradius-1
             if (       aerosol_radi(k,i,j,iaero) <= radmode(iptype,ir+1) &
                  .AND. aerosol_radi(k,i,j,iaero) >  radmode(iptype,ir  ) ) then ! interpolation
                indexR(k,iaero,i,j) = ir
                factR (k,iaero,i,j) = ( aerosol_radi(k,i,j,iaero) - radmode(iptype,ir) ) &
                                    / ( radmode(iptype,ir+1)      - radmode(iptype,ir) )
                exit
             endif
          enddo
       endif
    enddo
    enddo
    LOOP_END_INNER
    !$acc end parallel

    ! initialize
    !$acc parallel loop collapse(5)
    do icloud = 1, MSTRN_ncloud
    do idir = 1, 2
    do irgn = 1, 2
    LOOP_INNER
    do k = 1, rd_kmax+1
      rflux(k,i,j,irgn,idir,icloud) = 0.0_RP
    end do
    LOOP_END_INNER
    end do
    end do
    end do
    !$acc end parallel
    !$acc kernels
    do irgn_alb = 1, N_RAD_RGN
    do idir = 1, N_RAD_DIR
    LOOP_INNER
       rflux_sfc_dn(i,j,idir,irgn_alb) = 0.0_RP
    LOOP_END_INNER
    end do
    end do
    !$acc end kernels

    do iw = 1, MSTRN_nband
       if ( iflgb(I_SWLW,iw) == 0 ) then ! Long wave region
          irgn     = I_LW
          irgn_alb = I_R_IR ! IR
       else ! Short wave region
          irgn = I_SW
          if ( waveh(iw) >= 10000.0_RP / 0.7_RP ) then
             irgn_alb = I_R_VIS ! Visible
          else
             irgn_alb = I_R_NIR ! Near-IR
          endif
       endif

       chmax = nch(iw)

       !---< interpolation of gas parameters (P-T fitting) >---
       !$acc parallel loop collapse(2)
       LOOP_INNER
       !$acc loop collapse(2)
       do ich = 1, chmax
       do k = 1, rd_kmax
          tauGAS(k,ich,i,j) = 0.0_RP
       enddo
       enddo
       LOOP_END_INNER
       !$acc end parallel

       !--- Gas line absorption
       !$acc parallel loop collapse(2)
       LOOP_INNER
       !$acc loop seq
       do igas = 1, ngasabs(iw)
          gasno = igasabs(igas,iw)
          !$acc loop collapse(2)
          do ich = 1, chmax
          do k = 1, rd_kmax
             ip = indexP(k,i,j)
             A1 = AKD(ip-1,ich,1,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                + AKD(ip  ,ich,1,gasno,iw) * (          factP(k,i,j) )
             A2 = AKD(ip-1,ich,2,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                + AKD(ip  ,ich,2,gasno,iw) * (          factP(k,i,j) )
             A3 = AKD(ip-1,ich,3,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                + AKD(ip  ,ich,3,gasno,iw) * (          factP(k,i,j) )
             factPT = factT32(k,i,j)*(A3-A2) + A2 + factT21(k,i,j)*(A2-A1)

             length = gas(k,i,j,igasabs(igas,iw)) * PPM * dz_std(k,i,j)
             tauGAS(k,ich,i,j) = tauGAS(k,ich,i,j) + 10.0_RP**factPT * length
          end do
          end do
       end do
       LOOP_END_INNER
       !$acc end parallel

       !--- Gas broad absorption
       if ( iflgb(I_H2O_continuum,iw) == 1 ) then
          !$acc parallel loop collapse(2)
          LOOP_INNER
          !$acc loop collapse(2)
          do ich = 1, chmax
          do k = 1, rd_kmax
             ip = indexP(k,i,j)
             A1 = SKD(ip-1,ich,1,iw) * ( 1.0_RP-factP(k,i,j) )&
                + SKD(ip  ,ich,1,iw) * (        factP(k,i,j) )
             A2 = SKD(ip-1,ich,2,iw) * ( 1.0_RP-factP(k,i,j) )&
                + SKD(ip  ,ich,2,iw) * (        factP(k,i,j) )
             A3 = SKD(ip-1,ich,3,iw) * ( 1.0_RP-factP(k,i,j) )&
                + SKD(ip  ,ich,3,iw) * (        factP(k,i,j) )
             factPT = factT32(k,i,j)*(A3-A2) + A2 + factT21(k,i,j)*(A2-A1)

             qv = gas(k,i,j,1) * PPM * dz_std(k,i,j)
             length = qv*qv / ( qv + dz_std(k,i,j) )
             tauGAS(k,ich,i,j) = tauGAS(k,ich,i,j) + 10.0_RP**factPT * length
          enddo
          enddo
          LOOP_END_INNER
          !$acc end parallel
       endif

       if ( iflgb(I_CFC_continuum,iw) == 1 ) then
          !$acc parallel loop collapse(2) private(valsum)
          LOOP_INNER
          do k = 1, rd_kmax
             valsum(k) = 0.0_RP
          end do
          !$acc loop seq
          do icfc = 1, ncfc
          do k = 1, rd_kmax
             valsum(k) = valsum(k) + acfc_pow(icfc,iw) * cfc(k,i,j,icfc)
          end do
          end do
          do k = 1, rd_kmax
             valsum(k) = valsum(k) * PPM * dz_std(k,i,j)
          end do

          do ich = 1, chmax
          do k = 1, rd_kmax
             tauGAS(k,ich,i,j) = tauGAS(k,ich,i,j) + valsum(k)
          enddo
          enddo
          LOOP_END_INNER
          !$acc end parallel
       endif

       !---< particle (Rayleigh/Cloud/Aerosol) >---

       ! optical thickness, phase function
       ! im=1: extinction coefficient
       ! im=2,3,4: moments of the volume scattering phase function

       !--- Rayleigh scattering
       !$acc parallel loop collapse(2)
       LOOP_INNER
       !$acc loop collapse(2)
       do im = 1, MSTRN_nstream*2+2
       do k = 1, rd_kmax
          dp_P = rhodz(k,i,j) * GRAV / Pstd
          length = rayleigh(iw) * dp_P

          optparam(k,im,I_Cloud   ,i,j) = qmol(im,iw) * length
          optparam(k,im,I_ClearSky,i,j) = qmol(im,iw) * length
       enddo
       enddo
       LOOP_END_INNER
       !$acc end parallel

       !--- Cloud scattering
       !$acc parallel loop collapse(2)
       LOOP_INNER
       do k = 1, rd_kmax
          tauCLD_067u(k,i,j) = 0.0_RP
       enddo
       LOOP_END_INNER
       !$acc end parallel
       !$acc kernels
       LOOP_INNER
       !$acc loop seq
       do iaero = hydro_str, hydro_end
          iptype = aero2ptype(iaero)

          do im = 1, MSTRN_nstream*2+2
          do k = 1, rd_kmax
             ir = indexR(k,iaero,i,j)
             length = aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)

             q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(k,iaero,i,j) ) &
                   + q(ir+1,iptype,im,iw) * (        factR(k,iaero,i,j) )
             q_fit = max( q_fit, 0.0_RP )

             optparam(k,im,I_Cloud,i,j) = optparam(k,im,I_Cloud,i,j) + q_fit * length
          enddo
          enddo

          if ( waveh(iw) <= 1.493E+4_RP .AND. 1.493E+4_RP < waveh(iw+1) ) then ! 0.67 micron
             do k = 1, rd_kmax
                ir = indexR(k,iaero,i,j)
                length = aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)

                q_fit = q(ir  ,iptype,1,iw) * ( 1.0_RP-factR(k,iaero,i,j) ) &
                      + q(ir+1,iptype,1,iw) * (        factR(k,iaero,i,j) )
                q_fit = max( q_fit, 0.0_RP )

                tauCLD_067u(k,i,j) = tauCLD_067u(k,i,j) + q_fit * length
             end do
          end if

       enddo
       LOOP_END_INNER
       !$acc end kernels

       !--- Aerosol scattering
       !$acc kernels
       LOOP_INNER
       !$acc loop seq
       do iaero = aero_str, aero_end
          iptype = aero2ptype(iaero)

          do im = 1, MSTRN_nstream*2+2
          do k = 1, rd_kmax
             ir = indexR(k,iaero,i,j)
             length = aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)

             q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(k,iaero,i,j) ) &
                   + q(ir+1,iptype,im,iw) * (        factR(k,iaero,i,j) )
             q_fit = max( q_fit, 0.0_RP )

             optparam(k,im,I_Cloud   ,i,j) = optparam(k,im,I_Cloud   ,i,j) + q_fit * length
             optparam(k,im,I_ClearSky,i,j) = optparam(k,im,I_ClearSky,i,j) + q_fit * length
          enddo
          enddo
       enddo
       LOOP_END_INNER
       !$acc end kernels

       !$acc kernels
       !$acc loop independent collapse(2) &
       !$acc private(tauPR,omgPR,g,b,b_sfc,fsol_rgn,bbar,bbarh,tau,omg, &
       !$acc         flux,flux_direct,emisCLD, &
       !$acc         Tdir0,R0,T0,Em_LW,Em_SW,Ep_LW,Ep_SW,cf,tau_bar_sol,R,T,Em,Ep,Em0, &
       !$acc         R12mns,R12pls,E12mns,E12pls )
       LOOP_INNER

       !$acc loop collapse(2)
       do icloud = 1, MSTRN_ncloud
       do k = 1, rd_kmax
          tauPR(k,icloud) = optparam(k,1,icloud,i,j)
          omgPR(k,icloud) = optparam(k,1,icloud,i,j) - optparam(k,2,icloud,i,j)

          !--- g
          zerosw = 0.5_RP - sign(0.5_RP,omgPR(k,icloud)-RD_EPS)

          g(k,0,icloud) = 1.0_RP
          g(k,1,icloud) = optparam(k,3,icloud,i,j) * ( 1.0_RP-zerosw ) / ( omgPR(k,icloud)+zerosw )
          g(k,2,icloud) = optparam(k,4,icloud,i,j) * ( 1.0_RP-zerosw ) / ( omgPR(k,icloud)+zerosw )
          !g(k,1,icloud) = max( optparam(k,3,icloud,i,j) * ( 1.0_RP-zerosw ) / ( omgPR(k,icloud)+zerosw ), 0.0_RP )
          !g(k,2,icloud) = max( optparam(k,4,icloud,i,j) * ( 1.0_RP-zerosw ) / ( omgPR(k,icloud)+zerosw ), 0.0_RP )
       enddo
       enddo

       if ( irgn == I_SW ) then ! solar
          b(:,:,:,:) = 0.0_RP

          b_sfc = 0.0_RP
          fsol_rgn = fsol(iw) / fsol_tot * solins(i,j)

       elseif( irgn == I_LW ) then ! IR
          !--- set planck functions
          wl = 10000.0_RP / sqrt( waveh(iw) * waveh(iw+1) )

          ! from temp at cell center
          do k = 1, rd_kmax
             beta = 0.0_RP
             !$acc loop seq
             do iplk = MSTRN_nfitPLK, 1, -1
                beta = beta / ( wl*temp(k,i,j) ) + fitPLK(iplk,iw)
             enddo
             bbar(k) = exp(-beta) * temp(k,i,j) / (wl*wl)
          enddo

          ! from temp at cell wall
          do k = 1, rd_kmax+1
             beta = 0.0_RP
             !$acc loop seq
             do iplk = MSTRN_nfitPLK, 1, -1
                beta = beta / ( wl*temph(k,i,j) ) + fitPLK(iplk,iw)
             enddo
             bbarh(k) = exp(-beta) * temph(k,i,j) / (wl*wl)
          enddo

          ! from temp_sfc
          beta = 0.0_RP
          !$acc loop seq
          do iplk = MSTRN_nfitPLK, 1, -1
             beta = beta / ( wl*temp_sfc(i,j) ) + fitPLK(iplk,iw)
          enddo

          b_sfc = exp(-beta) * temp_sfc(i,j) / (wl*wl)
          fsol_rgn = 0.0_RP

       endif

       !--- total tau & omega
       !$acc loop collapse(3)
       do ich = 1, chmax
       do icloud = 1, MSTRN_ncloud
       do k = 1, rd_kmax
          tau(k,icloud,ich) = tauGAS(k,ich,i,j) + tauPR(k,icloud)
          !tau(k,icloud,ich) = max( tauGAS(k,ich,i,j) + tauPR(k,icloud), 0.0_RP )

          zerosw = 0.5_RP - sign( 0.5_RP, tau(k,icloud,ich)-RD_EPS ) ! if tau < EPS, zerosw = 1

          omg(k,icloud,ich) = ( 1.0_RP-zerosw ) * omgPR(k,icloud) / ( tau(k,icloud,ich)-zerosw ) &
                                + (        zerosw ) * 1.0_RP

          !omg(k,icloud,ich) = min( max( omg(k,icloud,ich), 0.0_RP ), 1.0_RP )
       enddo
       enddo
       enddo

       if( irgn == I_LW ) then ! IR
          !$acc loop collapse(3)
          do ich = 1, chmax
          do icloud = 1, MSTRN_ncloud
          do k = 1, rd_kmax
             zerosw = 0.5_RP - sign( 0.5_RP, tau(k,icloud,ich)-RD_EPS ) ! if tau < EPS, zerosw = 1

             b(k,0,icloud,ich) = bbarh(k)
             b(k,1,icloud,ich) = ( 1.0_RP-zerosw )       &
                               * ( -          bbarh(k+1) &
                                   + 4.0_RP * bbar (k  ) &
                                   - 3.0_RP * bbarh(k  ) &
                                 ) / ( tau(k,icloud,ich)-zerosw )
             b(k,2,icloud,ich) = ( 1.0_RP-zerosw )       &
                               * ( +          bbarh(k+1) &
                                   - 2.0_RP * bbar (k  ) &
                                   +          bbarh(k  ) &
                                 ) / ( tau(k,icloud,ich)**2-zerosw ) * 2.0_RP
          enddo
          enddo
          enddo

       endif ! solar/IR switch


       ! two-stream transfer
       !call PROF_rapstart('RD_MSTRN_twst', 3)
       call RD_MSTRN_two_stream( IA, IS, IE, IA, JS, JE, &
                                 RD_KMAX, chmax,                        & ! [IN]
                                 cosSZA(i,j),                           & ! [IN]
                                 fsol_rgn, irgn,                        & ! [IN]
                                 tau(:,:,:), omg(:,:,:),                & ! [IN]
                                 g(:,:,:), b(:,:,:,:),  b_sfc,          & ! [IN]
                                 albedo_sfc (i,j,I_R_direct, irgn_alb), & ! [IN]
                                 albedo_sfc (i,j,I_R_diffuse,irgn_alb), & ! [IN]
                                 cldfrac(:,i,j),                        & ! [IN]
                                 flux(:,:,:,:), flux_direct(:,:,:),     & ! [OUT]
                                 emisCLD(:,:),                          & ! [OUT]
                                 Tdir0(:,:,:), R0(:,:,:), T0(:,:,:),    & ! [WORK]
                                 Em_LW(:,:,:), Em_SW(:,:,:),            & ! [WORK]
                                 Ep_LW(:,:,:), Ep_SW(:,:,:),            & ! [WORK]
                                 cf(:,:), tau_bar_sol(:,:),             & ! [WORK]
                                 R(:,:), T(:,:),                        & ! [WORK]
                                 Em(:,:), Ep(:,:), Em0(:),              & ! [WORK]
                                 R12mns(:,:), R12pls(:,:),              & ! [WORK]
                                 E12mns(:,:), E12pls(:,:)               ) ! [WORK]


       !call PROF_rapend  ('RD_MSTRN_twst', 3)

       !$acc loop seq
       do ich = 1, chmax
          !$acc loop collapse(2)
          do icloud = 1, MSTRN_ncloud
             do k = 1, rd_kmax+1
                rflux(k,i,j,irgn,I_up,icloud) = rflux(k,i,j,irgn,I_up,icloud) + flux(k,I_up,icloud,ich) * wgtch(ich,iw)
                rflux(k,i,j,irgn,I_dn,icloud) = rflux(k,i,j,irgn,I_dn,icloud) + flux(k,I_dn,icloud,ich) * wgtch(ich,iw)
             enddo
          enddo
       enddo

       !$acc loop seq
       do ich = 1, chmax
          rflux_sfc_dn(i,j,I_R_direct ,irgn_alb) = rflux_sfc_dn(i,j,I_R_direct ,irgn_alb) &
                                                 + (                                    flux_direct(rd_kmax+1,I_Cloud,ich) ) * wgtch(ich,iw)
          rflux_sfc_dn(i,j,I_R_diffuse,irgn_alb) = rflux_sfc_dn(i,j,I_R_diffuse,irgn_alb) &
                                                 + ( flux(rd_kmax+1,I_dn,I_Cloud,ich) - flux_direct(rd_kmax+1,I_Cloud,ich) ) * wgtch(ich,iw)
       end do

       if ( waveh(iw) <= 952.0_RP .AND. 952.0_RP < waveh(iw+1) ) then ! 10.5 micron
          ! 10.5 micron emissivity for resolved clouds

          do k = 1, rd_kmax
             emisCLD_105u(k,i,j) = 0.0_RP
          end do
          !$acc loop seq
          do ich = 1, chmax
             do k = 1, rd_kmax
                emisCLD_105u(k,i,j) = emisCLD_105u(k,i,j) + emisCLD(k,ich) * wgtch(ich,iw)
             enddo
          end do
       endif

       LOOP_END_INNER
       !$acc end kernels

    enddo ! IW loop

    !$acc end data

    return
  end subroutine RD_MSTRN_DTRN3

  !-----------------------------------------------------------------------------
  !> Two stream calculation CORE
!OCL SERIAL
  subroutine RD_MSTRN_two_stream( &
       IA, IS, IE, JA, JS, JE,  &
       RD_KMAX, chmax,     &
       cosSZA,             &
       fsol, irgn,         &
       tau, omg,           &
       g, b, b_sfc,        &
       albedo_sfc_direct,  &
       albedo_sfc_diffuse, &
       cldfrac,            &
       flux, flux_direct,  &
       emisCLD,            &
       Tdir0, R0, T0,      &
       Em_LW, Em_SW,       &
       Ep_LW, Ep_SW,       &
       cf, tau_bar_sol,    &
       R, T,               &
       Em, Ep, Em0,        &
       R12mns, R12pls,     &
       E12mns, E12pls      )
    !$acc routine vector
    use scale_const, only: &
       PI   => CONST_PI,   &
       EPS  => CONST_EPS,  &
       EPS1 => CONST_EPS1
    implicit none

    integer,  intent(in)  :: IA, IS, IE, JA, JS, JE
    integer,  intent(in)  :: RD_KMAX
    integer,  intent(in), value :: chmax
    real(RP), intent(in)        :: cosSZA                                         ! cos(SZA) = mu0
    real(RP), intent(in), value :: fsol                                           ! solar radiation intensity
    integer,  intent(in), value  :: irgn                                          ! 1:LW 2:SW
    real(RP), intent(in)  :: tau        (rd_kmax,    MSTRN_ncloud,MSTRN_ch_limit) ! total optical thickness          (clear-sky/cloud)
    real(RP), intent(in)  :: omg        (rd_kmax,    MSTRN_ncloud,MSTRN_ch_limit) ! single scattering albedo         (clear-sky/cloud)
    real(RP), intent(in)  :: g          (rd_kmax,0:2,MSTRN_ncloud               ) ! two-stream approximation factors (clear-sky/cloud)
    real(RP), intent(in)  :: b          (rd_kmax,0:2,MSTRN_ncloud,MSTRN_ch_limit) ! planck expansion coefficients    (clear-sky/cloud)
    real(RP), intent(in), value  :: b_sfc                                         ! planck function at surface
    real(RP), intent(in)  :: albedo_sfc_direct                                    ! surface albedo (DIRECT)
    real(RP), intent(in)  :: albedo_sfc_diffuse                                   ! surface albedo (DIFFUSE)
    real(RP), intent(in)  :: cldfrac    (rd_kmax)                                 ! cloud fraction
    real(RP), intent(out) :: flux       (rd_kmax+1,2,MSTRN_ncloud,MSTRN_ch_limit) ! upward(sfc->TOA)/downward(TOA->sfc) flux (clear-sky/cloud)
    real(RP), intent(out) :: flux_direct(rd_kmax+1,  MSTRN_ncloud,MSTRN_ch_limit) ! downward(TOA->sfc) flux, solar direct    (clear-sky/cloud)
    real(RP), intent(out) :: emisCLD    (rd_kmax,                 MSTRN_ch_limit) ! cloud emissivity factor (cloud)

    ! work
    ! main factors
    real(RP), intent(out) :: Tdir0(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! transmission factor for solar direct (clear-sky/cloud)
    real(RP), intent(out) :: R0   (rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! reflection   factor                  (clear-sky/cloud)
    real(RP), intent(out) :: T0   (rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! transmission factor                  (clear-sky/cloud)
    real(RP), intent(out) :: Em_LW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! thermal source (sfc->TOA)            (clear-sky/cloud)
    real(RP), intent(out) :: Em_SW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! solar   source (sfc->TOA)            (clear-sky/cloud)
    real(RP), intent(out) :: Ep_LW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! thermal source (TOA->sfc)            (clear-sky/cloud)
    real(RP), intent(out) :: Ep_SW(rd_kmax,MSTRN_ncloud,MSTRN_ch_limit) ! solar   source (TOA->sfc)            (clear-sky/cloud)
    ! Averaged factors, considering cloud overwrap
    real(RP), intent(out) :: cf         (rd_kmax  ,MSTRN_ncloud*MSTRN_ch_limit)       ! cloud fraction
    real(RP), intent(out) :: tau_bar_sol(rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit)       ! solar insolation through accumulated optical thickness at each layer
    real(RP), intent(out) :: R          (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! reflection   factor
    real(RP), intent(out) :: T          (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! transmission factor
    real(RP), intent(out) :: Em         (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! source (sfc->TOA)
    real(RP), intent(out) :: Ep         (rd_kmax+1,MSTRN_ncloud*MSTRN_ch_limit) ! source (TOA->sfc)
    real(RP), intent(out) :: Em0        (MSTRN_ncloud)
    ! Doubling-Adding
    real(RP), intent(out) :: R12mns(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! reflection factor in doubling method
    real(RP), intent(out) :: R12pls(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! reflection factor in doubling method
    real(RP), intent(out) :: E12mns(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! source function   in doubling method
    real(RP), intent(out) :: E12pls(MSTRN_ncloud*MSTRN_ch_limit,rd_kmax+1) ! source function   in doubling method


    ! parameters with two-stream truncation
    real(RP) :: tau_new    ! optical thickness        : two-stream truncation
    real(RP) :: omg_new    ! single scattering albedo : two-stream truncation
    real(RP) :: g_new      ! asymmetric factor        : two-stream truncation
    real(RP) :: b_new0     ! planck function          : two-stream truncation
    real(RP) :: b_new1
    real(RP) :: b_new2
    real(RP) :: c0
    real(RP) :: c1
    real(RP) :: c2
    real(RP) :: Pmns, Ppls ! Phase  function          : two-stream truncation
    real(RP) :: Smns, Spls ! Source function          : two-stream truncation

    ! working
    real(RP) :: X, Y                 ! X-, X+
    real(RP) :: lamda                ! eigenvalue of X-, X+
    real(RP) :: E
    real(RP) :: Apls_mns, Bpls_mns   ! A+/A-, B+/B-
    real(DP) :: V0mns, V0pls         ! V0-, V0+
    real(DP) :: V1mns, V1pls         ! V1-, V1+
    real(DP) :: Dmns0, Dmns1, Dmns2  ! D0-, D1-, D2-
    real(DP) :: Dpls0, Dpls1, Dpls2  ! D0+, D1+, D2+
    real(RP) :: SIGmns, SIGpls       ! sigma-, sigma+
    real(RP) :: Qgamma               ! Q * gamma
    real(RP) :: zerosw, tmp

    real(RP) :: Tdir                                      ! transmission factor for solar direct

    real(RP) :: Umns, Upls ! flux intensity

    real(RP) :: factor
    real(RP) :: Wmns_irgn, M_irgn, W_irgn, Wpls_irgn, Wscale_irgn

    integer, parameter :: I_SFC2TOA = 1
    integer, parameter :: I_TOA2SFC = 2

    real(RP) :: sw
    integer  :: k, kk, icloud, ich, ic
    !---------------------------------------------------------------------------

    M_irgn      = M(irgn)
    W_irgn      = W(irgn)
    Wmns_irgn   = Wmns(irgn)
    Wpls_irgn   = Wpls(irgn)
    Wscale_irgn = Wscale(irgn)

    !$acc loop collapse(3)
    do ich = 1, chmax
    do icloud = 1, MSTRN_ncloud
!OCL LOOP_FISSION_TARGET(LS)
    do k = 1, rd_kmax

       !---< two-stream truncation >---
       tau_new = ( 1.0_RP - omg(k,icloud,ich)*g(k,2,icloud) ) * tau(k,icloud,ich)

       omg_new = ( 1.0_RP - g(k,2,icloud) ) / ( 1.0_RP - omg(k,icloud,ich)*g(k,2,icloud) ) * omg(k,icloud,ich)
       omg_new = min( omg_new, EPS1 )

       g_new   = ( g(k,1,icloud) - g(k,2,icloud) ) / ( 1.0_RP - g(k,2,icloud) )

#if defined(NVIDIA) || defined(SX)
       Tdir0(k,icloud,ich) = exp( -min( tau_new/cosSZA, 1.E+3_RP ) ) ! apply exp limiter
#else
       Tdir0(k,icloud,ich) = exp(-tau_new/cosSZA)
#endif

       factor   = ( 1.0_RP - omg(k,icloud,ich)*g(k,2,icloud) )
       b_new0 = b(k,0,icloud,ich)
       b_new1 = b(k,1,icloud,ich) / factor
       b_new2 = b(k,2,icloud,ich) / (factor*factor)
       c0     = Wmns_irgn * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new0
       c1     = Wmns_irgn * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new1
       c2     = Wmns_irgn * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new2

       !--- P+, P-
       Pmns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M_irgn*M_irgn )
       Ppls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M_irgn*M_irgn )

       !--- S+, S-
       Smns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M_irgn*cosSZA )
       Spls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M_irgn*cosSZA )

       !---< calculate R, T, e+, e- >---
       sw = 0.5_RP + sign(0.5_RP,tau_new-RD_EPS)

       !--- X, Y
       X     =  ( 1.0_RP - W_irgn * ( Ppls - Pmns ) ) / M_irgn
       Y     =  ( 1.0_RP - W_irgn * ( Ppls + Pmns ) ) / M_irgn
       !X     =  max( ( 1.0_RP - W_irgn * ( Ppls - Pmns ) ) / M_irgn, 1.E-30 )
       !Y     =  max( ( 1.0_RP - W_irgn * ( Ppls + Pmns ) ) / M_irgn, 1.E-30 )
       lamda = sqrt(X*Y)
#if defined(NVIDIA) || defined(SX)
       E     = exp( -min( lamda*tau_new, 1.E+3_RP ) ) ! apply exp limiter
#else
       E     = exp(-lamda*tau_new)
#endif

       !--- A+/A-, B+/B-
       Apls_mns = ( X * ( 1.0_DP+E ) - lamda * ( 1.0_DP-E ) ) &
                / ( X * ( 1.0_DP+E ) + lamda * ( 1.0_DP-E ) )
       Bpls_mns = ( X * ( 1.0_DP-E ) - lamda * ( 1.0_DP+E ) ) &
                / ( X * ( 1.0_DP-E ) + lamda * ( 1.0_DP+E ) )

       !--- R, T
       R0(k,icloud,ich) = (        sw ) * 0.5_RP * ( Apls_mns + Bpls_mns ) &
                        + ( 1.0_RP-sw ) * (          tau_new * (          Pmns ) / M_irgn )
       T0(k,icloud,ich) = (        sw ) * 0.5_RP * ( Apls_mns - Bpls_mns ) &
                        + ( 1.0_RP-sw ) * ( 1.0_RP - tau_new * ( 1.0_RP - Ppls ) / M_irgn )

       !--- thermal source
       Dmns0 = c0 / Y + 2.0_RP * c2 / (X*Y*Y) + c1 / (X*Y)
       Dpls0 = c0 / Y + 2.0_RP * c2 / (X*Y*Y) - c1 / (X*Y)
       Dmns1 = c1 / Y + 2.0_RP * c2 / (X*Y)
       Dpls1 = c1 / Y - 2.0_RP * c2 / (X*Y)
       Dmns2 = c2 / Y
       Dpls2 = c2 / Y

       V0mns = Dmns0
       V0pls = Dpls0
       V1mns = Dmns0 + Dmns1*tau_new + Dmns2*tau_new*tau_new
       V1pls = Dpls0 + Dpls1*tau_new + Dpls2*tau_new*tau_new

       Em_LW(k,icloud,ich) = (        sw ) * ( V0mns - R0(k,icloud,ich) * V0pls - T0(k,icloud,ich) * V1mns ) &
                           + ( 1.0_RP-sw ) * 0.5_RP * tau_new * ( 2.0_RP*c0 + c1*tau_new + c2*tau_new*tau_new )
       Ep_LW(k,icloud,ich) = (        sw ) * ( V1pls - T0(k,icloud,ich) * V0pls - R0(k,icloud,ich) * V1mns ) &
                           + ( 1.0_RP-sw ) * 0.5_RP * tau_new * ( 2.0_RP*c0 + c1*tau_new + c2*tau_new*tau_new )

       !--- solar source
       SIGmns = Wmns_irgn * ( Spls - Smns )
       SIGpls = Wmns_irgn * ( Spls + Smns )

       tmp    = X*Y*cosSZA - 1.0/cosSZA
       zerosw = 1.0_RP - sign(1.0_RP,abs(tmp)-EPS) ! if abs(tmp)<EPS then 2, otherwise 0
       Qgamma = ( SIGpls*X*cosSZA + SIGmns ) / ( tmp + zerosw*EPS )

       V0pls = 0.5_RP * ( ( 1.0_RP + 1.0_RP/(X*cosSZA) ) * Qgamma + SIGmns / X )
       V0mns = 0.5_RP * ( ( 1.0_RP - 1.0_RP/(X*cosSZA) ) * Qgamma - SIGmns / X )

       V1pls = V0pls * Tdir0(k,icloud,ich)
       V1mns = V0mns * Tdir0(k,icloud,ich)

       Em_SW(k,icloud,ich) = (        sw ) * ( V0mns - R0(k,icloud,ich) * V0pls - T0(k,icloud,ich) * V1mns ) &
                           + ( 1.0_RP-sw ) * Wmns_irgn * Smns * tau_new * sqrt( Tdir0(k,icloud,ich) )
       Ep_SW(k,icloud,ich) = (        sw ) * ( V1pls - T0(k,icloud,ich) * V0pls - R0(k,icloud,ich) * V1mns ) &
                           + ( 1.0_RP-sw ) * Wmns_irgn * Spls * tau_new * sqrt( Tdir0(k,icloud,ich) )
    enddo
    enddo
    enddo


    !---< consider partial cloud layer: semi-random over-wrapping >---

    do icloud = 1, MSTRN_ncloud
       if ( icloud == I_ClearSky ) then
          !$acc loop collapse(2)
          do ich = 1, chmax
          do k = 1, rd_kmax
             ic = icloud + ( ich - 1 ) * MSTRN_ncloud
             cf(k,ic) = 0.0_RP
          enddo
          enddo
       else
          !$acc loop independent collapse(2)
          do ich = 1, chmax
          do k = 1, rd_kmax
             ic = icloud + ( ich - 1 ) * MSTRN_ncloud
             cf(k,ic) = cldfrac(k)
          enddo
          enddo
       endif
    enddo


    !$acc loop collapse(2)
    do ich = 1, chmax
    do icloud = 1, MSTRN_ncloud
       ic = icloud + ( ich - 1 ) * MSTRN_ncloud
       tau_bar_sol(1,ic) = fsol
       !$acc loop seq
       do k = 2, rd_kmax+1 ! k-recurrence
          Tdir = (        cf(k-1,ic) ) * Tdir0(k-1,I_Cloud   ,ich) &
               + ( 1.0_RP-cf(k-1,ic) ) * Tdir0(k-1,I_ClearSky,ich)
          tau_bar_sol(k,ic) = tau_bar_sol(k-1,ic) * Tdir
       enddo
    enddo
    enddo

    !$acc loop collapse(3)
    do ich = 1, chmax
    do icloud = 1, MSTRN_ncloud
    do k = 1, rd_kmax+1
       ic = icloud + ( ich - 1 ) * MSTRN_ncloud
       flux_direct(k,icloud,ich) = cosSZA * tau_bar_sol(k,ic)
    enddo
    enddo
    enddo

    !$acc loop collapse(3)
    do ich = 1, chmax
    do icloud = 1, MSTRN_ncloud
    do k = 1, rd_kmax
       ic = icloud + ( ich - 1 ) * MSTRN_ncloud
       Em(k,ic) = (        cf(k,ic) ) * ( Em_LW(k,I_Cloud,   ich) &
                                        + Em_SW(k,I_Cloud,   ich) * tau_bar_sol(k,ic) ) &
                + ( 1.0_RP-cf(k,ic) ) * ( Em_LW(k,I_ClearSky,ich) &
                                        + Em_SW(k,I_ClearSky,ich) * tau_bar_sol(k,ic) )

       Ep(k,ic) = (        cf(k,ic) ) * ( Ep_LW(k,I_Cloud,   ich) &
                                        + Ep_SW(k,I_Cloud,   ich) * tau_bar_sol(k,ic) ) &
                + ( 1.0_RP-cf(k,ic) ) * ( Ep_LW(k,I_ClearSky,ich) &
                                        + Ep_SW(k,I_ClearSky,ich) * tau_bar_sol(k,ic) )

       R(k,ic) = (        cf(k,ic) ) * R0(k,I_Cloud   ,ich) &
               + ( 1.0_RP-cf(k,ic) ) * R0(k,I_ClearSky,ich)
       T(k,ic) = (        cf(k,ic) ) * T0(k,I_Cloud   ,ich) &
               + ( 1.0_RP-cf(k,ic) ) * T0(k,I_ClearSky,ich)
    enddo
    enddo
    enddo

    !---< Adding-Doubling method >---
    ! [note] TOA->Surface is positive direction. "pls" means upper to lower altitude.
    !$acc loop collapse(2)
    do ich = 1, chmax
    do icloud = 1, MSTRN_ncloud
       ! at lambert surface
       ic = icloud + ( ich - 1 ) * MSTRN_ncloud
       R(rd_kmax+1,ic) = (        cf(rd_kmax,ic) ) * albedo_sfc_diffuse &
                       + ( 1.0_RP-cf(rd_kmax,ic) ) * albedo_sfc_diffuse
       T(rd_kmax+1,ic) = 0.0_RP

       ! currently, Em0(I_Cloud) and Em0(I_ClearSky) is the same
       Em0(I_Cloud   ) = Wpls_irgn * ( flux_direct(rd_kmax+1,icloud,ich) * albedo_sfc_direct / (W_irgn*M_irgn) &
                       + 2.0_RP * PI * ( 1.0_RP-albedo_sfc_diffuse ) * b_sfc )
       Em0(I_ClearSky) = Wpls_irgn * ( flux_direct(rd_kmax+1,icloud,ich) * albedo_sfc_direct / (W_irgn*M_irgn) &
                       + 2.0_RP * PI * ( 1.0_RP-albedo_sfc_diffuse ) * b_sfc )

       Em(rd_kmax+1,ic) = (        cf(rd_kmax,ic) ) * Em0(I_Cloud   ) &
                        + ( 1.0_RP-cf(rd_kmax,ic) ) * Em0(I_ClearSky)
       Ep(rd_kmax+1,ic) = 0.0_RP


       R12pls(ic,rd_kmax+1) = R (rd_kmax+1,ic)
       E12mns(ic,rd_kmax+1) = Em(rd_kmax+1,ic)
       R12mns(ic,1) = R (1,ic)
       E12pls(ic,1) = Ep(1,ic)
    enddo
    enddo

    ! cloud emissivity
    ich = chmax / 2 + 1
    icloud = 2
    !$acc loop collapse(2)
    do ich = 1, chmax
    do k = 1, rd_kmax
       ic = icloud + ( ich - 1 ) * MSTRN_ncloud
       emisCLD(k,ich) = 1.0_RP - R(k,ic) - T(k,ic)
    enddo
    enddo

    !$acc loop seq
    do kk = 1, rd_kmax
!OCL ITERATIONS(max=20)
       do ic = 1, chmax * MSTRN_ncloud
          ! adding: surface to TOA
          k = rd_kmax - kk + 1
          R12pls(ic,k) = R (k,ic) + T(k,ic) / ( 1.0_RP - R12pls(ic,k+1) * R(k,ic) ) &
                       * ( R12pls(ic,k+1) * T (k,ic) )
          E12mns(ic,k) = Em(k,ic) + T(k,ic) / ( 1.0_RP - R12pls(ic,k+1) * R(k,ic) ) &
                       * ( R12pls(ic,k+1) * Ep(k,ic) + E12mns(ic,k+1) )
          ! adding: TOA to surface
          k = kk + 1
          R12mns(ic,k) = R (k,ic) + T(k,ic) / ( 1.0_RP - R12mns(ic,k-1) * R(k,ic) ) &
                       * ( R12mns(ic,k-1) * T (k,ic) )
          E12pls(ic,k) = Ep(k,ic) + T(k,ic) / ( 1.0_RP - R12mns(ic,k-1) * R(k,ic) ) &
                       * ( R12mns(ic,k-1) * Em(k,ic) + E12pls(ic,k-1) )
       enddo
    enddo

    !--- radiative flux at cell wall
    !$acc loop collapse(2)
    do ich = 1, chmax
    do icloud = 1, MSTRN_ncloud
       ic = icloud + ( ich - 1 ) * MSTRN_ncloud
       ! TOA boundary
       Upls = 0.0_RP
       Umns = E12mns(ic,1) + R12pls(ic,1) * Upls
       flux(1,I_up,icloud,ich) = Wscale_irgn * Umns
       flux(1,I_dn,icloud,ich) = Wscale_irgn * Upls + flux_direct(1,icloud,ich)
    enddo
    enddo

    !$acc loop collapse(3)
    do ich = 1, chmax
    do icloud = 1, MSTRN_ncloud
    do k = 2, rd_kmax+1
       ic = icloud + ( ich - 1 ) * MSTRN_ncloud
       Upls = ( E12pls(ic,k-1) + R12mns(ic,k-1)*E12mns(ic,k) ) / ( 1.0_RP - R12mns(ic,k-1)*R12pls(ic,k) )
       Umns = E12mns(ic,k) + R12pls(ic,k) * Upls
       flux(k,I_up,icloud,ich) = Wscale_irgn * Umns
       flux(k,I_dn,icloud,ich) = Wscale_irgn * Upls + flux_direct(k,icloud,ich)
    enddo
    enddo
    enddo

    return
  end subroutine RD_MSTRN_two_stream

end module scale_atmos_phy_rd_mstrnx
