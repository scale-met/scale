!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          mstrnX
!!          Ref: Nakajima and Tanaka(1986)
!!               Nakajima et al.(2000)
!!               Sekiguchi and Nakajima(2008)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-02-06 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_rd_mstrnx
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

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
  public :: ATMOS_PHY_RD_mstrnx

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
  private :: RD_albedo_ocean

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

  integer,  private :: I_MPAE2RD      (RD_naero)   ! look-up table between input aerosol category and MSTRN particle type
  data I_MPAE2RD(1      :N_HYD     ) / 1, 1, 2, 2, 2, 2 /
  data I_MPAE2RD(N_HYD+1:N_HYD+N_AE) / 1, 2, 3, 4, 5, 6, 7, 8, 9 /

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
  integer,  private, save :: MSTRN_nptype   = 11 !< # of particle species
                                                      !  1: water cloud
                                                      !  2: ice cloud
                                                      !  3: dust-like
                                                      !  4: soot
                                                      !  5: volcanic-ash
                                                      !  6: H2SO4
                                                      !  7: rural
                                                      !  8: sea salt
                                                      !  9: urban
                                                      ! 10: tropo.
                                                      ! 11: yellow dust
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
  integer,  private, save      :: MSTRN_nradius  =  6 !< # of radius mode for hygroscopic parameter
  integer,  private, parameter :: MSTRN_ncloud   =  2 !< # of cloud types [ClearSky/Cloud]

  logical,  private, save :: ATMOS_PHY_RD_MSTRN_ONLY_QCI = .false.


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
  real(RP), private, allocatable :: acfc    (:,:)       ! absorption coefficient table for CFC

  real(RP), private, allocatable :: fitPLK  (:,:)       ! fitting point for planck function
  real(RP), private, allocatable :: fsol    (:)         ! solar insolation    in each band
  real(RP), private              :: fsol_tot            ! total solar insolation
  real(RP), private, allocatable :: sfc     (:,:)       ! surface condition   in each band
  real(RP), private, allocatable :: rayleigh(:)         ! rayleigh scattering in each band
  real(RP), private, allocatable :: qmol    (:,:)       ! moments for rayleigh scattering phase function
  real(RP), private, allocatable :: q       (:,:,:,:)   ! moments for aerosol  scattering phase function

  integer,  private, allocatable :: hygro_flag(:)       ! flag for hygroscopic enlargement
  real(RP), private, allocatable :: radmode   (:,:)     ! radius mode for hygroscopic parameter



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

  ! for ocean albedo
  real(RP), private :: c_ocean_albedo(5,3)
  data c_ocean_albedo / -2.8108_RP   , -1.3651_RP,  2.9210E1_RP, -4.3907E1_RP,  1.8125E1_RP, &
                         6.5626E-1_RP, -8.7253_RP, -2.7749E1_RP,  4.9486E1_RP, -1.8345E1_RP, &
                        -6.5423E-1_RP,  9.9967_RP,  2.7769_RP  , -1.7620E1_RP,  7.0838_RP    /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_mstrnx_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_time, only: &
       TIME_NOWDATE
    use scale_grid_real, only: &
       REAL_BASEPOINT_LAT
    use scale_atmos_phy_rd_profile, only: &
       RD_PROFILE_setup       => ATMOS_PHY_RD_PROFILE_setup,       &
       RD_PROFILE_setup_zgrid => ATMOS_PHY_RD_PROFILE_setup_zgrid, &
       RD_PROFILE_read        => ATMOS_PHY_RD_PROFILE_read
    use scale_atmos_hydrometeor, only: &
       N_HYD
    use scale_atmos_aerosol, only: &
       N_AE
    implicit none

    character(len=*), intent(in) :: RD_TYPE

    real(RP)              :: ATMOS_PHY_RD_MSTRN_TOA
    integer               :: ATMOS_PHY_RD_MSTRN_KADD
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME
    integer               :: ATMOS_PHY_RD_MSTRN_nband
    integer               :: ATMOS_PHY_RD_MSTRN_nptype
    integer               :: ATMOS_PHY_RD_MSTRN_nradius

    namelist / PARAM_ATMOS_PHY_RD_MSTRN / &
       ATMOS_PHY_RD_MSTRN_TOA,                   &
       ATMOS_PHY_RD_MSTRN_KADD,                  &
       ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME,   &
       ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME,  &
       ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME, &
       ATMOS_PHY_RD_MSTRN_nband, &
       ATMOS_PHY_RD_MSTRN_nptype, &
       ATMOS_PHY_RD_MSTRN_nradius, &
       ATMOS_PHY_RD_MSTRN_ONLY_QCI

    integer :: ngas, ncfc
    integer :: ihydro, iaero
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[RADIATION] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Sekiguchi and Nakajima (2008) mstrnX radiation process'

    if ( RD_TYPE /= 'MSTRNX' ) then
       write(*,*) 'xxx RD_TYPE is not MSTRNX. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    ATMOS_PHY_RD_MSTRN_TOA                   = RD_TOA
    ATMOS_PHY_RD_MSTRN_KADD                  = RD_KADD
    ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME   = MSTRN_GASPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME  = MSTRN_AEROPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = MSTRN_HYGROPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_nband                 = MSTRN_nband
    ATMOS_PHY_RD_MSTRN_nptype                = MSTRN_nptype
    ATMOS_PHY_RD_MSTRN_nradius               = MSTRN_nradius

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_MSTRN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD_MSTRN. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_RD_MSTRN)

    RD_TOA                    = ATMOS_PHY_RD_MSTRN_TOA
    RD_KADD                   = ATMOS_PHY_RD_MSTRN_KADD
    MSTRN_GASPARA_INPUTFILE   = ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME
    MSTRN_AEROPARA_INPUTFILE  = ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME
    MSTRN_HYGROPARA_INPUTFILE = ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME
    MSTRN_nband               = ATMOS_PHY_RD_MSTRN_nband
    MSTRN_nptype              = ATMOS_PHY_RD_MSTRN_nptype
    MSTRN_nradius             = ATMOS_PHY_RD_MSTRN_nradius

    !--- setup MSTRN parameter
    call RD_MSTRN_setup( ngas, & ! [OUT]
                         ncfc  ) ! [OUT]

    !--- setup climatological profile
    call RD_PROFILE_setup

    RD_KMAX      = KMAX + RD_KADD

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
    call RD_PROFILE_setup_zgrid( RD_TOA, RD_KMAX, RD_KADD, & ! [IN]
                                 RD_zh(:), RD_z(:)         ) ! [INOUT]

    !--- read climatological profile
    call RD_PROFILE_read( RD_KMAX,                & ! [IN]
                          ngas,                   & ! [IN]
                          ncfc,                   & ! [IN]
                          RD_naero,               & ! [IN]
                          REAL_BASEPOINT_LAT,     & ! [IN]
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

    return
  end subroutine ATMOS_PHY_RD_mstrnx_setup

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_mstrnx( &
       DENS, RHOT, QTRC,      &
       CZ, FZ,                &
       fact_ocean,            &
       fact_land,             &
       fact_urban,            &
       temp_sfc, albedo_land, &
       solins, cosSZA,        &
       flux_rad,              &
       flux_rad_top,          &
       flux_rad_sfc_dn        )
!       Jval                            )
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       EPS  => CONST_EPS, &
       Mdry => CONST_Mdry, &
       Mvap => CONST_Mvap, &
       PPM  => CONST_PPM
    use scale_time, only: &
       TIME_NOWDATE
    use scale_grid_real, only: &
       REAL_BASEPOINT_LAT
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    use scale_atmos_phy_mp, only: &
       MP_EffectiveRadius => ATMOS_PHY_MP_EffectiveRadius, &
       MP_CloudFraction   => ATMOS_PHY_MP_CloudFraction,   &
       MP_Mixingratio     => ATMOS_PHY_MP_Mixingratio,     &
       MP_DENS            => ATMOS_PHY_MP_DENS
    use scale_atmos_phy_ae, only: &
       AE_EffectiveRadius => ATMOS_PHY_AE_EffectiveRadius, &
       AE_DENS            => ATMOS_PHY_AE_DENS, &
       QA_AE, &
       QS_AE
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_QV, &
       I_HC, &
       I_HI
    use scale_atmos_aerosol, only: &
       N_AE
    use scale_atmos_phy_rd_profile, only: &
       RD_PROFILE_read            => ATMOS_PHY_RD_PROFILE_read, &
       RD_PROFILE_use_climatology => ATMOS_PHY_RD_PROFILE_use_climatology
    implicit none

    real(RP), intent(in)  :: DENS           (KA,IA,JA)
    real(RP), intent(in)  :: RHOT           (KA,IA,JA)
    real(RP), intent(in)  :: QTRC           (KA,IA,JA,QA)
    real(RP), intent(in)  :: CZ             (  KA,IA,JA)    ! UNUSED
    real(RP), intent(in)  :: FZ             (0:KA,IA,JA)
    real(RP), intent(in)  :: fact_ocean     (IA,JA)
    real(RP), intent(in)  :: fact_land      (IA,JA)
    real(RP), intent(in)  :: fact_urban     (IA,JA)
    real(RP), intent(in)  :: temp_sfc       (IA,JA)
    real(RP), intent(in)  :: albedo_land    (IA,JA,2)
    real(RP), intent(in)  :: solins         (IA,JA)
    real(RP), intent(in)  :: cosSZA         (IA,JA)
    real(RP), intent(out) :: flux_rad       (KA,IA,JA,2,2,2)
    real(RP), intent(out) :: flux_rad_top   (IA,JA,2,2,2)
    real(RP), intent(out) :: flux_rad_sfc_dn(IA,JA,2,2)
!    real(RP), intent(out) :: Jval        (KA,IA,JA,CH_QA_photo)

    real(RP) :: temp   (KA,IA,JA)
    real(RP) :: pres   (KA,IA,JA)
    real(RP) :: qsat   (KA,IA,JA)
    real(RP) :: rh     (KA,IA,JA)
    real(RP) :: cldfrac(KA,IA,JA)
    real(RP) :: MP_Re  (KA,IA,JA,N_HYD)
    real(RP) :: MP_Qe  (KA,IA,JA,N_HYD)
    real(RP) :: AE_Re  (KA,IA,JA,N_AE)

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

    real(RP) :: zerosw

    integer :: ihydro, iaero, iq
    integer :: RD_k, k, i, j, v, ic
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Radiation(mstrnX)'

    call PROF_rapstart('RD_Profile', 3)

    call THERMODYN_temp_pres( temp(:,:,:),   & ! [OUT]
                              pres(:,:,:),   & ! [OUT]
                              DENS(:,:,:),   & ! [IN]
                              RHOT(:,:,:),   & ! [IN]
                              QTRC(:,:,:,:), & ! [IN]
                              TRACER_CV(:),  & ! [IN]
                              TRACER_R(:),   & ! [IN]
                              TRACER_MASS(:) ) ! [IN]

    call SATURATION_dens2qsat_liq( qsat(:,:,:), & ! [OUT]
                                   TEMP(:,:,:), & ! [IN]
                                   DENS(:,:,:)  ) ! [IN]

!OCL XFILL
    if ( I_QV > 0 ) then
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          rh(k,i,j) = QTRC(k,i,j,I_QV) / qsat(k,i,j)
       enddo
       enddo
       enddo
    else
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          rh(k,i,j) = 0.0_RP
       enddo
       enddo
       enddo
    endif

    call MP_CloudFraction( cldfrac(:,:,:), & ! [OUT]
                           QTRC(:,:,:,:),  & ! [IN]
                           EPS             ) ! [IN]

    call MP_EffectiveRadius( MP_Re(:,:,:,:), & ! [OUT]
                             QTRC (:,:,:,:), & ! [IN]
                             DENS (:,:,:)  , & ! [IN]
                             TEMP (:,:,:)    ) ! [IN]

    call AE_EffectiveRadius( AE_Re(:,:,:,:), & ! [OUT]
                             QTRC (:,:,:,:), & ! [IN]
                             rh   (:,:,:)    ) ! [IN]

    call MP_Mixingratio( MP_Qe(:,:,:,:), & ! [OUT]
                         QTRC (:,:,:,:)  ) ! [IN]

!    call AE_Mixingratio( AE_Qe(:,:,:,:), & ! [OUT]
!                         QTRC (:,:,:,:)  ) ! [IN]

    if ( ATMOS_PHY_RD_MSTRN_ONLY_QCI ) then
       do ihydro = 1, N_HYD
          if ( ihydro /= I_HC .and. ihydro /= I_HI ) then
             MP_Qe(:,:,:,ihydro) = 0.0_RP
          end if
       end do
    end if

    ! marge basic profile and value in model domain

    if ( RD_PROFILE_use_climatology ) then
       call RD_PROFILE_read( RD_KMAX,                & ! [IN]
                             MSTRN_ngas,             & ! [IN]
                             MSTRN_ncfc,             & ! [IN]
                             RD_naero,               & ! [IN]
                             REAL_BASEPOINT_LAT,     & ! [IN]
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
    endif

!OCL XFILL
    !$omp parallel do default(none)                                           &
    !$omp shared(JS,JE,IS,IE,RD_KADD,temph_merge,RD_temph,KE,RD_KMAX,KS,temp) &
    !$omp private(i,j,k,RD_k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KADD
          temph_merge(RD_k,i,j) = RD_temph(RD_k)
       enddo

       temp(KE+1,i,j) = temp(KE,i,j)
       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis

          temph_merge(RD_k,i,j) = 0.5_RP * ( temp(k,i,j) + temp(k+1,i,j) )
       enddo
       temph_merge(RD_KMAX+1,i,j) = temp(KS,i,j)
    enddo
    enddo

!OCL XFILL
    !$omp parallel do default(none) private(i,j,RD_k,k) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,RD_KADD,rhodz_merge,RD_rhodz,pres_merge,RD_pres,temp_merge) &
    !$omp shared(RD_temp,RD_KMAX,KS,dens,FZ,pres,temp)
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

!OCL XFILL
!OCL SERIAL
    do v = 1,  MSTRN_ngas
!OCL PARALLEL
    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KMAX
          gas_merge(RD_k,i,j,v) = RD_gas(RD_k,v)
       enddo
    enddo
    enddo
    enddo

    if ( I_QV > 0 ) then
       do j = JS, JE
       do i = IS, IE
          do RD_k = RD_KADD+1, RD_KMAX
             k = KS + RD_KMAX - RD_k ! reverse axis
             zerosw = sign(0.5_RP, QTRC(k,i,j,I_QV)-EPS) + 0.5_RP
             gas_merge(RD_k,i,j,1) = QTRC(k,i,j,I_QV) / Mvap * Mdry / PPM * zerosw ! [PPM]
          enddo
       enddo
       enddo
    endif

!OCL XFILL
!OCL SERIAL
    do v = 1,  MSTRN_ncfc
!OCL PARALLEL
    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KMAX
          cfc_merge(RD_k,i,j,v) = RD_cfc(RD_k,v)
       enddo
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KADD
          cldfrac_merge(RD_k,i,j) = RD_cldfrac(RD_k)
       enddo

       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis

          cldfrac_merge(RD_k,i,j) = 0.5_RP + sign( 0.5_RP, cldfrac(k,i,j)-min_cldfrac )
       enddo
    enddo
    enddo

!OCL XFILL
!OCL SERIAL
    do v = 1,  RD_naero
!OCL PARALLEL
    do j = JS, JE
    do i = IS, IE
    do RD_k = 1, RD_KADD
       aerosol_conc_merge(RD_k,i,j,v) = RD_aerosol_conc(RD_k,v)
       aerosol_radi_merge(RD_k,i,j,v) = RD_aerosol_radi(RD_k,v)
    enddo
    enddo
    enddo
    enddo

!OCL SERIAL
    do ihydro = 1, N_HYD
!OCL PARALLEL
    do j = JS, JE
    do i = IS, IE
    do RD_k = RD_KADD+1, RD_KMAX
       k = KS + RD_KMAX - RD_k ! reverse axis

       aerosol_conc_merge(RD_k,i,j,ihydro) = max( MP_Qe(k,i,j,ihydro), 0.0_RP ) &
                                           / MP_DENS(ihydro) * RHO_std / PPM ! [PPM to standard air]
       aerosol_radi_merge(RD_k,i,j,ihydro) = MP_Re(k,i,j,ihydro)
    enddo
    enddo
    enddo
    enddo

!OCL SERIAL
    do iaero = 1, N_AE

!!$       do j = JS, JE
!!$       do i = IS, IE
!!$       do RD_k = RD_KADD+1, RD_KMAX
!!$          k = KS + RD_KMAX - RD_k ! reverse axis
!!$
!!$          aerosol_conc_merge(RD_k,i,j,N_HYD+iaero) = max( AE_Qe(k,i,j,iaero), 0.0_RP ) &
!!$                                                   / AE_DENS(iaero) * RHO_std / PPM ! [PPM to standard air]
!!$          aerosol_radi_merge(RD_k,i,j,N_HYD+iaero) = AE_Re(k,i,j,iaero)
!!$       enddo
!!$       enddo
!!$       enddo

!OCL PARALLEL
       do j = JS, JE
       do i = IS, IE
       do RD_k = RD_KADD+1, RD_KMAX
          aerosol_conc_merge(RD_k,i,j,N_HYD+iaero) = RD_aerosol_conc(RD_k,N_HYD+iaero)
          aerosol_radi_merge(RD_k,i,j,N_HYD+iaero) = RD_aerosol_radi(RD_k,N_HYD+iaero)
       enddo
       enddo
       enddo

    enddo

    call PROF_rapend  ('RD_Profile', 3)
    call PROF_rapstart('RD_MSTRN_DTRN3', 3)

    ! calc radiative transfer
    call RD_MSTRN_DTRN3( RD_KMAX,                         & ! [IN]
                         MSTRN_ngas,                      & ! [IN]
                         MSTRN_ncfc,                      & ! [IN]
                         RD_naero,                        & ! [IN]
                         RD_hydro_str,                    & ! [IN]
                         RD_hydro_end,                    & ! [IN]
                         RD_aero_str,                     & ! [IN]
                         RD_aero_end,                     & ! [IN]
                         solins            (:,:),         & ! [IN]
                         cosSZA            (:,:),         & ! [IN]
                         rhodz_merge       (:,:,:),       & ! [IN]
                         pres_merge        (:,:,:),       & ! [IN]
                         temp_merge        (:,:,:),       & ! [IN]
                         temph_merge       (:,:,:),       & ! [IN]
                         temp_sfc          (:,:),         & ! [IN]
                         gas_merge         (:,:,:,:),     & ! [IN]
                         cfc_merge         (:,:,:,:),     & ! [IN]
                         aerosol_conc_merge(:,:,:,:),     & ! [IN]
                         aerosol_radi_merge(:,:,:,:),     & ! [IN]
                         I_MPAE2RD         (:),           & ! [IN]
                         cldfrac_merge     (:,:,:),       & ! [IN]
                         albedo_land       (:,:,:),       & ! [IN]
                         fact_ocean        (:,:),         & ! [IN]
                         fact_land         (:,:),         & ! [IN]
                         fact_urban        (:,:),         & ! [IN]
                         flux_rad_merge    (:,:,:,:,:,:), & ! [OUT]
                         flux_rad_sfc_dn   (:,:,:,:)      ) ! [OUT]

    call PROF_rapend  ('RD_MSTRN_DTRN3', 3)

    ! return to grid coordinate of model domain
!OCL SERIAL
    do ic = 1, 2
!OCL PARALLEL
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
    enddo

!OCL XFILL
!OCL SERIAL
    do ic = 1, 2
!OCL PARALLEL
    do j  = JS, JE
    do i  = IS, IE
       flux_rad_top(i,j,I_LW,I_up,ic) = flux_rad_merge(1,i,j,I_LW,I_up,ic)
       flux_rad_top(i,j,I_LW,I_dn,ic) = flux_rad_merge(1,i,j,I_LW,I_dn,ic)
       flux_rad_top(i,j,I_SW,I_up,ic) = flux_rad_merge(1,i,j,I_SW,I_up,ic)
       flux_rad_top(i,j,I_SW,I_dn,ic) = flux_rad_merge(1,i,j,I_SW,I_dn,ic)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_RD_mstrnx

  !-----------------------------------------------------------------------------
  !> Setup MSTRN parameter table
  subroutine RD_MSTRN_setup( &
       ngas, &
       ncfc  )
    use scale_const, only: &
       GRAV  => CONST_GRAV, &
       Rdry  => CONST_Rdry, &
       Pstd  => CONST_Pstd, &
       TEM00 => CONST_TEM00
    implicit none

    integer, intent(out) :: ngas
    integer, intent(out) :: ncfc

    integer :: nband, nstream, nfitP, nfitT, nflag !< gas             parameters for check
    integer :: nsfc, nptype, nplkord, nfitPLK      !< aerosol/surface parameters for check
    integer :: nradius                             !< hygroscopic     parameters for check

    character(len=H_LONG) :: dummy

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

    allocate( akd (MSTRN_ch_limit,MSTRN_nfitP,MSTRN_nfitT,MSTRN_ngas,MSTRN_nband) )
    allocate( skd (MSTRN_ch_limit,MSTRN_nfitP,MSTRN_nfitT,           MSTRN_nband) )
    allocate( acfc(MSTRN_ncfc,MSTRN_nband) )

    fid = IO_get_available_fid()
    open( fid,                                    &
          file   = trim(MSTRN_GASPARA_INPUTFILE), &
          form   = 'formatted',                   &
          status = 'old',                         &
          iostat = ierr                           )

       if ( ierr /= 0 ) then
          write(*,*) 'xxx Input data file does not found! ', trim(MSTRN_GASPARA_INPUTFILE)
          stop
       endif

       ! read gas parameters for check
       read(fid,*) nband, nstream, nfitP, nfitT, nflag, ncfc

       if ( nband /= MSTRN_nband ) then
          write(*,*) 'xxx Inconsistent parameter value! nband(given,file)=', MSTRN_nband, nband
          stop
       endif
       if ( nstream /= MSTRN_nstream ) then
          write(*,*) 'xxx Inconsistent parameter value! nstream(given,file)=', MSTRN_nstream, nstream
          stop
       endif
       if ( nfitP /= MSTRN_nfitP ) then
          write(*,*) 'xxx Inconsistent parameter value! nfitP(given,file)=', MSTRN_nfitP, nfitP
          stop
       endif
       if ( nfitT /= MSTRN_nfitT ) then
          write(*,*) 'xxx Inconsistent parameter value! nfitT(given,file)=', MSTRN_nfitT, nfitT
          stop
       endif
       if ( nflag /= MSTRN_nflag ) then
          write(*,*) 'xxx Inconsistent parameter value! nflag(given,file)=', MSTRN_nflag, nflag
          stop
       endif
       if ( ncfc /= MSTRN_ncfc ) then
          write(*,*) 'xxx Inconsistent parameter value! ncfc(given,file)=', MSTRN_ncfc, ncfc
          stop
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
                   read(fid,*) (akd(ich,ip,it,igasabs(igas,iw),iw),ich=1,nch(iw))
                enddo
                enddo
             enddo
          endif

          ! H2O continuum
          if ( iflgb(I_H2O_continuum,iw) > 0 ) then
             read(fid,*) dummy
             do it = 1, MSTRN_nfitT
             do ip = 1, MSTRN_nfitP
                read(fid,*) (skd(ich,ip,it,iw),ich=1,nch(iw))
             enddo
             enddo
          endif

          ! CFC absorption
          if ( iflgb(I_CFC_continuum,iw) > 0 ) then
             read(fid,*) dummy
             read(fid,*) (acfc(icfc,iw),icfc=1,MSTRN_ncfc)
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
    q(-1:0           ,:,:,:) = 0.D0 ! dummy for NaN
    q(MSTRN_nradius+1,:,:,:) = 0.D0 ! dummy for extrapolation

    open( fid,                                     &
          file   = trim(MSTRN_AEROPARA_INPUTFILE), &
          form   = 'formatted',                    &
          status = 'old',                          &
          iostat = ierr                            )

       if ( ierr /= 0 ) then
          write(*,*) 'xxx Input data file does not found! ', trim(MSTRN_AEROPARA_INPUTFILE)
          stop
       endif

       ! read aerosol/surface parameters for check
       read(fid,*) nband, nsfc, nptype, nstream, nplkord, nfitPLK

       if ( nband /= MSTRN_nband ) then
          write(*,*) 'xxx Inconsistent parameter value! nband(given,file)=', MSTRN_nband, nband
          stop
       endif
       if ( nsfc /= MSTRN_nsfc ) then
          write(*,*) 'xxx Inconsistent parameter value! nsfc(given,file)=', MSTRN_nsfc, nsfc
          stop
       endif
       if ( nptype /= MSTRN_nptype ) then
          write(*,*) 'xxx Inconsistent parameter value! nptype(given,file)=', MSTRN_nptype, nptype
          stop
       endif
       if ( nstream /= MSTRN_nstream ) then
          write(*,*) 'xxx Inconsistent parameter value! nstream(given,file)=', MSTRN_nstream, nstream
          stop
       endif
       if ( nplkord /= MSTRN_nplkord ) then
          write(*,*) 'xxx Inconsistent parameter value! nplkord(given,file)=', MSTRN_nplkord, nplkord
          stop
       endif
       if ( nfitPLK /= MSTRN_nfitPLK ) then
          write(*,*) 'xxx Inconsistent parameter value! nfitPLK(given,file)=', MSTRN_nfitPLK, nfitPLK
          stop
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
             do iptype = 1, nptype
                read(fid,*) q(1:MSTRN_nradius,iptype,im,iw)
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

    open( fid,                                      &
          file   = trim(MSTRN_HYGROPARA_INPUTFILE), &
          form   = 'formatted',                     &
          status = 'old',                           &
          iostat = ierr                             )

       if ( ierr /= 0 ) then
          write(*,*) 'xxx Input data file does not found! ', trim(MSTRN_HYGROPARA_INPUTFILE)
          stop
       endif

       read(fid,*) nptype

       if ( nptype /= MSTRN_nptype ) then
          write(*,*) 'xxx Inconsistent parameter value! nptype(given,file)=', MSTRN_nptype, nptype
          stop
       endif

       do iptype = 1, nptype
          read(fid,*) dummy
          read(fid,*) hygro_flag(iptype), nradius

          if ( nradius /= MSTRN_nradius ) then
             write(*,*) 'xxx Inconsistent parameter value! nradius(given,file)=', MSTRN_nradius, nradius
             stop
          endif

          read(fid,*) radmode(iptype,:)
       enddo

    close(fid)

    !----- report data -----
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.7)') '*** Baseline of total solar insolation : ', fsol_tot

    !---< constant parameter for main scheme >---
    RHO_std = Pstd / ( Rdry * TEM00 ) ! [kg/m3]

    M   (I_SW) = 1.0_RP / sqrt(3.0_RP)
    W   (I_SW) = 1.0_RP
    Wbar(I_SW) = 1.0_RP

    M   (I_LW) = 1.0_RP / 1.66_RP
    W   (I_LW) = 1.0_RP
    Wbar(I_LW) = 2.0_RP * M(I_LW)

    Wmns  (:) = sqrt( W(:) / M(:) )
    Wpls  (:) = sqrt( W(:) * M(:) )
    Wscale(:) = Wpls(:) / Wbar(:)

    !$acc enter data &
    !$acc& pcopyin(wgtch, fitPLK, logfitP, logfitT, fitT) &
    !$acc& pcopyin(radmode, ngasabs, igasabs) &
    !$acc& pcopyin(fsol, q, qmol, rayleigh, acfc, nch, AKD, SKD) &
    !$acc& pcopyin(Wmns, Wpls, Wscale, W, M) &
    !$acc& pcopyin(c_ocean_albedo)

    return
  end subroutine RD_MSTRN_setup

  !-----------------------------------------------------------------------------
  !> DTRN v3.2
  subroutine RD_MSTRN_DTRN3( &
       rd_kmax,      &
       ngas,         &
       ncfc,         &
       naero,        &
       hydro_str,    &
       hydro_end,    &
       aero_str,     &
       aero_end,     &
       solins,       &
       cosSZA,       &
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
       albedo_land,  &
       fact_ocean,   &
       fact_land,    &
       fact_urban,   &
       rflux,        &
       rflux_sfc_dn  )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       GRAV => CONST_GRAV, &
       Pstd => CONST_Pstd, &
       PPM  => CONST_PPM
    implicit none

    integer,  intent(in)  :: rd_kmax
    integer,  intent(in)  :: ngas
    integer,  intent(in)  :: ncfc
    integer,  intent(in)  :: naero
    integer,  intent(in)  :: hydro_str
    integer,  intent(in)  :: hydro_end
    integer,  intent(in)  :: aero_str
    integer,  intent(in)  :: aero_end
    real(RP), intent(in)  :: solins      (IA,JA)
    real(RP), intent(in)  :: cosSZA      (IA,JA)
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
    real(RP), intent(in)  :: albedo_land (IA,JA,2)
    real(RP), intent(in)  :: fact_ocean  (IA,JA)
    real(RP), intent(in)  :: fact_land   (IA,JA)
    real(RP), intent(in)  :: fact_urban  (IA,JA)
    real(RP), intent(out) :: rflux       (rd_kmax+1,IA,JA,2,2,MSTRN_ncloud)
    real(RP), intent(out) :: rflux_sfc_dn(IA,JA,2,2)                        ! surface downward radiation flux (LW/SW,direct/diffuse)

    ! for P-T fitting
    real(RP) :: dz_std (rd_kmax,IA,JA)       ! layer thickness at 0C, 1atm [cm]
    real(RP) :: logP   (rd_kmax,IA,JA)       ! log10(pres)
    real(RP) :: logT   (rd_kmax,IA,JA)       ! log10(temp)
    integer  :: indexP (rd_kmax,IA,JA)       ! index for interpolation in P-fitting
    real(RP) :: factP  (rd_kmax,IA,JA)       ! interpolation factor    in P-fitting
    real(RP) :: factT32(rd_kmax,IA,JA)       ! interpolation factor    in T-fitting
    real(RP) :: factT21(rd_kmax,IA,JA)       ! interpolation factor    in T-fitting
    integer  :: indexR (rd_kmax,IA,JA,naero) ! index for interpolation in R-fitting
    real(RP) :: factR  (rd_kmax,IA,JA,naero) ! interpolation factor    in R-fitting

    ! for optical thickness by gas
    real(RP) :: tauGAS(rd_kmax,IA,JA,MSTRN_ch_limit) ! optical thickness by gas absorption (total)
    real(RP) :: A1, A2, A3, factPT
    real(RP) :: qv, length
    integer  :: gasno

    ! for optical thickness by particles
    real(RP) :: tauPR   (rd_kmax,IA,JA,MSTRN_ncloud)               ! optical thickness        by Rayleigh/cloud/aerosol
    real(RP) :: omgPR   (rd_kmax,IA,JA,MSTRN_ncloud)               ! single scattering albedo by Rayleigh/cloud/aerosol
    real(RP) :: optparam(rd_kmax,IA,JA,MSTRN_nmoment,MSTRN_ncloud) ! optical parameters
    real(RP) :: q_fit, dp_P

    ! for albedo
    real(RP) :: albedo_sfc  (IA,JA,MSTRN_ncloud) ! surface albedo
!     real(RP) :: albedo_ocean(IA,JA,2, MSTRN_ncloud) ! surface albedo
!     real(RP) :: tau_column  (IA,JA, MSTRN_ncloud)

    ! for planck functions
    real(RP) :: bbar (rd_kmax  ,IA,JA) ! planck functions for thermal source at the interface
    real(RP) :: bbarh(rd_kmax+1,IA,JA) ! planck functions for thermal source at the center
    real(RP) :: b_sfc(IA,JA)           ! planck functions for thermal source at the surface
    real(RP) :: wl, beta

    ! for two-stream
    real(RP) :: tau(rd_kmax,IA,JA,    MSTRN_ncloud) ! total optical thickness
    real(RP) :: omg(rd_kmax,IA,JA,    MSTRN_ncloud) ! single scattering albedo
    real(RP) :: g  (rd_kmax,IA,JA,0:2,MSTRN_ncloud) ! two-stream approximation factors
                                                    ! 0: always 1
                                                    ! 1: asymmetry factor
                                                    ! 2: truncation factor
    real(RP) :: b  (rd_kmax,IA,JA,0:2,MSTRN_ncloud) ! planck expansion coefficients (zero if SW)
    real(RP) :: fsol_rgn(IA,JA)                     ! solar insolation              (zero if LW)

    real(RP) :: flux       (rd_kmax+1,IA,JA,2,MSTRN_ncloud) ! upward/downward flux
    real(RP) :: flux_direct(rd_kmax+1,IA,JA  ,MSTRN_ncloud) ! downward flux (direct solar)

    real(RP) :: zerosw
    real(RP) :: valsum
    integer  :: chmax
    integer  :: ip, ir, irgn
    integer  :: igas, icfc, iaero, iptype
    integer  :: iw, ich, iplk, icloud, im
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    !$acc data pcopy(rflux) &
    !$acc& pcopyin(solins, cosSZA, rhodz, pres, temp, temph, temp_sfc, gas, cfc) &
    !$acc& pcopyin(aerosol_conc, aerosol_radi, aero2ptype, cldfrac, albedo_land, fact_ocean, fact_land, fact_urban) &
    !$acc& create(dz_std, logP, logT, indexP, factP, factT32, factT21, indexR, factR) &
    !$acc& create(tauGAS, tauPR, omgPR, optparam, albedo_sfc, albedo_ocean, tau_column) &
    !$acc& create(bbar, bbarh, b_sfc) &
    !$acc& create(tau, omg, g, b, fsol_rgn, flux, flux_direct)

!OCL XFILL
    !$acc kernels pcopy(dz_std) pcopyin(rhodz)
    !$acc loop gang
    do j = JS, JE
    !$acc loop gang vector(8)
    do i = IS, IE
    !$acc loop gang vector(32)
    do k = 1, rd_kmax
       dz_std(k,i,j) = rhodz(k,i,j) / RHO_std * 100.0_RP ! [cm]
    enddo
    enddo
    enddo
    !$acc end kernels

!OCL XFILL
    !$acc kernels pcopy(logP, logT) pcopyin(pres, temp)
    !$acc loop gang
    do j = JS, JE
    !$acc loop gang vector(8)
    do i = IS, IE
    !$acc loop gang vector(32)
    do k = 1, rd_kmax
       logP(k,i,j) = log10( pres(k,i,j) )
       logT(k,i,j) = log10( temp(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels pcopy(indexP) pcopyin(logP, logfitP)
    !$acc loop gang
    do j = JS, JE
    !$acc loop gang vector(8)
    do i = IS, IE
    !$acc loop gang vector(32)
    do k = 1, rd_kmax
       indexP(k,i,j) = MSTRN_nfitP
       !$acc loop seq
       do ip = MSTRN_nfitP, 2, -1
          if( logP(k,i,j) >= logfitP(ip) ) indexP(k,i,j) = ip
       enddo
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels pcopy(factP, factT32, factT21) &
    !$acc& pcopyin(indexP, logP, logfitP, logT, logfitT, temp, fitt)
    !$acc loop gang
    !$omp parallel do default(none) private(i,j,k,ip) OMP_SCHEDULE_ collapse(2) &
    !$omp shared(JS,JE,IS,IE,rd_kmax,indexP,factP,factT32,logP,logfitP,logT,logfitT,temp,fitT) &
    !$omp shared(factT21)
    do j = JS, JE
    !$acc loop gang vector(8)
    do i = IS, IE
    !$acc loop gang vector(32) private(ip)
    do k = 1, rd_kmax
       ip = indexP(k,i,j)

       factP(k,i,j) = ( logP(k,i,j) - logfitP(ip-1) ) / ( logfitP(ip) - logfitP(ip-1) )

       factT32(k,i,j) = ( logT(k,i,j) - logfitT(2)  ) / ( logfitT(3) - logfitT(2) ) &
                      * ( temp(k,i,j) - fitT(1)     ) / ( fitT(3)    - fitT(1)    )
       factT21(k,i,j) = ( logT(k,i,j) - logfitT(2)  ) / ( logfitT(2) - logfitT(1) ) &
                      * ( fitT(3)     - temp(k,i,j) ) / ( fitT(3)    - fitT(1)    )
    enddo
    enddo
    enddo
    !$acc end kernels

    !---< interpolation of mode radius & hygroscopic parameter (R-fitting) >---
!OCL SERIAL
    do iaero = 1, naero
       iptype = aero2ptype(iaero)

       !$acc kernels pcopy(indexR, factR) pcopyin(aero2ptype, aerosol_radi, radmode)
       !$acc loop gang
!OCL PARALLEL
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,j,k,ir) &
       !$omp shared(JS,JE,IS,IE,rd_kmax,aerosol_radi,iaero,radmode,iptype,indexR,factR,MSTRN_nradius)
       do j = JS, JE
       !$acc loop gang vector(8)
       do i = IS, IE
       !$acc loop gang vector(32)
       do k = 1, rd_kmax
          if ( aerosol_radi(k,i,j,iaero) <= radmode(iptype,1) ) then ! extrapolation

             ir = 1
             indexR(k,i,j,iaero) = ir
             factR (k,i,j,iaero) = ( aerosol_radi(k,i,j,iaero) - radmode(iptype,ir) ) &
                                 / ( radmode(iptype,ir+1)      - radmode(iptype,ir) )

          elseif( aerosol_radi(k,i,j,iaero) > radmode(iptype,MSTRN_nradius) ) then ! extrapolation
             ! [Note] Extrapolation sometimes makes unexpected value
             ! optical thickness is set to zero. This treatment is Ad Hoc.

             ir = MSTRN_nradius
             indexR(k,i,j,iaero) = ir
             factR (k,i,j,iaero) = 1.0_RP

          else
             !$acc loop seq
             indexR(k,i,j,iaero) = -1
             do ir = 1, MSTRN_nradius-1
                if (       aerosol_radi(k,i,j,iaero) <= radmode(iptype,ir+1) &
                     .AND. aerosol_radi(k,i,j,iaero) >  radmode(iptype,ir  ) ) then ! interpolation

                   indexR(k,i,j,iaero) = ir
                   factR (k,i,j,iaero) = ( aerosol_radi(k,i,j,iaero) - radmode(iptype,ir) ) &
                                       / ( radmode(iptype,ir+1)      - radmode(iptype,ir) )

                   exit
                endif
             enddo
             ! indexR == -1 if some variables have NaN value.
!             write operation prevents optimization (auto parallelization)
!             if ( indexR(k,i,j,iaero) == -1 ) then
!                write(*,*) 'xxx invalid index', k,i,j, iaero, aerosol_radi(k,i,j,iaero)
!                call PRC_MPIstop
!             end if
          endif
       enddo
       enddo
       enddo
       !$acc end kernels
    enddo

    ! initialize
    !$acc kernels pcopy(rflux)
    rflux       (:,:,:,:,:,:) = 0.0_RP
    rflux_sfc_dn(:,:,:,:)     = 0.0_RP
    !$acc end kernels

    !$acc wait

    do iw = 1, MSTRN_nband
       irgn = iflgb(I_SWLW,iw) + 1
       chmax = nch(iw)

       !---< interpolation of gas parameters (P-T fitting) >---
!OCL XFILL
!OCL SERIAL
       do ich = 1, chmax
          !$acc kernels pcopy(tauGAS) async(0)
          !$acc loop gang
!OCL PARALLEL
          do j = JS, JE
          !$acc loop gang vector(8)
          do i = IS, IE
          !$acc loop gang vector(32)
          do k = 1, rd_kmax
             tauGAS(k,i,j,ich) = 0.0_RP
          enddo
          enddo
          enddo
          !$acc end kernels
       enddo

       !--- Gas line absorption
!OCL SERIAL
       do igas = 1, ngasabs(iw)
          gasno = igasabs(igas,iw)

          !$acc kernels pcopy(tauGAS) &
          !$acc& pcopyin(indexP, AKD, factP, factT32, factT21, gas, dz_std, ngasabs, igasabs) async(0)
          !$acc loop gang
!OCL PARALLEL
          !$omp parallel do default(none)                                                     &
          !$omp shared(JS,JE,IS,IE,rd_kmax,indexP,gas,igasabs,igas,iw,dz_std,chmax,AKD,gasno) &
          !$omp shared(factP,factT32,factT21,tauGAS)                                          &
          !$omp private(i,j,k,ip,length,A1,A2,A3,ich,factPT) OMP_SCHEDULE_
          do j = JS, JE
          !$acc loop gang vector(8)
          do i = IS, IE
          !$acc loop gang vector(32)
          do k = 1, rd_kmax
             ip = indexP(k,i,j)

             length = gas(k,i,j,igasabs(igas,iw)) * PPM * dz_std(k,i,j)

             !$acc loop seq
             do ich = 1, chmax
                A1 = AKD(ich,ip-1,1,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                   + AKD(ich,ip  ,1,gasno,iw) * (          factP(k,i,j) )
                A2 = AKD(ich,ip-1,2,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                   + AKD(ich,ip  ,2,gasno,iw) * (          factP(k,i,j) )
                A3 = AKD(ich,ip-1,3,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                   + AKD(ich,ip  ,3,gasno,iw) * (          factP(k,i,j) )

                factPT = factT32(k,i,j)*(A3-A2) + A2 + factT21(k,i,j)*(A2-A1)

                tauGAS(k,i,j,ich) = tauGAS(k,i,j,ich) + 10.0_RP**factPT * length
             enddo ! channel loop
          enddo
          enddo
          enddo
          !$acc end kernels
       enddo ! gas loop

       !--- Gas broad absorption
       if ( iflgb(I_H2O_continuum,iw) == 1 ) then
          !$acc kernels pcopy(tauGAS) &
          !$acc& pcopyin(indexP, SKD, factP, factT32, factT21, gas, dz_std) async(0)
          !$acc loop gang
          !$omp parallel do default(none) private(i,j,k,ich,ip,qv,length,A1,A2,A3,factPT) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IE,rd_kmax,indexP,gas,dz_std,chmax,SKD,factP,factT32) &
          !$omp shared(factT21,tauGAS,iw)
          do j = JS, JE
          !$acc loop gang vector(8)
          do i = IS, IE
          !$acc loop gang vector(32)
          do k = 1, rd_kmax
             ip = indexP(k,i,j)

             qv = gas(k,i,j,1) * PPM * dz_std(k,i,j)
             length = qv*qv / ( qv + dz_std(k,i,j) )

             !$acc loop seq
             do ich = 1, chmax
                A1 = SKD(ich,ip-1,1,iw) * ( 1.0_RP-factP(k,i,j) )&
                   + SKD(ich,ip  ,1,iw) * (        factP(k,i,j) )
                A2 = SKD(ich,ip-1,2,iw) * ( 1.0_RP-factP(k,i,j) )&
                   + SKD(ich,ip  ,2,iw) * (        factP(k,i,j) )
                A3 = SKD(ich,ip-1,3,iw) * ( 1.0_RP-factP(k,i,j) )&
                   + SKD(ich,ip  ,3,iw) * (        factP(k,i,j) )

                factPT = factT32(k,i,j)*(A3-A2) + A2 + factT21(k,i,j)*(A2-A1)

                tauGAS(k,i,j,ich) = tauGAS(k,i,j,ich) + 10.0_RP**factPT * length
             enddo ! channel loop
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       if ( iflgb(I_CFC_continuum,iw) == 1 ) then
          !$acc kernels pcopy(tauGAS) pcopyin(acfc, cfc, dz_std, nch) async(0)
          !$acc loop gang
          !$omp parallel do default(none) private(i,j,k,icfc,ich,valsum) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IE,rd_kmax,ncfc,acfc,iw,cfc,dz_std,chmax,tauGAS)
          do j = JS, JE
          !$acc loop gang vector(4)
          do i = IS, IE
          !$acc loop gang vector(32)
          do k = 1, rd_kmax
             valsum = 0.0_RP
             !$acc loop seq
             do icfc = 1, ncfc
                valsum = valsum + 10.0_RP**acfc(icfc,iw) * cfc(k,i,j,icfc)
             enddo
             valsum = valsum * PPM * dz_std(k,i,j)

             do ich = 1, chmax
             !$acc loop seq
                tauGAS(k,i,j,ich) = tauGAS(k,i,j,ich) + valsum
             enddo
          enddo
          enddo
          enddo
          !$acc end kernels
       endif

       !---< particle (Rayleigh/Cloud/Aerosol) >---

       ! optical thickness, phase function
       ! im=1: extinction coefficient
       ! im=2,3,4: moments of the volume scattering phase function

       !--- Rayleigh scattering
       !$acc kernels pcopy(optparam) pcopyin(rhodz, rayleigh, qmol) async(0)
       !$acc loop gang
       !$omp parallel do default(none) private(i,j,k,im,dp_P,length) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,rd_kmax,rhodz,GRAV,Pstd,rayleigh,iw,optparam) &
       !$omp shared(qmol)
       do j = JS, JE
       !$acc loop gang vector(8)
       do i = IS, IE
       !$acc loop gang vector(32)
       do k = 1, rd_kmax
          dp_P = rhodz(k,i,j) * GRAV / Pstd
          length = rayleigh(iw) * dp_P

          !$acc loop seq
          do im = 1, MSTRN_nstream*2+2
             optparam(k,i,j,im,I_Cloud   ) = qmol(im,iw) * length
             optparam(k,i,j,im,I_ClearSky) = qmol(im,iw) * length
          enddo
       enddo
       enddo
       enddo
       !$acc end kernels

       !--- Cloud scattering
       do iaero = hydro_str, hydro_end
          iptype = aero2ptype(iaero)

          !$acc kernels pcopy(optparam) pcopyin(indexR, aero2ptype, q, factR, dz_std, aerosol_conc) async(0)
          !$acc loop gang
          !$omp parallel do default(none) private(i,j,k,im,ir,length,q_fit) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IE,rd_kmax,indexR,iaero,aerosol_conc,dz_std) &
          !$omp shared(q,iptype,iw,factR,optparam)
          do j = JS, JE
          !$acc loop gang vector(8)
          do i = IS, IE
          !$acc loop gang vector(32)
          do k = 1, rd_kmax
             ir = indexR(k,i,j,iaero)

             length = aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)

             !$acc loop seq
             do im = 1, MSTRN_nstream*2+2
                q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(k,i,j,iaero) ) &
                      + q(ir+1,iptype,im,iw) * (        factR(k,i,j,iaero) )

                optparam(k,i,j,im,I_Cloud) = optparam(k,i,j,im,I_Cloud) + q_fit * length
             enddo
          enddo
          enddo
          enddo
          !$acc end kernels
       enddo

       !--- Aerosol scattering
       do iaero = aero_str, aero_end
          iptype = aero2ptype(iaero)

          !$acc kernels pcopy(optparam) pcopyin(indexR, aero2ptype, q, factR, dz_std, aerosol_conc) async(0)
          !$acc loop gang
          !$omp parallel do default(none) private(i,j,k,im,ir,length,q_fit) OMP_SCHEDULE_ collapse(2) &
          !$omp shared(JS,JE,IS,IE,rd_kmax,indexR,iaero,aerosol_conc,dz_std) &
          !$omp shared(iptype,iw,factR,optparam,q)
          do j = JS, JE
          !$acc loop gang vector(8)
          do i = IS, IE
          !$acc loop gang vector(32)
          do k = 1, rd_kmax
             ir = indexR(k,i,j,iaero)

             length = aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)

             !$acc loop seq
             do im = 1, MSTRN_nstream*2+2
                q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(k,i,j,iaero) ) &
                      + q(ir+1,iptype,im,iw) * (        factR(k,i,j,iaero) )

                optparam(k,i,j,im,I_Cloud   ) = optparam(k,i,j,im,I_Cloud   ) + q_fit * length
                optparam(k,i,j,im,I_ClearSky) = optparam(k,i,j,im,I_ClearSky) + q_fit * length
             enddo
          enddo
          enddo
          enddo
          !$acc end kernels
       enddo

!OCL SERIAL
       do icloud = 1, MSTRN_ncloud
          !$acc kernels pcopy(tauPR, omgPR, g) pcopyin(optparam) async(0)
          !$acc loop gang
!OCL PARALLEL
          !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
          !$omp private(i,j,k,zerosw) &
          !$omp shared(JS,JE,IS,IE,rd_kmax,tauPR,optparam,icloud,omgPR,g)
          do j = JS, JE
          !$acc loop gang vector(8)
          do i = IS, IE
          !$acc loop gang vector(32)
          do k = 1, rd_kmax
             tauPR(k,i,j,icloud) = optparam(k,i,j,1,icloud)
             omgPR(k,i,j,icloud) = optparam(k,i,j,1,icloud) - optparam(k,i,j,2,icloud)

             !--- g
             zerosw = 0.5_RP - sign(0.5_RP,omgPR(k,i,j,icloud)-RD_EPS)

             g(k,i,j,0,icloud) = 1.0_RP
             g(k,i,j,1,icloud) = optparam(k,i,j,3,icloud) * ( 1.0_RP-zerosw ) / ( omgPR(k,i,j,icloud)+zerosw )
             g(k,i,j,2,icloud) = optparam(k,i,j,4,icloud) * ( 1.0_RP-zerosw ) / ( omgPR(k,i,j,icloud)+zerosw )
             !g(k,i,j,1,icloud) = max( optparam(k,i,j,3,icloud) * ( 1.0_RP-zerosw ) / ( omgPR(k,i,j,icloud)+zerosw ), 0.0_RP )
             !g(k,i,j,2,icloud) = max( optparam(k,i,j,4,icloud) * ( 1.0_RP-zerosw ) / ( omgPR(k,i,j,icloud)+zerosw ), 0.0_RP )
          enddo
          enddo
          enddo
          !$acc end kernels
       enddo

       !--- Albedo
       ! [NOTE] mstrn has look-up table for albedo.
       !        Original scheme calculates albedo by using land-use index (and surface wetness).
       !        In the atmospheric model, albedo is calculated by surface model.
!OCL SERIAL
       do icloud = 1, MSTRN_ncloud
!           !$acc kernels pcopy(tau_column) pcopyin(tauPR) async(0)
!           !$acc loop gang
!!OCL PARALLEL
!           do j = JS, JE
!           !$acc loop gang private(valsum)
!           do i = IS, IE
!              valsum = 0.0_RP
!              !$acc loop gang vector(32) reduction(+:valsum)
!              do k = 1, rd_kmax
!                 valsum = valsum + tauPR(k,i,j,icloud) ! layer-total(for ocean albedo)
!              enddo
!              tau_column(i,j,icloud) = valsum
!           enddo
!           enddo
!           !$acc end kernels

!          call RD_albedo_ocean( cosSZA      (:,:),          & ! [IN]
!                                tau_column  (:,:,icloud),   & ! [IN]
!                                albedo_ocean(:,:,:,icloud) )  ! [OUT]

!          !$acc kernels pcopy(albedo_sfc) pcopyin(fact_ocean, fact_land, fact_urban, albedo_ocean, albedo_land) async(0)
          !$acc kernels pcopy(albedo_sfc) pcopyin(albedo_land) async(0)
          !$acc loop gang vector(4)
!OCL PARALLEL
          do j = JS, JE
          !$acc loop gang vector(32)
          do i = IS, IE
             albedo_sfc(i,j,icloud) = albedo_land(i,j,irgn)
!             albedo_sfc(i,j,icloud) = fact_ocean(i,j) * albedo_ocean(i,j,irgn,icloud) &
!                                    + fact_land (i,j) * albedo_land (i,j,irgn) &
!                                    + fact_urban(i,j) * albedo_land (i,j,irgn) ! tentative
          enddo
          enddo
          !$acc end kernels
       enddo

       ! sub-channel loop
       do ich = 1, chmax

          !--- total tau & omega
!OCL SERIAL
          do icloud = 1, 2
             !$acc kernels pcopy(tau, omg) pcopyin(tauGAS, tauPR, omgPR) async(0)
             !$acc loop gang
!OCL PARALLEL
             !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
             !$omp private(i,j,k,zerosw) &
             !$omp shared(JS,JE,IS,IE,rd_kmax,tau,icloud,tauGAS,ich,tauPR,omg,omgPR)
             do j = JS, JE
             !$acc loop gang vector(8)
             do i = IS, IE
             !$acc loop gang vector(32)
             do k = 1, rd_kmax
                tau(k,i,j,icloud) = tauGAS(k,i,j,ich) + tauPR(k,i,j,icloud)
                !tau(k,i,j,icloud) = max( tauGAS(k,i,j,ich) + tauPR(k,i,j,icloud), 0.0_RP )

                zerosw = 0.5_RP - sign( 0.5_RP, tau(k,i,j,icloud)-RD_EPS ) ! if tau < EPS, zerosw = 1

                omg(k,i,j,icloud) = ( 1.0_RP-zerosw ) * omgPR(k,i,j,icloud) / ( tau(k,i,j,icloud)-zerosw ) &
                                  + (        zerosw ) * 1.0_RP

                !omg(k,i,j,icloud) = min( max( omg(k,i,j,icloud), 0.0_RP ), 1.0_RP )
             enddo
             enddo
             enddo
             !$acc end kernels
          enddo

          !--- bn
          if ( irgn == I_SW ) then ! solar

!OCL XFILL
!OCL SERIAL
             do icloud = 1, 2
                !$acc kernels pcopy(b) async(0)
                !$acc loop gang
!OCL PARALLEL
                !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
                !$omp private(i,j,k) &
                !$omp shared(JS,JE,IS,IE,rd_kmax,b,icloud)
                do j = JS, JE
                !$acc loop gang vector(8)
                do i = IS, IE
                !$acc loop gang vector(32)
                do k = 1, rd_kmax
                   b(k,i,j,0,icloud) = 0.0_RP
                   b(k,i,j,1,icloud) = 0.0_RP
                   b(k,i,j,2,icloud) = 0.0_RP
                enddo
                enddo
                enddo
                !$acc end kernels
             enddo

             !$acc kernels pcopy(b_sfc, fsol_rgn) pcopyin(fsol, solins) async(0)
             !$acc loop gang vector(4)
             do j = JS, JE
             !$acc loop gang vector(32)
             do i = IS, IE
                b_sfc(i,j) = 0.0_RP
                fsol_rgn(i,j) = fsol(iw) / fsol_tot * solins(i,j)
             enddo
             enddo
             !$acc end kernels

          elseif( irgn == I_LW ) then ! IR
             !--- set planck functions
             wl = 10000.0_RP / sqrt( waveh(iw) * waveh(iw+1) )

             ! from temp at cell center
             !$acc kernels pcopy(bbar) pcopyin(temp, fitPLK) async(0)
             !$acc loop gang
             !$omp parallel do default(none) private(i,j,k,iplk,beta) OMP_SCHEDULE_ collapse(2) &
             !$omp shared(JS,JE,IS,IE,rd_kmax,wl,temp,fitPLK,iw,bbar)
             do j = JS, JE
             !$acc loop gang vector(8)
             do i = IS, IE
             !$acc loop gang vector(32)
             do k = 1, rd_kmax
                beta = 0.0_RP
                !$acc loop seq
                do iplk = MSTRN_nfitPLK, 1, -1
                   beta = beta / ( wl*temp(k,i,j) ) + fitPLK(iplk,iw)
                enddo

                bbar(k,i,j) = exp(-beta) * temp(k,i,j) / (wl*wl)
             enddo
             enddo
             enddo
             !$acc end kernels

             ! from temp at cell wall
             !$acc kernels pcopy(bbarh) pcopyin(temph, fitPLK) async(0)
             !$acc loop gang
             !$omp parallel do default(none)                            &
             !$omp shared(JS,JE,IS,IE,rd_kmax,wl,temph,fitPLK,iw,bbarh) &
             !$omp private(i,j,k,beta,iplk) OMP_SCHEDULE_ collapse(2)
             do j = JS, JE
             !$acc loop gang vector(8)
             do i = IS, IE
             !$acc loop gang vector(32)
             do k = 1, rd_kmax+1
                beta = 0.0_RP
                !$acc loop seq
                do iplk = MSTRN_nfitPLK, 1, -1
                   beta = beta / ( wl*temph(k,i,j) ) + fitPLK(iplk,iw)
                enddo

                bbarh(k,i,j) = exp(-beta) * temph(k,i,j) / (wl*wl)
             enddo
             enddo
             enddo
             !$acc end kernels

!OCL SERIAL
             do icloud = 1, 2
                !$acc kernels pcopy(b) pcopyin(tau, bbarh, bbar) async(0)
                !$acc loop gang
!OCL PARALLEL
                !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
                !$omp private(i,j,k,zerosw) &
                !$omp shared(JS,JE,IS,IE,rd_kmax,tau,icloud,b,bbarh,bbar)
                do j = JS, JE
                !$acc loop gang vector(8)
                do i = IS, IE
                !$acc loop gang vector(32)
                do k = 1, rd_kmax
                   zerosw = 0.5_RP - sign( 0.5_RP, tau(k,i,j,icloud)-RD_EPS ) ! if tau < EPS, zerosw = 1

                   b(k,i,j,0,icloud) = bbarh(k,i,j)
                   b(k,i,j,1,icloud) = ( 1.0_RP-zerosw )           &
                                     * ( -          bbarh(k+1,i,j) &
                                         + 4.0_RP * bbar (k  ,i,j) &
                                         - 3.0_RP * bbarh(k  ,i,j) &
                                       ) / ( tau(k,i,j,icloud)-zerosw )
                   b(k,i,j,2,icloud) = ( 1.0_RP-zerosw )           &
                                     * ( +          bbarh(k+1,i,j) &
                                         - 2.0_RP * bbar (k  ,i,j) &
                                         +          bbarh(k  ,i,j) &
                                       ) / ( tau(k,i,j,icloud)*tau(k,i,j,icloud)-zerosw ) * 2.0_RP
                enddo
                enddo
                enddo
                !$acc end kernels
             enddo

             ! from temp_sfc
             !$acc kernels pcopy(b_sfc) pcopyin(fitPLK, temp_sfc) async(0)
             !$acc loop gang vector(4)
             !$omp parallel do default(none) private(i,j,iplk,beta) OMP_SCHEDULE_ &
             !$omp shared(JS,JE,IS,IE,wl,temp_sfc,fitPLK,iw,b_sfc)
             do j = JS, JE
             !$acc loop gang vector(32)
             do i = IS, IE
                beta = 0.0_RP
                !$acc loop seq
                do iplk = MSTRN_nfitPLK, 1, -1
                   beta = beta / ( wl*temp_sfc(i,j) ) + fitPLK(iplk,iw)
                enddo

                b_sfc(i,j) = exp(-beta) * temp_sfc(i,j) / (wl*wl)
             enddo
             enddo
             !$acc end kernels

!OCL XFILL
             !$acc kernels pcopy(fsol_rgn) async(0)
             !$acc loop gang vector(4)
             do j = JS, JE
             !$acc loop gang vector(32)
             do i = IS, IE
                fsol_rgn(i,j) = 0.0_RP
             enddo
             enddo
             !$acc end kernels

          endif ! solar/IR switch

          !if( IO_L ) write(IO_FID_LOG,*) "tau sum", iw, ich, sum   (tau(:,IS:IE,JS:JE,1)), sum   (tau(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "tau max", iw, ich, maxval(tau(:,IS:IE,JS:JE,1)), maxval(tau(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "tau min", iw, ich, minval(tau(:,IS:IE,JS:JE,1)), minval(tau(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "omg sum", iw, ich, sum   (omg(:,IS:IE,JS:JE,1)), sum   (omg(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "omg max", iw, ich, maxval(omg(:,IS:IE,JS:JE,1)), maxval(omg(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "omg min", iw, ich, minval(omg(:,IS:IE,JS:JE,1)), minval(omg(:,IS:IE,JS:JE,2))

          ! two-stream transfer
          call PROF_rapstart('RD_MSTRN_twst', 3)
          call RD_MSTRN_two_stream( rd_kmax,                & ! [IN]
                                    iw, ich,                & ! [IN]
                                    cosSZA     (:,:),       & ! [IN]
                                    fsol_rgn   (:,:),       & ! [IN]
                                    irgn,                   & ! [IN]
                                    tau        (:,:,:,:),   & ! [IN]
                                    omg        (:,:,:,:),   & ! [IN]
                                    g          (:,:,:,:,:), & ! [IN]
                                    b          (:,:,:,:,:), & ! [IN]
                                    b_sfc      (:,:),       & ! [IN]
                                    albedo_sfc (:,:,:),     & ! [IN]
                                    cldfrac    (:,:,:),     & ! [IN]
                                    flux       (:,:,:,:,:), & ! [OUT]
                                    flux_direct(:,:,:,:)    ) ! [OUT]
          call PROF_rapend  ('RD_MSTRN_twst', 3)

!OCL SERIAL
          do icloud = 1, 2
             !$acc kernels pcopy(rflux) pcopyin(flux, wgtch) async(0)
             !$acc loop gang
!OCL PARALLEL
             !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
             !$omp shared(JS,JE,IS,IE,rd_kmax,rflux,irgn,icloud,wgtch,ich,iw,flux)
             do j = JS, JE
             !$acc loop gang vector(8)
             do i = IS, IE
             !$acc loop gang vector(32)
             do k = 1, rd_kmax+1
                rflux(k,i,j,irgn,I_up,icloud) = rflux(k,i,j,irgn,I_up,icloud) + flux(k,i,j,I_up,icloud) * wgtch(ich,iw)
                rflux(k,i,j,irgn,I_dn,icloud) = rflux(k,i,j,irgn,I_dn,icloud) + flux(k,i,j,I_dn,icloud) * wgtch(ich,iw)
             enddo
             enddo
             enddo
             !$acc end kernels
          enddo

          do j = JS, JE
          do i = IS, IE
             rflux_sfc_dn(i,j,irgn,1) = rflux_sfc_dn(i,j,irgn,1) + flux_direct(rd_kmax+1,i,j,I_Cloud) * wgtch(ich,iw)
             rflux_sfc_dn(i,j,irgn,2) = rflux_sfc_dn(i,j,irgn,2) &
                                      + ( flux(rd_kmax+1,i,j,I_dn,I_Cloud) - flux_direct(rd_kmax+1,i,j,I_Cloud) ) * wgtch(ich,iw)
          enddo
          enddo

          !if( IO_L ) write(IO_FID_LOG,*) "flux sum", iw, ich, sum   (flux(:,IS:IE,JS:JE,1)), sum   (flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux max", iw, ich, maxval(flux(:,IS:IE,JS:JE,1)), maxval(flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux min", iw, ich, minval(flux(:,IS:IE,JS:JE,1)), minval(flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux mal", iw, ich, maxloc(flux(:,IS:IE,JS:JE,1)), maxloc(flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux mil", iw, ich, minloc(flux(:,IS:IE,JS:JE,1)), minloc(flux(:,IS:IE,JS:JE,2))

       enddo ! ICH loop
    enddo ! IW loop

    !$acc wait

    !$acc end data

    return
  end subroutine RD_MSTRN_DTRN3

  !-----------------------------------------------------------------------------
  !> Two stream calculation CORE
  subroutine RD_MSTRN_two_stream( &
       rd_kmax,         &
       iw, ich,      &
       cosSZA0,      &
       fsol,         &
       irgn,         &
       tau,          &
       omg,          &
       g,            &
       b,            &
       b_sfc,        &
       albedo_sfc,   &
       cldfrac,      &
       flux,         &
       flux_direct   )
    use scale_const, only: &
       PI   => CONST_PI,   &
       HUGE => CONST_HUGE, &
       EPS  => CONST_EPS,  &
       EPS1 => CONST_EPS1
    implicit none

    integer,  intent(in)  :: rd_kmax
    integer,  intent(in)  :: iw, ich
    real(RP), intent(in)  :: cosSZA0    (IA,JA)                          ! cos(SZA) = mu0
    real(RP), intent(in)  :: fsol       (IA,JA)                          ! solar radiation intensity
    integer,  intent(in)  :: irgn                                        ! 1:LW 2:SW
    real(RP), intent(in)  :: tau        (rd_kmax,IA,JA,    MSTRN_ncloud) ! total optical thickness          (clear-sky/cloud)
    real(RP), intent(in)  :: omg        (rd_kmax,IA,JA,    MSTRN_ncloud) ! single scattering albedo         (clear-sky/cloud)
    real(RP), intent(in)  :: g          (rd_kmax,IA,JA,0:2,MSTRN_ncloud) ! two-stream approximation factors (clear-sky/cloud)
    real(RP), intent(in)  :: b          (rd_kmax,IA,JA,0:2,MSTRN_ncloud) ! planck expansion coefficients    (clear-sky/cloud)
    real(RP), intent(in)  :: b_sfc      (IA,JA)                          ! planck function at surface
    real(RP), intent(in)  :: albedo_sfc (IA,JA,MSTRN_ncloud)             ! surface albedo                   (clear-sky/cloud)
    real(RP), intent(in)  :: cldfrac    (rd_kmax,IA,JA)                  ! cloud fraction

    real(RP), intent(out) :: flux       (rd_kmax+1,IA,JA,2,MSTRN_ncloud) ! upward(sfc->TOA)/downward(TOA->sfc) flux (clear-sky/cloud)
    real(RP), intent(out) :: flux_direct(rd_kmax+1,IA,JA,  MSTRN_ncloud) ! downward(TOA->sfc) flux, solar direct    (clear-sky/cloud)

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
    real(RP) :: cosSZA(IA,JA)
    real(RP) :: X, Y                 ! X-, X+
    real(RP) :: lamda                ! eigenvalue of X-, X+
    real(RP) :: E
    real(RP) :: Apls_mns, Bpls_mns   ! A+/A-, B+/B-
    real(RP) :: V0mns, V0pls         ! V0-, V0+
    real(RP) :: V1mns, V1pls         ! V1-, V1+
    real(RP) :: Dmns0, Dmns1, Dmns2  ! D0-, D1-, D2-
    real(RP) :: Dpls0, Dpls1, Dpls2  ! D0+, D1+, D2+
    real(RP) :: SIGmns, SIGpls       ! sigma-, sigma+
    real(RP) :: Qgamma               ! Q * gamma
    real(RP) :: zerosw, tmp

    ! main factors
    real(RP) :: Tdir0(rd_kmax,IA,JA,MSTRN_ncloud) ! transmission factor for solar direct (clear-sky/cloud)
    real(RP) :: R0   (rd_kmax,IA,JA,MSTRN_ncloud) ! reflection   factor                  (clear-sky/cloud)
    real(RP) :: T0   (rd_kmax,IA,JA,MSTRN_ncloud) ! transmission factor                  (clear-sky/cloud)
    real(RP) :: Em_LW(rd_kmax,IA,JA,MSTRN_ncloud) ! thermal source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Em_SW(rd_kmax,IA,JA,MSTRN_ncloud) ! solar   source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Ep_LW(rd_kmax,IA,JA,MSTRN_ncloud) ! thermal source (TOA->sfc)            (clear-sky/cloud)
    real(RP) :: Ep_SW(rd_kmax,IA,JA,MSTRN_ncloud) ! solar   source (TOA->sfc)            (clear-sky/cloud)

    ! Averaged factors, considering cloud overwrap
    real(RP) :: cf         (rd_kmax  ,IA,JA) ! cloud fraction
    real(RP) :: tau_bar_sol(rd_kmax+1,IA,JA) ! solar insolation through accumulated optical thickness at each layer
    real(RP) :: Tdir       (rd_kmax+1,IA,JA) ! transmission factor for solar direct
    real(RP) :: R          (rd_kmax+1,IA,JA) ! reflection   factor
    real(RP) :: T          (rd_kmax+1,IA,JA) ! transmission factor
    real(RP) :: Em         (rd_kmax+1,IA,JA) ! source (sfc->TOA)
    real(RP) :: Ep         (rd_kmax+1,IA,JA) ! source (TOA->sfc)

    ! Doubling-Adding
    real(RP) :: R12mns(rd_kmax+1,IA,JA) ! reflection factor in doubling method
    real(RP) :: R12pls(rd_kmax+1,IA,JA) ! reflection factor in doubling method
    real(RP) :: E12mns(rd_kmax+1,IA,JA) ! source function   in doubling method
    real(RP) :: E12pls(rd_kmax+1,IA,JA) ! source function   in doubling method
    real(RP) :: Umns, Upls               ! flux intensity

    real(RP) :: Em0(MSTRN_ncloud)
    real(RP) :: factor
    real(RP) :: Wmns_irgn, M_irgn, W_irgn, Wpls_irgn, Wscale_irgn

    integer, parameter :: I_SFC2TOA = 1
    integer, parameter :: I_TOA2SFC = 2
    integer            :: direction

    real(RP) :: sw
    integer  :: k, i, j, icloud
    integer  :: kij
    !---------------------------------------------------------------------------

    M_irgn      = M(irgn)
    W_irgn      = W(irgn)
    Wmns_irgn   = Wmns(irgn)
    Wpls_irgn   = Wpls(irgn)
    Wscale_irgn = Wscale(irgn)

    !$acc data pcopy(flux, flux_direct) &
    !$acc& pcopyin(cosSZA0, fsol, tau, omg, g, b, b_sfc, albedo_sfc, cldfrac) &
    !$acc& create(cosSZA, Tdir0, R0, T0, Em_LW, Em_SW, Ep_LW, Ep_SW) &
    !$acc& create(tau_bar_sol, R, T, Em, Ep, R12mns, R12pls, E12mns, E12pls) &
    !$acc& create(Tdir, x_R, x_T, x_Em, x_Ep)

    !$acc kernels pcopyin(cosSZA0) pcopy(cosSZA) async(0)
    cosSZA(IS:IE,JS:JE) = max( cosSZA0(IS:IE,JS:JE), RD_cosSZA_min )
    !$acc end kernels

!OCL SERIAL
    do icloud = 1, 2
!OCL PARALLEL,NORECURRENCE(Tdir0,R0,T0,Em_LW,Ep_LW,Em_SW,Ep_SW),MFUNC
       !$acc kernels pcopy(tdir0, r0, t0, em_lw, ep_lw, em_sw, ep_sw) &
       !$acc& pcopyin(wmns, tau, g, omg, cossza, b, m, w) async(0)
       !$omp parallel do default(none) OMP_SCHEDULE_ collapse(2) &
       !$omp private(i,j,k) &
       !$omp private(tau_new,omg_new,g_new,factor,b_new0,b_new1,b_new2,c0,c1,c2,Pmns,Ppls,Smns) &
       !$omp private(Spls,sw,X,Y,lamda,E,Apls_mns,Bpls_mns,Dmns0,Dpls0,Dmns1,Dpls1) &
       !$omp private(Dmns2,Dpls2,V0mns,V0pls,V1mns,V1pls,SIGmns,SIGpls,tmp,zerosw) &
       !$omp private(Qgamma) &
       !$omp shared(JS,JE,IS,IE,rd_kmax,omg,icloud,g,tau,EPS1,Tdir0,cosSZA,b,Wmns_irgn,PI,M_irgn) &
       !$omp shared(W_irgn,Em_LW,Ep_LW,R0,T0,Ep_SW,Em_SW,EPS)
       do j = JS, JE
       do i = IS, IE
       do k = 1, rd_kmax

          !---< two-stream truncation >---
          tau_new = ( 1.0_RP - omg(k,i,j,icloud)*g(k,i,j,2,icloud) ) * tau(k,i,j,icloud)

          omg_new = ( 1.0_RP - g(k,i,j,2,icloud) ) / ( 1.0_RP - omg(k,i,j,icloud)*g(k,i,j,2,icloud) ) * omg(k,i,j,icloud)
          omg_new = min( omg_new, EPS1 )

          g_new   = ( g(k,i,j,1,icloud) - g(k,i,j,2,icloud) ) / ( 1.0_RP - g(k,i,j,2,icloud) )

#if defined(__PGI) || defined(__ES2)
          Tdir0(k,i,j,icloud) = exp( -min( tau_new/cosSZA(i,j), 1.E+3_RP ) ) ! apply exp limiter
#else
          Tdir0(k,i,j,icloud) = exp(-tau_new/cosSZA(i,j))
#endif

          factor   = ( 1.0_RP - omg(k,i,j,icloud)*g(k,i,j,2,icloud) )
          b_new0 = b(k,i,j,0,icloud)
          b_new1 = b(k,i,j,1,icloud) / factor
          b_new2 = b(k,i,j,2,icloud) / (factor*factor)
          c0     = Wmns_irgn * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new0
          c1     = Wmns_irgn * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new1
          c2     = Wmns_irgn * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new2

          !--- P+, P-
          Pmns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M_irgn*M_irgn )
          Ppls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M_irgn*M_irgn )

          !--- S+, S-
          Smns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M_irgn*cosSZA(i,j) )
          Spls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M_irgn*cosSZA(i,j) )

          !---< calculate R, T, e+, e- >---
          sw = 0.5_RP + sign(0.5_RP,tau_new-RD_EPS)
          tau_new = max( tau_new, RD_EPS )

          !--- X, Y
          X     =  ( 1.0_RP - W_irgn * ( Ppls - Pmns ) ) / M_irgn
          Y     =  ( 1.0_RP - W_irgn * ( Ppls + Pmns ) ) / M_irgn
          !X     =  max( ( 1.0_RP - W_irgn * ( Ppls - Pmns ) ) / M_irgn, 1.E-30 )
          !Y     =  max( ( 1.0_RP - W_irgn * ( Ppls + Pmns ) ) / M_irgn, 1.E-30 )
          lamda = sqrt(X*Y)
#if defined(__PGI) || defined(__ES2)
          E     = exp( -min( lamda*tau_new, 1.E+3_RP ) ) ! apply exp limiter
#else
          E     = exp(-lamda*tau_new)
#endif

          !--- A+/A-, B+/B-
          Apls_mns = ( X * ( 1.0_RP+E ) - lamda * ( 1.0_RP-E ) ) &
                   / ( X * ( 1.0_RP+E ) + lamda * ( 1.0_RP-E ) )
          Bpls_mns = ( X * ( 1.0_RP-E ) - lamda * ( 1.0_RP+E ) ) &
                   / ( X * ( 1.0_RP-E ) + lamda * ( 1.0_RP+E ) )

          !--- R, T
          R0(k,i,j,icloud) = (        sw ) * 0.5_RP * ( Apls_mns + Bpls_mns ) &
                           + ( 1.0_RP-sw ) * (          tau_new * (          Pmns ) / M_irgn )
          T0(k,i,j,icloud) = (        sw ) * 0.5_RP * ( Apls_mns - Bpls_mns ) &
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

          Em_LW(k,i,j,icloud) = (        sw ) * ( V0mns - R0(k,i,j,icloud) * V0pls - T0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * 0.5_RP * tau_new * ( 2.0_RP*c0 + c1*tau_new + c2*tau_new*tau_new )
          Ep_LW(k,i,j,icloud) = (        sw ) * ( V1pls - T0(k,i,j,icloud) * V0pls - R0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * 0.5_RP * tau_new * ( 2.0_RP*c0 + c1*tau_new + c2*tau_new*tau_new )

          !--- solar source
          SIGmns = Wmns_irgn * ( Spls - Smns )
          SIGpls = Wmns_irgn * ( Spls + Smns )

          tmp    = X*Y*cosSZA(i,j)-1.0/cosSZA(i,j)
          zerosw = 1.0_RP - sign(1.0_RP,abs(tmp)-EPS) ! if abs(tmp)<EPS then 2, otherwise 0
          Qgamma = ( SIGpls*X*cosSZA(i,j) + SIGmns ) / ( tmp + zerosw*EPS )

          V0pls = 0.5_RP * ( ( 1.0_RP + 1.0_RP/(X*cosSZA(i,j)) ) * Qgamma + SIGmns / X )
          V0mns = 0.5_RP * ( ( 1.0_RP - 1.0_RP/(X*cosSZA(i,j)) ) * Qgamma - SIGmns / X )

          V1pls = V0pls * Tdir0(k,i,j,icloud)
          V1mns = V0mns * Tdir0(k,i,j,icloud)

          Em_SW(k,i,j,icloud) = (        sw ) * ( V0mns - R0(k,i,j,icloud) * V0pls - T0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * Wmns_irgn * Smns * tau_new * sqrt( Tdir0(k,i,j,icloud) )
          Ep_SW(k,i,j,icloud) = (        sw ) * ( V1pls - T0(k,i,j,icloud) * V0pls - R0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * Wmns_irgn * Spls * tau_new * sqrt( Tdir0(k,i,j,icloud) )

       enddo
       enddo
       enddo
       !$acc end kernels
    enddo ! cloud loop

    !---< consider partial cloud layer: semi-random over-wrapping >---

    do icloud = 1, 2
       if ( icloud == 1 ) then
!OCL XFILL
          do j = JS, JE
          do i = IS, IE
          do k = 1, rd_kmax
             cf(k,i,j) = 0.0_RP
          enddo
          enddo
          enddo
       else
!OCL XFILL
          do j = JS, JE
          do i = IS, IE
          do k = 1, rd_kmax
             cf(k,i,j) = cldfrac(k,i,j)
          enddo
          enddo
          enddo
       endif

       !$acc kernels pcopy(Tdir, R, T, Em, Ep, flux_direct, tau_bar_sol) &
       !$acc& pcopyin(cldfrac, Tdir0, R0, T0, Em_LW, Em_SW, Ep_LW, Ep_SW, fsol, cosSZA, albedo_sfc, b_sfc, wpls, w, m) async(0)

       !$acc loop gang
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,rd_kmax,Tdir,cf,Tdir0)
       do j = JS, JE
       !$acc loop gang vector(8)
       do i = IS, IE
       !$acc loop gang vector(32)
       do k = 1, rd_kmax
          Tdir(k,i,j) = (        cf(k,i,j) ) * Tdir0(k,i,j,I_Cloud   ) &
                      + ( 1.0_RP-cf(k,i,j) ) * Tdir0(k,i,j,I_ClearSky)
       enddo
       enddo
       enddo

       !$acc loop gang vector(4)
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,tau_bar_sol,fsol,rd_kmax,Tdir)
       do j = JS, JE
       !$acc loop gang vector(32)
       do i = IS, IE
          tau_bar_sol(1,i,j) = fsol(i,j) ! k-recurrence
          !$acc loop seq
          do k = 2, rd_kmax+1
             tau_bar_sol(k,i,j) = tau_bar_sol(k-1,i,j) * Tdir(k-1,i,j)
          enddo
       enddo
       enddo

       !$acc loop gang
       !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,rd_kmax,Em,cf,Em_LW,Em_SW,tau_bar_sol,Ep) &
       !$omp shared(Ep_LW,Ep_SW,flux_direct,cosSZA,icloud)
       do j = JS, JE
       !$acc loop gang vector(8)
       do i = IS, IE
       !$acc loop gang vector(32)
       do k = 1, rd_kmax
          Em(k,i,j) = (        cf(k,i,j) ) * ( Em_LW(k,i,j,I_Cloud   ) &
                                             + Em_SW(k,i,j,I_Cloud   ) * tau_bar_sol(k,i,j) ) &
                    + ( 1.0_RP-cf(k,i,j) ) * ( Em_LW(k,i,j,I_ClearSky) &
                                             + Em_SW(k,i,j,I_ClearSky) * tau_bar_sol(k,i,j) )

          Ep(k,i,j) = (        cf(k,i,j) ) * ( Ep_LW(k,i,j,I_Cloud   ) &
                                             + Ep_SW(k,i,j,I_Cloud   ) * tau_bar_sol(k,i,j) ) &
                    + ( 1.0_RP-cf(k,i,j) ) * ( Ep_LW(k,i,j,I_ClearSky) &
                                             + Ep_SW(k,i,j,I_ClearSky) * tau_bar_sol(k,i,j) )

          flux_direct(k,i,j,icloud) = cosSZA(i,j) * tau_bar_sol(k,i,j)
       enddo
       enddo
       enddo

       !$acc loop gang vector(4)
       !$omp parallel do default(none) private(i,j,Em0) OMP_SCHEDULE_ &
       !$omp shared(JS,JE,IS,IE,rd_kmax,cf,albedo_sfc,T,flux_direct,cosSZA) &
       !$omp shared(tau_bar_sol,Wpls_irgn,W_irgn,M_irgn,PI,b_sfc,Em,Ep,icloud,R)
       do j = JS, JE
       !$acc loop gang vector(32) private(Em0)
       do i = IS, IE
          ! at lambert surface
          R(rd_kmax+1,i,j) = (        cf(rd_kmax,i,j) ) * albedo_sfc(i,j,I_Cloud   ) &
                           + ( 1.0_RP-cf(rd_kmax,i,j) ) * albedo_sfc(i,j,I_ClearSky)
          T(rd_kmax+1,i,j) = 0.0_RP

          flux_direct(rd_kmax+1,i,j,icloud) = cosSZA(i,j) * tau_bar_sol(rd_kmax+1,i,j)

          Em0(I_Cloud   ) = Wpls_irgn * ( flux_direct(rd_kmax+1,i,j,icloud) * albedo_sfc(i,j,I_Cloud   ) / (W_irgn*M_irgn) &
                          + 2.0_RP * PI * ( 1.0_RP-albedo_sfc(i,j,I_Cloud   ) ) * b_sfc(i,j) )
          Em0(I_ClearSky) = Wpls_irgn * ( flux_direct(rd_kmax+1,i,j,icloud) * albedo_sfc(i,j,I_ClearSky) / (W_irgn*M_irgn) &
                          + 2.0_RP * PI * ( 1.0_RP-albedo_sfc(i,j,I_ClearSky) ) * b_sfc(i,j) )

          Em(rd_kmax+1,i,j) = (        cf(rd_kmax,i,j) ) * Em0(I_Cloud   ) &
                            + ( 1.0_RP-cf(rd_kmax,i,j) ) * Em0(I_ClearSky)
          Ep(rd_kmax+1,i,j) = 0.0_RP
       enddo
       enddo

       !$acc end kernels

       !---< Adding-Doubling method >---
       ! [note] TOA->Surface is positive direction. "pls" means upper to lower altitude.

       !$acc kernels pcopy(R12mns, R12pls, E12mns, E12pls) pcopyin(R, T, Em, Ep) async(0)

       !$acc loop gang independent
       do direction = I_SFC2TOA, I_TOA2SFC

          if ( direction == I_SFC2TOA ) then ! adding: surface to TOA

!OCL LOOP_NOFUSION,XFILL
             !$acc loop gang vector(4)
             do j = JS, JE
             !$acc loop gang vector(32)
             do i = IS, IE
                R12pls(rd_kmax+1,i,j) = R (rd_kmax+1,i,j)
                E12mns(rd_kmax+1,i,j) = Em(rd_kmax+1,i,j)
             enddo
             enddo

             !$acc loop gang vector(4)
             !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
             !$omp shared(JS,JE,IS,IE,rd_kmax,R12pls,E12mns,R,R0,T,T0,cf,Ep,Em)
             do j = JS, JE
             !$acc loop gang vector(32)
             do i = IS, IE
                !$acc loop seq
                do k = rd_kmax, 1, -1
                   R(k,i,j) = (        cf(k,i,j) ) * R0(k,i,j,I_Cloud   ) &
                            + ( 1.0_RP-cf(k,i,j) ) * R0(k,i,j,I_ClearSky)
                   T(k,i,j) = (        cf(k,i,j) ) * T0(k,i,j,I_Cloud   ) &
                            + ( 1.0_RP-cf(k,i,j) ) * T0(k,i,j,I_ClearSky)

                   R12pls(k,i,j) = R (k,i,j) + T(k,i,j) / ( 1.0_RP - R12pls(k+1,i,j) * R(k,i,j)           ) &
                                                        * ( R12pls(k+1,i,j) * T (k,i,j)                   )
                   E12mns(k,i,j) = Em(k,i,j) + T(k,i,j) / ( 1.0_RP - R12pls(k+1,i,j) * R(k,i,j)           ) &
                                                        * ( R12pls(k+1,i,j) * Ep(k,i,j) + E12mns(k+1,i,j) )
                enddo
             enddo
             enddo
          else ! adding: TOA to surface

!OCL LOOP_NOFUSION,XFILL
             !$acc loop gang vector(4)
             do j = JS, JE
             !$acc loop gang vector(32)
             do i = IS, IE
                R12mns(1,i,j) = R (1,i,j)
                E12pls(1,i,j) = Ep(1,i,j)
             enddo
             enddo

             !$acc loop gang vector(4)
             !$omp parallel do default(none) private(i,j,k) OMP_SCHEDULE_ collapse(2) &
             !$omp shared(JS,JE,IS,IE,rd_kmax,E12pls,R12mns,R,T,Ep,Em)
             do j = JS, JE
             !$acc loop gang vector(32)
             do i = IS, IE
                !$acc loop seq
                do k = 2, rd_kmax+1
                   R12mns(k,i,j) = R (k,i,j) + T(k,i,j) / ( 1.0_RP - R12mns(k-1,i,j) * R(k,i,j)         ) &
                                                        * ( R12mns(k-1,i,j) *T (k,i,j)                  )
                   E12pls(k,i,j) = Ep(k,i,j) + T(k,i,j) / ( 1.0_RP - R12mns(k-1,i,j) * R(k,i,j)         ) &
                                                        * ( R12mns(k-1,i,j)*Em(k,i,j) + E12pls(k-1,i,j) )
                enddo
             enddo
             enddo

          endif

       enddo

       !$acc end kernels

       !--- radiative flux at cell wall:
       ! [note] "d" means upper to lower altitude.

       !$acc kernels pcopy(flux) pcopyin(E12mns, E12pls, R12mns, R12pls, flux_direct, wscale) async(0)

       !$acc loop gang
       !$omp parallel do default(none) private(i,j,Upls,Umns) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,rd_kmax,E12pls,R12mns,E12mns,R12pls,flux,icloud,Wscale_irgn) &
       !$omp shared(flux_direct)
       do j = JS, JE
       !$acc loop gang vector(8)
       do i = IS, IE
          ! TOA boundary
          Upls = 0.0_RP
          Umns = E12mns(1,i,j) + R12pls(1,i,j) * Upls

          flux(1,i,j,I_up,icloud) = Wscale_irgn * Umns
          flux(1,i,j,I_dn,icloud) = Wscale_irgn * Upls + flux_direct(1,i,j,icloud)
       enddo
       enddo
       !$acc end kernels

       !$acc kernels pcopy(flux) pcopyin(E12mns, E12pls, R12mns, R12pls, flux_direct, wscale) async(0)
       !$acc loop gang
       !$omp parallel do default(none) private(i,j,k,Upls,Umns) OMP_SCHEDULE_ collapse(2) &
       !$omp shared(JS,JE,IS,IE,rd_kmax,E12pls,R12mns,E12mns,R12pls,flux,icloud,Wscale_irgn) &
       !$omp shared(flux_direct)
       do j = JS, JE
       !$acc loop gang vector(8)
       do i = IS, IE
       !$acc loop gang vector(32)
       do k = 2, rd_kmax+1
          Upls = ( E12pls(k-1,i,j) + R12mns(k-1,i,j)*E12mns(k,i,j) ) / ( 1.0_RP - R12mns(k-1,i,j)*R12pls(k,i,j) )
          Umns = E12mns(k,i,j) + R12pls(k,i,j) * Upls

          flux(k,i,j,I_up,icloud) = Wscale_irgn * Umns
          flux(k,i,j,I_dn,icloud) = Wscale_irgn * Upls + flux_direct(k,i,j,icloud)
       enddo
       enddo
       enddo

       !$acc end kernels

    enddo ! cloud loop

    !$acc end data

    return
  end subroutine RD_MSTRN_two_stream

  !-----------------------------------------------------------------------------
  ! Sea surface reflectance by Payne
  subroutine RD_albedo_ocean( &
       cosSZA,       &
       tau,          &
       albedo_ocean )
    implicit none

    real(RP), intent(in)  :: cosSZA      (IA,JA)
    real(RP), intent(in)  :: tau         (IA,JA)
    real(RP), intent(out) :: albedo_ocean(IA,JA,2)

    real(RP) :: am1, tr1, s
    real(RP) :: sw

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    !$acc kernels pcopy(albedo_ocean) pcopyin(cosSZA, tau, c_ocean_albedo) async(0)
    !$acc loop gang vector(4)
    do j = JS, JE
    !$acc loop gang vector(32)
    do i = IS, IE
       am1 = max( min( cosSZA(i,j), 0.961_RP ), 0.0349_RP )

       sw = 0.5_RP + sign(0.5_RP,tau(i,j))

       tr1 = max( min( am1 / ( 4.0_RP * tau(i,j) ), 1.0_RP ), 0.05_RP )

       s = 0.0_RP
       !$acc loop seq
       do n = 1, 5
          s = s + c_ocean_albedo(n,1) * tr1**(n-1)           &
                + c_ocean_albedo(n,2) * tr1**(n-1) * am1     &
                + c_ocean_albedo(n,3) * tr1**(n-1) * am1*am1
       enddo

       albedo_ocean(i,j,I_SW) = ( 1.0_RP-sw ) * 0.05_RP &
                              + (        sw ) * exp(s)

       albedo_ocean(i,j,I_LW) = 0.05_RP
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine RD_albedo_ocean

end module scale_atmos_phy_rd_mstrnx
