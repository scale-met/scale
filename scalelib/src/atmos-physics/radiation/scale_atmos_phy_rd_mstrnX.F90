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
module scale_atmos_phy_rd_mstrnX
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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_mstrnX_setup
  public :: ATMOS_PHY_RD_mstrnX

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
  private :: albedo_sea

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: RD_cosSZA_min = 0.017_RP ! minimum SZA (>89.0)

  integer,  private, save :: RD_KADD = 0 !< RD_KMAX = KMAX + RD_KADD

  integer,  private, save :: RD_KMAX      ! # of computational cells: z for radiation scheme
  integer,  private, save :: RD_naero     ! # of cloud/aerosol species
  integer,  private, save :: RD_hydro_str ! start index for cloud
  integer,  private, save :: RD_hydro_end ! end   index for cloud
  integer,  private, save :: RD_aero_str  ! start index for aerosol
  integer,  private, save :: RD_aero_end  ! end   index for aerosol

  real(RP), private, allocatable, save :: RD_zh          (:)   ! altitude    at the interface [km]
  real(RP), private, allocatable, save :: RD_z           (:)   ! altitude    at the center    [km]
  real(RP), private, allocatable, save :: RD_rhodz       (:)   ! density * delta z            [kg/m2]
  real(RP), private, allocatable, save :: RD_pres        (:)   ! pressure    at the center    [hPa]
  real(RP), private, allocatable, save :: RD_presh       (:)   ! pressure    at the interface [hPa]
  real(RP), private, allocatable, save :: RD_temp        (:)   ! temperature at the center    [K]
  real(RP), private, allocatable, save :: RD_temph       (:)   ! temperature at the interface [K]
  real(RP), private, allocatable, save :: RD_gas         (:,:) ! gas species   volume mixing ratio [ppmv]
  real(RP), private, allocatable, save :: RD_cfc         (:,:) ! CFCs          volume mixing ratio [ppmv]
  real(RP), private, allocatable, save :: RD_aerosol_conc(:,:) ! cloud/aerosol volume mixing ratio [ppmv]
  real(RP), private, allocatable, save :: RD_aerosol_radi(:,:) ! cloud/aerosol effective radius    [cm]
  real(RP), private, allocatable, save :: RD_cldfrac     (:)   ! cloud fraction    [0-1]

  integer,  private, allocatable, save :: I_MPAE2RD      (:)   ! look-up table between input aerosol category and MSTRN particle type

  character(len=H_LONG), private :: MSTRN_GASPARA_INPUTFILE   = 'PARAG.29'     !< input file (gas parameter)
  character(len=H_LONG), private :: MSTRN_AEROPARA_INPUTFILE  = 'PARAPC.29'    !< input file (particle parameter)
  character(len=H_LONG), private :: MSTRN_HYGROPARA_INPUTFILE = 'VARDATA.RM29' !< input file (hygroscopic parameter)

  integer,  private, save      :: MSTRN_nband    = 29 !< # of wave bands

  logical,  private, save      :: MSTRN_single   = .false. !< # single radiation

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
  integer,  private, parameter :: MSTRN_nptype   = 11 !< # of particle species
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
  integer,  private, parameter :: MSTRN_nradius  =  6 !< # of radius mode for hygroscopic parameter
  integer,  private, parameter :: MSTRN_ncloud   =  2 !< # of cloud types [ClearSky/Cloud]



  real(RP), private, allocatable, save :: waveh   (:)         ! wavenumbers at band boundary [1/cm]

  real(RP), private, allocatable, save :: logfitP (:)         ! fitting point for log10(pressure)
  real(RP), private, allocatable, save :: fitT    (:)         ! fitting point for temperature
  real(RP), private, allocatable, save :: logfitT (:)         ! fitting point for log10(temperature)
  integer,  private, allocatable, save :: iflgb   (:,:)       ! optical properties flag   in each band
  integer,  private, allocatable, save :: nch     (:)         ! number  of subintervals   in each band
  real(RP), private, allocatable, save :: wgtch   (:,:)       ! weights of subintervals   in each band
  integer,  private, allocatable, save :: ngasabs (:)         ! number  of absorbers(gas) in each band
  integer,  private, allocatable, save :: igasabs (:,:)       ! index   of absorbers(gas) in each band

  real(RP), private, allocatable, save :: akd     (:,:,:,:,:) ! absorption coefficient table
  real(RP), private, allocatable, save :: skd     (:,:,:,:)   ! absorption coefficient table for H2O self broadening
  real(RP), private, allocatable, save :: acfc    (:,:)       ! absorption coefficient table for CFC

  real(RP), private, allocatable, save :: fitPLK  (:,:)       ! fitting point for planck function
  real(RP), private, allocatable, save :: fsol    (:)         ! solar insolation    in each band
  real(RP), private,              save :: fsol_tot            ! total solar insolation
  real(RP), private, allocatable, save :: sfc     (:,:)       ! surface condition   in each band
  real(RP), private, allocatable, save :: rayleigh(:)         ! rayleigh scattering in each band
  real(RP), private, allocatable, save :: qmol    (:,:)       ! moments for rayleigh scattering phase function
  real(RP), private, allocatable, save :: q       (:,:,:,:)   ! moments for aerosol  scattering phase function

  integer,  private, allocatable, save :: hygro_flag(:)       ! flag for hygroscopic enlargement
  real(RP), private, allocatable, save :: radmode   (:,:)     ! radius mode for hygroscopic parameter



  ! index for optical flag iflgb
  integer,  private, parameter :: I_SWLW          = 4
  integer,  private, parameter :: I_H2O_continuum = 5
  integer,  private, parameter :: I_CFC_continuum = 7
  ! for cloud type
  integer,  private, parameter :: I_ClearSky = 1
  integer,  private, parameter :: I_Cloud    = 2

  ! pre-calc
  real(RP), private, save :: RRHO_std         ! 1 / rho(0C,1atm) * 100 [cm*m2/kg]

  real(RP), private, save :: M(2)             ! discrete quadrature mu for two-stream approximation
  real(RP), private, save :: W(2)             ! discrete quadrature w  for two-stream approximation
  real(RP), private, save :: Wmns(2), Wpls(2) ! W-, W+
  real(RP), private, save :: Wbar(2), Wscale(2)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_mstrnX_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_time, only: &
       TIME_NOWDATE
    use scale_grid_real, only: &
       REAL_lat
    use scale_atmos_solarins, only: &
       ATMOS_SOLARINS_setup
    use scale_atmos_phy_rd_profile, only: &
       RD_PROFILE_setup            => ATMOS_PHY_RD_PROFILE_setup,            &
       RD_PROFILE_setup_zgrid      => ATMOS_PHY_RD_PROFILE_setup_zgrid,      &
       RD_PROFILE_read_climatorogy => ATMOS_PHY_RD_PROFILE_read_climatology
    implicit none

    character(len=H_SHORT), intent(in) :: RD_TYPE

    integer               :: ATMOS_PHY_RD_MSTRN_KADD
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME
    character(len=H_LONG) :: ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME
    integer               :: ATMOS_PHY_RD_MSTRN_nband
    logical               :: ATMOS_PHY_RD_MSTRN_single

    namelist / PARAM_ATMOS_PHY_RD_MSTRN / &
       ATMOS_PHY_RD_MSTRN_KADD,                  &
       ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME,   &
       ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME,  &
       ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME, &
       ATMOS_PHY_RD_MSTRN_nband                , &
       ATMOS_PHY_RD_MSTRN_single

    integer :: ngas, ncfc
    integer :: ihydro, iaero
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ MstrnX radiation process'

    if ( RD_TYPE /= 'MSTRNX' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx RD_TYPE is not MSTRNX. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    ATMOS_PHY_RD_MSTRN_KADD                  = RD_KADD
    ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME   = MSTRN_GASPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME  = MSTRN_AEROPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = MSTRN_HYGROPARA_INPUTFILE
    ATMOS_PHY_RD_MSTRN_nband                 = MSTRN_nband
    ATMOS_PHY_RD_MSTRN_single                = MSTRN_single

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_MSTRN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD_MSTRN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_RD_MSTRN)

    RD_KADD                   = ATMOS_PHY_RD_MSTRN_KADD
    MSTRN_GASPARA_INPUTFILE   = ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME
    MSTRN_AEROPARA_INPUTFILE  = ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME
    MSTRN_HYGROPARA_INPUTFILE = ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME
    MSTRN_nband               = ATMOS_PHY_RD_MSTRN_nband
    MSTRN_single              = ATMOS_PHY_RD_MSTRN_single

    !--- setup solar insolation
    call ATMOS_SOLARINS_setup( TIME_NOWDATE(1) )

    !--- setup MSTRN parameter
    call RD_MSTRN_setup( ngas, & ! [OUT]
                         ncfc  ) ! [OUT]

    !--- setup climatological profile
    call RD_PROFILE_setup

    RD_KMAX  = KMAX + RD_KADD
    RD_naero = MP_QA + AE_QA
    RD_hydro_str = 1
    RD_hydro_end = MP_QA
    RD_aero_str  = MP_QA + 1
    RD_aero_end  = MP_QA + AE_QA

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

    allocate( I_MPAE2RD(RD_naero) )

    ! make look-up table between hydrometeor tracer and mstrn particle type
    do ihydro = 1, MP_QA
       I_MPAE2RD(ihydro) = I_MP2RD(ihydro)
    enddo
    ! make look-up table between aerosol     tracer and mstrn particle type
    do iaero = 1, AE_QA
       I_MPAE2RD(MP_QA+iaero) = I_AE2RD(iaero)
    enddo

    !--- setup vartical grid for radiation (larger TOA than LES domain)
    call RD_PROFILE_setup_zgrid( RD_KMAX,  RD_KADD, & ! [IN]
                                 RD_zh(:), RD_z(:)  ) ! [INOUT]

    !--- read climatological profile
    call RD_PROFILE_read_climatorogy( RD_KMAX,                & ! [IN]
                                      ngas,                   & ! [IN]
                                      ncfc,                   & ! [IN]
                                      RD_naero,               & ! [IN]
                                      REAL_lat       (IS,JS), & ! [IN], tentative treatment
                                      TIME_NOWDATE   (:),     & ! [IN]
                                      RD_zh          (:),     & ! [IN]
                                      RD_z           (:),     & ! [IN]
                                      RD_rhodz       (:),     & ! [INOUT]
                                      RD_pres        (:),     & ! [INOUT]
                                      RD_presh       (:),     & ! [INOUT]
                                      RD_temp        (:),     & ! [INOUT]
                                      RD_temph       (:),     & ! [INOUT]
                                      RD_gas         (:,:),   & ! [INOUT]
                                      RD_cfc         (:,:),   & ! [INOUT]
                                      RD_aerosol_conc(:,:),   & ! [INOUT]
                                      RD_aerosol_radi(:,:),   & ! [INOUT]
                                      RD_cldfrac     (:)      ) ! [INOUT]

    return
  end subroutine ATMOS_PHY_RD_mstrnX_setup

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_mstrnX( &
       flux_rad, flux_rad_top,    & ! [out]
       solins, cosSZA,            & ! [out]
       DENS, RHOT, QTRC,          & ! [in]
       temp_sfc, param_sfc,       & ! [in]
       CZ, FZ, CDZ, RCDZ,         & ! [in]
       REAL_lon, REAL_lat,        & ! [in]
       TIME_NOWDATE               ) ! [in]
    use scale_const, only: &
       Mdry => CONST_Mdry, &
       Mvap => CONST_Mvap, &
       PPM  => CONST_PPM
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    use scale_atmos_solarins, only: &
       ATMOS_SOLARINS_insolation
    use scale_atmos_phy_mp, only: &
       MP_EffectiveRadius => ATMOS_PHY_MP_EffectiveRadius, &
       MP_CloudFraction   => ATMOS_PHY_MP_CloudFraction,   &
       MP_Mixingratio     => ATMOS_PHY_MP_Mixingratio,   &
       MP_DENS            => ATMOS_PHY_MP_DENS
    use scale_atmos_phy_ae, only: &
       AE_EffectiveRadius => ATMOS_PHY_AE_EffectiveRadius, &
       AE_DENS
    implicit none

    real(RP), intent(out) :: flux_rad(KA,IA,JA,2,2)
    real(RP), intent(out) :: flux_rad_top(IA,JA,2)
    real(RP), intent(out) :: solins(IA,JA)
    real(RP), intent(out) :: cosSZA(IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)  :: temp_sfc(IA,JA)
    real(RP), intent(in)  :: param_sfc(5)
    real(RP), intent(in)  :: CZ(KA)
    real(RP), intent(in)  :: FZ(KA-1)
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: REAL_lon(IA,JA)
    real(RP), intent(in)  :: REAL_lat(IA,JA)
    integer , intent(in)  :: TIME_NOWDATE(6)

    real(RP) :: temp   (KA,IA,JA)
    real(RP) :: pres   (KA,IA,JA)
    real(RP) :: qsat   (KA,IA,JA)
    real(RP) :: rh     (KA,IA,JA)
    real(RP) :: cldfrac(KA,IA,JA)
    real(RP) :: MP_Re  (KA,IA,JA,MP_QA)
    real(RP) :: MP_Qe  (KA,IA,JA,MP_QA)
    real(RP) :: AE_Re  (KA,IA,JA,AE_QA)

    real(RP), parameter :: min_cldfrac = 1.E-8_RP

    ! from solar insolation
    real(RP) :: Re_factor       ! The sun-Earth distance factor

    real(RP) :: rhodz_merge       (RD_KMAX)
    real(RP) :: pres_merge        (RD_KMAX)
    real(RP) :: temp_merge        (RD_KMAX)
    real(RP) :: temph_merge       (RD_KMAX+1)

    real(RP) :: gas_merge         (RD_KMAX,MSTRN_ngas)
    real(RP) :: cfc_merge         (RD_KMAX,MSTRN_ncfc)
    real(RP) :: aerosol_conc_merge(RD_KMAX,RD_naero  )
    real(RP) :: aerosol_radi_merge(RD_KMAX,RD_naero  )
    real(RP) :: cldfrac_merge     (RD_KMAX)

    ! output
    real(RP) :: flux_rad_merge(RD_KMAX+1,2,2)

    integer :: ihydro, iaero
    integer :: RD_k, k, i, j
    !---------------------------------------------------------------------------

     if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Radiation(mstrnX)'

     call THERMODYN_temp_pres( temp(:,:,:),  & ! [OUT]
                               pres(:,:,:),  & ! [OUT]
                               DENS(:,:,:),  & ! [IN]
                               RHOT(:,:,:),  & ! [IN]
                               QTRC(:,:,:,:) ) ! [IN]

     call SATURATION_dens2qsat_liq( qsat(:,:,:), & ! [OUT]
                                    TEMP(:,:,:), & ! [IN]
                                    DENS(:,:,:)  ) ! [IN]

     do j  = JS, JE
     do i  = IS, IE
     do k  = KS, KE
        rh(k,i,j) = QTRC(k,i,j,I_QV) / qsat(k,i,j)
     enddo
     enddo
     enddo

     call MP_CloudFraction( cldfrac(:,:,:), & ! [OUT]
                            QTRC(:,:,:,:)   ) ! [IN]

     call MP_EffectiveRadius( MP_Re(:,:,:,:), & ! [OUT]
                              QTRC(:,:,:,:),  & ! [IN]
                              DENS(:,:,:)     ) ! [IN]

     call AE_EffectiveRadius( AE_Re(:,:,:,:), & ! [OUT]
                              QTRC(:,:,:,:),  & ! [IN]
                              rh  (:,:,:)     ) ! [IN]

     call MP_Mixingratio( MP_Qe(:,:,:,:),    &  ! [OUT]
                          QTRC(:,:,:,:)      )  ! [IN]


     do j = JS, JE
     do i = IS, IE

       call ATMOS_SOLARINS_insolation( solins(i,j),  & ! [OUT]
                                       cosSZA(i,j),  & ! [OUT]
                                       Re_factor,      & ! [OUT]
                                       REAL_lon(i,j),  & ! [IN]
                                       REAL_lat(i,j),  & ! [IN]
                                       TIME_NOWDATE(:) ) ! [IN]

       ! marge basic profile and value in LES domain
       rhodz_merge(:) = RD_rhodz(:)
       pres_merge (:) = RD_pres (:)
       temp_merge (:) = RD_temp (:)
       temph_merge(:) = RD_temph(:)

       gas_merge         (:,:) = RD_gas         (:,:)
       cfc_merge         (:,:) = RD_cfc         (:,:)
       aerosol_conc_merge(:,:) = RD_aerosol_conc(:,:)
       aerosol_radi_merge(:,:) = RD_aerosol_radi(:,:)
       cldfrac_merge     (:)   = RD_cldfrac     (:)

       do k = KS, KE-1
          RD_k = RD_KMAX - ( k - KS ) ! reverse axis

          temph_merge(RD_k) = 0.5_RP * ( temp(k,i,j) + temp(k+1,i,j) )
       enddo
       temph_merge(RD_KMAX+1) = temp(KS,i,j)

       do k = KS, KE
          RD_k = RD_KMAX - ( k - KS ) ! reverse axis

          rhodz_merge(RD_k)   = dens(k,i,j) * CDZ(k)   ! [kg/m2]
          pres_merge (RD_k)   = pres(k,i,j) * 1.E-2_RP ! [hPa]
          temp_merge (RD_k)   = temp(k,i,j)

          gas_merge  (RD_k,1) = QTRC(k,i,j,I_QV) / Mvap * Mdry / PPM ! [PPM]
       enddo

       do k = KS, KE
          RD_k = RD_KMAX - ( k - KS ) ! reverse axis

          cldfrac_merge(RD_k) = 0.5_RP + sign( 0.5_RP, cldfrac(k,i,j)-min_cldfrac )
       enddo

       do ihydro = 1,  MP_QA
!          if ( I_MP2ALL(ihydro) > 0 ) then
             do k = KS, KE
                RD_k = RD_KMAX - ( k - KS ) ! reverse axis
                aerosol_conc_merge(RD_k,ihydro) = MP_Qe(k,i,j,ihydro) &
                                                / MP_DENS(ihydro) * DENS(k,i,j) / PPM ! [PPM]
                aerosol_radi_merge(RD_k,ihydro) = MP_Re(k,i,j,ihydro)
             enddo
!          endif
       enddo

       do iaero = 1,  AE_QA
          if ( I_AE2ALL(iaero) > 0 ) then
             do k = KS, KE
                RD_k = RD_KMAX - ( k - KS ) ! reverse axis

                aerosol_conc_merge(RD_k,MP_QA+iaero) = QTRC(k,i,j,I_AE2ALL(iaero)) &
                                                     / AE_DENS(iaero) * DENS(k,i,j) / PPM ! [PPM]
                aerosol_radi_merge(RD_k,MP_QA+iaero) = AE_Re(k,i,j,iaero)
             enddo
          endif
       enddo

       ! calc radiative transfer
       call RD_MSTRN_DTRN3( RD_KMAX,                  & ! [IN]
                            MSTRN_ngas,               & ! [IN]
                            MSTRN_ncfc,               & ! [IN]
                            RD_naero,                 & ! [IN]
                            RD_hydro_str,             & ! [IN]
                            RD_hydro_end,             & ! [IN]
                            RD_aero_str,              & ! [IN]
                            RD_aero_end,              & ! [IN]
                            solins(i,j),              & ! [IN]
                            cosSZA(i,j),              & ! [IN]
                            rhodz_merge       (:),    & ! [IN]
                            pres_merge        (:),    & ! [IN]
                            temp_merge        (:),    & ! [IN]
                            temph_merge       (:),    & ! [IN]
                            temp_sfc(i,j),            & ! [IN]
                            gas_merge         (:,:),  & ! [IN]
                            cfc_merge         (:,:),  & ! [IN]
                            aerosol_conc_merge(:,:),  & ! [IN]
                            aerosol_radi_merge(:,:),  & ! [IN]
                            I_MPAE2RD         (:),    & ! [IN]
                            cldfrac_merge     (:),    & ! [IN]
                            param_sfc         (:),    & ! [IN]
                            flux_rad_merge    (:,:,:) ) ! [OUT]

       ! return to grid coordinate of LES domain
       do k = KS-1, KE
          RD_k = RD_KMAX - ( k - KS ) ! reverse axis

          flux_rad(k,i,j,I_LW,I_up) = flux_rad_merge(RD_k,I_LW,I_up)
          flux_rad(k,i,j,I_LW,I_dn) = flux_rad_merge(RD_k,I_LW,I_dn)
          flux_rad(k,i,j,I_SW,I_up) = flux_rad_merge(RD_k,I_SW,I_up)
          flux_rad(k,i,j,I_SW,I_dn) = flux_rad_merge(RD_k,I_SW,I_dn)
       enddo

       flux_rad_top(i,j,I_LW) = flux_rad_merge(1,I_LW,I_up)-flux_rad_merge(1,I_LW,I_dn)
       flux_rad_top(i,j,I_SW) = flux_rad_merge(1,I_SW,I_up)-flux_rad_merge(1,I_SW,I_dn)

     enddo
     enddo

     ! single radiation
     if( MSTRN_single ) then
        do k = KS-1, KE
           flux_rad(k,:,:,I_LW,I_up) = flux_rad(k,IS,JS,I_LW,I_up)
           flux_rad(k,:,:,I_LW,I_dn) = flux_rad(k,IS,JS,I_LW,I_dn)
           flux_rad(k,:,:,I_SW,I_up) = flux_rad(k,IS,JS,I_SW,I_up)
           flux_rad(k,:,:,I_SW,I_dn) = flux_rad(k,IS,JS,I_SW,I_dn)
        end do

        flux_rad_top(:,:,I_LW) = flux_rad_top(IS,JS,I_LW)
        flux_rad_top(:,:,I_SW) = flux_rad_top(IS,JS,I_SW)

        solins(:,:) = solins(IS,JS)
        cosSZA(:,:) = cosSZA(IS,JS)
     endif

    return
  end subroutine ATMOS_PHY_RD_mstrnX

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
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Input data file does not found! ', trim(MSTRN_GASPARA_INPUTFILE)
          stop
       endif

       ! read gas parameters for check
       read(fid,*) nband, nstream, nfitP, nfitT, nflag, ncfc

       if ( nband /= MSTRN_nband ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nband(given,file)=', MSTRN_nband, nband
          stop
       endif
       if ( nstream /= MSTRN_nstream ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nstream(given,file)=', MSTRN_nstream, nstream
          stop
       endif
       if ( nfitP /= MSTRN_nfitP ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nfitP(given,file)=', MSTRN_nfitP, nfitP
          stop
       endif
       if ( nfitT /= MSTRN_nfitT ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nfitT(given,file)=', MSTRN_nfitT, nfitT
          stop
       endif
       if ( nflag /= MSTRN_nflag ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nflag(given,file)=', MSTRN_nflag, nflag
          stop
       endif
       if ( ncfc /= MSTRN_ncfc ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! ncfc(given,file)=', MSTRN_ncfc, ncfc
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

    allocate( qmol    (                           MSTRN_nmoment,MSTRN_nband) )
    allocate( q       (MSTRN_nradius,MSTRN_nptype,MSTRN_nmoment,MSTRN_nband) )

    open( fid,                                     &
          file   = trim(MSTRN_AEROPARA_INPUTFILE), &
          form   = 'formatted',                    &
          status = 'old',                          &
          iostat = ierr                            )

       if ( ierr /= 0 ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Input data file does not found! ', trim(MSTRN_AEROPARA_INPUTFILE)
          stop
       endif

       ! read aerosol/surface parameters for check
       read(fid,*) nband, nsfc, nptype, nstream, nplkord, nfitPLK

       if ( nband /= MSTRN_nband ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nband(given,file)=', MSTRN_nband, nband
          stop
       endif
       if ( nsfc /= MSTRN_nsfc ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nsfc(given,file)=', MSTRN_nsfc, nsfc
          stop
       endif
       if ( nptype /= MSTRN_nptype ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nptype(given,file)=', MSTRN_nptype, nptype
          stop
       endif
       if ( nstream /= MSTRN_nstream ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nstream(given,file)=', MSTRN_nstream, nstream
          stop
       endif
       if ( nplkord /= MSTRN_nplkord ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nplkord(given,file)=', MSTRN_nplkord, nplkord
          stop
       endif
       if ( nfitPLK /= MSTRN_nfitPLK ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nfitPLK(given,file)=', MSTRN_nfitPLK, nfitPLK
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
                read(fid,*) q(:,iptype,im,iw)
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
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Input data file does not found! ', trim(MSTRN_HYGROPARA_INPUTFILE)
          stop
       endif

       read(fid,*) nptype

       if ( nptype /= MSTRN_nptype ) then
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nptype(given,file)=', MSTRN_nptype, nptype
          stop
       endif

       do iptype = 1, nptype
          read(fid,*) dummy
          read(fid,*) hygro_flag(iptype), nradius

          if ( nradius /= MSTRN_nradius ) then
             if( IO_L ) write(IO_FID_LOG,*) 'xxx Inconsistent parameter value! nradius(given,file)=', MSTRN_nradius, nradius
             stop
          endif

          read(fid,*) radmode(iptype,:)
       enddo

    close(fid)

    !----- report data -----
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A)')       '*** insolation parameters (mstrn)'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.7)') '*** total solar insolation : ', fsol_tot

    !---< constant parameter for main scheme >---
    RRHO_std = Rdry * TEM00 / Pstd * 100.0_RP ! [cm*m2/kg]

    M   (I_SW) = 1.0_RP / sqrt(3.0_RP)
    W   (I_SW) = 1.0_RP
    Wbar(I_SW) = 1.0_RP

    M   (I_LW) = 1.0_RP / 1.66_RP
    W   (I_LW) = 1.0_RP
    Wbar(I_LW) = 2.0_RP * M(I_LW)

    Wmns  (:) = sqrt( W(:) / M(:) )
    Wpls  (:) = sqrt( W(:) * M(:) )
    Wscale(:) = Wpls(:) / Wbar(:)

    return
  end subroutine RD_MSTRN_setup

  !-----------------------------------------------------------------------------
  !> DTRN v3.2
  subroutine RD_MSTRN_DTRN3( &
       kmax,         &
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
       param_sfc,    &
       rflux         )
    use scale_const, only: &
         EPS  => CONST_EPS, &
         GRAV => CONST_GRAV, &
         Pstd => CONST_Pstd, &
         PPM  => CONST_PPM
    implicit none

    integer,  intent(in)  :: kmax
    integer,  intent(in)  :: ngas
    integer,  intent(in)  :: ncfc
    integer,  intent(in)  :: naero
    integer,  intent(in)  :: hydro_str
    integer,  intent(in)  :: hydro_end
    integer,  intent(in)  :: aero_str
    integer,  intent(in)  :: aero_end
    real(RP), intent(in)  :: solins
    real(RP), intent(in)  :: cosSZA
    real(RP), intent(in)  :: rhodz       (kmax)
    real(RP), intent(in)  :: pres        (kmax)
    real(RP), intent(in)  :: temp        (kmax)
    real(RP), intent(in)  :: temph       (kmax+1)
    real(RP), intent(in)  :: temp_sfc
    real(RP), intent(in)  :: gas         (kmax,ngas )
    real(RP), intent(in)  :: cfc         (kmax,ncfc )
    real(RP), intent(in)  :: aerosol_conc(kmax,naero)
    real(RP), intent(in)  :: aerosol_radi(kmax,naero)
    integer,  intent(in)  :: aero2ptype  (naero)
    real(RP), intent(in)  :: cldfrac     (kmax)
    real(RP), intent(in)  :: param_sfc   (5)
    real(RP), intent(out) :: rflux       (kmax+1,2,2)

    ! for planck functions
    real(RP) :: bbar (kmax,  MSTRN_nband) ! planck functions for thermal source at the interface
    real(RP) :: bbarh(kmax+1,MSTRN_nband) ! planck functions for thermal source at the center
    real(RP) :: b_sfc(       MSTRN_nband) ! planck functions for thermal source at the surface
    real(RP) :: wl, beta

    ! for P-T fitting
    real(RP) :: logP(kmax)   ! log10(pres)
    real(RP) :: logT(kmax)   ! log10(temp)
    real(RP) :: dz_std(kmax) ! layer thickness at 0C, 1atm [cm]

    integer  :: indexP (kmax)       ! index for interpolation in P-fitting
    real(RP) :: factP  (kmax)       ! interpolation factor    in P-fitting
    real(RP) :: factT32(kmax)       ! interpolation factor    in T-fitting
    real(RP) :: factT21(kmax)       ! interpolation factor    in T-fitting
    integer  :: indexR (naero,kmax) ! index for interpolation in R-fitting
    real(RP) :: factR  (naero,kmax) ! interpolation factor    in R-fitting

    ! for optical thickness by gas
    real(RP) :: tauGAS(kmax,MSTRN_ch_limit,MSTRN_nband) ! optical thickness by gas absorption (total)
    real(RP) :: tauKD (MSTRN_ch_limit)                  ! optical thickness by gas absorption (gases)
    real(RP) :: tauCON(MSTRN_ch_limit)                  ! optical thickness by gas absorption (H2O continuum)
    real(RP) :: tauCFC                                  ! optical thickness by gas absorption (CFC)
    real(RP) :: A1, A2, A3, factPT
    real(RP) :: qv
    integer  :: gasno

    ! for optical thickness by particles
    real(RP) :: tauPR   (kmax,MSTRN_ncloud)          ! optical thickness        by Rayleigh/cloud/aerosol
    real(RP) :: omgPR   (kmax,MSTRN_ncloud)          ! single scattering albedo by Rayleigh/cloud/aerosol
    real(RP) :: optparam(MSTRN_nmoment,MSTRN_ncloud) ! optical parameters
    real(RP) :: q_fit, dp_P

    ! for albedo
    real(RP) :: albedo_sfc(MSTRN_ncloud,MSTRN_nband) ! surface albedo
    real(RP) :: tau_column
    integer  :: luindex

    ! for two-stream
    real(RP) :: tau(    kmax,MSTRN_ncloud) ! total optical thickness
    real(RP) :: omg(    kmax,MSTRN_ncloud) ! single scattering albedo
    real(RP) :: g  (0:2,kmax,MSTRN_ncloud) ! two-stream approximation factors
                                           ! 0: always 1
                                           ! 1: asymmetry factor
                                           ! 2: truncation factor
    real(RP) :: b  (0:2,kmax,MSTRN_ncloud) ! planck expansion coefficients (zero if SW)
    real(RP) :: fsol_rgn                   ! solar insolation              (zero if LW)

    real(RP) :: flux       (kmax+1,2)      ! upward/downward flux
    real(RP) :: flux_direct(kmax+1)        ! downward flux (direct solar)

    real(RP) :: zerosw
    integer  :: ip, ir, irgn
    integer  :: igas, icfc, iaero, iptype
    integer  :: iw, ich, k, iplk, icloud, im
    !---------------------------------------------------------------------------

    !---< set planck functions >---
    do iw = 1, MSTRN_nband
       if ( iflgb(I_SWLW,iw)+1 == I_SW ) then

          bbar (:,iw) = 0.0_RP
          bbarh(:,iw) = 0.0_RP
          b_sfc(iw)   = 0.0_RP

       elseif( iflgb(I_SWLW,iw)+1 == I_LW ) then

          wl = 10000.0_RP / sqrt( waveh(iw) * waveh(iw+1) )

          ! from temp at cell center
          do k = 1, kmax
             beta = 0.0_RP
             do iplk = MSTRN_nfitPLK, 1, -1
                beta = beta / ( wl*temp(k) ) + fitPLK(iplk,iw)
             enddo
             bbar(k,iw) = exp(-beta) * temp(k) / (wl*wl)
          enddo

          ! from temp at cell wall
          do k = 1, kmax+1
             beta = 0.0_RP
             do iplk = MSTRN_nfitPLK, 1, -1
                beta = beta / ( wl*temph(k) ) + fitPLK(iplk,iw)
             enddo
             bbarh(k,iw) = exp(-beta) * temph(k) / (wl*wl)
          enddo

          ! from temp_sfc
          beta = 0.0_RP
          do iplk = MSTRN_nfitPLK, 1, -1
             beta = beta / ( wl*temp_sfc ) + fitPLK(iplk,iw)
          enddo
          b_sfc(iw) = exp(-beta) * temp_sfc / (wl*wl)

       endif
    enddo

    !---< interpolation of gas parameters (P-T fitting) >---
    do k = 1, kmax
       logP(k) = log10( pres(k) )
       logT(k) = log10( temp(k) )

       dz_std(k) = rhodz(k) * RRHO_std ! [cm]
    enddo

    do k = 1, kmax
       indexP(k) = MSTRN_nfitP
       do ip = MSTRN_nfitP, 2, -1
          if( logP(k) >= logfitP(ip) ) indexP(k) = ip
       enddo
    enddo

    do k = 1, kmax
       ip = indexP(k)
       factP(k) = ( logP(k) - logfitP(ip-1) ) / ( logfitP(ip) - logfitP(ip-1) )

       factT32(k) = ( logT(k) - logfitT(2) ) / ( logfitT(3) - logfitT(2) ) &
                  * ( temp(k) - fitT(1)    ) / ( fitT(3)    - fitT(1)    )
       factT21(k) = ( logT(k) - logfitT(2) ) / ( logfitT(2) - logfitT(1) ) &
                  * ( fitT(3) - temp(k)    ) / ( fitT(3)    - fitT(1)    )
    enddo

    tauGAS(:,:,:) = 0.0_RP

    do iw = 1, MSTRN_nband
       do k = 1, kmax
          ip = indexP(k)

          tauKD (:) = 0.0_RP
          tauCON(:) = 0.0_RP
          tauCFC    = 0.0_RP

          do igas = 1, ngasabs(iw)
             gasno = igasabs(igas,iw)

             do ich  = 1, nch(iw)
                A1 = AKD(ich,ip-1,1,gasno,iw) * ( 1.0_RP - factP(k) )&
                   + AKD(ich,ip  ,1,gasno,iw) * (          factP(k) )
                A2 = AKD(ich,ip-1,2,gasno,iw) * ( 1.0_RP - factP(k) )&
                   + AKD(ich,ip  ,2,gasno,iw) * (          factP(k) )
                A3 = AKD(ich,ip-1,3,gasno,iw) * ( 1.0_RP - factP(k) )&
                   + AKD(ich,ip  ,3,gasno,iw) * (          factP(k) )

                ! is this necessary?
                if (       A2 <= -20.0_RP &
                     .AND. A1 >  -20.0_RP &
                     .AND. A3 >  -20.0_RP ) then
                   print *, "found!", A1, A2, A3, (A1+A3)*0.5_RP
                endif

                factPT = factT32(k)*(A3-A2) + A2 + factT21(k)*(A2-A1)

                tauKD(ich) = tauKD(ich) &
                           + 10.0_RP**factPT * gas(k,igasabs(igas,iw)) * PPM * dz_std(k)
             enddo
          enddo

          if ( iflgb(I_H2O_continuum,iw) == 1 ) then
             qv = gas(k,1) * PPM * dz_std(k)

             do ich = 1, nch(iw)
                A1 = SKD(ich,ip-1,1,iw) * ( 1.0_RP-factP(k) )&
                   + SKD(ich,ip  ,1,iw) * (        factP(k) )
                A2 = SKD(ich,ip-1,2,iw) * ( 1.0_RP-factP(k) )&
                   + SKD(ich,ip  ,2,iw) * (        factP(k) )
                A3 = SKD(ich,ip-1,3,iw) * ( 1.0_RP-factP(k) )&
                   + SKD(ich,ip  ,3,iw) * (        factP(k) )

                factPT = factT32(k)*(A3-A2) + A2 + factT21(k)*(A2-A1)

                tauCON(ich) = 10.0_RP**factPT * qv*qv / ( qv + dz_std(k) )
             enddo
          endif

          if ( iflgb(I_CFC_continuum,iw) == 1 ) then
             tauCFC = 0.0_RP
             do icfc = 1, ncfc
                tauCFC = tauCFC + 10.0_RP**acfc(icfc,iw) * cfc(k,icfc) * dz_std(k)
             enddo
          endif

          do ich = 1, nch(iw)
             tauGAS(k,ich,iw) = tauKD(ich) + tauCON(ich) + tauCFC
          enddo

       enddo
    enddo

    !---< interpolation of mode radius & hygroscopic parameter (R-fitting) >---
    do iaero = 1, naero
       iptype = aero2ptype(iaero)
       do k = 1, kmax
          if ( aerosol_radi(k,iaero) <= radmode(iptype,1) ) then ! extrapolation

             ir = 1
             indexR(iaero,k) = ir
             factR (iaero,k) = ( radmode(iptype,ir  ) - aerosol_radi(k,iaero) ) &
                             / ( radmode(iptype,ir+1) - radmode(iptype,ir)    )

          elseif( aerosol_radi(k,iaero) > radmode(iptype,MSTRN_nradius) ) then ! extrapolation

             ir = MSTRN_nradius
             indexR(iaero,k) = ir
             factR (iaero,k) = ( radmode(iptype,ir) - aerosol_radi(k,iaero) ) &
                             / ( radmode(iptype,ir) - radmode(iptype,ir-1)  )

          else
             do ir = 1, MSTRN_nradius-1
                if (       aerosol_radi(k,iaero) <= radmode(iptype,ir+1) &
                     .AND. aerosol_radi(k,iaero) >  radmode(iptype,ir  ) ) then ! interpolation

                   indexR(iaero,k) = ir
                   factR (iaero,k) = ( aerosol_radi(k,iptype) - radmode(iptype,ir) ) &
                                   / ( radmode(iptype,ir+1)   - radmode(iptype,ir) )

                endif
             enddo
          endif
       enddo
    enddo

    ! initialize
    rflux(:,:,:) = 0.0_RP

    do iw = 1, MSTRN_nband

       !--- particle (Rayleigh/Cloud/Aerosol)
       do k = 1, kmax
          dp_P = rhodz(k) * GRAV / Pstd

          ! optical thickness, phase function
          ! im=1: extinction coefficient
          ! im=2,3,4: moments of the volume scattering phase function
          do im = 1, MSTRN_nstream*2+2

             !--- Rayleigh scattering

             optparam(im,I_Cloud   ) = rayleigh(iw) * qmol(im,iw) * dp_P
             optparam(im,I_ClearSky) = rayleigh(iw) * qmol(im,iw) * dp_P

             !--- Cloud scattering
             do iaero = hydro_str, hydro_end
                iptype = aero2ptype(iaero)
                ir     = indexR(iaero,k)

                q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(iaero,k) ) &
                      + q(ir+1,iptype,im,iw) * (        factR(iaero,k) )

                optparam(im,I_Cloud) = optparam(im,I_Cloud) + q_fit * aerosol_conc(k,iaero) * PPM * dz_std(k)
             enddo

             !--- Aerosol scattering
             do iaero = aero_str, aero_end
                iptype = aero2ptype(iaero)
                ir     = indexR(iaero,k)

                q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(iaero,k) ) &
                      + q(ir+1,iptype,im,iw) * (        factR(iaero,k) )

                optparam(im,I_Cloud   ) = optparam(im,I_Cloud   ) + q_fit * aerosol_conc(k,iaero) * PPM * dz_std(k)
                optparam(im,I_ClearSky) = optparam(im,I_ClearSky) + q_fit * aerosol_conc(k,iaero) * PPM * dz_std(k)
             enddo

          enddo

          do icloud = 1, MSTRN_ncloud
             tauPR(k,icloud) = optparam(1,icloud)
             omgPR(k,icloud) = optparam(1,icloud) - optparam(2,icloud)

             !--- g
             g(0,k,icloud) = 1.0_RP
             if ( omgPR(k,icloud) == 0.0_RP ) then
                g(1,k,icloud) = 0.0_RP
                g(2,k,icloud) = 0.0_RP
             else
                g(1,k,icloud) = optparam(3,icloud) / ( optparam(1,icloud) - optparam(2,icloud) )
                g(2,k,icloud) = optparam(4,icloud) / ( optparam(1,icloud) - optparam(2,icloud) )
             endif
          enddo

       enddo

       !--- Albedo
       luindex = int( param_sfc(1) )
       if ( luindex == 1 ) then ! ocean

          do icloud = 1, MSTRN_ncloud

             tau_column = 0.0_RP
             do k = 1, kmax
                tau_column = tau_column + tauPR(k,icloud) ! layer-total(for ocean albedo)
             enddo

             if ( tau_column > 0.0_RP .AND. iflgb(I_SWLW,iw)+1 == I_SW ) then
                albedo_sfc(icloud,iw) = albedo_sea(cosSZA,tau_column)
             else
                albedo_sfc(icloud,iw) = 0.05_RP
             endif

          enddo

       elseif( luindex == 2 ) then ! land surface
          albedo_sfc(:,iw) = (        param_sfc(2) ) * sfc(iw,2) & ! wet land
                           + ( 1.0_RP-param_sfc(2) ) * sfc(iw,3)   ! dry land
       else
          albedo_sfc(:,iw) = sfc(luindex,iw)
       endif

       ! sub-channel loop
       do ich = 1, nch(iw)

          !--- total tau
          do icloud = 1, 2
          do k      = 1, kmax
             tau(k,icloud) = tauGAS(k,ich,iw) + tauPR(k,icloud)
          enddo
          enddo

          !--- omega
          do icloud = 1, 2
          do k      = 1, kmax
             zerosw = 0.5_RP - sign( 0.5_RP, tau(k,icloud)-EPS ) ! if tau < EPS, zerosw = 1

             omg(k,icloud) = 1.0_RP ! initialize
             omg(k,icloud) = omgPR(k,icloud) / ( tau(k,icloud)-zerosw ) * ( 1.0_RP-zerosw )
          enddo
          enddo

          !--- bn
          if ( iflgb(I_SWLW,iw)+1 == I_SW ) then

             b(:,:,:) = 0.0_RP

             fsol_rgn = fsol(iw) / fsol_tot * solins

          elseif( iflgb(I_SWLW,iw)+1 == I_LW ) then

             do icloud = 1, 2
             do k      = 1, kmax
                b(0,k,icloud) = bbarh(k,iw)
                if ( tau(k,icloud) > 0.0_RP ) then
                   b(1,k,icloud) = ( -          bbarh(k+1,iw) &
                                     + 4.0_RP * bbar (k,  iw) &
                                     - 3.0_RP * bbarh(k,  iw) ) / tau(k,icloud)
                   b(2,k,icloud) = ( +          bbarh(k+1,iw) &
                                     - 2.0_RP * bbar (k,  iw) &
                                     +          bbarh(k,  iw) ) / (tau(k,icloud)*tau(k,icloud)) * 2.0_RP
                else
                   b(1,k,icloud) = 0.0_RP
                   b(2,k,icloud) = 0.0_RP
                endif
             enddo
             enddo

             fsol_rgn = 0.0_RP

          endif

          ! two-stream transfer
          irgn = iflgb(I_SWLW,iw) + 1

          !print *, iw, ich, sum(tau(:,1)), sum(omg(:,1)), sum(b(0,:,1)),  sum(b(1,:,1)),  sum(b(2,:,1))

          call RD_MSTRN_two_stream( kmax,               & ! [IN]
                                    cosSZA,             & ! [IN]
                                    fsol_rgn,           & ! [IN]
                                    irgn,               & ! [IN]
                                    tau        (:,:),   & ! [IN]
                                    omg        (:,:),   & ! [IN]
                                    g          (:,:,:), & ! [IN]
                                    b          (:,:,:), & ! [IN]
                                    b_sfc      (iw),    & ! [IN]
                                    albedo_sfc (:,iw),  & ! [IN]
                                    cldfrac    (:),     & ! [IN]
                                    flux       (:,:),   & ! [OUT]
                                    flux_direct(:)      ) ! [OUT]

          do k = 1, kmax+1
             rflux(k,irgn,I_up) = rflux(k,irgn,I_up) + flux(k,I_up) * wgtch(ich,iw)
             rflux(k,irgn,I_dn) = rflux(k,irgn,I_dn) + flux(k,I_dn) * wgtch(ich,iw)
          enddo

          !print *, iw, ich, sum(flux(:,1)), sum(flux(:,2)), sum(flux_direct(:)),wgtch(ich,iw)

       enddo ! ICH loop
    enddo ! IW loop

    return
  end subroutine RD_MSTRN_DTRN3

  !-----------------------------------------------------------------------------
  !> Two stream calculation CORE
  subroutine RD_MSTRN_two_stream( &
       kmax,       &
       cosSZA0,    &
       fsol,       &
       irgn,       &
       tau,        &
       omg,        &
       g,          &
       b,          &
       b_sfc,      &
       albedo_sfc, &
       cldfrac,    &
       flux,       &
       flux_direct )
    use scale_const, only: &
       PI   => CONST_PI,  &
       EPS1 => CONST_EPS1
    implicit none

    integer,  intent(in)  :: kmax
    real(RP), intent(in)  :: cosSZA0                    ! cos(SZA) = mu0
    real(RP), intent(in)  :: fsol                       ! solar radiation intensity
    integer,  intent(in)  :: irgn                       ! 1:SW 2:LW
    real(RP), intent(in)  :: tau(    kmax,MSTRN_ncloud) ! total optical thickness          (clear-sky/cloud)
    real(RP), intent(in)  :: omg(    kmax,MSTRN_ncloud) ! single scattering albedo         (clear-sky/cloud)
    real(RP), intent(in)  :: g  (0:2,kmax,MSTRN_ncloud) ! two-stream approximation factors (clear-sky/cloud)
    real(RP), intent(in)  :: b  (0:2,kmax,MSTRN_ncloud) ! planck expansion coefficients    (clear-sky/cloud)
    real(RP), intent(in)  :: b_sfc                      ! planck function at surface
    real(RP), intent(in)  :: albedo_sfc(MSTRN_ncloud)   ! surface albedo                   (clear-sky/cloud)
    real(RP), intent(in)  :: cldfrac(kmax)              ! cloud fraction

    real(RP), intent(out) :: flux       (kmax+1,2)      ! upward(sfc->TOA)/downward(TOA->sfc) flux
    real(RP), intent(out) :: flux_direct(kmax+1)        ! downward(TOA->sfc) flux, solar direct

    ! parameters with two-stream truncation
    real(RP) :: tau_new    ! optical thickness        : two-stream truncation
    real(RP) :: omg_new    ! single scattering albedo : two-stream truncation
    real(RP) :: g_new      ! asymmetric factor        : two-stream truncation
    real(RP) :: b_new(0:2) ! planck function          : two-stream truncation
    real(RP) :: c(0:2)
    real(RP) :: Pmns, Ppls ! Phase  function          : two-stream truncation
    real(RP) :: Smns, Spls ! Source function          : two-stream truncation

    ! working
    real(RP) :: cosSZA
    real(RP) :: X, Y                 ! X-, X+
    real(RP) :: lamda                ! eigenvalue of X-, X+
    real(RP) :: E
    real(RP) :: Apls_mns, Bpls_mns   ! A+/A-, B+/B-
    real(RP) :: V0mns, V0pls         ! V0-, V0+
    real(RP) :: V1mns, V1pls         ! V1-, V1+
    real(RP) :: Dmns(0:2), Dpls(0:2) ! D0-, D1-, D2-, D0+, D1+, D2+
    real(RP) :: SIGmns, SIGpls       ! sigma-, sigma+
    real(RP) :: Qgamma               ! Q * gamma

    ! main factors
    real(RP) :: Tdir0(kmax,MSTRN_ncloud) ! transmission factor for solar direct (clear-sky/cloud)
    real(RP) :: R0   (kmax,MSTRN_ncloud) ! reflection   factor                  (clear-sky/cloud)
    real(RP) :: T0   (kmax,MSTRN_ncloud) ! transmission factor                  (clear-sky/cloud)
    real(RP) :: Em_LW(kmax,MSTRN_ncloud) ! thermal source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Em_SW(kmax,MSTRN_ncloud) ! solar   source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Ep_LW(kmax,MSTRN_ncloud) ! thermal source (TOA->sfc)            (clear-sky/cloud)
    real(RP) :: Ep_SW(kmax,MSTRN_ncloud) ! solar   source (TOA->sfc)            (clear-sky/cloud)

    ! Averaged factors, considering cloud overwrap
    real(RP) :: tau_bar    ! accumulated optical thickness at each layer
    real(RP) :: R (kmax+1) ! reflection   factor
    real(RP) :: T (kmax+1) ! transmission factor
    real(RP) :: Em(kmax+1) ! source (sfc->TOA)
    real(RP) :: Ep(kmax+1) ! source (TOA->sfc)

    ! Doubling-Adding
    real(RP) :: R12mns(kmax+1), R12pls(kmax+1) ! reflection factor in doubling method
    real(RP) :: E12mns(kmax+1), E12pls(kmax+1) ! source function   in doubling method
    real(RP) :: Umns, Upls                     ! flux intensity

    real(RP) :: factor
    real(RP) :: Em0(MSTRN_ncloud)

    integer  :: k, icloud
    !---------------------------------------------------------------------------

    cosSZA = max( cosSZA0, RD_cosSZA_min )

    do icloud = 1, 2
    do k      = 1, kmax

       !---< two-stream truncation >---
       tau_new = ( 1.0_RP - omg(k,icloud)*g(2,k,icloud) ) * tau(k,icloud)

       omg_new = ( 1.0_RP - g(2,k,icloud) ) / ( 1.0_RP - omg(k,icloud)*g(2,k,icloud) ) * omg(k,icloud)
       omg_new = min( omg_new, EPS1 )

       g_new   = ( g(1,k,icloud) - g(2,k,icloud) ) / ( 1.0_RP - g(2,k,icloud) )

       Tdir0(k,icloud) = exp(-tau_new/cosSZA)

       factor   = ( 1.0_RP - omg(k,icloud)*g(2,k,icloud) )
       b_new(0) = b(0,k,icloud)
       b_new(1) = b(1,k,icloud) / factor
       b_new(2) = b(2,k,icloud) / (factor*factor)
       c(:)     = Wmns(irgn) * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new(:)

       !--- P+, P-
       Pmns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M(irgn)*M(irgn) )
       Ppls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M(irgn)*M(irgn) )

       !--- S+, S-
       Smns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M(irgn)*cosSZA )
       Spls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M(irgn)*cosSZA )

       !---< calculate R, T, e+, e- >---
       if ( tau_new <= 1.E-4_RP ) then

          !--- R, T
          R0(k,icloud) =          tau_new * (          Pmns ) / M(irgn)
          T0(k,icloud) = 1.0_RP - tau_new * ( 1.0_RP - Ppls ) / M(irgn)

          !--- thermal source
          Em_LW(k,icloud) = 0.5_RP * tau_new * ( 2.0_RP*c(0) + c(1)*tau_new + c(2)*tau_new*tau_new )
          Ep_LW(k,icloud) = 0.5_RP * tau_new * ( 2.0_RP*c(0) + c(1)*tau_new + c(2)*tau_new*tau_new )

          !--- solar source
          Em_SW(k,icloud) = Wmns(irgn) * Smns * tau_new * sqrt( Tdir0(k,icloud) )
          Ep_SW(k,icloud) = Wmns(irgn) * Spls * tau_new * sqrt( Tdir0(k,icloud) )

       else

          !--- X, Y
          X     =  ( 1.0_RP - W(irgn) * ( Ppls - Pmns ) ) / M(irgn)
          Y     =  ( 1.0_RP - W(irgn) * ( Ppls + Pmns ) ) / M(irgn)
          lamda = sqrt(X*Y)
          E     = exp(-lamda*tau_new)

          !--- A+/A-, B+/B-
          Apls_mns = ( X * ( 1.0_RP+E ) - lamda * ( 1.0_RP-E ) ) &
                   / ( X * ( 1.0_RP+E ) + lamda * ( 1.0_RP-E ) )
          Bpls_mns = ( X * ( 1.0_RP-E ) - lamda * ( 1.0_RP+E ) ) &
                   / ( X * ( 1.0_RP-E ) + lamda * ( 1.0_RP+E ) )

          !--- R, T
          R0(k,icloud) = 0.5_RP * ( Apls_mns + Bpls_mns )
          T0(k,icloud) = 0.5_RP * ( Apls_mns - Bpls_mns )

          !--- thermal source
          Dmns(0) = c(0) / Y + 2.0_RP * c(2) / (X*Y*Y) + c(1) / (X*Y)
          Dpls(0) = c(0) / Y + 2.0_RP * c(2) / (X*Y*Y) - c(1) / (X*Y)
          Dmns(1) = c(1) / Y + 2.0_RP * c(2) / (X*Y)
          Dpls(1) = c(1) / Y - 2.0_RP * c(2) / (X*Y)
          Dmns(2) = c(2) / Y
          Dpls(2) = c(2) / Y

          V0mns = Dmns(0)
          V0pls = Dpls(0)
          V1mns = Dmns(0) + Dmns(1)*tau_new + Dmns(2)*tau_new*tau_new
          V1pls = Dpls(0) + Dpls(1)*tau_new + Dpls(2)*tau_new*tau_new

          Em_LW(k,icloud) = V0mns - R0(k,icloud) * V0pls - T0(k,icloud) * V1mns
          Ep_LW(k,icloud) = V1pls - T0(k,icloud) * V0pls - R0(k,icloud) * V1mns

          !--- solar source
          SIGmns = Wmns(irgn) * ( Spls - Smns )
          SIGpls = Wmns(irgn) * ( Spls + Smns )

          Qgamma = ( SIGpls*X*cosSZA + SIGmns ) / ( X*Y*cosSZA - 1.0/cosSZA )

          V0pls = 0.5_RP * ( ( 1.0_RP + 1.0_RP/(X*cosSZA) ) * Qgamma + SIGmns / X )
          V0mns = 0.5_RP * ( ( 1.0_RP - 1.0_RP/(X*cosSZA) ) * Qgamma - SIGmns / X )

          V1pls = V0pls * Tdir0(k,icloud)
          V1mns = V0mns * Tdir0(k,icloud)

          Em_SW(k,icloud) = V0mns - R0(k,icloud) * V0pls - T0(k,icloud) * V1mns
          Ep_SW(k,icloud) = V1pls - T0(k,icloud) * V0pls - R0(k,icloud) * V1mns

       endif

    enddo ! k loop
    enddo ! cloud loop

    !---< consider partial cloud layer: semi-random over-wrapping >---
    tau_bar = 1.0_RP
    do k = 1, kmax

       R (k) = (        cldfrac(k) ) * R0(k,I_Cloud   ) &
             + ( 1.0_RP-cldfrac(k) ) * R0(k,I_ClearSky)

       T (k) = (        cldfrac(k) ) * T0(k,I_Cloud   ) &
             + ( 1.0_RP-cldfrac(k) ) * T0(k,I_ClearSky)

       Em(k) = (        cldfrac(k) ) * ( Em_LW(k,I_Cloud   ) + Em_SW(k,I_Cloud   ) * tau_bar * fsol ) &
             + ( 1.0_RP-cldfrac(k) ) * ( Em_LW(k,I_ClearSky) + Em_SW(k,I_ClearSky) * tau_bar * fsol )

       Ep(k) = (        cldfrac(k) ) * ( Ep_LW(k,I_Cloud   ) + Ep_SW(k,I_Cloud   ) * tau_bar * fsol ) &
             + ( 1.0_RP-cldfrac(k) ) * ( Ep_LW(k,I_ClearSky) + Ep_SW(k,I_ClearSky) * tau_bar * fsol )

       flux_direct(k) = cosSZA * tau_bar * fsol

       ! update tau_bar
       tau_bar = tau_bar * ( (        cldfrac(k) ) * Tdir0(k,I_Cloud   ) &
                           + ( 1.0_RP-cldfrac(k) ) * Tdir0(k,I_ClearSky) )
    enddo

    ! at lambert surface
    R (kmax+1) = (        cldfrac(kmax) ) * albedo_sfc(I_Cloud   ) &
               + ( 1.0_RP-cldfrac(kmax) ) * albedo_sfc(I_ClearSky)

    T (kmax+1) = 0.0_RP

    flux_direct(kmax+1) = cosSZA * tau_bar * fsol

    Em0(:) = Wpls(irgn) * ( flux_direct(kmax+1) * albedo_sfc(:) / (W(irgn)*M(irgn)) &
                          + 2.0_RP * PI * ( 1.0_RP-albedo_sfc(:) ) * b_sfc          )

    Em(kmax+1) = (        cldfrac(kmax) ) * Em0(I_Cloud   ) &
               + ( 1.0_RP-cldfrac(kmax) ) * Em0(I_ClearSky)

    Ep(kmax+1) = 0.0_RP

    !---< Adding-Doubling method >---
    ! [note] TOA->Surface is positive direction. "pls" means upper to lower altitude.

    ! adding: surface to TOA
    R12pls(kmax+1) = R (kmax+1)
    E12mns(kmax+1) = Em(kmax+1)

    do k = kmax, 1, -1
       R12pls(k) = R (k) + T(k) / ( 1.0_RP - R12pls(k+1) * R(k) ) * ( R12pls(k+1) * T (k)               )
       E12mns(k) = Em(k) + T(k) / ( 1.0_RP - R12pls(k+1) * R(k) ) * ( R12pls(k+1) * Ep(k) + E12mns(k+1) )
    enddo

    ! adding: TOA to surface
    R12mns(1) = R (1)
    E12pls(1) = Ep(1)

    do k = 2, kmax+1
       R12mns(k) = R (k) + T(k) / ( 1.0_RP - R12mns(k-1) * R(k) ) * ( R12mns(k-1)*T (k)               )
       E12pls(k) = Ep(k) + T(k) / ( 1.0_RP - R12mns(k-1) * R(k) ) * ( R12mns(k-1)*Em(k) + E12pls(k-1) )
    enddo

    !--- radiative flux at cell wall:
    ! [note] "d" means upper to lower altitude.

    ! TOA boundary
    flux(1,I_up) = Wscale(irgn) * E12mns(1)
    flux(1,I_dn) = flux_direct(1)

    do k = 2, kmax+1
       Upls = ( E12pls(k-1) + R12mns(k-1)*E12mns(k) ) / ( 1.0_RP - R12mns(k-1)*R12pls(k) )
       Umns =  E12mns(k) + R12pls(k) * Upls

       flux(k,I_up) = Wscale(irgn) * Umns
       flux(k,I_dn) = Wscale(irgn) * Upls + flux_direct(k)
    enddo

    return
  end subroutine RD_MSTRN_two_stream

  !-----------------------------------------------------------------------------
  ! Sea surface reflectance by Payne
  function albedo_sea( cosSZA, tau )
    implicit none

    real(RP), intent(in) :: cosSZA, tau
    real(RP)             :: albedo_sea

    real(RP) :: c(5,3)
    data c / -2.8108_RP   , -1.3651_RP,  2.9210E1_RP, -4.3907E1_RP,  1.8125E1_RP, &
              6.5626E-1_RP, -8.7253_RP, -2.7749E1_RP,  4.9486E1_RP, -1.8345E1_RP, &
             -6.5423E-1_RP,  9.9967_RP,  2.7769_RP  , -1.7620E1_RP,  7.0838_RP    /

    real(RP) :: am1, tr1, s

    integer  :: i
    !---------------------------------------------------------------------------

    am1 = max( min( cosSZA, 0.961_RP ), 0.0349_RP )

    tr1 = max( min( am1 / ( 4.0_RP * tau ), 1.0_RP ), 0.05_RP )

    s = 0.0_RP
    do i = 1, 5
       s = s + c(i,1) * tr1**(i-1)           &
             + c(i,2) * tr1**(i-1) * am1     &
             + c(i,3) * tr1**(i-1) * am1*am1
    enddo

    albedo_sea = exp(s)

    return
  end function albedo_sea

end module scale_atmos_phy_rd_mstrnX
