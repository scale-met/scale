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
  real(RP), private, parameter :: RD_EPS        = 1.E-4_RP ! minimum SZA (>89.0)

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
  subroutine ATMOS_PHY_RD_mstrnx_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_time, only: &
       TIME_NOWDATE
    use scale_grid_real, only: &
       REAL_LAT
    use scale_atmos_phy_rd_profile, only: &
       RD_PROFILE_setup       => ATMOS_PHY_RD_PROFILE_setup,       &
       RD_PROFILE_setup_zgrid => ATMOS_PHY_RD_PROFILE_setup_zgrid, &
       RD_PROFILE_read        => ATMOS_PHY_RD_PROFILE_read
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

    !--- setup MSTRN parameter
    call RD_MSTRN_setup( ngas, & ! [OUT]
                         ncfc  ) ! [OUT]

    !--- setup climatological profile
    call RD_PROFILE_setup

    RD_KMAX      = KMAX + RD_KADD
    RD_naero     = MP_QA + AE_QA
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
    call RD_PROFILE_read( RD_KMAX,                & ! [IN]
                          ngas,                   & ! [IN]
                          ncfc,                   & ! [IN]
                          RD_naero,               & ! [IN]
                          REAL_LAT       (IS,JS), & ! [IN], tentative treatment
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
       oceanfrc,              &
       temp_sfc, albedo_land, &
       solins, cosSZA,        &
       flux_rad,              &
       flux_rad_top           )
!       Jval                            )
    use scale_const, only: &
       Mdry => CONST_Mdry, &
       Mvap => CONST_Mvap, &
       PPM  => CONST_PPM
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
       AE_DENS
    implicit none

    real(RP), intent(in)  :: DENS        (KA,IA,JA)
    real(RP), intent(in)  :: RHOT        (KA,IA,JA)
    real(RP), intent(in)  :: QTRC        (KA,IA,JA,QA)
    real(RP), intent(in)  :: CZ          (  KA,IA,JA)    ! UNUSED
    real(RP), intent(in)  :: FZ          (0:KA,IA,JA)
    real(RP), intent(in)  :: oceanfrc    (IA,JA)
    real(RP), intent(in)  :: temp_sfc    (IA,JA)
    real(RP), intent(in)  :: albedo_land (IA,JA,2)
    real(RP), intent(in)  :: solins      (IA,JA)
    real(RP), intent(in)  :: cosSZA      (IA,JA)
    real(RP), intent(out) :: flux_rad    (KA,IA,JA,2,2)
    real(RP), intent(out) :: flux_rad_top(IA,JA,2)
!    real(RP), intent(out) :: Jval        (KA,IA,JA,CH_QA_photo)

    real(RP) :: temp   (KA,IA,JA)
    real(RP) :: pres   (KA,IA,JA)
    real(RP) :: qsat   (KA,IA,JA)
    real(RP) :: rh     (KA,IA,JA)
    real(RP) :: cldfrac(KA,IA,JA)
    real(RP) :: MP_Re  (KA,IA,JA,MP_QA)
    real(RP) :: MP_Qe  (KA,IA,JA,MP_QA)
    real(RP) :: AE_Re  (KA,IA,JA,AE_QA)

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
    real(RP) :: flux_rad_merge(RD_KMAX+1,IA,JA,2,2)

    integer :: ihydro, iaero
    integer :: RD_k, k, i, j, v
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Radiation(mstrnx)'

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

    ! marge basic profile and value in LES domain

    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KADD
          temph_merge(RD_k,i,j) = RD_temph(RD_k)
       enddo

       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis

          temph_merge(RD_k,i,j) = 0.5_RP * ( temp(k,i,j) + temp(k+1,i,j) )
       enddo
       temph_merge(RD_KMAX+1,i,j) = temp(KS,i,j)
    enddo
    enddo

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

    do v = 1,  MSTRN_ngas
    do j = JS, JE
    do i = IS, IE
       do RD_k = 1, RD_KMAX
          gas_merge(RD_k,i,j,v) = RD_gas(RD_k,v)
       enddo
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do RD_k = RD_KADD+1, RD_KMAX
          k = KS + RD_KMAX - RD_k ! reverse axis

          gas_merge(RD_k,i,j,1) = QTRC(k,i,j,I_QV) / Mvap * Mdry / PPM ! [PPM]
       enddo
    enddo
    enddo

    do v = 1,  MSTRN_ncfc
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

    do v = 1,  RD_naero
    do j = JS, JE
    do i = IS, IE
    do RD_k = 1, RD_KADD
       aerosol_conc_merge(RD_k,i,j,v) = RD_aerosol_conc(RD_k,v)
       aerosol_radi_merge(RD_k,i,j,v) = RD_aerosol_radi(RD_k,v)
    enddo
    enddo
    enddo
    enddo

    do ihydro = 1, MP_QA
    do j = JS, JE
    do i = IS, IE
    do RD_k = RD_KADD+1, RD_KMAX
       k = KS + RD_KMAX - RD_k ! reverse axis

       aerosol_conc_merge(RD_k,i,j,ihydro) = max( MP_Qe(k,i,j,ihydro), 0.0_RP ) &
                                           / MP_DENS(ihydro) * DENS(k,i,j) / PPM ! [PPM]
       aerosol_radi_merge(RD_k,i,j,ihydro) = MP_Re(k,i,j,ihydro)
    enddo
    enddo
    enddo
    enddo

    do iaero = 1, AE_QA
       if ( I_AE2ALL(iaero) > 0 ) then
          do j = JS, JE
          do i = IS, IE
          do RD_k = RD_KADD+1, RD_KMAX
             k = KS + RD_KMAX - RD_k ! reverse axis

             aerosol_conc_merge(RD_k,i,j,MP_QA+iaero) = max( QTRC(k,i,j,I_AE2ALL(iaero)), 0.0_RP ) &
                                                      / AE_DENS(iaero) * DENS(k,i,j) / PPM ! [PPM]
             aerosol_radi_merge(RD_k,i,j,MP_QA+iaero) = AE_Re(k,i,j,iaero)
          enddo
          enddo
          enddo
       else
          do j = JS, JE
          do i = IS, IE
          do RD_k = RD_KADD+1, RD_KMAX
             aerosol_conc_merge(RD_k,i,j,MP_QA+iaero) = RD_aerosol_conc(RD_k,MP_QA+iaero)
             aerosol_radi_merge(RD_k,i,j,MP_QA+iaero) = RD_aerosol_radi(RD_k,MP_QA+iaero)
          enddo
          enddo
          enddo
       endif
    enddo

    ! calc radiative transfer
    call RD_MSTRN_DTRN3( RD_KMAX,                      & ! [IN]
                         IA, IS, IE,                   & ! [IN]
                         JA, JS, JE,                   & ! [IN]
                         MSTRN_ngas,                   & ! [IN]
                         MSTRN_ncfc,                   & ! [IN]
                         RD_naero,                     & ! [IN]
                         RD_hydro_str,                 & ! [IN]
                         RD_hydro_end,                 & ! [IN]
                         RD_aero_str,                  & ! [IN]
                         RD_aero_end,                  & ! [IN]
                         solins            (:,:),      & ! [IN]
                         cosSZA            (:,:),      & ! [IN]
                         rhodz_merge       (:,:,:),    & ! [IN]
                         pres_merge        (:,:,:),    & ! [IN]
                         temp_merge        (:,:,:),    & ! [IN]
                         temph_merge       (:,:,:),    & ! [IN]
                         temp_sfc          (:,:),      & ! [IN]
                         gas_merge         (:,:,:,:),  & ! [IN]
                         cfc_merge         (:,:,:,:),  & ! [IN]
                         aerosol_conc_merge(:,:,:,:),  & ! [IN]
                         aerosol_radi_merge(:,:,:,:),  & ! [IN]
                         I_MPAE2RD         (:),        & ! [IN]
                         cldfrac_merge     (:,:,:),    & ! [IN]
                         albedo_land       (:,:,:),    & ! [IN]
                         oceanfrc          (:,:),      & ! [IN]
                         flux_rad_merge    (:,:,:,:,:) ) ! [OUT]

    ! return to grid coordinate of LES domain
    do j = JS, JE
    do i = IS, IE
    do RD_k = RD_KADD+1, RD_KMAX+1
       k = KS + RD_KMAX - RD_k ! reverse axis

       flux_rad(k,i,j,I_LW,I_up) = flux_rad_merge(RD_k,i,j,I_LW,I_up)
       flux_rad(k,i,j,I_LW,I_dn) = flux_rad_merge(RD_k,i,j,I_LW,I_dn)
       flux_rad(k,i,j,I_SW,I_up) = flux_rad_merge(RD_k,i,j,I_SW,I_up)
       flux_rad(k,i,j,I_SW,I_dn) = flux_rad_merge(RD_k,i,j,I_SW,I_dn)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       flux_rad_top(i,j,I_LW) = flux_rad_merge(1,i,j,I_LW,I_up)-flux_rad_merge(1,i,j,I_LW,I_dn)
       flux_rad_top(i,j,I_SW) = flux_rad_merge(1,i,j,I_SW,I_up)-flux_rad_merge(1,i,j,I_SW,I_dn)
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
       imax, IS, IE, &
       jmax, JS, JE, &
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
       oceanfrc,     &
       rflux         )
    use scale_const, only: &
       GRAV => CONST_GRAV, &
       Pstd => CONST_Pstd, &
       PPM  => CONST_PPM
    implicit none

    integer,  intent(in)  :: kmax
    integer,  intent(in)  :: imax, IS, IE
    integer,  intent(in)  :: jmax, JS, JE
    integer,  intent(in)  :: ngas
    integer,  intent(in)  :: ncfc
    integer,  intent(in)  :: naero
    integer,  intent(in)  :: hydro_str
    integer,  intent(in)  :: hydro_end
    integer,  intent(in)  :: aero_str
    integer,  intent(in)  :: aero_end
    real(RP), intent(in)  :: solins      (imax,jmax)
    real(RP), intent(in)  :: cosSZA      (imax,jmax)
    real(RP), intent(in)  :: rhodz       (kmax  ,imax,jmax)
    real(RP), intent(in)  :: pres        (kmax  ,imax,jmax)
    real(RP), intent(in)  :: temp        (kmax  ,imax,jmax)
    real(RP), intent(in)  :: temph       (kmax+1,imax,jmax)
    real(RP), intent(in)  :: temp_sfc    (imax,jmax)
    real(RP), intent(in)  :: gas         (kmax,imax,jmax,ngas )
    real(RP), intent(in)  :: cfc         (kmax,imax,jmax,ncfc )
    real(RP), intent(in)  :: aerosol_conc(kmax,imax,jmax,naero)
    real(RP), intent(in)  :: aerosol_radi(kmax,imax,jmax,naero)
    integer,  intent(in)  :: aero2ptype  (naero)
    real(RP), intent(in)  :: cldfrac     (kmax,imax,jmax)
    real(RP), intent(in)  :: albedo_land (imax,jmax,2)
    real(RP), intent(in)  :: oceanfrc    (imax,jmax)
    real(RP), intent(out) :: rflux       (kmax+1,imax,jmax,2,2)

    ! for P-T fitting
    real(RP) :: dz_std (kmax,imax,jmax)       ! layer thickness at 0C, 1atm [cm]
    real(RP) :: logP   (kmax,imax,jmax)       ! log10(pres)
    real(RP) :: logT   (kmax,imax,jmax)       ! log10(temp)
    integer  :: indexP (kmax,imax,jmax)       ! index for interpolation in P-fitting
    real(RP) :: factP  (kmax,imax,jmax)       ! interpolation factor    in P-fitting
    real(RP) :: factT32(kmax,imax,jmax)       ! interpolation factor    in T-fitting
    real(RP) :: factT21(kmax,imax,jmax)       ! interpolation factor    in T-fitting
    integer  :: indexR (kmax,imax,jmax,naero) ! index for interpolation in R-fitting
    real(RP) :: factR  (kmax,imax,jmax,naero) ! interpolation factor    in R-fitting

    ! for optical thickness by gas
    real(RP) :: tauGAS(kmax,imax,jmax,MSTRN_ch_limit) ! optical thickness by gas absorption (total)
    real(RP) :: A1, A2, A3, factPT
    real(RP) :: qv
    integer  :: gasno

    ! for optical thickness by particles
    real(RP) :: tauPR   (kmax,imax,jmax,MSTRN_ncloud)               ! optical thickness        by Rayleigh/cloud/aerosol
    real(RP) :: omgPR   (kmax,imax,jmax,MSTRN_ncloud)               ! single scattering albedo by Rayleigh/cloud/aerosol
    real(RP) :: optparam(kmax,imax,jmax,MSTRN_nmoment,MSTRN_ncloud) ! optical parameters
    real(RP) :: q_fit, dp_P

    ! for albedo
    real(RP) :: albedo_sfc  (imax,jmax,MSTRN_ncloud) ! surface albedo
    real(RP) :: albedo_ocean(imax,jmax,2)            ! surface albedo
    real(RP) :: tau_column  (imax,jmax)

    ! for planck functions
    real(RP) :: bbar (kmax  ,imax,jmax) ! planck functions for thermal source at the interface
    real(RP) :: bbarh(kmax+1,imax,jmax) ! planck functions for thermal source at the center
    real(RP) :: b_sfc(imax,jmax)        ! planck functions for thermal source at the surface
    real(RP) :: wl, beta

    ! for two-stream
    real(RP) :: tau(kmax,imax,jmax,    MSTRN_ncloud) ! total optical thickness
    real(RP) :: omg(kmax,imax,jmax,    MSTRN_ncloud) ! single scattering albedo
    real(RP) :: g  (kmax,imax,jmax,0:2,MSTRN_ncloud) ! two-stream approximation factors
                                                     ! 0: always 1
                                                     ! 1: asymmetry factor
                                                     ! 2: truncation factor
    real(RP) :: b  (kmax,imax,jmax,0:2,MSTRN_ncloud) ! planck expansion coefficients (zero if SW)
    real(RP) :: fsol_rgn(imax,jmax)                  ! solar insolation              (zero if LW)

    real(RP) :: flux       (kmax+1,imax,jmax,2)      ! upward/downward flux
    real(RP) :: flux_direct(kmax+1,imax,jmax)        ! downward flux (direct solar)

    real(RP) :: zerosw
    integer  :: ip, ir, irgn
    integer  :: igas, icfc, iaero, iptype
    integer  :: iw, ich, iplk, icloud, im
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = 1, kmax
       dz_std(k,i,j) = rhodz(k,i,j) * RRHO_std ! [cm]
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 1, kmax
       logP(k,i,j) = log10( pres(k,i,j) )
       logT(k,i,j) = log10( temp(k,i,j) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 1, kmax
       indexP(k,i,j) = MSTRN_nfitP
       do ip = MSTRN_nfitP, 2, -1
          if( logP(k,i,j) >= logfitP(ip) ) indexP(k,i,j) = ip
       enddo
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 1, kmax
       ip = indexP(k,i,j)

       factP(k,i,j) = ( logP(k,i,j) - logfitP(ip-1) ) / ( logfitP(ip) - logfitP(ip-1) )

       factT32(k,i,j) = ( logT(k,i,j) - logfitT(2)  ) / ( logfitT(3) - logfitT(2) ) &
                      * ( temp(k,i,j) - fitT(1)     ) / ( fitT(3)    - fitT(1)    )
       factT21(k,i,j) = ( logT(k,i,j) - logfitT(2)  ) / ( logfitT(2) - logfitT(1) ) &
                      * ( fitT(3)     - temp(k,i,j) ) / ( fitT(3)    - fitT(1)    )
    enddo
    enddo
    enddo

    !---< interpolation of mode radius & hygroscopic parameter (R-fitting) >---
    do iaero = 1, naero
       iptype = aero2ptype(iaero)

       do j = JS, JE
       do i = IS, IE
       do k = 1, kmax
          if ( aerosol_radi(k,i,j,iaero) <= radmode(iptype,1) ) then ! extrapolation

             ir = 1
             indexR(k,i,j,iaero) = ir
             factR (k,i,j,iaero) = ( aerosol_radi(k,i,j,iaero) - radmode(iptype,ir) ) &
                                 / ( radmode(iptype,ir+1)      - radmode(iptype,ir) )

          elseif( aerosol_radi(k,i,j,iaero) > radmode(iptype,MSTRN_nradius) ) then ! extrapolation

             ir = MSTRN_nradius - 1
             indexR(k,i,j,iaero) = ir
             factR (k,i,j,iaero) = ( radmode(iptype,ir+1) - aerosol_radi(k,i,j,iaero) ) &
                                 / ( radmode(iptype,ir+1) - radmode(iptype,ir)        ) &
                                 + 1.0_RP

          else
             do ir = 1, MSTRN_nradius-1
                if (       aerosol_radi(k,i,j,iaero) <= radmode(iptype,ir+1) &
                     .AND. aerosol_radi(k,i,j,iaero) >  radmode(iptype,ir  ) ) then ! interpolation

                   indexR(k,i,j,iaero) = ir
                   factR (k,i,j,iaero) = ( aerosol_radi(k,i,j,iaero) - radmode(iptype,ir) ) &
                                       / ( radmode(iptype,ir+1)      - radmode(iptype,ir) )

                endif
             enddo
          endif
       enddo
       enddo
       enddo
    enddo

    ! initialize
    rflux(:,:,:,:,:) = 0.0_RP

    do iw = 1, MSTRN_nband
       irgn = iflgb(I_SWLW,iw) + 1

       !---< interpolation of gas parameters (P-T fitting) >---
       do ich = 1, nch(iw)
       do j = JS, JE
       do i = IS, IE
       do k = 1, kmax
          tauGAS(k,i,j,ich) = 0.D0
       enddo
       enddo
       enddo
       enddo

       !--- Gas line absorption
       do igas = 1, ngasabs(iw)
          gasno = igasabs(igas,iw)

          do ich = 1, nch(iw)
          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax
             ip = indexP(k,i,j)

             A1 = AKD(ich,ip-1,1,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                + AKD(ich,ip  ,1,gasno,iw) * (          factP(k,i,j) )
             A2 = AKD(ich,ip-1,2,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                + AKD(ich,ip  ,2,gasno,iw) * (          factP(k,i,j) )
             A3 = AKD(ich,ip-1,3,gasno,iw) * ( 1.0_RP - factP(k,i,j) )&
                + AKD(ich,ip  ,3,gasno,iw) * (          factP(k,i,j) )

             factPT = factT32(k,i,j)*(A3-A2) + A2 + factT21(k,i,j)*(A2-A1)

             tauGAS(k,i,j,ich) = tauGAS(k,i,j,ich) &
                               + 10.0_RP**factPT * gas(k,i,j,igasabs(igas,iw)) * PPM * dz_std(k,i,j)
          enddo
          enddo
          enddo
          enddo ! channel loop
       enddo ! gas loop

       !--- Gas broad absorption
       if ( iflgb(I_H2O_continuum,iw) == 1 ) then
          do ich = 1, nch(iw)
          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax
             qv = gas(k,i,j,1) * PPM * dz_std(k,i,j)

             A1 = SKD(ich,ip-1,1,iw) * ( 1.0_RP-factP(k,i,j) )&
                + SKD(ich,ip  ,1,iw) * (        factP(k,i,j) )
             A2 = SKD(ich,ip-1,2,iw) * ( 1.0_RP-factP(k,i,j) )&
                + SKD(ich,ip  ,2,iw) * (        factP(k,i,j) )
             A3 = SKD(ich,ip-1,3,iw) * ( 1.0_RP-factP(k,i,j) )&
                + SKD(ich,ip  ,3,iw) * (        factP(k,i,j) )

             factPT = factT32(k,i,j)*(A3-A2) + A2 + factT21(k,i,j)*(A2-A1)

             tauGAS(k,i,j,ich) = tauGAS(k,i,j,ich) &
                               + 10.0_RP**factPT * qv*qv / ( qv + dz_std(k,i,j) )
          enddo
          enddo
          enddo
          enddo ! channel loop
       endif

       if ( iflgb(I_CFC_continuum,iw) == 1 ) then
          do icfc = 1, ncfc
          do ich = 1, nch(iw)
          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax
             tauGAS(k,i,j,ich) = tauGAS(k,i,j,ich) &
                               + 10.0_RP**acfc(icfc,iw) * cfc(k,i,j,icfc) * dz_std(k,i,j)
          enddo
          enddo
          enddo
          enddo
          enddo
       endif

       !---< particle (Rayleigh/Cloud/Aerosol) >---

       ! optical thickness, phase function
       ! im=1: extinction coefficient
       ! im=2,3,4: moments of the volume scattering phase function

       !--- Rayleigh scattering
       do im = 1, MSTRN_nstream*2+2
       do j = JS, JE
       do i = IS, IE
       do k = 1, kmax
          dp_P = rhodz(k,i,j) * GRAV / Pstd

          optparam(k,i,j,im,I_Cloud   ) = rayleigh(iw) * qmol(im,iw) * dp_P
          optparam(k,i,j,im,I_ClearSky) = rayleigh(iw) * qmol(im,iw) * dp_P
       enddo
       enddo
       enddo
       enddo

       !--- Cloud scattering
       do iaero = hydro_str, hydro_end
          iptype = aero2ptype(iaero)

          do im = 1, MSTRN_nstream*2+2
          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax
             ir = indexR(k,i,j,iaero)

             q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(k,i,j,iaero) ) &
                   + q(ir+1,iptype,im,iw) * (        factR(k,i,j,iaero) )

             optparam(k,i,j,im,I_Cloud   ) = optparam(k,i,j,im,I_Cloud   ) &
                                           + q_fit * aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)
          enddo
          enddo
          enddo
          enddo
       enddo

       !--- Aerosol scattering
       do iaero = aero_str, aero_end
          iptype = aero2ptype(iaero)

          do im = 1, MSTRN_nstream*2+2
          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax
             ir = indexR(k,i,j,iaero)

             q_fit = q(ir  ,iptype,im,iw) * ( 1.0_RP-factR(k,i,j,iaero) ) &
                   + q(ir+1,iptype,im,iw) * (        factR(k,i,j,iaero) )

             optparam(k,i,j,im,I_Cloud   ) = optparam(k,i,j,im,I_Cloud   ) &
                                           + q_fit * aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)
             optparam(k,i,j,im,I_ClearSky) = optparam(k,i,j,im,I_ClearSky) &
                                           + q_fit * aerosol_conc(k,i,j,iaero) * PPM * dz_std(k,i,j)
          enddo
          enddo
          enddo
          enddo
       enddo

       do icloud = 1, MSTRN_ncloud
       do j = JS, JE
       do i = IS, IE
       do k = 1, kmax
          tauPR(k,i,j,icloud) = optparam(k,i,j,1,icloud)
          omgPR(k,i,j,icloud) = optparam(k,i,j,1,icloud) - optparam(k,i,j,2,icloud)

          !--- g
          zerosw = 0.5_RP - sign(0.5_RP,omgPR(k,i,j,icloud)-RD_EPS)

          g(k,i,j,0,icloud) = 1.0_RP
          g(k,i,j,1,icloud) = optparam(k,i,j,3,icloud) * ( 1.0_RP-zerosw ) / ( omgPR(k,i,j,icloud)+zerosw )
          g(k,i,j,2,icloud) = optparam(k,i,j,4,icloud) * ( 1.0_RP-zerosw ) / ( omgPR(k,i,j,icloud)+zerosw )
       enddo
       enddo
       enddo
       enddo

       !--- Albedo
       ! [NOTE] mstrn has look-up table for albedo.
       !        Original scheme calculates albedo by using land-use index (and surface wetness).
       !        In the atmospheric model, albedo is calculated by surface model.
       do icloud = 1, MSTRN_ncloud
          do j = JS, JE
          do i = IS, IE
             tau_column(i,j) = 0.0_RP
             do k = 1, kmax
                tau_column(i,j) = tau_column(i,j) + tauPR(k,i,j,icloud) ! layer-total(for ocean albedo)
             enddo
          enddo
          enddo

          call RD_albedo_ocean( imax, IS, IE,       & ! [IN]
                                jmax, JS, JE,       & ! [IN]
                                cosSZA      (:,:),  & ! [IN]
                                tau_column  (:,:),  & ! [IN]
                                albedo_ocean(:,:,:) ) ! [OUT]

          do j = JS, JE
          do i = IS, IE
             albedo_sfc(i,j,icloud) = (        oceanfrc(i,j) ) * albedo_ocean(i,j,irgn) &
                                    + ( 1.0_RP-oceanfrc(i,j) ) * albedo_land (i,j,irgn)
          enddo
          enddo
       enddo

       ! sub-channel loop
       do ich = 1, nch(iw)

          !--- total tau
          do icloud = 1, 2
          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax
             tau(k,i,j,icloud) = tauGAS(k,i,j,ich) + tauPR(k,i,j,icloud)
          enddo
          enddo
          enddo
          enddo

          !--- omega
          do icloud = 1, 2
          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax
             zerosw = 0.5_RP - sign( 0.5_RP, tau(k,i,j,icloud)-RD_EPS ) ! if tau < EPS, zerosw = 1

             omg(k,i,j,icloud) = ( 1.0_RP-zerosw ) * omgPR(k,i,j,icloud) / ( tau(k,i,j,icloud)-zerosw ) &
                               + (        zerosw ) * 1.0_RP
          enddo
          enddo
          enddo
          enddo

          !--- bn
          if ( irgn == I_SW ) then ! solar

             do icloud = 1, 2
             do j = JS, JE
             do i = IS, IE
             do k = 1, kmax
                b(k,i,j,0,icloud) = 0.0_RP
                b(k,i,j,1,icloud) = 0.0_RP
                b(k,i,j,2,icloud) = 0.0_RP
             enddo
             enddo
             enddo
             enddo

             do j = JS, JE
             do i = IS, IE
                b_sfc(i,j) = 0.0_RP
                fsol_rgn(i,j) = fsol(iw) / fsol_tot * solins(i,j)
             enddo
             enddo

          elseif( irgn == I_LW ) then ! IR
             !--- set planck functions
             wl = 10000.0_RP / sqrt( waveh(iw) * waveh(iw+1) )

             ! from temp at cell center
             do j = JS, JE
             do i = IS, IE
             do k = 1, kmax
                beta = 0.0_RP
                do iplk = MSTRN_nfitPLK, 1, -1
                   beta = beta / ( wl*temp(k,i,j) ) + fitPLK(iplk,iw)
                enddo

                bbar(k,i,j) = exp(-beta) * temp(k,i,j) / (wl*wl)
             enddo
             enddo
             enddo

             ! from temp at cell wall
             do j = JS, JE
             do i = IS, IE
             do k = 1, kmax+1
                beta = 0.0_RP
                do iplk = MSTRN_nfitPLK, 1, -1
                   beta = beta / ( wl*temph(k,i,j) ) + fitPLK(iplk,iw)
                enddo

                bbarh(k,i,j) = exp(-beta) * temph(k,i,j) / (wl*wl)
             enddo
             enddo
             enddo

             do icloud = 1, 2
             do j = JS, JE
             do i = IS, IE
             do k = 1, kmax
                zerosw = 0.5_RP - sign( 0.5_RP, tau(k,i,j,icloud)-RD_EPS ) ! if tau < EPS, zerosw = 1

                b(k,i,j,0,icloud) = bbarh(k,i,j)
                b(k,i,j,1,icloud) = ( 1.0_RP-zerosw ) &
                                  * ( -          bbarh(k+1,i,j) &
                                      + 4.0_RP * bbar (k  ,i,j) &
                                      - 3.0_RP * bbarh(k  ,i,j) &
                                    ) / ( tau(k,i,j,icloud)-zerosw )
                b(k,i,j,2,icloud) = ( 1.0_RP-zerosw ) &
                                  * ( +          bbarh(k+1,i,j) &
                                      - 2.0_RP * bbar (k  ,i,j) &
                                      +          bbarh(k  ,i,j) &
                                    ) / ( tau(k,i,j,icloud)*tau(k,i,j,icloud)-zerosw ) * 2.0_RP
             enddo
             enddo
             enddo
             enddo

             ! from temp_sfc
             do j = JS, JE
             do i = IS, IE
                beta = 0.0_RP
                do iplk = MSTRN_nfitPLK, 1, -1
                   beta = beta / ( wl*temp_sfc(i,j) ) + fitPLK(iplk,iw)
                enddo

                b_sfc(i,j) = exp(-beta) * temp_sfc(i,j) / (wl*wl)
             enddo
             enddo

             do j = JS, JE
             do i = IS, IE
                fsol_rgn(i,j) = 0.0_RP
             enddo
             enddo

          endif ! solar/IR switch

          !if( IO_L ) write(IO_FID_LOG,*) "tau sum", iw, ich, sum   (tau(:,IS:IE,JS:JE,1)), sum   (tau(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "tau max", iw, ich, maxval(tau(:,IS:IE,JS:JE,1)), maxval(tau(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "tau min", iw, ich, minval(tau(:,IS:IE,JS:JE,1)), minval(tau(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "omg sum", iw, ich, sum   (omg(:,IS:IE,JS:JE,1)), sum   (omg(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "omg max", iw, ich, maxval(omg(:,IS:IE,JS:JE,1)), maxval(omg(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "omg min", iw, ich, minval(omg(:,IS:IE,JS:JE,1)), minval(omg(:,IS:IE,JS:JE,2))

          ! two-stream transfer
          call RD_MSTRN_two_stream( kmax,                   & ! [IN]
                                    imax, IS, IE,           & ! [IN]
                                    jmax, JS, JE,           & ! [IN]
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
                                    flux       (:,:,:,:),   & ! [OUT]
                                    flux_direct(:,:,:)      ) ! [OUT]

          do j = JS, JE
          do i = IS, IE
          do k = 1, kmax+1
             rflux(k,i,j,irgn,I_up) = rflux(k,i,j,irgn,I_up) + flux(k,i,j,I_up) * wgtch(ich,iw)
             rflux(k,i,j,irgn,I_dn) = rflux(k,i,j,irgn,I_dn) + flux(k,i,j,I_dn) * wgtch(ich,iw)
          enddo
          enddo
          enddo

          !if( IO_L ) write(IO_FID_LOG,*) "flux sum", iw, ich, sum   (flux(:,IS:IE,JS:JE,1)), sum   (flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux max", iw, ich, maxval(flux(:,IS:IE,JS:JE,1)), maxval(flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux min", iw, ich, minval(flux(:,IS:IE,JS:JE,1)), minval(flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux mal", iw, ich, maxloc(flux(:,IS:IE,JS:JE,1)), maxloc(flux(:,IS:IE,JS:JE,2))
          !if( IO_L ) write(IO_FID_LOG,*) "flux mil", iw, ich, minloc(flux(:,IS:IE,JS:JE,1)), minloc(flux(:,IS:IE,JS:JE,2))

       enddo ! ICH loop
    enddo ! IW loop

    return
  end subroutine RD_MSTRN_DTRN3

  !-----------------------------------------------------------------------------
  !> Two stream calculation CORE
  subroutine RD_MSTRN_two_stream( &
       kmax,         &
       imax, IS, IE, &
       jmax, JS, JE, &
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
       PI   => CONST_PI,  &
       EPS1 => CONST_EPS1
    implicit none

    integer,  intent(in)  :: kmax
    integer,  intent(in)  :: imax, IS, IE
    integer,  intent(in)  :: jmax, JS, JE
    integer,  intent(in)  :: iw, ich
    real(RP), intent(in)  :: cosSZA0    (imax,jmax)                       ! cos(SZA) = mu0
    real(RP), intent(in)  :: fsol       (imax,jmax)                       ! solar radiation intensity
    integer,  intent(in)  :: irgn                                         ! 1:SW 2:LW
    real(RP), intent(in)  :: tau        (kmax,imax,jmax,    MSTRN_ncloud) ! total optical thickness          (clear-sky/cloud)
    real(RP), intent(in)  :: omg        (kmax,imax,jmax,    MSTRN_ncloud) ! single scattering albedo         (clear-sky/cloud)
    real(RP), intent(in)  :: g          (kmax,imax,jmax,0:2,MSTRN_ncloud) ! two-stream approximation factors (clear-sky/cloud)
    real(RP), intent(in)  :: b          (kmax,imax,jmax,0:2,MSTRN_ncloud) ! planck expansion coefficients    (clear-sky/cloud)
    real(RP), intent(in)  :: b_sfc      (imax,jmax)                       ! planck function at surface
    real(RP), intent(in)  :: albedo_sfc (imax,jmax,MSTRN_ncloud)          ! surface albedo                   (clear-sky/cloud)
    real(RP), intent(in)  :: cldfrac    (kmax,imax,jmax)                  ! cloud fraction

    real(RP), intent(out) :: flux       (kmax+1,imax,jmax,2)              ! upward(sfc->TOA)/downward(TOA->sfc) flux
    real(RP), intent(out) :: flux_direct(kmax+1,imax,jmax)                ! downward(TOA->sfc) flux, solar direct

    ! parameters with two-stream truncation
    real(RP) :: tau_new    ! optical thickness        : two-stream truncation
    real(RP) :: omg_new    ! single scattering albedo : two-stream truncation
    real(RP) :: g_new      ! asymmetric factor        : two-stream truncation
    real(RP) :: b_new(0:2) ! planck function          : two-stream truncation
    real(RP) :: c(0:2)
    real(RP) :: Pmns, Ppls ! Phase  function          : two-stream truncation
    real(RP) :: Smns, Spls ! Source function          : two-stream truncation

    ! working
    real(RP) :: cosSZA(imax,jmax)
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
    real(RP) :: Tdir0(kmax,imax,jmax,MSTRN_ncloud) ! transmission factor for solar direct (clear-sky/cloud)
    real(RP) :: R0   (kmax,imax,jmax,MSTRN_ncloud) ! reflection   factor                  (clear-sky/cloud)
    real(RP) :: T0   (kmax,imax,jmax,MSTRN_ncloud) ! transmission factor                  (clear-sky/cloud)
    real(RP) :: Em_LW(kmax,imax,jmax,MSTRN_ncloud) ! thermal source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Em_SW(kmax,imax,jmax,MSTRN_ncloud) ! solar   source (sfc->TOA)            (clear-sky/cloud)
    real(RP) :: Ep_LW(kmax,imax,jmax,MSTRN_ncloud) ! thermal source (TOA->sfc)            (clear-sky/cloud)
    real(RP) :: Ep_SW(kmax,imax,jmax,MSTRN_ncloud) ! solar   source (TOA->sfc)            (clear-sky/cloud)

    ! Averaged factors, considering cloud overwrap
    real(RP) :: tau_bar(imax,jmax)   ! accumulated optical thickness at each layer
    real(RP) :: R (kmax+1,imax,jmax) ! reflection   factor
    real(RP) :: T (kmax+1,imax,jmax) ! transmission factor
    real(RP) :: Em(kmax+1,imax,jmax) ! source (sfc->TOA)
    real(RP) :: Ep(kmax+1,imax,jmax) ! source (TOA->sfc)

    ! Doubling-Adding
    real(RP) :: R12mns(kmax+1,imax,jmax), R12pls(kmax+1,imax,jmax) ! reflection factor in doubling method
    real(RP) :: E12mns(kmax+1,imax,jmax), E12pls(kmax+1,imax,jmax) ! source function   in doubling method
    real(RP) :: Umns, Upls                                         ! flux intensity

    real(RP) :: factor
    real(RP) :: Em0(imax,jmax,MSTRN_ncloud)

    real(RP) :: sw
    integer  :: k, i, j, icloud
    !---------------------------------------------------------------------------

    cosSZA(:,:) = max( cosSZA0(:,:), RD_cosSZA_min )

!OCL SERIAL
    do icloud = 1, 2
!OCL PARALLEL
       do j = JS, JE
       do i = IS, IE
       do k = 1, kmax

          !---< two-stream truncation >---
          tau_new = ( 1.0_RP - omg(k,i,j,icloud)*g(k,i,j,2,icloud) ) * tau(k,i,j,icloud)

          omg_new = ( 1.0_RP - g(k,i,j,2,icloud) ) / ( 1.0_RP - omg(k,i,j,icloud)*g(k,i,j,2,icloud) ) * omg(k,i,j,icloud)
          omg_new = min( omg_new, EPS1 )

          g_new   = ( g(k,i,j,1,icloud) - g(k,i,j,2,icloud) ) / ( 1.0_RP - g(k,i,j,2,icloud) )

          Tdir0(k,i,j,icloud) = exp(-tau_new/cosSZA(i,j))

          factor   = ( 1.0_RP - omg(k,i,j,icloud)*g(k,i,j,2,icloud) )
          b_new(0) = b(k,i,j,0,icloud)
          b_new(1) = b(k,i,j,1,icloud) / factor
          b_new(2) = b(k,i,j,2,icloud) / (factor*factor)
          c(:)     = Wmns(irgn) * 2.0_RP * PI * ( 1.0_RP - omg_new ) * b_new(:)

          !--- P+, P-
          Pmns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M(irgn)*M(irgn) )
          Ppls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M(irgn)*M(irgn) )

          !--- S+, S-
          Smns = omg_new * 0.5_RP * ( 1.0_RP - 3.0_RP * g_new * M(irgn)*cosSZA(i,j) )
          Spls = omg_new * 0.5_RP * ( 1.0_RP + 3.0_RP * g_new * M(irgn)*cosSZA(i,j) )

          !---< calculate R, T, e+, e- >---
          sw = 0.5_RP + sign(0.5_RP,tau_new-RD_EPS)
          tau_new = max( tau_new, RD_EPS )

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
          R0(k,i,j,icloud) = (        sw ) * 0.5_RP * ( Apls_mns + Bpls_mns ) &
                           + ( 1.0_RP-sw ) * (          tau_new * (          Pmns ) / M(irgn) )
          T0(k,i,j,icloud) = (        sw ) * 0.5_RP * ( Apls_mns - Bpls_mns ) &
                           + ( 1.0_RP-sw ) * ( 1.0_RP - tau_new * ( 1.0_RP - Ppls ) / M(irgn) )

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

          Em_LW(k,i,j,icloud) = (        sw ) * ( V0mns - R0(k,i,j,icloud) * V0pls - T0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * 0.5_RP * tau_new * ( 2.0_RP*c(0) + c(1)*tau_new + c(2)*tau_new*tau_new )
          Ep_LW(k,i,j,icloud) = (        sw ) * ( V1pls - T0(k,i,j,icloud) * V0pls - R0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * 0.5_RP * tau_new * ( 2.0_RP*c(0) + c(1)*tau_new + c(2)*tau_new*tau_new )

          !--- solar source
          SIGmns = Wmns(irgn) * ( Spls - Smns )
          SIGpls = Wmns(irgn) * ( Spls + Smns )

          Qgamma = ( SIGpls*X*cosSZA(i,j) + SIGmns ) / ( X*Y*cosSZA(i,j) - 1.0/cosSZA(i,j) )

          V0pls = 0.5_RP * ( ( 1.0_RP + 1.0_RP/(X*cosSZA(i,j)) ) * Qgamma + SIGmns / X )
          V0mns = 0.5_RP * ( ( 1.0_RP - 1.0_RP/(X*cosSZA(i,j)) ) * Qgamma - SIGmns / X )

          V1pls = V0pls * Tdir0(k,i,j,icloud)
          V1mns = V0mns * Tdir0(k,i,j,icloud)

          Em_SW(k,i,j,icloud) = (        sw ) * ( V0mns - R0(k,i,j,icloud) * V0pls - T0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * Wmns(irgn) * Smns * tau_new * sqrt( Tdir0(k,i,j,icloud) )
          Ep_SW(k,i,j,icloud) = (        sw ) * ( V1pls - T0(k,i,j,icloud) * V0pls - R0(k,i,j,icloud) * V1mns ) &
                              + ( 1.0_RP-sw ) * Wmns(irgn) * Spls * tau_new * sqrt( Tdir0(k,i,j,icloud) )


       enddo ! k loop
       enddo ! i loop
       enddo ! j loop
    enddo ! cloud loop

    !---< consider partial cloud layer: semi-random over-wrapping >---
    do j = JS, JE
    do i = IS, IE
       tau_bar(i,j) = 1.0_RP ! k-recurrence

       do k = 1, kmax

          R (k,i,j) = (        cldfrac(k,i,j) ) * R0(k,i,j,I_Cloud   ) &
                    + ( 1.0_RP-cldfrac(k,i,j) ) * R0(k,i,j,I_ClearSky)

          T (k,i,j) = (        cldfrac(k,i,j) ) * T0(k,i,j,I_Cloud   ) &
                    + ( 1.0_RP-cldfrac(k,i,j) ) * T0(k,i,j,I_ClearSky)

          Em(k,i,j) = (        cldfrac(k,i,j) ) * ( Em_LW(k,i,j,I_Cloud   ) &
                                                  + Em_SW(k,i,j,I_Cloud   ) * tau_bar(i,j) * fsol(i,j) ) &
                    + ( 1.0_RP-cldfrac(k,i,j) ) * ( Em_LW(k,i,j,I_ClearSky) &
                                                  + Em_SW(k,i,j,I_ClearSky) * tau_bar(i,j) * fsol(i,j) )

          Ep(k,i,j) = (        cldfrac(k,i,j) ) * ( Ep_LW(k,i,j,I_Cloud   ) &
                                                  + Ep_SW(k,i,j,I_Cloud   ) * tau_bar(i,j) * fsol(i,j) ) &
                    + ( 1.0_RP-cldfrac(k,i,j) ) * ( Ep_LW(k,i,j,I_ClearSky) &
                                                  + Ep_SW(k,i,j,I_ClearSky) * tau_bar(i,j) * fsol(i,j) )

          flux_direct(k,i,j) = cosSZA(i,j) * tau_bar(i,j) * fsol(i,j)

          ! update tau_bar
          tau_bar(i,j) = tau_bar(i,j) * ( (        cldfrac(k,i,j) ) * Tdir0(k,i,j,I_Cloud   ) &
                                        + ( 1.0_RP-cldfrac(k,i,j) ) * Tdir0(k,i,j,I_ClearSky) )
       enddo ! k loop
    enddo ! i loop
    enddo ! j loop

    do j = JS, JE
    do i = IS, IE
       ! at lambert surface
       R (kmax+1,i,j) = (        cldfrac(kmax,i,j) ) * albedo_sfc(i,j,I_Cloud   ) &
                      + ( 1.0_RP-cldfrac(kmax,i,j) ) * albedo_sfc(i,j,I_ClearSky)

       T (kmax+1,i,j) = 0.0_RP

       flux_direct(kmax+1,i,j) = cosSZA(i,j) * tau_bar(i,j) * fsol(i,j)

       Em0(i,j,I_Cloud   ) = Wpls(irgn) * ( flux_direct(kmax+1,i,j) * albedo_sfc(i,j,I_Cloud   ) / (W(irgn)*M(irgn)) &
                           + 2.0_RP * PI * ( 1.0_RP-albedo_sfc(i,j,I_Cloud   ) ) * b_sfc(i,j) )
       Em0(i,j,I_ClearSky) = Wpls(irgn) * ( flux_direct(kmax+1,i,j) * albedo_sfc(i,j,I_ClearSky) / (W(irgn)*M(irgn)) &
                           + 2.0_RP * PI * ( 1.0_RP-albedo_sfc(i,j,I_ClearSky) ) * b_sfc(i,j) )

       Em(kmax+1,i,j) = (        cldfrac(kmax,i,j) ) * Em0(i,j,I_Cloud   ) &
                      + ( 1.0_RP-cldfrac(kmax,i,j) ) * Em0(i,j,I_ClearSky)

       Ep(kmax+1,i,j) = 0.0_RP
    enddo ! i loop
    enddo ! j loop

    !---< Adding-Doubling method >---
    ! [note] TOA->Surface is positive direction. "pls" means upper to lower altitude.

    ! adding: surface to TOA
    do j = JS, JE
    do i = IS, IE
       R12pls(kmax+1,i,j) = R (kmax+1,i,j)
       E12mns(kmax+1,i,j) = Em(kmax+1,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = kmax, 1, -1
       R12pls(k,i,j) = R (k,i,j) &
                     + T(k,i,j) / ( 1.0_RP - R12pls(k+1,i,j) * R(k,i,j) ) * ( R12pls(k+1,i,j) * T (k,i,j)                   )
       E12mns(k,i,j) = Em(k,i,j) &
                     + T(k,i,j) / ( 1.0_RP - R12pls(k+1,i,j) * R(k,i,j) ) * ( R12pls(k+1,i,j) * Ep(k,i,j) + E12mns(k+1,i,j) )
    enddo
    enddo
    enddo

    ! adding: TOA to surface
    do j = JS, JE
    do i = IS, IE
       R12mns(1,i,j) = R (1,i,j)
       E12pls(1,i,j) = Ep(1,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 2, kmax+1
       R12mns(k,i,j) = R (k,i,j) &
                     + T(k,i,j) / ( 1.0_RP - R12mns(k-1,i,j) * R(k,i,j) ) * ( R12mns(k-1,i,j) *T (k,i,j) )
       E12pls(k,i,j) = Ep(k,i,j) &
                     + T(k,i,j) / ( 1.0_RP - R12mns(k-1,i,j) * R(k,i,j) ) * ( R12mns(k-1,i,j)*Em(k,i,j) + E12pls(k-1,i,j) )
    enddo
    enddo
    enddo

    !--- radiative flux at cell wall:
    ! [note] "d" means upper to lower altitude.

    ! TOA boundary
    do j = JS, JE
    do i = IS, IE
       flux(1,i,j,I_up) = Wscale(irgn) * E12mns(1,i,j)
       flux(1,i,j,I_dn) = flux_direct(1,i,j)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 2, kmax+1
       Upls = ( E12pls(k-1,i,j) + R12mns(k-1,i,j)*E12mns(k,i,j) ) / ( 1.0_RP - R12mns(k-1,i,j)*R12pls(k,i,j) )
       Umns =  E12mns(k,i,j) + R12pls(k,i,j) * Upls

       flux(k,i,j,I_up) = Wscale(irgn) * Umns
       flux(k,i,j,I_dn) = Wscale(irgn) * Upls + flux_direct(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine RD_MSTRN_two_stream

  !-----------------------------------------------------------------------------
  ! Sea surface reflectance by Payne
  subroutine RD_albedo_ocean( &
       imax, IS, IE, &
       jmax, JS, JE, &
       cosSZA,       &
       tau,          &
       albedo_ocean  )
    implicit none

    integer,  intent(in)  :: imax, IS, IE
    integer,  intent(in)  :: jmax, JS, JE
    real(RP), intent(in)  :: cosSZA      (imax,jmax)
    real(RP), intent(in)  :: tau         (imax,jmax)
    real(RP), intent(out) :: albedo_ocean(imax,jmax,2)

    real(RP) :: c(5,3)
    data c / -2.8108_RP   , -1.3651_RP,  2.9210E1_RP, -4.3907E1_RP,  1.8125E1_RP, &
              6.5626E-1_RP, -8.7253_RP, -2.7749E1_RP,  4.9486E1_RP, -1.8345E1_RP, &
             -6.5423E-1_RP,  9.9967_RP,  2.7769_RP  , -1.7620E1_RP,  7.0838_RP    /

    real(RP) :: am1, tr1, s
    real(RP) :: sw

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       albedo_ocean(i,j,I_SW) = 0.05_RP
    enddo
    enddo


    do j = JS, JE
    do i = IS, IE
       am1 = max( min( cosSZA(i,j), 0.961_RP ), 0.0349_RP )

       sw = 0.5_RP + sign(0.5_RP,tau(i,j))

       tr1 = max( min( am1 / ( 4.0_RP * tau(i,j) ), 1.0_RP ), 0.05_RP )

       s = 0.0_RP
       do n = 1, 5
          s = s + c(n,1) * tr1**(n-1)           &
                + c(n,2) * tr1**(n-1) * am1     &
                + c(n,3) * tr1**(n-1) * am1*am1
       enddo

       albedo_ocean(i,j,I_LW) = ( 1.0_RP-sw ) * 0.05_RP &
                              + (        sw ) * exp(s)
    enddo
    enddo

    return
  end subroutine RD_albedo_ocean

end module scale_atmos_phy_rd_mstrnx
