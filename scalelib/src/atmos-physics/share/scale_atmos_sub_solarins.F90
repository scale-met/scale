!-------------------------------------------------------------------------------
!> module SOLARINS
!!
!! @par Description
!!          calculate solar insolation
!!          based on berger(1978), berger et al.(1993)
!!          The original algorithm is valid to +/- 1,000,000 years from AD1950
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-01-29 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_solarins
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_SOLARINS_setup
  public :: ATMOS_SOLARINS_orbit
  public :: ATMOS_SOLARINS_insolation

  interface ATMOS_SOLARINS_insolation
     module procedure ATMOS_SOLARINS_insolation_0D
     module procedure ATMOS_SOLARINS_insolation_2D
  end interface ATMOS_SOLARINS_insolation

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, save :: ATMOS_SOLARINS_constant = 1360.250117_RP ! Solar constant [W/m2]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, save :: obliquity ! obliquity [rad]
  real(RP), private, save :: E         ! eccentricity
  real(RP), private, save :: omega     ! longitude of perigee [rad]
  real(RP), private, save :: lambda_m0 ! longitude at the vernal equinox [rad]

  integer,  private, parameter :: year_ref = 1950              ! reference year [year]
  integer,  private, save      :: ve_date(6)                   ! reference date of vernal equinox
  data ve_date / 1950, 3, 21, 0, 0, 0 /

  real(RP), private, parameter :: obliquity_ref = 23.320556_RP ! initial condition of obliquity (epsilon_star)
  real(RP), private, parameter :: psi_bar       = 50.439273_RP ! parameter for general precession [second of arc]
  real(RP), private, parameter :: zeta          =  3.392506_RP ! parameter for general precession [degree]

  !-----< Parameter tables from Berger(1978b) >-----
  integer,  private, parameter :: nObliq = 47         ! # of terms of the series expansion of epsilon
  real(RP), private, save      :: Obliq_amp  (nObliq) ! amplitude [second of arc]
  real(RP), private, save      :: Obliq_rate (nObliq) ! mean rate [second of arc/year]
  real(RP), private, save      :: Obliq_phase(nObliq) ! phase     [degree of arc]

  data Obliq_amp / &
        -2462.2214466_RP, & !  1
         -857.3232075_RP, & !  2
         -629.3231835_RP, & !  3
         -414.2804924_RP, & !  4
         -311.7632587_RP, & !  5
          308.9408604_RP, & !  6
         -162.5533601_RP, & !  7
         -116.1077911_RP, & !  8
          101.1189923_RP, & !  9
          -67.6856209_RP, & ! 10
           24.9079067_RP, & ! 11
           22.5811241_RP, & ! 12
          -21.1648355_RP, & ! 13
          -15.6549876_RP, & ! 14
           15.3936813_RP, & ! 15
           14.6660938_RP, & ! 16
          -11.7273029_RP, & ! 17
           10.2742696_RP, & ! 18
            6.4914588_RP, & ! 19
            5.8539148_RP, & ! 20
           -5.4872205_RP, & ! 21
           -5.4290191_RP, & ! 22
            5.1609570_RP, & ! 23
            5.0786314_RP, & ! 24
           -4.0735782_RP, & ! 25
            3.7227167_RP, & ! 26
            3.3971932_RP, & ! 27
           -2.8347004_RP, & ! 28
           -2.6550721_RP, & ! 29
           -2.5717867_RP, & ! 30
           -2.4712188_RP, & ! 31
            2.4625410_RP, & ! 32
            2.2464112_RP, & ! 33
           -2.0755511_RP, & ! 34
           -1.9713669_RP, & ! 35
           -1.8813061_RP, & ! 36
           -1.8468785_RP, & ! 37
            1.8186742_RP, & ! 38
            1.7601888_RP, & ! 39
           -1.5428851_RP, & ! 40
            1.4738838_RP, & ! 41
           -1.4593669_RP, & ! 42
            1.4192259_RP, & ! 43
           -1.1818980_RP, & ! 44
            1.1756474_RP, & ! 45
           -1.1316126_RP, & ! 46
            1.0896928_RP  / ! 47

  data Obliq_rate / &
        31.609974_RP, & !  1
        32.620504_RP, & !  2
        24.172203_RP, & !  3
        31.983787_RP, & !  4
        44.828336_RP, & !  5
        30.973257_RP, & !  6
        43.668246_RP, & !  7
        32.246691_RP, & !  8
        30.599444_RP, & !  9
        42.681324_RP, & ! 10
        43.836462_RP, & ! 11!
        47.439436_RP, & ! 12
        63.219948_RP, & ! 13
        64.230478_RP, & ! 14
         1.010530_RP, & ! 15
         7.437771_RP, & ! 16
        55.782177_RP, & ! 17
         0.373813_RP, & ! 18
        13.218362_RP, & ! 19
        62.583231_RP, & ! 20
        63.593761_RP, & ! 21
        76.438310_RP, & ! 22
        45.815258_RP, & ! 23
         8.448301_RP, & ! 24
        56.792707_RP, & ! 25
        49.747842_RP, & ! 26
        12.058272_RP, & ! 27
        75.278220_RP, & ! 28
        65.241008_RP, & ! 29
        64.604291_RP, & ! 30
         1.647247_RP, & ! 31
         7.811584_RP, & ! 32
        12.207832_RP, & ! 33
        63.856665_RP, & ! 34
        56.155990_RP, & ! 35
        77.448840_RP, & ! 36
         6.801054_RP, & ! 37
        62.209418_RP, & ! 38
        20.656133_RP, & ! 39
        48.344406_RP, & ! 40
        55.145460_RP, & ! 41
        69.000539_RP, & ! 42
        11.071350_RP, & ! 43
        74.291298_RP, & ! 44
        11.047742_RP, & ! 45
         0.636717_RP, & ! 46
        12.844549_RP  / ! 47

  data Obliq_phase / &
        251.9025_RP, & !  1
        280.8325_RP, & !  2
        128.3057_RP, & !  3
        292.7252_RP, & !  4
         15.3747_RP, & !  5
        263.7951_RP, & !  6
        308.4258_RP, & !  7
        240.0099_RP, & !  8
        222.9725_RP, & !  9
        268.7809_RP, & ! 10
        316.7998_RP, & ! 11
        319.6024_RP, & ! 12
        143.8050_RP, & ! 13
        172.7351_RP, & ! 14
         28.9300_RP, & ! 15
        123.5968_RP, & ! 16
         20.2082_RP, & ! 17
         40.8226_RP, & ! 18
        123.4722_RP, & ! 19
        155.6977_RP, & ! 20
        184.6277_RP, & ! 21
        267.2772_RP, & ! 22
         55.0196_RP, & ! 23
        152.5268_RP, & ! 24
         49.1382_RP, & ! 25
        204.6609_RP, & ! 26
         56.5233_RP, & ! 27
        200.3284_RP, & ! 28
        201.6651_RP, & ! 29
        213.5577_RP, & ! 30
         17.0374_RP, & ! 31
        164.4194_RP, & ! 32
         94.5422_RP, & ! 33
        131.9124_RP, & ! 34
         61.0309_RP, & ! 35
        296.2073_RP, & ! 36
        135.4894_RP, & ! 37
        114.8750_RP, & ! 38
        247.0691_RP, & ! 39
        256.6114_RP, & ! 40
         32.1008_RP, & ! 41
        143.6804_RP, & ! 42
         16.8784_RP, & ! 43
        160.6835_RP, & ! 44
         27.5932_RP, & ! 45
        348.1074_RP, & ! 46
         82.6496_RP  / ! 47

  integer,  private, parameter :: nEclip = 19         ! # of terms of the series expansion of ecliptic
  real(RP), private, save      :: Eclip_amp  (nEclip) ! amplitude
  real(RP), private, save      :: Eclip_rate (nEclip) ! mean rate [second of arc/year]
  real(RP), private, save      :: Eclip_phase(nEclip) ! phase     [degree of arc]

  data Eclip_amp / &
         0.01860798_RP, & !  1
         0.01627522_RP, & !  2
        -0.01300660_RP, & !  3
         0.00988829_RP, & !  4
        -0.00336700_RP, & !  5
         0.00333077_RP, & !  6
        -0.00235400_RP, & !  7
         0.00140015_RP, & !  8
         0.00100700_RP, & !  9
         0.00085700_RP, & ! 10
         0.00064990_RP, & ! 11
         0.00059900_RP, & ! 12
         0.00037800_RP, & ! 13
        -0.00033700_RP, & ! 14
         0.00027600_RP, & ! 15
         0.00018200_RP, & ! 16
        -0.00017400_RP, & ! 17
        -0.00012400_RP, & ! 18
         0.00001250_RP  / ! 19

  data Eclip_rate / &
         4.2072050_RP, & !  1
         7.3460910_RP, & !  2
        17.8572630_RP, & !  3
        17.2205460_RP, & !  4
        16.8467330_RP, & !  5
         5.1990790_RP, & !  6
        18.2310760_RP, & !  7
        26.2167580_RP, & !  8
         6.3591690_RP, & !  9
        16.2100160_RP, & ! 10
         3.0651810_RP, & ! 11
        16.5838290_RP, & ! 12
        18.4939800_RP, & ! 13
         6.1909530_RP, & ! 14
        18.8677930_RP, & ! 15
        17.4255670_RP, & ! 16
         6.1860010_RP, & ! 17
        18.4174410_RP, & ! 18
         0.6678630_RP  / ! 19

  data Eclip_phase / &
         28.620089_RP, & !  1
        193.788772_RP, & !  2
        308.307024_RP, & !  3
        320.199637_RP, & !  4
        279.376984_RP, & !  5
         87.195000_RP, & !  6
        349.129677_RP, & !  7
        128.443387_RP, & !  8
        154.143880_RP, & !  9
        291.269597_RP, & ! 10
        114.860583_RP, & ! 11
        332.092251_RP, & ! 12
        296.414411_RP, & ! 13
        145.769910_RP, & ! 14
        337.237063_RP, & ! 15
        152.092288_RP, & ! 16
        126.839891_RP, & ! 17
        210.667199_RP, & ! 18
         72.108838_RP  / ! 19

  integer,  private, parameter :: nPrece = 78         ! # of terms of the series expansion of general precession
  real(RP), private, save      :: Prece_amp  (nPrece) ! amplitude
  real(RP), private, save      :: Prece_rate (nPrece) ! mean rate [second of arc/year]
  real(RP), private, save      :: Prece_phase(nPrece) ! phase     [degree of arc]

  data Prece_amp / &
         7391.0225890_RP, & !  1
         2555.1526947_RP, & !  2
         2022.7629188_RP, & !  3
        -1973.6517951_RP, & !  4
         1240.2321818_RP, & !  5
          953.8679112_RP, & !  6
         -931.7537108_RP, & !  7
          872.3795383_RP, & !  8
          606.3544732_RP, & !  9
         -496.0274038_RP, & ! 10
          456.9608039_RP, & ! 11
          346.9462320_RP, & ! 12
         -305.8412902_RP, & ! 13
          249.6173246_RP, & ! 14
         -199.1027200_RP, & ! 15
          191.0560889_RP, & ! 16
         -175.2936572_RP, & ! 17
          165.9068833_RP, & ! 18
          161.1285917_RP, & ! 19
          139.7878093_RP, & ! 20
         -133.5228399_RP, & ! 21
          117.0673811_RP, & ! 22
          104.6907281_RP, & ! 23
           95.3227476_RP, & ! 24
           86.7824524_RP, & ! 25
           86.0857729_RP, & ! 26
           70.5893698_RP, & ! 27
          -69.9719343_RP, & ! 28
          -62.5817473_RP, & ! 29
           61.5450059_RP, & ! 30
          -57.9364011_RP, & ! 31
           57.1899832_RP, & ! 32
          -57.0236109_RP, & ! 33
          -54.2119253_RP, & ! 34
           53.2834147_RP, & ! 35
           52.1223575_RP, & ! 36
          -49.0059908_RP, & ! 37
          -48.3118757_RP, & ! 38
          -45.4191685_RP, & ! 39
          -42.2357920_RP, & ! 40
          -34.7971099_RP, & ! 41
           34.4623613_RP, & ! 42
          -33.8356643_RP, & ! 43
           33.6689362_RP, & ! 44
          -31.2521586_RP, & ! 45
          -30.8798701_RP, & ! 46
           28.4640769_RP, & ! 47
          -27.1960802_RP, & ! 48
           27.0860736_RP, & ! 49
          -26.3437456_RP, & ! 50
           24.7253740_RP, & ! 51
           24.6732126_RP, & ! 52
           24.4272733_RP, & ! 53
           24.0127327_RP, & ! 54
           21.7150294_RP, & ! 55
          -21.5375347_RP, & ! 56
           18.1148363_RP, & ! 57
          -16.9603104_RP, & ! 58
          -16.1765215_RP, & ! 59
           15.5567653_RP, & ! 60
           15.4846529_RP, & ! 61
           15.2150632_RP, & ! 62
           14.5047426_RP, & ! 63
          -14.3873316_RP, & ! 64
           13.1351419_RP, & ! 65
           12.8776311_RP, & ! 66
           11.9867234_RP, & ! 67
           11.9385578_RP, & ! 68
           11.7030822_RP, & ! 69
           11.6018181_RP, & ! 70
          -11.2617293_RP, & ! 71
          -10.4664199_RP, & ! 72
           10.4333970_RP, & ! 73
          -10.2377466_RP, & ! 74
           10.1934446_RP, & ! 75
          -10.1280191_RP, & ! 76
           10.0289441_RP, & ! 77
          -10.0034259_RP  / ! 78

  data Prece_rate / &
        31.609974_RP, & !  1
        32.620504_RP, & !  2
        24.172203_RP, & !  3
         0.636717_RP, & !  4
        31.983787_RP, & !  5
         3.138886_RP, & !  6
        30.973257_RP, & !  7
        44.828336_RP, & !  8
         0.991874_RP, & !  9
         0.373813_RP, & ! 10
        43.668246_RP, & ! 11
        32.246691_RP, & ! 12
        30.599444_RP, & ! 13
         2.147012_RP, & ! 14
        10.511172_RP, & ! 15
        42.681324_RP, & ! 16
        13.650058_RP, & ! 17
         0.986922_RP, & ! 18
         9.874455_RP, & ! 19
        13.013341_RP, & ! 20
         0.262904_RP, & ! 21
         0.004952_RP, & ! 22
         1.142024_RP, & ! 23
        63.219948_RP, & ! 24
         0.205021_RP, & ! 25
         2.151964_RP, & ! 26
        64.230478_RP, & ! 27
        43.836462_RP, & ! 28
        47.439436_RP, & ! 29
         1.384343_RP, & ! 30
         7.437771_RP, & ! 31
        18.829299_RP, & ! 32
         9.500642_RP, & ! 33
         0.431696_RP, & ! 34
         1.160090_RP, & ! 35
        55.782177_RP, & ! 36
        12.639528_RP, & ! 37
         1.155138_RP, & ! 38
         0.168216_RP, & ! 39
         1.647247_RP, & ! 40
        10.884985_RP, & ! 41
         5.610937_RP, & ! 42
        12.658184_RP, & ! 43
         1.010530_RP, & ! 44
         1.983748_RP, & ! 45
        14.023871_RP, & ! 46
         0.560178_RP, & ! 47
         1.273434_RP, & ! 48
        12.021467_RP, & ! 49
        62.583231_RP, & ! 50
        63.593761_RP, & ! 51
        76.438310_RP, & ! 52
         4.280910_RP, & ! 53
        13.218362_RP, & ! 54
        17.818769_RP, & ! 55
         8.359495_RP, & ! 56
        56.792707_RP, & ! 57
         8.448301_RP, & ! 58
         1.978796_RP, & ! 59
         8.863925_RP, & ! 60
         0.186365_RP, & ! 61
         8.996212_RP, & ! 62
         6.771027_RP, & ! 63
        45.815258_RP, & ! 64
        12.002811_RP, & ! 65
        75.278220_RP, & ! 66
        65.241008_RP, & ! 67
        18.870667_RP, & ! 68
        22.009553_RP, & ! 69
        64.604291_RP, & ! 70
        11.498094_RP, & ! 71
         0.578834_RP, & ! 72
         9.237738_RP, & ! 73
        49.747842_RP, & ! 74
         2.147012_RP, & ! 75
         1.196895_RP, & ! 76
         2.133898_RP, & ! 77
         0.173168_RP  / ! 78

  data Prece_phase / &
        251.9025_RP, & !  1
        280.8325_RP, & !  2
        128.3057_RP, & !  3
        348.1074_RP, & !  4
        292.7252_RP, & !  5
        165.1686_RP, & !  6
        263.7951_RP, & !  7
         15.3747_RP, & !  8
         58.5749_RP, & !  9
         40.8226_RP, & ! 10
        308.4258_RP, & ! 11
        240.0099_RP, & ! 12
        222.9725_RP, & ! 13
        106.5937_RP, & ! 14
        114.5182_RP, & ! 15
        268.7809_RP, & ! 16
        279.6869_RP, & ! 17
         39.6448_RP, & ! 18
        126.4108_RP, & ! 19
        291.5795_RP, & ! 20
        307.2848_RP, & ! 21
         18.9300_RP, & ! 22
        273.7596_RP, & ! 23
        143.8050_RP, & ! 24
        191.8927_RP, & ! 25
        125.5237_RP, & ! 26
        172.7351_RP, & ! 27
        316.7998_RP, & ! 28
        319.6024_RP, & ! 29
         69.7526_RP, & ! 30
        123.5968_RP, & ! 31
        217.6432_RP, & ! 32
         85.5882_RP, & ! 33
        156.2147_RP, & ! 34
         66.9489_RP, & ! 35
         20.2082_RP, & ! 36
        250.7568_RP, & ! 37
         48.0188_RP, & ! 38
          8.3739_RP, & ! 39
         17.0374_RP, & ! 40
        155.3409_RP, & ! 41
         94.1709_RP, & ! 42
        221.1120_RP, & ! 43
         28.9300_RP, & ! 44
        117.1498_RP, & ! 45
        320.5095_RP, & ! 46
        262.3602_RP, & ! 47
        336.2148_RP, & ! 48
        233.0046_RP, & ! 49
        155.6977_RP, & ! 50
        184.6277_RP, & ! 51
        267.2772_RP, & ! 52
         78.9281_RP, & ! 53
        123.4722_RP, & ! 54
        188.7132_RP, & ! 55
        180.1364_RP, & ! 56
         49.1382_RP, & ! 57
        152.5268_RP, & ! 58
         98.2198_RP, & ! 59
         97.4808_RP, & ! 60
        221.5376_RP, & ! 61
        168.2438_RP, & ! 62
        161.1199_RP, & ! 63
         55.0196_RP, & ! 64
        262.6495_RP, & ! 65
        200.3284_RP, & ! 66
        201.6651_RP, & ! 67
        294.6547_RP, & ! 68
         99.8233_RP, & ! 69
        213.5577_RP, & ! 70
        154.1631_RP, & ! 71
        232.7153_RP, & ! 72
        138.3034_RP, & ! 73
        204.6609_RP, & ! 74
        106.5938_RP, & ! 75
        250.4676_RP, & ! 76
        332.3345_RP, & ! 77
         27.3039_RP  / ! 78

  real(RP), private, save :: arcsec2d, arcsec2r ! unit converter

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> setup solar incidence module
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SOLARINS_setup( &
       iyear )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_D2R
    implicit none

    integer, intent(in) :: iyear ! year at setup

    namelist / PARAM_ATMOS_SOLARINS / &
       ATMOS_SOLARINS_constant

    real(RP) :: dyear ! delta t [year]

    integer  :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[SOLARINS]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ solar insolation'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_SOLARINS,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_SOLARINS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_SOLARINS)

    arcsec2d = 1.0_RP / (60.0_RP*60.0_RP)
    arcsec2r = arcsec2d * CONST_D2R

    call ATMOS_SOLARINS_orbit(iyear)

    dyear = real( iyear-year_ref, kind=RP )

    !----- report data -----
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A)')       '*** insolation parameters'
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)')    '*** reference year      : ', year_ref
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)')    '*** current   year      : ', iyear
    if( IO_L ) write(IO_FID_LOG,'(1x,A,I7)')    '*** difference from ref.: ', int(dyear)
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.7)') '*** solar constant                 [W/m2]: ', ATMOS_SOLARINS_constant
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.7)') '*** obliquity                       [deg]: ', obliquity / CONST_D2R
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.7)') '*** eccentricity                         : ', E
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.7)') '*** longitude of perihelion         [deg]: ', omega     / CONST_D2R
    if( IO_L ) write(IO_FID_LOG,'(1x,A,F12.7)') '*** longitude at the vernal equinox [deg]: ', lambda_m0 / CONST_D2R

    return
  end subroutine ATMOS_SOLARINS_setup

  !-----------------------------------------------------------------------------
  !> setup solar incidence module
  !-----------------------------------------------------------------------------
  subroutine ATMOS_SOLARINS_orbit( &
       iyear )
    use scale_const, only: &
         PI => CONST_PI, &
         CONST_D2R
    implicit none

    integer, intent(in) :: iyear ! year at setup

    real(RP) :: dyear ! delta t [year]

    real(RP) :: EsinOMG     ! e * sin(w) : ecliptic parameter
    real(RP) :: EcosOMG     ! e * cos(w) : ecliptic parameter
    real(RP) :: Perih_pi    ! pi         : longitude of fixed perihelion [rad]
    real(RP) :: Perih_psi   ! psi        : general precession [degree]
    real(RP) :: Perih_omega ! omega_bar  : longitude of perihelion [degree]

    real(RP) :: temp, beta, EcosOMG_mod
    integer  :: i
    !---------------------------------------------------------------------------

    ! time from reference year(1950.0 AD)
    dyear = real( iyear-year_ref, kind=RP )

    ! obliquity
    temp = 0.0_RP
    do i = 1, nObliq
       temp = temp + Obliq_amp(i)*arcsec2d * cos( Obliq_rate(i)*arcsec2r*dyear + Obliq_phase(i)*CONST_D2R )
    enddo
    obliquity = ( obliquity_ref + temp ) * CONST_D2R

    ! eccentricity
    EsinOMG = 0.0_RP
    EcosOMG = 0.0_RP
    do i = 1, nEclip
       EsinOMG = EsinOMG + Eclip_amp(i) * sin( Eclip_rate(i)*arcsec2r*dyear + Eclip_phase(i)*CONST_D2R )
       EcosOMG = EcosOMG + Eclip_amp(i) * cos( Eclip_rate(i)*arcsec2r*dyear + Eclip_phase(i)*CONST_D2R )
    enddo
    E = sqrt( EsinOMG*EsinOMG + EcosOMG*EcosOMG )

    ! longitude of fixed perihelion
    if ( EcosOMG == 0.0_RP ) then
       EcosOMG_mod  = 1.E-30_RP
    else
       EcosOMG_mod  = EcosOMG
    endif
    Perih_pi = atan(EsinOMG/EcosOMG_mod) + PI * ( 0.5_RP - sign(0.5_RP, EcosOMG) )

!    if( abs(EcosOMG) < 1.E-8_RP ) then
!       if    ( EsinOMG == 0.0_RP ) then
!          Perih_pi = 0.0_RP
!       elseif( EsinOMG <  0.0_RP ) then
!          Perih_pi = 1.5_RP * PI
!       elseif( EsinOMG >  0.0_RP ) then
!          Perih_pi = 0.5_RP * PI
!       endif
!
!    elseif( EcosOMG < 0.0_RP ) then
!
!       Perih_pi = atan(EsinOMG/EcosOMG) + PI
!
!    elseif( EcosOMG > 0.0_RP ) then
!
!       if ( EsinOMG <  0.0_RP ) then
!          Perih_pi = atan(EsinOMG/EcosOMG) + 2.0_RP * PI
!       else
!          Perih_pi = atan(EsinOMG/EcosOMG)
!       endif
!
!    endif
    Perih_pi = Perih_pi / CONST_D2R ! [rad]->[degree]

    ! general precession
    temp = 0.0_RP
    do i = 1, nPrece
       temp = temp + Prece_amp(i)*arcsec2d * sin( Prece_rate(i)*arcsec2r*dyear + Prece_phase(i)*CONST_D2R )
    enddo
    Perih_psi = psi_bar * arcsec2d * dyear + zeta + temp

    ! longitude of perihelion
    Perih_omega = Perih_pi + Perih_psi

    do i = 1, 1000
       if    ( Perih_omega + 180.0_RP <    0.0_RP ) then
          Perih_omega = Perih_omega + 360.0_RP
       elseif( Perih_omega + 180.0_RP >= 360.0_RP ) then
          Perih_omega = Perih_omega - 360.0_RP
       else
          exit
       endif
    enddo

    ! longitude of perigee (see berger et al.(1993))
    omega = ( Perih_omega + 180.0_RP ) * CONST_D2R

    beta = sqrt( 1.0_RP - E*E )

    ! longitude at the vernal equinox
    lambda_m0 = 2.0_RP * ( ( 1.0_RP/2.0_RP*E + 1.0_RP/8.0_RP*E*E*E ) * ( 1.0_RP        + beta ) * sin(       omega) &
                         - (                   1.0_RP/4.0_RP*E*E   ) * ( 1.0_RP/2.0_RP + beta ) * sin(2.0_RP*omega) &
                         + (                   1.0_RP/8.0_RP*E*E*E ) * ( 1.0_RP/3.0_RP + beta ) * sin(3.0_RP*omega) )

    return
  end subroutine ATMOS_SOLARINS_orbit

  !-----------------------------------------------------------------------------
  !> calc factor of Earths solar insolation
  subroutine ATMOS_SOLARINS_insolation_0D( &
      solins,    &
      cosSZA,    &
      Re_factor, &
      lon,       &
      lat,       &
      now_date   )
    use scale_const, only: &
         PI => CONST_PI
    use scale_calendar, only: &
         CALENDAR_getDayOfYear,  &
         CALENDAR_ymd2absday,    &
         CALENDAR_hms2abssec,    &
         I_year, I_month, I_day, &
         I_hour, I_min, I_sec
    implicit none

    real(RP), intent(out) :: solins      ! solar insolation
    real(RP), intent(out) :: cosSZA      ! cos(Solar Zenith Angle)
    real(RP), intent(out) :: Re_factor   ! factor of the distance of Earth from the sun (1/rho2)
    real(RP), intent(in)  :: lon         ! longitude
    real(RP), intent(in)  :: lat         ! latitude
    integer,  intent(in)  :: now_date(6) ! date(yyyy,mm,dd,hh,mm,ss)

    real(RP) :: lambda_m       ! mean longitude from vernal equinox
    real(RP) :: lambda         !
    real(RP) :: sinDEC, cosDEC ! sin/cos(solar declination)
    real(RP) :: hourangle      ! hour angle: relative longitude of subsolar point

    integer  :: absday, absday_ve
    real(RP) :: DayOfYear, abssec
    real(RP) :: nu
    !---------------------------------------------------------------------------

    call CALENDAR_getDayOfYear( DayOfYear, now_date(I_year) )

    call CALENDAR_ymd2absday( absday,            & ! [OUT]
                              now_date(I_year),  & ! [IN]
                              now_date(I_month), & ! [IN]
                              now_date(I_day)    ) ! [IN]

    call CALENDAR_ymd2absday( absday_ve,         & ! [OUT]
                              now_date(I_year),  & ! [IN]
                              ve_date (I_month), & ! [IN]
                              ve_date (I_day)    ) ! [IN]

    call CALENDAR_hms2abssec( abssec,            & ! [OUT]
                              now_date(I_hour),  & ! [IN]
                              now_date(I_min),   & ! [IN]
                              now_date(I_sec),   & ! [IN]
                              0.0_RP             ) ! [IN]

    lambda_m = lambda_m0 + 2.0_RP * PI * real(absday-absday_ve,kind=RP) / DayOfYear

    nu = lambda_m - omega

    ! actual longitude from vernal equinox
    lambda = lambda_m &
           + ( 2.0_RP*E -  1.0_RP/ 4.0_RP*E*E*E ) * sin(        nu ) &
           + (             5.0_RP/ 4.0_RP*E*E   ) * sin( 2.0_RP*nu ) &
           + (            13.0_RP/12.0_RP*E*E*E ) * sin( 3.0_RP*nu )

    ! solar declination
    sinDEC = sin(lambda) * sin(obliquity)
    cosDEC = sqrt( 1.0_RP - sinDEC*sinDEC )

    ! hour angle
    hourangle = lon + 2.0_RP * PI * abssec / (24.0_RP*60.0_RP*60.0_RP)

    ! cos(Solar Zenith Angle)
    cosSZA = sin(lat)*sinDEC - cos(lat)*cosDEC*cos(hourangle)

    ! 1 / (rho*rho)
    Re_factor = 1.0_RP                          &
              + 2.0_RP*E     * cos(        nu ) &
              + 0.5_RP*E*E   * cos( 2.0_RP*nu ) &
              + 2.5_RP*E*E                      &
              + 4.0_RP*E*E*E * cos( 3.0_RP*nu )

    solins = ATMOS_SOLARINS_constant * Re_factor * max(cosSZA,0.0_RP)

    return
  end subroutine ATMOS_SOLARINS_insolation_0D

  !-----------------------------------------------------------------------------
  !> calc factor of Earths solar insolation
  subroutine ATMOS_SOLARINS_insolation_2D( &
      solins,    &
      cosSZA,    &
      lon,       &
      lat,       &
      now_date   )
    use scale_grid_index
    use scale_const, only: &
         PI => CONST_PI
    use scale_calendar, only: &
         CALENDAR_getDayOfYear,  &
         CALENDAR_ymd2absday,    &
         CALENDAR_hms2abssec,    &
         I_year, I_month, I_day, &
         I_hour, I_min, I_sec
    implicit none

    real(RP), intent(out) :: solins(IA,JA) ! solar insolation
    real(RP), intent(out) :: cosSZA(IA,JA) ! cos(Solar Zenith Angle)
    real(RP), intent(in)  :: lon   (IA,JA) ! longitude
    real(RP), intent(in)  :: lat   (IA,JA) ! latitude
    integer,  intent(in)  :: now_date(6)   ! date(yyyy,mm,dd,hh,mm,ss)

    real(RP) :: lambda_m       ! mean longitude from vernal equinox
    real(RP) :: lambda         !
    real(RP) :: sinDEC, cosDEC ! sin/cos(solar declination)
    real(RP) :: hourangle(IA,JA) ! hour angle: relative longitude of subsolar point

    integer  :: absday, absday_ve
    real(RP) :: DayOfYear, abssec
    real(RP) :: nu
    real(RP) :: Re_factor   ! factor of the distance of Earth from the sun (1/rho2)

    integer  :: i, j
    !---------------------------------------------------------------------------

    call CALENDAR_getDayOfYear( DayOfYear, now_date(I_year) )

    call CALENDAR_ymd2absday( absday,            & ! [OUT]
                              now_date(I_year),  & ! [IN]
                              now_date(I_month), & ! [IN]
                              now_date(I_day)    ) ! [IN]

    call CALENDAR_ymd2absday( absday_ve,         & ! [OUT]
                              now_date(I_year),  & ! [IN]
                              ve_date (I_month), & ! [IN]
                              ve_date (I_day)    ) ! [IN]

    call CALENDAR_hms2abssec( abssec,            & ! [OUT]
                              now_date(I_hour),  & ! [IN]
                              now_date(I_min),   & ! [IN]
                              now_date(I_sec),   & ! [IN]
                              0.0_RP             ) ! [IN]

    lambda_m = lambda_m0 + 2.0_RP * PI * real(absday-absday_ve,kind=RP) / DayOfYear

    nu = lambda_m - omega

    ! 1 / (rho*rho)
    Re_factor = 1.0_RP                          &
              + 2.0_RP*E     * cos(        nu ) &
              + 0.5_RP*E*E   * cos( 2.0_RP*nu ) &
              + 2.5_RP*E*E                      &
              + 4.0_RP*E*E*E * cos( 3.0_RP*nu )

    ! actual longitude from vernal equinox
    lambda = lambda_m &
           + ( 2.0_RP*E -  1.0_RP/ 4.0_RP*E*E*E ) * sin(        nu ) &
           + (             5.0_RP/ 4.0_RP*E*E   ) * sin( 2.0_RP*nu ) &
           + (            13.0_RP/12.0_RP*E*E*E ) * sin( 3.0_RP*nu )

    ! solar declination
    sinDEC = sin(lambda) * sin(obliquity)
    cosDEC = sqrt( 1.0_RP - sinDEC*sinDEC )

    do j = JS, JE
    do i = IS, IE
       ! hour angle
       hourangle(i,j) = lon(i,j) + 2.0_RP * PI * abssec / (24.0_RP*60.0_RP*60.0_RP)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       ! cos(Solar Zenith Angle)
       cosSZA(i,j) = sin(lat(i,j))*sinDEC - cos(lat(i,j))*cosDEC*cos(hourangle(i,j))
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       solins(i,j) = ATMOS_SOLARINS_constant * Re_factor * max(cosSZA(i,j),0.0_RP)
    enddo
    enddo

    return
  end subroutine ATMOS_SOLARINS_insolation_2D

end module scale_atmos_solarins
!-------------------------------------------------------------------------------
