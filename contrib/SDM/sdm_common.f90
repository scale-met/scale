!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM common variables
!!
!! @par Description
!!          Common variables used in the Super Droplet Method (SDM)
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-06-13 (S.Shima) [new] Separated common variables from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-14 (S.Shima) [rev] Removed unused variables: KMIN,KPLS,...,
!! @li      2014-07-18 (Y.Sato)  [add] add QTRC_sdm
!! @li      2014-12-12 (Y.Sato)  [mod] modify LatHet and DNS_RL as those used in SCALE Library
!! @li      2014-12-12 (Y.Sato)  [mod] modify characteristics of aeorosl from ammonium sulfate to ammonium bisulfate
!! @li      2014-12-19 (Y.Sato)  [mod] modify the location for defining LatHet and DNS_RL, and modify some typo
!! @li      2014-12-25 (Y.Sato)  [mod] modify LatHet to LH0
!! @li      2015-06-27 (S.Shima) [add] Add num_threads to store the threads number of auto parallelization on K/FX10
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio

  use scale_const, only: &
     ONE_PI => CONST_PI, &
     rrst   => CONST_R, &       ! Gas constant [J/(K*mol)]
     mass_air => CONST_Mdry, &  ! Molecular mass of air [g/mol]
     GasV_C => CONST_Rvap, &    ! Gas Constant of vapor [J/K/kg]
     LH0    => CONST_LH0, &
     CPvap  => CONST_CPvap, &
     CL     => CONST_CL, &
     TEM00  => CONST_TEM00, &
     DWATR  => CONST_DWATR

  use rng_uniform_mt, only: c_rng_uniform_mt

  !-----------------------------------------------------------------------------
  implicit none
  public
  !-----------------------------------------------------------------------------
  !
  !++ Parameters & variables
  !
  !-----------------------------------------------------------------------------
  character(len=H_LONG), save :: SD_IN_BASENAME = ''
  character(len=H_LONG), save :: SD_OUT_BASENAME = ''
  character(len=H_LONG), save :: RANDOM_IN_BASENAME = ''
  character(len=H_LONG), save :: RANDOM_OUT_BASENAME = ''
  type(c_rng_uniform_mt), save :: rng_s2c
  integer, save :: fid_sd_i, fid_sd_o
  integer, save :: fid_random_i, fid_random_o
  logical, save :: sd_rest_flg_in = .false. ! restart flg of Super Droplet
  logical, save :: sd_first = .true. 
  integer, save :: num_threads ! number of threads of auto parallelization on K/FX10
  !
  !++ Basic variable for SDM
  !
  !------------------------------------------------------------------------------
  integer(DP), allocatable, save :: sdn_s2c(:)   ! multipilicity
  real(RP), allocatable, save :: sdri_s2c(:)     ! index-i(real) of s.d.
  real(RP), allocatable, save :: sdrj_s2c(:)     ! index-j(real) of s.d.
  real(RP), allocatable, save :: sdrk_s2c(:)     ! index-k(real) of s.d.
  real(RP), allocatable, save :: sdx_s2c(:)      ! x-cordinate of s.d.
  real(RP), allocatable, save :: sdy_s2c(:)      ! y-cordinate of s.d.
  real(RP), allocatable, save :: sdz_s2c(:)      ! z-cordinate of s.d.
  real(RP), allocatable, save :: sdr_s2c(:)      ! y-cordinate of s.d.
  real(RP), allocatable, save :: sdu_s2c(:)      ! x-components velocity of s.d.
  real(RP), allocatable, save :: sdv_s2c(:)      ! y-components velocity of s.d.
  real(RP), allocatable, save :: sdvz_s2c(:)     ! z-components velocity of s.d. and terminal veloicty of s.d.
  real(RP), allocatable, save :: sdasl_s2c(:,:)  ! aeosol mass of s.d.
  real(RP), allocatable, save :: sdrkl_s2c(:,:)  ! index-k(real) at 'sdm_zlower'
  real(RP), allocatable, save :: sdrku_s2c(:,:)  ! index-k(real) at 'sdm_upper'
  integer(DP), allocatable, save :: sdn_tmp(:)   ! multiplicity of super-droplets
  real(RP), allocatable, save :: sdrk_tmp(:)     ! index-k(real) of super-droplets
  real(RP), allocatable, save :: sdx_tmp(:)      ! x-coordinate of super-droplets
  real(RP), allocatable, save :: sdy_tmp(:)      ! y-coordinate of super-droplets
  real(RP), allocatable, save :: sdz_tmp(:)      ! z-coordinate of super-droplets
  real(RP), allocatable, save :: sdr_tmp(:)      ! equivalent radius of super-droplets
  real(RP), allocatable, save :: sdu_tmp(:)      ! x-components velocity of super-droplets
  real(RP), allocatable, save :: sdv_tmp(:)      ! y-components velocity of super-droplets
  real(RP), allocatable, save :: sdvz_tmp(:)     ! terminal velocity of super-droplets zeta components of contravariant velocity
  real(RP), allocatable, save :: sdasl_tmp(:,:)  ! aerosol mass of super-droplets
  ! SDM for aerosol formation
  integer(DP), allocatable, save :: sdn_fm(:)    ! multiplicity of super-droplets
  real(RP), allocatable, save :: sdri_fm(:)      ! index-i(real) of super-droplets
  real(RP), allocatable, save :: sdrj_fm(:)      ! index-j(real) of super-droplets
  real(RP), allocatable, save :: sdrk_fm(:)      ! index-k(real) of super-droplets
  real(RP), allocatable, save :: sdx_fm(:)       ! x-coordinate of super-droplets
  real(RP), allocatable, save :: sdy_fm(:)       ! y-coordinate of super-droplets
  real(RP), allocatable, save :: sdz_fm(:)       ! z-coordinate of super-droplets
  real(RP), allocatable, save :: sdr_fm(:)       ! equivalent radius of super-droplets
  real(RP), allocatable, save :: sdvz_fm(:)      ! terminal velocity of super-droplets zeta components of contravariant velocity
  real(RP), allocatable, save :: sdasl_fm(:,:)   ! aerosol mass of super-droplets formed by gas-to-particle conversion
!!$  real(RP), allocatable, save :: rhod_crs(:,:,:) ! dry air density
!!$  real(RP), allocatable, save :: rhoc_sdm(:,:,:) ! density of cloud water
!!$  real(RP), allocatable, save :: rhor_sdm(:,:,:) ! density of rain water
!!$  real(RP), allocatable, save :: rhoa_sdm(:,:,:) ! density of aerosol
  real(RP), allocatable, save :: rand_s2c(:)
  integer, allocatable, save :: sortid_s2c(:)
  integer, allocatable, save :: sortfreq_s2c(:)
  integer, allocatable, save :: sortkey_s2c(:)
  integer, allocatable, save :: sorttag_s2c(:)
  real(DP), allocatable, save :: rbuf(:,:,:) ! We should prepare different rbuf and sbuf for multiplicity(long int) and other variables(real RP)
  real(DP), allocatable, save :: sbuf(:,:,:) ! Will be fixed in the near future
  integer, allocatable, save :: sdm_itmp1(:)
  integer, allocatable, save :: sdm_itmp2(:)
  integer, allocatable, save :: sdm_itmp3(:)
  integer, allocatable, save :: sd_itmp1(:,:)
  integer, allocatable, save :: sd_itmp2(:,:)
  integer, allocatable, save :: sd_itmp3(:,:)
  real(RP), allocatable, save :: sd_dtmp1(:)
  real(RP), allocatable, save :: QTRC_sdm(:,:,:,:)
  !------------------------------------------------------------------------------
  !
  !++ Basic variable for SDM
  !
  !------------------------------------------------------------------------------
  real(RP) :: xmax_sdm, ymax_sdm      ! Distance in X, Y-direction
  real(RP) :: rkumax, rkumin, rklmax  ! Maxinum and Minuimum real number of sd_rku and sd_rkl
  real(RP) :: minzph                  ! Minimum of zph(i,j,2)
  integer :: knum_sdm                ! Maximum number of SDM grid in vertical
!  integer :: ni=IE-IS+1, nj=JE-JS+1, nz=KE-KS+1     ! Number of grid for each dirrection with halo
  integer :: stat
  integer :: dstw_sub, dste_sub, dsts_sub, dstn_sub
  integer :: srcw_sub, srce_sub, srcs_sub, srcn_sub
  !------------------------------------------------------------------------------
  !
  !++ Basic parameter for SDM
  !
  !------------------------------------------------------------------------------
!!$  real(RP), parameter :: LatHet = 2.453E+6_RP   ! Latent heat of water at 293K [J/kg]
!!$  real(RP), parameter :: LatHet = LH0 - ( CPvap - CL )*TEM00   ! Latent heat of water at 0K [J/kg]
  real(RP), parameter :: LatHet = LH0           ! Latent heat of water at 273.15K [J/kg] 
  real(RP), parameter :: Heat_C = 2.550E-2_RP   ! thermal conductivity at 293K, 100kPa [J/m s K]
  !--- Y.Sato
!!$  real(RP), parameter :: DNS_RL = 998.203_RP ! Density of liquid water at 293K [kg m-3]
  real(RP), parameter :: DNS_RL = DWATR         ! Density of liquid water [kg m-3]
  !--- Y.Sato
  real(RP), parameter :: Diff_C = 2.52E-5_RP    ! Diffusion Constant
  !--- Y.Sato
  real(RP), parameter :: LatGas = LatHet / GasV_C
  real(RP), parameter :: L_RL_K = LatHet * DNS_RL / Heat_C
  !--- Y.Sato
  real(RP), parameter :: RLRv_D = DNS_RL * GasV_C / Diff_C
!  real(RP), parameter :: exp_K = GasD_C / cp_C    ! exponent k used in adiabatic process
  real(RP), parameter :: m2micro = 1.0E+6_RP    ! Convert unit [m] -> [micron]
  real(RP), parameter :: micro2m = 1.0E-6_RP    ! Convert unit [micron] -> [m]
  real(RP), parameter :: TWO_PI = 2.0_RP * 3.14159265358979_RP
  real(RP), parameter :: O_THRD = 1.0_RP / 3.0_RP
  real(RP), parameter :: F_THRD = 4.0_RP / 3.0_RP
  real(RP), parameter :: O_SIX  = 1.0_RP / 6.0_RP
  real(RP), parameter :: T_EIGH = 3.0_RP / 8.0_RP
  ! thereshold value between valid super-droplet and invalid super-droplets
  real(RP), parameter :: VALID2INVALID = -999.0_RP 
  real(RP), parameter :: PRECIPI = -999.9_RP  ! value indicated as super-droplets is precipitation
  real(RP), parameter :: PREC2INVALID = -999.99_RP ! thereshold value between precipitation and invalid super-droplets
  real(RP), parameter :: INVALID = -999.999_RP ! value indicated as invalid super-droplets
  !------------------------------------------------------------------------------
  !
  !++ Parameter relating to microphysics
  !
  !------------------------------------------------------------------------------
  real(RP), parameter :: boltz = 1.38066E-23_RP  ! Boltzmann constant
  ! Molecular mass of sea salt (NaCl) contained as water-solble aerosol in S.D. [g]
  real(RP), parameter :: mass_nacl = 58.44277_RP 
!!  ! Molecular mass of ammonium sulfate contained as water-solble aerosol in S.D. [g]
!!  real(RP), parameter :: mass_amsul = 132.14_RP  
  ! Molecular mass of ammonium bisulfate contained as water-solble aerosol in S.D. [g]
  real(RP), parameter :: mass_amsul = 115.11_RP  
  ! Degree of ion dissociation of ammonium sulfate and sea salt contained 
  ! as water-soluble aerosol in super droplets
  real(RP), parameter :: ion_amsul = 2.0_RP, ion_nacl = 2.0_RP 
!!  ! Density of ammonium sulfate and sea salt [g/m3] contained as water-soluble aerosol in super droplets
!!  real(RP), parameter :: rho_amsul = 1.769E+6_RP, rho_nacl = 2.165E+6_RP
  ! Density of ammonium bisulfate and sea salt [g/m3] contained as water-soluble aerosol in super droplets
  real(RP), parameter :: rho_amsul = 1.78E+6_RP, rho_nacl = 2.165E+6_RP
  ! parameter for 3mode log-normal distribution of aerosol
  real(RP)  :: n1_amsul                                  !! [1/m3]
  real(RP)  :: n2_amsul
  real(RP)  :: n3_nacl
  real(RP)  :: rb1_amsul                                 !! [m]
  real(RP)  :: rb2_amsul
  real(RP)  :: rb3_nacl
  real(RP)  :: sgm1_amsul                                !! [-]
  real(RP)  :: sgm2_amsul
  real(RP)  :: sgm3_nacl
  real(RP)  :: rmax_amsul                                !! [m]
  real(RP)  :: rmax_nacl
  real(RP)  :: rmin_amsul
  real(RP)  :: rmin_nacl
  ! parameter for 3mode log-normal distribution under RICO-observation : Derksen(2009)
  real(RP), parameter :: n1_amsul_derksn   = 118.E+6_RP   !! [1/m3]
  real(RP), parameter :: n2_amsul_derksn   = 11.E+6_RP
  real(RP), parameter :: n3_nacl_derksn    = 0.1E+6_RP
  real(RP), parameter :: rb1_amsul_derksn  = 1.9E-08_RP   !! [m]
  real(RP), parameter :: rb2_amsul_derksn  = 5.6E-08_RP
  real(RP), parameter :: rb3_nacl_derksn   = 8.0E-07_RP
  real(RP), parameter :: sgm1_amsul_derksn = 3.3_RP       !! [-]
  real(RP), parameter :: sgm2_amsul_derksn = 1.6_RP
  real(RP), parameter :: sgm3_nacl_derksn  = 2.2_RP
  real(RP), parameter :: rmax_amsul_derksn = 5.0E-06_RP   !! [m]
  real(RP), parameter :: rmax_nacl_derksn  = 5.0E-06_RP
  real(RP), parameter :: rmin_amsul_derksn = 2.0E-09_RP
  real(RP), parameter :: rmin_nacl_derksn  = 1.0E-07_RP
  ! parameter for 3mode log-nomiral distribution under RICO-observation : vanZanten(2010)
  real(RP), parameter :: n1_amsul_zanten   = 90.E+6_RP    !! [1/m3]
  real(RP), parameter :: n2_amsul_zanten   = 15.E+6_RP
  real(RP), parameter :: rb1_amsul_zanten  = 3.0E-08_RP   !! [m]
  real(RP), parameter :: rb2_amsul_zanten  = 1.4E-07_RP
  real(RP), parameter :: sgm1_amsul_zanten = 1.28_RP      !! [-]
  real(RP), parameter :: sgm2_amsul_zanten = 1.75_RP
  real(RP), parameter :: rmax_amsul_zanten = 5.0E-06_RP   !! [m]
  real(RP), parameter :: rmin_amsul_zanten = 1.0E-08_RP
  ! parameter for 3mode log-nomiral distribution under RICO-observation : modified vanZanten(2010)
  real(RP), parameter :: rb_amsul_zanten_m   = 3.0E-08_RP   !! [m]
  real(RP), parameter :: sgm_amsul_zanten_m  = 1.28_RP      !! [-]
  real(RP), parameter :: rmax_amsul_zanten_m = 1.0E-07_RP   !! [m]
  real(RP), parameter :: rmin_amsul_zanten_m = 1.0E-08_RP
  real(RP), parameter :: Es_T_A = 2.53E+11_RP  ! A[kg/m s2] : Es = A * exp(-B/T)
  real(RP), parameter :: Es_T_B = 5.42E+3_RP   ! B[K] : Es = A * exp(-B/T)
  real(RP), parameter :: CurveF = 3.3E-7_RP    !  Curvature term Const. a=CoruveF/T [m K]
  real(RP), parameter :: ASL_FF = 4.3E-6_RP    !  Cloud Condensation term Const. b=ASL_FF*(Ion*M)/ms [m3]
  !------------------------------------------------------------------------------
  !
  !++ Terminal velocity
  !
  !------------------------------------------------------------------------------
  real(RP) :: vz_b(7), vz_c(6)
  data vz_b / -0.318657E+1_RP,  0.9926960_RP,    -0.153193E-2_RP, -0.987059E-3_RP, &
              -0.578878e-3_RP,  0.855176E-4_RP,  -0.327815E-5_RP /
  data vz_c / -0.500015E+1_RP,  0.523778E+1_RP,  -0.204914E+1_RP,              &
               0.47529400_RP,  -0.542819E-1_RP,   0.238449E-2_RP /
  !------------------------------------------------------------------------------
  !
  !++ Collision/Coalescence
  !
  !------------------------------------------------------------------------------
  real(RP) :: ratcol(21), r0col(15), ecoll(15,21)
  data r0col /   6.0_RP,  8.0_RP, 10.0_RP, 15.0_RP,  20.0_RP,  25.0_RP,  30.0_RP,  &
                40.0_RP, 50.0_RP, 60.0_RP, 70.0_RP, 100.0_RP, 150.0_RP, 200.0_RP, &
                300.0_RP /
  data ratcol / 0.00_RP, 0.050_RP, 0.10_RP, 0.150_RP, 0.20_RP, 0.250_RP,        &
                0.30_RP, 0.350_RP, 0.40_RP, 0.450_RP, 0.50_RP, 0.550_RP,        &
                0.60_RP, 0.650_RP, 0.70_RP, 0.750_RP, 0.80_RP, 0.850_RP,        &
                0.90_RP, 0.950_RP, 1.00_RP /
  data ecoll /                                                      &
      0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP  &
     ,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0010_RP,0.0030_RP,0.0030_RP,0.0030_RP,0.0040_RP,0.0050_RP  &
     ,0.0050_RP,0.0050_RP,0.0100_RP,0.1000_RP,0.0500_RP,0.2000_RP,0.5000_RP,0.7700_RP,0.8700_RP,0.9700_RP  &
     ,0.0070_RP,0.0070_RP,0.0070_RP,0.0080_RP,0.0090_RP,0.0100_RP,0.0100_RP,0.0700_RP,0.4000_RP,0.4300_RP  &
     ,0.5800_RP,0.7900_RP,0.9300_RP,0.9600_RP,1.0000_RP,0.0090_RP,0.0090_RP,0.0090_RP,0.0120_RP,0.0150_RP  &
     ,0.0100_RP,0.0200_RP,0.2800_RP,0.6000_RP,0.6400_RP,0.7500_RP,0.9100_RP,0.9700_RP,0.9800_RP,1.0000_RP  &
     ,0.0140_RP,0.0140_RP,0.0140_RP,0.0150_RP,0.0160_RP,0.0300_RP,0.0600_RP,0.5000_RP,0.7000_RP,0.7700_RP  &
     ,0.8400_RP,0.9500_RP,0.9700_RP,1.0000_RP,1.0000_RP,0.0170_RP,0.0170_RP,0.0170_RP,0.0200_RP,0.0220_RP  &
     ,0.0600_RP,0.1000_RP,0.6200_RP,0.7800_RP,0.8400_RP,0.8800_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0300_RP,0.0300_RP,0.0240_RP,0.0220_RP,0.0320_RP,0.0620_RP,0.2000_RP,0.6800_RP,0.8300_RP,0.8700_RP  &
     ,0.9000_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0250_RP,0.0250_RP,0.0250_RP,0.0360_RP,0.0430_RP  &
     ,0.1300_RP,0.2700_RP,0.7400_RP,0.8600_RP,0.8900_RP,0.9200_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0270_RP,0.0270_RP,0.0270_RP,0.0400_RP,0.0520_RP,0.2000_RP,0.4000_RP,0.7800_RP,0.8800_RP,0.9000_RP  &
     ,0.9400_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0300_RP,0.0300_RP,0.0300_RP,0.0470_RP,0.0640_RP  &
     ,0.2500_RP,0.5000_RP,0.8000_RP,0.9000_RP,0.9100_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0400_RP,0.0400_RP,0.0330_RP,0.0370_RP,0.0680_RP,0.2400_RP,0.5500_RP,0.8000_RP,0.9000_RP,0.9100_RP  &
     ,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0350_RP,0.0350_RP,0.0350_RP,0.0550_RP,0.0790_RP  &
     ,0.2900_RP,0.5800_RP,0.8000_RP,0.9000_RP,0.9100_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0370_RP,0.0370_RP,0.0370_RP,0.0620_RP,0.0820_RP,0.2900_RP,0.5900_RP,0.7800_RP,0.9000_RP,0.9100_RP  &
     ,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0370_RP,0.0370_RP,0.0370_RP,0.0600_RP,0.0800_RP  &
     ,0.2900_RP,0.5800_RP,0.7700_RP,0.8900_RP,0.9100_RP,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0370_RP,0.0370_RP,0.0370_RP,0.0410_RP,0.0750_RP,0.2500_RP,0.5400_RP,0.7600_RP,0.8800_RP,0.9200_RP  &
     ,0.9500_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0370_RP,0.0370_RP,0.0370_RP,0.0520_RP,0.0670_RP  &
     ,0.2500_RP,0.5100_RP,0.7700_RP,0.8800_RP,0.9300_RP,0.9700_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP  &
     ,0.0370_RP,0.0370_RP,0.0370_RP,0.0470_RP,0.0570_RP,0.2500_RP,0.4900_RP,0.7700_RP,0.8900_RP,0.9500_RP  &
     ,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,1.0000_RP,0.0360_RP,0.0360_RP,0.0360_RP,0.0420_RP,0.0480_RP  &
     ,0.2300_RP,0.4700_RP,0.7800_RP,0.9200_RP,1.0000_RP,1.0200_RP,1.0200_RP,1.0200_RP,1.0200_RP,1.0200_RP  &
     ,0.0400_RP,0.0400_RP,0.0350_RP,0.0330_RP,0.0400_RP,0.1120_RP,0.4500_RP,0.7900_RP,1.0100_RP,1.0300_RP  &
     ,1.0400_RP,1.0400_RP,1.0400_RP,1.0400_RP,1.0400_RP,0.0330_RP,0.0330_RP,0.0330_RP,0.0330_RP,0.0330_RP  &
     ,0.1190_RP,0.4700_RP,0.9500_RP,1.3000_RP,1.7000_RP,2.3000_RP,2.3000_RP,2.3000_RP,2.3000_RP,2.3000_RP  &
     ,0.0270_RP,0.0270_RP,0.0270_RP,0.0270_RP,0.0270_RP,0.1250_RP,0.5200_RP,1.4000_RP,2.3000_RP,3.0000_RP  &
     ,4.0000_RP,4.0000_RP,4.0000_RP,4.0000_RP,4.0000_RP /
  !------------------------------------------------------------------------------
  !
  !++ For Namelist
  !
  !------------------------------------------------------------------------------
  real(RP), save :: sdm_dtcmph(3)              ! Time interval [sec] of condensation, coalescence, and droplet motion
  logical, save :: sdm_calvar(3)              ! flag for calculation of condensation, coalecscence, and droplet motion
  real(RP), save :: sdm_rdnc     = 1000.0_RP   ! Number of real droplet
  real(RP), save :: sdm_sdnmlvol = 1000.0_RP   ! Normal volume for number concentration of S.D. [m3]
  real(RP), save :: sdm_inisdnc  = 1000.0_RP   ! Number of super droplets per sdm_nmlvol at initial [1/m^3]
  real(RP), save :: sdm_zlower   = 5.0_RP      ! Lowest limitation of super droplet [m]
  real(RP), save :: sdm_zupper   = 2.E+3_RP    ! Hightest limitation of super droplet [m]
  integer, save :: sdm_extbuf   = 30          ! Rate to buffer size of super droplets for extra [%]
  real(RP), save :: sdm_rqc2qr   = 1.E-5_RP    ! Threshold of radius between qc and qr [m]
  integer, save :: sdm_colkrnl  = 2           ! Kernel type for coalescence process (0:Golovin, 1:Long, 2:Halls,3:No col_effi hydrodynamic)
  integer, save :: sdm_colbrwn  = 0           ! Flag of Brownian Coagulation and Scavenge process
  integer, save :: sdm_mvexchg  = 0           ! flag of exchange momentum betweeen super-droplets and fluid 0:No, 1:Yes
  integer, save :: sdm_aslset   = 1.0_RP      ! Conrol flag to set the species and way of chemical material as water-soluble aerosol
  real(RP), save :: sdm_aslfmdt  = 0.1_RP      ! time interval [sec] of aerosol nucleation
  real(RP), save :: sdm_aslfmsdnc = 1000.0_RP  ! Number of S.D. per sdnmvol as aeroosl(s)
  real(RP), save :: sdm_aslfmrate = 0.0_RP ! Formation rate of super droplets as aerosol [1/(m^3 s)]
  real(RP), save :: sdm_aslmw(1:20)      ! User specified molecular mass of chemical material contained as water-soluble aerosol in S.D.
  real(RP), save :: sdm_aslion(1:20)     ! User specified ion of chemical material contained as water-soluble aerosol in S.D.
  real(RP), save :: sdm_aslrho(1:20)     ! User specified density of chemical material contained as water-soluble aerosol in S.D.
  real(RP), save :: sdm_nadjdt = 0.1_RP  ! Time interval of adjust S.D. number [s]
  integer, save :: sdm_nadjvar = 0      ! Control flag of adjustment super-droplet number in each grid
                                            ! 0:No adjust number of droplet
                                            ! 1:adjust number of droplet by adding
                                            ! 2:adjust number of droplet by removal 
                                            ! 3:adjust number of droplet by adding and removal
  integer, save :: sdm_dmpvar  = 000   ! Control flag to output Super Droplets
                                       ! 1st digit (10^0) corresponds to text output
                                       ! 2nd digit (10^1) corresponds to binary output
                                       ! 3rd digit (10^2) corresponds to binary output of large droplets 
                                       ! 0: off, 1: on, 2: output with sort data (only for binary output)
                                       ! currently only 001 is supported   
  real(RP), save :: sdm_dmpitva  = 0.e0 ! Time interval of text output [s]
  integer, save :: sdm_dmpnskip = 1     ! Base skip to store super droplets in text format
  real(RP), save :: sdm_dmpitvb  = 0.e0  ! Time interval of binary output of all droplets [s]
  real(RP), save :: sdm_dmpitvl  = 0.e0  ! Time interval of binary output of large droplets [s]
  real(RP), save :: sdm_dmpsdsiz = 0.e0  ! Threshold radius to store large super droplets in binary format [m]

  data sdm_dtcmph / 0.1_RP,0.1_RP,0.1_RP /
  data sdm_calvar / .false.,.false.,.false. /
  data sdm_aslmw  / 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP  /
  data sdm_aslion / 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP  /
  data sdm_aslrho / 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, &
                    0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP  /
  real(RP), allocatable, save  :: prr_crs(:,:,:) ! prr_crs(IA,JA,1:2) Precipitation and accumulation for rain
  real(RP), allocatable, save  :: zph_crs(:,:,:), dxiv_sdm(:), dyiv_sdm(:)!, dziv_sdm(:)
  real(RP), allocatable, save  :: dx_sdm(:), dy_sdm(:)!, dz_sdm(:)   ! Dx, Dy, Dz for SDM (normally they are equal to those of SCALE)
!  integer, parameter :: nqw = QQA
  integer, save :: sdfmnum_s2c
  real(RP), save:: sdininum_s2c
  integer, save :: sdnumasl_s2c, sdnum_s2c
  integer, save :: ni_s2c, nj_s2c, nk_s2c
  integer, save :: tag
  integer, parameter :: nomp = 1
  integer, save :: wbc=1, ebc=1, sbc=1, nbc=1  ! only periodic boundary is applied
  integer, save :: nsub   ! Number of sub domain in group domain
  integer, save :: nisub, njsub
  integer, save :: sthopt=1, trnopt=0
!  real(RP), save :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA)
!  real(RP), save :: crs_dtmp3(KA,IA,JA), crs_dtmp4(KA,IA,JA)
!  real(RP), save :: crs_dtmp5(KA,IA,JA), crs_dtmp6(KA,IA,JA)
  real(RP), save :: sdm_dtevl  ! time step of {condensation/evaporation} process
  real(RP), save :: sdm_dtcol  ! time step of {stochastic coalescence} process
  real(RP), save :: sdm_dtadv  ! time step of {motion of super-droplets} process
!  real(RP), allocatable, save :: qwtr_crs(:,:,:,:)
!  real(RP), allocatable, save :: j31(:,:,:) ! z-x components of Jacobian
!  real(RP), allocatable, save :: j32(:,:,:) ! z-y components of Jacobian
!  real(RP), allocatable, save :: jcb(:,:,:) ! Jacobian at scalar points
!  real(RP), allocatable, save :: jcb8w(:,:,:) ! Jacobian at w points
!  real(RP), allocatable, save :: mf(:,:)  ! Map scale factors
!!$  integer,  allocatable, save :: KMIN1(:), IMIN1(:), JMIN1(:)
!!$  integer,  allocatable, save :: KPLS1(:), IPLS1(:), JPLS1(:)
  integer, save :: nclstp(0:3)

  logical, save :: docondensation = .true.
  logical, save :: doautoconversion = .true.
  logical, save :: doprecipitation  = .true.
  logical, save :: domovement = .true.
  logical, save :: donegative_fixer  = .true.  ! apply negative fixer?
end module m_sdm_common
