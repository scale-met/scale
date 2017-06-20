!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by 6 water category, double moment bulk scheme
!!          Seiki and Nakajima(2014) J. Atmos. Sci., 71, 833–853
!!
!!          Reference:  -- Journals
!!                       Seifert and Beheng(2006)  : Meteorol.Atmos.Phys.,vol.92,pp.45-66
!!                       Seifert and Beheng(2001)  : Atmos.Res.,vol.59-60,pp.265-281
!!                       Seifert(2008)             : J.Atmos.Sci.,vol.65,pp.3608-3619
!!                       Lin et al.(1983)          : J.Appl.Meteor.,vol.22,pp.1065-1092
!!                       Ruttledge and Hobbs(1983) : J.Atmos.Sci.,vol.40,pp.1185-1206
!!                       Ruttledge and Hobbs(1984) : J.Atmos.Sci.,vol.40,pp.2949-2977
!!                       Cotton etal.(1986)        : J.C.Appl.Meteor.,25,pp.1658-1680
!!                       Cotton and Field (2002)   : QJRMS.,vol.128,pp2417-pp2437
!!                       Beard(1980)               : J.Atmos.Sci.,vol.37,pp.1363-1374 [Add] 10/08/03
!!                       Berry and Reinhardt(1974a): J.Atmos.Sci.,vol.31,pp.1814-1824
!!                       Berry and Reinhardt(1974b): J.Atmos.Sci.,vol.31,pp.1825-1831
!!                       Fu(1996)                  : J.Climate, vol.9, pp.2058-2082   [Add] 10/08/03
!!                       Fu etal(1998)             : J.Climate, vol.11, pp.2223-2237  [Add] 10/08/03
!!                       Ghan etal.(1997)          : J.Geophys.Res.,vol.102,pp.21777-21794, [Add] 09/08/18
!!                       Hong et al.(2004)         : Mon.Wea.Rev.,pp.103-120
!!                       Heymsfeild and Iaquinta(2000): J.Atmos.Sci., vol.57, pp.916-938 [Add] 10/08/03
!!                       Heymsfield and Kajikawa(1987): J.Atmos.Sci., vol.44, pp.1088-1099
!!                       Johnson(1981)             : J.Atmos.Sci., vol.38, pp.215-218 [Add] 09/08/18
!!                       McFarquhar and Heymsfield(1996): J.Atmos.Sci.,vol.53,pp.2401-2423
!!                       Mitchell(1996)            : J.Atmos.Sci., vol.53, pp.1710-1723. [Add] 10/08/03
!!                       Morrison etal.(2005)      : Mon.Wea.Rev.,vol.62,pp.1665-1677, [Add] 09/08/18
!!                       Locatelli and Hobbs (1974): J.Geophys.Res., vol.70, pp.2185-2197
!!                       Lohmann(2002)             : J.Atmos.Sci.,vol.59,pp.647-656
!!                       Takano and Liou(1989)     : J.Atmos.Sci.,vol.46,pp.3-19
!!                       Takano and Liou(1994)     : J.Atmos.Sci.,vol.52,pp.818-837
!!                       Auer and Veal(1970)       : J.Atmos.Sci.,vol.27,pp.919-926
!!                       Ikawa et al.(1991)        : J.M.S.J.,vol.69,pp.641-667
!!                       Murakami(1990)            : J.M.S.J.,vol.68,pp.107-128
!!                      -- Books
!!                       Pruppacher and Klett(1997): Kluwer Academic Publishers
!!                          Microphysics of Clouds and Precipitation, 2nd.edit.
!!                       Seinfeld and Pandis(1998) : Wiley Interscience
!!                          Atmospheric Chemistry and Physics.
!!                       Jacobson(2005)            : Cambridge press
!!                          Fundamentals of Atmospheric Modeling, 2nd.edit.
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-10-24 (T.Seiki)    [new] import from NICAM(11/08/30 ver.)
!! @li      2015-09-08 (Y.Sato)     [add] Add evaporated cloud number concentration
!!
!<
!-------------------------------------------------------------------------------

#ifdef PROFILE_FAPP
#define PROFILE_START(name) call fapp_start(name, 1, 1)
#define PROFILE_STOP(name)  call fapp_stop (name, 1, 1)
#elif defined(PROFILE_FINEPA)
#define PROFILE_START(name) call start_collection(name)
#define PROFILE_STOP(name)  call stop_collection (name)
#else
#define PROFILE_START(name)
#define PROFILE_STOP(name)
#endif

#include "macro_thermodyn.h"
module scale_atmos_phy_mp_sn14
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_const, only: &
     GRAV   => CONST_GRAV,   &
     PI     => CONST_PI,     &
     UNDEF8 => CONST_UNDEF8, &
     Rdry   => CONST_Rdry,   &
     CPdry  => CONST_CPdry,  &
     CVdry  => CONST_CVdry,  &
     P00    => CONST_PRE00,  &
     T00    => CONST_TEM00,  &
     Rvap   => CONST_Rvap,   &
     CPvap  => CONST_CPvap,  &
     CVvap  => CONST_CVvap,  &
     CL     => CONST_CL,     &
     CI     => CONST_CI,     &
     LHV    => CONST_LHV00,  &
     LHF    => CONST_LHF00,  &
     LHV0   => CONST_LHV0,   &
     LHF0   => CONST_LHF0,   &
     LHS0   => CONST_LHS0,   &
     LHV00  => CONST_LHV00,  &
     LHF00  => CONST_LHF00,  &
     PSAT0  => CONST_PSAT0,  &
     EMELT  => CONST_EMELT,  &
     DWATR  => CONST_DWATR

  use scale_atmos_hydrometeor, only: &
     N_HYD, &
     I_QV,  &
     I_QC,  &
     I_QR,  &
     I_QI,  &
     I_QS,  &
     I_QG,  &
     I_NC,  &
     I_NR,  &
     I_NI,  &
     I_NS,  &
     I_NG,  &
     I_HC,  &
     I_HR,  &
     I_HI,  &
     I_HS,  &
     I_HG
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_sn14_config
  public :: ATMOS_PHY_MP_sn14_setup
  public :: ATMOS_PHY_MP_sn14
  public :: ATMOS_PHY_MP_sn14_CloudFraction
  public :: ATMOS_PHY_MP_sn14_EffectiveRadius
  public :: ATMOS_PHY_MP_sn14_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: QA_MP  = 11

  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_sn14_NAME(QA_MP)
  character(len=H_MID)  , public, target :: ATMOS_PHY_MP_sn14_DESC(QA_MP)
  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_sn14_UNIT(QA_MP)

  real(RP), public, target :: ATMOS_PHY_MP_sn14_DENS(N_HYD) ! hydrometeor density [kg/m3]=[g/L]

  data ATMOS_PHY_MP_sn14_NAME / &
                 'QV', &
                 'QC', &
                 'QR', &
                 'QI', &
                 'QS', &
                 'QG', &
                 'NC', &
                 'NR', &
                 'NI', &
                 'NS', &
                 'NG'  /

  data ATMOS_PHY_MP_sn14_DESC / &
                 'Ratio of Water Vapor mass to total mass (Specific humidity)',   &
                 'Ratio of Cloud Water mass to total mass',   &
                 'Ratio of Rain Water mass to total mass',    &
                 'Ratio of Cloud Ice mass to total mass',     &
                 'Ratio of Snow mass to total mass',          &
                 'Ratio of Graupel mass to total mass',       &
                 'Cloud Water Number Density', &
                 'Rain Water Number Density',  &
                 'Cloud Ice Number Density',   &
                 'Snow Number Density',        &
                 'Graupel Number Density'      /

  data ATMOS_PHY_MP_sn14_UNIT / &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'kg/kg',  &
                 'num/kg', &
                 'num/kg', &
                 'num/kg', &
                 'num/kg', &
                 'num/kg'  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: mp_sn14_init
  private :: mp_sn14
  private :: MP_terminal_velocity

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !

  integer, private, parameter :: HYDRO_MAX = 5

  integer, private, parameter :: I_mp_QC = 1
  integer, private, parameter :: I_mp_QR = 2
  integer, private, parameter :: I_mp_QI = 3
  integer, private, parameter :: I_mp_QS = 4
  integer, private, parameter :: I_mp_QG = 5
  integer, private, parameter :: I_mp_NC = 6
  integer, private, parameter :: I_mp_NR = 7
  integer, private, parameter :: I_mp_NI = 8
  integer, private, parameter :: I_mp_NS = 9
  integer, private, parameter :: I_mp_NG = 10

  integer, private :: QS_MP
  integer, private :: QE_MP

  ! production rate
  !   nucleation
  integer, parameter :: I_LCccn =  1
  integer, parameter :: I_NCccn =  2
  integer, parameter :: I_LIccn =  3
  integer, parameter :: I_NIccn =  4
  !   freezing
  integer, parameter :: I_LChom =  5
  integer, parameter :: I_NChom =  6
  integer, parameter :: I_LChet =  7
  integer, parameter :: I_NChet =  8
  integer, parameter :: I_LRhet =  9
  integer, parameter :: I_NRhet = 10
  !   melting
  integer, parameter :: I_LImlt = 11
  integer, parameter :: I_NImlt = 12
  integer, parameter :: I_LSmlt = 13
  integer, parameter :: I_NSmlt = 14
  integer, parameter :: I_LGmlt = 15
  integer, parameter :: I_NGmlt = 16
  !   vapor deposition
  integer, parameter :: I_LRdep = 17
  integer, parameter :: I_NRdep = 18
  integer, parameter :: I_LIdep = 19
  integer, parameter :: I_NIdep = 20
  integer, parameter :: I_LSdep = 21
  integer, parameter :: I_NSdep = 22
  integer, parameter :: I_LGdep = 23
  integer, parameter :: I_NGdep = 24
  integer, parameter :: I_LCdep = 25
!    integer, parameter :: I_NCdep = 26
  !   warm collection process
  !     auto-conversion
  integer, parameter :: I_LCaut = 26
  integer, parameter :: I_NCaut = 27
  integer, parameter :: I_NRaut = 28
  !     accretion
  integer, parameter :: I_LCacc = 29
  integer, parameter :: I_NCacc = 30
  !     self-colletion, break-up
  integer, parameter :: I_NRslc = 31
  integer, parameter :: I_NRbrk = 32

  !   partial conversion(ice, snow => graupel)
  integer, parameter :: I_LIcon = 33
  integer, parameter :: I_NIcon = 34
  integer, parameter :: I_LScon = 35
  integer, parameter :: I_NScon = 36
  !   enhanced melting due to
  integer, parameter :: I_LIacm = 37 ! ice-cloud
  integer, parameter :: I_NIacm = 38
  integer, parameter :: I_LIarm = 39 ! ice-rain
  integer, parameter :: I_NIarm = 40
  integer, parameter :: I_LSacm = 41 ! snow-cloud
  integer, parameter :: I_NSacm = 42
  integer, parameter :: I_LSarm = 43 ! snow-rain
  integer, parameter :: I_NSarm = 44
  integer, parameter :: I_LGacm = 45 ! graupel-cloud
  integer, parameter :: I_NGacm = 46
  integer, parameter :: I_LGarm = 47 ! graupel-rain
  integer, parameter :: I_NGarm = 48
  !   ice multiplication by splintering
  integer, parameter :: I_LGspl = 49
  integer, parameter :: I_LSspl = 50
  integer, parameter :: I_NIspl = 51

  integer, parameter :: PQ_MAX = 51

  ! production rate of mixed-phase collection process
  ! PXXacYY2ZZ means XX collect YY produce ZZ
  integer, parameter :: I_LIacLC2LI   =  1 ! cloud-ice
  integer, parameter :: I_NIacNC2NI   =  2
  integer, parameter :: I_LSacLC2LS   =  3 ! cloud-snow(cloud change)
  integer, parameter :: I_NSacNC2NS   =  4
  integer, parameter :: I_LGacLC2LG   =  5 ! cloud-graupel
  integer, parameter :: I_NGacNC2NG   =  6
  integer, parameter :: I_LRacLI2LG_I =  7 ! rain-ice(ice change)
  integer, parameter :: I_NRacNI2NG_I =  8
  integer, parameter :: I_LRacLI2LG_R =  9 ! rain-ice(rain change)
  integer, parameter :: I_NRacNI2NG_R = 10
  integer, parameter :: I_LRacLS2LG_S = 11 ! rain-snow(snow change)
  integer, parameter :: I_NRacNS2NG_S = 12
  integer, parameter :: I_LRacLS2LG_R = 13 ! rain-snow(rain change)
  integer, parameter :: I_NRacNS2NG_R = 14
  integer, parameter :: I_LRacLG2LG   = 15 ! rain-graupel(rain change)
  integer, parameter :: I_NRacNG2NG   = 16
  integer, parameter :: I_LIacLI2LS   = 17 ! ice-ice
  integer, parameter :: I_NIacNI2NS   = 18
  integer, parameter :: I_LIacLS2LS   = 19 ! ice-snow(ice change)
  integer, parameter :: I_NIacNS2NS   = 20
  integer, parameter :: I_NSacNS2NS   = 21 ! snow-snow
  integer, parameter :: I_NGacNG2NG   = 22 ! graupel-graupel
  integer, parameter :: I_LGacLS2LG   = 23 ! snow-graupel
  integer, parameter :: I_NGacNS2NG   = 24

  integer, parameter :: Pac_MAX       = 24



  character(len=H_SHORT), save :: WLABEL(HYDRO_MAX)

  ! empirical value from Meyers etal.(1991), 1[/liter] = 1.d3[/m3]
  real(RP), private, parameter :: nqmin(HYDRO_MAX) = (/ 1.E+4_RP, 1.0_RP, 1.0_RP, 1.E-4_RP, 1.E-4_RP /) ! [1/m3]
  ! refer to Seifert(2002) (Dr. Thesis, Table.5.1)
  ! max mass, for D_min=79um, 2mm, 5mm, 1cm, 1cm
  real(RP), private, parameter :: xqmax(HYDRO_MAX) = (/ 2.6E-10_RP, 5.0E-6_RP, 1.377E-6_RP, 7.519E-6_RP, 4.90E-5_RP /)! [kg]
  ! SB06, Table 1.
  ! min mass, for D_min=2um, 79um, 10um, 20um, 100um
  real(RP), private, parameter :: xqmin(HYDRO_MAX) = (/ 4.20E-15_RP, 2.60E-10_RP, 3.382E-13_RP, 1.847E-12_RP, 1.230E-10_RP /)! [kg]



  ! for all processes
  ! SB06, Table 1.
  real(RP), private, parameter :: xc_min = 4.20E-15_RP   ! [kg] : min mass, D_min=2um
  real(RP), private, parameter :: xr_min = 2.60E-10_RP   ! [kg] : min mass, D_min=79um
  real(RP), private, parameter :: xi_min = 3.382E-13_RP  ! [kg] : min mass, D_min=10um
  real(RP), private, parameter :: xs_min = 1.847E-12_RP  ! [kg] : min mass, D_min=20um
  real(RP), private, parameter :: xg_min = 1.230E-10_RP  ! [kg] : min mass, D_min=100um
  ! refer to Seifert(2002) (Dr. Thesis, Table.5.1)
  real(RP), private, parameter :: xc_max = 2.6E-10_RP    ! [kg] : max, D_max=79um
  real(RP), private, parameter :: xr_max = 5.00E-6_RP    ! [kg] : max, D_max=2mm
  real(RP), private, parameter :: xi_max = 1.377E-6_RP   ! [kg] : max, D_max=5mm
  real(RP), private, parameter :: xs_max = 7.519E-6_RP   ! [kg] : max, D_max=1cm
  real(RP), private, parameter :: xg_max = 4.900E-5_RP   ! [kg] : max, D_max=1cm
  ! filter similar to Ikawa et al.(1991) sec.3.5
  real(RP), private, parameter :: xmin_filter= xc_min
  ! filter of effective radius(1 micron)
  real(RP), private, parameter :: rmin_re= 1.E-6_RP
  !
  ! SB06(95),(96)
  real(RP), private, parameter :: n0r_min= 2.5E+5_RP    ! [m-4]: min intercept parameter of rain
  real(RP), private, parameter :: n0r_max= 2.0E+7_RP    ! [m-4]: max
  real(RP), private, parameter :: lambdar_min= 1.E+3_RP ! [m-1]: min slope parameter of rain
  real(RP), private, parameter :: lambdar_max= 1.E+4_RP ! [m-1]: max
  ! empirical value from Meyers etal.(1991), 1[/liter] = 1.d3[/m3]
  real(RP), private, parameter :: nc_min = 1.E+4_RP     ! [m-3] empirical T.Mitsui
  real(RP), private, parameter :: nr_min = 1.0_RP     ! [m-3] 1/1000 [/liter]
  real(RP), private, parameter :: ni_min = 1.0_RP     ! [m-3]
  real(RP), private, parameter :: ns_min = 1.E-4_RP    ! [m-3]
  real(RP), private, parameter :: ng_min = 1.E-4_RP    ! [m-3]
  ! empirical filter
  real(RP), private, parameter :: lc_min = xc_min*nc_min
  real(RP), private, parameter :: lr_min = xr_min*nr_min
  real(RP), private, parameter :: li_min = xi_min*ni_min
  real(RP), private, parameter :: ls_min = xs_min*ns_min
  real(RP), private, parameter :: lg_min = xg_min*ng_min
  !
  real(RP), private, parameter :: x_sep   = 2.6E-10_RP ! boundary mass between cloud and rain
  !
  real(RP), private, parameter :: tem_min=100.0_RP
  real(RP), private, parameter :: rho_min=1.E-5_RP     ! 3.e-3 is lower limit recognized in many experiments.
  real(RP), private, parameter :: rhoi   = 916.70_RP
  !
  integer, private, save :: ntmax_phase_change = 1
  integer, private, save :: ntmax_collection   = 1
  integer, private, save :: MP_ntmax_sedimentation= 1 ! 10/08/03 [Add] T.Mitsui
  !
  !--- standard density
  real(RP), private, parameter :: rho_0 = 1.280_RP
  !--- max number of Nc( activatable aerosol number concentration )
  real(RP), allocatable, private, save :: nc_uplim_d(:,:,:)
  !
  !--- thermal conductivity of air
  real(RP), private, parameter :: Ka0  = 2.428E-2_RP
  !<--- J/m/s/K : 0C/1atm
  real(RP), private, parameter :: dKa_dT = 7.47E-5_RP
  !<--- J/m/s/K/K : dependency of temperature
  !====== Ka = Ka0 + temc*dKa_dT
  !
  !--- Dynamic viscosity
  real(RP), private, parameter :: mua0 = 1.718E-5_RP
  !<--- Kg/m/s : 0C/1atm
  real(RP), private, parameter :: dmua_dT = 5.28E-8_RP
  !<--- Kg/m/s/K : dependency of temperature
  !======  mua = mua0 + temc*dmua_dT
  !
  real(RP), private, save :: xc_ccn = 1.E-12_RP ! [kg]
  real(RP), private, save :: xi_ccn = 1.E-12_RP ! [kg] ! [move] 11/08/30 T.Mitsui
  !
  ! capacity of diffusional growth
  ! ( dependent of their geometries )
  real(RP), private, save :: cap(HYDRO_MAX)
  !
  ! constants for Diameter-Mass relation
  ! D = a * x^b
  real(RP), private, save :: a_m(HYDRO_MAX)
  real(RP), private, save :: b_m(HYDRO_MAX)
  ! constants for Terminal velocity-Mass relation
  ! vt = alpha * x^beta * f
  real(RP), private, save :: alpha_v(HYDRO_MAX,2)
  real(RP), private, save :: beta_v(HYDRO_MAX,2)
  real(RP), private, save :: alpha_vn(HYDRO_MAX,2)   !
  real(RP), private, save :: beta_vn(HYDRO_MAX,2)    !
  real(RP), private, save :: gamma_v(HYDRO_MAX)
  !  Aerodynamical factor for correction of terminal velocity.(Heymsfield and Iaquinta, 2000)
  !  vt(tem,pre) = vt0 * (pre/pre0)**a_pre0 * (tem/tem0)**a_tem0
  real(RP), private, parameter :: pre0_vt   = 300.E+2_RP  ! 300hPa
  real(RP), private, parameter :: tem0_vt   = 233.0_RP  ! -40degC
  real(RP), private, parameter :: a_pre0_vt = -0.1780_RP
  real(RP), private, parameter :: a_tem0_vt = -0.3940_RP
  ! Parameters to determine Droplet Size Distribution
  ! as a General Gamma Distribution
  ! f(x) = A x^nu exp(-lambda x^mu )
  ! for Marshall Palmer Distribution ( popular for rain )
  !     mu=1/3, nu=-2/3
  ! for Gamma Distribution ( popular for cloud )
  !     mu=1
  real(RP), private, save :: nu(HYDRO_MAX)
  real(RP), private, save :: mu(HYDRO_MAX)
  ! Mitchell(1996), JAS, vol.53, No.12, pp.1710-1723
  !  area = a_area*D^b_area
  !  area = ax_area*x^bx_area
  ! Auer and Veal(1970), JAS, vol.27, pp.919-pp.926
  ! height = a_h*x^b_h( based on h=a_ar*D^b_ar,  ar:aspect ratio)
  real(RP), private, save :: a_area(HYDRO_MAX)       !
  real(RP), private, save :: b_area(HYDRO_MAX)       !
  real(RP), private, save :: ax_area(HYDRO_MAX)      !
  real(RP), private, save :: bx_area(HYDRO_MAX)      !
  ! parameters for radius of equivalent area
  ! r_ea = a_rea*x**b_rea
  real(RP), private, save :: a_rea(HYDRO_MAX)        !
  real(RP), private, save :: b_rea(HYDRO_MAX)        !
  real(RP), private, save :: a_rea2(HYDRO_MAX)       !
  real(RP), private, save :: b_rea2(HYDRO_MAX)       !
  real(RP), private, save :: a_rea3(HYDRO_MAX)       !
  real(RP), private, save :: b_rea3(HYDRO_MAX)       !
  !
  real(RP), private, save :: a_d2vt(HYDRO_MAX)       !
  real(RP), private, save :: b_d2vt(HYDRO_MAX)       !
  ! coefficient of x^2 moment of DSD
  ! Z = integral x*x*f(x) dx
  !   = coef_m2*N*(L/N)^2
  real(RP), private, save :: coef_m2(HYDRO_MAX)
  ! radar reflectivity coefficient defined by diameter
  real(RP), private, save :: coef_d6(HYDRO_MAX)      !
  ! volume coefficient defined by diameter
  real(RP), private, save :: coef_d3(HYDRO_MAX)      !
  ! coefficient of weighted mean diameter
  real(RP), private, save :: coef_d(HYDRO_MAX)
  ! coefficient of weighted mean d*d*v
  real(RP), private, save :: coef_d2v(HYDRO_MAX)     !
  ! coefficient of moment of d*d*v
  real(RP), private, save :: coef_md2v(HYDRO_MAX)    !
  !
  ! for effective radius(spherical particle)
  real(RP), private, save :: coef_r2(HYDRO_MAX)
  real(RP), private, save :: coef_r3(HYDRO_MAX)
  real(RP), private, save :: coef_re(HYDRO_MAX)
  ! for effective radius(hexagonal plate)
  real(RP), private, save :: coef_rea2(HYDRO_MAX)    !
  real(RP), private, save :: coef_rea3(HYDRO_MAX)    !
  logical, private, save :: opt_M96_ice=.true.               !
  logical, private, save :: opt_M96_column_ice=.false.       !
  !
  ! coefficeint of weighted mean terminal velocity
  ! vt0 is number weighted and
  ! vt1 is mass   weighted.
  real(RP), private, save :: coef_vt0(HYDRO_MAX,2)
  real(RP), private, save :: coef_vt1(HYDRO_MAX,2)
  real(RP), private, save :: coef_deplc
  real(RP), private, save :: coef_dave_N(HYDRO_MAX) !
  real(RP), private, save :: coef_dave_L(HYDRO_MAX) !
  ! diameter of terminal velocity branch
  !
  real(RP), private, save :: d0_ni=261.76E-6_RP   !
  real(RP), private, save :: d0_li=398.54E-6_RP
  real(RP), private, parameter     :: d0_ns=270.03E-6_RP   !
  real(RP), private, parameter     :: d0_ls=397.47E-6_RP   !
  real(RP), private, parameter     :: d0_ng=269.08E-6_RP   !
  real(RP), private, parameter     :: d0_lg=376.36E-6_RP   !
  !
  real(RP), private, parameter :: coef_vtr_ar1=9.65_RP    ! coef. for large branch
  ! original parameter of Rogers etal.(1993)
  real(RP), private, parameter :: coef_vtr_br1=10.43_RP    ! ...
  real(RP), private, parameter :: coef_vtr_cr1=600.0_RP    ! ...
  real(RP), private, parameter :: coef_vtr_ar2=4.E+3_RP     ! coef. for small branch
  real(RP), private, parameter :: coef_vtr_br2=12.E+3_RP     ! ...
  real(RP), private, parameter :: d_vtr_branch=0.745E-3_RP  ! 0.745 mm (diameter dividing 2-branches)
  ! equilibrium diameter of rain break-up
  real(RP), private, parameter :: dr_eq   = 1.10E-3_RP ! eqilibrium diameter, Seifert 2008(36)
  ! coefficient of General Gamma.
  !  f(x)  = A x^nu exp(-lambda x^mu )
  ! lambda = coef_lambda * (L/N)^{-mu}
  !     A  = coef_A*N*lambda^slope_A
  real(RP), private, save :: coef_A(HYDRO_MAX)
  real(RP), private, save :: coef_lambda(HYDRO_MAX)
!  real(RP), private, save :: slope_A(HYDRO_MAX)
  ! coefficeint of weighted ventilation effect.
  ! large, and small branch is by PK97(13-60),(13-61),(13-88),(13-89)
  real(RP), private, save :: ah_vent  (HYDRO_MAX,2) !
  real(RP), private, save :: bh_vent  (HYDRO_MAX,2) !
  real(RP), private, save :: ah_vent0 (HYDRO_MAX,2) !
  real(RP), private, save :: bh_vent0 (HYDRO_MAX,2) !
  real(RP), private, save :: ah_vent1 (HYDRO_MAX,2) !
  real(RP), private, save :: bh_vent1 (HYDRO_MAX,2) !
  ! coefficient of collision growth
  real(RP), private, save :: delta_b0 (HYDRO_MAX)
  real(RP), private, save :: delta_b1 (HYDRO_MAX)
  real(RP), private, save :: delta_ab0(HYDRO_MAX,HYDRO_MAX)
  real(RP), private, save :: delta_ab1(HYDRO_MAX,HYDRO_MAX)
  !
  real(RP), private, save :: theta_b0 (HYDRO_MAX)
  real(RP), private, save :: theta_b1 (HYDRO_MAX)
  real(RP), private, save :: theta_ab0(HYDRO_MAX,HYDRO_MAX)
  real(RP), private, save :: theta_ab1(HYDRO_MAX,HYDRO_MAX)
  !
  logical, private, save :: opt_debug=.false.
  !
  logical, private, save :: opt_debug_tem=.false.
  logical, private, save :: opt_debug_inc=.true.
  logical, private, save :: opt_debug_act=.true.
  logical, private, save :: opt_debug_ree=.true.
  logical, private, save :: opt_debug_bcs=.true.

  integer, private, save :: MP_NSTEP_SEDIMENTATION
  real(RP), private, save :: MP_RNSTEP_SEDIMENTATION
  real(DP), private, save :: MP_DTSEC_SEDIMENTATION

  !
  ! metrics of vertical coordinate
  ! not used in SCALE-RM
  !
  real(RP), private, allocatable, save :: gsgam2_d (:,:,:)
  real(RP), private, allocatable, save :: gsgam2h_d(:,:,:)
  real(RP), private, allocatable, save :: gam2_d   (:,:,:)
  real(RP), private, allocatable, save :: gam2h_d  (:,:,:)
  real(RP), private, allocatable, save :: rgsgam2_d(:,:,:)
  real(RP), private, allocatable, save :: rgs_d    (:,:,:)
  real(RP), private, allocatable, save :: rgsh_d   (:,:,:)

  logical, private, save :: MP_doautoconversion = .true.
  logical, private, save :: MP_doprecipitation  = .true.
  logical, private, save :: MP_couple_aerosol   = .false. ! apply CCN effect?
  real(RP), private, save :: MP_ssw_lim = 1.E+1_RP


  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Configure
  subroutine ATMOS_PHY_MP_sn14_config( &
       MP_TYPE, &
       QA, QS   )
    use scale_process, only: &
       PRC_MPIstop
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    character(len=*), intent(in)  :: MP_TYPE
    integer,          intent(out) :: QA
    integer,          intent(out) :: QS

    integer :: QS2
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Cloud Microphysics Tracer] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Tracers for Seiki and Nakajima (2014) 2-moment bulk 6 category'

    if ( MP_TYPE /= 'SN14' ) then
       write(*,*) 'xxx ATMOS_PHY_MP_TYPE is not SN14. Check!'
       call PRC_MPIstop
    end if

    call ATMOS_HYDROMETEOR_regist( QS_MP,                       & ! [OUT]
                                   1, 2, 3,                     & ! [IN]
                                   ATMOS_PHY_MP_sn14_NAME(1:6), & ! [IN]
                                   ATMOS_PHY_MP_sn14_DESC(1:6), & ! [IN]
                                   ATMOS_PHY_MP_sn14_UNIT(1:6)  ) ! [IN]

    call TRACER_regist( QS2,                          & ! [OUT]
                        5,                            & ! [IN]
                        ATMOS_PHY_MP_sn14_NAME(7:11), & ! [IN]
                        ATMOS_PHY_MP_sn14_DESC(7:11), & ! [IN]
                        ATMOS_PHY_MP_sn14_UNIT(7:11)  ) ! [IN]

    QA = QA_MP
    QS = QS_MP
    QE_MP = QS_MP + QA_MP - 1

    I_QV = QS
    I_QC = QS + I_mp_QC
    I_QR = QS + I_mp_QR
    I_QI = QS + I_mp_QI
    I_QS = QS + I_mp_QS
    I_QG = QS + I_mp_QG
    I_NC = QS + I_mp_NC
    I_NR = QS + I_mp_NR
    I_NI = QS + I_mp_NI
    I_NS = QS + I_mp_NS
    I_NG = QS + I_mp_NG

    return
  end subroutine ATMOS_PHY_MP_sn14_config

  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  subroutine ATMOS_PHY_MP_sn14_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_grid, only: &
       CDZ => GRID_CDZ
    use scale_const, only: &
       CONST_UNDEF, &
       CONST_DWATR, &
       CONST_DICE
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doautoconversion, &
       MP_doprecipitation,  &
       MP_ssw_lim,          &
       MP_couple_aerosol,   &
       MP_ntmax_sedimentation

    real(RP), parameter :: max_term_vel = 10.0_RP  !-- terminal velocity for calculate dt of sedimentation
    integer :: nstep_max
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Cloud Microphysics] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Seiki and Nakajima (2014) 2-moment bulk 6 category'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_ATMOS_PHY_MP)

    ATMOS_PHY_MP_sn14_DENS(:) = CONST_UNDEF
    ATMOS_PHY_MP_sn14_DENS(I_HC) = CONST_DWATR
    ATMOS_PHY_MP_sn14_DENS(I_HR) = CONST_DWATR
    ATMOS_PHY_MP_sn14_DENS(I_HI) = CONST_DICE
    ATMOS_PHY_MP_sn14_DENS(I_HS) = CONST_DICE
    ATMOS_PHY_MP_sn14_DENS(I_HG) = CONST_DICE

    WLABEL(1) = "CLOUD"
    WLABEL(2) = "RAIN"
    WLABEL(3) = "ICE"
    WLABEL(4) = "SNOW"
    WLABEL(5) = "GRAUPEL"

    call mp_sn14_init

    allocate(nc_uplim_d(1,IA,JA))
    nc_uplim_d(:,:,:) = 150.E6_RP

    nstep_max = int( ( TIME_DTSEC_ATMOS_PHY_MP * max_term_vel ) / minval( CDZ ) )
    MP_ntmax_sedimentation = max( MP_ntmax_sedimentation, nstep_max )

    MP_NSTEP_SEDIMENTATION  = MP_ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(MP_ntmax_sedimentation,kind=RP)
    MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Timestep of sedimentation is divided into : ', MP_ntmax_sedimentation, ' step(s)'
    if( IO_L ) write(IO_FID_LOG,*) '*** DT of sedimentation [s]                   : ', MP_DTSEC_SEDIMENTATION

    !--- For kij
    allocate( gsgam2_d (KA,IA,JA) )
    allocate( gsgam2h_d(KA,IA,JA) )
    allocate( gam2_d   (KA,IA,JA) )
    allocate( gam2h_d  (KA,IA,JA) )
    allocate( rgsgam2_d(KA,IA,JA) )
    allocate( rgs_d    (KA,IA,JA) )
    allocate( rgsh_d   (KA,IA,JA) )
    gsgam2_d (:,:,:) = 1.0_RP
    gsgam2h_d(:,:,:) = 1.0_RP
    gam2_d   (:,:,:) = 1.0_RP
    gam2h_d  (:,:,:) = 1.0_RP
    rgsgam2_d(:,:,:) = 1.0_RP
    rgs_d    (:,:,:) = 1.0_RP
    rgsh_d   (:,:,:) = 1.0_RP

    return
  end subroutine ATMOS_PHY_MP_sn14_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sn14( &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN,       &
       EVAPORATE, &
       SFLX_rain, &
       SFLX_snow  )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)    :: CCN(KA,IA,JA)
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Cloud microphysics(SN14)'

#ifdef PROFILE_FIPP
    call fipp_start()
#endif

    call MP_negativefilter( DENS, QTRC )

    call mp_sn14( DENS,      & ! [INOUT]
                  MOMZ,      & ! [INOUT]
                  MOMX,      & ! [INOUT]
                  MOMY,      & ! [INOUT]
                  RHOT,      & ! [INOUT]
                  QTRC,      & ! [INOUT]
                  CCN,       & ! [IN]
                  EVAPORATE, & ! [OUT]
                  SFLX_rain, & ! [OUT]
                  SFLX_snow  ) ! [OUT]

    call MP_negativefilter( DENS, QTRC )

#ifdef PROFILE_FIPP
    call fipp_stop()
#endif

    return
  end subroutine ATMOS_PHY_MP_sn14

  !-----------------------------------------------------------------------------
  subroutine mp_sn14_init
    use scale_process, only: &
       PRC_MPIstop
    use scale_specfunc, only: &
        gammafunc => SF_gamma
    implicit none

    real(RP) :: w1(HYDRO_MAX)
    real(RP) :: w2(HYDRO_MAX)
    real(RP) :: w3(HYDRO_MAX)
    real(RP) :: w4(HYDRO_MAX)
    real(RP) :: w5(HYDRO_MAX)
    real(RP) :: w6(HYDRO_MAX)
    real(RP) :: w7(HYDRO_MAX)
    real(RP) :: w8(HYDRO_MAX)

    ! work for calculation of capacity, Mitchell and Arnott (1994) , eq.(9)
    real(RP) :: ar_ice_fix = 0.7_RP
    real(RP) :: wcap1, wcap2
    ! work for ventilation coefficient
    logical :: flag_vent0(HYDRO_MAX), flag_vent1(HYDRO_MAX)
    integer :: ierr
    integer :: iw, ia, ib
    integer :: n
    !
    namelist /nm_mp_sn14_init/       &
         opt_debug,                  &
         opt_debug_tem,              &
         opt_debug_inc,              &
         opt_debug_act,              &
         opt_debug_ree,              &
         opt_debug_bcs,              &
         ntmax_phase_change,         &
         ntmax_collection
    !
    namelist /nm_mp_sn14_particles/ &
         a_m, b_m, alpha_v, beta_v, gamma_v, &
         alpha_vn, beta_vn,    &
         a_area, b_area, cap,  &
         nu, mu,               &
         opt_M96_column_ice,   &
         opt_M96_ice,          &
         ar_ice_fix
    real(RP), parameter :: eps_gamma=1.E-30_RP

    a_m(:)         = UNDEF8
    b_m(:)         = UNDEF8
    alpha_v(:,:)   = UNDEF8
    beta_v(:,:)    = UNDEF8
    alpha_vn(:,:)  = UNDEF8
    beta_vn(:,:)   = UNDEF8
    gamma_v(:)     = UNDEF8
    a_d2vt(:)      = UNDEF8
    b_d2vt(:)      = UNDEF8
    a_area(:)      = UNDEF8
    b_area(:)      = UNDEF8
    ax_area(:)     = UNDEF8
    bx_area(:)     = UNDEF8
    a_rea(:)       = UNDEF8
    b_rea(:)       = UNDEF8
    a_rea2(:)      = UNDEF8
    b_rea2(:)      = UNDEF8
    a_rea3(:)      = UNDEF8
    b_rea3(:)      = UNDEF8
    nu(:)          = UNDEF8
    mu(:)          = UNDEF8
    cap(:)         = UNDEF8
    coef_m2(:)     = UNDEF8
    coef_dave_N(:) = UNDEF8
    coef_dave_L(:) = UNDEF8
    coef_d(:)      = UNDEF8
    coef_d3(:)     = UNDEF8
    coef_d6(:)     = UNDEF8
    coef_d2v(:)    = UNDEF8
    coef_md2v(:)   = UNDEF8
    coef_r2(:)     = UNDEF8
    coef_r3(:)     = UNDEF8
    coef_re(:)     = UNDEF8
    coef_rea2(:)   = UNDEF8
    coef_rea3(:)   = UNDEF8
    coef_A(:)      = UNDEF8
!    slope_A(:)     = UNDEF8
    coef_lambda(:) = UNDEF8
    coef_vt0(:,:)  = UNDEF8
    coef_vt1(:,:)  = UNDEF8
    delta_b0(:)    = UNDEF8
    delta_b1(:)    = UNDEF8
    delta_ab0(:,:) = UNDEF8
    delta_ab1(:,:) = UNDEF8
    theta_b0(:)    = UNDEF8
    theta_b1(:)    = UNDEF8
    theta_ab0(:,:) = UNDEF8
    theta_ab1(:,:) = UNDEF8
    !
    ah_vent(:,:)   = UNDEF8
    ah_vent0(:,:)  = UNDEF8
    ah_vent1(:,:)  = UNDEF8
    bh_vent(:,:)   = UNDEF8
    bh_vent0(:,:)  = UNDEF8
    bh_vent1(:,:)  = UNDEF8

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_sn14_init,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist nm_mp_sn14_init. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=nm_mp_sn14_init)

    !
    ! default setting
    !
    ! Area parameters with mks unit originated by Mitchell(1996)
    a_area(I_mp_QC) = PI/4.0_RP ! sphere
    a_area(I_mp_QR) = PI/4.0_RP ! sphere
    a_area(I_mp_QI) = 0.65_RP*1.E-4_RP*100.0_RP**(2.00_RP)   ! Mitchell(1996), Hexagonal Plate
    a_area(I_mp_QS) = 0.2285_RP*1.E-4_RP*100.0_RP**(1.88_RP) ! Mitchell(1996), Aggregates
    a_area(I_mp_QG) = 0.50_RP*1.E-4_RP*100.0_RP**(2.0_RP)    ! Mitchell(1996), Lump Graupel
    b_area(I_mp_QC) = 2.0_RP
    b_area(I_mp_QR) = 2.0_RP
    b_area(I_mp_QI) = 2.0_RP
    b_area(I_mp_QS) = 1.88_RP
    b_area(I_mp_QG) = 2.0_RP
    !
    ! Seifert and Beheng(2006), Table. 1 or List of symbols
    !----------------------------------------------------------
    ! Diameter-Mass relationship
    ! D = a * x^b
    a_m(I_mp_QC) = 0.124_RP
    a_m(I_mp_QR) = 0.124_RP
    a_m(I_mp_QI) = 0.217_RP
    a_m(I_mp_QS) = 8.156_RP
    a_m(I_mp_QG) = 0.190_RP
    b_m(I_mp_QC) = 1.0_RP/3.0_RP
    b_m(I_mp_QR) = 1.0_RP/3.0_RP
    b_m(I_mp_QI) = 0.302_RP
    b_m(I_mp_QS) = 0.526_RP
    b_m(I_mp_QG) = 0.323_RP
    !----------------------------------------------------------
    ! Terminal velocity-Mass relationship
    ! vt = alpha * x^beta * (rho0/rho)^gamma
    alpha_v(I_mp_QC,:)= 3.75E+5_RP
    alpha_v(I_mp_QR,:)= 159.0_RP ! not for sedimantation
    alpha_v(I_mp_QI,:)= 317.0_RP
    alpha_v(I_mp_QS,:)= 27.70_RP
    alpha_v(I_mp_QG,:)= 40.0_RP
    beta_v(I_mp_QC,:) = 2.0_RP/3.0_RP
    beta_v(I_mp_QR,:) = 0.266_RP ! not for sedimantation
    beta_v(I_mp_QI,:) = 0.363_RP
    beta_v(I_mp_QS,:) = 0.216_RP
    beta_v(I_mp_QG,:) = 0.230_RP
    gamma_v(I_mp_QC)  = 1.0_RP
    ! This is high Reynolds number limit(Beard 1980)
    gamma_v(I_mp_QR)  = 1.0_RP/2.0_RP
    gamma_v(I_mp_QI)  = 1.0_RP/2.0_RP
    gamma_v(I_mp_QS)  = 1.0_RP/2.0_RP
    gamma_v(I_mp_QG)  = 1.0_RP/2.0_RP
    !----------------------------------------------------------
    ! DSD parameters
    ! f(x) = A x^nu exp( -lambda x^mu )
    ! Gamma Disribution           : mu=1  , nu:arbitrary
    ! Marshall-Palmer Distribution: mu=1/3, nu:-2/3
    ! In the case of MP, f(D) dD = f(x)dx
    !                    f(x)    = c * f(D)/D^2 (c:coefficient)
    nu(I_mp_QC) =  1.0_RP          ! arbitrary for Gamma
    nu(I_mp_QR) = -1.0_RP/3.0_RP     ! nu(diameter)=1, equilibrium condition.
    nu(I_mp_QI) =  1.0_RP          !
    nu(I_mp_QS) =  1.0_RP          !
    nu(I_mp_QG) =  1.0_RP          !
    !
    mu(I_mp_QC) = 1.0_RP           ! Gamma
    mu(I_mp_QR) = 1.0_RP/3.0_RP      ! Marshall Palmer
    mu(I_mp_QI) = 1.0_RP/3.0_RP      !
    mu(I_mp_QS) = 1.0_RP/3.0_RP      !
    mu(I_mp_QG) = 1.0_RP/3.0_RP      !
    !----------------------------------------------------------
    ! Geomeries for diffusion growth
    ! Pruppacher and Klett(1997), (13-77)-(13-80) and
    ! originally derived by McDonald(1963b)
    ! sphere: cap=2
    ! plate : cap=pi
    ! needle with aspect ratio a/b
    !       : cap=log(2*a/b)
    cap(I_mp_QC) = 2.0_RP    ! sphere
    cap(I_mp_QR) = 2.0_RP    ! sphere
    cap(I_mp_QI) = PI ! hexagonal plate
    cap(I_mp_QS) = 2.0_RP    ! mix aggregates
    cap(I_mp_QG) = 2.0_RP    ! lump
    !
    alpha_vn(:,:) = alpha_v(:,:)
    beta_vn(:,:) = beta_v(:,:)
    !------------------------------------------------------------------------
    !
    ! additional setting
    !

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_sn14_particles,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist nm_mp_sn14_particles. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=nm_mp_sn14_particles)

    ! [Add] 10/08/03 T.Mitsui
    ! particles shapes are
    if( opt_M96_ice ) then
       ! ice is randomly oriented Hexagonal plate (Auer and Veal 1970, Takano and Liou 1995, Mitchell 1996)
       ! snow is assemblages of planar polycrystals(Mitchell 1996)
       ! graupel is Lump graupel(R4b) is assumed(Mitchell 1996)
       a_area(I_mp_QI)   = 0.120284936_RP
       a_area(I_mp_QS)   = 0.131488_RP
       a_area(I_mp_QG)   = 0.5_RP
       b_area(I_mp_QI)   = 1.850000_RP
       b_area(I_mp_QS)   = 1.880000_RP
       b_area(I_mp_QG)   = 2.0_RP
       a_m(I_mp_QI)      = 1.23655360084766_RP
       a_m(I_mp_QS)      = a_m(I_mp_QI)
       a_m(I_mp_QG)      = 0.346111225718402_RP
       b_m(I_mp_QI)      = 0.408329930583912_RP
       b_m(I_mp_QS)      = b_m(I_mp_QI)
       b_m(I_mp_QG)      = 0.357142857142857_RP
       !
       if( opt_M96_column_ice )then
          d0_ni=240.49E-6_RP   ! this is column
          d0_li=330.09E-6_RP   ! this is column
          a_area(I_mp_QI)= (0.684_RP*1.E-4_RP)*10.0_RP**(2.0_RP*2.00_RP)
          b_area(I_mp_QI)= 2.0_RP
          a_m(I_mp_QI)   = 0.19834046116844_RP
          b_m(I_mp_QI)   = 0.343642611683849_RP
          ! [Add] 11/08/30 T.Mitsui
          ! approximated by the capacity for prolate spheroid with constant aspect ratio
          wcap1        = sqrt(1.0_RP-ar_ice_fix**2)
          wcap2        = log( (1.0_RP+wcap1)/ar_ice_fix )
          cap(I_mp_QI) = 2.0_RP*wcap2/wcap1
          !
       end if
       !
       ! These value are derived by least-square fitting in the range
       ! qi [100um:1000um] in diameter
       ! qs [100um:1000um] in diameter
       ! qg [200um:2000um] in diameter
       !                      small branch          , large branch
       alpha_v (I_mp_QI,:) =(/ 5798.60107421875_RP, 167.347076416016_RP/)
       alpha_vn(I_mp_QI,:) =(/ 12408.177734375_RP, 421.799865722656_RP/)
       if( opt_M96_column_ice )then
          alpha_v (I_mp_QI,:) = (/2901.0_RP, 32.20_RP/)
          alpha_vn(I_mp_QI,:) = (/9675.2_RP, 64.16_RP/)
       end if
       alpha_v (I_mp_QS,:)    =(/ 15173.3916015625_RP, 305.678619384766_RP/)
       alpha_vn(I_mp_QS,:)    =(/ 29257.1601562500_RP, 817.985717773438_RP/)
       alpha_v (I_mp_QG,:)    =(/ 15481.6904296875_RP, 311.642242431641_RP/)
       alpha_vn(I_mp_QG,:)    =(/ 27574.6562500000_RP, 697.536132812500_RP/)
       !
       beta_v (I_mp_QI,:)  =(/ 0.504873454570770_RP, 0.324817866086960_RP/)
       beta_vn(I_mp_QI,:)  =(/ 0.548495233058929_RP, 0.385287821292877_RP/)
       if( opt_M96_column_ice )then
          beta_v (I_mp_QI,:)  =(/ 0.465552181005478_RP, 0.223826110363007_RP/)
          beta_vn(I_mp_QI,:)  =(/ 0.530453503131866_RP, 0.273761242628098_RP/)
       end if
       beta_v (I_mp_QS,:)     =(/ 0.528109610080719_RP, 0.329863965511322_RP/)
       beta_vn(I_mp_QS,:)     =(/ 0.567154467105865_RP, 0.393876969814301_RP/)
       beta_v (I_mp_QG,:)     =(/ 0.534656763076782_RP, 0.330253750085831_RP/)
       beta_vn(I_mp_QG,:)     =(/ 0.570551633834839_RP, 0.387124240398407_RP/)
    end if
    !
    ! area-diameter relation => area-mass relation
    ax_area(:) = a_area(:)*a_m(:)**b_area(:)
    bx_area(:) = b_area(:)*b_m(:)
    !
    ! radius of equivalent area - m ass relation
    ! pi*rea**2 = ax_area*x**bx_area
    a_rea(:)   = sqrt(ax_area(:)/PI)
    b_rea(:)   = bx_area(:)/2.0_RP
    a_rea2(:)  = a_rea(:)**2
    b_rea2(:)  = b_rea(:)*2.0_RP
    a_rea3(:)  = a_rea(:)**3
    b_rea3(:)  = b_rea(:)*3.0_RP
    !
    a_d2vt(:)=alpha_v(:,2)*(1.0_RP/alpha_v(:,2))**(beta_v(:,2)/b_m(:))
    b_d2vt(:)=(beta_v(:,2)/b_m(:))
    !
    ! Calculation of Moment Coefficient
    !
    w1(:) = 0.0_RP
    w2(:) = 0.0_RP
    w3(:) = 0.0_RP
    w4(:) = 0.0_RP
    w5(:) = 0.0_RP
    w6(:) = 0.0_RP
    w7(:) = 0.0_RP
    w8(:) = 0.0_RP
    !-------------------------------------------------------
    ! moment coefficient
    ! SB06 (82)
    ! M^n /= coef_mn * N * (L/N)**n
    ! M^2 = Z = coef_m2 * N *(L/N)**2
    ! a*M^b = a*integral x^b f(x) dx = ave D
    do iw=1, HYDRO_MAX
       n = 2
       w1(iw) = gammafunc( (n+nu(iw)+1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw)+1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw)+2.0_RP)/mu(iw) )
       coef_m2(iw) = w1(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**n
       !
       w4(iw) = gammafunc( (b_m(iw)+nu(iw)+1.0_RP)/mu(iw) )
       coef_d(iw) = a_m(iw) * w4(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**b_m(iw)
       w5(iw) = gammafunc( (2.0_RP*b_m(iw)+beta_v(iw,2)+nu(iw)+1.0_RP)/mu(iw) )
       w6(iw) = gammafunc( (3.0_RP*b_m(iw)+beta_v(iw,2)+nu(iw)+1.0_RP)/mu(iw) )
       coef_d2v(iw) = a_m(iw) * w6(iw)/w5(iw)* ( w2(iw)/w3(iw) )**b_m(iw)
       coef_md2v(iw)=           w5(iw)/w2(iw)* ( w2(iw)/w3(iw) )**(2.0_RP*b_m(iw)+beta_v(iw,2))
       ! 09/04/14 [Add] T.Mitsui, volume and radar reflectivity
       w7(iw) = gammafunc( (3.0_RP*b_m(iw)+nu(iw)+1.0_RP)/mu(iw) )
       coef_d3(iw)  = a_m(iw)**3 * w7(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**(3.0_RP*b_m(iw))
       w8(iw) = gammafunc( (6.0_RP*b_m(iw)+nu(iw)+1.0_RP)/mu(iw) )
       coef_d6(iw)  = a_m(iw)**6 * w8(iw)/w2(iw)*(  w2(iw)/w3(iw)  )**(6.0_RP*b_m(iw))
    end do
    !
    coef_deplc = coef_d(I_mp_QC)/a_m(I_mp_QC)
    !-------------------------------------------------------
    ! coefficient of 2nd and 3rd moments for effective radius
    ! for spherical particle
    do iw=1, HYDRO_MAX
       ! integ r^2 f(x)dx
       w1(iw) = gammafunc( (2.0_RP*b_m(iw)+nu(iw)+1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw)+1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw)+2.0_RP)/mu(iw) )
       ! integ r^3 f(x)dx
       w4(iw) = gammafunc( (3.0_RP*b_m(iw)+nu(iw)+1.0_RP)/mu(iw) )
       !
       coef_r2(iw)=w1(iw)/w2(iw)*( w2(iw)/w3(iw) )**(2.0_RP*b_m(iw))
       coef_r3(iw)=w4(iw)/w2(iw)*( w2(iw)/w3(iw) )**(3.0_RP*b_m(iw))
       coef_re(iw)=coef_r3(iw)/coef_r2(iw)
       !
    end do
    !-------------------------------------------------------
    ! coefficient for effective radius of equivalent area and
    ! coefficient for volume of equivalent area
    do iw=1, HYDRO_MAX
       w1(iw) = gammafunc( (nu(iw)+1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw)+2.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (b_rea2(iw)+nu(iw)+1.0_RP)/mu(iw) )
       w4(iw) = gammafunc( (b_rea3(iw)+nu(iw)+1.0_RP)/mu(iw) )
       !
       coef_rea2(iw) = w3(iw)/w1(iw)*( w1(iw)/w2(iw) )**b_rea2(iw)
       coef_rea3(iw) = w4(iw)/w1(iw)*( w1(iw)/w2(iw) )**b_rea3(iw)
    end do
    !-------------------------------------------------------
    ! coefficient of gamma-distribution
    ! SB06(80)
    do iw=1, HYDRO_MAX
       w1(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       coef_A(iw)      = mu(iw)/w1(iw)
!       slope_A(iw)     = w1(iw)
       coef_lambda(iw) = (w1(iw)/w2(iw))**(-mu(iw))
    end do
    !-------------------------------------------------------
    ! coefficient for terminal velocity in sedimentation
    ! SB06(78)
    do ia=1,2
       do iw=1, HYDRO_MAX
          n = 0
          w1(iw) = gammafunc( (beta_vn(iw,ia) + nu(iw) + 1.0_RP + n)/mu(iw) )
          w2(iw) = gammafunc( (                 nu(iw) + 1.0_RP + n)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
          w4(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
          ! coefficient of terminal velocity for number
          coef_vt0(iw,ia)=alpha_vn(iw,ia)*w1(iw)/w2(iw)*(w3(iw)/w4(iw))**beta_vn(iw,ia)
          n = 1
          w1(iw) = gammafunc( (beta_v(iw,ia) + nu(iw) + 1.0_RP + n)/mu(iw) )
          w2(iw) = gammafunc( (                nu(iw) + 1.0_RP + n)/mu(iw) )
          ! coefficient of terminal velocity for mass
          coef_vt1(iw,ia)=alpha_v(iw,ia)*w1(iw)/w2(iw)*(w3(iw)/w4(iw))**beta_v(iw,ia)
       end do
    end do
    ! coefficient for weighted diameter used in calculation of terminal velocity
    do iw=1, HYDRO_MAX
       w1(iw) = gammafunc( (       b_m(iw) + nu(iw) + 1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (1.0_RP + b_m(iw) + nu(iw) + 1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w4(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       coef_dave_N(iw) =  (w1(iw)/w3(iw))*(w3(iw)/w4(iw))**        b_m(iw)
       coef_dave_L(iw) =  (w2(iw)/w3(iw))*(w3(iw)/w4(iw))**(1.0_RP+b_m(iw))
    end do
    !-------------------------------------------------------
    !
    ah_vent(I_mp_QC,1:2) = (/1.0000_RP,1.0000_RP/) ! no effect
    ah_vent(I_mp_QR,1:2) = (/1.0000_RP,0.780_RP/)
    ah_vent(I_mp_QI,1:2) = (/1.0000_RP,0.860_RP/)
    ah_vent(I_mp_QS,1:2) = (/1.0000_RP,0.780_RP/)
    ah_vent(I_mp_QG,1:2) = (/1.0000_RP,0.780_RP/)
    bh_vent(I_mp_QC,1:2) = (/0.0000_RP,0.0000_RP/)
    bh_vent(I_mp_QR,1:2) = (/0.108_RP,0.308_RP/)
    bh_vent(I_mp_QI,1:2) = (/0.140_RP,0.280_RP/)
    bh_vent(I_mp_QS,1:2) = (/0.108_RP,0.308_RP/)
    bh_vent(I_mp_QG,1:2) = (/0.108_RP,0.308_RP/)
    !
    do iw=1, HYDRO_MAX
       n = 0
       if( (nu(iw) + b_m(iw) + n) > eps_gamma  )then
          w1(iw) = gammafunc( (nu(iw) + b_m(iw) + n)/mu(iw) )
          w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
          ah_vent0(iw,1)= ah_vent(iw,1)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.0_RP)
          ah_vent0(iw,2)= ah_vent(iw,2)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.0_RP)
          flag_vent0(iw)=.true.
       else
          ah_vent0(iw,1)= 1.0_RP
          ah_vent0(iw,2)= 1.0_RP
          flag_vent0(iw)=.false.
       end if
       n = 1
       if( (nu(iw) + b_m(iw) + n) > eps_gamma  )then
          w1(iw) = gammafunc( (nu(iw) + b_m(iw) + n)/mu(iw) )
          w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
          ah_vent1(iw,1)= ah_vent(iw,1)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.0_RP)
          ah_vent1(iw,2)= ah_vent(iw,2)*(w1(iw)/w2(iw))*(w2(iw)/w3(iw))**(b_m(iw)+n-1.0_RP)
          flag_vent1(iw)=.true.
       else
          ah_vent1(iw,1)= 1.0_RP
          ah_vent1(iw,2)= 1.0_RP
          flag_vent1(iw)=.true.
       end if
    end do
    do iw=1, HYDRO_MAX
       n = 0
       if( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,1) + n) < eps_gamma )then
          flag_vent0(iw)=.false.
       end if
       if(flag_vent0(iw))then
          w1(iw) = gammafunc( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,1) + n)/mu(iw) )
          w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
          ! [Add] 11/08/30 T.Mitsui
          w4(iw) = gammafunc( (nu(iw) + 2.0_RP*b_m(iw) + beta_v(iw,1) + n)/mu(iw) )
          bh_vent0(iw,1)=bh_vent(iw,1)*(w4(iw)/w2(iw))*(w2(iw)/w3(iw))**(2.00_RP*b_m(iw)+beta_v(iw,1)+n-1.0_RP)
          w5(iw) = gammafunc( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,2) + n)/mu(iw) )
          bh_vent0(iw,2)=bh_vent(iw,2)*(w5(iw)/w2(iw))*(w2(iw)/w3(iw))**(1.5_RP*b_m(iw)+0.5_RP*beta_v(iw,2)+n-1.0_RP)
       else
          bh_vent0(iw,1) = 0.0_RP
          bh_vent0(iw,2) = 0.0_RP
       end if
       !
       n = 1
       if( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,1) + n) < eps_gamma )then
          flag_vent1(iw)=.false.
       end if
       if(flag_vent1(iw))then
          w1(iw) = gammafunc( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,1) + n)/mu(iw) )
          w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
          w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
          ! [Add] 11/08/30 T.Mitsui
          w4(iw) = gammafunc( (nu(iw) + 2.0_RP*b_m(iw) + beta_v(iw,1) + n)/mu(iw) )
          bh_vent1(iw,1)=bh_vent(iw,1)*(w4(iw)/w2(iw))*(w2(iw)/w3(iw))**(2.00_RP*b_m(iw)+beta_v(iw,1)+n-1.0_RP)
          !
          w5(iw) = gammafunc( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,2) + n)/mu(iw) )
          bh_vent1(iw,2)=bh_vent(iw,2)*(w5(iw)/w2(iw))*(w2(iw)/w3(iw))**(1.5_RP*b_m(iw)+0.5_RP*beta_v(iw,2)+n-1.0_RP)
       else
          bh_vent1(iw,1) = 0.0_RP
          bh_vent1(iw,2) = 0.0_RP
       end if
    end do
    !-------------------------------------------------------
    ! coefficient for collision process
    ! stochastic coefficient for collision cross section
    ! sb06 (90) -- self collection
    do iw=1, HYDRO_MAX
       n = 0
       w1(iw) = gammafunc( (2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       delta_b0(iw) = w1(iw)/w2(iw) &
            *( w2(iw)/w3(iw) )**(2.0_RP*b_rea(iw) + n)
       n = 1
       w1(iw) = gammafunc( (2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       delta_b1(iw) = w1(iw)/w2(iw) &
            *( w2(iw)/w3(iw) )**(2.0_RP*b_rea(iw) + n)
    end do
    ! stochastic coefficient for collision cross section
    ! sb06(91) -- riming( collide with others )
    do iw=1, HYDRO_MAX
       n = 0
       w1(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       w4(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP    )/mu(iw) )
       n = 1
       w5(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
    end do
    ! ia > ib ( larger particles "a" catch smaller particles "b" )
    do ia=1, HYDRO_MAX
       do ib=1, HYDRO_MAX
          n=0 !
          ! NOTE, collected  particle has a moment of n.
          !       collecting particle has only number(n=0).
          delta_ab0(ia,ib) = 2.0_RP*(w1(ib)/w2(ib))*(w4(ia)/w2(ia)) &
               * ( w2(ib)/w3(ib) )**(b_rea(ib)+n) &
               * ( w2(ia)/w3(ia) )**(b_rea(ia)  )
          n=1 !
          delta_ab1(ia,ib) = 2.0_RP*(w5(ib)/w2(ib))*(w4(ia)/w2(ia)) &
               * ( w2(ib)/w3(ib) )**(b_rea(ib)+n) &
               * ( w2(ia)/w3(ia) )**(b_rea(ia)  )
       end do
    end do
    ! stochastic coefficient for terminal velocity
    ! sb06(92) -- self collection
    ! assuming equivalent area circle.
    do iw=1, HYDRO_MAX
       n = 0
       w1(iw) = gammafunc( (2.0_RP*beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w2(iw) = gammafunc( (                      2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w4(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       theta_b0(iw) = w1(iw)/w2(iw) * ( w3(iw)/w4(iw) )**(2.0_RP*beta_v(iw,2))
       n = 1
       w1(iw) = gammafunc( (2.0_RP*beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w2(iw) = gammafunc( (                        2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       theta_b1(iw) = w1(iw)/w2(iw) * ( w3(iw)/w4(iw) )**(2.0_RP*beta_v(iw,2))
    end do
    !
    ! stochastic coefficient for terminal velocity
    ! sb06(93) -- riming( collide with others )
    do iw=1, HYDRO_MAX
       n = 0
       w1(iw) = gammafunc( (beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w2(iw) = gammafunc( (                   2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w3(iw) = gammafunc( (beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP    )/mu(iw) )
       w4(iw) = gammafunc( (                   2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP    )/mu(iw) )
       !
       w5(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w6(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       n = 1
       w7(iw) = gammafunc( (beta_v(iw,2) + b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w8(iw) = gammafunc( (                   b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
    end do
    ! ia > ib ( larger particles "a" catch smaller particles "b" )
    do ia=1, HYDRO_MAX
       do ib=1, HYDRO_MAX
          theta_ab0(ia,ib) = 2.0_RP * (w1(ib)/w2(ib))*(w3(ia)/w4(ia)) &
               * (w5(ia)/w6(ia))**beta_v(ia,2) &
               * (w5(ib)/w6(ib))**beta_v(ib,2)
          theta_ab1(ia,ib) = 2.0_RP * (w7(ib)/w8(ib))*(w3(ia)/w4(ia)) &
               * (w5(ia)/w6(ia))**beta_v(ia,2) &
               * (w5(ib)/w6(ib))**beta_v(ib,2)
       end do
    end do

    if( IO_L ) write(IO_FID_LOG,'(100a16)')      "LABEL       ",WLABEL(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "capacity    ",cap(:) ! [Add] 11/08/30 T.Mitsui
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_m2     ",coef_m2(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_d      ",coef_d(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_d3     ",coef_d3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_d6     ",coef_d6(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_d2v    ",coef_d2v(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_md2v   ",coef_md2v(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "a_d2vt      ",a_d2vt(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "b_d2vt      ",b_d2vt(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_r2     ",coef_r2(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_r3     ",coef_r3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_re     ",coef_re(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "a_area      ",a_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "b_area      ",b_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "ax_area     ",ax_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "bx_area     ",bx_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "a_rea       ",a_rea(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "b_rea       ",b_rea(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "a_rea3      ",a_rea3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "b_rea3      ",b_rea3(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_rea2   ",coef_rea2(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_rea3   ",coef_rea3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_vt0    ",coef_vt0(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_vt1    ",coef_vt1(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_A      ",coef_A(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "coef_lambda ",coef_lambda(:)

    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "ah_vent0 sml",ah_vent0(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "ah_vent0 lrg",ah_vent0(:,2)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "ah_vent1 sml",ah_vent1(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "ah_vent1 lrg",ah_vent1(:,2)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "bh_vent0 sml",bh_vent0(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "bh_vent0 lrg",bh_vent0(:,2)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "bh_vent1 sml",bh_vent1(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "bh_vent1 lrg",bh_vent1(:,2)

    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "delta_b0    ",delta_b0(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "delta_b1    ",delta_b1(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "theta_b0    ",theta_b0(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100ES16.6)') "theta_b1    ",theta_b1(:)

    do ia=1, HYDRO_MAX
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100ES16.6)') "delta0(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab0(ia,ib),ib=1,HYDRO_MAX)
    enddo
    do ia=1, HYDRO_MAX
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100ES16.6)') "delta1(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab1(ia,ib),ib=1,HYDRO_MAX)
    enddo
    do ia=1, HYDRO_MAX
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100ES16.6)') "theta0(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab0(ia,ib),ib=1,HYDRO_MAX)
    enddo
    do ia=1, HYDRO_MAX
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100ES16.6)') "theta1(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab1(ia,ib),ib=1,HYDRO_MAX)
    enddo

    return
  end subroutine mp_sn14_init
  !-----------------------------------------------------------------------------
  subroutine mp_sn14 ( &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN,       &
       EVAPORATE, &
       SFLX_rain, &
       SFLX_snow  )
    use scale_time, only: &
       dt_DP => TIME_DTSEC_ATMOS_PHY_MP
    use scale_grid, only: &
       z    => GRID_CZ, &
       dz   => GRID_CDZ
    use scale_atmos_phy_mp_common, only: &
         MP_precipitation => ATMOS_PHY_MP_precipitation
    use scale_tracer, only: &
       QA, &
       TRACER_R, &
       TRACER_CV, &
       TRACER_MASS
    use scale_atmos_saturation, only: &
       moist_psat_liq      => ATMOS_SATURATION_psat_liq,   &
       moist_psat_ice      => ATMOS_SATURATION_psat_ice
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)    :: CCN(KA,IA,JA)
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)

    !
    ! primary variables
    !
    real(RP) :: rhoq(I_QV:I_NG,KA,IA,JA)
    real(RP) :: pott          (KA,IA,JA)
    !
    ! diagnostic variables
    !
    real(RP) :: qdry(KA,IA,JA)
    real(RP) :: temp(KA,IA,JA)
    real(RP) :: pres(KA,IA,JA)
    real(RP) :: rhoe(KA,IA,JA)
    real(RP) :: velz(KA,IA,JA)
    !
    real(RP) :: rhoq2(I_QV:I_NG,KA,IA,JA)
    !
    real(RP) :: xq(HYDRO_MAX,KA,IA,JA)
    !
    real(RP) :: dq_xa(HYDRO_MAX,KA,IA,JA)
    real(RP) :: vt_xa(HYDRO_MAX,2,KA,IA,JA) ! terminal velocity

    real(RP) :: wtemp(KA,IA,JA) ! filtered temperature
    real(RP) :: esw(KA,IA,JA) ! saturated vapor pressure(water)
    real(RP) :: esi(KA,IA,JA) ! saturated vapor pressure(ice)
    !
    real(RP) :: rho_fac
    real(RP) :: rho_fac_q(HYDRO_MAX,KA,IA,JA) ! factor for tracers, 1:cloud, 2:rain, 3:ice, 4: snow, 5:graupel
    real(RP) :: cva(KA,IA,JA)    !
    real(RP) :: cpa(KA,IA,JA)       ! [Add] 09/08/18 T.Mitsui
    !
    real(RP) :: drhogqv               ! d (rho*qv*gsgam2)
    real(RP) :: drhogqc, drhognc      !        qc, nc
    real(RP) :: drhogqr, drhognr      !        qr, nr
    real(RP) :: drhogqi, drhogni      !        qi, ni
    real(RP) :: drhogqs, drhogns      !        qs, ns
    real(RP) :: drhogqg, drhogng      !        qg, ng

    ! production rate
    real(RP) :: PQ(PQ_MAX,KA,IA,JA)

    real(RP) :: wrm_dqc, wrm_dnc
    real(RP) :: wrm_dqr, wrm_dnr

    ! production rate of mixed-phase collection process
    real(RP) :: Pac(Pac_MAX,KA,IA,JA)

    real(RP) :: gc_dqc, gc_dnc
    real(RP) :: sc_dqc, sc_dnc
    real(RP) :: ic_dqc, ic_dnc
    real(RP) :: rg_dqg, rg_dng
    real(RP) :: rg_dqr, rg_dnr
    real(RP) :: rs_dqr, rs_dnr, rs_dqs, rs_dns
    real(RP) :: ri_dqr, ri_dnr
    real(RP) :: ri_dqi, ri_dni
    real(RP) :: ii_dqi, ii_dni
    real(RP) :: is_dqi, is_dni, ss_dns
    real(RP) :: gs_dqs, gs_dns, gg_dng
    ! mixed-phase collection process total plus(clp_), total minus(clm_)
    real(RP) :: clp_dqc, clp_dnc, clm_dqc, clm_dnc
    real(RP) :: clp_dqr, clp_dnr, clm_dqr, clm_dnr
    real(RP) :: clp_dqi, clp_dni, clm_dqi, clm_dni
    real(RP) :: clp_dqs, clp_dns, clm_dqs, clm_dns
    real(RP) :: clp_dqg, clp_dng, clm_dqg, clm_dng
    real(RP) :: fac1, fac3, fac4, fac6, fac7, fac9, fac10
    ! production rate of partial conversion(ice, snow => graupel)
    real(RP) :: pco_dqi, pco_dni
    real(RP) :: pco_dqs, pco_dns
    real(RP) :: pco_dqg, pco_dng
    ! production rate of enhanced melting due to
    real(RP) :: eml_dqc, eml_dnc
    real(RP) :: eml_dqr, eml_dnr
    real(RP) :: eml_dqi, eml_dni
    real(RP) :: eml_dqs, eml_dns
    real(RP) :: eml_dqg, eml_dng
    ! production rate of ice multiplication by splintering
    real(RP) :: spl_dqi, spl_dni
    real(RP) :: spl_dqg, spl_dqs

    real(RP) :: rrho(KA,IA,JA)

    !-----------------------------------------------
    ! work for explicit supersaturation modeling
    !-----------------------------------------------
    real(RP) :: dTdt_equiv_d(KA,IA,JA)    !
    !--------------------------------------------------
    !
    ! variables for output
    !
    !--------------------------------------------------
    ! work for column production term
    real(RP) :: sl_PLCdep(IA,JA)
    real(RP) :: sl_PLRdep(IA,JA), sl_PNRdep(IA,JA) !
    !--------------------------------------------------
    real(RP) :: qke_d(KA,IA,JA)

    real(RP), parameter :: eps       = 1.E-19_RP
    real(RP), parameter :: eps_qv    = 1.E-19_RP
    real(RP), parameter :: eps_rhoge = 1.E-19_RP
    real(RP), parameter :: eps_rhog  = 1.E-19_RP
    integer :: ntdiv

    real(RP) :: Rmoist

    real(RP) :: velw(KA,IA,JA,QA_MP-1)
    real(RP) :: FLX_rain (KA,IA,JA)
    real(RP) :: FLX_snow (KA,IA,JA)
    real(RP) :: FLX_tot  (KA,IA,JA)
    real(RP) :: wflux_rain(KA,IA,JA)
    real(RP) :: wflux_snow(KA,IA,JA)
    integer  :: step

    real(RP) :: sw
    real(RP) :: dt

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    !============================================================================
    !
    !--  Each process is integrated sequentially.
    !    1. Nucleation and filter
    !    2. Phase change
    !    3. Collection
    !    4. Saturation adjustment( only for qc, nc )
    !    5. Sedimentation
    !    6. filter( keep non-negative value for radiation scheme )
    !    0. calculation of optical moments
    !
    !============================================================================

    dt = real(dt_DP,kind=RP)

    !----------------------------------------------------------------------------
    !
    ! 1.Nucleation of cloud water and cloud ice
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP_Preprocess', 3)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       do iq = I_QV, I_NG
          rhoq(iq,k,i,j) = DENS(k,i,j) * QTRC(k,i,j,iq)
       enddo
       rhoq2(I_QV,k,i,j) = DENS(k,i,j)*QTRC(k,i,j,I_QV)
       rhoq2(I_NI,k,i,j) = max( 0.0_RP, DENS(k,i,j)*QTRC(k,i,j,I_NI) )
       rhoq2(I_NC,k,i,j) = max( 0.0_RP, DENS(k,i,j)*QTRC(k,i,j,I_NC) )
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       velz(KS-1,i,j) = 0.0_RP
       do k = KS, KE-1
          velz(k,i,j) = MOMZ(k,i,j) / ( DENS(k,i,j) + DENS(k+1,i,j) ) * 2.0_RP
       enddo
       velz(KE,i,j) = 0.0_RP
    end do
    end do

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rrho(k,i,j) = 1.0_RP / DENS(k,i,j)
       pott(k,i,j) = RHOT(k,i,j) * rrho(k,i,j)
       CALC_QDRY( qdry(k,i,j), QTRC, TRACER_MASS, k, i, j, iq )
       CALC_CV( cva(k,i,j), qdry(k,i,j), QTRC, k, i, j, iq, CVdry, TRACER_CV )
       CALC_R( Rmoist, qdry(k,i,j), QTRC, k, i, j, iq, Rdry, TRACER_R )
       cpa(k,i,j) = cva(k,i,j) + Rmoist
       CALC_PRE( pres(k,i,j), DENS(k,i,j), pott(k,i,j), Rmoist, cpa(k,i,j), P00 )
       temp(k,i,j) = pres(k,i,j) / ( DENS(k,i,j) * Rmoist )
       rhoe(k,i,j) = DENS(k,i,j) * temp(k,i,j) * cva(k,i,j)
       wtemp(k,i,j) = max(temp(k,i,j), tem_min)
    enddo
    enddo
    enddo

    if( opt_debug_tem ) call debug_tem_kij( 1, temp(:,:,:), DENS(:,:,:), pres(:,:,:), QTRC(:,:,:,I_QV) )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rho_fac              = rho_0 / max(DENS(k,i,j),rho_min)
       rho_fac_q(I_mp_QC,k,i,j) = rho_fac**gamma_v(I_mp_QC)
       rho_fac_q(I_mp_QR,k,i,j) = rho_fac**gamma_v(I_mp_QR)
       rho_fac_q(I_mp_QI,k,i,j) = (pres(k,i,j)/pre0_vt)**a_pre0_vt * (temp(k,i,j)/tem0_vt)**a_tem0_vt
       rho_fac_q(I_mp_QS,k,i,j) = rho_fac_q(I_mp_QI,k,i,j)
       rho_fac_q(I_mp_QG,k,i,j) = rho_fac_q(I_mp_QI,k,i,j)
    enddo
    enddo
    enddo

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
       sl_PLCdep(i,j) = 0.0_RP
       sl_PLRdep(i,j) = 0.0_RP
       sl_PNRdep(i,j) = 0.0_RP
    end do
    end do

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       qke_d(k,i,j) = 0.0_RP ! 2*TKE
    enddo
    enddo
    enddo

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dTdt_equiv_d(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    call PROF_rapend  ('MP_Preprocess', 3)

    call PROF_rapstart('MP_Nucleation', 3)

    call nucleation_kij(    &
         z, velz,           & ! in
         DENS, wtemp, pres, & ! in
         rhoq2,             & ! (in)
         PQ,                & ! out
         cpa,               & ! in
         dTdt_equiv_d,      & ! in
         qke_d,             & ! in
         CCN,               & ! in
         dt                 ) ! in

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! nucleation
       drhogqc = dt * PQ(I_LCccn,k,i,j)
       drhognc = dt * PQ(I_NCccn,k,i,j)
       drhogqi = dt * PQ(I_LIccn,k,i,j)
       drhogni = dt * PQ(I_NIccn,k,i,j)
       drhogqv = max( -rhoq(I_QV,k,i,j), -drhogqc-drhogqi )
       fac1    = drhogqv / min( -drhogqc-drhogqi, -eps ) ! limiting coefficient

       rhoq(I_QV,k,i,j) = rhoq(I_QV,k,i,j) + drhogqv
       rhoq(I_QC,k,i,j) = max(0.0_RP, rhoq(I_QC,k,i,j) + drhogqc*fac1)
       rhoq(I_QI,k,i,j) = max(0.0_RP, rhoq(I_QI,k,i,j) + drhogqi*fac1)
       rhoq(I_NC,k,i,j) = max(0.0_RP, rhoq(I_NC,k,i,j) + drhognc)
       rhoq(I_NI,k,i,j) = max(0.0_RP, rhoq(I_NI,k,i,j) + drhogni)

       ! cloud number concentration filter
       rhoq(I_NC,k,i,j) = min( rhoq(I_NC,k,i,j), nc_uplim_d(1,i,j) )

       rhoe(k,i,j) = rhoe(k,i,j) - LHV * drhogqv + LHF * drhogqi*fac1

       QTRC(k,i,j,I_QV) = rhoq(I_QV,k,i,j) * rrho(k,i,j)
       QTRC(k,i,j,I_QC) = rhoq(I_QC,k,i,j) * rrho(k,i,j)
       QTRC(k,i,j,I_QI) = rhoq(I_QI,k,i,j) * rrho(k,i,j)
       QTRC(k,i,j,I_NC) = rhoq(I_NC,k,i,j) * rrho(k,i,j)
       QTRC(k,i,j,I_NI) = rhoq(I_NI,k,i,j) * rrho(k,i,j)

       CALC_QDRY( qdry(k,i,j), QTRC, TRACER_MASS, k, i, j, iq )
       CALC_CV( cva(k,i,j), qdry(k,i,j), QTRC, k, i, j, iq, CVdry, TRACER_CV )
       CALC_R( Rmoist, qdry(k,i,j), QTRC, k, i, j, iq, Rdry, TRACER_R )
       temp(k,i,j) = rhoe(k,i,j) / ( DENS(k,i,j) * cva(k,i,j) )
       pres(k,i,j) = DENS(k,i,j) * Rmoist * temp(k,i,j)
       wtemp(k,i,j) = max( temp(k,i,j), tem_min )
    enddo
    enddo
    enddo

!    if( opt_debug )     call debugreport_nucleation
    if( opt_debug_tem ) call debug_tem_kij( 2, temp(:,:,:), DENS(:,:,:), pres(:,:,:), QTRC(:,:,:,I_QV) )


    call PROF_rapend  ('MP_Nucleation', 3)
    !----------------------------------------------------------------------------
    !
    ! 2.Phase change: Freezing, Melting, Vapor deposition
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP_Phase_change', 3)

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rhoq2(I_QR,k,i,j)     = rhoq(I_QR,k,i,j)
       rhoq2(I_NR,k,i,j)     = rhoq(I_NR,k,i,j)
       xq(I_mp_QR,k,i,j)     = max(xr_min,  min(xr_max, rhoq2(I_QR,k,i,j)/(rhoq2(I_NR,k,i,j)+nr_min) ))

       dq_xa(I_mp_QR,k,i,j)  = a_m(I_mp_QR)*xq(I_mp_QR,k,i,j)**b_m(I_mp_QR)
       vt_xa(I_mp_QR,1,k,i,j) = alpha_v(I_mp_QR,1)*(xq(I_mp_QR,k,i,j)**beta_v(I_mp_QR,1))*rho_fac_q(I_mp_QR,k,i,j)
       vt_xa(I_mp_QR,2,k,i,j) = vt_xa(I_mp_QR,1,k,i,j)

       !! Following values shoud be already filtered to be non-zero before sbroutine was called.
       ! Mass concentration [kg/m3]
       rhoq2(I_QV,k,i,j) = rhoq(I_QV,k,i,j)
       rhoq2(I_QC,k,i,j) = rhoq(I_QC,k,i,j)
       rhoq2(I_QI,k,i,j) = rhoq(I_QI,k,i,j)
       rhoq2(I_QS,k,i,j) = rhoq(I_QS,k,i,j)
       rhoq2(I_QG,k,i,j) = rhoq(I_QG,k,i,j)
       ! Number concentration[/m3] (should be filtered to avoid zero division.)
       rhoq2(I_NC,k,i,j) = rhoq(I_NC,k,i,j)
       rhoq2(I_NI,k,i,j) = rhoq(I_NI,k,i,j)
       rhoq2(I_NS,k,i,j) = rhoq(I_NS,k,i,j)
       rhoq2(I_NG,k,i,j) = rhoq(I_NG,k,i,j)

       ! Mass of mean particle [kg]
       ! SB06(94)
       !
       xq(I_mp_QC,k,i,j)     = min(xc_max, max(xc_min, rhoq2(I_QC,k,i,j)/(rhoq2(I_NC,k,i,j)+nc_min) ))
       xq(I_mp_QI,k,i,j)     = min(xi_max, max(xi_min, rhoq2(I_QI,k,i,j)/(rhoq2(I_NI,k,i,j)+ni_min) ))
       xq(I_mp_QS,k,i,j)     = min(xs_max, max(xs_min, rhoq2(I_QS,k,i,j)/(rhoq2(I_NS,k,i,j)+ns_min) ))
       xq(I_mp_QG,k,i,j)     = min(xg_max, max(xg_min, rhoq2(I_QG,k,i,j)/(rhoq2(I_NG,k,i,j)+ng_min) ))
       ! diamter of average mass
       ! SB06(32)
       dq_xa(I_mp_QC,k,i,j)  = a_m(I_mp_QC)*xq(I_mp_QC,k,i,j)**b_m(I_mp_QC)
       dq_xa(I_mp_QI,k,i,j)  = a_m(I_mp_QI)*xq(I_mp_QI,k,i,j)**b_m(I_mp_QI)
       dq_xa(I_mp_QS,k,i,j)  = a_m(I_mp_QS)*xq(I_mp_QS,k,i,j)**b_m(I_mp_QS)
       dq_xa(I_mp_QG,k,i,j)  = a_m(I_mp_QG)*xq(I_mp_QG,k,i,j)**b_m(I_mp_QG)

       ! terminal velocity of average mass
       vt_xa(I_mp_QC,1,k,i,j) = alpha_v(I_mp_QC,1)*(xq(I_mp_QC,k,i,j)**beta_v(I_mp_QC,1))*rho_fac_q(I_mp_QC,k,i,j)
       vt_xa(I_mp_QI,1,k,i,j) = alpha_v(I_mp_QI,1)*(xq(I_mp_QI,k,i,j)**beta_v(I_mp_QI,1))*rho_fac_q(I_mp_QI,k,i,j)
       vt_xa(I_mp_QS,1,k,i,j) = alpha_v(I_mp_QS,1)*(xq(I_mp_QS,k,i,j)**beta_v(I_mp_QS,1))*rho_fac_q(I_mp_QS,k,i,j)
       vt_xa(I_mp_QG,1,k,i,j) = alpha_v(I_mp_QG,1)*(xq(I_mp_QG,k,i,j)**beta_v(I_mp_QG,1))*rho_fac_q(I_mp_QG,k,i,j)
       vt_xa(I_mp_QC,2,k,i,j) = alpha_v(I_mp_QC,2)*(xq(I_mp_QC,k,i,j)**beta_v(I_mp_QC,2))*rho_fac_q(I_mp_QC,k,i,j)
       vt_xa(I_mp_QI,2,k,i,j) = alpha_v(I_mp_QI,2)*(xq(I_mp_QI,k,i,j)**beta_v(I_mp_QI,2))*rho_fac_q(I_mp_QI,k,i,j)
       vt_xa(I_mp_QS,2,k,i,j) = alpha_v(I_mp_QS,2)*(xq(I_mp_QS,k,i,j)**beta_v(I_mp_QS,2))*rho_fac_q(I_mp_QS,k,i,j)
       vt_xa(I_mp_QG,2,k,i,j) = alpha_v(I_mp_QG,2)*(xq(I_mp_QG,k,i,j)**beta_v(I_mp_QG,2))*rho_fac_q(I_mp_QG,k,i,j)

    end do
    end do
    end do

    call moist_psat_liq( esw, wtemp )
    call moist_psat_ice( esi, wtemp )

    call freezing_water_kij( &
         dt,             & ! in
         PQ,             & ! inout
         rhoq2, xq, temp ) ! in

    call dep_vapor_melt_ice_kij( &
         PQ,                 & ! inout
         DENS, wtemp, pres, qdry, & ! in
         rhoq2,               & ! in
         esw, esi,           & ! in
         xq,                 & ! in
         vt_xa,              & ! in
         dq_xa               ) ! in

    !
    ! update subroutine
    !
    call update_by_phase_change_kij( &
         ntdiv, ntmax_phase_change,  & ! in
         dt,                         & ! in
         gsgam2_d,                   & ! in
         z,                          & ! in
         dz,                         & ! in
         velz,                       & ! in
         dTdt_equiv_d,               & ! in
         DENS,                       & ! in
         rhoe,                       & ! inout
         rhoq, QTRC,                 & ! inout
         temp, pres,                 & ! inout
         cva,                        & ! out
         esw, esi, rhoq2,            & ! in
         PQ,                         & ! inout
         EVAPORATE,                  & ! out
         sl_PLCdep,                  & ! inout
         sl_PLRdep, sl_PNRdep        ) ! inout

!    if( opt_debug )     call debugreport_phasechange
    if( opt_debug_tem ) call debug_tem_kij( 3, temp(:,:,:), DENS(:,:,:), pres(:,:,:), QTRC(:,:,:,I_QV) )

    call PROF_rapend  ('MP_Phase_change', 3)

    !---------------------------------------------------------------------------
    !
    ! 3.Collection process
    !
    !---------------------------------------------------------------------------
    call PROF_rapstart('MP_Collection', 3)

    ! parameter setting
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! Mass concentration [kg/m3]
       rhoq2(I_QV,k,i,j) = rhoq(I_QV,k,i,j)
       rhoq2(I_QC,k,i,j) = rhoq(I_QC,k,i,j)
       rhoq2(I_QR,k,i,j) = rhoq(I_QR,k,i,j)
       rhoq2(I_QI,k,i,j) = rhoq(I_QI,k,i,j)
       rhoq2(I_QS,k,i,j) = rhoq(I_QS,k,i,j)
       rhoq2(I_QG,k,i,j) = rhoq(I_QG,k,i,j)
       ! Number concentration[/m3]
       rhoq2(I_NC,k,i,j) = rhoq(I_NC,k,i,j)
       rhoq2(I_NR,k,i,j) = rhoq(I_NR,k,i,j)
       rhoq2(I_NI,k,i,j) = rhoq(I_NI,k,i,j)
       rhoq2(I_NS,k,i,j) = rhoq(I_NS,k,i,j)
       rhoq2(I_NG,k,i,j) = rhoq(I_NG,k,i,j)

       ! Mass of mean particle [kg]
       xq(I_mp_QC,k,i,j) = min(xc_max, max(xc_min, rhoq2(I_QC,k,i,j)/(rhoq2(I_NC,k,i,j)+nc_min) ) )
       xq(I_mp_QR,k,i,j) = min(xr_max, max(xr_min, rhoq2(I_QR,k,i,j)/(rhoq2(I_NR,k,i,j)+nr_min) ) )
       xq(I_mp_QI,k,i,j) = min(xi_max, max(xi_min, rhoq2(I_QI,k,i,j)/(rhoq2(I_NI,k,i,j)+ni_min) ) )
       xq(I_mp_QS,k,i,j) = min(xs_max, max(xs_min, rhoq2(I_QS,k,i,j)/(rhoq2(I_NS,k,i,j)+ns_min) ) )
       xq(I_mp_QG,k,i,j) = min(xg_max, max(xg_min, rhoq2(I_QG,k,i,j)/(rhoq2(I_NG,k,i,j)+ng_min) ) )

       ! effective cross section is assume as area equivalent circle
       dq_xa(I_mp_QC,k,i,j) = 2.0_RP*a_rea(I_mp_QC)*xq(I_mp_QC,k,i,j)**b_rea(I_mp_QC)
       dq_xa(I_mp_QR,k,i,j) = 2.0_RP*a_rea(I_mp_QR)*xq(I_mp_QR,k,i,j)**b_rea(I_mp_QR)
       dq_xa(I_mp_QI,k,i,j) = 2.0_RP*a_rea(I_mp_QI)*xq(I_mp_QI,k,i,j)**b_rea(I_mp_QI)
       dq_xa(I_mp_QS,k,i,j) = 2.0_RP*a_rea(I_mp_QS)*xq(I_mp_QS,k,i,j)**b_rea(I_mp_QS)
       dq_xa(I_mp_QG,k,i,j) = 2.0_RP*a_rea(I_mp_QG)*xq(I_mp_QG,k,i,j)**b_rea(I_mp_QG)

       ! terminal velocity of average mass
       ! SB06(33)
       vt_xa(I_mp_QC,2,k,i,j) = alpha_v(I_mp_QC,2)*(xq(I_mp_QC,k,i,j)**beta_v(I_mp_QC,2))*rho_fac_q(I_mp_QC,k,i,j)
       vt_xa(I_mp_QR,2,k,i,j) = alpha_v(I_mp_QR,2)*(xq(I_mp_QR,k,i,j)**beta_v(I_mp_QR,2))*rho_fac_q(I_mp_QR,k,i,j)
       vt_xa(I_mp_QI,2,k,i,j) = alpha_v(I_mp_QI,2)*(xq(I_mp_QI,k,i,j)**beta_v(I_mp_QI,2))*rho_fac_q(I_mp_QI,k,i,j)
       vt_xa(I_mp_QS,2,k,i,j) = alpha_v(I_mp_QS,2)*(xq(I_mp_QS,k,i,j)**beta_v(I_mp_QS,2))*rho_fac_q(I_mp_QS,k,i,j)
       vt_xa(I_mp_QG,2,k,i,j) = alpha_v(I_mp_QG,2)*(xq(I_mp_QG,k,i,j)**beta_v(I_mp_QG,2))*rho_fac_q(I_mp_QG,k,i,j)
    enddo
    enddo
    enddo

    ! Auto-conversion, Accretion, Self-collection, Break-up
    ! [Mod] T.Seiki
    if ( MP_doautoconversion ) then
       call aut_acc_slc_brk_kij(  &
            PQ, & ! inout
            rhoq2, xq, dq_xa, &
            DENS               )
    else
!OCL XFILL
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PQ(I_LCaut,k,i,j) = 0.0_RP
          PQ(I_NCaut,k,i,j) = 0.0_RP
          PQ(I_NRaut,k,i,j) = 0.0_RP
          PQ(I_LCacc,k,i,j) = 0.0_RP
          PQ(I_NCacc,k,i,j) = 0.0_RP
          PQ(I_NRslc,k,i,j) = 0.0_RP
          PQ(I_NRbrk,k,i,j) = 0.0_RP
       end do
       end do
       end do
    endif

    call mixed_phase_collection_kij( &
         ! collection process
         Pac, PQ,                    & ! inout
         temp, rhoq2,                & ! in
         xq, dq_xa, vt_xa            ) ! in
!         DENS(:,:,:),                ) ! in

    call ice_multiplication_kij( &
         PQ,                     & ! inout
         Pac,                    & ! in
         temp, rhoq2, xq         ) ! in

    !
    ! update
    ! rhogq = l*gsgam
    !
    PROFILE_START("sn14_update_rhoq")
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       ! warm collection process
       wrm_dqc = max( dt*( PQ(I_LCaut,k,i,j)+PQ(I_LCacc,k,i,j) ), -rhoq2(I_QC,k,i,j)  )
       wrm_dnc = max( dt*( PQ(I_NCaut,k,i,j)+PQ(I_NCacc,k,i,j) ), -rhoq2(I_NC,k,i,j)  )
       wrm_dnr = max( dt*( PQ(I_NRaut,k,i,j)+PQ(I_NRslc,k,i,j)+PQ(I_NRbrk,k,i,j) ), -rhoq2(I_NR,k,i,j) )
       wrm_dqr = -wrm_dqc
       ! mixed phase collection
       ! Pxxacyy2zz xx and yy decrease and zz increase .
       !
       ! At first fixer is applied to decreasing particles.
       ! order of fixer: graupel-cloud, snow-cloud, ice-cloud, graupel-rain, snow-rain, ice-rain,
       !                 snow-ice,  ice-ice, graupel-snow, snow-snow
       ! cloud mass decrease
       gc_dqc  = max( dt*Pac(I_LGacLC2LG,k,i,j)  , min(0.0_RP, -rhoq2(I_QC,k,i,j)-wrm_dqc               )) ! => dqg
       sc_dqc  = max( dt*Pac(I_LSacLC2LS,k,i,j)  , min(0.0_RP, -rhoq2(I_QC,k,i,j)-wrm_dqc-gc_dqc        )) ! => dqs
       ic_dqc  = max( dt*Pac(I_LIacLC2LI,k,i,j)  , min(0.0_RP, -rhoq2(I_QC,k,i,j)-wrm_dqc-gc_dqc-sc_dqc )) ! => dqi
       ! cloud num. decrease
       gc_dnc  = max( dt*Pac(I_NGacNC2NG,k,i,j)  , min(0.0_RP, -rhoq2(I_NC,k,i,j)-wrm_dnc               )) ! => dnc:minus
       sc_dnc  = max( dt*Pac(I_NSacNC2NS,k,i,j)  , min(0.0_RP, -rhoq2(I_NC,k,i,j)-wrm_dnc-gc_dnc        )) ! => dnc:minus
       ic_dnc  = max( dt*Pac(I_NIacNC2NI,k,i,j)  , min(0.0_RP, -rhoq2(I_NC,k,i,j)-wrm_dnc-gc_dnc-sc_dnc )) ! => dnc:minus

       ! rain mass decrease ( tem < 273.15K)
       sw = sign(0.5_RP, T00-temp(k,i,j)) + 0.5_RP ! if( temp(k,i,j) <= T00 )then sw=1, else sw=0
       rg_dqr  = max( dt*Pac(I_LRacLG2LG,  k,i,j), min(0.0_RP, -rhoq2(I_QR,k,i,j)-wrm_dqr               )) * sw
       rg_dqg  = max( dt*Pac(I_LRacLG2LG,  k,i,j), min(0.0_RP, -rhoq2(I_QG,k,i,j)                       )) * ( 1.0_RP - sw )
       rs_dqr  = max( dt*Pac(I_LRacLS2LG_R,k,i,j), min(0.0_RP, -rhoq2(I_QR,k,i,j)-wrm_dqr-rg_dqr        )) * sw
       ri_dqr  = max( dt*Pac(I_LRacLI2LG_R,k,i,j), min(0.0_RP, -rhoq2(I_QR,k,i,j)-wrm_dqr-rg_dqr-rs_dqr )) * sw
       ! rain num. decrease
       rg_dnr  = max( dt*Pac(I_NRacNG2NG,  k,i,j), min(0.0_RP, -rhoq2(I_NR,k,i,j)-wrm_dnr               )) * sw
       rg_dng  = max( dt*Pac(I_NRacNG2NG,  k,i,j), min(0.0_RP, -rhoq2(I_NG,k,i,j)                       )) * ( 1.0_RP - sw )
       rs_dnr  = max( dt*Pac(I_NRacNS2NG_R,k,i,j), min(0.0_RP, -rhoq2(I_NR,k,i,j)-wrm_dnr-rg_dnr        )) * sw
       ri_dnr  = max( dt*Pac(I_NRacNI2NG_R,k,i,j), min(0.0_RP, -rhoq2(I_NR,k,i,j)-wrm_dnr-rg_dnr-rs_dnr )) * sw

       ! ice mass decrease
       fac1    = (ri_dqr-eps)/ (dt*Pac(I_LRacLI2LG_R,k,i,j)-eps) ! suppress factor by filter of rain
       ri_dqi  = max( dt*Pac(I_LRacLI2LG_I,k,i,j)*fac1, min(0.0_RP, -rhoq2(I_QI,k,i,j)+ic_dqc               )) ! => dqg
       ii_dqi  = max( dt*Pac(I_LIacLI2LS,k,i,j)       , min(0.0_RP, -rhoq2(I_QI,k,i,j)+ic_dqc-ri_dqi        )) ! => dqs
       is_dqi  = max( dt*Pac(I_LIacLS2LS,k,i,j)       , min(0.0_RP, -rhoq2(I_QI,k,i,j)+ic_dqc-ri_dqi-ii_dqi )) ! => dqs
       ! ice num. decrease
       fac4    = (ri_dnr-eps)/ (dt*Pac(I_NRacNI2NG_R,k,i,j)-eps) ! suppress factor by filter of rain
       ri_dni  = max( dt*Pac(I_NRacNI2NG_I,k,i,j)*fac4, min(0.0_RP, -rhoq2(I_NI,k,i,j)               )) ! => dni:minus
       ii_dni  = max( dt*Pac(I_NIacNI2NS,k,i,j)       , min(0.0_RP, -rhoq2(I_NI,k,i,j)-ri_dni        )) ! => dni:minus,dns:plus(*0.5)
       is_dni  = max( dt*Pac(I_NIacNS2NS,k,i,j)       , min(0.0_RP, -rhoq2(I_NI,k,i,j)-ri_dni-ii_dni )) ! => dni:minus,dns:plus
       ! snow mass decrease
       fac3    = (rs_dqr-eps)/(dt*Pac(I_LRacLS2LG_R,k,i,j)-eps) ! suppress factor by filter of rain
       rs_dqs  = max( dt*Pac(I_LRacLS2LG_S,k,i,j)*fac3, min(0.0_RP, -rhoq2(I_QS,k,i,j)+sc_dqc+ii_dqi+is_dqi        )) ! => dqg
       gs_dqs  = max( dt*Pac(I_LGacLS2LG,k,i,j)       , min(0.0_RP, -rhoq2(I_QS,k,i,j)+sc_dqc+ii_dqi+is_dqi-rs_dqs )) ! => dqg
       ! snow num. decrease
       fac6    = (rs_dnr-eps)/(dt*Pac(I_NRacNS2NG_R,k,i,j)-eps) ! suppress factor by filter of rain
!       fac7    = (is_dni-eps)/(dt*Pac(I_NIacNS2NS,  k,i,j)-eps) ! suppress factor by filter of ice
       rs_dns  = max( dt*Pac(I_NRacNS2NG_S,k,i,j)*fac6, min(0.0_RP, -rhoq2(I_NS,k,i,j)+0.50_RP*ii_dni+is_dni       )) ! => dns:minus
       gs_dns  = max( dt*Pac(I_NGacNS2NG,k,i,j)       , min(0.0_RP, -rhoq2(I_NS,k,i,j)+0.50_RP*ii_dni+is_dni-rs_dns )) ! => dns:minus
       ss_dns  = max( dt*Pac(I_NSacNS2NS,k,i,j)       , min(0.0_RP, -rhoq2(I_NS,k,i,j)+0.50_RP*ii_dni+is_dni-rs_dns-gs_dns ))
       !
       gg_dng  = max( dt*Pac(I_NGacNG2NG,k,i,j)       , min(0.0_RP, -rhoq2(I_NG,k,i,j) ))
       !
       ! total plus in mixed phase collection(clp_)
       ! mass
       ! if( temp(k,i,j) <= T00 )then sw=1, else sw=0
       clp_dqc =  0.0_RP
       clp_dqr = (-rg_dqg-rs_dqs-ri_dqi) * (1.0_RP-sw)
       clp_dqi = -ic_dqc
       clp_dqs = -sc_dqc-ii_dqi-is_dqi
       clp_dqg = -gc_dqc -gs_dqs  + (-rg_dqr-rs_dqr-rs_dqs-ri_dqr-ri_dqi) * sw
       ! num.( number only increase when a+b=>c,  dnc=-dna)
       clp_dnc = 0.0_RP
       clp_dnr = 0.0_RP
       clp_dni = 0.0_RP
       clp_dns = -ii_dni*0.5_RP
       clp_dng = (-rs_dnr-ri_dnr) * sw
       ! total minus in mixed phase collection(clm_)
       ! mass
       clm_dqc = gc_dqc+sc_dqc+ic_dqc
       clm_dqr = (rg_dqr+rs_dqr+ri_dqr) * sw
       clm_dqi = ri_dqi+ii_dqi+is_dqi
       clm_dqs = rs_dqs+gs_dqs
       clm_dqg = rg_dqg * (1.0_RP-sw)
       ! num.
       clm_dnc = gc_dnc+sc_dnc+ic_dnc
       clm_dnr = (rg_dnr+rs_dnr+ri_dnr) * sw
       clm_dni = ri_dni+ii_dni+is_dni
       clm_dns = rs_dns+ss_dns+gs_dns
       clm_dng = gg_dng + rg_dng * (1.0_RP-sw)

       ! partial conversion
       ! 08/05/08 [Mod] T.Mitsui
       pco_dqi = max( dt*PQ(I_LIcon,k,i,j), -clp_dqi )
       pco_dqs = max( dt*PQ(I_LScon,k,i,j), -clp_dqs )
       pco_dqg = -pco_dqi-pco_dqs
       ! 08/05/08 [Mod] T.Mitsui
       pco_dni = max( dt*PQ(I_NIcon,k,i,j), -clp_dni )
       pco_dns = max( dt*PQ(I_NScon,k,i,j), -clp_dns )
       pco_dng = -pco_dni-pco_dns
       ! enhanced melting ( always negative value )
       ! ice-cloud melting produces cloud, others produce rain
       eml_dqi =  max( dt*PQ(I_LIacm,k,i,j), min(0.0_RP, -rhoq2(I_QI,k,i,j)-(clp_dqi+clm_dqi)-pco_dqi ))
       eml_dqs =  max( dt*PQ(I_LSacm,k,i,j), min(0.0_RP, -rhoq2(I_QS,k,i,j)-(clp_dqs+clm_dqs)-pco_dqs ))
       eml_dqg =  max( dt*(PQ(I_LGacm,k,i,j)+PQ(I_LGarm,k,i,j)+PQ(I_LSarm,k,i,j)+PQ(I_LIarm,k,i,j)), &
                  min(0.0_RP, -rhoq2(I_QG,k,i,j)-(clp_dqg+clm_dqg)-pco_dqg ))
       eml_dqc = -eml_dqi
       eml_dqr = -eml_dqs-eml_dqg
       !
       eml_dni =  max( dt*PQ(I_NIacm,k,i,j), min(0.0_RP, -rhoq2(I_NI,k,i,j)-(clp_dni+clm_dni)-pco_dni ))
       eml_dns =  max( dt*PQ(I_NSacm,k,i,j), min(0.0_RP, -rhoq2(I_NS,k,i,j)-(clp_dns+clm_dns)-pco_dns ))
       eml_dng =  max( dt*(PQ(I_NGacm,k,i,j)+PQ(I_NGarm,k,i,j)+PQ(I_NSarm,k,i,j)+PQ(I_NIarm,k,i,j)), &
                  min(0.0_RP, -rhoq2(I_NG,k,i,j)-(clp_dng+clm_dng)-pco_dng ))
       eml_dnc = -eml_dni
       eml_dnr = -eml_dns-eml_dng
       !
       ! ice multiplication
       spl_dqg = max( dt*PQ(I_LGspl,k,i,j), min(0.0_RP, -rhoq2(I_QG,k,i,j)-(clp_dqg+clm_dqg)-pco_dqg-eml_dqg ))
       spl_dqs = max( dt*PQ(I_LSspl,k,i,j), min(0.0_RP, -rhoq2(I_QS,k,i,j)-(clp_dqs+clm_dqs)-pco_dqs-eml_dqs ))
       spl_dqi = -spl_dqg-spl_dqs
       fac9    = (spl_dqg-eps)/(dt*PQ(I_LGspl,k,i,j)-eps)
       fac10   = (spl_dqs-eps)/(dt*PQ(I_LSspl,k,i,j)-eps)
       spl_dni = dt*PQ(I_NIspl,k,i,j)*fac9*fac10
       !
       ! total cloud change
       drhogqc = (wrm_dqc + clp_dqc + clm_dqc           + eml_dqc )
       drhognc = (wrm_dnc + clp_dnc + clm_dnc           + eml_dnc )
       ! total rain change
       drhogqr = (wrm_dqr + clp_dqr + clm_dqr           + eml_dqr )
       drhognr = (wrm_dnr + clp_dnr + clm_dnr           + eml_dnr )
       ! total ice change
       drhogqi = (          clp_dqi + clm_dqi + pco_dqi + eml_dqi + spl_dqi)
       drhogni = (          clp_dni + clm_dni + pco_dni + eml_dni + spl_dni)
       ! total snow change
       drhogqs = (          clp_dqs + clm_dqs + pco_dqs + eml_dqs + spl_dqs)
       drhogns = (          clp_dns + clm_dns + pco_dns + eml_dns )
       ! total graupel change
       drhogqg = (          clp_dqg + clm_dqg + pco_dqg + eml_dqg + spl_dqg)
       drhogng = (          clp_dng + clm_dng + pco_dng + eml_dng )
       !
       !--- update
       !
       rhoq(I_QC,k,i,j) = max(0.0_RP, rhoq(I_QC,k,i,j) + drhogqc )
       rhoq(I_NC,k,i,j) = max(0.0_RP, rhoq(I_NC,k,i,j) + drhognc )
       rhoq(I_QR,k,i,j) = max(0.0_RP, rhoq(I_QR,k,i,j) + drhogqr )
       rhoq(I_NR,k,i,j) = max(0.0_RP, rhoq(I_NR,k,i,j) + drhognr )
       rhoq(I_QI,k,i,j) = max(0.0_RP, rhoq(I_QI,k,i,j) + drhogqi )
       rhoq(I_NI,k,i,j) = max(0.0_RP, rhoq(I_NI,k,i,j) + drhogni )
       rhoq(I_QS,k,i,j) = max(0.0_RP, rhoq(I_QS,k,i,j) + drhogqs )
       rhoq(I_NS,k,i,j) = max(0.0_RP, rhoq(I_NS,k,i,j) + drhogns )
       rhoq(I_QG,k,i,j) = max(0.0_RP, rhoq(I_QG,k,i,j) + drhogqg )
       rhoq(I_NG,k,i,j) = max(0.0_RP, rhoq(I_NG,k,i,j) + drhogng )
       !
       ! update
       ! rhogq = l*gsgam
       rhoe(k,i,j) = rhoe(k,i,j) + LHF * ( drhogqi + drhogqs + drhogqg )
    enddo
    enddo
    enddo
    PROFILE_STOP("sn14_update_rhoq")

    call PROF_rapend  ('MP_Collection', 3)

    call PROF_rapstart('MP_Postprocess', 3)

    !--- update mixing ratio
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       do iq = I_QV, I_NG
          QTRC(k,i,j,iq) = rhoq(iq,k,i,j) * rrho(k,i,j)
       enddo

       CALC_QDRY( qdry(k,i,j), QTRC, TRACER_MASS, k, i, j, iq )
       CALC_CV( cva(k,i,j), qdry(k,i,j), QTRC, k, i, j, iq, CVdry, TRACER_CV )
       CALC_R( Rmoist, qdry(k,i,j), QTRC, k, i, j, iq, Rdry, TRACER_R )
       cpa(k,i,j) = cva(k,i,j) + Rmoist
       temp(k,i,j) = rhoe(k,i,j) / ( DENS(k,i,j) * cva(k,i,j) )
       pres(k,i,j) = DENS(k,i,j) * Rmoist * temp(k,i,j)
       RHOT(k,i,j) = temp(k,i,j) * ( P00 / pres(k,i,j) )**(Rmoist/cpa(k,i,j)) &
               * DENS(k,i,j)
    enddo
    enddo
    enddo

!    if( opt_debug )     call debugreport_collection
    if( opt_debug_tem ) call debug_tem_kij( 4, temp(:,:,:), DENS(:,:,:), pres(:,:,:), QTRC(:,:,:,I_QV) )

    call PROF_rapend  ('MP_Postprocess', 3)

    !----------------------------------------------------------------------------
    !
    ! 4.Saturation adjustment
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP_Saturation_adjustment', 3)
    ! nothing to do
    call PROF_rapend  ('MP_Saturation_adjustment', 3)
    !----------------------------------------------------------------------------
    !
    ! 5. Sedimentation ( terminal velocity must be negative )
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP_Sedimentation', 3)

    if ( MP_doprecipitation ) then

    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE
       FLX_rain(k,i,j) = 0.0_RP
       FLX_snow(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    velw(:,:,:,:) = 0.0_RP

    do step = 1, MP_NSTEP_SEDIMENTATION

       call MP_terminal_velocity( velw(:,:,:,:), & ! [OUT]
                                  rhoq(:,:,:,:), & ! [IN]
                                  DENS(:,:,:),   & ! [IN]
                                  temp(:,:,:),   & ! [IN]
                                  pres(:,:,:)    ) ! [IN]

       call MP_precipitation( wflux_rain(:,:,:),     & ! [OUT]
                              wflux_snow(:,:,:),     & ! [OUT]
                              DENS      (:,:,:),     & ! [INOUT]
                              MOMZ      (:,:,:),     & ! [INOUT]
                              MOMX      (:,:,:),     & ! [INOUT]
                              MOMY      (:,:,:),     & ! [INOUT]
                              rhoe      (:,:,:),     & ! [INOUT]
                              QTRC      (:,:,:,:),   & ! [INOUT]
                              QA_MP,                 & ! [IN]
                              QS_MP,                 & ! [IN]
                              velw      (:,:,:,:),   & ! [IN]
                              temp      (:,:,:),     & ! [IN]
                              TRACER_CV(:),          & ! [IN]
                              MP_DTSEC_SEDIMENTATION ) ! [IN]

       do j = JS, JE
       do i = IS, IE
       do k = KS-1, KE
          FLX_rain(k,i,j) = FLX_rain(k,i,j) + wflux_rain(k,i,j) * MP_RNSTEP_SEDIMENTATION
          FLX_snow(k,i,j) = FLX_snow(k,i,j) + wflux_snow(k,i,j) * MP_RNSTEP_SEDIMENTATION
       enddo
       enddo
       enddo

    enddo

    endif

    do j = JS, JE
    do i = IS, IE
       SFLX_rain(i,j) = FLX_rain(KS-1,i,j)
       SFLX_snow(i,j) = FLX_snow(KS-1,i,j)
    end do
    end do

    call PROF_rapend  ('MP_Sedimentation', 3)

    return
  end subroutine mp_sn14

  !-----------------------------------------------------------------------------
  subroutine debug_tem_kij( &
      point,    &
      tem,      &
      rho,      &
      pre,      &
      qv        )
    use scale_process, only: &
       PRC_myrank
    implicit none

    integer, intent(in) :: point
    real(RP), intent(in) :: tem(KA,IA,JA)
    real(RP), intent(in) :: rho(KA,IA,JA)
    real(RP), intent(in) :: pre(KA,IA,JA)
    real(RP), intent(in) :: qv (KA,IA,JA)

    integer :: k ,i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if (      tem(k,i,j) < tem_min &
            .OR. rho(k,i,j) < rho_min &
            .OR. pre(k,i,j) < 1.0_RP    ) then

          if( IO_L ) write(IO_FID_LOG,'(A,I3,A,4(F16.5),3(I6))') &
          "*** point: ", point, " low tem,rho,pre:", tem(k,i,j), rho(k,i,j), pre(k,i,j), qv(k,i,j), k, i, j, PRC_myrank
       endif
    enddo
    enddo
    enddo

    return
  end subroutine debug_tem_kij

  subroutine nucleation_kij( &
       z, velz,          &
       rho, tem, pre,    &
       rhoq,             &
       PQ,               &
       cpa,              & ! in
       dTdt_rad,         & ! in
       qke,              & ! in
       CCN,               & ! in
       dt                ) ! in
    use scale_process, only: &
       PRC_MPIstop
    use scale_tracer, only: &
       QA
    use scale_atmos_saturation, only: &
       moist_psat_liq       => ATMOS_SATURATION_psat_liq, &
       moist_psat_ice       => ATMOS_SATURATION_psat_ice,   &
       moist_pres2qsat_liq  => ATMOS_SATURATION_pres2qsat_liq, &
       moist_pres2qsat_ice  => ATMOS_SATURATION_pres2qsat_ice,   &
       moist_dqsi_dtem_rho  => ATMOS_SATURATION_dqsi_dtem_rho
    implicit none

    real(RP), intent(in)  :: z(KA)      !
    real(RP), intent(in)  :: velz(KA,IA,JA)   ! w of half point
    real(RP), intent(in)  :: rho(KA,IA,JA)    ! [Add] 09/08/18 T.Mitsui
    real(RP), intent(in)  :: tem(KA,IA,JA)    ! [Add] 09/08/18 T.Mitsui
    real(RP), intent(in)  :: pre(KA,IA,JA)    ! [Add] 09/08/18 T.Mitsui
    !
    real(RP), intent(in)  :: rhoq(I_QV:I_NG,KA,IA,JA)     !
    real(RP), intent(out) :: PQ(PQ_MAX,KA,IA,JA)
    !
    real(RP), intent(in) ::  cpa(KA,IA,JA)      ! in  09/08/18 [Add] T.Mitsui
    real(RP), intent(in)  :: dTdt_rad(KA,IA,JA) ! 09/08/18 T.Mitsui
    real(RP), intent(in)  :: qke(KA,IA,JA)      ! 09/08/18 T.Mitsui
    real(RP), intent(in)  :: dt
    real(RP), intent(in)  :: CCN(KA,IA,JA)
    !
    ! namelist variables
    !
    ! total aerosol number concentration [/m3]
    real(RP), parameter :: c_ccn_ocean= 1.00E+8_RP
    real(RP), parameter :: c_ccn_land = 1.26E+9_RP
    real(RP), save      :: c_ccn      = 1.00E+8_RP
    ! aerosol activation factor
    real(RP), parameter :: kappa_ocean= 0.462_RP
    real(RP), parameter :: kappa_land = 0.308_RP
    real(RP), save      :: kappa      = 0.462_RP
    real(RP), save      :: c_in       = 1.0_RP
    ! SB06 (36)
    real(RP), save :: nm_M92 = 1.E+3_RP
    real(RP), save :: am_M92 = -0.639_RP
    real(RP), save :: bm_M92 = 12.96_RP
    !
    real(RP), save :: in_max = 1000.E+3_RP ! max num. of Ice-Nuclei [num/m3]
    real(RP), save :: ssi_max= 0.60_RP
    real(RP), save :: ssw_max= 1.1_RP  ! [%]
    !
    logical, save :: flag_first = .true.
    real(RP), save :: qke_min = 0.03_RP ! sigma=0.1[m/s], 09/08/18 T.Mitsui
    real(RP), save :: tem_ccn_low=233.150_RP  ! = -40 degC  ! [Add] 10/08/03 T.Mitsui
    real(RP), save :: tem_in_low =173.150_RP  ! = -100 degC ! [Add] 10/08/03 T.Mitsui
    logical, save :: nucl_twomey = .false.
    logical, save :: inucl_w     = .false.
    !
    namelist /nm_mp_sn14_nucleation/ &
         in_max,                     & !
         c_ccn, kappa,               & ! cloud nucleation
         nm_M92, am_M92, bm_M92,     & ! ice nucleation
         xc_ccn, xi_ccn,             &
         tem_ccn_low,                & ! [Add] 10/08/03 T.Mitsui
         tem_in_low,                 & ! [Add] 10/08/03 T.Mitsui
         ssw_max, ssi_max,           &
         nucl_twomey, inucl_w        ! [Add] 13/01/30 Y.Sato
    !
!    real(RP) :: c_ccn_map(1,IA,JA)   ! c_ccn horizontal distribution
!    real(RP) :: kappa_map(1,IA,JA)   ! kappa horizontal distribution
!    real(RP) :: c_in_map(1,IA,JA)    ! c_in  horizontal distribution ! [Add] 11/08/30 T.Mitsui
    real(RP) :: esw(KA,IA,JA)      ! saturation vapor pressure, water
    real(RP) :: esi(KA,IA,JA)      !                            ice
    real(RP) :: ssw(KA,IA,JA)      ! super saturation (water)
    real(RP) :: ssi(KA,IA,JA)      ! super saturation (ice)
!    real(RP) :: w_dsswdz(KA,IA,JA) ! w*(d_ssw/ d_z) super saturation(water) flux
    real(RP) :: w_dssidz(KA,IA,JA) ! w*(d_ssi/ d_z), 09/04/14 T.Mitsui
!    real(RP) :: ssw_below(KA,IA,JA)! ssw(k-1)
    real(RP) :: ssi_below(KA,IA,JA)! ssi(k-1), 09/04/14 T.Mitsui
    real(RP) :: z_below(KA,IA,JA)  ! z(k-1)
    real(RP) :: dz                   ! z(k)-z(k-1)
    real(RP) :: pv                   ! vapor pressure
    ! work variables for Twomey Equation.
    real(RP) :: qsw(KA,IA,JA)
    real(RP) :: qsi(KA,IA,JA)
    real(RP) :: dqsidtem_rho(KA,IA,JA)
    real(RP) :: dssidt_rad(KA,IA,JA)
!    real(RP) :: dni_ratio(KA,IA,JA)
    real(RP) :: wssi, wdssi
    !
!    real(RP) :: xi_nuc(1,IA,JA)    ! xi use the value @ cloud base
!    real(RP) :: alpha_nuc(1,IA,JA) ! alpha_nuc
!    real(RP) :: eta_nuc(1,IA,JA)   ! xi use the value @ cloud base
    !
    real(RP) :: sigma_w(KA,IA,JA)
    real(RP) :: weff(KA,IA,JA)
    real(RP) :: weff_max(KA,IA,JA)
    !
    real(RP) :: coef_ccn(IA,JA)
    real(RP) :: slope_ccn(IA,JA)
    real(RP) :: nc_new(KA,IA,JA)
    real(RP) :: nc_new_below(KA,IA,JA)
    real(RP) :: dnc_new
    real(RP) :: nc_new_max   ! Lohmann (2002),
    real(RP) :: a_max
    real(RP) :: b_max
    logical :: flag_nucleation(KA,IA,JA)
    !
    real(RP) :: r_gravity
    real(RP), parameter :: r_sqrt3=0.577350269_RP ! = sqrt(1.d0/3.d0)
    real(RP), parameter :: eps=1.E-30_RP
    !====> ! 09/08/18
    !
    real(RP) :: dlcdt_max, dli_max ! defined by supersaturation
    real(RP) :: dncdt_max, dni_max ! defined by supersaturation
    real(RP) :: rdt
    !
    integer :: i, j, k
    !
    if( flag_first )then
       rewind(IO_FID_CONF)
       read(IO_FID_CONF, nml=nm_mp_sn14_nucleation, end=100)
100    if( IO_NML ) write(IO_FID_NML,nml=nm_mp_sn14_nucleation)
       flag_first=.false.

       if ( MP_couple_aerosol .AND. nucl_twomey ) then
          write(*,*) "xxx [mp_sn14/nucleation] nucl_twomey should be false when MP_couple_aerosol is true, stop"
          call PRC_MPIstop
       endif
    endif
    !
!    c_ccn_map(1,:,:) = c_ccn
!    kappa_map(1,:,:) = kappa
!    c_in_map(1,:,:)  = c_in
    !
!    nc_uplim_d(1,:,:)  = c_ccn_map(1,:,:)*1.5_RP
    do j = JS, JE
    do i = IS, IE
       nc_uplim_d(1,i,j)  = c_ccn*1.5_RP
    end do
    end do
    !
    rdt            = 1.0_RP/dt
    r_gravity      = 1.0_RP/GRAV
    !
    call moist_psat_liq     ( esw, tem )
    call moist_psat_ice     ( esi, tem )
    call moist_pres2qsat_liq( qsw, tem, pre )
    call moist_pres2qsat_ice( qsi, tem, pre )
    call moist_dqsi_dtem_rho( dqsidtem_rho, tem, rho )
    !
    ! Lohmann (2002),JAS, eq.(1) but changing unit [cm-3] => [m-3]
    a_max = 1.E+6_RP*0.1_RP*(1.E-6_RP)**1.27_RP
    b_max = 1.27_RP
    !
    ssi_max = 1.0_RP

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          pv        = rhoq(I_QV,k,i,j)*Rvap*tem(k,i,j)
          ssw(k,i,j) = min( MP_ssw_lim, ( pv/esw(k,i,j)-1.0_RP ) )*100.0_RP
          ssi(k,i,j) = (pv/esi(k,i,j) - 1.00_RP)
!          ssw_below(k+1,i,j) = ssw(k,i,j)
          ssi_below(k+1,i,j) = ssi(k,i,j)
          z_below(k+1,i,j)   = z(k)
       end do
!       ssw_below(KS,i,j) = ssw(KS,i,j)
       ssi_below(KS,i,j) = ssi(KS,i,j)
       z_below(KS,i,j)   = z(KS-1)

       ! dS/dz is evaluated by first order upstream difference
       !***  Solution for Twomey Equation ***
!       coef_ccn(i,j)  = 1.E+6_RP*0.88_RP*(c_ccn_map(1,i,j)*1.E-6_RP)**(2.0_RP/(kappa_map(1,i,j) + 2.0_RP)) * &
       coef_ccn(i,j)  = 1.E+6_RP*0.88_RP*(c_ccn*1.E-6_RP)**(2.0_RP/(kappa + 2.0_RP)) &
!           * (70.0_RP)**(kappa_map(1,i,j)/(kappa_map(1,i,j) + 2.0_RP))
            * (70.0_RP)**(kappa/(kappa + 2.0_RP))
!       slope_ccn(i,j) = 1.5_RP*kappa_map(1,i,j)/(kappa_map(1,i,j) + 2.0_RP)
       slope_ccn(i,j) = 1.5_RP*kappa/(kappa + 2.0_RP)
       !
       do k=KS, KE
          sigma_w(k,i,j) = r_sqrt3*sqrt(max(qke(k,i,j),qke_min))
       end do
       sigma_w(KS-1,i,j) = sigma_w(KS,i,j)
       sigma_w(KE+1,i,j) = sigma_w(KE,i,j)
       ! effective vertical velocity
       do k=KS, KE
          weff(k,i,j) = 0.5_RP*(velz(k-1,i,j) + velz(k,i,j)) - cpa(k,i,j)*r_gravity*dTdt_rad(k,i,j)
       end do

    end do
    end do
    !
    if( MP_couple_aerosol ) then

         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            if( ssw(k,i,j) > 1.e-10_RP .AND. pre(k,i,j) > 300.E+2_RP ) then
               nc_new(k,i,j) = max( CCN(k,i,j), c_ccn )
            else
               nc_new(k,i,j) = 0.0_RP
            endif
         enddo
         enddo
         enddo

    else

      if( nucl_twomey ) then
        ! diagnose cloud condensation nuclei

        do j = JS, JE
        do i = IS, IE
        do k = KS, KE
           ! effective vertical velocity (maximum vertical velocity in turbulent flow)
           weff_max(k,i,j) = weff(k,i,j) + sigma_w(k,i,j)
           ! large scale upward motion region and saturated
           if( (weff(k,i,j) > 1.E-8_RP) .AND. (ssw(k,i,j) > 1.E-10_RP)  .AND. pre(k,i,j) > 300.E+2_RP )then
              ! Lohmann (2002), eq.(1)
              nc_new_max   = coef_ccn(i,j)*weff_max(k,i,j)**slope_ccn(i,j)
              nc_new(k,i,j) = a_max*nc_new_max**b_max
           else
              nc_new(k,i,j) = 0.0_RP
           end if
        end do
        end do
        end do
      else
        ! calculate cloud condensation nuclei
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if( ssw(k,i,j) > 1.e-10_RP .AND. pre(k,i,j) > 300.E+2_RP ) then
                nc_new(k,i,j) = c_ccn*ssw(k,i,j)**kappa
             else
                nc_new(k,i,j) = 0.0_RP
             endif
          enddo
          enddo
          enddo
      endif

    endif

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          ! nc_new is bound by upper limit
          if( nc_new(k,i,j) > nc_uplim_d(1,i,j) )then ! no more CCN
             flag_nucleation(k,i,j) = .false.
             nc_new_below(k+1,i,j)  = 1.E+30_RP
          else if( nc_new(k,i,j) > eps )then ! nucleation can occur
             flag_nucleation(k,i,j) = .true.
             nc_new_below(k+1,i,j)  = nc_new(k,i,j)
          else ! nucleation cannot occur(unsaturated or negative w)
             flag_nucleation(k,i,j) = .false.
             nc_new_below(k+1,i,j)  = 0.0_RP
          end if
       end do
       nc_new_below(KS,i,j) = 0.0_RP
!       do k=KS, KE
           ! search maximum value of nc_new
!          if(  ( nc_new(k,i,j) < nc_new_below(k,i,j) ) .OR. &
!               ( nc_new_below(k,i,j) > c_ccn_map(1,i,j)*0.05_RP ) )then ! 5% of c_ccn
!               ( nc_new_below(k,i,j) > c_ccn*0.05_RP ) )then ! 5% of c_ccn
!             flag_nucleation(k,i,j) = .false.
!          end if
!       end do

    end do
    end do

    if( MP_couple_aerosol ) then

         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            ! nucleation occurs at only cloud base.
            ! if CCN is more than below parcel, nucleation newly occurs
            ! effective vertical velocity
            if ( flag_nucleation(k,i,j) .AND. & ! large scale upward motion region and saturated
                 tem(k,i,j) > tem_ccn_low ) then
               dlcdt_max    = (rhoq(I_QV,k,i,j) - esw(k,i,j)/(Rvap*tem(k,i,j)))*rdt
               dncdt_max    = dlcdt_max/xc_min
!               dnc_new      = nc_new(k,i,j)-rhoq(I_NC,k,i,j)
               dnc_new      = nc_new(k,i,j)
               PQ(I_NCccn,k,i,j) = min( dncdt_max, dnc_new*rdt )
               PQ(I_LCccn,k,i,j) = min( dlcdt_max, xc_min*PQ(I_NCccn,k,i,j) )
            else
               PQ(I_NCccn,k,i,j) = 0.0_RP
               PQ(I_LCccn,k,i,j) = 0.0_RP
            end if
         end do
         end do
         end do

    else

      if( nucl_twomey ) then
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            ! nucleation occurs at only cloud base.
            ! if CCN is more than below parcel, nucleation newly occurs
            ! effective vertical velocity
            if ( flag_nucleation(k,i,j) .AND. & ! large scale upward motion region and saturated
                 tem(k,i,j) > tem_ccn_low .AND. &
                 nc_new(k,i,j) > rhoq(I_NC,k,i,j) ) then
               dlcdt_max    = (rhoq(I_QV,k,i,j) - esw(k,i,j)/(Rvap*tem(k,i,j)))*rdt
               dncdt_max    = dlcdt_max/xc_min
               dnc_new      = nc_new(k,i,j)-rhoq(I_NC,k,i,j)
               PQ(I_NCccn,k,i,j) = min( dncdt_max, dnc_new*rdt )
               PQ(I_LCccn,k,i,j) = min( dlcdt_max, xc_min*PQ(I_NCccn,k,i,j) )
            else
               PQ(I_NCccn,k,i,j) = 0.0_RP
               PQ(I_LCccn,k,i,j) = 0.0_RP
            end if
         end do
         end do
         end do
      else
         do j = JS, JE
         do i = IS, IE
         do k = KS, KE
            ! effective vertical velocity
            if(  tem(k,i,j) > tem_ccn_low .AND. &
                 nc_new(k,i,j) > rhoq(I_NC,k,i,j) ) then
               dlcdt_max    = (rhoq(I_QV,k,i,j) - esw(k,i,j)/(Rvap*tem(k,i,j)))*rdt
               dncdt_max    = dlcdt_max/xc_min
               dnc_new      = nc_new(k,i,j)-rhoq(I_NC,k,i,j)
               PQ(I_NCccn,k,i,j) = min( dncdt_max, dnc_new*rdt )
               PQ(I_LCccn,k,i,j) = min( dlcdt_max, xc_min*PQ(I_NCccn,k,i,j) )
            else
               PQ(I_NCccn,k,i,j) = 0.0_RP
               PQ(I_LCccn,k,i,j) = 0.0_RP
            end if
         end do
         end do
         end do
      endif
    endif

    !
    ! ice nucleation
    !
    ! +++ NOTE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Based on Phillips etal.(2006).
    ! However this approach doesn't diagnose Ni itself but diagnose tendency.
    ! Original approach adjust Ni instantaneously .
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dz             = z(k) - z_below(k,i,j)
       w_dssidz(k,i,j) = velz(k,i,j)*(ssi(k,i,j) - ssi_below(k,i,j))/dz ! 09/04/14 [Add] T.Mitsui
       dssidt_rad(k,i,j) = -rhoq(I_QV,k,i,j)/(rho(k,i,j)*qsi(k,i,j)*qsi(k,i,j))*dqsidtem_rho(k,i,j)*dTdt_rad(k,i,j)
       dli_max        = (rhoq(I_QV,k,i,j) - esi(k,i,j)/(Rvap*tem(k,i,j)))*rdt
       dni_max        = min( dli_max/xi_ccn, (in_max-rhoq(I_NI,k,i,j))*rdt )
       wdssi          = min( w_dssidz(k,i,j)+dssidt_rad(k,i,j), 0.01_RP)
       wssi           = min( ssi(k,i,j), ssi_max)
       ! SB06(34),(35)
       if(  ( wdssi    > eps        ) .AND. & !
            (tem(k,i,j) < 273.15_RP   ) .AND. & !
            (rhoq(I_NI,k,i,j)  < in_max     ) .AND. &
            (wssi      >= eps       ) )then   !
!             PNIccn(k,i,j) = min(dni_max, c_in_map(1,i,j)*bm_M92*nm_M92*0.3_RP*exp(0.3_RP*bm_M92*(wssi-0.1_RP))*wdssi)
          if( inucl_w ) then
             PQ(I_NIccn,k,i,j) = min(dni_max, c_in*bm_M92*nm_M92*0.3_RP*exp(0.3_RP*bm_M92*(wssi-0.1_RP))*wdssi)
          else
             PQ(I_NIccn,k,i,j) = min(dni_max, max(c_in*nm_M92*exp(0.3_RP*bm_M92*(wssi-0.1_RP) )-rhoq(I_NI,k,i,j),0.0_RP )*rdt )
          endif
          PQ(I_LIccn,k,i,j) = min(dli_max, PQ(I_NIccn,k,i,j)*xi_ccn )
          ! only for output
!             dni_ratio(k,i,j) = dssidt_rad(k,i,j)/( w_dssidz(k,i,j)+dssidt_rad(k,i,j) )
       else
          PQ(I_NIccn,k,i,j) = 0.0_RP
          PQ(I_LIccn,k,i,j) = 0.0_RP
       end if
    end do
    end do
    end do

    return
  end subroutine nucleation_kij
  !----------------------------
  subroutine ice_multiplication_kij( &
       PQ,           & ! out
       Pac,          & ! in
       tem, rhoq, xq ) ! in

    ! ice multiplication by splintering
    ! we consider Hallet-Mossop process
    use scale_specfunc, only: &
       gammafunc => SF_gamma
    use scale_tracer, only: &
       QA
    implicit none
    real(RP), intent(inout):: PQ(PQ_MAX,KA,IA,JA)
    !
    real(RP), intent(in) :: Pac(Pac_MAX,KA,IA,JA)
    real(RP), intent(in) :: tem(KA,IA,JA)
    real(RP), intent(in) :: rhoq(I_QV:I_NG,KA,IA,JA)
    real(RP), intent(in) :: xq(HYDRO_MAX,KA,IA,JA)
    !
    logical, save :: flag_first      = .true.
    ! production of (350.d3 ice particles)/(cloud riming [g]) => 350*d6 [/kg]
    real(RP), parameter :: pice = 350.0E+6_RP
    ! production of (1 ice particle)/(250 cloud particles riming)
    real(RP), parameter :: pnc  = 250.0_RP
    real(RP), parameter :: rc_cr= 12.E-6_RP ! critical size[micron]
    ! temperature factor
    real(RP) :: fp
    ! work for incomplete gamma function
    real(RP), save :: xc_cr    ! mass[kg] of cloud with r=critical size[micron]
    real(RP), save :: alpha    ! slope parameter of gamma function
    real(RP), save :: gm, lgm  ! gamma(alpha), log(gamma(alpha))
    real(RP) :: igm            ! in complete gamma(x,alpha)
    real(RP) :: x
    ! coefficient of expansion using in calculation of igm
    real(RP) :: a0,a1,a2,a3,a4,a5
    real(RP) :: a6,a7,a8,a9,a10
    real(RP) :: an1,an2,b0,b1,b2,c0,c1,c2
    real(RP) :: d0,d1,d2,e1,e2,h0,h1,h2
    real(RP), parameter :: eps=1.0E-30_RP
    ! number of cloud droplets larger than 12 micron(radius).
    real(RP) :: n12
    !
    real(RP) :: wn, wni, wns, wng
    integer :: i, j, k
    !
    if( flag_first )then
       flag_first = .false.
       ! work for Incomplete Gamma function
       xc_cr = (2.0_RP*rc_cr/a_m(I_mp_QC))**(1.0_RP/b_m(I_mp_QC))
       alpha = (nu(I_mp_QC)+1.0_RP)/mu(I_mp_QC)
       gm    = gammafunc(alpha)
       lgm   = log(gm)
    end if
    !
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! Here we assume particle temperature is same as environment temperature.
       ! If you want to treat in a better manner,
       ! you can diagnose with eq.(64) in CT(86)
       if     (tem(k,i,j) >  270.16_RP)then
          fp   = 0.0_RP
       else if(tem(k,i,j) >= 268.16_RP)then
          fp   = (270.16_RP-tem(k,i,j))*0.5_RP
       else if(tem(k,i,j) >= 265.16_RP)then
          fp   = (tem(k,i,j)-265.16_RP)*0.333333333_RP
       else
          fp   = 0.0_RP
       end if
       ! Approximation of Incomplete Gamma function
       ! Here we calculate with algorithm by Numerical Recipes.
       ! This approach is based on eq.(78) in Cotton etal.(1986),
       ! but more accurate and expanded for Generalized Gamma Distribution.
       x       = coef_lambda(I_mp_QC)*(xc_cr/xq(I_mp_QC,k,i,j))**mu(I_mp_QC)
       !
       if(x<1.E-2_RP*alpha)then        ! negligible
          igm  = 0.0_RP
       else if(x<alpha+1.0_RP)then    ! series expansion
          ! 10th-truncation is enough for cloud droplet.
          a0   = 1.0_RP/alpha         ! n=0
          a1   = a0*x/(alpha+1.0_RP)  ! n=1
          a2   = a1*x/(alpha+2.0_RP)  ! n=2
          a3   = a2*x/(alpha+3.0_RP)  ! n=3
          a4   = a3*x/(alpha+4.0_RP)  ! n=4
          a5   = a4*x/(alpha+5.0_RP)  ! n=5
          a6   = a5*x/(alpha+6.0_RP)  ! n=6
          a7   = a6*x/(alpha+7.0_RP)  ! n=7
          a8   = a7*x/(alpha+8.0_RP)  ! n=8
          a9   = a8*x/(alpha+9.0_RP)  ! n=9
          a10  = a9*x/(alpha+10.0_RP) ! n=10
          igm  = (a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10)*exp( -x + alpha*log(x) - lgm )
       else if(x<alpha*1.d2) then ! continued fraction expansion
          ! 2nd-truncation is enough for cloud droplet.
          ! setup
          b0   = x+1.0_RP-alpha
          c0   = 1.0_RP/eps
          d0   = 1.0_RP/b0
          h0   = d0
          ! n=1
          an1  = -(1.0_RP-alpha)
          b1   = b0 + 2.0_RP
          d1   = 1.0_RP/(an1*d0+b1)
          c1   = b1+an1/c0
          e1   = d1*c1
          h1   = h0*e1
          ! n=2
          an2  = -2.0_RP*(2.0_RP-alpha)
          b2   = b1 + 2.0_RP
          d2   = 1.0_RP/(an2*d1+b2)
          c2   = b2+an2/c1
          e2   = d2*c2
          h2   = h1*e2
          !
          igm  = 1.0_RP - exp( -x + alpha*log(x) - lgm )*h2
       else                       ! negligible
          igm  = 1.0_RP
       end if
       ! n12 is number of cloud droplets larger than 12 micron.
       n12 = rhoq(I_NC,k,i,j)*(1.0_RP-igm)
       ! eq.(82) CT(86)
       wn           = (pice + n12/((rhoq(I_QC,k,i,j)+xc_min)*pnc) )*fp ! filtered by xc_min
       wni          = wn*(-Pac(I_LIacLC2LI,k,i,j)  ) ! riming production rate is all negative
       wns          = wn*(-Pac(I_LSacLC2LS,k,i,j)  )
       wng          = wn*(-Pac(I_LGacLC2LG,k,i,j)  )
       PQ(I_NIspl,k,i,j) = wni+wns+wng
       !
       PQ(I_LSspl,k,i,j) = - wns*xq(I_mp_QI,k,i,j) ! snow    => ice
       PQ(I_LGspl,k,i,j) = - wng*xq(I_mp_QI,k,i,j) ! graupel => ice
       !
    end do
    end do
    end do
    !
    return
  end subroutine ice_multiplication_kij
  !----------------------------
  subroutine mixed_phase_collection_kij( &
       ! collection process
       Pac, PQ,                          & ! out
       wtem, rhoq,                       & ! in
       xq, dq_xave,  vt_xave             ) ! in
!       rho                               ) ! [Add] 11/08/30
    use scale_tracer, only: &
       QA
    use scale_atmos_saturation, only: &
       moist_psat_ice => ATMOS_SATURATION_psat_ice
    implicit none

    !--- mixed-phase collection process
    !                  And all we set all production term as a negative sign to avoid confusion.
    !
    real(RP), intent(out):: Pac(Pac_MAX,KA,IA,JA)
    !--- partial conversion
    real(RP), intent(inout):: PQ(PQ_MAX,KA,IA,JA)
    !
    real(RP), intent(in) :: wtem(KA,IA,JA)
    !--- mass/number concentration[kg/m3]
    real(RP), intent(in) :: rhoq(I_QV:I_NG,KA,IA,JA)
    ! necessary ?
    real(RP), intent(in) :: xq(HYDRO_MAX,KA,IA,JA)
    !--- diameter of averaged mass( D(ave x) )
    real(RP), intent(in) :: dq_xave(HYDRO_MAX,KA,IA,JA)
    !--- terminal velocity of averaged mass( vt(ave x) )
    real(RP), intent(in) :: vt_xave(HYDRO_MAX,2,KA,IA,JA)
    ! [Add] 11/08/30 T.Mitsui, for autoconversion of ice
!    real(RP), intent(in) :: rho(KA,IA,JA)
    !
    ! namelist variables
    !=== for collection
    !--- threshold of diameters to collide with others
    real(RP), save :: dc0 =  15.0E-6_RP       ! lower threshold of cloud
    real(RP), save :: dc1 =  40.0E-6_RP       ! upper threshold of cloud
    real(RP), save :: di0 = 150.0E-6_RP       ! lower threshold of cloud
    real(RP), save :: ds0 = 150.0E-6_RP       ! lower threshold of cloud
    real(RP), save :: dg0 = 150.0E-6_RP       ! lower threshold of cloud
    !--- standard deviation of terminal velocity[m/s]
    real(RP), save :: sigma_c=0.00_RP        ! cloud
    real(RP), save :: sigma_r=0.00_RP        ! rain
    real(RP), save :: sigma_i=0.2_RP        ! ice
    real(RP), save :: sigma_s=0.2_RP        ! snow
    real(RP), save :: sigma_g=0.00_RP        ! graupel
    !--- max collection efficiency for cloud
    real(RP), save :: E_im = 0.80_RP        ! ice max
    real(RP), save :: E_sm = 0.80_RP        ! snow max
    real(RP), save :: E_gm = 1.00_RP        ! graupel max
    !--- collection efficiency between 2 species
    real(RP), save :: E_ir=1.0_RP            ! ice     x rain
    real(RP), save :: E_sr=1.0_RP            ! snow    x rain
    real(RP), save :: E_gr=1.0_RP            ! graupel x rain
    real(RP), save :: E_ii=1.0_RP            ! ice     x ice
    real(RP), save :: E_si=1.0_RP            ! snow    x ice
    real(RP), save :: E_gi=1.0_RP            ! graupel x ice
    real(RP), save :: E_ss=1.0_RP            ! snow    x snow
    real(RP), save :: E_gs=1.0_RP            ! graupel x snow
    real(RP), save :: E_gg=1.0_RP            ! graupel x graupel
    !=== for partial conversion
    !--- flag: 1=> partial conversion to graupel, 0=> no conversion
    integer, save :: i_iconv2g=1          ! ice  => graupel
    integer, save :: i_sconv2g=1          ! snow => graupel
    !--- bulk density of graupel
    real(RP), save :: rho_g   = 900.0_RP     ! [kg/m3]
    !--- space filling coefficient [%]
    real(RP), save :: cfill_i = 0.68_RP     ! ice
    real(RP), save :: cfill_s = 0.01_RP     ! snow
    !--- critical diameter for ice conversion
    real(RP), save :: di_cri  = 500.E-6_RP    ! [m]
    ! [Add] 10/08/03 T.Mitsui
    logical, save :: opt_stick_KS96=.false.
    logical, save :: opt_stick_CO86=.false.
    real(RP), parameter :: a_dec = 0.883_RP
    real(RP), parameter :: b_dec = 0.093_RP
    real(RP), parameter :: c_dec = 0.00348_RP
    real(RP), parameter :: d_dec = 4.5185E-5_RP
    !
    logical, save :: flag_first = .true.
    namelist /nm_mp_sn14_collection/ &
         dc0, dc1, di0, ds0, dg0,    &
         sigma_c, sigma_r, sigma_i, sigma_s, sigma_g, &
         opt_stick_KS96,   &
         opt_stick_CO86,   &
         E_im, E_sm, E_gm, &
         E_ir, E_sr, E_gr, E_ii, E_si, E_gi, E_ss, E_gs, E_gg, &
         i_iconv2g, i_sconv2g, rho_g, cfill_i, cfill_s, di_cri
    !
    real(RP) :: tem(KA,IA,JA)
    !
    !--- collection efficency of each specie
    real(RP) :: E_c, E_r, E_i, E_s, E_g    !
    real(RP) :: E_ic, E_sc, E_gc           !
    !--- sticking efficiency
    real(RP) :: E_stick(KA,IA,JA)
    ! [Add] 10/08/03 T.Mitsui
    real(RP) :: temc, temc2, temc3
    real(RP) :: E_dec
    real(RP) :: esi_rat
    real(RP) :: esi(KA,IA,JA)
    !
    real(RP) :: temc_p, temc_m             ! celcius tem.
    ! [Add] 11/08/30 T.Mitsui, estimation of autoconversion time
!    real(RP) :: ci_aut(KA,IA,JA)
!    real(RP) :: taui_aut(KA,IA,JA)
!    real(RP) :: tau_sce(KA,IA,JA)
    !--- DSD averaged diameter for each species
    real(RP) :: ave_dc                     ! cloud
!    real(RP) :: ave_dr                     ! rain
    real(RP) :: ave_di                     ! ice
    real(RP) :: ave_ds                     ! snow
    real(RP) :: ave_dg                     ! graupel
    !--- coefficient of collection equations(L:mass, N:number)
    real(RP) :: coef_acc_LCI,   coef_acc_NCI   ! cloud     - cloud ice
    real(RP) :: coef_acc_LCS,   coef_acc_NCS   ! cloud     - snow
    !
    real(RP) :: coef_acc_LCG,   coef_acc_NCG   ! cloud     - graupel
    real(RP) :: coef_acc_LRI_I, coef_acc_NRI_I ! rain      - cloud ice
    real(RP) :: coef_acc_LRI_R, coef_acc_NRI_R ! rain      - cloud ice
    real(RP) :: coef_acc_LRS_S, coef_acc_NRS_S ! rain      - snow
    real(RP) :: coef_acc_LRS_R, coef_acc_NRS_R ! rain      - snow
    real(RP) :: coef_acc_LRG,   coef_acc_NRG   ! rain      - graupel
    real(RP) :: coef_acc_LII,   coef_acc_NII   ! cloud ice - cloud ice
    real(RP) :: coef_acc_LIS,   coef_acc_NIS   ! cloud ice - snow
    real(RP) ::                 coef_acc_NSS   ! snow      - snow
    real(RP) ::                 coef_acc_NGG   ! grauepl   - graupel
    real(RP) :: coef_acc_LSG,   coef_acc_NSG   ! snow      - graupel
    !--- (diameter) x (diameter)
    real(RP) :: dcdc, dcdi, dcds, dcdg
    real(RP) :: drdr, drdi, drds, drdg
    real(RP) :: didi, dids, didg
    real(RP) :: dsds, dsdg
    real(RP) :: dgdg
    !--- (terminal velocity) x (terminal velocity)
    real(RP) :: vcvc, vcvi, vcvs, vcvg
    real(RP) :: vrvr, vrvi, vrvs, vrvg
    real(RP) :: vivi, vivs, vivg
    real(RP) :: vsvs, vsvg
    real(RP) :: vgvg
    !
    real(RP) :: wx_cri, wx_crs
    real(RP) :: coef_emelt
    real(RP) :: w1

    real(RP) :: sw
    !
    integer :: i, j, k
    !
    if( flag_first )then
       rewind( IO_FID_CONF )
       read( IO_FID_CONF, nml=nm_mp_sn14_collection, end=100 )
100    if( IO_NML ) write(IO_FID_NML,nml=nm_mp_sn14_collection)
       flag_first = .false.
    end if
    !
    ! [Add] 10/08/03 T.Mitsui
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       tem(k,i,j) = max( wtem(k,i,j), tem_min ) ! 11/08/30 T.Mitsui
    end do
    end do
    end do

    call moist_psat_ice( esi, tem )

    if( opt_stick_KS96 )then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ! Khain and Sednev (1996), eq.(3.15)
          temc          = tem(k,i,j) - T00
          temc2         = temc*temc
          temc3         = temc2*temc
          E_dec         = max(0.0_RP, a_dec + b_dec*temc + c_dec*temc2 + d_dec*temc3 )
          esi_rat       = rhoq(I_QV,k,i,j)*Rvap*tem(k,i,j)/esi(k,i,j)
          E_stick(k,i,j) = min(1.0_RP, E_dec*esi_rat)
       end do
       end do
       end do
    else if( opt_stick_CO86 )then
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ! [Add] 11/08/30 T.Mitsui, Cotton et al. (1986)
          temc          = min(tem(k,i,j) - T00,0.0_RP)
          w1            = 0.035_RP*temc-0.7_RP
          E_stick(k,i,j) = 10._RP**w1
       end do
       end do
       end do
    else
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ! Lin et al. (1983)
          temc_m        = min(tem(k,i,j) - T00,0.0_RP) ! T < 273.15
          E_stick(k,i,j) = exp(0.09_RP*temc_m)
       end do
       end do
       end do
    end if
    !
    PROFILE_START("sn14_collection")
!OCL NOSIMD
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       !
!          temc_m = min(tem(k,i,j) - T00,0.0_RP) ! T < 273.15
       temc_p = max(tem(k,i,j) - T00,0.0_RP) ! T > 273.15
       ! averaged diameter using SB06(82)
       ave_dc = coef_d(I_mp_QC)*xq(I_mp_QC,k,i,j)**b_m(I_mp_QC)
       ave_di = coef_d(I_mp_QI)*xq(I_mp_QI,k,i,j)**b_m(I_mp_QI)
       ave_ds = coef_d(I_mp_QS)*xq(I_mp_QS,k,i,j)**b_m(I_mp_QS)
       ave_dg = coef_d(I_mp_QG)*xq(I_mp_QG,k,i,j)**b_m(I_mp_QG)
       !------------------------------------------------------------------------
       ! coellection efficiency are given as follows
       E_c = max(0.0_RP, min(1.0_RP, (ave_dc-dc0)/(dc1-dc0) ))
       sw = 0.5_RP - sign(0.5_RP, di0-ave_di) ! if(ave_di>di0)then sw=1
       E_i = E_im * sw
       sw = 0.5_RP - sign(0.5_RP, ds0-ave_ds) ! if(ave_ds>ds0)then sw=1
       E_s = E_sm * sw
       sw = 0.5_RP - sign(0.5_RP, dg0-ave_dg) ! if(ave_dg>dg0)then sw=1
       E_g = E_gm * sw
       E_ic = E_i*E_c
       E_sc = E_s*E_c
       E_gc = E_g*E_c
       !------------------------------------------------------------------------
       ! Collection:  a collects b ( assuming particle size a>b )
       dcdc = dq_xave(I_mp_QC,k,i,j) * dq_xave(I_mp_QC,k,i,j)
       drdr = dq_xave(I_mp_QR,k,i,j) * dq_xave(I_mp_QR,k,i,j)
       didi = dq_xave(I_mp_QI,k,i,j) * dq_xave(I_mp_QI,k,i,j)
       dsds = dq_xave(I_mp_QS,k,i,j) * dq_xave(I_mp_QS,k,i,j)
       dgdg = dq_xave(I_mp_QG,k,i,j) * dq_xave(I_mp_QG,k,i,j)
       dcdi = dq_xave(I_mp_QC,k,i,j) * dq_xave(I_mp_QI,k,i,j)
       dcds = dq_xave(I_mp_QC,k,i,j) * dq_xave(I_mp_QS,k,i,j)
       dcdg = dq_xave(I_mp_QC,k,i,j) * dq_xave(I_mp_QG,k,i,j)
       drdi = dq_xave(I_mp_QR,k,i,j) * dq_xave(I_mp_QI,k,i,j)
       drds = dq_xave(I_mp_QR,k,i,j) * dq_xave(I_mp_QS,k,i,j)
       drdg = dq_xave(I_mp_QR,k,i,j) * dq_xave(I_mp_QG,k,i,j)
       dids = dq_xave(I_mp_QI,k,i,j) * dq_xave(I_mp_QS,k,i,j)
!       didg = dq_xave(I_mp_QI,k,i,j) * dq_xave(I_mp_QG,k,i,j)
       dsdg = dq_xave(I_mp_QS,k,i,j) * dq_xave(I_mp_QG,k,i,j)
       !
       vcvc = vt_xave(I_mp_QC,2,k,i,j)* vt_xave(I_mp_QC,2,k,i,j)
       vrvr = vt_xave(I_mp_QR,2,k,i,j)* vt_xave(I_mp_QR,2,k,i,j)
       vivi = vt_xave(I_mp_QI,2,k,i,j)* vt_xave(I_mp_QI,2,k,i,j)
       vsvs = vt_xave(I_mp_QS,2,k,i,j)* vt_xave(I_mp_QS,2,k,i,j)
       vgvg = vt_xave(I_mp_QG,2,k,i,j)* vt_xave(I_mp_QG,2,k,i,j)
       vcvi = vt_xave(I_mp_QC,2,k,i,j)* vt_xave(I_mp_QI,2,k,i,j)
       vcvs = vt_xave(I_mp_QC,2,k,i,j)* vt_xave(I_mp_QS,2,k,i,j)
       vcvg = vt_xave(I_mp_QC,2,k,i,j)* vt_xave(I_mp_QG,2,k,i,j)
       vrvi = vt_xave(I_mp_QR,2,k,i,j)* vt_xave(I_mp_QI,2,k,i,j)
       vrvs = vt_xave(I_mp_QR,2,k,i,j)* vt_xave(I_mp_QS,2,k,i,j)
       vrvg = vt_xave(I_mp_QR,2,k,i,j)* vt_xave(I_mp_QG,2,k,i,j)
       vivs = vt_xave(I_mp_QI,2,k,i,j)* vt_xave(I_mp_QS,2,k,i,j)
!       vivg = vt_xave(I_mp_QI,2,k,i,j)* vt_xave(I_mp_QG,2,k,i,j)
       vsvg = vt_xave(I_mp_QS,2,k,i,j)* vt_xave(I_mp_QG,2,k,i,j)
       !------------------------------------------------------------------------
       !
       !+++ pattern 1: a + b => a  (a>b)
       !                           (i-c, s-c, g-c, s-i, g-r, s-g)
       !------------------------------------------------------------------------
       ! cloud-ice => ice
       ! reduction term of cloud
       coef_acc_LCI = &
              ( delta_b1(I_mp_QC)*dcdc + delta_ab1(I_mp_QI,I_mp_QC)*dcdi + delta_b0(I_mp_QI)*didi ) &
            * ( theta_b1(I_mp_QC)*vcvc - theta_ab1(I_mp_QI,I_mp_QC)*vcvi + theta_b0(I_mp_QI)*vivi   &
            +  sigma_i + sigma_c )**0.5_RP
       coef_acc_NCI = &
              ( delta_b0(I_mp_QC)*dcdc + delta_ab0(I_mp_QI,I_mp_QC)*dcdi + delta_b0(I_mp_QI)*didi ) &
            * ( theta_b0(I_mp_QC)*vcvc - theta_ab0(I_mp_QI,I_mp_QC)*vcvi + theta_b0(I_mp_QI)*vivi   &
            +  sigma_i + sigma_c )**0.5_RP
       Pac(I_LIacLC2LI,k,i,j)= -0.25_RP*pi*E_ic*rhoq(I_NI,k,i,j)*rhoq(I_QC,k,i,j)*coef_acc_LCI
       Pac(I_NIacNC2NI,k,i,j)= -0.25_RP*pi*E_ic*rhoq(I_NI,k,i,j)*rhoq(I_NC,k,i,j)*coef_acc_NCI
       ! cloud-snow => snow
       ! reduction term of cloud
       coef_acc_LCS = &
              ( delta_b1(I_mp_QC)*dcdc + delta_ab1(I_mp_QS,I_mp_QC)*dcds + delta_b0(I_mp_QS)*dsds ) &
            * ( theta_b1(I_mp_QC)*vcvc - theta_ab1(I_mp_QS,I_mp_QC)*vcvs + theta_b0(I_mp_QS)*vsvs   &
            +  sigma_s + sigma_c )**0.5_RP
       coef_acc_NCS = &
              ( delta_b0(I_mp_QC)*dcdc + delta_ab0(I_mp_QS,I_mp_QC)*dcds + delta_b0(I_mp_QS)*dsds ) &
            * ( theta_b0(I_mp_QC)*vcvc - theta_ab0(I_mp_QS,I_mp_QC)*vcvs + theta_b0(I_mp_QS)*vsvs   &
            +  sigma_s + sigma_c )**0.5_RP
       Pac(I_LSacLC2LS,k,i,j)= -0.25_RP*pi*E_sc*rhoq(I_NS,k,i,j)*rhoq(I_QC,k,i,j)*coef_acc_LCS
       Pac(I_NSacNC2NS,k,i,j)= -0.25_RP*pi*E_sc*rhoq(I_NS,k,i,j)*rhoq(I_NC,k,i,j)*coef_acc_NCS
       ! cloud-graupel => graupel
       ! reduction term of cloud
       coef_acc_LCG = &
              ( delta_b1(I_mp_QC)*dcdc + delta_ab1(I_mp_QG,I_mp_QC)*dcdg + delta_b0(I_mp_QG)*dgdg ) &
            * ( theta_b1(I_mp_QC)*vcvc - theta_ab1(I_mp_QG,I_mp_QC)*vcvg + theta_b0(I_mp_QG)*vgvg   &
            +  sigma_g + sigma_c )**0.5_RP
       coef_acc_NCG = &
              ( delta_b0(I_mp_QC)*dcdc + delta_ab0(I_mp_QG,I_mp_QC)*dcdg + delta_b0(I_mp_QG)*dgdg ) &
            * ( theta_b0(I_mp_QC)*vcvc - theta_ab0(I_mp_QG,I_mp_QC)*vcvg + theta_b0(I_mp_QG)*vgvg   &
            +  sigma_g + sigma_c )**0.5_RP
       Pac(I_LGacLC2LG,k,i,j)= -0.25_RP*pi*E_gc*rhoq(I_NG,k,i,j)*rhoq(I_QC,k,i,j)*coef_acc_LCG
       Pac(I_NGacNC2NG,k,i,j)= -0.25_RP*pi*E_gc*rhoq(I_NG,k,i,j)*rhoq(I_NC,k,i,j)*coef_acc_NCG
       ! snow-graupel => graupel
       coef_acc_LSG = &
              ( delta_b1(I_mp_QS)*dsds + delta_ab1(I_mp_QG,I_mp_QS)*dsdg + delta_b0(I_mp_QG)*dgdg ) &
            * ( theta_b1(I_mp_QS)*vsvs - theta_ab1(I_mp_QG,I_mp_QS)*vsvg + theta_b0(I_mp_QG)*vgvg   &
            +  sigma_g + sigma_s )**0.5_RP
       coef_acc_NSG = &
              ( delta_b0(I_mp_QS)*dsds + delta_ab0(I_mp_QG,I_mp_QS)*dsdg + delta_b0(I_mp_QG)*dgdg ) &
            ! [fix] T.Mitsui 08/05/08
            * ( theta_b0(I_mp_QS)*vsvs - theta_ab0(I_mp_QG,I_mp_QS)*vsvg + theta_b0(I_mp_QG)*vgvg   &
            +  sigma_g + sigma_s )**0.5_RP
       Pac(I_LGacLS2LG,k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_gs*rhoq(I_NG,k,i,j)*rhoq(I_QS,k,i,j)*coef_acc_LSG
       Pac(I_NGacNS2NG,k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_gs*rhoq(I_NG,k,i,j)*rhoq(I_NS,k,i,j)*coef_acc_NSG
       !------------------------------------------------------------------------
       ! ice-snow => snow
       ! reduction term of ice
       coef_acc_LIS = &
              ( delta_b1(I_mp_QI)*didi + delta_ab1(I_mp_QS,I_mp_QI)*dids + delta_b0(I_mp_QS)*dsds ) &
            * ( theta_b1(I_mp_QI)*vivi - theta_ab1(I_mp_QS,I_mp_QI)*vivs + theta_b0(I_mp_QS)*vsvs   &
            +  sigma_i + sigma_s )**0.5_RP
       coef_acc_NIS = &
              ( delta_b0(I_mp_QI)*didi + delta_ab0(I_mp_QS,I_mp_QI)*dids + delta_b0(I_mp_QS)*dsds ) &
            * ( theta_b0(I_mp_QI)*vivi - theta_ab0(I_mp_QS,I_mp_QI)*vivs + theta_b0(I_mp_QS)*vsvs   &
            +  sigma_i + sigma_s )**0.5_RP
       Pac(I_LIacLS2LS,k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_si*rhoq(I_NS,k,i,j)*rhoq(I_QI,k,i,j)*coef_acc_LIS
       Pac(I_NIacNS2NS,k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_si*rhoq(I_NS,k,i,j)*rhoq(I_NI,k,i,j)*coef_acc_NIS
       !
       sw = sign(0.5_RP, T00-tem(k,i,j)) + 0.5_RP
       ! if ( tem(k,i,j) <= T00 )then
          ! rain-graupel => graupel
          ! reduction term of rain
          ! sw = 1
       ! else
          ! rain-graupel => rain
          ! reduction term of graupel
          ! sw = 0
       coef_acc_LRG = &
              ( ( delta_b1(I_mp_QR)*drdr + delta_ab1(I_mp_QG,I_mp_QR)*drdg + delta_b0(I_mp_QG)*dgdg ) * sw &
              + ( delta_b1(I_mp_QG)*dgdg + delta_ab1(I_mp_QR,I_mp_QG)*drdg + delta_b0(I_mp_QR)*drdr ) * (1.0_RP-sw) ) &
            * sqrt( ( theta_b1(I_mp_QR)*vrvr - theta_ab1(I_mp_QG,I_mp_QR)*vrvg + theta_b0(I_mp_QG)*vgvg ) * sw &
                  + ( theta_b1(I_mp_QG)*vgvg - theta_ab1(I_mp_QR,I_mp_QG)*vrvg + theta_b0(I_mp_QR)*vrvr ) * (1.0_RP-sw) &
                  + sigma_r + sigma_g )
       Pac(I_LRacLG2LG,k,i,j) = -0.25_RP*pi*E_gr*coef_acc_LRG &
            * ( rhoq(I_NG,k,i,j)*rhoq(I_QR,k,i,j) * sw &
              + rhoq(I_NR,k,i,j)*rhoq(I_QG,k,i,j) * (1.0_RP-sw) )
       coef_acc_NRG = &
              ( delta_b0(I_mp_QR)*drdr + delta_ab0(I_mp_QG,I_mp_QR)*drdg + delta_b0(I_mp_QG)*dgdg ) &
            * ( theta_b0(I_mp_QR)*vrvr - theta_ab0(I_mp_QG,I_mp_QR)*vrvg + theta_b0(I_mp_QG)*vgvg   &
            +  sigma_r + sigma_g )**0.5_RP
       Pac(I_NRacNG2NG,k,i,j) = -0.25_RP*pi*E_gr*rhoq(I_NG,k,i,j)*rhoq(I_NR,k,i,j)*coef_acc_NRG
       !
       !------------------------------------------------------------------------
       !
       !+++ pattern 2: a + b => c  (a>b)
       !                           (r-i,r-s)
       !------------------------------------------------------------------------
       ! rain-ice => graupel
       ! reduction term of ice
       coef_acc_LRI_I = &
              ( delta_b1(I_mp_QI)*didi + delta_ab1(I_mp_QR,I_mp_QI)*drdi + delta_b0(I_mp_QR)*drdr ) &
            * ( theta_b1(I_mp_QI)*vivi - theta_ab1(I_mp_QR,I_mp_QI)*vrvi + theta_b0(I_mp_QR)*vrvr   &
            +  sigma_r + sigma_i )**0.5_RP
       coef_acc_NRI_I = &
              ( delta_b0(I_mp_QI)*didi + delta_ab0(I_mp_QR,I_mp_QI)*drdi + delta_b0(I_mp_QR)*drdr ) &
            * ( theta_b0(I_mp_QI)*vivi - theta_ab0(I_mp_QR,I_mp_QI)*vrvi + theta_b0(I_mp_QR)*vrvr   &
            +  sigma_r + sigma_i )**0.5_RP
       Pac(I_LRacLI2LG_I,k,i,j)= -0.25_RP*pi*E_ir*rhoq(I_NR,k,i,j)*rhoq(I_QI,k,i,j)*coef_acc_LRI_I
       Pac(I_NRacNI2NG_I,k,i,j)= -0.25_RP*pi*E_ir*rhoq(I_NR,k,i,j)*rhoq(I_NI,k,i,j)*coef_acc_NRI_I
       ! reduction term of rain
       coef_acc_LRI_R = &
              ( delta_b1(I_mp_QR)*drdr + delta_ab1(I_mp_QI,I_mp_QR)*drdi + delta_b0(I_mp_QI)*didi ) &
            * ( theta_b1(I_mp_QR)*vrvr - theta_ab1(I_mp_QI,I_mp_QR)*vrvi + theta_b0(I_mp_QI)*vivi   &
            +  sigma_r + sigma_i )**0.5_RP
       coef_acc_NRI_R = &
              ( delta_b0(I_mp_QR)*drdr + delta_ab0(I_mp_QI,I_mp_QR)*drdi + delta_b0(I_mp_QI)*didi ) &
            * ( theta_b0(I_mp_QR)*vrvr - theta_ab0(I_mp_QI,I_mp_QR)*vrvi + theta_b0(I_mp_QI)*vivi   &
            +  sigma_r + sigma_i )**0.5_RP
       Pac(I_LRacLI2LG_R,k,i,j)= -0.25_RP*pi*E_ir*rhoq(I_NI,k,i,j)*rhoq(I_QR,k,i,j)*coef_acc_LRI_R
       Pac(I_NRacNI2NG_R,k,i,j)= -0.25_RP*pi*E_ir*rhoq(I_NI,k,i,j)*rhoq(I_NR,k,i,j)*coef_acc_NRI_R
       ! rain-snow => graupel
       ! reduction term of snow
       coef_acc_LRS_S = &
              ( delta_b1(I_mp_QS)*dsds + delta_ab1(I_mp_QR,I_mp_QS)*drds + delta_b0(I_mp_QR)*drdr ) &
            * ( theta_b1(I_mp_QS)*vsvs - theta_ab1(I_mp_QR,I_mp_QS)*vrvs + theta_b0(I_mp_QR)*vrvr   &
            +  sigma_r + sigma_s )**0.5_RP
       coef_acc_NRS_S = &
              ( delta_b0(I_mp_QS)*dsds + delta_ab0(I_mp_QR,I_mp_QS)*drds + delta_b0(I_mp_QR)*drdr ) &
            * ( theta_b0(I_mp_QS)*vsvs - theta_ab0(I_mp_QR,I_mp_QS)*vrvs + theta_b0(I_mp_QR)*vrvr   &
            +  sigma_r + sigma_s )**0.5_RP
       Pac(I_LRacLS2LG_S,k,i,j)= -0.25_RP*pi*E_sr*rhoq(I_NR,k,i,j)*rhoq(I_QS,k,i,j)*coef_acc_LRS_S
       Pac(I_NRacNS2NG_S,k,i,j)= -0.25_RP*pi*E_sr*rhoq(I_NR,k,i,j)*rhoq(I_NS,k,i,j)*coef_acc_NRS_S
       ! reduction term of rain
       coef_acc_LRS_R = &
              ( delta_b1(I_mp_QR)*drdr + delta_ab1(I_mp_QS,I_mp_QR)*drds + delta_b0(I_mp_QS)*dsds ) &
            * ( theta_b1(I_mp_QR)*vrvr - theta_ab1(I_mp_QS,I_mp_QR)*vrvs + theta_b0(I_mp_QS)*vsvs   &
            +  sigma_r + sigma_s )**0.5_RP
       coef_acc_NRS_R = &
              ( delta_b0(I_mp_QR)*drdr + delta_ab0(I_mp_QS,I_mp_QR)*drds + delta_b0(I_mp_QS)*dsds ) &
            * ( theta_b0(I_mp_QR)*vrvr - theta_ab0(I_mp_QS,I_mp_QR)*vrvs + theta_b0(I_mp_QS)*vsvs   &
            +  sigma_r + sigma_s )**0.5_RP
       Pac(I_LRacLS2LG_R,k,i,j)= -0.25_RP*pi*E_sr*rhoq(I_NS,k,i,j)*rhoq(I_QR,k,i,j)*coef_acc_LRS_R
       Pac(I_NRacNS2NG_R,k,i,j)= -0.25_RP*pi*E_sr*rhoq(I_NS,k,i,j)*rhoq(I_NR,k,i,j)*coef_acc_NRS_R
       !------------------------------------------------------------------------
       !
       !+++ pattern 3: a + a => b  (i-i)
       !
       !------------------------------------------------------------------------
       ! ice-ice ( reduction is double, but includes double-count)
       coef_acc_LII = &
              ( delta_b0(I_mp_QI)*didi + delta_ab1(I_mp_QI,I_mp_QI)*didi + delta_b1(I_mp_QI)*didi ) &
            * ( theta_b0(I_mp_QI)*vivi - theta_ab1(I_mp_QI,I_mp_QI)*vivi + theta_b1(I_mp_QI)*vivi   &
            +  sigma_i + sigma_i )**0.5_RP
       coef_acc_NII = &
              ( delta_b0(I_mp_QI)*didi + delta_ab0(I_mp_QI,I_mp_QI)*didi + delta_b0(I_mp_QI)*didi ) &
            * ( theta_b0(I_mp_QI)*vivi - theta_ab0(I_mp_QI,I_mp_QI)*vivi + theta_b0(I_mp_QI)*vivi   &
            +  sigma_i + sigma_i )**0.5_RP
       Pac(I_LIacLI2LS,k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_ii*rhoq(I_NI,k,i,j)*rhoq(I_QI,k,i,j)*coef_acc_LII
       Pac(I_NIacNI2NS,k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_ii*rhoq(I_NI,k,i,j)*rhoq(I_NI,k,i,j)*coef_acc_NII
       !
!          ci_aut(k,i,j)   =  0.25_RP*pi*E_ii*rhoq(I_NI,k,i,j)*coef_acc_LII
!          taui_aut(k,i,j) = 1._RP/max(E_stick(k,i,j)*ci_aut(k,i,j),1.E-10_RP)
!          tau_sce(k,i,j)  = rhoq(I_QI,k,i,j)/max(rhoq(I_QI,k,i,j)+rhoq(I_QS,k,i,j),1.E-10_RP)
       !------------------------------------------------------------------------
       !
       !+++ pattern 4: a + a => a  (s-s)
       !
       !------------------------------------------------------------------------
       ! snow-snow => snow
       coef_acc_NSS = &
              ( delta_b0(I_mp_QS)*dsds + delta_ab0(I_mp_QS,I_mp_QS)*dsds + delta_b0(I_mp_QS)*dsds ) &
            * ( theta_b0(I_mp_QS)*vsvs - theta_ab0(I_mp_QS,I_mp_QS)*vsvs + theta_b0(I_mp_QS)*vsvs   &
            +  sigma_s + sigma_s )**0.5_RP
       Pac(I_NSacNS2NS,k,i,j)= -0.125_RP*pi*E_stick(k,i,j)*E_ss*rhoq(I_NS,k,i,j)*rhoq(I_NS,k,i,j)*coef_acc_NSS
       !
       ! graupel-grauple => graupel
       coef_acc_NGG = &
              ( delta_b0(I_mp_QG)*dgdg + delta_ab0(I_mp_QG,I_mp_QG)*dgdg + delta_b0(I_mp_QG)*dgdg ) &
            * ( theta_b0(I_mp_QG)*vgvg - theta_ab0(I_mp_QG,I_mp_QG)*vgvg + theta_b0(I_mp_QG)*vgvg   &
            +  sigma_g + sigma_g )**0.5_RP
       Pac(I_NGacNG2NG,k,i,j)= -0.125_RP*pi*E_stick(k,i,j)*E_gg*rhoq(I_NG,k,i,j)*rhoq(I_NG,k,i,j)*coef_acc_NGG
       !
       !------------------------------------------------------------------------
       !--- Partial conversion
       ! SB06(70),(71)
       ! i_iconv2g: option whether partial conversions work or not
       ! ice-cloud => graupel
       sw = 0.5_RP - sign(0.5_RP,di_cri-ave_di) ! if( ave_di > di_cri )then sw=1
       wx_cri = cfill_i*DWATR/rho_g*( pi/6.0_RP*rho_g*ave_di*ave_di*ave_di/xq(I_mp_QI,k,i,j) - 1.0_RP ) * sw
       PQ(I_LIcon,k,i,j) = i_iconv2g * Pac(I_LIacLC2LI,k,i,j)/max(1.0_RP, wx_cri) * sw
       PQ(I_NIcon,k,i,j) = i_iconv2g * PQ(I_LIcon,k,i,j)/xq(I_mp_QI,k,i,j) * sw

       ! snow-cloud => graupel
       wx_crs = cfill_s*DWATR/rho_g*( pi/6.0_RP*rho_g*ave_ds*ave_ds*ave_ds/xq(I_mp_QS,k,i,j) - 1.0_RP )
       PQ(I_LScon,k,i,j) = i_sconv2g * (Pac(I_LSacLC2LS,k,i,j))/max(1.0_RP, wx_crs)
       PQ(I_NScon,k,i,j) = i_sconv2g * PQ(I_LScon,k,i,j)/xq(I_mp_QS,k,i,j)
       !------------------------------------------------------------------------
       !--- enhanced melting( due to collection-freezing of water droplets )
       !    originally from Rutledge and Hobbs(1984). eq.(A.21)
       ! if T > 273.15 then temc_p=T-273.15, else temc_p=0
       ! 08/05/08 [fix] T.Mitsui LHF00 => LHF0
       ! melting occurs around T=273K, so LHF0 is suitable both SIMPLE and EXACT,
       ! otherwise LHF can have sign both negative(EXACT) and positive(SIMPLE).
!!$       coef_emelt   = -CL/LHF00*temc_p
       coef_emelt   =  CL/LHF0*temc_p
       ! cloud-graupel
       PQ(I_LGacm,k,i,j) =  coef_emelt*Pac(I_LGacLC2LG,k,i,j)
       PQ(I_NGacm,k,i,j) =  PQ(I_LGacm,k,i,j)/xq(I_mp_QG,k,i,j)
       ! rain-graupel
       PQ(I_LGarm,k,i,j) =  coef_emelt*Pac(I_LRacLG2LG,k,i,j)
       PQ(I_NGarm,k,i,j) =  PQ(I_LGarm,k,i,j)/xq(I_mp_QG,k,i,j)
       ! cloud-snow
       PQ(I_LSacm,k,i,j) =  coef_emelt*(Pac(I_LSacLC2LS,k,i,j))
       PQ(I_NSacm,k,i,j) =  PQ(I_LSacm,k,i,j)/xq(I_mp_QS,k,i,j)
       ! rain-snow
       PQ(I_LSarm,k,i,j) =  coef_emelt*(Pac(I_LRacLS2LG_R,k,i,j)+Pac(I_LRacLS2LG_S,k,i,j))
       PQ(I_NSarm,k,i,j) =  PQ(I_LSarm,k,i,j)/xq(I_mp_QG,k,i,j) ! collect? might be I_mp_QS
       ! cloud-ice
       PQ(I_LIacm,k,i,j) =  coef_emelt*Pac(I_LIacLC2LI,k,i,j)
       PQ(I_NIacm,k,i,j) =  PQ(I_LIacm,k,i,j)/xq(I_mp_QI,k,i,j)
       ! rain-ice
       PQ(I_LIarm,k,i,j) =  coef_emelt*(Pac(I_LRacLI2LG_R,k,i,j)+Pac(I_LRacLI2LG_I,k,i,j))
       PQ(I_NIarm,k,i,j) =  PQ(I_LIarm,k,i,j)/xq(I_mp_QG,k,i,j) ! collect? might be I_mp_QI
    end do
    end do
    end do
    PROFILE_STOP("sn14_collection")

    !
    return
  end subroutine mixed_phase_collection_kij
  !----------------------------
  ! Auto-conversion, Accretion, Self-collection, Break-up
  subroutine aut_acc_slc_brk_kij(  &
       PQ,                     &
       rhoq, xq, dq_xave,      &
       rho                     )
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(inout) :: PQ(PQ_MAX,KA,IA,JA)
    !
    real(RP), intent(in)  :: rhoq(I_QV:I_NG,KA,IA,JA)
    real(RP), intent(in)  :: xq(HYDRO_MAX,KA,IA,JA)
    real(RP), intent(in)  :: dq_xave(HYDRO_MAX,KA,IA,JA)
    real(RP), intent(in)  :: rho(KA,IA,JA)
    !
    ! parameter for autoconversion
    real(RP), parameter :: kcc     = 4.44E+9_RP  ! collision efficiency [m3/kg2/sec]
    real(RP), parameter :: tau_min = 1.E-20_RP  ! empirical filter by T.Mitsui
    real(RP), parameter :: rx_sep  = 1.0_RP/x_sep ! 1/x_sep, 10/08/03 [Add] T.Mitsui
    !
    ! parameter for accretion
    real(RP), parameter :: kcr     = 5.8_RP     ! collision efficiency [m3/kg2/sec]
    real(RP), parameter :: thr_acc = 5.E-5_RP   ! threshold for universal function original
    !
    ! parameter for self collection and collison break-up
    real(RP), parameter :: krr     = 4.33_RP ! k_rr,      S08 (35)
    real(RP), parameter :: kaprr   = 60.7_RP  ! kappa_rr,  SB06(11)
    real(RP), parameter :: kbr     = 1000._RP ! k_br,      SB06(14)
    real(RP), parameter :: kapbr   = 2.3E+3_RP   ! kappa_br,  SB06(15)
    real(RP), parameter :: dr_min  = 0.35E-3_RP ! minimum diameter, SB06(13)-(15)
    !
    ! work variables
    real(RP) :: coef_nuc0 ! coefficient of number for Auto-conversion
    real(RP) :: coef_nuc1 !                mass
    real(RP) :: coef_aut0 !                number
    real(RP) :: coef_aut1 !                mass
    real(RP) :: lwc       ! lc+lr
    real(RP) :: tau       ! conversion ratio: qr/(qc+qr) ranges [0:1]
    real(RP) :: rho_fac   ! factor of air density
    real(RP) :: psi_aut   ! Universal function of Auto-conversion
    real(RP) :: psi_acc   ! Universal function of Accretion
    real(RP) :: psi_brk   ! Universal function of Breakup
    real(RP) :: ddr       ! diameter difference from equilibrium
    !
    integer :: i, j, k
    !
    coef_nuc0 = (nu(I_mp_QC)+2.0_RP)/(nu(I_mp_QC)+1.0_RP)
    coef_nuc1 = (nu(I_mp_QC)+2.0_RP)*(nu(I_mp_QC)+4.0_RP)/(nu(I_mp_QC)+1.0_RP)/(nu(I_mp_QC)+1.0_RP)
    coef_aut0 =  -kcc*coef_nuc0
    coef_aut1 =  -kcc/x_sep/20._RP*coef_nuc1
    !
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       lwc = rhoq(I_QR,k,i,j)+rhoq(I_QC,k,i,j)
       if( lwc > xc_min )then
          tau  = max(tau_min, rhoq(I_QR,k,i,j)/lwc)
       else
          tau  = tau_min
       end if
       rho_fac = sqrt(rho_0/max(rho(k,i,j),rho_min))
       !
       ! Auto-conversion ( cloud-cloud => rain )
       psi_aut       = 400._RP*(tau**0.7_RP)*(1.0_RP - (tau**0.7_RP))**3   ! (6) SB06
       PQ(I_NCaut,k,i,j)  = coef_aut0*rhoq(I_QC,k,i,j)*rhoq(I_QC,k,i,j)*rho_fac*rho_fac    ! (9) SB06 sc+aut
       ! lc = lwc*(1-tau), lr = lwc*tau
       PQ(I_LCaut,k,i,j)  = coef_aut1*lwc*lwc*xq(I_mp_QC,k,i,j)*xq(I_mp_QC,k,i,j) &          ! (4) SB06
            *((1.0_RP-tau)*(1.0_RP-tau) + psi_aut)*rho_fac*rho_fac        !
       PQ(I_NRaut,k,i,j)  = -rx_sep*PQ(I_LCaut,k,i,j)                           ! (A7) SB01
       !
       ! Accretion ( cloud-rain => rain )
       psi_acc       =(tau/(tau+thr_acc))**4                          ! (8) SB06
       PQ(I_LCacc,k,i,j)  = -kcr*rhoq(I_QC,k,i,j)*rhoq(I_QR,k,i,j)*rho_fac*psi_acc         ! (7) SB06
       PQ(I_NCacc,k,i,j)  = -kcr*rhoq(I_NC,k,i,j)*rhoq(I_QR,k,i,j)*rho_fac*psi_acc         ! (A6) SB01
       !
       ! Self-collection ( rain-rain => rain )
       PQ(I_NRslc,k,i,j)  = -krr*rhoq(I_NR,k,i,j)*rhoq(I_QR,k,i,j)*rho_fac                 ! (A.8) SB01
       !
       ! Collisional breakup of rain
       ddr           = min(1.E-3_RP, dq_xave(I_mp_QR,k,i,j) - dr_eq )
       if      (dq_xave(I_mp_QR,k,i,j) < dr_min )then                          ! negligible
          psi_brk      = -1.0_RP
          PQ(I_NRbrk,k,i,j) = 0.0_RP
       else if (dq_xave(I_mp_QR,k,i,j) <= dr_eq  )then
          psi_brk      = kbr*ddr + 1.0_RP                               ! (14) SB06 (+1 is necessary)
          PQ(I_NRbrk,k,i,j) = - (psi_brk + 1.0_RP)*PQ(I_NRslc,k,i,j)              ! (13) SB06
       else
          psi_brk      = 2.0_RP*exp(kapbr*ddr) - 1.0_RP                   ! (15) SB06
          PQ(I_NRbrk,k,i,j) = - (psi_brk + 1.0_RP)*PQ(I_NRslc,k,i,j)              ! (13) SB06
       end if
       !
    end do
    end do
    end do
    !
    return
  end subroutine aut_acc_slc_brk_kij
  ! Vapor Deposition, Ice Melting
  subroutine dep_vapor_melt_ice_kij( &
       PQ,                   & ! out
       rho, tem, pre, qd,  & ! in
       rhoq,                   & ! in
       esw, esi,             & ! in
       xq, vt_xave, dq_xave  ) ! in
    use scale_const, only: &
       eps => CONST_EPS
    use scale_tracer, only: &
       QA
    implicit none

    ! Diffusion growth or Evaporation, Sublimation
    real(RP), intent(inout) :: PQ(PQ_MAX,KA,IA,JA)  ! mass change   for cloud, [Add]  09/08/18 T.Mitsui

    real(RP), intent(in)  :: rho(KA,IA,JA)     ! air density
    real(RP), intent(in)  :: tem(KA,IA,JA)     ! air temperature
    real(RP), intent(in)  :: pre(KA,IA,JA)     ! air pressure
    real(RP), intent(in)  :: qd(KA,IA,JA)      ! mixing ratio of dry air
    real(RP), intent(in)  :: esw(KA,IA,JA)     ! saturation vapor pressure(liquid water)
    real(RP), intent(in)  :: esi(KA,IA,JA)     ! saturation vapor pressure(solid water)
    real(RP), intent(in)  :: rhoq(I_QV:I_NG,KA,IA,JA)
    real(RP), intent(in)  :: xq(HYDRO_MAX,KA,IA,JA)    ! mean mass
    ! Notice following values differ from mean terminal velocity or diameter.
    ! mean(vt(x)) /= vt(mean(x)) and mean(D(x)) /= D(mean(x))
    ! Following ones are vt(mean(x)) and D(mean(x)).
    real(RP), intent(in)  :: vt_xave(HYDRO_MAX,2,KA,IA,JA) ! terminal velocity of mean cloud 09/08/18 [Add], T.Mitsui
    !
    real(RP), intent(in)  :: dq_xave(HYDRO_MAX,KA,IA,JA) ! diameter
    !
    real(RP) :: rho_lim            ! limited density              09/08/18 T.Mitsui
    real(RP) :: temc_lim           ! limited temperature[celsius] 09/08/18 T.Mitsui
    real(RP) :: pre_lim            ! limited density              09/08/18 T.Mitsui
    real(RP) :: temc               ! temperature[celsius]
!    real(RP) :: pv                 ! vapor pressure
    real(RP) :: qv                 ! mixing ratio of water vapor [Add] 09/08/18
!    real(RP) :: ssw                ! super saturation ratio(liquid water)
!    real(RP) :: ssi                ! super saturation ratio(ice water)
    real(RP) :: nua, r_nua         ! kinematic viscosity of air
    real(RP) :: mua                ! viscosity of air
    real(RP) :: Kalfa              ! thermal conductance
    real(RP) :: Dw                 ! diffusivity of water vapor
    real(RP) :: Dt                 ! diffusivity of heat
    real(RP) :: Gw, Gi             ! diffusion factor by balance between heat and vapor
    real(RP) :: Gwr, Gii, Gis, Gig ! for rain, ice, snow and graupel.
    real(RP) :: Gm                 ! melting factor by balance between heat and vapor
    real(RP) :: Nsc_r3             !
    ! [Mod] 11/08/30 T.Mitsui, considering large and small branches
!    real(RP) :: Nrecs_r2            ! 09/08/18 [Add] T.Mitsui
    real(RP) :: Nrers_r2, Nreis_r2  !
    real(RP) :: Nress_r2, Nregs_r2  !
!    real(RP) :: Nrecl_r2            ! 09/08/18 [Add] T.Mitsui
    real(RP) :: Nrerl_r2, Nreil_r2  !
    real(RP) :: Nresl_r2, Nregl_r2  !
    real(RP) :: NscNrer_s, NscNrer_l
    real(RP) :: NscNrei_s, NscNrei_l
    real(RP) :: NscNres_s, NscNres_l
    real(RP) :: NscNreg_s, NscNreg_l
    real(RP) :: ventLR_s, ventLR_l
    real(RP) :: ventNI_s, ventNI_l, ventLI_s, ventLI_l
    real(RP) :: ventNS_s, ventNS_l, ventLS_s, ventLS_l
    real(RP) :: ventNG_s, ventNG_l, ventLG_s, ventLG_l
    !
    real(RP) :: wtr, wti, wts, wtg
    real(RP), parameter :: r_14=1.0_RP/1.4_RP
    real(RP), parameter :: r_15=1.0_RP/1.5_RP
    !
    real(RP) ::         ventLR
    real(RP) :: ventNI, ventLI
    real(RP) :: ventNS, ventLS
    real(RP) :: ventNG, ventLG
    !
    real(RP), parameter :: Re_max=1.E+3_RP
    real(RP), parameter :: Re_min=1.E-4_RP

    real(RP) :: sw
    !
    integer :: i, j, k
    !
    ! Notice,T.Mitsui
    ! Vapor deposition and melting would not be solved iteratively to reach equilibrium.
    ! Because following phenomena are not adjustment but transition.
    ! Just time-scales differ among them.
    ! If we would treat more appropreately, there would be time-splitting method to solve each ones.

    PROFILE_START("sn14_dep_vapor")
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       temc    = tem(k,i,j) - T00   ! degC
       temc_lim= max(temc, -40._RP )       ! [Add] 09/08/18 T.Mitsui, Pruppacher and Klett(1997),(13-3)
       rho_lim = max(rho(k,i,j),rho_min)   ! [Add] 09/08/18 T.Mitsui
       qv      = rhoq(I_QV,k,i,j)/rho_lim
       pre_lim = rho_lim*(qd(k,i,j)*Rdry + qv*Rvap)*(temc_lim+T00) ![Add] 09/08/18 T.Mitsui
       !--------------------------------------------------------------------
       ! Diffusion growth part is described in detail
       ! by Pruppacher and Klett (1997) Sec. 13.2(liquid) and 13.3(solid)
       !
       ! G:factor of thermal diffusion(1st.term) and vapor diffusion(2nd. term)
       ! SB06(23),(38), Lin et al(31),(52) or others
       ! Dw is introduced by Pruppacher and Klett(1997),(13-3)
       Dw      = 0.211E-4_RP* (((temc_lim+T00)/T00)**1.94_RP) *(P00/pre_lim)
       Kalfa      = Ka0  + temc_lim*dKa_dT
       mua     = mua0 + temc_lim*dmua_dT
       nua     = mua/rho_lim
       r_nua   = 1.0_RP/nua
       Gw      = (LHV0/Kalfa/tem(k,i,j))*(LHV0/Rvap/tem(k,i,j)-1.0_RP)+(Rvap*tem(k,i,j)/Dw/esw(k,i,j))
       Gi      = (LHS0/Kalfa/tem(k,i,j))*(LHS0/Rvap/tem(k,i,j)-1.0_RP)+(Rvap*tem(k,i,j)/Dw/esi(k,i,j))
       ! capacities account for their surface geometries
       Gwr     = 4.0_RP*PI/cap(I_mp_QR)/Gw
       Gii     = 4.0_RP*PI/cap(I_mp_QI)/Gi
       Gis     = 4.0_RP*PI/cap(I_mp_QS)/Gi
       Gig     = 4.0_RP*PI/cap(I_mp_QG)/Gi
       ! vent: ventilation effect( asymmetry vapor field around particles due to aerodynamic )
       ! SB06 (30),(31) and each coefficient is by (88),(89)
       Nsc_r3  = (nua/Dw)**(0.33333333_RP)                    ! (Schmidt number )^(1/3)
       !
!       Nrecs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QC,1,k,i,j)*dq_xave(I_mp_QC,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud
       Nrers_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QR,1,k,i,j)*dq_xave(I_mp_QR,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) rain
       Nreis_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QI,1,k,i,j)*dq_xave(I_mp_QI,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
       Nress_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QS,1,k,i,j)*dq_xave(I_mp_QS,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) snow
       Nregs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QG,1,k,i,j)*dq_xave(I_mp_QG,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) graupel
       !
!       Nrecl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QC,2,k,i,j)*dq_xave(I_mp_QC,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud
       Nrerl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QR,2,k,i,j)*dq_xave(I_mp_QR,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) rain
       Nreil_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QI,2,k,i,j)*dq_xave(I_mp_QI,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
       Nresl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QS,2,k,i,j)*dq_xave(I_mp_QS,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) snow
       Nregl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(I_mp_QG,2,k,i,j)*dq_xave(I_mp_QG,k,i,j)*r_nua))) ! (Reynolds number)^(1/2) graupel
       NscNrer_s=Nsc_r3*Nrers_r2 ! small rain
       NscNrer_l=Nsc_r3*Nrerl_r2 ! large rain
       !
       NscNrei_s=Nsc_r3*Nreis_r2 ! small ice
       NscNrei_l=Nsc_r3*Nreil_r2 ! large ice
       !
       NscNres_s=Nsc_r3*Nress_r2 ! small snow
       NscNres_l=Nsc_r3*Nresl_r2 ! large snow
       !
       NscNreg_s=Nsc_r3*Nregs_r2 ! small snow
       NscNreg_l=Nsc_r3*Nregl_r2 ! large snow
       !
       ventLR_s = ah_vent1(I_mp_QR,1) + bh_vent1(I_mp_QR,1)*NscNrer_s
       ventLR_l = ah_vent1(I_mp_QR,2) + bh_vent1(I_mp_QR,2)*NscNrer_l
       !
       ventNI_s = ah_vent0(I_mp_QI,1) + bh_vent0(I_mp_QI,1)*NscNrei_s
       ventNI_l = ah_vent0(I_mp_QI,2) + bh_vent0(I_mp_QI,2)*NscNrei_l
       ventLI_s = ah_vent1(I_mp_QI,1) + bh_vent1(I_mp_QI,1)*NscNrei_s
       ventLI_l = ah_vent1(I_mp_QI,2) + bh_vent1(I_mp_QI,2)*NscNrei_l
       !
       ventNS_s = ah_vent0(I_mp_QS,1) + bh_vent0(I_mp_QS,1)*NscNres_s
       ventNS_l = ah_vent0(I_mp_QS,2) + bh_vent0(I_mp_QS,2)*NscNres_l
       ventLS_s = ah_vent1(I_mp_QS,1) + bh_vent1(I_mp_QS,1)*NscNres_s
       ventLS_l = ah_vent1(I_mp_QS,2) + bh_vent1(I_mp_QS,2)*NscNres_l
       !
       ventNG_s = ah_vent0(I_mp_QG,1) + bh_vent0(I_mp_QG,1)*NscNreg_s
       ventNG_l = ah_vent0(I_mp_QG,2) + bh_vent0(I_mp_QG,2)*NscNreg_l
       ventLG_s = ah_vent1(I_mp_QG,1) + bh_vent1(I_mp_QG,1)*NscNreg_s
       ventLG_l = ah_vent1(I_mp_QG,2) + bh_vent1(I_mp_QG,2)*NscNreg_l
       !
       ! branch is 1.4 for rain, snow, graupel; is 1.0 for ice (PK97, 13-60,-61,-88,-89).
       !
       wtr     = ( min(max( NscNrer_s*r_14, 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.4*0.5 and 1.4*2
       wti     = ( min(max( NscNrei_s     , 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.0*0.5 and 1.0*2
       wts     = ( min(max( NscNres_s*r_14, 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.4*0.5 and 1.4*2
       wtg     = ( min(max( NscNreg_s*r_14, 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.4*0.5 and 1.4*2
       ! interpolation between two branches
       ventNI  = (1.0_RP-wti)*ventNI_s + wti*ventNI_l
       ventNS  = (1.0_RP-wts)*ventNS_s + wts*ventNS_l
       ventNG  = (1.0_RP-wtg)*ventNG_s + wtg*ventNG_l
       !
       ventLR  = (1.0_RP-wtr)*ventLR_s + wtr*ventLR_l
       ventLI  = (1.0_RP-wti)*ventLI_s + wti*ventLI_l
       ventLS  = (1.0_RP-wts)*ventLS_s + wts*ventLS_l
       ventLG  = (1.0_RP-wtg)*ventLG_s + wtg*ventLG_l
       !
       ! SB06(29)
       ! [Mod] 08/05/08 T.Mitsui, recover PNXdep, and rain is only evaporation.
       ! Ni, Ns, Ng should decrease in nature so we add this term.
       ! And vapor deposition never occur unless number exist.
       ! [Add comment]  09/08/18
       ! recover condensation/evaporation of rain,
       ! and ventilation effects are not taken into account for cloud.
       !
!!$***************************************************************************
!!$        NOTICE:
!!$         09/08/18 [Mod] Hereafter PLxdep means inverse of timescale.
!!$***************************************************************************
!!$     PQ(I_LCdep,k,i,j) = Gwr*ssw*rhoq(I_NC,k,i,j)*dq_xave(I_mp_QC,k,i,j)*coef_deplc
!!$     PQ(I_LRdep,k,i,j) = Gwr*ssw*rhoq(I_NR,k,i,j)*dq_xave(I_mp_QR,k,i,j)*ventLR
!!$     PQ(I_LIdep,k,i,j) = Gii*ssi*rhoq(I_NI,k,i,j)*dq_xave(I_mp_QI,k,i,j)*ventLI
!!$     PQ(I_LSdep,k,i,j) = Gis*ssi*rhoq(I_NS,k,i,j)*dq_xave(I_mp_QS,k,i,j)*ventLS
!!$     PQ(I_LGdep,k,i,j) = Gig*ssi*rhoq(I_NG,k,i,j)*dq_xave(I_mp_QG,k,i,j)*ventLG
       PQ(I_LCdep,k,i,j) = Gwr*rhoq(I_NC,k,i,j)*dq_xave(I_mp_QC,k,i,j)*coef_deplc
       PQ(I_LRdep,k,i,j) = Gwr*rhoq(I_NR,k,i,j)*dq_xave(I_mp_QR,k,i,j)*ventLR
       PQ(I_LIdep,k,i,j) = Gii*rhoq(I_NI,k,i,j)*dq_xave(I_mp_QI,k,i,j)*ventLI
       PQ(I_LSdep,k,i,j) = Gis*rhoq(I_NS,k,i,j)*dq_xave(I_mp_QS,k,i,j)*ventLS
       PQ(I_LGdep,k,i,j) = Gig*rhoq(I_NG,k,i,j)*dq_xave(I_mp_QG,k,i,j)*ventLG
       PQ(I_NRdep,k,i,j) = PQ(I_LRdep,k,i,j)/xq(I_mp_QR,k,i,j)
       PQ(I_NIdep,k,i,j) = 0.0_RP
       PQ(I_NSdep,k,i,j) = PQ(I_LSdep,k,i,j)/xq(I_mp_QS,k,i,j)
       PQ(I_NGdep,k,i,j) = PQ(I_LGdep,k,i,j)/xq(I_mp_QG,k,i,j)
       !
       !------------------------------------------------------------------------
       ! Melting part is described by Pruppacher and Klett (1997) Sec.16.3.1
       ! Here we omit "Shedding" of snow-flakes and ice-particles.
       ! "Shedding" may be applicative if you refer
       ! eq.(38) in Cotton etal.(1986) Jour. Clim. Appl. Meteor. p.1658-1680.
       ! SB06(73)
       Dt      = Kalfa/(CPvap*rho_0)
       ! Gm: factor caused by balance between
       !     "water evaporation cooling(1st.)" and "fusion heating(2nd.)"
       ! SB06(76)
       ! [fix] 08/05/08 T.Mitsui  LHF00 => EMELT  and  esw => PSAT0
       ! LHS0 is more suitable than LHS because melting occurs around 273.15 K.
       Gm      = 2.0_RP*PI/EMELT&
               * ( (Kalfa*Dt/Dw)*(temc) + (Dw*LHS0/Rvap)*(esi(k,i,j)/tem(k,i,j)-PSAT0/T00) )
       ! SB06(76)
       ! Notice! melting only occurs where T > 273.15 K else doesn't.
       ! [fix] 08/05/08 T.Mitsui, Gm could be both positive and negative value.
       !       See Pruppacher and Klett(1997) eq.(16-79) or Rasmussen and Pruppacher(1982)
       sw = ( sign(0.5_RP,temc) + 0.5_RP ) * ( sign(0.5_RP,Gm-eps) + 0.5_RP ) ! sw = 1 if( (temc>=0.0_RP) .AND. (Gm>0.0_RP) ), otherwise sw = 0
       !  if Gm==0 then rh and tem is critical value for melting process.
       ! 08/05/16 [Mod] T.Mitsui, change term of PLimlt. N_i => L_i/ (limited x_i)
       ! because melting never occur when N_i=0.
       PQ(I_LImlt,k,i,j) = - Gm * rhoq(I_QI,k,i,j)*dq_xave(I_mp_QI,k,i,j)*ventLI/xq(I_mp_QI,k,i,j) * sw
       ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
       PQ(I_NImlt,k,i,j) = - Gm * rhoq(I_NI,k,i,j)*dq_xave(I_mp_QI,k,i,j)*ventNI/xq(I_mp_QI,k,i,j) * sw ! 09/08/18 [Mod] recover, T.Mitsui
       PQ(I_LSmlt,k,i,j) = - Gm * rhoq(I_QS,k,i,j)*dq_xave(I_mp_QS,k,i,j)*ventLS/xq(I_mp_QS,k,i,j) * sw
       ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
       PQ(I_NSmlt,k,i,j) = - Gm * rhoq(I_NS,k,i,j)*dq_xave(I_mp_QS,k,i,j)*ventNS/xq(I_mp_QS,k,i,j) * sw ! 09/08/18 [Mod] recover, T.Mitsui
       PQ(I_LGmlt,k,i,j) = - Gm * rhoq(I_QG,k,i,j)*dq_xave(I_mp_QG,k,i,j)*ventLG/xq(I_mp_QG,k,i,j) * sw
       ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
       PQ(I_NGmlt,k,i,j) = - Gm * rhoq(I_NG,k,i,j)*dq_xave(I_mp_QG,k,i,j)*ventNG/xq(I_mp_QG,k,i,j) * sw ! 09/08/18 [Mod] recover, T.Mitsui

    end do
    end do
    end do
    PROFILE_STOP("sn14_dep_vapor")
    !
    return
  end subroutine dep_vapor_melt_ice_kij
  !-----------------------------------------------------------------------------
  subroutine freezing_water_kij( &
       dt,           &
       PQ,           &
       rhoq, xq, tem )
    use scale_tracer, only: &
       QA
    implicit none
    !
    ! In this subroutine,
    ! We assumed surface temperature of droplets are same as environment.

    real(RP), intent(in) :: dt
    real(RP), intent(inout):: PQ(PQ_MAX,KA,IA,JA)
    !
    real(RP), intent(in) :: tem(KA,IA,JA)
    !
    real(RP), intent(in) :: rhoq(I_QV:I_NG,KA,IA,JA)
    real(RP), intent(in) :: xq(HYDRO_MAX,KA,IA,JA)
    !
    real(RP), parameter :: temc_min = -65.0_RP
    real(RP), parameter :: a_het = 0.2_RP  ! SB06 (44)
    real(RP), parameter :: b_het = 0.65_RP ! SB06 (44)
    !
    real(RP) :: coef_m2_c
    real(RP) :: coef_m2_r
    ! temperature [celsius]
    real(RP) :: temc, temc2, temc3, temc4
    ! temperature function of homegenous/heterogenous freezing
    real(RP) :: Jhom, Jhet
    real(RP) :: rdt
    real(RP) :: tmp
    !
    integer :: i,j,k
    !
    rdt = 1.0_RP/dt
    !
    coef_m2_c =   coef_m2(I_mp_QC)
    coef_m2_r =   coef_m2(I_mp_QR)
    !
    PROFILE_START("sn14_freezing")
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          temc = max( tem(k,i,j) - T00, temc_min )
          ! These cause from aerosol-droplet interaction.
          ! Bigg(1953) formula, Khain etal.(2000) eq.(4.5), Pruppacher and Klett(1997) eq.(9-48)
          Jhet =  a_het*exp( -b_het*temc - 1.0_RP )
          ! These cause in nature.
          ! Cotton and Field 2002, QJRMS. (12)
          if( temc < -65.0_RP )then
             jhom = 10.0_RP**(24.37236_RP)*1.E+3_RP
             Jhet =  a_het*exp( 65.0_RP*b_het - 1.0_RP ) ! 09/04/14 [Add], fixer T.Mitsui
          else if( temc < -30.0_RP ) then
             temc2 = temc*temc
             temc3 = temc*temc2
             temc4 = temc2*temc2
             jhom = 10.0_RP**(&
                  - 243.40_RP - 14.75_RP*temc - 0.307_RP*temc2 &
                  - 0.00287_RP*temc3 - 0.0000102_RP*temc4 ) *1.E+3_RP
          else if( temc < 0.0_RP) then
             jhom = 10._RP**(-7.63_RP-2.996_RP*(temc+30.0_RP))*1.E+3_RP
          else
             Jhom = 0.0_RP
             Jhet = 0.0_RP
          end if
          ! Note, xc should be limited in range[xc_min:xc_max].
          ! and PNChom need to be calculated by NC
          ! because reduction rate of Nc need to be bound by NC.
          ! For the same reason PLChom also bound by LC and xc.
          ! Basically L and N should be independent
          ! but average particle mass x should be in suitable range.
          ! Homogenous Freezing
          PQ(I_LChom,k,i,j) = 0.0_RP
          PQ(I_NChom,k,i,j) = 0.0_RP
          ! Heterogenous Freezing
#if defined(__PGI) || defined(__ES2)
          tmp = min( xq(I_mp_QC,k,i,j)*(Jhet+Jhom)*dt, 1.E+3_RP) ! apply exp limiter
          PQ(I_LChet,k,i,j) = -rdt*rhoq(I_QC,k,i,j)*( 1.0_RP - exp( -coef_m2_c*tmp ) )
          PQ(I_NChet,k,i,j) = -rdt*rhoq(I_NC,k,i,j)*( 1.0_RP - exp( -          tmp ) )

          tmp = min( xq(I_mp_QR,k,i,j)*(Jhet+Jhom)*dt, 1.E+3_RP) ! apply exp limiter
          PQ(I_LRhet,k,i,j) = -rdt*rhoq(I_QR,k,i,j)*( 1.0_RP - exp( -coef_m2_r*tmp ) )
          PQ(I_NRhet,k,i,j) = -rdt*rhoq(I_NR,k,i,j)*( 1.0_RP - exp( -          tmp ) )
#else
          PQ(I_LChet,k,i,j) = -rdt*rhoq(I_QC,k,i,j)*( 1.0_RP - exp( -coef_m2_c*xq(I_mp_QC,k,i,j)*(Jhet+Jhom)*dt ) )
          PQ(I_NChet,k,i,j) = -rdt*rhoq(I_NC,k,i,j)*( 1.0_RP - exp( -          xq(I_mp_QC,k,i,j)*(Jhet+Jhom)*dt ) )
          PQ(I_LRhet,k,i,j) = -rdt*rhoq(I_QR,k,i,j)*( 1.0_RP - exp( -coef_m2_r*xq(I_mp_QR,k,i,j)*(Jhet+Jhom)*dt ) )
          PQ(I_NRhet,k,i,j) = -rdt*rhoq(I_NR,k,i,j)*( 1.0_RP - exp( -          xq(I_mp_QR,k,i,j)*(Jhet+Jhom)*dt ) )
#endif
       end do
    end do
    end do
    PROFILE_STOP("sn14_freezing")
    !
    return
  end subroutine freezing_water_kij
  !-----------------------------------------------------------------------------
  subroutine MP_terminal_velocity( &
      velw, &
      rhoq, &
      DENS, &
      temp, &
      pres  )
    use scale_const, only: &
       CONST_UNDEF
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: velw(KA,IA,JA,QA_MP-1) ! terminal velocity of cloud mass
    real(RP), intent(in)  :: rhoq(I_QV:I_NG,KA,IA,JA) ! rho * q
    real(RP), intent(in)  :: DENS(KA,IA,JA)    ! rho
    real(RP), intent(in)  :: temp(KA,IA,JA)    ! temperature
    real(RP), intent(in)  :: pres(KA,IA,JA)    ! pressure

    real(RP) :: xq       ! average mass of 1 particle( mass/number )

    real(RP) :: rhofac   ! density factor for terminal velocity( air friction )
    real(RP) :: rhofac_q(HYDRO_MAX)

    real(RP) :: rlambdar ! work for diagnosis of Rain DSD ( Seifert, 2008 )
    real(RP) :: mud_r
    real(RP) :: dq, dql  ! weigthed diameter.   Improved Rogers etal. (1993) formula by T.Mitsui


    real(RP) :: weight ! weighting coefficient for 2-branches is determined by ratio between 0.745mm and weighted diameter.  SB06 Table.1
    real(RP) :: velq_s ! terminal velocity for small branch of Rogers formula
    real(RP) :: velq_l ! terminal velocity for large branch of Rogers formula

    real(RP) :: tmp
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    PROFILE_START("sn14_terminal_vel")

    mud_r = 3.0_RP * nu(I_mp_QR) + 2.0_RP


    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rhofac = rho_0 / max( DENS(k,i,j), rho_min )

       ! QC
       rhofac_q(I_mp_QC) = rhofac ** gamma_v(I_mp_QC)
       xq = max( xqmin(I_mp_QC), min( xqmax(I_mp_QC), rhoq(I_QC,k,i,j) / ( rhoq(I_NC,k,i,j) + nqmin(I_mp_QC) ) ) )

       velw(k,i,j,I_mp_QC) = -rhofac_q(I_mp_QC) * coef_vt1(I_mp_QC,1) * xq**beta_v(I_mp_QC,1)
       ! NC
       velw(k,i,j,I_mp_NC) = -rhofac_q(I_mp_QC) * coef_vt0(I_mp_QC,1) * xq**beta_vn(I_mp_QC,1)

       ! QR
       rhofac_q(I_mp_QR) = rhofac ** gamma_v(I_mp_QR)
       xq = max( xqmin(I_mp_QR), min( xqmax(I_mp_QR), rhoq(I_QR,k,i,j) / ( rhoq(I_NR,k,i,j) + nqmin(I_mp_QR) ) ) )

       rlambdar = a_m(I_mp_QR) * xq**b_m(I_mp_QR) &
                   * ( (mud_r+3.0_RP) * (mud_r+2.0_RP) * (mud_r+1.0_RP) )**(-0.333333333_RP)
       dq = ( 4.0_RP + mud_r ) * rlambdar ! D^(3)+mu weighted mean diameter
       dql = dq
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + tanh( PI * log( dq/d_vtr_branch ) ) ) ) )

       velq_s = coef_vtr_ar2 * dq &
            * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar )**(-5.0_RP-mud_r) )
       velq_l = coef_vtr_ar1 - coef_vtr_br1 &
            * ( 1.0_RP + coef_vtr_cr1*rlambdar )**(-4.0_RP-mud_r)
       velw(k,i,j,I_mp_QR) = -rhofac_q(I_mp_QR) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )
       ! NR
       dq = ( 1.0_RP + mud_r ) * rlambdar
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + tanh( PI * log( dq/d_vtr_branch ) ) ) ) )

       velq_s = coef_vtr_ar2 * dql &
            * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar )**(-2.0_RP-mud_r) )
       velq_l = coef_vtr_ar1 - coef_vtr_br1 &
            * ( 1.0_RP + coef_vtr_cr1*rlambdar )**(-1.0_RP-mud_r)
       velw(k,i,j,I_mp_NR) = -rhofac_q(I_mp_QR) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )

       ! QI
       rhofac_q(I_mp_QI) = ( pres(k,i,j)/pre0_vt )**a_pre0_vt * ( temp(k,i,j)/tem0_vt )**a_tem0_vt
       xq = max( xqmin(I_mp_QI), min( xqmax(I_mp_QI), rhoq(I_QI,k,i,j) / ( rhoq(I_NI,k,i,j) + nqmin(I_mp_QI) ) ) )

       tmp = a_m(I_mp_QI) * xq**b_m(I_mp_QI)
       dq = coef_dave_L(I_mp_QI) * tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_li ) ) ) )

       velq_s = coef_vt1(I_mp_QI,1) * xq**beta_v (I_mp_QI,1)
       velq_l = coef_vt1(I_mp_QI,2) * xq**beta_v (I_mp_QI,2)
       velw(k,i,j,I_mp_QI) = -rhofac_q(I_mp_QI) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )
       ! NI
       dq = coef_dave_N(I_mp_QI) * tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ni ) ) ) )

       velq_s = coef_vt0(I_mp_QI,1) * xq**beta_vn(I_mp_QI,1)
       velq_l = coef_vt0(I_mp_QI,2) * xq**beta_vn(I_mp_QI,2)
       velw(k,i,j,I_mp_NI) = -rhofac_q(I_mp_QI) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )

       ! QS
       rhofac_q(I_mp_QS) = rhofac_q(I_mp_QI)
       xq = max( xqmin(I_mp_QS), min( xqmax(I_mp_QS), rhoq(I_QS,k,i,j) / ( rhoq(I_NS,k,i,j) + nqmin(I_mp_QS) ) ) )

       tmp = a_m(I_mp_QS) * xq**b_m(I_mp_QS)
       dq = coef_dave_L(I_mp_QS) * tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ls ) ) ) )

       velq_s = coef_vt1(I_mp_QS,1) * xq**beta_v (I_mp_QS,1)
       velq_l = coef_vt1(I_mp_QS,2) * xq**beta_v (I_mp_QS,2)
       velw(k,i,j,I_mp_QS) = -rhofac_q(I_mp_QS) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )
       ! NS
       dq = coef_dave_N(I_mp_QS) * tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ns ) ) ) )

       velq_s = coef_vt0(I_mp_QS,1) * xq**beta_vn(I_mp_QS,1)
       velq_l = coef_vt0(I_mp_QS,2) * xq**beta_vn(I_mp_QS,2)
       velw(k,i,j,I_mp_NS) = -rhofac_q(I_mp_QS) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )

       ! QG
       rhofac_q(I_mp_QG) = rhofac_q(I_mp_QI)
       xq = max( xqmin(I_mp_QG), min( xqmax(I_mp_QG), rhoq(I_QG,k,i,j) / ( rhoq(I_NG,k,i,j) + nqmin(I_mp_QG) ) ) )

       tmp = a_m(I_mp_QG) * xq**b_m(I_mp_QG)
       dq = coef_dave_L(I_mp_QG) * tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_lg ) ) ) )

       velq_s = coef_vt1(I_mp_QG,1) * xq**beta_v (I_mp_QG,1)
       velq_l = coef_vt1(I_mp_QG,2) * xq**beta_v (I_mp_QG,2)
       velw(k,i,j,I_mp_QG) = -rhofac_q(I_mp_QG) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )
       ! NG
       dq = coef_dave_N(I_mp_QG) * tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ng ) ) ) )

       velq_s = coef_vt0(I_mp_QG,1) * xq**beta_vn(I_mp_QG,1)
       velq_l = coef_vt0(I_mp_QG,2) * xq**beta_vn(I_mp_QG,2)
       velw(k,i,j,I_mp_NG) = -rhofac_q(I_mp_QG) &
            * ( velq_l * (          weight ) &
              + velq_s * ( 1.0_RP - weight ) )
    enddo
    enddo
    enddo

    do iq = 1, QA_MP-1
    do j = JS, JE
    do i = IS, IE
       velw(1:KS-2,i,j,iq) = CONST_UNDEF
       velw(KS-1,i,j,iq) = velw(KS,i,j,iq)
       velw(KE+1:KA,i,j,iq) = CONST_UNDEF
    enddo
    enddo
    enddo

    PROFILE_STOP("sn14_terminal_vel")

    return
  end subroutine MP_terminal_velocity
  !----------------------------------------------------------------
  subroutine update_by_phase_change_kij(   &
       ntdiv, ntmax,        & ! in [Add] 10/08/03
       dt,                  & ! in
       gsgam2,              & ! in
       z,                   & ! in
       dz,                  & ! in
       velz,                & ! in
       dTdt_rad,            & ! in
       rho,                 & ! in
       rhoe,                & ! inout
       rhoq, q,             & ! inout
       tem, pre,            & ! inout
       cva,                 & ! out
       esw, esi, rhoq2,     & ! in
       PQ,                  & ! in
       qc_evaporate,        & ! in
       sl_PLCdep,           &
       sl_PLRdep, sl_PNRdep )
    use scale_tracer, only: &
       QA, &
       TRACER_R,  &
       TRACER_CV, &
       TRACER_CP, &
       TRACER_MASS
    use scale_atmos_saturation, only: &
       moist_pres2qsat_liq  => ATMOS_SATURATION_pres2qsat_liq,  &
       moist_pres2qsat_ice  => ATMOS_SATURATION_pres2qsat_ice,  &
       moist_dqsw_dtem_rho  => ATMOS_SATURATION_dqsw_dtem_rho,  &
       moist_dqsi_dtem_rho  => ATMOS_SATURATION_dqsi_dtem_rho,  &
       moist_dqsw_dtem_dpre => ATMOS_SATURATION_dqsw_dtem_dpre, &
       moist_dqsi_dtem_dpre => ATMOS_SATURATION_dqsi_dtem_dpre
    implicit none

    integer, intent(in)    :: ntdiv               ! [Add] 10/08/03
    integer, intent(in)    :: ntmax               ! [Add] 10/08/03
    !
    real(RP), intent(in)    :: dt                 ! time step[s]
    real(RP), intent(in)    :: gsgam2(KA,IA,JA)   ! metric
    real(RP), intent(in)    :: z(KA)              ! altitude [m]
    real(RP), intent(in)    :: dz(KA)             ! altitude [m]
    real(RP), intent(in)    :: velz(KA,IA,JA)     ! vertical velocity @ half point[m/s]
    real(RP), intent(in)    :: dTdt_rad(KA,IA,JA) ! temperture tendency by radiation[K/s]
    real(RP), intent(in)    :: rho(KA,IA,JA)      ! density[kg/m3]
    real(RP), intent(inout) :: rhoe(KA,IA,JA)     ! internal energy[J/m3]
    real(RP), intent(inout) :: rhoq(I_QV:I_NG,KA,IA,JA)  ! tracers[kg/m3]
    real(RP), intent(inout) :: q(KA,IA,JA,QA)     ! tracers mixing ratio[kg/kg]
    real(RP), intent(inout) :: tem(KA,IA,JA)      ! temperature[K]
    real(RP), intent(inout) :: pre(KA,IA,JA)      ! pressure[Pa]
    real(RP), intent(out)   :: cva(KA,IA,JA)      ! specific heat at constant volume
    real(RP), intent(in)    :: esw(KA,IA,JA)      ! saturated vapor pressure for liquid
    real(RP), intent(in)    :: esi(KA,IA,JA)      !                          for ice
    real(RP), intent(in)    :: rhoq2(I_QV:I_NG,KA,IA,JA)
    !+++ tendency[kg/m3/s]
    real(RP), intent(inout) :: PQ(PQ_MAX,KA,IA,JA)
    real(RP), intent(out)   :: qc_evaporate(KA,IA,JA)
    !+++ Column integrated tendency[kg/m2/s]
    real(RP), intent(inout) :: sl_PLCdep(IA,JA)
    real(RP), intent(inout) :: sl_PLRdep(IA,JA), sl_PNRdep(IA,JA)
    !
    real(RP) :: Rmoist
    !
    real(RP) :: xi                     ! mean mass of ice particles
    real(RP) :: rrho                   ! 1/rho
    real(RP) :: wtem(KA,IA,JA)         ! temperature[K]
    real(RP) :: qdry                   ! mixing ratio of dry air
    !
    real(RP) :: r_cva                  ! specific heat at constant volume
    real(RP) :: cpa                    ! specific heat at constant pressure
    real(RP) :: r_cpa                  ! specific heat at constant pressure
    real(RP) :: qsw(KA,IA,JA)          ! saturated mixing ratio for liquid
    real(RP) :: qsi(KA,IA,JA)          ! saturated mixing ratio for solid
    real(RP) :: dqswdtem_rho(KA,IA,JA) ! (dqsw/dtem)_rho
    real(RP) :: dqsidtem_rho(KA,IA,JA) ! (dqsi/dtem)_rho
    real(RP) :: dqswdtem_pre(KA,IA,JA) ! (dqsw/dtem)_pre
    real(RP) :: dqsidtem_pre(KA,IA,JA) ! (dqsi/dtem)_pre
    real(RP) :: dqswdpre_tem(KA,IA,JA) ! (dqsw/dpre)_tem
    real(RP) :: dqsidpre_tem(KA,IA,JA) ! (dqsi/dpre)_tem
    !
    real(RP) :: w                      ! vetical velocity[m/s]
    real(RP) :: Acnd                   ! Pdynliq + Bergeron-Findeisen
    real(RP) :: Adep                   ! Pdyndep + Bergeron-Findeisen
    real(RP) :: aliqliq, asolliq
    real(RP) :: aliqsol, asolsol
    real(RP) :: Pdynliq                ! production term of ssw by vertical motion
    real(RP) :: Pdynsol                ! production term of ssi by vertical motion
    real(RP) :: Pradliq                ! production term of ssw by radiation
    real(RP) :: Pradsol                ! production term of ssi by radiation
    real(RP) :: taucnd, r_taucnd       ! time scale of ssw change by MP
    real(RP) :: taudep, r_taudep       ! time scale of ssi change by MP
    real(RP) :: taucnd_c(KA,IA,JA), r_taucnd_c  ! by cloud
    real(RP) :: taucnd_r(KA,IA,JA), r_taucnd_r  ! by rain
    real(RP) :: taudep_i(KA,IA,JA), r_taudep_i  ! by ice
    real(RP) :: taudep_s(KA,IA,JA), r_taudep_s  ! by snow
    real(RP) :: taudep_g(KA,IA,JA), r_taudep_g  ! by graupel
    ! alternative tendency through changing ssw and ssi
    real(RP) :: PNCdep ! [Add] 11/08/30 T.Mitsui
    real(RP) :: PLR2NR, PLI2NI, PLS2NS, PLG2NG
    real(RP) :: coef_a_cnd, coef_b_cnd
    real(RP) :: coef_a_dep, coef_b_dep
    !
    real(RP) :: frz_dqc
    real(RP) :: frz_dnc
    real(RP) :: frz_dqr
    real(RP) :: frz_dnr
    real(RP) :: mlt_dqi
    real(RP) :: mlt_dni
    real(RP) :: mlt_dqs
    real(RP) :: mlt_dns
    real(RP) :: mlt_dqg
    real(RP) :: mlt_dng
    real(RP) :: dep_dqi
    real(RP) :: dep_dni
    real(RP) :: dep_dqs
    real(RP) :: dep_dns
    real(RP) :: dep_dqg
    real(RP) :: dep_dng
    real(RP) :: dep_dqr
    real(RP) :: dep_dnr
    real(RP) :: dep_dqc
    real(RP) :: dep_dnc   ! 11/08/30 [Add] T.Mitsui, dep_dnc
    real(RP) :: r_xc_ccn, r_xi_ccn ! 11/08/30 [Add] T.Mitsui
    !
    real(RP) :: drhoqv
    real(RP) :: drhoqc, drhoqr, drhoqi, drhoqs, drhoqg
    real(RP) :: drhonc, drhonr, drhoni, drhons, drhong
    !
    real(RP) :: fac1, fac2, fac3, fac4, fac5, fac6
    real(RP) :: r_rvaptem        ! 1/(Rvap*tem)
    real(RP) :: pv               ! vapor pressure
    real(RP) :: lvsw, lvsi       ! saturated vapor density
    real(RP) :: dlvsw, dlvsi     !
    ! [Add] 11/08/30 T.Mitsui
    real(RP) :: dcnd, ddep       ! total cndensation/deposition
    real(RP) :: uplim_cnd        ! upper limit of condensation
    real(RP) :: lowlim_cnd       ! lower limit of evaporation
    ! [Add] 11/08/30 T.Mitsui
    real(RP) :: uplim_dep        ! upper limit of condensation
    real(RP) :: lowlim_dep       ! lower limit of evaporation
    real(RP) :: ssw, ssi         ! supersaturation ratio
    real(RP) :: r_esw, r_esi     ! 1/esw, 1/esi
    real(RP) :: r_lvsw, r_lvsi   ! 1/(lvsw*ssw), 1/(lvsi*ssi)
    real(RP) :: r_dt             ! 1/dt
    real(RP) :: ssw_o, ssi_o
!    real(RP) :: dt_dyn
!    real(RP) :: dt_mp
    !
!    real(RP)       :: tem_lh(KA,IA,JA)
!    real(RP)       :: dtemdt_lh(KA,IA,JA)
    real(RP), save :: fac_cndc        = 1.0_RP
    logical, save :: opt_fix_taucnd_c=.false.
    logical, save :: flag_first      =.true.
    !
    namelist /nm_mp_sn14_condensation/ &
         opt_fix_taucnd_c, fac_cndc

    real(RP) :: fac_cndc_wrk
    !
    real(RP), parameter :: tau100day   = 1.E+7_RP
    real(RP), parameter :: r_tau100day = 1.E-7_RP
    real(RP), parameter :: eps=1.E-30_RP
    !
    integer :: i,j,k,iqw
    real(RP) :: sw
    !

    ! [Add] 11/08/30 T.Mitsui
    if( flag_first )then
       flag_first = .false.
       rewind(IO_FID_CONF)
       read  (IO_FID_CONF,nml=nm_mp_sn14_condensation, end=100)
100    if( IO_NML ) write(IO_FID_NML,nml=nm_mp_sn14_condensation)
    end if
    !
!    dt_dyn     = dt*ntmax
!    dt_mp      = dt*(ntdiv-1)
    !
    r_dt       = 1.0_RP/dt
    !
    r_xc_ccn=1.0_RP/xc_ccn
!    r_xi_ccn=1.0_RP/xi_ccn
    !
    if( opt_fix_taucnd_c )then
       fac_cndc_wrk = fac_cndc**(1.0_RP-b_m(I_mp_QC))
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PQ(I_LCdep,k,i,j)  = PQ(I_LCdep,k,i,j)*fac_cndc_wrk
       end do
       end do
       end do
       if( IO_L ) write(IO_FID_LOG,*) "taucnd:fac_cndc_wrk=",fac_cndc_wrk
    end if

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! Temperature lower limit is only used for saturation condition.
       ! On the other hand original "tem" is used for calculation of latent heat or energy equation.
       wtem(k,i,j)  = max( tem(k,i,j), tem_min )
    end do
    end do
    end do

    call moist_pres2qsat_liq ( qsw, wtem, pre )
    call moist_pres2qsat_ice ( qsi, wtem, pre )
    call moist_dqsw_dtem_rho ( dqswdtem_rho, wtem, rho )
    call moist_dqsi_dtem_rho ( dqsidtem_rho, wtem, rho )
    call moist_dqsw_dtem_dpre( dqswdtem_pre, dqswdpre_tem, wtem, pre )
    call moist_dqsi_dtem_dpre( dqsidtem_pre, dqsidpre_tem, wtem, pre )

    PROFILE_START("sn14_update")
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if( z(k) <= 25000.0_RP )then
          w = 0.5_RP*(velz(k-1,i,j) + velz(k,i,j))
       else
          w = 0.0_RP
       end if
       if( pre(k,i,j) < esw(k,i,j)+1.E-10_RP )then
          qsw(k,i,j) = 1.0_RP
          dqswdtem_rho(k,i,j) = 0.0_RP
          dqswdtem_pre(k,i,j) = 0.0_RP
          dqswdpre_tem(k,i,j) = 0.0_RP
       end if
       if( pre(k,i,j) < esi(k,i,j)+1.E-10_RP )then
          qsi(k,i,j) = 1.0_RP
          dqsidtem_rho(k,i,j) = 0.0_RP
          dqsidtem_pre(k,i,j) = 0.0_RP
          dqsidpre_tem(k,i,j) = 0.0_RP
       end if

       r_rvaptem        = 1.0_RP/(Rvap*wtem(k,i,j))
       lvsw             = esw(k,i,j)*r_rvaptem        ! rho=p/(Rv*T)
       lvsi             = esi(k,i,j)*r_rvaptem        !
       pv               = rhoq2(I_QV,k,i,j)*Rvap*tem(k,i,j)
       r_esw            = 1.0_RP/esw(k,i,j)
       r_esi            = 1.0_RP/esi(k,i,j)
       ssw              = min( MP_ssw_lim, ( pv*r_esw-1.0_RP ) )
       ssi              = pv*r_esi - 1.0_RP
       r_lvsw           = 1.0_RP/lvsw
       r_lvsi           = 1.0_RP/lvsi
       r_taucnd_c       = PQ(I_LCdep,k,i,j)*r_lvsw
       r_taucnd_r       = PQ(I_LRdep,k,i,j)*r_lvsw
       r_taudep_i       = PQ(I_LIdep,k,i,j)*r_lvsi
       r_taudep_s       = PQ(I_LSdep,k,i,j)*r_lvsi
       r_taudep_g       = PQ(I_LGdep,k,i,j)*r_lvsi
!       taucnd_c(k,i,j)   = 1.0_RP/(r_taucnd_c+r_tau100day)
!       taucnd_r(k,i,j)   = 1.0_RP/(r_taucnd_r+r_tau100day)
!       taudep_i(k,i,j)   = 1.0_RP/(r_taudep_i+r_tau100day)
!       taudep_s(k,i,j)   = 1.0_RP/(r_taudep_s+r_tau100day)
!       taudep_g(k,i,j)   = 1.0_RP/(r_taudep_g+r_tau100day)

       CALC_QDRY( qdry, q, TRACER_MASS, k, i, j, iqw )
       CALC_CV( cva(k,i,j), qdry, q, k, i, j, iqw, CVdry, TRACER_CV )
       CALC_CP( cpa, qdry, q, k, i, j, iqw, CPdry, TRACER_CP )
       CALC_R( Rmoist, qdry, q, k, i, j, iqw, Rdry, TRACER_R )
       r_cva = 1.0_RP / cva(k,i,j)
       r_cpa = 1.0_RP / cpa

       ! Coefficient of latent heat release for ssw change by PLCdep and PLRdep
       aliqliq          = 1.0_RP &
               + r_cva*( LHV00              + (CVvap-CL)*tem(k,i,j) )*dqswdtem_rho(k,i,j)
       ! Coefficient of latent heat release for ssw change by PLIdep, PLSdep and PLGdep
       asolliq          = 1.0_RP &
               + r_cva*( LHV00 + LHF00 + (CVvap-CI)*tem(k,i,j) )*dqswdtem_rho(k,i,j)
       ! Coefficient of latent heat release for ssi change by PLCdep and PLRdep
       aliqsol          = 1.0_RP &
               + r_cva*( LHV00              + (CVvap-CL)*tem(k,i,j) )*dqsidtem_rho(k,i,j)
       ! Coefficient of latent heat release for ssi change by PLIdep, PLSdep and PLGdep
       asolsol          = 1.0_RP &
               + r_cva*( LHV00 + LHF00 + (CVvap-CI)*tem(k,i,j) )*dqsidtem_rho(k,i,j)
       Pdynliq          = w * GRAV * ( r_cpa*dqswdtem_pre(k,i,j) + rho(k,i,j)*dqswdpre_tem(k,i,j) )
       Pdynsol          = w * GRAV * ( r_cpa*dqsidtem_pre(k,i,j) + rho(k,i,j)*dqsidpre_tem(k,i,j) )
       Pradliq          = -dTdt_rad(k,i,j)    * dqswdtem_rho(k,i,j)
       Pradsol          = -dTdt_rad(k,i,j)    * dqsidtem_rho(k,i,j)

       ssw_o            = ssw
       ssi_o            = ssi
!       ssw_o            = ssw - Pdynliq*(dt_dyn-dt_mp)/qsw(k,i,j) + Pradliq*r_qsw*dt_mp
!       ssi_o            = ssi - Pdynsol*(dt_dyn-dt_mp)/qsi(k,i,j) + Pradsol*r_qsi*dt_mp

       Acnd             = Pdynliq + Pradliq &
               - ( r_taudep_i+r_taudep_s+r_taudep_g ) * ( qsw(k,i,j) - qsi(k,i,j) )
       Adep             = Pdynsol + Pradsol &
               + ( r_taucnd_c+r_taucnd_r )            * ( qsw(k,i,j) - qsi(k,i,j) )
       r_taucnd         = &
               + aliqliq*( r_taucnd_c+r_taucnd_r ) &
               + asolliq*( r_taudep_i+r_taudep_s+r_taudep_g )
       r_taudep         = &
               + aliqsol*( r_taucnd_c+r_taucnd_r )&
               + asolsol*( r_taudep_i+r_taudep_s+r_taudep_g )

       uplim_cnd        = max( rho(k,i,j)*ssw_o*qsw(k,i,j)*r_dt, 0.0_RP )
       lowlim_cnd       = min( rho(k,i,j)*ssw_o*qsw(k,i,j)*r_dt, 0.0_RP )
       if( r_taucnd < r_tau100day )then
!          taucnd            = tau100day
          PQ(I_LCdep,k,i,j) = max(lowlim_cnd, min(uplim_cnd, PQ(I_LCdep,k,i,j)*ssw_o ))
          PQ(I_LRdep,k,i,j) = max(lowlim_cnd, min(uplim_cnd, PQ(I_LRdep,k,i,j)*ssw_o ))
          PQ(I_NRdep,k,i,j) = min(0.0_RP, PQ(I_NRdep,k,i,j)*ssw_o )
!          PLR2NR           = 0.0_RP
       else
          taucnd     = 1.0_RP/r_taucnd
          ! Production term for liquid water content
          coef_a_cnd = rho(k,i,j)*Acnd*taucnd
          coef_b_cnd = rho(k,i,j)*taucnd*r_dt*(ssw_o*qsw(k,i,j)-Acnd*taucnd) * ( exp(-dt*r_taucnd) - 1.0_RP )
          PQ(I_LCdep,k,i,j) = coef_a_cnd*r_taucnd_c - coef_b_cnd*r_taucnd_c
          PLR2NR            = PQ(I_NRdep,k,i,j)/(PQ(I_LRdep,k,i,j)+1.E-30_RP)
          PQ(I_LRdep,k,i,j) = coef_a_cnd*r_taucnd_r - coef_b_cnd*r_taucnd_r
          PQ(I_NRdep,k,i,j) = min(0.0_RP, PQ(I_LRdep,k,i,j)*PLR2NR )
       end if

       uplim_dep        = max( rho(k,i,j)*ssi_o*qsi(k,i,j)*r_dt, 0.0_RP )
       lowlim_dep       = min( rho(k,i,j)*ssi_o*qsi(k,i,j)*r_dt, 0.0_RP )
       if( r_taudep < r_tau100day )then
!          taudep            = tau100day
          PQ(I_LIdep,k,i,j) = max(lowlim_dep, min(uplim_dep, PQ(I_LIdep,k,i,j)*ssi_o ))
          PQ(I_LSdep,k,i,j) = max(lowlim_dep, min(uplim_dep, PQ(I_LSdep,k,i,j)*ssi_o ))
          PQ(I_LGdep,k,i,j) = max(lowlim_dep, min(uplim_dep, PQ(I_LGdep,k,i,j)*ssi_o ))
          PQ(I_NIdep,k,i,j) = min(0.0_RP, PQ(I_NIdep,k,i,j)*ssi_o )
          PQ(I_NSdep,k,i,j) = min(0.0_RP, PQ(I_NSdep,k,i,j)*ssi_o )
          PQ(I_NGdep,k,i,j) = min(0.0_RP, PQ(I_NGdep,k,i,j)*ssi_o )
       else
          taudep     = 1.0_RP/r_taudep
          ! Production term for ice water content
          coef_a_dep = rho(k,i,j)*Adep*taudep
          coef_b_dep = rho(k,i,j)*taudep*r_dt*(ssi_o*qsi(k,i,j)-Adep*taudep) * ( exp(-dt*r_taudep) - 1.0_RP )
          PLI2NI           = PQ(I_NIdep,k,i,j)/max(PQ(I_LIdep,k,i,j),1.E-30_RP)
          PLS2NS           = PQ(I_NSdep,k,i,j)/max(PQ(I_LSdep,k,i,j),1.E-30_RP)
          PLG2NG           = PQ(I_NGdep,k,i,j)/max(PQ(I_LGdep,k,i,j),1.E-30_RP)
          PQ(I_LIdep,k,i,j) =  coef_a_dep*r_taudep_i - coef_b_dep*r_taudep_i
          PQ(I_LSdep,k,i,j) =  coef_a_dep*r_taudep_s - coef_b_dep*r_taudep_s
          PQ(I_LGdep,k,i,j) =  coef_a_dep*r_taudep_g - coef_b_dep*r_taudep_g
          PQ(I_NIdep,k,i,j) = min(0.0_RP, PQ(I_LIdep,k,i,j)*PLI2NI )
          PQ(I_NSdep,k,i,j) = min(0.0_RP, PQ(I_LSdep,k,i,j)*PLS2NS )
          PQ(I_NGdep,k,i,j) = min(0.0_RP, PQ(I_LGdep,k,i,j)*PLG2NG )
       end if

       sw = 0.5_RP - sign(0.5_RP, PQ(I_LCdep,k,i,j)+eps) != 1 for PLCdep<=-eps
       PNCdep = min(0.0_RP, ((rhoq2(I_QC,k,i,j)+PQ(I_LCdep,k,i,j)*dt)*r_xc_ccn - rhoq2(I_NC,k,i,j))*r_dt ) * sw
!       if( PQ(I_LCdep,k,i,j) < -eps )then
!          PNCdep = min(0.0_RP, ((rhoq2(I_QC,k,i,j)+PQ(I_LCdep,k,i,j)*dt)*r_xc_ccn - rhoq2(I_NC,k,i,j))*r_dt )
!       else
!          PNCdep = 0.0_RP
!       end if
!       if( PQ(I_LIdep,k,i,j) < -eps )then
!          PQ(I_NIdep,k,i,j) = min(0.0_RP, ((li(k,i,j)+PQ(I_LIdep,k,i,j)*dt)*r_xi_ccn - rhoq2(I_NI,k,i,j))*r_dt )
!       else
!          PQ(I_NIdep,k,i,j) = 0.0_RP
!       end if

       !--- evaporation/condensation
       r_rvaptem = 1.0_RP/(Rvap*wtem(k,i,j))
       lvsw    = esw(k,i,j)*r_rvaptem
       dlvsw   = rhoq2(I_QV,k,i,j)-lvsw
       dcnd    = dt*(PQ(I_LCdep,k,i,j)+PQ(I_LRdep,k,i,j))

       sw = ( sign(0.5_RP,dcnd) + sign(0.5_RP,dlvsw) ) &
          * ( 0.5_RP + sign(0.5_RP,abs(dcnd)-eps) ) ! to avoid zero division
       ! sw= 1: always supersaturated
       ! sw=-1: always unsaturated
       ! sw= 0: partially unsaturated during timestep
       fac1 = min(dlvsw*sw,dcnd*sw)*sw / (abs(sw)-1.0_RP+dcnd) & ! sw=1,-1
            + 1.0_RP - abs(sw)                                   ! sw=0
       dep_dqc = max( dt*PQ(I_LCdep,k,i,j)*fac1, &
                     -rhoq2(I_QC,k,i,j) - 1e30_RP*(sw+1.0_RP) )*abs(sw) != -lc for sw=-1, -inf for sw=1
       dep_dqr = max( dt*PQ(I_LRdep,k,i,j)*fac1, &
                     -rhoq2(I_QR,k,i,j) - 1e30_RP*(sw+1.0_RP) )*abs(sw) != -lr for sw=-1, -inf for sw=1
!       if     ( (dcnd >  eps) .AND. (dlvsw > eps) )then
!          ! always supersaturated
!          fac1    = min(dlvsw,dcnd)/dcnd
!          dep_dqc =  dt*PQ(I_LCdep,k,i,j)*fac1
!          dep_dqr =  dt*PQ(I_LRdep,k,i,j)*fac1
!       else if( (dcnd < -eps) .AND. (dlvsw < -eps) )then
!          ! always unsaturated
!          fac1    = max( dlvsw,dcnd )/dcnd
!          dep_dqc = max( dt*PQ(I_LCdep,k,i,j)*fac1, -rhoq2(I_QC,k,i,j) )
!          dep_dqr = max( dt*PQ(I_LRdep,k,i,j)*fac1, -rhoq2(I_QR,k,i,j) )
!       else
!          ! partially unsaturated during timestep
!          fac1    = 1.0_RP
!          dep_dqc = 0.0_RP
!          dep_dqr = 0.0_RP
!       end if

       ! evaporation always lose number(always negative).
       dep_dnc = max( dt*PNCdep*fac1, -rhoq2(I_NC,k,i,j) ) ! ss>0 dep=0, ss<0 dep<0 ! [Add] 11/08/30 T.Mitsui
       dep_dnr = max( dt*PQ(I_NRdep,k,i,j)*fac1, -rhoq2(I_NR,k,i,j) ) ! ss>0 dep=0, ss<0 dep<0

       qc_evaporate(k,i,j) = - dep_dnc ! [Add] Y.Sato 15/09/08

       !--- deposition/sublimation
       lvsi    = esi(k,i,j)*r_rvaptem
       ddep    = dt*(PQ(I_LIdep,k,i,j)+PQ(I_LSdep,k,i,j)+PQ(I_LGdep,k,i,j))
       dlvsi   = rhoq2(I_QV,k,i,j)-lvsi  ! limiter for esi>1.d0

       sw = ( sign(0.5_RP,ddep) + sign(0.5_RP,dlvsi) ) &
          * ( 0.5_RP + sign(0.5_RP,abs(ddep)-eps) ) ! to avoid zero division
       ! sw= 1: always supersaturated
       ! sw=-1: always unsaturated
       ! sw= 0: partially unsaturated during timestep
       fac2 = min(dlvsi*sw,ddep*sw)*sw / (abs(sw)-1.0_RP+ddep) & ! sw=1,-1
            + 1.0_RP - abs(sw)                                   ! sw=0
       dep_dqi = max( dt*PQ(I_LIdep,k,i,j) &
                      * ( 1.0_RP-abs(sw) + fac2*abs(sw) ), & != fac2 for sw=-1,1, 1 for sw=0
                     -rhoq2(I_QI,k,i,j) - 1e30_RP*(sw+1.0_RP) )      != -li for sw=-1, -inf for sw=0,1
       dep_dqs = max( dt*PQ(I_LSdep,k,i,j) &
                      * ( 1.0_RP-abs(sw) + fac2*abs(sw) ), & != fac2 for sw=-1,1, 1 for sw=0
                     -rhoq2(I_QS,k,i,j) - 1e30_RP*(sw+1.0_RP) )      != -ls for sw=-1, -inf for sw=0,1
       dep_dqg = max( dt*PQ(I_LGdep,k,i,j) &
                      * ( 1.0_RP-abs(sw) + fac2*abs(sw) ), & != fac2 for sw=-1,1, 1 for sw=0
                     -rhoq2(I_QG,k,i,j) - 1e30_RP*(sw+1.0_RP) )      != -lg for sw=-1, -inf for sw=0,1
!       if      ( (ddep >  eps) .AND. (dlvsi > eps) )then
!          ! always supersaturated
!          fac2    = min(dlvsi,ddep)/ddep
!          dep_dqi = dt*PQ(I_LIdep,k,i,j)*fac2
!          dep_dqs = dt*PQ(I_LSdep,k,i,j)*fac2
!          dep_dqg = dt*PQ(I_LGdep,k,i,j)*fac2
!       else if ( (ddep < -eps) .AND. (dlvsi < -eps) )then
!          ! always unsaturated
!          fac2    = max(dlvsi,ddep)/ddep
!          dep_dqi = max(dt*PQ(I_LIdep,k,i,j)*fac2, -rhoq2(I_QI,k,i,j) )
!          dep_dqs = max(dt*PQ(I_LSdep,k,i,j)*fac2, -rhoq2(I_QS,k,i,j) )
!          dep_dqg = max(dt*PQ(I_LGdep,k,i,j)*fac2, -rhoq2(I_QG,k,i,j) )
!       else
!          ! partially unsaturated during timestep
!          fac2    = 1.0_RP
!          dep_dqi = dt*PQ(I_LIdep,k,i,j)
!          dep_dqs = dt*PQ(I_LSdep,k,i,j)
!          dep_dqg = dt*PQ(I_LGdep,k,i,j)
!       end if

       ! evaporation always lose number(always negative).
       dep_dni = max( dt*PQ(I_NIdep,k,i,j)*fac2, -rhoq2(I_NI,k,i,j) ) ! ss>0 dep=0, ss<0 dep<0
       dep_dns = max( dt*PQ(I_NSdep,k,i,j)*fac2, -rhoq2(I_NS,k,i,j) ) ! ss>0 dep=0, ss<0 dep<0
       dep_dng = max( dt*PQ(I_NGdep,k,i,j)*fac2, -rhoq2(I_NG,k,i,j) ) ! ss>0 dep=0, ss<0 dep<0

       !--- freezing of cloud drop
       frz_dqc = max( dt*(PQ(I_LChom,k,i,j)+PQ(I_LChet,k,i,j)), -rhoq2(I_QC,k,i,j)-dep_dqc ) ! negative value
       frz_dnc = max( dt*(PQ(I_NChom,k,i,j)+PQ(I_NChet,k,i,j)), -rhoq2(I_NC,k,i,j)-dep_dnc ) ! negative value
       fac3    = ( frz_dqc-eps )/( dt*(PQ(I_LChom,k,i,j)+PQ(I_LChet,k,i,j))-eps )
       fac4    = ( frz_dnc-eps )/( dt*(PQ(I_NChom,k,i,j)+PQ(I_NChet,k,i,j))-eps )
       PQ(I_LChom,k,i,j) = fac3*PQ(I_LChom,k,i,j)
       PQ(I_LChet,k,i,j) = fac3*PQ(I_LChet,k,i,j)
       PQ(I_NChom,k,i,j) = fac4*PQ(I_NChom,k,i,j)
       PQ(I_NChet,k,i,j) = fac4*PQ(I_NChet,k,i,j)

       !--- melting
       ! ice change
       mlt_dqi = max( dt*PQ(I_LImlt,k,i,j), -rhoq2(I_QI,k,i,j)-dep_dqi )  ! negative value
       mlt_dni = max( dt*PQ(I_NImlt,k,i,j), -rhoq2(I_NI,k,i,j)-dep_dni )  ! negative value
       ! snow change
       mlt_dqs = max( dt*PQ(I_LSmlt,k,i,j), -rhoq2(I_QS,k,i,j)-dep_dqs )  ! negative value
       mlt_dns = max( dt*PQ(I_NSmlt,k,i,j), -rhoq2(I_NS,k,i,j)-dep_dns )  ! negative value
       ! graupel change
       mlt_dqg = max( dt*PQ(I_LGmlt,k,i,j), -rhoq2(I_QG,k,i,j)-dep_dqg )  ! negative value
       mlt_dng = max( dt*PQ(I_NGmlt,k,i,j), -rhoq2(I_NG,k,i,j)-dep_dng )  ! negative value

       !--- freezing of larger droplets
       frz_dqr = max( dt*(PQ(I_LRhet,k,i,j)), min(0.0_RP, -rhoq2(I_QR,k,i,j)-dep_dqr) ) ! negative value
       frz_dnr = max( dt*(PQ(I_NRhet,k,i,j)), min(0.0_RP, -rhoq2(I_NR,k,i,j)-dep_dnr) ) ! negative value

       fac5         = ( frz_dqr-eps )/( dt*PQ(I_LRhet,k,i,j)-eps )
       PQ(I_LRhet,k,i,j) = fac5*PQ(I_LRhet,k,i,j)
       fac6         = ( frz_dnr-eps )/( dt*PQ(I_NRhet,k,i,j)-eps )
       PQ(I_NRhet,k,i,j) = fac6*PQ(I_NRhet,k,i,j)

       ! water vapor change
       drhoqv = -(dep_dqc+dep_dqi+dep_dqs+dep_dqg+dep_dqr)

       rhoq(I_QV,k,i,j) = max(0.0_RP, rhoq(I_QV,k,i,j) + drhoqv )

       rhoe(k,i,j) = rhoe(k,i,j) - LHV * drhoqv

       xi = min(xi_max, max(xi_min, rhoq2(I_QI,k,i,j)/(rhoq2(I_NI,k,i,j)+ni_min) ))
       sw = 0.5_RP + sign(0.5_RP,xi-x_sep) ! if (xi>=x_sep) then sw=1 else sw=0
                                                  ! sw=1: large ice crystals turn into rain by melting

       ! total cloud change
       drhoqc = ( frz_dqc - mlt_dqi*(1.0_RP-sw) + dep_dqc )
       drhonc = ( frz_dnc - mlt_dni*(1.0_RP-sw) + dep_dnc )
       ! total rain change
       drhoqr = ( frz_dqr - mlt_dqg - mlt_dqs - mlt_dqi*sw + dep_dqr )
       drhonr = ( frz_dnr - mlt_dng - mlt_dns - mlt_dni*sw + dep_dnr )

       rhoq(I_QC,k,i,j) = max(0.0_RP, rhoq(I_QC,k,i,j) + drhoqc )
       rhoq(I_NC,k,i,j) = max(0.0_RP, rhoq(I_NC,k,i,j) + drhonc )
       rhoq(I_QR,k,i,j) = max(0.0_RP, rhoq(I_QR,k,i,j) + drhoqr )
       rhoq(I_NR,k,i,j) = max(0.0_RP, rhoq(I_NR,k,i,j) + drhonr )

       ! total ice change
       drhoqi = (-frz_dqc + mlt_dqi             + dep_dqi )
       drhoni = (-frz_dnc + mlt_dni             + dep_dni )

       rhoq(I_QI,k,i,j) = max(0.0_RP, rhoq(I_QI,k,i,j) + drhoqi )
       rhoq(I_NI,k,i,j) = max(0.0_RP, rhoq(I_NI,k,i,j) + drhoni )

       rhoe(k,i,j) = rhoe(k,i,j) + LHF * drhoqi

       ! total snow change
       drhoqs = (           mlt_dqs             + dep_dqs )
       drhons = (           mlt_dns             + dep_dns )

       rhoq(I_QS,k,i,j) = max(0.0_RP, rhoq(I_QS,k,i,j) + drhoqs )
       rhoq(I_NS,k,i,j) = max(0.0_RP, rhoq(I_NS,k,i,j) + drhons )

       rhoe(k,i,j) = rhoe(k,i,j) + LHF * drhoqs

       ! total graupel change
       drhoqg = (-frz_dqr + mlt_dqg             + dep_dqg )
       drhong = (-frz_dnr + mlt_dng             + dep_dng )

       rhoq(I_QG,k,i,j) = max(0.0_RP, rhoq(I_QG,k,i,j) + drhoqg )
       rhoq(I_NG,k,i,j) = max(0.0_RP, rhoq(I_NG,k,i,j) + drhong )

       rhoe(k,i,j) = rhoe(k,i,j) + LHF * drhoqg

       !--- update mixing ratio
       rrho = 1.0_RP/rho(k,i,j)

       q(k,i,j,I_QV) = rhoq(I_QV,k,i,j) * rrho
       q(k,i,j,I_QC) = rhoq(I_QC,k,i,j) * rrho
       q(k,i,j,I_QR) = rhoq(I_QR,k,i,j) * rrho
       q(k,i,j,I_QI) = rhoq(I_QI,k,i,j) * rrho
       q(k,i,j,I_QS) = rhoq(I_QS,k,i,j) * rrho
       q(k,i,j,I_QG) = rhoq(I_QG,k,i,j) * rrho
       q(k,i,j,I_NC) = rhoq(I_NC,k,i,j) * rrho
       q(k,i,j,I_NR) = rhoq(I_NR,k,i,j) * rrho
       q(k,i,j,I_NI) = rhoq(I_NI,k,i,j) * rrho
       q(k,i,j,I_NS) = rhoq(I_NS,k,i,j) * rrho
       q(k,i,j,I_NG) = rhoq(I_NG,k,i,j) * rrho

       CALC_QDRY( qdry, q, TRACER_MASS, k, i, j, iqw )
       CALC_CV( cva(k,i,j), qdry, q, k, i, j, iqw, CVdry, TRACER_CV )
       CALC_R( Rmoist, qdry, q, k, i, j, iqw, Rdry, TRACER_R )
       tem(k,i,j) = rhoe(k,i,j) / ( rho(k,i,j) * cva(k,i,j) )
       pre(k,i,j) = rho(k,i,j) * Rmoist * tem(k,i,j)

       sl_PLCdep(i,j) = sl_PLCdep(i,j) + dep_dqc*Dz(k)*gsgam2(k,i,j)
       sl_PLRdep(i,j) = sl_PLRdep(i,j) + dep_dqr*Dz(k)*gsgam2(k,i,j)
       sl_PNRdep(i,j) = sl_PNRdep(i,j) + dep_dnr*Dz(k)*gsgam2(k,i,j)
    end do
    end do
    end do
    PROFILE_STOP("sn14_update")

    return
  end subroutine update_by_phase_change_kij
  !-------------------------------------------------------------------------------
  subroutine MP_negativefilter( &
       DENS, &
       QTRC  )
    use scale_tracer, only: &
       QA
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP) :: diffq(KA,IA,JA)
    real(RP) :: r_xmin

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_filter', 3)

    r_xmin = 1.0_RP / xmin_filter

    ! total hydrometeor (before correction)
    do j = JS, JE
    do i = IS, IE

    do k = KS, KE
       diffq(k,i,j) = QTRC(k,i,j,I_QV) &
                    + QTRC(k,i,j,I_QC) &
                    + QTRC(k,i,j,I_QR) &
                    + QTRC(k,i,j,I_QI) &
                    + QTRC(k,i,j,I_QS) &
                    + QTRC(k,i,j,I_QG)
    enddo

    ! remove negative value of hydrometeor (mass, number)
    do iq = I_QC, I_NG
    do k  = KS, KE
       QTRC(k,i,j,iq) = max(0.0_RP, QTRC(k,i,j,iq))
    enddo
    enddo

    ! apply correction of hydrometeor to total density
    ! [note] mass conservation is broken here to fill rounding error.
    do k  = KS, KE
       DENS(k,i,j) = DENS(k,i,j)        &
                   * ( 1.0_RP           &
                     + QTRC(k,i,j,I_QV) &
                     + QTRC(k,i,j,I_QC) &
                     + QTRC(k,i,j,I_QR) &
                     + QTRC(k,i,j,I_QI) &
                     + QTRC(k,i,j,I_QS) &
                     + QTRC(k,i,j,I_QG) &
                     - diffq(k,i,j)      ) ! after-before
    enddo

    ! avoid unrealistical value of number concentration
    ! due to numerical diffusion in advection

    do k  = KS, KE
       if ( QTRC(k,i,j,I_NC) > QTRC(k,i,j,I_QC)*r_xmin ) then
          QTRC(k,i,j,I_NC) = QTRC(k,i,j,I_QC)*r_xmin
       endif
    enddo
    do k  = KS, KE
       if ( QTRC(k,i,j,I_NR) > QTRC(k,i,j,I_QR)*r_xmin ) then
          QTRC(k,i,j,I_NR) = QTRC(k,i,j,I_QR)*r_xmin
       endif
    enddo
    do k  = KS, KE
       if ( QTRC(k,i,j,I_NI) > QTRC(k,i,j,I_QI)*r_xmin ) then
          QTRC(k,i,j,I_NI) = QTRC(k,i,j,I_QI)*r_xmin
       endif
    enddo
    do k  = KS, KE
       if ( QTRC(k,i,j,I_NS) > QTRC(k,i,j,I_QS)*r_xmin ) then
          QTRC(k,i,j,I_NS) = QTRC(k,i,j,I_QS)*r_xmin
       endif
    enddo
    do k  = KS, KE
       if ( QTRC(k,i,j,I_NG) > QTRC(k,i,j,I_QG)*r_xmin ) then
          QTRC(k,i,j,I_NG) = QTRC(k,i,j,I_QG)*r_xmin
       endif
    enddo

    enddo
    enddo

    call PROF_rapend('MP_filter', 3)

    return
  end subroutine MP_negativefilter
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sn14_CloudFraction( &
       cldfrac,       &
       QTRC,          &
       mask_criterion )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA)
    real(RP), intent(in)  :: mask_criterion

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do iq = I_QC, I_QG
          qhydro = qhydro + QTRC(k,i,j,iq)
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-mask_criterion)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_sn14_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_sn14_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0, &
       TEMP0  )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius          [cm]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)     ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)        ! density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)        ! temperature               [K]

    ! mass concentration[kg/m3] and mean particle mass[kg]
    real(RP) :: xc(KA,IA,JA)
    real(RP) :: xr(KA,IA,JA)
    real(RP) :: xi(KA,IA,JA)
    real(RP) :: xs(KA,IA,JA)
    real(RP) :: xg(KA,IA,JA)
    ! diameter of average mass[kg/m3]
    real(RP) :: dc_ave(KA,IA,JA)
    real(RP) :: dr_ave(KA,IA,JA)
    ! radius of average mass
    real(RP) :: rc, rr
    ! 2nd. and 3rd. order moment of DSD
    real(RP) :: ri2m(KA,IA,JA), ri3m(KA,IA,JA)
    real(RP) :: rs2m(KA,IA,JA), rs3m(KA,IA,JA)
    real(RP) :: rg2m(KA,IA,JA), rg3m(KA,IA,JA)

    real(RP) :: coef_Fuetal1998
    ! r2m_min is minimum value(moment of 1 particle with 1 micron)
    real(RP), parameter :: r2m_min=1.E-12_RP
    real(RP), parameter :: um2cm = 100.0_RP

    real(RP) :: limitsw, zerosw
    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! mean particle mass[kg]
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       xc(k,i,j) = min( xc_max, max( xc_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QC)/(QTRC0(k,i,j,I_NC)+nc_min) ) )
       xr(k,i,j) = min( xr_max, max( xr_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QR)/(QTRC0(k,i,j,I_NR)+nr_min) ) )
       xi(k,i,j) = min( xi_max, max( xi_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QI)/(QTRC0(k,i,j,I_NI)+ni_min) ) )
       xs(k,i,j) = min( xs_max, max( xs_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QS)/(QTRC0(k,i,j,I_NS)+ns_min) ) )
       xg(k,i,j) = min( xg_max, max( xg_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QG)/(QTRC0(k,i,j,I_NG)+ng_min) ) )
    enddo
    enddo
    enddo

    ! diameter of average mass : SB06 eq.(32)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dc_ave(k,i,j) = a_m(I_mp_QC) * xc(k,i,j)**b_m(I_mp_QC)
       dr_ave(k,i,j) = a_m(I_mp_QR) * xr(k,i,j)**b_m(I_mp_QR)
    enddo
    enddo
    enddo

    ! cloud effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rc = 0.5_RP * dc_ave(k,i,j)
       limitsw = 0.5_RP + sign(0.5_RP, rc-rmin_re )
       Re(k,i,j,I_HC) = coef_re(I_mp_QC) * rc * limitsw * um2cm
    enddo
    enddo
    enddo

    ! rain effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rr = 0.5_RP * dr_ave(k,i,j)
       limitsw = 0.5_RP + sign(0.5_RP, rr-rmin_re )
       Re(k,i,j,I_HR) = coef_re(I_mp_QR) * rr * limitsw * um2cm
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ri2m(k,i,j) = PI * coef_rea2(I_mp_QI) * QTRC0(k,i,j,I_NI) * a_rea2(I_mp_QI) * xi(k,i,j)**b_rea2(I_mp_QI)
       rs2m(k,i,j) = PI * coef_rea2(I_mp_QS) * QTRC0(k,i,j,I_NS) * a_rea2(I_mp_QS) * xs(k,i,j)**b_rea2(I_mp_QS)
       rg2m(k,i,j) = PI * coef_rea2(I_mp_QG) * QTRC0(k,i,j,I_NG) * a_rea2(I_mp_QG) * xg(k,i,j)**b_rea2(I_mp_QG)
    enddo
    enddo
    enddo

    ! Fu(1996), eq.(3.11) or Fu et al.(1998), eq.(2.5)
    coef_Fuetal1998 = 3.0_RP / (4.0_RP*rhoi)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ri3m(k,i,j) = coef_Fuetal1998 * QTRC0(k,i,j,I_NI) * xi(k,i,j)
       rs3m(k,i,j) = coef_Fuetal1998 * QTRC0(k,i,j,I_NS) * xs(k,i,j)
       rg3m(k,i,j) = coef_Fuetal1998 * QTRC0(k,i,j,I_NG) * xg(k,i,j)
    enddo
    enddo
    enddo

    ! ice effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       zerosw = 0.5_RP - sign(0.5_RP, ri2m(k,i,j) - r2m_min )
       Re(k,i,j,I_HI) = ri3m(k,i,j) / ( ri2m(k,i,j) + zerosw ) * ( 1.0_RP - zerosw ) * um2cm
    enddo
    enddo
    enddo

    ! snow effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       zerosw = 0.5_RP - sign(0.5_RP, rs2m(k,i,j) - r2m_min )
       Re(k,i,j,I_HS) = rs3m(k,i,j) / ( rs2m(k,i,j) + zerosw ) * ( 1.0_RP - zerosw ) * um2cm
    enddo
    enddo
    enddo

    ! graupel effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       zerosw = 0.5_RP - sign(0.5_RP, rg2m(k,i,j) - r2m_min )
       Re(k,i,j,I_HG) = rg3m(k,i,j) / ( rg2m(k,i,j) + zerosw ) * ( 1.0_RP - zerosw ) * um2cm
    enddo
    enddo
    enddo

    Re(:,:,:,I_HG+1:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_sn14_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sn14_Mixingratio( &
       Qe,   &
       QTRC0 )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometeor, only: &
       N_HYD
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,N_HYD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)     ! tracer mass concentration [kg/kg]

    integer  :: ihydro, iqa
    !---------------------------------------------------------------------------


    Qe(:,:,:,I_HC) = QTRC0(:,:,:,I_QC)
    Qe(:,:,:,I_HR) = QTRC0(:,:,:,I_QR)
    Qe(:,:,:,I_HI) = QTRC0(:,:,:,I_QI)
    Qe(:,:,:,I_HS) = QTRC0(:,:,:,I_QS)
    Qe(:,:,:,I_HG) = QTRC0(:,:,:,I_QG)
    Qe(:,:,:,I_HG+1:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_sn14_Mixingratio

end module scale_atmos_phy_mp_sn14
