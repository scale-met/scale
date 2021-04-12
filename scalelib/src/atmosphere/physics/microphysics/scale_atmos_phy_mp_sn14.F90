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
!<
!-------------------------------------------------------------------------------

#include "scalelib.h"
module scale_atmos_phy_mp_sn14
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

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
     DWATR  => CONST_DWATR,  &
     SMALL  => CONST_EPS

  use scale_atmos_hydrometeor, only: &
     N_HYD

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_sn14_setup
  public :: ATMOS_PHY_MP_sn14_finalize
  public :: ATMOS_PHY_MP_sn14_tendency
  public :: ATMOS_PHY_MP_sn14_terminal_velocity
  public :: ATMOS_PHY_MP_sn14_effective_radius
  public :: ATMOS_PHY_MP_sn14_cloud_fraction
  public :: ATMOS_PHY_MP_sn14_qtrc2qhyd
  public :: ATMOS_PHY_MP_sn14_qtrc2nhyd
  public :: ATMOS_PHY_MP_sn14_qhyd2qtrc

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: QA_MP  = 11

  integer,                parameter, public :: ATMOS_PHY_MP_SN14_ntracers = QA_MP
  integer,                parameter, public :: ATMOS_PHY_MP_SN14_nwaters = 2
  integer,                parameter, public :: ATMOS_PHY_MP_SN14_nices = 3
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_SN14_tracer_names(QA_MP) = (/ &
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
       'NG'  /)
  character(len=H_MID)  , parameter, public :: ATMOS_PHY_MP_SN14_tracer_descriptions(QA_MP) = (/ &
       'Ratio of Water Vapor mass to total mass (Specific humidity)', &
       'Ratio of Cloud Water mass to total mass                    ', &
       'Ratio of Rain Water mass to total mass                     ', &
       'Ratio of Cloud Ice mass ratio to total mass                ', &
       'Ratio of Snow miass ratio to total mass                    ', &
       'Ratio of Graupel mass ratio to total mass                  ', &
       'Cloud Water Number Density                                 ', &
       'Rain Water Number Density                                  ', &
       'Cloud Ice Number Density                                   ', &
       'Snow Number Density                                        ', &
       'Graupel Number Density                                     '/)
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_SN14_tracer_units(QA_MP) = (/ &
       'kg/kg ', &
       'kg/kg ', &
       'kg/kg ', &
       'kg/kg ', &
       'kg/kg ', &
       'kg/kg ', &
       'num/kg', &
       'num/kg', &
       'num/kg', &
       'num/kg', &
       'num/kg'  /)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: mp_sn14_init
  private :: mp_sn14

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  integer,  private, parameter   :: I_QV = 1
  integer,  private, parameter   :: I_QC = 2
  integer,  private, parameter   :: I_QR = 3
  integer,  private, parameter   :: I_QI = 4
  integer,  private, parameter   :: I_QS = 5
  integer,  private, parameter   :: I_QG = 6
  integer,  private, parameter   :: I_NC = 7
  integer,  private, parameter   :: I_NR = 8
  integer,  private, parameter   :: I_NI = 9
  integer,  private, parameter   :: I_NS = 10
  integer,  private, parameter   :: I_NG = 11

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

  ! for charge density
  integer, parameter :: I_CGNGacNS2NG = 25
  integer, parameter :: I_CGNGacNI2NG = 26
  integer, parameter :: I_NGspl = 49
  integer, parameter :: I_NSspl = 50
  integer, parameter :: Pcrg_MAX = 26

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
  real(RP), private, save :: a_m(HYDRO_MAX), log_a_m(HYDRO_MAX)
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
  real(RP), private, save :: coef_vt0(HYDRO_MAX,2), log_coef_vt0(HYDRO_MAX,2)
  real(RP), private, save :: coef_vt1(HYDRO_MAX,2), log_coef_vt1(HYDRO_MAX,2)
  real(RP), private, save :: coef_deplc
  real(RP), private, save :: coef_dave_N(HYDRO_MAX), log_coef_dave_N(HYDRO_MAX)
  real(RP), private, save :: coef_dave_L(HYDRO_MAX), log_coef_dave_L(HYDRO_MAX)
  ! diameter of terminal velocity branch
  !
  real(RP), private, save :: d0_ni=261.76E-6_RP, log_d0_ni
  real(RP), private, save :: d0_li=398.54E-6_RP, log_d0_li
  real(RP), private, parameter :: d0_ns=270.03E-6_RP, log_d0_ns = log(d0_ns)
  real(RP), private, parameter :: d0_ls=397.47E-6_RP, log_d0_ls = log(d0_ls)
  real(RP), private, parameter :: d0_ng=269.08E-6_RP, log_d0_ng = log(d0_ng)
  real(RP), private, parameter :: d0_lg=376.36E-6_RP, log_d0_lg = log(d0_lg)
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
  logical, private, save :: opt_debug_inc=.true.
  logical, private, save :: opt_debug_act=.true.
  logical, private, save :: opt_debug_ree=.true.
  logical, private, save :: opt_debug_bcs=.true.

  logical, private, save :: MP_doautoconversion = .true.
  logical, private, save :: MP_couple_aerosol   = .false. ! apply CCN effect?
  real(RP), private, save :: MP_ssw_lim = 1.E+1_RP

  !
  ! namelist variables for nucleation
  !
  ! total aerosol number concentration [/m3]
  real(RP), private, parameter :: c_ccn_ocean= 1.00E+8_RP
  real(RP), private, parameter :: c_ccn_land = 1.26E+9_RP
  real(RP), private, save      :: c_ccn      = 1.00E+8_RP
  ! aerosol activation factor
  real(RP), private, parameter :: kappa_ocean= 0.462_RP
  real(RP), private, parameter :: kappa_land = 0.308_RP
  real(RP), private, save      :: kappa      = 0.462_RP
  real(RP), private, save      :: c_in       = 1.0_RP
  ! SB06 (36)
  real(RP), private, save :: nm_M92 = 1.E+3_RP
  real(RP), private, save :: am_M92 = -0.639_RP
  real(RP), private, save :: bm_M92 = 12.96_RP
  !
  real(RP), private, save :: in_max = 1000.E+3_RP ! max num. of Ice-Nuclei [num/m3]
  real(RP), private, save :: ssi_max= 0.60_RP
  real(RP), private, save :: ssw_max= 1.1_RP  ! [%]
  !
  real(RP), private, save :: qke_min = 0.03_RP ! sigma=0.1[m/s], 09/08/18 T.Mitsui
  real(RP), private, save :: tem_ccn_low=233.150_RP  ! = -40 degC  ! [Add] 10/08/03 T.Mitsui
  real(RP), private, save :: tem_in_low =173.150_RP  ! = -100 degC ! [Add] 10/08/03 T.Mitsui
  logical,  private, save :: nucl_twomey = .false.
  logical,  private, save :: inucl_w     = .false.


  ! for incomplete gamma function
  real(RP), private, parameter :: rc_cr= 12.E-6_RP ! critical size[micron]
  real(RP), private, save :: xc_cr    ! mass[kg] of cloud with r=critical size[micron]
  real(RP), private, save :: alpha    ! slope parameter of gamma function
  real(RP), private, save :: gm, lgm  ! gamma(alpha), log(gamma(alpha))


  !=== for collection
  !--- threshold of diameters to collide with others
  real(RP), private, save :: dc0 =  15.0E-6_RP       ! lower threshold of cloud
  real(RP), private, save :: dc1 =  40.0E-6_RP       ! upper threshold of cloud
  real(RP), private, save :: di0 = 150.0E-6_RP       ! lower threshold of cloud
  real(RP), private, save :: ds0 = 150.0E-6_RP       ! lower threshold of cloud
  real(RP), private, save :: dg0 = 150.0E-6_RP       ! lower threshold of cloud
  !--- standard deviation of terminal velocity[m/s]
  real(RP), private, save :: sigma_c=0.00_RP        ! cloud
  real(RP), private, save :: sigma_r=0.00_RP        ! rain
  real(RP), private, save :: sigma_i=0.2_RP        ! ice
  real(RP), private, save :: sigma_s=0.2_RP        ! snow
  real(RP), private, save :: sigma_g=0.00_RP        ! graupel
  !--- max collection efficiency for cloud
  real(RP), private, save :: E_im = 0.80_RP        ! ice max
  real(RP), private, save :: E_sm = 0.80_RP        ! snow max
  real(RP), private, save :: E_gm = 1.00_RP        ! graupel max
  !--- collection efficiency between 2 species
  real(RP), private, save :: E_ir=1.0_RP            ! ice     x rain
  real(RP), private, save :: E_sr=1.0_RP            ! snow    x rain
  real(RP), private, save :: E_gr=1.0_RP            ! graupel x rain
  real(RP), private, save :: E_ii=1.0_RP            ! ice     x ice
  real(RP), private, save :: E_si=1.0_RP            ! snow    x ice
  real(RP), private, save :: E_gi=1.0_RP            ! graupel x ice
  real(RP), private, save :: E_ss=1.0_RP            ! snow    x snow
  real(RP), private, save :: E_gs=1.0_RP            ! graupel x snow
  real(RP), private, save :: E_gg=1.0_RP            ! graupel x graupel
  !=== for partial conversion
  !--- flag: 1=> partial conversion to graupel, 0=> no conversion
  integer,  private, save :: i_iconv2g=1          ! ice  => graupel
  integer,  private, save :: i_sconv2g=1          ! snow => graupel
  !--- bulk density of graupel
  real(RP), private, save :: rho_g   = 900.0_RP     ! [kg/m3]
  !--- space filling coefficient [%]
  real(RP), private, save :: cfill_i = 0.68_RP     ! ice
  real(RP), private, save :: cfill_s = 0.01_RP     ! snow
  !--- critical diameter for ice conversion
  real(RP), private, save :: di_cri  = 500.E-6_RP    ! [m]
  logical,  private, save :: opt_stick_KS96=.false.
  logical,  private, save :: opt_stick_CO86=.false.


  real(RP), private, save :: fac_cndc        = 1.0_RP
  logical,  private, save :: opt_fix_taucnd_c=.false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_sn14_setup
  !! setup
  !<
  subroutine ATMOS_PHY_MP_sn14_setup( &
    KA, IA, JA )
    use scale_prc, only: &
       PRC_abort
    implicit none

    integer, intent(in) :: KA
    integer, intent(in) :: IA
    integer, intent(in) :: JA

    namelist / PARAM_ATMOS_PHY_MP_SN14 / &
       MP_doautoconversion, &
       MP_ssw_lim,          &
       MP_couple_aerosol

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_sn14_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_MP_sn14_setup",*) 'Seiki and Nakajima (2014) 2-moment bulk 6 category'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SN14,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sn14_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sn14_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SN14. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_SN14)

    WLABEL(1) = "CLOUD"
    WLABEL(2) = "RAIN"
    WLABEL(3) = "ICE"
    WLABEL(4) = "SNOW"
    WLABEL(5) = "GRAUPEL"

    call mp_sn14_init

    allocate(nc_uplim_d(1,IA,JA))
    nc_uplim_d(:,:,:) = 150.E6_RP

    return
  end subroutine ATMOS_PHY_MP_sn14_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_MP_sn14_finalize

    deallocate(nc_uplim_d)

    return
  end subroutine ATMOS_PHY_MP_sn14_finalize

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_sn14_tendency
  !! calculate tendency
  !<
  subroutine ATMOS_PHY_MP_sn14_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, &
       W,    &
       QTRC, &
       PRES, &
       TEMP, &
       Qdry, &
       CPtot, &
       CVtot, &
       CCN, &
       dt, &
       cz, &
       fz, &
       RHOQ_t, &
       RHOE_t, &
       CPtot_t, &
       CVtot_t, &
       EVAPORATE, &
       flg_lt, &
       d0_crg, v0_crg, &
       dqcrg, &
       beta_crg, &
       QTRC_crg, &
       QSPLT_in, Sarea, &
       RHOQcrg_t      )
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS     (KA,IA,JA)
    real(RP), intent(in) :: W        (KA,IA,JA)
    real(RP), intent(in) :: QTRC     (KA,IA,JA,QA_MP)
    real(RP), intent(in) :: PRES(KA,IA,JA)
    real(RP), intent(in) :: TEMP(KA,IA,JA)
    real(RP), intent(in) :: Qdry(KA,IA,JA)
    real(RP), intent(in) :: CPtot(KA,IA,JA)
    real(RP), intent(in) :: CVtot(KA,IA,JA)
    real(RP), intent(in) :: CCN      (KA,IA,JA)
    real(DP), intent(in) :: dt
    real(RP), intent(in) :: cz(  KA,IA,JA)
    real(RP), intent(in) :: fz(0:KA,IA,JA)

    real(RP), intent(out) :: RHOQ_t   (KA,IA,JA,QA_MP)
    real(RP), intent(out) :: RHOE_t   (KA,IA,JA)
    real(RP), intent(out) :: CPtot_t(KA,IA,JA)
    real(RP), intent(out) :: CVtot_t(KA,IA,JA)
    real(RP), intent(out) :: EVAPORATE(KA,IA,JA)   !--- number of evaporated cloud [/m3]

    ! Optional for Lightning
    logical,  intent(in), optional :: flg_lt
    real(RP), intent(in), optional :: d0_crg, v0_crg
    real(RP), intent(in), optional :: dqcrg(KA,IA,JA)
    real(RP), intent(in), optional :: beta_crg(KA,IA,JA)
    real(RP), intent(in), optional :: QTRC_crg(KA,IA,JA,HYDRO_MAX)
    real(RP), intent(out), optional :: QSPLT_in(KA,IA,JA,3)
    real(RP), intent(out), optional :: Sarea(KA,IA,JA,HYDRO_MAX)
    real(RP), intent(out), optional :: RHOQcrg_t(KA,IA,JA,HYDRO_MAX)
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / microphysics / SN14'

#ifdef PROFILE_FIPP
    call fipp_start()
#endif

    !##### MP Main #####
    call MP_sn14 ( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS(:,:,:), W(:,:,:), QTRC(:,:,:,:), PRES(:,:,:), TEMP(:,:,:), & ! (in)
       Qdry(:,:,:), CPtot(:,:,:), CVtot(:,:,:), CCN(:,:,:),            & ! (in)
       real(dt,RP), cz(:,:,:), fz(:,:,:),                              & ! (in)
       RHOQ_t(:,:,:,:), RHOE_t(:,:,:), CPtot_t(:,:,:), CVtot_t(:,:,:), & ! (out)
       EVAPORATE(:,:,:),                                               & ! (out)
       flg_lt, d0_crg, v0_crg, dqcrg(:,:,:), beta_crg(:,:,:),          & ! (optional in)
       QTRC_crg(:,:,:,:),                                              & ! (optional in)
       QSPLT_in(:,:,:,:), Sarea(:,:,:,:), RHOQcrg_t(:,:,:,:)           ) ! (optional out)

#ifdef PROFILE_FIPP
    call fipp_stop()
#endif

    return
  end subroutine ATMOS_PHY_MP_sn14_tendency

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_sn14_cloud_fraction
  !! Calculate Cloud Fraction
  !<
  subroutine ATMOS_PHY_MP_sn14_cloud_fraction( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC,           &
       mask_criterion, &
       cldfrac         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA_MP-1)
    real(RP), intent(in)  :: mask_criterion

    real(RP), intent(out) :: cldfrac(KA,IA,JA)

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    !$omp parallel do &
    !$omp private(qhydro)
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do iq = 1, QA_MP-1
          qhydro = qhydro + QTRC(k,i,j,iq)
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-mask_criterion)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_sn14_cloud_fraction
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_sn14_effective_radius
  !! Calculate Effective Radius
  !<
  subroutine ATMOS_PHY_MP_sn14_effective_radius( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS0, TEMP0, QTRC0, &
       Re                   )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: DENS0(KA,IA,JA)        ! density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)        ! temperature               [K]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,I_QC:I_NG)     ! tracer mass concentration [kg/kg]

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius          [cm]

    ! mass concentration[kg/m3] and mean particle mass[kg]
    real(RP) :: xc(KA)
    real(RP) :: xr(KA)
    real(RP) :: xi(KA)
    real(RP) :: xs(KA)
    real(RP) :: xg(KA)
    ! diameter of average mass[kg/m3]
    real(RP) :: dc_ave(KA)
    real(RP) :: dr_ave(KA)
    ! radius of average mass
    real(RP) :: rc, rr
    ! 2nd. and 3rd. order moment of DSD
    real(RP) :: ri2m(KA), ri3m(KA)
    real(RP) :: rs2m(KA), rs3m(KA)
    real(RP) :: rg2m(KA), rg3m(KA)

    real(RP), parameter :: coef_Fuetal1998 = 3.0_RP / (4.0_RP*rhoi)

    ! r2m_min is minimum value(moment of 1 particle with 1 micron)
    real(RP), parameter :: r2m_min=1.E-12_RP
    real(RP), parameter :: um2cm = 100.0_RP

    real(RP) :: limitsw, zerosw
    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do &
    !$omp private(xc,xr,xi,xs,xg,dc_ave,dr_ave,rc,rr,ri2m,ri3m,rs2m,rs3m,rg2m,rg3m, &
    !$omp         limitsw,zerosw)
    do j  = JS, JE
    do i  = IS, IE

       ! mean particle mass[kg]
       do k  = KS, KE
          xc(k) = min( xc_max, max( xc_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QC)/(QTRC0(k,i,j,I_NC)+nc_min) ) )
          xr(k) = min( xr_max, max( xr_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QR)/(QTRC0(k,i,j,I_NR)+nr_min) ) )
          xi(k) = min( xi_max, max( xi_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QI)/(QTRC0(k,i,j,I_NI)+ni_min) ) )
          xs(k) = min( xs_max, max( xs_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QS)/(QTRC0(k,i,j,I_NS)+ns_min) ) )
          xg(k) = min( xg_max, max( xg_min, DENS0(k,i,j)*QTRC0(k,i,j,I_QG)/(QTRC0(k,i,j,I_NG)+ng_min) ) )
       enddo

       ! diameter of average mass : SB06 eq.(32)
       do k = KS, KE
          dc_ave(k) = a_m(I_mp_QC) * xc(k)**b_m(I_mp_QC)
          dr_ave(k) = a_m(I_mp_QR) * xr(k)**b_m(I_mp_QR)
       enddo

       ! cloud effective radius
       do k = KS, KE
          rc = 0.5_RP * dc_ave(k)
          limitsw = 0.5_RP + sign(0.5_RP, rc-rmin_re )
          Re(k,i,j,I_HC) = coef_re(I_mp_QC) * rc * limitsw * um2cm
       enddo

       ! rain effective radius
       do k = KS, KE
          rr = 0.5_RP * dr_ave(k)
          limitsw = 0.5_RP + sign(0.5_RP, rr-rmin_re )
          Re(k,i,j,I_HR) = coef_re(I_mp_QR) * rr * limitsw * um2cm
       enddo

       do k = KS, KE
          ri2m(k) = PI * coef_rea2(I_mp_QI) * QTRC0(k,i,j,I_NI) * a_rea2(I_mp_QI) * xi(k)**b_rea2(I_mp_QI)
          rs2m(k) = PI * coef_rea2(I_mp_QS) * QTRC0(k,i,j,I_NS) * a_rea2(I_mp_QS) * xs(k)**b_rea2(I_mp_QS)
          rg2m(k) = PI * coef_rea2(I_mp_QG) * QTRC0(k,i,j,I_NG) * a_rea2(I_mp_QG) * xg(k)**b_rea2(I_mp_QG)
       enddo

       ! Fu(1996), eq.(3.11) or Fu et al.(1998), eq.(2.5)
       do k = KS, KE
          ri3m(k) = coef_Fuetal1998 * QTRC0(k,i,j,I_NI) * xi(k)
          rs3m(k) = coef_Fuetal1998 * QTRC0(k,i,j,I_NS) * xs(k)
          rg3m(k) = coef_Fuetal1998 * QTRC0(k,i,j,I_NG) * xg(k)
       enddo

       ! ice effective radius
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, ri2m(k) - r2m_min )
          Re(k,i,j,I_HI) = ri3m(k) / ( ri2m(k) + zerosw ) * ( 1.0_RP - zerosw ) * um2cm
       enddo

       ! snow effective radius
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, rs2m(k) - r2m_min )
          Re(k,i,j,I_HS) = rs3m(k) / ( rs2m(k) + zerosw ) * ( 1.0_RP - zerosw ) * um2cm
       enddo

       ! graupel effective radius
       do k = KS, KE
          zerosw = 0.5_RP - sign(0.5_RP, rg2m(k) - r2m_min )
          Re(k,i,j,I_HG) = rg3m(k) / ( rg2m(k) + zerosw ) * ( 1.0_RP - zerosw ) * um2cm
       enddo

       do k = KS, KE
          Re(k,i,j,I_HH) = 0.0_RP
       end do

    enddo
    enddo


    return
  end subroutine ATMOS_PHY_MP_sn14_effective_radius
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_sn14_qtrc2qhyd
  !! Calculate mass ratio of each category
  !<
  subroutine ATMOS_PHY_MP_sn14_qtrc2qhyd( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC0, &
       Qe     )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA_MP-1) ! tracer mass concentration [kg/kg]

    real(RP), intent(out) :: Qe   (KA,IA,JA,N_HYD)   ! mixing ratio of each cateory [kg/kg]

    integer :: k, i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Qe(k,i,j,I_HC) = QTRC0(k,i,j,I_mp_QC)
       Qe(k,i,j,I_HR) = QTRC0(k,i,j,I_mp_QR)
       Qe(k,i,j,I_HI) = QTRC0(k,i,j,I_mp_QI)
       Qe(k,i,j,I_HS) = QTRC0(k,i,j,I_mp_QS)
       Qe(k,i,j,I_HG) = QTRC0(k,i,j,I_mp_QG)
       Qe(k,i,j,I_HH) = 0.0_RP
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_MP_sn14_qtrc2qhyd

  !> Calculate number concentration of each category
  subroutine ATMOS_PHY_MP_sn14_qtrc2nhyd( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC0, &
       Ne     )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA_MP-1) ! tracer mass concentration [kg/kg]

    real(RP), intent(out) :: Ne   (KA,IA,JA,N_HYD)   ! number density of each cateory [1/m3]

    integer :: k, i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       Ne(k,i,j,I_HC) = QTRC0(k,i,j,I_mp_NC)
       Ne(k,i,j,I_HR) = QTRC0(k,i,j,I_mp_NR)
       Ne(k,i,j,I_HI) = QTRC0(k,i,j,I_mp_NI)
       Ne(k,i,j,I_HS) = QTRC0(k,i,j,I_mp_NS)
       Ne(k,i,j,I_HG) = QTRC0(k,i,j,I_mp_NG)
       Ne(k,i,j,I_HH) = 0.0_RP
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_MP_sn14_qtrc2nhyd

  subroutine ATMOS_PHY_MP_sn14_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Qe, &
       QTRC,     &
       QNUM      )
    use scale_const, only: &
       PI    => CONST_PI, &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: Qe(KA,IA,JA,N_HYD) ! mass ratio of each cateory [kg/kg]

    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP-1) ! tracer mass concentration [kg/kg]

    real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD)

    real(RP), parameter :: Dc   =  20.E-6_RP ! typical particle diameter for cloud  [m]
    real(RP), parameter :: Dr   = 200.E-6_RP ! typical particle diameter for rain   [m]
    real(RP), parameter :: Di   =  80.E-6_RP ! typical particle diameter for ice    [m]
    real(RP), parameter :: Ds   =  80.E-6_RP ! typical particle diameter for snow   [m]
    real(RP), parameter :: Dg   = 200.E-6_RP ! typical particle diameter for grapel [m]
    real(RP), parameter :: RHOw =  1000.0_RP ! typical density for warm particles   [kg/m3]
    real(RP), parameter :: RHOf =   100.0_RP ! typical density for frozen particles [kg/m3]
    real(RP), parameter :: RHOg =   400.0_RP ! typical density for grapel particles [kg/m3]
    real(RP), parameter :: b    =     3.0_RP ! assume spherical form

    real(RP) :: piov6

    integer :: k, i, j
    !---------------------------------------------------------------------------


!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_mp_QC) = Qe(k,i,j,I_HC)
    end do
    end do
    end do

!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_mp_QR) = Qe(k,i,j,I_HR)
    end do
    end do
    end do

!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_mp_QI) = Qe(k,i,j,I_HI)
    end do
    end do
    end do

!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_mp_QS) = Qe(k,i,j,I_HS)
    end do
    end do
    end do

!OCL XFILL
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC(k,i,j,I_mp_QG) = Qe(k,i,j,I_HG) + Qe(k,i,j,I_HH)
    end do
    end do
    end do

    piov6 = PI / 6.0_RP

    if ( present(QNUM) ) then

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( QNUM(k,i,j,I_HC) .ne. UNDEF ) then
             QTRC(k,i,j,I_mp_NC) = QNUM(k,i,j,I_HC)
          else
             QTRC(k,i,j,I_mp_NC) = QTRC(k,i,j,I_mp_QC) / ( (piov6*RHOw) * Dc**b )
          end if
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( QNUM(k,i,j,I_HR) .ne. UNDEF ) then
             QTRC(k,i,j,I_mp_NR) = QNUM(k,i,j,I_HR)
          else
             QTRC(k,i,j,I_mp_NR) = QTRC(k,i,j,I_mp_QR) / ( (piov6*RHOw) * Dr**b )
          end if
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( QNUM(k,i,j,I_HI) .ne. UNDEF ) then
             QTRC(k,i,j,I_mp_NI) = QNUM(k,i,j,I_HI)
          else
             QTRC(k,i,j,I_mp_NI) = QTRC(k,i,j,I_mp_QI) / ( (piov6*RHOf) * Di**b )
          end if
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( QNUM(k,i,j,I_HS) .ne. UNDEF ) then
             QTRC(k,i,j,I_mp_NS) = QNUM(k,i,j,I_HS)
          else
             QTRC(k,i,j,I_mp_NS) = QTRC(k,i,j,I_mp_QS) / ( (piov6*RHOf) * Ds**b )
          end if
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          if ( QNUM(k,i,j,I_HG) .ne. UNDEF ) then
             if ( QNUM(k,i,j,I_HH) .ne. UNDEF ) then
                QTRC(k,i,j,I_mp_NG) = QNUM(k,i,j,I_HG) + QNUM(k,i,j,I_HH)
             else
                QTRC(k,i,j,I_mp_NG) = QNUM(k,i,j,I_HG)
             end if
          else
             QTRC(k,i,j,I_mp_NG) = QTRC(k,i,j,I_mp_QG) / ( (piov6*RHOg) * Dg**b )
          end if
       end do
       end do
       end do

    else

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_mp_NC) = QTRC(k,i,j,I_mp_QC) / ( (piov6*RHOw) * Dc**b )
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_mp_NR) = QTRC(k,i,j,I_mp_QR) / ( (piov6*RHOw) * Dr**b )
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_mp_NI) = QTRC(k,i,j,I_mp_QI) / ( (piov6*RHOf) * Di**b )
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_mp_NS) = QTRC(k,i,j,I_mp_QS) / ( (piov6*RHOf) * Ds**b )
       end do
       end do
       end do

!OCL XFILL
       !$omp parallel do
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,I_mp_NG) = QTRC(k,i,j,I_mp_QG) / ( (piov6*RHOg) * Dg**b )
       end do
       end do
       end do

    end if

    return
  end subroutine ATMOS_PHY_MP_sn14_qhyd2qtrc

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_sn14_terminal_velocity
  !! Calculate terminal velocity
  !<
!OCL SERIAL
  subroutine ATMOS_PHY_MP_sn14_terminal_velocity( &
       KA, KS, KE, &
      DENS, &
      TEMP, &
      RHOQ, &
      PRES, &
      vterm )
    use scale_const, only: &
       CONST_UNDEF
    implicit none

    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: RHOQ(KA,I_QC:I_NG) ! rho * q
    real(RP), intent(in)  :: DENS(KA)    ! rho
    real(RP), intent(in)  :: TEMP(KA)    ! temperature
    real(RP), intent(in)  :: PRES(KA)    ! pressure

    real(RP), intent(out) :: vterm(KA,QA_MP-1) ! terminal velocity of cloud mass

    real(RP) :: xq, log_xq ! average mass of 1 particle( mass/number )

    real(RP) :: rhofac   ! density factor for terminal velocity
    real(RP) :: rhofac_q(KA), log_rhofac_q(KA)

    real(RP) :: rlambdar(KA) ! work for diagnosis of Rain DSD ( Seifert, 2008 )
    real(RP) :: mud_r
    real(RP) :: dq, log_dq   ! weigthed diameter.   Improved Rogers etal. (1993) formula by T.Mitsui


    real(RP) :: weight, weightk(KA) ! weighting coefficient for 2-branches is determined by ratio between 0.745mm and weighted diameter.  SB06 Table.1
    real(RP) :: velq_s ! terminal velocity for small branch of Rogers formula
    real(RP) :: velq_l ! terminal velocity for large branch of Rogers formula

    real(RP) :: tmp
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    ! QC, NC
    do k = KS, KE
       rhofac = rho_0 / max( DENS(k), rho_min )

!       rhofac_q(k) = rhofac ** gamma_v(I_mp_QC)
!       xq = max( xqmin(I_mp_QC), min( xqmax(I_mp_QC), rhoq(k,I_QC) / ( rhoq(k,I_NC) + nqmin(I_mp_QC) ) ) )
       log_rhofac_q(k) = log(rhofac) * gamma_v(I_mp_QC)
       log_xq = log( max( xqmin(I_mp_QC), min( xqmax(I_mp_QC), rhoq(k,I_QC) / ( rhoq(k,I_NC) + nqmin(I_mp_QC) ) ) ) )

!       vterm(k,I_mp_QC) = -rhofac_q(k) * coef_vt1(I_mp_QC,1) * xq**beta_v(I_mp_QC,1)
       vterm(k,I_mp_QC) = - exp( log_rhofac_q(k) + log_coef_vt1(I_mp_QC,1) + log_xq * beta_v(I_mp_QC,1) )

!       vterm(k,I_mp_NC) = -rhofac_q(k) * coef_vt0(I_mp_QC,1) * xq**beta_vn(I_mp_QC,1)
       vterm(k,I_mp_NC) = - exp( log_rhofac_q(k) + log_coef_vt0(I_mp_QC,1) + log_xq * beta_vn(I_mp_QC,1) )
    end do

    ! QR, NR
    do k = KS, KE
       rhofac = rho_0 / max( DENS(k), rho_min )
       rhofac_q(k) = rhofac**gamma_v(I_mp_QR)
    end do
    mud_r = 3.0_RP * nu(I_mp_QR) + 2.0_RP
    do k = KS, KE
       xq = max( xqmin(I_mp_QR), min( xqmax(I_mp_QR), rhoq(k,I_QR) / ( rhoq(k,I_NR) + nqmin(I_mp_QR) ) ) )

       rlambdar(k) = a_m(I_mp_QR) * xq**b_m(I_mp_QR) &
                   * ( (mud_r+3.0_RP) * (mud_r+2.0_RP) * (mud_r+1.0_RP) )**(-0.333333333_RP)
    end do
    do k = KS, KE
       dq = ( 4.0_RP + mud_r ) * rlambdar(k) ! D^(3)+mu weighted mean diameter
       weightk(k) = 0.5_RP * ( 1.0_RP + tanh( PI * log( dq/d_vtr_branch ) ) )
    end do
    do k = KS, KE
       velq_s = coef_vtr_ar2 * dq &
              * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar(k) )**(-5.0_RP-mud_r) )
       velq_l = coef_vtr_ar1 &
              - coef_vtr_br1 * ( 1.0_RP + coef_vtr_cr1*rlambdar(k) )**(-4.0_RP-mud_r)
       weight = min( 1.0_RP, max( 0.0_RP, weightk(k) ) )
       vterm(k,I_mp_QR) = -rhofac_q(k) * ( velq_l * (          weight ) &
                                         + velq_s * ( 1.0_RP - weight ) )

    end do
    do k = KS, KE
       dq = ( 1.0_RP + mud_r ) * rlambdar(k)
       weightk(k) = 0.5_RP * ( 1.0_RP + tanh( PI * log( dq/d_vtr_branch ) ) )
    end do
    do k = KS, KE
       dq = ( 4.0_RP + mud_r ) * rlambdar(k)
       velq_s = coef_vtr_ar2 * dq &
              * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar(k) )**(-2.0_RP-mud_r) )
       velq_l = coef_vtr_ar1 - coef_vtr_br1 &
              * ( 1.0_RP + coef_vtr_cr1*rlambdar(k) )**(-1.0_RP-mud_r)
       weight = min( 1.0_RP, max( 0.0_RP, weightk(k) ) )
       vterm(k,I_mp_NR) = -rhofac_q(k) * ( velq_l * (          weight ) &
                                         + velq_s * ( 1.0_RP - weight ) )
    end do

    ! QI, NI
!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
!       rhofac_q(k) = ( PRES(k)/pre0_vt )**a_pre0_vt * ( TEMP(k)/tem0_vt )**a_tem0_vt
!       xq = max( xqmin(I_mp_QI), min( xqmax(I_mp_QI), rhoq(k,I_QI) / ( rhoq(k,I_NI) + nqmin(I_mp_QI) ) ) )
       log_rhofac_q(k) = log( PRES(k)/pre0_vt ) * a_pre0_vt + log( TEMP(k)/tem0_vt ) * a_tem0_vt
       log_xq = log( max( xqmin(I_mp_QI), min( xqmax(I_mp_QI), rhoq(k,I_QI) / ( rhoq(k,I_NI) + nqmin(I_mp_QI) ) ) ) )

!       tmp = a_m(I_mp_QI) * xq**b_m(I_mp_QI)
!       dq = coef_dave_L(I_mp_QI) * tmp
!       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_li ) ) ) )
       tmp = log_a_m(I_mp_QI) + log_xq * b_m(I_mp_QI)
       log_dq = log_coef_dave_L(I_mp_QI) + tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log_dq - log_d0_li ) ) )

!       velq_s = coef_vt1(I_mp_QI,1) * xq**beta_v (I_mp_QI,1)
!       velq_l = coef_vt1(I_mp_QI,2) * xq**beta_v (I_mp_QI,2)
!       vterm(k,I_mp_QI) = -rhofac_q(k) * ( velq_l * (          weight ) &
!                                         + velq_s * ( 1.0_RP - weight ) )
       velq_s = log_coef_vt1(I_mp_QI,1) + log_xq * beta_v(I_mp_QI,1)
       velq_l = log_coef_vt1(I_mp_QI,2) + log_xq * beta_v(I_mp_QI,2)
       vterm(k,I_mp_QI) = - exp( log_rhofac_q(k) + velq_l * (          weight ) &
                                                 + velq_s * ( 1.0_RP - weight ) )

!       dq = coef_dave_N(I_mp_QI) * tmp
!       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ni ) ) ) )
       log_dq = log_coef_dave_N(I_mp_QI) + tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log_dq - log_d0_ni ) ) )

!       velq_s = coef_vt0(I_mp_QI,1) * xq**beta_vn(I_mp_QI,1)
!       velq_l = coef_vt0(I_mp_QI,2) * xq**beta_vn(I_mp_QI,2)
!       vterm(k,I_mp_NI) = -rhofac_q(k) * ( velq_l * (          weight ) &
!                                         + velq_s * ( 1.0_RP - weight ) )
       velq_s = log_coef_vt0(I_mp_QI,1) + log_xq * beta_vn(I_mp_QI,1)
       velq_l = log_coef_vt0(I_mp_QI,2) + log_xq * beta_vn(I_mp_QI,2)
       vterm(k,I_mp_NI) = - exp( log_rhofac_q(k) + velq_l * (          weight ) &
                                                 + velq_s * ( 1.0_RP - weight ) )
    end do

    ! QS, NS
    do k = KS, KE
!       xq = max( xqmin(I_mp_QS), min( xqmax(I_mp_QS), rhoq(k,I_QS) / ( rhoq(k,I_NS) + nqmin(I_mp_QS) ) ) )
       log_xq = log( max( xqmin(I_mp_QS), min( xqmax(I_mp_QS), rhoq(k,I_QS) / ( rhoq(k,I_NS) + nqmin(I_mp_QS) ) ) ) )

!       tmp = a_m(I_mp_QS) * xq**b_m(I_mp_QS)
!       dq = coef_dave_L(I_mp_QS) * tmp
!       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ls ) ) ) )
       tmp = log_a_m(I_mp_QS) + log_xq * b_m(I_mp_QS)
       log_dq = log_coef_dave_L(I_mp_QS) + tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log_dq - log_d0_ls ) ) )

!       velq_s = coef_vt1(I_mp_QS,1) * xq**beta_v (I_mp_QS,1)
!       velq_l = coef_vt1(I_mp_QS,2) * xq**beta_v (I_mp_QS,2)
!       vterm(k,I_mp_QS) = -rhofac_q(k) * ( velq_l * (          weight ) &
!                                         + velq_s * ( 1.0_RP - weight ) )
       velq_s = log_coef_vt1(I_mp_QS,1) + log_xq * beta_v(I_mp_QS,1)
       velq_l = log_coef_vt1(I_mp_QS,2) + log_xq * beta_v(I_mp_QS,2)
       vterm(k,I_mp_QS) = - exp( log_rhofac_q(k) + velq_l * (          weight ) &
                                                 + velq_s * ( 1.0_RP - weight ) )

!       dq = coef_dave_N(I_mp_QS) * tmp
!       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ns ) ) ) )
       log_dq = log_coef_dave_N(I_mp_QS) + tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log_dq - log_d0_ns ) ) )

!       velq_s = coef_vt0(I_mp_QS,1) * xq**beta_vn(I_mp_QS,1)
!       velq_l = coef_vt0(I_mp_QS,2) * xq**beta_vn(I_mp_QS,2)
!       vterm(k,I_mp_NS) = -rhofac_q(k) * ( velq_l * (          weight ) &
!                                         + velq_s * ( 1.0_RP - weight ) )
       velq_s = log_coef_vt0(I_mp_QS,1) + log_xq * beta_vn(I_mp_QS,1)
       velq_l = log_coef_vt0(I_mp_QS,2) + log_xq * beta_vn(I_mp_QS,2)
       vterm(k,I_mp_NS) = - exp( log_rhofac_q(k) + velq_l * (          weight ) &
                                                 + velq_s * ( 1.0_RP - weight ) )
    end do

    ! QG, NG
    do k = KS, KE
!       xq = max( xqmin(I_mp_QG), min( xqmax(I_mp_QG), rhoq(k,I_QG) / ( rhoq(k,I_NG) + nqmin(I_mp_QG) ) ) )
       log_xq = log( max( xqmin(I_mp_QG), min( xqmax(I_mp_QG), rhoq(k,I_QG) / ( rhoq(k,I_NG) + nqmin(I_mp_QG) ) ) ) )

!       tmp = a_m(I_mp_QG) * xq**b_m(I_mp_QG)
!       dq = coef_dave_L(I_mp_QG) * tmp
!       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_lg ) ) ) )
       tmp = log_a_m(I_mp_QG) + log_xq * b_m(I_mp_QG)
       log_dq = log_coef_dave_L(I_mp_QG) + tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log_dq - log_d0_lg ) ) )

!       velq_s = coef_vt1(I_mp_QG,1) * xq**beta_v (I_mp_QG,1)
!       velq_l = coef_vt1(I_mp_QG,2) * xq**beta_v (I_mp_QG,2)
!       vterm(k,I_mp_QG) = -rhofac_q(k) * ( velq_l * (          weight ) &
!                                         + velq_s * ( 1.0_RP - weight ) )
       velq_s = log_coef_vt1(I_mp_QG,1) + log_xq * beta_v(I_mp_QG,1)
       velq_l = log_coef_vt1(I_mp_QG,2) + log_xq * beta_v(I_mp_QG,2)
       vterm(k,I_mp_QG) = - exp( log_rhofac_q(k) + velq_l * (          weight ) &
                                                 + velq_s * ( 1.0_RP - weight ) )

!       dq = coef_dave_N(I_mp_QG) * tmp
!       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log( dq/d0_ng ) ) ) )
       log_dq = log_coef_dave_N(I_mp_QG) + tmp
       weight = min( 1.0_RP, max( 0.0_RP, 0.5_RP * ( 1.0_RP + log_dq - log_d0_ng ) ) )

!       velq_s = coef_vt0(I_mp_QG,1) * xq**beta_vn(I_mp_QG,1)
!       velq_l = coef_vt0(I_mp_QG,2) * xq**beta_vn(I_mp_QG,2)
!       vterm(k,I_mp_NG) = -rhofac_q(k) * ( velq_l * (          weight ) &
!                                         + velq_s * ( 1.0_RP - weight ) )
       velq_s = log_coef_vt0(I_mp_QG,1) + log_xq * beta_vn(I_mp_QG,1)
       velq_l = log_coef_vt0(I_mp_QG,2) + log_xq * beta_vn(I_mp_QG,2)
       vterm(k,I_mp_NG) = - exp( log_rhofac_q(k) + velq_l * (          weight ) &
                                                 + velq_s * ( 1.0_RP - weight ) )
    enddo

    do iq = 1, QA_MP-1
       vterm(1:KS-2 ,iq) = CONST_UNDEF
       vterm(KS-1   ,iq) = vterm(KS,iq)
       vterm(KE+1:KA,iq) = CONST_UNDEF
    enddo

    return
  end subroutine ATMOS_PHY_MP_sn14_terminal_velocity


  ! private
  !-----------------------------------------------------------------------------
  subroutine mp_sn14_init
    use scale_prc, only: &
       PRC_abort
    use scale_specfunc, only: &
        gammafunc => SF_gamma
    implicit none

    real(RP), parameter :: eps_gamma=1.E-30_RP

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
    namelist / PARAM_ATMOS_PHY_MP_SN14_init / &
         opt_debug,                  &
         opt_debug_inc,              &
         opt_debug_act,              &
         opt_debug_ree,              &
         opt_debug_bcs,              &
         ntmax_phase_change,         &
         ntmax_collection
    !
    namelist / PARAM_ATMOS_PHY_MP_SN14_particles / &
         a_m, b_m, alpha_v, beta_v, gamma_v, &
         alpha_vn, beta_vn,    &
         a_area, b_area, cap,  &
         nu, mu,               &
         opt_M96_column_ice,   &
         opt_M96_ice,          &
         ar_ice_fix

    namelist / PARAM_ATMOS_PHY_MP_SN14_nucleation / &
         in_max,                     & !
         c_ccn, kappa,               & ! cloud nucleation
         nm_M92, am_M92, bm_M92,     & ! ice nucleation
         xc_ccn, xi_ccn,             &
         tem_ccn_low,                & ! [Add] 10/08/03 T.Mitsui
         tem_in_low,                 & ! [Add] 10/08/03 T.Mitsui
         ssw_max, ssi_max,           &
         nucl_twomey, inucl_w        ! [Add] 13/01/30 Y.Sato

    namelist / PARAM_ATMOS_PHY_MP_SN14_collection / &
         dc0, dc1, di0, ds0, dg0,    &
         sigma_c, sigma_r, sigma_i, sigma_s, sigma_g, &
         opt_stick_KS96,   &
         opt_stick_CO86,   &
         E_im, E_sm, E_gm, &
         E_ir, E_sr, E_gr, E_ii, E_si, E_gi, E_ss, E_gs, E_gg, &
         i_iconv2g, i_sconv2g, rho_g, cfill_i, cfill_s, di_cri

    !
    namelist / PARAM_ATMOS_PHY_MP_SN14_condensation / &
         opt_fix_taucnd_c, fac_cndc


    a_m(:)         = UNDEF8
    log_a_m(:)     = UNDEF8
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
    log_coef_dave_N(:) = UNDEF8
    log_coef_dave_L(:) = UNDEF8
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
    log_coef_vt0(:,:) = UNDEF8
    log_coef_vt1(:,:) = UNDEF8
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
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SN14_init,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sn14_init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sn14_init",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SN14_init. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_SN14_init)

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
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SN14_particles,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sn14_init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sn14_init",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SN14_particles. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_SN14_particles)

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
          coef_vt0(iw,ia) = alpha_vn(iw,ia) * w1(iw) / w2(iw) * ( w3(iw) / w4(iw) )**beta_vn(iw,ia)
          log_coef_vt0(iw,ia) = log( alpha_vn(iw,ia) * w1(iw) / w2(iw) ) + log( w3(iw) / w4(iw) ) * beta_vn(iw,ia)
          n = 1
          w1(iw) = gammafunc( (beta_v(iw,ia) + nu(iw) + 1.0_RP + n)/mu(iw) )
          w2(iw) = gammafunc( (                nu(iw) + 1.0_RP + n)/mu(iw) )
          ! coefficient of terminal velocity for mass
          coef_vt1(iw,ia) = alpha_v(iw,ia) * w1(iw) / w2(iw) * ( w3(iw) / w4(iw) )**beta_v(iw,ia)
          log_coef_vt1(iw,ia) = log( alpha_v(iw,ia) * w1(iw) / w2(iw) ) + log( w3(iw) / w4(iw) ) * beta_v(iw,ia)
       end do
    end do

    ! coefficient for weighted diameter used in calculation of terminal velocity
    do iw=1, HYDRO_MAX
       w1(iw) = gammafunc( (         b_m(iw) + nu(iw) + 1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (1.0_RP + b_m(iw) + nu(iw) + 1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w4(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       coef_dave_N(iw) = ( w1(iw) / w3(iw) ) * ( w3(iw) / w4(iw) )**(       b_m(iw))
       coef_dave_L(iw) = ( w2(iw) / w3(iw) ) * ( w3(iw) / w4(iw) )**(1.0_RP+b_m(iw))
       log_coef_dave_N(iw) = log( w1(iw) / w3(iw) ) + log( w3(iw) / w4(iw) ) * (       b_m(iw))
       log_coef_dave_L(iw) = log( w2(iw) / w3(iw) ) + log( w3(iw) / w4(iw) ) * (1.0_RP+b_m(iw))
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
!          w1(iw) = gammafunc( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,1) + n)/mu(iw) )
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

    rewind(IO_FID_CONF)
    read(IO_FID_CONF, nml=PARAM_ATMOS_PHY_MP_SN14_nucleation, iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sn14_init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sn14_init",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SN14_nucleation. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_SN14_nucleation)
    if ( MP_couple_aerosol .AND. nucl_twomey ) then
       LOG_ERROR("ATMOS_PHY_MP_SN14_nucleation_kij",*) "nucl_twomey should be false when MP_couple_aerosol is true, stop"
       call PRC_abort
    endif


    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_MP_SN14_collection, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sn14_init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sn14_init",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SN14_collection. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_SN14_collection)


    ! [Add] 11/08/30 T.Mitsui
    rewind(IO_FID_CONF)
    read  (IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SN14_condensation, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sn14_init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sn14_init",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SN14_condensation. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_SN14_condensation)



    ! work for Incomplete Gamma function
    xc_cr = (2.0_RP*rc_cr/a_m(I_mp_QC))**(1.0_RP/b_m(I_mp_QC))
    alpha = (nu(I_mp_QC)+1.0_RP)/mu(I_mp_QC)
    gm    = gammafunc(alpha)
    lgm   = log(gm)


    ! log
    do ia=1, HYDRO_MAX
       log_a_m(ia) = log( a_m(ia) )
    end do
    log_d0_li = log( d0_li )
    log_d0_ni = log( d0_ni )


    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(100a16)')      "LABEL       ",WLABEL(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "capacity    ",cap(:) ! [Add] 11/08/30 T.Mitsui
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_m2     ",coef_m2(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_d      ",coef_d(:)
    !
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_d3     ",coef_d3(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_d6     ",coef_d6(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_d2v    ",coef_d2v(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_md2v   ",coef_md2v(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "a_d2vt      ",a_d2vt(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "b_d2vt      ",b_d2vt(:)
    !
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_r2     ",coef_r2(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_r3     ",coef_r3(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_re     ",coef_re(:)
    !
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "a_area      ",a_area(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "b_area      ",b_area(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "ax_area     ",ax_area(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "bx_area     ",bx_area(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "a_rea       ",a_rea(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "b_rea       ",b_rea(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "a_rea3      ",a_rea3(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "b_rea3      ",b_rea3(:)
    !
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_rea2   ",coef_rea2(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_rea3   ",coef_rea3(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_vt0    ",coef_vt0(:,1)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_vt1    ",coef_vt1(:,1)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_A      ",coef_A(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "coef_lambda ",coef_lambda(:)

    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "ah_vent0 sml",ah_vent0(:,1)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "ah_vent0 lrg",ah_vent0(:,2)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "ah_vent1 sml",ah_vent1(:,1)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "ah_vent1 lrg",ah_vent1(:,2)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "bh_vent0 sml",bh_vent0(:,1)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "bh_vent0 lrg",bh_vent0(:,2)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "bh_vent1 sml",bh_vent1(:,1)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "bh_vent1 lrg",bh_vent1(:,2)

    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "delta_b0    ",delta_b0(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "delta_b1    ",delta_b1(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "theta_b0    ",theta_b0(:)
    LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,100ES16.6)') "theta_b1    ",theta_b1(:)

    do ia=1, HYDRO_MAX
       LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,a10,a,100ES16.6)') "delta0(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab0(ia,ib),ib=1,HYDRO_MAX)
    enddo
    do ia=1, HYDRO_MAX
       LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,a10,a,100ES16.6)') "delta1(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab1(ia,ib),ib=1,HYDRO_MAX)
    enddo
    do ia=1, HYDRO_MAX
       LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,a10,a,100ES16.6)') "theta0(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab0(ia,ib),ib=1,HYDRO_MAX)
    enddo
    do ia=1, HYDRO_MAX
       LOG_INFO("ATMOS_PHY_MP_sn14_init",'(a,a10,a,100ES16.6)') "theta1(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab1(ia,ib),ib=1,HYDRO_MAX)
    enddo

    return
  end subroutine mp_sn14_init
  !-----------------------------------------------------------------------------
  subroutine mp_sn14 ( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS,  &
       W,     &
       QTRC,  &
       PRES0, &
       TEMP0, &
       Qdry,  &
       CPtot0, &
       CVtot0, &
       CCN,    &
       dt,   &
       cz,   &
       fz,   &
       RHOQ_t, &
       RHOE_t, &
       CPtot_t, &
       CVtot_t, &
       EVAPORATE, &
       flg_lt, &
       d0_crg, &
       v0_crg, &
       dqcrg, &
       beta_crg, &
       QTRC_crg, &
       QSPLT_in, &
       Sarea, &
       RHOQcrg_t_mp )

    use scale_atmos_hydrometeor, only: &
       CP_VAPOR, &
       CP_WATER, &
       CP_ICE,   &
       CV_VAPOR, &
       CV_WATER, &
       CV_ICE
    use scale_atmos_saturation, only: &
       moist_psat_liq      => ATMOS_SATURATION_psat_liq,   &
       moist_psat_ice      => ATMOS_SATURATION_psat_ice
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS  (KA,IA,JA)
    real(RP), intent(in) :: W     (KA,IA,JA)
    real(RP), intent(in) :: QTRC  (KA,IA,JA,QA_MP)
    real(RP), intent(in) :: PRES0 (KA,IA,JA)
    real(RP), intent(in) :: TEMP0 (KA,IA,JA)
    real(RP), intent(in) :: Qdry  (KA,IA,JA)
    real(RP), intent(in) :: CPtot0(KA, IA, JA)
    real(RP), intent(in) :: CVtot0(KA, IA, JA)
    real(RP), intent(in) :: CCN   (KA,IA,JA)
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: cz(  KA,IA,JA)
    real(RP), intent(in) :: fz(0:KA,IA,JA)

    real(RP),intent(out) :: RHOQ_t(KA,IA,JA,QA_MP)
    real(RP),intent(out) :: RHOE_t(KA,IA,JA)
    real(RP),intent(out) :: CPtot_t(KA,IA,JA)
    real(RP),intent(out) :: CVtot_t(KA,IA,JA)

    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)   ! number of evaporated cloud [/m3/s]

    !--- for lightning
    logical,  intent(in), optional :: flg_lt
    real(RP), intent(in), optional :: d0_crg, v0_crg
    real(RP), intent(in), optional :: dqcrg(KA,IA,JA)
    real(RP), intent(in), optional :: beta_crg(KA,IA,JA)
    real(RP), intent(in), optional :: QTRC_crg(KA,IA,JA,HYDRO_MAX)
    real(RP), intent(out), optional :: QSPLT_in(KA,IA,JA,3)
    real(RP), intent(out), optional :: Sarea(KA,IA,JA,HYDRO_MAX)
    real(RP), intent(out), optional :: RHOQcrg_t_mp(KA,IA,JA,HYDRO_MAX)

    real(RP) :: pres (KA)
    real(RP) :: temp (KA)
    real(RP) :: cva  (KA)
    real(RP) :: cpa  (KA)
    real(RP) :: rrho (KA) !> 1/DENS
    real(RP) :: rhoe (KA)
    real(RP) :: rhoq (KA,I_QV:I_NG)
    real(RP) :: rhoq2(KA,I_QV:I_NG)
    !
    real(RP) :: RHOQ0_t (KA,QA_MP)
    real(RP) :: RHOE0_t (KA)
    real(RP) :: CPtot0_t(KA)
    real(RP) :: CVtot0_t(KA)
    !
    real(RP) :: xq(KA,HYDRO_MAX)
    !
    real(RP) :: dq_xa(KA,HYDRO_MAX)
    real(RP) :: vt_xa(KA,HYDRO_MAX,2) ! terminal velocity

    real(RP) :: wtemp(KA) ! filtered temperature
    real(RP) :: esw  (KA) ! saturated vapor pressure(water)
    real(RP) :: esi  (KA) ! saturated vapor pressure(ice)
    !
    real(RP) :: rho_fac
    real(RP) :: rho_fac_q(KA,HYDRO_MAX) ! factor for tracers, 1:cloud, 2:rain, 3:ice, 4: snow, 5:graupel
    !
    real(RP) :: drhoqv         ! d (rho*qv)
    real(RP) :: drhoqc, drhonc !        qc, nc
    real(RP) :: drhoqr, drhonr !        qr, nr
    real(RP) :: drhoqi, drhoni !        qi, ni
    real(RP) :: drhoqs, drhons !        qs, ns
    real(RP) :: drhoqg, drhong !        qg, ng

    ! production rate
    real(RP) :: PQ(KA,PQ_MAX)
    real(RP) :: wrm_dqc, wrm_dnc
    real(RP) :: wrm_dqr, wrm_dnr

    ! production rate of mixed-phase collection process
    real(RP) :: Pac(KA,Pac_MAX)
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

    !-----------------------------------------------
    ! work for explicit supersaturation modeling
    !-----------------------------------------------
    real(RP) :: dTdt_equiv_d(KA)
    !--------------------------------------------------
    !
    ! variables for output
    !
    !--------------------------------------------------
    ! work for column production term
    real(RP) :: sl_PLCdep
    real(RP) :: sl_PLRdep, sl_PNRdep
    !--------------------------------------------------
    real(RP) :: qke_d(KA)

    real(RP), parameter :: eps      = 1.E-19_RP
    real(RP), parameter :: eps_qv   = 1.E-19_RP
    real(RP), parameter :: eps_rhoe = 1.E-19_RP
    real(RP), parameter :: eps_rho  = 1.E-19_RP

    ! for limitter
    real(RP) :: di2l, dtem
    real(RP) :: fact

    real(RP) :: sw

    integer  :: k, i, j, iq

    real(RP) :: dqv, dql, dqi
    real(RP) :: dcv, dcp

    !---- for Lightning component
!    logical, private, save :: MP_doice_graupel_collection = .false.
!    real(RP), private :: flg_igcol = 0.0_RP
    !--- for Charge split
    real(RP) :: v0_crg_l=0.0_RP, d0_crg_l=0.0_RP
    real(RP) :: facq(I_QC:I_QG), f_crg
    integer :: grid(2), pp, qq
    real(RP) :: drhoqcrg_c, drhoqcrg_r
    real(RP) :: drhoqcrg_i, drhoqcrg_s, drhoqcrg_g

    ! production rate of charge density
    real(RP) :: Pcrg1(KA,PQ_MAX)
    real(RP) :: Pcrg2(KA,Pcrg_MAX)
    real(RP) :: rhoq_crg(KA,I_QC:I_QG)
    real(RP) :: rhoq2_crg(KA,I_QC:I_QG)
    real(RP) :: QTRC0(KA,QA_MP)
    real(RP) :: Crs(KA,HYDRO_MAX)
    real(RP),allocatable :: RHOQcrg0_t(:,:,:,:)

    real(RP) :: crg_split_s
    real(RP) :: crg_split_g
    real(RP) :: crg_split_i
    real(RP) :: wrm_dnc_crg
    real(RP) :: wrm_dnr_crg
    real(RP) :: gc_dnc_crg
    real(RP) :: sc_dnc_crg
    real(RP) :: ic_dnc_crg
    real(RP) :: rg_dng_crg
    real(RP) :: rg_dnr_crg
    real(RP) :: rs_dnr_crg
    real(RP) :: rs_dns_crg
    real(RP) :: ri_dnr_crg
    real(RP) :: ri_dni_crg
    real(RP) :: ii_dni_crg
    real(RP) :: is_dni_crg
    real(RP) :: ss_dns_crg
    real(RP) :: gs_dns_crg
    real(RP) :: gi_dni_crg
    real(RP) :: gg_dng_crg
    ! mixed-phase collection process total plus(clp_), total minus(clm_)
    real(RP) :: clp_dnc_crg, clm_dnc_crg
    real(RP) :: clp_dnr_crg, clm_dnr_crg
    real(RP) :: clp_dni_crg, clm_dni_crg
    real(RP) :: clp_dns_crg, clm_dns_crg
    real(RP) :: clp_dng_crg, clm_dng_crg
    ! production rate of partial conversion(ice, snow => graupel)
    real(RP) :: pco_dni_crg
    real(RP) :: pco_dns_crg
    real(RP) :: pco_dng_crg
    ! production rate of enhanced melting due to
    real(RP) :: eml_dnc_crg
    real(RP) :: eml_dnr_crg
    real(RP) :: eml_dni_crg
    real(RP) :: eml_dns_crg
    real(RP) :: eml_dng_crg
    ! production rate of ice multiplication by splintering
    real(RP) :: spl_dni_crg
    real(RP) :: spl_dns_crg
    real(RP) :: spl_dng_crg
    real(RP) :: rate1
    logical  :: flg_lt_l

    real(RP) :: sw1, sw2
    real(RP) :: tmp
    integer :: ip
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------

    if ( present(flg_lt) ) then
       flg_lt_l = flg_lt
    else
       flg_lt_l = .false.
    end if


    !--- Lightning component is on
    if( flg_lt_l ) then
!       flg_igcol = 0.0_RP
       d0_crg_l = d0_crg
       v0_crg_l = v0_crg

!OCL ZFILL
       !$omp workshare
       Qsplt_in(:,:,:,:) = 0.0_RP
       !$omp end workshare

       allocate(RHOQcrg0_t(KA,IA,JA,HYDRO_MAX))
!OCL ZFILL
       !$omp workshare
       RHOQcrg_t_mp(:,:,:,:) = 0.0_RP
       !$omp end workshare

!      if( MP_doice_graupel_collection ) then
!        flg_igcol = 1.0_RP
!      else
!        write(*,*) 'xxx MP_doice_graupel_collection should be true for TK78 Stop!'
!        call PRC_MPIstop
!        flg_igcol = 0.0_RP
!      endif
    endif

    !$omp parallel do default(none) &
    !$omp shared(KA,KS,KE,IS,IE,JS,JE, &
    !$omp        CP_VAPOR,CP_WATER,CP_ICE,CV_VAPOR,CV_WATER,CV_ICE,LHV,LHF,LHF0, &
    !$omp        MP_doautoconversion, &
    !$omp        DENS,W,QTRC,TEMP0,PRES0,QDRY,CPtot0,CVtot0,CCN, &
    !$omp        cz,fz,dt, &
    !$omp        RHOQ_t,RHOE_t,CPtot_t,CVtot_t,EVAPORATE, &
    !$omp        c_ccn,gamma_v,nc_uplim_d,a_m,b_m,alpha_v,beta_v,a_rea,b_rea,ntmax_phase_change, &
    !$omp        QTRC_crg,QSPLT_in,RHOQcrg_t_mp,RHOQcrg0_t,Sarea, &
    !$omp        d0_crg_l,v0_crg_l,beta_crg,dqcrg,flg_lt_l) &
    !$omp private (pres,temp,rrho,rhoe,rhoq,rhoq2,cva,cpa,rhoq0_t,rhoe0_t,cptot0_t,cvtot0_t, &
    !$omp          xq,dq_xa,vt_xa,wtemp,esw,esi,rho_fac,rho_fac_q, &
    !$omp          drhoqv,drhoqc,drhonc,drhoqr,drhonr,drhoqi,drhoni,drhoqs,drhons,drhoqg,drhong, &
    !$omp          PQ,Pac,wrm_dqc,wrm_dnc,wrm_dqr,wrm_dnr, &
    !$omp          gc_dqc,gc_dnc,sc_dqc,sc_dnc,ic_dqc,ic_dnc, &
    !$omp          rg_dqg,rg_dng,rg_dqr,rg_dnr,rs_dqr,rs_dnr,rs_dqs,rs_dns,ri_dqr,ri_dnr,ri_dqi,ri_dni, &
    !$omp          ii_dqi,ii_dni,is_dqi,is_dni,ss_dns,gs_dqs,gs_dns,gg_dng, &
    !$omp          clp_dqc,clp_dnc,clm_dqc,clm_dnc,clp_dqr,clp_dnr,clm_dqr,clm_dnr, &
    !$omp          clp_dqi,clp_dni,clm_dqi,clm_dni,clp_dqs,clp_dns,clm_dqs,clm_dns,clp_dqg,clp_dng,clm_dqg,clm_dng, &
    !$omp          rhoq_crg,rhoq2_crg,Pcrg1,Pcrg2, &
    !$omp          drhoqcrg_c,drhoqcrg_r,drhoqcrg_i,drhoqcrg_s,drhoqcrg_g, &
    !$omp          crg_split_s,crg_split_g,crg_split_i,wrm_dnc_crg,wrm_dnr_crg,QTRC0,crs, &
    !$omp          gc_dnc_crg,sc_dnc_crg,ic_dnc_crg,rg_dng_crg,rg_dnr_crg,rs_dnr_crg,rs_dns_crg,ri_dnr_crg,ri_dni_crg, &
    !$omp          ii_dni_crg,is_dni_crg,ss_dns_crg,gs_dns_crg,gi_dni_crg,gg_dng_crg, &
    !$omp          clp_dnc_crg,clm_dnc_crg,clp_dnr_crg,clm_dnr_crg,clp_dni_crg,clm_dni_crg,clp_dns_crg,clm_dns_crg,clp_dng_crg,clm_dng_crg, &
    !$omp          pco_dni_crg,pco_dns_crg,pco_dng_crg,eml_dnc_crg,eml_dnr_crg,eml_dni_crg,eml_dns_crg,eml_dng_crg,spl_dni_crg,spl_dns_crg,spl_dng_crg, &
    !$omp          fac1,fac3,fac4,fac6,fac7,fac9,fac10, &
    !$omp          pco_dqi,pco_dni,pco_dqs,pco_dns,pco_dqg,pco_dng, &
    !$omp          eml_dqc,eml_dnc,eml_dqr,eml_dnr,eml_dqi,eml_dni,eml_dqs,eml_dns,eml_dqg,eml_dng, &
    !$omp          spl_dqi,spl_dni,spl_dqg,spl_dqs, &
    !$omp          dTdt_equiv_d,sl_PLCdep,sl_PLRdep,sl_PNRdep,qke_d, &
    !$omp          dqv,dql,dqi,dcv,dcp, &
    !$omp          di2l,dtem,fact,sw,sw1,sw2,tmp )
    do j = JS, JE
    do i = IS, IE

       ! total tendency
       do k = KS, KE
          RHOQ_t(k,i,j,:) = 0.0_RP
          RHOE_t(k,i,j)   = 0.0_RP
          CPtot_t(k,i,j)  = 0.0_RP
          CVtot_t(k,i,j)  = 0.0_RP
       end do

       ! intermidiate variable
       do k = KS, KE
          cpa (k) = CPtot0(k,i,j)
          cva (k) = CVtot0(k,i,j)
          pres(k) = PRES0 (k,i,j)
          temp(k) = TEMP0 (k,i,j)
       enddo

       !============================================================================
       !
       !--  Each process is integrated sequentially.
       !    1. Nucleation and filter
       !    2. Phase change
       !    3. Collection
       !
       !============================================================================

       !----------------------------------------------------------------------------
       !
       ! 1.Nucleation of cloud water and cloud ice
       !
       !----------------------------------------------------------------------------
       do iq = I_QV, I_NG
       do k = KS, KE
          rhoq(k,iq) = DENS(k,i,j) * QTRC(k,i,j,iq)
       enddo
       enddo
       do k = KS, KE
          rhoq2(k,I_QV) = DENS(k,i,j)*QTRC(k,i,j,I_QV)
          rhoq2(k,I_NI) = max( 0.0_RP, DENS(k,i,j)*QTRC(k,i,j,I_NI) )
          rhoq2(k,I_NC) = max( 0.0_RP, DENS(k,i,j)*QTRC(k,i,j,I_NC) )
       enddo

       do k = KS, KE
          rrho (k) = 1.0_RP / DENS(k,i,j)
          rhoe (k) = DENS(k,i,j) * temp(k) * cva(k)
          wtemp(k) = max(temp(k), tem_min)
       enddo

#ifdef DEBUG
       call debug_tem( KA, KS, KE, &
            1, i, j, temp(:), DENS(:,i,j), pres(:), QTRC(:,i,j,I_QV) )
#endif

       do k = KS, KE
          rho_fac              = rho_0 / max(DENS(k,i,j),rho_min)
          rho_fac_q(k,I_mp_QC) = rho_fac**gamma_v(I_mp_QC)
          rho_fac_q(k,I_mp_QR) = rho_fac**gamma_v(I_mp_QR)
          rho_fac_q(k,I_mp_QI) = (pres(k)/pre0_vt)**a_pre0_vt * (temp(k)/tem0_vt)**a_tem0_vt
          rho_fac_q(k,I_mp_QS) = rho_fac_q(k,I_mp_QI)
          rho_fac_q(k,I_mp_QG) = rho_fac_q(k,I_mp_QI)
       enddo


       sl_PLCdep = 0.0_RP
       sl_PLRdep = 0.0_RP
       sl_PNRdep = 0.0_RP

!OCL XFILL
       do k = KS, KE
          qke_d(k) = 0.0_RP ! 2*TKE
       enddo

!OCL XFILL
       do k = KS, KE
          dTdt_equiv_d(k) = 0.0_RP
       enddo

!       nc_uplim_d(1)  = c_ccn_map(1,i,j)*1.5_RP
       nc_uplim_d(1,i,j)  = c_ccn*1.5_RP

       call nucleation( &
            KA, KS, KE, &
            cz(:,i,j), fz(:,i,j),           & ! (in)
            w(:,i,j), DENS(:,i,j),          & ! (in)
            wtemp(:), pres(:), qdry(:,i,j), & ! (in)
            rhoq2(:,:), cpa(:),             & ! (in)
            dTdt_equiv_d(:),                & ! (in)
            qke_d(:),                       & ! (in)
            CCN(:,i,j), nc_uplim_d(1,i,j),  & ! (in)
            dt,                             & ! (in)
            PQ(:,:)                         ) ! (out)

       ! tendency

       ! nucleation
       do k = KS, KE
          drhoqc = PQ(k,I_LCccn)
          drhoqi = PQ(k,I_LIccn)
          tmp = - drhoqc - drhoqi
          drhoqv = max( - rhoq(k,I_QV) / dt , tmp )
          fac1 = drhoqv / min( tmp, -eps ) ! drhoqc and drhoqi must be >= 0, otherwise fac1 can be artificially huge value.

          drhoqc = drhoqc * fac1
          drhoqi = drhoqi * fac1

          RHOQ0_t(k,I_QV) = drhoqv
          RHOQ0_t(k,I_QC) = drhoqc
          RHOQ0_t(k,I_QI) = drhoqi

          RHOE0_t(k) = - LHV * drhoqv + LHF * drhoqi

          dqv = rrho(k) * drhoqv
          dql = rrho(k) * drhoqc
          dqi = rrho(k) * drhoqi

          dcv = CV_VAPOR * dqv + CV_WATER * dql + CV_ICE * dqi
          dcp = CP_VAPOR * dqv + CP_WATER * dql + CP_ICE * dqi

          CVtot0_t(k) = dcv
          CPtot0_t(k) = dcp

          drhonc = PQ(k,I_NCccn) * fac1
          drhoni = PQ(k,I_NIccn) * fac1
          RHOQ0_t(k,I_NC) = drhonc
          RHOQ0_t(k,I_NI) = drhoni
       end do

       ! total tendency
       do k = KS, KE
          RHOE_t (k,i,j) = RHOE_t (k,i,j) + RHOE0_t (k)
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CVtot0_t(k)
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CPtot0_t(k)
       enddo

       ! update intermidiate variable
       do k = KS, KE
          rhoq(k,I_QV) = rhoq(k,I_QV) + RHOQ0_t(k,I_QV)*dt
          rhoq(k,I_QC) = max(0.0_RP, rhoq(k,I_QC) + RHOQ0_t(k,I_QC)*dt )
          rhoq(k,I_QI) = max(0.0_RP, rhoq(k,I_QI) + RHOQ0_t(k,I_QI)*dt )
          rhoq(k,I_NC) = max(0.0_RP, rhoq(k,I_NC) + RHOQ0_t(k,I_NC)*dt )
          rhoq(k,I_NI) = max(0.0_RP, rhoq(k,I_NI) + RHOQ0_t(k,I_NI)*dt )

          ! cloud number concentration filter
          rhoq(k,I_NC) = min( rhoq(k,I_NC), nc_uplim_d(1,i,j) )
       end do

       do k = KS, KE
          rhoe(k) = rhoe(k) + RHOE0_t (k)*dt
          cva (k) = cva (k) + CVtot0_t(k)*dt
          cpa (k) = cpa (k) + CPtot0_t(k)*dt

          temp (k) = rhoe(k) / ( DENS(k,i,j) * cva(k) )
          pres (k) = DENS(k,i,j) * (cpa(k)-cva(k)) * temp(k)
          wtemp(k) = max( temp(k), tem_min )
       enddo


#ifdef DEBUG
       call debug_tem( KA, KS, KE, &
            2, i, j, temp(:), DENS(:,i,j), pres(:), QTRC(:,i,j,I_QV) )
#endif


       !----------------------------------------------------------------------------
       !
       ! 2.Phase change: Freezing, Melting, Vapor deposition
       !
       !----------------------------------------------------------------------------
!OCL LOOP_FISSION_TARGET(LS)
       do k = KS, KE
          rhoq2(k,I_QR)     = rhoq(k,I_QR)
          rhoq2(k,I_NR)     = rhoq(k,I_NR)
          xq(k,I_mp_QR)     = max(xr_min,  min(xr_max, rhoq2(k,I_QR)/(rhoq2(k,I_NR)+nr_min) ))

          dq_xa(k,I_mp_QR)  = a_m(I_mp_QR)*xq(k,I_mp_QR)**b_m(I_mp_QR)
          vt_xa(k,I_mp_QR,1) = alpha_v(I_mp_QR,1)*(xq(k,I_mp_QR)**beta_v(I_mp_QR,1))*rho_fac_q(k,I_mp_QR)
          vt_xa(k,I_mp_QR,2) = vt_xa(k,I_mp_QR,1)

          !! Following values shoud be already filtered to be non-zero before sbroutine was called.
          ! Mass concentration [kg/m3]
          rhoq2(k,I_QV) = rhoq(k,I_QV)
          rhoq2(k,I_QC) = rhoq(k,I_QC)
          rhoq2(k,I_QI) = rhoq(k,I_QI)
          rhoq2(k,I_QS) = rhoq(k,I_QS)
          rhoq2(k,I_QG) = rhoq(k,I_QG)
          ! Number concentration[/m3] (should be filtered to avoid zero division.)
          rhoq2(k,I_NC) = rhoq(k,I_NC)
          rhoq2(k,I_NI) = rhoq(k,I_NI)
          rhoq2(k,I_NS) = rhoq(k,I_NS)
          rhoq2(k,I_NG) = rhoq(k,I_NG)

          ! Mass of mean particle [kg]
          ! SB06(94)
          !
          xq(k,I_mp_QC)     = min(xc_max, max(xc_min, rhoq2(k,I_QC)/(rhoq2(k,I_NC)+nc_min) ))
          xq(k,I_mp_QI)     = min(xi_max, max(xi_min, rhoq2(k,I_QI)/(rhoq2(k,I_NI)+ni_min) ))
          xq(k,I_mp_QS)     = min(xs_max, max(xs_min, rhoq2(k,I_QS)/(rhoq2(k,I_NS)+ns_min) ))
          xq(k,I_mp_QG)     = min(xg_max, max(xg_min, rhoq2(k,I_QG)/(rhoq2(k,I_NG)+ng_min) ))
          ! diamter of average mass
          ! SB06(32)
          dq_xa(k,I_mp_QC)  = a_m(I_mp_QC)*xq(k,I_mp_QC)**b_m(I_mp_QC)
          dq_xa(k,I_mp_QI)  = a_m(I_mp_QI)*xq(k,I_mp_QI)**b_m(I_mp_QI)
          dq_xa(k,I_mp_QS)  = a_m(I_mp_QS)*xq(k,I_mp_QS)**b_m(I_mp_QS)
          dq_xa(k,I_mp_QG)  = a_m(I_mp_QG)*xq(k,I_mp_QG)**b_m(I_mp_QG)

          ! terminal velocity of average mass
          vt_xa(k,I_mp_QC,1) = alpha_v(I_mp_QC,1)*(xq(k,I_mp_QC)**beta_v(I_mp_QC,1))*rho_fac_q(k,I_mp_QC)
          vt_xa(k,I_mp_QI,1) = alpha_v(I_mp_QI,1)*(xq(k,I_mp_QI)**beta_v(I_mp_QI,1))*rho_fac_q(k,I_mp_QI)
          vt_xa(k,I_mp_QS,1) = alpha_v(I_mp_QS,1)*(xq(k,I_mp_QS)**beta_v(I_mp_QS,1))*rho_fac_q(k,I_mp_QS)
          vt_xa(k,I_mp_QG,1) = alpha_v(I_mp_QG,1)*(xq(k,I_mp_QG)**beta_v(I_mp_QG,1))*rho_fac_q(k,I_mp_QG)
          vt_xa(k,I_mp_QC,2) = alpha_v(I_mp_QC,2)*(xq(k,I_mp_QC)**beta_v(I_mp_QC,2))*rho_fac_q(k,I_mp_QC)
          vt_xa(k,I_mp_QI,2) = alpha_v(I_mp_QI,2)*(xq(k,I_mp_QI)**beta_v(I_mp_QI,2))*rho_fac_q(k,I_mp_QI)
          vt_xa(k,I_mp_QS,2) = alpha_v(I_mp_QS,2)*(xq(k,I_mp_QS)**beta_v(I_mp_QS,2))*rho_fac_q(k,I_mp_QS)
          vt_xa(k,I_mp_QG,2) = alpha_v(I_mp_QG,2)*(xq(k,I_mp_QG)**beta_v(I_mp_QG,2))*rho_fac_q(k,I_mp_QG)
       end do

       rhoq_crg(:,:) = 0.0_RP
       rhoq2_crg(:,:) = 0.0_RP
       if (flg_lt_l) then
          do k = KS, KE
             do iq = I_QC, I_QG
                rhoq_crg(k,iq) = DENS(k,i,j) * QTRC_crg(k,i,j,iq-1)
                rhoq2_crg(k,iq) = rhoq_crg(k,iq)
             enddo
          enddo
       end if

       call moist_psat_liq( KA, KS, KE, &
                            wtemp(:), esw(:) ) ! [IN], [OUT]
       call moist_psat_ice( KA, KS, KE, &
                            wtemp(:), esi(:) ) ! [IN], [OUT]

       call freezing_water( &
            KA, KS, KE, &
            dt,                           & ! (in)
            rhoq2(:,:), xq(:,:), temp(:), & ! (in)
            PQ(:,:)                       ) ! (inout)

       call dep_vapor_melt_ice( &
            KA, KS, KE, &
            DENS(:,i,j), wtemp(:), pres(:), qdry(:,i,j), rhoq2(:,:), & ! (in)
            esw(:), esi(:), xq(:,:), vt_xa(:,:,:), dq_xa(:,:),       & ! (in)
            PQ(:,:)                                                  ) ! (inout)

       !
       ! update subroutine
       !
       if( flg_lt_l ) then
          call update_by_phase_change( &
               KA, KS, KE, &
               ntmax_phase_change, dt,          & ! (in)
               cz(:,i,j), fz(:,i,j),            & ! (in)
               w(:,i,j),                        & ! (in)
               dTdt_equiv_d(:),                 & ! (in)
               DENS(:,i,j), qdry(:,i,j),        & ! (in)
               esw(:), esi(:),                  & ! (in)
               rhoq2(:,:), pres(:), temp(:),    & ! (in)
               cpa(:), cva(:),                  & ! (in)
               flg_lt_l,                        & ! in
               PQ(:,:),                         & ! (inout)
               sl_PLCdep, sl_PLRdep, sl_PNRdep, & ! (inout)
               RHOQ0_t(:,:), RHOE0_t(:),        & ! (out)
               CPtot0_t(:), CVtot0_t(:),        & ! (out)
               EVAPORATE(:,i,j),                & ! (out)
               rhoq2_crg(:,I_QC:I_QG),          & ! (in:optional)
               RHOQcrg0_t(:,i,j,:)              ) ! (inout:optional)
       else
          call update_by_phase_change( &
               KA, KS, KE, &
               ntmax_phase_change, dt,          & ! (in)
               cz(:,i,j), fz(:,i,j),            & ! (in)
               w(:,i,j),                        & ! (in)
               dTdt_equiv_d(:),                 & ! (in)
               DENS(:,i,j), qdry(:,i,j),        & ! (in)
               esw(:), esi(:),                  & ! (in)
               rhoq2(:,:), pres(:), temp(:),    & ! (in)
               cpa(:), cva(:),                  & ! (in)
               flg_lt_l,                        & ! in
               PQ(:,:),                         & ! (inout)
               sl_PLCdep, sl_PLRdep, sl_PNRdep, & ! (inout)
               RHOQ0_t(:,:), RHOE0_t(:),        & ! (out)
               CPtot0_t(:), CVtot0_t(:),        & ! (out)
               EVAPORATE(:,i,j)                 ) ! (out)
       endif

       ! total tendency
       do k = KS, KE
          RHOE_t (k,i,j) = RHOE_t (k,i,j) + RHOE0_t (k)
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CVtot0_t(k)
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CPtot0_t(k)
       enddo

       ! update intermidiate variable
       do iq = 1, QA_MP
       do k = KS, KE
          rhoq(k,iq) = max(0.0_RP, rhoq(k,iq) + RHOQ0_t(k,iq)*dt )
       enddo
       enddo

       do k = KS, KE
          rhoe(k) = rhoe(k) + RHOE0_t (k)*dt
          cva (k) = cva (k) + CVtot0_t(k)*dt
          cpa (k) = cpa (k) + CPtot0_t(k)*dt
          temp(k) = rhoe(k) / ( DENS(k,i,j) * cva(k) )
          pres(k) = DENS(k,i,j) * ( cpa(k) - cva(k) ) * temp(k)
       enddo

       if (flg_lt_l) then
          do iq = I_QC, I_QG
          do k = KS, KE
             rhoq_crg(k,iq) = rhoq_crg(k,iq) + RHOQcrg0_t(k,i,j,iq-1)*dt  ! need limiter?
          enddo
          enddo
       end if

#ifdef DEBUG
       call debug_tem( KA, KS, KE, &
            3, i, j, temp(:), DENS(:,i,j), pres(:), QTRC(:,i,j,I_QV) )
#endif

       !---------------------------------------------------------------------------
       !
       ! 3.Collection process
       !
       !---------------------------------------------------------------------------

       ! parameter setting

       ! Mass concentration [kg/m3]
       do k = KS, KE
          rhoq2(k,I_QV) = rhoq(k,I_QV)
          rhoq2(k,I_QC) = rhoq(k,I_QC)
          rhoq2(k,I_QR) = rhoq(k,I_QR)
          rhoq2(k,I_QI) = rhoq(k,I_QI)
          rhoq2(k,I_QS) = rhoq(k,I_QS)
          rhoq2(k,I_QG) = rhoq(k,I_QG)
       end do
       ! Number concentration[/m3]
       do k = KS, KE
          rhoq2(k,I_NC) = rhoq(k,I_NC)
          rhoq2(k,I_NR) = rhoq(k,I_NR)
          rhoq2(k,I_NI) = rhoq(k,I_NI)
          rhoq2(k,I_NS) = rhoq(k,I_NS)
          rhoq2(k,I_NG) = rhoq(k,I_NG)
       end do

       ! Mass of mean particle [kg]
       do k = KS, KE
          xq(k,I_mp_QC) = min(xc_max, max(xc_min, rhoq2(k,I_QC)/(rhoq2(k,I_NC)+nc_min) ) )
          xq(k,I_mp_QR) = min(xr_max, max(xr_min, rhoq2(k,I_QR)/(rhoq2(k,I_NR)+nr_min) ) )
          xq(k,I_mp_QI) = min(xi_max, max(xi_min, rhoq2(k,I_QI)/(rhoq2(k,I_NI)+ni_min) ) )
          xq(k,I_mp_QS) = min(xs_max, max(xs_min, rhoq2(k,I_QS)/(rhoq2(k,I_NS)+ns_min) ) )
          xq(k,I_mp_QG) = min(xg_max, max(xg_min, rhoq2(k,I_QG)/(rhoq2(k,I_NG)+ng_min) ) )
       end do

       ! effective cross section is assume as area equivalent circle
       do k = KS, KE
          dq_xa(k,I_mp_QC) = 2.0_RP*a_rea(I_mp_QC)*xq(k,I_mp_QC)**b_rea(I_mp_QC)
          dq_xa(k,I_mp_QR) = 2.0_RP*a_rea(I_mp_QR)*xq(k,I_mp_QR)**b_rea(I_mp_QR)
          dq_xa(k,I_mp_QI) = 2.0_RP*a_rea(I_mp_QI)*xq(k,I_mp_QI)**b_rea(I_mp_QI)
          dq_xa(k,I_mp_QS) = 2.0_RP*a_rea(I_mp_QS)*xq(k,I_mp_QS)**b_rea(I_mp_QS)
          dq_xa(k,I_mp_QG) = 2.0_RP*a_rea(I_mp_QG)*xq(k,I_mp_QG)**b_rea(I_mp_QG)
       end do

       ! terminal velocity of average mass
       ! SB06(33)
       do k = KS, KE
          vt_xa(k,I_mp_QC,2) = alpha_v(I_mp_QC,2)*(xq(k,I_mp_QC)**beta_v(I_mp_QC,2))*rho_fac_q(k,I_mp_QC)
          vt_xa(k,I_mp_QR,2) = alpha_v(I_mp_QR,2)*(xq(k,I_mp_QR)**beta_v(I_mp_QR,2))*rho_fac_q(k,I_mp_QR)
          vt_xa(k,I_mp_QI,2) = alpha_v(I_mp_QI,2)*(xq(k,I_mp_QI)**beta_v(I_mp_QI,2))*rho_fac_q(k,I_mp_QI)
          vt_xa(k,I_mp_QS,2) = alpha_v(I_mp_QS,2)*(xq(k,I_mp_QS)**beta_v(I_mp_QS,2))*rho_fac_q(k,I_mp_QS)
          vt_xa(k,I_mp_QG,2) = alpha_v(I_mp_QG,2)*(xq(k,I_mp_QG)**beta_v(I_mp_QG,2))*rho_fac_q(k,I_mp_QG)
       enddo

       Pcrg1(:,:) = 0.0_RP
       Pcrg2(:,:) = 0.0_RP
       if (flg_lt_l) then
          do iq = I_QC, I_QG
          do k = KS, KE
             rhoq2_crg(k,iq) = rhoq_crg(k,iq)
          enddo
          enddo
       end if

       ! Auto-conversion, Accretion, Self-collection, Break-up
       if ( MP_doautoconversion ) then
          call aut_acc_slc_brk(        &
               KA, KS, KE,             &
               flg_lt_l,               & ! (in)
               rhoq2(:,:),             & ! (in)
               rhoq2_crg(:,I_QC:I_QG), & ! (in)
               xq(:,:), dq_xa(:,:),    & ! (in)
               DENS(:,i,j),            & ! (in)
               PQ(:,:),                & ! (in)
               Pcrg1(:,:)              ) ! (inout)

       else
!OCL XFILL
          do k = KS, KE
             PQ(k,I_LCaut) = 0.0_RP
             PQ(k,I_NCaut) = 0.0_RP
             PQ(k,I_NRaut) = 0.0_RP
             PQ(k,I_LCacc) = 0.0_RP
             PQ(k,I_NCacc) = 0.0_RP
             PQ(k,I_NRslc) = 0.0_RP
             PQ(k,I_NRbrk) = 0.0_RP
             !--- for lightning
             Pcrg1(k,I_LCaut) = 0.0_RP
             Pcrg1(k,I_NCaut) = 0.0_RP
             Pcrg1(k,I_NRaut) = 0.0_RP
             Pcrg1(k,I_LCacc) = 0.0_RP
             Pcrg1(k,I_NCacc) = 0.0_RP
             Pcrg1(k,I_NRslc) = 0.0_RP
             Pcrg1(k,I_NRbrk) = 0.0_RP
          end do
       endif

       call mixed_phase_collection(            &
      ! collection process
            KA, KS, KE,                        & ! (in)
            flg_lt_l,                          & ! (in)
            d0_crg_l, v0_crg_l,                & ! (in)
            beta_crg(:,i,j),                   & ! (in)
            dqcrg(:,i,j),                      & ! (in)
            temp(:), rhoq2(:,:),               & ! (in)
            rhoq2_crg(:,I_QC:I_QG),            & ! (in)
            xq(:,:), dq_xa(:,:), vt_xa(:,:,:), & ! (in)
            PQ(:,:),                           & ! (inout)
            Pcrg1(:,:),                        & ! (inout)
            Pcrg2(:,:),                        & ! (inout)
            Pac(:,:)                           ) ! (out)

       call ice_multiplication(     &
            KA, KS, KE,             & ! (in)
            flg_lt_l,               & ! (in)
            Pac(:,:),               & ! (in)
            temp(:), rhoq2(:,:),    & ! (in)
            rhoq2_crg(:,I_QC:I_QG), & ! (in)
            xq(:,:),                & ! (in)
            PQ(:,:),                & ! (inout)
            Pcrg1(:,:)              ) ! (inout)

       !
       ! update
       !
!OCL LOOP_FISSION_TARGET(LS)
       do k = KS, KE
          ! warm collection process
          wrm_dqc = max( dt*( PQ(k,I_LCaut)+PQ(k,I_LCacc)               ), -rhoq2(k,I_QC) )
          wrm_dnc = max( dt*( PQ(k,I_NCaut)+PQ(k,I_NCacc)               ), -rhoq2(k,I_NC) )
          wrm_dnr = max( dt*( PQ(k,I_NRaut)+PQ(k,I_NRslc)+PQ(k,I_NRbrk) ), -rhoq2(k,I_NR) )
          wrm_dqr = -wrm_dqc
          ! for charge density
          if (flg_lt_l) then
             wrm_dnc_crg = dt*( Pcrg1(k,I_NCaut)+Pcrg1(k,I_NCacc) )   ! C + C -> R
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QC)),abs(wrm_dnc_crg) )
             wrm_dnc_crg = sign( sw1,wrm_dnc_crg )
             wrm_dnr_crg = - wrm_dnc_crg
          end if
          ! mixed phase collection
          ! Pxxacyy2zz xx and yy decrease and zz increase .
          !
          ! At first fixer is applied to decreasing particles.
          ! order of fixer: graupel-cloud, snow-cloud, ice-cloud, graupel-rain, snow-rain, ice-rain,
          !                 snow-ice,  ice-ice, graupel-snow, snow-snow
          ! cloud mass decrease
          gc_dqc  = max( dt*Pac(k,I_LGacLC2LG), min(0.0_RP, -rhoq2(k,I_QC)-wrm_dqc               )) ! => dqg
          sc_dqc  = max( dt*Pac(k,I_LSacLC2LS), min(0.0_RP, -rhoq2(k,I_QC)-wrm_dqc-gc_dqc        )) ! => dqs
          ic_dqc  = max( dt*Pac(k,I_LIacLC2LI), min(0.0_RP, -rhoq2(k,I_QC)-wrm_dqc-gc_dqc-sc_dqc )) ! => dqi
          ! cloud num. decrease
          gc_dnc  = max( dt*Pac(k,I_NGacNC2NG), min(0.0_RP, -rhoq2(k,I_NC)-wrm_dnc               )) ! => dnc:minus
          sc_dnc  = max( dt*Pac(k,I_NSacNC2NS), min(0.0_RP, -rhoq2(k,I_NC)-wrm_dnc-gc_dnc        )) ! => dnc:minus
          ic_dnc  = max( dt*Pac(k,I_NIacNC2NI), min(0.0_RP, -rhoq2(k,I_NC)-wrm_dnc-gc_dnc-sc_dnc )) ! => dnc:minus
          ! Decrease of absolute value of cloud charge density
          if (flg_lt_l) then
             gc_dnc_crg  = dt*Pcrg2(k,I_NGacNC2NG)  ! C + G -> G  ( move from c to g )
             sc_dnc_crg  = dt*Pcrg2(k,I_NSacNC2NS)  ! C + S -> S  ( move from c to s )
             ic_dnc_crg  = dt*Pcrg2(k,I_NIacNC2NI)  ! C + I -> I  ( move from c to i )
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QC)+wrm_dnc_crg                      ),abs(gc_dnc_crg) )
             gc_dnc_crg = sign( sw1,gc_dnc_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QC)+wrm_dnc_crg+gc_dnc_crg           ),abs(sc_dnc_crg) )
             sc_dnc_crg = sign( sw1,sc_dnc_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QC)+wrm_dnc_crg+gc_dnc_crg+sc_dnc_crg),abs(ic_dnc_crg) )
             ic_dnc_crg = sign( sw1,ic_dnc_crg )
          end if

          ! rain mass decrease ( tem < 273.15K)
          sw = sign(0.5_RP, T00-temp(k)) + 0.5_RP ! if( temp(k,i,j) <= T00 )then sw=1, else sw=0
          rg_dqr  = max( dt*Pac(k,I_LRacLG2LG  ), min(0.0_RP, -rhoq2(k,I_QR)-wrm_dqr               )) * sw
          rg_dqg  = max( dt*Pac(k,I_LRacLG2LG  ), min(0.0_RP, -rhoq2(k,I_QG)                       )) * ( 1.0_RP - sw )
          rs_dqr  = max( dt*Pac(k,I_LRacLS2LG_R), min(0.0_RP, -rhoq2(k,I_QR)-wrm_dqr-rg_dqr        )) * sw
          ri_dqr  = max( dt*Pac(k,I_LRacLI2LG_R), min(0.0_RP, -rhoq2(k,I_QR)-wrm_dqr-rg_dqr-rs_dqr )) * sw
          ! rain num. decrease
          rg_dnr  = max( dt*Pac(k,I_NRacNG2NG  ), min(0.0_RP, -rhoq2(k,I_NR)-wrm_dnr               )) * sw
          rg_dng  = max( dt*Pac(k,I_NRacNG2NG  ), min(0.0_RP, -rhoq2(k,I_NG)                       )) * ( 1.0_RP - sw )
          rs_dnr  = max( dt*Pac(k,I_NRacNS2NG_R), min(0.0_RP, -rhoq2(k,I_NR)-wrm_dnr-rg_dnr        )) * sw
          ri_dnr  = max( dt*Pac(k,I_NRacNI2NG_R), min(0.0_RP, -rhoq2(k,I_NR)-wrm_dnr-rg_dnr-rs_dnr )) * sw
          ! Decrease of absolute value of rain charge density
          if (flg_lt_l) then
             rg_dnr_crg  = dt*Pcrg2(k,I_NRacNG2NG  )* sw                 ! R + G -> G
             rg_dng_crg  = dt*Pcrg2(k,I_NRacNG2NG  )* ( 1.0_RP - sw )    ! R + G -> R
             rs_dnr_crg  = dt*Pcrg2(k,I_NRacNS2NG_R)* sw                 ! R + S -> G
             ri_dnr_crg  = dt*Pcrg2(k,I_NRacNI2NG_R)* sw                 ! R + I -> G
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QR)+wrm_dnr_crg),abs(rg_dnr_crg) )
             rg_dnr_crg = sign( sw1,rg_dnr_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QG)            ),abs(rg_dng_crg) )
             rg_dng_crg = sign( sw1,rg_dng_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QR)+wrm_dnr_crg),abs(rs_dnr_crg) )
             rs_dnr_crg = sign( sw1,rs_dnr_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QR)+wrm_dnr_crg),abs(ri_dnr_crg) )
             ri_dnr_crg = sign( sw1,ri_dnr_crg )
          end if
          ! ice mass decrease
          fac1    = (ri_dqr-eps)/ (dt*Pac(k,I_LRacLI2LG_R)-eps) ! suppress factor by filter of rain
          ri_dqi  = max( dt*Pac(k,I_LRacLI2LG_I)*fac1, min(0.0_RP, -rhoq2(k,I_QI)+ic_dqc               )) ! => dqg
          ii_dqi  = max( dt*Pac(k,I_LIacLI2LS  )     , min(0.0_RP, -rhoq2(k,I_QI)+ic_dqc-ri_dqi        )) ! => dqs
          is_dqi  = max( dt*Pac(k,I_LIacLS2LS  )     , min(0.0_RP, -rhoq2(k,I_QI)+ic_dqc-ri_dqi-ii_dqi )) ! => dqs
!!         !-- Y.Sato added(2018/8/31)
!!         gi_dqi  = max( dt*Pac(I_LGacLI2LG,k)       , min(0.0_RP, -rhoq2(I_QI,k)+ic_dqc-ri_dqi-ii_dqi-is_dqi )) ! => dqg
          ! ice num. decrease
          fac4    = (ri_dnr-eps)/ (dt*Pac(k,I_NRacNI2NG_R)-eps) ! suppress factor by filter of rain
          ri_dni  = max( dt*Pac(k,I_NRacNI2NG_I)*fac4, min(0.0_RP, -rhoq2(k,I_NI)               )) ! => dni:minus
          ii_dni  = max( dt*Pac(k,I_NIacNI2NS  )     , min(0.0_RP, -rhoq2(k,I_NI)-ri_dni        )) ! => dni:minus,dns:plus(*0.5)
          is_dni  = max( dt*Pac(k,I_NIacNS2NS  )     , min(0.0_RP, -rhoq2(k,I_NI)-ri_dni-ii_dni )) ! => dni:minus,dns:plus
!!         !-- Y.Sato added(2018/8/31)
!!         gi_dni  = max( dt*Pac(I_NGacNI2NG,k)       , min(0.0_RP, -rhoq2(I_NI,k)-ri_dni-ii_dni-is_dni )) ! => dns:minus
          ! Decrease of absolute value of ice charge density
          if (flg_lt_l) then
             ri_dni_crg  = dt*Pcrg2(k,I_NRacNI2NG_I)*fac4   !  I + R -> G
             ii_dni_crg  = dt*Pcrg2(k,I_NIacNI2NS)          !  I + I -> S
             is_dni_crg  = dt*Pcrg2(k,I_NIacNS2NS)          !  I + S -> S
!!            !-- Y.Sato added(2018/8/31)
!!            gi_dni_crg  = dt*Pcrg2(k,I_NGacNI2NG)          ! G + S -> G
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QI)-ic_dnc_crg)                      ,abs(ri_dni_crg) )
             ri_dni_crg = sign( sw1,ri_dni_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QI)-ic_dnc_crg+ri_dni_crg)           ,abs(ii_dni_crg) )
             ii_dni_crg = sign( sw1,ii_dni_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QI)-ic_dnc_crg+ri_dni_crg+ii_dni_crg),abs(is_dni_crg) )
             is_dni_crg = sign( sw1,is_dni_crg )
!!            !-- Y.Sato added(2018/8/31)
!!            sw1 = min( abs(rhoq2_crg(k,I_QI)-ic_dnc_crg+ri_dni_crg+ii_dni_crg+is_dni_crg),abs(gi_dni_crg) )
!!            gi_dni_crg = sign( sw1,gi_dni_crg )
          end if
          ! snow mass decrease
          fac3    = (rs_dqr-eps)/(dt*Pac(k,I_LRacLS2LG_R)-eps) ! suppress factor by filter of rain
          rs_dqs  = max( dt*Pac(k,I_LRacLS2LG_S)*fac3, min(0.0_RP, -rhoq2(k,I_QS)+sc_dqc+ii_dqi+is_dqi        )) ! => dqg
          gs_dqs  = max( dt*Pac(k,I_LGacLS2LG  )     , min(0.0_RP, -rhoq2(k,I_QS)+sc_dqc+ii_dqi+is_dqi-rs_dqs )) ! => dqg
          ! snow num. decrease
          fac6    = (rs_dnr-eps)/(dt*Pac(k,I_NRacNS2NG_R)-eps) ! suppress factor by filter of rain
          !       fac7    = (is_dni-eps)/(dt*Pac(I_NIacNS2NS,  k,i,j)-eps) ! suppress factor by filter of ice
          rs_dns  = max( dt*Pac(k,I_NRacNS2NG_S)*fac6, min(0.0_RP, -rhoq2(k,I_NS)+0.50_RP*ii_dni+is_dni       )) ! => dns:minus
          gs_dns  = max( dt*Pac(k,I_NGacNS2NG  )     , min(0.0_RP, -rhoq2(k,I_NS)+0.50_RP*ii_dni+is_dni-rs_dns )) ! => dns:minus
          ss_dns  = max( dt*Pac(k,I_NSacNS2NS  )     , min(0.0_RP, -rhoq2(k,I_NS)+0.50_RP*ii_dni+is_dni-rs_dns-gs_dns ))
          !
          gg_dng  = max( dt*Pac(k,I_NGacNG2NG)       , min(0.0_RP, -rhoq2(k,I_NG) ))
          if (flg_lt_l) then
             ! Decrease of absolute value of snow charge density
             rs_dns_crg  = dt*Pcrg2(k,I_NRacNS2NG_S)*fac6  ! R + S -> G
             gs_dns_crg  = dt*Pcrg2(k,I_NGacNS2NG)         ! G + S -> G
             ss_dns_crg  = 0.0_RP                              ! S + S -> S (No charge transfer)
             gg_dng_crg  = 0.0_RP                              ! G + G -> G (No charge transfer)
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QS)-sc_dnc_crg-ii_dni-is_dni_crg),           abs(rs_dns_crg) )
             rs_dns_crg = sign( sw1,rs_dns_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QS)-sc_dnc_crg-ii_dni-is_dni_crg+rs_dns_crg),abs(gs_dns_crg) )
             gs_dns_crg = sign( sw1,gs_dns_crg )
             !--- Charge split
             sw1 = sign(0.5_RP, abs( Pcrg2(k,I_CGNGacNS2NG) )-EPS ) + 0.5_RP ! if abs Pcrg2 is smaller than EPS, sw=1, else sw=0
             sw2 = sign(0.5_RP, abs( Pcrg2(k,I_CGNGacNI2NG) )-EPS ) + 0.5_RP ! if abs Pcrg2 is smaller than EPS, sw=1, else sw=0
             crg_split_g =  dt*Pcrg2(k,I_CGNGacNS2NG)*sw1 &
                         +  dt*Pcrg2(k,I_CGNGacNI2NG)*sw2
             crg_split_s = -dt*Pcrg2(k,I_CGNGacNS2NG)*sw1
             crg_split_i = 0.0_RP
!!            crg_split_i = -dt*Pcrg2(k,I_CGNGacNI2NG)*sw2  !Y.Sato (2018/8/31)
             QSPLT_in(k,i,j,1) = crg_split_g / dt    ! fC/s
             QSPLT_in(k,i,j,3) = crg_split_s / dt    ! fC/s
             QSPLT_in(k,i,j,2) = crg_split_i / dt    ! fC/s
          end if
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
          ! Decrease of absolute value of graupel charge density
          if (flg_lt_l) then
             clp_dnc_crg = 0.0_RP
             clp_dnr_crg = -rg_dng_crg*(1.0_RP-sw)
             clp_dni_crg = -ic_dnc_crg !&
!!                          +crg_split_i ! Y.Sato added (2018/8/31)
             clp_dns_crg = -sc_dnc_crg-ii_dni_crg-is_dni_crg-ss_dns_crg &
                           +crg_split_s
             clp_dng_crg = -gc_dnc_crg+(-rg_dnr_crg-rs_dnr_crg-ri_dnr_crg)*sw &
                           -ri_dni_crg-rs_dns_crg-gs_dns_crg-gg_dng_crg &
!!                          -gi_dni_crg & !--- Y.Sato (2018/8/31)
                           +crg_split_g
          end if

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
          ! charge density
          if (flg_lt_l) then
             clm_dnc_crg = gc_dnc_crg+sc_dnc_crg+ic_dnc_crg
             clm_dnr_crg = (rg_dnr_crg+rs_dnr_crg+ri_dnr_crg) * sw
             clm_dni_crg = ri_dni_crg+ii_dni_crg+is_dni_crg  !!&
!!                        + gi_dni_crg  !Y.Sato (2018/8/31)
             clm_dns_crg = rs_dns_crg+gs_dns_crg+ss_dns_crg
             clm_dng_crg = gg_dng_crg+rg_dng_crg*(1.0_RP-sw)
          end if

          ! partial conversion
          ! 08/05/08 [Mod] T.Mitsui
          pco_dqi = max( dt*PQ(k,I_LIcon), -clp_dqi )
          pco_dqs = max( dt*PQ(k,I_LScon), -clp_dqs )
          pco_dqg = -pco_dqi-pco_dqs
          ! 08/05/08 [Mod] T.Mitsui
          pco_dni = max( dt*PQ(k,I_NIcon), -clp_dni )
          pco_dns = max( dt*PQ(k,I_NScon), -clp_dns )
          pco_dng = -pco_dni-pco_dns
          !-- for charge density
          if (flg_lt_l) then
             pco_dni_crg = dt*Pcrg1(k,I_NIcon)
             pco_dns_crg = dt*Pcrg1(k,I_NScon)
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QI)+clp_dni_crg  ),abs(pco_dni_crg) )
             pco_dni_crg = sign( sw1,pco_dni_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QS)+clp_dns_crg  ),abs(pco_dns_crg) )
             pco_dns_crg = sign( sw1,pco_dns_crg )
             pco_dng_crg = -pco_dni_crg-pco_dns_crg
          end if
          ! enhanced melting ( always negative value )
          ! ice-cloud melting produces cloud, others produce rain
          eml_dqi =  max( dt*PQ(k,I_LIacm), min(0.0_RP, -rhoq2(k,I_QI)-(clp_dqi+clm_dqi)-pco_dqi ))
          eml_dqs =  max( dt*PQ(k,I_LSacm), min(0.0_RP, -rhoq2(k,I_QS)-(clp_dqs+clm_dqs)-pco_dqs ))
          eml_dqg =  max( dt*(PQ(k,I_LGacm)+PQ(k,I_LGarm)+PQ(k,I_LSarm)+PQ(k,I_LIarm)), &
                     min(0.0_RP, -rhoq2(k,I_QG)-(clp_dqg+clm_dqg)-pco_dqg ))
          eml_dqc = -eml_dqi
          eml_dqr = -eml_dqs-eml_dqg
          !
          eml_dni =  max( dt*PQ(k,I_NIacm), min(0.0_RP, -rhoq2(k,I_NI)-(clp_dni+clm_dni)-pco_dni ))
          eml_dns =  max( dt*PQ(k,I_NSacm), min(0.0_RP, -rhoq2(k,I_NS)-(clp_dns+clm_dns)-pco_dns ))
          eml_dng =  max( dt*(PQ(k,I_NGacm)+PQ(k,I_NGarm)+PQ(k,I_NSarm)+PQ(k,I_NIarm)), &
                     min(0.0_RP, -rhoq2(k,I_NG)-(clp_dng+clm_dng)-pco_dng ))
          eml_dnc = -eml_dni
          eml_dnr = -eml_dns-eml_dng
          !-- for charge density
          if (flg_lt_l) then
             eml_dni_crg =  dt*Pcrg1(k,I_NIacm)   ! I+C->C
             eml_dns_crg =  dt*Pcrg1(k,I_NSacm)   ! S+C->R
             eml_dng_crg =  dt*(Pcrg1(k,I_NGacm)+Pcrg1(k,I_NGarm)+Pcrg1(k,I_NSarm)+Pcrg1(k,I_NIarm)) ! G+R->R, G+C->R
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QI)+clp_dni_crg+clm_dni_crg+pco_dni_crg  ),abs(eml_dni_crg) )
             eml_dni_crg = sign( sw1,eml_dni_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QS)+clp_dns_crg+clm_dns_crg+pco_dns_crg  ),abs(eml_dns_crg) )
             eml_dns_crg = sign( sw1,eml_dns_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QG)+clp_dng_crg+clm_dng_crg+pco_dng_crg  ),abs(eml_dng_crg) )
             eml_dng_crg = sign( sw1,eml_dng_crg )

             eml_dnc_crg = -eml_dni_crg
             eml_dnr_crg = -eml_dns_crg-eml_dng_crg
          end if
          !
          ! ice multiplication
          spl_dqg = max( dt*PQ(k,I_LGspl), min(0.0_RP, -rhoq2(k,I_QG)-(clp_dqg+clm_dqg)-pco_dqg-eml_dqg ))
          spl_dqs = max( dt*PQ(k,I_LSspl), min(0.0_RP, -rhoq2(k,I_QS)-(clp_dqs+clm_dqs)-pco_dqs-eml_dqs ))
          spl_dqi = -spl_dqg-spl_dqs
          fac9    = (spl_dqg-eps)/(dt*PQ(k,I_LGspl)-eps)
          fac10   = (spl_dqs-eps)/(dt*PQ(k,I_LSspl)-eps)
          spl_dni = dt*PQ(k,I_NIspl)*fac9*fac10
          !-- for charge density
          if (flg_lt_l) then
             spl_dns_crg = dt*Pcrg1(k,I_NSspl)*fac9*fac10
             spl_dng_crg = dt*Pcrg1(k,I_NGspl)*fac9*fac10
             !--- limiter
             sw1 = min( abs(rhoq2_crg(k,I_QS)+clp_dns_crg+pco_dns_crg+eml_dns_crg  ),abs(spl_dns_crg) )
             spl_dns_crg = sign( sw1,spl_dns_crg )
             sw1 = min( abs(rhoq2_crg(k,I_QG)+clp_dng_crg+pco_dng_crg+eml_dng_crg  ),abs(spl_dng_crg) )
             spl_dng_crg = sign( sw1,spl_dng_crg )
             spl_dni_crg = -spl_dns_crg-spl_dng_crg
          end if


          !
          ! melting and freezing limiter
          di2l = clp_dqc + clp_dqr + clm_dqc + clm_dqr + eml_dqc + eml_dqr ! = - ( clp_dqi + clp_dqs + clp_dqg + clm_dqi + clm_dqs + clm_dqg + eml_dqi + eml_dqs + eml_dqg )
          dtem = - di2l * LHF0 /  ( cva(k) * DENS(k,i,j) )
          if ( abs(dtem) < EPS ) then
             fact = 1.0_RP
          else
             fact = min( 1.0_RP, max( 0.0_RP, ( T00 - temp(k) ) / dtem ) )
          end if

          !
          ! total cloud change
          drhoqc = wrm_dqc + ( clp_dqc + clm_dqc + eml_dqc ) * fact
          drhonc = wrm_dnc + ( clp_dnc + clm_dnc + eml_dnc ) * fact
          ! total rain change
          drhoqr = wrm_dqr + ( clp_dqr + clm_dqr + eml_dqr ) * fact
          drhonr = wrm_dnr + ( clp_dnr + clm_dnr + eml_dnr ) * fact
          ! total ice change
          drhoqi =           ( clp_dqi + clm_dqi + eml_dqi ) * fact + pco_dqi + spl_dqi
          drhoni =           ( clp_dni + clm_dni + eml_dni ) * fact + pco_dni + spl_dni
          ! total snow change
          drhoqs =           ( clp_dqs + clm_dqs + eml_dqs ) * fact + pco_dqs + spl_dqs
          drhons =           ( clp_dns + clm_dns + eml_dns ) * fact + pco_dns
          ! total graupel change
          drhoqg =           ( clp_dqg + clm_dqg + eml_dqg ) * fact + pco_dqg + spl_dqg
          drhong =           ( clp_dng + clm_dng + eml_dng ) * fact + pco_dng
          !-- for charge density
          if (flg_lt_l) then
             drhoqcrg_c = wrm_dnc_crg + ( clp_dnc_crg + clm_dnc_crg + eml_dnc_crg ) * fact
             drhoqcrg_r = wrm_dnr_crg + ( clp_dnr_crg + clm_dnr_crg + eml_dnr_crg ) * fact
             drhoqcrg_i = ( clp_dni_crg + clm_dni_crg + eml_dni_crg ) * fact + pco_dni_crg + spl_dni_crg
             drhoqcrg_s = ( clp_dns_crg + clm_dns_crg + eml_dns_crg ) * fact + pco_dns_crg + spl_dns_crg
             drhoqcrg_g = ( clp_dng_crg + clm_dng_crg + eml_dng_crg ) * fact + pco_dng_crg + spl_dng_crg
          end if
          !

          ! tendency
          RHOQ0_t(k,I_QC) = drhoqc / dt
          RHOQ0_t(k,I_NC) = drhonc / dt
          RHOQ0_t(k,I_QR) = drhoqr / dt
          RHOQ0_t(k,I_NR) = drhonr / dt
          RHOQ0_t(k,I_QI) = drhoqi / dt
          RHOQ0_t(k,I_NI) = drhoni / dt
          RHOQ0_t(k,I_QS) = drhoqs / dt
          RHOQ0_t(k,I_NS) = drhons / dt
          RHOQ0_t(k,I_QG) = drhoqg / dt
          RHOQ0_t(k,I_NG) = drhong / dt

          if (flg_lt_l) then
             RHOQcrg0_t(k,i,j,I_QC-1) = drhoqcrg_c / dt
             RHOQcrg0_t(k,i,j,I_QR-1) = drhoqcrg_r / dt
             RHOQcrg0_t(k,i,j,I_QI-1) = drhoqcrg_i / dt
             RHOQcrg0_t(k,i,j,I_QS-1) = drhoqcrg_s / dt
             RHOQcrg0_t(k,i,j,I_QG-1) = drhoqcrg_g / dt
          end if

          RHOE0_t(k) = LHF * ( drhoqi + drhoqs + drhoqg ) / dt

          dql = rrho(k) * ( drhoqc + drhoqr )
          dqi = rrho(k) * ( drhoqi + drhoqs + drhoqg )

          dcv = CV_WATER * dql + CV_ICE * dqi
          dcp = CP_WATER * dql + CP_ICE * dqi

          CVtot0_t(k) = dcv / dt
          CPtot0_t(k) = dcp / dt

       enddo

       ! total tendency
       do k = KS, KE
          RHOE_t (k,i,j) = RHOE_t (k,i,j) + RHOE0_t (k)
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CVtot0_t(k)
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CPtot0_t(k)
       enddo

       !--- update
       do iq = I_QC, I_NG
       do k = KS, KE
          rhoq(k,iq) = max(0.0_RP, rhoq(k,iq) + RHOQ0_t(k,iq)*dt )
       enddo
       enddo

       !--- for lithgning component
       if (flg_lt_l) then
          do iq = I_QC, I_QG
          do k = KS, KE
             rhoq_crg(k,iq) = rhoq_crg(k,iq) + RHOQcrg0_t(k,i,j,iq-1)*dt  !-- need limiter?
          enddo
          enddo
          do iq = I_QC, I_NG
          do k = KS, KE
             QTRC0(k,iq) = rhoq(k,iq) / DENS(k,i,j)
          enddo
          enddo

          call Cross_Section( KA, KS, KE,  & ! [IN]
                              QA_MP,       & ! [IN]
                              QTRC0(:,:),  & ! [IN]
                              DENS(:,i,j), & ! [IN]
                              Crs(:,:)     ) ! [OUT]

          do k = KS, KE
             Sarea(k,i,j,I_mp_QC) = Crs(k,I_mp_QC)
             Sarea(k,i,j,I_mp_QR) = Crs(k,I_mp_QR)
             Sarea(k,i,j,I_mp_QI) = Crs(k,I_mp_QI)
             Sarea(k,i,j,I_mp_QS) = Crs(k,I_mp_QS)
             Sarea(k,i,j,I_mp_QG) = Crs(k,I_mp_QG)
          enddo
       end if

       ! total tendency
       do iq = I_QV, I_NG
       do k = KS, KE
          RHOQ_t(k,i,j,iq) = ( rhoq(k,iq) - DENS(k,i,j)*QTRC(k,i,j,iq) )/dt
       enddo
       enddo

       if (flg_lt_l) then
          do iq = I_QC, I_QG
          do k = KS, KE
             RHOQcrg_t_mp(k,i,j,iq-1) = ( rhoq_crg(k,iq) - DENS(k,i,j)*QTRC_crg(k,i,j,iq-1) )/dt
          enddo
          enddo
       end if

#ifdef DEBUG
       call debug_tem( KA, KS, KE, &
            4, i, j, temp(:), DENS(:,i,j), pres(:), QTRC(:,i,j,I_QV) )
#endif

    end do
    end do

    return
  end subroutine mp_sn14

  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine debug_tem( &
       KA, KS, KE, &
       point, i, j, &
       tem, rho, pre, qv )
    use scale_prc, only: &
       PRC_myrank
    implicit none
    integer, intent(in) :: KA, KS, KE

    integer,  intent(in) :: point
    integer,  intent(in) :: i, j
    real(RP), intent(in) :: tem(KA)
    real(RP), intent(in) :: rho(KA)
    real(RP), intent(in) :: pre(KA)
    real(RP), intent(in) :: qv (KA)

    integer :: k
    !---------------------------------------------------------------------------

    do k = KS, KE
       if (      tem(k) < tem_min &
            .OR. rho(k) < rho_min &
            .OR. pre(k) < 1.0_RP    ) then

          LOG_INFO("ATMOS_PHY_MP_SN14_debug_tem_kij",'(A,I3,A,4(F16.5),3(I6))') &
          "point: ", point, " low tem,rho,pre:", tem(k), rho(k), pre(k), qv(k), k, i, j, PRC_myrank
       endif
    enddo

    return
  end subroutine debug_tem

!OCL SERIAL
  subroutine nucleation( &
       KA, KS, KE, &
       cz, fz, w,           &
       rho, tem, pre, qdry, &
       rhoq,                &
       cpa,                 & ! in
       dTdt_rad,            & ! in
       qke,                 & ! in
       CCN, nc_uplim_d,     & ! in
       dt,                  & ! in
       PQ                   )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_saturation, only: &
       moist_psat_liq       => ATMOS_SATURATION_psat_liq, &
       moist_psat_ice       => ATMOS_SATURATION_psat_ice,   &
       moist_pres2qsat_liq  => ATMOS_SATURATION_pres2qsat_liq, &
       moist_pres2qsat_ice  => ATMOS_SATURATION_pres2qsat_ice,   &
       moist_dqsi_dtem_dens => ATMOS_SATURATION_dqs_dtem_dens_liq
    implicit none

    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in)  :: cz(  KA)
    real(RP), intent(in)  :: fz(0:KA)
    real(RP), intent(in)  :: w   (KA)    ! w of full level
    real(RP), intent(in)  :: rho (KA)
    real(RP), intent(in)  :: tem (KA)
    real(RP), intent(in)  :: pre (KA)
    real(RP), intent(in)  :: qdry(KA)
    !
    real(RP), intent(in)  :: rhoq(KA,I_QV:I_NG)
    !
    real(RP), intent(in) ::  cpa(KA)
    real(RP), intent(in) :: dTdt_rad(KA)
    real(RP), intent(in) :: qke(KA)
    real(RP), intent(in) :: dt
    real(RP), intent(in) :: CCN(KA)
    real(RP), intent(in) :: nc_uplim_d
    !
    real(RP), intent(out) :: PQ(KA,PQ_MAX)
    !
    !
!    real(RP) :: c_ccn_map(1)   ! c_ccn horizontal distribution
!    real(RP) :: kappa_map(1)   ! kappa horizontal distribution
!    real(RP) :: c_in_map(1)    ! c_in  horizontal distribution ! [Add] 11/08/30 T.Mitsui
    real(RP) :: esw(KA)       ! saturation vapor pressure, water
    real(RP) :: esi(KA)       !                            ice
    real(RP) :: ssw(KA)       ! super saturation (water)
    real(RP) :: ssi(KA)       ! super saturation (ice)
!    real(RP) :: w_dsswdz(KA)  ! w*(d_ssw/ d_z) super saturation(water) flux
    real(RP) :: w_dssidz(KA)  ! w*(d_ssi/ d_z), 09/04/14 T.Mitsui
!    real(RP) :: ssw_below(KA) ! ssw(k-1)
    real(RP) :: ssi_below(KA) ! ssi(k-1), 09/04/14 T.Mitsui
    real(RP) :: z_below(KA)   ! z(k-1)
    real(RP) :: dzh           ! z(k)-z(k-1)
    real(RP) :: pv            ! vapor pressure
    ! work variables for Twomey Equation.
    real(RP) :: qsw(KA)
    real(RP) :: qsi(KA)
    real(RP) :: dqsidtem_rho(KA)
    real(RP) :: dssidt_rad(KA)
    real(RP) :: wssi, wdssi
    !
!    real(RP) :: xi_nuc(1)    ! xi use the value @ cloud base
!    real(RP) :: alpha_nuc(1) ! alpha_nuc
!    real(RP) :: eta_nuc(1)   ! xi use the value @ cloud base
    !
    real(RP) :: sigma_w(KA)
    real(RP) :: weff(KA)
    real(RP) :: weff_max(KA)
    real(RP) :: velz(KA)
    !
    real(RP) :: coef_ccn
    real(RP) :: slope_ccn
    real(RP) :: nc_new(KA)
    real(RP) :: nc_new_below(KA)
    real(RP) :: dnc_new
    real(RP) :: nc_new_max   ! Lohmann (2002),
    real(RP) :: a_max
    real(RP) :: b_max
    logical :: flag_nucleation(KA)
    !
    real(RP) :: r_gravity
    real(RP), parameter :: r_sqrt3=0.577350269_RP ! = sqrt(1.d0/3.d0)
    real(RP), parameter :: eps=1.E-30_RP
    !====> ! 09/08/18
    !
    real(RP) :: dlcdt_max, dli_max ! defined by supersaturation
    real(RP) :: dncdt_max, dni_max ! defined by supersaturation
    real(RP) :: rdt

    real(RP) :: tmp
    !
    integer :: k
    !
    !
!    c_ccn_map(1) = c_ccn
!    kappa_map(1) = kappa
!    c_in_map(1)  = c_in
    !
    !
    rdt            = 1.0_RP/dt
    r_gravity      = 1.0_RP/GRAV
    !
    call moist_psat_liq      ( KA, KS, KE, &
                              tem(:), esw(:) ) ! [IN], [OUT]
    call moist_psat_ice      ( KA, KS, KE, &
                              tem(:), esi(:) ) ! [IN], [OUT]
    call moist_pres2qsat_liq ( KA, KS, KE, &
                              tem(:), pre(:), qdry(:), & ! [IN]
                              qsw(:)                   ) ! [OUT]
    call moist_pres2qsat_ice ( KA, KS, KE, &
                              tem(:), pre(:), qdry(:), & ! [IN]
                              qsi(:)                   ) ! [OUT]
    call moist_dqsi_dtem_dens( KA, KS, KE, &
                               tem(:), rho(:), & ! [IN]
                               dqsidtem_rho(:) ) ! [OUT]
    !
    ! Lohmann (2002),JAS, eq.(1) but changing unit [cm-3] => [m-3]
    a_max = 1.E+6_RP*0.1_RP*(1.E-6_RP)**1.27_RP
    b_max = 1.27_RP
    !
    ssi_max = 1.0_RP

    do k = KS, KE
       pv     = rhoq(k,I_QV)*Rvap*tem(k)
       ssw(k) = min( MP_ssw_lim, ( pv/esw(k)-1.0_RP ) )*100.0_RP
       ssi(k) = ( pv/esi(k) - 1.00_RP )
!       ssw_below(k+1) = ssw(k)
       ssi_below(k+1) = ssi(k)
       z_below(k+1)   = cz(k)
    end do
!    ssw_below(KS) = ssw(KS)
    ssi_below(KS) = ssi(KS)
    z_below(KS)   = cz(KS-1)

    ! dS/dz is evaluated by first order upstream difference
    !***  Solution for Twomey Equation ***
!    coef_ccn  = 1.E+6_RP*0.88_RP*(c_ccn_map(1)*1.E-6_RP)**(2.0_RP/(kappa_map(1) + 2.0_RP)) * &
    coef_ccn  = 1.E+6_RP*0.88_RP*(c_ccn*1.E-6_RP)**(2.0_RP/(kappa + 2.0_RP)) &
!           * (70.0_RP)**(kappa_map(1)/(kappa_map(1) + 2.0_RP))
            * (70.0_RP)**(kappa/(kappa + 2.0_RP))
!    slope_ccn = 1.5_RP*kappa_map(1)/(kappa_map(1) + 2.0_RP)
    slope_ccn = 1.5_RP*kappa/(kappa + 2.0_RP)
    !
    do k=KS, KE
       sigma_w(k) = r_sqrt3*sqrt(max(qke(k),qke_min))
    end do
    sigma_w(KS-1) = sigma_w(KS)
    sigma_w(KE+1) = sigma_w(KE)
    ! effective vertical velocity
    do k=KS, KE
       weff(k) = w(k) - cpa(k)*r_gravity*dTdt_rad(k)
    end do
    !
    if( MP_couple_aerosol ) then

       do k = KS, KE
          if( ssw(k) > 1.e-10_RP .AND. pre(k) > 300.E+2_RP ) then
             nc_new(k) = max( CCN(k), c_ccn )
          else
             nc_new(k) = 0.0_RP
          endif
       enddo

    else

       if( nucl_twomey ) then
          ! diagnose cloud condensation nuclei
          do k = KS, KE
             ! effective vertical velocity (maximum vertical velocity in turbulent flow)
             weff_max(k) = weff(k) + sigma_w(k)
             ! large scale upward motion region and saturated
             if( (weff(k) > 1.E-8_RP) .AND. (ssw(k) > 1.E-10_RP)  .AND. pre(k) > 300.E+2_RP )then
                ! Lohmann (2002), eq.(1)
                nc_new_max   = coef_ccn*weff_max(k)**slope_ccn
                nc_new(k) = a_max*nc_new_max**b_max
             else
                nc_new(k) = 0.0_RP
             end if
          end do

       else
          ! calculate cloud condensation nuclei
          do k = KS, KE
             if( ssw(k) > 1.e-10_RP .AND. pre(k) > 300.E+2_RP ) then
                nc_new(k) = c_ccn*ssw(k)**kappa
             else
                nc_new(k) = 0.0_RP
             endif
          enddo
      endif

    endif

    do k = KS, KE
       ! nc_new is bound by upper limit
       if( nc_new(k) > nc_uplim_d ) then ! no more CCN
          flag_nucleation(k) = .false.
          nc_new_below(k+1)  = 1.E+30_RP
       else if( nc_new(k) > eps ) then ! nucleation can occur
          flag_nucleation(k) = .true.
          nc_new_below(k+1)  = nc_new(k)
       else ! nucleation cannot occur(unsaturated or negative w)
          flag_nucleation(k) = .false.
          nc_new_below(k+1)  = 0.0_RP
       end if
    end do
    nc_new_below(KS) = 0.0_RP
!    do k=KS, KE
        ! search maximum value of nc_new
!       if(  ( nc_new(k) < nc_new_below(k) ) .OR. &
!            ( nc_new_below(k) > c_ccn_map(1)*0.05_RP ) )then ! 5% of c_ccn
!            ( nc_new_below(k) > c_ccn*0.05_RP ) )then ! 5% of c_ccn
!          flag_nucleation(k) = .false.
!       end if
!    end do

    if( MP_couple_aerosol ) then
       do k = KS, KE
          ! nucleation occurs at only cloud base.
          ! if CCN is more than below parcel, nucleation newly occurs
          ! effective vertical velocity
          if ( flag_nucleation(k) .AND. & ! large scale upward motion region and saturated
               tem(k) > tem_ccn_low ) then
             dlcdt_max     = ( rhoq(k,I_QV) - esw(k) / ( Rvap * tem(k) ) ) * rdt
             dlcdt_max     = max( dlcdt_max, 0.0_RP ) ! dlcdt_max can be artificially negative due to truncation error in floating point operation
             dncdt_max     = dlcdt_max/xc_min
!             dnc_new       = nc_new(k)-rhoq(k,I_NC)
             dnc_new       = nc_new(k)
             PQ(k,I_NCccn) = min( dncdt_max, dnc_new*rdt )
             PQ(k,I_LCccn) = min( dlcdt_max, xc_min*PQ(k,I_NCccn) )
          else
             PQ(k,I_NCccn) = 0.0_RP
             PQ(k,I_LCccn) = 0.0_RP
          end if
       end do
    else

       if( nucl_twomey ) then
          do k = KS, KE
             ! nucleation occurs at only cloud base.
             ! if CCN is more than below parcel, nucleation newly occurs
             ! effective vertical velocity
             if ( flag_nucleation(k) .AND. & ! large scale upward motion region and saturated
                  tem(k) > tem_ccn_low .AND. &
                  nc_new(k) > rhoq(k,I_NC) ) then
                dlcdt_max     = ( rhoq(k,I_QV) - esw(k) / ( Rvap * tem(k) ) ) * rdt
                dlcdt_max     = max( dlcdt_max, 0.0_RP ) ! dlcdt_max can be artificially negative due to truncation error in floating point operation
                dncdt_max     = dlcdt_max/xc_min
                dnc_new       = nc_new(k)-rhoq(k,I_NC)
                PQ(k,I_NCccn) = min( dncdt_max, dnc_new*rdt )
                PQ(k,I_LCccn) = min( dlcdt_max, xc_min*PQ(k,I_NCccn) )
             else
                PQ(k,I_NCccn) = 0.0_RP
                PQ(k,I_LCccn) = 0.0_RP
             end if
          end do
       else
          do k = KS, KE
             ! effective vertical velocity
             if(  tem(k) > tem_ccn_low .AND. &
                  nc_new(k) > rhoq(k,I_NC) ) then
                dlcdt_max     = ( rhoq(k,I_QV) - esw(k) / ( Rvap * tem(k) ) ) * rdt
                dlcdt_max     = max( dlcdt_max, 0.0_RP ) ! dlcdt_max can be artificially negative due to truncation error in floating point operation
                dncdt_max     = dlcdt_max/xc_min
                dnc_new       = nc_new(k)-rhoq(k,I_NC)
                PQ(k,I_NCccn) = min( dncdt_max, dnc_new*rdt )
                PQ(k,I_LCccn) = min( dlcdt_max, xc_min*PQ(k,I_NCccn) )
             else
                PQ(k,I_NCccn) = 0.0_RP
                PQ(k,I_LCccn) = 0.0_RP
             end if
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
    do k = KS, KE-1
       velz(k)           = ( w(k) * ( cz(k+1) - fz(k) ) + w(k+1) * ( fz(k) - cz(k) ) ) / ( cz(k+1) - cz(k) ) ! @ half level
    end do
    velz(KE) = 0.0_RP
    do k = KS, KE
       dzh           = cz(k) - z_below(k)
       w_dssidz(k)   = velz(k) * (ssi(k) - ssi_below(k))/dzh
       dssidt_rad(k) = -rhoq(k,I_QV)/(rho(k)*qsi(k)*qsi(k))*dqsidtem_rho(k)*dTdt_rad(k)
       dli_max       = ( rhoq(k,I_QV) - esi(k) / ( Rvap * tem(k) ) ) * rdt
       dli_max       = max( dli_max, 0.0_RP ) ! dli_max can be artificially negative due to truncation error in floating point operation
       dni_max       = min( dli_max/xi_ccn, (in_max-rhoq(k,I_NI))*rdt )
       wdssi         = min( w_dssidz(k)+dssidt_rad(k), 0.01_RP)
       wssi          = min( ssi(k), ssi_max)
       ! SB06(34),(35)
       if(  ( wdssi       > eps         ) .AND. & !
            (tem(k)       < 273.15_RP   ) .AND. & !
            (rhoq(k,I_NI) < in_max      ) .AND. &
            (wssi      >= eps       ) )then   !
          tmp = c_in * nm_M92 * exp( 0.3_RP * bm_M92 * ( wssi - 0.1_RP ) )
          if( inucl_w ) then
             tmp = bm_M92 * 0.3_RP * tmp * wdssi
          else
             tmp = max( tmp - rhoq(k,I_NI), 0.0_RP ) * rdt
          endif
          PQ(k,I_NIccn) = min(dni_max, tmp)
          PQ(k,I_LIccn) = min(dli_max, PQ(k,I_NIccn)*xi_ccn )
       else
          PQ(k,I_NIccn) = 0.0_RP
          PQ(k,I_LIccn) = 0.0_RP
       end if
    end do

    return
  end subroutine nucleation
  !----------------------------
!OCL SERIAL
  subroutine ice_multiplication( &
    KA, KS, KE,    & ! in
    flg_lt,        & ! in
    Pac,           & ! in
    tem, rhoq,     & ! in
    rhoq_crg, xq,  & ! in
    PQ, Pcrg1      ) ! inout

    ! ice multiplication by splintering
    ! we consider Hallet-Mossop process
    use scale_specfunc, only: &
       gammafunc => SF_gamma
    implicit none

    integer, intent(in) :: KA, KS, KE
    !
    real(RP), intent(in) :: Pac(KA,Pac_MAX)
    real(RP), intent(in) :: tem(KA)
    real(RP), intent(in) :: rhoq(KA,I_QV:I_NG)
    real(RP), intent(in) :: xq(KA,HYDRO_MAX)
    !
    real(RP), intent(inout):: PQ(KA,PQ_MAX)
    ! for lightning
    logical,  intent(in) :: flg_lt
    real(RP), intent(in) :: rhoq_crg(KA,I_QC:I_QG)
    real(RP), intent(inout):: Pcrg1(KA,PQ_MAX)
    !
    ! production of (350.d3 ice particles)/(cloud riming [g]) => 350*d6 [/kg]
    real(RP), parameter :: pice = 350.0E+6_RP
    ! production of (1 ice particle)/(250 cloud particles riming)
    real(RP), parameter :: pnc  = 250.0_RP
    ! temperature factor
    real(RP) :: fp

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
    integer :: k, iq
    real(RP) :: sw1
    !
    !
!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       ! Here we assume particle temperature is same as environment temperature.
       ! If you want to treat in a better manner,
       ! you can diagnose with eq.(64) in CT(86)
       if     (tem(k) >  270.16_RP)then
          fp   = 0.0_RP
       else if(tem(k) >= 268.16_RP)then
          fp   = (270.16_RP-tem(k))*0.5_RP
       else if(tem(k) >= 265.16_RP)then
          fp   = (tem(k)-265.16_RP)*0.333333333_RP
       else
          fp   = 0.0_RP
       end if
       ! Approximation of Incomplete Gamma function
       ! Here we calculate with algorithm by Numerical Recipes.
       ! This approach is based on eq.(78) in Cotton etal.(1986),
       ! but more accurate and expanded for Generalized Gamma Distribution.
       x       = coef_lambda(I_mp_QC)*(xc_cr/xq(k,I_mp_QC))**mu(I_mp_QC)
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
       n12 = rhoq(k,I_NC)*(1.0_RP-igm)
       ! eq.(82) CT(86)
       wn            = (pice + n12/((rhoq(k,I_QC)+xc_min)*pnc) )*fp ! filtered by xc_min
       wni           = wn*(-Pac(k,I_LIacLC2LI)  ) ! riming production rate is all negative
       wns           = wn*(-Pac(k,I_LSacLC2LS)  )
       wng           = wn*(-Pac(k,I_LGacLC2LG)  )
       PQ(k,I_NIspl) = wni+wns+wng
       !
       PQ(k,I_LSspl) = - wns*xq(k,I_mp_QI) ! snow    => ice
       PQ(k,I_LGspl) = - wng*xq(k,I_mp_QI) ! graupel => ice
       if (flg_lt) then
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NS)-SMALL ) !--- if NC is small,ignore charge transfer
          Pcrg1(k,I_NSspl) = - wns*(1.0_RP-sw1) &
                                / (rhoq(k,I_NS)+sw1)*rhoq_crg(k,I_QS)
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NG)-SMALL ) !--- if NC is small,ignore charge transfer
          Pcrg1(k,I_NGspl) = - wng*(1.0_RP-sw1) &
                                / (rhoq(k,I_NG)+sw1)*rhoq_crg(k,I_QG)
          Pcrg1(k,I_NIspl) = - ( Pcrg1(k,I_NSspl) + Pcrg1(k,I_NGspl) )
       end if
       !
    end do
    !
    return
  end subroutine ice_multiplication
  !----------------------------
!OCL SERIAL
  subroutine mixed_phase_collection( &
    ! collection process
       KA, KS, KE,            & ! in
       flg_lt,                & ! in
       d0_crg, v0_crg,        & ! in
       beta_crg, dqcrg,       & ! in
       wtem, rhoq, rhoq_crg,  & ! in
       xq, dq_xave,  vt_xave, & ! in
    !       rho                            ! [Add] 11/08/30
       PQ,                    & ! inout
       Pcrg1, Pcrg2,          & ! inout
       Pac                    ) ! out
    use scale_atmos_saturation, only: &
       moist_psat_ice => ATMOS_SATURATION_psat_ice
    implicit none

    integer, intent(in) :: KA, KS, KE

    !--- mixed-phase collection process
    !                  And all we set all production term as a negative sign to avoid confusion.
    !
    real(RP), intent(in) :: wtem(KA)
    !--- mass/number concentration[kg/m3]
    real(RP), intent(in) :: rhoq(KA,I_QV:I_NG)
    ! necessary ?
    real(RP), intent(in) :: xq(KA,HYDRO_MAX)
    !--- diameter of averaged mass( D(ave x) )
    real(RP), intent(in) :: dq_xave(KA,HYDRO_MAX)
    !--- terminal velocity of averaged mass( vt(ave x) )
    real(RP), intent(in) :: vt_xave(KA,HYDRO_MAX,2)
    ! [Add] 11/08/30 T.Mitsui, for autoconversion of ice
    !    real(RP), intent(in) :: rho(KA)
    !--- partial conversion
    real(RP), intent(inout):: PQ(KA,PQ_MAX)
    !
    real(RP), intent(out):: Pac(KA,Pac_MAX)
    !--- for lightning component
    logical,  intent(in) :: flg_lt
    real(RP), intent(in) :: beta_crg(KA)
    real(RP), intent(in) :: dqcrg(KA)
    real(RP), intent(in) :: d0_crg, v0_crg
    real(RP), intent(in) :: rhoq_crg(KA,I_QC:I_QG)
    real(RP), intent(inout):: Pcrg1(KA,PQ_MAX)
    real(RP), intent(inout):: Pcrg2(KA,Pcrg_MAX)

    real(RP), parameter :: a_dec = 0.883_RP
    real(RP), parameter :: b_dec = 0.093_RP
    real(RP), parameter :: c_dec = 0.00348_RP
    real(RP), parameter :: d_dec = 4.5185E-5_RP
    !
    !
    real(RP) :: tem(KA)
    !
    !--- collection efficency of each specie
    real(RP) :: E_c(KA), E_r, E_i, E_s, E_g
    real(RP) :: E_ic, E_sc, E_gc
    !--- sticking efficiency
    real(RP) :: E_stick(KA)
    ! [Add] 10/08/03 T.Mitsui
    real(RP) :: temc, temc2, temc3
    real(RP) :: E_dec
    real(RP) :: esi_rat
    real(RP) :: esi(KA)
    !
    real(RP) :: temc_p, temc_m             ! celcius tem.
!    real(RP) :: ci_aut(KA)
!    real(RP) :: taui_aut(KA)
!    real(RP) :: tau_sce(KA)
    !--- DSD averaged diameter for each species
    real(RP) :: ave_dc                     ! cloud
!    real(RP) :: ave_dr                     ! rain
    real(RP) :: ave_di(KA)                 ! ice
    real(RP) :: ave_ds(KA)                 ! snow
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
    real(RP) :: coef_acc_LSG,   coef_acc_NSG(KA) ! snow      - graupel
    !--- (diameter) x (diameter)
    real(RP) :: dcdc(KA), dcdi, dcds, dcdg
    real(RP) :: drdr(KA), drdi(KA), drds(KA), drdg
    real(RP) :: didi(KA), dids, didg
    real(RP) :: dsds(KA), dsdg
    real(RP) :: dgdg(KA)
    !--- (terminal velocity) x (terminal velocity)
    real(RP) :: vcvc(KA), vcvi, vcvs, vcvg
    real(RP) :: vrvr(KA), vrvi(KA), vrvs(KA), vrvg
    real(RP) :: vivi(KA), vivs, vivg
    real(RP) :: vsvs(KA), vsvg
    real(RP) :: vgvg(KA)
    !
    real(RP) :: wx_cri, wx_crs
    real(RP) :: coef_emelt
    real(RP) :: w1

    real(RP) :: sw, sw1, sw2
    real(RP) :: alpha_lt
    !
    integer :: k, iqw
    !
    !
    do k = KS, KE
       tem(k) = max( wtem(k), tem_min ) ! 11/08/30 T.Mitsui
    end do

    call moist_psat_ice( KA, KS, KE, &
                         tem(:), esi(:)  ) ! [IN], [OUT]

    if( opt_stick_KS96 )then
       do k = KS, KE
          ! Khain and Sednev (1996), eq.(3.15)
          temc       = tem(k) - T00
          temc2      = temc*temc
          temc3      = temc2*temc
          E_dec      = max(0.0_RP, a_dec + b_dec*temc + c_dec*temc2 + d_dec*temc3 )
          esi_rat    = rhoq(k,I_QV)*Rvap*tem(k)/esi(k)
          E_stick(k) = min(1.0_RP, E_dec*esi_rat)
       end do
    else if( opt_stick_CO86 )then
       do k = KS, KE
          ! [Add] 11/08/30 T.Mitsui, Cotton et al. (1986)
          temc       = min(tem(k) - T00,0.0_RP)
          w1         = 0.035_RP*temc-0.7_RP
          E_stick(k) = 10._RP**w1
       end do
    else
       do k = KS, KE
          ! Lin et al. (1983)
          temc_m     = min(tem(k) - T00,0.0_RP) ! T < 273.15
          E_stick(k) = exp(0.09_RP*temc_m)
       end do
    end if

    do k = KS, KE
       ! averaged diameter using SB06(82)
       ave_dc = coef_d(I_mp_QC)*xq(k,I_mp_QC)**b_m(I_mp_QC)
       !------------------------------------------------------------------------
       ! coellection efficiency are given as follows
       E_c(k) = max(0.0_RP, min(1.0_RP, (ave_dc-dc0)/(dc1-dc0) ))
    end do

    !------------------------------------------------------------------------
    ! Collection:  a collects b ( assuming particle size a>b )
    do k = KS, KE
       dcdc(k) = dq_xave(k,I_mp_QC) * dq_xave(k,I_mp_QC)
       drdr(k) = dq_xave(k,I_mp_QR) * dq_xave(k,I_mp_QR)
       didi(k) = dq_xave(k,I_mp_QI) * dq_xave(k,I_mp_QI)
       dsds(k) = dq_xave(k,I_mp_QS) * dq_xave(k,I_mp_QS)
       dgdg(k) = dq_xave(k,I_mp_QG) * dq_xave(k,I_mp_QG)
       drdi(k) = dq_xave(k,I_mp_QR) * dq_xave(k,I_mp_QI)
       drds(k) = dq_xave(k,I_mp_QR) * dq_xave(k,I_mp_QS)
    end do
    do k = KS, KE
       vcvc(k) = vt_xave(k,I_mp_QC,2) * vt_xave(k,I_mp_QC,2)
       vrvr(k) = vt_xave(k,I_mp_QR,2) * vt_xave(k,I_mp_QR,2)
       vivi(k) = vt_xave(k,I_mp_QI,2) * vt_xave(k,I_mp_QI,2)
       vsvs(k) = vt_xave(k,I_mp_QS,2) * vt_xave(k,I_mp_QS,2)
       vgvg(k) = vt_xave(k,I_mp_QG,2) * vt_xave(k,I_mp_QG,2)
       vrvi(k) = vt_xave(k,I_mp_QR,2) * vt_xave(k,I_mp_QI,2)
       vrvs(k) = vt_xave(k,I_mp_QR,2) * vt_xave(k,I_mp_QS,2)
    end do

    !------------------------------------------------------------------------
    !
    !+++ pattern 1: a + b => a  (a>b)
    !                           (i-c, s-c, g-c, s-i, g-r, s-g)
    !------------------------------------------------------------------------

    ! cloud-ice => ice
    ! reduction term of cloud
    do k = KS, KE
       ave_di(k) = coef_d(I_mp_QI)*xq(k,I_mp_QI)**b_m(I_mp_QI)

       sw = 0.5_RP - sign(0.5_RP, di0-ave_di(k)) ! if(ave_di>di0)then sw=1
       E_i = E_im * sw
       E_ic = E_i*E_c(k)

       dcdi = dq_xave(k,I_mp_QC) * dq_xave(k,I_mp_QI)
       vcvi = vt_xave(k,I_mp_QC,2) * vt_xave(k,I_mp_QI,2)

       coef_acc_LCI = &
              ( delta_b1(I_mp_QC)*dcdc(k) + delta_ab1(I_mp_QI,I_mp_QC)*dcdi + delta_b0(I_mp_QI)*didi(k) ) &
        * sqrt( theta_b1(I_mp_QC)*vcvc(k) - theta_ab1(I_mp_QI,I_mp_QC)*vcvi + theta_b0(I_mp_QI)*vivi(k) &
            +  sigma_i + sigma_c )
       coef_acc_NCI = &
              ( delta_b0(I_mp_QC)*dcdc(k) + delta_ab0(I_mp_QI,I_mp_QC)*dcdi + delta_b0(I_mp_QI)*didi(k) ) &
        * sqrt( theta_b0(I_mp_QC)*vcvc(k) - theta_ab0(I_mp_QI,I_mp_QC)*vcvi + theta_b0(I_mp_QI)*vivi(k) &
            +  sigma_i + sigma_c )
       Pac(k,I_LIacLC2LI)= -0.25_RP*pi*E_ic*rhoq(k,I_NI)*rhoq(k,I_QC)*coef_acc_LCI
       Pac(k,I_NIacNC2NI)= -0.25_RP*pi*E_ic*rhoq(k,I_NI)*rhoq(k,I_NC)*coef_acc_NCI
    end do

    ! cloud-snow => snow
    ! reduction term of cloud
    do k = KS, KE
       ave_ds(k) = coef_d(I_mp_QS)*xq(k,I_mp_QS)**b_m(I_mp_QS)
       sw = 0.5_RP - sign(0.5_RP, ds0-ave_ds(k)) ! if(ave_ds>ds0)then sw=1
       E_s = E_sm * sw
       E_sc = E_s*E_c(k)

       dcds = dq_xave(k,I_mp_QC) * dq_xave(k,I_mp_QS)
       vcvs = vt_xave(k,I_mp_QC,2) * vt_xave(k,I_mp_QS,2)

       coef_acc_LCS = &
              ( delta_b1(I_mp_QC)*dcdc(k) + delta_ab1(I_mp_QS,I_mp_QC)*dcds + delta_b0(I_mp_QS)*dsds(k) ) &
        * sqrt( theta_b1(I_mp_QC)*vcvc(k) - theta_ab1(I_mp_QS,I_mp_QC)*vcvs + theta_b0(I_mp_QS)*vsvs(k) &
            +  sigma_s + sigma_c )
       coef_acc_NCS = &
              ( delta_b0(I_mp_QC)*dcdc(k) + delta_ab0(I_mp_QS,I_mp_QC)*dcds + delta_b0(I_mp_QS)*dsds(k) ) &
        * sqrt( theta_b0(I_mp_QC)*vcvc(k) - theta_ab0(I_mp_QS,I_mp_QC)*vcvs + theta_b0(I_mp_QS)*vsvs(k) &
            +  sigma_s + sigma_c )
       Pac(k,I_LSacLC2LS)= -0.25_RP*pi*E_sc*rhoq(k,I_NS)*rhoq(k,I_QC)*coef_acc_LCS
       Pac(k,I_NSacNC2NS)= -0.25_RP*pi*E_sc*rhoq(k,I_NS)*rhoq(k,I_NC)*coef_acc_NCS
    end do

    ! cloud-graupel => graupel
    ! reduction term of cloud
    do k = KS, KE
       ave_dg = coef_d(I_mp_QG)*xq(k,I_mp_QG)**b_m(I_mp_QG)
       sw = 0.5_RP - sign(0.5_RP, dg0-ave_dg) ! if(ave_dg>dg0)then sw=1
       E_g = E_gm * sw
       E_gc = E_g*E_c(k)

       dcdg = dq_xave(k,I_mp_QC) * dq_xave(k,I_mp_QG)
       vcvg = vt_xave(k,I_mp_QC,2) * vt_xave(k,I_mp_QG,2)

       coef_acc_LCG = &
              ( delta_b1(I_mp_QC)*dcdc(k) + delta_ab1(I_mp_QG,I_mp_QC)*dcdg + delta_b0(I_mp_QG)*dgdg(k) ) &
        * sqrt( theta_b1(I_mp_QC)*vcvc(k) - theta_ab1(I_mp_QG,I_mp_QC)*vcvg + theta_b0(I_mp_QG)*vgvg(k) &
            +  sigma_g + sigma_c )
       coef_acc_NCG = &
              ( delta_b0(I_mp_QC)*dcdc(k) + delta_ab0(I_mp_QG,I_mp_QC)*dcdg + delta_b0(I_mp_QG)*dgdg(k) ) &
        * sqrt( theta_b0(I_mp_QC)*vcvc(k) - theta_ab0(I_mp_QG,I_mp_QC)*vcvg + theta_b0(I_mp_QG)*vgvg(k) &
            +  sigma_g + sigma_c )
       Pac(k,I_LGacLC2LG)= -0.25_RP*pi*E_gc*rhoq(k,I_NG)*rhoq(k,I_QC)*coef_acc_LCG
       Pac(k,I_NGacNC2NG)= -0.25_RP*pi*E_gc*rhoq(k,I_NG)*rhoq(k,I_NC)*coef_acc_NCG
    end do

    ! snow-graupel => graupel
    do k = KS, KE
       dsdg = dq_xave(k,I_mp_QS) * dq_xave(k,I_mp_QG)
       vsvg = vt_xave(k,I_mp_QS,2) * vt_xave(k,I_mp_QG,2)

       coef_acc_LSG = &
              ( delta_b1(I_mp_QS)*dsds(k) + delta_ab1(I_mp_QG,I_mp_QS)*dsdg + delta_b0(I_mp_QG)*dgdg(k) ) &
        * sqrt( theta_b1(I_mp_QS)*vsvs(k) - theta_ab1(I_mp_QG,I_mp_QS)*vsvg + theta_b0(I_mp_QG)*vgvg(k) &
            +  sigma_g + sigma_s )
       coef_acc_NSG(k) = &
              ( delta_b0(I_mp_QS)*dsds(k) + delta_ab0(I_mp_QG,I_mp_QS)*dsdg + delta_b0(I_mp_QG)*dgdg(k) ) &
        * sqrt( theta_b0(I_mp_QS)*vsvs(k) - theta_ab0(I_mp_QG,I_mp_QS)*vsvg + theta_b0(I_mp_QG)*vgvg(k) &
            +  sigma_g + sigma_s )
       Pac(k,I_LGacLS2LG)= -0.25_RP*pi*E_stick(k)*E_gs*rhoq(k,I_NG)*rhoq(k,I_QS)*coef_acc_LSG
       Pac(k,I_NGacNS2NG)= -0.25_RP*pi*E_stick(k)*E_gs*rhoq(k,I_NG)*rhoq(k,I_NS)*coef_acc_NSG(k)
    end do

!!    !-----------------
!!    !  (start) Y.Sato added on 2018/8/31
!!    !--- ice-graupel => graupel
!!    do k = KS, KE
!!       didg = dq_xave(k,I_mp_QI) * dq_xave(k,I_mp_QG)
!!       vivg = vt_xave(k,I_mp_QI,2)* vt_xave(k,I_mp_QG,2)
!!       coef_acc_LIG = &
!!              ( delta_b1(I_QI)*didi(k) + delta_ab1(I_QG,I_QI)*didg + delta_b0(I_QG)*dgdg(k) ) &
!!        * sqrt( theta_b1(I_QI)*vivi(k) - theta_ab1(I_QG,I_QI)*vivg + theta_b0(I_QG)*vgvg(k) &
!!            +  sigma_g + sigma_i )
!!       coef_acc_NIG(k) = &
!!              ( delta_b0(I_QI)*didi(k) + delta_ab0(I_QG,I_QI)*didg + delta_b0(I_QG)*dgdg(k) ) &
!!            ! [fix] T.Mitsui 08/05/08
!!        * sqrt( theta_b0(I_QI)*vivi(k) - theta_ab0(I_QG,I_QI)*vivg + theta_b0(I_QG)*vgvg(k) &
!!            +  sigma_g + sigma_i )
!!       Pac(k,I_LGacLI2LG)= -0.25_RP*pi*E_stick(k)*E_gi*rhoq(k,I_NG)*rhoq(k,I_QI)*coef_acc_LIG*flg_igcol
!!       Pac(k,I_NGacNI2NG)= -0.25_RP*pi*E_stick(k)*E_gi*rhoq(k,I_NG)*rhoq(k,I_NI)*coef_acc_NIG(k)*flg_igcol
!!       !  (end) Y.Sato added on 2018/8/31
!!       !------------------
!!    end do

    !------------------------------------------------------------------------
    ! ice-snow => snow
    ! reduction term of ice
    do k = KS, KE
       dids = dq_xave(k,I_mp_QI) * dq_xave(k,I_mp_QS)
       vivs = vt_xave(k,I_mp_QI,2) * vt_xave(k,I_mp_QS,2)

       coef_acc_LIS = &
              ( delta_b1(I_mp_QI)*didi(k) + delta_ab1(I_mp_QS,I_mp_QI)*dids + delta_b0(I_mp_QS)*dsds(k) ) &
        * sqrt( theta_b1(I_mp_QI)*vivi(k) - theta_ab1(I_mp_QS,I_mp_QI)*vivs + theta_b0(I_mp_QS)*vsvs(k) &
            +  sigma_i + sigma_s )
       coef_acc_NIS = &
              ( delta_b0(I_mp_QI)*didi(k) + delta_ab0(I_mp_QS,I_mp_QI)*dids + delta_b0(I_mp_QS)*dsds(k) ) &
        * sqrt( theta_b0(I_mp_QI)*vivi(k) - theta_ab0(I_mp_QS,I_mp_QI)*vivs + theta_b0(I_mp_QS)*vsvs(k) &
            +  sigma_i + sigma_s )
       Pac(k,I_LIacLS2LS)= -0.25_RP*pi*E_stick(k)*E_si*rhoq(k,I_NS)*rhoq(k,I_QI)*coef_acc_LIS
       Pac(k,I_NIacNS2NS)= -0.25_RP*pi*E_stick(k)*E_si*rhoq(k,I_NS)*rhoq(k,I_NI)*coef_acc_NIS
    end do

    do k = KS, KE
       sw = sign(0.5_RP, T00-tem(k)) + 0.5_RP
       ! if ( tem(k) <= T00 )then
          ! rain-graupel => graupel
          ! reduction term of rain
          ! sw = 1
       ! else
          ! rain-graupel => rain
          ! reduction term of graupel
          ! sw = 0

       drdg = dq_xave(k,I_mp_QR) * dq_xave(k,I_mp_QG)
       vrvg = vt_xave(k,I_mp_QR,2) * vt_xave(k,I_mp_QG,2)

       coef_acc_LRG = &
              ( ( delta_b1(I_mp_QR)*drdr(k) + delta_ab1(I_mp_QG,I_mp_QR)*drdg + delta_b0(I_mp_QG)*dgdg(k) ) * sw &
              + ( delta_b1(I_mp_QG)*dgdg(k) + delta_ab1(I_mp_QR,I_mp_QG)*drdg + delta_b0(I_mp_QR)*drdr(k) ) * (1.0_RP-sw) ) &
            * sqrt( ( theta_b1(I_mp_QR)*vrvr(k) - theta_ab1(I_mp_QG,I_mp_QR)*vrvg + theta_b0(I_mp_QG)*vgvg(k) ) * sw &
                  + ( theta_b1(I_mp_QG)*vgvg(k) - theta_ab1(I_mp_QR,I_mp_QG)*vrvg + theta_b0(I_mp_QR)*vrvr(k) ) * (1.0_RP-sw) &
                  + sigma_r + sigma_g )
       Pac(k,I_LRacLG2LG) = -0.25_RP*pi*E_gr*coef_acc_LRG &
            * ( rhoq(k,I_NG)*rhoq(k,I_QR) * sw &
              + rhoq(k,I_NR)*rhoq(k,I_QG) * (1.0_RP-sw) )
       coef_acc_NRG = &
              ( delta_b0(I_mp_QR)*drdr(k) + delta_ab0(I_mp_QG,I_mp_QR)*drdg + delta_b0(I_mp_QG)*dgdg(k) ) &
        * sqrt( theta_b0(I_mp_QR)*vrvr(k) - theta_ab0(I_mp_QG,I_mp_QR)*vrvg + theta_b0(I_mp_QG)*vgvg(k) &
            +  sigma_r + sigma_g )
       Pac(k,I_NRacNG2NG) = -0.25_RP*pi*E_gr*rhoq(k,I_NG)*rhoq(k,I_NR)*coef_acc_NRG
    end do

    !------------------------------------------------------------------------
    !
    !+++ pattern 2: a + b => c  (a>b)
    !                           (r-i,r-s)
    !------------------------------------------------------------------------

    ! rain-ice => graupel
    ! reduction term of ice
    do k = KS, KE
       coef_acc_LRI_I = &
              ( delta_b1(I_mp_QI)*didi(k) + delta_ab1(I_mp_QR,I_mp_QI)*drdi(k) + delta_b0(I_mp_QR)*drdr(k) ) &
        * sqrt( theta_b1(I_mp_QI)*vivi(k) - theta_ab1(I_mp_QR,I_mp_QI)*vrvi(k) + theta_b0(I_mp_QR)*vrvr(k) &
            +  sigma_r + sigma_i )
       coef_acc_NRI_I = &
              ( delta_b0(I_mp_QI)*didi(k) + delta_ab0(I_mp_QR,I_mp_QI)*drdi(k) + delta_b0(I_mp_QR)*drdr(k) ) &
        * sqrt( theta_b0(I_mp_QI)*vivi(k) - theta_ab0(I_mp_QR,I_mp_QI)*vrvi(k) + theta_b0(I_mp_QR)*vrvr(k) &
            +  sigma_r + sigma_i )
       Pac(k,I_LRacLI2LG_I)= -0.25_RP*pi*E_ir*rhoq(k,I_NR)*rhoq(k,I_QI)*coef_acc_LRI_I
       Pac(k,I_NRacNI2NG_I)= -0.25_RP*pi*E_ir*rhoq(k,I_NR)*rhoq(k,I_NI)*coef_acc_NRI_I
    end do

    ! reduction term of rain
    do k = KS, KE
       coef_acc_LRI_R = &
              ( delta_b1(I_mp_QR)*drdr(k) + delta_ab1(I_mp_QI,I_mp_QR)*drdi(k) + delta_b0(I_mp_QI)*didi(k) ) &
        * sqrt( theta_b1(I_mp_QR)*vrvr(k) - theta_ab1(I_mp_QI,I_mp_QR)*vrvi(k) + theta_b0(I_mp_QI)*vivi(k) &
            +  sigma_r + sigma_i )
       coef_acc_NRI_R = &
              ( delta_b0(I_mp_QR)*drdr(k) + delta_ab0(I_mp_QI,I_mp_QR)*drdi(k) + delta_b0(I_mp_QI)*didi(k) ) &
        * sqrt( theta_b0(I_mp_QR)*vrvr(k) - theta_ab0(I_mp_QI,I_mp_QR)*vrvi(k) + theta_b0(I_mp_QI)*vivi(k) &
            +  sigma_r + sigma_i )
       Pac(k,I_LRacLI2LG_R)= -0.25_RP*pi*E_ir*rhoq(k,I_NI)*rhoq(k,I_QR)*coef_acc_LRI_R
       Pac(k,I_NRacNI2NG_R)= -0.25_RP*pi*E_ir*rhoq(k,I_NI)*rhoq(k,I_NR)*coef_acc_NRI_R
    end do

    ! rain-snow => graupel
    ! reduction term of snow
    do k = KS, KE
       coef_acc_LRS_S = &
              ( delta_b1(I_mp_QS)*dsds(k) + delta_ab1(I_mp_QR,I_mp_QS)*drds(k) + delta_b0(I_mp_QR)*drdr(k) ) &
        * sqrt( theta_b1(I_mp_QS)*vsvs(k) - theta_ab1(I_mp_QR,I_mp_QS)*vrvs(k) + theta_b0(I_mp_QR)*vrvr(k) &
            +  sigma_r + sigma_s )
       coef_acc_NRS_S = &
              ( delta_b0(I_mp_QS)*dsds(k) + delta_ab0(I_mp_QR,I_mp_QS)*drds(k) + delta_b0(I_mp_QR)*drdr(k) ) &
        * sqrt( theta_b0(I_mp_QS)*vsvs(k) - theta_ab0(I_mp_QR,I_mp_QS)*vrvs(k) + theta_b0(I_mp_QR)*vrvr(k) &
            +  sigma_r + sigma_s )
       Pac(k,I_LRacLS2LG_S)= -0.25_RP*pi*E_sr*rhoq(k,I_NR)*rhoq(k,I_QS)*coef_acc_LRS_S
       Pac(k,I_NRacNS2NG_S)= -0.25_RP*pi*E_sr*rhoq(k,I_NR)*rhoq(k,I_NS)*coef_acc_NRS_S
    end do

    ! reduction term of rain
    do k = KS, KE
       coef_acc_LRS_R = &
              ( delta_b1(I_mp_QR)*drdr(k) + delta_ab1(I_mp_QS,I_mp_QR)*drds(k) + delta_b0(I_mp_QS)*dsds(k) ) &
        * sqrt( theta_b1(I_mp_QR)*vrvr(k) - theta_ab1(I_mp_QS,I_mp_QR)*vrvs(k) + theta_b0(I_mp_QS)*vsvs(k) &
            +  sigma_r + sigma_s )
       coef_acc_NRS_R = &
              ( delta_b0(I_mp_QR)*drdr(k) + delta_ab0(I_mp_QS,I_mp_QR)*drds(k) + delta_b0(I_mp_QS)*dsds(k) ) &
        * sqrt( theta_b0(I_mp_QR)*vrvr(k) - theta_ab0(I_mp_QS,I_mp_QR)*vrvs(k) + theta_b0(I_mp_QS)*vsvs(k) &
            +  sigma_r + sigma_s )
       Pac(k,I_LRacLS2LG_R)= -0.25_RP*pi*E_sr*rhoq(k,I_NS)*rhoq(k,I_QR)*coef_acc_LRS_R
       Pac(k,I_NRacNS2NG_R)= -0.25_RP*pi*E_sr*rhoq(k,I_NS)*rhoq(k,I_NR)*coef_acc_NRS_R
    end do

    !------------------------------------------------------------------------
    !
    !+++ pattern 3: a + a => b  (i-i)
    !
    !------------------------------------------------------------------------

    ! ice-ice ( reduction is double, but includes double-count)
    do k = KS, KE
       coef_acc_LII = &
              ( delta_b0(I_mp_QI)*didi(k) + delta_ab1(I_mp_QI,I_mp_QI)*didi(k) + delta_b1(I_mp_QI)*didi(k) ) &
        * sqrt( theta_b0(I_mp_QI)*vivi(k) - theta_ab1(I_mp_QI,I_mp_QI)*vivi(k) + theta_b1(I_mp_QI)*vivi(k) &
            +  sigma_i + sigma_i )
       coef_acc_NII = &
              ( delta_b0(I_mp_QI)*didi(k) + delta_ab0(I_mp_QI,I_mp_QI)*didi(k) + delta_b0(I_mp_QI)*didi(k) ) &
        * sqrt( theta_b0(I_mp_QI)*vivi(k) - theta_ab0(I_mp_QI,I_mp_QI)*vivi(k) + theta_b0(I_mp_QI)*vivi(k) &
            +  sigma_i + sigma_i )
       Pac(k,I_LIacLI2LS)= -0.25_RP*pi*E_stick(k)*E_ii*rhoq(k,I_NI)*rhoq(k,I_QI)*coef_acc_LII
       Pac(k,I_NIacNI2NS)= -0.25_RP*pi*E_stick(k)*E_ii*rhoq(k,I_NI)*rhoq(k,I_NI)*coef_acc_NII
       !
!          ci_aut(k)   =  0.25_RP*pi*E_ii*rhoq(k,I_NI)*coef_acc_LII
!          taui_aut(k) = 1._RP/max(E_stick(k)*ci_aut(k),1.E-10_RP)
!          tau_sce(k)  = rhoq(k,I_QI)/max(rhoq(k,I_QIj)+rhoq(k,I_QS),1.E-10_RP)
    end do

    !------------------------------------------------------------------------
    !
    !+++ pattern 4: a + a => a  (s-s)
    !
    !------------------------------------------------------------------------

    ! snow-snow => snow
    do k = KS, KE
       coef_acc_NSS = &
              ( delta_b0(I_mp_QS)*dsds(k) + delta_ab0(I_mp_QS,I_mp_QS)*dsds(k) + delta_b0(I_mp_QS)*dsds(k) ) &
        * sqrt( theta_b0(I_mp_QS)*vsvs(k) - theta_ab0(I_mp_QS,I_mp_QS)*vsvs(k) + theta_b0(I_mp_QS)*vsvs(k) &
            +  sigma_s + sigma_s )
       Pac(k,I_NSacNS2NS)= -0.125_RP*pi*E_stick(k)*E_ss*rhoq(k,I_NS)*rhoq(k,I_NS)*coef_acc_NSS
    end do

    ! graupel-grauple => graupel
    do k = KS, KE
       coef_acc_NGG = &
              ( delta_b0(I_mp_QG)*dgdg(k) + delta_ab0(I_mp_QG,I_mp_QG)*dgdg(k) + delta_b0(I_mp_QG)*dgdg(k) ) &
        * sqrt( theta_b0(I_mp_QG)*vgvg(k) - theta_ab0(I_mp_QG,I_mp_QG)*vgvg(k) + theta_b0(I_mp_QG)*vgvg(k) &
            +  sigma_g + sigma_g )
       Pac(k,I_NGacNG2NG)= -0.125_RP*pi*E_stick(k)*E_gg*rhoq(k,I_NG)*rhoq(k,I_NG)*coef_acc_NGG
    end do

    !------------------------------------------------------------------------
    !--- Partial conversion
    ! SB06(70),(71)
    ! i_iconv2g: option whether partial conversions work or not
    ! ice-cloud => graupel
    do k = KS, KE
       sw = 0.5_RP - sign(0.5_RP,di_cri-ave_di(k)) ! if( ave_di > di_cri )then sw=1
       wx_cri = cfill_i*DWATR/rho_g*( pi/6.0_RP*rho_g*ave_di(k)**3/xq(k,I_mp_QI) - 1.0_RP ) * sw
       PQ(k,I_LIcon) = i_iconv2g * Pac(k,I_LIacLC2LI)/max(1.0_RP, wx_cri) * sw
       PQ(k,I_NIcon) = i_iconv2g * PQ(k,I_LIcon)/xq(k,I_mp_QI) * sw
    end do

    ! snow-cloud => graupel
    do k = KS, KE
       wx_crs = cfill_s*DWATR/rho_g*( pi/6.0_RP*rho_g*ave_ds(k)**3/xq(k,I_mp_QS) - 1.0_RP )
       PQ(k,I_LScon) = i_sconv2g * (Pac(k,I_LSacLC2LS))/max(1.0_RP, wx_crs)
       PQ(k,I_NScon) = i_sconv2g * PQ(k,I_LScon)/xq(k,I_mp_QS)
    end do

    !------------------------------------------------------------------------
    !--- enhanced melting( due to collection-freezing of water droplets )
    !    originally from Rutledge and Hobbs(1984). eq.(A.21)
    ! if T > 273.15 then temc_p=T-273.15, else temc_p=0
    ! 08/05/08 [fix] T.Mitsui LHF00 => LHF0
    ! melting occurs around T=273K, so LHF0 is suitable both SIMPLE and EXACT,
    ! otherwise LHF can have sign both negative(EXACT) and positive(SIMPLE).
    do k = KS, KE
!          temc_m = min(tem(k) - T00,0.0_RP) ! T < 273.15
       temc_p = max(tem(k) - T00,0.0_RP) ! T > 273.15
!!$       coef_emelt   = -CL/LHF00*temc_p
       coef_emelt   =  CL/LHF0*temc_p
       ! cloud-graupel
       PQ(k,I_LGacm) =  coef_emelt*Pac(k,I_LGacLC2LG)
       PQ(k,I_NGacm) =  PQ(k,I_LGacm)/xq(k,I_mp_QG)
       ! rain-graupel
       PQ(k,I_LGarm) =  coef_emelt*Pac(k,I_LRacLG2LG)
       PQ(k,I_NGarm) =  PQ(k,I_LGarm)/xq(k,I_mp_QG)
       ! cloud-snow
       PQ(k,I_LSacm) =  coef_emelt*(Pac(k,I_LSacLC2LS))
       PQ(k,I_NSacm) =  PQ(k,I_LSacm)/xq(k,I_mp_QS)
       ! rain-snow
       PQ(k,I_LSarm) =  coef_emelt*(Pac(k,I_LRacLS2LG_R)+Pac(k,I_LRacLS2LG_S))
       PQ(k,I_NSarm) =  PQ(k,I_LSarm)/xq(k,I_mp_QG)
       ! cloud-ice
       PQ(k,I_LIacm) =  coef_emelt*Pac(k,I_LIacLC2LI)
       PQ(k,I_NIacm) =  PQ(k,I_LIacm)/xq(k,I_mp_QI)
       ! rain-ice
       PQ(k,I_LIarm) =  coef_emelt*(Pac(k,I_LRacLI2LG_R)+Pac(k,I_LRacLI2LG_I))
       PQ(k,I_NIarm) =  PQ(k,I_LIarm)/xq(k,I_mp_QG)
    end do


    if ( flg_lt ) then
       !--- C + I -> I (decrease from cloud chgarge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NC)-SMALL ) !--- if NC is small,  ignore charge transfer
          Pcrg2(k,I_NIacNC2NI) = Pac(k,I_NIacNC2NI)*(1.0_RP-sw1) / (rhoq(k,I_NC)+sw1)*rhoq_crg(k,I_QC)
       end do

       !--- C + S -> S (decrease from cloud charge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NC)-SMALL ) !--- if NC is small,  ignore charge transfer
          Pcrg2(k,I_NSacNC2NS) = Pac(k,I_NSacNC2NS)*(1.0_RP-sw1) / (rhoq(k,I_NC)+sw1)*rhoq_crg(k,I_QC)
       end do

       !--- C + G -> G (decrease from cloud charge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NC)-SMALL ) !--- if NC is small,  ignore charge transfer
          Pcrg2(k,I_NGacNC2NG) = Pac(k,I_NGacNC2NG)*(1.0_RP-sw1) / (rhoq(k,I_NC)+sw1)*rhoq_crg(k,I_QC)
       end do

       !--- S + G -> G (decrease from snow charge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NS)-SMALL ) !--- if NS is small,  ignore charge transfer
          Pcrg2(k,I_NGacNS2NG) = Pac(k,I_NGacNS2NG)*(1.0_RP-sw1) / (rhoq(k,I_NS)+sw1)*rhoq_crg(k,I_QS)
       end do

       !--- Charge split by Snow-Graupel rebound--------------------------------
       do k = KS, KE
          alpha_lt = 5.0_RP * ( dq_xave(k,I_mp_QS) / d0_crg )**2 * vt_xave(k,I_mp_QG,2) / v0_crg
          alpha_lt = min( alpha_lt, 10.0_RP )
          Pcrg2(k,I_CGNGacNS2NG)= 0.25_RP*pi*( 1.0_RP - E_stick(k) )*E_gs &
                                * rhoq(k,I_NG)*rhoq(k,I_NS)*coef_acc_NSG(k) &
                                * ( dqcrg(k)*alpha_lt ) &
                                * beta_crg(k)
       end do

!!      !--- I + G -> G (decrease from snow charge)
!!       do k = KS, KE
!!         sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NI)-SMALL ) !--- if NS is small,  ignore charge transfer
!!         Pcrg2(k,I_NGacNI2NG) = Pac(k,I_NGacNI2NG)*(1.0_RP-sw1) / (rhoq(k,I_NI)+sw1)*rhoq_crg(k,I_QI)*flg_igcol
!!          !--- Charge split by Ice-Graupel rebound--------------------------------
!!          alpha_lt = 5.0_RP * ( dq_xave(k,I_mp_QI) / d0_crg )**2 * vt_xave(k,I_mp_QG,2) / v0_crg
!!          alpha_lt = min( alpha_lt, 10.0_RP )
!!          Pcrg2(k,I_CGNGacNI2NG)= 0.25_RP*pi*( 1.0_RP - E_stick(k) )*E_gi &
!!                                * rhoq(k,I_NG)*rhoq(k,I_NI)*coef_acc_NIG(k) &
!!                                * ( dqcrg(k)*alpha_lt ) &
!!                                * beta_crg(k) * flg_igcol
!!       end do

       !--- I + S -> S (decrease from ice charge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NI)-SMALL ) !--- if NI is small,  ignore charge transfer
          Pcrg2(k,I_NIacNS2NS) = Pac(k,I_NIacNS2NS)*(1.0_RP-sw1) / (rhoq(k,I_NI)+sw1)*rhoq_crg(k,I_QI)
       end do

       !--- R+G->R (T>T00 sw=0, decrease from graupel charge), ->G(T<=T00 sw=1, dcrerase from rain charge)
       do k = KS, KE
          sw = sign(0.5_RP, T00-tem(k)) + 0.5_RP
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NR)-SMALL ) !--- if NR is small,  ignore charge transfer
          sw2 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NG)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg2(k,I_NRacNG2NG) = Pac(k,I_NRacNG2NG)*(1.0_RP-sw1)/(rhoq(k,I_NR)+sw1)*rhoq_crg(k,I_QR) * sw &
                               + Pac(k,I_NRacNG2NG)*(1.0_RP-sw2)/(rhoq(k,I_NG)+sw2)*rhoq_crg(k,I_QG) * (1.0_RP-sw)
       end do

       !--- R + I -> G (decrease from both ice and rain charge, but only ice charge at this part)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NI)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg2(k,I_NRacNI2NG_I) = Pac(k,I_NRacNI2NG_I)*(1.0_RP-sw1) / (rhoq(k,I_NI)+sw1)*rhoq_crg(k,I_QI)
       end do

       !--- R + I -> G (decrease from both ice and rain charge, but only rain charge at this part)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NR)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg2(k,I_NRacNI2NG_R) = Pac(k,I_NRacNI2NG_R)*(1.0_RP-sw1) / (rhoq(k,I_NR)+sw1)*rhoq_crg(k,I_QR)
       end do

       !--- R + S -> G (decrease from both snow and rain charge, but only snow charge at this part)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NS)-SMALL ) !--- if NS is small,  ignore charge transfer
          Pcrg2(k,I_NRacNS2NG_S) = Pac(k,I_NRacNS2NG_S)*(1.0_RP-sw1) / (rhoq(k,I_NS)+sw1)*rhoq_crg(k,I_QS)
       end do

       !--- R + S -> G (decrease from both snow and rain charge, but only rain charge at this part)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NR)-SMALL ) !--- if NR is small,  ignore charge transfer
          Pcrg2(k,I_NRacNS2NG_R) = Pac(k,I_NRacNS2NG_R)*(1.0_RP-sw1) / (rhoq(k,I_NR)+sw1)*rhoq_crg(k,I_QR)
       end do

       !--- I + I -> S (decrease from ice charge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NI)-SMALL ) !--- if NI is small,  ignore charge transfer
          Pcrg2(k,I_NIacNI2NS) = Pac(k,I_NIacNI2NS)*(1.0_RP-sw1) / (rhoq(k,I_NI)+sw1)*rhoq_crg(k,I_QI)
       end do

       !--- I + C -> G (decrease from ice charge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NI)-SMALL ) !--- if NI is small,  ignore charge transfer
          Pcrg1(k,I_NIcon) = i_iconv2g * PQ(k,I_NIcon)*(1.0_RP-sw) / (rhoq(k,I_NI)+sw1)*rhoq_crg(k,I_QI)
       end do

       !--- S + C -> G (decrease from snow charge)
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NS)-SMALL ) !--- if NS is small,  ignore charge transfer
          Pcrg1(k,I_NScon) = i_sconv2g * PQ(k,I_NScon)*(1.0_RP-sw) / (rhoq(k,I_NS)+sw1)*rhoq_crg(k,I_QS)
       end do

       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NG)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg1(k,I_NGacm) = PQ(k,I_NGacm)*(1.0_RP-sw1) / (rhoq(k,I_NG)+sw1)*rhoq_crg(k,I_QG)
       end do
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NG)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg1(k,I_NGarm) = PQ(k,I_NGarm)*(1.0_RP-sw1) / (rhoq(k,I_NG)+sw1)*rhoq_crg(k,I_QG)
       end do
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NS)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg1(k,I_NSacm) = PQ(k,I_NSacm)*(1.0_RP-sw1) / (rhoq(k,I_NS)+sw1)*rhoq_crg(k,I_QS)
       end do
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NS)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg1(k,I_NSarm) = PQ(k,I_NSarm)*(1.0_RP-sw1) / (rhoq(k,I_NS)+sw1)*rhoq_crg(k,I_QS)
       end do
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NI)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg1(k,I_NIacm) = PQ(k,I_NIacm)*(1.0_RP-sw1) / (rhoq(k,I_NI)+sw1)*rhoq_crg(k,I_QI)
       end do
       do k = KS, KE
          sw1 = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NR)-SMALL ) !--- if NG is small,  ignore charge transfer
          Pcrg1(k,I_NIarm) = PQ(k,I_NIarm)*(1.0_RP-sw1) / (rhoq(k,I_NR)+sw1)*rhoq_crg(k,I_QI)
       end do
    end if

    return
  end subroutine mixed_phase_collection
  !----------------------------
  ! Auto-conversion, Accretion, Self-collection, Break-up
!OCL SERIAL
  subroutine aut_acc_slc_brk(  &
       KA, KS, KE,             &
       flg_lt,                 &
       rhoq, rhoq_crg,         &
       xq, dq_xave,            &
       rho,                    &
       PQ, Pcrg                )
    implicit none

    integer, intent(in) :: KA, KS, KE
    !
    real(RP), intent(in)  :: rhoq(KA,I_QV:I_NG)
    real(RP), intent(in)  :: rhoq_crg(KA,I_QC:I_QG)
    logical,  intent(in)  :: flg_lt
    real(RP), intent(in)  :: xq(KA,HYDRO_MAX)
    real(RP), intent(in)  :: dq_xave(KA,HYDRO_MAX)
    real(RP), intent(in)  :: rho(KA)
    !
    real(RP), intent(inout) :: PQ(KA,PQ_MAX)
    real(RP), intent(inout) :: Pcrg(KA,PQ_MAX)
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
    integer :: k, iqw
    real(RP) :: sw
    !
    coef_nuc0 = (nu(I_mp_QC)+2.0_RP)/(nu(I_mp_QC)+1.0_RP)
    coef_nuc1 = (nu(I_mp_QC)+2.0_RP)*(nu(I_mp_QC)+4.0_RP)/(nu(I_mp_QC)+1.0_RP)/(nu(I_mp_QC)+1.0_RP)
    coef_aut0 =  -kcc*coef_nuc0
    coef_aut1 =  -kcc/x_sep/20._RP*coef_nuc1
    !
    do k = KS, KE
       lwc = rhoq(k,I_QR) + rhoq(k,I_QC)
       if( lwc > xc_min )then
          tau  = max(tau_min, rhoq(k,I_QR)/lwc)
       else
          tau  = tau_min
       end if
       rho_fac = sqrt(rho_0/max(rho(k),rho_min))
       !
       ! Auto-conversion ( cloud-cloud => rain )
       psi_aut       = 400._RP*(tau**0.7_RP)*(1.0_RP - (tau**0.7_RP))**3   ! (6) SB06
       PQ(k,I_NCaut)  = coef_aut0*rhoq(k,I_QC)*rhoq(k,I_QC)*rho_fac*rho_fac    ! (9) SB06 sc+aut
       ! lc = lwc*(1-tau), lr = lwc*tau
       PQ(k,I_LCaut)  = coef_aut1*lwc*lwc*xq(k,I_mp_QC)*xq(k,I_mp_QC) &          ! (4) SB06
            *((1.0_RP-tau)*(1.0_RP-tau) + psi_aut)*rho_fac*rho_fac        !
       PQ(k,I_NRaut)  = -rx_sep*PQ(k,I_LCaut)                           ! (A7) SB01
       !--- for Charge density
       if (flg_lt) then
          sw = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NC)-SMALL )
          Pcrg(k,I_NCaut) = PQ(k,I_NCaut)*(1.0_RP-sw)/(rhoq(k,I_NC)+sw)*rhoq_crg(k,I_QC)
          Pcrg(k,I_NRaut) = -Pcrg(k,I_NCaut)
       end if
       !
       ! Accretion ( cloud-rain => rain )
       psi_acc       =(tau/(tau+thr_acc))**4                          ! (8) SB06
       PQ(k,I_LCacc)  = -kcr*rhoq(k,I_QC)*rhoq(k,I_QR)*rho_fac*psi_acc         ! (7) SB06
       PQ(k,I_NCacc)  = -kcr*rhoq(k,I_NC)*rhoq(k,I_QR)*rho_fac*psi_acc         ! (A6) SB01
       !--- for Charge density
       if (flg_lt) then
          sw = 0.5_RP - sign( 0.5_RP, rhoq(k,I_NC)-SMALL )
          Pcrg(k,I_NCacc)  = PQ(k,I_NCacc)*(1.0_RP-sw)/(rhoq(k,I_NC)+sw)*rhoq(k,I_QC)
       end if
       !
       ! Self-collection ( rain-rain => rain )
       PQ(k,I_NRslc)  = -krr*rhoq(k,I_NR)*rhoq(k,I_QR)*rho_fac                 ! (A.8) SB01
       !
       ! Collisional breakup of rain
       ddr           = min(1.E-3_RP, dq_xave(k,I_mp_QR) - dr_eq )
       if      ( dq_xave(k,I_mp_QR) < dr_min ) then       ! negligible
          psi_brk      = -1.0_RP
       else if ( dq_xave(k,I_mp_QR) <= dr_eq  ) then
          psi_brk      = kbr*ddr                          ! (14) SB06
       else
          psi_brk      = exp(kapbr*ddr) - 1.0_RP          ! (15) SB06 (SB06 has a typo)
       end if
       PQ(k,I_NRbrk) = - (psi_brk + 1.0_RP)*PQ(k,I_NRslc) ! (13) SB06
       !
    end do
    !
    return
  end subroutine aut_acc_slc_brk
  ! Vapor Deposition, Ice Melting
!OCL SERIAL
  subroutine dep_vapor_melt_ice( &
       KA, KS, KE, &
       rho, tem, pre, qd,    & ! in
       rhoq,                 & ! in
       esw, esi,             & ! in
       xq, vt_xave, dq_xave, & ! in
       PQ                    ) ! inout
    use scale_const, only: &
       eps => CONST_EPS
    implicit none

    integer, intent(in) :: KA, KS, KE

    ! Diffusion growth or Evaporation, Sublimation
    real(RP), intent(inout) :: PQ(KA,PQ_MAX)  ! mass change   for cloud

    real(RP), intent(in)  :: rho(KA)     ! air density
    real(RP), intent(in)  :: tem(KA)     ! air temperature
    real(RP), intent(in)  :: pre(KA)     ! air pressure
    real(RP), intent(in)  :: qd (KA)      ! mixing ratio of dry air
    real(RP), intent(in)  :: esw(KA)     ! saturation vapor pressure(liquid water)
    real(RP), intent(in)  :: esi(KA)     ! saturation vapor pressure(solid water)
    real(RP), intent(in)  :: rhoq(KA,I_QV:I_NG)
    real(RP), intent(in)  :: xq(KA,HYDRO_MAX)    ! mean mass
    ! Notice following values differ from mean terminal velocity or diameter.
    ! mean(vt(x)) /= vt(mean(x)) and mean(D(x)) /= D(mean(x))
    ! Following ones are vt(mean(x)) and D(mean(x)).
    real(RP), intent(in)  :: vt_xave(KA,HYDRO_MAX,2) ! terminal velocity of mean cloud
    !
    real(RP), intent(in)  :: dq_xave(KA,HYDRO_MAX) ! diameter
    !
    real(RP) :: rho_lim            ! limited density
    real(RP) :: temc_lim           ! limited temperature[celsius]
    real(RP) :: pre_lim            ! limited density
    real(RP) :: temc               ! temperature[celsius]
!    real(RP) :: pv                 ! vapor pressure
    real(RP) :: qv                 ! mixing ratio of water vapor
!    real(RP) :: ssw                ! super saturation ratio(liquid water)
!    real(RP) :: ssi                ! super saturation ratio(ice water)
    real(RP) :: nua, r_nua         ! kinematic viscosity of air
    real(RP) :: mua                ! viscosity of air
    real(RP) :: Kalfa(KA)          ! thermal conductance
    real(RP) :: Dw                 ! diffusivity of water vapor
    real(RP) :: Dt                 ! diffusivity of heat
    real(RP) :: Gw, Gi             ! diffusion factor by balance between heat and vapor
    real(RP) :: Gwr, Gii, Gis, Gig ! for rain, ice, snow and graupel.
    real(RP) :: Gm                 ! melting factor by balance between heat and vapor
    real(RP) :: Nsc_r3             !
    ! [Mod] 11/08/30 T.Mitsui, considering large and small branches
!    real(RP) :: Nrecs_r2
    real(RP) :: Nrers_r2, Nreis_r2  !
    real(RP) :: Nress_r2, Nregs_r2  !
!    real(RP) :: Nrecl_r2
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
    real(RP) :: ventNI(KA), ventLI(KA)
    real(RP) :: ventNS(KA), ventLS(KA)
    real(RP) :: ventNG(KA), ventLG(KA)
    !
    real(RP), parameter :: Re_max=1.E+3_RP
    real(RP), parameter :: Re_min=1.E-4_RP

    real(RP) :: sw
    !
    integer :: k
    !
    ! Notice,T.Mitsui
    ! Vapor deposition and melting would not be solved iteratively to reach equilibrium.
    ! Because following phenomena are not adjustment but transition.
    ! Just time-scales differ among them.
    ! If we would treat more appropreately, there would be time-splitting method to solve each ones.

!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       temc    = tem(k) - T00   ! degC
       temc_lim= max(temc, -40._RP )       ! [Add] 09/08/18 T.Mitsui, Pruppacher and Klett(1997),(13-3)
       rho_lim = max(rho(k),rho_min)   ! [Add] 09/08/18 T.Mitsui
       qv      = rhoq(k,I_QV)/rho_lim
       pre_lim = rho_lim*(qd(k)*Rdry + qv*Rvap)*(temc_lim+T00) ![Add] 09/08/18 T.Mitsui
       !--------------------------------------------------------------------
       ! Diffusion growth part is described in detail
       ! by Pruppacher and Klett (1997) Sec. 13.2(liquid) and 13.3(solid)
       !
       ! G:factor of thermal diffusion(1st.term) and vapor diffusion(2nd. term)
       ! SB06(23),(38), Lin et al(31),(52) or others
       ! Dw is introduced by Pruppacher and Klett(1997),(13-3)
       Dw      = 0.211E-4_RP* (((temc_lim+T00)/T00)**1.94_RP) *(P00/pre_lim)
       Kalfa(k) = Ka0  + temc_lim*dKa_dT
       mua     = mua0 + temc_lim*dmua_dT
       nua     = mua/rho_lim
       r_nua   = 1.0_RP/nua
       Gw      = (LHV0/Kalfa(k)/tem(k))*(LHV0/Rvap/tem(k)-1.0_RP)+(Rvap*tem(k)/Dw/esw(k))
       Gi      = (LHS0/Kalfa(k)/tem(k))*(LHS0/Rvap/tem(k)-1.0_RP)+(Rvap*tem(k)/Dw/esi(k))
       ! capacities account for their surface geometries
       Gwr     = 4.0_RP*PI/cap(I_mp_QR)/Gw
       Gii     = 4.0_RP*PI/cap(I_mp_QI)/Gi
       Gis     = 4.0_RP*PI/cap(I_mp_QS)/Gi
       Gig     = 4.0_RP*PI/cap(I_mp_QG)/Gi
       ! vent: ventilation effect( asymmetry vapor field around particles due to aerodynamic )
       ! SB06 (30),(31) and each coefficient is by (88),(89)
       Nsc_r3  = (nua/Dw)**(0.33333333_RP)                    ! (Schmidt number )^(1/3)
       !
!       Nrecs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QC,1)*dq_xave(k,I_mp_QC)*r_nua))) ! (Reynolds number)^(1/2) cloud
       Nrers_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QR,1)*dq_xave(k,I_mp_QR)*r_nua))) ! (Reynolds number)^(1/2) rain
       Nreis_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QI,1)*dq_xave(k,I_mp_QI)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
       Nress_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QS,1)*dq_xave(k,I_mp_QS)*r_nua))) ! (Reynolds number)^(1/2) snow
       Nregs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QG,1)*dq_xave(k,I_mp_QG)*r_nua))) ! (Reynolds number)^(1/2) graupel
       !
!       Nrecl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QC,2)*dq_xave(k,I_mp_QC)*r_nua))) ! (Reynolds number)^(1/2) cloud
       Nrerl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QR,2)*dq_xave(k,I_mp_QR)*r_nua))) ! (Reynolds number)^(1/2) rain
       Nreil_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QI,2)*dq_xave(k,I_mp_QI)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
       Nresl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QS,2)*dq_xave(k,I_mp_QS)*r_nua))) ! (Reynolds number)^(1/2) snow
       Nregl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,I_mp_QG,2)*dq_xave(k,I_mp_QG)*r_nua))) ! (Reynolds number)^(1/2) graupel
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
       ventNI(k)  = (1.0_RP-wti)*ventNI_s + wti*ventNI_l
       ventNS(k)  = (1.0_RP-wts)*ventNS_s + wts*ventNS_l
       ventNG(k)  = (1.0_RP-wtg)*ventNG_s + wtg*ventNG_l
       !
       ventLR     = (1.0_RP-wtr)*ventLR_s + wtr*ventLR_l
       ventLI(k)  = (1.0_RP-wti)*ventLI_s + wti*ventLI_l
       ventLS(k)  = (1.0_RP-wts)*ventLS_s + wts*ventLS_l
       ventLG(k)  = (1.0_RP-wtg)*ventLG_s + wtg*ventLG_l
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
!!$         Hereafter PLxdep means inverse of timescale.
!!$***************************************************************************
!!$     PQ(k,I_LCdep) = Gwr*ssw*rhoq(k,I_NC)*dq_xave(k,I_mp_QC)*coef_deplc
!!$     PQ(k,I_LRdep) = Gwr*ssw*rhoq(k,I_NR)*dq_xave(k,I_mp_QR)*ventLR(k)
!!$     PQ(k,I_LIdep) = Gii*ssi*rhoq(k,I_NI)*dq_xave(k,I_mp_QI)*ventLI(k)
!!$     PQ(k,I_LSdep) = Gis*ssi*rhoq(k,I_NS)*dq_xave(k,I_mp_QS)*ventLS(k)
!!$     PQ(k,I_LGdep) = Gig*ssi*rhoq(k,I_NG)*dq_xave(k,I_mp_QG)*ventLG(k)
       PQ(k,I_LCdep) = Gwr*rhoq(k,I_NC)*dq_xave(k,I_mp_QC)*coef_deplc
       PQ(k,I_LRdep) = Gwr*rhoq(k,I_NR)*dq_xave(k,I_mp_QR)*ventLR
       PQ(k,I_LIdep) = Gii*rhoq(k,I_NI)*dq_xave(k,I_mp_QI)*ventLI(k)
       PQ(k,I_LSdep) = Gis*rhoq(k,I_NS)*dq_xave(k,I_mp_QS)*ventLS(k)
       PQ(k,I_LGdep) = Gig*rhoq(k,I_NG)*dq_xave(k,I_mp_QG)*ventLG(k)
       PQ(k,I_NRdep) = PQ(k,I_LRdep)/xq(k,I_mp_QR)
       PQ(k,I_NIdep) = 0.0_RP
       PQ(k,I_NSdep) = PQ(k,I_LSdep)/xq(k,I_mp_QS)
       PQ(k,I_NGdep) = PQ(k,I_LGdep)/xq(k,I_mp_QG)
    end do

    do k = KS, KE
       temc    = tem(k) - T00   ! degC
       !------------------------------------------------------------------------
       ! Melting part is described by Pruppacher and Klett (1997) Sec.16.3.1
       ! Here we omit "Shedding" of snow-flakes and ice-particles.
       ! "Shedding" may be applicative if you refer
       ! eq.(38) in Cotton etal.(1986) Jour. Clim. Appl. Meteor. p.1658-1680.
       ! SB06(73)
       Dt      = Kalfa(k)/(CPvap*rho_0)
       ! Gm: factor caused by balance between
       !     "water evaporation cooling(1st.)" and "fusion heating(2nd.)"
       ! SB06(76)
       ! [fix] 08/05/08 T.Mitsui  LHF00 => EMELT  and  esw => PSAT0
       ! LHS0 is more suitable than LHS because melting occurs around 273.15 K.
       Gm      = 2.0_RP*PI/EMELT&
               * ( (Kalfa(k)*Dt/Dw)*(temc) + (Dw*LHS0/Rvap)*(esi(k)/tem(k)-PSAT0/T00) )
       ! SB06(76)
       ! Notice! melting only occurs where T > 273.15 K else doesn't.
       ! [fix] 08/05/08 T.Mitsui, Gm could be both positive and negative value.
       !       See Pruppacher and Klett(1997) eq.(16-79) or Rasmussen and Pruppacher(1982)
       sw = ( sign(0.5_RP,temc) + 0.5_RP ) * ( sign(0.5_RP,Gm-eps) + 0.5_RP ) ! sw = 1 if( (temc>=0.0_RP) .AND. (Gm>0.0_RP) ), otherwise sw = 0
       !  if Gm==0 then rh and tem is critical value for melting process.
       ! 08/05/16 [Mod] T.Mitsui, change term of PLimlt. N_i => L_i/ (limited x_i)
       ! because melting never occur when N_i=0.
       PQ(k,I_LImlt) = - Gm * rhoq(k,I_QI)*dq_xave(k,I_mp_QI)*ventLI(k)/xq(k,I_mp_QI) * sw
       PQ(k,I_NImlt) = - Gm * rhoq(k,I_NI)*dq_xave(k,I_mp_QI)*ventNI(k)/xq(k,I_mp_QI) * sw
       PQ(k,I_LSmlt) = - Gm * rhoq(k,I_QS)*dq_xave(k,I_mp_QS)*ventLS(k)/xq(k,I_mp_QS) * sw
       PQ(k,I_NSmlt) = - Gm * rhoq(k,I_NS)*dq_xave(k,I_mp_QS)*ventNS(k)/xq(k,I_mp_QS) * sw
       PQ(k,I_LGmlt) = - Gm * rhoq(k,I_QG)*dq_xave(k,I_mp_QG)*ventLG(k)/xq(k,I_mp_QG) * sw
       PQ(k,I_NGmlt) = - Gm * rhoq(k,I_NG)*dq_xave(k,I_mp_QG)*ventNG(k)/xq(k,I_mp_QG) * sw
    end do
    !
    return
  end subroutine dep_vapor_melt_ice
  !-----------------------------------------------------------------------------
!OCL SERIAL
  subroutine freezing_water( &
       KA, KS, KE, &
       dt,            &
       rhoq, xq, tem, &
       PQ             )
    implicit none
    !
    ! In this subroutine,
    ! We assumed surface temperature of droplets are same as environment.

    integer, intent(in) :: KA, KS, KE

    real(RP), intent(in) :: dt
    !
    real(RP), intent(in) :: tem(KA)
    !
    real(RP), intent(in) :: rhoq(KA,I_QV:I_NG)
    real(RP), intent(in) :: xq(KA,HYDRO_MAX)
    !
    real(RP), intent(inout):: PQ(KA,PQ_MAX)
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
    real(RP) :: Jhom, Jhet, Jh(KA)
    real(RP) :: rdt
    real(RP) :: tmp
    !
    integer :: k
    !
    rdt = 1.0_RP/dt
    !
    coef_m2_c =   coef_m2(I_mp_QC)
    coef_m2_r =   coef_m2(I_mp_QR)
    !

    ! Note, xc should be limited in range[xc_min:xc_max].
    ! and PNChom need to be calculated by NC
    ! because reduction rate of Nc need to be bound by NC.
    ! For the same reason PLChom also bound by LC and xc.
    ! Basically L and N should be independent
    ! but average particle mass x should be in suitable range.

    ! Homogenous Freezing
    do k = KS, KE
       PQ(k,I_LChom) = 0.0_RP
       PQ(k,I_NChom) = 0.0_RP
    end do

    ! Heterogenous Freezing
    do k = KS, KE
       temc = max( tem(k) - T00, temc_min )
       ! These cause from aerosol-droplet interaction.
       ! Bigg(1953) formula, Khain etal.(2000) eq.(4.5), Pruppacher and Klett(1997) eq.(9-48)
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
          Jhet =  a_het*exp( -b_het*temc - 1.0_RP )
       else if( temc < 0.0_RP) then
          jhom = 10._RP**(-7.63_RP-2.996_RP*(temc+30.0_RP))*1.E+3_RP
          Jhet =  a_het*exp( -b_het*temc - 1.0_RP )
       else
          Jhom = 0.0_RP
          Jhet = 0.0_RP
       end if
       Jh(k) = ( Jhet + Jhom ) * dt
    end do
    do k = KS, KE
#if defined(PGI) || defined(SX)
       tmp = min( xq(k,I_mp_QC)*Jh(k), 1.E+3_RP) ! apply exp limiter
       PQ(k,I_LChet) = -rdt*rhoq(k,I_QC)*( 1.0_RP - exp( -coef_m2_c*tmp ) )
       PQ(k,I_NChet) = -rdt*rhoq(k,I_NC)*( 1.0_RP - exp( -          tmp ) )

       tmp = min( xq(k,I_mp_QR)*Jh(k), 1.E+3_RP) ! apply exp limiter
       PQ(k,I_LRhet) = -rdt*rhoq(k,I_QR)*( 1.0_RP - exp( -coef_m2_r*tmp ) )
       PQ(k,I_NRhet) = -rdt*rhoq(k,I_NR)*( 1.0_RP - exp( -          tmp ) )
#else
       PQ(k,I_LChet) = -rdt*rhoq(k,I_QC)*( 1.0_RP - exp( -coef_m2_c*xq(k,I_mp_QC)*Jh(k) ) )
       PQ(k,I_NChet) = -rdt*rhoq(k,I_NC)*( 1.0_RP - exp( -          xq(k,I_mp_QC)*Jh(k) ) )
       PQ(k,I_LRhet) = -rdt*rhoq(k,I_QR)*( 1.0_RP - exp( -coef_m2_r*xq(k,I_mp_QR)*Jh(k) ) )
       PQ(k,I_NRhet) = -rdt*rhoq(k,I_NR)*( 1.0_RP - exp( -          xq(k,I_mp_QR)*Jh(k) ) )
#endif
    end do

    !
    return
  end subroutine freezing_water

  !----------------------------------------------------------------
!OCL SERIAL
  subroutine update_by_phase_change(   &
       KA, KS, KE, &
       ntmax,                & ! in
       dt,                   & ! in
       cz,                   & ! in
       fz,                   & ! in
       w,                    & ! in
       dTdt_rad,             & ! in
       rho,                  & ! in
       qdry,                 & ! in
       esw, esi, rhoq2,      & ! in
       pre, tem,             & ! in
       cpa,cva,              & ! in
       flg_lt,               & ! in
       PQ,                   & ! inout
       sl_PLCdep,            & ! inout
       sl_PLRdep, sl_PNRdep, & ! inout
       RHOQ_t,               & ! out
       RHOE_t,               & ! out
       CPtot_t,              & ! out
       CVtot_t,              & ! out
       qc_evaporate,         & ! out
       rhoq2_crg,            & ! in:optional
       RHOQcrg_t             ) ! out:optional

    use scale_atmos_hydrometeor, only: &
       CP_VAPOR, &
       CP_WATER, &
       CP_ICE,   &
       CV_VAPOR, &
       CV_WATER, &
       CV_ICE
    use scale_atmos_saturation, only: &
       moist_pres2qsat_liq  => ATMOS_SATURATION_pres2qsat_liq,  &
       moist_pres2qsat_ice  => ATMOS_SATURATION_pres2qsat_ice,  &
       moist_dqs_dtem_dens_liq  => ATMOS_SATURATION_dqs_dtem_dens_liq,  &
       moist_dqs_dtem_dens_ice  => ATMOS_SATURATION_dqs_dtem_dens_ice,  &
       moist_dqs_dtem_dpre_liq => ATMOS_SATURATION_dqs_dtem_dpre_liq, &
       moist_dqs_dtem_dpre_ice => ATMOS_SATURATION_dqs_dtem_dpre_ice
    implicit none

    integer, intent(in) :: KA, KS, KE

    integer, intent(in)    :: ntmax
    !
    real(RP), intent(in)   :: dt           ! time step[s]
    real(RP), intent(in)   :: cz(KA)       ! altitude [m]
    real(RP), intent(in)   :: fz(0:KA)       ! altitude difference [m]
    real(RP), intent(in)   :: w (KA)       ! vertical velocity @ full level [m/s]
    real(RP), intent(in)   :: dTdt_rad(KA) ! temperture tendency by radiation[K/s]
    real(RP), intent(in)   :: rho (KA)     ! density[kg/m3]
    real(RP), intent(in)   :: qdry(KA)     ! dry air mass ratio [kg/kg]
    real(RP), intent(in)   :: esw (KA)     ! saturated vapor pressure for liquid
    real(RP), intent(in)   :: esi (KA)     !                          for ice
    real(RP), intent(in)   :: rhoq2(KA,I_QV:I_NG)

    real(RP), intent(in)   :: tem(KA)      ! temperature[K]
    real(RP), intent(in)   :: pre(KA)      ! pressure[Pa]
    real(RP), intent(in)   :: cpa(KA)      !
    real(RP), intent(in)   :: cva(KA)      ! specific heat at constant volume

    !+++ tendency[kg/m3/s]
    real(RP), intent(inout) :: PQ(KA,PQ_MAX)
    !+++ Column integrated tendency[kg/m2/s]
    real(RP), intent(inout) :: sl_PLCdep
    real(RP), intent(inout) :: sl_PLRdep, sl_PNRdep

    real(RP),intent(out) :: RHOQ_t(KA,QA_MP)
    real(RP),intent(out) :: RHOE_t(KA)
    real(RP),intent(out) :: CPtot_t(KA)
    real(RP),intent(out) :: CVtot_t(KA)

    !+++ tendency[kg/m3/s]
    real(RP), intent(out)   :: qc_evaporate(KA)

    !--- for lightning component
    logical, intent(in)   :: flg_lt ! false -> without lightning, true-> with lightning
    real(RP), intent(in), optional  :: rhoq2_crg(KA,I_QC:I_QG)
    real(RP), intent(inout), optional :: RHOQcrg_t(KA,HYDRO_MAX)

    real(RP) :: xi               ! mean mass of ice particles
    real(RP) :: rrho             ! 1/rho
    real(RP) :: wtem(KA)         ! temperature[K]
    !
    real(RP) :: r_cva            ! specific heat at constant volume
    !real(RP) :: cpa             ! specific heat at constant pressure
    real(RP) :: r_cpa            ! specific heat at constant pressure
    real(RP) :: qsw(KA)          ! saturated mixing ratio for liquid
    real(RP) :: qsi(KA)          ! saturated mixing ratio for solid
    real(RP) :: dqswdtem_rho(KA) ! (dqsw/dtem)_rho
    real(RP) :: dqsidtem_rho(KA) ! (dqsi/dtem)_rho
    real(RP) :: dqswdtem_pre(KA) ! (dqsw/dtem)_pre
    real(RP) :: dqsidtem_pre(KA) ! (dqsi/dtem)_pre
    real(RP) :: dqswdpre_tem(KA) ! (dqsw/dpre)_tem
    real(RP) :: dqsidpre_tem(KA) ! (dqsi/dpre)_tem
    !
    real(RP) :: w2(KA)                 ! vetical velocity[m/s]
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
    real(RP) :: taucnd_c(KA), r_taucnd_c  ! by cloud
    real(RP) :: taucnd_r(KA), r_taucnd_r  ! by rain
    real(RP) :: taudep_i(KA), r_taudep_i  ! by ice
    real(RP) :: taudep_s(KA), r_taudep_s  ! by snow
    real(RP) :: taudep_g(KA), r_taudep_g  ! by graupel
    ! alternative tendency through changing ssw and ssi
    real(RP) :: PNCdep
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
    real(RP) :: dep_dqi(KA)
    real(RP) :: dep_dni(KA)
    real(RP) :: dep_dqs(KA)
    real(RP) :: dep_dns(KA)
    real(RP) :: dep_dqg(KA)
    real(RP) :: dep_dng(KA)
    real(RP) :: dep_dqr(KA)
    real(RP) :: dep_dnr(KA)
    real(RP) :: dep_dqc(KA)
    real(RP) :: dep_dnc(KA)   ! 11/08/30 [Add] T.Mitsui, dep_dnc
    real(RP) :: r_xc_ccn, r_xi_ccn ! 11/08/30 [Add] T.Mitsui
    !
    real(RP) :: drhoqv(KA)
    real(RP) :: drhoqc(KA), drhoqr(KA), drhoqi(KA), drhoqs(KA), drhoqg(KA)
    real(RP) :: drhonc(KA), drhonr(KA), drhoni(KA), drhons(KA), drhong(KA)
    !-- for Charge densicty
    real(RP) :: drhoqcrg_c(KA), drhoqcrg_r(KA)
    real(RP) :: drhoqcrg_i(KA), drhoqcrg_s(KA), drhoqcrg_g(KA)
    real(RP) :: frz_dnc_crg
    real(RP) :: frz_dnr_crg
    real(RP) :: mlt_dni_crg
    real(RP) :: mlt_dns_crg
    real(RP) :: mlt_dng_crg
    real(RP) :: dep_dni_crg
    real(RP) :: dep_dns_crg
    real(RP) :: dep_dng_crg
    real(RP) :: dep_dnr_crg
    real(RP) :: dep_dnc_crg
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

    real(RP) :: fac_cndc_wrk
    !
    real(RP), parameter :: tau100day   = 1.E+7_RP
    real(RP), parameter :: r_tau100day = 1.E-7_RP
    real(RP), parameter :: eps=1.E-30_RP
    !
    real(RP) :: PLCdep(KA), PLRdep(KA), PNRdep(KA)
    real(RP) :: dz
    !
    integer :: k,iqw
    real(RP) :: sw, sw2
    real(RP) :: dqv, dql, dqi
    real(RP) :: dcv, dcp
    real(RP) :: dqc_crg, dqr_crg, dqi_crg, dqs_crg, dqg_crg

    !
    real(RP) :: fact

    !
!    dt_dyn     = dt*ntmax
    !
    r_dt       = 1.0_RP/dt
    !
    r_xc_ccn=1.0_RP/xc_ccn
!    r_xi_ccn=1.0_RP/xi_ccn
    !
    if( opt_fix_taucnd_c )then
       fac_cndc_wrk = fac_cndc**(1.0_RP-b_m(I_mp_QC))
       do k = KS, KE
          PQ(k,I_LCdep)  = PQ(k,I_LCdep)*fac_cndc_wrk
       end do
       LOG_INFO("ATMOS_PHY_MP_SN14_update_by_phase_change",*) "taucnd:fac_cndc_wrk=",fac_cndc_wrk
    end if

!OCL XFILL
    do k = KS, KE
       ! Temperature lower limit is only used for saturation condition.
       ! On the other hand original "tem" is used for calculation of latent heat or energy equation.
       wtem(k)  = max( tem(k), tem_min )
    end do

    call moist_pres2qsat_liq( KA, KS, KE, &
                              wtem(:), pre(:), qdry(:), & ! [IN]
                              qsw(:)                    ) ! [OUT]
    call moist_pres2qsat_ice( KA, KS, KE, &
                              wtem(:), pre(:), qdry(:), & ! [IN]
                              qsi(:)                    ) ! [OUT]
    call moist_dqs_dtem_dens_liq( KA, KS, KE, &
                                  wtem(:), rho(:), & ! [IN]
                                  dqswdtem_rho(:)  ) ! [OUT]
    call moist_dqs_dtem_dens_ice( KA, KS, KE, &
                                  wtem(:), rho(:), & ! [IN]
                                  dqsidtem_rho(:)  ) ! [OUT]
    call moist_dqs_dtem_dpre_liq( KA, KS, KE, &
                                  wtem(:), pre(:), qdry(:),        & ! [IN]
                                  dqswdtem_pre(:), dqswdpre_tem(:) ) ! [OUT]
    call moist_dqs_dtem_dpre_ice( KA, KS, KE, &
                                  wtem(:), pre(:), qdry(:),        & ! [IN]
                                  dqsidtem_pre(:), dqsidpre_tem(:) ) ! [OUT]

    do k = KS, KE
       if( cz(k) <= 25000.0_RP )then
          w2(k) = w(k)
       else
          w2(k) = 0.0_RP
       end if
       if( pre(k) < esw(k)+1.E-10_RP )then
          qsw(k) = 1.0_RP
          dqswdtem_rho(k) = 0.0_RP
          dqswdtem_pre(k) = 0.0_RP
          dqswdpre_tem(k) = 0.0_RP
       end if
       if( pre(k) < esi(k)+1.E-10_RP )then
          qsi(k) = 1.0_RP
          dqsidtem_rho(k) = 0.0_RP
          dqsidtem_pre(k) = 0.0_RP
          dqsidpre_tem(k) = 0.0_RP
       end if
    end do

!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       r_rvaptem        = 1.0_RP/(Rvap*wtem(k))
       lvsw             = esw(k)*r_rvaptem        ! rho=p/(Rv*T)
       lvsi             = esi(k)*r_rvaptem        !
       pv               = rhoq2(k,I_QV)*Rvap*tem(k)
       r_esw            = 1.0_RP/esw(k)
       r_esi            = 1.0_RP/esi(k)
       ssw              = min( MP_ssw_lim, ( pv*r_esw-1.0_RP ) )
       ssi              = pv*r_esi - 1.0_RP
       r_lvsw           = 1.0_RP/lvsw
       r_lvsi           = 1.0_RP/lvsi
       r_taucnd_c       = PQ(k,I_LCdep)*r_lvsw
       r_taucnd_r       = PQ(k,I_LRdep)*r_lvsw
       r_taudep_i       = PQ(k,I_LIdep)*r_lvsi
       r_taudep_s       = PQ(k,I_LSdep)*r_lvsi
       r_taudep_g       = PQ(k,I_LGdep)*r_lvsi
!       taucnd_c(k)   = 1.0_RP/(r_taucnd_c+r_tau100day)
!       taucnd_r(k)   = 1.0_RP/(r_taucnd_r+r_tau100day)
!       taudep_i(k)   = 1.0_RP/(r_taudep_i+r_tau100day)
!       taudep_s(k)   = 1.0_RP/(r_taudep_s+r_tau100day)
!       taudep_g(k)   = 1.0_RP/(r_taudep_g+r_tau100day)

       r_cva = 1.0_RP / cva(k)
       r_cpa = 1.0_RP / cpa(k)

       ! Coefficient of latent heat release for ssw change by PLCdep and PLRdep
       aliqliq = 1.0_RP &
               + r_cva*( LHV00              + (CVvap-CL)*tem(k) )*dqswdtem_rho(k)
       ! Coefficient of latent heat release for ssw change by PLIdep, PLSdep and PLGdep
       asolliq = 1.0_RP &
               + r_cva*( LHV00 + LHF00 + (CVvap-CI)*tem(k) )*dqswdtem_rho(k)
       ! Coefficient of latent heat release for ssi change by PLCdep and PLRdep
       aliqsol = 1.0_RP &
               + r_cva*( LHV00              + (CVvap-CL)*tem(k) )*dqsidtem_rho(k)
       ! Coefficient of latent heat release for ssi change by PLIdep, PLSdep and PLGdep
       asolsol = 1.0_RP &
               + r_cva*( LHV00 + LHF00 + (CVvap-CI)*tem(k) )*dqsidtem_rho(k)
       Pdynliq = w2(k) * GRAV * ( r_cpa*dqswdtem_pre(k) + rho(k)*dqswdpre_tem(k) )
       Pdynsol = w2(k) * GRAV * ( r_cpa*dqsidtem_pre(k) + rho(k)*dqsidpre_tem(k) )
       Pradliq = -dTdt_rad(k)    * dqswdtem_rho(k)
       Pradsol = -dTdt_rad(k)    * dqsidtem_rho(k)

       ssw_o   = ssw
       ssi_o   = ssi
!       ssw_o   = ssw - Pdynliq*(dt_dyn-dt_mp)/qsw(k) + Pradliq*r_qsw*dt_mp
!       ssi_o   = ssi - Pdynsol*(dt_dyn-dt_mp)/qsi(k) + Pradsol*r_qsi*dt_mp

       r_taucnd         = &
               + aliqliq*( r_taucnd_c+r_taucnd_r ) &
               + asolliq*( r_taudep_i+r_taudep_s+r_taudep_g )
       r_taudep         = &
               + aliqsol*( r_taucnd_c+r_taucnd_r )&
               + asolsol*( r_taudep_i+r_taudep_s+r_taudep_g )

       if( r_taucnd < r_tau100day )then
          uplim_cnd        = max( rho(k)*ssw_o*qsw(k)*r_dt, 0.0_RP )
          lowlim_cnd       = min( rho(k)*ssw_o*qsw(k)*r_dt, 0.0_RP )
!          taucnd            = tau100day
          PQ(k,I_LCdep) = max(lowlim_cnd, min(uplim_cnd, PQ(k,I_LCdep)*ssw_o ))
          PQ(k,I_LRdep) = max(lowlim_cnd, min(uplim_cnd, PQ(k,I_LRdep)*ssw_o ))
          PQ(k,I_NRdep) = min(0.0_RP, PQ(k,I_NRdep)*ssw_o )
!          PLR2NR           = 0.0_RP
       else
          Acnd = Pdynliq + Pradliq &
               - ( r_taudep_i+r_taudep_s+r_taudep_g ) * ( qsw(k) - qsi(k) )

          taucnd     = 1.0_RP/r_taucnd
          ! Production term for liquid water content
          coef_a_cnd = rho(k)*Acnd*taucnd
          coef_b_cnd = rho(k)*taucnd*r_dt*(ssw_o*qsw(k)-Acnd*taucnd) * ( exp(-dt*r_taucnd) - 1.0_RP )
          PQ(k,I_LCdep) = coef_a_cnd*r_taucnd_c - coef_b_cnd*r_taucnd_c
          PLR2NR            = PQ(k,I_NRdep)/(PQ(k,I_LRdep)+1.E-30_RP)
          PQ(k,I_LRdep) = coef_a_cnd*r_taucnd_r - coef_b_cnd*r_taucnd_r
          PQ(k,I_NRdep) = min(0.0_RP, PQ(k,I_LRdep)*PLR2NR )
       end if

       if( r_taudep < r_tau100day )then
          uplim_dep        = max( rho(k)*ssi_o*qsi(k)*r_dt, 0.0_RP )
          lowlim_dep       = min( rho(k)*ssi_o*qsi(k)*r_dt, 0.0_RP )
!          taudep            = tau100day
          PQ(k,I_LIdep) = max(lowlim_dep, min(uplim_dep, PQ(k,I_LIdep)*ssi_o ))
          PQ(k,I_LSdep) = max(lowlim_dep, min(uplim_dep, PQ(k,I_LSdep)*ssi_o ))
          PQ(k,I_LGdep) = max(lowlim_dep, min(uplim_dep, PQ(k,I_LGdep)*ssi_o ))
          PQ(k,I_NIdep) = min(0.0_RP, PQ(k,I_NIdep)*ssi_o )
          PQ(k,I_NSdep) = min(0.0_RP, PQ(k,I_NSdep)*ssi_o )
          PQ(k,I_NGdep) = min(0.0_RP, PQ(k,I_NGdep)*ssi_o )
       else
          Adep = Pdynsol + Pradsol &
               + ( r_taucnd_c+r_taucnd_r )            * ( qsw(k) - qsi(k) )

          taudep     = 1.0_RP/r_taudep
          ! Production term for ice water content
          coef_a_dep = rho(k)*Adep*taudep
          coef_b_dep = rho(k)*taudep*r_dt*(ssi_o*qsi(k)-Adep*taudep) * ( exp(-dt*r_taudep) - 1.0_RP )
          PLI2NI           = PQ(k,I_NIdep)/max(PQ(k,I_LIdep),1.E-30_RP)
          PLS2NS           = PQ(k,I_NSdep)/max(PQ(k,I_LSdep),1.E-30_RP)
          PLG2NG           = PQ(k,I_NGdep)/max(PQ(k,I_LGdep),1.E-30_RP)
          PQ(k,I_LIdep) =  coef_a_dep*r_taudep_i - coef_b_dep*r_taudep_i
          PQ(k,I_LSdep) =  coef_a_dep*r_taudep_s - coef_b_dep*r_taudep_s
          PQ(k,I_LGdep) =  coef_a_dep*r_taudep_g - coef_b_dep*r_taudep_g
          PQ(k,I_NIdep) = min(0.0_RP, PQ(k,I_LIdep)*PLI2NI )
          PQ(k,I_NSdep) = min(0.0_RP, PQ(k,I_LSdep)*PLS2NS )
          PQ(k,I_NGdep) = min(0.0_RP, PQ(k,I_LGdep)*PLG2NG )
       end if
    end do

    !--- evaporation/condensation
!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       sw = 0.5_RP - sign(0.5_RP, PQ(k,I_LCdep)+eps) != 1 for PLCdep<=-eps
       PNCdep = min(0.0_RP, ((rhoq2(k,I_QC)+PQ(k,I_LCdep)*dt)*r_xc_ccn - rhoq2(k,I_NC))*r_dt ) * sw
!       if( PQ(k,I_LCdep) < -eps )then
!          PNCdep = min(0.0_RP, ((rhoq2(k,I_QC)+PQ(k,I_LCdep)*dt)*r_xc_ccn - rhoq2(k,I_NC))*r_dt )
!       else
!          PNCdep = 0.0_RP
!       end if
!       if( PQ(k,I_LIdep) < -eps )then
!          PQ(k,I_NIdep) = min(0.0_RP, ((li(k)+PQ(k,I_LIdep)*dt)*r_xi_ccn - rhoq2(k,I_NI))*r_dt )
!       else
!          PQ(k,I_NIdep) = 0.0_RP
!       end if

       r_rvaptem = 1.0_RP/(Rvap*wtem(k))
       lvsw    = esw(k)*r_rvaptem
       dlvsw   = rhoq2(k,I_QV)-lvsw
       dcnd    = dt*(PQ(k,I_LCdep)+PQ(k,I_LRdep))

       sw = ( sign(0.5_RP,dcnd) + sign(0.5_RP,dlvsw) ) &
          * ( 0.5_RP + sign(0.5_RP,abs(dcnd)-eps) ) ! to avoid zero division
       ! sw= 1: always supersaturated
       ! sw=-1: always unsaturated
       ! sw= 0: partially unsaturated during timestep
       fac1 = min(dlvsw*sw,dcnd*sw)*sw / (abs(sw)-1.0_RP+dcnd) & ! sw=1,-1
            + 1.0_RP - abs(sw)                                   ! sw=0
       dep_dqc(k) = max( dt*PQ(k,I_LCdep)*fac1, &
                        -rhoq2(k,I_QC) - 1e30_RP*(sw+1.0_RP) )*abs(sw) != -lc for sw=-1, -inf for sw=1
       dep_dqr(k) = max( dt*PQ(k,I_LRdep)*fac1, &
                        -rhoq2(k,I_QR) - 1e30_RP*(sw+1.0_RP) )*abs(sw) != -lr for sw=-1, -inf for sw=1
!       if     ( (dcnd >  eps) .AND. (dlvsw > eps) )then
!          ! always supersaturated
!          fac1    = min(dlvsw,dcnd)/dcnd
!          dep_dqc(k) =  dt*PQ(k,I_LCdep)*fac1
!          dep_dqr(k) =  dt*PQ(k,I_LRdep)*fac1
!       else if( (dcnd < -eps) .AND. (dlvsw < -eps) )then
!          ! always unsaturated
!          fac1    = max( dlvsw,dcnd )/dcnd
!          dep_dqc(k) = max( dt*PQ(k,I_LCdep)*fac1, -rhoq2(k,I_QC) )
!          dep_dqr(k) = max( dt*PQ(k,I_LRdep)*fac1, -rhoq2(k,I_QR) )
!       else
!          ! partially unsaturated during timestep
!          fac1    = 1.0_RP
!          dep_dqc(k) = 0.0_RP
!          dep_dqr(k) = 0.0_RP
!       end if

       ! evaporation always lose number(always negative).
       dep_dnc(k) = max( dt*PNCdep*fac1, -rhoq2(k,I_NC) ) ! ss>0 dep=0, ss<0 dep<0 ! [Add] 11/08/30 T.Mitsui
       dep_dnr(k) = max( dt*PQ(k,I_NRdep)*fac1, -rhoq2(k,I_NR) ) ! ss>0 dep=0, ss<0 dep<0

       qc_evaporate(k) = - dep_dnc(k) ! [Add] Y.Sato 15/09/08

       !--- reduce charge density of cloud and rain by evaporation
       if (flg_lt) then
          sw = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NC)-SMALL ) !--- if NC is small,  ignore charge transfer
          dep_dnc_crg = dep_dnc(k)*( 1.0_RP-sw )/( rhoq2(k,I_NC)+sw )*rhoq2_crg(k,I_QC)
          sw = min( abs(rhoq2_crg(k,I_QC)),abs(dep_dnc_crg) )
          dep_dnc_crg = sign( sw,dep_dnc_crg )

          drhoqcrg_c(k) = dep_dnc_crg


          sw = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NR)-SMALL ) !--- if NR is small,  ignore charge transfer
          dep_dnr_crg = dep_dnr(k) * ( 1.0_RP-sw )/( rhoq2(k,I_NR)+sw )*rhoq2_crg(k,I_QR)
          sw = min( abs(rhoq2_crg(k,I_QR)),abs(dep_dnr_crg) )
          dep_dnr_crg = sign( sw,dep_dnr_crg )

          drhoqcrg_r(k) = dep_dnr_crg
       end if
    end do

       !--- deposition/sublimation
!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       lvsi    = esi(k)*r_rvaptem
       ddep    = dt*(PQ(k,I_LIdep)+PQ(k,I_LSdep)+PQ(k,I_LGdep))
       dlvsi   = rhoq2(k,I_QV)-lvsi  ! limiter for esi>1.d0

       sw = ( sign(0.5_RP,ddep) + sign(0.5_RP,dlvsi) ) &
          * ( 0.5_RP + sign(0.5_RP,abs(ddep)-eps) ) ! to avoid zero division
       ! sw= 1: always supersaturated
       ! sw=-1: always unsaturated
       ! sw= 0: partially unsaturated during timestep
       fac2 = min(dlvsi*sw,ddep*sw)*sw / (abs(sw)-1.0_RP+ddep) & ! sw=1,-1
            + 1.0_RP - abs(sw)                                   ! sw=0
       dep_dqi(k) = max( dt*PQ(k,I_LIdep) &
                         * ( 1.0_RP-abs(sw) + fac2*abs(sw) ), & != fac2 for sw=-1,1, 1 for sw=0
                        -rhoq2(k,I_QI) - 1e30_RP*(sw+1.0_RP) )  != -li for sw=-1, -inf for sw=0,1
       dep_dqs(k) = max( dt*PQ(k,I_LSdep) &
                         * ( 1.0_RP-abs(sw) + fac2*abs(sw) ), & != fac2 for sw=-1,1, 1 for sw=0
                        -rhoq2(k,I_QS) - 1e30_RP*(sw+1.0_RP) )  != -ls for sw=-1, -inf for sw=0,1
       dep_dqg(k) = max( dt*PQ(k,I_LGdep) &
                         * ( 1.0_RP-abs(sw) + fac2*abs(sw) ), & != fac2 for sw=-1,1, 1 for sw=0
                        -rhoq2(k,I_QG) - 1e30_RP*(sw+1.0_RP) )  != -lg for sw=-1, -inf for sw=0,1
!       if      ( (ddep >  eps) .AND. (dlvsi > eps) )then
!          ! always supersaturated
!          fac2    = min(dlvsi,ddep)/ddep
!          dep_dqi(k) = dt*PQ(k,I_LIdep)*fac2
!          dep_dqs(k) = dt*PQ(k,I_LSdep)*fac2
!          dep_dqg(k) = dt*PQ(k,I_LGdep)*fac2
!       else if ( (ddep < -eps) .AND. (dlvsi < -eps) )then
!          ! always unsaturated
!          fac2    = max(dlvsi,ddep)/ddep
!          dep_dqi(k) = max(dt*PQ(k,I_LIdep)*fac2, -rhoq2(k,I_QI) )
!          dep_dqs(k) = max(dt*PQ(k,I_LSdep)*fac2, -rhoq2(k,I_QS) )
!          dep_dqg(k) = max(dt*PQ(k,I_LGdep)*fac2, -rhoq2(k,I_QG) )
!       else
!          ! partially unsaturated during timestep
!          fac2    = 1.0_RP
!          dep_dqi(k) = dt*PQ(k,I_LIdep)
!          dep_dqs(k) = dt*PQ(k,I_LSdep)
!          dep_dqg(k) = dt*PQ(k,I_LGdep)
!       end if

       ! evaporation always lose number(always negative).
       dep_dni(k) = max( dt*PQ(k,I_NIdep)*fac2, -rhoq2(k,I_NI) ) ! ss>0 dep=0, ss<0 dep<0
       dep_dns(k) = max( dt*PQ(k,I_NSdep)*fac2, -rhoq2(k,I_NS) ) ! ss>0 dep=0, ss<0 dep<0
       dep_dng(k) = max( dt*PQ(k,I_NGdep)*fac2, -rhoq2(k,I_NG) ) ! ss>0 dep=0, ss<0 dep<0

       !--- for Charge density
       if (flg_lt) then
          sw = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NI)-SMALL ) !--- if NI is small,  ignore charge transfer
          dep_dni_crg = dep_dni(k)*( 1.0_RP-sw )/( rhoq2(k,I_NI)+sw ) * rhoq2_crg(k,I_QI)
          sw = min( abs(rhoq2_crg(k,I_QI)),abs(dep_dni_crg) )
          dep_dni_crg = sign( sw,dep_dni_crg )

          drhoqcrg_i(k) = dep_dni_crg


          sw = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NS)-SMALL ) !--- if NS is small,  ignore charge transfer
          dep_dns_crg = dep_dns(k)*( 1.0_RP-sw )/( rhoq2(k,I_NS)+sw ) * rhoq2_crg(k,I_QS)
          sw = min( abs(rhoq2_crg(k,I_QS)),abs(dep_dns_crg) )
          dep_dns_crg = sign( sw,dep_dns_crg )

          drhoqcrg_s(k) = dep_dns_crg


          sw = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NG)-SMALL ) !--- if NG is small,  ignore charge transfer
          dep_dng_crg = dep_dng(k)*( 1.0_RP-sw )/( rhoq2(k,I_NG)+sw ) * rhoq2_crg(k,I_QG)
          sw = min( abs(rhoq2_crg(k,I_QG)),abs(dep_dng_crg) )
          dep_dng_crg = sign( sw,dep_dng_crg )

          drhoqcrg_g(k) = dep_dng_crg

       end if
    end do

    !--- freezing of cloud drop
!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       frz_dqc = max( dt*(PQ(k,I_LChom)+PQ(k,I_LChet)), -rhoq2(k,I_QC)-dep_dqc(k) ) ! negative value
       frz_dnc = max( dt*(PQ(k,I_NChom)+PQ(k,I_NChet)), -rhoq2(k,I_NC)-dep_dnc(k) ) ! negative value

       drhoqc(k) =   frz_dqc
       drhonc(k) =   frz_dnc
       drhoqi(k) = - frz_dqc
       drhoni(k) = - frz_dnc

       fac3    = ( frz_dqc-eps )/( dt*(PQ(k,I_LChom)+PQ(k,I_LChet))-eps )
       fac4    = ( frz_dnc-eps )/( dt*(PQ(k,I_NChom)+PQ(k,I_NChet))-eps )
       PQ(k,I_LChom) = fac3*PQ(k,I_LChom)
       PQ(k,I_LChet) = fac3*PQ(k,I_LChet)
       PQ(k,I_NChom) = fac4*PQ(k,I_NChom)
       PQ(k,I_NChet) = fac4*PQ(k,I_NChet)
       if (flg_lt) then
          sw = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NC)-SMALL ) !--- if NI is small,  ignore charge transfer
          frz_dnc_crg = frz_dnc*( 1.0_RP-sw )/( rhoq2(k,I_NC)+sw ) * rhoq2_crg(k,I_QC)
          !--- limiter
          sw = min( abs(rhoq2_crg(k,I_QC)+dep_dnc_crg),abs(frz_dnc_crg) )
          frz_dnc_crg = sign( sw,frz_dnc_crg )

          drhoqcrg_c(k) = drhoqcrg_c(k) + frz_dnc_crg
          drhoqcrg_i(k) = drhoqcrg_i(k) - frz_dnc_crg
       end if
    end do

    !--- melting
!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       ! ice change
       mlt_dqi = max( dt*PQ(k,I_LImlt), -rhoq2(k,I_QI)-dep_dqi(k) )  ! negative value
       mlt_dni = max( dt*PQ(k,I_NImlt), -rhoq2(k,I_NI)-dep_dni(k) )  ! negative value

       xi = min(xi_max, max(xi_min, rhoq2(k,I_QI)/(rhoq2(k,I_NI)+ni_min) ))
       sw = 0.5_RP + sign(0.5_RP,xi-x_sep) ! if (xi>=x_sep) then sw=1 else sw=0
                                           ! sw=1: large ice crystals turn into rain by melting

       ! snow change
       mlt_dqs = max( dt*PQ(k,I_LSmlt), -rhoq2(k,I_QS)-dep_dqs(k) )  ! negative value
       mlt_dns = max( dt*PQ(k,I_NSmlt), -rhoq2(k,I_NS)-dep_dns(k) )  ! negative value

       ! graupel change
       mlt_dqg = max( dt*PQ(k,I_LGmlt), -rhoq2(k,I_QG)-dep_dqg(k) )  ! negative value
       mlt_dng = max( dt*PQ(k,I_NGmlt), -rhoq2(k,I_NG)-dep_dng(k) )  ! negative value

       drhoqc(k) = drhoqc(k) - mlt_dqi * (1.0_RP-sw)
       drhonc(k) = drhonc(k) - mlt_dni * (1.0_RP-sw)

       drhoqr(k) =           - mlt_dqi * sw          - mlt_dqs - mlt_dqg
       drhonr(k) =           - mlt_dni * sw          - mlt_dns - mlt_dng

       drhoqi(k) = drhoqi(k) + mlt_dqi
       drhoni(k) = drhoni(k) + mlt_dni

       drhoqs(k) =                                     mlt_dqs
       drhons(k) =                                     mlt_dns

       drhoqg(k) =                                               mlt_dqg
       drhong(k) =                                               mlt_dng

       !--- for charge density
       if (flg_lt) then
          sw2 = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NI)-SMALL ) !--- if NI is small,  ignore charge transfer   ! I -> C
          mlt_dni_crg = mlt_dni*( 1.0_RP-sw2 ) / ( rhoq2(k,I_NI)+sw2 ) * rhoq2_crg(k,I_QI)
          !--- limiter (|rhoq2(NC)| is already reduced by deposition (dep_dni_crg and -frz_dnc_crg) )
          !-- Charge abs(frz_dnc_crg) is already moved from cloud to ice
          sw2 = min( abs(rhoq2_crg(k,I_QI)+dep_dni_crg-frz_dnc_crg),abs(mlt_dni_crg) )
          mlt_dni_crg = sign( sw2, mlt_dni_crg )

          drhoqcrg_c(k) = drhoqcrg_c(k) - mlt_dni_crg * (1.0_RP-sw)
          drhoqcrg_r(k) = drhoqcrg_r(k) - mlt_dni_crg * sw
          drhoqcrg_i(k) = drhoqcrg_i(k) + mlt_dni_crg


          sw2 = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NS)-SMALL ) !--- if NS is small,  ignore charge transfer   ! S -> C
          mlt_dns_crg = mlt_dns*( 1.0_RP-sw2 ) / ( rhoq2(k,I_NS)+sw2 ) * rhoq2_crg(k,I_QS)

          sw2 = min( abs(rhoq2_crg(k,I_QS)+dep_dns_crg            ),abs(mlt_dns_crg) )
          mlt_dns_crg = sign( sw2, mlt_dns_crg )

          drhoqcrg_r(k) = drhoqcrg_r(k) - mlt_dns_crg
          drhoqcrg_s(k) = drhoqcrg_s(k) + mlt_dns_crg


          sw2 = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NG)-SMALL ) !--- if NG is small,  ignore charge transfer   ! G -> C
          mlt_dng_crg = mlt_dng*( 1.0_RP-sw2 ) / ( rhoq2(k,I_NG)+sw2 ) * rhoq2_crg(k,I_QG)
          sw2 = min( abs(rhoq2_crg(k,I_QG)+dep_dng_crg            ),abs(mlt_dng_crg) )
          mlt_dng_crg = sign( sw2, mlt_dng_crg )

          drhoqcrg_r(k) = drhoqcrg_r(k) - mlt_dng_crg
          drhoqcrg_g(k) = drhoqcrg_g(k) + mlt_dng_crg
       end if
    end do

    !--- freezing of larger droplets
!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       frz_dqr = max( dt*(PQ(k,I_LRhet)), min(0.0_RP, -rhoq2(k,I_QR)-dep_dqr(k)) ) ! negative value
       frz_dnr = max( dt*(PQ(k,I_NRhet)), min(0.0_RP, -rhoq2(k,I_NR)-dep_dnr(k)) ) ! negative value

       drhoqr(k) = drhoqr(k) + frz_dqr
       drhonr(k) = drhonr(k) + frz_dnr
       drhoqg(k) = drhoqg(k) - frz_dqr
       drhong(k) = drhong(k) - frz_dnr

       !--- for charge density
       if (flg_lt) then
          sw = 0.5_RP - sign( 0.5_RP, rhoq2(k,I_NR)-SMALL ) !--- if NR is small,  ignore charge transfer
          frz_dnr_crg = frz_dnr*( 1.0_RP-sw ) /( rhoq2(k,I_NR)+sw ) * rhoq2_crg(k,I_QR)
          !--- limiter
          sw = min( abs(rhoq2_crg(k,I_QR)+dep_dnr_crg            ),abs(frz_dnr_crg) )
          frz_dnr_crg = sign( sw,frz_dnr_crg )

          drhoqcrg_r(k) = drhoqcrg_r(k) + frz_dnr_crg
          drhoqcrg_g(k) = drhoqcrg_g(k) - frz_dnr_crg
       end if

       fac5         = ( frz_dqr-eps )/( dt*PQ(k,I_LRhet)-eps )
       PQ(k,I_LRhet) = fac5*PQ(k,I_LRhet)
       fac6         = ( frz_dnr-eps )/( dt*PQ(k,I_NRhet)-eps )
       PQ(k,I_NRhet) = fac6*PQ(k,I_NRhet)
    end do

!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       ! water vapor change
       dqv = -( dep_dqc(k) + dep_dqi(k) + dep_dqs(k) + dep_dqg(k) + dep_dqr(k) )

       ! limiter
       sw = 0.5_RP - sign(0.5_RP, abs(dqv) - eps) ! if |dqv| < eps then sw = 1
       fact = ( max( rhoq2(k,I_QV) + dqv * dt, 0.0_RP ) - rhoq2(k,I_QV) ) / dt / ( dqv + sw ) * ( 1.0_RP - sw ) &
            + 1.0_RP * sw
       fact = min( 1.0_RP, max( 0.0_RP, fact ) )

       drhoqv(k) = dqv * fact

       drhoqc(k) = drhoqc(k) + dep_dqc(k) * fact
       drhonc(k) = drhonc(k) + dep_dnc(k) * fact
       drhoqr(k) = drhoqr(k) + dep_dqr(k) * fact
       drhonr(k) = drhonr(k) + dep_dnr(k) * fact
       drhoqi(k) = drhoqi(k) + dep_dqi(k) * fact
       drhoni(k) = drhoni(k) + dep_dni(k) * fact
       drhoqs(k) = drhoqs(k) + dep_dqs(k) * fact
       drhons(k) = drhons(k) + dep_dns(k) * fact
       drhoqg(k) = drhoqg(k) + dep_dqg(k) * fact
       drhong(k) = drhong(k) + dep_dng(k) * fact

       dz = fz(k) - fz(k-1)
       sl_PLCdep = sl_PLCdep + dep_dqc(k) * dz * fact
       sl_PLRdep = sl_PLRdep + dep_dqr(k) * dz * fact
       sl_PNRdep = sl_PNRdep + dep_dnr(k) * dz * fact
    end do

!OCL LOOP_FISSION_TARGET(LS)
    do k = KS, KE
       ! tendency
       RHOQ_t(k,I_QV) = drhoqv(k) / dt
       RHOQ_t(k,I_QC) = drhoqc(k) / dt
       RHOQ_t(k,I_NC) = drhonc(k) / dt
       RHOQ_t(k,I_QR) = drhoqr(k) / dt
       RHOQ_t(k,I_NR) = drhonr(k) / dt
       RHOQ_t(k,I_QI) = drhoqi(k) / dt
       RHOQ_t(k,I_NI) = drhoni(k) / dt
       RHOQ_t(k,I_QS) = drhoqs(k) / dt
       RHOQ_t(k,I_NS) = drhons(k) / dt
       RHOQ_t(k,I_QG) = drhoqg(k) / dt
       RHOQ_t(k,I_NG) = drhong(k) / dt

       RHOE_t(k) = ( - LHV * drhoqv(k) + LHF * ( drhoqi(k) + drhoqs(k) + drhoqg(k) ) ) / dt

       rrho = 1.0_RP/rho(k)
       dqv = rrho * drhoqv(k)
       dql = rrho * ( drhoqc(k) + drhoqr(k) )
       dqi = rrho * ( drhoqi(k) + drhoqs(k) + drhoqg(k) )

       dcv = CV_VAPOR * dqv + CV_WATER * dql + CV_ICE * dqi
       dcp = CP_VAPOR * dqv + CP_WATER * dql + CP_ICE * dqi

       CVtot_t(k) = dcv/dt
       CPtot_t(k) = dcp/dt

    end do

    ! tendency of charge density
    if (flg_lt) then
       do k = KS, KE
          RHOQcrg_t(k,I_mp_QC) = drhoqcrg_c(k) / dt
          RHOQcrg_t(k,I_mp_QR) = drhoqcrg_r(k) / dt
          RHOQcrg_t(k,I_mp_QI) = drhoqcrg_i(k) / dt
          RHOQcrg_t(k,I_mp_QS) = drhoqcrg_s(k) / dt
          RHOQcrg_t(k,I_mp_QG) = drhoqcrg_g(k) / dt
       end do
    end if


    return
  end subroutine update_by_phase_change
  !-----------------------------------------------------------------------------
  !> Calculate Cross Section
  subroutine Cross_Section( &
       KA, KS, KE, &
       QA_MP,      &
       QTRC0,      &
       DENS0,      &
       Crs         )
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: QA_MP
    real(RP), intent(in)  :: QTRC0(KA,QA_MP)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA)          ! density                   [kg/m3]
    real(RP), intent(out) :: Crs  (KA,HYDRO_MAX)! Cross section             [cm]

    ! mass concentration[kg/m3] and mean particle mass[kg]
    real(RP) :: xc(KA)
    real(RP) :: xr(KA)
    real(RP) :: xi(KA)
    real(RP) :: xs(KA)
    real(RP) :: xg(KA)
    ! diameter of average mass[kg/m3]
!    real(RP) :: dc_ave(KA)
!    real(RP) :: dr_ave(KA)
    ! radius of average mass
    real(RP) :: rc, rr

    real(RP) :: coef_Fuetal1998
    ! r2m_min is minimum value(moment of 1 particle with 1 micron)
    real(RP), parameter :: r2m_min=1.E-12_RP
    real(RP), parameter :: um2cm = 100.0_RP

    real(RP) :: limitsw, zerosw
    integer :: k
    !---------------------------------------------------------------------------

    ! mean particle mass[kg]
    do k  = KS, KE
       xc(k) = min( xc_max, max( xc_min, DENS0(k)*QTRC0(k,I_QC)/(QTRC0(k,I_NC)+nc_min) ) )
       xr(k) = min( xr_max, max( xr_min, DENS0(k)*QTRC0(k,I_QR)/(QTRC0(k,I_NR)+nr_min) ) )
       xi(k) = min( xi_max, max( xi_min, DENS0(k)*QTRC0(k,I_QI)/(QTRC0(k,I_NI)+ni_min) ) )
       xs(k) = min( xs_max, max( xs_min, DENS0(k)*QTRC0(k,I_QS)/(QTRC0(k,I_NS)+ns_min) ) )
       xg(k) = min( xg_max, max( xg_min, DENS0(k)*QTRC0(k,I_QG)/(QTRC0(k,I_NG)+ng_min) ) )
    enddo

    ! diameter of average mass : SB06 eq.(32)
!    do k = KS, KE
!       dc_ave(k) = a_m(I_QC) * xc(k)**b_m(I_QC)
!       dr_ave(k) = a_m(I_QR) * xr(k)**b_m(I_QR)
!    enddo

    do k = KS, KE
       Crs(k,I_mp_QC) = PI * coef_r2(I_mp_QC) * QTRC0(k,I_NC) * a_rea2(I_mp_QC) * xc(k)**b_rea2(I_mp_QC)
    enddo

    do k = KS, KE
       Crs(k,I_mp_QR) = PI * coef_r2(I_mp_QR) * QTRC0(k,I_NR) * a_rea2(I_mp_QR) * xr(k)**b_rea2(I_mp_QR)
    enddo

    do k = KS, KE
       Crs(k,I_mp_QI) = PI * coef_rea2(I_mp_QI) * QTRC0(k,I_NI) * a_rea2(I_mp_QI) * xi(k)**b_rea2(I_mp_QI)
       Crs(k,I_mp_QS) = PI * coef_rea2(I_mp_QS) * QTRC0(k,I_NS) * a_rea2(I_mp_QS) * xs(k)**b_rea2(I_mp_QS)
       Crs(k,I_mp_QG) = PI * coef_rea2(I_mp_QG) * QTRC0(k,I_NG) * a_rea2(I_mp_QG) * xg(k)**b_rea2(I_mp_QG)
    enddo

    return
  end subroutine Cross_Section

end module scale_atmos_phy_mp_sn14
