!-------------------------------------------------------------------------------
!
!+  Double moment Water 6 scheme
!
!-------------------------------------------------------------------------------
#include "macro_thermodyn.h"
module scale_atmos_phy_mp_sn14
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module contains subroutines for the sn14 parametrization.
  !
  !
  !++ Current Corresponding Author : T.Seiki
  !
  !++ History: SN14
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !        0.00   11/10/24 T.Seiki, import from NICAM(11/08/30 ver.)
  !
  !      -----------------------------------------------------------------------
  !      Reference:  -- Journals
  !                   Seifert and Beheng(2006)  : Meteorol.Atmos.Phys.,vol.92,pp.45-66
  !                   Seifert and Beheng(2001)  : Atmos.Res.,vol.59-60,pp.265-281
  !                   Seifert(2008)             : J.Atmos.Sci.,vol.65,pp.3608-3619
  !                   Lin et al.(1983)          : J.Appl.Meteor.,vol.22,pp.1065-1092
  !                   Ruttledge and Hobbs(1983) : J.Atmos.Sci.,vol.40,pp.1185-1206
  !                   Ruttledge and Hobbs(1984) : J.Atmos.Sci.,vol.40,pp.2949-2977
  !                   Cotton etal.(1986)        : J.C.Appl.Meteor.,25,pp.1658-1680
  !                   Cotton and Field (2002)   : QJRMS.,vol.128,pp2417-pp2437
  !                   Beard(1980)               : J.Atmos.Sci.,vol.37,pp.1363-1374 [Add] 10/08/03
  !                   Berry and Reinhardt(1974a): J.Atmos.Sci.,vol.31,pp.1814-1824
  !                   Berry and Reinhardt(1974b): J.Atmos.Sci.,vol.31,pp.1825-1831
  !                   Fu(1996)                  : J.Climate, vol.9, pp.2058-2082   [Add] 10/08/03
  !                   Fu etal(1998)             : J.Climate, vol.11, pp.2223-2237  [Add] 10/08/03
  !                   Ghan etal.(1997)          : J.Geophys.Res.,vol.102,pp.21777-21794, [Add] 09/08/18
  !                   Hong et al.(2004)         : Mon.Wea.Rev.,pp.103-120
  !                   Heymsfeild and Iaquinta(2000): J.Atmos.Sci., vol.57, pp.916-938 [Add] 10/08/03
  !                   Heymsfield and Kajikawa(1987): J.Atmos.Sci., vol.44, pp.1088-1099
  !                   Johnson(1981)             : J.Atmos.Sci., vol.38, pp.215-218 [Add] 09/08/18
  !                   McFarquhar and Heymsfield(1996): J.Atmos.Sci.,vol.53,pp.2401-2423
  !                   Mitchell(1996)            : J.Atmos.Sci., vol.53, pp.1710-1723. [Add] 10/08/03
  !                   Morrison etal.(2005)      : Mon.Wea.Rev.,vol.62,pp.1665-1677, [Add] 09/08/18
  !                   Locatelli and Hobbs (1974): J.Geophys.Res., vol.70, pp.2185-2197
  !                   Lohmann(2002)             : J.Atmos.Sci.,vol.59,pp.647-656
  !                   Takano and Liou(1989)     : J.Atmos.Sci.,vol.46,pp.3-19
  !                   Takano and Liou(1994)     : J.Atmos.Sci.,vol.52,pp.818-837
  !                   Auer and Veal(1970)       : J.Atmos.Sci.,vol.27,pp.919-926
  !                   Ikawa et al.(1991)        : J.M.S.J.,vol.69,pp.641-667
  !                   Murakami(1990)            : J.M.S.J.,vol.68,pp.107-128
  !                  -- Books
  !                   Pruppacher and Klett(1997): Kluwer Academic Publishers
  !                      Microphysics of Clouds and Precipitation, 2nd.edit.
  !                   Seinfeld and Pandis(1998) : Wiley Interscience
  !                      Atmospheric Chemistry and Physics.
  !                   Jacobson(2005)            : Cambridge press
  !                      Fundamentals of Atmospheric Modeling, 2nd.edit.
  !                  -- Source code
  !                   scale_mp_nsw6.f90 in NICAM
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_tracer_sn14
  use scale_const, only: &
     GRAV   => CONST_GRAV,    &
     PI     => CONST_PI,      &
     UNDEF8 => CONST_UNDEF8,  &
     Rdry   => CONST_Rdry,    &
     CPdry  => CONST_CPdry,   &
     CVdry  => CONST_CVdry,   &
     P00    => CONST_PRE00,   &
     T00    => CONST_TEM00,   &
     Rvap   => CONST_Rvap,    &
     CPvap  => CONST_CPvap,   &
     CVvap  => CONST_CVvap,   &
     CL     => CONST_CL,      &
     CI     => CONST_CI,      &
     LHV    => CONST_LH00,    &
     LHF    => CONST_LHF00,   &
     LHV0    => CONST_LH0,    &
     LHF0    => CONST_LHF0,   &
     LHS0    => CONST_LHS0,   &
     LHV00    => CONST_LH00,  &
     LHF00    => CONST_LHF00, &
     PSAT0  => CONST_PSAT0,   &
     EMELT  => CONST_EMELT,   &
     DWATR  => CONST_DWATR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_sn14_setup
  public :: ATMOS_PHY_MP_sn14
  public :: ATMOS_PHY_MP_sn14_CloudFraction
  public :: ATMOS_PHY_MP_sn14_EffectiveRadius
  public :: ATMOS_PHY_MP_sn14_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, target :: ATMOS_PHY_MP_DENS(MP_QA) ! hydrometeor density [kg/m3]=[g/L]

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
  integer, parameter :: HYDRO_MAX  = 6    ! total number of mixing ratio of water
  character(len=H_SHORT), save :: WLABEL(11)

  ! empirical value from Meyers etal.(1991), 1[/liter] = 1.d3[/m3]
  real(RP), private, parameter :: nqmin(6) = (/ 0.0_RP, 1.E+4_RP, 1.0_RP, 1.0_RP, 1.E-4_RP, 1.E-4_RP /) ! [1/m3]
  ! refer to Seifert(2002) (Dr. Thesis, Table.5.1)
  ! max mass, for D_min=79um, 2mm, 5mm, 1cm, 1cm
  real(RP), private, parameter :: xqmax(6) = (/ 0.0_RP, 2.6E-10_RP, 5.0E-6_RP, 1.377E-6_RP, 7.519E-6_RP, 4.90E-5_RP /)! [kg]
  ! SB06, Table 1.
  ! min mass, for D_min=2um, 79um, 10um, 20um, 100um
  real(RP), private, parameter :: xqmin(6) = (/ 0.0_RP, 4.20E-15_RP, 2.60E-10_RP, 3.382E-13_RP, 1.847E-12_RP, 1.230E-10_RP /)! [kg]



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
  ! for Seifert(2008)
  ! work parameter for gamma function, imported from MISC_gammafunc
  real(RP), private, parameter :: gfac_coef(6)=(/&
       +76.180091729471460_RP, -86.505320329416770_RP, +24.014098240830910_RP,&
       -1.2317395724501550_RP, +0.1208650973866179E-2_RP, -0.5395239384953E-5_RP /)
  real(RP), private, parameter :: gfac_ser0=1.000000000190015_RP
  !
  integer, private, save :: ntmax_phase_change = 1
  integer, private, save :: ntmax_collection   = 1
  integer, private, save :: ntmax_sedimentation= 1 ! 10/08/03 [Add] T.Mitsui
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
  real(RP), private, save :: MP_DTSEC_SEDIMENTATION

  !
  ! metrics of vertical coordinate
  ! not used in SCALE-LES
  !
  real(RP), private, allocatable, save :: gsgam2_d (:,:,:)
  real(RP), private, allocatable, save :: gsgam2h_d(:,:,:)
  real(RP), private, allocatable, save :: gam2_d   (:,:,:)
  real(RP), private, allocatable, save :: gam2h_d  (:,:,:)
  real(RP), private, allocatable, save :: rgsgam2_d(:,:,:)
  real(RP), private, allocatable, save :: rgs_d    (:,:,:)
  real(RP), private, allocatable, save :: rgsh_d   (:,:,:)

  logical, private, save :: doautoconversion = .true.
  logical, private, save :: doprecipitation  = .true.
  real(RP), private, save :: MP_ssw_lim = 1.E+1_RP

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sn14_setup( MP_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_DWATR, &
       CONST_DICE
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    implicit none
    character(len=H_SHORT), intent(in) :: MP_TYPE
    NAMELIST / PARAM_ATMOS_PHY_MP / &
       doautoconversion, &
       doprecipitation,  &
       MP_ssw_lim

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for SN14'

    if ( MP_TYPE .ne. 'SN14' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_MP_TYPE is not SN14. Check!'
       call PRC_MPIstop
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    ATMOS_PHY_MP_DENS(I_mp_QC) = CONST_DWATR
    ATMOS_PHY_MP_DENS(I_mp_QR) = CONST_DWATR
    ATMOS_PHY_MP_DENS(I_mp_QI) = CONST_DICE
    ATMOS_PHY_MP_DENS(I_mp_QS) = CONST_DICE
    ATMOS_PHY_MP_DENS(I_mp_QG) = CONST_DICE

    WLABEL( 1) = "VAPOR"
    WLABEL( 2) = "CLOUD"
    WLABEL( 3) = "RAIN"
    WLABEL( 4) = "ICE"
    WLABEL( 5) = "SNOW"
    WLABEL( 6) = "GRAUPEL"
    WLABEL( 7) = "CLOUD_NUM"
    WLABEL( 8) = "RAIN_NUM"
    WLABEL( 9) = "ICE_NUM"
    WLABEL(10) = "SNOW_NUM"
    WLABEL(11) = "GRAUPEL_NUM"

    call mp_sn14_init( IA, JA )


    MP_NSTEP_SEDIMENTATION  = ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(ntmax_sedimentation,kind=RP)
    MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

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
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC  )
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QAD)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics'

    call PROF_rapstart('MP0 Setup')
    call MP_negativefilter( DENS, QTRC )
    call PROF_rapend  ('MP0 Setup')

    call mp_sn14( DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )


    call PROF_rapstart('MP6 Filter')
    call MP_negativefilter( DENS, QTRC )
    call PROF_rapend  ('MP6 Filter')

    return
  end subroutine ATMOS_PHY_MP_sn14

  !-----------------------------------------------------------------------------
  subroutine mp_sn14_init ( IAA, JA )
    use scale_process, only: &
       PRC_MPIstop
    use scale_specfunc, only: &
        gammafunc => SF_gamma
    implicit none

    integer, intent(in) :: IAA, JA

    real(RP), allocatable :: w1(:),w2(:),w3(:),w4(:),w5(:),w6(:),w7(:),w8(:)
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
         ntmax_collection,           &
         ntmax_sedimentation
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
    !

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[SN14]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_sn14_init,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist nm_mp_sn14_init. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=nm_mp_sn14_init)

    !
    ! default setting
    !
    ! Area parameters with mks unit originated by Mitchell(1996)
    a_area(I_QC) = PI/4.0_RP ! sphere
    a_area(I_QR) = PI/4.0_RP ! sphere
    a_area(I_QI) = 0.65_RP*1.E-4_RP*100.0_RP**(2.00_RP)   ! Mitchell(1996), Hexagonal Plate
    a_area(I_QS) = 0.2285_RP*1.E-4_RP*100.0_RP**(1.88_RP) ! Mitchell(1996), Aggregates
    a_area(I_QG) = 0.50_RP*1.E-4_RP*100.0_RP**(2.0_RP)    ! Mitchell(1996), Lump Graupel
    b_area(I_QC) = 2.0_RP
    b_area(I_QR) = 2.0_RP
    b_area(I_QI) = 2.0_RP
    b_area(I_QS) = 1.88_RP
    b_area(I_QG) = 2.0_RP
    !
    ! Seifert and Beheng(2006), Table. 1 or List of symbols
    !----------------------------------------------------------
    ! Diameter-Mass relationship
    ! D = a * x^b
    a_m(I_QC) = 0.124_RP
    a_m(I_QR) = 0.124_RP
    a_m(I_QI) = 0.217_RP
    a_m(I_QS) = 8.156_RP
    a_m(I_QG) = 0.190_RP
    b_m(I_QC) = 1.0_RP/3.0_RP
    b_m(I_QR) = 1.0_RP/3.0_RP
    b_m(I_QI) = 0.302_RP
    b_m(I_QS) = 0.526_RP
    b_m(I_QG) = 0.323_RP
    !----------------------------------------------------------
    ! Terminal velocity-Mass relationship
    ! vt = alpha * x^beta * (rho0/rho)^gamma
    alpha_v(I_QC,:)= 3.75E+5_RP
    alpha_v(I_QR,:)= 159.0_RP ! not for sedimantation
    alpha_v(I_QI,:)= 317.0_RP
    alpha_v(I_QS,:)= 27.70_RP
    alpha_v(I_QG,:)= 40.0_RP
    beta_v(I_QC,:) = 2.0_RP/3.0_RP
    beta_v(I_QR,:) = 0.266_RP ! not for sedimantation
    beta_v(I_QI,:) = 0.363_RP
    beta_v(I_QS,:) = 0.216_RP
    beta_v(I_QG,:) = 0.230_RP
    gamma_v(I_QC)  = 1.0_RP
    ! This is high Reynolds number limit(Beard 1980)
    gamma_v(I_QR)  = 1.0_RP/2.0_RP
    gamma_v(I_QI)  = 1.0_RP/2.0_RP
    gamma_v(I_QS)  = 1.0_RP/2.0_RP
    gamma_v(I_QG)  = 1.0_RP/2.0_RP
    !----------------------------------------------------------
    ! DSD parameters
    ! f(x) = A x^nu exp( -lambda x^mu )
    ! Gamma Disribution           : mu=1  , nu:arbitrary
    ! Marshall-Palmer Distribution: mu=1/3, nu:-2/3
    ! In the case of MP, f(D) dD = f(x)dx
    !                    f(x)    = c * f(D)/D^2 (c:coefficient)
    nu(I_QC) =  1.0_RP          ! arbitrary for Gamma
    nu(I_QR) = -1.0_RP/3.0_RP     ! nu(diameter)=1, equilibrium condition.
    nu(I_QI) =  1.0_RP          !
    nu(I_QS) =  1.0_RP          !
    nu(I_QG) =  1.0_RP          !
    !
    mu(I_QC) = 1.0_RP           ! Gamma
    mu(I_QR) = 1.0_RP/3.0_RP      ! Marshall Palmer
    mu(I_QI) = 1.0_RP/3.0_RP      !
    mu(I_QS) = 1.0_RP/3.0_RP      !
    mu(I_QG) = 1.0_RP/3.0_RP      !
    !----------------------------------------------------------
    ! Geomeries for diffusion growth
    ! Pruppacher and Klett(1997), (13-77)-(13-80) and
    ! originally derived by McDonald(1963b)
    ! sphere: cap=2
    ! plate : cap=pi
    ! needle with aspect ratio a/b
    !       : cap=log(2*a/b)
    cap(I_QC) = 2.0_RP    ! sphere
    cap(I_QR) = 2.0_RP    ! sphere
    cap(I_QI) = PI ! hexagonal plate
    cap(I_QS) = 2.0_RP    ! mix aggregates
    cap(I_QG) = 2.0_RP    ! lump
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
    if( IO_L ) write(IO_FID_LOG,nml=nm_mp_sn14_particles)

    ! [Add] 10/08/03 T.Mitsui
    ! particles shapes are
    if( opt_M96_ice ) then
       ! ice is randomly oriented Hexagonal plate (Auer and Veal 1970, Takano and Liou 1995, Mitchell 1996)
       ! snow is assemblages of planar polycrystals(Mitchell 1996)
       ! graupel is Lump graupel(R4b) is assumed(Mitchell 1996)
       a_area(I_QI)   = 0.120284936_RP
       a_area(I_QS)   = 0.131488_RP
       a_area(I_QG)   = 0.5_RP
       b_area(I_QI)   = 1.850000_RP
       b_area(I_QS)   = 1.880000_RP
       b_area(I_QG)   = 2.0_RP
       a_m(I_QI)      = 1.23655360084766_RP
       a_m(I_QS)      = a_m(I_QI)
       a_m(I_QG)      = 0.346111225718402_RP
       b_m(I_QI)      = 0.408329930583912_RP
       b_m(I_QS)      = b_m(I_QI)
       b_m(I_QG)      = 0.357142857142857_RP
       !
       if( opt_M96_column_ice )then
          d0_ni=240.49E-6_RP   ! this is column
          d0_li=330.09E-6_RP   ! this is column
          a_area(I_QI)= (0.684_RP*1.E-4_RP)*10.0_RP**(2.0_RP*2.00_RP)
          b_area(I_QI)= 2.0_RP
          a_m(I_QI)   = 0.19834046116844_RP
          b_m(I_QI)   = 0.343642611683849_RP
          ! [Add] 11/08/30 T.Mitsui
          ! approximated by the capacity for prolate spheroid with constant aspect ratio
          wcap1       = sqrt(1.0_RP-ar_ice_fix**2)
          wcap2       = log( (1.0_RP+wcap1)/ar_ice_fix )
          cap(I_QI)   = 2.0_RP*wcap2/wcap1
          !
       end if
       !
       ! These value are derived by least-square fitting in the range
       ! qi [100um:1000um] in diameter
       ! qs [100um:1000um] in diameter
       ! qg [200um:2000um] in diameter
       !                      small branch          , large branch
       alpha_v (I_QI,:) =(/ 5798.60107421875_RP, 167.347076416016_RP/)
       alpha_vn(I_QI,:) =(/ 12408.177734375_RP, 421.799865722656_RP/)
       if( opt_M96_column_ice )then
          alpha_v (I_QI,:) = (/2901.0_RP, 32.20_RP/)
          alpha_vn(I_QI,:) = (/9675.2_RP, 64.16_RP/)
       end if
       alpha_v (I_QS,:)    =(/ 15173.3916015625_RP, 305.678619384766_RP/)
       alpha_vn(I_QS,:)    =(/ 29257.1601562500_RP, 817.985717773438_RP/)
       alpha_v (I_QG,:)    =(/ 15481.6904296875_RP, 311.642242431641_RP/)
       alpha_vn(I_QG,:)    =(/ 27574.6562500000_RP, 697.536132812500_RP/)
       !
       beta_v (I_QI,:)  =(/ 0.504873454570770_RP, 0.324817866086960_RP/)
       beta_vn(I_QI,:)  =(/ 0.548495233058929_RP, 0.385287821292877_RP/)
       if( opt_M96_column_ice )then
          beta_v (I_QI,:)  =(/ 0.465552181005478_RP, 0.223826110363007_RP/)
          beta_vn(I_QI,:)  =(/ 0.530453503131866_RP, 0.273761242628098_RP/)
       end if
       beta_v (I_QS,:)     =(/ 0.528109610080719_RP, 0.329863965511322_RP/)
       beta_vn(I_QS,:)     =(/ 0.567154467105865_RP, 0.393876969814301_RP/)
       beta_v (I_QG,:)     =(/ 0.534656763076782_RP, 0.330253750085831_RP/)
       beta_vn(I_QG,:)     =(/ 0.570551633834839_RP, 0.387124240398407_RP/)
    end if
    !
    ! area-diameter relation => area-mass relation
    ax_area(I_QC:I_QG) = a_area(I_QC:I_QG)*a_m(I_QC:I_QG)**b_area(I_QC:I_QG)
    bx_area(I_QC:I_QG) = b_area(I_QC:I_QG)*b_m(I_QC:I_QG)
    !
    ! radius of equivalent area - m ass relation
    ! pi*rea**2 = ax_area*x**bx_area
    a_rea(I_QC:I_QG)   = sqrt(ax_area(I_QC:I_QG)/PI)
    b_rea(I_QC:I_QG)   = bx_area(I_QC:I_QG)/2.0_RP
    a_rea2(I_QC:I_QG)  = a_rea(I_QC:I_QG)**2
    b_rea2(I_QC:I_QG)  = b_rea(I_QC:I_QG)*2.0_RP
    a_rea3(I_QC:I_QG)  = a_rea(I_QC:I_QG)**3
    b_rea3(I_QC:I_QG)  = b_rea(I_QC:I_QG)*3.0_RP
    !
    a_d2vt(I_QC:I_QG)=alpha_v(I_QC:I_QG,2)*(1.0_RP/alpha_v(I_QC:I_QG,2))**(beta_v(I_QC:I_QG,2)/b_m(I_QC:I_QG))
    b_d2vt(I_QC:I_QG)=(beta_v(I_QC:I_QG,2)/b_m(I_QC:I_QG))
    !
    ! Calculation of Moment Coefficient
    !
    allocate( w1(I_QC:I_QG), w2(I_QC:I_QG), w3(I_QC:I_QG), w4(I_QC:I_QG) )
    allocate( w5(I_QC:I_QG), w6(I_QC:I_QG), w7(I_QC:I_QG), w8(I_QC:I_QG) )
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
    ! M^n = coef_mn * N * (L/N)**n
    ! M^2 = Z = coef_m2 * N *(L/N)**2
    ! a*M^b = a*integral x^b f(x) dx = ave D
    do iw=I_QC,I_QG
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
    coef_deplc = coef_d(I_QC)/a_m(I_QC)
    !-------------------------------------------------------
    ! coefficient of 2nd and 3rd moments for effective radius
    ! for spherical particle
    do iw=I_QC,I_QG
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
    do iw=I_QC,I_QG
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
    do iw=I_QC,I_QG
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
       do iw=I_QC,I_QG
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
    do iw=I_QC,I_QG
       w1(iw) = gammafunc( (       b_m(iw) + nu(iw) + 1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (1.0_RP + b_m(iw) + nu(iw) + 1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w4(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       coef_dave_N(iw) =  (w1(iw)/w3(iw))*(w3(iw)/w4(iw))**      b_m(iw)
       coef_dave_L(iw) =  (w2(iw)/w3(iw))*(w3(iw)/w4(iw))**(1.0_RP+b_m(iw))
    end do
    !-------------------------------------------------------
    !
    ah_vent(I_QC,1:2) = (/1.0000_RP,1.0000_RP/) ! no effect
    ah_vent(I_QR,1:2) = (/1.0000_RP,0.780_RP/)
    ah_vent(I_QI,1:2) = (/1.0000_RP,0.860_RP/)
    ah_vent(I_QS,1:2) = (/1.0000_RP,0.780_RP/)
    ah_vent(I_QG,1:2) = (/1.0000_RP,0.780_RP/)
    bh_vent(I_QC,1:2) = (/0.0000_RP,0.0000_RP/)
    bh_vent(I_QR,1:2) = (/0.108_RP,0.308_RP/)
    bh_vent(I_QI,1:2) = (/0.140_RP,0.280_RP/)
    bh_vent(I_QS,1:2) = (/0.108_RP,0.308_RP/)
    bh_vent(I_QG,1:2) = (/0.108_RP,0.308_RP/)
    !
    do iw=I_QC,I_QG
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
    do iw=I_QC,I_QG
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
    do iw=I_QC,I_QG
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
    do iw=I_QC,I_QG
       n = 0
       w1(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       w4(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP    )/mu(iw) )
       n = 1
       w5(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
    end do
    ! ia > ib ( larger particles "a" catch smaller particles "b" )
    do ia=I_QC,I_QG
       do ib=I_QC,I_QG
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
    do iw=I_QC,I_QG
       n = 0
       w1(iw) = gammafunc( (2.0_RP*beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
       w2(iw) = gammafunc( (                        2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
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
    do iw=I_QC,I_QG
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
    do ia=I_QC,I_QG
       do ib=I_QC,I_QG
          theta_ab0(ia,ib) = 2.0_RP * (w1(ib)/w2(ib))*(w3(ia)/w4(ia)) &
               * (w5(ia)/w6(ia))**beta_v(ia,2) &
               * (w5(ib)/w6(ib))**beta_v(ib,2)
          theta_ab1(ia,ib) = 2.0_RP * (w7(ib)/w8(ib))*(w3(ia)/w4(ia)) &
               * (w5(ia)/w6(ia))**beta_v(ia,2) &
               * (w5(ib)/w6(ib))**beta_v(ib,2)
       end do
    end do

    deallocate(w1,w2,w3,w4,w5,w6,w7,w8)

    if( IO_L ) write(IO_FID_LOG,'(100a16)')     "LABEL       ",WLABEL(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "capacity    ",cap(:) ! [Add] 11/08/30 T.Mitsui
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_m2     ",coef_m2(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_d      ",coef_d(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_d3     ",coef_d3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_d6     ",coef_d6(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_d2v    ",coef_d2v(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_md2v   ",coef_md2v(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "a_d2vt      ",a_d2vt(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "b_d2vt      ",b_d2vt(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_r2     ",coef_r2(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_r3     ",coef_r3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_re     ",coef_re(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "a_area      ",a_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "b_area      ",b_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "ax_area     ",ax_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "bx_area     ",bx_area(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "a_rea       ",a_rea(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "b_rea       ",b_rea(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "a_rea3      ",a_rea3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "b_rea3      ",b_rea3(:)
    !
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_rea2   ",coef_rea2(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_rea3   ",coef_rea3(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_vt0    ",coef_vt0(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_vt1    ",coef_vt1(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_A      ",coef_A(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "coef_lambda ",coef_lambda(:)

    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "ah_vent0 sml",ah_vent0(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "ah_vent0 lrg",ah_vent0(:,2)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "ah_vent1 sml",ah_vent1(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "ah_vent1 lrg",ah_vent1(:,2)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "bh_vent0 sml",bh_vent0(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "bh_vent0 lrg",bh_vent0(:,2)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "bh_vent1 sml",bh_vent1(:,1)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "bh_vent1 lrg",bh_vent1(:,2)

    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "delta_b0    ",delta_b0(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "delta_b1    ",delta_b1(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "theta_b0    ",theta_b0(:)
    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "theta_b1    ",theta_b1(:)

    do ia=QQS,QQE
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "delta0(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab0(ia,ib),ib=QQS,QQE)
    enddo
    do ia=QQS,QQE
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "delta1(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab1(ia,ib),ib=QQS,QQE)
    enddo
    do ia=QQS,QQE
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "theta0(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab0(ia,ib),ib=QQS,QQE)
    enddo
    do ia=QQS,QQE
       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "theta1(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab1(ia,ib),ib=QQS,QQE)
    enddo

    allocate(nc_uplim_d(1,IAA,JA))
    nc_uplim_d(:,:,:) = 150.d6

    return
  end subroutine mp_sn14_init
  !-----------------------------------------------------------------------------
  subroutine mp_sn14 ( &
    DENS, MOMZ, MOMX, MOMY, RHOT, QTRC )
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_grid, only: &
       z    => GRID_CZ, &
       dz   => GRID_CDZ
    use scale_atmos_phy_mp_common, only: &
         MP_precipitation => ATMOS_PHY_MP_precipitation
    use scale_atmos_thermodyn, only: &
       AQ_CV, &
       AQ_CP
    use scale_atmos_saturation, only: &
       moist_psat_liq      => ATMOS_SATURATION_psat_liq,   &
       moist_psat_ice      => ATMOS_SATURATION_psat_ice
    use scale_history, only: &
       HIST_in
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    !
    ! primary variables
    !
    real(RP) :: rhogq_d (KA,IA,JA,QA)
    real(RP) :: th_d    (KA,IA,JA)
    !
    ! diagnostic variables
    !
    real(RP) :: qd_d(KA,IA,JA)
    real(RP) :: q_d(KA,IA,JA,QA)
    real(RP) :: tem_d(KA,IA,JA)
    real(RP) :: pre_d(KA,IA,JA)
    real(RP) :: rhoge_d(KA,IA,JA)
    real(RP) :: w_d(KA,IA,JA)
    !
    real(RP) :: lv_d(KA,IA,JA)
    real(RP) :: lc_d(KA,IA,JA), nc_d(KA,IA,JA)
    real(RP) :: lr_d(KA,IA,JA), nr_d(KA,IA,JA)
    real(RP) :: li_d(KA,IA,JA), ni_d(KA,IA,JA)
    real(RP) :: ls_d(KA,IA,JA), ns_d(KA,IA,JA)
    real(RP) :: lg_d(KA,IA,JA), ng_d(KA,IA,JA)
    !
    real(RP) :: xc_d(KA,IA,JA) ! lc/nc
    real(RP) :: xr_d(KA,IA,JA) ! lr/nr
    real(RP) :: xi_d(KA,IA,JA) ! li/ni
    real(RP) :: xs_d(KA,IA,JA) ! ls/ns
    real(RP) :: xg_d(KA,IA,JA) ! lg/ng
    !
    real(RP) :: dc_xa_d(KA,IA,JA) !
    real(RP) :: dr_xa_d(KA,IA,JA) !
    real(RP) :: di_xa_d(KA,IA,JA) !
    real(RP) :: ds_xa_d(KA,IA,JA) !
    real(RP) :: dg_xa_d(KA,IA,JA) !
    real(RP) :: vt_xa_d(KA,IA,JA,HYDRO_MAX,2) !

    real(RP) :: wtem_d(KA,IA,JA)      ! filtered temperature
    real(RP) :: esw_d(KA,IA,JA)       ! saturated vapor pressure(water)
    real(RP) :: esi_d(KA,IA,JA)       ! saturated vapor pressure(ice)
    !
    real(RP) :: rho_fac
    real(RP) :: rho_fac_c_d(KA,IA,JA) ! factor for cloud
    real(RP) :: rho_fac_r_d(KA,IA,JA) !            rain
    real(RP) :: rho_fac_i_d(KA,IA,JA) !            ice
    real(RP) :: rho_fac_s_d(KA,IA,JA) !            snow
    real(RP) :: rho_fac_g_d(KA,IA,JA) !            graupel
    real(RP) :: cva_d(KA,IA,JA)    !
    real(RP) :: cpa_d(KA,IA,JA)       ! [Add] 09/08/18 T.Mitsui
    !
    real(RP) :: drhogqv               ! d (rho*qv*gsgam2)
    real(RP) :: drhogqc, drhognc      !        qc, nc
    real(RP) :: drhogqr, drhognr      !        qr, nr
    real(RP) :: drhogqi, drhogni      !        qi, ni
    real(RP) :: drhogqs, drhogns      !        qs, ns
    real(RP) :: drhogqg, drhogng      !        qg, ng

    ! production rate of nucleation
    real(RP) :: PLCccn_d(KA,IA,JA), PNCccn_d(KA,IA,JA)
    real(RP) :: PLIccn_d(KA,IA,JA), PNIccn_d(KA,IA,JA)
    ! production rate of freezing
    real(RP) :: PLChom_d(KA,IA,JA), PNChom_d(KA,IA,JA)
    real(RP) :: PLChet_d(KA,IA,JA), PNChet_d(KA,IA,JA)
    real(RP) :: PLRhet_d(KA,IA,JA), PNRhet_d(KA,IA,JA)
    ! production rate of melting
    real(RP) :: PLImlt_d(KA,IA,JA), PNImlt_d(KA,IA,JA)
    real(RP) :: PLSmlt_d(KA,IA,JA), PNSmlt_d(KA,IA,JA)
    real(RP) :: PLGmlt_d(KA,IA,JA), PNGmlt_d(KA,IA,JA)
    ! production rate of vapor deposition
    real(RP) :: PLRdep_d(KA,IA,JA), PNRdep_d(KA,IA,JA)
    real(RP) :: PLIdep_d(KA,IA,JA), PNIdep_d(KA,IA,JA)
    real(RP) :: PLSdep_d(KA,IA,JA), PNSdep_d(KA,IA,JA)
    real(RP) :: PLGdep_d(KA,IA,JA), PNGdep_d(KA,IA,JA)
    real(RP) :: PLCdep_d(KA,IA,JA) !, PNCdep_d(KA,IA,JA)

    ! production rate of warm collection process
    ! auto-conversion
    real(RP) :: PLCaut_d(KA,IA,JA), PNCaut_d(KA,IA,JA)
    real(RP) :: PNRaut_d(KA,IA,JA)
    ! accretion
    real(RP) :: PLCacc_d(KA,IA,JA), PNCacc_d(KA,IA,JA)
    ! self-colletion, break-up
    real(RP) :: PNRslc_d(KA,IA,JA), PNRbrk_d(KA,IA,JA)
    real(RP) :: wrm_dqc, wrm_dnc
    real(RP) :: wrm_dqr, wrm_dnr
    ! production rate of mixed-phase collection process
    ! PXXacYY2ZZ means XX collect YY produce ZZ
    real(RP) :: PLIacLC2LI_d(KA,IA,JA)  , PNIacNC2NI_d(KA,IA,JA)   ! cloud-ice
    !
    real(RP) :: PLSacLC2LS_d(KA,IA,JA),   PNSacNC2NS_d(KA,IA,JA)   ! cloud-snow(cloud change)
    !
    real(RP) :: PLGacLC2LG_d(KA,IA,JA),   PNGacNC2NG_d(KA,IA,JA)   ! cloud-graupel
    real(RP) :: PLRacLI2LG_I_d(KA,IA,JA), PNRacNI2NG_I_d(KA,IA,JA) ! rain-ice(ice change)
    real(RP) :: PLRacLI2LG_R_d(KA,IA,JA), PNRacNI2NG_R_d(KA,IA,JA) ! rain-ice(rain change)
    real(RP) :: PLRacLS2LG_S_d(KA,IA,JA), PNRacNS2NG_S_d(KA,IA,JA) ! rain-snow(snow change)
    real(RP) :: PLRacLS2LG_R_d(KA,IA,JA), PNRacNS2NG_R_d(KA,IA,JA) ! rain-snow(rain change)
    real(RP) :: PLRacLG2LG_d(KA,IA,JA),   PNRacNG2NG_d(KA,IA,JA)   ! rain-graupel(rain change)
    real(RP) :: PLIacLI2LS_d(KA,IA,JA),   PNIacNI2NS_d(KA,IA,JA)   ! ice-ice
    real(RP) :: PLIacLS2LS_d(KA,IA,JA),   PNIacNS2NS_d(KA,IA,JA)   ! ice-snow(ice change)
    real(RP) :: PNSacNS2NS_d(KA,IA,JA)                             ! snow-snow
    real(RP) :: PNGacNG2NG_d(KA,IA,JA)                             ! graupel-graupel
    real(RP) :: PLGacLS2LG_d(KA,IA,JA),   PNGacNS2NG_d(KA,IA,JA)   ! snow-graupel

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
    real(RP) :: PLIcon_d(KA,IA,JA), PNIcon_d(KA,IA,JA)
    real(RP) :: PLScon_d(KA,IA,JA), PNScon_d(KA,IA,JA)
    real(RP) :: pco_dqi, pco_dni
    real(RP) :: pco_dqs, pco_dns
    real(RP) :: pco_dqg, pco_dng
    ! production rate of enhanced melting due to
    real(RP) :: PLIacm_d(KA,IA,JA), PNIacm_d(KA,IA,JA) ! ice-cloud
    real(RP) :: PLIarm_d(KA,IA,JA), PNIarm_d(KA,IA,JA) ! ice-rain
    real(RP) :: PLSacm_d(KA,IA,JA), PNSacm_d(KA,IA,JA) ! snow-cloud
    real(RP) :: PLSarm_d(KA,IA,JA), PNSarm_d(KA,IA,JA) ! snow-rain
    real(RP) :: PLGacm_d(KA,IA,JA), PNGacm_d(KA,IA,JA) ! graupel-cloud
    real(RP) :: PLGarm_d(KA,IA,JA), PNGarm_d(KA,IA,JA) ! graupel-rain
    real(RP) :: eml_dqc, eml_dnc
    real(RP) :: eml_dqr, eml_dnr
    real(RP) :: eml_dqi, eml_dni
    real(RP) :: eml_dqs, eml_dns
    real(RP) :: eml_dqg, eml_dng
    ! production rate of ice multiplication by splintering
    real(RP) :: PLGspl_d(KA,IA,JA)
    real(RP) :: PLSspl_d(KA,IA,JA)
    real(RP) :: PNIspl_d(KA,IA,JA)
    real(RP) :: spl_dqi, spl_dni
    real(RP) :: spl_dqg, spl_dqs

    real(RP) :: rrhog_d(KA,IA,JA)

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
    real(RP) :: sl_PLCdep_d(1,IA,JA)
    real(RP) :: sl_PLRdep_d(1,IA,JA), sl_PNRdep_d(1,IA,JA) !
    !--------------------------------------------------
    logical :: flag_history_in
    real(RP) :: qke_d(KA,IA,JA)

    real(RP), parameter :: eps       = 1.E-30_RP
    real(RP), parameter :: eps_qv    = 1.E-50_RP
    real(RP), parameter :: eps_rhoge = 1.E-50_RP
    real(RP), parameter :: eps_rhog  = 1.E-50_RP
    real(RP) :: wdt
    integer :: ntdiv

    real(RP) :: Rmoist
    real(RP) :: cpa

    real(RP) :: velw(KA,IA,JA,QA)
    real(RP) :: flux_rain (KA,IA,JA)
    real(RP) :: flux_snow (KA,IA,JA)
    real(RP) :: flux_prec (IA,JA)
    real(RP) :: wflux_rain(KA,IA,JA)
    real(RP) :: wflux_snow(KA,IA,JA)
    integer :: step

    integer :: k, i, j, iq
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
    !----------------------------------------------------------------------------
    !
    ! 1.Nucleation of cloud water and cloud ice
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MPX ijkconvert')

    do iq = 1, QA
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             q_d(k,i,j,iq) = QTRC(k,i,j,iq)
          enddo
          do k = KS, KE
             rhogq_d(k,i,j,iq) = DENS(k,i,j) * q_d(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       sl_PLCdep_d(1,i,j) = 0.0_RP
       sl_PLRdep_d(1,i,j) = 0.0_RP
       sl_PNRdep_d(1,i,j) = 0.0_RP

       do k = KS-1, KE+1
          rrhog_d(k,i,j) = 1.0_RP / DENS(k,i,j)
       enddo
!!       do k = KS-1, KE+1
!!          w_d(k,i,j) = ( MOMZ(k,i,j) + MOMZ(k-1,i,j) ) * rrhog_d(k,i,j)
       w_d(KS-1,i,j) = 0.0_RP
       do k = KS, KE-1
          w_d(k,i,j) = MOMZ(k,i,j) / ( DENS(k,i,j) + DENS(k+1,i,j) ) * 2.0_RP
       enddo
       w_d(KE:KE+1,i,j) = 0.0_RP
       do k = KS, KE
          th_d(k,i,j) = RHOT(k,i,j) * rrhog_d(k,i,j)
       enddo
       do k = KS, KE
          CALC_QDRY( qd_d(k,i,j), q_d, k, i, j, iq )
       enddo
       do k = KS, KE
          CALC_CV( cva_d(k,i,j), qd_d(k,i,j), q_d, k, i, j, iq, CVdry, AQ_CV )
       enddo
       do k = KS, KE
          CALC_R( Rmoist, q_d(k,i,j,I_QV), qd_d(k,i,j), Rdry, Rvap )
          cpa = cva_d(k,i,j) + Rmoist
          CALC_PRE( pre_d(k,i,j), DENS(k,i,j), th_d(k,i,j), Rmoist, cpa, P00 )
          tem_d(k,i,j) = pre_d(k,i,j) / ( DENS(k,i,j) * Rmoist )
       enddo
       do k = KS, KE
          rhoge_d(k,i,j) = DENS(k,i,j) * tem_d(k,i,j) * cva_d(k,i,j)
       enddo
    enddo
    enddo

    if( opt_debug_tem ) call debug_tem_kij( 1, tem_d(:,:,:), DENS(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call PROF_rapend  ('MPX ijkconvert')

    call PROF_rapstart('MP1 Nucleation')

    do j = JS,  JE
    do i = IS,  IE
       do k = KS, KE
          rho_fac         = rho_0 / max(DENS(k,i,j),rho_min)
          rho_fac_c_d(k,i,j) = rho_fac**gamma_v(I_QC)
          rho_fac_r_d(k,i,j) = rho_fac**gamma_v(I_QR)
          rho_fac_i_d(k,i,j) = (pre_d(k,i,j)/pre0_vt)**a_pre0_vt * (tem_d(k,i,j)/tem0_vt)**a_tem0_vt
          rho_fac_s_d(k,i,j) = rho_fac_i_d(k,i,j)
          rho_fac_g_d(k,i,j) = rho_fac_i_d(k,i,j)
       enddo
       do k = KS, KE
          CALC_CP( cpa_d(k,i,j), qd_d(k,i,j), q_d, k, i, j, iq, CVdry, AQ_CP )
          wtem_d(k,i,j) = max(tem_d(k,i,j), tem_min)
       enddo
       do k = KS, KE
          qke_d       (k,i,j) = 0.0_RP ! 2*TKE
          dTdt_equiv_d(k,i,j) = 0.0_RP
       enddo
       do k = KS, KE
          lv_d(k,i,j) = DENS(k,i,j)*q_d(k,i,j,I_QV)
          ni_d(k,i,j) = max( 0.0_RP, DENS(k,i,j)*q_d(k,i,j,I_NI) )
          nc_d(k,i,j) = max( 0.0_RP, DENS(k,i,j)*q_d(k,i,j,I_NC) )
       enddo
    enddo
    enddo

    call nucleation_kij(     &
         IA, JA, KA,     & ! in
         IS, IE,         & ! in
         JS, JE,         & ! in
         KS, KE,         & ! in
         z, w_d,         & ! in
         DENS, wtem_d, pre_d, & ! in
         lv_d, nc_d, ni_d,     & ! in
         PNCccn_d, PLCccn_d, & ! out
         PNIccn_d, PLIccn_d, & ! out
         cpa_d,            & ! in
         dTdt_equiv_d,     & ! in
         qke_d,            & ! in
         flag_history_in,& ! in
         dt              ) ! in

    do j = JS,  JE
    do i = IS,  IE
       do k = KS, KE
          ! nucleation
          drhogqc = dt * PLCccn_d(k,i,j)
          drhognc = dt * PNCccn_d(k,i,j)
          drhogqi = dt * PLIccn_d(k,i,j)
          drhogni = dt * PNIccn_d(k,i,j)
          drhogqv = max( -rhogq_d(k,i,j,I_QV), -drhogqc-drhogqi )
          fac1    = drhogqv / min( -drhogqc-drhogqi, -eps ) ! limiting coefficient

          rhogq_d(k,i,j,I_QV) = rhogq_d(k,i,j,I_QV) + drhogqv
          rhogq_d(k,i,j,I_QC) = max(0.0_RP, rhogq_d(k,i,j,I_QC) + drhogqc*fac1)
          rhogq_d(k,i,j,I_QI) = max(0.0_RP, rhogq_d(k,i,j,I_QI) + drhogqi*fac1)
          rhogq_d(k,i,j,I_NC) = max(0.0_RP, rhogq_d(k,i,j,I_NC) + drhognc)
          rhogq_d(k,i,j,I_NI) = max(0.0_RP, rhogq_d(k,i,j,I_NI) + drhogni)

          ! cloud number concentration filter
          rhogq_d(k,i,j,I_NC) = min( rhogq_d(k,i,j,I_NC), nc_uplim_d(1,i,j) )

          rhoge_d(k,i,j) = rhoge_d(k,i,j) - LHV * drhogqv + LHF * drhogqi*fac1
       enddo
       !
       do k = KS, KE
          q_d(k,i,j,I_QV) = rhogq_d(k,i,j,I_QV) * rrhog_d(k,i,j)
          q_d(k,i,j,I_QC) = rhogq_d(k,i,j,I_QC) * rrhog_d(k,i,j)
          q_d(k,i,j,I_QI) = rhogq_d(k,i,j,I_QI) * rrhog_d(k,i,j)
          q_d(k,i,j,I_NC) = rhogq_d(k,i,j,I_NC) * rrhog_d(k,i,j)
          q_d(k,i,j,I_NI) = rhogq_d(k,i,j,I_NI) * rrhog_d(k,i,j)

          CALC_QDRY( qd_d(k,i,j), q_d, k, i, j, iq )
          CALC_CV( cva_d(k,i,j), qd_d(k,i,j), q_d, k, i, j, iq, CVdry, AQ_CV )
          CALC_R ( Rmoist, q_d(k,i,j,I_QV), qd_d(k,i,j), Rdry, Rvap )

          tem_d(k,i,j) = rhoge_d(k,i,j) / ( DENS(k,i,j) * cva_d(k,i,j) )
          pre_d(k,i,j) = DENS(k,i,j) * Rmoist * tem_d(k,i,j)
       enddo
    enddo
    enddo

!    if( opt_debug )     call debugreport_nucleation
    if( opt_debug_tem ) call debug_tem_kij( 2, tem_d(:,:,:), DENS(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call PROF_rapend  ('MP1 Nucleation')
    !----------------------------------------------------------------------------
    !
    ! 2.Phase change: Freezing, Melting, Vapor deposition
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP2 Phase change')

    ! parameter setting
    wdt=dt

!       if(  ntdiv     == ntmax_phase_change  )then
          flag_history_in=.true.
!       end if
       !
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             lr_d(k,i,j)     = rhogq_d(k,i,j,I_QR)
             nr_d(k,i,j)     = rhogq_d(k,i,j,I_NR)
             xr_d(k,i,j)     = max(xr_min,  min(xr_max, lr_d(k,i,j)/(nr_d(k,i,j)+nr_min) ))
             dr_xa_d(k,i,j)  = a_m(I_QR)*xr_d(k,i,j)**b_m(I_QR)
             vt_xa_d(k,i,j,I_QR,1) = alpha_v(I_QR,1)*(xr_d(k,i,j)**beta_v(I_QR,1))*rho_fac_r_d(k,i,j)
             vt_xa_d(k,i,j,I_QR,2) = vt_xa_d(k,i,j,I_QR,1)
             !! Following values shoud be already filtered to be non-zero before sbroutine was called.
             ! Mass concentration [kg/m3]
             lv_d(k,i,j)     = rhogq_d(k,i,j,I_QV)
             !
             lc_d(k,i,j)     = rhogq_d(k,i,j,I_QC)
             li_d(k,i,j)     = rhogq_d(k,i,j,I_QI)
             ls_d(k,i,j)     = rhogq_d(k,i,j,I_QS)
             lg_d(k,i,j)     = rhogq_d(k,i,j,I_QG)
             ! Number concentration[/m3] (should be filtered to avoid zero division.)
             nc_d(k,i,j)     = rhogq_d(k,i,j,I_NC)
             ni_d(k,i,j)     = rhogq_d(k,i,j,I_NI)
             ns_d(k,i,j)     = rhogq_d(k,i,j,I_NS)
             ng_d(k,i,j)     = rhogq_d(k,i,j,I_NG)
             ! Mass of mean particle [kg]
             ! SB06(94)
             !
             xc_d(k,i,j)     = min(xc_max, max(xc_min, lc_d(k,i,j)/(nc_d(k,i,j)+nc_min) ))
             xi_d(k,i,j)     = min(xi_max, max(xi_min, li_d(k,i,j)/(ni_d(k,i,j)+ni_min) ))
             xs_d(k,i,j)     = min(xs_max, max(xs_min, ls_d(k,i,j)/(ns_d(k,i,j)+ns_min) ))
             xg_d(k,i,j)     = min(xg_max, max(xg_min, lg_d(k,i,j)/(ng_d(k,i,j)+ng_min) ))
             ! diamter of average mass
             ! SB06(32)
             dc_xa_d(k,i,j)  = a_m(I_QC)*xc_d(k,i,j)**b_m(I_QC)
             di_xa_d(k,i,j)  = a_m(I_QI)*xi_d(k,i,j)**b_m(I_QI)
             ds_xa_d(k,i,j)  = a_m(I_QS)*xs_d(k,i,j)**b_m(I_QS)
             dg_xa_d(k,i,j)  = a_m(I_QG)*xg_d(k,i,j)**b_m(I_QG)
             ! terminal velocity of average mass
             vt_xa_d(k,i,j,I_QC,1) = alpha_v(I_QC,1)*(xc_d(k,i,j)**beta_v(I_QC,1))*rho_fac_c_d(k,i,j)
             vt_xa_d(k,i,j,I_QI,1) = alpha_v(I_QI,1)*(xi_d(k,i,j)**beta_v(I_QI,1))*rho_fac_i_d(k,i,j)
             vt_xa_d(k,i,j,I_QS,1) = alpha_v(I_QS,1)*(xs_d(k,i,j)**beta_v(I_QS,1))*rho_fac_s_d(k,i,j)
             vt_xa_d(k,i,j,I_QG,1) = alpha_v(I_QG,1)*(xg_d(k,i,j)**beta_v(I_QG,1))*rho_fac_g_d(k,i,j)
             vt_xa_d(k,i,j,I_QC,2) = alpha_v(I_QC,2)*(xc_d(k,i,j)**beta_v(I_QC,2))*rho_fac_c_d(k,i,j)
             vt_xa_d(k,i,j,I_QI,2) = alpha_v(I_QI,2)*(xi_d(k,i,j)**beta_v(I_QI,2))*rho_fac_i_d(k,i,j)
             vt_xa_d(k,i,j,I_QS,2) = alpha_v(I_QS,2)*(xs_d(k,i,j)**beta_v(I_QS,2))*rho_fac_s_d(k,i,j)
             vt_xa_d(k,i,j,I_QG,2) = alpha_v(I_QG,2)*(xg_d(k,i,j)**beta_v(I_QG,2))*rho_fac_g_d(k,i,j)
             !
             wtem_d(k,i,j) = max( tem_d(k,i,j), tem_min )
             !
             CALC_QDRY( qd_d(k,i,j), q_d, k, i, j, iq )
          end do
       end do
       end do
       !
       call moist_psat_liq( esw_d, wtem_d )
       call moist_psat_ice( esi_d, wtem_d )
       !
       call freezing_water_kij( &
            IA, JA, KA,    & ! in
            IS, IE,     & ! in
            JS, JE,     & ! in
            KS, KE,     & ! in
            wdt,            & ! in
            PLChom_d, PNChom_d, & ! out
            PLChet_d, PNChet_d, & ! out
            PLRhet_d, PNRhet_d, & ! out
            lc_d, lr_d, nc_d, nr_d, xc_d, xr_d, tem_d ) ! in
       !
       call dep_vapor_melt_ice_kij( &
            IA, JA, KA,           & ! in
            IS, IE,               & ! in
            JS, JE,               & ! in
            KS, KE,               & ! in
            PLCdep_d,             & ! out
            PLRdep_d, PNRdep_d,   & ! out
            PLIdep_d, PNIdep_d,   & ! out
            PLSdep_d, PNSdep_d,   & ! out
            PLGdep_d, PNGdep_d,   & ! out
            PLImlt_d, PNImlt_d,   & ! out
            PLSmlt_d, PNSmlt_d,   & ! out
            PLGmlt_d, PNGmlt_d,   & ! out
            DENS, wtem_d, pre_d, lv_d,& ! in
            qd_d,                 & ! in
            esw_d, esi_d,         & ! in
            nc_d,                 & ! in
            nr_d, ni_d, ns_d, ng_d,   & ! in
            li_d, ls_d, lg_d,     & ! in
            xc_d,                 & ! in
            xr_d, xi_d, xs_d, xg_d,   & ! in
            vt_xa_d,              & ! in
            dc_xa_d,              & ! in
            dr_xa_d, di_xa_d,     & ! in
            ds_xa_d, dg_xa_d      ) ! in
       !
       ! update subroutine
       !
       call update_by_phase_change_kij(    &
            ntdiv, ntmax_phase_change,     & ! in
            IA, JA, KA,                    & ! in
            IS, IE, JS, JE, KS, KE,        & ! in
            QA,                            & ! in
            wdt,                           & ! in
            gsgam2_d,                      & ! in
            z,                             & ! in
            dz,                            & ! in
            w_d,                           & ! in
            dTdt_equiv_d,                  & ! in
            DENS,                          & ! in
            rhoge_d,                       & ! inout
            rhogq_d, q_d,                  & ! inout
            tem_d, pre_d,                  & ! inout
            cva_d,                         & ! out
            esw_d, esi_d, LV_d,            & ! in
            LC_d, LR_d, LI_d, LS_d, LG_d,  & ! in
            NC_d, NR_d, NI_d, NS_d, NG_d,  & ! in
            !+++ tendency terms
            ! homogeneous freezing
            PLChom_d, PNChom_d,            & ! in
            ! heterogeneous freezing
            PLChet_d, PNChet_d, PLRhet_d, PNRhet_d, &
            ! condensation/evaporation, deposition/sublimation
            PLCdep_d,           & ! inout
            PLRdep_d, PNRdep_d, & ! inout
            PLIdep_d, PNIdep_d, & ! inout
            PLSdep_d, PNSdep_d, & ! inout
            PLGdep_d, PNGdep_d, & ! inout
            ! melting term
            PLImlt_d, PNImlt_d, & ! in
            PLSmlt_d, PNSmlt_d, & ! in
            PLGmlt_d, PNGmlt_d, & ! in
            flag_history_in,    & ! in
            sl_PLCdep_d,        & ! out
            sl_PLRdep_d, sl_PNRdep_d ) ! out

!       if( opt_debug )     call debugreport_phasechange
       if( opt_debug_tem ) call debug_tem_kij( 3, tem_d(:,:,:), DENS(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call PROF_rapend  ('MP2 Phase change')
    !----------------------------------------------------------------------------
    !
    ! 3.Collection process
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP3 Collection')

    wdt = dt
    flag_history_in=.true.

    ! parameter setting
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          ! Mass concentration [kg/m3]
          lv_d(k,i,j) = rhogq_d(k,i,j,I_QV)
          lc_d(k,i,j) = rhogq_d(k,i,j,I_QC)
          lr_d(k,i,j) = rhogq_d(k,i,j,I_QR)
          li_d(k,i,j) = rhogq_d(k,i,j,I_QI)
          ls_d(k,i,j) = rhogq_d(k,i,j,I_QS)
          lg_d(k,i,j) = rhogq_d(k,i,j,I_QG)
       enddo
       do k = KS, KE
          ! Number concentration[/m3]
          nc_d(k,i,j) = rhogq_d(k,i,j,I_NC)
          nr_d(k,i,j) = rhogq_d(k,i,j,I_NR)
          ni_d(k,i,j) = rhogq_d(k,i,j,I_NI)
          ns_d(k,i,j) = rhogq_d(k,i,j,I_NS)
          ng_d(k,i,j) = rhogq_d(k,i,j,I_NG)
       enddo
       do k = KS, KE
          ! Mass of mean particle [kg]
          xc_d(k,i,j) = min(xc_max, max(xc_min, lc_d(k,i,j)/(nc_d(k,i,j)+nc_min) ) )
          xr_d(k,i,j) = min(xr_max, max(xr_min, lr_d(k,i,j)/(nr_d(k,i,j)+nr_min) ) )
          xi_d(k,i,j) = min(xi_max, max(xi_min, li_d(k,i,j)/(ni_d(k,i,j)+ni_min) ) )
          xs_d(k,i,j) = min(xs_max, max(xs_min, ls_d(k,i,j)/(ns_d(k,i,j)+ns_min) ) )
          xg_d(k,i,j) = min(xg_max, max(xg_min, lg_d(k,i,j)/(ng_d(k,i,j)+ng_min) ) )
       enddo
       do k = KS, KE
          ! effective cross section is assume as area equivalent circle
          dc_xa_d(k,i,j)  = 2.0_RP*a_rea(I_QC)*xc_d(k,i,j)**b_rea(I_QC)
          dr_xa_d(k,i,j)  = 2.0_RP*a_rea(I_QR)*xr_d(k,i,j)**b_rea(I_QR)
          di_xa_d(k,i,j)  = 2.0_RP*a_rea(I_QI)*xi_d(k,i,j)**b_rea(I_QI)
          ds_xa_d(k,i,j)  = 2.0_RP*a_rea(I_QS)*xs_d(k,i,j)**b_rea(I_QS)
          dg_xa_d(k,i,j)  = 2.0_RP*a_rea(I_QG)*xg_d(k,i,j)**b_rea(I_QG)
       enddo
       do k = KS, KE
          ! terminal velocity of average mass
          ! SB06(33)
          vt_xa_d(k,i,j,I_QC,2) = alpha_v(I_QC,2)*(xc_d(k,i,j)**beta_v(I_QC,2))*rho_fac_c_d(k,i,j)
          vt_xa_d(k,i,j,I_QR,2) = alpha_v(I_QR,2)*(xr_d(k,i,j)**beta_v(I_QR,2))*rho_fac_r_d(k,i,j)
          vt_xa_d(k,i,j,I_QI,2) = alpha_v(I_QI,2)*(xi_d(k,i,j)**beta_v(I_QI,2))*rho_fac_i_d(k,i,j)
          vt_xa_d(k,i,j,I_QS,2) = alpha_v(I_QS,2)*(xs_d(k,i,j)**beta_v(I_QS,2))*rho_fac_s_d(k,i,j)
          vt_xa_d(k,i,j,I_QG,2) = alpha_v(I_QG,2)*(xg_d(k,i,j)**beta_v(I_QG,2))*rho_fac_g_d(k,i,j)
       enddo
    enddo
    enddo

    ! Auto-conversion, Accretion, Self-collection, Break-up
    ! [Mod] T.Seiki
    if ( doautoconversion ) then
       call aut_acc_slc_brk_kij(  &
            IA, JA, KA,      &
            IS, IE,       &
            JS, JE,       &
            KS, KE,       &
            PLCaut_d, PNCaut_d,   &
            PNRaut_d,           &
            PLCacc_d, PNCacc_d,   &
            PNRslc_d, PNRbrk_d,   &
            lc_d, lr_d, nc_d, nr_d,xc_d,&
            dr_xa_d,            &
            DENS                )
    else
       PLCaut_d = 0.0_RP
       PNCaut_d = 0.0_RP
       PNRaut_d = 0.0_RP
       PLCacc_d = 0.0_RP
       PNCacc_d = 0.0_RP
       PNRslc_d = 0.0_RP
       PNRbrk_d = 0.0_RP
    endif
       !
       call mixed_phase_collection_kij(         &
            IA, JA, KA,                         & ! in
            IS, IE,                             & ! in
            JS, JE,                             & ! in
            KS, KE,                             & ! in
            ! collection process
            PLIacLC2LI_d,   PNIacNC2NI_d,       & ! out
            PLSacLC2LS_d,   PNSacNC2NS_d,       & ! out
            !
            PLGacLC2LG_d,   PNGacNC2NG_d,       & ! out
            PLRacLI2LG_I_d, PNRacNI2NG_I_d,     & ! out
            PLRacLI2LG_R_d, PNRacNI2NG_R_d,     & ! out
            PLRacLS2LG_S_d, PNRacNS2NG_S_d,     & ! out
            PLRacLS2LG_R_d, PNRacNS2NG_R_d,     & ! out
            PLRacLG2LG_d,   PNRacNG2NG_d,       & ! out
            PLIacLI2LS_d,   PNIacNI2NS_d,       & ! out
            PLIacLS2LS_d,   PNIacNS2NS_d,       & ! out
            PNSacNS2NS_d,                       & ! out
            PNGacNG2NG_d,                       & ! out
            PLGacLS2LG_d,   PNGacNS2NG_d,       & ! out
            ! partial conversion (ice, snow => graupel)
            PLIcon_d, PNIcon_d, PLScon_d, PNScon_d, & ! out
            ! enhanced melting (latent heat effect )
            PLIacm_d, PNIacm_d, PLIarm_d, PNIarm_d, & ! out
            PLSacm_d, PNSacm_d, PLSarm_d, PNSarm_d, & ! out
            PLGacm_d, PNGacm_d, PLGarm_d, PNGarm_d, & ! out
            tem_d,                            & ! in
            lv_d,                             & ! in
            lc_d, lr_d, li_d, ls_d, lg_d,             & ! in
            nc_d, nr_d, ni_d, ns_d, ng_d,             & ! in
            xc_d, xr_d, xi_d, xs_d, xg_d,             & ! in
            dc_xa_d, dr_xa_d, di_xa_d, ds_xa_d, dg_xa_d,& ! in
            vt_xa_d(:,:,:,I_QC,2), vt_xa_d(:,:,:,I_QR,2), &
            vt_xa_d(:,:,:,I_QI,2), vt_xa_d(:,:,:,I_QS,2), vt_xa_d(:,:,:,I_QG,2), &
!            DENS(:,:,:),                    & ! in
            flag_history_in                 ) ! in

       !
       call ice_multiplication_kij( &
            IA, JA, KA,        & ! in
            IS, IE,         & ! in
            JS, JE,         & ! in
            KS, KE,         & ! in
            PLGspl_d, PLSspl_d,     & ! out
            PNIspl_d,             & ! out
            PLIacLC2LI_d,         & ! in
            PLSacLC2LS_d,         & ! in
            PLGacLC2LG_d,         & ! in
            tem_d, lc_d, nc_d, xi_d, xc_d ) ! in
       !
       ! update
       ! rhogq = l*gsgam
       !
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             !
             ! warm collection process
             wrm_dqc = max( wdt*( PLCaut_d(k,i,j)+PLCacc_d(k,i,j) ), -lc_d(k,i,j)  )
             wrm_dnc = max( wdt*( PNCaut_d(k,i,j)+PNCacc_d(k,i,j) ), -nc_d(k,i,j)  )
             wrm_dnr = max( wdt*( PNRaut_d(k,i,j)+PNRslc_d(k,i,j)+PNRbrk_d(k,i,j) ), -nr_d(k,i,j) )
             wrm_dqr = -wrm_dqc
             ! mixed phase collection
             ! Pxxacyy2zz xx and yy decrease and zz increase .
             !
             ! At first fixer is applied to decreasing particles.
             ! order of fixer: graupel-cloud, snow-cloud, ice-cloud, graupel-rain, snow-rain, ice-rain,
             !                 snow-ice,  ice-ice, graupel-snow, snow-snow
             ! cloud mass decrease
             gc_dqc  = max( wdt*PLGacLC2LG_d(k,i,j)  , min(0.0_RP, -lc_d(k,i,j)-wrm_dqc               )) ! => dqg
             sc_dqc  = max( wdt*PLSacLC2LS_d(k,i,j)  , min(0.0_RP, -lc_d(k,i,j)-wrm_dqc-gc_dqc        )) ! => dqs
             ic_dqc  = max( wdt*PLIacLC2LI_d(k,i,j)  , min(0.0_RP, -lc_d(k,i,j)-wrm_dqc-gc_dqc-sc_dqc )) ! => dqi
             ! cloud num. decrease
             gc_dnc  = max( wdt*PNGacNC2NG_d(k,i,j)  , min(0.0_RP, -nc_d(k,i,j)-wrm_dnc               )) ! => dnc:minus
             sc_dnc  = max( wdt*PNSacNC2NS_d(k,i,j)  , min(0.0_RP, -nc_d(k,i,j)-wrm_dnc-gc_dnc        )) ! => dnc:minus
             ic_dnc  = max( wdt*PNIacNC2NI_d(k,i,j)  , min(0.0_RP, -nc_d(k,i,j)-wrm_dnc-gc_dnc-sc_dnc )) ! => dnc:minus
             ! rain mass decrease ( tem < 273.15K)
             if( tem_d(k,i,j) <= T00 )then
                rg_dqr  = max( wdt*PLRacLG2LG_d  (k,i,j), min(0.0_RP, -lr_d(k,i,j)-wrm_dqr               )) ! => dqg
                rg_dqg  = 0.0_RP
                rs_dqr  = max( wdt*PLRacLS2LG_R_d(k,i,j), min(0.0_RP, -lr_d(k,i,j)-wrm_dqr-rg_dqr        )) ! => dqg
                ri_dqr  = max( wdt*PLRacLI2LG_R_d(k,i,j), min(0.0_RP, -lr_d(k,i,j)-wrm_dqr-rg_dqr-rs_dqr )) ! => dqg
                ! rain num. decrease
                rg_dnr  = max( wdt*PNRacNG2NG_d  (k,i,j), min(0.0_RP, -nr_d(k,i,j)-wrm_dnr               )) ! => dnr:minus,dng:plus
                rg_dng  = 0.0_RP
                rs_dnr  = max( wdt*PNRacNS2NG_R_d(k,i,j), min(0.0_RP, -nr_d(k,i,j)-wrm_dnr-rg_dnr        )) ! => dnr:minus,dng:plus
                ri_dnr  = max( wdt*PNRacNI2NG_R_d(k,i,j), min(0.0_RP, -nr_d(k,i,j)-wrm_dnr-rg_dnr-rs_dnr )) ! => dnr:minus,dng:plus
             else
                rg_dqg  = max( wdt*PLRacLG2LG_d  (k,i,j), min(0.0_RP, -lg_d(k,i,j)                       )) ! => dqg
                rg_dqr  = 0.0_RP ! r+g -> r
                rs_dqr  = 0.0_RP ! r+s -> r
                ri_dqr  = 0.0_RP ! r+i -> r
                ! rain num. decrease
                rg_dng  = max( wdt*PNRacNG2NG_d  (k,i,j), min(0.0_RP, -ng_d(k,i,j)                       )) ! => dnr:minus,dng:plus
                rg_dnr  = 0.0_RP ! r+g -> r
                rs_dnr  = 0.0_RP ! r+s -> r
                ri_dnr  = 0.0_RP ! r+i -> r
             end if
             ! ice mass decrease
             fac1    = (ri_dqr-eps)/ (wdt*PLRacLI2LG_R_d(k,i,j)-eps) ! suppress factor by filter of rain
             ri_dqi  = max( wdt*PLRacLI2LG_I_d(k,i,j)*fac1, min(0.0_RP, -li_d(k,i,j)+ic_dqc               )) ! => dqg
             ii_dqi  = max( wdt*PLIacLI2LS_d(k,i,j)       , min(0.0_RP, -li_d(k,i,j)+ic_dqc-ri_dqi        )) ! => dqs
             is_dqi  = max( wdt*PLIacLS2LS_d(k,i,j)       , min(0.0_RP, -li_d(k,i,j)+ic_dqc-ri_dqi-ii_dqi )) ! => dqs
             ! ice num. decrease
             fac4    = (ri_dnr-eps)/ (wdt*PNRacNI2NG_R_d(k,i,j)-eps) ! suppress factor by filter of rain
             ri_dni  = max( wdt*PNRacNI2NG_I_d(k,i,j)*fac4, min(0.0_RP, -ni_d(k,i,j)               )) ! => dni:minus
             ii_dni  = max( wdt*PNIacNI2NS_d(k,i,j)       , min(0.0_RP, -ni_d(k,i,j)-ri_dni        )) ! => dni:minus,dns:plus(*0.5)
             is_dni  = max( wdt*PNIacNS2NS_d(k,i,j)       , min(0.0_RP, -ni_d(k,i,j)-ri_dni-ii_dni )) ! => dni:minus,dns:plus
             ! snow mass decrease
             fac3    = (rs_dqr-eps)/(wdt*PLRacLS2LG_R_d(k,i,j)-eps) ! suppress factor by filter of rain
             rs_dqs  = max( wdt*PLRacLS2LG_S_d(k,i,j)*fac3, min(0.0_RP, -ls_d(k,i,j)+sc_dqc+ii_dqi+is_dqi        )) ! => dqg
             gs_dqs  = max( wdt*PLGacLS2LG_d(k,i,j)       , min(0.0_RP, -ls_d(k,i,j)+sc_dqc+ii_dqi+is_dqi-rs_dqs )) ! => dqg
             ! snow num. decrease
             fac6    = (rs_dnr-eps)/(wdt*PNRacNS2NG_R_d(k,i,j)-eps) ! suppress factor by filter of rain
!             fac7    = (is_dni-eps)/(wdt*PNIacNS2NS_d  (k,i,j)-eps) ! suppress factor by filter of ice
             rs_dns  = max( wdt*PNRacNS2NG_S_d(k,i,j)*fac6, min(0.0_RP, -ns_d(k,i,j)+0.50_RP*ii_dni+is_dni       )) ! => dns:minus
             gs_dns  = max( wdt*PNGacNS2NG_d(k,i,j)       , min(0.0_RP, -ns_d(k,i,j)+0.50_RP*ii_dni+is_dni-rs_dns )) ! => dns:minus
             ss_dns  = max( wdt*PNSacNS2NS_d(k,i,j)       , min(0.0_RP, -ns_d(k,i,j)+0.50_RP*ii_dni+is_dni-rs_dns-gs_dns ))
             !
             gg_dng  = max( wdt*PNGacNG2NG_d(k,i,j)       , min(0.0_RP, -ng_d(k,i,j) ))
             !
             ! total plus in mixed phase collection(clp_)
             ! mass
             if( tem_d(k,i,j) <= T00 )then
                clp_dqc =  0.0_RP
                clp_dqr =  0.0_RP
                clp_dqi = -ic_dqc
                clp_dqs = -sc_dqc-ii_dqi-is_dqi
                clp_dqg = -gc_dqc-rg_dqr-rs_dqr-rs_dqs-ri_dqr-ri_dqi-gs_dqs
                ! num.( number only increase when a+b=>c,  dnc=-dna)
                clp_dnc = 0.0_RP
                clp_dnr = 0.0_RP
                clp_dni = 0.0_RP
                clp_dns = -ii_dni*0.5_RP
                clp_dng =       -rs_dnr-ri_dnr
                ! total minus in mixed phase collection(clm_)
                ! mass
                clm_dqc = gc_dqc+sc_dqc+ic_dqc
                clm_dqr = rg_dqr+rs_dqr+ri_dqr
                clm_dqi = ri_dqi+ii_dqi+is_dqi
                clm_dqs = rs_dqs+gs_dqs
                clm_dqg = 0.0_RP
                ! num.
                clm_dnc = gc_dnc+sc_dnc+ic_dnc
                clm_dnr = rg_dnr+rs_dnr+ri_dnr
                clm_dni = ri_dni+ii_dni+is_dni
                clm_dns = rs_dns+ss_dns+gs_dns
!!$             clm_dng = 0.D0
                clm_dng = gg_dng ! [Mod] 11/08/30 T.Mitsui
             else
                clp_dqc =  0.0_RP
                clp_dqr = -rg_dqg-rs_dqs-ri_dqi
                clp_dqi = -ic_dqc
                clp_dqs = -sc_dqc-ii_dqi-is_dqi
                clp_dqg = -gc_dqc-gs_dqs
                ! num.( number only increase when a+b=>c,  dnc=-dna)
                clp_dnc = 0.0_RP
                clp_dnr = 0.0_RP
                clp_dni = 0.0_RP
                clp_dns = -ii_dni*0.5_RP
                clp_dng = 0.0_RP
                ! total minus in mixed phase collection(clm_)
                ! mass
                clm_dqc = gc_dqc+sc_dqc+ic_dqc
                clm_dqr = 0.0_RP
                clm_dqi = ri_dqi+ii_dqi+is_dqi
                clm_dqs = rs_dqs+gs_dqs
                clm_dqg = rg_dqg
                ! num.
                clm_dnc = gc_dnc+sc_dnc+ic_dnc
                clm_dnr = 0.0_RP
                clm_dni = ri_dni+ii_dni+is_dni
                clm_dns = rs_dns+ss_dns+gs_dns
!!$             clm_dng = rg_dng
                clm_dng = rg_dng+gg_dng ! [Mod] 11/08/30 T.Mitsui
             end if
             !
             ! partial conversion
             ! 08/05/08 [Mod] T.Mitsui
             pco_dqi = max( wdt*PLIcon_d(k,i,j), -clp_dqi )
             pco_dqs = max( wdt*PLScon_d(k,i,j), -clp_dqs )
             pco_dqg = -pco_dqi-pco_dqs
             ! 08/05/08 [Mod] T.Mitsui
             pco_dni = max( wdt*PNIcon_d(k,i,j), -clp_dni )
             pco_dns = max( wdt*PNScon_d(k,i,j), -clp_dns )
             pco_dng = -pco_dni-pco_dns
             ! enhanced melting ( always negative value )
             ! ice-cloud melting produces cloud, others produce rain
             eml_dqi =  max( wdt*PLIacm_d(k,i,j), min(0.0_RP, -li_d(k,i,j)-(clp_dqi+clm_dqi)-pco_dqi ))
             eml_dqs =  max( wdt*PLSacm_d(k,i,j), min(0.0_RP, -ls_d(k,i,j)-(clp_dqs+clm_dqs)-pco_dqs ))
             eml_dqg =  max( wdt*(PLGacm_d(k,i,j)+PLGarm_d(k,i,j)+PLSarm_d(k,i,j)+PLIarm_d(k,i,j)), &
                  min(0.0_RP, -lg_d(k,i,j)-(clp_dqg+clm_dqg)-pco_dqg ))
             eml_dqc = -eml_dqi
             eml_dqr = -eml_dqs-eml_dqg
             !
             eml_dni =  max( wdt*PNIacm_d(k,i,j), min(0.0_RP, -ni_d(k,i,j)-(clp_dni+clm_dni)-pco_dni ))
             eml_dns =  max( wdt*PNSacm_d(k,i,j), min(0.0_RP, -ns_d(k,i,j)-(clp_dns+clm_dns)-pco_dns ))
             eml_dng =  max( wdt*(PNGacm_d(k,i,j)+PNGarm_d(k,i,j)+PNSarm_d(k,i,j)+PNIarm_d(k,i,j)), &
                  min(0.0_RP, -ng_d(k,i,j)-(clp_dng+clm_dng)-pco_dng ))
             eml_dnc = -eml_dni
             eml_dnr = -eml_dns-eml_dng
             !
             ! ice multiplication
             spl_dqg = max( wdt*PLGspl_d(k,i,j), min(0.0_RP, -lg_d(k,i,j)-(clp_dqg+clm_dqg)-pco_dqg-eml_dqg ))
             spl_dqs = max( wdt*PLSspl_d(k,i,j), min(0.0_RP, -ls_d(k,i,j)-(clp_dqs+clm_dqs)-pco_dqs-eml_dqs ))
             spl_dqi = -spl_dqg-spl_dqs
             fac9    = (spl_dqg-eps)/(wdt*PLGspl_d(k,i,j)-eps)
             fac10   = (spl_dqs-eps)/(wdt*PLSspl_d(k,i,j)-eps)
             spl_dni = wdt*PNIspl_d(k,i,j)*fac9*fac10
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
             rhogq_d(k,i,j,I_QC) = max(0.0_RP, rhogq_d(k,i,j,I_QC) + drhogqc )
             rhogq_d(k,i,j,I_NC) = max(0.0_RP, rhogq_d(k,i,j,I_NC) + drhognc )
             rhogq_d(k,i,j,I_QR) = max(0.0_RP, rhogq_d(k,i,j,I_QR) + drhogqr )
             rhogq_d(k,i,j,I_NR) = max(0.0_RP, rhogq_d(k,i,j,I_NR) + drhognr )
             rhogq_d(k,i,j,I_QI) = max(0.0_RP, rhogq_d(k,i,j,I_QI) + drhogqi )
             rhogq_d(k,i,j,I_NI) = max(0.0_RP, rhogq_d(k,i,j,I_NI) + drhogni )
             rhogq_d(k,i,j,I_QS) = max(0.0_RP, rhogq_d(k,i,j,I_QS) + drhogqs )
             rhogq_d(k,i,j,I_NS) = max(0.0_RP, rhogq_d(k,i,j,I_NS) + drhogns )
             rhogq_d(k,i,j,I_QG) = max(0.0_RP, rhogq_d(k,i,j,I_QG) + drhogqg )
             rhogq_d(k,i,j,I_NG) = max(0.0_RP, rhogq_d(k,i,j,I_NG) + drhogng )
             !
             ! update
             ! rhogq = l*gsgam
             rhoge_d(k,i,j) = rhoge_d(k,i,j) + LHF * ( drhogqi + drhogqs + drhogqg )
          enddo
       enddo
       enddo

    !--- update mixing ratio
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       do iq = 1,  QA
          q_d(k,i,j,iq) = rhogq_d(k,i,j,iq) * rrhog_d(k,i,j)
       enddo

       CALC_QDRY( qd_d(k,i,j), q_d, k, i, j, iq )
       CALC_CV( cva_d(k,i,j), qd_d(k,i,j), q_d, k, i, j, iq, CVdry, AQ_CV )
       CALC_R( Rmoist, q_d(k,i,j,I_QV), qd_d(k,i,j), Rdry, Rvap )

       tem_d(k,i,j) = rhoge_d(k,i,j) / ( DENS(k,i,j) * cva_d(k,i,j) )
       pre_d(k,i,j) = DENS(k,i,j) * Rmoist * tem_d(k,i,j)
    enddo
    enddo
    enddo

!    if( opt_debug )     call debugreport_collection
    if( opt_debug_tem ) call debug_tem_kij( 4, tem_d(:,:,:), DENS(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call PROF_rapend  ('MP3 Collection')

    call PROF_rapstart('MPX ijkconvert')

    do j = JS, JE
    do i = IS, IE

       do k  = KS, KE
          CALC_CP( cpa_d(k,i,j), qd_d(k,i,j), q_d, k, i, j, iq, CPdry, AQ_CP )
          CALC_R( Rmoist, q_d(k,i,j,I_QV), qd_d(k,i,j), Rdry, Rvap )
          RHOT(k,i,j) = tem_d(k,i,j) * ( P00 / pre_d(k,i,j) )**(Rmoist/cpa_d(k,i,j)) &
               * DENS(k,i,j)
       enddo
    enddo
    enddo

    do iq = 1, QA
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = rhogq_d(k,i,j,iq) * rrhog_d(k,i,j)
       enddo
       enddo
       enddo
    enddo

    call PROF_rapend  ('MPX ijkconvert')

    !----------------------------------------------------------------------------
    !
    ! 4.Saturation adjustment
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP4 Saturation adjustment')
    ! nothing to do
    call PROF_rapend  ('MP4 Saturation adjustment')
    !----------------------------------------------------------------------------
    !
    ! 5. Sedimentation ( terminal velocity must be negative )
    !
    !----------------------------------------------------------------------------
    call PROF_rapstart('MP5 Sedimentation')

    if ( doprecipitation ) then

    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE
       flux_rain(k,i,j) = 0.0_RP
       flux_snow(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do step = 1, MP_NSTEP_SEDIMENTATION

       call MP_terminal_velocity( velw(:,:,:,:),    &
                                  rhogq_d(:,:,:,:), &
                                  DENS(:,:,:),      &
                                  tem_d(:,:,:),     &
                                  pre_d(:,:,:)      )

!       call precipitation( wflux_rain(:,:,:),     &
!                           wflux_snow(:,:,:),     &
!                           velw(:,:,:,:),         &
!                           rhogq_d(:,:,:,:),      &
!                           rhoge_d(:,:,:),        &
!                           tem_d(:,:,:),          &
!                           MP_DTSEC_SEDIMENTATION )
       call MP_precipitation( &
            wflux_rain, wflux_snow, &
            DENS, MOMZ, MOMX, MOMY, &
            rhoge_d, QTRC, &
            velw, tem_d, &
            MP_DTSEC_SEDIMENTATION )

       do j = JS, JE
       do i = IS, IE
          do k = KS-1, KE
             flux_rain(k,i,j) = flux_rain(k,i,j) + wflux_rain(k,i,j) * MP_RNSTEP_SEDIMENTATION
             flux_snow(k,i,j) = flux_snow(k,i,j) + wflux_snow(k,i,j) * MP_RNSTEP_SEDIMENTATION
          enddo
          flux_prec(i,j) = flux_rain(KS-1,i,j) + flux_snow(KS-1,i,j)
       enddo
       enddo

!       if( opt_debug ) call debugreport_sedimentation

    enddo

    endif

    call HIST_in( flux_rain(KS-1,:,:), 'RAIN', 'surface rain rate', 'kg/m2/s', dt)
    call HIST_in( flux_snow(KS-1,:,:), 'SNOW', 'surface snow rate', 'kg/m2/s', dt)
    call HIST_in( flux_prec(:,:),      'PREC', 'surface precipitaion rate', 'kg/m2/s', dt)

    call PROF_rapend  ('MP5 Sedimentation')

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
       IA, JA, KA,      &
       IS, IE,       &
       JS, JE,       &
       KS, KE,       &
       z, w,             &
       rho, tem, pre,    &
       LV, NC, NI,       &
       PNCccn, PLCccn,   &
       PNIccn, PLIccn,   &
       cpa,              & ! in
       dTdt_rad,         & ! in
       qke,              & ! in
       flag_history_in,  & ! in
       dt                ) ! in
    use scale_atmos_saturation, only: &
       moist_psat_liq       => ATMOS_SATURATION_psat_liq, &
       moist_psat_ice       => ATMOS_SATURATION_psat_ice,   &
       moist_pres2qsat_liq  => ATMOS_SATURATION_pres2qsat_liq, &
       moist_pres2qsat_ice  => ATMOS_SATURATION_pres2qsat_ice,   &
       moist_dqsi_dtem_rho  => ATMOS_SATURATION_dqsi_dtem_rho
    implicit none

    integer, intent(in)  :: KA, IA, JA
    integer, intent(in)  :: KS, IS, JS
    integer, intent(in)  :: KE, IE, JE
    !
    real(RP), intent(in)  :: z(KA)      !
    real(RP), intent(in)  :: rho(KA,IA,JA)    ! [Add] 09/08/18 T.Mitsui
    real(RP), intent(in)  :: tem(KA,IA,JA)    ! [Add] 09/08/18 T.Mitsui
    real(RP), intent(in)  :: pre(KA,IA,JA)    ! [Add] 09/08/18 T.Mitsui
    real(RP), intent(in)  :: w(KA,IA,JA)      ! w of half point
    !
    real(RP), intent(in)  :: LV(KA,IA,JA)     !
    real(RP), intent(in)  :: NC(KA,IA,JA)     ! [Add] 09/04/14 T.Mitsui
    real(RP), intent(in)  :: NI(KA,IA,JA)     !
    real(RP), intent(out) :: PNCccn(KA,IA,JA) !
    real(RP), intent(out) :: PLCccn(KA,IA,JA) !
    real(RP), intent(out) :: PNIccn(KA,IA,JA) !
    real(RP), intent(out) :: PLIccn(KA,IA,JA)
    !
    real(RP), intent(in) ::  cpa(KA,IA,JA)      ! in  09/08/18 [Add] T.Mitsui
    real(RP), intent(in)  :: dTdt_rad(KA,IA,JA) ! 09/08/18 T.Mitsui
    real(RP), intent(in)  :: qke(KA,IA,JA)      ! 09/08/18 T.Mitsui
    real(RP), intent(in)  :: dt
    logical, intent(in)  :: flag_history_in      ! in 10/08/03 [Add] T.Mitsui
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
100    if( IO_L ) write(IO_FID_LOG, nml=nm_mp_sn14_nucleation)
       flag_first=.false.
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
    !
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          pv        = LV(k,i,j)*Rvap*tem(k,i,j)
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
       coef_ccn(i,j)  = 1.E+6_RP*0.88_RP*(c_ccn*1.E-6_RP)**(2.0_RP/(kappa + 2.0_RP)) * &
!            (70.0_RP)**(kappa_map(1,i,j)/(kappa_map(1,i,j) + 2.0_RP))
            (70.0_RP)**(kappa/(kappa + 2.0_RP))
!       slope_ccn(i,j) = 1.5_RP*kappa_map(1,i,j)/(kappa_map(1,i,j) + 2.0_RP)
       slope_ccn(i,j) = 1.5_RP*kappa/(kappa + 2.0_RP)
       !
       do k=KS, KE
          sigma_w(k,i,j) = r_sqrt3*sqrt(max(qke(k,i,j),qke_min))
       end do
       sigma_w(KS-1,i,j) = sigma_w(KS,i,j)
       sigma_w(KE+1,i,j) = sigma_w(KE,i,j)
       ! effective vertical velocity
       do k=KS, KE-1
          weff(k,i,j) = 0.5_RP*(w(k,i,j) + w(k+1,i,j)) - cpa(k,i,j)*r_gravity*dTdt_rad(k,i,j)
       end do
       weff(KS-1,i,j) = weff(KS,i,j)
       weff(KE,i,j)   = weff(KE-1,i,j)
       !
       if( nucl_twomey ) then
       ! diagnose cloud condensation nuclei
        do k=KS, KE
           ! effective vertical velocity (maximum vertical velocity in turbulent flow)
           weff_max(k,i,j) = weff(k,i,j) + sigma_w(k,i,j)
           ! large scale upward motion region and saturated
           if( (weff(k,i,j) > 1.E-8_RP) .and. (ssw(k,i,j) > 1.E-10_RP)  .and. pre(k,i,j) > 300.E+2_RP )then
              ! Lohmann (2002), eq.(1)
              nc_new_max   = coef_ccn(i,j)*weff_max(k,i,j)**slope_ccn(i,j)
              nc_new(k,i,j) = a_max*nc_new_max**b_max
           else
              nc_new(k,i,j) = 0.0_RP
           end if
        end do
       else
        ! calculate cloud condensation nuclei
        do k=KS, KE
         if( ssw(k,i,j) > 1.e-10_RP .and. pre(k,i,j) > 300.E+2_RP ) then
           nc_new(k,i,j) = c_ccn*ssw(k,i,j)**kappa
         else
           nc_new(k,i,j) = 0.0_RP
         endif
        enddo
       endif
       !
       do k=KS, KE
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
       ! search maximum value of nc_new
!       do k=KS, KE
!          if(  ( nc_new(k,i,j) < nc_new_below(k,i,j) ) .or. &
!               ( nc_new_below(k,i,j) > c_ccn_map(1,i,j)*0.05_RP ) )then ! 5% of c_ccn
!               ( nc_new_below(k,i,j) > c_ccn*0.05_RP ) )then ! 5% of c_ccn
!             flag_nucleation(k,i,j) = .false.
!          end if
!       end do
       if( nucl_twomey ) then
        ! nucleation occurs at only cloud base.
        ! if CCN is more than below parcel, nucleation newly occurs
        do k=KS, KE
          ! effective vertical velocity
          if(   flag_nucleation(k,i,j)               .and. & ! large scale upward motion region and saturated
               ( tem(k,i,j)    > tem_ccn_low       ) .and. &
               ( nc_new(k,i,j) > NC(k,i,j) )                )then
             dlcdt_max    = (LV(k,i,j) - esw(k,i,j)/(Rvap*tem(k,i,j)))*rdt
             dncdt_max    = dlcdt_max/xc_min
             dnc_new      = nc_new(k,i,j)-NC(k,i,j)
             PNCccn(k,i,j) = min( dncdt_max, dnc_new*rdt )
             PLCccn(k,i,j) = min( dlcdt_max, xc_min*PNCccn(k,i,j) )
          else
             PNCccn(k,i,j) = 0.0_RP
             PLCccn(k,i,j) = 0.0_RP
          end if
        end do
       else
        do k=KS, KE
          ! effective vertical velocity
          if(  ( tem(k,i,j)    > tem_ccn_low       ) .and. &
               ( nc_new(k,i,j) > NC(k,i,j) )                )then
             dlcdt_max    = (LV(k,i,j) - esw(k,i,j)/(Rvap*tem(k,i,j)))*rdt
             dncdt_max    = dlcdt_max/xc_min
             dnc_new      = nc_new(k,i,j)-NC(k,i,j)
             PNCccn(k,i,j) = min( dncdt_max, dnc_new*rdt )
             PLCccn(k,i,j) = min( dlcdt_max, xc_min*PNCccn(k,i,j) )
          else
             PNCccn(k,i,j) = 0.0_RP
             PLCccn(k,i,j) = 0.0_RP
          end if
        end do
       endif
       !
       ! ice nucleation
       !
       ! +++ NOTE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! Based on Phillips etal.(2006).
       ! However this approach doesn't diagnose Ni itself but diagnose tendency.
       ! Original approach adjust Ni instantaneously .
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       do k=KS, KE
          dz             = z(k) - z_below(k,i,j)
          w_dssidz(k,i,j) = w(k,i,j)*(ssi(k,i,j) - ssi_below(k,i,j))/dz ! 09/04/14 [Add] T.Mitsui
          dssidt_rad(k,i,j) = -LV(k,i,j)/(rho(k,i,j)*qsi(k,i,j)*qsi(k,i,j))*dqsidtem_rho(k,i,j)*dTdt_rad(k,i,j)
          dli_max        = (LV(k,i,j) - esi(k,i,j)/(Rvap*tem(k,i,j)))*rdt
          dni_max        = min( dli_max/xi_ccn, (in_max-NI(k,i,j))*rdt )
          wdssi          = min( w_dssidz(k,i,j)+dssidt_rad(k,i,j), 0.01_RP)
          wssi           = min( ssi(k,i,j), ssi_max)
          ! SB06(34),(35)
          if(  ( wdssi    > eps        ) .and. & !
               (tem(k,i,j) < 273.15_RP   ) .and. & !
               (NI(k,i,j)  < in_max     ) .and. &
               (wssi      >= eps       ) )then   !
!             PNIccn(k,i,j) = min(dni_max, c_in_map(1,i,j)*bm_M92*nm_M92*0.3_RP*exp(0.3_RP*bm_M92*(wssi-0.1_RP))*wdssi)
             if( inucl_w ) then
              PNIccn(k,i,j) = min(dni_max, c_in*bm_M92*nm_M92*0.3_RP*exp(0.3_RP*bm_M92*(wssi-0.1_RP))*wdssi)
             else
                PNIccn(k,i,j) = min(dni_max, max(c_in*nm_M92*exp(0.3_RP*bm_M92*(wssi-0.1_RP) )-NI(k,i,j),0.0_RP )*rdt )
             endif
             PLIccn(k,i,j) = min(dli_max, PNIccn(k,i,j)*xi_ccn )
             ! only for output
!             dni_ratio(k,i,j) = dssidt_rad(k,i,j)/( w_dssidz(k,i,j)+dssidt_rad(k,i,j) )
          else
             PNIccn(k,i,j) = 0.0_RP
             PLIccn(k,i,j) = 0.0_RP
          end if
       end do

    end do
    end do
    !
    return
  end subroutine nucleation_kij
  !----------------------------
  subroutine ice_multiplication_kij( &
       IA, JA, KA,              & ! in
       IS, IE,               & ! in
       JS, JE,               & ! in
       KS, KE,               & ! in
       PLGspl, PLSspl, PNIspl,   & ! out
       PLIacLC2LI,               & ! in
       PLSacLC2LS,               & ! in
       PLGacLC2LG,               & ! in
       tem, LC, nc, xi, xc       ) ! in

    ! ice multiplication by splintering
    ! we consider Hallet-Mossop process
    use scale_specfunc, only: &
       gammafunc => SF_gamma
    implicit none
    !
    integer, intent(in) :: IA, JA, KA
    integer, intent(in) :: IS, JS, KS
    integer, intent(in) :: IE, JE, KE
    real(RP), intent(out):: PLGspl(KA,IA,JA)
    real(RP), intent(out):: PLSspl(KA,IA,JA) ! [Add]
    real(RP), intent(out):: PNIspl(KA,IA,JA)
    !
    real(RP), intent(in) :: PLIacLC2LI(KA,IA,JA)
    real(RP), intent(in) :: PLSacLC2LS(KA,IA,JA)
    real(RP), intent(in) :: PLGacLC2LG(KA,IA,JA)
    real(RP), intent(in) :: tem(KA,IA,JA)
    real(RP), intent(in) :: LC(KA,IA,JA)
    real(RP), intent(in) :: NC(KA,IA,JA)
    real(RP), intent(in) :: xi(KA,IA,JA)
    real(RP), intent(in) :: xc(KA,IA,JA)
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
       xc_cr = (2.0_RP*rc_cr/a_m(I_QC))**(1.0_RP/b_m(I_QC))
       alpha = (nu(I_QC)+1.0_RP)/mu(I_QC)
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
          x       = coef_lambda(I_QC)*(xc_cr/xc(k,i,j))**mu(I_QC)
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
          n12 = NC(k,i,j)*(1.0_RP-igm)
          ! eq.(82) CT(86)
          wn           = (pice + n12/((LC(k,i,j)+xc_min)*pnc) )*fp ! filtered by xc_min
          wni          = wn*(-PLIacLC2LI(k,i,j)  ) ! riming production rate is all negative
          wns          = wn*(-PLSacLC2LS(k,i,j)  )
          wng          = wn*(-PLGacLC2LG(k,i,j)  )
          PNIspl(k,i,j) = wni+wns+wng
          !
          PLSspl(k,i,j) = - wns*xi(k,i,j) ! snow    => ice
          PLGspl(k,i,j) = - wng*xi(k,i,j) ! graupel => ice
          !
       end do
    end do
    end do
    !
    return
  end subroutine ice_multiplication_kij
  !----------------------------
  subroutine mixed_phase_collection_kij(   &
       IA, JA, KA,                    & ! in
       IS, IE,                     & ! in
       JS, JE,                     & ! in
       KS, KE,                     & ! in
       ! collection process
       PLIacLC2LI,   PNIacNC2NI,       & ! out
       PLSacLC2LS,   PNSacNC2NS,       & ! out
       !
       PLGacLC2LG,   PNGacNC2NG,       & ! out
       PLRacLI2LG_I, PNRacNI2NG_I,     & ! out
       PLRacLI2LG_R, PNRacNI2NG_R,     & ! out
       PLRacLS2LG_S, PNRacNS2NG_S,     & ! out
       PLRacLS2LG_R, PNRacNS2NG_R,     & ! out
       PLRacLG2LG,   PNRacNG2NG,       & ! out
       PLIacLI2LS,   PNIacNI2NS,       & ! out
       PLIacLS2LS,   PNIacNS2NS,       & ! out
       PNSacNS2NS,                     & ! out
       PNGacNG2NG,                     & ! out [ADD] 11/08/30 T.Mitsui
       PLGacLS2LG,   PNGacNS2NG,       & ! out
       ! partial conversion (ice, snow => graupel)
       PLIcon, PNIcon, PLScon, PNScon, & ! out
       ! enhanced melting
       PLIacm, PNIacm, PLIarm, PNIarm, & ! out
       PLSacm, PNSacm, PLSarm, PNSarm, & ! out
       PLGacm, PNGacm, PLGarm, PNGarm, & ! out
       wtem,                           & ! in
       LV,                             & ! in [Add] 10/08/03 T.Mitsui
       LC, LR, LI, LS, LG,             & ! in
       NC, NR, NI, NS, NG,             & ! in
       xc, xr, xi, xs, xg,             & ! in
       dc_xave, dr_xave, di_xave, ds_xave, dg_xave,& ! in
       vtc_xave, vtr_xave, vti_xave, vts_xave, vtg_xave, &
!       rho, & ! [Add] 11/08/30
       flag_history_in                 ) ! in [Add] 11/08/30
    use scale_atmos_saturation, only: &
       moist_psat_ice => ATMOS_SATURATION_psat_ice
    implicit none

    integer, intent(in) :: IA,JA,KA
    integer, intent(in) :: IS,JS,KS
    integer, intent(in) :: IE,JE,KE
    !--- mixed-phase collection process
    !    PXXacYY2ZZ_x: XX collecting YY to form ZZ.
    !                  _x means the contributions of XX or YY.
    !                  And all we set all production term as a negative sign to avoid confusion.
    !
    !--- ice-cloud     => ice
    real(RP), intent(out):: PLIacLC2LI(KA,IA,JA) ! mass
    real(RP), intent(out):: PNIacNC2NI(KA,IA,JA) ! number
    !--- snow-cloud    => snow
    real(RP), intent(out):: PLSacLC2LS(KA,IA,JA) ! reduction of cloud
    real(RP), intent(out):: PNSacNC2NS(KA,IA,JA) !
    !--- graupel-cloud => graupel
    real(RP), intent(out):: PLGacLC2LG(KA,IA,JA)
    real(RP), intent(out):: PNGacNC2NG(KA,IA,JA)
    !--- rain-ice      => graupel
    real(RP), intent(out):: PLRacLI2LG_R(KA,IA,JA) ! reduction of rain
    real(RP), intent(out):: PNRacNI2NG_R(KA,IA,JA) !
    real(RP), intent(out):: PLRacLI2LG_I(KA,IA,JA) ! reduction of ice
    real(RP), intent(out):: PNRacNI2NG_I(KA,IA,JA) !
    !--- rain-snow     => graupel
    real(RP), intent(out):: PLRacLS2LG_R(KA,IA,JA) ! reduction of rain
    real(RP), intent(out):: PNRacNS2NG_R(KA,IA,JA) !
    real(RP), intent(out):: PLRacLS2LG_S(KA,IA,JA) ! reduction of snow
    real(RP), intent(out):: PNRacNS2NG_S(KA,IA,JA) !
    !--- rain-graupel  => graupel
    real(RP), intent(out):: PLRacLG2LG(KA,IA,JA) ! reduction of rain
    real(RP), intent(out):: PNRacNG2NG(KA,IA,JA) ! reduction of graupel
    !--- ice-ice     => snow
    real(RP), intent(out):: PLIacLI2LS(KA,IA,JA)
    real(RP), intent(out):: PNIacNI2NS(KA,IA,JA)
    !--- ice-snow     => snow
    real(RP), intent(out):: PLIacLS2LS(KA,IA,JA)
    real(RP), intent(out):: PNIacNS2NS(KA,IA,JA)
    !--- snow-snow     => snow
    real(RP), intent(out):: PNSacNS2NS(KA,IA,JA)
    !--- graupel-graupel=> graupel
    real(RP), intent(out):: PNGacNG2NG(KA,IA,JA)
    !--- graupel-snow     => graupel
    real(RP), intent(out):: PLGacLS2LG(KA,IA,JA)
    real(RP), intent(out):: PNGacNS2NG(KA,IA,JA)
    !--- partial conversion
    !--- ice-cloud => graupel
    real(RP), intent(out):: PLIcon(KA,IA,JA)
    real(RP), intent(out):: PNIcon(KA,IA,JA)
    !--- snow-cloud => graupel
    real(RP), intent(out):: PLScon(KA,IA,JA)
    real(RP), intent(out):: PNScon(KA,IA,JA)
    !--- enhanced melting
    !--- graupel-cloud melting => rain
    real(RP), intent(out):: PLGacm(KA,IA,JA)
    real(RP), intent(out):: PNGacm(KA,IA,JA)
    !--- graupel-rain melting => rain
    real(RP), intent(out):: PLGarm(KA,IA,JA)
    real(RP), intent(out):: PNGarm(KA,IA,JA)
    !--- snow-cloud melting => rain
    real(RP), intent(out):: PLSacm(KA,IA,JA)
    real(RP), intent(out):: PNSacm(KA,IA,JA)
    !--- snow-rain melting => rain
    real(RP), intent(out):: PLSarm(KA,IA,JA)
    real(RP), intent(out):: PNSarm(KA,IA,JA)
    !--- ice-cloud melting => cloud ?
    real(RP), intent(out):: PLIacm(KA,IA,JA)
    real(RP), intent(out):: PNIacm(KA,IA,JA)
    !--- ice-rain melting => rain
    real(RP), intent(out):: PLIarm(KA,IA,JA)
    real(RP), intent(out):: PNIarm(KA,IA,JA)
    !
    real(RP), intent(in) :: wtem(KA,IA,JA)
    !--- mass concentration[kg/m3]
    real(RP), intent(in) :: LV(KA,IA,JA)
    real(RP), intent(in) :: LC(KA,IA,JA)
    real(RP), intent(in) :: LR(KA,IA,JA)
    real(RP), intent(in) :: LI(KA,IA,JA)
    real(RP), intent(in) :: LS(KA,IA,JA)
    real(RP), intent(in) :: LG(KA,IA,JA)
    !--- number concentration[/m3]
    real(RP), intent(in) :: NC(KA,IA,JA)
    real(RP), intent(in) :: NR(KA,IA,JA)
    real(RP), intent(in) :: NI(KA,IA,JA)
    real(RP), intent(in) :: NS(KA,IA,JA)
    real(RP), intent(in) :: NG(KA,IA,JA)
    ! necessary ?
    real(RP), intent(in) :: xc(KA,IA,JA) ! LC/NC
    real(RP), intent(in) :: xr(KA,IA,JA) ! LR/NR
    real(RP), intent(in) :: xi(KA,IA,JA) ! LI/NI
    real(RP), intent(in) :: xs(KA,IA,JA) ! LS/NS
    real(RP), intent(in) :: xg(KA,IA,JA) ! LG/NG
    !--- diameter of averaged mass( D(ave x) )
    real(RP), intent(in) :: dc_xave(KA,IA,JA)
    real(RP), intent(in) :: dr_xave(KA,IA,JA)
    real(RP), intent(in) :: di_xave(KA,IA,JA)
    real(RP), intent(in) :: ds_xave(KA,IA,JA)
    real(RP), intent(in) :: dg_xave(KA,IA,JA)
    !--- terminal velocity of averaged mass( vt(ave x) )
    real(RP), intent(in) :: vtc_xave(KA,IA,JA)
    real(RP), intent(in) :: vtr_xave(KA,IA,JA)
    real(RP), intent(in) :: vti_xave(KA,IA,JA)
    real(RP), intent(in) :: vts_xave(KA,IA,JA)
    real(RP), intent(in) :: vtg_xave(KA,IA,JA)
    ! [Add] 11/08/30 T.Mitsui, for autoconversion of ice
!    real(RP), intent(in) :: rho(KA,IA,JA)
    logical, intent(in) :: flag_history_in
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
    !
    integer :: i, j, k
    !
    if( flag_first )then
       rewind( IO_FID_CONF )
       read( IO_FID_CONF, nml=nm_mp_sn14_collection, end=100 )
100    if( IO_L ) write( IO_FID_LOG, nml=nm_mp_sn14_collection )
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
             esi_rat       = LV(k,i,j)*Rvap*tem(k,i,j)/esi(k,i,j)
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
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          !
!          temc_m = min(tem(k,i,j) - T00,0.0_RP) ! T < 273.15
          temc_p = max(tem(k,i,j) - T00,0.0_RP) ! T > 273.15
          ! averaged diameter using SB06(82)
          ave_dc = coef_d(I_QC)*xc(k,i,j)**b_m(I_QC)
          ave_di = coef_d(I_QI)*xi(k,i,j)**b_m(I_QI)
          ave_ds = coef_d(I_QS)*xs(k,i,j)**b_m(I_QS)
          ave_dg = coef_d(I_QG)*xg(k,i,j)**b_m(I_QG)
          !------------------------------------------------------------------------
          ! coellection efficiency are given as follows
          E_c = max(0.0_RP, min(1.0_RP, (ave_dc-dc0)/(dc1-dc0) ))
          if(ave_di>di0)then
             E_i = E_im
          else
             E_i = 0.0_RP
          end if
          if(ave_ds>ds0)then
             E_s = E_sm
          else
             E_s = 0.0_RP
          end if
          if(ave_dg>dg0)then
             E_g = E_gm
          else
             E_g = 0.0_RP
          end if
          E_ic = E_i*E_c
          E_sc = E_s*E_c
          E_gc = E_g*E_c
          !------------------------------------------------------------------------
          ! Collection:  a collects b ( assuming particle size a>b )
          dcdc = dc_xave(k,i,j) * dc_xave(k,i,j)
          drdr = dr_xave(k,i,j) * dr_xave(k,i,j)
          didi = di_xave(k,i,j) * di_xave(k,i,j)
          dsds = ds_xave(k,i,j) * ds_xave(k,i,j)
          dgdg = dg_xave(k,i,j) * dg_xave(k,i,j)
          dcdi = dc_xave(k,i,j) * di_xave(k,i,j)
          dcds = dc_xave(k,i,j) * ds_xave(k,i,j)
          dcdg = dc_xave(k,i,j) * dg_xave(k,i,j)
          drdi = dr_xave(k,i,j) * di_xave(k,i,j)
          drds = dr_xave(k,i,j) * ds_xave(k,i,j)
          drdg = dr_xave(k,i,j) * dg_xave(k,i,j)
          dids = di_xave(k,i,j) * ds_xave(k,i,j)
!          didg = di_xave(k,i,j) * dg_xave(k,i,j)
          dsdg = ds_xave(k,i,j) * dg_xave(k,i,j)
          !
          vcvc = vtc_xave(k,i,j)* vtc_xave(k,i,j)
          vrvr = vtr_xave(k,i,j)* vtr_xave(k,i,j)
          vivi = vti_xave(k,i,j)* vti_xave(k,i,j)
          vsvs = vts_xave(k,i,j)* vts_xave(k,i,j)
          vgvg = vtg_xave(k,i,j)* vtg_xave(k,i,j)
          vcvi = vtc_xave(k,i,j)* vti_xave(k,i,j)
          vcvs = vtc_xave(k,i,j)* vts_xave(k,i,j)
          vcvg = vtc_xave(k,i,j)* vtg_xave(k,i,j)
          vrvi = vtr_xave(k,i,j)* vti_xave(k,i,j)
          vrvs = vtr_xave(k,i,j)* vts_xave(k,i,j)
          vrvg = vtr_xave(k,i,j)* vtg_xave(k,i,j)
          vivs = vti_xave(k,i,j)* vts_xave(k,i,j)
!          vivg = vti_xave(k,i,j)* vtg_xave(k,i,j)
          vsvg = vts_xave(k,i,j)* vtg_xave(k,i,j)
          !------------------------------------------------------------------------
          !
          !+++ pattern 1: a + b => a  (a>b)
          !                           (i-c, s-c, g-c, s-i, g-r, s-g)
          !------------------------------------------------------------------------
          ! cloud-ice => ice
          ! reduction term of cloud
          coef_acc_LCI = &
               (  delta_b1(I_QC)*dcdc + delta_ab1(I_QI,I_QC)*dcdi + delta_b0(I_QI)*didi ) &
               *( theta_b1(I_QC)*vcvc - theta_ab1(I_QI,I_QC)*vcvi + theta_b0(I_QI)*vivi   &
               +  sigma_i + sigma_c )**0.5_RP
          coef_acc_NCI = &
               (  delta_b0(I_QC)*dcdc + delta_ab0(I_QI,I_QC)*dcdi + delta_b0(I_QI)*didi ) &
               *( theta_b0(I_QC)*vcvc - theta_ab0(I_QI,I_QC)*vcvi + theta_b0(I_QI)*vivi   &
               +  sigma_i + sigma_c )**0.5_RP
          PLIacLC2LI(k,i,j)= -0.25_RP*pi*E_ic*NI(k,i,j)*LC(k,i,j)*coef_acc_LCI
          PNIacNC2NI(k,i,j)= -0.25_RP*pi*E_ic*NI(k,i,j)*NC(k,i,j)*coef_acc_NCI
          ! cloud-snow => snow
          ! reduction term of cloud
          coef_acc_LCS = &
               (  delta_b1(I_QC)*dcdc + delta_ab1(I_QS,I_QC)*dcds + delta_b0(I_QS)*dsds ) &
               *( theta_b1(I_QC)*vcvc - theta_ab1(I_QS,I_QC)*vcvs + theta_b0(I_QS)*vsvs   &
               +  sigma_s + sigma_c )**0.5_RP
          coef_acc_NCS = &
               (  delta_b0(I_QC)*dcdc + delta_ab0(I_QS,I_QC)*dcds + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QC)*vcvc - theta_ab0(I_QS,I_QC)*vcvs + theta_b0(I_QS)*vsvs   &
               +  sigma_s + sigma_c )**0.5_RP
          PLSacLC2LS(k,i,j)= -0.25_RP*pi*E_sc*NS(k,i,j)*LC(k,i,j)*coef_acc_LCS
          PNSacNC2NS(k,i,j)= -0.25_RP*pi*E_sc*NS(k,i,j)*NC(k,i,j)*coef_acc_NCS
          ! cloud-graupel => graupel
          ! reduction term of cloud
          coef_acc_LCG = &
               (  delta_b1(I_QC)*dcdc + delta_ab1(I_QG,I_QC)*dcdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b1(I_QC)*vcvc - theta_ab1(I_QG,I_QC)*vcvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_c )**0.5_RP
          coef_acc_NCG = &
               (  delta_b0(I_QC)*dcdc + delta_ab0(I_QG,I_QC)*dcdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b0(I_QC)*vcvc - theta_ab0(I_QG,I_QC)*vcvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_c )**0.5_RP
          PLGacLC2LG(k,i,j)= -0.25_RP*pi*E_gc*NG(k,i,j)*LC(k,i,j)*coef_acc_LCG
          PNGacNC2NG(k,i,j)= -0.25_RP*pi*E_gc*NG(k,i,j)*NC(k,i,j)*coef_acc_NCG
          ! snow-graupel => graupel
          coef_acc_LSG = &
               (  delta_b1(I_QS)*dsds + delta_ab1(I_QG,I_QS)*dsdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b1(I_QS)*vsvs - theta_ab1(I_QG,I_QS)*vsvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_s )**0.5_RP
          coef_acc_NSG = &
               (  delta_b0(I_QS)*dsds + delta_ab0(I_QG,I_QS)*dsdg + delta_b0(I_QG)*dgdg ) &
               ! [fix] T.Mitsui 08/05/08
               *( theta_b0(I_QS)*vsvs - theta_ab0(I_QG,I_QS)*vsvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_s )**0.5_RP
          PLGacLS2LG(k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_gs*NG(k,i,j)*LS(k,i,j)*coef_acc_LSG
          PNGacNS2NG(k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_gs*NG(k,i,j)*NS(k,i,j)*coef_acc_NSG
          !------------------------------------------------------------------------
          ! ice-snow => snow
          ! reduction term of ice
          coef_acc_LIS = &
               (  delta_b1(I_QI)*didi + delta_ab1(I_QS,I_QI)*dids + delta_b0(I_QS)*dsds ) &
               *( theta_b1(I_QI)*vivi - theta_ab1(I_QS,I_QI)*vivs + theta_b0(I_QS)*vsvs   &
               +  sigma_i + sigma_s )**0.5_RP
          coef_acc_NIS = &
               (  delta_b0(I_QI)*didi + delta_ab0(I_QS,I_QI)*dids + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QI)*vivi - theta_ab0(I_QS,I_QI)*vivs + theta_b0(I_QS)*vsvs   &
               +  sigma_i + sigma_s )**0.5_RP
          PLIacLS2LS(k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_si*NS(k,i,j)*LI(k,i,j)*coef_acc_LIS
          PNIacNS2NS(k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_si*NS(k,i,j)*NI(k,i,j)*coef_acc_NIS
          !
          if ( tem(k,i,j) <= T00 )then
             ! rain-graupel => graupel
             ! reduction term of rain
             coef_acc_LRG = &
                  (  delta_b1(I_QR)*drdr + delta_ab1(I_QG,I_QR)*drdg + delta_b0(I_QG)*dgdg ) &
                  *( theta_b1(I_QR)*vrvr - theta_ab1(I_QG,I_QR)*vrvg + theta_b0(I_QG)*vgvg   &
                  +  sigma_r + sigma_g )**0.5_RP
             PLRacLG2LG(k,i,j) = -0.25_RP*pi*E_gr*NG(k,i,j)*LR(k,i,j)*coef_acc_LRG
          else
             ! rain-graupel => rain
             ! reduction term of graupel
             coef_acc_LRG = &
                  (  delta_b1(I_QG)*dgdg + delta_ab1(I_QR,I_QG)*drdg + delta_b0(I_QR)*drdr ) &
                  *( theta_b1(I_QG)*vgvg - theta_ab1(I_QR,I_QG)*vrvg + theta_b0(I_QR)*vrvr   &
                  +  sigma_r + sigma_g )**0.5_RP
             PLRacLG2LG(k,i,j) = -0.25_RP*pi*E_gr*NR(k,i,j)*LG(k,i,j)*coef_acc_LRG
          end if
          coef_acc_NRG = &
               (  delta_b0(I_QR)*drdr + delta_ab0(I_QG,I_QR)*drdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b0(I_QR)*vrvr - theta_ab0(I_QG,I_QR)*vrvg + theta_b0(I_QG)*vgvg   &
               +  sigma_r + sigma_g )**0.5_RP
          PNRacNG2NG(k,i,j) = -0.25_RP*pi*E_gr*NG(k,i,j)*NR(k,i,j)*coef_acc_NRG
          !
          !------------------------------------------------------------------------
          !
          !+++ pattern 2: a + b => c  (a>b)
          !                           (r-i,r-s)
          !------------------------------------------------------------------------
          ! rain-ice => graupel
          ! reduction term of ice
          coef_acc_LRI_I = &
               (  delta_b1(I_QI)*didi + delta_ab1(I_QR,I_QI)*drdi + delta_b0(I_QR)*drdr ) &
               *( theta_b1(I_QI)*vivi - theta_ab1(I_QR,I_QI)*vrvi + theta_b0(I_QR)*vrvr   &
               +  sigma_r + sigma_i )**0.5_RP
          coef_acc_NRI_I = &
               (  delta_b0(I_QI)*didi + delta_ab0(I_QR,I_QI)*drdi + delta_b0(I_QR)*drdr ) &
               *( theta_b0(I_QI)*vivi - theta_ab0(I_QR,I_QI)*vrvi + theta_b0(I_QR)*vrvr   &
               +  sigma_r + sigma_i )**0.5_RP
          PLRacLI2LG_I(k,i,j)= -0.25_RP*pi*E_ir*NR(k,i,j)*LI(k,i,j)*coef_acc_LRI_I
          PNRacNI2NG_I(k,i,j)= -0.25_RP*pi*E_ir*NR(k,i,j)*NI(k,i,j)*coef_acc_NRI_I
          ! reduction term of rain
          coef_acc_LRI_R = &
               (  delta_b1(I_QR)*drdr + delta_ab1(I_QI,I_QR)*drdi + delta_b0(I_QI)*didi ) &
               *( theta_b1(I_QR)*vrvr - theta_ab1(I_QI,I_QR)*vrvi + theta_b0(I_QI)*vivi   &
               +  sigma_r + sigma_i )**0.5_RP
          coef_acc_NRI_R = &
               (  delta_b0(I_QR)*drdr + delta_ab0(I_QI,I_QR)*drdi + delta_b0(I_QI)*didi ) &
               *( theta_b0(I_QR)*vrvr - theta_ab0(I_QI,I_QR)*vrvi + theta_b0(I_QI)*vivi   &
               +  sigma_r + sigma_i )**0.5_RP
          PLRacLI2LG_R(k,i,j)= -0.25_RP*pi*E_ir*NI(k,i,j)*LR(k,i,j)*coef_acc_LRI_R
          PNRacNI2NG_R(k,i,j)= -0.25_RP*pi*E_ir*NI(k,i,j)*NR(k,i,j)*coef_acc_NRI_R
          ! rain-snow => graupel
          ! reduction term of snow
          coef_acc_LRS_S = &
               (  delta_b1(I_QS)*dsds + delta_ab1(I_QR,I_QS)*drds + delta_b0(I_QR)*drdr ) &
               *( theta_b1(I_QS)*vsvs - theta_ab1(I_QR,I_QS)*vrvs + theta_b0(I_QR)*vrvr   &
               +  sigma_r + sigma_s )**0.5_RP
          coef_acc_NRS_S = &
               (  delta_b0(I_QS)*dsds + delta_ab0(I_QR,I_QS)*drds + delta_b0(I_QR)*drdr ) &
               *( theta_b0(I_QS)*vsvs - theta_ab0(I_QR,I_QS)*vrvs + theta_b0(I_QR)*vrvr   &
               +  sigma_r + sigma_s )**0.5_RP
          PLRacLS2LG_S(k,i,j)= -0.25_RP*pi*E_sr*NR(k,i,j)*LS(k,i,j)*coef_acc_LRS_S
          PNRacNS2NG_S(k,i,j)= -0.25_RP*pi*E_sr*NR(k,i,j)*NS(k,i,j)*coef_acc_NRS_S
          ! reduction term of rain
          coef_acc_LRS_R = &
               (  delta_b1(I_QR)*drdr + delta_ab1(I_QS,I_QR)*drds + delta_b0(I_QS)*dsds ) &
               *( theta_b1(I_QR)*vrvr - theta_ab1(I_QS,I_QR)*vrvs + theta_b0(I_QS)*vsvs   &
               +  sigma_r + sigma_s )**0.5_RP
          coef_acc_NRS_R = &
               (  delta_b0(I_QR)*drdr + delta_ab0(I_QS,I_QR)*drds + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QR)*vrvr - theta_ab0(I_QS,I_QR)*vrvs + theta_b0(I_QS)*vsvs   &
               +  sigma_r + sigma_s )**0.5_RP
          PLRacLS2LG_R(k,i,j)= -0.25_RP*pi*E_sr*NS(k,i,j)*LR(k,i,j)*coef_acc_LRS_R
          PNRacNS2NG_R(k,i,j)= -0.25_RP*pi*E_sr*NS(k,i,j)*NR(k,i,j)*coef_acc_NRS_R
          !------------------------------------------------------------------------
          !
          !+++ pattern 3: a + a => b  (i-i)
          !
          !------------------------------------------------------------------------
          ! ice-ice ( reduction is double, but includes double-count)
          coef_acc_LII = &
               (  delta_b0(I_QI)*didi + delta_ab1(I_QI,I_QI)*didi + delta_b1(I_QI)*didi ) &
               *( theta_b0(I_QI)*vivi - theta_ab1(I_QI,I_QI)*vivi + theta_b1(I_QI)*vivi   &
               +  sigma_i + sigma_i )**0.5_RP
          coef_acc_NII = &
               (  delta_b0(I_QI)*didi + delta_ab0(I_QI,I_QI)*didi + delta_b0(I_QI)*didi ) &
               *( theta_b0(I_QI)*vivi - theta_ab0(I_QI,I_QI)*vivi + theta_b0(I_QI)*vivi   &
               +  sigma_i + sigma_i )**0.5_RP
          PLIacLI2LS(k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_ii*NI(k,i,j)*LI(k,i,j)*coef_acc_LII
          PNIacNI2NS(k,i,j)= -0.25_RP*pi*E_stick(k,i,j)*E_ii*NI(k,i,j)*NI(k,i,j)*coef_acc_NII
          !
!          ci_aut(k,i,j)   =  0.25_RP*pi*E_ii*NI(k,i,j)*coef_acc_LII
!          taui_aut(k,i,j) = 1._RP/max(E_stick(k,i,j)*ci_aut(k,i,j),1.E-10_RP)
!          tau_sce(k,i,j)  = LI(k,i,j)/max(LI(k,i,j)+LS(k,i,j),1.E-10_RP)
          !------------------------------------------------------------------------
          !
          !+++ pattern 4: a + a => a  (s-s)
          !
          !------------------------------------------------------------------------
          ! snow-snow => snow
          coef_acc_NSS = &
               (  delta_b0(I_QS)*dsds + delta_ab0(I_QS,I_QS)*dsds + delta_b0(I_QS)*dsds ) &
               *( theta_b0(I_QS)*vsvs - theta_ab0(I_QS,I_QS)*vsvs + theta_b0(I_QS)*vsvs   &
               +  sigma_s + sigma_s )**0.5_RP
          PNSacNS2NS(k,i,j)= -0.125_RP*pi*E_stick(k,i,j)*E_ss*NS(k,i,j)*NS(k,i,j)*coef_acc_NSS
          !
          ! graupel-grauple => graupel
          coef_acc_NGG = &
               (  delta_b0(I_QG)*dgdg + delta_ab0(I_QG,I_QG)*dgdg + delta_b0(I_QG)*dgdg ) &
               *( theta_b0(I_QG)*vgvg - theta_ab0(I_QG,I_QG)*vgvg + theta_b0(I_QG)*vgvg   &
               +  sigma_g + sigma_g )**0.5_RP
          PNGacNG2NG(k,i,j)= -0.125_RP*pi*E_stick(k,i,j)*E_gg*NG(k,i,j)*NG(k,i,j)*coef_acc_NGG
          !
          !------------------------------------------------------------------------
          !--- Partial conversion
          ! SB06(70),(71)
          ! i_iconv2g: option whether partial conversions work or not
          ! ice-cloud => graupel
          if( ave_di > di_cri )then
             wx_cri = cfill_i*DWATR/rho_g*( pi/6.0_RP*rho_g*ave_di*ave_di*ave_di/xi(k,i,j) - 1.0_RP )
             PLIcon(k,i,j) = i_iconv2g*  PLIacLC2LI(k,i,j)/max(1.0_RP, wx_cri)
             PNIcon(k,i,j) = i_iconv2g*  PLIcon(k,i,j)/xi(k,i,j)
          else
             wx_cri       = 0.0_RP
             PLIcon(k,i,j) = 0.0_RP
             PNIcon(k,i,j) = 0.0_RP
          end if
          ! snow-cloud => graupel
          wx_crs = cfill_s*DWATR/rho_g*( pi/6.0_RP*rho_g*ave_ds*ave_ds*ave_ds/xs(k,i,j) - 1.0_RP )
          PLScon(k,i,j) = i_sconv2g*  (PLSacLC2LS(k,i,j))/max(1.0_RP, wx_crs)
          PNScon(k,i,j) = i_sconv2g*  PLScon(k,i,j)/xs(k,i,j)
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
          PLGacm(k,i,j) =  coef_emelt*PLGacLC2LG(k,i,j)
          PNGacm(k,i,j) =  PLGacm(k,i,j)/xg(k,i,j)
          ! rain-graupel
          PLGarm(k,i,j) =  coef_emelt*PLRacLG2LG(k,i,j)
          PNGarm(k,i,j) =  PLGarm(k,i,j)/xg(k,i,j)
          ! cloud-snow
          PLSacm(k,i,j) =  coef_emelt*(PLSacLC2LS(k,i,j))
          PNSacm(k,i,j) =  PLSacm(k,i,j)/xs(k,i,j)
          ! rain-snow
          PLSarm(k,i,j) =  coef_emelt*(PLRacLS2LG_R(k,i,j)+PLRacLS2LG_S(k,i,j))
          PNSarm(k,i,j) =  PLSarm(k,i,j)/xg(k,i,j)
          ! cloud-ice
          PLIacm(k,i,j) =  coef_emelt*PLIacLC2LI(k,i,j)
          PNIacm(k,i,j) =  PLIacm(k,i,j)/xi(k,i,j)
          ! rain-ice
          PLIarm(k,i,j) =  coef_emelt*(PLRacLI2LG_R(k,i,j)+PLRacLI2LG_I(k,i,j))
          PNIarm(k,i,j) =  PLIarm(k,i,j)/xg(k,i,j)
       end do
    end do
    end do
    !
!!$    if ( flag_history_in )then
!!$       call history_in( 'ml_PLIaut' , PLIacLI2LS(:,:) )
!!$       call history_in( 'ml_PNIaut' , PNIacNI2NS(:,:) )
!!$       call history_in( 'ml_PLIacc' , PLIacLS2LS(:,:) )
!!$       call history_in( 'ml_PNIacc' , PNIacNS2NS(:,:) )
!!$       call history_in( 'ml_ci_aut' , ci_aut(:,:) )
!!$       call history_in( 'ml_ti_aut' , taui_aut(:,:) )
!!$       call history_in( 'ml_qi_aut' , LI(:,:)/rho(:,:) )
!!$       call history_in( 'ml_qs_aut' , LS(:,:)/rho(:,:) )
!!$       call history_in( 'ml_ni_aut' , NI(:,:) )
!!$       call history_in( 'ml_ns_aut' , NS(:,:) )
!!$       call history_in( 'ml_tem_aut' , tem(:,:) )
!!$       call history_in( 'ml_rho_aut' , rho(:,:) )
!!$       call history_in( 'ml_tau_sce' , tau_sce(:,:) )
!!$    end if
    !
    return
  end subroutine mixed_phase_collection_kij
  !----------------------------
  ! Auto-conversion, Accretion, Self-collection, Break-up
  subroutine aut_acc_slc_brk_kij(  &
       IA, JA, KA,            &
       IS, IE,             &
       JS, JE,             &
       KS, KE,             &
       PLCaut, PNCaut,         &
       PNRaut,                 &
       PLCacc, PNCacc,         &
       PNRslc,                 &
       PNRbrk,                 &
       LC, LR, NC, NR, xc,     &
       dr_xave,                &
       rho                     )
    implicit none

    integer, intent(in)  :: IA, JA, KA
    integer, intent(in)  :: IS, JS, KS
    integer, intent(in)  :: IE, JE, KE
    !
    real(RP), intent(out) :: PLCaut(KA,IA,JA) ! Lc change for Auto-conversion
    real(RP), intent(out) :: PNCaut(KA,IA,JA) ! Nc
    real(RP), intent(out) :: PNRaut(KA,IA,JA) ! Nr
    real(RP), intent(out) :: PLCacc(KA,IA,JA) ! Lc change for Accretion
    real(RP), intent(out) :: PNCacc(KA,IA,JA) ! Nc
    real(RP), intent(out) :: PNRslc(KA,IA,JA) ! Nr change for Self-collection
    real(RP), intent(out) :: PNRbrk(KA,IA,JA) ! Nr change for Breakup
    !
    real(RP), intent(in)  :: LC(KA,IA,JA)
    real(RP), intent(in)  :: LR(KA,IA,JA)
    real(RP), intent(in)  :: NC(KA,IA,JA)
    real(RP), intent(in)  :: NR(KA,IA,JA)
    real(RP), intent(in)  :: xc(KA,IA,JA)
    real(RP), intent(in)  :: dr_xave(KA,IA,JA)
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
    coef_nuc0 = (nu(I_QC)+2.0_RP)/(nu(I_QC)+1.0_RP)
    coef_nuc1 = (nu(I_QC)+2.0_RP)*(nu(I_QC)+4.0_RP)/(nu(I_QC)+1.0_RP)/(nu(I_QC)+1.0_RP)
    coef_aut0 =  -kcc*coef_nuc0
    coef_aut1 =  -kcc/x_sep/20._RP*coef_nuc1
    !
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          lwc = LR(k,i,j)+LC(k,i,j)
          if( lwc > xc_min )then
             tau  = max(tau_min, LR(k,i,j)/lwc)
          else
             tau  = tau_min
          end if
          rho_fac = sqrt(rho_0/max(rho(k,i,j),rho_min))
          !
          ! Auto-conversion ( cloud-cloud => rain )
          psi_aut       = 400._RP*(tau**0.7_RP)*(1.0_RP - (tau**0.7_RP))**3   ! (6) SB06
          PNCaut(k,i,j)  = coef_aut0*LC(k,i,j)*LC(k,i,j)*rho_fac*rho_fac    ! (9) SB06 sc+aut
          ! lc = lwc*(1-tau), lr = lwc*tau
          PLCaut(k,i,j)  = coef_aut1*lwc*lwc*xc(k,i,j)*xc(k,i,j) &          ! (4) SB06
               *((1.0_RP-tau)*(1.0_RP-tau) + psi_aut)*rho_fac*rho_fac        !
          PNRaut(k,i,j)  = -rx_sep*PLCaut(k,i,j)                           ! (A7) SB01
          !
          ! Accretion ( cloud-rain => rain )
          psi_acc       =(tau/(tau+thr_acc))**4                          ! (8) SB06
          PLCacc(k,i,j)  = -kcr*LC(k,i,j)*LR(k,i,j)*rho_fac*psi_acc         ! (7) SB06
          PNCacc(k,i,j)  = -kcr*NC(k,i,j)*LR(k,i,j)*rho_fac*psi_acc         ! (A6) SB01
          !
          ! Self-collection ( rain-rain => rain )
          PNRslc(k,i,j)  = -krr*NR(k,i,j)*LR(k,i,j)*rho_fac                 ! (A.8) SB01
          !
          ! Collisional breakup of rain
          ddr           = min(1.E-3_RP, dr_xave(k,i,j) - dr_eq )
          if      (dr_xave(k,i,j) < dr_min )then                          ! negligible
             psi_brk      = -1.0_RP
             PNRbrk(k,i,j) = 0.0_RP
          else if (dr_xave(k,i,j) <= dr_eq  )then
             psi_brk      = kbr*ddr + 1.0_RP                               ! (14) SB06 (+1 is necessary)
             PNRbrk(k,i,j) = - (psi_brk + 1.0_RP)*PNRslc(k,i,j)              ! (13) SB06
          else
             psi_brk      = 2.0_RP*exp(kapbr*ddr) - 1.0_RP                   ! (15) SB06
             PNRbrk(k,i,j) = - (psi_brk + 1.0_RP)*PNRslc(k,i,j)              ! (13) SB06
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
       IA, JA, KA,          & ! in
       IS, IE,           & ! in
       JS, JE,           & ! in
       KS, KE,           & ! in
       PLCdep,               & ! out
       PLRdep, PNRdep,       & ! out
       PLIdep, PNIdep,       & ! out
       PLSdep, PNSdep,       & ! out
       PLGdep, PNGdep,       & ! out
       PLImlt, PNImlt,       & ! out
       PLSmlt, PNsmlt,       & ! out
       PLGmlt, PNGmlt,       & ! out
       rho, tem, pre, LV,    & ! in
       qd,                   & ! in
       esw, esi,             & ! in
       NC, NR, NI, NS, NG,   & ! in
       LI, LS, LG,           & ! in
       xc, xr, xi, xs, xg,   & ! in
       vt_xave,              &
       dc_xave, dr_xave, di_xave, ds_xave, dg_xave ) ! in
    implicit none

    integer, intent(in)  :: IA,JA,KA
    integer, intent(in)  :: IS,JS,KS
    integer, intent(in)  :: IE,JE,KE
    ! Diffusion growth or Evaporation, Sublimation
    real(RP), intent(out) :: PLCdep(KA,IA,JA)  ! mass change   for cloud, [Add]  09/08/18 T.Mitsui
    real(RP), intent(out) :: PLRdep(KA,IA,JA)  ! mass change   for rain deposion
    real(RP), intent(out) :: PNRdep(KA,IA,JA)  ! number change
    real(RP), intent(out) :: PLIdep(KA,IA,JA)  ! mass          for cloud ice
    real(RP), intent(out) :: PNIdep(KA,IA,JA)  ! number
    real(RP), intent(out) :: PLSdep(KA,IA,JA)  ! mass          for snow
    real(RP), intent(out) :: PNSdep(KA,IA,JA)  ! number
    real(RP), intent(out) :: PLGdep(KA,IA,JA)  ! mass          for graupel
    real(RP), intent(out) :: PNGdep(KA,IA,JA)  ! number
    ! Melting under condition(T > 273.15K and Gm > 0.0 )
    real(RP), intent(out) :: PLImlt(KA,IA,JA)  ! mass          for cloud ice melting
    real(RP), intent(out) :: PNImlt(KA,IA,JA)  ! number
    real(RP), intent(out) :: PLSmlt(KA,IA,JA)  ! mass          for snow
    real(RP), intent(out) :: PNSmlt(KA,IA,JA)  ! number
    real(RP), intent(out) :: PLGmlt(KA,IA,JA)  ! mass          for graupel
    real(RP), intent(out) :: PNGmlt(KA,IA,JA)  ! number

    real(RP), intent(in)  :: rho(KA,IA,JA)     ! air density
    real(RP), intent(in)  :: tem(KA,IA,JA)     ! air temperature
    real(RP), intent(in)  :: pre(KA,IA,JA)     ! air pressure
    real(RP), intent(in)  :: qd(KA,IA,JA)      ! mixing ratio of dry air
    real(RP), intent(in)  :: esw(KA,IA,JA)     ! saturation vapor pressure(liquid water)
    real(RP), intent(in)  :: esi(KA,IA,JA)     ! saturation vapor pressure(solid water)
    real(RP), intent(in)  :: LV(KA,IA,JA)      ! mass   of vapor
    real(RP), intent(in)  :: NC(KA,IA,JA)      ! number of cloud  09/08/18 [Add] T.Mitsui
    real(RP), intent(in)  :: NR(KA,IA,JA)      ! number of rain
    real(RP), intent(in)  :: NI(KA,IA,JA)      !        of cloud ice
    real(RP), intent(in)  :: NS(KA,IA,JA)      !        of snow
    real(RP), intent(in)  :: NG(KA,IA,JA)      !        of graupel
    real(RP), intent(in)  :: LI(KA,IA,JA)      ! mass   of cloud ice
    real(RP), intent(in)  :: LS(KA,IA,JA)      ! mass   of snow
    real(RP), intent(in)  :: LG(KA,IA,JA)      ! mass   of graupel
    real(RP), intent(in)  :: xc(KA,IA,JA)      ! mean mass of cloud(filtered) [Add] 09/08/18 T.Mitsui
    real(RP), intent(in)  :: xr(KA,IA,JA)      ! mean mass of rain(filtered)
    real(RP), intent(in)  :: xi(KA,IA,JA)      !           of cloud ice(filtered)
    real(RP), intent(in)  :: xs(KA,IA,JA)      !           of snow(filtered)
    real(RP), intent(in)  :: xg(KA,IA,JA)      !           of graupel(filtered)
    ! Notice following values differ from mean terminal velocity or diameter.
    ! mean(vt(x)) /= vt(mean(x)) and mean(D(x)) /= D(mean(x))
    ! Following ones are vt(mean(x)) and D(mean(x)).
    real(RP), intent(in)  :: vt_xave(KA,IA,JA,HYDRO_MAX,2) ! terminal velocity of mean cloud 09/08/18 [Add], T.Mitsui
    !
    real(RP), intent(in)  :: dc_xave(KA,IA,JA) ! diameter of mean cloud 09/08/18 [Add] T.Mitsui
    real(RP), intent(in)  :: dr_xave(KA,IA,JA) ! diameter of mean rain
    real(RP), intent(in)  :: di_xave(KA,IA,JA) !                  ice
    real(RP), intent(in)  :: ds_xave(KA,IA,JA) !                  snow
    real(RP), intent(in)  :: dg_xave(KA,IA,JA) !                  graupel
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
    real(RP) ::                   ventLR(KA,IA,JA)     !
    real(RP) :: ventNI(KA,IA,JA), ventLI(KA,IA,JA)     !
    real(RP) :: ventNS(KA,IA,JA), ventLS(KA,IA,JA)     !
    real(RP) :: ventNG(KA,IA,JA), ventLG(KA,IA,JA)     !
    !
    real(RP), parameter :: Re_max=1.E+3_RP
    real(RP), parameter :: Re_min=1.E-4_RP
    !
    integer :: i, j, k
    !
    ! Notice,T.Mitsui
    ! Vapor deposition and melting would not be solved iteratively to reach equilibrium.
    ! Because following phenomena are not adjustment but transition.
    ! Just time-scales differ among them.
    ! If we would treat more appropreately, there would be time-splitting method to solve each ones.
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          temc    = tem(k,i,j) - T00   ! degC
          temc_lim= max(temc, -40._RP )       ! [Add] 09/08/18 T.Mitsui, Pruppacher and Klett(1997),(13-3)
          rho_lim = max(rho(k,i,j),rho_min)   ! [Add] 09/08/18 T.Mitsui
          qv      = LV(k,i,j)/rho_lim
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
          Gwr     = 4.0_RP*PI/cap(I_QR)/Gw
          Gii     = 4.0_RP*PI/cap(I_QI)/Gi
          Gis     = 4.0_RP*PI/cap(I_QS)/Gi
          Gig     = 4.0_RP*PI/cap(I_QG)/Gi
          ! vent: ventilation effect( asymmetry vapor field around particles due to aerodynamic )
          ! SB06 (30),(31) and each coefficient is by (88),(89)
          Nsc_r3  = (nua/Dw)**(0.33333333_RP)                    ! (Schmidt number )^(1/3)
          !
!          Nrecs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QC,1)*dc_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud
          Nrers_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QR,1)*dr_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) rain
          Nreis_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QI,1)*di_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
          Nress_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QS,1)*ds_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) snow
          Nregs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QG,1)*dg_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) graupel
          !
!          Nrecl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QC,2)*dc_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud
          Nrerl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QR,2)*dr_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) rain
          Nreil_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QI,2)*di_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud ice
          Nresl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QS,2)*ds_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) snow
          Nregl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QG,2)*dg_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) graupel
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
          ventLR_s = ah_vent1(I_QR,1) + bh_vent1(I_QR,1)*NscNrer_s
          ventLR_l = ah_vent1(I_QR,2) + bh_vent1(I_QR,2)*NscNrer_l
          !
          ventNI_s = ah_vent0(I_QI,1) + bh_vent0(I_QI,1)*NscNrei_s
          ventNI_l = ah_vent0(I_QI,2) + bh_vent0(I_QI,2)*NscNrei_l
          ventLI_s = ah_vent1(I_QI,1) + bh_vent1(I_QI,1)*NscNrei_s
          ventLI_l = ah_vent1(I_QI,2) + bh_vent1(I_QI,2)*NscNrei_l
          !
          ventNS_s = ah_vent0(I_QS,1) + bh_vent0(I_QS,1)*NscNres_s
          ventNS_l = ah_vent0(I_QS,2) + bh_vent0(I_QS,2)*NscNres_l
          ventLS_s = ah_vent1(I_QS,1) + bh_vent1(I_QS,1)*NscNres_s
          ventLS_l = ah_vent1(I_QS,2) + bh_vent1(I_QS,2)*NscNres_l
          !
          ventNG_s = ah_vent0(I_QG,1) + bh_vent0(I_QG,1)*NscNreg_s
          ventNG_l = ah_vent0(I_QG,2) + bh_vent0(I_QG,2)*NscNreg_l
          ventLG_s = ah_vent1(I_QG,1) + bh_vent1(I_QG,1)*NscNreg_s
          ventLG_l = ah_vent1(I_QG,2) + bh_vent1(I_QG,2)*NscNreg_l
          !
          ! branch is 1.4 for rain, snow, graupel; is 1.0 for ice (PK97, 13-60,-61,-88,-89).
          !
          wtr     = ( min(max( NscNrer_s*r_14, 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.4*0.5 and 1.4*2
          wti     = ( min(max( NscNrei_s     , 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.0*0.5 and 1.0*2
          wts     = ( min(max( NscNres_s*r_14, 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.4*0.5 and 1.4*2
          wtg     = ( min(max( NscNreg_s*r_14, 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.4*0.5 and 1.4*2
          ! interpolation between two branches
          ventNI(k,i,j)  = (1.0_RP-wti)*ventNI_s + wti*ventNI_l
          ventNS(k,i,j)  = (1.0_RP-wts)*ventNS_s + wts*ventNS_l
          ventNG(k,i,j)  = (1.0_RP-wtg)*ventNG_s + wtg*ventNG_l
          !
          ventLR(k,i,j)  = (1.0_RP-wtr)*ventLR_s + wtr*ventLR_l
          ventLI(k,i,j)  = (1.0_RP-wti)*ventLI_s + wti*ventLI_l
          ventLS(k,i,j)  = (1.0_RP-wts)*ventLS_s + wts*ventLS_l
          ventLG(k,i,j)  = (1.0_RP-wtg)*ventLG_s + wtg*ventLG_l
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
!!$       PLCdep(k,i,j) = Gwr*ssw*NC(k,i,j)*dc_xave(k,i,j)*coef_deplc
!!$       PLRdep(k,i,j) = Gwr*ssw*NR(k,i,j)*dr_xave(k,i,j)*ventLR
!!$       PLIdep(k,i,j) = Gii*ssi*NI(k,i,j)*di_xave(k,i,j)*ventLI
!!$       PLSdep(k,i,j) = Gis*ssi*NS(k,i,j)*ds_xave(k,i,j)*ventLS
!!$       PLGdep(k,i,j) = Gig*ssi*NG(k,i,j)*dg_xave(k,i,j)*ventLG
          PLCdep(k,i,j) = Gwr*NC(k,i,j)*dc_xave(k,i,j)*coef_deplc
          PLRdep(k,i,j) = Gwr*NR(k,i,j)*dr_xave(k,i,j)*ventLR(k,i,j)
          PLIdep(k,i,j) = Gii*NI(k,i,j)*di_xave(k,i,j)*ventLI(k,i,j)
          PLSdep(k,i,j) = Gis*NS(k,i,j)*ds_xave(k,i,j)*ventLS(k,i,j)
          PLGdep(k,i,j) = Gig*NG(k,i,j)*dg_xave(k,i,j)*ventLG(k,i,j)
          PNRdep(k,i,j) = PLRdep(k,i,j)/xr(k,i,j)
          PNIdep(k,i,j) = 0.0_RP
          PNSdep(k,i,j) = PLSdep(k,i,j)/xs(k,i,j)
          PNGdep(k,i,j) = PLGdep(k,i,j)/xg(k,i,j)
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
          if( (temc>=0.0_RP) .and. (Gm>0.0_RP) )then !  if Gm==0 then rh and tem is critical value for melting process.
             ! 08/05/16 [Mod] T.Mitsui, change term of PLimlt. N_i => L_i/ (limited x_i)
             ! because melting never occur when N_i=0.
             PLImlt(k,i,j) = - Gm * LI(k,i,j)*di_xave(k,i,j)*ventLI(k,i,j)/xi(k,i,j)
             ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
             PNImlt(k,i,j) = - Gm * NI(k,i,j)*di_xave(k,i,j)*ventNI(k,i,j)/xi(k,i,j) ! 09/08/18 [Mod] recover, T.Mitsui
             PLSmlt(k,i,j) = - Gm * LS(k,i,j)*ds_xave(k,i,j)*ventLS(k,i,j)/xs(k,i,j)
             ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
             PNSmlt(k,i,j) = - Gm * NS(k,i,j)*ds_xave(k,i,j)*ventNS(k,i,j)/xs(k,i,j) ! 09/08/18 [Mod] recover, T.Mitsui
             PLGmlt(k,i,j) = - Gm * LG(k,i,j)*dg_xave(k,i,j)*ventLG(k,i,j)/xg(k,i,j)
             ! [Mod] 08/08/23 T.Mitsui for Seifert(2008)
             PNGmlt(k,i,j) = - Gm * NG(k,i,j)*dg_xave(k,i,j)*ventNG(k,i,j)/xg(k,i,j) ! 09/08/18 [Mod] recover, T.Mitsui
          else
             PLImlt(k,i,j) = 0.0_RP
             PNImlt(k,i,j) = 0.0_RP
             PLSmlt(k,i,j) = 0.0_RP
             PNSmlt(k,i,j) = 0.0_RP
             PLGmlt(k,i,j) = 0.0_RP
             PNGmlt(k,i,j) = 0.0_RP
          end if
          !
       end do
    end do
    end do
    !
    return
  end subroutine dep_vapor_melt_ice_kij
  !-----------------------------------------------------------------------------
  subroutine freezing_water_kij( &
       IA, JA, KA,          &
       IS, IE,           &
       JS, JE,           &
       KS, KE,           &
       dt,                   &
       PLChom, PNChom,       &
       PLChet, PNChet,       &
       PLRhet, PNRhet,       &
       LC, LR, NC, NR, xc, xr, tem   )
    !
    ! In this subroutine,
    ! We assumed surface temperature of droplets are same as environment.
    implicit none

    integer, intent(in) :: IA, JA, KA
    integer, intent(in) :: IS, JS, KS
    integer, intent(in) :: IE, JE, KE
    !
    real(RP), intent(in) :: dt
    ! freezing in nature (homogenous)
    real(RP), intent(out):: PLChom(KA,IA,JA) ! cloud water => cloud ice
    real(RP), intent(out):: PNChom(KA,IA,JA) ! cloud water => cloud ice
    ! freezing via aerosols (heterogenous)
    real(RP), intent(out):: PLChet(KA,IA,JA) ! cloud water => cloud ice
    real(RP), intent(out):: PLRhet(KA,IA,JA) ! rain        => graupel
    real(RP), intent(out):: PNChet(KA,IA,JA) ! cloud water => cloud ice
    real(RP), intent(out):: PNRhet(KA,IA,JA) ! rain        => graupel
    !
    real(RP), intent(in) :: tem(KA,IA,JA)
    !
    real(RP), intent(in) :: LC(KA,IA,JA)
    real(RP), intent(in) :: LR(KA,IA,JA)
    real(RP), intent(in) :: NC(KA,IA,JA)
    real(RP), intent(in) :: NR(KA,IA,JA)
    real(RP), intent(in) :: xc(KA,IA,JA)
    real(RP), intent(in) :: xr(KA,IA,JA)
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
    !
    integer :: i,j,k
    !
    rdt = 1.0_RP/dt
    !
    coef_m2_c =   coef_m2(I_QC)
    coef_m2_r =   coef_m2(I_QR)
    !
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
          PLChom(k,i,j) = 0.0_RP
          PNChom(k,i,j) = 0.0_RP
          ! Heterogenous Freezing
          PLChet(k,i,j) = -rdt*LC(k,i,j)*( 1.0_RP - exp( -coef_m2_c*xc(k,i,j)*(Jhet+Jhom)*dt ) )
          PNChet(k,i,j) = -rdt*NC(k,i,j)*( 1.0_RP - exp( -          xc(k,i,j)*(Jhet+Jhom)*dt ) )
          PLRhet(k,i,j) = -rdt*LR(k,i,j)*( 1.0_RP - exp( -coef_m2_r*xr(k,i,j)*(Jhet+Jhom)*dt ) )
          PNRhet(k,i,j) = -rdt*NR(k,i,j)*( 1.0_RP - exp( -          xr(k,i,j)*(Jhet+Jhom)*dt ) )
       end do
    end do
    end do
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
    implicit none

    real(RP), intent(out) :: velw(KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(RP), intent(in)  :: rhoq(KA,IA,JA,QA) ! rho * q
    real(RP), intent(in)  :: DENS(KA,IA,JA)    ! rho
    real(RP), intent(in)  :: temp(KA,IA,JA)    ! temperature
    real(RP), intent(in)  :: pres(KA,IA,JA)    ! pressure

    real(RP) :: xq      (KA,6)  ! average mass of 1 particle( mass/number )

    real(RP) :: rhofac  (KA)    ! density factor for terminal velocity( air friction )
    real(RP) :: rhofac_q(KA,6)

    real(RP) :: rlambdar(KA)    ! work for diagnosis of Rain DSD ( Seifert, 2008 )
    real(RP) :: mud_r           !
    real(RP) :: dq      (KA,QA) ! work for Rogers etal. (1993)
    real(RP) :: weight  (KA,QA) !
    real(RP) :: velq_s  (KA,QA) ! terminal velocity for small branch of Rogers formula
    real(RP) :: velq_l  (KA,QA) ! terminal velocity for large branch of Rogers formula

    integer :: k, i, j
    !---------------------------------------------------------------------------

    mud_r = 3.0_RP * nu(I_QR) + 2.0_RP

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          rhofac(k) = rho_0 / max( DENS(k,i,j), rho_min )
       enddo

       ! average mass
       do k = KS, KE
          xq(k,I_QC) = rhoq(k,i,j,I_QC) / ( rhoq(k,i,j,I_NC) + nqmin(I_QC) )
          xq(k,I_QR) = rhoq(k,i,j,I_QR) / ( rhoq(k,i,j,I_NR) + nqmin(I_QR) )
          xq(k,I_QI) = rhoq(k,i,j,I_QI) / ( rhoq(k,i,j,I_NI) + nqmin(I_QI) )
          xq(k,I_QS) = rhoq(k,i,j,I_QS) / ( rhoq(k,i,j,I_NS) + nqmin(I_QS) )
          xq(k,I_QG) = rhoq(k,i,j,I_QG) / ( rhoq(k,i,j,I_NG) + nqmin(I_QG) )
       enddo

       ! limiter
       do k  = KS, KE
          xq(k,I_QC) = max( min( xq(k,I_QC), xqmax(I_QC) ), xqmin(I_QC) )
          xq(k,I_QR) = max( min( xq(k,I_QR), xqmax(I_QR) ), xqmin(I_QR) )
          xq(k,I_QI) = max( min( xq(k,I_QI), xqmax(I_QI) ), xqmin(I_QI) )
          xq(k,I_QS) = max( min( xq(k,I_QS), xqmax(I_QS) ), xqmin(I_QS) )
          xq(k,I_QG) = max( min( xq(k,I_QG), xqmax(I_QG) ), xqmin(I_QG) )
       enddo

       do k = KS, KE
          rlambdar(k) = a_m(I_QR) * xq(k,I_QR)**b_m(I_QR) &
                      * ( (mud_r+3.0_RP) * (mud_r+2.0_RP) * (mud_r+1.0_RP) )**(-0.333333333_RP)
       enddo

       ! Improved Rogers formula by T.Mitsui
       ! weigthed diameter
       do k = KS, KE
          dq(k,I_QR) = ( 4.0_RP + mud_r ) * rlambdar(k) ! D^(3)+mu weighted mean diameter
          dq(k,I_NR) = ( 1.0_RP + mud_r ) * rlambdar(k)
       enddo
       do k = KS, KE
          dq(k,I_QI) = coef_dave_L(I_QI) * a_m(I_QI) * xq(k,I_QI)**b_m(I_QI)
          dq(k,I_NI) = coef_dave_N(I_QI) * a_m(I_QI) * xq(k,I_QI)**b_m(I_QI)
       enddo
       do k = KS, KE
          dq(k,I_QS) = coef_dave_L(I_QS) * a_m(I_QS) * xq(k,I_QS)**b_m(I_QS)
          dq(k,I_NS) = coef_dave_N(I_QS) * a_m(I_QS) * xq(k,I_QS)**b_m(I_QS)
       enddo
       do k = KS, KE
          dq(k,I_QG) = coef_dave_L(I_QG) * a_m(I_QG) * xq(k,I_QG)**b_m(I_QG)
          dq(k,I_NG) = coef_dave_N(I_QG) * a_m(I_QG) * xq(k,I_QG)**b_m(I_QG)
       enddo

       ! small/large branch terminal velocity
       do k = KS, KE
          velq_s(k,I_QR) = coef_vtr_ar2 * dq(k,I_QR) * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar(k) )**(-5-mud_r) )
          velq_l(k,I_QR) = coef_vtr_ar1 - coef_vtr_br1 * ( 1.0_RP + coef_vtr_cr1*rlambdar(k) )**(-4-mud_r)
          velq_s(k,I_NR) = coef_vtr_ar2 * dq(k,I_QR) * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar(k) )**(-2-mud_r) )
          velq_l(k,I_NR) = coef_vtr_ar1 - coef_vtr_br1 * ( 1.0_RP + coef_vtr_cr1*rlambdar(k) )**(-1-mud_r)
       enddo
       do k = KS, KE
          velq_s(k,I_QI) = coef_vt1(I_QI,1) * xq(k,I_QI)**beta_v (I_QI,1)
          velq_l(k,I_QI) = coef_vt1(I_QI,2) * xq(k,I_QI)**beta_v (I_QI,2)
          velq_s(k,I_NI) = coef_vt0(I_QI,1) * xq(k,I_QI)**beta_vn(I_QI,1)
          velq_l(k,I_NI) = coef_vt0(I_QI,2) * xq(k,I_QI)**beta_vn(I_QI,2)
       enddo
       do k = KS, KE
          velq_s(k,I_QS) = coef_vt1(I_QS,1) * xq(k,I_QS)**beta_v (I_QS,1)
          velq_l(k,I_QS) = coef_vt1(I_QS,2) * xq(k,I_QS)**beta_v (I_QS,2)
          velq_s(k,I_NS) = coef_vt0(I_QS,1) * xq(k,I_QS)**beta_vn(I_QS,1)
          velq_l(k,I_NS) = coef_vt0(I_QS,2) * xq(k,I_QS)**beta_vn(I_QS,2)
       enddo
       do k = KS, KE
          velq_s(k,I_QG) = coef_vt1(I_QG,1) * xq(k,I_QG)**beta_v (I_QG,1)
          velq_l(k,I_QG) = coef_vt1(I_QG,2) * xq(k,I_QG)**beta_v (I_QG,2)
          velq_s(k,I_NG) = coef_vt0(I_QG,1) * xq(k,I_QG)**beta_vn(I_QG,1)
          velq_l(k,I_NG) = coef_vt0(I_QG,2) * xq(k,I_QG)**beta_vn(I_QG,2)
       enddo

       ! weighting coefficient for 2-branches is determined by ratio between 0.745mm and weighted diameter
       ! SB06 Table.1
       do k = KS, KE
          weight(k,I_QR) = 0.5_RP * ( 1.0_RP + tanh( PI * log( dq(k,I_QR)/d_vtr_branch ) ) ) ! Lr
          weight(k,I_QI) = 0.5_RP * ( 1.0_RP + log( dq(k,I_QI)/d0_li ) )
          weight(k,I_QS) = 0.5_RP * ( 1.0_RP + log( dq(k,I_QS)/d0_ls ) )
          weight(k,I_QG) = 0.5_RP * ( 1.0_RP + log( dq(k,I_QG)/d0_lg ) )
          weight(k,I_NR) = 0.5_RP * ( 1.0_RP + tanh( PI * log( dq(k,I_NR)/d_vtr_branch ) ) ) ! Nr
          weight(k,I_NI) = 0.5_RP * ( 1.0_RP + log( dq(k,I_NI)/d0_ni ) )
          weight(k,I_NS) = 0.5_RP * ( 1.0_RP + log( dq(k,I_NS)/d0_ns ) )
          weight(k,I_NG) = 0.5_RP * ( 1.0_RP + log( dq(k,I_NG)/d0_ng ) )
       enddo

          ! filter is used to avoid unrealistic terminal velocity( when ni,ns,ng are too small )
       do k = KS, KE
          weight(k,I_QI) = min( max( weight(k,I_QI), 0.0_RP ), 1.0_RP )
          weight(k,I_QS) = min( max( weight(k,I_QS), 0.0_RP ), 1.0_RP )
          weight(k,I_QG) = min( max( weight(k,I_QG), 0.0_RP ), 1.0_RP )
          weight(k,I_NI) = min( max( weight(k,I_NI), 0.0_RP ), 1.0_RP )
          weight(k,I_NS) = min( max( weight(k,I_NS), 0.0_RP ), 1.0_RP )
          weight(k,I_NG) = min( max( weight(k,I_NG), 0.0_RP ), 1.0_RP )
       enddo

       do k = KS, KE
          rhofac_q(k,I_QC) = rhofac(k) ** gamma_v(I_QC)
          rhofac_q(k,I_QR) = rhofac(k) ** gamma_v(I_QR)
          rhofac_q(k,I_QI) = ( pres(k,i,j)/pre0_vt )**a_pre0_vt * ( temp(k,i,j)/tem0_vt )**a_tem0_vt
          rhofac_q(k,I_QS) = rhofac_q(k,I_QI)
          rhofac_q(k,I_QG) = rhofac_q(k,I_QI)
       enddo

       ! interpolated terminal velocity
       ! SB06(78) these are defined as negative value
       velw(:,i,j,I_QV) = CONST_UNDEF

       do k = KS, KE
          velw(k,i,j,I_QC) = -rhofac_q(k,I_QC) * coef_vt1(I_QC,1) * xq(k,I_QC)**beta_v(I_QC,1)
          velw(k,i,j,I_NC) = -rhofac_q(k,I_QC) * coef_vt0(I_QC,1) * xq(k,I_QC)**beta_vn(I_QC,1)
       enddo
       velw(KS-1,i,j,I_QC) = velw(KS,i,j,I_QC)
       velw(KS-1,i,j,I_NC) = velw(KS,i,j,I_NC)

       do k = KS, KE
          velw(k,i,j,I_QR) = -rhofac_q(k,I_QR) * ( velq_l(k,I_QR) * (          weight(k,I_QR) ) &
                                                 + velq_s(k,I_QR) * ( 1.0_RP - weight(k,I_QR) ) )
          velw(k,i,j,I_NR) = -rhofac_q(k,I_QR) * ( velq_l(k,I_NR) * (          weight(k,I_NR) ) &
                                                 + velq_s(k,I_NR) * ( 1.0_RP - weight(k,I_NR) ) )
       enddo
       velw(KS-1,i,j,I_QR) = velw(KS,i,j,I_QR)
       velw(KS-1,i,j,I_NR) = velw(KS,i,j,I_NR)
       do k = KS, KE
          velw(k,i,j,I_QI) = -rhofac_q(k,I_QI) * ( velq_l(k,I_QI) * (          weight(k,I_QI) ) &
                                                 + velq_s(k,I_QI) * ( 1.0_RP - weight(k,I_QI) ) )
          velw(k,i,j,I_NI) = -rhofac_q(k,I_QI) * ( velq_l(k,I_NI) * (          weight(k,I_NI) ) &
                                                 + velq_s(k,I_NI) * ( 1.0_RP - weight(k,I_NI) ) )
       enddo
       velw(KS-1,i,j,I_QI) = velw(KS,i,j,I_QI)
       velw(KS-1,i,j,I_NI) = velw(KS,i,j,I_NI)
       do k = KS, KE
          velw(k,i,j,I_QS) = -rhofac_q(k,I_QS) * ( velq_l(k,I_QS) * (          weight(k,I_QS) ) &
                                                 + velq_s(k,I_QS) * ( 1.0_RP - weight(k,I_QS) ) )
          velw(k,i,j,I_NS) = -rhofac_q(k,I_QS) * ( velq_l(k,I_NS) * (          weight(k,I_NS) ) &
                                                 + velq_s(k,I_NS) * ( 1.0_RP - weight(k,I_NS) ) )
       enddo
       velw(KS-1,i,j,I_QS) = velw(KS,i,j,I_QS)
       velw(KS-1,i,j,I_NS) = velw(KS,i,j,I_NS)
       do k = KS, KE
          velw(k,i,j,I_QG) = -rhofac_q(k,I_QG) * ( velq_l(k,I_QG) * (          weight(k,I_QG) ) &
                                                 + velq_s(k,I_QG) * ( 1.0_RP - weight(k,I_QG) ) )
          velw(k,i,j,I_NG) = -rhofac_q(k,I_QG) * ( velq_l(k,I_NG) * (          weight(k,I_NG) ) &
                                                 + velq_s(k,I_NG) * ( 1.0_RP - weight(k,I_NG) ) )
       enddo
       velw(KS-1,i,j,I_QG) = velw(KS,i,j,I_QG)
       velw(KS-1,i,j,I_NG) = velw(KS,i,j,I_NG)

    enddo
    enddo

    velw(   1:KS-2,:,:,:) = CONST_UNDEF
    velw(KE+1:KA  ,:,:,:) = CONST_UNDEF

    return
  end subroutine MP_terminal_velocity
  !----------------------------------------------------------------
  subroutine update_by_phase_change_kij(   &
       ntdiv    , ntmax,         & ! in [Add] 10/08/03
       IA, JA, KA,               &
       IS, IE,                   & ! in
       JS, JE,                   & ! in
       KS, KE,                   & ! in
       QA,                       & ! in
       dt,                       & ! in
       gsgam2,                   & ! in
       z,                        & ! in
       dz,                       & ! in
       wh,                       & ! in
       dTdt_rad,                 & ! in
       rhog,                     & ! in
       rhoge,                    & ! inout
       rhogq, q,                 & ! inout
       tem, pre,                 & ! inout
       cva,                      & ! out
       esw, esi, LV,             & ! in
       LC, LR, LI, LS, LG,       & ! in
       NC, NR, NI, NS, NG,       & ! in
       !+++ tendency terms
       ! homogeneous freezing
       PLChom, PNChom,           & ! in
       ! heterogeneous freezing
       PLChet, PNChet, PLRhet, PNRhet, &
       ! condensation/evaporation, deposition/sublimation
       PLCdep,         &
       PLRdep, PNRdep, &
       PLIdep, PNIdep, &
       PLSdep, PNSdep, &
       PLGdep, PNGdep, &
       ! melting term
       PLImlt, PNImlt, &
       PLSmlt, PNSmlt, &
       PLGmlt, PNGmlt, &
       !
       flag_history_in,&
       ! column integrated tendency terms
       sl_PLCdep,      &      !
       sl_PLRdep, sl_PNRdep ) !
    use scale_atmos_thermodyn, only: &
       AQ_CV, &
       AQ_CP
    use scale_atmos_saturation, only: &
       moist_pres2qsat_liq  => ATMOS_SATURATION_pres2qsat_liq,  &
       moist_pres2qsat_ice  => ATMOS_SATURATION_pres2qsat_ice,  &
       moist_dqsw_dtem_rho  => ATMOS_SATURATION_dqsw_dtem_rho,  &
       moist_dqsi_dtem_rho  => ATMOS_SATURATION_dqsi_dtem_rho,  &
       moist_dqsw_dtem_dpre => ATMOS_SATURATION_dqsw_dtem_dpre, &
       moist_dqsi_dtem_dpre => ATMOS_SATURATION_dqsi_dtem_dpre
    implicit none

    integer, intent(in)    :: ntdiv                   ! [Add] 10/08/03
    integer, intent(in)    :: ntmax                   ! [Add] 10/08/03
    !
    integer, intent(in)    :: IA, JA, KA                   !
    integer, intent(in)    :: IS, IE        !
    integer, intent(in)    :: JS, JE        !
    integer, intent(in)    :: KS, KE        !
    integer, intent(in)    :: QA                   ! tracer number
    real(RP), intent(in)    :: dt                      ! time step[s]
    real(RP), intent(in)    :: gsgam2(KA,IA,JA)      ! metric
    real(RP), intent(in)    :: z(KA)           ! altitude [m]
    real(RP), intent(in)    :: dz(KA)          ! altitude [m]
    real(RP), intent(in)    :: wh(KA,IA,JA)          ! vertical velocity @ half point[m/s]
    real(RP), intent(in)    :: dTdt_rad(KA,IA,JA)    ! temperture tendency by radiation[K/s]
    real(RP), intent(in)    :: rhog(KA,IA,JA)        ! density[kg/m3]
    real(RP), intent(inout) :: rhoge(KA,IA,JA)       ! internal energy[J/m3]
    real(RP), intent(inout) :: rhogq(KA,IA,JA,QA) ! tracers[kg/m3]
    real(RP), intent(inout) :: q(KA,IA,JA,QA)     ! tracers mixing ratio[kg/kg]
    real(RP), intent(inout) :: tem(KA,IA,JA)         ! temperature[K]
    real(RP), intent(inout) :: pre(KA,IA,JA)         ! pressure[Pa]
    real(RP), intent(out)   :: cva(KA,IA,JA)         ! specific heat at constant volume
    real(RP), intent(in)    :: esw(KA,IA,JA)         ! saturated vapor pressure for liquid
    real(RP), intent(in)    :: esi(KA,IA,JA)         !                          for ice
    real(RP), intent(in)    :: lv(KA,IA,JA)          ! vapor mass [kg/m3]
    real(RP), intent(in)    :: lc(KA,IA,JA), nc(KA,IA,JA) ! cloud mass [kg/m3], number[/m3]
    real(RP), intent(in)    :: lr(KA,IA,JA), nr(KA,IA,JA) ! rain
    real(RP), intent(in)    :: li(KA,IA,JA), ni(KA,IA,JA) ! ice
    real(RP), intent(in)    :: ls(KA,IA,JA), ns(KA,IA,JA) ! snow
    real(RP), intent(in)    :: lg(KA,IA,JA), ng(KA,IA,JA) ! graupel
    !+++ Freezing tendency[kg/m3/s]
    real(RP), intent(inout) :: PLChom(KA,IA,JA), PNChom(KA,IA,JA)
    real(RP), intent(inout) :: PLChet(KA,IA,JA), PNChet(KA,IA,JA)
    real(RP), intent(inout) :: PLRhet(KA,IA,JA), PNRhet(KA,IA,JA)
    !+++ Condensation/Evaporation, Deposition/Sublimaion tendency[kg/m3/s]
    real(RP), intent(inout) :: PLCdep(KA,IA,JA)
    real(RP), intent(inout) :: PLRdep(KA,IA,JA), PNRdep(KA,IA,JA)
    real(RP), intent(inout) :: PLIdep(KA,IA,JA), PNIdep(KA,IA,JA)
    real(RP), intent(inout) :: PLSdep(KA,IA,JA), PNSdep(KA,IA,JA)
    real(RP), intent(inout) :: PLGdep(KA,IA,JA), PNGdep(KA,IA,JA)
    !+++ Melting Tendency[kg/m3/s]
    real(RP), intent(in)    :: PLImlt(KA,IA,JA), PNImlt(KA,IA,JA)
    real(RP), intent(in)    :: PLSmlt(KA,IA,JA), PNSmlt(KA,IA,JA)
    real(RP), intent(in)    :: PLGmlt(KA,IA,JA), PNGmlt(KA,IA,JA)
    !+++
    logical, intent(in)    :: flag_history_in
    !+++ Column integrated tendency[kg/m2/s]
    real(RP), intent(inout) :: sl_PLCdep(1,IA,JA)
    real(RP), intent(inout) :: sl_PLRdep(1,IA,JA), sl_PNRdep(1,IA,JA)
    !
    real(RP) :: Rmoist
    !
    real(RP) :: xi(KA,IA,JA)                     ! mean mass of ice particles
    real(RP) :: rrhog(KA,IA,JA)                  ! 1/rhog
    real(RP) :: wtem(KA,IA,JA)                   ! temperature[K]
    real(RP) :: qd(KA,IA,JA)                     ! mixing ratio of dry air
    !
    real(RP) :: r_cva                    ! specific heat at constant volume
    real(RP) :: cpa(KA,IA,JA), r_cpa   ! specific heat at constant pressure
    real(RP) :: qsw(KA,IA,JA)          ! saturated mixing ratio for liquid
    real(RP) :: qsi(KA,IA,JA)          ! saturated mixing ratio for solid
    real(RP) :: dqswdtem_rho(KA,IA,JA) ! (dqsw/dtem)_rho
    real(RP) :: dqsidtem_rho(KA,IA,JA) ! (dqsi/dtem)_rho
    real(RP) :: dqswdtem_pre(KA,IA,JA) ! (dqsw/dtem)_pre
    real(RP) :: dqsidtem_pre(KA,IA,JA) ! (dqsi/dtem)_pre
    real(RP) :: dqswdpre_tem(KA,IA,JA) ! (dqsw/dpre)_tem
    real(RP) :: dqsidpre_tem(KA,IA,JA) ! (dqsi/dpre)_tem
    !
    real(RP) :: w(KA,IA,JA)                     ! vetical velocity[m/s]
    real(RP) :: Acnd                              ! Pdynliq + Bergeron-Findeisen
    real(RP) :: Adep                              ! Pdyndep + Bergeron-Findeisen
    real(RP) :: aliqliq, asolliq
    real(RP) :: aliqsol, asolsol
    real(RP) :: Pdynliq                           ! production term of ssw by vertical motion
    real(RP) :: Pdynsol                           ! production term of ssi by vertical motion
    real(RP) :: Pradliq                           ! production term of ssw by radiation
    real(RP) :: Pradsol                           ! production term of ssi by radiation
    real(RP) :: taucnd(KA,IA,JA),   r_taucnd    ! time scale of ssw change by MP
    real(RP) :: taudep(KA,IA,JA),   r_taudep    ! time scale of ssi change by MP
    real(RP) :: taucnd_c(KA,IA,JA), r_taucnd_c  ! by cloud
    real(RP) :: taucnd_r(KA,IA,JA), r_taucnd_r  ! by rain
    real(RP) :: taudep_i(KA,IA,JA), r_taudep_i  ! by ice
    real(RP) :: taudep_s(KA,IA,JA), r_taudep_s  ! by snow
    real(RP) :: taudep_g(KA,IA,JA), r_taudep_g  ! by graupel
    ! alternative tendency through changing ssw and ssi
    real(RP) :: PNCdep(KA,IA,JA) ! [Add] 11/08/30 T.Mitsui
    real(RP) :: PLR2NR, PLI2NI, PLS2NS, PLG2NG
    real(RP) :: coef_a_cnd, coef_b_cnd
    real(RP) :: coef_a_dep, coef_b_dep
    !
    real(RP) :: frz_dqc, frz_dnc
    real(RP) :: frz_dqr, frz_dnr
    real(RP) :: mlt_dqi, mlt_dni
    real(RP) :: mlt_dqs, mlt_dns
    real(RP) :: mlt_dqg, mlt_dng
    real(RP) :: dep_dqi, dep_dni
    real(RP) :: dep_dqs, dep_dns
    real(RP) :: dep_dqg, dep_dng
    real(RP) :: dep_dqr, dep_dnr
    real(RP) :: dep_dqc, dep_dnc   ! 11/08/30 [Add] T.Mitsui, dep_dnc
    real(RP) :: r_xc_ccn, r_xi_ccn ! 11/08/30 [Add] T.Mitsui
    !
    real(RP) :: drhogqv
    real(RP) :: drhogqc, drhogqr, drhogqi, drhogqs, drhogqg
    real(RP) :: drhognc, drhognr, drhogni, drhogns, drhogng
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
    !

    ! [Add] 11/08/30 T.Mitsui
    if( flag_first )then
       flag_first = .false.
       rewind(IO_FID_CONF)
       read  (IO_FID_CONF,nml=nm_mp_sn14_condensation, end=100)
100    if( IO_L ) write (IO_FID_LOG,nml=nm_mp_sn14_condensation)
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
       fac_cndc_wrk = fac_cndc**(1.0_RP-b_m(I_QC))
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          PLCdep(k,i,j)  = PLCdep(k,i,j)*fac_cndc_wrk
       end do
       end do
       end do
       if( IO_L ) write(IO_FID_LOG,*) "taucnd:fac_cndc_wrk=",fac_cndc_wrk
    end if
    !
    do j = JS, JE
    do i = IS, IE
       do k = KS,KE
          rrhog(k,i,j) = 1.0_RP/rhog(k,i,j)
          ! Temperature lower limit is only used for saturation condition.
          ! On the other hand original "tem" is used for calculation of latent heat or energy equation.
          wtem(k,i,j)  = max( tem(k,i,j), tem_min )
          !
          if( z(k) <= 25000.0_RP )then
             w(k,i,j) = 0.5_RP*(wh(k,i,j) + wh(k+1,i,j))
          else
             w(k,i,j) = 0.0_RP
          end if
          !
          CALC_QDRY( qd(k,i,j), q, k, i, j, iqw )
          CALC_CV( cva(k,i,j), qd(k,i,j), q, k, i, j, iqw, CVdry, AQ_CV )
          CALC_CP( cpa(k,i,j), qd(k,i,j), q, k, i, j, iqw, CPdry, AQ_CP )
       end do
    end do
    end do
    !
    call moist_pres2qsat_liq ( qsw, wtem, pre )
    call moist_pres2qsat_ice ( qsi, wtem, pre )
    call moist_dqsw_dtem_rho ( dqswdtem_rho, wtem, rhog )
    call moist_dqsi_dtem_rho ( dqsidtem_rho, wtem, rhog )
    call moist_dqsw_dtem_dpre( dqswdtem_pre, dqswdpre_tem, wtem, pre )
    call moist_dqsi_dtem_dpre( dqsidtem_pre, dqsidpre_tem, wtem, pre )
    !
    do j=JS, JE
    do i=IS, IE
       do k=KS, KE
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
       end do
       !
       do k = KS, KE
          r_rvaptem        = 1.0_RP/(Rvap*wtem(k,i,j))
          lvsw             = esw(k,i,j)*r_rvaptem        ! rho=p/(Rv*T)
          lvsi             = esi(k,i,j)*r_rvaptem        !
          pv               = lv(k,i,j)*Rvap*tem(k,i,j)
          r_esw            = 1.0_RP/esw(k,i,j)
          r_esi            = 1.0_RP/esi(k,i,j)
          ssw              = min( MP_ssw_lim, ( pv*r_esw-1.0_RP ) )
          ssi              = pv*r_esi - 1.0_RP
          r_lvsw           = 1.0_RP/lvsw
          r_lvsi           = 1.0_RP/lvsi
          r_taucnd_c       = PLCdep(k,i,j)*r_lvsw
          r_taucnd_r       = PLRdep(k,i,j)*r_lvsw
          r_taudep_i       = PLIdep(k,i,j)*r_lvsi
          r_taudep_s       = PLSdep(k,i,j)*r_lvsi
          r_taudep_g       = PLGdep(k,i,j)*r_lvsi
!          taucnd_c(k,i,j)   = 1.0_RP/(r_taucnd_c+r_tau100day)
!          taucnd_r(k,i,j)   = 1.0_RP/(r_taucnd_r+r_tau100day)
!          taudep_i(k,i,j)   = 1.0_RP/(r_taudep_i+r_tau100day)
!          taudep_s(k,i,j)   = 1.0_RP/(r_taudep_s+r_tau100day)
!          taudep_g(k,i,j)   = 1.0_RP/(r_taudep_g+r_tau100day)
          !
          r_cva            = 1.0_RP/cva(k,i,j)
          r_cpa            = 1.0_RP/cpa(k,i,j)
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
          Pdynliq          = w(k,i,j)*GRAV * ( r_cpa*dqswdtem_pre(k,i,j) + rhog(k,i,j)*dqswdpre_tem(k,i,j) )
          Pdynsol          = w(k,i,j)*GRAV * ( r_cpa*dqsidtem_pre(k,i,j) + rhog(k,i,j)*dqsidpre_tem(k,i,j) )
          Pradliq          = -dTdt_rad(k,i,j)    * dqswdtem_rho(k,i,j)
          Pradsol          = -dTdt_rad(k,i,j)    * dqsidtem_rho(k,i,j)
          !
          ssw_o            = ssw
          ssi_o            = ssi
!             ssw_o            = ssw - Pdynliq*(dt_dyn-dt_mp)/qsw(k,i,j) + Pradliq*r_qsw*dt_mp
!             ssi_o            = ssi - Pdynsol*(dt_dyn-dt_mp)/qsi(k,i,j) + Pradsol*r_qsi*dt_mp
          !
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
          !
          uplim_cnd        = max( rhog(k,i,j)*ssw_o*qsw(k,i,j)*r_dt, 0.0_RP )
          lowlim_cnd       = min( rhog(k,i,j)*ssw_o*qsw(k,i,j)*r_dt, 0.0_RP )
          if( r_taucnd < r_tau100day )then
             taucnd(k,i,j)     = tau100day
             PLCdep(k,i,j) = max(lowlim_cnd, min(uplim_cnd, PLCdep(k,i,j)*ssw_o ))
             PLRdep(k,i,j) = max(lowlim_cnd, min(uplim_cnd, PLRdep(k,i,j)*ssw_o ))
             PNRdep(k,i,j) = min(0.0_RP, PNRdep(k,i,j)*ssw_o )
             PLR2NR           = 0.0_RP
          else
             taucnd(k,i,j)     = 1.0_RP/r_taucnd
             ! Production term for liquid water content
             coef_a_cnd = rhog(k,i,j)*Acnd*taucnd(k,i,j)
             coef_b_cnd = rhog(k,i,j)*taucnd(k,i,j)*r_dt*(ssw_o*qsw(k,i,j)-Acnd*taucnd(k,i,j)) * ( exp(-dt*r_taucnd) - 1.0_RP )
             PLCdep(k,i,j) = coef_a_cnd*r_taucnd_c - coef_b_cnd*r_taucnd_c
             PLR2NR           = PNRdep(k,i,j)/(PLRdep(k,i,j)+1.E-30_RP)
             PLRdep(k,i,j) = coef_a_cnd*r_taucnd_r - coef_b_cnd*r_taucnd_r
             PNRdep(k,i,j) = min(0.0_RP, PLRdep(k,i,j)*PLR2NR )
          end if
          !
          uplim_dep        = max( rhog(k,i,j)*ssi_o*qsi(k,i,j)*r_dt, 0.0_RP )
          lowlim_dep       = min( rhog(k,i,j)*ssi_o*qsi(k,i,j)*r_dt, 0.0_RP )
          if( r_taudep < r_tau100day )then
             taudep(k,i,j)     = tau100day
             PLIdep(k,i,j) = max(lowlim_dep, min(uplim_dep, PLIdep(k,i,j)*ssi_o ))
             PLSdep(k,i,j) = max(lowlim_dep, min(uplim_dep, PLSdep(k,i,j)*ssi_o ))
             PLGdep(k,i,j) = max(lowlim_dep, min(uplim_dep, PLGdep(k,i,j)*ssi_o ))
             PNIdep(k,i,j) = min(0.0_RP, PNIdep(k,i,j)*ssi_o )
             PNSdep(k,i,j) = min(0.0_RP, PNSdep(k,i,j)*ssi_o )
             PNGdep(k,i,j) = min(0.0_RP, PNGdep(k,i,j)*ssi_o )
          else
             taudep(k,i,j)     = 1.0_RP/r_taudep
             ! Production term for ice water content
             coef_a_dep = rhog(k,i,j)*Adep*taudep(k,i,j)
             coef_b_dep = rhog(k,i,j)*taudep(k,i,j)*r_dt*(ssi_o*qsi(k,i,j)-Adep*taudep(k,i,j)) * ( exp(-dt*r_taudep) - 1.0_RP )
             PLI2NI           = PNIdep(k,i,j)/max(PLIdep(k,i,j),1.E-30_RP)
             PLS2NS           = PNSdep(k,i,j)/max(PLSdep(k,i,j),1.E-30_RP)
             PLG2NG           = PNGdep(k,i,j)/max(PLGdep(k,i,j),1.E-30_RP)
             PLIdep(k,i,j) =  coef_a_dep*r_taudep_i - coef_b_dep*r_taudep_i
             PLSdep(k,i,j) =  coef_a_dep*r_taudep_s - coef_b_dep*r_taudep_s
             PLGdep(k,i,j) =  coef_a_dep*r_taudep_g - coef_b_dep*r_taudep_g
             PNIdep(k,i,j) = min(0.0_RP, PLIdep(k,i,j)*PLI2NI )
             PNSdep(k,i,j) = min(0.0_RP, PLSdep(k,i,j)*PLS2NS )
             PNGdep(k,i,j) = min(0.0_RP, PLGdep(k,i,j)*PLG2NG )
          end if
          !
       end do
       !
       do k = KS, KE
          if( PLCdep(k,i,j) < -eps )then
             PNCdep(k,i,j) = min(0.0_RP, ((lc(k,i,j)+PLCdep(k,i,j)*dt)*r_xc_ccn - nc(k,i,j))*r_dt )
          else
             PNCdep(k,i,j) = 0.0_RP
          end if
!          if( PLIdep(k,i,j) < -eps )then
!             PNIdep(k,i,j) = min(0.0_RP, ((li(k,i,j)+PLIdep(k,i,j)*dt)*r_xi_ccn - ni(k,i,j))*r_dt )
!          else
!             PNIdep(k,i,j) = 0.0_RP
!          end if
       end do
       !
       do k = KS, KE
          xi(k,i,j) = min(xi_max, max(xi_min, li(k,i,j)/(ni(k,i,j)+ni_min) ))
       end do
       !
       do k = KS, KE
          !
          !--- evaporation/condensation, deposition/sublimation
          !
          r_rvaptem = 1.0_RP/(Rvap*wtem(k,i,j))
          lvsw    = esw(k,i,j)*r_rvaptem
          lvsi    = esi(k,i,j)*r_rvaptem
          dlvsw   = lv(k,i,j)-lvsw
          dlvsi   = lv(k,i,j)-lvsi  ! limiter for esi>1.d0
          dcnd    = dt*(PLCdep(k,i,j)+PLRdep(k,i,j))
          ddep    = dt*(PLIdep(k,i,j)+PLSdep(k,i,j)+PLGdep(k,i,j))
          dep_dqc = 0.0_RP
          dep_dqr = 0.0_RP
          dep_dqi = 0.0_RP
          dep_dqs = 0.0_RP
          dep_dqg = 0.0_RP
          !    always supersaturated
          if     ( (dcnd >  eps) .and. (dlvsw > eps) )then
             fac1    = min(dlvsw,dcnd)/dcnd
             dep_dqc =  dt*PLCdep(k,i,j)*fac1
             dep_dqr =  dt*PLRdep(k,i,j)*fac1
             ! always unsaturated
          else if( (dcnd < -eps) .and. (dlvsw < -eps) )then
             fac1    = max( dlvsw,dcnd )/dcnd
             dep_dqc = max( dt*PLCdep(k,i,j)*fac1, -lc(k,i,j) )
             dep_dqr = max( dt*PLRdep(k,i,j)*fac1, -lr(k,i,j) )
          else
             ! partially unsaturated during timestep
             fac1    = 1.0_RP
             dep_dqc = 0.0_RP
             dep_dqr = 0.0_RP
          end if
          !
          fac2    = 1.0_RP
          !    always supersaturated
          if      ( (ddep >  eps) .and. (dlvsi > eps) )then
             fac2    = min(dlvsi,ddep)/ddep
             dep_dqi = dt*PLIdep(k,i,j)*fac2
             dep_dqs = dt*PLSdep(k,i,j)*fac2
             dep_dqg = dt*PLGdep(k,i,j)*fac2
             ! always unsaturated
          else if ( (ddep < -eps) .and. (dlvsi < -eps) )then
             fac2    = max(dlvsi,ddep)/ddep
             dep_dqi = max(dt*PLIdep(k,i,j)*fac2, -li(k,i,j) )
             dep_dqs = max(dt*PLSdep(k,i,j)*fac2, -ls(k,i,j) )
             dep_dqg = max(dt*PLGdep(k,i,j)*fac2, -lg(k,i,j) )
          else
             ! partially unsaturated during timestep
             fac2    = 1.0_RP
             dep_dqi = dt*PLIdep(k,i,j)
             dep_dqs = dt*PLSdep(k,i,j)
             dep_dqg = dt*PLGdep(k,i,j)
          end if
          ! evaporation always lose number(always negative).
          dep_dnc = max( dt*PNCdep(k,i,j)*fac1, -nc(k,i,j) ) ! ss>0 dep=0, ss<0 dep<0 ! [Add] 11/08/30 T.Mitsui
          dep_dnr = max( dt*PNRdep(k,i,j)*fac1, -nr(k,i,j) ) ! ss>0 dep=0, ss<0 dep<0
          dep_dni = max( dt*PNIdep(k,i,j)*fac2, -ni(k,i,j) ) ! ss>0 dep=0, ss<0 dep<0
          dep_dns = max( dt*PNSdep(k,i,j)*fac2, -ns(k,i,j) ) ! ss>0 dep=0, ss<0 dep<0
          dep_dng = max( dt*PNGdep(k,i,j)*fac2, -ng(k,i,j) ) ! ss>0 dep=0, ss<0 dep<0
          !
          !--- freezing of cloud drop
          !
          frz_dqc = max( dt*(PLChom(k,i,j)+PLChet(k,i,j)), -lc(k,i,j)-dep_dqc ) ! negative value
          frz_dnc = max( dt*(PNChom(k,i,j)+PNChet(k,i,j)), -nc(k,i,j)-dep_dnc ) ! negative value
          fac3    = ( frz_dqc-eps )/( dt*(PLChom(k,i,j)+PLChet(k,i,j))-eps )
          fac4    = ( frz_dnc-eps )/( dt*(PNChom(k,i,j)+PNChet(k,i,j))-eps )
          PLChom(k,i,j) = fac3*PLChom(k,i,j)
          PLChet(k,i,j) = fac3*PLChet(k,i,j)
          PNChom(k,i,j) = fac4*PNChom(k,i,j)
          PNChet(k,i,j) = fac4*PNChet(k,i,j)
          !
          !--- melting
          !
          ! ice change
          mlt_dqi = max( dt*PLImlt(k,i,j), -li(k,i,j)-dep_dqi )  ! negative value
          mlt_dni = max( dt*PNImlt(k,i,j), -ni(k,i,j)-dep_dni )  ! negative value
          ! snow change
          mlt_dqs = max( dt*PLSmlt(k,i,j), -ls(k,i,j)-dep_dqs )  ! negative value
          mlt_dns = max( dt*PNSmlt(k,i,j), -ns(k,i,j)-dep_dns )  ! negative value
          ! graupel change
          mlt_dqg = max( dt*PLGmlt(k,i,j), -lg(k,i,j)-dep_dqg )  ! negative value
          mlt_dng = max( dt*PNGmlt(k,i,j), -ng(k,i,j)-dep_dng )  ! negative value
          !
          !--- freezing of larger droplets
          !
          frz_dqr = max( dt*(PLRhet(k,i,j)), min(0.0_RP, -lr(k,i,j)-dep_dqr) ) ! negative value
          frz_dnr = max( dt*(PNRhet(k,i,j)), min(0.0_RP, -nr(k,i,j)-dep_dnr) ) ! negative value
          !
          fac5         = ( frz_dqr-eps )/( dt*PLRhet(k,i,j)-eps )
          PLRhet(k,i,j) = fac5*PLRhet(k,i,j)
          fac6         = ( frz_dnr-eps )/( dt*PNRhet(k,i,j)-eps )
          PNRhet(k,i,j) = fac6*PNRhet(k,i,j)
          !
          ! water vapor change
          drhogqv = -(dep_dqc+dep_dqi+dep_dqs+dep_dqg+dep_dqr)
          if ( xi(k,i,j) > x_sep )then ! large ice crystals turn into rain by melting
             ! total cloud change
             drhogqc = ( frz_dqc                               + dep_dqc )
             drhognc = ( frz_dnc                               + dep_dnc )
             ! total rain change
             drhogqr = ( frz_dqr - mlt_dqg - mlt_dqs - mlt_dqi + dep_dqr )
             drhognr = ( frz_dnr - mlt_dng - mlt_dns - mlt_dni + dep_dnr )
             !
          else
             ! total cloud change
             drhogqc = ( frz_dqc - mlt_dqi           + dep_dqc )
             drhognc = ( frz_dnc - mlt_dni           + dep_dnc )
             ! total rain change
             drhogqr = ( frz_dqr - mlt_dqg - mlt_dqs + dep_dqr )
             drhognr = ( frz_dnr - mlt_dng - mlt_dns + dep_dnr )
          end if
          ! total ice change
          drhogqi = (-frz_dqc + mlt_dqi           + dep_dqi )
          drhogni = (-frz_dnc + mlt_dni           + dep_dni )
          ! total snow change
          drhogqs = (           mlt_dqs           + dep_dqs )
          drhogns = (           mlt_dns           + dep_dns )
          ! total graupel change
          drhogqg = (-frz_dqr + mlt_dqg           + dep_dqg )
          drhogng = (-frz_dnr + mlt_dng           + dep_dng )
          !
          !--- update
          ! filter for rounding error
          rhogq(k,i,j,I_QV) = max(0.0_RP, rhogq(k,i,j,I_QV) + drhogqv )
          !
          rhogq(k,i,j,I_QC) = max(0.0_RP, rhogq(k,i,j,I_QC) + drhogqc )
          rhogq(k,i,j,I_NC) = max(0.0_RP, rhogq(k,i,j,I_NC) + drhognc )
          rhogq(k,i,j,I_QR) = max(0.0_RP, rhogq(k,i,j,I_QR) + drhogqr )
          rhogq(k,i,j,I_NR) = max(0.0_RP, rhogq(k,i,j,I_NR) + drhognr )
          rhogq(k,i,j,I_QI) = max(0.0_RP, rhogq(k,i,j,I_QI) + drhogqi )
          rhogq(k,i,j,I_NI) = max(0.0_RP, rhogq(k,i,j,I_NI) + drhogni )
          rhogq(k,i,j,I_QS) = max(0.0_RP, rhogq(k,i,j,I_QS) + drhogqs )
          rhogq(k,i,j,I_NS) = max(0.0_RP, rhogq(k,i,j,I_NS) + drhogns )
          rhogq(k,i,j,I_QG) = max(0.0_RP, rhogq(k,i,j,I_QG) + drhogqg )
          rhogq(k,i,j,I_NG) = max(0.0_RP, rhogq(k,i,j,I_NG) + drhogng )
          !
          rhoge(k,i,j) = rhoge(k,i,j) &
               - LHV * drhogqv &
               + LHF * ( drhogqi + drhogqs + drhogqg )
       end do
       !
       do k = KS, KE
          !--- update mixing ratio
          q(k,i,j,I_QV) = rhogq(k,i,j,I_QV) * rrhog(k,i,j)
          q(k,i,j,I_QC) = rhogq(k,i,j,I_QC) * rrhog(k,i,j)
          q(k,i,j,I_QR) = rhogq(k,i,j,I_QR) * rrhog(k,i,j)
          q(k,i,j,I_QI) = rhogq(k,i,j,I_QI) * rrhog(k,i,j)
          q(k,i,j,I_QS) = rhogq(k,i,j,I_QS) * rrhog(k,i,j)
          q(k,i,j,I_QG) = rhogq(k,i,j,I_QG) * rrhog(k,i,j)
          !
          q(k,i,j,I_NC) = rhogq(k,i,j,I_NC) * rrhog(k,i,j)
          q(k,i,j,I_NR) = rhogq(k,i,j,I_NR) * rrhog(k,i,j)
          q(k,i,j,I_NI) = rhogq(k,i,j,I_NI) * rrhog(k,i,j)
          q(k,i,j,I_NS) = rhogq(k,i,j,I_NS) * rrhog(k,i,j)
          q(k,i,j,I_NG) = rhogq(k,i,j,I_NG) * rrhog(k,i,j)
       end do
       !
       do k = KS, KE
          CALC_QDRY( qd(k,i,j), q, k, i, j, iqw )
          CALC_CV( cva(k,i,j), qd(k,i,j), q, k, i, j, iqw, CVdry, AQ_CV )
          CALC_R( Rmoist, q(k,i,j,I_QV), qd(k,i,j), Rdry, Rvap )
          tem(k,i,j) = rhoge(k,i,j) / ( rhog(k,i,j) * cva(k,i,j) )
          pre(k,i,j) = rhog(k,i,j) * Rmoist * tem(k,i,j)
       end do
       !
       do k = KS, KE
          sl_PLCdep(1,i,j) = sl_PLCdep(1,i,j) + dep_dqc*Dz(k)*gsgam2(k,i,j)
          sl_PLRdep(1,i,j) = sl_PLRdep(1,i,j) + dep_dqr*Dz(k)*gsgam2(k,i,j)
          sl_PNRdep(1,i,j) = sl_PNRdep(1,i,j) + dep_dnr*Dz(k)*gsgam2(k,i,j)
       end do
    end do
    end do
    !
    return

  end subroutine update_by_phase_change_kij
  !-------------------------------------------------------------------------------
  subroutine MP_negativefilter( &
       DENS, &
       QTRC  )
    implicit none
    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)

    real(RP) :: diffq(KA,IA,JA)
    real(RP) :: r_xmin

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

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
    do iq = 1, QA
    do k  = KS, KE
       if ( QTRC(k,i,j,iq) < 0.0_RP ) then
          QTRC(k,i,j,iq) = 0.0_RP
       endif
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

    return
  end subroutine MP_negativefilter
  !-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sn14_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QAD)

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do iq = 1, MP_QA
          qhydro = qhydro + QTRC(k,i,j,I_MP2ALL(iq))
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
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
       DENS0  )
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Re   (KA,IA,JA,MP_QAD) ! effective radius
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]

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
       dc_ave(k,i,j) = a_m(I_QC) * xc(k,i,j)**b_m(I_QC)
       dr_ave(k,i,j) = a_m(I_QR) * xr(k,i,j)**b_m(I_QR)
    enddo
    enddo
    enddo

    ! cloud effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rc = 0.5_RP * dc_ave(k,i,j)
       limitsw = 0.5_RP + sign(0.5_RP, rc-rmin_re )
       Re(k,i,j,I_mp_QC) = coef_re(I_QC) * rc * limitsw
    enddo
    enddo
    enddo

    ! rain effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rr = 0.5_RP * dr_ave(k,i,j)
       limitsw = 0.5_RP + sign(0.5_RP, rr-rmin_re )
       Re(k,i,j,I_mp_QR) = coef_re(I_QR) * rr * limitsw
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ri2m(k,i,j) = PI * coef_rea2(I_QI) * QTRC0(k,i,j,I_NI) * a_rea2(I_QI) * xi(k,i,j)**b_rea2(I_QI)
       rs2m(k,i,j) = PI * coef_rea2(I_QS) * QTRC0(k,i,j,I_NS) * a_rea2(I_QS) * xs(k,i,j)**b_rea2(I_QS)
       rg2m(k,i,j) = PI * coef_rea2(I_QG) * QTRC0(k,i,j,I_NG) * a_rea2(I_QG) * xg(k,i,j)**b_rea2(I_QG)
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
       Re(k,i,j,I_mp_QI) = ri3m(k,i,j) / ( ri2m(k,i,j) + zerosw ) * ( 1.0_RP - zerosw )
    enddo
    enddo
    enddo

    ! snow effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       zerosw = 0.5_RP - sign(0.5_RP, rs2m(k,i,j) - r2m_min )
       Re(k,i,j,I_mp_QS) = rs3m(k,i,j) / ( rs2m(k,i,j) + zerosw ) * ( 1.0_RP - zerosw )
    enddo
    enddo
    enddo

    ! graupel effective radius
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       zerosw = 0.5_RP - sign(0.5_RP, rg2m(k,i,j) - r2m_min )
       Re(k,i,j,I_mp_QG) = rg3m(k,i,j) / ( rg2m(k,i,j) + zerosw ) * ( 1.0_RP - zerosw )
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_sn14_EffectiveRadius
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sn14_Mixingratio( &
       Qe,    &
       QTRC0  )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_tracer, only: &
       QAD => QA, &
       MP_QAD => MP_QA
    implicit none

    real(RP), intent(out) :: Qe   (KA,IA,JA,MP_QAD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QAD)    ! tracer mass concentration [kg/kg]

    integer  :: ihydro
    !---------------------------------------------------------------------------

    do ihydro = 1, MP_QA
       Qe(:,:,:,ihydro) = QTRC0(:,:,:,I_MP2ALL(ihydro))
    enddo

    return
  end subroutine ATMOS_PHY_MP_sn14_Mixingratio
  !-----------------------------------------------------------------------------
end module scale_atmos_phy_mp_sn14
!-------------------------------------------------------------------------------
