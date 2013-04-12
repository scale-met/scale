!-------------------------------------------------------------------------------
!
!+  NICAM Double moment Water 6 scheme
!
!-------------------------------------------------------------------------------
module mod_atmos_phy_mp
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines for the sn13w parametrization.
  !
  !       
  !++ Current Corresponding Author : T.Seiki
  ! 
  !++ History: NDW6
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !        0.00   11/10/24 T.Seiki, import from NICAM(11/08/30 ver.)
  !        0.01   12/11/09 Y.Sato,  arrange for NDW2 version(warm rain only)      
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
  !                   mod_mp_nsw6.f90 in NICAM
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_const, only: &
     GRAV   => CONST_GRAV,    &
     PI     => CONST_PI,      &
     UNDEF8 => CONST_UNDEF8,  &
     Rdry   => CONST_Rdry,    &
     CPdry  => CONST_CPdry,   &
     CVdry  => CONST_CVdry,   &
     RovCP  => CONST_RovCP,   &
     CPovCV => CONST_CPovCV,  &
     P00    => CONST_PRE00,   &
     T00    => CONST_TEM00,   &
     Rvap   => CONST_Rvap,    &
     CPvap  => CONST_CPvap,   &
     CVvap  => CONST_CVvap,   &
     CL     => CONST_CL,      &
     CI     => CONST_CI,      &
     EPSvap => CONST_EPSvap,  &
     LHV    => CONST_LH00,    &
     LHF    => CONST_LHF00,   &
     LHS    => CONST_LHS00,   &
     LHV0    => CONST_LH0,    &
     LHF0    => CONST_LHF0,   &
     LHS0    => CONST_LHS0,   &
     LHV00    => CONST_LH00,  &
     LHF00    => CONST_LHF00, &
     PSAT0  => CONST_PSAT0,   &
     EMELT  => CONST_EMELT,   &
     DWATR  => CONST_DWATR
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_setup
  public :: ATMOS_PHY_MP

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: mp_sn13w_init
  private :: mp_sn13w
  private :: MP_terminal_velocity

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  integer, parameter :: HYDRO_MAX  = 3    ! total number of mixing ratio of water
  character(len=32), save :: WLABEL(11)

  ! empirical value from Meyers etal.(1991), 1[/liter] = 1.d3[/m3]
  real(RP), private, parameter :: nqmin(3) = (/ 0.0_RP, 1.E+4_RP, 1.0_RP /) ! [1/m3]
  ! refer to Seifert(2002) (Dr. Thesis, Table.5.1)
  ! max mass, for D_min=79um, 2mm, 5mm, 1cm, 1cm
  real(RP), private, parameter :: xqmax(3) = (/ 0.0_RP, 2.6E-10_RP, 5.0E-6_RP /)! [kg]
  ! SB06, Table 1.
  ! min mass, for D_min=2um, 79um, 10um, 20um, 100um
  real(RP), private, parameter :: xqmin(3) = (/ 0.0_RP, 4.20E-15_RP, 2.60E-10_RP /)! [kg]



  ! for all processes
  ! SB06, Table 1.
  real(RP), private, parameter :: xc_min = 4.20E-15_RP   ! [kg] : min mass, D_min=2um
  real(RP), private, parameter :: xr_min = 2.60E-10_RP   ! [kg] : min mass, D_min=79um
  ! refer to Seifert(2002) (Dr. Thesis, Table.5.1)
  real(RP), private, parameter :: xc_max = 2.6E-10_RP    ! [kg] : max, D_max=79um
  real(RP), private, parameter :: xr_max = 5.00E-6_RP    ! [kg] : max, D_max=2mm
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
  ! empirical filter
  real(RP), private, parameter :: lc_min = xc_min*nc_min
  real(RP), private, parameter :: lr_min = xr_min*nr_min
  !
  real(RP), private, parameter :: x_sep   = 2.6E-10_RP ! boundary mass between cloud and rain
  !
  real(RP), private, parameter :: tem_min=100.0_RP
  real(RP), private, parameter :: rho_min=1.E-5_RP     ! 3.e-3 is lower limit recognized in many experiments.
  ! for Seifert(2008)
  ! work parameter for gamma function, imported from MISC_gammafunc
!  real(RP), private, parameter :: gfac_coef(6)=(/&
!       +76.180091729471460_RP, -86.505320329416770_RP, +24.014098240830910_RP,&
!       -1.2317395724501550_RP, +0.1208650973866179E-2_RP, -0.5395239384953E-5_RP /)
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
  real(RP), private, save :: slope_A(HYDRO_MAX)
  ! coefficeint of weighted ventilation effect.
  ! large, and small branch is by PK97(13-60),(13-61),(13-88),(13-89)
  real(RP), private, save :: ah_vent  (HYDRO_MAX,2) !
  real(RP), private, save :: bh_vent  (HYDRO_MAX,2) !
  real(RP), private, save :: ah_vent0 (HYDRO_MAX,2) !
  real(RP), private, save :: bh_vent0 (HYDRO_MAX,2) !
  real(RP), private, save :: ah_vent1 (HYDRO_MAX,2) !
  real(RP), private, save :: bh_vent1 (HYDRO_MAX,2) !
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$  ! coefficient of collision growth
!!$  real(RP), private, save :: delta_b0 (HYDRO_MAX)
!!$  real(RP), private, save :: delta_b1 (HYDRO_MAX)
!!$  real(RP), private, save :: delta_ab0(HYDRO_MAX,HYDRO_MAX)
!!$  real(RP), private, save :: delta_ab1(HYDRO_MAX,HYDRO_MAX)
!!$  !
!!$  real(RP), private, save :: theta_b0 (HYDRO_MAX)
!!$  real(RP), private, save :: theta_b1 (HYDRO_MAX)
!!$  real(RP), private, save :: theta_ab0(HYDRO_MAX,HYDRO_MAX)
!!$  real(RP), private, save :: theta_ab1(HYDRO_MAX,HYDRO_MAX)

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
  subroutine ATMOS_PHY_MP_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_MP
    implicit none

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       doautoconversion, &
       doprecipitation,  &
       MP_ssw_lim

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Cloud Microphisics]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Wrapper for NDW6'

    if ( ATMOS_TYPE_PHY_MP /= 'NDW3' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_MP is not NDW6. Check!'
       call PRC_MPIstop
    endif

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

    WLABEL( 1) = "VAPOR"
    WLABEL( 2) = "CLOUD"
    WLABEL( 3) = "RAIN"
    WLABEL( 4) = "CLOUD_NUM"
    WLABEL( 5) = "RAIN_NUM"

    call mp_sn13w_init( IA, JA )

    MP_NSTEP_SEDIMENTATION  = ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(ntmax_sedimentation,kind=8)
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
  end subroutine ATMOS_PHY_MP_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Microphysics'

    call TIME_rapstart('MP0 Setup')
    call MP_negativefilter
    call TIME_rapend  ('MP0 Setup')

    call mp_sn13w

    call TIME_rapstart('MP6 Filter')
    call MP_negativefilter
    call TIME_rapend  ('MP6 Filter')

    call ATMOS_vars_fillhalo

    call ATMOS_vars_total

    return
  end subroutine ATMOS_PHY_MP

  !-----------------------------------------------------------------------------
  subroutine mp_sn13w_init ( IAA, JA )
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
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
    namelist /nm_mp_sn13w_init/       &
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
    namelist /nm_mp_sn13w_particles/ &
         a_m, b_m, alpha_v, beta_v, gamma_v, &
         alpha_vn, beta_vn,    &
         a_area, b_area, cap,  &
         nu, mu,               &
!         opt_M96_column_ice,   &
!         opt_M96_ice,          &
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
    slope_A(:)     = UNDEF8
    coef_lambda(:) = UNDEF8
    coef_vt0(:,:)  = UNDEF8
    coef_vt1(:,:)  = UNDEF8
!    delta_b0(:)    = UNDEF8
!    delta_b1(:)    = UNDEF8
!    delta_ab0(:,:) = UNDEF8
!    delta_ab1(:,:) = UNDEF8
!    theta_b0(:)    = UNDEF8
!    theta_b1(:)    = UNDEF8
!    theta_ab0(:,:) = UNDEF8
!    theta_ab1(:,:) = UNDEF8
    !
    ah_vent(:,:)   = UNDEF8
    ah_vent0(:,:)  = UNDEF8
    ah_vent1(:,:)  = UNDEF8
    bh_vent(:,:)   = UNDEF8
    bh_vent0(:,:)  = UNDEF8
    bh_vent1(:,:)  = UNDEF8
    !

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[NDW6]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_sn13w_init,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist nm_mp_sn13w_init. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=nm_mp_sn13w_init)

    !
    ! default setting
    !
    ! Area parameters with mks unit originated by Mitchell(1996)
    a_area(I_QC) = PI/4.0_RP ! sphere
    a_area(I_QR) = PI/4.0_RP ! sphere
    b_area(I_QC) = 2.0_RP
    b_area(I_QR) = 2.0_RP
    !
    ! Seifert and Beheng(2006), Table. 1 or List of symbols
    !----------------------------------------------------------
    ! Diameter-Mass relationship
    ! D = a * x^b
    a_m(I_QC) = 0.124_RP
    a_m(I_QR) = 0.124_RP
    b_m(I_QC) = 1.0_RP/3.0_RP
    b_m(I_QR) = 1.0_RP/3.0_RP
    !----------------------------------------------------------
    ! Terminal velocity-Mass relationship
    ! vt = alpha * x^beta * (rho0/rho)^gamma
    alpha_v(I_QC,:)= 3.75E+5_RP
    alpha_v(I_QR,:)= 159.0_RP ! not for sedimantation
    beta_v(I_QC,:) = 2.0_RP/3.0_RP
    beta_v(I_QR,:) = 0.266_RP ! not for sedimantation
    gamma_v(I_QC)  = 1.0_RP
    ! This is high Reynolds number limit(Beard 1980)
    gamma_v(I_QR)  = 1.0_RP/2.0_RP
    !----------------------------------------------------------
    ! DSD parameters
    ! f(x) = A x^nu exp( -lambda x^mu )
    ! Gamma Disribution           : mu=1  , nu:arbitrary
    ! Marshall-Palmer Distribution: mu=1/3, nu:-2/3
    ! In the case of MP, f(D) dD = f(x)dx
    !                    f(x)    = c * f(D)/D^2 (c:coefficient)
    nu(I_QC) =  1.0_RP          ! arbitrary for Gamma
    nu(I_QR) = -1.0_RP/3.0_RP     ! nu(diameter)=1, equilibrium condition.
    !
    mu(I_QC) = 1.0_RP           ! Gamma
    mu(I_QR) = 1.0_RP/3.0_RP      ! Marshall Palmer
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
    !
    alpha_vn(:,:) = alpha_v(:,:)
    beta_vn(:,:) = beta_v(:,:)
    !------------------------------------------------------------------------
    !
    ! additional setting
    !

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=nm_mp_sn13w_particles,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist nm_mp_sn13w_particles. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=nm_mp_sn13w_particles)

    !
    ! area-diameter relation => area-mass relation
    ax_area(I_QC:I_QR) = a_area(I_QC:I_QR)*a_m(I_QC:I_QR)**b_area(I_QC:I_QR)
    bx_area(I_QC:I_QR) = b_area(I_QC:I_QR)*b_m(I_QC:I_QR)
    !
    ! radius of equivalent area - m ass relation
    ! pi*rea**2 = ax_area*x**bx_area
    a_rea(I_QC:I_QR)   = sqrt(ax_area(I_QC:I_QR)/PI)
    b_rea(I_QC:I_QR)   = bx_area(I_QC:I_QR)/2.0_RP
    a_rea2(I_QC:I_QR)  = a_rea(I_QC:I_QR)**2
    b_rea2(I_QC:I_QR)  = b_rea(I_QC:I_QR)*2.0_RP
    a_rea3(I_QC:I_QR)  = a_rea(I_QC:I_QR)**3
    b_rea3(I_QC:I_QR)  = b_rea(I_QC:I_QR)*3.0_RP
    !
    a_d2vt(I_QC:I_QR)=alpha_v(I_QC:I_QR,2)*(1.0_RP/alpha_v(I_QC:I_QR,2))**(beta_v(I_QC:I_QR,2)/b_m(I_QC:I_QR))
    b_d2vt(I_QC:I_QR)=(beta_v(I_QC:I_QR,2)/b_m(I_QC:I_QR))
    !
    ! Calculation of Moment Coefficient
    !
    allocate( w1(I_QC:I_QR), w2(I_QC:I_QR), w3(I_QC:I_QR), w4(I_QC:I_QR) )
    allocate( w5(I_QC:I_QR), w6(I_QC:I_QR), w7(I_QC:I_QR), w8(I_QC:I_QR) )
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
    do iw=I_QC,I_QR
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
    do iw=I_QC,I_QR
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
    do iw=I_QC,I_QR
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
    do iw=I_QC,I_QR
       w1(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
       w2(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
       coef_A(iw)      = mu(iw)/w1(iw)
       slope_A(iw)     = w1(iw)
       coef_lambda(iw) = (w1(iw)/w2(iw))**(-mu(iw))
    end do
    !-------------------------------------------------------
    ! coefficient for terminal velocity in sedimentation
    ! SB06(78)
    do ia=1,2
       do iw=I_QC,I_QR
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
    do iw=I_QC,I_QR
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
    bh_vent(I_QC,1:2) = (/0.0000_RP,0.0000_RP/)
    bh_vent(I_QR,1:2) = (/0.108_RP,0.308_RP/)
    !
    do iw=I_QC,I_QR
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
       endif
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
       endif
    end do
    do iw=I_QC,I_QR
       n = 0
       if( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,1) + n) < eps_gamma )then
          flag_vent0(iw)=.false.
       endif
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
       endif
       !
       n = 1
       if( (nu(iw) + 1.5_RP*b_m(iw) + 0.5_RP*beta_v(iw,1) + n) < eps_gamma )then
          flag_vent1(iw)=.false.
       endif
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
       endif
    end do
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    !-------------------------------------------------------
!!$    ! coefficient for collision process
!!$    ! stochastic coefficient for collision cross section
!!$    ! sb06 (90) -- self collection
!!$    do iw=I_QC,I_QR
!!$       n = 0
!!$       w1(iw) = gammafunc( (2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
!!$       w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
!!$       delta_b0(iw) = w1(iw)/w2(iw) &
!!$            *( w2(iw)/w3(iw) )**(2.0_RP*b_rea(iw) + n)
!!$       n = 1
!!$       w1(iw) = gammafunc( (2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       delta_b1(iw) = w1(iw)/w2(iw) &
!!$            *( w2(iw)/w3(iw) )**(2.0_RP*b_rea(iw) + n)
!!$    end do
!!$    ! stochastic coefficient for collision cross section
!!$    ! sb06(91) -- riming( collide with others )
!!$    do iw=I_QC,I_QR
!!$       n = 0
!!$       w1(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w2(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
!!$       w3(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
!!$       w4(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP    )/mu(iw) )
!!$       n = 1
!!$       w5(iw) = gammafunc( (b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$    end do
!!$    ! ia > ib ( larger particles "a" catch smaller particles "b" )
!!$    do ia=I_QC,I_QR
!!$       do ib=I_QC,I_QR
!!$          n=0 !
!!$          ! NOTE, collected  particle has a moment of n.
!!$          !       collecting particle has only number(n=0).
!!$          delta_ab0(ia,ib) = 2.0_RP*(w1(ib)/w2(ib))*(w4(ia)/w2(ia)) &
!!$               * ( w2(ib)/w3(ib) )**(b_rea(ib)+n) &
!!$               * ( w2(ia)/w3(ia) )**(b_rea(ia)  )
!!$          n=1 !
!!$          delta_ab1(ia,ib) = 2.0_RP*(w5(ib)/w2(ib))*(w4(ia)/w2(ia)) &
!!$               * ( w2(ib)/w3(ib) )**(b_rea(ib)+n) &
!!$               * ( w2(ia)/w3(ia) )**(b_rea(ia)  )
!!$       end do
!!$    end do
!!$    ! stochastic coefficient for terminal velocity
!!$    ! sb06(92) -- self collection
!!$    ! assuming equivalent area circle.
!!$    do iw=I_QC,I_QR
!!$       n = 0
!!$       w1(iw) = gammafunc( (2.0_RP*beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w2(iw) = gammafunc( (                        2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w3(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
!!$       w4(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
!!$       theta_b0(iw) = w1(iw)/w2(iw) * ( w3(iw)/w4(iw) )**(2.0_RP*beta_v(iw,2))
!!$       n = 1
!!$       w1(iw) = gammafunc( (2.0_RP*beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w2(iw) = gammafunc( (                        2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       theta_b1(iw) = w1(iw)/w2(iw) * ( w3(iw)/w4(iw) )**(2.0_RP*beta_v(iw,2))
!!$    end do
!!$    !
!!$    ! stochastic coefficient for terminal velocity
!!$    ! sb06(93) -- riming( collide with others )
!!$    do iw=I_QC,I_QR
!!$       n = 0
!!$       w1(iw) = gammafunc( (beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w2(iw) = gammafunc( (                   2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w3(iw) = gammafunc( (beta_v(iw,2) + 2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP    )/mu(iw) )
!!$       w4(iw) = gammafunc( (                   2.0_RP*b_rea(iw) + nu(iw) + 1.0_RP    )/mu(iw) )
!!$       !
!!$       w5(iw) = gammafunc( (nu(iw) + 1.0_RP)/mu(iw) )
!!$       w6(iw) = gammafunc( (nu(iw) + 2.0_RP)/mu(iw) )
!!$       n = 1
!!$       w7(iw) = gammafunc( (beta_v(iw,2) + b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$       w8(iw) = gammafunc( (                   b_rea(iw) + nu(iw) + 1.0_RP + n)/mu(iw) )
!!$    end do
!!$    ! ia > ib ( larger particles "a" catch smaller particles "b" )
!!$    do ia=I_QC,I_QR
!!$       do ib=I_QC,I_QR
!!$          theta_ab0(ia,ib) = 2.0_RP * (w1(ib)/w2(ib))*(w3(ia)/w4(ia)) &
!!$               * (w5(ia)/w6(ia))**beta_v(ia,2) &
!!$               * (w5(ib)/w6(ib))**beta_v(ib,2)
!!$          theta_ab1(ia,ib) = 2.0_RP * (w7(ib)/w8(ib))*(w3(ia)/w4(ia)) &
!!$               * (w5(ia)/w6(ia))**beta_v(ia,2) &
!!$               * (w5(ib)/w6(ib))**beta_v(ib,2)
!!$       end do
!!$    end do
!!$
!!$    deallocate(w1,w2,w3,w4,w5,w6,w7,w8)

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

!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "delta_b0    ",delta_b0(:)
!!$    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "delta_b1    ",delta_b1(:)
!!$    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "theta_b0    ",theta_b0(:)
!!$    if( IO_L ) write(IO_FID_LOG,'(a,100e16.6)') "theta_b1    ",theta_b1(:)
!!$
!!$    do ia=QQS,QQE
!!$       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "delta0(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab0(ia,ib),ib=QQS,QQE)
!!$    enddo
!!$    do ia=QQS,QQE
!!$       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "delta1(a,b)=(",trim(WLABEL(ia)),",b)=",(delta_ab1(ia,ib),ib=QQS,QQE)
!!$    enddo
!!$    do ia=QQS,QQE
!!$       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "theta0(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab0(ia,ib),ib=QQS,QQE)
!!$    enddo
!!$    do ia=QQS,QQE
!!$       if( IO_L ) write(IO_FID_LOG,'(a,a10,a,100e16.6)') "theta1(a,b)=(",trim(WLABEL(ia)),",b)=",(theta_ab1(ia,ib),ib=QQS,QQE)
!!$    enddo

    allocate(nc_uplim_d(1,IAA,JA))
    nc_uplim_d(:,:,:) = 150.d6

    return
  end subroutine mp_sn13w_init
  !-----------------------------------------------------------------------------
  subroutine mp_sn13w
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP, &
       ct => TIME_NOWDAYSEC
    use mod_grid, only: &
       z    => GRID_CZ, &
       dz   => GRID_CDZ
    use mod_atmos_precipitation, only: &
       precipitation => ATMOS_PRECIPITATION
    use mod_atmos_thermodyn, only: &
       !-- For kji
       thrmdyn_qd      => ATMOS_THERMODYN_qd, &
       thrmdyn_cv      => ATMOS_THERMODYN_cv, &
       thrmdyn_cp      => ATMOS_THERMODYN_cp, &
       thrmdyn_tempre  => ATMOS_THERMODYN_tempre, &
       thrmdyn_tempre2 => ATMOS_THERMODYN_tempre2, &
       CVw => AQ_CV
    use mod_atmos_saturation, only: &
       moist_psat_water    => ATMOS_SATURATION_psat_liq,     &
       moist_psat_ice      => ATMOS_SATURATION_psat_ice,     &
       moist_dqsw_dtem_rho => ATMOS_SATURATION_dqsw_dtem_rho
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none

    !
    ! primary variables
    !
    real(RP) :: rhogvx_d(KA,IA,JA)
    real(RP) :: rhogvy_d(KA,IA,JA)
    real(RP) :: rhogw_d (KA,IA,JA)
    real(RP) :: rhogq_d (KA,IA,JA,QA)
    real(RP) :: rhog_d  (KA,IA,JA)
    real(RP) :: th_d    (KA,IA,JA)
    !
    ! diagnostic variables
    !
    real(RP) :: qd_d(KA,IA,JA)
    real(RP) :: q_d(KA,IA,JA,QA)
    real(RP) :: rho_d(KA,IA,JA)
    real(RP) :: tem_d(KA,IA,JA)
    real(RP) :: pre_d(KA,IA,JA)
    real(RP) :: rhoge_d(KA,IA,JA)
    real(RP) :: w_d(KA,IA,JA)
    !
    real(RP) :: lv_d(KA,IA,JA)
    real(RP) :: lc_d(KA,IA,JA), nc_d(KA,IA,JA)
    real(RP) :: lr_d(KA,IA,JA), nr_d(KA,IA,JA)
    !
    real(RP) :: xc_d(KA,IA,JA) ! lc/nc
    real(RP) :: xr_d(KA,IA,JA) ! lr/nr
    !
    real(RP) :: dc_xa_d(KA,IA,JA) !
    real(RP) :: dr_xa_d(KA,IA,JA) !
    real(RP) :: vt_xa_d(KA,IA,JA,HYDRO_MAX,2) !

    real(RP) :: wtem_d(KA,IA,JA)      ! filtered temperature
    real(RP) :: esw_d(KA,IA,JA)       ! saturated vapor pressure(water)
    real(RP) :: esi_d(KA,IA,JA)       ! saturated vapor pressure(ice)
    !
    real(RP) :: rho_fac
    real(RP) :: rho_fac_c_d(KA,IA,JA) ! factor for cloud
    real(RP) :: rho_fac_r_d(KA,IA,JA) !            rain
    real(RP) :: cva_d(KA,IA,JA)    !
    real(RP) :: cpa_d(KA,IA,JA)       ! [Add] 09/08/18 T.Mitsui
    !
    real(RP) :: drhogqv               ! d (rho*qv*gsgam2)
    real(RP) :: drhogqc, drhognc      !        qc, nc
    real(RP) :: drhogqr, drhognr      !        qr, nr

    ! production rate of nucleation
    real(RP) :: PLCccn_d(KA,IA,JA), PNCccn_d(KA,IA,JA)
    ! production rate of vapor deposition
    real(RP) :: PLRdep_d(KA,IA,JA), PNRdep_d(KA,IA,JA)
    real(RP) :: PLCdep_d(KA,IA,JA), PNCdep_d(KA,IA,JA)

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

    real(RP) :: fac1
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    real(RP) :: gc_dqc, gc_dnc
!!$    real(RP) :: sc_dqc, sc_dnc
!!$    real(RP) :: ic_dqc, ic_dnc
!!$    real(RP) :: rg_dqg, rg_dng
!!$    real(RP) :: rg_dqr, rg_dnr
!!$    real(RP) :: rs_dqr, rs_dnr, rs_dqs, rs_dns
!!$    real(RP) :: ri_dqr, ri_dnr
!!$    real(RP) :: ri_dqi, ri_dni
!!$    real(RP) :: ii_dqi, ii_dni
!!$    real(RP) :: is_dqi, is_dni, ss_dns
!!$    real(RP) :: gs_dqs, gs_dns, gg_dng
!!$    ! mixed-phase collection process total plus(clp_), total minus(clm_)
!!$    real(RP) :: clp_dqc, clp_dnc, clm_dqc, clm_dnc
!!$    real(RP) :: clp_dqr, clp_dnr, clm_dqr, clm_dnr
!!$    real(RP) :: clp_dqi, clp_dni, clm_dqi, clm_dni
!!$    real(RP) :: clp_dqs, clp_dns, clm_dqs, clm_dns
!!$    real(RP) :: clp_dqg, clp_dng, clm_dqg, clm_dng
!!$    real(RP) :: fac1, fac3, fac4, fac6, fac7, fac9, fac10
!!$    ! production rate of partial conversion(ice, snow => graupel)
!!$    real(RP) :: PLIcon_d(KA,IA,JA), PNIcon_d(KA,IA,JA)
!!$    real(RP) :: PLScon_d(KA,IA,JA), PNScon_d(KA,IA,JA)
!!$    real(RP) :: pco_dqi, pco_dni
!!$    real(RP) :: pco_dqs, pco_dns
!!$    real(RP) :: pco_dqg, pco_dng
!!$    real(RP) :: eml_dqc, eml_dnc
!!$    real(RP) :: eml_dqr, eml_dnr
!!$    real(RP) :: eml_dqi, eml_dni
!!$    real(RP) :: eml_dqs, eml_dns
!!$    real(RP) :: eml_dqg, eml_dng
!!$    ! production rate of ice multiplication by splintering
!!$    real(RP) :: PLGspl_d(KA,IA,JA)
!!$    real(RP) :: PLSspl_d(KA,IA,JA)
!!$    real(RP) :: PNIspl_d(KA,IA,JA)
!!$    real(RP) :: spl_dqi, spl_dni
!!$    real(RP) :: spl_dqg, spl_dqs

    real(RP) :: rrhog_d(KA,IA,JA)

    !-----------------------------------------------
    ! work for explicit supersaturation modeling
    !-----------------------------------------------
    real(RP) :: dTdt_equiv_d(KA,IA,JA)    !
    real(RP) :: dlvsi_d(KA,IA,JA)
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
    real(RP) :: qke(IJA,KA)
    real(RP) :: qke_d(KA,IA,JA)

    real(RP), parameter :: eps       = 1.E-30_RP
    real(RP), parameter :: eps_qv    = 1.E-50_RP
    real(RP), parameter :: eps_rhoge = 1.E-50_RP
    real(RP), parameter :: eps_rhog  = 1.E-50_RP
    real(RP) :: r_dt
    real(RP) :: wdt, r_wdt
    real(RP) :: r_ntmax
    integer :: ntdiv

    real(RP) :: rhoe(KA,IA,JA)
    real(RP) :: temp(KA,IA,JA)
    real(RP) :: pres(KA,IA,JA)
    real(RP) :: rhoq(KA,IA,JA,QA)

    real(RP) :: velw(KA,IA,JA,QA)
    real(RP) :: flux_rain (KA,IA,JA)
    real(RP) :: flux_snow (KA,IA,JA)
    real(RP) :: wflux_rain(KA,IA,JA)
    real(RP) :: wflux_snow(KA,IA,JA)
    integer :: step

    integer :: k, i, j, iq, ij
    !---------------------------------------------------------------------------

    r_dt = 1.0_RP / dt

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
    call TIME_rapstart('MPX ijkconvert')

    do j = 1, JA
    do i = 1, IA
       sl_PLCdep_d(1,i,j) = 0.0_RP
       sl_PLRdep_d(1,i,j) = 0.0_RP
       sl_PNRdep_d(1,i,j) = 0.0_RP

!OCL XFILL
       do k = KS, KE
          rhog_d  (k,i,j) = DENS(k,i,j)
       enddo
!OCL XFILL
       do k = KS, KE
          rhogw_d (k,i,j) = MOMZ(k,i,j)
       enddo
       do k = KS, KE
          th_d    (k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       enddo
       rhog_d  (1:KS-1,i,j) = rhog_d  (KS,i,j)
       rhogw_d (1:KS-1,i,j) = rhogw_d (KS,i,j)
       th_d    (1:KS-1,i,j) = th_d    (KS,i,j)

       rhog_d  (KE+1:KA,i,j) = rhog_d  (KE,i,j)
       rhogw_d (KE+1:KA,i,j) = rhogw_d (KE,i,j)
       th_d    (KE+1:KA,i,j) = th_d    (KE,i,j)
    enddo
    enddo

    do iq = 1, QA
       do j = 1, JA
       do i = 1, IA
          do k = KS, KE
             q_d(k,i,j,iq) = QTRC(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo

    do iq = 1, QA
       do j = 1, JA
       do i = 1, IA
          q_d(1:KS-1, i,j  ,iq) = q_d(KS,i,j,iq)
          q_d(KE+1:KA,i,j  ,iq) = q_d(KE,i,j,iq)
       enddo
       enddo
    enddo

    do iq = 1, QA
    do k  = 1, KA
    do j = 1, JA
    do i = 1, IA
       rhogq_d(k,i,j,iq) = rhog_d(k,i,j) * q_d(k,i,j,iq)
    enddo
    enddo
    enddo
    enddo

    do k  = 2, KA
    do j = 1, JA
    do i = 1, IA
       rho_d  (k,i,j) = rhog_d(k,i,j)
       rrhog_d(k,i,j) = 1.0_RP / rhog_d(k,i,j)
       w_d    (k,i,j) = ( rhogw_d(k,i,j) + rhogw_d(k-1,i,j) ) * rrhog_d(k,i,j)
    enddo
    enddo
    enddo

    call thrmdyn_qd( qd_d, q_d )
    call thrmdyn_cv( cva_d, q_d, qd_d )
    call thrmdyn_tempre2( tem_d, pre_d, rho_d, th_d, qd_d, q_d )

    do k = 1, KA
    do j = 1, JA
    do i = 1, IA
       rhoge_d(k,i,j) = rhog_d(k,i,j) * tem_d(k,i,j) * cva_d(k,i,j)
    enddo
    enddo
    enddo

    if( opt_debug_tem ) call debug_tem_kij( 1, tem_d(:,:,:), rho_d(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call TIME_rapend  ('MPX ijkconvert')

    call TIME_rapstart('MP1 Nucleation')

    do k = KS, KE
    do j = 1,  JA
    do i = 1,  IA
       rho_fac         = rho_0 / max(rho_d(k,i,j),rho_min)
       rho_fac_c_d(k,i,j) = rho_fac**gamma_v(I_QC)
       rho_fac_r_d(k,i,j) = rho_fac**gamma_v(I_QR)
    enddo
    enddo
    enddo
    do j=1, JA
    do i=1, IA
       rho_fac_c_d(1:KS-1,i,j)    = 1.0_RP
       rho_fac_r_d(1:KS-1,i,j)    = 1.0_RP

       rho_fac_c_d(KE+1:KA,i,j) = rho_fac_c_d(KE,i,j)
       rho_fac_r_d(KE+1:KA,i,j) = rho_fac_r_d(KE,i,j)
    end do
    end do

    call thrmdyn_cp( cpa_d, q_d, qd_d )
    do k = 1, KA
    do j = 1, JA
    do i = 1, IA
       wtem_d(k,i,j) = max(tem_d(k,i,j), tem_min)
    enddo
    enddo
    enddo

    do k = 1, KA
    do j = 1, JA
    do i = 1, IA
       qke_d       (k,i,j) = 0.0_RP ! 2*TKE
       dTdt_equiv_d(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do k = 1, KA
    do j = 1, JA
    do i = 1, IA
       lv_d(k,i,j) = rho_d(k,i,j)*q_d(k,i,j,I_QV)
       nc_d(k,i,j) = max( 0.0_RP, rho_d(k,i,j)*q_d(k,i,j,I_NC) )
    enddo
    enddo
    enddo

    call nucleation_kij(     &
         IA, JA, KA,     & ! in
         IS, IE,         & ! in
         JS, JE,         & ! in
         KS, KE,         & ! in
         z, w_d,         & ! in
         rho_d, wtem_d, pre_d, & ! in
         lv_d, nc_d,     & ! in
         PNCccn_d, PLCccn_d, & ! out
         cva_d, cpa_d,       & ! in
         dTdt_equiv_d,     & ! in
         qke_d,            & ! in
         flag_history_in,& ! in
         ct, dt         ) ! in

    do k = KS, KE
    do j = 1,  JA
    do i = 1,  IA
       ! nucleation
       drhogqc = dt * PLCccn_d(k,i,j)
       drhognc = dt * PNCccn_d(k,i,j)
       drhogqv = max( -rhogq_d(k,i,j,I_QV), -drhogqc )
       fac1    = drhogqv / min( -drhogqc, -eps ) ! limiting coefficient

       rhogq_d(k,i,j,I_QV) = rhogq_d(k,i,j,I_QV) + drhogqv
       rhogq_d(k,i,j,I_QC) = max(0.0_RP, rhogq_d(k,i,j,I_QC) + drhogqc*fac1)
       rhogq_d(k,i,j,I_NC) = max(0.0_RP, rhogq_d(k,i,j,I_NC) + drhognc)

       ! cloud number concentration filter
       rhoge_d(k,i,j) = rhoge_d(k,i,j) - LHV * drhogqv 
    enddo
    enddo
    enddo

    do k  = KS, KE
    do j = 1, JA
    do i = 1, IA
       q_d(k,i,j,I_QV) = rhogq_d(k,i,j,I_QV) * rrhog_d(k,i,j)
       q_d(k,i,j,I_QC) = rhogq_d(k,i,j,I_QC) * rrhog_d(k,i,j)
       q_d(k,i,j,I_NC) = rhogq_d(k,i,j,I_NC) * rrhog_d(k,i,j)
    enddo
    enddo
    enddo

    call thrmdyn_qd( qd_d, q_d )
    call thrmdyn_cv( cva_d, q_d, qd_d )

    do k  = KS, KE
    do j = 1,    JA
    do i = 1,    IA
       tem_d(k,i,j) = rhoge_d(k,i,j) / ( rhog_d(k,i,j) * cva_d(k,i,j) )
       pre_d(k,i,j) = rho_d(k,i,j)*( qd_d(k,i,j)*Rdry+q_d(k,i,j,I_QV)*Rvap )*tem_d(k,i,j)
    enddo
    enddo
    enddo

!    if( opt_debug )     call debugreport_nucleation
    if( opt_debug_tem ) call debug_tem_kij( 2, tem_d(:,:,:), rho_d(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call TIME_rapend  ('MP1 Nucleation')
    !----------------------------------------------------------------------------
    !
    ! 2.Phase change: Freezing, Melting, Vapor deposition
    !
    !----------------------------------------------------------------------------
    call TIME_rapstart('MP2 Phase change')

    ! parameter setting
    wdt=dt
    r_wdt=1.0_RP/wdt

!       if(  ntdiv     == ntmax_phase_change  )then
          flag_history_in=.true.
!       endif
       !
       do k=1, KA
          do j=1, JA
          do i=1, IA
             lr_d(k,i,j)     = rhogq_d(k,i,j,I_QR)
             nr_d(k,i,j)     = rhogq_d(k,i,j,I_NR)
             xr_d(k,i,j)     = max(xr_min,  min(xr_max, lr_d(k,i,j)/(nr_d(k,i,j)+nr_min) ))
             dr_xa_d(k,i,j)  = a_m(I_QR)*xr_d(k,i,j)**b_m(I_QR)
             vt_xa_d(k,i,j,I_QR,1) = alpha_v(I_QR,1)*(xr_d(k,i,j)**beta_v(I_QR,1))*rho_fac_r_d(k,i,j)
             vt_xa_d(k,i,j,I_QR,2) = vt_xa_d(k,i,j,I_QR,1)
          end do
          end do
       end do
       !
       do k=1, KA
          do j=1, JA
          do i=1, IA
             !! Following values shoud be already filtered to be non-zero before sbroutine was called.
             ! Mass concentration [kg/m3]
             lv_d(k,i,j)     = rhogq_d(k,i,j,I_QV)
             !
             lc_d(k,i,j)     = rhogq_d(k,i,j,I_QC)
             ! Number concentration[/m3] (should be filtered to avoid zero division.)
             nc_d(k,i,j)     = rhogq_d(k,i,j,I_NC)
             ! Mass of mean particle [kg]
             ! SB06(94)
             !
             xc_d(k,i,j)     = min(xc_max, max(xc_min, lc_d(k,i,j)/(nc_d(k,i,j)+nc_min) ))
             ! diamter of average mass
             ! SB06(32)
             dc_xa_d(k,i,j)  = a_m(I_QC)*xc_d(k,i,j)**b_m(I_QC)
             ! terminal velocity of average mass
             vt_xa_d(k,i,j,I_QC,1) = alpha_v(I_QC,1)*(xc_d(k,i,j)**beta_v(I_QC,1))*rho_fac_c_d(k,i,j)
             vt_xa_d(k,i,j,I_QC,2) = alpha_v(I_QC,2)*(xc_d(k,i,j)**beta_v(I_QC,2))*rho_fac_c_d(k,i,j)
          end do
          end do
       end do
       !
       wtem_d(:,:,:) = max(tem_d(:,:,:), tem_min)
       call moist_psat_water( esw_d, wtem_d )
       !
       call thrmdyn_qd( qd_d, q_d )
       !
       call dep_vapor_melt_ice_kij( &
            IA, JA, KA,           & ! in
            IS, IE,               & ! in
            JS, JE,               & ! in
            KS, KE,               & ! in
            PLCdep_d,             & ! out
            PLRdep_d, PNRdep_d,   & ! out
            rho_d, wtem_d, pre_d, lv_d,& ! in
            qd_d,                 & ! in
!!!!!!!! mod xxxxxxxxxxxxxxxxxxxx T.Seiki             
!!$         esw_d,                & ! in
            esw_d,                & ! in
            nc_d,                 & ! in
            nr_d,                 & ! in
            xc_d,                 & ! in
            xr_d,                 & ! in
            vt_xa_d,              & ! in
            dc_xa_d,              & ! in
            dr_xa_d               ) ! in
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
            rhog_d,                        & ! in
            rhoge_d,                       & ! inout
            rhogq_d, q_d,                  & ! inout
            tem_d, pre_d, rho_d,           & ! inout
            cva_d,                         & ! out
!!!!!!!! mod xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$         esw_d, esi_d, LV_d,            & ! in
            esw_d, LV_d,            & ! in
            LC_d, LR_d,                    & ! in
            NC_d, NR_d,                    & ! in
            !+++ tendency terms
            ! condensation/evaporation, deposition/sublimation
            PLCdep_d,           & ! inout
            PLRdep_d, PNRdep_d, & ! inout
            ! melting term
            flag_history_in,    & ! in
            sl_PLCdep_d,        & ! out
            sl_PLRdep_d, sl_PNRdep_d ) ! out

!       if( opt_debug )     call debugreport_phasechange
       if( opt_debug_tem ) call debug_tem_kij( 3, tem_d(:,:,:), rho_d(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call TIME_rapend  ('MP2 Phase change')
    !----------------------------------------------------------------------------
    !
    ! 3.Collection process
    !
    !----------------------------------------------------------------------------
    call TIME_rapstart('MP3 Collection')

    wdt = dt
    flag_history_in=.true.
    do k = 1, KA
    do j = 1, JA
    do i = 1, IA
       ! Mass concentration [kg/m3]
       lv_d(k,i,j) = rhogq_d(k,i,j,I_QV)
       lc_d(k,i,j) = rhogq_d(k,i,j,I_QC)
       lr_d(k,i,j) = rhogq_d(k,i,j,I_QR)
       ! Number concentration[/m3]
       nc_d(k,i,j) = rhogq_d(k,i,j,I_NC)
       nr_d(k,i,j) = rhogq_d(k,i,j,I_NR)
       ! Mass of mean particle [kg]
       xc_d(k,i,j) = min(xc_max, max(xc_min, lc_d(k,i,j)/(nc_d(k,i,j)+nc_min) ) )
       xr_d(k,i,j) = min(xr_max, max(xr_min, lr_d(k,i,j)/(nr_d(k,i,j)+nr_min) ) )
       ! effective cross section is assume as area equivalent circle
       dr_xa_d(k,i,j)  = 2.0_RP*a_rea(I_QR)*xr_d(k,i,j)**b_rea(I_QR)
    enddo
    enddo
    enddo
    ! Auto-conversion, Accretion, Self-collection, Break-up
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
            rho_d, tem_d          )
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
    !
    ! update
    ! rhogq = l*gsgam
    !
    do k=KS, KE
    do j=1, JA
    do i=1, IA
       !
       ! warm collection process
       wrm_dqc = max( wdt*( PLCaut_d(k,i,j)+PLCacc_d(k,i,j) ), -lc_d(k,i,j)  )
       wrm_dnc = max( wdt*( PNCaut_d(k,i,j)+PNCacc_d(k,i,j) ), -nc_d(k,i,j)  )
       wrm_dnr = max( wdt*( PNRaut_d(k,i,j)+PNRslc_d(k,i,j)+PNRbrk_d(k,i,j) ), -nr_d(k,i,j) )
       !--- update
       !
       rhogq_d(k,i,j,I_QC) = max(0.0_RP, rhogq_d(k,i,j,I_QC) + wrm_dqc )
       rhogq_d(k,i,j,I_NC) = max(0.0_RP, rhogq_d(k,i,j,I_NC) + wrm_dnc )
       rhogq_d(k,i,j,I_QR) = max(0.0_RP, rhogq_d(k,i,j,I_QR) - wrm_dqc )
       rhogq_d(k,i,j,I_NR) = max(0.0_RP, rhogq_d(k,i,j,I_NR) + wrm_dnr )
    end do
    end do
    end do
    !--- update mixing ratio
    do iq = 1,  QA
    do k  = KS, KE
    do j  = 1,  JA
    do i  = 1,  IA
       q_d(k,i,j,iq) = rhogq_d(k,i,j,iq) * rrhog_d(k,i,j)
    enddo
    enddo
    enddo
    enddo

    call thrmdyn_qd( qd_d, q_d )
    call thrmdyn_cv( cva_d, q_d, qd_d )

    do k  = KS, KE
    do j = 1,    JA
    do i = 1,    IA
       tem_d(k,i,j) = rhoge_d(k,i,j) / ( rhog_d(k,i,j) * cva_d(k,i,j) )
       pre_d(k,i,j) = rho_d(k,i,j)*( qd_d(k,i,j)*Rdry+q_d(k,i,j,I_QV)*Rvap )*tem_d(k,i,j)
    enddo
    enddo
    enddo

!    if( opt_debug )     call debugreport_collection
    if( opt_debug_tem ) call debug_tem_kij( 4, tem_d(:,:,:), rho_d(:,:,:), pre_d(:,:,:), q_d(:,:,:,I_QV) )

    call TIME_rapend  ('MP3 Collection')

    call TIME_rapstart('MPX ijkconvert')

    do k  = KS, KE
    do j = 1,    JA
    do i = 1,    IA
       th_d(k,i,j) = tem_d(k,i,j) * ( P00 / pre_d(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          RHOT(k,i,j) = th_d(k,i,j) * rhog_d(k,i,j)
       enddo
!OCL XFILL
       do k = KS, KE
          rhoe(k,i,j) = rhoge_d(k,i,j)
       enddo
!OCL XFILL
       do k = KS, KE
          temp(k,i,j) = tem_d(k,i,j)
       enddo
!OCL XFILL
       do k = KS, KE
          pres(k,i,j) = pre_d(k,i,j)
       enddo
    enddo
    enddo

    do iq = 1, QA
       do j  = JS, JE
       do i  = IS, IE
!OCL XFILL
          do k = KS, KE
             rhoq(k,i,j,iq) = rhogq_d(k,i,j,iq)
          enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,iq) = rhoq(k,i,j,iq) / DENS(k,i,j)
       enddo
       enddo
       enddo
    enddo

    call TIME_rapend  ('MPX ijkconvert')

    !----------------------------------------------------------------------------
    !
    ! 4.Saturation adjustment
    !
    !----------------------------------------------------------------------------
    call TIME_rapstart('MP4 Saturation adjustment')
    ! nothing to do
    call TIME_rapend  ('MP4 Saturation adjustment')
    !----------------------------------------------------------------------------
    !
    ! 5. Sedimentation ( terminal velocity must be negative )
    !
    !----------------------------------------------------------------------------
    call TIME_rapstart('MP5 Sedimentation')

    if ( doprecipitation ) then

!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       flux_rain(k,i,j) = 0.0_RP
       flux_snow(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo

    do step = 1, MP_NSTEP_SEDIMENTATION

       call MP_terminal_velocity( velw(:,:,:,:), &
                                  rhoq(:,:,:,:), &
                                  temp(:,:,:),   &
                                  pres(:,:,:)    )

       call precipitation( wflux_rain(:,:,:),     &
                           wflux_snow(:,:,:),     &
                           velw(:,:,:,:),         &
                           rhoq(:,:,:,:),         &
                           rhoe(:,:,:),           &
                           temp(:,:,:),           &
                           MP_DTSEC_SEDIMENTATION )

       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          flux_rain(k,i,j) = flux_rain(k,i,j) + wflux_rain(k,i,j) * MP_RNSTEP_SEDIMENTATION
       enddo
       enddo
       enddo

!       if( opt_debug ) call debugreport_sedimentation

    enddo

    endif

    call TIME_rapend  ('MP5 Sedimentation')

    return
  end subroutine mp_sn13w

  !-----------------------------------------------------------------------------
  subroutine debug_tem_kij( &
      point,    &
      tem,      &
      rho,      &
      pre,      &
      qv        )
    use mod_process, only: &
       PRC_myrank
    implicit none

    integer, intent(in) :: point
    real(RP), intent(in) :: tem(KA,IA,JA)
    real(RP), intent(in) :: rho(KA,IA,JA)
    real(RP), intent(in) :: pre(KA,IA,JA)
    real(RP), intent(in) :: qv (KA,IA,JA)

    integer :: k ,i, j
    !---------------------------------------------------------------------------

    do k  = 1,  KA
    do j = JS, JE
    do i = IS, IE
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
  !
  subroutine nucleation_kij( &
       IA, JA, KA,      &
       IS, IE,       &
       JS, JE,       &
       KS, KE,       &
       z, w,             &
       rho, tem, pre,    &
       LV, NC,       &
       PNCccn, PLCccn,   &
       cva, cpa,         & ! in
       dTdt_rad,         & ! in
       qke,              & ! in
       flag_history_in,  & ! in
       ct, dt            ) ! in
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_const, only: &
       GRAV   => CONST_GRAV, &
       UNDEF8 => CONST_UNDEF8
    use mod_atmos_saturation, only: &
       moist_psat_water     => ATMOS_SATURATION_psat_liq,      &
       moist_psat_ice       => ATMOS_SATURATION_psat_ice,      &
       moist_qsat_water     => ATMOS_SATURATION_pres2qsat_liq, &
       moist_qsat_ice       => ATMOS_SATURATION_pres2qsat_ice, &
       moist_dqsw_dtem_rho  => ATMOS_SATURATION_dqsw_dtem_rho, &
       moist_dqsi_dtem_rho  => ATMOS_SATURATION_dqsi_dtem_rho, &
       moist_dqsw_dtem_dpre => ATMOS_SATURATION_dqsw_dtem_dpre
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
    real(RP), intent(out) :: PNCccn(KA,IA,JA) !
    real(RP), intent(out) :: PLCccn(KA,IA,JA) !
    !
    real(RP), intent(in) ::  cva(KA,IA,JA)      ! in  09/08/18 [Add] T.Mitsui
    real(RP), intent(in) ::  cpa(KA,IA,JA)      ! in  09/08/18 [Add] T.Mitsui
    real(RP), intent(in)  :: dTdt_rad(KA,IA,JA) ! 09/08/18 T.Mitsui
    real(RP), intent(in)  :: qke(KA,IA,JA)      ! 09/08/18 T.Mitsui
    real(RP), intent(in)  :: ct ! current time[sec.],   09/04/14x T.Mitsui
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
    real(RP), save :: ssw_max= 1.1_RP  ! [%]
    !
    logical, save :: flag_first = .true.
    real(RP), save :: qke_min = 0.03_RP ! sigma=0.1[m/s], 09/08/18 T.Mitsui
    real(RP), save :: tem_ccn_low=233.150_RP  ! = -40 degC  ! [Add] 10/08/03 T.Mitsui
    !
    namelist /nm_mp_sn13w_nucleation/ &
         c_ccn, kappa,               & ! cloud nucleation
         xc_ccn,                     &
         tem_ccn_low,                & ! [Add] 10/08/03 T.Mitsui
         ssw_max
    !
    real(RP) :: c_ccn_map(1,IA,JA)   ! c_ccn horizontal distribution
    real(RP) :: kappa_map(1,IA,JA)   ! kappa horizontal distribution
    real(RP) :: esw(KA,IA,JA)      ! saturation vapor pressure, water
    real(RP) :: ssw(KA,IA,JA)      ! super saturation (water)
    real(RP) :: ssw_below(KA,IA,JA)! ssw(k-1)
    real(RP) :: z_below(KA,IA,JA)  ! z(k-1)
    real(RP) :: dz                   ! z(k)-z(k-1)
    real(RP) :: pv                   ! vapor pressure
    real(RP) :: n_in                 ! number of ice nuclei
    ! work variables for Twomey Equation.
    real(RP) :: qsw(KA,IA,JA)
    real(RP) :: r_qsw(KA,IA,JA)
    !
    real(RP) :: alpha_nuc(1,IA,JA) ! alpha_nuc
    real(RP) :: eta_nuc(1,IA,JA)   ! xi use the value @ cloud base
    real(RP) :: Dw, Q1, Q2
    !
    real(RP) :: sigma_w(KA,IA,JA)
    real(RP) :: weff(KA,IA,JA)
    real(RP) :: weff_max(KA,IA,JA)
    real(RP) :: w_m(KA,IA,JA)
    !
    real(RP) :: coef_ccn(IA,JA)
    real(RP) :: slope_ccn(IA,JA)
    real(RP) :: nc_new(KA,IA,JA)
    real(RP) :: nc_new_below(KA,IA,JA)
    real(RP) :: dnc_new
    real(RP) :: nc_new_max   ! Lohmann (2002),
    real(RP) :: a_max
    real(RP) :: b_max
    logical  :: flag_nucleation(KA,IA,JA)
    !
    real(RP) :: r_gravity
    real(RP), parameter :: r_sqrt3=0.577350269_RP ! = sqrt(1.d0/3.d0)
    real(RP), parameter :: eps=1.E-30_RP
    !====> ! 09/08/18
    !
    real(RP) :: dlcdt_max ! defined by supersaturation
    real(RP) :: dncdt_max ! defined by supersaturation
    real(RP) :: rdt
    real(RP) :: work
    !
    integer :: i, j, k, kk
    !
    if( flag_first )then
       rewind(IO_FID_CONF)
       read(IO_FID_CONF, nml=nm_mp_sn13w_nucleation, end=100)
100    if( IO_L ) write(IO_FID_LOG, nml=nm_mp_sn13w_nucleation)
       flag_first=.false.
    endif
    !
    c_ccn_map(1,:,:) = c_ccn
    kappa_map(1,:,:) = kappa
    !
    nc_uplim_d(1,:,:)  = c_ccn_map(1,:,:)*1.5_RP
    !
    rdt            = 1.0_RP/dt
    r_gravity      = 1.0_RP/GRAV
    PNCccn(:,:,:)    = 0.0_RP
    PLCccn(:,:,:)    = 0.0_RP
    ssw(:,:,:)       = 0.0_RP
    ssw_below(:,:,:) = 0.0_RP
    z_below(:,:,:)   = 0.0_RP
    weff(:,:,:)      = 0.0_RP
    work           = r_sqrt3*sqrt(qke_min)
    sigma_w(:,:,:)   = work
    nc_new(:,:,:)      = 0.0_RP
    nc_new_below(:,:,:)= 0.0_RP
    !
    call moist_psat_water( esw, tem )
    call moist_qsat_water( qsw, tem, pre )
    r_qsw(:,:,:) = 1.0_RP/qsw
    !
    do k=KS, KE
       do j=1, JA
       do i=1, IA
          pv        = LV(k,i,j)*Rvap*tem(k,i,j)
          ssw(k,i,j) = min( MP_ssw_lim, ( pv/esw(k,i,j)-1.0_RP ) )*100.0_RP
          ssw_below(k+1,i,j) = ssw(k,i,j)
          z_below(k+1,i,j)   = z(k)
       end do
       end do
    end do
    ssw_below(KS,:,:) = ssw(KS,:,:)
    z_below(KS,:,:)   = z(KS-1)
    !
    !
    ! dS/dz is evaluated by first order upstream difference
    !***  Solution for Twomey Equation ***
    do j=1, JA
    do i=1, IA
       coef_ccn(i,j)  = 1.E+6_RP*0.88_RP*(c_ccn_map(1,i,j)*1.E-6_RP)**(2.0_RP/(kappa_map(1,i,j) + 2.0_RP)) * &
            (70.0_RP)**(kappa_map(1,i,j)/(kappa_map(1,i,j) + 2.0_RP))
       slope_ccn(i,j) = 1.5_RP*kappa_map(1,i,j)/(kappa_map(1,i,j) + 2.0_RP)
    end do
    end do
    !
    do k=KS, KE
       sigma_w(k,:,:) = r_sqrt3*sqrt(max(qke(k,:,:),qke_min))
    end do
    sigma_w(KS-1,:,:) = sigma_w(KS,:,:)
    sigma_w(KE+1,:,:) = sigma_w(KE,:,:)
    ! effective vertical velocity
    do k=KS, KE-1
       do j=1, JA
       do i=1, IA
          weff(k,i,j) = 0.5_RP*(w(k,i,j) + w(k+1,i,j)) - cpa(k,i,j)*r_gravity*dTdt_rad(k,i,j)
       end do
       end do
    end do
    weff(KS-1,:,:) = weff(KS,:,:)
    weff(KE,:,:)   = weff(KE-1,:,:)
    ! Lohmann (2002),JAS, eq.(1) but changing unit [cm-3] => [m-3]
    a_max = 1.E+6_RP*0.1_RP*(1.E-6_RP)**1.27_RP
    b_max = 1.27_RP
    ! diagnose cloud condensation nuclei
    do k=KS, KE
       do j=1, JA
       do i=1, IA
          ! effective vertical velocity (maximum vertical velocity in turbulent flow)
          weff_max(k,i,j) = weff(k,i,j) + sigma_w(k,i,j)
          ! large scale upward motion region and saturated
          if( (weff(k,i,j) > 1.E-8_RP) .AND. (ssw(k,i,j) > 1.E-10_RP)  .AND. pre(k,i,j) > 300.E+2_RP )then
             ! Lohmann (2002), eq.(1)
             nc_new_max   = coef_ccn(i,j)*weff_max(k,i,j)**slope_ccn(i,j)
             nc_new(k,i,j) = a_max*nc_new_max**b_max
          endif
       end do
       end do
    end do
    !
    flag_nucleation(:,:,:)=.false.
    do k=KS, KE
       do j=1, JA
       do i=1, IA
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
          endif
       end do
       end do
    end do
    nc_new_below(KS,:,:) = 0.0_RP
    ! search maximum value of nc_new
    do k=KS, KE
       do j=1, JA
       do i=1, IA
          if(  ( nc_new(k,i,j) < nc_new_below(k,i,j) ) .or. &
               ( nc_new_below(k,i,j) > c_ccn_map(1,i,j)*0.05_RP ) )then ! 5% of c_ccn
             flag_nucleation(k,i,j) = .false.
          endif
       end do
       end do
    end do
    ! nucleation occurs at only cloud base.
    ! if CCN is more than below parcel, nucleation newly occurs
    do k=KS, KE
       do j=1, JA
       do i=1, IA
          ! effective vertical velocity
          if(   flag_nucleation(k,i,j)               .AND. & ! large scale upward motion region and saturated
               ( tem(k,i,j)    > tem_ccn_low       ) .AND. &
               ( nc_new(k,i,j) > NC(k,i,j) )                )then
             dlcdt_max    = (LV(k,i,j) - esw(k,i,j)/(Rvap*tem(k,i,j)))*rdt
             dncdt_max    = dlcdt_max/xc_min
             dnc_new      = nc_new(k,i,j)-NC(k,i,j)
             PNCccn(k,i,j) = min( dncdt_max, dnc_new*rdt )
             PLCccn(k,i,j) = min( dlcdt_max, xc_min*PNCccn(k,i,j) )
          endif
       end do
       end do
    end do
    !
    return
  end subroutine nucleation_kij
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
       rho, tem                )
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
    real(RP), intent(in)  :: tem(KA,IA,JA)
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
    PLCaut(1:KS,:,:)=0.0_RP
    PNCaut(1:KS,:,:)=0.0_RP
    PNRaut(1:KS,:,:)=0.0_RP
    PLCacc(1:KS,:,:)=0.0_RP
    PNCacc(1:KS,:,:)=0.0_RP
    PNRslc(1:KS,:,:)=0.0_RP
    PNRbrk(1:KS,:,:)=0.0_RP
    !
    PLCaut(KE:KA,:,:)=0.0_RP
    PNCaut(KE:KA,:,:)=0.0_RP
    PNRaut(KE:KA,:,:)=0.0_RP
    PLCacc(KE:KA,:,:)=0.0_RP
    PNCacc(KE:KA,:,:)=0.0_RP
    PNRslc(KE:KA,:,:)=0.0_RP
    PNRbrk(KE:KA,:,:)=0.0_RP
    !
    coef_nuc0 = (nu(I_QC)+2.0_RP)/(nu(I_QC)+1.0_RP)
    coef_nuc1 = (nu(I_QC)+2.0_RP)*(nu(I_QC)+4.0_RP)/(nu(I_QC)+1.0_RP)/(nu(I_QC)+1.0_RP)
    coef_aut0 =  -kcc*coef_nuc0
    coef_aut1 =  -kcc/x_sep/20._RP*coef_nuc1
    !
    do k=KS, KE
       do j=1, JA
       do i=1, IA
          lwc = LR(k,i,j)+LC(k,i,j)
          if( lwc > xc_min )then
             tau  = max(tau_min, LR(k,i,j)/lwc)
          else
             tau  = tau_min
          endif
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
          endif
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
       rho, tem, pre, LV,    & ! in
       qd,                   & ! in
!!!!!!!! mod xxxxxxxxxxxxxxxxxxxx T.Seiki             
!!$    esw, esi,             & ! in
       esw,                  & ! in
       NC, NR,               & ! in
       xc, xr,               & ! in
       vt_xave,              & ! in
       dc_xave, dr_xave      ) ! in
    implicit none

    integer, intent(in)  :: IA,JA,KA
    integer, intent(in)  :: IS,JS,KS
    integer, intent(in)  :: IE,JE,KE
    ! Diffusion growth or Evaporation, Sublimation
    real(RP), intent(out) :: PLCdep(KA,IA,JA)  ! mass change   for cloud, [Add]  09/08/18 T.Mitsui
    real(RP), intent(out) :: PLRdep(KA,IA,JA)  ! mass change   for rain deposion
    real(RP), intent(out) :: PNRdep(KA,IA,JA)  ! number change
 
    real(RP), intent(in)  :: rho(KA,IA,JA)     ! air density
    real(RP), intent(in)  :: tem(KA,IA,JA)     ! air temperature
    real(RP), intent(in)  :: pre(KA,IA,JA)     ! air pressure
    real(RP), intent(in)  :: qd(KA,IA,JA)      ! mixing ratio of dry air
    real(RP), intent(in)  :: esw(KA,IA,JA)     ! saturation vapor pressure(liquid water)
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$ real(RP), intent(in)  :: esi(KA,IA,JA)     ! saturation vapor pressure(solid water)
    real(RP), intent(in)  :: LV(KA,IA,JA)      ! mass   of vapor
    real(RP), intent(in)  :: NC(KA,IA,JA)      ! number of cloud  09/08/18 [Add] T.Mitsui
    real(RP), intent(in)  :: NR(KA,IA,JA)      ! number of rain
    real(RP), intent(in)  :: xc(KA,IA,JA)      ! mean mass of cloud(filtered) [Add] 09/08/18 T.Mitsui
    real(RP), intent(in)  :: xr(KA,IA,JA)      ! mean mass of rain(filtered)
    ! Notice following values differ from mean terminal velocity or diameter.
    ! mean(vt(x)) /= vt(mean(x)) and mean(D(x)) /= D(mean(x))
    ! Following ones are vt(mean(x)) and D(mean(x)).
    real(RP), intent(in)  :: vt_xave(KA,IA,JA,HYDRO_MAX,2) ! terminal velocity of mean cloud 09/08/18 [Add], T.Mitsui
    !
    real(RP), intent(in)  :: dc_xave(KA,IA,JA) ! diameter of mean cloud 09/08/18 [Add] T.Mitsui
    real(RP), intent(in)  :: dr_xave(KA,IA,JA) ! diameter of mean rain
    !
    real(RP) :: rho_lim            ! limited density              09/08/18 T.Mitsui
    real(RP) :: temc_lim           ! limited temperature[celsius] 09/08/18 T.Mitsui
    real(RP) :: pre_lim            ! limited density              09/08/18 T.Mitsui
    real(RP) :: temc               ! temperature[celsius]
    real(RP) :: pv                 ! vapor pressure
    real(RP) :: qv                 ! mixing ratio of water vapor [Add] 09/08/18
    real(RP) :: ssw                ! super saturation ratio(liquid water)
    real(RP) :: nua, r_nua         ! kinematic viscosity of air
    real(RP) :: mua                ! viscosity of air
    real(RP) :: Kalfa              ! thermal conductance
    real(RP) :: Dw                 ! diffusivity of water vapor
    real(RP) :: Dt                 ! diffusivity of heat
    real(RP) :: Gw                 ! diffusion factor by balance between heat and vapor
    real(RP) :: Gwr                ! for rain, ice, snow and graupel.
    real(RP) :: Nsc_r3             !
    ! [Mod] 11/08/30 T.Mitsui, considering large and small branches
    real(RP) :: Nrecs_r2            ! 09/08/18 [Add] T.Mitsui
    real(RP) :: Nrers_r2, Nreis_r2  !
    real(RP) :: Nress_r2, Nregs_r2  !
    real(RP) :: Nrecl_r2            ! 09/08/18 [Add] T.Mitsui
    real(RP) :: Nrerl_r2, Nreil_r2  !
    real(RP) :: Nresl_r2, Nregl_r2  !
    real(RP) :: NscNrer_s, NscNrer_l
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    real(RP) :: NscNrei_s, NscNrei_l
!!$    real(RP) :: NscNres_s, NscNres_l
!!$    real(RP) :: NscNreg_s, NscNreg_l
    real(RP) :: ventLR_s, ventLR_l
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    real(RP) :: ventNI_s, ventNI_l, ventLI_s, ventLI_l
!!$    real(RP) :: ventNS_s, ventNS_l, ventLS_s, ventLS_l
!!$    real(RP) :: ventNG_s, ventNG_l, ventLG_s, ventLG_l
    !
    real(RP) :: wtr, wti, wts, wtg
    real(RP), parameter :: r_14=1.0_RP/1.4_RP
    real(RP), parameter :: r_15=1.0_RP/1.5_RP
    !
    real(RP) :: ventNR, ventLR(KA,IA,JA)     !
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    real(RP) :: ventNI(KA,IA,JA), ventLI(KA,IA,JA)     !
!!$    real(RP) :: ventNS(KA,IA,JA), ventLS(KA,IA,JA)     !
!!$    real(RP) :: ventNG(KA,IA,JA), ventLG(KA,IA,JA)     !
    real(RP) :: ah_vent1_rs(KA,IA,JA)
    real(RP) :: ah_vent1_rl(KA,IA,JA)
    real(RP) :: bh_vent1_rs(KA,IA,JA)
    real(RP) :: bh_vent1_rl(KA,IA,JA)
    !
    real(RP), parameter :: Re_max=1.E+3_RP
    real(RP), parameter :: Re_min=1.E-4_RP
    !
    integer :: i, j, k
    !
    PLCdep(1:KS,:,:)=0.0_RP
    PLRdep(1:KS,:,:)=0.0_RP
    PNRdep(1:KS,:,:)=0.0_RP
    !
    PLCdep(KE:KA,:,:)=0.0_RP
    PLRdep(KE:KA,:,:)=0.0_RP
    PNRdep(KE:KA,:,:)=0.0_RP
    !
    ah_vent1_rs(:,:,:)  = ah_vent1(I_QR,1)
    ah_vent1_rl(:,:,:)  = ah_vent1(I_QR,2)
    bh_vent1_rs(:,:,:)  = bh_vent1(I_QR,1)
    bh_vent1_rl(:,:,:)  = bh_vent1(I_QR,2)
    !
    ventNR=0.0_RP
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    ventNI=0.0_RP
!!$    ventNS=0.0_RP
!!$    ventNG=0.0_RP
    ventLR=0.0_RP
!!!!!!!! delete xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    ventLI=0.0_RP
!!$    ventLS=0.0_RP
!!$    ventLG=0.0_RP
    !
    ! Notice,T.Mitsui
    ! Vapor deposition and melting would not be solved iteratively to reach equilibrium.
    ! Because following phenomena are not adjustment but transition.
    ! Just time-scales differ among them.
    ! If we would treat more appropreately, there would be time-splitting method to solve each ones.
    do k=KS, KE
       do j=1, JA
       do i=1, IA
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
          ! capacities account for their surface geometries
          Gwr     = 4.0_RP*PI/cap(I_QR)/Gw
          ! vent: ventilation effect( asymmetry vapor field around particles due to aerodynamic )
          ! SB06 (30),(31) and each coefficient is by (88),(89)
          Nsc_r3  = (nua/Dw)**(0.33333333_RP)                    ! (Schmidt number )^(1/3)
          !
          Nrecs_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QC,1)*dc_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud
          Nrers_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QR,1)*dr_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) rain
          !
          Nrecl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QC,2)*dc_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) cloud
          Nrerl_r2 = sqrt(max(Re_min,min(Re_max,vt_xave(k,i,j,I_QR,2)*dr_xave(k,i,j)*r_nua))) ! (Reynolds number)^(1/2) rain
          NscNrer_s=Nsc_r3*Nrers_r2 ! small rain
          NscNrer_l=Nsc_r3*Nrerl_r2 ! large rain
          !
          ventLR_s = ah_vent1(I_QR,1) + bh_vent1(I_QR,1)*NscNrer_s
          ventLR_l = ah_vent1(I_QR,2) + bh_vent1(I_QR,2)*NscNrer_l
          !
          ! branch is 1.4 for rain, snow, graupel; is 1.0 for ice (PK97, 13-60,-61,-88,-89).
          !
          wtr     = ( min(max( NscNrer_s*r_14, 0.5_RP), 2.0_RP) -0.5_RP )*r_15 ! weighting between 1.4*0.5 and 1.4*2
          !
          ventLR(k,i,j)  = (1.0_RP-wtr)*ventLR_s + wtr*ventLR_l
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
          PLCdep(k,i,j) = Gwr*NC(k,i,j)*dc_xave(k,i,j)*coef_deplc
          PLRdep(k,i,j) = Gwr*NR(k,i,j)*dr_xave(k,i,j)*ventLR(k,i,j)
          PNRdep(k,i,j) = PLRdep(k,i,j)/xr(k,i,j)
          !
          !------------------------------------------------------------------------
          ! Melting part is described by Pruppacher and Klett (1997) Sec.16.3.1
          ! Here we omit "Shedding" of snow-flakes and ice-particles.
          ! "Shedding" may be applicative if you refer
          ! eq.(38) in Cotton etal.(1986) Jour. Clim. Appl. Meteor. p.1658-1680.
          ! SB06(73)
          Dt      = Kalfa/(CPvap*rho_0)
       end do
       end do
    end do
    !
    return
  end subroutine dep_vapor_melt_ice_kij
  !-----------------------------------------------------------------------------
  function gammafunc( xx ) result(f)
    implicit none
    real(RP), intent(in) :: xx
    real(RP) :: f
    real(RP) :: coef(6)=(/&
         +76.18009172947146_RP,&
         -86.50532032941677_RP,&
         +24.01409824083091_RP,&
         -1.231739572450155_RP,&
         +0.1208650973866179E-2_RP,&
         -0.5395239384953E-5_RP&
         /)
    integer :: j
    real(RP) :: x,y,tmp,ser

    x=xx
    y=x
    tmp=x+5.5_RP
    tmp = tmp - (x+0.5_RP)*log(tmp)
    ser=1.000000000190015_RP
    do j=1,6
       y=y+1
       ser = ser+coef(j)/y
    end do
    f = exp(-tmp+log(2.5066282746310005_RP*ser/x))
  end function gammafunc
  !-----------------------------------------------------------------------------
  function gammafunc_3d( IJA, KA, x )
    implicit none
    integer, intent(in) :: IJA
    integer, intent(in) :: KA
    real(RP), intent(in) :: x(IJA,KA)
    real(RP)             :: gammafunc_3d(IJA,KA)
    real(RP), parameter  :: coef(6)=(/&
         +76.18009172947146_RP,&
         -86.50532032941677_RP,&
         +24.01409824083091_RP,&
         -1.231739572450155_RP,&
         +0.1208650973866179E-2_RP,&
         -0.5395239384953E-5_RP&
         /)
    real(RP), parameter :: ser0=1.000000000190015_RP
    real(RP) :: tmp(IJA,KA)
    real(RP) :: ser(IJA,KA)
    integer :: ij,k,iter
    !
    ser(:,:) = ser0 &
         + coef(1)/(x(:,:)+1.0_RP) + coef(2)/(x(:,:)+2.0_RP) + coef(3)/(x(:,:)+3.0_RP) &
         + coef(4)/(x(:,:)+4.0_RP) + coef(5)/(x(:,:)+5.0_RP) + coef(6)/(x(:,:)+6.0_RP)
    tmp(:,:) = x(:,:)+5.5_RP - (x(:,:)+0.5_RP)*log(x(:,:)+5.5_RP)
    gammafunc_3d(:,:)   = exp(-tmp(:,:)+log(2.5066282746310005_RP*ser(:,:)/x(:,:)))
    return
  end function gammafunc_3d
  !-----------------------------------------------------------------------------
  function igammafunc_3d( IJA, KA, x, alpha, gm ) result(igm)
    ! Section 6.2 in "Numerical Recipes in C"
    ! function is represented by(6.2.1)
    ! g(x,alpha)=1/gamma(alpha)*integral_0^x { exp(-t) * t^(alpha-1) }dt
    ! g(0)=0 and g(+infinity)=1
    implicit none
    !
    integer, intent(in) :: IJA
    integer, intent(in) :: KA
    real(RP), intent(in) :: x(IJA,KA)     !
    real(RP), intent(in) :: alpha(IJA,KA) !
    real(RP), intent(in) :: gm(IJA,KA)    ! gamma function
    ! incomplete gamma function
    real(RP) :: igm(IJA,KA)
    !
    real(RP) :: lx(IJA,KA)
    real(RP) :: lgm(IJA,KA)
    ! work
    ! coefficient of expansion using in calculation of igm
    real(RP) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15
    real(RP) :: an1,an2,an3,an4,an5
    real(RP) :: b0,b1,b2,b3,b4,b5
    real(RP) :: c0,c1,c2,c3,c4,c5
    real(RP) :: d0,d1,d2,d3,d4,d5
    real(RP) :: e0,e1,e2,e3,e4,e5
    real(RP) :: h0,h1,h2,h3,h4,h5
    real(RP), parameter :: eps=1.E-30_RP
    !
    integer :: ij,k
    !
    lgm(:,:)=log(gm(:,:))
    lx(:,:) =log( x(:,:))
    !
    ! Incomplete Gamma Function
    !
    do k=1,KA
       do ij=1,IJA
          if     ( x(ij,k) < 1.E-2_RP*alpha(ij,k) )then ! negligible
             igm(ij,k)=0.0_RP
          else if( x(ij,k) < alpha(ij,k)+1.0_RP )then ! series expansion (6.2.5)
             !
             ! 10th-truncation is enough for cloud droplet.
             a0   = 1.0_RP/alpha(ij,k)         ! n=0
             a1   = a0*x(ij,k)/(alpha(ij,k)+1.0_RP)  ! n=1
             a2   = a1*x(ij,k)/(alpha(ij,k)+2.0_RP)  ! n=2
             a3   = a2*x(ij,k)/(alpha(ij,k)+3.0_RP)  ! n=3
             a4   = a3*x(ij,k)/(alpha(ij,k)+4.0_RP)  ! n=4
             a5   = a4*x(ij,k)/(alpha(ij,k)+5.0_RP)  ! n=5
             !
             a6   = a5*x(ij,k)/(alpha(ij,k)+6.0_RP)  ! n=6
             a7   = a6*x(ij,k)/(alpha(ij,k)+7.0_RP)  ! n=7
             a8   = a7*x(ij,k)/(alpha(ij,k)+8.0_RP)  ! n=8
             a9   = a8*x(ij,k)/(alpha(ij,k)+9.0_RP)  ! n=9
             a10  = a9*x(ij,k)/(alpha(ij,k)+10.0_RP) ! n=10
             !
             a11  = a10*x(ij,k)/(alpha(ij,k)+11.0_RP) ! n=11
             a12  = a11*x(ij,k)/(alpha(ij,k)+12.0_RP) ! n=12
             a13  = a12*x(ij,k)/(alpha(ij,k)+13.0_RP) ! n=13
             a14  = a13*x(ij,k)/(alpha(ij,k)+14.0_RP) ! n=14
             a15  = a14*x(ij,k)/(alpha(ij,k)+15.0_RP) ! n=15
             !
             igm(ij,k) = (a0+a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13+a14+a15)&
                  * exp( -x(ij,k) + alpha(ij,k)*lx(ij,k) - lgm(ij,k) )
          else if( x(ij,k) < alpha(ij,k)*1.E+2_RP ) then ! continued fraction expansion (6.2.6)
             !
             ! 3rd-truncation is enough for rain droplet.
             ! setup
             b0   = x(ij,k)+1.0_RP-alpha(ij,k)
             c0   = 1.0_RP/eps
             d0   = 1.0_RP/b0
             h0   = d0
             ! n=1
             an1  = -     (1.0_RP-alpha(ij,k))
             b1   = b0 + 2.0_RP
             d1   = 1.0_RP/(an1*d0+b1)
             c1   = b1+an1/c0
             e1   = d1*c1
             h1   = h0*e1
             ! n=2
             an2  = -2.0_RP*(2.0_RP-alpha(ij,k))
             b2   = b1 + 2.0_RP
             d2   = 1.0_RP/(an2*d1+b2)
             c2   = b2+an2/c1
             e2   = d2*c2
             h2   = h1*e2
             ! n=3
             an3  = -3.0_RP*(3.0_RP-alpha(ij,k))
             b3   = b2 + 2.0_RP
             d3   = 1.0_RP/(an3*d2+b3)
             c3   = b3+an3/c2
             e3   = d3*c3
             h3   = h2*e3
             ! n=4
             an4  = -4.0_RP*(4.0_RP-alpha(ij,k))
             b4   = b3 + 2.0_RP
             d4   = 1.0_RP/(an4*d3+b4)
             c4   = b4+an4/c3
             e4   = d4*c4
             h4   = h3*e4
             ! n=5
             an5  = -5.0_RP*(5.0_RP-alpha(ij,k))
             b5   = b4 + 2.0_RP
             d5   = 1.0_RP/(an5*d4+b5)
             c5   = b5+an5/c4
             e5   = d5*c5
             h5   = h4*e5
             !
             igm(ij,k)  = 1.0_RP - exp( -x(ij,k) + alpha(ij,k)*lx(ij,k) - lgm(ij,k) )*h5
          else   ! negligible
             igm(ij,k)  = 1.0_RP
          endif
       end do
    end do
    !
    return
  end function igammafunc_3d
  !-------------------------------------------------------------------------------
  function betafunc_3d( IJA, KA, x, w )
    implicit none
    integer, intent(in) :: IJA
    integer, intent(in) :: KA
    real(RP), intent(in) :: x(IJA,KA)
    real(RP), intent(in) :: w(IJA,KA)
    real(RP)             :: betafunc_3d(IJA,KA)
    real(RP), parameter  :: coef(6)=(/&
         +76.18009172947146_RP,&
         -86.50532032941677_RP,&
         +24.01409824083091_RP,&
         -1.231739572450155_RP,&
         +0.1208650973866179E-2_RP,&
         -0.5395239384953E-5_RP&
         /)
    real(RP), parameter :: ser0=1.000000000190015_RP
    real(RP) :: y(IJA,KA)
    real(RP) :: tmp_x(IJA,KA)
    real(RP) :: tmp_w(IJA,KA)
    real(RP) :: tmp_xw(IJA,KA)
    real(RP) :: ser_x(IJA,KA)
    real(RP) :: ser_w(IJA,KA)
    real(RP) :: ser_xw(IJA,KA)
    real(RP) :: lg_x(IJA,KA)
    real(RP) :: lg_w(IJA,KA)
    real(RP) :: lg_xw(IJA,KA)
    integer :: ij,k,iter
    !
    ! log(gamma(x))
    !
    ser_x(:,:) = ser0 &
         + coef(1)/(x(:,:)+1.0_RP) + coef(2)/(x(:,:)+2.0_RP) + coef(3)/(x(:,:)+3.0_RP) &
         + coef(4)/(x(:,:)+4.0_RP) + coef(5)/(x(:,:)+5.0_RP) + coef(6)/(x(:,:)+6.0_RP)
    tmp_x(:,:) = x(:,:)+5.5_RP - (x(:,:)+0.50_RP)*log(x(:,:)+5.5_RP)
    lg_x(:,:)  = -tmp_x(:,:)+log(2.5066282746310005_RP*ser_x(:,:)/x(:,:))
    !
    ! log(gamma(w))
    !
    ser_w(:,:) = ser0 &
         + coef(1)/(w(:,:)+1.0_RP) + coef(2)/(w(:,:)+2.0_RP) + coef(3)/(w(:,:)+3.0_RP) &
         + coef(4)/(w(:,:)+4.0_RP) + coef(5)/(w(:,:)+5.0_RP) + coef(6)/(w(:,:)+6.0_RP)
    tmp_w(:,:) = w(:,:)+5.5_RP - (w(:,:)+0.5_RP)*log(w(:,:)+5.5_RP)
    lg_w(:,:)  = -tmp_w(:,:)+log(2.5066282746310005_RP*ser_w(:,:)/w(:,:))
    !
    ! log(gamma(x+w))
    !
    y(:,:) = x(:,:) + w(:,:)
    ser_xw(:,:) = ser0 &
         + coef(1)/(y(:,:)+1.0_RP) + coef(2)/(y(:,:)+2.0_RP) + coef(3)/(y(:,:)+3.0_RP) &
         + coef(4)/(y(:,:)+4.0_RP) + coef(5)/(y(:,:)+5.0_RP) + coef(6)/(y(:,:)+6.0_RP)
    tmp_xw(:,:) = y(:,:)+5.5_RP - (y(:,:)+0.50_RP)*log(y(:,:)+5.5_RP)
    lg_xw(:,:)  = -tmp_xw(:,:)+log(2.5066282746310005_RP*ser_xw(:,:)/y(:,:))
    !
    betafunc_3d(:,:) = exp( lg_x(:,:) + lg_w(:,:) - lg_xw(:,:))
    !
    return
  end function betafunc_3d

  !-----------------------------------------------------------------------------
  subroutine MP_terminal_velocity( &
      velw, &
      rhoq, &
      temp, &
      pres  )
    use mod_const, only: &
       PI => CONST_PI
    use mod_atmos_vars, only: &
       DENS
    implicit none

    real(RP), intent(out) :: velw(KA,IA,JA,QA) ! terminal velocity of cloud mass
    real(RP), intent(in)  :: rhoq(KA,IA,JA,QA) ! rho * q
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
       enddo

       ! limiter
       do k  = KS, KE
          xq(k,I_QC) = max( min( xq(k,I_QC), xqmax(I_QC) ), xqmin(I_QC) )
          xq(k,I_QR) = max( min( xq(k,I_QR), xqmax(I_QR) ), xqmin(I_QR) )
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

       ! small/large branch terminal velocity
       do k = KS, KE
          velq_s(k,I_QR) = coef_vtr_ar2 * dq(k,I_QR) * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar(k) )**(-5-mud_r) )
          velq_l(k,I_QR) = coef_vtr_ar1 - coef_vtr_br1 * ( 1.0_RP + coef_vtr_cr1*rlambdar(k) )**(-4-mud_r)
          velq_s(k,I_NR) = coef_vtr_ar2 * dq(k,I_QR) * ( 1.0_RP - ( 1.0_RP + coef_vtr_br2*rlambdar(k) )**(-2-mud_r) )
          velq_l(k,I_NR) = coef_vtr_ar1 - coef_vtr_br1 * ( 1.0_RP + coef_vtr_cr1*rlambdar(k) )**(-1-mud_r)
       enddo

       ! weighting coefficient for 2-branches is determined by ratio between 0.745mm and weighted diameter
       ! SB06 Table.1
       do k = KS, KE
          weight(k,I_QR) = 0.5_RP * ( 1.0_RP + tanh( PI * log( dq(k,I_QR)/d_vtr_branch ) ) ) ! Lr
          weight(k,I_NR) = 0.5_RP * ( 1.0_RP + tanh( PI * log( dq(k,I_NR)/d_vtr_branch ) ) ) ! Nr
       enddo

       do k = KS, KE
          rhofac_q(k,I_QC) = rhofac(k) ** gamma_v(I_QC)
          rhofac_q(k,I_QR) = rhofac(k) ** gamma_v(I_QR)
       enddo

       ! interpolated terminal velocity
       ! SB06(78) these are defined as negative value
       do k = KS, KE
          velw(k,i,j,I_QC) = -rhofac_q(k,I_QC) * coef_vt1(I_QC,1) * xq(k,I_QC)**beta_v(I_QC,1)
          velw(k,i,j,I_NC) = -rhofac_q(k,I_QC) * coef_vt0(I_QC,1) * xq(k,I_QC)**beta_v(I_QC,1)
       enddo
       do k = KS, KE
          velw(k,i,j,I_QR) = -rhofac_q(k,I_QR) * ( velq_l(k,I_QR) * (        weight(k,I_QR) ) &
                                                 + velq_s(k,I_QR) * ( 1.0_RP - weight(k,I_QR) ) )
          velw(k,i,j,I_NR) = -rhofac_q(k,I_QR) * ( velq_l(k,I_NR) * (        weight(k,I_NR) ) &
                                                 + velq_s(k,I_NR) * ( 1.0_RP - weight(k,I_NR) ) )
       enddo
    enddo
    enddo

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
       tem, pre, rho,            & ! inout
       cva,                      & ! out
!!!!!!!! mod xxxxxxxxxxxxxxxxxxxx T.Seiki 
!!$    esw, esi, LV,             & ! in
       esw, LV,             & ! in
       LC, LR,                   & ! in
       NC, NR,                   & ! in
       !+++ tendency terms
       ! condensation/evaporation, deposition/sublimation
       PLCdep,         &
       PLRdep, PNRdep, &
       !
       flag_history_in,&
       ! column integrated tendency terms
       sl_PLCdep,      &      !
       sl_PLRdep, sl_PNRdep ) !
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_atmos_thermodyn, only: &
       thrmdyn_qd      => ATMOS_THERMODYN_qd, &
       thrmdyn_cv      => ATMOS_THERMODYN_cv, &
       thrmdyn_cp      => ATMOS_THERMODYN_cp
    use mod_atmos_saturation, only: &
       moist_qsat_water     => ATMOS_SATURATION_qsat_liq,      &
       moist_dqsw_dtem_rho  => ATMOS_SATURATION_dqsw_dtem_rho, &
       moist_dqsw_dtem_dpre => ATMOS_SATURATION_dqsw_dtem_dpre
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
    real(RP), intent(in)    :: rho(KA,IA,JA)         ! air density[kg/m3]
    real(RP), intent(out)   :: cva(KA,IA,JA)         ! specific heat at constant volume
    real(RP), intent(in)    :: esw(KA,IA,JA)         ! saturated vapor pressure for liquid
    real(RP), intent(in)    :: lv(KA,IA,JA)          ! vapor mass [kg/m3]
    real(RP), intent(in)    :: lc(KA,IA,JA), nc(KA,IA,JA) ! cloud mass [kg/m3], number[/m3]
    real(RP), intent(in)    :: lr(KA,IA,JA), nr(KA,IA,JA) ! rain
    !+++ Condensation/Evaporation, Deposition/Sublimaion tendency[kg/m3/s]
    real(RP), intent(inout) :: PLCdep(KA,IA,JA)
    real(RP), intent(inout) :: PLRdep(KA,IA,JA), PNRdep(KA,IA,JA)
    !+++
    logical, intent(in)    :: flag_history_in
    !+++ Column integrated tendency[kg/m2/s]
    real(RP), intent(inout) :: sl_PLCdep(1,IA,JA)
    real(RP), intent(inout) :: sl_PLRdep(1,IA,JA), sl_PNRdep(1,IA,JA)
    !
    real(RP) :: xi(KA,IA,JA)                     ! mean mass of ice particles
    real(RP) :: rrhog(KA,IA,JA)                  ! 1/rhog
    real(RP) :: wtem(KA,IA,JA)                   ! temperature[K]
    real(RP) :: qd(KA,IA,JA)                     ! mixing ratio of dry air
    !
    real(RP) :: r_cva                    ! specific heat at constant volume
    real(RP) :: cpa(KA,IA,JA), r_cpa   ! specific heat at constant pressure
    real(RP) :: qsw(KA,IA,JA), r_qsw   ! saturated mixing ratio for liquid
    real(RP) :: dqswdtem_rho(KA,IA,JA) ! (dqsw/dtem)_rho
    real(RP) :: dqswdtem_pre(KA,IA,JA) ! (dqsw/dtem)_pre
    real(RP) :: dqswdpre_tem(KA,IA,JA) ! (dqsw/dpre)_tem
    !
    real(RP) :: w(KA,IA,JA)                     ! vetical velocity[m/s]
    real(RP) :: Acnd                              ! Pdynliq + Bergeron-Findeisen
    real(RP) :: aliqliq
    real(RP) :: Pdynliq                           ! production term of ssw by vertical motion
    real(RP) :: Pradliq                           ! production term of ssw by radiation
    real(RP) :: taucnd(KA,IA,JA),   r_taucnd    ! time scale of ssw change by MP
    real(RP) :: taucnd_c(KA,IA,JA), r_taucnd_c  ! by cloud
    real(RP) :: taucnd_r(KA,IA,JA), r_taucnd_r  ! by rain
    ! alternative tendency through changing ssw and ssi
    real(RP) :: PNCdep(KA,IA,JA) ! [Add] 11/08/30 T.Mitsui
    real(RP) :: PLR2NR
    real(RP) :: coef_a_cnd, coef_b_cnd
    !
    real(RP) :: dep_dqr, dep_dnr
    real(RP) :: dep_dqc, dep_dnc   ! 11/08/30 [Add] T.Mitsui, dep_dnc
    real(RP) :: r_xc_ccn           ! 11/08/30 [Add] T.Mitsui
    !
    real(RP) :: drhogqv
    real(RP) :: drhogqc, drhogqr
    real(RP) :: drhognc, drhognr
    !
    real(RP) :: fac1
    real(RP) :: r_rvaptem        ! 1/(Rvap*tem)
    real(RP) :: pv               ! vapor pressure
    real(RP) :: lvsw             ! saturated vapor density
    real(RP) :: dlvsw            !
    ! [Add] 11/08/30 T.Mitsui
    real(RP) :: dcnd             ! total cndensation/deposition
    real(RP) :: uplim_cnd        ! upper limit of condensation
    real(RP) :: lowlim_cnd       ! lower limit of evaporation
    real(RP) :: ssw              ! supersaturation ratio
    real(RP) :: r_esw            ! 1/esw
    real(RP) :: r_lvsw           ! 1/(lvsw*ssw)
    real(RP) :: evap_max         !
    real(RP) :: r_dt             ! 1/dt
    real(RP) :: ssw_o
    real(RP) :: dt_dyn
!    real(RP) :: dt_mp
    !
    real(RP), parameter :: tau100day   = 1.E+7_RP
    real(RP), parameter :: r_tau100day = 1.E-7_RP
    real(RP), parameter :: eps=1.E-30_RP
    !
    integer :: i,j,k
    !
    dt_dyn     = dt*ntmax
!    dt_mp      = dt*(ntdiv-1)
    !
    r_dt       = 1.0_RP/dt
    rrhog(:,:,:) = 1.0_RP/rhog(:,:,:)
    ! Temperature lower limit is only used for saturation condition.
    ! On the other hand original "tem" is used for calculation of latent heat or energy equation.
    wtem(:,:,:)  = max( tem(:,:,:), tem_min )
    !
    w(1:KS,:,:)  = 0.0_RP
    w(KE:KA,:,:) = 0.0_RP
    do k=KS,KE
       do j=1,JA
       do i=1,IA
          ! [Add] 11/08/30 T.Mitsui
          if( z(k) <= 25000.0_RP )then
             w(k,i,j) = 0.5_RP*(wh(k,i,j) + wh(k+1,i,j))
          else
             w(k,i,j) = 0.0_RP
          endif
       end do
       end do
    end do
    !
    call thrmdyn_qd(   &
         qd,           & !--- out
         q )             !--- in
    call thrmdyn_cv(   &
         cva,          & !--- out
         q,            & !--- in
         qd )            !--- in
    call thrmdyn_cp(   &
         cpa,          & !--- out
         q,            & !--- in
         qd )            !--- in
    call moist_qsat_water    ( qsw, wtem, pre )
    call moist_dqsw_dtem_rho ( dqswdtem_rho, wtem, rho )
    call moist_dqsw_dtem_dpre( dqswdtem_pre, dqswdpre_tem, wtem, pre )
    !
    do k=1, KA
       do j=1, JA
       do i=1, IA
          if( pre(k,i,j) < esw(k,i,j)+1.E-10_RP )then
             qsw(k,i,j) = 1.0_RP
             dqswdtem_rho(k,i,j) = 0.0_RP
             dqswdtem_pre(k,i,j) = 0.0_RP
             dqswdpre_tem(k,i,j) = 0.0_RP
          endif
       end do
       end do
    end do
    !
    ! taucnd, taudep
    taucnd_c(:,:,:)   = tau100day
    taucnd_r(:,:,:)   = tau100day
    taucnd(:,:,:)     = tau100day
    !
    r_xc_ccn=1.0_RP/xc_ccn
    do k=KS, KE
       do j=1, JA
          do i=1, IA
             r_rvaptem        = 1.0_RP/(Rvap*wtem(k,i,j))
             lvsw             = esw(k,i,j)*r_rvaptem        ! rho=p/(Rv*T)
             pv               = lv(k,i,j)*Rvap*tem(k,i,j)
             r_esw            = 1.0_RP/esw(k,i,j)
             ssw              = min( MP_ssw_lim, ( pv*r_esw-1.0_RP ) )
             r_lvsw           = 1.0_RP/lvsw
             r_taucnd_c       = PLCdep(k,i,j)*r_lvsw
             r_taucnd_r       = PLRdep(k,i,j)*r_lvsw
             taucnd_c(k,i,j)   = 1.0_RP/(r_taucnd_c+r_tau100day)
             taucnd_r(k,i,j)   = 1.0_RP/(r_taucnd_r+r_tau100day)
             !
             r_cva            = 1.0_RP/cva(k,i,j)
             r_cpa            = 1.0_RP/cpa(k,i,j)
             ! Coefficient of latent heat release for ssw change by PLCdep and PLRdep
             aliqliq          = 1.0_RP &
                  + r_cva*( LHV00              + (CVvap-CL)*tem(k,i,j) )*dqswdtem_rho(k,i,j)
             Pdynliq          = w(k,i,j)*GRAV * ( r_cpa*dqswdtem_pre(k,i,j) + rho(k,i,j)*dqswdpre_tem(k,i,j) )
             Pradliq          = -dTdt_rad(k,i,j)    * dqswdtem_rho(k,i,j)
             !
             r_qsw            = 1.0_RP/qsw(k,i,j)
             ssw_o            = ssw
             !
             Acnd             = Pdynliq + Pradliq 
             r_taucnd         = &
                  + aliqliq*( r_taucnd_c+r_taucnd_r ) 
             !
             uplim_cnd        = max( rho(k,i,j)*ssw_o*qsw(k,i,j)*r_dt, 0.0_RP )
             lowlim_cnd       = min( rho(k,i,j)*ssw_o*qsw(k,i,j)*r_dt, 0.0_RP )
             ! [fix] xxxxxx T.Seiki
!!$          if( r_taudep < r_tau100day )then
             if( r_taucnd < r_tau100day )then
                taucnd(k,i,j)     = tau100day
                PLCdep(k,i,j) = max(lowlim_cnd, min(uplim_cnd, PLCdep(k,i,j)*ssw_o ))
                PLRdep(k,i,j) = max(lowlim_cnd, min(uplim_cnd, PLRdep(k,i,j)*ssw_o ))
                PNRdep(k,i,j) = min(0.0_RP, PNRdep(k,i,j)*ssw_o )
                PLR2NR           = 0.0_RP
             else
                taucnd(k,i,j)     = 1.0_RP/r_taucnd
                ! Production term for liquid water content
                coef_a_cnd = rho(k,i,j)*Acnd*taucnd(k,i,j)
                coef_b_cnd = rho(k,i,j)*taucnd(k,i,j)*r_dt*(ssw_o*qsw(k,i,j)-Acnd*taucnd(k,i,j)) * ( exp(-dt*r_taucnd) - 1.0_RP )
                PLCdep(k,i,j) = coef_a_cnd*r_taucnd_c - coef_b_cnd*r_taucnd_c
                PLRdep(k,i,j) = coef_a_cnd*r_taucnd_r - coef_b_cnd*r_taucnd_r
                PLR2NR           = PNRdep(k,i,j)/(PLRdep(k,i,j)+1.d-30)
                PNRdep(k,i,j) = min(0.0_RP, PLRdep(k,i,j)*PLR2NR )
             endif
             !
             if( PLCdep(k,i,j) < -eps )then
                PNCdep(k,i,j) = min(0.0_RP, ((lc(k,i,j)+PLCdep(k,i,j)*dt)*r_xc_ccn - nc(k,i,j))*r_dt )
             else
                PNCdep=0.0_RP
             endif
             !
          end do
       end do
    end do
    !

    do k=KS, KE
       do j=1, JA
          do i=1, IA
             !
             !--- evaporation/condensation, deposition/sublimation
             !
             r_rvaptem = 1.0_RP/(Rvap*wtem(k,i,j))
             lvsw    = esw(k,i,j)*r_rvaptem
             dlvsw   = lv(k,i,j)-lvsw
             dcnd    = dt*(PLCdep(k,i,j)+PLRdep(k,i,j))
             dep_dqc = 0.0_RP
             dep_dqr = 0.0_RP
             !    always supersaturated
             if     ( (dcnd >  eps) .AND. (dlvsw > eps) )then
                fac1    = min(dlvsw,dcnd)/dcnd
                dep_dqc =  dt*PLCdep(k,i,j)*fac1
                dep_dqr =  dt*PLRdep(k,i,j)*fac1
                ! always unsaturated
             else if( (dcnd < -eps) .AND. (dlvsw < -eps) )then
                fac1    = max( dlvsw,dcnd )/dcnd
                dep_dqc = max( dt*PLCdep(k,i,j)*fac1, -lc(k,i,j) )
                dep_dqr = max( dt*PLRdep(k,i,j)*fac1, -lr(k,i,j) )
             else
                ! partially unsaturated during timestep
                fac1    = 1.0_RP
                dep_dqc = 0.0_RP
                dep_dqr = 0.0_RP
             endif
             !
             ! evaporation always lose number(always negative).
             dep_dnc = max( dt*PNCdep(k,i,j)*fac1, -nc(k,i,j) ) ! ss>0 dep=0, ss<0 dep<0 ! [Add] 11/08/30 T.Mitsui
             dep_dnr = max( dt*PNRdep(k,i,j)*fac1, -nr(k,i,j) ) ! ss>0 dep=0, ss<0 dep<0
             !
             !--- melting
             !
             ! water vapor change
             drhogqv = -(dep_dqc+dep_dqr)
             ! total cloud change
             drhogqc = ( dep_dqc )
             ! [Mod] 11/08/30 T.Mitsui
             drhognc = ( dep_dnc )
             ! total rain change
             drhogqr = ( dep_dqr )
             drhognr = ( dep_dnr )
             !
             !--- update
             ! filter for rounding error
             rhogq(k,i,j,I_QV) = max(0.0_RP, rhogq(k,i,j,I_QV) + drhogqv )
             !
             rhogq(k,i,j,I_QC) = max(0.0_RP, rhogq(k,i,j,I_QC) + drhogqc )
             rhogq(k,i,j,I_NC) = max(0.0_RP, rhogq(k,i,j,I_NC) + drhognc )
             rhogq(k,i,j,I_QR) = max(0.0_RP, rhogq(k,i,j,I_QR) + drhogqr )
             rhogq(k,i,j,I_NR) = max(0.0_RP, rhogq(k,i,j,I_NR) + drhognr )
             !
             rhoge(k,i,j) = rhoge(k,i,j) - LHV * drhogqv 
             !
             !--- update mixing ratio
             q(k,i,j,I_QV) = rhogq(k,i,j,I_QV) * rrhog(k,i,j)
             q(k,i,j,I_QC) = rhogq(k,i,j,I_QC) * rrhog(k,i,j)
             q(k,i,j,I_QR) = rhogq(k,i,j,I_QR) * rrhog(k,i,j)
             !
             q(k,i,j,I_NC) = rhogq(k,i,j,I_NC) * rrhog(k,i,j)
             q(k,i,j,I_NR) = rhogq(k,i,j,I_NR) * rrhog(k,i,j)
             !
!             if( i >= IS .AND. i <= IE .AND. j >= JS .AND. j <= JE ) then
             sl_PLCdep(1,i,j) = sl_PLCdep(1,i,j) + dep_dqc*Dz(k)*gsgam2(k,i,j)
             sl_PLRdep(1,i,j) = sl_PLRdep(1,i,j) + dep_dqr*Dz(k)*gsgam2(k,i,j)
             sl_PNRdep(1,i,j) = sl_PNRdep(1,i,j) + dep_dnr*Dz(k)*gsgam2(k,i,j)
!             endif
          end do
       end do
    end do
    call thrmdyn_qd(   &
         qd,           & !--- out
         q )             !--- in
    call thrmdyn_cv(   &
         cva,          & !--- out
         q,            & !--- in
         qd )            !--- in
    tem(:,:,:) = rhoge(:,:,:) / ( rhog(:,:,:) * cva(:,:,:) )
    pre(:,:,:) = rho(:,:,:)*( qd(:,:,:)*Rdry+q(:,:,:,I_QV)*Rvap )*tem(:,:,:)
    !
    return
    
  end subroutine update_by_phase_change_kij
  !-------------------------------------------------------------------------------
  subroutine MP_negativefilter
     use mod_atmos_vars, only: &
       DENS, &
       QTRC
    implicit none

    real(RP) :: diffq(KA,IA,JA)
    real(RP) :: r_xmin

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    ! total hydrometeor (before correction)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       diffq(k,i,j) = QTRC(k,i,j,I_QV) &
                    + QTRC(k,i,j,I_QC) &
                    + QTRC(k,i,j,I_QR) 
    enddo
    enddo
    enddo

    ! remove negative value of hydrometeor (mass, number)
    do iq = 1, QA
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( QTRC(k,i,j,iq) < 0.0_RP ) then
          QTRC(k,i,j,iq) = 0.0_RP
       endif
    enddo
    enddo
    enddo
    enddo

    ! apply correction of hydrometeor to total density
    ! [note] mass conservation is broken here to fill rounding error.
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       DENS(k,i,j) = DENS(k,i,j)        &
                   * ( 1.0_RP           &
                     + QTRC(k,i,j,I_QV) &
                     + QTRC(k,i,j,I_QC) &
                     + QTRC(k,i,j,I_QR) &
                     - diffq(k,i,j)      ) ! after-before
    enddo
    enddo
    enddo

    ! avoid unrealistical value of number concentration
    ! due to numerical diffusion in advection
    r_xmin = 1.0_RP / xmin_filter
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( QTRC(k,i,j,I_NC) > QTRC(k,i,j,I_QC)*r_xmin ) then
          QTRC(k,i,j,I_NC) = QTRC(k,i,j,I_QC)*r_xmin
       endif
    enddo
    enddo
    enddo
!
    do j  = 1, JA
    do i  = 1, IA
    do k  = 1, KA
       if ( QTRC(k,i,j,I_NR) > QTRC(k,i,j,I_QR)*r_xmin ) then
          QTRC(k,i,j,I_NR) = QTRC(k,i,j,I_QR)*r_xmin
       endif
    enddo
    enddo
    enddo

    return
  end subroutine MP_negativefilter
  !-------------------------------------------------------------------------------
end module mod_atmos_phy_mp
!-------------------------------------------------------------------------------
